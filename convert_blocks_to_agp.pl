#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

use Parallel::ForkManager;
use Term::ExtendedColor qw/:all/;
use List::Util qw/min max/;

$OUTPUT_AUTOFLUSH = 1;

my %mstheader = (
    population_type              => "RIL2",
    population_name              => "HeliconiusWGS",
    distance_function            => "kosambi",
    cut_off_p_value              => "0.000001",
    no_map_dist                  => "15",
    no_map_size                  => "2",
    missing_threshold            => "0.25",
    estimation_before_clustering => "yes",
    detect_bad_data              => "yes",
    objective_function           => "ML",
);

my %args;

$args{vcf_file} = "";
$args{blocks}   = "";
$args{genetics} = "";
$args{output}   = "test";
my $options_okay = GetOptions(
    'vcffile=s'  => \$args{vcf_file},
    'blocks=s'   => \$args{blocks},
    'genetics=s' => \$args{genetics},
    'output=s'   => \$args{output},

);

my %hom = (
    "A" => 0,
    "H" => 0,
);

croak "No VCF file! Please specify with -v"      if $args{vcf_file} eq "";
croak "No blocks file! Please specify with -b"   if $args{blocks} eq "";
croak "No genetics file! Please specify with -g" if $args{genetics} eq "";

my $genetics = load_genetics( $args{genetics} );

my ( $scfl, $samples ) = load_vcf_header( $args{vcf_file}, $genetics );

my $markers = load_markers( $args{blocks} );

my $marker_lookup = run_mstmap( $markers, $samples, $args{output} );

my $map = integrity_check( $markers, $marker_lookup, $args{output}, $scfl );

#write_agp( $markers, $marker_lookup, $args{output}, $scfl );

sub load_markers {
    my ($blocks) = @_;
    my %markers;
    my %header;

    open my $blocksfile, '<', $blocks
      or croak "Can't open $blocks: $OS_ERROR\n";
    while ( my $blocksline = <$blocksfile> ) {
        my @f = split /\t/, uncolor $blocksline;

        if ( $blocksline =~ /^Scaffold/ ) {
            next if keys %header ne 0;
            for my $i ( 0 .. $#f ) {
                $f[$i] =~ s/(\s+)//;
                $header{ $f[$i] } = $i;
            }
            next;
        }

        my $maternal = transform( $f[ $header{"Maternal-AHAH"} ] );
        $markers{maternal}{$maternal}{ $f[0] }{ $f[1] } = $f[2] if $maternal ne "";

        for my $pat ($f[$header{"Paternal-AHAH"}], $f[$header{"Intercross-ABHABH_HHA"}]) {
            next if ($pat =~ "[~\.\ ]");
            $pat = check_mirror($pat);
            $markers{intercross}{$pat}{$f[0]}{$f[1]} = $f[2];
        }
    }
    close $blocksfile;

    my @mats = sort { keys %{ $markers{maternal}{$b} } <=> keys %{ $markers{maternal}{$a} } }
        keys %{$markers{maternal}};
    for my $mata (@mats)
    {
        next if !defined $markers{maternal}{$mata};
        my $counta = keys %{ $markers{maternal}{$mata} };
        for my $matb (@mats) {
            next if !defined $markers{maternal}{$matb};
            next if ($mata eq $matb);
            my $countb = keys %{ $markers{maternal}{$matb} };
            my $hamming = get_dist($mata, $matb);
            if ($hamming <= 4) {
                for my $scf (keys %{$markers{maternal}{$matb}}) {
                    for my $start (keys %{$markers{maternal}{$matb}{$scf}}) {
                        $markers{maternal}{$mata}{$scf}{$start} = $markers{maternal}{$matb}{$scf}{$start};
                    }
                }
                delete $markers{maternal}{$matb};
            }
        }
    }

    for my $p (sort keys %{$markers{maternal}}) {
        my $scf = keys %{$markers{maternal}{$p}};
        print "$p\t$scf\n";
    }
    my %lengths;
    for my $pattern ( keys %{$markers{intercross}} ) {
        my $length = 0;
        for my $scf ( keys %{ $markers{intercross}{$pattern} } ) {
            for my $start ( keys %{ $markers{intercross}{$pattern}{$scf} } ) {
                $length += $markers{intercross}{$pattern}{$scf}{$start} - $start + 1;
            }
        }
        delete $markers{intercross}{$pattern} if $length == 1;
    }

    return \%markers;
}

sub get_dist {
    my ($mata, $matb) = @_;
    my @matas = split //, $mata;
    my @matbs = split //, $matb;
    my $hamming = 0;
    for my $i (0..$#matas) {
        next if $matas[$i] eq ' ' or $matbs[$i] eq ' ';
        $hamming++ if $matas[$i] ne $matbs[$i];
    }
    return $hamming;
}

sub check_mirror {
    my ($pattern, $markers) = @_;
    my $mirror = $pattern;
    $mirror =~ tr/AB/BA/;
    return defined $markers->{$mirror} ? $mirror : $pattern;
}

sub transform {
    my ($pattern) = @_;

    return "" if $pattern =~ /^( +)$/;
    my @s = split //, $pattern;
    my $trans = "";
    for my $i ( 0 .. $#s - 1 ) {
        if ( defined $hom{ $s[$i] } && defined $hom{ $s[ $i + 1 ] } ) {
            $trans .= $s[$i] eq $s[ $i + 1 ] ? '-' : 'x';
        }
        else {
            $trans .= ' ';
        }
    }
    return $trans;
}

sub load_genetics {
    my ($geneticsfilename) = @_;

    my %genetics;

    open my $geneticsfile, "<", $geneticsfilename
      or croak "Can't open marker type file $geneticsfilename: $OS_ERROR\n";

    my $infoline;
    while ( $infoline = <$geneticsfile> ) {
        chomp $infoline;
        last if ( $infoline =~ /Type/ );

        if ( $infoline =~ /^Female/ or $infoline =~ /^Male/ ) {
            my ( $sex, $samples ) = split /\t/, $infoline;
            map { $genetics{samples}{$_}++ } split /,/, $samples;
        }
    }
    return \%genetics;
}

sub load_vcf_header {
    my ( $vcffilename, $genetics ) = @_;

    open my $vcffile, '<', $vcffilename
      or croak "Can't open $vcffilename: $OS_ERROR!\n";

    my %scfl;
    my @samples;
    while ( my $vcfline = <$vcffile> ) {

        if ( $vcfline =~ /^##contig=<ID=(.+),length=(.+)>$/ ) {
            $scfl{$1} = $2;
        }
        if ( $vcfline =~ /^#CHROM/ ) {
            chomp $vcfline;
            my @header = split /\t/, $vcfline;
            map { push @samples, $_ if defined $genetics->{samples}{$_} }
              @header;
            last;
        }
    }
    close $vcffile;

    return ( \%scfl, \@samples );
}

sub run_mstmap {

    my ( $markers, $samples, $output ) = @_;

    my @mstheader = (
        "population_type",   "population_name",
        "distance_function", "cut_off_p_value",
        "no_map_dist",       "no_map_size",
        "missing_threshold", "estimation_before_clustering",
        "detect_bad_data",   "objective_function",
        "number_of_loci",    "number_of_individual",
    );

    my %marker_lookup;

    open my $mstmapin, ">", "$output.mstmap.markers"
      or croak "Can't open $output.mstmap.markers: $OS_ERROR\n";

    open my $codein, ">", "$output.mstmap.marker.code"
      or croak "Can't open $output.mstmap.marker.code: $OS_ERROR\n";

    $mstheader{"number_of_loci"}       = keys %{$markers->{intercross}};
    $mstheader{"number_of_individual"} = @{$samples};
    map { print $mstmapin "$_ $mstheader{$_}\n"; } @mstheader;

    print $mstmapin "locus_name";
    map { print $mstmapin "\t$_" } ( 1 .. @{$samples} );

    #    map { print $mstmapin "\t$_" } @{$samples};

    print $mstmapin "\n";

    print $codein "ID\tPattern\n";
    my $id = 1;
    while ( my ($marker) = each %{$markers->{intercross}} ) {
        print $mstmapin "$id";
        print $codein "$id\t$marker\n";
        my @gt = split //, $marker;
        map {
            print $mstmapin "\t";
            print $mstmapin $_ eq 'H' ? 'X' : $_;
        } @gt;
        print $mstmapin "\n";
        $marker_lookup{$id} = $marker;
        $id++;
    }

    close $codein;
    close $mstmapin;

    system(
"MSTMap.exe $output.mstmap.markers $output.mstmap.map > $output.mstmap.log"
    );

    return \%marker_lookup;
}

sub integrity_check {
    my ( $markers, $marker_lookup, $output, $scfl ) = @_;

    my $lg = load_map($output);

    my $scflg = match_scaffolds_to_map( $markers, $marker_lookup, $lg );

    for my $scf ( sort keys %{$scflg} ) {
        for my $lg ( sort keys %{ $scflg->{$scf} } ) {
            my $blockstart;
            my $blockend;
            my $blocklen;
            my $blockcm = 0;
            for my $start ( sort { $a <=> $b } keys %{ $scflg->{$scf}{$lg} } ) {
                if ( $blockcm ne $scflg->{$scf}{$lg}{$start}{cm} ) {
                    print
                      "$scf\t$lg\t$blockcm\t$blockstart\t$blockend\t$blocklen\n"
                      if defined $blockstart;
                    $blockcm    = $scflg->{$scf}{$lg}{$start}{cm};
                    $blockstart = $start;
                    $blockend   = $scflg->{$scf}{$lg}{$start}{end};
                    $blocklen   = $scflg->{$scf}{$lg}{$start}{len};
                }
                else {
                    $blockend = $scflg->{$scf}{$lg}{$start}{end};
                    $blocklen += $scflg->{$scf}{$lg}{$start}{len};
                }
            }
            print "$scf\t$lg\t$blockcm\t$blockstart\t$blockend\t$blocklen\n";
        }
    }
}

sub match_scaffolds_to_map {
    my ( $markers, $marker_lookup, $lg ) = @_;

    my %scf;

    for my $lgnum ( keys %{$lg} ) {
        for my $cm ( keys %{ $lg->{$lgnum} } ) {
            for my $marker ( keys %{ $lg->{$lgnum}{$cm} } ) {
                my $pattern = $marker_lookup->{$marker};
                for my $scaffold ( keys %{ $markers->{intercross}{$pattern} } ) {
                    for my $start ( keys %{ $markers->{intercross}{$pattern}{$scaffold} } )
                    {
                        $scf{$scaffold}{$lgnum}{$start}{cm} = $cm;
                        $scf{$scaffold}{$lgnum}{$start}{end} =
                          $markers->{intercross}{$pattern}{$scaffold}{$start};
                        $scf{$scaffold}{$lgnum}{$start}{len} =
                          $markers->{intercross}{$pattern}{$scaffold}{$start} - $start + 1;
                    }
                }
            }
        }
    }

    return \%scf;
}

sub load_map {
    my ($output) = @_;
    open my $mstout, '<', "$output.mstmap.map"
      or croak "Can't open MST output! $OS_ERROR\n";
    my %lg;
    my $lgnum;
    my $ingroup = 0;
    while ( my $mstline = <$mstout> ) {
        chomp $mstline;

        if ( $mstline =~ /size of the linkage groups/ ) {
            my $lgsizeline = <$mstout>;
            my @lgsizes = split /\s+/, $lgsizeline;
            shift @lgsizes;    # Semicolon at start
            my %lgsize;
            my $lgcount = 0;
            map { $lgsize{$_}++; $lgcount++ } @lgsizes;
            print "LG size\tLGs\n";
            for my $lgmarknum ( sort { $a <=> $b } keys %lgsize ) {
                print "$lgmarknum\t$lgsize{$lgmarknum}\n";
            }
            print "Found $lgcount linkage groups\n";
        }

        $lgnum = $1 if ( $mstline =~ /^group (.+)$/ );
        if ( $mstline eq ";BEGINOFGROUP" ) {
            $ingroup = 1;
            next;
        }
        $ingroup = 0 if ( $mstline eq ";ENDOFGROUP" );
        if ($ingroup) {
            if ( $mstline =~ /^(.+)\t([\d\.]+)$/ ) {
                my $marker = $1;
                my $cm     = $2;
                $lg{$lgnum}{$cm}{$marker}++;
            }
        }
    }
    close $mstout;
    return \%lg;
}

sub write_agp {
    my ( $markers, $marker_lookup, $output, $scfl ) = @_;

    open my $mstout, '<', "$output.mstmap.map"
      or croak "Can't open MST output! $OS_ERROR\n";

    my $chr     = "";
    my $ingroup = 0;
    open my $chragp, '>', "$output.chrom.agp"
      or croak "Can't open AGP file! $OS_ERROR\n";

    my %lg;
    while ( my $mstline = <$mstout> ) {
        chomp $mstline;
        $chr = $1 if ( $mstline =~ /^group (.+)$/ );
        if ( $mstline eq ";BEGINOFGROUP" ) {
            $ingroup = 1;
            next;
        }
        if ( $mstline eq ";ENDOFGROUP" ) {
            $ingroup = 0;
            output_lg( $chr, \%lg, $markers, $marker_lookup, $scfl, $chragp );
            %lg = ();
        }
        if ($ingroup) {
            if ( $mstline =~ /^(.+)\t([\d\.]+)$/ ) {
                my $marker = $1;
                my $cm     = $2;
                $lg{$cm}{$marker}++;
            }
        }
    }
    close $chragp;
    close $mstout;
}

sub output_lg {
    my ( $chr, $lg, $markers, $marker_lookup, $scfl, $chragp ) = @_;
    my $pos  = 1;
    my $part = 1;
    my @agplines;
    foreach my $cm ( sort { $a <=> $b } keys %{$lg} ) {
        foreach my $marker ( keys %{ $lg->{$cm} } ) {
            my $pattern = $marker_lookup->{$marker};
            foreach my $scaffold ( sort keys %{ $markers->{intercross}{$pattern} } ) {
                my $chrend   = $pos + $scfl->{$scaffold} - 1;
                my $gapstart = $chrend + 1;
                my $gapend   = $gapstart + 99;
                push @agplines,
"$chr\t$pos\t$chrend\t$part\tD\t$scaffold\t1\t$scfl->{$scaffold}\t+\t# Contigs | $cm cM\n";
                $part++;
                push @agplines,
                  "$chr\t$gapstart\t$gapend\t$part\tN\t100\tfragment\tno\t\n";
                $part++;
                $pos += $scfl->{$scaffold} + 100;
            }
        }
    }
    pop @agplines;    # Lose final gap
    print $chragp @agplines;
}
