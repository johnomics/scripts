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

$OUTPUT_AUTOFLUSH = 1;

my %mstheader = (
    population_type              => "DH",
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

croak "No VCF file! Please specify with -v"      if $args{vcf_file} eq "";
croak "No blocks file! Please specify with -b"   if $args{blocks} eq "";
croak "No genetics file! Please specify with -g" if $args{genetics} eq "";

my $genetics = load_genetics( $args{genetics} );

my ( $scfl, $samples ) = load_vcf_header( $args{vcf_file}, $genetics );

my $markers = load_markers( $args{blocks} );

my $marker_lookup = run_mstmap( $markers, $samples, $args{output} );

write_agp( $markers, $marker_lookup, $args{output}, $scfl );

sub load_markers {
    my ($blocks) = @_;
    my %markers;

    open my $blocksfile, '<', $blocks
      or croak "Can't open $blocks: $OS_ERROR\n";
    while ( my $blocksline = <$blocksfile> ) {
        next if ( $blocksline =~ /^Scaffold/ );
        my @f = split /\t/, uncolor $blocksline;
        my $paternal = $f[12];
        next if ( $paternal =~ "[ ~\.]" );
        $markers{$paternal}{ $f[0] }{ $f[1] } = $f[2];
    }
    close $blocksfile;
    return \%markers;
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

    $mstheader{"number_of_loci"}       = keys %{$markers};
    $mstheader{"number_of_individual"} = @{$samples};
    map { print $mstmapin "$_ $mstheader{$_}\n"; } @mstheader;

    print $mstmapin "locus_name";
    map { print $mstmapin "\t$_" } @{$samples};
    print $mstmapin "\n";

    my $id = 1;
    while ( my ($marker) = each %{$markers} ) {
        print $mstmapin "$id";
        my @gt = split //, $marker;
        map {
            print $mstmapin "\t";
            print $mstmapin $_ eq 'H' ? 'B' : 'A';
        } @gt;
        print $mstmapin "\n";
        $marker_lookup{$id} = $marker;
        $id++;
    }

    close $mstmapin;

    system(
"MSTMap.exe $output.mstmap.markers $output.mstmap.map > $output.mstmap.log"
    );

    return \%marker_lookup;
}

sub write_agp {
    my ( $markers, $marker_lookup, $output, $scfl ) = @_;

    open my $mstout, '<', "$output.mstmap.map"
      or croak "Can't open MST output! $OS_ERROR\n";

    my $chr     = "";
    my $ingroup = 0;
    open my $chragp, '>', "$output.chrom.agp" or croak "Can't open AGP file! $OS_ERROR\n";
    
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
    my $pos = 1;
    my $part = 1;
    my @agplines;
    foreach my $cm ( sort { $a <=> $b } keys %{$lg} ) {
        foreach my $marker ( keys %{ $lg->{$cm} } ) {
            my $pattern = $marker_lookup->{$marker};
            foreach my $scaffold ( sort keys %{ $markers->{$pattern} } ) {
                my $chrend = $pos + $scfl->{$scaffold} - 1;
                my $gapstart = $chrend+1;
                my $gapend = $gapstart + 99;
                push @agplines, "$chr\t$pos\t$chrend\t$part\tD\t$scaffold\t1\t$scfl->{$scaffold}\t+\t# Contigs | $cm cM\n";
                $part++;
                push @agplines, "$chr\t$gapstart\t$gapend\t$part\tN\t100\tfragment\tno\t\n";
                $part++;
                $pos += $scfl->{$scaffold}+100;
            }
        }
    }
    pop @agplines; # Lose final gap
    print $chragp @agplines;
}
