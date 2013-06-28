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

my $samplenum = 69;
my $empty     = ' ' x $samplenum;

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
$args{blocks} = "";
$args{output} = "test";

my $options_okay =
  GetOptions( 'blocks=s' => \$args{blocks}, 'output=s' => \$args{output} );

my ( $header, $types, $blocklist ) = load_blocks( $args{blocks} );

print "After LOAD BLOCKS\n";
for my $block (@{$blocklist}) {
    print Dumper $block if ($block->{'Start'} == 53270);
}

fill_blocks( $header, $types, $blocklist );
print "After FILL BLOCKS\n";

for my $block (@{$blocklist}) {
    print Dumper $block if ($block->{'Start'} == 53270);
}

correct_maternal( "Maternal-AHAH", $blocklist );
print "After CORRECT MATERNAL\n";

for my $block (@{$blocklist}) {
    print Dumper $block if ($block->{'Start'} == 53270);
}

collapse( $blocklist, $types );
print "After COLLAPSE\n";

for my $block (@{$blocklist}) {
    print Dumper $block if ($block->{'Start'} == 53270);
}

make_maps( $blocklist, $args{output} );

my %presents;
for my $block ( @{$blocklist} ) {
    map { print STDERR "$block->{$_}\t" } @{$header};
    print STDERR "\n";
}

sub make_maps {
    my ( $blocklist, $output ) = @_;

    my %matpat;
    my %patmat;
    my $glength = 0;
    my %pattern_block;
    for my $block ( @{$blocklist} ) {
        my $mat = $block->{"Maternal-AHAH"};
        my $pat = $block->{"Paternal-AHAH"};
        next if $mat =~ /[ \-]/;
        if ( $pat =~ /^( +)$/ ) {
            $pat = convert_intercross($block);
        }
        next if $pat =~ /[ \-]/;
        $matpat{$mat}{$pat}{length} += $block->{'Length'};
        $matpat{$mat}{$pat}{blocks}++;
        $patmat{$pat}{$mat}{length} += $block->{'Length'};
        $patmat{$pat}{$mat}{blocks}++;
        $pattern_block{$pat}{ $block->{'Scaffold'} }{ $block->{'Start'} } =
          $block->{'End'};
    }
    exit;

    for my $pat (sort keys %patmat) {
        if (keys %{$patmat{$pat}} > 1) {
            for my $mat (sort keys %{$patmat{$pat}}) {
                print "$pat\t$mat\t$patmat{$pat}{$mat}{blocks}\t$patmat{$pat}{$mat}{length}\n";
            }
            for my $scf (keys %{$pattern_block{$pat}}) {
                for my $start (keys %{$pattern_block{$pat}{$scf}}) {
                    my $length = $pattern_block{$pat}{$scf}{$start} - $start + 1;
                    print "\t$scf\t$start\t$pattern_block{$pat}{$scf}{$start}\t$length\n";
                }
            }
        }
    }
    my %scfmap;
    my %genome;
    for my $mat ( keys %matpat ) {
        for my $pat ( keys %{ $matpat{$mat} } ) {
            delete $matpat{$mat}{$pat}
              if (  $matpat{$mat}{$pat}{length} < 1000
                and $matpat{$mat}{$pat}{blocks} <= 2 );
        }
        next if ( keys %{ $matpat{$mat} } == 0 );
        my $markercode = run_mstmap( $mat, $matpat{$mat}, $output );
        $genome{$mat} = load_map( $output, $mat );
        for my $lg ( sort keys %{ $genome{$mat} } ) {
            for my $cm ( sort { $a <=> $b } keys %{ $genome{$mat}{$lg} } ) {
                for my $marker (
                    sort { $a <=> $b }
                    keys %{ $genome{$mat}{$lg}{$cm} }
                  )
                {
                    my $pattern = $markercode->{$marker};
                    for my $scf ( sort keys %{ $pattern_block{$pattern} } ) {
                        for my $start (
                            sort { $a <=> $b }
                            keys %{ $pattern_block{$pattern}{$scf} }
                          )
                        {
                            my $length =
                              $pattern_block{$pattern}{$scf}{$start} -
                              $start + 1;
                            $scfmap{$scf}{$start}{mat} = $mat;
                            $scfmap{$scf}{$start}{lg}  = $lg;
                            $scfmap{$scf}{$start}{cm}  = $cm;
                            if ($scf eq 'HE670564') {
                                print "$mat\t$lg\t$cm\t$marker\t$pattern\t$scf\t$start\n";
                                print Dumper $scfmap{$scf};
                            } 
                        }
                    }
                }
            }
        }
    }
exit;
    my $curscf = "";
    my @scfblocks;
    my %scfstats;
    for my $block ( @{$blocklist} ) {
        if ( $block->{'Scaffold'} ne $curscf ) {
            if ( $curscf ne "" ) {
                validate_scaffold( \@scfblocks, \%scfmap, \%genome,
                    \%scfstats );
            }
            $curscf    = $block->{'Scaffold'};
            @scfblocks = ();
        }
        push @scfblocks, $block;
    }
    validate_scaffold( \@scfblocks, \%scfmap, \%genome, \%scfstats );

    my %genomestat;
    foreach my $scf ( sort keys %scfstats ) {

        my $stat = $scfstats{$scf};

        print
"$scf\t$stat->{markerblocks}\t$stat->{allblocks}\t$stat->{length}\t$stat->{chromosomes}\t$stat->{lgs}\t$stat->{gaps}\n";
        if ( $stat->{chromosomes} == 0 ) {
            if ( $stat->{lgs} == 0 and $stat->{gaps} == 0 ) {
                $genomestat{"Unassigned"}{scf}++;
                $genomestat{"Unassigned"}{len} += $stat->{length};
            }
        }
        elsif ( $stat->{chromosomes} == 1 ) {
            if ( $stat->{lgs} == 1 ) {
                if ( $stat->{gaps} == 0 ) {
                    $genomestat{"Assigned"}{scf}++;
                    $genomestat{"Assigned"}{len} += $stat->{length};
                }
                else {
                    $genomestat{"Gaps"}{scf}++;
                    $genomestat{"Gaps"}{len} += $stat->{length};
                }
            }
            else {
                $genomestat{"Multiple LGs"}{scf}++;
                $genomestat{"Multiple LGs"}{len} += $stat->{length};
            }
        }
        else {
            $genomestat{"Multiple Chrs"}{scf}++;
            $genomestat{"Multiple Chrs"}{len} += $stat->{length};
        }
    }
    my $genomesize = 0;
    my $genomescf  = 0;
    for my $stat ( sort keys %genomestat ) {
        printf STDERR "%16s\t%4d\t%9d\n", $stat, $genomestat{$stat}{scf},
          $genomestat{$stat}{len};
        $genomesize += $genomestat{$stat}{len};
        $genomescf  += $genomestat{$stat}{scf};
    }
    printf STDERR "%16s\t%4d\t%9d\n", 'Genome', $genomescf, $genomesize;
}

sub validate_scaffold {
    my ( $scfblocks, $scfmap, $genome, $stats ) = @_;
    my $blocks_with_marker = 0;
    my $scf                = $scfblocks->[0]{'Scaffold'};
    my %scfcms;
    my $scflen = 0;
    for my $block ( @{$scfblocks} ) {
        $scflen += $block->{'Length'};
        if ( defined $scfmap->{ $block->{'Scaffold'} }{ $block->{'Start'} } ) {
            my $mappos = $scfmap->{$scf}{ $block->{'Start'} };
            $blocks_with_marker++;
            $scfcms{ $mappos->{mat} }{ $mappos->{lg} }{ $mappos->{cm} }
              {blocks}++;
            $scfcms{ $mappos->{mat} }{ $mappos->{lg} }{ $mappos->{cm} }{length}
              += $block->{'Length'};
        }
    }
    my $scfchroms = keys %scfcms // 0;
    my $scflgs    = 0;
    my $scfgaps   = 0;
    for my $mat ( sort keys %scfcms ) {
        for my $lg ( sort keys %{ $scfcms{$mat} } ) {
            $scflgs++;
            my @scfcm       = sort { $a <=> $b } keys %{ $scfcms{$mat}{$lg} };
            my $start_check = 0;
            my $gap_cms     = 0;
            for my $lgcm ( sort { $a <=> $b } keys %{ $genome->{$mat}{$lg} } ) {
                last if @scfcm == 0;
                if ($start_check) {
                    if ( $scfcm[0] ne $lgcm ) {
                        $gap_cms++;
                    }
                }
                if ( $scfcm[0] eq $lgcm ) {
                    print
"$scf\t$mat\t$lg\t$lgcm\t$scfcms{$mat}{$lg}{$lgcm}{length}\t$scfcms{$mat}{$lg}{$lgcm}{blocks}\n";
                    $start_check = 1;
                    shift @scfcm;
                }
            }
            $scfgaps++ if $gap_cms > 0;
        }
    }
    $stats->{$scf}{markerblocks} = $blocks_with_marker;
    $stats->{$scf}{allblocks}    = scalar @{$scfblocks};
    $stats->{$scf}{length}       = $scflen;
    $stats->{$scf}{chromosomes}  = $scfchroms;
    $stats->{$scf}{lgs}          = $scflgs;
    $stats->{$scf}{gaps}         = $scfgaps;

}

sub convert_intercross {
    my ($block) = @_;

    for my $ic ( "Intercross-ABHABH_HHA", "Intercross-ABHABH_HHH" ) {
        next if $block->{$ic} =~ /[ \-]/;
        my $mat = $block->{'Maternal-AHAH'};
        my @int = split //, phase( $block->{$ic} );
        my @mat = split //, $mat;
        my @pat =
          map { $int[$_] eq 'H' ? ( $mat[$_] eq 'A' ? 'B' : 'A' ) : $mat[$_] }
          0 .. $#int;
        my $pat = join '', @pat;
        return $pat;
    }
    return " ";
}

sub phase {
    my ($pat) = @_;
    my @p = split //, $pat;
    my $out = "";
    my %trans;
    my @checkbases;
    if ( $pat =~ 'A' && $pat =~ 'B' && $pat =~ 'H' ) {
        @checkbases = ( 'A', 'B' );
    }
    elsif ( $pat =~ 'A' && $pat =~ 'H' && $pat !~ 'B' ) {
        @checkbases = ( 'A', 'H' );
    }
    elsif ( $pat =~ 'B' && $pat =~ 'H' && $pat !~ 'A' ) {
        @checkbases = ( 'B', 'H' );
    }
    return " " x length($pat) if !@checkbases;
    for my $gt (@p) {
        my $added = 0;
        for my $bi (0,1) {
            if ( $gt eq $checkbases[$bi] ) {
                if ( !defined $trans{$gt} ) {
                    $trans{$checkbases[$bi]} = 'A';
                    my $other = $bi == 0 ? 1 : 0;
                    $trans{$checkbases[$other]} = 'B';
                }
                $out .= $trans{$gt};
                $added++;
            }
        }
        $out .= $gt =~ /[ ~\.]/ ? '-' : $gt if !$added;
    }
    return $out;
}

sub load_map {
    my ( $output, $pattern ) = @_;
    open my $mstout, '<', "$output.$pattern.mstmap.map"
      or croak "Can't open map for $pattern! $OS_ERROR\n";
    my %lg;
    my $lgnum;
    my $ingroup = 0;
    while ( my $mstline = <$mstout> ) {
        chomp $mstline;

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

sub make_transpat {
    my ($pat) = @_;
    my @gt = split //, $pat;
    my $transpat = "";
    for my $i ( 0 .. $#gt - 1 ) {
        $transpat .=
            ( $gt[$i] !~ /[ABH]/ or $gt[ $i + 1 ] !~ /[ABH]/ ) ? ' '
          : $gt[$i] eq $gt[ $i + 1 ] ? '-'
          :                            'x';

    }
    return $transpat;
}

sub run_mstmap {

    my ( $pattern, $markers, $output ) = @_;

    my @mstheader = (
        "population_type",   "population_name",
        "distance_function", "cut_off_p_value",
        "no_map_dist",       "no_map_size",
        "missing_threshold", "estimation_before_clustering",
        "detect_bad_data",   "objective_function",
        "number_of_loci",    "number_of_individual",
    );

    my %marker_lookup;

    open my $mstmapin, ">", "$output.$pattern.mstmap.markers"
      or croak "Can't open $output.$pattern.mstmap.markers: $OS_ERROR\n";

    open my $codein, ">", "$output.$pattern.mstmap.marker.code"
      or croak "Can't open $output.$pattern.mstmap.marker.code: $OS_ERROR\n";

    my $samples = split //, ( keys %{$markers} )[0];
    $mstheader{"number_of_loci"}       = keys %{$markers};
    $mstheader{"number_of_individual"} = $samples;
    map { print $mstmapin "$_ $mstheader{$_}\n"; } @mstheader;

    print $mstmapin "locus_name";
    map { print $mstmapin "\t$_" } ( 1 .. $samples );
    print $mstmapin "\n";

    print $codein "ID\tPattern\n";
    my $id = 1;
    while ( my ($marker) = each %{$markers} ) {
        print $mstmapin "$id";
        print $codein "$id\t$marker\n";
        my @gt = split //, $marker;
        map { print $mstmapin "\t$_"; } @gt;
        print $mstmapin "\n";
        $marker_lookup{$id} = $marker;
        $id++;
    }

    close $codein;
    close $mstmapin;

    system(
"MSTMap.exe $output.$pattern.mstmap.markers $output.$pattern.mstmap.map > $output.$pattern.mstmap.log"
    );

    return \%marker_lookup;
}

sub collapse {
    my ( $blocklist, $types ) = @_;
    my $i = 0;
    my $j = 1;
    while ( $j <= $#{$blocklist} ) {
        my $same = 1;
        $same = 0
          if ( $blocklist->[$i]{'Scaffold'} ne $blocklist->[$j]{'Scaffold'} );
        map {
            if ( $blocklist->[$i]{$_} ne $blocklist->[$j]{$_} ) {
                $same = 0;
            }
        } keys %{$types};

        if ($same) {
            $blocklist->[$i]{'End'}    = $blocklist->[$j]{'End'};
            $blocklist->[$i]{'Length'} = sprintf "%8s",
              $blocklist->[$i]{'Length'} + $blocklist->[$j]{'Length'};
            splice( @{$blocklist}, $j, 1 );
        }
        else {
            $i = $j;
            $j++;
        }
    }
}

sub fill_blocks {
    my ( $header, $types, $blocklist ) = @_;
    fill_type_pair( "Intercross-ABHABH_HHA", "Intercross-ABHABH_HHH",
        $blocklist );
    print "After Intercross-ABHABH_HHA Intercross-ABHABH_HHH\n";

        for my $block (@{$blocklist}) {
            print Dumper $block if ($block->{'Start'} == 53270);
        }
    
    fill_type_pair( "Intercross-ABHABH_HHH", "Intercross-ABHABH_HHA",
        $blocklist );
        print "After Intercross-ABHABH_HHH Intercross-ABHABH_HHA\n";

            for my $block (@{$blocklist}) {
                print Dumper $block if ($block->{'Start'} == 53270);
            }

    fill_type_pair( "Intercross-ABHABH_HHA", "Paternal-AHAH", $blocklist );

    print "After Intercross-ABHABH_HHA Paternal-AHAH\n";

        for my $block (@{$blocklist}) {
            print Dumper $block if ($block->{'Start'} == 53270);
        }

    fill_type_pair( "Paternal-AHAH", "Intercross-ABHABH_HHA", $blocklist );

    print "After Paternal-AHAH Intercross-ABHABH_HHA\n";

        for my $block (@{$blocklist}) {
            print Dumper $block if ($block->{'Start'} == 53270);
        }

    fill_type_pair( "Intercross-ABHABH_HHA", "Maternal-AHAH", $blocklist );

    print "After Intercross-ABHABH_HHA Maternal-AHAH\n";

        for my $block (@{$blocklist}) {
            print Dumper $block if ($block->{'Start'} == 53270);
        }


    fill_type_pair( "Intercross-ABHABH_HHH", "Maternal-AHAH", $blocklist );
    print "After Intercross-ABHABH_HHH Maternal-AHAH\n";

        for my $block (@{$blocklist}) {
            print Dumper $block if ($block->{'Start'} == 53270);
        }


    fill_type_pair( "Paternal-AHAH",         "Maternal-AHAH", $blocklist );
    print "After Paternal-AHAH Maternal-AHAH\n";

        for my $block (@{$blocklist}) {
            print Dumper $block if ($block->{'Start'} == 53270);
        }



    fill_type_pair( "Paternal-AHAB_AHB",     "Sex-HB",        $blocklist );

    print "After Paternal-AHAB_AHB Sex-HB\n";

        for my $block (@{$blocklist}) {
            print Dumper $block if ($block->{'Start'} == 53270);
        }




}

sub correct_maternal {
    my ( $mattype, $blocklist ) = @_;

    fix_maternal_errors( $mattype, $blocklist );

    my %mat;
    for my $i ( 0 .. $#{$blocklist} ) {
        $mat{ $blocklist->[$i]{$mattype} }{length} +=
          $blocklist->[$i]{'Length'};
        push @{ $mat{ $blocklist->[$i]{$mattype} }{blocks} }, $i;
    }
    my @matl = sort { $mat{$b}{length} <=> $mat{$a}{length} } keys %mat;
    my %merged;
    for my $mata (@matl) {
        next if defined $merged{$mata};
        next if $mata eq $empty;
        next if $mata =~ /\-/;
        for my $matb (@matl) {
            next if $matb eq $empty;
            next if defined $merged{$matb};
            next if $mata eq $matb;
            next if $mat{$mata}{length} < $mat{$matb}{length};
            my $fix = 0;
            $fix = 1 if hamming( $mata,         $matb );
            $fix = 1 if hamming( mirror($mata), $matb );

            if ($fix) {
                $merged{$matb} = $mata;
                for my $b ( @{ $mat{$matb}{blocks} } ) {
                    $blocklist->[$b]{$mattype} = $mata;
                }
            }
        }
    }
}

sub fix_maternal_errors {
    my ( $mattype, $blocklist ) = @_;

    my %mat;
    for my $i ( 0 .. $#{$blocklist} ) {
        next if $blocklist->[$i]{$mattype} eq $empty;
        $mat{ $blocklist->[$i]{'Scaffold'} }{ $blocklist->[$i]{'Start'} }
          {pattern} = $blocklist->[$i]{$mattype};
        $mat{ $blocklist->[$i]{'Scaffold'} }{ $blocklist->[$i]{'Start'} }{block}
          = $i;
    }

    for my $scf ( keys %mat ) {
        next if keys %{ $mat{$scf} } < 3;
        my @starts = sort { $a <=> $b } keys %{ $mat{$scf} };
        for my $i ( 1 .. $#starts - 1 ) {
            next
              if $mat{$scf}{ $starts[ $i - 1 ] }{pattern} eq
              $mat{$scf}{ $starts[$i] }{pattern};
            my @prev = split //, $mat{$scf}{ $starts[ $i - 1 ] }{pattern};
            my @cur  = split //, $mat{$scf}{ $starts[$i] }{pattern};
            my @next = split //, $mat{$scf}{ $starts[ $i + 1 ] }{pattern};
            for my $c ( 0 .. $#prev ) {
                $cur[$c] = $prev[$c]
                  if $prev[$c] =~ /[ABH]/
                  && $prev[$c] eq $next[$c]
                  && $cur[$c] ne $prev[$c];
            }
            my $newcur = join '', @cur;
            $blocklist->[ $mat{$scf}{ $starts[$i] }{block} ]{$mattype} =
              $newcur;
            $mat{$scf}{ $starts[$i] }{pattern} = $newcur;
        }
    }
}

sub fill_type_pair {
    my ( $typea, $typeb, $blocklist ) = @_;

    my %ab;

    for my $i ( 0 .. $#{$blocklist} ) {
        $ab{ $blocklist->[$i]{$typea} }{ $blocklist->[$i]{$typeb} }{length} +=
          $blocklist->[$i]{'Length'};
        push
          @{ $ab{ $blocklist->[$i]{$typea} }{ $blocklist->[$i]{$typeb} }{blocks}
          }, $i;
    }

    for my $apat ( keys %ab ) {
        next if $apat eq $empty;
        next if $apat =~ /\-/;
        next if keys %{ $ab{$apat} } == 1;
        my $maxbpat   = "";
        my $maxbcount = 0;
        for my $bpat ( keys %{ $ab{$apat} } ) {
            if ( $ab{$apat}{$bpat}{length} > $maxbcount ) {
                $maxbpat   = $bpat;
                $maxbcount = $ab{$apat}{$bpat}{length};
            }
        }
        next if $maxbpat eq $empty;
        next if $maxbpat =~ /\-/;

        for my $bpat ( keys %{ $ab{$apat} } ) {
            next if $bpat eq $maxbpat;

            my $fix = 0;
            $fix = 1 if $bpat eq $empty;
            $fix = 1 if hamming( $maxbpat, $bpat );
            $fix = 1 if hamming( mirror($maxbpat), $bpat );

            if ($fix) {
                for my $b ( @{ $ab{$apat}{$bpat}{blocks} } ) {
                    $blocklist->[$b]{$typeb} = $maxbpat;
                }
            }
        }
    }
}

sub mirror {
    my $pat = shift;
    if (   ( $pat =~ 'A' and $pat =~ 'B' and $pat =~ 'H' )
        or ( $pat =~ 'A' and $pat =~ 'B' and $pat !~ 'H' ) )
    {
        $pat =~ tr/AB/BA/;
    }
    if ( $pat =~ 'A' and $pat =~ 'H' and $pat !~ 'B' ) {
        $pat =~ tr /AH/HA/;
    }
    if ( $pat =~ 'B' and $pat =~ 'H' and $pat !~ 'A' ) {
        $pat =~ tr /HB/BH/;
    }
    return $pat;
}

sub hamming {
    my ( $a, $b ) = @_;
    my $hamming = ( $a ^ $b ) =~ tr/\001-\255//;
    return $hamming <= 6;
}

sub load_blocks {
    my $blockfilename = shift;

    my @blocklist;
    my %types;
    my @header;
    my %matpat;

    open my $blocksfile, '<', $blockfilename
      or croak "Can't open $blockfilename: $OS_ERROR\n";
    while ( my $blocksline = <$blocksfile> ) {
        chomp $blocksline;
        my @f = split /\t/, uncolor $blocksline;

        if ( $blocksline =~ /^Scaffold/ ) {
            next if @header ne 0;
            for my $i ( 0 .. $#f ) {
                $f[$i] =~ s/(\s+)//;
                push @header, $f[$i];
                $types{ $f[$i] } = $i if $i > 3;
            }
        }
        else {
            my %block;
            for my $i ( 0 .. $#f ) {
                $block{ $header[$i] } = $f[$i];
                $block{ $header[$i] } = phase( $block{ $header[$i] } )
                  if defined $types{ $header[$i] } && $f[$i] !~ /^( +)$/;
            }
            push @blocklist, \%block;
        }
    }
    close $blocksfile;
    return ( \@header, \%types, \@blocklist );
}