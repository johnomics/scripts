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

fill_blocks( $header, $types, $blocklist );

correct_maternal( "Maternal-AHAH", $blocklist );

collapse( $blocklist, $types );

make_maps( $blocklist, $args{output} );
exit;
my %presents;
for my $block ( @{$blocklist} ) {
    map { print "$block->{$_}\t" } @{$header};
    print "\n";
    my $present = "";
    for my $f ( @{$header} ) {
        if ( $f =~ /[ABH]{2}/ ) {
            $present .= $block->{$f} !~ /^( +)$/ ? '1' : '.';
        }
    }
    $presents{$present}{length} += $block->{'Length'};
    $presents{$present}{blocks}++;
}

for my $present ( sort keys %presents ) {
    print
      "$present\t$presents{$present}{blocks}\t$presents{$present}{length}\n";
}

sub make_maps {
    my ( $blocklist, $output ) = @_;

    my %matpat;
    my $glength = 0;
    my %pattern_block;
    for my $block ( @{$blocklist} ) {
        my $mattrans = make_transpat( $block->{"Maternal-AHAH"} );
        my $pattrans = make_transpat( $block->{"Paternal-AHAH"} );
        next if $mattrans =~ /( )/;
        if ( $pattrans =~ /^( +)$/ ) {
            $pattrans = convert_intercross($block);
        }
        next if $pattrans =~ /( )/;
        $matpat{$mattrans}{$pattrans}{length} += $block->{'Length'};
        $matpat{$mattrans}{$pattrans}{blocks}++;
        $pattern_block{$pattrans}{ $block->{'Scaffold'} }{ $block->{'Start'} }
          = $block->{'End'};
    }

    for my $mat ( keys %matpat ) {
        for my $pat ( keys %{ $matpat{$mat} } ) {
            delete $matpat{$mat}{$pat}
              if (  $matpat{$mat}{$pat}{length} < 1000
                and $matpat{$mat}{$pat}{blocks} <= 2 );
        }
        next if ( keys %{ $matpat{$mat} } == 0 );
        my $markercode = run_mstmap( $mat, $matpat{$mat}, $output );
        my $map = load_map( $output, $mat );
        for my $lg ( sort keys %{$map} ) {
            for my $cm ( sort { $a <=> $b } keys %{ $map->{$lg} } ) {
                for my $marker ( sort { $a <=> $b } keys %{ $map->{$lg}{$cm} } )
                {
                    my $pattern = $markercode->{$marker};
                    print
"$mat\t$lg\t$cm\t$marker\t$pattern\t$matpat{$mat}{$pattern}{blocks}\t$matpat{$mat}{$pattern}{length}\n";
                    for my $scf ( keys %{ $pattern_block{$pattern} } ) {
                        for
                          my $start ( keys %{ $pattern_block{$pattern}{$scf} } )
                        {
                            my $length =
                              $pattern_block{$pattern}{$scf}{$start} -
                              $start + 1;
                            print "\t$scf\t$start\t$pattern_block{$pattern}{$scf}{$start}\t$length\n";
                        }
                    }
                }
            }
        }

        print "\n";
    }

    print "\n";

}

sub convert_intercross {
    my ($block) = @_;

    for my $ic ( "Intercross-ABHABH_HHA", "Intercross-ABHABH_HHH" ) {
        next if $block->{$ic} =~ /^( +)$/;
        my $mat = $block->{'Maternal-AHAH'};
        my @int = split //, $block->{$ic};
        my @mat = split //, $mat;
        my @pat =
          map { $int[$_] eq 'H' ? ( $mat[$_] eq 'A' ? 'H' : 'A' ) : $mat[$_] }
          0 .. $#int;
        my $pat = join '', @pat;
        return make_transpat($pat);
    }
    return " ";
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
        map {
            print $mstmapin "\t";
            print $mstmapin $_ eq 'x' ? 'A' : 'B';
        } @gt;
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
    fill_type_pair( "Intercross-ABHABH_HHH", "Intercross-ABHABH_HHA",
        $blocklist );
    fill_type_pair( "Intercross-ABHABH_HHA", "Paternal-AHAH", $blocklist );
    fill_type_pair( "Paternal-AHAH", "Intercross-ABHABH_HHA", $blocklist );
    fill_type_pair( "Intercross-ABHABH_HHA", "Maternal-AHAH", $blocklist );
    fill_type_pair( "Intercross-ABHABH_HHH", "Maternal-AHAH", $blocklist );
    fill_type_pair( "Paternal-AHAH",         "Maternal-AHAH", $blocklist );
    fill_type_pair( "Paternal-AHAB_AHB",     "Sex-HB",        $blocklist );
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
        next if $mata =~ /[~\.]/;
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
        next if $apat =~ /[~.]/;
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
        next if $maxbpat =~ /[~.]/;

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
            }
            push @blocklist, \%block;
        }
    }
    close $blocksfile;
    return ( \@header, \%types, \@blocklist );
}

sub check_error {
    my $p = shift;
    return $p =~ /[\~\.\ ]/ ? $empty : $p;
}
