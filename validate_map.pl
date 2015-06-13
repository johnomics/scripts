#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Memoize;
use DBD::SQLite;
use Parallel::ForkManager;
use Term::ExtendedColor qw/:all/;

use List::Util qw/min max/;

$OUTPUT_AUTOFLUSH = 1;

memoize "mirror";
memoize "phase";
memoize "consistent";
memoize "linked";
memoize "calc_LR";

my %args;
$args{input}   = "";
$args{verbose} = "";
$args{threads} = 1;

my $options_okay = GetOptions(
    'input=s'   => \$args{input},
    'threads=i' => \$args{threads},
    'verbose'   => \$args{verbose}
);

croak "No database!" if $args{input} eq "";

my $metadata = { input => $args{input}, verbose => $args{verbose}, threads => $args{threads} };

print STDERR "Loading map...\n";

my $markers = load_map( $args{input}, $metadata );

print STDERR "Loading blocks...\n";
my $blockdb = load_blocks( $args{input}, $markers, $metadata );

my $empty_blocks = scaffold_summary( $blockdb, $metadata );

fill_empty_blocks( $empty_blocks, $markers, $metadata );
print STDERR "Done\n";

exit;

sub load_map {
    my ( $input, $metadata ) = @_;

    my %markers;

    print STDERR "Load map from database...\n" if $metadata->{verbose};
    my $dbh = DBI->connect( "dbi:SQLite:dbname=$input", "", "" );

    my $sth = $dbh->prepare("SELECT * FROM chromosome_map ORDER BY chromosome, cm");
    my $maplist = $dbh->selectall_arrayref( $sth, { Slice => {} } );
    $sth->finish;
    $dbh->disconnect;
    
    for my $marker (@{$maplist}) {
        $marker->{clean} = uncolor $marker->{clean};

        my $i1 = create_intercross( $marker->{print}, $marker->{original} );
        my $i2 = create_intercross( $marker->{print}, mirror($marker->{original}) );

        map { s/B/H/g if /^[AB]+$/ } ( $marker->{print}, $marker->{clean} );

        my $ref = { chromosome => $marker->{chromosome}, cM => $marker->{cm}, maternal => $marker->{print} };
        for my $pattern ( $marker->{print}, $marker->{original}, $marker->{clean} ) {
            next if defined $markers{$pattern};
            $ref->{paternal} = $pattern eq $marker->{print} ? $metadata->{empty} : $marker->{clean};
            $markers{$pattern} = $ref;
            $markers{ mirror($pattern) } = $ref;
        }

        # Add intercross patterns
        $markers{$i1}           = $ref;
        $markers{ mirror($i1) } = $ref;
        $markers{$i2}           = $ref;
        $markers{ mirror($i2) } = $ref;

        add_presence_absence_patterns( $i1, \%markers, $ref );
        add_presence_absence_patterns( $i2, \%markers, $ref );
    }

    # Make empty pattern by checking first Maternal pattern (could be any pattern)
    my $pattern = ( keys %markers )[0];
    $metadata->{samples} = length $pattern;
    $metadata->{empty}   = ' ' x $metadata->{samples};

    \%markers;
}

sub add_presence_absence_patterns {
    my ( $pattern, $markers, $ref ) = @_;

    my $paa = $pattern;
    my $pab = $pattern;
    $paa =~ s/[AH]/0/g;
    $paa =~ s/B/\./g;
    $pab =~ s/[BH]/0/g;
    $pab =~ s/A/\./g;

    $markers->{$paa} = $ref;
    $markers->{$pab} = $ref;

    $paa =~ s/0/1/g;
    $pab =~ s/0/1/g;

    $markers->{$paa} = $ref;
    $markers->{$pab} = $ref;

    return;
}

sub load_blocks {
    my ( $input, $markers, $metadata ) = @_;

    my $blocklist = [];
    my $header;
    my $types;

    print STDERR "Load blocks from database...\n" if $metadata->{verbose};
    my $dbh = DBI->connect( "dbi:SQLite:dbname=$input", "", "" );

    my $sth = $dbh->prepare(
        "SELECT scaffold, start, end, length, Maternal, Paternal FROM cleanblocks ORDER BY scaffold, start");
    my $fileblocklist = $dbh->selectall_arrayref( $sth, { Slice => {} } );
    push @{$blocklist}, @{$fileblocklist};
    $sth->finish;

    $dbh->disconnect;

    my $blockdb = DBI->connect( "dbi:SQLite:dbname=:memory:", "", "" );
    $sth = $blockdb->prepare(
"CREATE TABLE mapblocks (scaffold text, start integer, end integer, length integer, Maternal text, Paternal text)"
    );
    $sth->execute;
    my $insert_handle = $blockdb->prepare_cached('INSERT INTO mapblocks VALUES (?,?,?,?,?,?)');

    for my $block ( @{$blocklist} ) {
        if ( defined $markers->{ $block->{'Paternal'} } ) {
            $block->{'Paternal'} = $markers->{ $block->{'Paternal'} }{paternal};
        }
        $insert_handle->execute(
            $block->{scaffold}, $block->{start},    $block->{end},
            $block->{length},   $block->{Maternal}, $block->{Paternal}
        );
    }

    $sth = $blockdb->prepare("CREATE INDEX mapblocks_maternal on mapblocks (Maternal)");
    $sth->execute;

    $sth = $blockdb->prepare("CREATE INDEX mapblocks_paternal on mapblocks (Paternal)");
    $sth->execute;

    $sth = $blockdb->prepare("CREATE INDEX mapblocks_scaffold on mapblocks (scaffold)");
    $sth->execute;

    $blockdb;
}

sub scaffold_summary {
    my ( $blockdb, $metadata ) = @_;

    my $empty = $metadata->{empty};

    my $scaffolds = $blockdb->selectcol_arrayref('select distinct scaffold from mapblocks');

    print STDERR "Generating summary...\n";
    my %block_stats;
    my %scaffold_stats;
    my @empty_blocks;
    for my $scaffold ( @{$scaffolds} ) {
        my $blocks =
          $blockdb->selectall_arrayref( "select * from mapblocks where scaffold=\"$scaffold\" order by start",
            { Slice => {} } );
        my $blocknum = @{$blocks};

        my $scaffold_length;
        my $maternal = 0;
        my $paternal = 0;

        if ( $blocknum == 1 ) {
            add_stat( 'single', $blocks->[0], \%block_stats );
            add_stat( 'single', $blocks->[0], \%scaffold_stats );
            $scaffold_length = $blocks->[0]{length};
            push @empty_blocks, $blocks->[0];
        }
        else {
            for my $i ( 0 .. $#{$blocks} ) {
                my $block = $blocks->[$i];
                my $blank = 0;

                $scaffold_length += $block->{length};
                if ( $i == 0 or $i == $#{$blocks} ) {
                    add_stat( 'end', $block, \%block_stats );
                    $blank++;
                }
                else {
                    if ( $block->{Maternal} eq $empty ) {
                        add_stat( 'middle', $block, \%block_stats );
                        $blank++;
                    }
                    else {
                        $maternal++;
                        if ( $block->{Paternal} eq $empty ) {
                            add_stat( 'maternal', $block, \%block_stats );
                        }
                        else {
                            add_stat( 'assigned', $block, \%block_stats );
                            $paternal++;
                        }
                    }
                }
                push @empty_blocks, $block if $blank;
            }
        }

        add_stat( 'mapped', { length => $scaffold_length }, \%scaffold_stats ) if $maternal or $paternal;
        add_stat( 'maternal', { length => $scaffold_length }, \%scaffold_stats ) if $maternal and !$paternal;
    }

    print "Scaffolds:\n";
    for my $stat ( sort { $scaffold_stats{$b}{length} <=> $scaffold_stats{$a}{length} } keys %scaffold_stats ) {
        printf "%10s\t%6d\t%10d\n", $stat, $scaffold_stats{$stat}{count}, $scaffold_stats{$stat}{length};
    }
    print "Blocks:\n";
    my $total_blocks = 0;
    my $total_length = 0;
    for my $stat ( sort { $block_stats{$b}{length} <=> $block_stats{$a}{length} } keys %block_stats ) {
        printf "%10s\t%6d\t%10d\n", $stat, $block_stats{$stat}{count}, $block_stats{$stat}{length};
        $total_blocks += $block_stats{$stat}{count};
        $total_length += $block_stats{$stat}{length};
    }
    printf "TOTAL     \t%6d\t%10d\n", $total_blocks, $total_length;

    \@empty_blocks;
}

sub add_stat {
    my ( $type, $block, $stats ) = @_;
    $stats->{$type}{count}++;
    $stats->{$type}{length} += $block->{length};
    $block->{status} = $type;
}

sub fill_empty_blocks {
    my ( $empty_blocks, $markers, $metadata ) = @_;

    my $partitions = get_partitions( $empty_blocks, $metadata->{threads} );

    my $pm = new Parallel::ForkManager( $metadata->{threads} );

    for my $part ( 1 .. $metadata->{threads} ) {
        $pm->start and next;

        print "$part: ", scalar @{ $partitions->{$part} }, "\n";
        fill_partition_blocks( $partitions->{$part}, $markers, $metadata );

        $pm->finish;
    }
    $pm->wait_all_children;
}

sub get_partitions {
    my ( $empty_blocks, $threads ) = @_;

    my %partitions;
    my $part = 1;

    my $length;
    for my $block ( @{$empty_blocks} ) {
        $length += $block->{length};
    }
    my $part_threshold = $length / $threads;

    my $part_length = 0;
    for my $block ( sort { $b->{length} <=> $a->{length} } @{$empty_blocks} ) {
        if ( $part_length > $part_threshold ) {
            $part_length = 0;
            $part++;
        }
        push @{ $partitions{$part} }, $block;
        $part_length += $block->{length};
    }

    \%partitions;
}

sub fill_partition_blocks {
    my ( $blocks, $markers, $metadata ) = @_;

    my $dbh = DBI->connect( "dbi:SQLite:dbname=$metadata->{input}", "", "" );

    my %coverage;
    my %blocks_covered;
    for my $block ( @{$blocks} ) {
        my $sth = $dbh->prepare(
"SELECT position, pattern FROM markers where scaffold=\"$block->{scaffold}\" and position >= $block->{start} and position <= $block->{end} and marker_type=\"Reject\" ORDER BY scaffold, position"
        );
        my $rejects = $dbh->selectall_arrayref( $sth, { Slice => {} } );
        $sth->finish;

        my %patterns;
        my %cMs;
        for my $reject ( @{$rejects} ) {
            my $pattern = $reject->{pattern};
            next if $pattern =~ /^0+$/ or $pattern =~ /^1+$/ or $pattern =~ /^\.+$/ or $pattern =~ /^H+$/;
            my $valid_pattern = check_reject( $pattern, $markers );

            $patterns{$pattern}{count}++;
            $patterns{$pattern}{min} = $patterns{$pattern}{min} // $reject->{position};
            $patterns{$pattern}{max} = $reject->{position};

            if ( defined $valid_pattern ) {
                $patterns{$pattern}{chromosome} = $markers->{$valid_pattern}{chromosome};
                $patterns{$pattern}{cM}         = $markers->{$valid_pattern}{cM};
                $patterns{$pattern}{valid}      = $valid_pattern;

                my $chr_cM = "$markers->{ $valid_pattern }{chromosome}:$markers->{ $valid_pattern }{cM}";
                $cMs{$chr_cM}{count}++;
                $cMs{$chr_cM}{min} = $cMs{$chr_cM}{min} // $reject->{position};
                $cMs{$chr_cM}{max} = $reject->{position};
            }
        }

        validate_block_markers( \%cMs );

        for my $cM (keys %cMs) {
            $coverage{$block->{status}} += $cMs{$cM}{length}
        }
        $blocks_covered{$block->{status}}++ if keys %cMs != 0;
        
#        output_patterns( $block, $rejects, \%patterns, \%cMs );

    }
    $dbh->disconnect;

    my $total_coverage = 0;
    my $total_blocks = 0;
    for my $status (sort keys %blocks_covered) {
        print "$status\t$blocks_covered{$status}\t$coverage{$status}\n";
        $total_coverage += $coverage{$status};
        $total_blocks += $blocks_covered{$status};
    }
    print "Total\t$total_blocks\t$total_coverage\n";
}

sub check_reject {
    my ( $reject, $markers ) = @_;

    return $reject if defined $markers->{$reject};

    return;
}

sub validate_block_markers {
    my ($cMs) = @_;

    for my $cM ( keys %{$cMs} ) {
        $cMs->{$cM}{length} = $cMs->{$cM}{max} - $cMs->{$cM}{min} + 1;
        $cMs->{$cM}{density} = sprintf "%6.2f", $cMs->{$cM}{length} / $cMs->{$cM}{count};
    }
    my @cMs = sort { $cMs->{$b}{count} <=> $cMs->{$a}{count} } keys %{$cMs};

    for my $cM (@cMs) {
        next if !defined $cMs->{$cM};
        for my $cM2 (@cMs) {
            next if !defined $cMs->{$cM2};
            next if $cM eq $cM2;
            if ( overlap( $cMs->{$cM}, $cMs->{$cM2} ) ) {
                delete $cMs->{$cM2};
            }
        }

    }
}

sub overlap {
    my ( $cM1, $cM2 ) = @_;
    return ( $cM1->{min} < $cM2->{min} and $cM1->{max} > $cM2->{min} )
      or ( $cM1->{min} < $cM2->{max}   and $cM1->{max} > $cM2->{max} );
}

sub output_patterns {
    my ( $block, $rejects, $patterns, $cMs ) = @_;
    my $snps = @{$rejects};
    my $output = sprintf "%-16s\t%8d\t%8d\t%8d\t%6d\t%8s\n", $block->{scaffold}, $block->{start}, $block->{end},
      $block->{length},
      $snps, $block->{status};

    my @sorted_patterns = sort { $patterns->{$b}{count} <=> $patterns->{$a}{count} } keys %{$patterns};

    for my $p ( 0 .. $#sorted_patterns ) {
        last if $p == 10;
        my $sp     = $patterns->{ $sorted_patterns[$p] };
        my $len    = $sp->{max} - $sp->{min} + 1;
        my $len_pc = sprintf "%5.2f", $len / $block->{length} * 100;
        $output .= "\t$sorted_patterns[$p]\t$sp->{min}\t$sp->{max}\t$len\t$len_pc\t$sp->{count}";
        $output .= "\t$sp->{chromosome}\t$sp->{cM}" if defined $sp->{cM};
        $output .= "\n";
    }

    for my $cM ( sort { $cMs->{$b}{count} <=> $cMs->{$a}{count} } keys %{$cMs} ) {
        my $c      = $cMs->{$cM};
        my $len    = $c->{max} - $c->{min} + 1;
        my $len_pc = sprintf "%5.2f", $len / $block->{length} * 100;
        $output .= "$cM\t$c->{count}\t$c->{min}\t$c->{max}\t$len\t$len_pc\n";
    }
    print $output;
    return;
}

## LIBRARY FUNCTIONS

sub create_intercross {
    my ( $maternal, $paternal ) = @_;
    my @m = split //, $maternal;
    my @p = split //, $paternal;

    my @i = map { ( $m[$_] eq $p[$_] ) ? $m[$_] : 'H' } 0 .. $#m;

    join '', @i;
}

sub intercross_in_phase {
    my ( $i1, $i2 ) = @_;
    my @i1 = split //, $i1;
    my @i2 = split //, $i2;

    my $h_match    = 0;
    my $h_mismatch = 0;
    for my $i ( 0 .. $#i1 ) {
        if ( $i1[$i] eq 'H' ) {
            if ( $i2[$i] eq 'H' ) {
                $h_match++;
            }
            else {
                $h_mismatch++;
            }
        }
    }
    $h_match > $h_mismatch;
}

sub merge_coupled_intercross {
    my ( $i1, $i2 ) = @_;
    my @i1 = split //, $i1;
    my @i2 = split //, $i2;

    my @out;
    for my $i ( 0 .. $#i1 ) {
        $out[$i] = $i1[$i] eq $i2[$i] ? $i1[$i] : '-';
    }
    my $out = join '', @out;

    $out;
}

sub merge_repulsion_intercross {
    my ( $i1, $i2 ) = @_;

    # I1 should be the pattern that starts with A; I2 starts with H
    if ( $i2 =~ /^A/ ) {
        my $t = $i1;
        $i1 = $i2;
        $i2 = $t;
    }

    my @i1 = split //, $i1;
    my @i2 = split //, $i2;

    my @out1;
    my @out2;
    for my $i ( 0 .. $#i1 ) {
        if (   ( $i1[$i] eq 'A' and $i2[$i] eq 'B' )
            or ( $i1[$i] eq 'B' and $i2[$i] eq 'A' ) )
        {
            $out1[$i] = '-';
            $out2[$i] = '-';
        }
        elsif ( $i1[$i] eq $i2[$i] or $i2[$i] eq 'H' ) {    # A,A; B,B; H,H; A,H; B,H
            $out1[$i] = $i1[$i];
            $out2[$i] = $i1[$i];
        }
        elsif ( $i1[$i] eq 'H' ) {
            $out1[$i] = $i2[$i];
            $out2[$i] = $i2[$i] eq 'A' ? 'B' : 'A';
        }
    }

    my $out1 = join '', @out1;
    my $out2 = join '', @out2;

    ( $out1, $out2 );
}

sub consistent {
    my ( $pattern1, $pattern2 ) = @_;

    my @pattern1 = split //, $pattern1;
    my @pattern2 = split //, $pattern2;
    my $distance = 0;
    for my $i ( 0 .. $#pattern1 ) {
        return 0
          if $pattern1[$i] ne $pattern2[$i]
          and $pattern1[$i] ne 'H'
          and $pattern2[$i] ne 'H'
          and $pattern1[$i] ne '.'
          and $pattern2[$i] ne '.';
    }
    return 1;
}

sub phase {
    my $pat = shift;
    my @checkbases;
    if ( $pat =~ 'A' && $pat =~ 'B' && $pat =~ 'H' ) {
        @checkbases = ( 'A', 'B' );
    }
    elsif ( $pat =~ 'A' && $pat =~ 'B' && $pat !~ 'H' ) {
        @checkbases = ( 'A', 'B' );
    }
    elsif ( $pat =~ 'B' && $pat =~ 'H' && $pat !~ 'A' ) {
        @checkbases = ( 'B', 'H' );
    }
    return " " x length($pat) if !@checkbases;

    my @p = split //, $pat;
    my $out = "";
    my %trans;

    for my $gt (@p) {
        my $added = 0;
        for my $bi ( 0, 1 ) {
            if ( $gt eq $checkbases[$bi] ) {
                if ( !defined $trans{$gt} ) {
                    $trans{ $checkbases[$bi] } = 'A';
                    my $other = $bi == 0 ? 1 : 0;
                    $trans{ $checkbases[$other] } = 'B';
                }
                $out .= $trans{$gt};
                $added++;
            }
        }
        $out .= $gt =~ /[ ~\.]/ ? '-' : $gt if !$added;
    }

    $out;
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

    $pat;
}

sub hamming {
    return ( $_[0] ^ $_[1] ) =~ tr/\001-\255//;
}

## Linkage tests

sub linked {
    my ( $a, $b ) = @_;
    ( $a =~ /H/ and $b =~ /H/ )
      ? linked_intercross( $a, $b )
      : linked_backcross( $a, $b );
}

sub linked_backcross {
    my ( $a, $b ) = @_;

    # If either pattern is Intercross, remove individuals where Intercross is H
    if ( $a =~ /H/ or $b =~ /H/ ) {
        my @a = split //, $a;
        my @b = split //, $b;
        my @new_a;
        my @new_b;
        for my $i ( 0 .. $#a ) {
            next if $a[$i] =~ /H/ or $b[$i] =~ /H/;
            push @new_a, $a[$i];
            push @new_b, $b[$i];
        }
        $a = join '', @new_a;
        $b = join '', @new_b;
    }
    my $N = length $a;
    my $R = hamming( $a, $b );
    my $r = $R / $N;
    my $LOD =
      ( ( $N - $R ) * ( log( 1 - $r ) / log(10) ) ) + ( $R * ( log($r) / log(10) ) ) + ( $N * ( log(2) / log(10) ) );
    return $LOD > 3;
}

sub linked_intercross {
    my ( $a, $b ) = @_;

    my $N      = length $a;
    my $haps   = get_haps( $a, $b );
    my $r      = f2_em( $haps, $N, 0.25 );
    my $p_re   = calc_p_re($r);
    my $R      = calc_R( $haps, $p_re );
    my $LR_r   = calc_LR( $r, $R, $N );
    my $LR_ind = calc_LR( 0.5, $R, $N );
    my $LOD    = log( $LR_r / $LR_ind ) / log(10);
    return $LOD > 3;
}

sub f2_em {
    my ( $haps, $N, $r ) = @_;
    for my $i ( 1 .. 100 ) {
        my $p_re  = calc_p_re($r);
        my $R     = calc_R( $haps, $p_re );
        my $S     = calc_S( $haps, $p_re );
        my $new_r = $R / ( 2 * $N );
        last if sprintf( "%.5f", $new_r ) eq sprintf( "%.5f", $r );
        $r = $new_r;
    }

    $r;
}

sub get_haps {
    my ( $a, $b ) = @_;

    my @a = split //, $a;
    my @b = split //, $b;

    my %haps;
    for my $i ( 'A', 'B', 'H' ) {
        for my $j ( 'A', 'B', 'H' ) {
            $haps{ $i . $j } = 0;
        }
    }
    for my $i ( 0 .. $#a ) {
        my $hap = $a[$i] . $b[$i];
        $haps{$hap}++;
    }
    return \%haps;
}

sub calc_p_re {
    my $r = shift;
    ( $r**2 ) / ( ( ( 1 - $r )**2 ) + ( $r**2 ) )

}

sub calc_R {
    my ( $haps, $p_re ) = @_;

    0 * $haps->{'AA'} +
      1 * $haps->{'AH'} +
      2 * $haps->{'AB'} +
      1 * $haps->{'HA'} +
      2 * $haps->{'HH'} * $p_re +
      1 * $haps->{'HB'} +
      2 * $haps->{'BA'} +
      1 * $haps->{'BH'} +
      0 * $haps->{'BB'};
}

sub calc_S {
    my ( $haps, $p_re ) = @_;

    2 * $haps->{'AA'} +
      1 * $haps->{'AH'} +
      0 * $haps->{'AB'} +
      1 * $haps->{'HA'} +
      2 * $haps->{'HH'} * ( 1 - $p_re ) +
      1 * $haps->{'HB'} +
      0 * $haps->{'BA'} +
      1 * $haps->{'BH'} +
      2 * $haps->{'BB'};
}

sub calc_LR {
    my ( $r, $R, $N ) = @_;
    ( 1 - $r )**( $N - $R ) * ( $r**$R );
}

sub int_hamming {

    # Calculate hamming distance between intercross and paternal pattern
    my ( $int, $pat ) = @_;
    my $hamming = 0;
    my @i       = split //, $int;
    my @p       = split //, $pat;
    for my $b ( 0 .. $#i ) {
        next
          if $i[$b] eq 'H'
          or $i[$b] eq '-'
          or $p[$b] eq 'H'
          or $p[$b] eq '-';
        $hamming++ if $i[$b] ne $p[$b];
    }

    $hamming;
}

sub int_match {

    # Check for complete match between intercross and other pattern,
    # ignoring missing bases in intercross pattern
    my ( $intercross, $pattern ) = @_;
    my @i = split //, $intercross;
    my @p = split //, $pattern;
    my $match = 1;
    for my $a ( 0 .. $#i ) {
        next if $i[$a] =~ /[H\-]/ or $p[$a] =~ /[H\-]/;
        if ( $i[$a] ne $p[$a] ) { $match = 0; last; }
    }

    $match;
}

sub separate_intercross {
    my ( $pattern, $intercross ) = @_;

    $intercross = mirror($intercross)
      if int_hamming( $pattern, mirror($intercross) ) < int_hamming( $pattern, $intercross );
    my @intercross = split //, $intercross;
    my @pattern    = split //, $pattern;
    my @complement =
      map {
            $intercross[$_] eq '-' ? '-'
          : $intercross[$_] eq 'H' ? ( $pattern[$_] eq 'A' ? 'B' : 'A' )
          : $intercross[$_]
      } 0 .. $#intercross;
    my $complement = join '', @complement;

    $complement;
}
