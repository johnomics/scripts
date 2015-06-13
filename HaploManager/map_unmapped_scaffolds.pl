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

sub permute (&@);

my $PERCENT_THRESHOLD = 25;
my $SNP_THRESHOLD     = 5;

my %args;
$args{markerdb} = "";
$args{cleandb}  = "";
$args{gff}      = "";
$args{output}   = "";
$args{verbose}  = "";
$args{threads}  = 1;

my $options_okay = GetOptions(
    'markerdb=s' => \$args{markerdb},
    'cleandb=s'  => \$args{cleandb},
    'gff=s'      => \$args{gff},
    'output=s'   => \$args{output},
    'threads=i'  => \$args{threads},
    'verbose'    => \$args{verbose}
);

croak "No marker database!" if $args{markerdb} eq "";
croak "Marker database doesn't exist!" if !-e $args{markerdb};

croak "No clean database!" if $args{cleandb} eq "";
croak "Clean database doesn't exist!" if !-e $args{cleandb};

croak "GFF file doesn't exist!" if $args{gff} and !-e $args{gff};

my $metadata = {
    markerdb => $args{markerdb},
    cleandb  => $args{cleandb},
    gff      => $args{gff},
    output   => $args{output},
    verbose  => $args{verbose},
    threads  => $args{threads}
};

print STDERR "Loading map...\n";

my $markers = load_map( $args{cleandb}, $metadata );

my $genes;
$genes = load_genes( $genes, $args{gff} ) if $args{gff};

my $unmapped_scaffolds = load_unmapped_scaffolds( $args{cleandb}, $metadata );

fill_unmapped_scaffolds( $unmapped_scaffolds, $markers, $genes, $metadata );
print STDERR "Done\n";

exit;

sub load_genes {
    my ( $genes, $gff ) = @_;
    print STDERR "Loading genes...\n";

    open my $gff_fh, '<', $gff or croak "Can't open GFF file $gff!\n";

    while ( my $feature = <$gff_fh> ) {
        chomp $feature;
        my @f = split /\t/, $feature;
        next if $f[2] ne 'gene';
        push @{ $genes->{ $f[0] } }, { start => $f[3], end => $f[4], attributes => $f[8] };
    }
    close $gff_fh;

    return $genes;
}

sub load_map {
    my ( $input, $metadata ) = @_;

    print STDERR "Load map from database...\n" if $metadata->{verbose};
    my $dbh = DBI->connect( "dbi:SQLite:dbname=$input", "", "" );

    my $sth = $dbh->prepare("SELECT * FROM chromosome_map ORDER BY chromosome, cm");
    my $maplist = $dbh->selectall_arrayref( $sth, { Slice => {} } );
    $sth->finish;
    $dbh->disconnect;

    # Make empty pattern by checking first Maternal pattern (could be any pattern)
    $metadata->{samples} = length $maplist->[0]{print};
    $metadata->{empty}   = ' ' x $metadata->{samples};

    map { $_->{clean} = uncolor $_->{clean} } @{$maplist};

    interpolate_missing_markers($maplist);

    my %markers;
    for my $marker ( @{$maplist} ) {

        my $i1 = create_intercross( $marker->{print}, $marker->{original} );
        my $i2 = create_intercross( $marker->{print}, mirror( $marker->{original} ) );

        map { s/B/H/g if /^[AB]+$/ } ( $marker->{print}, $marker->{clean} );

        my $ref = { chromosome => $marker->{chromosome}, cM => $marker->{cm}, maternal => $marker->{print} };
        for my $pattern ( $marker->{print}, $marker->{original}, $marker->{clean} ) {
            next if defined $markers{$pattern};
            $ref->{paternal} = $pattern eq $marker->{print} ? $metadata->{empty} : $marker->{clean};
            $markers{$pattern} = $ref;
            $markers{ mirror($pattern) } = $ref;
            if ( $pattern eq $marker->{print} ) {
                $markers{$pattern}{cM} = 'M';
                $markers{ mirror($pattern) }{cM} = 'M';
            }
        }

        # Add intercross patterns
        $markers{$i1}           = $ref;
        $markers{ mirror($i1) } = $ref;
        $markers{$i2}           = $ref;
        $markers{ mirror($i2) } = $ref;

        add_presence_absence_patterns( $i1, \%markers, $ref );
        add_presence_absence_patterns( $i2, \%markers, $ref );
    }

    \%markers;
}

sub interpolate_missing_markers {
    my ($maplist) = @_;

    my %clean;
    for my $marker ( @{$maplist} ) {
        $clean{ $marker->{chromosome} }{ $marker->{cm} }{print} = $marker->{print};
        $clean{ $marker->{chromosome} }{ $marker->{cm} }{clean} = $marker->{clean};
    }

    for my $chromosome ( sort { $a <=> $b } keys %clean ) {
        my @cms = sort { $a <=> $b } keys %{ $clean{$chromosome} };
        my $i = 0;
        while ( $i < @cms - 2 ) {
            my $cmi = $cms[$i];
            my $cmj = $cms[ $i + 1 ];

            my @cleani = split //, $clean{$chromosome}{$cmi}{clean};
            my @cleanj = split //, $clean{$chromosome}{$cmj}{clean};
            my @diffs;
            for my $c ( 0 .. @cleani - 1 ) {
                push @diffs, $c if $cleani[$c] ne $cleanj[$c];
            }

            if ( @diffs > 1 ) {

                #                print "$chromosome\t$cmi\t$clean{$chromosome}{$cmi}{clean}\n";
                my @orders;
                permute { push @orders, \@_; } @diffs;
                for my $order (@orders) {
                    add_variants( $maplist, $chromosome, $cmi, $clean{$chromosome}{$cmi}, $order );
                }

                #                print "$chromosome\t$cmj\t$clean{$chromosome}{$cmj}{clean}\n";
                #                print "\n";
            }
            $i++;
        }
    }
}

sub add_variants {
    my ( $maplist, $chromosome, $cm, $cmi, $order ) = @_;
    for my $sample ( @{$order} ) {
        next if $sample eq $order->[-1];
        $cm += 0.735;
        substr( $cmi->{clean}, $sample, 1 ) =~ tr/AB/BA/;
        my %marker = (
            chromosome => $chromosome,
            cm         => $cm,
            print      => $cmi->{print},
            clean      => $cmi->{clean},
            original   => $cmi->{clean},
            length     => 0
        );
        push @{$maplist}, \%marker;
    }
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

sub load_unmapped_scaffolds {
    my ( $input, $metadata ) = @_;
    my $dbh = DBI->connect( "dbi:SQLite:dbname=$input", "", "" );

    my $sth = $dbh->prepare("SELECT * from scaffold_map");
    my $scaffold_map = $dbh->selectall_arrayref( $sth, { Slice => {} } );
    $sth->finish;
    $dbh->disconnect;

    my %scaffolds;
    for my $part ( @{$scaffold_map} ) {
        if ( $part->{type} eq 'active' or $part->{type} eq 'retained' ) {
            push @{ $scaffolds{ $part->{scaffold} }{ $part->{chromosome} }{ $part->{cm} } }, $part;
        }
    }

    my %unmapped_scaffolds;
    my $unmapped_count;
    my $unmapped_length;
    my $unmapped_parts;
    for my $scaffold ( keys %scaffolds ) {
        if ( ( keys %{ $scaffolds{$scaffold} } == 1 ) and defined $scaffolds{$scaffold}{0} ) {
            $unmapped_count++;
            for my $chromosome ( keys %{ $scaffolds{$scaffold} } ) {
                for my $cm ( keys %{ $scaffolds{$scaffold}{$chromosome} } ) {
                    for my $part ( @{ $scaffolds{$scaffold}{$chromosome}{$cm} } ) {
                        push @{ $unmapped_scaffolds{$scaffold}{parts} }, $part;
                        $unmapped_scaffolds{$scaffold}{length} += $part->{length};
                        $unmapped_length += $part->{length};
                        $unmapped_parts++;
                    }
                }
            }
        }
    }
    print "Found $unmapped_count unmapped scaffolds in $unmapped_parts parts covering $unmapped_length bp from "
      . scalar( keys %scaffolds )
      . " active or retained scaffolds\n";

    \%unmapped_scaffolds;
}

sub fill_unmapped_scaffolds {
    my ( $unmapped_scaffolds, $markers, $genes, $metadata ) = @_;

    my $partitions = get_partitions( $unmapped_scaffolds, $metadata->{threads} );

    my $pm = new Parallel::ForkManager( $metadata->{threads} );

    for my $partition ( 1 .. $metadata->{threads} ) {
        $pm->start and next;

        $pm->finish if !defined $partitions->{$partition};
        print "$partition: ", scalar keys %{ $partitions->{$partition} }, "\n";
        fill_partition_parts( $partitions->{$partition}, $markers, $genes, $metadata );

        $pm->finish;
    }
    $pm->wait_all_children;
}

sub get_partitions {
    my ( $scaffolds, $threads ) = @_;

    my %partitions;
    my $partnum = 1;

    my $length;
    for my $scaffold ( keys %{$scaffolds} ) {
        $length += $scaffolds->{$scaffold}{length};
    }
    my $part_threshold = $length / $threads;

    my $part_length = 0;
    for my $scaffold ( sort { $scaffolds->{$b}{length} <=> $scaffolds->{$a}{length} } keys %{$scaffolds} ) {
        if ( $part_length > $part_threshold ) {
            $part_length = 0;
            $partnum++;
        }
        $partitions{$partnum}{$scaffold} = $scaffolds->{$scaffold};
        $part_length += $scaffolds->{$scaffold}{length};
    }
    \%partitions;
}

sub fill_partition_parts {
    my ( $scaffolds, $markers, $genes, $metadata ) = @_;

    my $dbh = DBI->connect( "dbi:SQLite:dbname=$metadata->{markerdb}", "", "" );

    open my $outfile, '>', $metadata->{output} or croak "Can't open output file $metadata->{output}! $OS_ERROR\n";

    my %coverage;
    my %parts_covered;

    for my $scaffold ( sort { $scaffolds->{$b}{length} <=> $scaffolds->{$a}{length} } keys %{$scaffolds} ) {
        print "\n\n\t$scaffold\t$scaffolds->{$scaffold}{length}\n";

        my @rejects;
        my $genecount = 0;
        my %rejectpos;
        my $parts = $scaffolds->{$scaffold}{parts};
        for my $part ( @{$parts} ) {
            print "\t$part->{scaffold}\t$part->{start}\t$part->{end}\t$part->{length}\t$part->{comments}\n";

            if ( defined $genes->{$scaffold} ) {
                for my $gene ( sort { $a->{start} <=> $b->{start} } @{ $genes->{$scaffold} } ) {
                    if ( $gene->{start} < $part->{end} and $gene->{end} > $part->{start} ) {
                        $genecount++;
                        print "\t\tGene\t$gene->{start}\t$gene->{end}\t$gene->{attributes}\n";
                    }
                }
            }

            my ( $oldname, $astart, $bend, $cstart, $dend, $oldstart, $oldend, $strand ) = split "\t",
              $part->{comments};
            $part->{oldname}   = $oldname;
            $part->{oldstart}  = $oldstart;
            $part->{oldend}    = $oldend;
            $part->{oldstrand} = $strand;
            my $sth = $dbh->prepare(
"SELECT position, pattern FROM markers where scaffold=\"$part->{oldname}\" and position >= $part->{oldstart} and position <= $part->{oldend} and marker_type=\"Reject\" ORDER BY scaffold, position"
            );
            my $rejects = $dbh->selectall_arrayref( $sth, { Slice => {} } );

            for my $reject ( @{$rejects} ) {
                my $newpos;
                if ( $strand == 1 ) {
                    $newpos = $reject->{position} - $oldstart + $part->{start};
                }
                else {
                    $newpos = $part->{end} - ( $reject->{position} - $oldstart );
                }
                $rejectpos{$newpos}{scaffold} = $part->{oldname};
                $rejectpos{$newpos}{position} = $reject->{position};
                $reject->{position}           = $newpos;

            }
            push @rejects, sort { $a->{position} <=> $b->{position} } @{$rejects};
            $sth->finish;
        }

        print "\tSNPs: " . scalar(@rejects) . "\n";
        my %patterns;
        my %cMs;
        my %chrs;
        for my $reject (@rejects) {
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

                my $chr = $markers->{$valid_pattern}{chromosome};
                $chrs{$chr}{count}++;
                $chrs{$chr}{min} = $chrs{$chr}{min} // $reject->{position};
                $chrs{$chr}{max} = $reject->{position};

                my $cM = $markers->{$valid_pattern}{cM};
                if ( $cM ne 'M' ) {
                    my $chrcM = "$chr:$cM";
                    $cMs{$chrcM}{count}++;
                    $cMs{$chrcM}{min} = $cMs{$chrcM}{min} // $reject->{position};
                    $cMs{$chrcM}{max} = $reject->{position};
                }
            }
        }

        validate_block_markers( \%cMs );

        output_patterns( $scaffold, $scaffolds->{$scaffold}{length}, $genecount, \%patterns, \%cMs, \%chrs );

        my $output_blocks =
          make_output_blocks( \%chrs, \%cMs, $scaffold, $scaffolds->{$scaffold}{length}, \%rejectpos, $parts );

        output_blocks( $output_blocks, $parts, $dbh, $outfile );

    }
    $dbh->disconnect;
    close $outfile;
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
    my ( $scaffold, $length, $genecount, $patterns, $cMs, $chrs ) = @_;

    my $output = "";
    my @sorted_patterns = sort { $patterns->{$b}{count} <=> $patterns->{$a}{count} } keys %{$patterns};

    for my $p ( 0 .. $#sorted_patterns ) {
        last if $p == 10;
        my $sp     = $patterns->{ $sorted_patterns[$p] };
        my $len    = $sp->{max} - $sp->{min} + 1;
        my $len_pc = sprintf "%5.2f", $len / $length * 100;
        $output .= "\t$sorted_patterns[$p]\t$sp->{min}\t$sp->{max}\t$len\t$len_pc\t$sp->{count}";
        $output .= "\t$sp->{chromosome}\t$sp->{cM}" if defined $sp->{cM};
        $output .= "\n";
    }

    $output .= "HOME\t$scaffold\t$length\t$genecount";
    $output = output_stats( $chrs, $output, $length );
    $output = output_stats( $cMs,  $output, $length );

    print $output;
    print "\n";

    return;
}

sub output_stats {
    my ( $stats, $output, $length ) = @_;

    for my $stat ( sort { $stats->{$b}{count} <=> $stats->{$a}{count} } keys %{$stats} ) {
        my $c      = $stats->{$stat};
        my $len    = $c->{max} - $c->{min} + 1;
        my $len_pc = sprintf "%5.2f", $len / $length * 100;
        $stats->{$stat}{pc} = $len_pc;
        $output .= "\t$stat\t$c->{count}\t$c->{min}\t$c->{max}\t$len\t$len_pc";
    }
    return $output;

}

sub make_output_blocks {
    my ( $chrs, $cMs, $scaffold, $length, $rejectpos, $parts ) = @_;
    return if keys %{$chrs} == 0 or keys %{$cMs} == 0;
    my $chr   = ( sort { $chrs->{$b}{count} <=> $chrs->{$a}{count} } keys %{$chrs} )[0];
    my $chrcM = ( sort { $cMs->{$b}{count} <=> $cMs->{$a}{count} } keys %{$cMs} )[0];
    my ( $cMchr, $cM ) = split /:/, $chrcM;
    return if $chr ne $cMchr;

    my $h = $chrs->{$chr};
    my $m = $cMs->{$chrcM};
    return if $h->{pc} < $PERCENT_THRESHOLD or $h->{count} < $SNP_THRESHOLD;

    my @borders = ( { pos => $h->{min}, chr => $chr, cM => -1 } );
    if ( $m->{pc} >= $PERCENT_THRESHOLD and $m->{count} >= $SNP_THRESHOLD ) {
        push @borders, { pos => $m->{min}, chr => $chr, cM => $cM };
        push @borders, { pos => $m->{max}, chr => $chr, cM => -1 };
    }

    push @borders, { pos => $h->{max}, chr => 0, cM => -1 };

    my @output_blocks;

    my $curchr    = 0;
    my $curcm     = -1;
    my $curborder = shift @borders;
    my $curoldpos = $rejectpos->{ $curborder->{pos} };
    for my $part ( @{$parts} ) {
        my $curdir = $part->{oldstrand};
        my $curstart = $curdir == 1 ? $part->{oldstart} : $part->{oldend};
        next if $part->{oldname} =~ /PB/;
        while ( $curoldpos->{scaffold} eq $part->{oldname}
            and $curoldpos->{position} >= $part->{oldstart}
            and $curoldpos->{position} <= $part->{oldend} )
        {
            my $curend = $curoldpos->{position} - $curdir;
            if ( $curdir == 1 and $curend > $curstart or $curdir == -1 and $curend < $curstart ) {
                ( $curstart, $curend ) = order( $curstart, $curend );
                push @output_blocks,
                  { name => $part->{oldname}, start => $curstart, end => $curend, chr => $curchr, cm => $curcm };
            }
            $curstart = $curoldpos->{position};
            $curchr   = $curborder->{chr};
            $curcm    = $curborder->{cM};
            if (@borders) {
                $curborder = shift @borders;
                $curoldpos = $rejectpos->{ $curborder->{pos} };
            }
            else {
                last;
            }
        }
        my $curend = $curdir == 1 ? $part->{oldend} : $part->{oldstart};
        ( $curstart, $curend ) = order( $curstart, $curend );
        push @output_blocks,
          { name => $part->{oldname}, start => $curstart, end => $curend, chr => $curchr, cm => $curcm };
    }
    return \@output_blocks;
}

sub output_blocks {
    my ( $output_blocks, $parts, $dbh, $outfile ) = @_;

    return if not defined $output_blocks;
    for my $part ( @{$parts} ) {
        print "$part->{oldname}\t$part->{oldstart}\t$part->{oldend}\t$part->{oldstrand}\n";
        my @partoutblocks;
        my $newchrfound = 0;
        for my $outblock ( sort {$a->{name} cmp $b->{name} or $a->{start} <=> $b->{start}} @{$output_blocks} ) {
            next
              if $outblock->{name} ne $part->{oldname}
              or $outblock->{start} > $part->{oldend}
              or $outblock->{end} < $part->{oldstart};
            push @partoutblocks, $outblock;
            $newchrfound++ if $outblock->{chr} != 0;
        }

        next if not $newchrfound;
        
        my $sth = $dbh->prepare("select * from scaffold_map where scaffold='$part->{oldname}' order by start");
        my $blocks = $dbh->selectall_arrayref( $sth, { Slice => {} } );
        for my $block ( @{$blocks} ) {
            next if $block->{start} > $part->{oldend} or $block->{end} < $part->{oldstart};
            croak
"Block not empty!\n\tPart\t$part->{oldname}\t$part->{oldstart}\t$part->{oldend}\t$part->{oldstrand}\nBlock\t$block->{start}\t$block->{end}\t$block->{length}\t$block->{chromosome}\t$block->{cm}\n"
              if $block->{chromosome} != 0 and $block->{cm} != -1;
            next if !@partoutblocks;
            print $outfile "$block->{scaffold}\t$block->{start}\tD\n";

            if ( $block->{start} < $partoutblocks[0]{start} ) {
                print $outfile "$block->{scaffold}\t$block->{start}\t". ($partoutblocks[0]{start}-1) .",0,-1\n";
            }
            for my $outblock (@partoutblocks) {
                print $outfile
                  "$outblock->{name}\t$outblock->{start}\t$outblock->{end},$outblock->{chr},$outblock->{cm}\n";
            }
            if ($block->{end} > $partoutblocks[-1]{end} ) {
                print $outfile "$block->{scaffold}\t" . ($partoutblocks[-1]{end}+1) . "\t$block->{end},0,-1\n";
            }
        }
    }
}

## LIBRARY FUNCTIONS

sub order {
    my ( $a, $b ) = @_;
    if ( $a > $b ) {
        my $tmp = $a;
        $a = $b;
        $b = $tmp;
    }
    return ( $a, $b );
}

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

# http://stackoverflow.com/questions/9122315/permutations-using-perl
sub permute (&@) {
    my $code = shift;
    my @idx  = 0 .. $#_;
    while ( $code->( @_[@idx] ) ) {
        my $p = $#idx;
        --$p while $idx[ $p - 1 ] > $idx[$p];
        my $q = $p or return;
        push @idx, reverse splice @idx, $p;
        ++$q while $idx[ $p - 1 ] > $idx[$q];
        @idx[ $p - 1, $q ] = @idx[ $q, $p - 1 ];
    }
}
