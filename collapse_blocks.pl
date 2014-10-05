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
use List::Util qw/min max/;

$OUTPUT_AUTOFLUSH = 1;

memoize "mirror";

my $samplenum = 69;
my $empty     = ' ' x $samplenum;

my %mstheader = (
    population_type              => "RIL2",
    population_name              => "HeliconiusWGS",
    distance_function            => "kosambi",
    cut_off_p_value              => "0.000001",
    no_map_dist                  => "15",
    no_map_size                  => "2",
    missing_threshold            => "0.75",
    estimation_before_clustering => "yes",
    detect_bad_data              => "yes",
    objective_function           => "ML",
);

my %sex = ( "Paternal-AHAB_AHB" => 0, "Paternal-AHAB_AHA" => 0 );

# Chroms only used for graphics and output! Patterns are inferred correctly without them
my %chroms = (
    "ABBBBBBABABAABAAAABBBBBAABABBABABAABAAABABBBABBABBABBAABBABBAABBAABBB" => "1",
    "ABBBABABBAABBBBBBABBAABBBBBBABBBABBABABAABBAAABBABBABBAAAAAAAABBAABAB" => "2",
    "ABABABAABAABBBBBBABBBAABBBBAABAABAABAABAAAAABAABBAAAAAAAABBBBABBABBBA" => "3",
    "ABBABBBBBABBABAAAAABAABBBBBBAAAAABABBAABABABAABABBAABABBAAABBBAABBABB" => "4",
    "AABABBABBBABABABBBBABBBBBAABBABABAABABBABAAAAABBBAAABBAABBBBABBBBBABA" => "5",
    "AABABBBBBABBAABABAABAABABABAAAAABAAAABAABBBBAABAABBBBAAABABAABAAAABAB" => "6",
    "AABBAABBABBBAAABBBBBBBBABAAABBAAABAAAAAAABBAAABAAAABAABBBBABABABAABBB" => "7",
    "AABBABBAABBBAAABAABABAAAAABBABBAAAABBABAABABABBBAABAABAAABABAAAAAAAAB" => "8",
    "AABBBABBABAABBABBABBBAABABBBABBABAAAAAABABABAAAAABBBBBAABBABBAABABBAA" => "9",
    "ABAAABABAAAABAABBBABABABABBABABABBABABBBBABBABBBBBBBBABAAAAABBBABBBBB" => "10",
    "ABAABBBABBAABAABBAAAAAABABBAABAAAAABBAAAAAAABABBABAAABBABAAAABABABBAB" => "11",
    "AAABBAABBBBABBBABABABBABBAAAAABAAABBABAABBBABAABAAABAAABBBBBBAABBBBAA" => "12",
    "ABAABBBBAAAABABAABBBAAAAAABBBBBBBBAAAAABAAABABABBBABABAABAAABBAABBBAA" => "13",
    "AAABABBBABBAABBABABABBBAAAABBAAAAABBAABBAAABAABBBBBAABBABBAABBAABBAAB" => "14",
    "ABBBBBABBAAAABBBBBBABBABBBAAABAABBABAAAAAAAABBABAABABBABBAAAAABAABBAA" => "15",
    "ABBAABBABBBABBBBBBAAAABABBAAABBBBAABABAABBAAABAAABBBBBBBBBABBBAAABABB" => "16",
    "ABBBAABAABBABAABABABABBAABAABBABBAABBAAAABAAABBBAAAABAABBBAAABBBAABAB" => "17",
    "ABAAAABAAAABAABAAAAABAABBAABABBAABAABBAABBBAAAABBAABABBBABBABBBAAABBB" => "18",
    "ABBAABBAABBBAABBBBBAAABBBBAABBBAABAABBBABABBAAABAABBAABBAAAAAABBAABAA" => "19",
    "ABBBABAAABBBBABAABBBBAABBABBABBABBABBBBBBBBBAAABBAAAABBBBBBBABBBBAAAA" => "20",
    "ABBABABBBBAABBBBBBBAAAAABABAAABBBBABABABBBABBABABBAAAAAABABBAABBAABAB" => "21",
);

my %args;
$args{input}  = "";
$args{output} = "test";

my $options_okay =
  GetOptions( 'input=s' => \$args{input}, 'output=s' => \$args{output} );

print STDERR "Loading blocks...\n";
my ( $header, $types, $blocklist ) = load_blocks( $args{input} );

# Print header for statistics summary
printf STDERR "%-60s", "Step";
for my $t ( 4 .. $#{$header} ) {
    printf STDERR "\t%-40s", $header->[$t];
}
print STDERR "\n";

get_block_stats( "After load_blocks", $header, $types, $blocklist );

fill_blocks($blocklist);

correct_maternal( "Maternal-AHAH", $blocklist );

get_block_stats( "After correct_maternal", $header, $types, $blocklist );

collapse( $blocklist, $types );

get_block_stats( "After collapse", $header, $types, $blocklist );

my ($matpat, $patmat, $pattern_block) = combine_maternal_and_paternal_patterns( $blocklist );

get_block_stats( "After filling maternal and paternal patterns", $header, $types, $blocklist );

print STDERR "Writing clean blocks...\n";
output_blocks( $args{input}, $blocklist );

print STDERR "Making linkage maps...\n";
clean_intercross_patterns($matpat, $patmat, $pattern_block);

my ($scfmap, $genome, $markerscf) = make_linkage_maps( $blocklist, $args{output}, $matpat, $pattern_block);

my ($scfstats) = validate_scaffolds( $blocklist, $scfmap, $genome);

print STDERR "Writing maps and stats...\n";
output_maps( $blocklist, $args{output}, $genome, $scfstats );

output_stats( $scfstats, $scfmap, $markerscf );

sub combine_maternal_and_paternal_patterns {
    my ($blocklist) = @_;
    
    my %matpat;
    my %patmat;
    my %pattern_block;
    my $sexmat;
    for my $block ( @{$blocklist} ) {
        # Find complete Maternal and Paternal patterns if they exist,
        # looking for sex patterns first
        my $mat = $block->{'Sex-HB'};
        $sexmat = $mat if $mat !~ /[ \-]/;
        my $pat = $block->{'Paternal-AHAB_AHA'};
        $pat = $block->{'Paternal-AHAB_AHB'} if ( $pat =~ /[ \-]/ );
        $mat = defined $sexmat ? $sexmat : 'S' x length $pat
          if ( $pat !~ /[ \-]/ and $mat =~ /[ \-]/ );

        $mat = $block->{"Maternal-AHAH"} if ( $mat =~ /[ \-]/ );
        $pat = $block->{"Paternal-AHAH"} if ( $pat =~ /[ \-]/ );

        # If maternal pattern is complete or absent
        # and paternal pattern is absent,
        # try to convert intercross patterns to paternal
        next if ( $mat =~ /\-/ );
        if ( $pat =~ /^( +)$/ ) {
            ( $mat, $pat ) = convert_intercross_block($block);
        }
        next if $mat =~ /[ \-]/;
        next if $pat =~ /[ \-]/;

        # If maternal and paternal patterns found, store them
        $block->{'Maternal-AHAH'} = $mat
          if ( $block->{'Maternal-AHAH'} eq $empty );
        $block->{'Paternal-AHAH'} = $pat
          if ( $block->{'Paternal-AHAH'} eq $empty );
        $matpat{$mat}{$pat}{length} += $block->{'length'};
        $matpat{$mat}{$pat}{blocks}++;
        $patmat{$pat}{$mat}{length} += $block->{'length'};
        $patmat{$pat}{$mat}{blocks}++;
        $pattern_block{$pat}{ $block->{'scaffold'} }{ $block->{'start'} } =
          $block->{'end'};
    }

    for my $mat ( keys %matpat ) {
        if ( $mat =~ /^(S+)$/ ) {
            for my $sexpat ( keys %{ $matpat{$mat} } ) {
                $patmat{$sexpat}{$sexmat}{length} +=
                  $patmat{$sexpat}{$mat}{length};
                $patmat{$sexpat}{$sexmat}{blocks} +=
                  $patmat{$sexpat}{$mat}{blocks};
                delete $patmat{$sexpat}{$mat};
                $matpat{$sexmat}{$sexpat}{length} +=
                  $matpat{$mat}{$sexpat}{length};
                $matpat{$sexmat}{$sexpat}{blocks} +=
                  $matpat{$mat}{$sexpat}{blocks};
            }
            delete $matpat{$mat};
        }
    }
    return (\%matpat, \%patmat, \%pattern_block);
    
}

sub clean_intercross_patterns {
    my ($matpat, $patmat, $pattern_block ) = @_;
    
    my @int;
    my @pat;

    for my $pat ( keys %{$pattern_block} ) {
        $pat =~ 'H' ? push @int, $pat : push @pat, $pat;
    }

    for my $int (@int) {
        my $imat         = "I" x length $int;
        my @i            = split //, $int;
        my $int_pmatched = 0;
        for my $pat (@pat) {
            if ( int_match( \@i, $pat ) ) {
                $int_pmatched++;
                my $pmat =
                  ( sort { $patmat->{$pat}{$b}{length} <=> $patmat->{$pat}{$a}{length} } keys %{ $patmat->{$pat} } )[0];
                $patmat->{$int}{$pmat}{length} = $patmat->{$int}{$imat}{length};
                $patmat->{$int}{$pmat}{blocks} = $patmat->{$int}{$imat}{blocks};
                delete $patmat->{$int}{$imat};
                $matpat->{$pmat}{$int}{length} = $matpat->{$imat}{$int}{length};
                $matpat->{$pmat}{$int}{blocks} = $matpat->{$imat}{$int}{blocks};
                delete $matpat->{$imat}{$int};
                last;
            }
        }
        if ( !$int_pmatched ) {
            my $minh     = length $int;
            my $minh_mat = $empty;
            for my $mat ( keys %{$matpat} ) {
                my $forh = int_hamming( \@i, $mat );
                my $revh = int_hamming( \@i, mirror($mat) );
                my $math = $forh < $revh ? $forh : $revh;
                if ( $math < $minh ) {
                    $minh     = $math;
                    $minh_mat = $mat;
                }
            }

            #            if ( $minh == 0 ) {
            if ( $minh_mat ne $empty ) {
                $patmat->{$int}{$minh_mat}{length} = $patmat->{$int}{$imat}{length};
                $patmat->{$int}{$minh_mat}{blocks} = $patmat->{$int}{$imat}{blocks};
                delete $patmat->{$int}{$imat};
                $matpat->{$minh_mat}{$int}{length} = $matpat->{$imat}{$int}{length};
                $matpat->{$minh_mat}{$int}{blocks} = $matpat->{$imat}{$int}{blocks};
                delete $matpat->{$imat}{$int};
            }
        }
    }

    return;
}

sub make_linkage_maps {
    my ($blocklist, $output, $matpat, $pattern_block) = @_;

    my %scfmap;
    my %genome;
    my %markerscf;
    for my $mat ( keys %{$matpat} ) {
        # Delete paternal patterns less than 10kb long or occurring in only 1 or 2 blocks
        for my $pat ( keys %{ $matpat->{$mat} } ) {
            delete $matpat->{$mat}{$pat}
              if (  $matpat->{$mat}{$pat}{length} < 1000
                and $matpat->{$mat}{$pat}{blocks} <= 2 );
        }

        # Skip maternal pattern if no paternal patterns left

        # Make map for this maternal pattern
        next if ( keys %{ $matpat->{$mat} } == 0 );
        my $markercode = run_mstmap( $mat, $matpat->{$mat}, $output );
        $genome{$mat} = load_map( $output, $mat );

        # If more than one linkage map returned for this chromosome,
        # attempt to rephase markers and remake the map
        if ( keys %{ $genome{$mat} } == 2 ) {
            my $mirlg = ( keys %{ $genome{$mat} } )[0];
            my %phased;
            for my $lg ( keys %{ $genome{$mat} } ) {
                for my $cm ( keys %{ $genome{$mat}{$lg} } ) {
                    for my $markernum ( keys %{ $genome{$mat}{$lg}{$cm} } ) {
                        my $pattern = $markercode->{$markernum}{orig};
                        my $fixed = $lg eq $mirlg ? mirror($pattern) : $pattern;
                        $phased{$fixed}++;
                    }
                }
            }
            $markercode = run_mstmap( $mat, \%phased, $output );
            $genome{$mat} = load_map( $output, $mat );
        }

        # Assign scaffolds to map position
        for my $lg ( sort keys %{ $genome{$mat} } ) {
            for my $cm ( sort { $a <=> $b } keys %{ $genome{$mat}{$lg} } ) {
                for my $marker (
                    sort { $a <=> $b }
                    keys %{ $genome{$mat}{$lg}{$cm} }
                  )
                {
                    my $pattern = $markercode->{$marker}{orig};
                    $pattern = mirror($pattern)
                      if !defined $pattern_block->{$pattern};
                    for my $scf ( sort keys %{ $pattern_block->{$pattern} } ) {
                        for my $start (
                            sort { $a <=> $b }
                            keys %{ $pattern_block->{$pattern}{$scf} }
                          )
                        {
                            my $length = $pattern_block->{$pattern}{$scf}{$start} - $start + 1;
                            $scfmap{$scf}{$start}{mat} = $mat;
                            $scfmap{$scf}{$start}{lg}  = $lg;
                            $scfmap{$scf}{$start}{cm}  = $cm;
                            $markerscf{"$mat:$lg:$cm"}{$scf}++;
                        }
                    }
                }
            }
        }
    }

    return (\%scfmap, \%genome, \%markerscf);
}

sub output_maps {
    my ( $blocklist, $output, $genome, $scfstats ) = @_;

    open my $chrommapfh, ">", "$output.chrommap.tsv"
      or croak "Can't open chromosome map!\n";
    open my $scfmapfh, ">", "$output.scfmap.tsv"
      or croak "Can't open scaffold map!\n";

    print $chrommapfh "Chromosome\tcM\tStart\tLength\n";

    print $scfmapfh "Chromosome\tScaffold\tScfStart\tScfEnd\tChrStart\tLength\tScfChroms\tScfGaps\tScfOriented\n";

    # Match inferred chromosomes to real chromosomes
    my %chromnums;
    my $newchromnum = 100;
    for my $mat (keys %{$genome}) {
        if (defined $chroms{$mat}) {
            $chromnums{$mat} = $chroms{$mat}
        }
        else {
            $chromnums{$mat} = $newchromnum;
            $newchromnum++;
        }
    }

    foreach my $mat ( sort { $chromnums{$a} <=> $chromnums{$b} } keys %{$genome} ) {
#        croak "More than one linkage group found for $chroms{$mat}!\n"
#          if ( keys %{ $genome->{$mat} } > 1 );
        foreach my $lg ( keys %{ $genome->{$mat} } ) {
            my $chrpos_cm  = 1;
            my $chrpos_scf = 1;

            my @cms = sort { $a <=> $b } keys %{ $genome->{$mat}{$lg} };
            my @chromscf;
            foreach my $cm_i ( 0 .. $#cms ) {
                print $chrommapfh "$chromnums{$mat}\t$cms[$cm_i]\t$chrpos_cm\t$genome->{$mat}{$lg}{$cms[$cm_i]}{len}\n";
                my @ordered_scfs =
                  order_scfs( $genome->{$mat}{$lg}{ $cms[$cm_i] }{scf}, $genome->{$mat}{$lg}, \@cms, $cm_i );
                foreach my $scf (@ordered_scfs) {
                    if ( defined( $chromscf[-1] )
                        and $chromscf[-1]{scf} eq $scf )
                    {
                        my $start =
                          $genome->{$mat}{$lg}{ $cms[$cm_i] }{scf}{$scf}{start};
                        my $end =
                          $genome->{$mat}{$lg}{ $cms[$cm_i] }{scf}{$scf}{end};
                        $chromscf[-1]{scfstart} = $start if ( $start < $chromscf[-1]{scfstart} );
                        $chromscf[-1]{scfend}   = $end   if ( $chromscf[-1]{scfend} < $end );
                        $chromscf[-1]{length} += $genome->{$mat}{$lg}{ $cms[$cm_i] }{scf}{$scf}{length};
                        $chromscf[-1]{oriented}++;
                    }
                    else {
                        push @chromscf,
                          {
                            scf      => $scf,
                            scfstart => $genome->{$mat}{$lg}{ $cms[$cm_i] }{scf}{$scf}{start},
                            scfend   => $genome->{$mat}{$lg}{ $cms[$cm_i] }{scf}{$scf}{end},
                            chrstart => $chrpos_scf,
                            length   => $genome->{$mat}{$lg}{ $cms[$cm_i] }{scf}{$scf}{length},
                            oriented => 0
                          };
                    }
                    $chrpos_scf +=
                      $genome->{$mat}{$lg}{ $cms[$cm_i] }{scf}{$scf}{length};
                }
                $chrpos_cm += $genome->{$mat}{$lg}{ $cms[$cm_i] }{len};
            }
            foreach my $i ( 0 .. $#chromscf ) {
                my $scfi = $chromscf[$i]{scf};
                print $scfmapfh
"$chromnums{$mat}\t$scfi\t$chromscf[$i]{scfstart}\t$chromscf[$i]{scfend}\t$chromscf[$i]{chrstart}\t$chromscf[$i]{length}";
                print $scfmapfh "\t$scfstats->{$scfi}{chromosomes}\t$scfstats->{$scfi}{gaps}\t$chromscf[$i]{oriented}";
                print $scfmapfh "\n";
            }
        }
    }

    close $scfmapfh;
    close $chrommapfh;

    #    check_unassigned($args{input}, \%scfstats, \%patmat);

}

sub output_stats {
    my ($scfstats, $scfmap, $markerscf) = @_;

    my %genomestat;
    my %local_assembly_markers;
    foreach my $scf ( sort keys %{$scfstats} ) {

        my $stat = $scfstats->{$scf};
        $stat->{ordered} = 0;
        if (    $stat->{chromosomes} == 1
            and $stat->{lgs} == 1
            and $stat->{gaps} == 0
            and !$stat->{oriented} )
        {
            my %markers;
            for my $pos ( keys %{ $scfmap->{$scf} } ) {
                $markers{"$scfmap->{$scf}{$pos}{mat}:$scfmap->{$scf}{$pos}{lg}:$scfmap->{$scf}{$pos}{cm}"}++;
            }
            croak "Should only be one marker at non-oriented scaffold $scf!"
              if ( keys %markers != 1 );
            my $marker                 = ( keys %markers )[0];
            my $unoriented_marker_scfs = 0;
            for my $markerscf ( keys %{ $markerscf->{$marker} } ) {
                next if $markerscf eq $scf;
                if (   $scfstats->{$markerscf}{chromosomes} != 1
                    or $scfstats->{$markerscf}{lgs} != 1
                    or $scfstats->{$markerscf}{gaps} > 0
                    or !$scfstats->{$markerscf}{oriented} )
                {
                    $unoriented_marker_scfs++;
                }
            }
            if ( !$unoriented_marker_scfs ) {
                $stat->{ordered}++;
            }

            $local_assembly_markers{$marker}++ if !$stat->{ordered};
        }

#        print
#"$scf\t$stat->{markerblocks}\t$stat->{allblocks}\t$stat->{length}\t$stat->{chromosomes}\t$stat->{lgs}\t$stat->{gaps}\n";
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
                    if ( $stat->{oriented} ) {
                        $genomestat{"Oriented"}{scf}++;
                        $genomestat{"Oriented"}{len} += $stat->{length};
                    }
                    if ( $stat->{ordered} ) {
                        $genomestat{"Ordered"}{scf}++;
                        $genomestat{"Ordered"}{len} += $stat->{length};
                    }
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
        printf STDERR "%16s\t%4d\t%9d\n", $stat, $genomestat{$stat}{scf}, $genomestat{$stat}{len};
        next if $stat =~ /Oriented/ or $stat =~ /Ordered/;
        $genomesize += $genomestat{$stat}{len};
        $genomescf  += $genomestat{$stat}{scf};
    }
    printf STDERR "%16s\t%4d\t%9d\n", 'Genome', $genomescf, $genomesize;

    print STDERR "Marker blocks requiring local assembly: ", scalar keys %local_assembly_markers, "\n";
    my %lam_block_scfs;
    my %lam_block_sizes;
    for my $lam ( sort keys %local_assembly_markers ) {
        my $scfnum = keys %{ $markerscf->{$lam} };
        $lam_block_scfs{$scfnum}++;
        for my $scf ( keys %{ $markerscf->{$lam} } ) {
            $lam_block_sizes{$lam} += $scfstats->{$scf}{length};
        }
    }

    #    for my $lam ( sort { $lam_block_sizes{$b} <=> $lam_block_sizes{$a} }
    #        keys %lam_block_sizes )
    #    {
    #        print "$lam\t", scalar keys %{ $markerscf{$lam} },
    #          "\t$lam_block_sizes{$lam}\n";
    #    }
    #    for my $scfnum ( sort { $a <=> $b } keys %lam_block_scfs ) {
    #        print "$scfnum\t$lam_block_scfs{$scfnum}\n";
    #    }
}

sub order_scfs {
    my ( $scf_ref, $lg_ref, $cms_ref, $cm_i ) = @_;

    my @scfs = keys %{$scf_ref};
    my @ordered_scfs;

    my $first = "";
    my $last  = "";
    map {
        $first = $_
          if (  ( $cm_i > 0 )
            and ( defined $lg_ref->{ $cms_ref->[ $cm_i - 1 ] }{scf}{$_} ) );
        $last = $_
          if (  ( $cm_i < $#{$cms_ref} )
            and ( defined $lg_ref->{ $cms_ref->[ $cm_i + 1 ] }{scf}{$_} ) );
    } @scfs;

    push @ordered_scfs, $first if $first ne "";
    map { push @ordered_scfs, $_ if ( ( $_ ne $first ) and ( $_ ne $last ) ); } @scfs;
    push @ordered_scfs, $last if $last ne "" and $first ne $last;
    return @ordered_scfs;
}

sub int_hamming {
    # Calculate hamming distance between intercross and paternal pattern
    my ( $int, $pat ) = @_;
    my $hamming = 0;
    my @p = split //, $pat;
    for my $b ( 0 .. $#{$int} ) {
        next if $int->[$b] eq 'H' or $int->[$b] eq '-';
        $hamming++ if $int->[$b] ne $p[$b];
    }
    return $hamming;
}

sub int_match {
    # Check for complete match between intercross and other pattern,
    # ignoring missing bases in intercross pattern
    my ( $inta, $pat ) = @_;
    my @p = split //, $pat;
    my $match = 1;
    for my $a ( 0 .. $#{$inta} ) {
        next if $inta->[$a] eq 'H' or $inta->[$a] eq '-';
        if ( $inta->[$a] ne $p[$a] ) { $match = 0; last; }
    }
    return $match;
}

sub check_unassigned {
    my ( $input, $scfstats, $patmat ) = @_;
    my %unassigned;
    for my $scf ( keys %{$scfstats} ) {
        $unassigned{$scf}++ if $scfstats->{$scf}{markerblocks} == 0;
    }

    open my $snps, '<', "$input.markers.out"
      or croak "Can't open SNP file! $OS_ERROR\n";
    my %scfsnps;
    my $scf;
    while ( my $snp = <$snps> ) {
        if ( $snp =~ /Reject/ ) {
            my @f = split /\t/, $snp;
            next if !( defined $unassigned{ $f[0] } );
            $scf = $f[0];
            my $pattern = $f[8];
            next if $pattern =~ /[01]/;
            $scfsnps{$pattern}++;
        }
        if ( $snp =~ /^\-/ and keys %scfsnps > 0 ) {
            foreach my $pattern (
                sort { $scfsnps{$b} <=> $scfsnps{$a} }
                keys %scfsnps
              )
            {
                next if $scfsnps{$pattern} == 1;

                #                print
                #"$scf\t$scfstats->{$scf}{length}\t$pattern\t$scfsnps{$pattern}\n";
            }
            %scfsnps = ();
        }
    }
    close $snps;
}

sub validate_scaffolds {
    my ($blocklist, $scfmap, $genome) = @_;

    my $curscf = "";
    my @scfblocks;
    my %scfstats;
    for my $block ( @{$blocklist} ) {
        if ( $block->{'scaffold'} ne $curscf ) {
            if ( $curscf ne "" ) {
                validate_scaffold( \@scfblocks, $scfmap, $genome, \%scfstats );
            }
            $curscf    = $block->{'scaffold'};
            @scfblocks = ();
        }
        push @scfblocks, $block;
    }
    validate_scaffold( \@scfblocks, $scfmap, $genome, \%scfstats );

    return \%scfstats;
}

sub validate_scaffold {
    my ( $scfblocks, $scfmap, $genome, $stats ) = @_;
    my $blocks_with_marker = 0;
    my $scf                = $scfblocks->[0]{'scaffold'};
    my %scfcms;
    my $scflen = 0;
    my %scfmarkers;
    for my $block ( @{$scfblocks} ) {
        $scflen += $block->{'length'};
        if ( defined $scfmap->{ $block->{'scaffold'} }{ $block->{'start'} } ) {
            my $mappos = $scfmap->{$scf}{ $block->{'start'} };
            $blocks_with_marker++;
            $scfcms{ $mappos->{mat} }{ $mappos->{lg} }{ $mappos->{cm} }{blocks}++;
            $scfcms{ $mappos->{mat} }{ $mappos->{lg} }{ $mappos->{cm} }{length} +=
              $block->{'length'};

            if (  !( defined $scfcms{ $mappos->{mat} }{ $mappos->{lg} }{ $mappos->{cm} }{start} )
                or ( $block->{'start'} < $scfcms{ $mappos->{mat} }{ $mappos->{lg} }{ $mappos->{cm} }{start} ) )
            {
                $scfcms{ $mappos->{mat} }{ $mappos->{lg} }{ $mappos->{cm} }{start} = $block->{'start'};
            }

            if (  !( defined $scfcms{ $mappos->{mat} }{ $mappos->{lg} }{ $mappos->{cm} }{end} )
                or ( $block->{'end'} > $scfcms{ $mappos->{mat} }{ $mappos->{lg} }{ $mappos->{cm} }{end} ) )
            {
                $scfcms{ $mappos->{mat} }{ $mappos->{lg} }{ $mappos->{cm} }{end} =
                  $block->{'end'};
            }

            $scfmarkers{"$mappos->{mat}:$mappos->{lg}:$mappos->{cm}"}{ $block->{'scaffold'} }++;
        }
    }
    my $scfchroms = keys %scfcms // 0;
    my $scflgs    = 0;
    my $scfgaps   = 0;
    $stats->{$scf}{oriented} = 0;

    for my $mat ( sort keys %scfcms ) {
        for my $lg ( sort keys %{ $scfcms{$mat} } ) {
            $scflgs++;
            my @scfcm = sort { $a <=> $b } keys %{ $scfcms{$mat}{$lg} };
            $stats->{$scf}{oriented}++ if @scfcm > 1;
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

              #                    my $chrom = $chroms{$mat} // "Unknown";
              #                    print
              #"$scf\t$chrom\t$mat\t$lg\t$lgcm\t$scfcms{$mat}{$lg}{$lgcm}{length}\t$scfcms{$mat}{$lg}{$lgcm}{blocks}\n";
                    $start_check = 1;
                    shift @scfcm;
                    $genome->{$mat}{$lg}{$lgcm}{scf}{$scf}{start} =
                      $scfcms{$mat}{$lg}{$lgcm}{start};
                    $genome->{$mat}{$lg}{$lgcm}{scf}{$scf}{end} =
                      $scfcms{$mat}{$lg}{$lgcm}{end};
                    $genome->{$mat}{$lg}{$lgcm}{scf}{$scf}{length} =
                      $scfcms{$mat}{$lg}{$lgcm}{length};
                    $genome->{$mat}{$lg}{$lgcm}{len} +=
                      $scfcms{$mat}{$lg}{$lgcm}{length};
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

sub convert_intercross_block {
    my ($block) = @_;

    for my $ic ( "Intercross-ABHABH_HHA", "Intercross-ABHABH_HHH" ) {
        next if $block->{$ic} =~ /[ \-]/;
        my $mat = $block->{'Maternal-AHAH'};
        if ( $mat =~ /^( +)$/ ) {
            return ( "I" x length( $block->{$ic} ), phase( $block->{$ic} ) );
        }
        my $pat = convert_intercross( $mat, $block->{$ic} );
        return ( $mat, $pat );
    }
    return ( " ", " " );
}

sub convert_intercross {
    my ( $mat, $int ) = @_;
    return $int if $mat =~ /I/;
    my @int = split //, phase($int);
    my @mat = split //, $mat;
    my @pat =
      map { $int[$_] eq '-' ? '-' : $int[$_] eq 'H' ? ( $mat[$_] eq 'A' ? 'B' : 'A' ) : $mat[$_] } 0 .. $#int;
    my $pat = join '', @pat;
    return $pat;
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
    map { delete $lg{$_} if keys %{ $lg{$_} } == 1 } keys %lg;
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
        "population_type", "population_name",    "distance_function", "cut_off_p_value",
        "no_map_dist",     "no_map_size",        "missing_threshold", "estimation_before_clustering",
        "detect_bad_data", "objective_function", "number_of_loci",    "number_of_individual",
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

    #    my @markerlist = keys %{$markers};
    #    my $phased_markers = phase_markers (\@markerlist, $pattern);

    for my $marker ( keys %{$markers} ) {
        my $outmarker = $marker;
        $outmarker = convert_intercross( $pattern, $marker )
          if $marker =~ /H/;
        $outmarker = check_mirror( $outmarker, $markers );
        $id = output_marker( $id, $outmarker, $marker, $mstmapin, $codein, \%marker_lookup );
    }

    close $codein;
    close $mstmapin;

    system("MSTMap.exe $output.$pattern.mstmap.markers $output.$pattern.mstmap.map > $output.$pattern.mstmap.log");

    return \%marker_lookup;
}

sub phase_markers {
    my ( $markers, $pattern ) = @_;
    my @sm =
      sort map { $_ =~ /H/ ? convert_intercross( $pattern, $_ ) : $_ } @{$markers};
    my $len    = length $sm[0];
    my $mirror = 0;
    my %phased;
    $phased{ $sm[0] } = $sm[0];
    for my $i ( 1 .. $#sm ) {
        $mirror = 1 if hamming( $sm[ $i - 1 ], $sm[$i] ) > ( $len / 2 );
        my $phase = $mirror ? mirror( $sm[$i] ) : $sm[$i];
        $phased{$phase} = $sm[$i];
    }

    return \%phased;
}

sub check_mirror {
    my ( $marker, $markers ) = @_;
    my $mirror     = mirror($marker);
    my @markerlist = keys %{$markers};
    my $markerh    = sum_hamming( $marker, \@markerlist );
    my $mirrorh    = sum_hamming( $mirror, \@markerlist );
    return $markerh < $mirrorh ? $marker : $mirror;

    #    my $mirror_minh = min_hamming( $mirror, \@markerlist );
    #    my $marker_minh = min_hamming( $marker, \@markerlist );
    #    return $mirror_minh < $marker_minh ? $mirror : $marker;
}

sub min_hamming {
    my ( $marker, $list ) = @_;
    my $minh = length $marker;
    for my $l ( @{$list} ) {
        next if $l eq $marker;
        my $h = hamming( $l, $marker );
        $minh = $h if $h < $minh;
    }
    return $minh;
}

sub sum_hamming {
    my ( $marker, $list ) = @_;
    my $h = 0;
    for my $l ( @{$list} ) {
        next if $l eq $marker;
        $h += hamming( $l, $marker );
    }
    return $h;
}

sub output_marker {
    my ( $id, $outmarker, $origmarker, $mstmapin, $codein, $marker_lookup ) = @_;

    print $codein "$id\t$outmarker\t$origmarker\n";
    print $mstmapin "$id";
    my @gt = split //, $outmarker;
    map { print $mstmapin "\t"; print $mstmapin $_ eq 'H' ? 'X' : $_; } @gt;
    print $mstmapin "\n";
    $marker_lookup->{$id}{out}  = $outmarker;
    $marker_lookup->{$id}{orig} = $origmarker;
    $id++;
    return $id;
}

sub collapse {
    my ( $blocklist, $types ) = @_;
    my $i = 0;
    my $j = 1;
    while ( $j <= $#{$blocklist} ) {
        my $same = 1;
        $same = 0
          if ( $blocklist->[$i]{'scaffold'} ne $blocklist->[$j]{'scaffold'} );
        map {
            if ( $blocklist->[$i]{$_} ne $blocklist->[$j]{$_} ) {
                $same = 0;
            }
        } keys %{$types};

        if ($same) {
            $blocklist->[$i]{'end'}    = $blocklist->[$j]{'end'};
            $blocklist->[$i]{'length'} = $blocklist->[$i]{'length'} + $blocklist->[$j]{'length'};
            splice( @{$blocklist}, $j, 1 );
        }
        else {
            $i = $j;
            $j++;
        }
    }
}

sub correct_maternal {
    my ( $mattype, $blocklist ) = @_;

    fix_maternal_errors( $mattype, $blocklist );

    # Extract maternal patterns and calculate lengths covered
    my %mat;
    for my $i ( 0 .. $#{$blocklist} ) {
        $mat{ $blocklist->[$i]{$mattype} }{length} +=
          $blocklist->[$i]{'length'};
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
            $fix = 1 if hamming( $mata,         $matb ) <= 6;
            $fix = 1 if hamming( mirror($mata), $matb ) <= 6;

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

    # Extract just maternal patterns from blocks and index by scaffold and start position
    my %mat;
    for my $i ( 0 .. $#{$blocklist} ) {
        next if $blocklist->[$i]{$mattype} eq $empty;
        $mat{ $blocklist->[$i]{'scaffold'} }{ $blocklist->[$i]{'start'} }{pattern} = $blocklist->[$i]{$mattype};
        $mat{ $blocklist->[$i]{'scaffold'} }{ $blocklist->[$i]{'start'} }{block}   = $i;
    }

    # Correct changes in maternal patterns
    for my $scf ( keys %mat ) {
        next if keys %{ $mat{$scf} } < 3;
        my @starts = sort { $a <=> $b } keys %{ $mat{$scf} };
        for my $i ( 1 .. $#starts - 1 ) {
            next
              if $mat{$scf}{ $starts[ $i - 1 ] }{pattern} eq $mat{$scf}{ $starts[$i] }{pattern};
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

sub fill_blocks {
    my ($blocklist) = @_;
    fill_type_pair( "Intercross-ABHABH_HHA", "Intercross-ABHABH_HHH", $blocklist );
    fill_type_pair( "Intercross-ABHABH_HHH", "Intercross-ABHABH_HHA", $blocklist );
    fill_type_pair( "Intercross-ABHABH_HHA", "Intercross-ABHABH_HHH", $blocklist );
    fill_type_pair( "Intercross-ABHABH_HHH", "Intercross-ABHABH_HHA", $blocklist );
    fill_type_pair( "Intercross-ABHABH_HHA", "Paternal-AHAH",         $blocklist );
    fill_type_pair( "Intercross-ABHABH_HHH", "Paternal-AHAH",         $blocklist );
    fill_type_pair( "Paternal-AHAH",         "Intercross-ABHABH_HHA", $blocklist );
    fill_type_pair( "Paternal-AHAH",         "Intercross-ABHABH_HHH", $blocklist );
    fill_type_pair( "Maternal-ABHABH",       "Maternal-AHAH",         $blocklist );
    fill_type_pair( "Intercross-ABHABH_HHA", "Maternal-AHAH",         $blocklist );
    fill_type_pair( "Intercross-ABHABH_HHH", "Maternal-AHAH",         $blocklist );
    fill_type_pair( "Paternal-AHAH",         "Maternal-AHAH",         $blocklist );
    fill_type_pair( "Paternal-AHAB_AHB",     "Paternal-AHAB_AHA",     $blocklist );
    fill_type_pair( "Paternal-AHAB_AHA",     "Paternal-AHAB_AHB",     $blocklist );
    fill_type_pair( "Paternal-AHAB_AHB",     "Sex-HB",                $blocklist );
    fill_type_pair( "Paternal-AHAB_AHA",     "Sex-HB",                $blocklist );
}

sub fill_type_pair {
    my ( $typea, $typeb, $blocklist ) = @_;

    my %ab;

    # Pull out type a and b fields from all blocks and calculate the total length covered by these types
    for my $b ( 0 .. $#{$blocklist} ) {
        $ab{ $blocklist->[$b]{$typea} }{ $blocklist->[$b]{$typeb} }{length} +=
          $blocklist->[$b]{'length'};
        push @{ $ab{ $blocklist->[$b]{$typea} }{ $blocklist->[$b]{$typeb} }{blocks} }, $b;
    }

    for my $apat ( keys %ab ) {

        next if $apat eq $empty;
        next if $apat =~ /\-/;
        next if keys %{ $ab{$apat} } == 1;

        # Find type b pattern most found with this type a pattern
        my $maxbpat    = "";
        my $maxblength = 0;
        for my $bpat ( keys %{ $ab{$apat} } ) {
            if ( $ab{$apat}{$bpat}{length} > $maxblength ) {
                $maxbpat    = $bpat;
                $maxblength = $ab{$apat}{$bpat}{length};
            }
        }
        next if $maxbpat eq $empty;
        next if $maxbpat =~ /\-/;

        # Correct all b patterns associated with this a pattern
        # to the most often found b pattern,
        # assuming the distance between patterns is less than 6bp
        for my $bpat ( keys %{ $ab{$apat} } ) {
            next if $bpat eq $maxbpat;

            my $fix = 0;
            $fix = 1 if $bpat eq $empty;
            $fix = 1 if hamming( $maxbpat, $bpat ) <= 6;
            $fix = 1 if hamming( mirror($maxbpat), $bpat ) <= 6;

            if ($fix) {
                for my $b ( @{ $ab{$apat}{$bpat}{blocks} } ) {
                    $blocklist->[$b]{$typeb} = $maxbpat;
                }
            }
        }
    }

    get_block_stats( "After fill $typea\-\>$typeb", $header, $types, $blocklist );
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
    return $hamming;
}

sub get_header_types {
    my $dbh = shift;

    my $sth = $dbh->prepare("SELECT * FROM blocks WHERE scaffold = 'getcolumnnames'");
    $sth->execute;
    my @header = @{ $sth->{NAME} };

    my %types;
    map { $types{$_} = "" } @header;
    delete $types{"scaffold"};
    delete $types{"start"};
    delete $types{"end"};
    delete $types{"length"};

    ( \@header, \%types );
}

sub load_blocks {
    my $input = shift;

    my $blocklist = [];
    my $header;
    my $types;

    my @inputfiles = split ',', $input;
    for my $inputfile (@inputfiles) {
        my $dbh = DBI->connect( "dbi:SQLite:dbname=$inputfile", "", "" );
        
        ( $header, $types ) = get_header_types($dbh);
        
        my $sth = $dbh->prepare_cached("SELECT * FROM blocks");
        my $fileblocklist = $dbh->selectall_arrayref( $sth, { Slice => {} } );
        push @{$blocklist}, @{$fileblocklist};
        $sth->finish;

        $dbh->disconnect;
    }
    
    # Make empty pattern
    my $empty;
    for my $block ( @{$blocklist} ) {
        if ( $block->{"Maternal-AHAH"} ne "" ) {
            $empty = ' ' x length $block->{"Maternal-AHAH"};
            last;
        }
    }
    
    for my $block ( @{$blocklist} ) {
        for my $type ( keys %{$block} ) {
            if ( $block->{$type} eq "" ) {
                $block->{$type} = $empty;
                next;
            }
            $block->{$type} = convert_sex( $block->{$type} ) if defined $sex{$type};
            $block->{$type} = phase( $block->{$type} )       if defined $types->{$type};
        }
    }
    return ( $header, $types, $blocklist );
}

sub output_blocks {
    my $input     = shift;
    my $blocklist = shift;

    my $inputfile = (split ',',$input)[0];
    my $dbh = DBI->connect( "dbi:SQLite:dbname=$inputfile", "", "", { AutoCommit => 0 } );

    my $sth = $dbh->prepare("DROP TABLE IF EXISTS cleanblocks");
    $sth->execute;
    $dbh->commit;

    my ( $header, $types ) = get_header_types($dbh);

    my $columns   = $#{$header};
    my $statement = "CREATE TABLE cleanblocks (scaffold text, start integer, end integer, length integer";
    map { $statement .= ", \"$_\" text" } @{$header}[ 4 .. $columns ];
    $statement .= ")";
    $sth = $dbh->prepare($statement);
    $sth->execute;
    $dbh->commit;

    $statement = "INSERT INTO cleanblocks VALUES (?,?,?,?";
    map { $statement .= ",?" } @{$header}[ 4 .. $columns ];
    $statement .= ")";
    my $insert_handle = $dbh->prepare_cached($statement);

    for my $block ( @{$blocklist} ) {
        my @values = map { $block->{$_} =~ /^( +)$/ ? '' : $block->{$_} } @{$header};
        $insert_handle->execute(@values);
    }
    $dbh->commit;
    $dbh->disconnect;
}

sub convert_sex {
    my ($pat) = @_;
    $pat =~ s/B/H/g;
    return $pat;
}

sub get_block_stats {
    my ( $step, $header, $types, $blocklist ) = @_;
    my %typeblocks;
    my %typepatterns;
    my %typebases;
    for my $block ( @{$blocklist} ) {
        for my $type ( keys %{$types} ) {
            next if $block->{$type} eq $empty;
            if ( $block->{$type} =~ /\-/ ) {
                $typeblocks{error}{$type}++;
                $typepatterns{error}{$type}{ $block->{$type} }++;
                $typebases{error}{$type} += $block->{'length'};
            }
            else {
                $typeblocks{ok}{$type}++;
                $typepatterns{ok}{$type}{ $block->{$type} }++;
                $typebases{ok}{$type} += $block->{'length'};
            }
        }
    }
    printf STDERR "%-60s", $step;
    for my $t ( 4 .. $#{$header} ) {
        my $typestats = "";
        my $type      = $header->[$t];
        $typestats .= $typeblocks{ok}{$type}    // '0';
        $typestats .= ':';
        $typestats .= $typeblocks{error}{$type} // '0';
        $typestats .= ';';
        $typestats .=
          defined $typepatterns{ok}{$type}
          ? scalar keys %{ $typepatterns{ok}{$type} }
          : '0';
        $typestats .= ':';
        $typestats .=
          defined $typepatterns{error}{$type}
          ? scalar keys %{ $typepatterns{error}{$type} }
          : '0';
        $typestats .= ';';
        $typestats .= $typebases{ok}{$type} // '0';
        $typestats .= ':';
        $typestats .= $typebases{error}{$type} // '0';

        printf STDERR "\t%-40s", $typestats;
    }
    print STDERR "\n";
}

