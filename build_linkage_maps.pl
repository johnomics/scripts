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

# Chroms only used for graphics and output! Patterns are inferred correctly without them
my %chromosomes = (
    "ABBBBBBABABAABAAAABBBBBAABABBABABAABAAABABBBABBABBABBAABBABBAABBAABBB" => "1",
    "ABBBABABBAABBBBBBABBAABBBBBBABBBABBABABAABBAAABBABBABBAAAAAAAABBAABAB" => "2",
    "ABABABAABAABBBBBBABBBAABBBBAABAABAABAABAAAAABAABBAAAAAAAABBBBABBABBBA" => "3",
    "ABBABBBBBABBABAAAAABAABBBBBBAAAAABABBAABABABAABABBAABABBAAABBBAABBABB" => "4",
    "AABABBABBBABABABBBBABBBBBAABBABABAABABBABAAAAABBBAAABBAABBBBABBBBBABA" => "5",
    "AABABBBBBABBAABABAABAABABABAAAAABAAAABAABBBBAABAABBBBAAABABAABAAAABAB" => "6",
    "AABBAABBABBBAAABBBBBBBBABAAABBAAABAAAAAAABBAAABAAAABAABBBBABABABAABBB" => "7",
    "AABBABBAABBBAAABAABABAAAAABBABBAAAABBABAABABABBBAABAABAAABABAAAAAAAAB" => "8",
    "AABBBABBABAABBABBABBBAABABBBABBABAAAAAABABABAAAAABBBBBAABBABBAABABBAA" => "9",
    "BAHHBABHBBBHABBAABBABHBABAABABABAHBABAAAABAABAAAAAAAABHBBHBBHAABHBAAA" => "10",
    "ABAABBBABBAABAABBAAAAAABABBAABAAAAABBAAAAAAABABBABAAABBABAAAABABABBAB" => "11",
    "AAABBAABBBBABBBABABABBABBAAAAABAAABBABAABBBABAABAAABAAABBBBBBAABBBBAA" => "12",
    "ABAABBBBAAAABABAABBBAAAAAABBBBBBBBAAAAABAAABABABBBABABAABAAABBAABBBAA" => "13",
    "AAABAHHBABHBAHHABBHHBBHAAAABBABBAHBBAAHBAAHBAHBBBHBAABHABBHBBBAABHAHB" => "14",
    "ABBBBBABBAAAABBBBBBABBABBBAAABAABBABAAAAAAAABBABAABABBABBAAAAABAABBAA" => "15",
    "ABBAABBABBBABBBBBBAAAABABBAAABBBBAABABAABBAAABAAABBBBBBBBBABBBAAABABB" => "16",
    "ABBBAABAABBABAABABABABBAABAABBABBAABBAAAABAAABBBAAAABAABBBAAABBBAABAB" => "17",
    "ABAAAABAAAABAABAAAAABAABBAABABBAABAABBAABBBAAAABBAABABBBABBABBBAAABBB" => "18",
    "ABBAABBAABBBAABBBBBAAABBBBAABBBAABAABBBABABBAAABAABBAABBAAAAAABBAABAA" => "19",
    "ABBBABAAABBBBABAABBBBAABBABBABBABBABBBBBBBBBAAABBAAAABBBBBBBABBBBAAAA" => "20",
    "ABBABABBBBAABBBBBBBAAAAABABAAABBBBABABABBBABBABABBAAAAAABABBAABBAABAB" => "21",
);

my %args;
$args{input}   = "";
$args{output}  = "test";
$args{verbose} = "";

my $options_okay = GetOptions(
    'input=s'  => \$args{input},
    'output=s' => \$args{output},
    'verbose'  => \$args{verbose}
);

my $metadata = { verbose => $args{verbose} };

print STDERR "INFER MARKERS\n";
my ( $blocklist, $linkage_groups, $marker_blocks ) =
  load_markers( $args{input}, $metadata );

print STDERR "MAKING LINKAGE MAPS\n";
make_linkage_maps( $blocklist, $args{output}, $linkage_groups, $marker_blocks );

output_marker_blocks( $linkage_groups, $marker_blocks, $args{output} );

print STDERR "Done\n";

exit;

## INFER MARKERS
sub load_markers {
    my ( $input, $metadata ) = @_;

    print STDERR "Loading blocks...\n";
    my $blocklist = load_blocks( $input, $metadata );

    correct_paternal( $blocklist, $metadata );

    collapse( $blocklist, $metadata );

    #    fill_empty_blocks( $blocklist, $valid_patterns, $metadata );

    my ( $linkage_groups, $marker_blocks ) = build_linkage_groups( $blocklist, $metadata );

    ( $blocklist, $linkage_groups, $marker_blocks );
}

sub load_blocks {
    my ( $input, $metadata ) = @_;

    my $blocklist = [];
    my $header;
    my $types;

    print STDERR "Load cleanblocks from database...\n" if $metadata->{verbose};
    my $dbh = DBI->connect( "dbi:SQLite:dbname=$input", "", "" );

    ( $header, $types ) = get_header_types($dbh);

    my $sth = $dbh->prepare("SELECT * FROM cleanblocks ORDER BY scaffold, start");
    my $fileblocklist = $dbh->selectall_arrayref( $sth, { Slice => {} } );
    push @{$blocklist}, @{$fileblocklist};
    $sth->finish;

    $dbh->disconnect;

    # Make empty pattern by checking first Maternal pattern (could be any pattern)
    my $empty;
    my $samplenum;
    for my $block ( @{$blocklist} ) {
        if ( $block->{Maternal} ne "" ) {
            $samplenum = length $block->{Maternal};
            $empty     = ' ' x $samplenum;
            last;
        }
    }

    $metadata->{input}   = $input;
    $metadata->{header}  = $header;
    $metadata->{samples} = $samplenum;
    $metadata->{types}   = $types;
    $metadata->{empty}   = $empty;

    get_block_stats( "After loading blocks", $blocklist, $metadata );

    $blocklist;
}

sub collapse {
    my ( $blocklist, $metadata ) = @_;
    my $i = 0;
    my $j = 1;
    while ( $j <= $#{$blocklist} ) {
        my $same = 1;
        $same = 0
          if ( $blocklist->[$i]{'scaffold'} ne $blocklist->[$j]{'scaffold'} );
        $same = 0
          if ( $blocklist->[$i]{'validity'} ne $blocklist->[$j]{'validity'} );

        map {
            if ( $blocklist->[$i]{$_} ne $blocklist->[$j]{$_} ) {
                my $phase_j = phase( $blocklist->[$j]{$_} );
                if ( $blocklist->[$i]{$_} eq $phase_j ) {
                    $blocklist->[$j]{$_} = $phase_j;
                }
                else {
                    $same = 0;
                }
            }
        } keys %{ $metadata->{types} };

        if ($same) {
            $blocklist->[$i]{'end'} = $blocklist->[$j]{'end'};
            $blocklist->[$i]{'length'} =
              $blocklist->[$i]{'length'} + $blocklist->[$j]{'length'};
            splice( @{$blocklist}, $j, 1 );
        }
        else {
            $i = $j;
            $j++;
        }
    }

    update_block_stats( "After collapse", $blocklist, $metadata );

    return;
}

sub correct_paternal {
    my ( $blocklist, $metadata ) = @_;

    my $patterns = get_paternal_patterns($blocklist);

    my $pattern_num = scalar keys %{$patterns};
    for my $pattern ( sort { $patterns->{$b}{length} <=> $patterns->{$a}{length} } keys %{$patterns} ) {

        #        if ( $patterns->{$pattern}{length} > 50000 ) {
        for my $upgrade (
            sort { $patterns->{$pattern}{validity}{$b}{length} <=> $patterns->{$pattern}{validity}{$a}{length} }
            keys %{ $patterns->{$pattern}{validity} }
          )
        {
            next if $upgrade =~ /I[ed]/ or $upgrade =~ /[X\-]/;
            for my $b ( @{ $patterns->{$pattern}{validity}{$upgrade}{blocks} } ) {
                $blocklist->[$b]{validity} = "   MVPVV";
            }
        }

        #        }
    }

    update_block_stats( "After correcting paternal", $blocklist, $metadata );

    return;

}

sub get_paternal_patterns {
    my ($blocklist) = @_;
    my %patterns;
    for my $b ( 0 .. $#{$blocklist} ) {
        my $block    = $blocklist->[$b];
        my $paternal = $block->{'Paternal'};
        $patterns{$paternal}{length} += $block->{length};
        $patterns{$paternal}{blocks}++;
        $patterns{$paternal}{validity}{ $block->{validity} }{length} +=
          $block->{length};
        push @{ $patterns{$paternal}{validity}{ $block->{validity} }{blocks} }, $b;
        $patterns{$paternal}{maternal}{ $block->{'Maternal'} }{length} += $block->{length};
        $patterns{$paternal}{maternal}{ $block->{'Maternal'} }{blocks}++;
    }

    # Collapse mirrors
    for my $paternal ( sort { $patterns{$a}{length} <=> $patterns{$b}{length} } keys %patterns ) {
        next if !defined $patterns{$paternal};
        my $mirror = mirror($paternal);
        if ( defined $patterns{$mirror} ) {
            $patterns{$mirror}{length} += $patterns{$paternal}{length};
            $patterns{$mirror}{blocks} += $patterns{$paternal}{blocks};
            for my $validity ( keys %{ $patterns{$paternal}{validity} } ) {
                $patterns{$mirror}{validity}{$validity}{length} += $patterns{$paternal}{validity}{$validity}{length};
                for my $b ( @{ $patterns{$paternal}{validity}{$validity}{blocks} } ) {
                    $blocklist->[$b]{'Paternal'} = $mirror;
                }
                push @{ $patterns{$mirror}{validity}{$validity}{blocks} },
                  @{ $patterns{$paternal}{validity}{$validity}{blocks} };
            }
            for my $maternal ( keys %{ $patterns{$paternal}{maternal} } ) {
                $patterns{$mirror}{maternal}{$maternal}{length} += $patterns{$paternal}{maternal}{$maternal}{length};
                $patterns{$mirror}{maternal}{$maternal}{blocks} += $patterns{$paternal}{maternal}{$maternal}{blocks};
            }
            delete $patterns{$mirror};
        }
    }

    \%patterns;
}

sub fill_empty_blocks {
    my ( $blocklist, $valid_patterns, $metadata ) = @_;

    my @empty_blocks;
    for my $block ( @{$blocklist} ) {
        next if $block->{'validity'} ne 'i  m p  ';
        push @empty_blocks, $block;
    }

    my $dbh = DBI->connect( "dbi:SQLite:dbname=$metadata->{input}", "", "" );

    for my $block ( sort { $b->{length} <=> $a->{length} } @empty_blocks ) {
        my $sth = $dbh->prepare(
"SELECT position, pattern FROM markers where scaffold=\"$block->{scaffold}\" and position >= $block->{start} and position <= $block->{end} and marker_type=\"Reject\" ORDER BY scaffold, position"
        );
        my $rejects = $dbh->selectall_arrayref( $sth, { Slice => {} } );
        $sth->finish;

        my $snps = @{$rejects};
        printf "%-16s\t%8d\t%8d\t%8d\t%6d\n", $block->{scaffold}, $block->{start}, $block->{end},
          $block->{length}, $snps;

        my %patterns;
        for my $reject ( @{$rejects} ) {

            #			next if $reject->{pattern} =~ /[01]/ or $reject->{pattern} =~ /^\.+$/;
            $patterns{ $reject->{pattern} }++;
        }

        my @sorted_patterns = sort { $patterns{$b} <=> $patterns{$a} } keys %patterns;
        for my $p ( 0 .. $#sorted_patterns ) {
            last if $p == 10;
            print "\t$sorted_patterns[$p]\t$patterns{$sorted_patterns[$p]}\n";
        }
    }

    $dbh->disconnect;

}

sub build_linkage_groups {
    my ( $blocklist, $metadata ) = @_;

    my %linkage_groups;
    my %marker_blocks;

    for my $b ( 0 .. $#{$blocklist} ) {
        my $block = $blocklist->[$b];

        next if $block->{validity} !~ /V$/;

        my $maternal = $block->{Maternal};
        my $paternal = $block->{Paternal};

        next
          if ( $maternal eq $metadata->{empty}
            or $paternal eq $metadata->{empty} );

        $linkage_groups{$maternal}{$paternal}{length} += $block->{length};
        $linkage_groups{$maternal}{$paternal}{blocks}++;
        $linkage_groups{$maternal}{$paternal}{output} = $paternal;
        push @{ $linkage_groups{$maternal}{$paternal}{validity}{ $block->{validity} } }, $b;
        $marker_blocks{$paternal}{ $block->{scaffold} }{ $block->{start} } =
          $block->{end};
    }

    ( \%linkage_groups, \%marker_blocks );
}

## MAKE LINKAGE MAPS

sub make_linkage_maps {
    my ( $blocklist, $output, $linkage_groups, $marker_block ) = @_;

    my %genome;

    my $pm = new Parallel::ForkManager( scalar keys %{$linkage_groups} );

    for my $chromosome_print ( sort { $chromosomes{$a} <=> $chromosomes{$b} } keys %{$linkage_groups} ) {
        my $chromosome = $chromosomes{$chromosome_print} // '-';
        $pm->start() and next;
        my $lg         = $linkage_groups->{$chromosome_print};
        my $marker_num = keys %{ $linkage_groups->{$chromosome_print} };
        printf STDERR
          "Building map for print $chromosome_print, chromosome %3s.\t%4d paternal markers to process\n",
          $chromosome, $marker_num;

        my %revised_lg;
        my $marker_count;
        my %rejects;
        for my $marker ( sort { $lg->{$b}{length} <=> $lg->{$a}{length} } keys %{$lg} ) {
            $marker_count++;
            printf STDERR "%2d: processed %3d of %3d markers, retained %3d\n", $chromosome, $marker_count, $marker_num,
              scalar keys %revised_lg
              if $marker_count % 10 == 0;
            my ( $fixed_marker, $outcome ) =
              integrate_marker( $marker, $lg->{$marker}{length}, \%revised_lg, $chromosome_print, $output );
            $lg->{$marker}{output} = $fixed_marker if $outcome eq 'Valid';
            $rejects{$outcome}{$marker} = $fixed_marker;
        }
        printf STDERR "%2d: ", $chromosome;
        for my $outcome (sort keys %rejects) {
            my $marker_num = keys %{$rejects{$outcome}};
            my $bases;
            for my $marker (keys %{$rejects{$outcome}}) {
                $bases += $lg->{$marker}{length};
                delete $lg->{$marker} if $outcome ne 'Valid';
            }
            printf STDERR "%8s:%3d:%8d ", $outcome, $marker_num, $bases;
        }
        print STDERR "\n";
        $genome{$chromosome_print} = make_linkage_map( $chromosome_print, $lg, $output );

        $pm->finish();
    }
    $pm->wait_all_children;

    \%genome;
}

sub integrate_marker {
    my ( $marker, $length, $rlg, $chromosome_print, $output ) = @_;

    $rlg->{$marker}{output} = $marker;
    $rlg->{$marker}{length} = $length;
    return ( $marker, 'Valid' ) if keys %{$rlg} < 3 or $length >= 100000;

    # Build map with marker in it
    my $marker_id = run_mstmap( $chromosome_print, $rlg, $output );
    my $map = load_map( $output, $chromosome_print );

    my $chromosome_lgs = keys %{$map};
    if ( $chromosome_lgs >= 2 ) {
        phase_linkage_groups( $map, $rlg, $marker_id );
        $marker_id = run_mstmap( $chromosome_print, $rlg, $output );
        $map = load_map( $output, $chromosome_print );
        my $chromosome_lgs = keys %{$map};
    }

    my ( $fixed_marker, $outcome ) = check_marker_on_map( $marker, $map, $rlg, $marker_id );

    delete $rlg->{$marker} if $outcome eq 'End' or $outcome eq 'Disorder';

    if ( $outcome eq 'Valid' and $fixed_marker ne $marker ) {
        if ( !defined $rlg->{$fixed_marker} ) {
            $rlg->{$fixed_marker}{output} = $fixed_marker;
            $rlg->{$fixed_marker}{length} = $rlg->{$marker}{length};
        }
        else {
            $rlg->{$fixed_marker}{length} += $rlg->{$marker}{length};
        }
        delete $rlg->{$marker};
    }
    else {
    }

    ( $fixed_marker, $outcome );
}

sub check_marker_on_map {
    my ( $marker, $map, $markers, $marker_id ) = @_;

    my $original_marker;
    my $new_marker = '';
    my $marker_status = '';
#    print "\n";
    for my $lg ( sort keys %{$map} ) {
        my @cMs = sort { $a <=> $b } keys %{ $map->{$lg} };
        for my $i ( 0 .. $#cMs ) {
            my $new = 0;
            my $cM           = $cMs[$i];
            my $cM_marker_id = ( keys %{ $map->{$lg}{$cM} } )[0];
            my $cM_marker    = $marker_id->{$cM_marker_id};
#            print "\n$lg\t$cM\t$cM_marker\t$markers->{$cM_marker}{output}\t$markers->{$cM_marker}{length}";

            if ($cM_marker eq $marker) {
                $new++;
#                print "\tNew";
            }

            if ( $i == 0 or $i == $#cMs ) {
#                print "\tEnd";
                $marker_status = 'End' if $new;
                next;
            }

            $original_marker = $markers->{$cM_marker}{output};

            my $this = $cM_marker;

            my $prev = $marker_id->{ ( keys %{ $map->{$lg}{ $cMs[ $i - 1 ] } } )[0] };
            my $next = $marker_id->{ ( keys %{ $map->{$lg}{ $cMs[ $i + 1 ] } } )[0] };

            my @prev = split //, $markers->{$prev}{output};
            my @this = split //, $markers->{$this}{output};
            my @next = split //, $markers->{$next}{output};

            my $fixed_marker = '';
            for my $i ( 0 .. $#this ) {
                if ( $prev[$i] eq $next[$i] and $prev[$i] ne $this[$i] ) {
                    $this[$i] = $prev[$i];
                }
            }
            $fixed_marker = join '', @this;
            $new_marker = $fixed_marker if $new;
            
            if ($fixed_marker ne $original_marker) {
                if ($new) {
#                    print "\tFixed"
                }
                else {
#                    print "\tDisordered";
                    $marker_status = 'Disorder';
                }
            }
        }
    }
    $marker_status = $marker_status ne '' ? $marker_status : $new_marker eq '' ? 'Missing' : 'Valid';
 #   print "\t$marker_status";
    ( $new_marker, $marker_status );
}

sub make_linkage_map {
    my ( $chromosome_print, $markers, $output ) = @_;

    # Make map for this maternal pattern
    my $marker_id = run_mstmap( $chromosome_print, $markers, $output );

    my $map = load_map( $output, $chromosome_print );
    my $chromosome_lgs = keys %{$map};

#    print STDERR "\tBuilt $chromosome_lgs linkage groups";

    # If more than one linkage map returned for this chromosome,
    # attempt to rephase markers and remake the map
    my $phased_markers;
    if ( $chromosome_lgs >= 2 ) {
        phase_linkage_groups( $map, $markers, $marker_id );
        $marker_id = run_mstmap( $chromosome_print, $markers, $output );
        $map = load_map( $output, $chromosome_print );
        my $chromosome_lgs = keys %{$map};

#        print STDERR "\tAfter phasing, $chromosome_lgs linkage groups";
    }

#    smooth_markers( $map, $markers, $marker_id );
#    $marker_id = run_mstmap( $chromosome_print, $markers, $output );
#    $map = load_map( $output, $chromosome_print );

    my $map_markers = write_map_markers( $map, $markers, $marker_id, $output, $chromosome_print );

#    print STDERR "\n";

    $map_markers;
}

sub run_mstmap {

    my ( $chromosome_print, $markers, $output ) = @_;

    my $mstmap_input = write_mstmap_header( $output, $chromosome_print, $markers );

    my $id = 1;
    my %marker_id;
    for my $marker ( keys %{$markers} ) {
        $marker_id{$id} = $marker;
        $markers->{$marker}{output} =
          check_mirror( $markers->{$marker}{output}, $markers );
        $id = write_marker( $id, $markers->{$marker}, $mstmap_input );
    }

    close $mstmap_input;

    system(
"MSTMap.exe $output.$chromosome_print.mstmap.input $output.$chromosome_print.mstmap.map > $output.$chromosome_print.mstmap.log"
    );

    \%marker_id;
}

sub write_mstmap_header {
    my ( $output, $chromosome_print, $markers ) = @_;

    my %mst_header = (
        population_type              => "RIL2",
        population_name              => "HeliconiusWGS",
        distance_function            => "kosambi",
        cut_off_p_value              => "0.000001",
        no_map_dist                  => "0",
        no_map_size                  => "0",
        missing_threshold            => "1",
        estimation_before_clustering => "yes",
        detect_bad_data              => "yes",
        objective_function           => "ML",
    );

    my @mst_header = (
        "population_type", "population_name",    "distance_function", "cut_off_p_value",
        "no_map_dist",     "no_map_size",        "missing_threshold", "estimation_before_clustering",
        "detect_bad_data", "objective_function", "number_of_loci",    "number_of_individual",
    );

    open my $mstmap_input, ">", "$output.$chromosome_print.mstmap.input"
      or croak "Can't open $output.$chromosome_print.mstmap.input: $OS_ERROR\n";

    my $samples = split //, ( keys %{$markers} )[0];
    $mst_header{"number_of_loci"}       = keys %{$markers};
    $mst_header{"number_of_individual"} = $samples;
    map { print $mstmap_input "$_ $mst_header{$_}\n"; } @mst_header;

    print $mstmap_input "locus_name";
    map { print $mstmap_input "\t$_" } ( 1 .. $samples );
    print $mstmap_input "\n";

    $mstmap_input;
}

sub write_marker {
    my ( $id, $marker, $mstmap_input, $marker_lookup ) = @_;

    print $mstmap_input "$id";
    my @gt = split //, $marker->{output};
    map { print $mstmap_input "\t"; print $mstmap_input $_ eq 'H' ? 'X' : $_; } @gt;
    print $mstmap_input "\n";
    $id++;

    $id;
}

sub load_map {
    my ( $output, $chromosome_print ) = @_;

    open my $mstmap_output, '<', "$output.$chromosome_print.mstmap.map"
      or croak "Can't open map for $chromosome_print! $OS_ERROR\n";

    my %lg;
    my $lg_num;
    my $in_group = 0;
    while ( my $mst_line = <$mstmap_output> ) {
        chomp $mst_line;

        $lg_num = $1 if ( $mst_line =~ /^group (.+)$/ );
        if ( $mst_line eq ";BEGINOFGROUP" ) {
            $in_group = 1;
            next;
        }
        $in_group = 0 if ( $mst_line eq ";ENDOFGROUP" );
        if ($in_group) {
            if ( $mst_line =~ /^(.+)\t([\d\.]+)$/ ) {
                my $id = $1;
                my $cM = $2;
                $lg{$lg_num}{$cM}{$id}++;
            }
        }
    }

    # Delete linkage groups with a single marker
    map { delete $lg{$_} if keys %{ $lg{$_} } == 1 } keys %lg;
    close $mstmap_output;

    \%lg;
}

sub check_mirror {
    my ( $marker, $markers ) = @_;
    my $mirror     = mirror($marker);
    my @markerlist = map { $markers->{$_}{output} } keys %{$markers};
    my $markerh    = sum_hamming( $marker, \@markerlist );
    my $mirrorh    = sum_hamming( $mirror, \@markerlist );

    $markerh < $mirrorh ? $marker : $mirror;
}

sub sum_hamming {
    my ( $marker, $list ) = @_;
    my $h = 0;
    for my $l ( @{$list} ) {
        next if $l eq $marker;
        $h += hamming( $l, $marker );
    }

    $h;
}

sub phase_linkage_groups {
    my ( $map, $markers, $marker_id ) = @_;

    my $first_lg = ( keys %{$map} )[0];
    for my $lg ( keys %{$map} ) {
        for my $cM ( keys %{ $map->{$lg} } ) {
            for my $id ( keys %{ $map->{$lg}{$cM} } ) {
                my $marker = $marker_id->{$id};
                my $output = $markers->{$marker}{output};
                my $phased = $lg eq $first_lg ? mirror($output) : $output;
                $markers->{$marker}{output} = $phased;
            }
        }
    }
}

sub smooth_markers {
    my ( $map, $markers, $marker_id ) = @_;

    my @markerlist;
    for my $lg ( sort keys %{$map} ) {
        for my $cM ( sort { $a <=> $b } keys %{ $map->{$lg} } ) {
            for my $id ( sort { $a <=> $b } keys %{ $map->{$lg}{$cM} } ) {

               #                print
               #"$marker_id->{$id}\t$markers->{ $marker_id->{$id} }{output}\t$markers->{ $marker_id->{$id} }{length}\n";
                push @markerlist,
                  {
                    pattern  => $markers->{ $marker_id->{$id} }{output},
                    original => $marker_id->{$id},
                    bases    => $markers->{ $marker_id->{$id} }{length}
                  };
            }
        }
    }

    my %output;
    my $samples;
    for my $mref (@markerlist) {
        $samples = length $mref->{pattern};
        my @out = split //, $mref->{pattern};
        $output{ $mref->{pattern} } = \@out;
    }

    my %haplotypes;
    for my $sample ( 0 .. $samples - 1 ) {
        my $haplotype_bases = 0;
        my $haplotype_len   = 0;
        my $cur_gt          = '';
        for my $m ( 0 .. $#markerlist ) {
            my $gt = $output{ $markerlist[$m]{pattern} }[$sample];
            if ( $cur_gt ne $gt ) {

                #                print "$sample\t$cur_gt\t$haplotype_len\t$haplotype_bases\n" if $cur_gt ne '';
                if ( $cur_gt eq '' ) {    #Skip the first haplotype, as this will often be short
                    $cur_gt = $gt;
                    next;
                }
                $cur_gt = $gt;
                if ( $haplotype_len <= 2 and $haplotype_bases <= 5000 ) {
                    for my $h ( 1 .. $haplotype_len ) {
                        $output{ $markerlist[ $m - $h ]{pattern} }[$sample] =
                          $cur_gt;
                    }
                }
                $haplotype_len++;
                $haplotype_bases = $markerlist[$m]{bases};
            }
            else {
                $haplotype_len = 1;
                $haplotype_bases += $markerlist[$m]{bases};
            }
        }

        #        print "$sample\t$cur_gt\t$haplotype_len\t$haplotype_bases\n";
    }

    map { $markers->{ $_->{original} }{output} = join '', @{ $output{ $_->{pattern} } }; } @markerlist;
}

sub write_map_markers {
    my ( $map, $markers, $marker_id, $output, $chromosome_print ) = @_;

    my %chromosome;
    open my $map_marker_file, ">", "$output.$chromosome_print.mstmap.markers"
      or croak "Can't open $output.$chromosome_print.mstmap.markers: $OS_ERROR\n";

    print $map_marker_file "ID\tLG\tCM\tOriginal\tOutput\tLength\n";
    for my $lg ( sort keys %{$map} ) {
        for my $cM ( sort { $a <=> $b } keys %{ $map->{$lg} } ) {
            for my $id ( sort { $a <=> $b } keys %{ $map->{$lg}{$cM} } ) {
                my $marker = $marker_id->{$id};
                $chromosome{$lg}{$cM}{$marker}{length} =
                  $markers->{$marker}{length};
                $chromosome{$lg}{$cM}{$marker}{output} =
                  $markers->{$marker}{output};
                print $map_marker_file "$id\t$lg\t$cM";
                print $map_marker_file "\t$marker\t$markers->{$marker}{output}";
                print $map_marker_file "\t$markers->{$marker}{length}";
                print $map_marker_file "\n";
            }
        }
    }
    close $map_marker_file;

    \%chromosome;
}

## MAP GENOME

sub output_marker_blocks {
    my ( $linkage_groups, $marker_blocks, $output ) = @_;

    open my $chromosome_map_file, '>', "$output.chromosome.map.tsv"
      or croak "Can't open chromosome map file! $OS_ERROR\n";
    print $chromosome_map_file "Chromosome\tPrint\tcM\tOriginalMarker\tCleanMarker\tLength\n";
    open my $scaffold_map_file, '>', "$output.scaffold.map.tsv"
      or croak "Can't open scaffold map file! $OS_ERROR\n";
    print $scaffold_map_file "Chromosome\tcM\tScaffold\tStart\tEnd\tLength\n";

    my %chromosome_numbers;
    my $new_chromosome_number = 100;
    for my $print ( keys %{$linkage_groups} ) {
        my $chromosome = $chromosomes{$print} // $new_chromosome_number;
        $new_chromosome_number++ if $chromosome == $new_chromosome_number;
        $chromosome_numbers{$chromosome} = $print;
    }

    my %genome;
    for my $chromosome ( sort { $a <=> $b } keys %chromosome_numbers ) {
        my $chromosome_print = $chromosome_numbers{$chromosome};

        open my $chrom_file, '<', "$output.$chromosome_print.mstmap.markers"
          or croak "Can't open file for chromosome $chromosome! $OS_ERROR\n";
        my $header = <$chrom_file>;
        while ( my $chrom_line = <$chrom_file> ) {
            chomp $chrom_line;
            my ( $id, $lg, $cM, $original, $output, $length ) = split /\t/, $chrom_line;
            $genome{$chromosome_print}{$lg}{$cM}{$original}{output} = $output;
            $genome{$chromosome_print}{$lg}{$cM}{$original}{length} = $length;
        }
        close $chrom_file;

        croak "More than one linkage group found for chromosome print $chromosome_print"
          if keys %{ $genome{$chromosome_print} } > 1;
        for my $lg ( sort keys %{ $genome{$chromosome_print} } ) {
            for my $cM (
                sort { $a <=> $b }
                keys %{ $genome{$chromosome_print}{$lg} }
              )
            {
                for my $marker ( keys %{ $genome{$chromosome_print}{$lg}{$cM} } ) {

                    print $chromosome_map_file "$chromosome\t$chromosome_print\t$cM\t$marker\t";

                    map {
                        my $colour = $_ eq 'A' ? 'red1' : $_ eq 'B' ? 'royalblue1' : $_ eq 'H' ? 'sandybrown' : 'bold';
                        print $chromosome_map_file fg $colour, $_;
                    } split //, $genome{$chromosome_print}{$lg}{$cM}{$marker}{output};
                    print $chromosome_map_file "\t$genome{$chromosome_print}{$lg}{$cM}{$marker}{length}";
                    print $chromosome_map_file "\n";

                    for my $scaffold ( sort keys %{ $marker_blocks->{$marker} } ) {
                        my $block_start = -1;
                        my $last_end    = -1;
                        for my $start (
                            sort { $a <=> $b }
                            keys %{ $marker_blocks->{$marker}{$scaffold} }
                          )
                        {
                            $block_start = $start if $block_start == -1;
                            if ( $last_end != -1 and $start > $last_end + 1 ) {
                                my $length = $last_end - $block_start + 1;
                                print $scaffold_map_file
                                  "$chromosome\t$cM\t$scaffold\t$block_start\t$last_end\t$length\n";
                                $block_start = $start;
                            }
                            $last_end =
                              $marker_blocks->{$marker}{$scaffold}{$start};
                        }
                        my $length = $last_end - $block_start + 1;
                        print $scaffold_map_file "$chromosome\t$cM\t$scaffold\t$block_start\t$last_end\t$length\n";
                    }
                }
            }
        }
        print $chromosome_map_file "\n";
    }
    close $chromosome_map_file;
    close $scaffold_map_file;
}

## STATISTICS

sub get_block_stats {
    my ( $step, $blocklist, $metadata ) = @_;

    my $stats = get_stats( $blocklist, $metadata );

    $metadata->{stats} = $stats;
    return if !$metadata->{verbose};

    for my $type ( @{ $metadata->{header} } ) {
        next if !defined $metadata->{types}{$type};
        output_type_stats( $step, $type, $stats ) if $metadata->{verbose};
    }

    return;
}

sub update_block_stats {
    my ( $step, $blocklist, $metadata ) = @_;

    my $stats = get_stats( $blocklist, $metadata );

    my $all_equal = 1;

    for my $type ( @{ $metadata->{header} } ) {
        next if !defined $metadata->{types}{$type};

        next if stats_equal( $metadata->{stats}{$type}, $stats->{$type} );

        $all_equal = 0;

        output_type_stats( $step, $type, $stats ) if $metadata->{verbose};
    }

    if ( $all_equal and $metadata->{verbose} ) {
        printf STDERR "%-60s", $step;
        print STDERR "No change\n";
    }

    $metadata->{stats} = $stats;

    return;
}

sub output_type_stats {
    my ( $step, $type, $stats ) = @_;

    printf STDERR "%-60s", $step;
    printf STDERR "%-25s", $type;
    for my $full ('ok') {
        printf STDERR "%12d", $stats->{$type}{$full}{patterns};
        printf STDERR "%12d", $stats->{$type}{$full}{blocks};
        printf STDERR "%12d", $stats->{$type}{$full}{bases};
    }
    print STDERR "\n";

    return;
}

sub get_stats {
    my ( $blocklist, $metadata ) = @_;

    my %stats;
    my $empty = $metadata->{empty};

    my %validities;
    for my $block ( @{$blocklist} ) {
        for my $type ( sort keys %{ $metadata->{types} } ) {

            next if $block->{$type} eq $empty;

            my $full = $block->{$type} =~ /\-/ ? 'no' : 'ok';

            $stats{$type}{$full}{blocks}++;
            $stats{$type}{$full}{bases} += $block->{'length'};
            $stats{$type}{$full}{patterns}{ $block->{$type} }++;
        }
        $validities{ $block->{validity} }{blocks}++;
        $validities{ $block->{validity} }{bases} += $block->{'length'};
        $validities{ $block->{validity} }{patterns}{ $block->{'Paternal'} }++;
    }
    for my $validity (
        sort { $validities{$b}{bases} <=> $validities{$a}{bases} }
        keys %validities
      )
    {
        printf STDERR "%s\t%5d\t%5d\t%10d\n", $validity, $validities{$validity}{blocks},
          scalar keys %{ $validities{$validity}{patterns} },
          $validities{$validity}{bases};
    }
    for my $type ( keys %{ $metadata->{types} } ) {
        for my $full ( 'ok', 'no' ) {
            $stats{$type}{$full}{blocks} = $stats{$type}{$full}{blocks} // 0;
            $stats{$type}{$full}{bases}  = $stats{$type}{$full}{bases}  // 0;
            $stats{$type}{$full}{patterns} =
              defined $stats{$type}{$full}{patterns}
              ? scalar keys %{ $stats{$type}{$full}{patterns} }
              : 0;
        }
    }

    \%stats;
}

sub stats_equal {
    my ( $oldstats, $newstats ) = @_;
    my $stats_equal = 1;

    for my $full ( 'ok', 'no' ) {
        $stats_equal = 0
          if $oldstats->{$full}{patterns} ne $newstats->{$full}{patterns};
        $stats_equal = 0
          if $oldstats->{$full}{blocks} ne $newstats->{$full}{blocks};
        $stats_equal = 0
          if $oldstats->{$full}{bases} ne $newstats->{$full}{bases};
    }

    $stats_equal;
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
          and $pattern2[$i] ne 'H';
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