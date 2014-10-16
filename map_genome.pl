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
memoize "phase";

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

my $options_okay =
  GetOptions( 'input=s' => \$args{input}, 'output=s' => \$args{output}, 'verbose' => \$args{verbose} );

my $metadata = { verbose => $args{verbose} };

print STDERR "INFER MARKERS\n";
my ( $blocklist, $linkage_groups, $marker_blocks ) = infer_markers( $args{input}, $metadata );

print STDERR "MAKING LINKAGE MAPS\n";
my $genome = make_linkage_maps( $blocklist, $args{output}, $linkage_groups, $marker_blocks );

output_marker_blocks( $genome, $marker_blocks, $args{output} );

print STDERR "Done\n";

exit;

## INFER MARKERS
sub infer_markers {
    my ( $input, $metadata ) = @_;

    print STDERR "Loading blocks...\n";
    my $blocklist = load_blocks( $input, $metadata );

    validate_blocks( $blocklist, $metadata );

    fill_blocks( $blocklist, $metadata );

    correct_maternal( "Maternal", $blocklist, $metadata );

    collapse( $blocklist, $metadata );

    my ( $linkage_groups, $marker_blocks ) = build_linkage_groups( $blocklist, $metadata );

    print STDERR "Writing clean blocks...\n";
    output_blocks( $input, $blocklist );

    ( $blocklist, $linkage_groups, $marker_blocks );
}

sub load_blocks {
    my ( $input, $metadata ) = @_;

    my $blocklist = [];
    my $header;
    my $types;

    print STDERR "Load blocks from database...\n" if $metadata->{verbose};
    my @inputfiles = split ',', $input;
    for my $inputfile (@inputfiles) {
        my $dbh = DBI->connect( "dbi:SQLite:dbname=$inputfile", "", "" );

        ( $header, $types ) = get_header_types($dbh);

        my $sth = $dbh->prepare("SELECT * FROM blocks");
        my $fileblocklist = $dbh->selectall_arrayref( $sth, { Slice => {} } );
        push @{$blocklist}, @{$fileblocklist};
        $sth->finish;

        $dbh->disconnect;
    }

    print STDERR "Phasing blocks...\n" if $metadata->{verbose};

    # Make empty pattern
    my $empty;
    my $samplenum;
    for my $block ( @{$blocklist} ) {
        if ( $block->{Maternal} ne "" ) {
            $samplenum = length $block->{Maternal};
            $empty     = ' ' x $samplenum;
            last;
        }
    }

    for my $block ( @{$blocklist} ) {
        for my $type ( keys %{$block} ) {
            if ( $block->{$type} eq "" ) {
                $block->{$type} = $empty;
                next;
            }
            $block->{$type} = phase( $block->{$type} ) if defined $types->{$type};
        }
    }
    $metadata->{header}  = $header;
    $metadata->{samples} = $samplenum;
    $metadata->{types}   = $types;
    $metadata->{empty}   = $empty;

    get_block_stats( "After loading blocks", $blocklist, $metadata );

    $blocklist;
}

sub validate_blocks {
    my ( $blocklist, $metadata ) = @_;

    my $empty = $metadata->{empty};

    for my $block ( @{$blocklist} ) {
        my $hom_consistent = check_consistent( $block, "IntercrossHom", $empty );
        my $het_consistent = check_consistent( $block, "IntercrossHet", $empty );

        $block->{IntercrossHet} = $empty if $block->{IntercrossHom} ne $empty and $hom_consistent and !$het_consistent;
        $block->{IntercrossHom} = $empty if !$hom_consistent and $block->{IntercrossHet} ne $empty and $het_consistent;
        empty_block( $block, $metadata ) if !$hom_consistent and !$het_consistent;
    }

    get_block_stats( "After validating blocks", $blocklist, $metadata );
}

sub check_consistent {
    my ( $block, $intercross_type, $empty ) = @_;

    my $maternal   = $block->{"Maternal"};
    my $paternal   = $block->{"Paternal"};
    my $intercross = $block->{$intercross_type};

    my $consistent = 1;

    if ( $maternal ne $empty and $paternal ne $empty and $intercross ne $empty ) {
        $consistent = 0
          if create_intercross( $maternal, $paternal ) ne $intercross
          and create_intercross( $maternal,         mirror($paternal) ) ne $intercross
          and create_intercross( mirror($maternal), $paternal ) ne $intercross
          and create_intercross( mirror($maternal), mirror($paternal) ) ne $intercross;
    }
    elsif ( $maternal eq $empty and $paternal ne $empty and $intercross ne $empty ) {
        my $ok = 0;
        if ( consistent( $paternal, $intercross ) ) {
            $block->{"Maternal"} = separate_intercross( $paternal, $intercross );
            $ok = 1;
        }
        elsif ( consistent( mirror($paternal), $intercross ) ) {
            $block->{"Maternal"} = separate_intercross( mirror($paternal), $intercross );
            $ok = 1;
        }
        $consistent = $ok;
    }

    return $consistent;
}

sub empty_block {
    my ( $block, $metadata ) = @_;
    for my $type ( keys %{ $metadata->{types} } ) {
        $block->{$type} = $metadata->{empty};
    }
    return;
}

sub correct_maternal {
    my ( $maternal_type, $blocklist, $metadata ) = @_;

    my $empty = $metadata->{empty};

    my ( $maternal, $maternal_patterns ) = get_type_blocks( $maternal_type, $blocklist, $metadata );

    my %merged;
    for my $a ( @{$maternal_patterns} ) {
        next if defined( $merged{$a} ) or $a eq $empty or $a =~ /\-/;
        for my $b ( @{$maternal_patterns} ) {
            next
              if $b eq $empty
              or defined( $merged{$b} )
              or $a eq $b
              or $maternal->{$a}{length} < $maternal->{$b}{length};

            if ( hamming( $a, $b ) <= 6 or hamming( mirror($a), $b ) <= 6 ) {
                $merged{$b} = $a;
                for my $i ( @{ $maternal->{$b}{blocks} } ) {
                    $blocklist->[$i]{$maternal_type} = $a;
                }
            }
        }
    }

    ( $maternal, $maternal_patterns ) = get_type_blocks( $maternal_type, $blocklist, $metadata );

    for my $a ( @{$maternal_patterns} ) {
        for my $i ( @{ $maternal->{$a}{blocks} } ) {
            if ( $maternal->{$a}{length} < 10000 ) {
                $blocklist->[$i]{$maternal_type} = $empty;
            }
            else {
                $blocklist->[$i]{$maternal_type} = phase( $blocklist->[$i]{$maternal_type} );
            }
        }
    }

    update_block_stats( "After correcting maternal patterns", $blocklist, $metadata );

    fill_maternal( $blocklist, $metadata );
    update_block_stats( "After filling maternal patterns", $blocklist, $metadata );

    ( $maternal, $maternal_patterns ) = get_type_blocks( $maternal_type, $blocklist, $metadata );

    infer_intercross_maternal( $maternal, $maternal_patterns, $blocklist, $metadata );

    update_block_stats( "After inferring intercross maternal patterns", $blocklist, $metadata );

    return;
}

sub infer_intercross_maternal {
    my ( $maternal, $maternal_patterns, $blocklist, $metadata ) = @_;

    my $empty = $metadata->{empty};

    my ( $intercross_blocks, $intercross_patterns ) = get_type_blocks( "Intercross", $blocklist, $metadata );

    my %maternal_candidates;
    for my $intercross ( @{$intercross_patterns} ) {
        next if $intercross eq $empty;

        # Check to see if Maternal patterns are empty for this Intercross pattern and skip if not
        my %maternal_filled;
        for my $block ( @{ $intercross_blocks->{$intercross}{blocks} } ) {
            if ( $blocklist->[$block]{"Maternal"} ne $empty ) {
                $maternal_filled{ $blocklist->[$block]{"Maternal"} }++;
            }
        }
        next if ( keys %maternal_filled );

        # Get any matching Maternal patterns for this Intercross marker
        my @matching_maternal;
        for my $maternal ( @{$maternal_patterns} ) {
            next if $maternal eq $empty;
            if ( int_hamming( $intercross, $maternal ) <= 6 or int_hamming( mirror($intercross), $maternal ) <= 6 ) {
                push @matching_maternal, $maternal;
            }
        }

        # If only one Maternal pattern found, set Intercross blocks to this Maternal pattern
        if ( @matching_maternal == 1 ) {
            for my $block ( @{ $intercross_blocks->{$intercross}{blocks} } ) {
                next if $blocklist->[$block]{"Maternal"} ne $empty;
                $blocklist->[$block]{"Maternal"} = $matching_maternal[0];
                $blocklist->[$block]{"Paternal"} = separate_intercross( $blocklist->[$block]{"Maternal"}, $intercross )
                  if $blocklist->[$block]{"Paternal"} eq $empty;
            }
        }

        # If no Maternal patterns found, they need to be inferred from the Intercross pattern
        elsif ( @matching_maternal == 0 ) {
            push @{ $maternal_candidates{$intercross}{patterns} }, $intercross;
            $maternal_candidates{$intercross}{length} = $intercross_blocks->{$intercross}{length};
        }
        else {
            # More than one matching Maternal pattern found
        }
    }

    update_block_stats( "After matching intercross to existing maternals", $blocklist, $metadata );

    my $merged = merge_maternal_candidates( \%maternal_candidates );

    # Update maternal patterns for inferred intercross chromosomes
    for my $m ( keys %{$merged} ) {
        for my $i ( @{ $merged->{$m}{patterns} } ) {
            my $pattern = $i->{original};
            for my $block ( @{ $intercross_blocks->{$pattern}{blocks} } ) {
                $blocklist->[$block]{"Maternal"} = $m if $blocklist->[$block]{"Maternal"} eq $empty;
                $blocklist->[$block]{"Paternal"} = separate_intercross( $m, $pattern )
                  if $blocklist->[$block]{"Paternal"} eq $empty;
            }
        }
    }

    return;
}

sub merge_maternal_candidates {
    my ($candidates) = @_;

    my @candidates = sort { $candidates->{$b}{length} <=> $candidates->{$a}{length} }
      keys %{$candidates};

    my $longest_candidate = shift @candidates;
    my %merged            = (
        $longest_candidate => {
            patterns => [
                {
                    pattern  => $longest_candidate,
                    original => $longest_candidate,
                    length   => $candidates->{$longest_candidate}{length}
                }
            ],
            length => $candidates->{$longest_candidate}{length}
        }
    );

    for my $c (@candidates) {
        my $patterns = @{ $candidates->{$c}{patterns} };

        my $update_pattern = $c;
        my $out_c          = $c;

        my @merged = sort { $merged{$b}{length} <=> $merged{$a}{length} } keys %merged;
        for my $m (@merged) {
            my $new_pattern = "";
            if ( int_hamming( $c, $m ) <= 6 ) {
                $new_pattern = get_intercross_consensus( $c, $candidates->{$c}{length}, $merged{$m}{patterns} );
            }

            if ( int_hamming( mirror($c), $m ) <= 6 ) {
                $new_pattern = get_intercross_consensus( mirror($c), $candidates->{$c}{length}, $merged{$m}{patterns} );
                $out_c = mirror($c);
            }

            # Update merged pattern
            if ( $new_pattern ne "" ) {
                $update_pattern = $new_pattern;
                if ( $m ne $update_pattern ) {
                    push @{ $merged{$update_pattern}{patterns} }, @{ $merged{$m}{patterns} };
                    $merged{$update_pattern}{length} +=
                      $merged{$m}{length};
                    delete $merged{$m};
                }
                last;
            }
        }

        # Add candidate to merged patterns
        push @{ $merged{$update_pattern}{patterns} },
          { original => $c, pattern => $out_c, length => $candidates->{$c}{length} };
        $merged{$update_pattern}{length} += $candidates->{$c}{length};
        $patterns = @{ $merged{$update_pattern}{patterns} };
    }
    \%merged;
}

sub get_intercross_consensus {
    my ( $c, $c_length, $merged ) = @_;
    my $cM = { pattern => $c, length => $c_length };
    my @consensus_calls;

    for my $m ( $cM, @{$merged} ) {
        my @m = split //, $m->{pattern};
        for my $i ( 0 .. $#m ) {
            $consensus_calls[$i]{ $m[$i] } += $m->{length};
        }
    }

    my $consensus = "";
    for my $con (@consensus_calls) {
        if ( keys %{$con} == 1 ) {
            $consensus .= ( keys %{$con} )[0];
        }
        else {
            map { $con->{$_} = 0 if !defined $con->{$_} } ( 'A', 'B', 'H' );
            my $max = ( $con->{'A'} > $con->{'B'} ) ? 'A' : 'B';
            $consensus .= ( $con->{$max} > ( $con->{'H'} * 0.1 ) ) ? $max : 'H';
        }
    }
    $consensus;
}

sub get_type_blocks {
    my ( $typename, $blocklist, $metadata ) = @_;

    my @types;
    for my $type ( keys %{ $metadata->{types} } ) {
        push @types, $type if $type =~ /$typename/;
    }
    my %type_blocks;
    for my $i ( 0 .. $#{$blocklist} ) {
        for my $type (@types) {
            my $block = $blocklist->[$i];
            $type_blocks{ $block->{$type} }{type} = $type;
            $type_blocks{ $block->{$type} }{length} += $block->{length};
            push @{ $type_blocks{ $block->{$type} }{blocks} }, $i;
        }
    }

    my @type_patterns = sort { $type_blocks{$b}{length} <=> $type_blocks{$a}{length} } keys %type_blocks;

    ( \%type_blocks, \@type_patterns );
}

sub fill_maternal {
    my ( $blocklist, $metadata ) = @_;

    for my $type ( keys %{ $metadata->{types} } ) {
        next if $type eq "Maternal";
        fill_type_pair( $type, "Maternal", $blocklist, $metadata );
    }
    return;
}

sub fill_blocks {
    my ( $blocklist, $metadata ) = @_;

    fill_type_pair( "Paternal",      "IntercrossHet", $blocklist, $metadata );
    fill_type_pair( "Paternal",      "IntercrossHom", $blocklist, $metadata );
    fill_type_pair( "IntercrossHet", "Paternal",      $blocklist, $metadata );
    fill_type_pair( "IntercrossHom", "Paternal",      $blocklist, $metadata );
    fill_type_pair( "IntercrossHet", "IntercrossHom", $blocklist, $metadata );
    fill_type_pair( "IntercrossHom", "IntercrossHet", $blocklist, $metadata );

    fill_maternal( $blocklist, $metadata );

    update_block_stats( "After filling blocks", $blocklist, $metadata );

    return;
}

sub fill_type_pair {
    my ( $typea, $typeb, $blocklist, $metadata ) = @_;

    my $empty = $metadata->{empty};
    my %ab;

    # Pull out type a and b fields from all blocks and calculate the total length covered by these types
    for my $b ( 0 .. $#{$blocklist} ) {
        $ab{ $blocklist->[$b]{$typea} }{ $blocklist->[$b]{$typeb} }{length} +=
          $blocklist->[$b]{'length'};
        push @{ $ab{ $blocklist->[$b]{$typea} }{ $blocklist->[$b]{$typeb} }{blocks} }, $b;
    }

    for my $a ( keys %ab ) {

        next if $a eq $empty or $a =~ /\-/ or keys %{ $ab{$a} } == 1;

        # Find type b pattern most found with this type a pattern
        my $maxb        = "";
        my $maxb_length = 0;
        for my $b ( keys %{ $ab{$a} } ) {
            if ( $ab{$a}{$b}{length} > $maxb_length ) {
                $maxb        = $b;
                $maxb_length = $ab{$a}{$b}{length};
            }
        }
        next if $maxb eq $empty;
        next if $maxb =~ /\-/;

        # Correct all b patterns associated with this a pattern
        # to the most often found b pattern,
        # assuming the distance between patterns is less than 6bp
        for my $b ( keys %{ $ab{$a} } ) {
            next if $b eq $maxb;

            my $fix = 0;
            $fix = 1 if $b eq $empty;
            $fix = 1 if hamming( $maxb, $b ) <= 6;
            $fix = 1 if hamming( mirror($maxb), $b ) <= 6;

            if ($fix) {
                for my $i ( @{ $ab{$a}{$b}{blocks} } ) {
                    $blocklist->[$i]{$typeb} = $maxb;
                }
            }
        }
    }

    return;
}

sub collapse {
    my ( $blocklist, $metadata ) = @_;
    my $i = 0;
    my $j = 1;
    while ( $j <= $#{$blocklist} ) {
        my $same = 1;
        $same = 0
          if ( $blocklist->[$i]{'scaffold'} ne $blocklist->[$j]{'scaffold'} );
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
            $blocklist->[$i]{'end'}    = $blocklist->[$j]{'end'};
            $blocklist->[$i]{'length'} = $blocklist->[$i]{'length'} + $blocklist->[$j]{'length'};
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

sub build_linkage_groups {
    my ( $blocklist, $metadata ) = @_;

    my %linkage_groups;
    my %marker_blocks;

    for my $block ( @{$blocklist} ) {
        my $maternal = $block->{Maternal};
        my $paternal = $block->{Paternal};

        next if ( $maternal eq $metadata->{empty} or $paternal eq $metadata->{empty} );

        $linkage_groups{$maternal}{$paternal}{length} += $block->{length};
        $linkage_groups{$maternal}{$paternal}{blocks}++;
        $linkage_groups{$maternal}{$paternal}{lookup} = $paternal;
        $marker_blocks{$paternal}{ $block->{scaffold} }{ $block->{start} } =
          $block->{end};
    }

    ( \%linkage_groups, \%marker_blocks );
}

sub int_hamming {

    # Calculate hamming distance between intercross and paternal pattern
    my ( $int, $pat ) = @_;
    my $hamming = 0;
    my @i       = split //, $int;
    my @p       = split //, $pat;
    for my $b ( 0 .. $#i ) {
        next if $i[$b] eq 'H' or $i[$b] eq '-' or $p[$b] eq 'H' or $p[$b] eq '-';
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

sub output_blocks {
    my $input     = shift;
    my $blocklist = shift;

    my $inputfile = ( split ',', $input )[0];
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
        my @values = map { $block->{$_} } @{$header};
        $insert_handle->execute(@values);
    }
    $dbh->commit;
    $dbh->disconnect;

    return;
}

## MAKE LINKAGE MAPS

sub make_linkage_maps {
    my ( $blocklist, $output, $linkage_groups, $marker_block ) = @_;

    my %genome;

    for my $chromosome_print ( keys %{$linkage_groups} ) {
        my $chromosome = $chromosomes{$chromosome_print} // '-';
        printf STDERR "Building map for print $chromosome_print chromosome %3s with %4d paternal markers",
          $chromosome, scalar keys %{ $linkage_groups->{$chromosome_print} };

        clean_print( $linkage_groups, $chromosome_print );
        next if !defined $linkage_groups->{$chromosome_print};

        $genome{$chromosome_print} =
          make_linkage_map( $chromosome_print, $linkage_groups->{$chromosome_print}, $output );

    }

    \%genome;
}

sub clean_print {
    my ( $linkage_groups, $chromosome_print ) = @_;

    # Delete paternal patterns less than 1kb long or occurring in only 1 or 2 blocks
    for my $marker ( keys %{ $linkage_groups->{$chromosome_print} } ) {
        $linkage_groups->{$chromosome_print}{$marker}{output} = $marker;
        delete $linkage_groups->{$chromosome_print}{$marker}
          if (  $linkage_groups->{$chromosome_print}{$marker}{length} < 1000
            and $linkage_groups->{$chromosome_print}{$marker}{blocks} <= 2 );
    }

    if ( keys %{ $linkage_groups->{$chromosome_print} } == 0 ) {
        delete $linkage_groups->{$chromosome_print};
    }
}

sub make_linkage_map {
    my ( $chromosome_print, $markers, $output ) = @_;

    # Make map for this maternal pattern
    my $marker_id = run_mstmap( $chromosome_print, $markers, $output );

    my $map = load_map( $output, $chromosome_print );
    my $chromosome_lgs = keys %{$map};
    print STDERR "\tBuilt $chromosome_lgs linkage groups";

    # If more than one linkage map returned for this chromosome,
    # attempt to rephase markers and remake the map
    my $phased_markers;
    if ( $chromosome_lgs == 2 ) {
        phase_linkage_groups( $map, $markers, $marker_id );
        $marker_id = run_mstmap( $chromosome_print, $markers, $output );
        $map = load_map( $output, $chromosome_print );
        my $chromosome_lgs = keys %{$map};
        print STDERR "\tAfter phasing, $chromosome_lgs linkage groups";
    }

    smooth_markers( $map, $markers, $marker_id );
    $marker_id = run_mstmap( $chromosome_print, $markers, $output );

    my $map_markers = write_map_markers( $map, $markers, $marker_id, $output, $chromosome_print );
    print STDERR "\n";

    $map_markers;
}

sub run_mstmap {

    my ( $chromosome_print, $markers, $output ) = @_;

    my $mstmap_input = write_mstmap_header( $output, $chromosome_print, $markers );

    my $id = 1;
    my %marker_id;
    for my $marker ( keys %{$markers} ) {
        $marker_id{$id} = $marker;
        $markers->{$marker}{output} = check_mirror( $markers->{$marker}{output}, $markers );
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
        no_map_dist                  => "15",
        no_map_size                  => "2",
        missing_threshold            => "0.75",
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
                push @markerlist, { pattern => $markers->{ $marker_id->{$id} }{output}, original => $marker_id->{$id} };
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

    for my $sample ( 0 .. $samples - 1 ) {
        my $haplotype_len = 0;
        my $cur_gt        = '';
        for my $m ( 0 .. $#markerlist ) {
            my $gt = $output{ $markerlist[$m]{pattern} }[$sample];
            if ( $cur_gt ne $gt ) {
                $cur_gt = $gt;
                if ( $haplotype_len <= 3 ) {
                    for my $h ( 1 .. $haplotype_len ) {
                        $output{ $markerlist[ $m - $h ]{pattern} }[$sample] = $cur_gt;
                    }
                }
                $haplotype_len = 1;
            }
            else {
                $haplotype_len++;
            }
        }
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
                $chromosome{$lg}{$cM}{$marker}{length} = $markers->{$marker}{length};
                $chromosome{$lg}{$cM}{$marker}{output} = $markers->{$marker}{output};
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
    my ( $genome, $marker_blocks, $output ) = @_;

    open my $chromosome_map_file, '>', "$output.chromosome.map.tsv"
      or croak "Can't open chromosome map file! $OS_ERROR\n";
    print $chromosome_map_file "Chromosome\tPrint\tcM\tOriginalMarker\tCleanMarker\tLength\n";
    open my $scaffold_map_file, '>', "$output.scaffold.map.tsv" or croak "Can't open scaffold map file! $OS_ERROR\n";
    print $scaffold_map_file "Chromosome\tcM\tScaffold\tStart\tEnd\tLength\n";

    my %chromosome_numbers;
    my $new_chromosome_number = 100;
    for my $print ( keys %{$genome} ) {
        my $chromosome = $chromosomes{$print} // $new_chromosome_number;
        $new_chromosome_number++ if $chromosome == $new_chromosome_number;
        $chromosome_numbers{$chromosome} = $print;
    }
    for my $chromosome ( sort { $a <=> $b } keys %chromosome_numbers ) {
        my $chromosome_print = $chromosome_numbers{$chromosome};
        croak "More than one linkage group found for chromosome print $chromosome_print"
          if keys %{ $genome->{$chromosome_print} } > 1;
        for my $lg ( sort keys %{ $genome->{$chromosome_print} } ) {
            for my $cM ( sort { $a <=> $b } keys %{ $genome->{$chromosome_print}{$lg} } ) {
                for my $marker ( keys %{ $genome->{$chromosome_print}{$lg}{$cM} } ) {
                    print $chromosome_map_file "$chromosome\t$chromosome_print\t$cM\t$marker";
                    print $chromosome_map_file "\t$genome->{$chromosome_print}{$lg}{$cM}{$marker}{output}";
                    print $chromosome_map_file "\t$genome->{$chromosome_print}{$lg}{$cM}{$marker}{length}";
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
                            $last_end = $marker_blocks->{$marker}{$scaffold}{$start};
                        }
                        my $length = $last_end - $block_start + 1;
                        print $scaffold_map_file "$chromosome\t$cM\t$scaffold\t$block_start\t$last_end\t$length\n";
                    }
                }
            }
        }
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
    for my $full ( 'ok', 'no' ) {
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
    for my $block ( @{$blocklist} ) {
        for my $type ( keys %{ $metadata->{types} } ) {
            next if $block->{$type} eq $empty;

            my $full = $block->{$type} =~ /\-/ ? 'no' : 'ok';

            $stats{$type}{$full}{blocks}++;
            $stats{$type}{$full}{bases} += $block->{'length'};
            $stats{$type}{$full}{patterns}{ $block->{$type} }++;
        }
    }
    for my $type ( keys %{ $metadata->{types} } ) {
        for my $full ( 'ok', 'no' ) {
            $stats{$type}{$full}{blocks} = $stats{$type}{$full}{blocks} // 0;
            $stats{$type}{$full}{bases}  = $stats{$type}{$full}{bases}  // 0;
            $stats{$type}{$full}{patterns} =
              defined $stats{$type}{$full}{patterns} ? scalar keys %{ $stats{$type}{$full}{patterns} } : 0;
        }
    }

    \%stats;
}

sub stats_equal {
    my ( $oldstats, $newstats ) = @_;
    my $stats_equal = 1;

    for my $full ( 'ok', 'no' ) {
        $stats_equal = 0 if $oldstats->{$full}{patterns} ne $newstats->{$full}{patterns};
        $stats_equal = 0 if $oldstats->{$full}{blocks} ne $newstats->{$full}{blocks};
        $stats_equal = 0 if $oldstats->{$full}{bases} ne $newstats->{$full}{bases};
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

sub consistent {
    my ( $pattern1, $pattern2 ) = @_;

    my @pattern1 = split //, $pattern1;
    my @pattern2 = split //, $pattern2;
    my $distance = 0;
    map { $distance++ if $pattern1[$_] ne $pattern2[$_] and $pattern1[$_] ne 'H' and $pattern2[$_] ne 'H' }
      0 .. $#pattern1;
    return ( $distance == 0 ) ? 1 : 0;
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
    my ( $a, $b ) = @_;
    my $hamming = ( $a ^ $b ) =~ tr/\001-\255//;

    $hamming;
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
