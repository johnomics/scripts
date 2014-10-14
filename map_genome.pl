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
my ( $scfmap, $genome, $markerscf ) = make_linkage_maps( $blocklist, $args{output}, $linkage_groups, $marker_blocks );

print STDERR "MAPPING GENOME\n";
map_genome( $blocklist, $scfmap, $genome, $markerscf, $args{output} );
print STDERR "Done\n";

exit;

## INFER MARKERS
sub infer_markers {
    my ( $input, $metadata ) = @_;

    print STDERR "Loading blocks...\n";
    my $blocklist = load_blocks( $input, $metadata );

    fill_blocks( $blocklist, $metadata );

    correct_maternal( "Maternal-AHAH", $blocklist, $metadata );

    collapse( $blocklist, $metadata );

    my ( $linkage_groups, $patmat, $pattern_blocks ) = build_linkage_groups( $blocklist, $metadata );

    #    print STDERR "Processing Intercross patterns...\n";
    #    clean_intercross_patterns( $matpat, $patmat, $pattern_block, $metadata );

    print STDERR "Writing clean blocks...\n";
    output_blocks( $input, $blocklist );

    ( $blocklist, $linkage_groups, $pattern_blocks );
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
        if ( $block->{"Maternal-AHAH"} ne "" ) {
            $samplenum = length $block->{"Maternal-AHAH"};
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
            $block->{$type} = convert_sex( $block->{$type} ) if defined $sex{$type};
            $block->{$type} = phase( $block->{$type} )       if defined $types->{$type};
        }
    }
    $metadata->{header}  = $header;
    $metadata->{samples} = $samplenum;
    $metadata->{types}   = $types;

    get_block_stats( "After loading blocks", $blocklist, $metadata );

    $blocklist;
}

sub convert_sex {
    my ($pat) = @_;
    $pat =~ s/B/H/g;
    $pat;
}

sub correct_maternal {
    my ( $maternal_type, $blocklist, $metadata ) = @_;

    my $empty = ' ' x $metadata->{samples};

    fill_maternal( $blocklist, $metadata );
    update_block_stats( "After first maternal fill", $blocklist, $metadata );

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
        if ( $maternal->{$a}{length} < 10000 ) {
            for my $i ( @{ $maternal->{$a}{blocks} } ) {
                $blocklist->[$i]{$maternal_type} = $empty;
            }
        }
    }

    update_block_stats( "After correcting maternal patterns", $blocklist, $metadata );

    fill_maternal( $blocklist, $metadata );
    update_block_stats( "After second maternal fill", $blocklist, $metadata );

    ( $maternal, $maternal_patterns ) = get_type_blocks( $maternal_type, $blocklist, $metadata );

    infer_intercross_maternal( $maternal, $maternal_patterns, $blocklist, $metadata );

    update_block_stats( "After inferring intercross maternal patterns", $blocklist, $metadata );

    return;
}

sub infer_intercross_maternal {
    my ( $maternal, $maternal_patterns, $blocklist, $metadata ) = @_;

    my $empty = ' ' x $metadata->{samples};

    my ( $intercross_blocks, $intercross_patterns ) = get_type_blocks( "Intercross", $blocklist, $metadata );

    my %maternal_candidates;
    for my $intercross ( @{$intercross_patterns} ) {
        next if $intercross eq $empty;

        # Check to see if Maternal patterns are empty for this Intercross pattern and skip if not
        my %maternal_filled;
        for my $block ( @{ $intercross_blocks->{$intercross}{blocks} } ) {
            if ( $blocklist->[$block]{"Maternal-AHAH"} ne $empty ) {
                $maternal_filled{ $blocklist->[$block]{"Maternal-AHAH"} }++;
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
                next if $blocklist->[$block]{"Maternal-AHAH"} ne $empty;
                $blocklist->[$block]{"Maternal-AHAH"} = $matching_maternal[0];
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
                next if $blocklist->[$block]{"Maternal-AHAH"} ne $empty;
                $blocklist->[$block]{"Maternal-AHAH"} = $m;
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
    my $cm = { pattern => $c, length => $c_length };
    my @consensus_calls;

    for my $m ( $cm, @{$merged} ) {
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
        next if $type eq "Maternal-AHAH";
        fill_type_pair( $type, "Maternal-AHAH", $blocklist, $metadata );
    }
    return;
}

sub fill_blocks {
    my ( $blocklist, $metadata ) = @_;

    # Get all Paternal and Intercross types
    my %ordertypes = ( "Intercross" => '' );

    my @ordering;
    for my $type ( keys %{ $metadata->{types} } ) {
        my $shorttype;
        if ( $type =~ /(.+)-(.+)/ ) {
            $shorttype = $1;
        }
        if ( defined $ordertypes{$shorttype} ) {
            push @ordering, $type;
        }
    }
    push @ordering, "Paternal-AHAH";
    my %done;
    for my $a (@ordering) {
        for my $b (@ordering) {
            next if $a eq $b;
            next if defined $done{$a}{$b};
            fill_type_pair( $a, $b, $blocklist, $metadata );
            $done{$a}{$b}++;
            fill_type_pair( $b, $a, $blocklist, $metadata );
            $done{$b}{$a}++;

        }
    }

    fill_type_pair( "Paternal-AHAB_AHB", "Paternal-AHAB_AHA", $blocklist, $metadata );
    fill_type_pair( "Paternal-AHAB_AHA", "Paternal-AHAB_AHB", $blocklist, $metadata );
    fill_type_pair( "Paternal-AHAB_AHB", "Sex-HB",            $blocklist, $metadata );
    fill_type_pair( "Paternal-AHAB_AHA", "Sex-HB",            $blocklist, $metadata );

    update_block_stats( "After filling blocks", $blocklist, $metadata );

    return;
}

sub fill_type_pair {
    my ( $typea, $typeb, $blocklist, $metadata ) = @_;

    my $empty = ' ' x $metadata->{samples};
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
                $same = 0;
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
    my %patmat;
    my %pattern_blocks;
    my $sexmat;

    my $empty = ' ' x $metadata->{samples};
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
            ( $mat, $pat ) = convert_intercross_block( $block, $metadata );
        }
        next if $mat =~ /[ \-]/;
        next if $pat =~ /[ \-]/;

        # If maternal and paternal patterns found, store them
        $block->{'Maternal-AHAH'} = $mat
          if ( $block->{'Maternal-AHAH'} eq $empty );
        $block->{'Paternal-AHAH'} = $pat
          if ( $block->{'Paternal-AHAH'} eq $empty );
        $linkage_groups{$mat}{$pat}{length} += $block->{'length'};
        $linkage_groups{$mat}{$pat}{blocks}++;
        $patmat{$pat}{$mat}{length} += $block->{'length'};
        $patmat{$pat}{$mat}{blocks}++;
        $pattern_blocks{$pat}{ $block->{'scaffold'} }{ $block->{'start'} } =
          $block->{'end'};
    }

    for my $mat ( keys %linkage_groups ) {
        if ( $mat =~ /^(S+)$/ ) {
            for my $sexpat ( keys %{ $linkage_groups{$mat} } ) {
                $patmat{$sexpat}{$sexmat}{length} +=
                  $patmat{$sexpat}{$mat}{length};
                $patmat{$sexpat}{$sexmat}{blocks} +=
                  $patmat{$sexpat}{$mat}{blocks};
                delete $patmat{$sexpat}{$mat};
                $linkage_groups{$sexmat}{$sexpat}{length} +=
                  $linkage_groups{$mat}{$sexpat}{length};
                $linkage_groups{$sexmat}{$sexpat}{blocks} +=
                  $linkage_groups{$mat}{$sexpat}{blocks};
            }
            delete $linkage_groups{$mat};
        }
    }

    update_block_stats( "After filling maternal and paternal patterns", $blocklist, $metadata );

    ( \%linkage_groups, \%patmat, \%pattern_blocks );
}

sub convert_intercross_block {
    my ( $block, $metadata ) = @_;

    my @intercross_types;
    for my $type ( keys %{ $metadata->{types} } ) {
        push @intercross_types, $type if $type =~ /Intercross/;
    }

    for my $ic (@intercross_types) {
        next if $block->{$ic} =~ /[ \-]/;
        my $mat = $block->{'Maternal-AHAH'};
        if ( $mat =~ /^( +)$/ ) {
            return ( "I" x length( $block->{$ic} ), phase( $block->{$ic} ) );
        }
        my $pat = convert_intercross( $mat, $block->{$ic} );
        return ( $mat, $pat );
    }

    ( " ", " " );
}

sub clean_intercross_patterns {
    my ( $matpat, $patmat, $pattern_block, $metadata ) = @_;

    my $empty = ' ' x $metadata->{samples};
    my @int;
    my @pat;

    my %int_lengths;
    for my $pat ( keys %{$pattern_block} ) {
        if ( $pat =~ 'H' ) {
            push @int, $pat;
            for my $scf ( keys %{ $pattern_block->{$pat} } ) {
                for my $start ( keys %{ $pattern_block->{$pat}{$scf} } ) {
                    $int_lengths{$pat} += $pattern_block->{$pat}{$scf}{$start} - $start + 1;
                }
            }
        }
        else {
            push @pat, $pat;
        }
    }

    for my $int ( sort { $int_lengths{$a} <=> $int_lengths{$b} } keys %int_lengths ) {
        my $imat         = "I" x length $int;
        my $int_pmatched = 0;
        for my $pat (@pat) {
            if ( int_match( $int, $pat ) ) {
                $int_pmatched++;
                my $pmat = (
                    sort { $patmat->{$pat}{$b}{length} <=> $patmat->{$pat}{$a}{length} }
                      keys %{ $patmat->{$pat} }
                )[0];
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
            my $minh     = 2;
            my $minh_mat = $empty;

            #            for my $mat ( keys %{$matpat} ) {
            #                if ( int_match( $int, $mat ) ) {
            #                    $minh     = 0;
            #                    $minh_mat = $mat;
            #                    last;
            #                }

            #                my $forh = int_hamming( $int, $mat );
            #                my $revh = int_hamming( $int, mirror($mat) );
            #                my $math = $forh < $revh ? $forh : $revh;
            #                if ( $math < $minh ) {
            #                    $minh     = $math;
            #                    $minh_mat = $mat;
            #                }
            #            }

            #            print "$int\t$int_lengths{$int}\t$minh";
            #            if ( defined $chroms{$minh_mat} ) {
            #                print "\t$chroms{$minh_mat}";
            #            }
            #            else {
            #                print "\t-";
            #            }
            #            print "\n";
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

    update_block_stats( "After processing intercross patterns", $blocklist, $metadata );
    return;
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
    my ( $blocklist, $output, $matpat, $pattern_block ) = @_;

    my %scfmap;
    my %genome;
    my %markerscf;
    for my $mat ( keys %{$matpat} ) {

        my $paternal_patterns = keys %{ $matpat->{$mat} };
        my $chromosome = $chroms{$mat} // '-';
        printf STDERR "Building map for maternal pattern $mat chromosome %3s with %4d paternal patterns",
          $chromosome, $paternal_patterns;

        # Delete paternal patterns less than 1kb long or occurring in only 1 or 2 blocks
        for my $pat ( keys %{ $matpat->{$mat} } ) {
            delete $matpat->{$mat}{$pat}
              if (  $matpat->{$mat}{$pat}{length} < 1000
                and $matpat->{$mat}{$pat}{blocks} <= 2 );
        }

        # Skip maternal pattern if no paternal patterns left
        if ( keys %{ $matpat->{$mat} } == 0 ) {
            print STDERR "\tSkipping\n";
            next;
        }

        # Make map for this maternal pattern
        my $markercode = run_mstmap( $mat, $matpat->{$mat}, $output );
        $genome{$mat} = load_map( $output, $mat );

        my $linkage_groups = keys %{ $genome{$mat} };
        print STDERR "\tBuilt $linkage_groups linkage groups";

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
            my $linkage_groups = keys %{ $genome{$mat} };
            print STDERR "\tAfter phasing, have $linkage_groups linkage groups";
        }

        print STDERR "\n";

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

    ( \%scfmap, \%genome, \%markerscf );
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

    # Delete linkage groups with a single marker
    map { delete $lg{$_} if keys %{ $lg{$_} } == 1 } keys %lg;
    close $mstout;

    \%lg;
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

    \%marker_lookup;
}

sub check_mirror {
    my ( $marker, $markers ) = @_;
    my $mirror     = mirror($marker);
    my @markerlist = keys %{$markers};
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

    $id;
}

## MAP GENOME

sub map_genome {
    my ( $blocklist, $scfmap, $genome, $markerscf, $output ) = @_;

    my $scfstats = validate_scaffolds( $blocklist, $scfmap, $genome );

    print STDERR "Writing maps and stats...\n";
    output_genome_maps( $blocklist, $output, $genome, $scfstats );

    output_genome_stats( $scfstats, $scfmap, $markerscf );

    return;
}

sub validate_scaffolds {
    my ( $blocklist, $scfmap, $genome ) = @_;

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

    \%scfstats;
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

    return;
}

sub output_genome_maps {
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
    for my $mat ( keys %{$genome} ) {
        if ( defined $chroms{$mat} ) {
            $chromnums{$mat} = $chroms{$mat};
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
    return;
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

    @ordered_scfs;
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

    $snps;
}

sub output_genome_stats {
    my ( $scfstats, $scfmap, $markerscf ) = @_;

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

#            croak "Should only be one marker at non-oriented scaffold $scf!"
#              if ( keys %markers != 1 );
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

    return;
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
    my $empty = ' ' x $metadata->{samples};
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
    my ($maternal, $paternal) = @_;
    my @m = split //, $maternal;
    my @p = split //, $paternal;
    
    my @i = map {($m[$_] eq $p[$_]) ? $m[$_] : 'H'} 0..$#m;
    
    join '', @i;
}

sub convert_intercross {
    my ( $mat, $int ) = @_;
    return $int if $mat =~ /I/;
    my @int = split //, phase($int);
    my @mat = split //, $mat;
    my @pat =
      map { $int[$_] eq '-' ? '-' : $int[$_] eq 'H' ? ( $mat[$_] eq 'A' ? 'B' : 'A' ) : $mat[$_] } 0 .. $#int;
    my $pat = join '', @pat;

    $pat;
}

sub phase {
    my $pat = shift;
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
