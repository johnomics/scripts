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

my %args;
$args{input} = "";
$args{db}    = "";
$args{agp}   = "";
my $options_okay = GetOptions( 'input=s' => \$args{input}, 'db=s' => \$args{db}, 'agp=s', \$args{agp} );

croak "Please supply a database name with -d" if $args{db} eq "";

my $agp = load_agp( $args{agp} );

my ( $chromosomes, $scaffold_blocks ) = load_blocks( $args{input} );

extend_blocks( $scaffold_blocks, $chromosomes, $args{db}, $agp );

exit;

sub load_agp {
    my $agp_filename = shift;

    open my $agp_file, "<", $agp_filename or croak "Can't open AGP file $agp_filename!\n";
    my %agp;
    while ( my $agp_line = <$agp_file> ) {
        my ( $scf, $start, $end, $part, $type, @f ) = split /\t/, $agp_line;
        $agp{$scf}{$part}{start} = $start;
        $agp{$scf}{$part}{end}   = $end;
        $agp{$scf}{$part}{type}  = $type;
    }
    close $agp_file;

    return \%agp;
}

sub load_blocks {
    my ($input) = @_;

    my $chromosomes = load_chromosomes($input);

    my $scaffold_blocks = load_scaffold_blocks($input, $chromosomes);

    ( $chromosomes, $scaffold_blocks );
}

sub load_chromosomes {
    my ($input) = @_;
    open my $chromosome_map, '<', "$input.chromosome.map.tsv" or croak "Can't open chromosome map! $OS_ERROR\n";

    my %chromosomes;
    my $header = <$chromosome_map>;
    while ( my $chromosome_line = <$chromosome_map> ) {
        chomp $chromosome_line;
        my ( $chromosome, $print, $cM, $original, $clean, $length ) = split /\t/, $chromosome_line;
        $chromosomes{$chromosome}{print}                  = $print;
        $chromosomes{$chromosome}{markers}{$cM}{original} = $original;
        $chromosomes{$chromosome}{markers}{$cM}{clean}    = $clean;
        $chromosomes{$chromosome}{markers}{$cM}{length}   = $length;
    }
    close $chromosome_map;

    \%chromosomes;
}

sub load_scaffold_blocks {
    my ($input, $chromosomes) = @_;

    open my $scaffold_map, '<', "$input.scaffold.map.tsv" or croak "Can't open scaffold map! $OS_ERROR\n";

    my %scaffold_blocks;
    my $header = <$scaffold_map>;
    while ( my $scaffold_line = <$scaffold_map> ) {
        chomp $scaffold_line;
        my ( $chromosome, $cM, $scaffold, $start, $end, $length ) = split /\t/, $scaffold_line;
        next if $scaffold =~ /sch/;
        $scaffold_blocks{$scaffold}{$start}{chromosome} = $chromosome;
        $scaffold_blocks{$scaffold}{$start}{cM}         = $cM;
        $scaffold_blocks{$scaffold}{$start}{end}        = $end;
        $scaffold_blocks{$scaffold}{$start}{length}     = $length;
        $scaffold_blocks{$scaffold}{$start}{type}       = 'Block';
        
        $chromosomes{$chromosome}{markers}{$cM}{scaffolds}{$scaffold} = $start;
    }
    close $scaffold_map;

    \%scaffold_blocks;
}

sub extend_blocks {
    my ( $scaffold_blocks, $chromosomes, $db, $agp ) = @_;

    for my $scaffold ( keys %{$scaffold_blocks} ) {

        load_gaps( $scaffold, $scaffold_blocks, $agp );

        fill_gaps( $scaffold, $scaffold_blocks, $chromosomes, $db );
    }
}

sub load_gaps {
    my ( $scaffold, $scaffold_blocks, $agp ) = @_;

    my $gap_start = 1;

    for my $block_start ( sort { $a <=> $b } keys %{ $scaffold_blocks->{$scaffold} } ) {
        my $gap_end = $block_start - 1;
        $scaffold_blocks->{$scaffold}{$gap_start}{end}    = $gap_end;
        $scaffold_blocks->{$scaffold}{$gap_start}{length} = $gap_end - $gap_start + 1;
        $scaffold_blocks->{$scaffold}{$gap_start}{type}   = "Gap";
        $gap_start                                        = $scaffold_blocks->{$scaffold}{$block_start}{end} + 1;
    }

    my $scaffold_length = 0;
    for my $part ( sort { $a <=> $b } keys %{ $agp->{$scaffold} } ) {
        $scaffold_length = $agp->{$scaffold}{$part}{end};
    }
    $scaffold_blocks->{$scaffold}{$gap_start}{end}    = $scaffold_length;
    $scaffold_blocks->{$scaffold}{$gap_start}{length} = $scaffold_length - $gap_start + 1;
    $scaffold_blocks->{$scaffold}{$gap_start}{type}   = "Gap";
}

sub fill_gaps {
    my ( $scaffold, $scaffold_blocks, $chromosomes, $db ) = @_;

    my $dbh = DBI->connect( "dbi:SQLite:dbname=$db", "", "" );

    my @block_starts = sort { $a <=> $b } keys %{ $scaffold_blocks->{$scaffold} };
    for my $start_i ( 0 .. $#block_starts ) {

        my $rejects = get_rejects( $block_starts[$start_i], $scaffold, $scaffold_blocks, $dbh );
        next if @{$rejects} == 0;

        my ( $previous, $next ) =
          get_neighbours( $start_i, $scaffold_blocks->{$scaffold}, \@block_starts, $chromosomes );

        for my $reject ( @{$rejects} ) {
            print "$reject->{position}\t$reject->{pattern}\t$previous\t$next";
            if ( match_reject( $reject->{pattern}, $previous ) ) {
                print " P";
            }
            else {
                print "  ";
            }
            if ( match_reject( $reject->{pattern}, $next ) ) {
                print " N";
            }
            else {
                print "  ";
            }

            print "\n";
        }
        print "\n";
    }

    $dbh->disconnect;
}

sub get_rejects {
    my ( $start, $scaffold, $scaffold_blocks, $dbh ) = @_;

    my $rejects = [];
    my $block   = $scaffold_blocks->{$scaffold}{$start};
    return $rejects if $block->{type} eq "Block";

    my $sth = $dbh->prepare(
"SELECT * FROM markers where scaffold=\"$scaffold\" and position >= $start and position <= $block->{end} and marker_type='Reject' order by position"
    );
    $rejects = $dbh->selectall_arrayref( $sth, { Slice => {} } );
    $sth->finish;

    print "$scaffold\t$block->{type}\t$start\t$block->{end}\t$block->{length}\n";

    $rejects;

}

sub get_neighbours {
    my ( $start_i, $scaffold_block, $block_starts, $chromosomes ) = @_;

    my $empty    = ' ' x 69;
    my $previous = $empty;
    my $next     = $empty;
    if ( $start_i > 0 ) {
        my $previous_block = $scaffold_block->{ $block_starts->[ $start_i - 1 ] };
        $previous = $chromosomes->{ $previous_block->{chromosome} }{markers}{ $previous_block->{cM} }{original};
    }
    if ( $start_i < $#{$block_starts} ) {
        my $next_block = $scaffold_block->{ $block_starts->[ $start_i + 1 ] };
        $next = $chromosomes->{ $next_block->{chromosome} }{markers}{ $next_block->{cM} }{original};
    }

    ( $previous, $next );
}

sub match_reject {
    my ( $reject, $pattern ) = @_;
    my $match = 0;
    return $match if $reject =~ /[01]/;
    return $match if $reject =~ /^[\.H]+$/;

    $reject =~ tr/H/B/ if ( $reject =~ /^[AH]+$/ );

    $match = 1 if match( $reject, $pattern );
    $match = 1 if match( $reject, mirror($pattern) );

    $match;
}

sub match {
    my ( $a, $b ) = @_;
    my @a = split //, $a;
    my @b = split //, $b;

    my $match = 1;
    for my $i ( 0 .. $#a ) {
        next if $a[$i] eq '.' or $b[$i] eq '.';
        next if $a[$i] eq 'H' or $b[$i] eq 'H';
        $match = 0 if $a[$i] ne $b[$i];
    }

    $match;
}

sub mirror {
    my $pat = shift;
    $pat =~ tr/AB/BA/;
    $pat;
}
