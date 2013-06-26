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

my %args;
$args{blocks} = "";

my $options_okay = GetOptions( 'blocks=s' => \$args{blocks}, );

my ( $header, $types, $blocklist ) = load_blocks( $args{blocks} );

my $markers = correct_blocks( $header, $types, $blocklist );

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

sub correct_blocks {
    my ( $header, $types, $blocklist ) = @_;
    correct_maternal( "Maternal-AHAH", $blocklist );
    correct_type_pair( "Intercross-ABHABH_HHA", "Paternal-AHAH", $blocklist );
    correct_type_pair( "Paternal-AHAH", "Intercross-ABHABH_HHA", $blocklist );
    correct_type_pair( "Intercross-ABHABH_HHA", "Intercross-ABHABH_HHH",
        $blocklist );
    correct_type_pair( "Intercross-ABHABH_HHH", "Intercross-ABHABH_HHA",
        $blocklist );
    correct_type_pair( "Intercross-ABHABH_HHA", "Maternal-AHAH", $blocklist );
    correct_type_pair( "Intercross-ABHABH_HHH", "Maternal-AHAH", $blocklist );
    correct_type_pair( "Paternal-AHAH",         "Maternal-AHAH", $blocklist );
    correct_type_pair( "Paternal-AHAB_AHB",     "Sex-HB",        $blocklist );
}

sub correct_maternal {
    my ( $mattype, $blocklist ) = @_;

    my %mat;
    for my $i ( 0 .. $#{$blocklist} ) {
        $mat{ $blocklist->[$i]{$mattype} }{length} +=
          $blocklist->[$i]{'Length'};
        push @{ $mat{ $blocklist->[$i]{$mattype} }{blocks} }, $i;
    }
    my @matl = reverse sort {$mat{$a}{length} <=> $mat{$b}{length}} keys %mat;
    my %merged;
    for my $mata ( @matl ) {
        next if defined $merged{$mata};
        next if $mata eq $empty;
        next if $mata =~ /[~\.]/;
        for my $matb ( @matl ) {
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

sub correct_type_pair {
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
    return $hamming <= 3;
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
