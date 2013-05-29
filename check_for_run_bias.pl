#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

use Term::ExtendedColor qw/:all/;

$OUTPUT_AUTOFLUSH = 1;

my %args;

$args{block_file} = "";
$args{maxscf} = 0;
my $options_okay = GetOptions( 'block=s' => \$args{block_file}, 'maxscf=i' => \$args{maxscf});
croak "No block file! Please specify -b $OS_ERROR\n"
  if ( $args{block_file} eq "" );

open my $blockfile, '<', $args{block_file}
  or croak "Can't open $args{block_file}! $OS_ERROR\n";

my $prev_scf;
my $scf_count = 0;
my @scflines;
my @block;
my $final_scf;
while ( my $marker = <$blockfile> ) {
    next if ( $marker =~ /^(-+)$/ );
    chomp $marker;
    my ( $scf, $snp, $type, $mq, $fs, $orig, $edge, $cons, $corrected ) =
      split /\t/,
      uncolor($marker);

    next if ( $type ne 'maternal' );

    if ( !defined $prev_scf ) {
        $prev_scf = $scf;
    }
    if ( $prev_scf ne $scf ) {
        process_blocks( $prev_scf, \@scflines, \@block );
        @scflines  = ();
        $prev_scf  = $scf;
        $final_scf = $scf;

        $scf_count++;
        print STDERR '.'            if ( $scf_count % 10 == 0 );
        print STDERR "$scf_count\n" if ( $scf_count % 100 == 0 );
        last if ($args{maxscf} > 0 and $args{maxscf} <= $scf_count);
    }
    push @scflines, { pos => $snp, calls => $cons };
}

process_blocks( $final_scf, \@scflines, \@block ) if (@scflines > 0);

close $blockfile;


my %blocksum;
foreach my $sample (0..$#block) {
    foreach my $scf ( keys %{ $block[$sample] } ) {
        foreach my $pos1 ( sort { $a <=> $b } keys %{ $block[$sample]{$scf} } )
        {
            foreach my $pos2 (
                sort { $a <=> $b }
                keys %{ $block[$sample]{$scf}{$pos1} }
              )
            {
                my $blocksize = $pos2 - $pos1 + 1;
                $blocksum{$blocksize}{$sample} +=
                  $block[$sample]{$scf}{$pos1}{$pos2};
            }
        }
    }
}

print "Size";
map { print "\t$_" } ( 0 .. 68 );
print "\n";

foreach my $blocksize ( sort { $a <=> $b } keys %blocksum ) {
    print $blocksize;
    foreach my $blocksamples ( 0 .. 68 ) {
        print "\t";
        print $blocksum{$blocksize}{$blocksamples} // "0";
    }
    print "\n";
}

sub process_blocks {
    my ( $scfname, $scf, $block ) = @_;

    my @start;
    my @prev_call;
    my $prev_pos;
    my @calls;
    for my $marker ( @{$scf} ) {
        @calls = split //, $marker->{calls};

        for my $i ( 0 .. $#calls ) {
            if ( !defined $start[$i] ) {
                $start[$i] = $marker->{pos};
            }

            if ( defined $prev_call[$i] and $prev_call[$i] ne $calls[$i] ) {
                $block->[$i]{$scfname}{ $start[$i] }{$prev_pos}++;
                $start[$i] = $marker->{pos};
            }
            $prev_call[$i] = $calls[$i];
        }
        $prev_pos = $marker->{pos};
    }
    for my $i ( 0 .. $#calls ) {
        $block->[$i]{$scfname}{ $start[$i] }{$prev_pos}++;
    }
}
