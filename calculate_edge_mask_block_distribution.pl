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

my $options_okay = GetOptions( 'block=s' => \$args{block_file}, );
croak "No block file! Please specify -b $OS_ERROR\n"
  if ( $args{block_file} eq "" );

open my $blockfile, '<', $args{block_file}
  or croak "Can't open $args{block_file}! $OS_ERROR\n";

my $prev_scf;
my $scf_count = 0;
my @scflines;
my %block;
my $final_scf;
while ( my $marker = <$blockfile> ) {
    next if ( $marker =~ /^(-+)$/ );
    chomp $marker;
    my ( $scf, $snp, $type, $orig, $edge, $cons, $corrected ) = split /\t/,
      uncolor($marker);

    next if ( $type eq 'paternal' );

    if ( !defined $prev_scf ) {
        $prev_scf = $scf;
    }
    if ( $prev_scf ne $scf ) {
        process_blocks( $prev_scf, \@scflines, \%block );
        @scflines = ();
        $prev_scf = $scf;
        $final_scf = $scf;

        $scf_count++;
        print STDERR '.'            if ( $scf_count % 10 == 0 );
        print STDERR "$scf_count\n" if ( $scf_count % 100 == 0 );
    }
    push @scflines, { pos => $snp, calls => $corrected };
}

process_blocks( $final_scf, \@scflines, \%block );

close $blockfile;

my %blocksum;
foreach my $scf (keys %block) {
    foreach my $pos1 (sort {$a<=>$b} keys %{$block{$scf}}) {
        foreach my $pos2 (sort {$a<=>$b} keys %{$block{$scf}{$pos1}}) {
            my $blocksize = $pos2-$pos1+1;
            $blocksum{$blocksize}{$block{$scf}{$pos1}{$pos2}}++;
        }
    }
}

print "Size";
map {print "\t$_"} (1..69);
print "\n";

foreach my $blocksize (sort {$a<=>$b} keys %blocksum) {
    print $blocksize;
    foreach my $blocksamples (1..69) {
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

            if (defined $prev_call[$i] and $prev_call[$i] ne $calls[$i]) {
                $block->{$scfname}{$start[$i]}{$prev_pos}++;
                $start[$i] = $marker->{pos};
            }
            $prev_call[$i] = $calls[$i];
        }
        $prev_pos = $marker->{pos};
    }
    for my $i (0..$#calls) {
        $block->{$scfname}{$start[$i]}{$prev_pos}++;
    }
}

sub output_block {
    my ($block) = @_;
    
    foreach my $scf (sort keys %{$block}) {
        foreach my $pos1 (sort {$a<=>$b} keys %{$block->{$scf}}) {
            foreach my $pos2 (sort {$a<=>$b} keys %{$block->{$scf}{$pos1}}) {
                print "$scf\t$pos1\t$pos2\t$block->{$scf}{$pos1}{$pos2}\n";
            }
        }
    }
}