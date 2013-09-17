#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Parallel::ForkManager;
use List::Util qw/sum/;
use Data::Dumper;

my $options_okay = GetOptions();

my @popnames = ("1", "2", "3","999");
my @pops = ( 4, 4, 4, 4 );

my @popnum;
for my $i ( 1 .. @pops ) {
    for my $j ( 1 .. $pops[ $i - 1 ] ) {
        push @popnum, $i;
    }
}

my $nsam    = sum @pops;
my $npops   = @pops;
my $howmany = 1000;
my $seqlen  = 5000;

my $mutation = 2.5e-9;
my $n = 1000000;
my $theta = 4 * $mutation * $n;

my $t_gf    = 0.2;      # time of modern admixture event
my $t_h     = 0.4;       # time of 12 split
my $t_n     = 6.8;       # time of 123 split
my $t_o     = 12;         # time of outgroup split
my $t_s     = 10;         # time of ancient population divergence
my $f       = 0.1;       # % of gene flow
my $fournm  = 0.4;       # modern migration between 1 and 2 M_12 = 4*N0*0.1
my $fournma = 0.4;       # ancient migration between 1 and 2

my %model;
my $iarg = "";
map { $iarg .= "$_ " } @pops;
$model{common} = "-t $theta ";
$model{common} .= "-I $npops $iarg ";    # Create populations
$model{common} .=
  "-m 1 2 $fournm -m 2 1 $fournm ";    # Populations 1 and 2 have gene flow
$model{common} .= "-ej $t_o 1 4 ";      # 1 splits from outgroup 4 at time t_o

$model{noflow} = "-ej $t_n 3 1 ";    # 1 and 3 split at t_n
$model{noflow} .= "-ej $t_h 2 1 ";    # 1 and 2 split at t_h

$model{admixture} = $model{noflow};
$model{admixture} .= "-es $t_gf 2 1-$f -ej $t_gf 5 3 "
  ;    # f% flows from pop 2 to pop 3 (via new pop 5)
$model{admixture} .= "-es $t_gf 3 1-$f -ej $t_gf 5 2 "
    ;    # f% flows from pop 3 to pop 2 (via new pop 5)


$model{structure} = "-ej $t_s 2 1 ";    # 1 and 2 split at t_s
$model{structure} .= "-ej $t_n 3 1 ";      # 1 and 3 split at t_n
$model{structure} .= "-em $t_h 2 1 $fournma -em $t_h 1 2 $fournma "
  ;                                      # ancient migration between 1 and 2

for my $h ( "noflow", "admixture", "structure" ) {
    print "Running $h model...\n";
    system "ms $nsam $howmany -T $model{common} $model{$h} | grep ^\\( > $h.ms";
    system "seq-gen -mHKY -l $seqlen $h.ms > $h.sg 2> $h.sg.log";

    open my $sgout, "<", "$h.sg"
      or croak "Can't open seq-gen output!\n";

    open my $windows, ">", "$h.windows.tsv"
      or croak "Can't open window output!\n";
    print $windows "Window\tPopulation\tSample\tSequence\n";

    my $window = 0;
    while ( my $header = <$sgout> ) {
        chomp $header;
        my ( $null, $seqs, $len ) = split /\s/, $header;
        $window++;
        my %windowseq;
        for my $seq ( 1 .. $seqs ) {
            my $seqline = <$sgout>;
            chomp $seqline;
            my ( $seqnum, $seqbases ) = split /\s+/, $seqline;
            $windowseq{$seqnum} = $seqbases;
        }

        for my $sam ( sort { $a <=> $b } keys %windowseq ) {
            print $windows "$window\t$popnum[$sam-1]\t$sam\t$windowseq{$sam}\n";
        }
    }
    close $windows;
    close $sgout;
}
