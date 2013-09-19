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

my @popnames = ( "1", "2", "3", "4", "4", "4", "4" );
my @pops     = ( 8,   8,   8,   2,   2,   2,   2 );

my @popnum;
for my $i ( 0 .. $#pops ) {
    for my $j ( 1 .. $pops[$i] ) {
        push @popnum, $popnames[$i];
    }
}

my $nsam    = sum @pops;
my $npops   = @pops;
my $howmany = 1000;
my $seqlen  = 5000;

my $mutation = 2.5e-9;
my $n        = 1000000;
my $theta    = 4 * $mutation * $n;

my $t_gf    = 0.02;     # time of modern admixture event
my $t_h     = 0.8;      # time of 12 split
my $t_n     = 1.7;      # time of 123 split
my $t_o     = 3;        # time of outgroup split
my $t_s     = 2;        # time of ancient population divergence
my $f       = 0.1;      # proportion of gene flow
my $fournm  = 40000;    # modern migration between 1 and 2
my $fournma = 40000;    # ancient migration between 1 and 2

my %model;
my $iarg = "";
map { $iarg .= "$_ " } @pops;
$model{common} = "-t $theta ";
$model{common} = "-I $npops $iarg ";    # Create populations
$model{common} .=
  "-m 1 2 $fournm -m 2 1 $fournm ";     # Populations 1 and 2 have gene flow
$model{common} .= "-ej 1 5 4 ";
$model{common} .= "-ej 1.4 6 4 ";
$model{common} .= "-ej $t_n 7 4 ";
$model{common} .= "-ej $t_o 4 1 ";      # 1 splits from outgroup 4 at time t_o

$model{noflow} = "-ej $t_n 3 1 ";       # 1 and 3 split at t_n
$model{noflow} .= "-ej $t_h 2 1 ";      # 1 and 2 split at t_h

$model{admixture} = $model{noflow};
$model{admixture} .=
    "-es $t_gf 2 1-$f -ej $t_gf "
  . ( @pops + 1 )
  . " 3 ";    # f% flows from pop 2 to pop 3 (via new pop)
$model{admixture} .=
    "-es $t_gf 3 1-$f -ej $t_gf "
  . ( @pops + 1 )
  . " 2 ";    # f% flows from pop 3 to pop 2 (via new pop)

$model{structure} = "-ej $t_s 2 1 ";    # 1 and 2 split at t_s
$model{structure} .= "-ej $t_n 3 1 ";   # 1 and 3 split at t_n
$model{structure} .= "-em $t_h 2 1 $fournma -em $t_h 1 2 $fournma "
  ;                                     # ancient migration between 1 and 2

my $model_pm = new Parallel::ForkManager(3);
$model_pm->set_max_procs(3);

my $modelnum = 0;
for my $h ( "noflow", "admixture", "structure" ) {
    $modelnum++;
    $model_pm->start($modelnum) and next;

    print "Running $h model...\n";
    system "ms $nsam $howmany -T $model{common} $model{$h} | grep ^\\( > $h.ms";

    for my $s ( 0.1, 0.3, 0.5, 0.7 ) {
        system "seq-gen -mHKY -l $seqlen -s $s $h.ms > $h.s$s.sg 2> $h.s$s.sg.log";

        open my $sgout, "<", "$h.s$s.sg"
          or croak "Can't open seq-gen output!\n";

        open my $windows, ">", "$h.s$s.windows.tsv"
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
                print $windows
                  "$window\t$popnum[$sam-1]\t$sam\t$windowseq{$sam}\n";
            }
        }
        close $windows;
        close $sgout;
    }
    $model_pm->finish();
}

$model_pm->wait_all_children;
