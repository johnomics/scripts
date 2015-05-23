#!/usr/bin/env perl

# batchhm.pl
# Run HaploMerger until pure haploid genome is produced

# -f FASTA file of genome, zipped or unzipped
# -d Directory containing HaploMerger scripts and config files:
#    - hm.batchA-G
#    - scoreMatrix.q
#    - all_lastz.ctl
# -o Name of output directory
# -p final scaffold prefix

# John Davey
# johnomics@gmail.com
# Begun 6 April 2015

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use File::Basename 'fileparse';
use File::Copy 'copy';
use File::Copy::Recursive 'dircopy';
use File::Path 'rmtree';
use IO::Uncompress::Gunzip qw($GunzipError);
use Cwd;

# Autoflush output so reporting on progress works
$| = 1;

my $input           = "";
my $configdir       = "";
my $outputdir       = "";
my $prefix          = "";
my $scaffold_prefix = "";
my $badmerges       = "";
my $tsv             = "";
my $rounds          = 0;
my $g;

my $options_okay = GetOptions(
    'input=s'           => \$input,
    'configdir=s'       => \$configdir,
    'outputdir=s'       => \$outputdir,
    'prefix=s'          => \$prefix,
    'scaffold_prefix=s' => \$scaffold_prefix,
    'badmerges=s'       => \$badmerges,
    'tsv=s'             => \$tsv,
    'rounds=i'          => \$rounds,
    'g'                 => \$g,
);

croak "Can't open badmerges file $badmerges!" if $badmerges and !-e $badmerges;
croak "No TSV file!" if $tsv eq "";

my $i = 1;
my $command = "runhm.pl -i $input -c $configdir -o $outputdir\_$i -p $prefix$i -g";
$command .= " -s $scaffold_prefix" if $scaffold_prefix;
$command .= " -b $badmerges" if $badmerges;
system $command;
system "map_merge.py -m $outputdir\_$i -p $prefix$i";
system "transfer_merge.py -d $tsv -n $outputdir\_$i/genome.genomex.result/$outputdir\_$i\_old.tsv -o $outputdir\_$i\_orig.tsv";

while (-e "$outputdir\_$i/genome.genomex.result/$outputdir\_$i\_refined.fa" and -s "$outputdir\_$i/genome.genomex.result/$outputdir\_$i\_refined.fa" > 0) {
    my $j = $i+1;
    my $filename = "$outputdir\_$i/genome.genomex.result/$outputdir\_$i";
    $filename .= $g ? "_refined.fa" : ".fa";

    if ($badmerges) {
        my $newbadmerges = "$outputdir\_$i/bad_merges.txt";
        system "transfer_bad_merges.pl -i $badmerges -o $newbadmerges -t $outputdir\_$i/genome.genomex.result/$outputdir\_$i\_old.tsv";
        $badmerges = $newbadmerges;
    }

    my $command = "runhm.pl -i $filename -c $configdir -o $outputdir\_$j -p $prefix$j -g";
    $command .= " -s $scaffold_prefix" if $scaffold_prefix;
    $command .= " -b $badmerges" if $badmerges;

    system $command;
    system "map_merge.py -m $outputdir\_$j -p $prefix$j";
    system "transfer_merge.py -d $outputdir\_$i\_orig.tsv -n $outputdir\_$j/genome.genomex.result/$outputdir\_$j\_old.tsv -o $outputdir\_$j\_orig.tsv";

    $i++;
    
    last if $rounds > 0 and $rounds < $i;
}


