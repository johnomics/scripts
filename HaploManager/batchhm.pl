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
my $g;

my $options_okay = GetOptions(
    'input=s'           => \$input,
    'configdir=s'       => \$configdir,
    'outputdir=s'       => \$outputdir,
    'prefix=s'          => \$prefix,
    'g'                 => \$g,
);

my $i = 1;
system("runhm.pl -i $input -c $configdir -o $outputdir\_$i -p $prefix$i -g");

while (-e "$outputdir\_$i/genome.genomex.result/$outputdir\_$i\_refined.fa" and -s "$outputdir\_$i/genome.genomex.result/$outputdir\_$i\_refined.fa" > 0) {
    my $j = $i+1;
    my $filename = "$outputdir\_$i/genome.genomex.result/$outputdir\_$i";
    $filename .= $g ? "_refined.fa" : ".fa";
    
    system("runhm.pl -i $filename -c $configdir -o $outputdir\_$j -p $prefix$j -g");
    $i++;
}