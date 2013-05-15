#!/usr/bin/env perl

# partition_reference_into_intervals.pl
#
# Takes a FASTA reference file and a number of jobs, or a unit size in Mbp,
# and writes out a set of intervals files 1.intervals, 2.intervals,...
# ready for use with GATK -L option on SGE
#
# Author: John Davey john.davey@ed.ac.uk
# Begun 23/07/11

#############################################################################
###                                                                       ###
###                                 CODE                                  ###
###                                                                       ###
#############################################################################

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Data::Dumper;
use File::Basename;

# Autoflush output so reporting on progress works
$| = 1;

my $ref_filename = "";
my $unit         = 20;
my $options_okay = GetOptions(
    'reference=s' => \$ref_filename,
    'unit=f'      => \$unit,
);

croak "\nUsage: perl partition_reference_into_intervals.pl -r reference -u unit (Mbp)\n"
  if !$options_okay;

croak "Please specify a reference with -r" if ( $ref_filename eq "" );

my ( $refbasename, $refdir, $refsuf ) =
  fileparse( $ref_filename, ".fasta", ".fas", ".fa", ".fna" );

# Create and clean up intervals directory for this reference
if ( !( -e "$refbasename\.intervals" ) ) {
    mkdir("$refbasename\.intervals")
      or croak "Can't create $refbasename\.intervals directory: $OS_ERROR\n";
}
unlink glob "$refbasename\.intervals/*.intervals";

open my $ref_file, '<', $ref_filename
  or croak "Can't open $ref_filename: $OS_ERROR!\n";

my $threshold_switch = 0;
my $unit_size        = 0;
my $scaffold_num     = 0;
my $scaffold_size    = 0;
my $full             = 0;
my $first            = 1;

my %job_file;
my $job_id = 1;
open $job_file{$job_id}, '>', "$refbasename\.intervals/$job_id.intervals"
  or croak
  "Can't open job file $refbasename\.intervals/$job_id.intervals: $OS_ERROR\n";

print STDERR "Job\tScfs\tBases\n";
while ( my $ref_line = <$ref_file> ) {
    chomp $ref_line;
    if ( $ref_line =~ /^>(.+)$/ ) {    # header line
        if (!$first) {print {$job_file{$job_id}} ":1-$scaffold_size\n";}
        $first = 0;
        $scaffold_size    = 0;
        if ($full) {
            print STDERR "$job_id\t$scaffold_num\t$unit_size\n";
            close $job_file{$job_id};
            $job_id++;
            open $job_file{$job_id}, '>',
              "$refbasename\.intervals/$job_id.intervals"
              or croak
"Can't open job file $refbasename\.intervals/$job_id.intervals: $OS_ERROR\n";
            $scaffold_num     = 0;
            $threshold_switch = 0;
            $unit_size        = 0;
            $full             = 0;
        }
        my $scaffold_name = ( split / /, $1 )[0];
        print { $job_file{$job_id} } "$scaffold_name";
        $scaffold_num++;
    }
    else {    # Sequence line
        $scaffold_size += length($ref_line);
        $unit_size += length($ref_line);
        next if ($full);
        $threshold_switch += length($ref_line);
        if ( $threshold_switch > ( $unit * 1_000_000 ) ) {
            $full = 1;
        }
    }
}
print {$job_file{$job_id}} ":1-$scaffold_size\n";
print STDERR "$job_id\t$scaffold_num\t$unit_size\n";
close $job_file{$job_id};


close $ref_file;
