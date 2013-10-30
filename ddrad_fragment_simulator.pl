#!/usr/bin/env perl

# ddrad_fragment_simulator.pl
#
# Input: reference genome and two restriction enzymes
# Output: restriction fragments with start, end, length and GC statistics
# Author: John Davey johnomics@gmail.com

# Based on simulate_rad_fragments.pl for Molecular Ecology RAD paper

# Modifications begun 30 October 2013

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
use Parallel::ForkManager;

# Autoflush output so reporting on progress works
$| = 1;

use constant ILLUMINA_MIN => 200;
use constant ILLUMINA_MAX => 700;

my $in_filename = "";
my $p1enzyme = "CTGCAG"; # PstI
my $p2enzyme = "TTAA";  # MseI
my $filep2   = "";
my $threads = 1;
my $options_okay = GetOptions(
    'input=s'      => \$in_filename,
    'left=s'       => \$p1enzyme,
    'right=s'      => \$p2enzyme,
    'filep2=s'     => \$filep2,
    'threads=i'    => \$threads,
);

croak
"\nUsage: perl simulate_rad_sites.pl -i fastq_file -l p1enzyme -r p2enzyme\n"
  if !$options_okay;

croak
"\nPlease specify an input file with -i\nUsage: perl simulate_rad_sites.pl -i fastq_file -l p1enzyme -r p2enzyme\n"
  if $in_filename eq "";

croak
"\nPlease specify restriction enzyme sites with -l and -r\nUsage: perl simulate_rad_sites.pl -i fastq_file -l p1enzyme -r p2enzyme\n"
  if $p1enzyme eq "" or $p2enzyme eq "";

my %enzymes;
if ($filep2 ne "") {
    open my $p2file, "<", $filep2 or croak "Can't open enzyme file! $OS_ERROR\n";

    while (my $enzyme_line = <$p2file>) {
        if ($enzyme_line =~ /^(.+?)(\s+)(.+)$/) {
            my $name = $1;
            my $site = $3;
            $site =~ s/\^//g;
            $enzymes{$name} = $site;
        }
    }
    close $p2file;
}
else {
    $enzymes{p2} = $p2enzyme;
}

my %genome;
open my $in_file, '<', $in_filename or croak "Can't open input file $in_filename: $OS_ERROR!\n";

my $seq = "";
my $chrom = "";
my %lengths;
while (my $in_line = <$in_file>) {
    chomp $in_line;
    if ($in_line =~ />(.+) /) {
        $genome{$chrom} = $seq if $seq ne "";
        $chrom = $1;
        $seq = "";
    }
    else {
        $seq .= uc $in_line;
    }
}
$genome{$chrom} = $seq;
close $in_file;


print "Genome\tP1Name\tP1Seq\tP1Sites\tP1Fragments\tP2Name\tP2Seq\tP1P2Fragments\tP1P2Pippin\tP1P2Illumina\tP1SitesP2Cut\tP1SitesP2CutProp\tP1SitesP2CutPippin\tP1SitesP2CutPippinProp\tP1SitesP2CutIllumina\tP1SitesP2CutIlluminaProp\n";
my $pm = Parallel::ForkManager->new($threads);

$pm->run_on_finish (
    sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
        print ${$data_ref} if (defined $data_ref);
    }
);

for my $enzyme (keys %enzymes) {
    next if $enzymes{$enzyme} =~ /$p1enzyme/;

    my $pid = $pm->start and next;
    my $frags = 0;
    my $illumina = 0;
    my $pippin = 0;

    my %lengths;
    my $p1count = 0;
    my $p1p2count = 0;
    my $p1p2illcount = 0;
    my $p1p2pippincount = 0;
    for my $chrom (keys %genome) {
        my ($chrp1count, $chrp1p2count, $chrp1p2illcount,$chrp1p2pippincount) = process_chrom($genome{$chrom},\%lengths,$p1enzyme,$enzymes{$enzyme});
        $p1count += $chrp1count;
        $p1p2count += $chrp1p2count;
        $p1p2illcount += $chrp1p2illcount;
        $p1p2pippincount += $chrp1p2pippincount;
    }
    
    my $p1p2prop = sprintf "%5.2f", $p1p2count/$p1count*100;
    my $p1p2illprop = sprintf "%5.2f", $p1p2illcount/$p1count*100;
    my $p1p2pippinprop = sprintf "%5.2f", $p1p2pippincount/$p1count*100;
    
    for my $len (sort {$a<=>$b} keys %lengths) {
        $frags += $lengths{$len};
        $illumina += $lengths{$len} if $len >= ILLUMINA_MIN && $len <= ILLUMINA_MAX;
        $pippin += $lengths{$len} if $len >= ILLUMINA_MIN;
    }
    my $p1frags = $p1count * 2;
    my $outstr = "$in_filename\tPstI\t$p1enzyme\t$p1count\t$p1frags\t$enzyme\t$enzymes{$enzyme}\t$frags\t$pippin\t$illumina\t$p1p2count\t$p1p2prop\t$p1p2pippincount\t$p1p2pippinprop\t$p1p2illcount\t$p1p2illprop\n";
    
    $pm->finish(0, \$outstr);
}

$pm->wait_all_children;


sub process_chrom {
    my ($seq, $len, $p1enzyme, $p2enzyme) = @_;
    my @frags = split /($p1enzyme)/, $seq;
    my $p1count = () = $seq =~ /$p1enzyme/g;
    my $p1p2count = 0;
    my $p1p2illcount = 0;
    my $p1p2pippincount = 0;
    for my $f (0..$#frags) {
        my $found = 0;
        my $ill_found = 0;
        my $pippin_found = 0;
        if ($frags[$f] eq $p1enzyme) {
            if ($frags[$f-1] =~ /(.*)$p2enzyme(.+)$/) {
                my $l = length $2;
                $len->{$l}++;
                $found++;
                $ill_found++ if $l >= ILLUMINA_MIN && $l <= ILLUMINA_MAX;
                $pippin_found++ if $l >= ILLUMINA_MIN;
            }
            if ($frags[$f+1] =~ /^(.+?)$p2enzyme/) {
                my $l = length $1;
                $len->{$l}++;
                $found++;
                $ill_found++ if $l >= ILLUMINA_MIN && $l <= ILLUMINA_MAX;
                $pippin_found++ if $l >= ILLUMINA_MIN;
            }
        }
        $p1p2count++ if $found;
        $p1p2illcount++ if $ill_found;
        $p1p2pippincount++ if $pippin_found;
    }
    ($p1count,$p1p2count,$p1p2illcount,$p1p2pippincount);
}