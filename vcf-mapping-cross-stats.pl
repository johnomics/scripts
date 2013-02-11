#!/usr/bin/env perl

# vcf-mapping-cross-stats.pl
#
# Input: one VCF file containing a whole mapping cross
# Output: statistics about numbers of missing individuals, mapping qualities etc
# Author: John Davey john.davey@ed.ac.uk
# Begun 01/02/13 based on convert_vcf_to_joinmap.pl

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
use Pod::Usage;
use Data::Dumper;

# Autoflush output so reporting on progress works
$| = 1;

my $vcf_filename     = "";
my $lengths_filename = "Hmel1-1_primaryScaffolds.lengths";
my $max_pos          = 0;

my $options_okay = GetOptions(
    'vcf=s'     => \$vcf_filename,
    'lengths=s' => \$lengths_filename,
    'max_pos=i' => \$max_pos,
);
croak "No VCF file! Please specify -v $OS_ERROR\n" if ( $vcf_filename eq "" );

my %scflen;
my $genome_length = 0;
open my $lengths_file, '<', $lengths_filename
  or croak "Can't open $lengths_filename $OS_ERROR!\n";
while (my $scf_line = <$lengths_file>) {
    chomp $scf_line;
    my ($scf, $length) = split "\t", $scf_line;
    $scflen{$scf}=$length;
    $genome_length += $length;
}
close $lengths_file;

open my $vcf_file, '<', $vcf_filename
  or croak "Can't open $vcf_filename $OS_ERROR!\n";

my @sample_names;

my $base_count;

my %position_missing;
my %sample_stat;
my %stat_range;
while ( my $vcf_line = <$vcf_file> ) {
    chomp $vcf_line;
    if ( $vcf_line =~ /^#CHROM/ ) {
        @sample_names = split /\t/, $vcf_line;
        for my $i ( 0 .. 8 ) { shift @sample_names; }
    }
    next if ( $vcf_line =~ /^#/ );

    $base_count++;
    if ( $base_count % 10000 == 0 )   { print STDERR "." }
    if ( $base_count % 100000 == 0 )  { printf STDERR "%8d", $base_count; }
    if ( $base_count % 1000000 == 0 ) { print STDERR "\n"; }
    last if ( ( $max_pos > 0 ) && ( $base_count > $max_pos ) );

    my @fields = split /\t/, $vcf_line;

    my %position;
    my @format_fields = split /:/, $fields[8];
    my $gt_field_num  = -1;
    my $gq_field_num  = -1;
    my $dp_field_num  = -1;
    for my $i ( 0 .. ( @format_fields - 1 ) ) {
        if ( $format_fields[$i] eq "GT" ) { $gt_field_num = $i; }
        if ( $format_fields[$i] eq "GQ" ) { $gq_field_num = $i; }
        if ( $format_fields[$i] eq "DP" ) { $dp_field_num = $i; }
    }

    my $missing = 0;
    for my $i ( 9 .. ( @sample_names + 8 ) ) {
        if ( $fields[$i] eq "./." ) {
            $missing++;
            $sample_stat{dp}{$sample_names[$i-9]}{0}++;
            $sample_stat{gq}{$sample_names[$i-9]}{0}++;
        }
        else {
            my @sample_f = split /:/, $fields[$i];
            $sample_stat{dp}{$sample_names[$i-9]}{$sample_f[$dp_field_num]}++;
            $sample_stat{gq}{$sample_names[$i-9]}{$sample_f[$gq_field_num]}++;
            $stat_range{dp}{$sample_f[$dp_field_num]}++;
            $stat_range{gq}{$sample_f[$gq_field_num]}++;
        }
    }
    $position_missing{$missing}{$fields[0]}++;
}

close $vcf_file;

print "\n";



foreach my $stat ("dp", "gq") {
    print "Ind\\$stat";
    foreach my $stat_val (sort {$a<=>$b} keys %{$stat_range{$stat}}) {
        print "\t$stat_val";
    }
    print "\n";

    foreach my $sample (sort {$a<=>$b} keys %{$sample_stat{$stat}}) {
        print "$sample";
        foreach my $stat_val (sort {$a<=>$b} keys %{$stat_range{$stat}}) {
            if (defined $sample_stat{$stat}{$sample}{$stat_val}) {
                print "\t$sample_stat{$stat}{$sample}{$stat_val}";
            }
            else {
                print "\t0";
            }
        }
        print "\n";
    }
    print "\n\n";
}

foreach my $missing (sort {$a<=>$b} keys %position_missing) {
    my $missing_bases = 0;
    foreach my $scf (keys %{$position_missing{$missing}}) {
        $missing_bases += $position_missing{$missing}{$scf};
    }
    my %cumul_scf;
    foreach my $cumul_miss (0..$missing) {
        foreach my $scf (keys %{$position_missing{$cumul_miss}}) {
            $cumul_scf{$scf}++;
        }
    }
    my $base_coverage = 0;
    foreach my $scf (keys %cumul_scf) {
        $base_coverage += $scflen{$scf};
    }
    my $pc_base_coverage = sprintf "%5.2f", $base_coverage / $genome_length * 100;
    my $cumul_scf_num = keys %cumul_scf;
    my $pc_scfs = sprintf "%5.2f", $cumul_scf_num / keys(%scflen) * 100;

    print "$missing\t$missing_bases\t$cumul_scf_num\t$pc_scfs\t$base_coverage\t$pc_base_coverage\n";
}
