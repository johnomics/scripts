#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Memoize;

use Parallel::ForkManager;
use Term::ExtendedColor qw/:all/;
use List::Util qw/min max/;

$OUTPUT_AUTOFLUSH = 1;

open my $vcf, '<', "/whale-data/jd626/edinburgh/isilon/heliconius/scaffolding/PstI/PstI.Hmel_primaryScaffs_NoMito_Yb-scaff.pass1.vcf" or croak "Can't open VCF file";

my %markers;
while (my $vcfline = <$vcf>) {
    if ($vcfline =~ /^scf7180001249824/) {
        chomp $vcfline;
        my @F = split /\t/, $vcfline;
        my $scf = $F[0];
        my $pos = $F[1];
        my $f0gm = $F[-3];
        my $f1mo = $F[-2];
        my $f1fa = $F[-1];
        my @pgts = map {(split /:/, $_)[0]} ($f0gm, $f1mo, $f1fa);
        my %pgts;
        map {$pgts{$_}++} @pgts;
        next if defined $pgts{'./.'};
        next if keys %pgts == 1;
        my $pgt = "@pgts";
        my @ogts = map {(split /:/, $_)[0]} @F[9..($#F-3)];
        $markers{$pgt}{$pos}="@ogts";
    }
}
close $vcf;

for my $pgt (keys %markers) {
    for my $pos (sort {$a<=>$b} keys %{$markers{$pgt}}) {
        print gts_to_pattern($pgt);
        print "\t$pos\t";
        print gts_to_pattern($markers{$pgt}{$pos});
        print "\n";
    }
}


sub gts_to_pattern {
    return map {$_ eq '0/0' ? 'A' : $_ eq '0/1' ? 'H' : $_ eq '1/1' ? 'B' : '-'} split / /, shift;
}