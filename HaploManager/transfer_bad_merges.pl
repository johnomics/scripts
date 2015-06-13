#!/usr/bin/env perl

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

my $input  = "";
my $output = "";
my $tsv    = "";

my $options_okay = GetOptions(
    'input=s'  => \$input,
    'output=s' => \$output,
    'tsv=s'    => \$tsv,
);

croak "No input filename!"        if $input eq "";
croak "No output filename!"       if $output eq "";
croak "No TSV filename!"          if $tsv eq "";
croak "TSV file doesn't exist!"   if !-e $tsv;
croak "Input file doesn't exist!" if !-e $input;

open my $tsvfh, '<', $tsv or croak "Can't open TSV file!\n";

my %transfers;
while ( my $tsvline = <$tsvfh> ) {
    chomp $tsvline;
    my ( $scf1, $start1, $end1, $scf2, $start2, $end2, $strand, $comment ) = split "\t", $tsvline;
    $transfers{$scf1}{$scf2}++;
}
close $tsvfh;

open my $newbmfh, '>', $output or croak "Can't open new bad merges file!\n";

open my $bmfh, '<', $input or croak "Can't open badmerges file $input!\n";

my %output_dict;
while ( my $bmline = <$bmfh> ) {
    chomp $bmline;
    next if $bmline =~ /^$/;
    my ( $bm1, $bm2 ) = split "\t", $bmline;
    my @bm1_transfers = defined $transfers{$bm1} ? keys %{$transfers{$bm1}} : ();
    my @bm2_transfers = defined $transfers{$bm2} ? keys %{$transfers{$bm2}} : ();

    for my $bm1t (@bm1_transfers) {
        for my $bm2t (@bm2_transfers) {
            next if $bm1t eq $bm2t;
            my $newmerge = "$bm1t\t$bm2t\n";
            next if defined $output_dict{$newmerge};
            print $newbmfh $newmerge;
            $output_dict{$newmerge}++;
        }
    }
}

close $bmfh;
close $newbmfh;