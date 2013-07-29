#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Parallel::ForkManager;
use Memoize;

memoize 'test_block';

my $alignment    = "";
my $reference    = "";
my $output_dir   = "aligner_test_out";
my $threads      = 1;
my $options_okay = GetOptions(
    'alignment=s' => \$alignment,
    'reference=s' => \$reference,
    'output=s'    => \$output_dir,
    'threads=i'   => \$threads,
);

croak
  "Please specify a file containing a gapped alignment in FASTA format with -a"
  if $alignment eq "";
croak
  "Please specify the name of the reference sequence in the alignment with -r"
  if $reference eq "";

my $seqs_ref = load_alignment($alignment);
croak "Reference must be present in alignment"
  if !defined $seqs_ref->{$reference};

-d $output_dir or mkdir $output_dir;

open my $refout, ">", "$output_dir/$reference\_ref.fa" or croak "Can't open output reference FASTA file!\n";
my $refoutseq = $seqs_ref->{$reference};
$refoutseq =~ s/-//g;
print $refout ">$reference\n$refoutseq\n";
close $refout;

my $pm = new Parallel::ForkManager($threads);
$pm->set_max_procs($threads);

for my $seqname ( keys %{$seqs_ref} ) {
    next if $seqname eq $reference;
    for my $readlen ( 50, 100 ) {
        for my $endtype ( 'single', 'paired' ) {
            next if $endtype eq 'paired';
            $pm->start and next;
            print "$seqname\t$readlen\t$endtype\n";
            generate_read_alignment(
                $seqname,  $readlen,   $endtype,
                $seqs_ref, $reference, $output_dir
            );
            $pm->finish;
        }
    }
}

$pm->wait_all_children;

sub generate_read_alignment {
    my ( $seqname, $readlen, $endtype, $seqs_ref, $reference, $output_dir ) =
      @_;
    open my $sam, ">", "$output_dir/$seqname.$readlen.$endtype.sam"
      or croak
      "Can't open SAM file for $seqname, $readlen bp reads, $endtype end\n";
    open my $fq, ">", "$output_dir/$seqname.$readlen.$endtype\_1.fq"
      or croak "Can't open read 1 FASTQ file for $seqname, $readlen bp reads\n";

    my $testseq = $seqs_ref->{$seqname};
    my $refseq  = $seqs_ref->{$reference};
    print STDERR
      "$seqname alignment is not the same length as reference alignment!\n"
      and next
      if length($testseq) != length($refseq);

    my ( $testbp, $refbp, $refpos ) = get_aligned_bases( $testseq, $refseq );

    for my $i ( 1 .. $#{$refbp} - $readlen + 1 ) {

        my @testread;
        my @refread;
        my $cigar = "";
        my $start;
        my @snps;

        my $ins       = 0;
        my $del       = 0;
        my $mat       = 0;
        my $write_ins = 0;
        my $write_del = 0;
        my $write_mat = 0;

        my $j = $i - 1;

        while ( @testread < $readlen and $j <= $#{$refbp} ) {
            if ( $testbp->[$j] eq '-' ) {
                $write_ins++ if $ins > 0;
                $write_mat++ if $mat > 0;
                $del++;
            }
            elsif ( $refbp->[$j] eq '-' ) {
                $write_del++ if $del > 0;
                $write_mat++ if $mat > 0;
                $ins++;
            }
            else {
                $write_ins++ if $ins > 0;
                $write_del++ if $del > 0;
                $mat++;

                push @snps, $refpos->[$j] if ( $testbp->[$j] ne $refbp->[$j] );
            }
            ( $write_ins, $ins, $cigar ) =
              test_block( $write_ins, $ins, 'I', $cigar );
            ( $write_del, $del, $cigar ) =
              test_block( $write_del, $del, 'D', $cigar );
            ( $write_mat, $mat, $cigar ) =
              test_block( $write_mat, $mat, 'M', $cigar );

            if ( $del == 0 ) {
                $start = $refpos->[$j] if not defined $start;
                push @testread, $testbp->[$j];
                push @refread,  $refbp->[$j];
            }
            $j++;
        }
        $cigar .= $ins . "I" if $ins > 0;
        $cigar .= $del . "D" if $del > 0;
        $cigar .= $mat . "M" if $mat > 0;

        my $testreadseq = join( '', @testread );

        next if $testreadseq =~ /N/ or length($testreadseq) < $readlen;

        my $snpstr = join ',', @snps;
        my $snpnum = @snps;
        print $sam "$start\t"
          . join( '', @testread )
          . "\t$snpnum\t$snpstr\t$cigar\n";
        print $fq '@'
          . "$start\_$snpnum\_$snpstr\_$cigar\n$testreadseq\n+\n"
          . 'I' x $readlen . "\n";
    }
    close $fq;
    close $sam;
}

sub test_block {
    my ( $write, $len, $let, $cigar ) = @_;
    if ($write) {
        $cigar .= "$len$let";
        $write = 0;
        $len   = 0;
    }
    return ( $write, $len, $cigar );
}

sub get_aligned_bases {
    my ( $testseq, $refseq ) = @_;

    my @testbases;
    my @refbases;
    my @refpos;

    my $refi = 0;
    for my $i ( 0 .. length($refseq) ) {
        my $ti = substr $testseq, $i, 1;
        my $ri = substr $refseq,  $i, 1;

        next if ( $ti eq '-' and $ri eq '-' );

        push @testbases, $ti;
        push @refbases,  $ri;

        if ($ri eq '-') {
            push @refpos, '-';
        }
        else {
            $refi++;
            push @refpos, $refi;
        }
    }

    return ( \@testbases, \@refbases, \@refpos );
}

sub load_alignment {
    my $alignment = shift;

    my %seqs;
    open my $fa, '<', $alignment
      or croak "Can't open alignment file! $OS_ERROR\n";
    my $seqname  = "";
    my $seqbases = "";
    while ( my $faline = <$fa> ) {
        chomp $faline;
        if ( $faline =~ /^>([^\s]+)/ ) {
            $seqs{$seqname} = $seqbases if $seqname ne "";
            $seqname        = $1;
            $seqbases       = "";
        }
        else {
            $seqbases .= $faline;
        }
    }
    $seqs{$seqname} = $seqbases;
    close $fa;

    return \%seqs;
}
