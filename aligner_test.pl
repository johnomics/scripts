#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Data::Dumper;
use Parallel::ForkManager;

my $alignment    = "";
my $reference    = "";
my $output_dir   = "aligner_test_out";
my $threads      = 1;
my $commands     = "";
my $maxseqs      = 0;
my $options_okay = GetOptions(
    'alignment=s' => \$alignment,
    'reference=s' => \$reference,
    'commands=s'  => \$commands,
    'output=s'    => \$output_dir,
    'threads=i'   => \$threads,
    'maxseqs=i'   => \$maxseqs,
);

croak
  "Please specify a file containing a gapped alignment in FASTA format with -a"
  if $alignment eq "";
croak
  "Please specify the name of the reference sequence in the alignment with -r"
  if $reference eq "";
croak "Please specify a file containing aligner commands with -c"
  if $commands eq "";

my $seqs = load_alignment($alignment);
croak "Reference must be present in alignment"
  if !defined $seqs->{$reference};

my $aligners = load_commands($commands);

-d $output_dir or mkdir $output_dir;

prepare_reference( $reference, $output_dir, $seqs, $aligners );

my $pm = new Parallel::ForkManager($threads);
$pm->set_max_procs($threads);

my $seqsdone = 0;
for my $seqname ( keys %{$seqs} ) {
    next if $seqname eq $reference;

    last if $maxseqs and $seqsdone == $maxseqs;
    my ( $testbp, $refbp, $refpos, $stats ) =
      get_aligned_bases( $seqname, $reference, $seqs );

    next if !defined $testbp;

    for my $readlen (50) {
        $pm->start and next;
        generate_read_alignment( $testbp, $refbp, $refpos, $seqname, $readlen,
            $seqs, $reference, $output_dir );

        my $refpath   = "$output_dir/$reference";
        my $fastqpath = "$output_dir/$seqname.$readlen";

        for my $aligner ( sort keys %{$aligners} ) {

            run_aligner( $aligners->{$aligner}, $refpath, $fastqpath );
            my $result = compare_alignment( $fastqpath, $aligner );

            print
"$stats->{bases}\t$stats->{matches}\t$stats->{ns}\t$stats->{snps}\t$stats->{ins}\t$stats->{del}\t$readlen\t$aligner\t$result\t$seqname\n";
        }

        $pm->finish;
    }
    $seqsdone++;
}

$pm->wait_all_children;

sub compare_alignment {
    my ( $fastqpath, $aligner ) = @_;
    open my $result, "<", "$fastqpath.$aligner.sam"
      or croak "Can't open results file $fastqpath.$aligner.sam!\n";

    my $matchpos   = 0;
    my $matchcigar = 0;
    my $reads      = 0;
    while ( my $samline = <$result> ) {
        chomp $samline;
        next if $samline =~ /^@/;
        $reads++;
        my ( $readname, $flag, $ref, $refpos, $mq, $cigar, @other ) =
          split /\t/, $samline;

        my ( $true_readpos, $snps, $snppos, $true_cigar ) = split /_/,
          $readname;

        $matchpos++   if ( $true_readpos eq $refpos );
        $matchcigar++ if ( $true_cigar eq $cigar );
    }
    close $result;

    return "$reads\t$matchpos\t$matchcigar";
}

sub run_aligner {
    my ( $aligner, $refpath, $fastqpath ) = @_;

    for my $cmd ( @{ $aligner->{align}{cmd} } ) {
        $cmd =~ s/REFERENCE/$refpath/g;
        $cmd =~ s/FASTQ/$fastqpath/g;
        system($cmd);
    }
    return;
}

sub generate_read_alignment {
    my ( $testbp, $refbp, $refpos, $seqname, $readlen,
        $seqs, $reference, $output_dir )
      = @_;
    open my $sam, ">", "$output_dir/$seqname.$readlen.sam"
      or croak "Can't open SAM file for $seqname, $readlen bp reads\n";
    open my $fq, ">", "$output_dir/$seqname.$readlen.fq"
      or croak "Can't open read 1 FASTQ file for $seqname, $readlen bp reads\n";

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

                push @snps, $refpos->[$j]
                  if (  $testbp->[$j] ne $refbp->[$j]
                    and $testbp->[$j] ne 'N'
                    and $refbp->[$j] ne 'N' );
            }
            if ($write_ins) {
                $cigar .= $ins . "I";
                $write_ins = $ins = 0;
            }
            elsif ($write_del) {
                $cigar .= $del . "D";
                $write_del = $del = 0;
            }
            elsif ($write_mat) {
                $cigar .= $mat . "M";
                $write_mat = $mat = 0;
            }

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

        #        next if $testreadseq =~ /N/;
        next if length($testreadseq) < $readlen;

        my $snpstr = join ',', @snps;
        my $snpnum = @snps;
        print $sam "$start\t"
          . join( '', @testread ) . "\t"
          . join( '', @refread )
          . "\t$snpnum\t$snpstr\t$cigar\n";
        print $fq '@'
          . "$start\_$snpnum\_$snpstr\_$cigar\n$testreadseq\n+\n"
          . 'I' x $readlen . "\n";
    }
    close $fq;
    close $sam;
}

sub get_aligned_bases {
    my ( $seqname, $reference, $seqs ) = @_;

    my @testbases;
    my @refbases;
    my @refpos;
    my %stats;

    $stats{bases}   = 0;
    $stats{ins}     = 0;
    $stats{del}     = 0;
    $stats{ns}      = 0;
    $stats{matches} = 0;

    my $testseq = $seqs->{$seqname};
    my $refseq  = $seqs->{$reference};
    print STDERR
      "$seqname alignment is not the same length as reference alignment!\n"
      and return
      if length($testseq) != length($refseq);

    my $refi = 0;
    for my $i ( 0 .. length($refseq) ) {
        my $ti = substr $testseq, $i, 1;
        my $ri = substr $refseq,  $i, 1;

        next if ( $ti eq '-' and $ri eq '-' );

        $stats{bases}++;
        push @testbases, $ti;
        push @refbases,  $ri;

        if ( $ri eq '-' ) {
            push @refpos, '-';
            $stats{ins}++;
        }
        else {
            $refi++;
            push @refpos, $refi;

            if ( $ti eq '-' ) {
                $stats{del}++;
            }
            elsif ( $ti eq 'N' or $ri eq 'N' ) {
                $stats{ns}++;
            }
            elsif ( $ti ne $ri ) {
                $stats{snps}++;
            }
            else {
                $stats{matches}++;
            }
        }
    }

    ( \@testbases, \@refbases, \@refpos, \%stats );
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

    \%seqs;
}

sub load_commands {
    my $commands = shift;

    open my $cmdfile, '<', $commands
      or croak "Can't open commands file $commands! $OS_ERROR\n";
    my $aligner = "";
    my $stage   = "";
    my %aligners;
    while ( my $cmdline = <$cmdfile> ) {
        chomp $cmdline;
        if ( $cmdline =~ /^#(.+)(\s+)(.+)$/ ) {
            $aligner = $1;
            $stage   = $3;
        }
        else {
            push @{ $aligners{$aligner}{$stage}{cmd} }, $cmdline;
        }
    }
    close $cmdfile;

    \%aligners;
}

sub prepare_reference {
    my ( $reference, $output_dir, $seqs, $aligners ) = @_;

    my $refpath = "$output_dir/$reference";
    open my $refout, ">", "$refpath.fa"
      or croak "Can't open output reference FASTA file!\n";
    my $refoutseq = $seqs->{$reference};
    $refoutseq =~ s/-//g;
    print $refout ">$reference\n$refoutseq\n";
    close $refout;

    for my $aligner ( keys %{$aligners} ) {
        for my $cmd ( @{ $aligners->{$aligner}{index}{cmd} } ) {
            $cmd =~ s/REFERENCE/$refpath/g;
            print "$aligner:$cmd\n";
            system($cmd);
        }
    }
}
