#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

use File::Basename;

$OUTPUT_AUTOFLUSH = 1;

my %args;

$args{evalue} = 1e-50;

my $options_okay =
  GetOptions( 'assembly=s' => \$args{assembly}, 'evalue=s' => \$args{evalue} );
croak "No assembly file! Please specify -a $OS_ERROR\n"
  if ( $args{assembly} eq "" );

my $blastfile = basename( $args{assembly} );
$blastfile = "dennis_refs." . $blastfile;
$blastfile =~ s/\.fa$/\.blastn\.out/;

open my $blast, '<', $blastfile
  or croak
  "Can't open BLAST output $blastfile in current directory! $OS_ERROR\n";

my %contig;
my %hit;
while ( my $hit = <$blast> ) {
    chomp $hit;
    my @f = split /\t/, $hit;
    next if $f[10] > $args{evalue};
    $contig{ $f[1] }++;
    $hit{ $f[0] }{ $f[1] }{qstart} = $f[6];
    $hit{ $f[0] }{ $f[1] }{qend}   = $f[7];
    $hit{ $f[0] }{ $f[1] }{sstart} = $f[8];
    $hit{ $f[0] }{ $f[1] }{send}   = $f[9];
}
close $blast;

open my $assembly, '<', $args{assembly}
  or croak "Can't open assembly file $args{assembly}! $OS_ERROR\n";

my $contig = 0;
my %conseq;
while ( my $fasta = <$assembly> ) {
    if ( $fasta =~ /^>(.+)$/ ) {
        $contig = defined $contig{$1} ? $1 : "";
    }
    else {
        if ($contig) {
            chomp $fasta;
            $conseq{$contig} .= $fasta;
        }
    }
}
close $assembly;


my $scaffoldfile = $blastfile;
$scaffoldfile =~ s/\.blastn\.out$/.scaffold.fa/;

open my $scaffold, '>', $scaffoldfile or croak "Can't open scaffold file $scaffoldfile! $OS_ERROR\n";

while ( my $ref = each %hit ) {
    my $header       = "$ref";
    my $ref_assembly = "";
    my $gapstart     = 0;
    foreach my $contig (
        sort { $hit{$ref}{$a}{qstart} <=> $hit{$ref}{$b}{qstart} }
        keys %{ $hit{$ref} }
      )
    {
        my $seq = $conseq{$contig};
        $header .= " $contig";
        if ( $hit{$ref}{$contig}{sstart} > $hit{$ref}{$contig}{send} ) {
            $seq = revcomp($seq);
            $header .= "-RC";
        }
        $ref_assembly .= $seq;
        if ( $gapstart != 0 ) {
            my $gapsize = $hit{$ref}{$contig}{qend} - $gapstart + 1;
            $gapsize = 1 if ( $gapsize < 1 );
            $ref_assembly .= 'N' x $gapsize;
        }
        $gapstart = $hit{$ref}{$contig}{qend};
    }
    print $scaffold ">$header\n$ref_assembly\n";
}

close $scaffold;

sub revcomp {
    my $dna = shift;
    my $rc  = reverse $dna;
    $rc =~ tr/ACGTacgt/TGCAtgca/;
    return $rc;
}
