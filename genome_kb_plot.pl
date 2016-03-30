#!/usr/bin/env perl

# genome_kb_plot.pl
# Make Kumar-Blaxter plots for multiple FASTA files

# -o Output PDF name
# -w Width of final PDF figure (default 10)
# -h Height of final PDF figure (default 10)
# List of files to process

# John Davey
# johnomics@gmail.com
# Updated Friday 20 March 2015

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use File::Basename 'fileparse';
use Data::Dumper;
use POSIX qw /ceil/;

my $output = "output";
my $width  = 10;
my $height = 10;
my $ysize  = 0;
my $xsize  = 0;
my $config = "";
my $facet = '';
my $eps = '';

my $options_okay = GetOptions(
    'output=s' => \$output,
    'width=f'  => \$width,
    'height=f' => \$height,
    'ysize=i'  => \$ysize,
    'xsize=i'  => \$xsize,
    'config=s' => \$config,
    'facet'    => \$facet,
    'eps'      => \$eps,
);

sub get_genome_files {
    my ($config) = @_;
    open my $configfh, '<', $config or croak "Can't open config file!\n";
    my %genomefiles;
    my @genomenames;
    my @genomecolours;
    while (my $configln = <$configfh>) {
        chomp $configln;
        my ($name, $file, $colour) = split "\t", $configln;
        push @genomenames, $name;
        push @genomecolours, $colour if $colour;
        $genomefiles{$file}{name} = $name;
        $genomefiles{$file}{colour} = $colour if $colour;
    }
    close $configfh;
    return \@genomenames, \@genomecolours, \%genomefiles;
}

my %genomefiles;
my @genomenames;
my @genomecolours;
if ($config) {
    my ($genomenamesref, $genomecoloursref, $genomefileref) = get_genome_files($config);
    %genomefiles = %{$genomefileref};
    @genomenames = @{$genomenamesref};
    @genomecolours = @{$genomecoloursref};
}
else {
    for my $file (@ARGV) {
        $genomefiles{$file}{colour} = 0;
        my ($filename, $dirs, $suffix) = fileparse ($file, qr/\.[^.]*/ );
        $genomefiles{$file}{name} = $filename;
        @genomenames = sort keys %genomefiles;
    }
}

my %gls;
my %gsizes;
for my $genomefile (sort keys %genomefiles) {

    print "$genomefile\n";
    open my $genome, "<", $genomefile or die "Can't open $genomefile\n";
    my $length     = 0;
    my $genomesize = 0;
    my $scaffolds  = 0;
    while ( my $gl = <$genome> ) {
        if ( $gl =~ /^>/ ) {
            if ( $length > 0 ) {
                push @{ $gls{$genomefile} }, $length;
                $genomesize += $length;
                $length = 0;
                $scaffolds++;
            }
        }
        else {
            chomp $gl;
            $length += length $gl;
        }
    }
    push @{ $gls{$genomefile} }, $length if $length > 0;
    $genomesize += $length;

    $ysize = $genomesize if $genomesize > $ysize;
    $xsize = $scaffolds  if $scaffolds > $xsize;
    close $genome;
    
    $gsizes{$genomefile} = $genomesize;

    @{ $gls{$genomefile} } = sort { $b <=> $a } @{ $gls{$genomefile} };
}

open my $tsv, '>', "$output.tsv" or croak "Can't open output TSV file $output.tsv! $OS_ERROR\n";

print $tsv "Genome\tScaffold\tGenomeLength\n";
for my $gf (sort keys %genomefiles) {
    my $cumsum = 0;
    my @milestones = (50, 90, 95);
    for my $i ( 0 .. $#{ $gls{$gf} } ) {
        $cumsum += $gls{$gf}[$i];
        
        if (@milestones and $cumsum > $gsizes{$gf} * ($milestones[0]/100)) {
            print "$gf\tN$milestones[0]=". ($i+1) . "\tL$milestones[0]=$gls{$gf}[$i]\n";
            shift @milestones
        }
        print $tsv $genomefiles{$gf}{name} . "\t$i\t$cumsum\n";
    }
}

close $tsv;

open my $rscript, '>', "$output.R" or croak "Can't open output R script $output.R! $OS_ERROR\n";
print $rscript "library(ggplot2)\n";
print $rscript "library(scales)\n";
print $rscript "genomecols<-c(" . join(",", @genomecolours) . ")\n" if @genomecolours;
print $rscript "kb.df<-read.delim(\"$output.tsv\")\n";
print $rscript "kb.df\$Genome<-factor(kb.df\$Genome, levels=c(\"". join("\",\"", @genomenames). "\"))\n";
if ($eps) {
    print $rscript "postscript(\"$output.eps\", paper=\"special\", width=$width, heigh=$height, horizontal=FALSE)\n";
}
else {
    print $rscript "pdf(\"$output.pdf\", width=$width, height=$height)\n";
}

my $ymultiple = 1 . 0 x (length($ysize) - 1);
my $ymax = ceil($ysize/$ymultiple) * $ymultiple;
if ($facet) {
    print $rscript "ggplot(kb.df, aes(Scaffold, GenomeLength)) + geom_line() + facet_wrap(~Genome)";
}
else {
    print $rscript "ggplot(kb.df, aes(Scaffold, GenomeLength, colour=Genome)) + geom_line()";
}
print $rscript " + scale_x_continuous(limits=c(0, $xsize)) + scale_y_continuous(limits=c(0, $ymax),breaks=seq(0,$ymax,$ymultiple),labels=seq(0,$ymax,$ymultiple)/1000000)";
print $rscript "+theme_bw(base_size = 10)+xlab(\"Number of scaffolds, ranked in order of size\")+ylab(\"Cumulative genome length (Mb)\")";
print $rscript "+scale_colour_manual(values=genomecols)" if @genomecolours;
print $rscript "\n";
print $rscript "dev.off()\n";

close $rscript;

system("Rscript $output.R");

unlink "$output.tsv";
unlink "$output.R";
