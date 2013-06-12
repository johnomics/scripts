#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

use Term::ExtendedColor qw/:all/;
use Number::Format qw/:subs :vars/;
use List::Util qw/min max/;

$OUTPUT_AUTOFLUSH = 1;

my %chrommarker = (
    "HHAAHHAAHAAAHHHAAAAAAAAHAHHHAAHHHAHHHHHHHAAHHHAHHHHAHHAAAAHAHAHAHHAAA" =>
      "7",
    "AAHHAAHHAHHHAAAHHHHHHHHAHAAAHHAAAHAAAAAAAHHAAAHAAAAHAAHHHHAHAHAHAAHHH" =>
      "7",
    "AHHAAHHAAHHHAAHHHHHAAAHHHHAAHHHAAHAAHHHAHAHHAAAHAAHHAAHHAAAAAAHHAAHAA" =>
      "19",
    "HAAHHAAHHAAAHHAAAAAHHHAAAAHHAAAHHAHHAAAHAHAAHHHAHHAAHHAAHHHHHHAAHHAHH" =>
      "19",
    "AHAAAAHAAAAHAAHAAAAAHAAHHAAHAHHAAHAAHHAAHHHAAAAHHAAHAHHHAHHAHHHAAAHHH" =>
      "18",
    "HAHHHHAHHHHAHHAHHHHHAHHAAHHAHAAHHAHHAAHHAAAHHHHAAHHAHAAAHAAHAAAHHHAAA" =>
      "18",
    "HHHAAHHAAAAHAAAHAHAHAAHAAHHHHHAHHHAAHAHHAAAHAHHAHHHAHHHAAAAAAHHAAAAHH" =>
      "12",
    "AAAHHAAHHHHAHHHAHAHAHHAHHAAAAAHAAAHHAHAAHHHAHAAHAAAHAAAHHHHHHAAHHHHAA" =>
      "12",
    "AHHHHHHAHAHAAHAAAAHHHHHAAHAHHAHAHAAHAAAHAHHHAHHAHHAHHAAHHAHHAAHHAAHHH" =>
      "1",
    "HAAAAAAHAHAHHAHHHHAAAAAHHAHAAHAHAHHAHHHAHAAAHAAHAAHAAHHAAHAAHHAAHHAAA" =>
      "1",
    "HAHHAAAAHHHHAHAHHAAAHHHHHHAAAAAAAAHHHHHAHHHAHAHAAAHAHAHHAHHHAAHHAAAHH" =>
      "13",
    "AHAAHHHHAAAAHAHAAHHHAAAAAAHHHHHHHHAAAAAHAAAHAHAHHHAHAHAAHAAAHHAAHHHAA" =>
      "13",
    "HHAHAAAAAHAAHHAHAHHAHHAHAHAHHHHHAHHHHAHHAAAAHHAHHAAAAHHHAHAHHAHHHHAHA" =>
      "6",
    "AAHAHHHHHAHHAAHAHAAHAAHAHAHAAAAAHAAAAHAAHHHHAAHAAHHHHAAAHAHAAHAAAAHAH" =>
      "6",
    "HAAAAAHAAHHHHAAAAAAHAAHAAAHHHAHHAAHAHHHHHHHHAAHAHHAHAAHAAHHHHHAHHAAHH" =>
      "15",
    "AHHHHHAHHAAAAHHHHHHAHHAHHHAAAHAAHHAHAAAAAAAAHHAHAAHAHHAHHAAAAAHAAHHAA" =>
      "15",
    "AHAHAHAAHAAHHHHHHAHHHAAHHHHAAHAAHAAHAAHAAAAAHAAHHAAAAAAAAHHHHAHHAHHHA" =>
      "3",
    "HAHAHAHHAHHAAAAAAHAAAHHAAAAHHAHHAHHAHHAHHHHHAHHAAHHHHHHHHAAAAHAAHAAAH" =>
      "3",
    "AHHAAHHAHHHAHHHHHHAAAAHAHHAAAHHHHAAHAHAAHHAAAHAAAHHHHHHHHHAHHHAAAHAHH" =>
      "16",
    "HAAHHAAHAAAHAAAAAAHHHHAHAAHHHAAAAHHAHAHHAAHHHAHHHAAAAAAAAAHAAAHHHAHAA" =>
      "16",
    "HAHHAAAHAAHHAHHAAHHHHHHAHAAHHAHHHHHAAHHHHHHHAHAAHAHHHAAHAHHHHAHAHAAHA" =>
      "11",
    "AHAAHHHAHHAAHAAHHAAAAAAHAHHAAHAAAAAHHAAAAAAAHAHHAHAAAHHAHAAAAHAHAHHAH" =>
      "11",
    "AAHHAHHAAHHHAAAHAAHAHAAAAAHHAHHAAAAHHAHAAHAHAHHHAAHAAHAAAHAHAAAAAAAAH" =>
      "8",
    "HHAAHAAHHAAAHHHAHHAHAHHHHHAAHAAHHHHAAHAHHAHAHAAAHHAHHAHHHAHAHHHHHHHHA" =>
      "8",
    "HAAAHHAHHAAHAHHAHAHAHAAHHAHHAAHAAHHAAHHHHAHHHAAAHHHHAHHAAAHHHAAAHHAHA" =>
      "17",
    "AHHHAAHAAHHAHAAHAHAHAHHAAHAAHHAHHAAHHAAAAHAAAHHHAAAAHAAHHHAAAHHHAAHAH" =>
      "17",
    "HHAHAAHAAAHAHAHAAAAHAAAAAHHAAHAHAHHAHAAHAHHHHHAAAHHHAAHHAAAAHAAAAAHAH" =>
      "5",
    "AAHAHHAHHHAHAHAHHHHAHHHHHAAHHAHAHAAHAHHAHAAAAAHHHAAAHHAAHHHHAHHHHHAHA" =>
      "5",
    "HHAAAHAAHAHHAAHAAHAAAHHAHAAAHAAHAHHHHHHAHAHAHHHHHAAAAAHHAAHAAHHAHAAHH" =>
      "9",
    "AAHHHAHHAHAAHHAHHAHHHAAHAHHHAHHAHAAAAAAHAHAHAAAAAHHHHHAAHHAHHAAHAHHAA" =>
      "9",
    "HAAHAAAAAHAAHAHHHHHAHHAAAAAAHHHHHAHAAHHAHAHAHHAHAAHHAHAAHHHAAAHHAAHAA" =>
      "4",
    "HAAHAHAAAAHHAAAAAAAHHHHHAHAHHHAAAAHAHAHAAAHAAHAHAAHHHHHHAHAAHHAAHHAHA" =>
      "Z",
    "AHHAHAHHHHAAHHHHHHHAAAAAHAHAAAHHHHAHAHAHHHAHHAHAHHAAAAAAHAHHAAHHAAHAH" =>
      "Z",
    "AHHHAHAHHAAHHHHHHAHHAAHHHHHHAHHHAHHAHAHAAHHAAAHHAHHAHHAAAAAAAAHHAAHAH" =>
      "2",
    "HAAAHAHAAHHAAAAAAHAAHHAAAAAAHAAAHAAHAHAHHAAHHHAAHAAHAAHHHHHHHHAAHHAHA" =>
      "2",
    "HAAAHAHHHAAAAHAHHAAAAHHAAHAAHAAHAAHAAAAAAAAAHHHAAHHHHAAAAAAAHAAAAHHHH" =>
      "20",
    "AHHHAHAAAHHHHAHAAHHHHAAHHAHHAHHAHHAHHHHHHHHHAAAHHAAAAHHHHHHHAHHHHAAAA" =>
      "20"
);

my %args;

$args{markerfile} = "";
$args{lengthfile} = "";

my $options_okay = GetOptions(
    'markers=s' => \$args{markerfile},
    'lengths=s' => \$args{lengthfile}
);
croak "No marker file! Please specify -m $OS_ERROR\n"
  if ( $args{markerfile} eq "" );

croak "No lengths file! Please specify -l $OS_ERROR\n"
  if ( $args{lengthfile} eq "" );

open my $lengthfile, '<', $args{lengthfile}
  or croak "Can't open $args{lengthfile}! $OS_ERROR\n";

my %scflen;
my $genomelen = 0;
while ( my $lengthline = <$lengthfile> ) {
    chomp $lengthline;
    my ( $scf, $len ) = split /\t/, $lengthline;
    $scflen{$scf} = $len;
    $genomelen += $len;
}
close $lengthfile;

my $scfnum = keys %scflen;
print "Genome has $scfnum scaffolds and is $genomelen bp long\n";

open my $markerfile, '<', $args{markerfile}
  or croak "Can't open $args{markerfile}! $OS_ERROR\n";

my @scflines;
my $prev_scf;
my %chrom;
my %scf;
my %marker;
my $scf_count    = 0;
my $marker_count = 0;
while ( my $marker = <$markerfile> ) {
    next
      if ( $marker =~ /^(-+)$/ );
    chomp $marker;
    my ( $scf, $pos, $type, $parents, $mq, $fs, $p, $rmsobs, $rmspat, $orig, $edge, $cons, $corrected ) =
      split /\t/,
      uncolor($marker);
    next if ( $type !~ "Maternal-A" );
    $marker_count++;
    $prev_scf = $scf if !defined $prev_scf;

    if ( $scf ne $prev_scf ) {
        $scf_count++;
        print STDERR "."            if ( $scf_count % 10 == 0 );
        print STDERR "$scf_count\n" if ( $scf_count % 100 == 0 );

        process_scf( $prev_scf, \@scflines );
        @scflines = ();
    }
    $prev_scf = $scf;
    push @scflines, "$pos:$corrected";
}
process_scf( $prev_scf, \@scflines );
close $markerfile;

print "\n";

my $genfound = 0;
foreach my $chr ( sort { $a <=> $b } keys %chrom ) {
    foreach my $scf ( keys %{ $chrom{$chr}{scf} } ) {
        $chrom{$chr}{scflen} += $scflen{$scf};
    }
    my $scfnum = keys %{ $chrom{$chr}{scf} };
    print "$chr\t";
    printf "%4d\t", $scfnum;
    print format_picture( $chrom{$chr}{markerlen}, '###,###,###' ) . "\t"
      . format_picture( $chrom{$chr}{scflen}, '###,###,###' ) . "\n";
    $genfound += $chrom{$chr}{markerlen};
}

my %scfmarkers;
foreach my $scf ( keys %scf ) {
    my $unassigned = 0;
    my $assigned   = 0;
    foreach my $chr ( keys %{ $scf{$scf} } ) {
        if ( $chr eq "0" ) {
            $unassigned += $scf{$scf}{$chr};
        }
        else {
            $assigned += $scf{$scf}{$chr};
        }
    }
    my $all    = $unassigned + $assigned;
    my $unique = keys %{ $scf{$scf} };
    $unique-- if ( defined $scf{$scf}{"0"} );
    $scfmarkers{$all}{$assigned}{$unassigned}{$unique}{count}++;
    $scfmarkers{$all}{$assigned}{$unassigned}{$unique}{length} += $scflen{$scf};
}

my %unassigned;
foreach my $scf ( keys %marker ) {
    my $min = min keys %{ $marker{$scf} };
    my $max = max keys %{ $marker{$scf} };

    foreach my $start ( keys %{ $marker{$scf} } ) {
        if ( $marker{$scf}{$start}{chr} eq "0" ) {
            my $length_bin =
              $marker{$scf}{$start}{length} < 10 ? $marker{$scf}{$start}{length}
              : $marker{$scf}{$start}{length} < 100
              ? int( $marker{$scf}{$start}{length} / 10 ) * 10
              : int( $marker{$scf}{$start}{length} / 100 ) * 100;
            if ( $start == $min or $start == $max ) {
                $unassigned{$length_bin}{end}++;
            }
            else {
                $unassigned{$length_bin}{middle}++;
            }
        }
    }
}

print "\n\nLen\tEnds\tMiddles\n";
my %unassigned_regions;
foreach my $length_bin ( sort { $a <=> $b } keys %unassigned ) {
    print "$length_bin";
    foreach my $type ( sort keys %{ $unassigned{$length_bin} } ) {
        $unassigned_regions{$type} += $unassigned{$length_bin}{$type};
        print "\t$unassigned{$length_bin}{$type}";
    }
    print "\n";
}

print
"Unassigned regions: $unassigned_regions{middle} middle, $unassigned_regions{end} ends\n";

print "\n\nAll\tChr\tNo chr\tUnique chr\tCount\tLength\n";
my $total_scf         = 0;
my $total_length      = 0;
my $unique_scf        = 0;
my $unique_length     = 0;
my $unassigned_scf    = 0;
my $unassigned_length = 0;
foreach my $markercount ( sort { $a <=> $b } keys %scfmarkers ) {
    foreach my $assigned (
        sort { $a <=> $b }
        keys %{ $scfmarkers{$markercount} }
      )
    {
        foreach my $unassigned (
            sort { $a <=> $b }
            keys %{ $scfmarkers{$markercount}{$assigned} }
          )
        {
            foreach my $unique (
                sort { $a <=> $b }
                keys %{ $scfmarkers{$markercount}{$assigned}{$unassigned} }
              )
            {
                print
"$markercount\t$assigned\t$unassigned\t$unique\t$scfmarkers{$markercount}{$assigned}{$unassigned}{$unique}{count}\t$scfmarkers{$markercount}{$assigned}{$unassigned}{$unique}{length}\n";
                $total_scf +=
                  $scfmarkers{$markercount}{$assigned}{$unassigned}
                  {$unique}{count};
                $total_length +=
                  $scfmarkers{$markercount}{$assigned}{$unassigned}
                  {$unique}{length};

                if ( $unique == 1 ) {
                    $unique_scf +=
                      $scfmarkers{$markercount}{$assigned}{$unassigned}
                      {$unique}{count};
                    $unique_length +=
                      $scfmarkers{$markercount}{$assigned}{$unassigned}
                      {$unique}{length};
                }

                if ( $unique == 0 ) {
                    $unassigned_scf +=
                      $scfmarkers{$markercount}{$assigned}{$unassigned}
                      {$unique}{count};
                    $unassigned_length +=
                      $scfmarkers{$markercount}{$assigned}{$unassigned}
                      {$unique}{length};
                }
            }
        }
    }
}

print "\nType\tScaffolds\tLength\tMarker length\n";
print "Total\t$total_scf\t"
  . format_picture( $total_length, "###,###,###" ) . "\t"
  . format_picture( $genfound,     "###,###,###" ) . "\n";
print "Unique\t$unique_scf\t"
  . format_picture( $unique_length, "###,###,###" ) . "\n";
print "Unassigned\t$unassigned_scf\t"
  . format_picture( $unassigned_length, "###,###,###" ) . "\n";
print "\n$marker_count markers found\n";

sub process_scf {
    my ( $scf, $scflines ) = @_;

    my $prev_pattern;
    my $start;
    my $prev_pos;

    my $snps = @{$scflines};
    foreach my $marker ( @{$scflines} ) {
        my ( $pos, $pattern ) = split /:/, $marker;
        $prev_pattern = $pattern if !defined $prev_pattern;
        $start        = $pos     if !defined $start;
        if ( $pattern ne $prev_pattern ) {
            process_marker( $scf, $start, $prev_pos, $prev_pattern );
            $start = $pos;
        }
        $prev_pattern = $pattern;
        $prev_pos     = $pos;
    }
    process_marker( $scf, $start, $prev_pos, $prev_pattern );
}

sub process_marker {
    my ( $scf, $start, $end, $pattern ) = @_;
    my $length = $end - $start + 1;
    my $chr = $chrommarker{$pattern} // "0";
    $chrom{$chr}{scf}{$scf}++;
    $chrom{$chr}{markerlen} += $length;
    $scf{$scf}{$chr}++;

    my $closest_chr = $chr;
    my $min_hamming = length $pattern;
    if ( $chr eq "0" ) {
        foreach my $chrprint ( keys %chrommarker ) {
            my $hamming = ( $pattern ^ $chrprint ) =~ tr/\001-\255//;
            if ( $hamming < $min_hamming ) {
                $min_hamming = $hamming;
                $closest_chr = $chrommarker{$chrprint};
            }
        }
    }
    $marker{$scf}{$start}{pattern} = $pattern;
    $marker{$scf}{$start}{chr}     = $chr;
    $marker{$scf}{$start}{length}  = $length;
    print
"$scf\t$start\t$end\t$pattern\t$length\t$chr\t$closest_chr\t$min_hamming\n";
}
