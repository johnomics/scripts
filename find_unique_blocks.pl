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

my %args;

$args{markerfile}   = "";
$args{lengthfile}   = "";
$args{geneticsfile} = "";

my $options_okay = GetOptions(
    'markers=s'  => \$args{markerfile},
    'lengths=s'  => \$args{lengthfile},
    'genetics=s' => \$args{geneticsfile},
);
croak "No marker file! Please specify -m $OS_ERROR\n"
  if ( $args{markerfile} eq "" );

croak "No lengths file! Please specify -l $OS_ERROR\n"
  if ( $args{lengthfile} eq "" );

croak "No genetics file! Please specify -g $OS_ERROR\n"
  if ( $args{geneticsfile} eq "" );

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

my %genetics;

open my $geneticsfile, '<', $args{geneticsfile}
  or croak "Can't open $args{geneticsfile}! $OS_ERROR\n";
my $infoline;
while ( my $infoline = <$geneticsfile> ) {
    last if ( $infoline =~ /Type/ );
}

while ( my $marker_type = <$geneticsfile> ) {
    chomp $marker_type;
    my ( $parents, $males, $females, $type, $corrections ) = split /\t/,
      $marker_type;

    #    next if $type !~ 'Maternal';
    $genetics{$type}{$parents}{males}   = $males;
    $genetics{$type}{$parents}{females} = $females;
}
close $geneticsfile;

open my $markerfile, '<', $args{markerfile}
  or croak "Can't open $args{markerfile}! $OS_ERROR\n";

my %scflines;
my $prev_scf;
my %scf;
my %marker;
my $scf_count    = 0;
my $marker_count = 0;

while ( my $marker = <$markerfile> ) {
    next
      if ( $marker =~ /^(-+)$/ );
    chomp $marker;
    my (
        $scf,  $pos,  $type,   $parents, $mq,
        $fs,   $p,    $rmsobs, $rmspat,  $orig,
        $edge, $cons, $corrected
      )
      = split /\t/,
      uncolor($marker);
    next if $type eq 'Reject';

    #    next if $type !~ 'Maternal';
    $marker_count++;
    $prev_scf = $scf if !defined $prev_scf;

    if ( $scf ne $prev_scf ) {
        $scf_count++;
        print STDERR "."            if ( $scf_count % 10 == 0 );
        print STDERR "$scf_count\n" if ( $scf_count % 100 == 0 );

        process_scf( $prev_scf, \%scflines, \%genetics );
        %scflines = ();
    }
    $prev_scf                = $scf;
    $scflines{$pos}{type}    = $type;
    $scflines{$pos}{pattern} = $corrected;
}
process_scf( $prev_scf, \%scflines );
close $markerfile;

print "\n";

sub process_scf {
    my ( $scf, $scflines, $genetics ) = @_;

    printf "%8s\t%8s\t%8s\t%8s", 'Scaffold', 'Start', 'End', 'Length';
    map { printf "\t%-69s", $_; } sort keys %genetics;
    print "\n";

    my %curpat;
    for my $type ( keys %genetics ) {
        $curpat{$type} = " " x 69;
    }
    my $prevpos = 1;
    for my $pos ( sort { $a <=> $b } keys %{$scflines} ) {
        my $type = $scflines->{$pos}{type};
        if ( $curpat{$type} ne $scflines->{$pos}{pattern} ) {
            printf "%8s\t%8d\t%8d\t%8d", $scf, $prevpos, $pos - 1,
              $pos - $prevpos;
            for my $printtype ( sort keys %curpat ) {
                if ( $type eq $printtype ) {
                    print "\t";
                    my @cur = split //, $curpat{$printtype};
                    my @new = split //, $scflines->{$pos}{pattern};
                    for my $i ( 0 .. $#cur ) {
                        my $col = $cur[$i] eq $new[$i] ? 'black' : 'red1';
                        print fg $col, $cur[$i];
                    }
                }
                else {
                    print "\t$curpat{$printtype}";
                }
            }
            print "\n";
            $curpat{$type} = $scflines->{$pos}{pattern};
            $prevpos = $pos;
        }
    }
    printf "%8s\t%8d\t%8d\t%8d\n", $scf, $prevpos, $scflen{$scf},
      $scflen{$scf} - $prevpos + 1;
}
