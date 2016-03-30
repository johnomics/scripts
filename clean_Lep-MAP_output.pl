#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use POSIX;
use English;
use Data::Dumper;
use Getopt::Long;
use Term::ExtendedColor qw/:all/;
use List::Util qw/min max/;

my $crossprefix = "";
my $chromosomes = 0;
my $species     = '';
my $verbose     = '';

my $options_okay = GetOptions(
    'prefix=s'      => \$crossprefix,
    'chromosomes=s' => \$chromosomes,
    'species=s'     => \$species,
    'verbose'       => \$verbose,
);

croak
  "Please specify a cross prefix with -p (to match existing files PREFIX_posteriors.gz, PREFIX_id.err, PREFIX_id.txt)"
  if $crossprefix eq "";
croak "Please specify number of chromosomes with -c" if $chromosomes <= 0;
croak "Please specify a species label with -l" if $species eq '';

croak "$crossprefix\_id.err does not exist" if !-e "$crossprefix\_id.err";
croak "$crossprefix\_id.txt does not exist" if !-e "$crossprefix\_id.txt";

my $start = time;
print STDERR "Start time: " . localtime($start) . "\n";

my @snps       = load_snps($crossprefix);
my @markers    = load_markers($crossprefix);
my @markersnps = load_markersnps($crossprefix);
my %errors     = load_errors($crossprefix);
my %reverses   = load_reverses($crossprefix);

my %paternals;
my %maternals;

for my $chr ( 1 .. $chromosomes ) {
    my %summary = get_markers( \%maternals, $crossprefix, $chr, \@snps, \@markers, \@markersnps, \%errors );

    my %patterns = get_patterns( \%summary, \%paternals, $chr );

    %patterns = reorder_map( \%patterns, defined $reverses{$chr}, $verbose );

    for my $cM ( sort { $a <=> $b } keys %patterns ) {
        $paternals{$chr}{$cM} = $patterns{$cM}{pattern};
    }
}

open my $markeroutput, ">", "$crossprefix.markers.tsv"
  or croak "Can't open $crossprefix.markers.tsv for writing\n";

write_maternals( \%maternals, $markeroutput );

write_paternals( \%paternals, $markeroutput, 1 + keys %maternals );

close $markeroutput;

my $end = time;
print STDERR "End time: " . localtime($end) . "\n";

my $runtime = $end - $start;
my $hour    = int( $runtime / 3600 );
my $min     = int( ( $runtime - $hour * 3600 ) / 60 );
my $sec     = $runtime - $hour * 3600 - $min * 60;
printf STDERR "Run time: %02d:%02d:%02d\n", $hour, $min, $sec;

sub load_snps {
    my $crossprefix = shift;
    my @snps;
    print STDERR "Loading SNPs file $crossprefix\_snps.txt\n";
    open my $snpfile, "<", "$crossprefix\_snps.txt" or croak "Can't open $crossprefix\_snps.txt\n";
    while ( my $line = <$snpfile> ) {
        next if $line =~ /^#/;
        chomp $line;
        $line =~ s/\*//g;
        push @snps, $line;
    }
    close $snpfile;

    @snps;
}

sub load_markers {
    my $crossprefix = shift;
    my @markers;
    open my $markerfile, "<", "$crossprefix\_id.err" or croak "Can't open $crossprefix\_id.err\n";
    print STDERR "Loading marker file $crossprefix\_id.err\n";
    while ( my $line = <$markerfile> ) {
        next if $line !~ /^\d/;
        chomp $line;
        my ( $id, $type, $count, $paternal, $maternal ) = split "\t", $line;
        $markers[$id] = {
            type     => $type,
            count    => $count,
            paternal => $paternal,
            maternal => $maternal
        };
    }
    close $markerfile;
    @markers;
}

sub load_markersnps {
    my $crossprefix = shift;
    my @markersnps;
    print STDERR "Loading marker SNPs file $crossprefix\_id.txt\n";
    open my $markersnpfile, "<", "$crossprefix\_id.txt" or croak "Can't open $crossprefix\_id.txt\n";
    while ( my $line = <$markersnpfile> ) {
        next if $line =~ /^#/;
        chomp $line;
        push @markersnps, $line;
    }
    close $markersnpfile;

    @markersnps;
}

sub load_errors {
    my $crossprefix = shift;
    my %errors;
    if ( -e "$crossprefix.errors.tsv" ) {
        print STDERR "Loading errors file $crossprefix.errors.tsv\n";
        open my $errorfile, '<', "$crossprefix.errors.tsv" or croak "Can't open $crossprefix.errors.tsv\n";
        while ( my $line = <$errorfile> ) {
            chomp $line;
            my ( $chromosome, $cM, $destination ) = split /\t/, $line;
            $errors{$chromosome}{$cM} = $destination;
        }
    }
    else {
        print STDERR "No errors file found\n";
    }
    %errors;
}

sub load_reverses {
    my $crossprefix = shift;
    my %reverses;
    if ( -e "$crossprefix.reverse.txt" ) {
        print STDERR "Loading reverses from $crossprefix.reverse.txt\n";
        open my $reversefile, '<', "$crossprefix.reverse.txt" or croak "Can't open $crossprefix.reverse.txt\n";
        while ( my $line = <$reversefile> ) {
            chomp $line;
            $reverses{$line}++;
        }
    }
    else {
        print STDERR "No reverses file found\n";
    }
    %reverses;
}

sub get_markers {
    my ( $maternals, $crossprefix, $chr, $snps, $markers, $markersnps, $errors ) = @_;
    my %summary;
    open my $orderfile, "<", "$crossprefix.order$chr.txt" or croak "can't open $crossprefix.order$chr.txt\n";
    my $phases    = '';
    my $unique_id = -1;
    while ( my $line = <$orderfile> ) {
        chomp $line;
        my @f = split /\t/, $line;

        my $marker_number = $f[0];
        next if $marker_number !~ /^\d+$/;

        my $paternal_cm = $f[1];
        if ( defined $errors->{$chr}{$paternal_cm} ) {
            next if $errors->{$chr}{$paternal_cm} eq '-';
            $paternal_cm = $errors->{$chr}{$paternal_cm};
        }

        my $error = "-";
        if ( $f[3] =~ /\( (.+) \)/ ) {
            $error = $1;
        }

        my $marker = $markers->[ $markersnps->[ $marker_number - 1 ] ];
        if ( $line !~ "duplicate" ) {
            $phases = $f[-1];
            if ( $marker->{type} == 2 ) {
                $maternals->{$chr}{ $marker->{maternal} } = { count => $marker->{count} }
                  if $marker->{type} == 2;
            }
            else {
                $unique_id++;
                $summary{$paternal_cm}{ $marker->{type} }{$unique_id} = {
                    count    => $marker->{count},
                    error    => $error,
                    phase    => $phases,
                    maternal => $marker->{maternal},
                    paternal => $marker->{paternal}
                };
            }
        }

        push @{ $summary{$paternal_cm}{ $marker->{type} }{$unique_id}{pos} }, $snps->[ $marker_number - 1 ]
          if $marker->{type} != 2;
    }
    close $orderfile;
    %summary;
}

sub get_patterns {
    my ( $summary, $paternals, $chr ) = @_;
    my %patterns;
    for my $cM ( sort { $a <=> $b } keys %{$summary} ) {
        next if not defined $summary->{$cM}{1};
        my $id = (
            sort {
                     $summary->{$cM}{1}{$a}{error} <=> $summary->{$cM}{1}{$b}{error}
                  or $summary->{$cM}{1}{$b}{count} <=> $summary->{$cM}{1}{$a}{count}
            } keys %{ $summary->{$cM}{1} }
        )[0];
        my $marker = $summary->{$cM}{1}{$id};
        next if $marker->{error} != 0 or $marker->{paternal} =~ /\?/;

        $patterns{$cM}{pattern} =
          ( $marker->{phase} eq '0-' ) ? $marker->{paternal} : phase( $marker->{paternal} );
        $patterns{$cM}{count} = $marker->{count};
        $patterns{$cM}{pos}   = $marker->{pos};
    }
    %patterns;
}

sub reorder_map {
    my ( $patterns, $reverse, $verbose ) = @_;

    my @cMs = sort { $a <=> $b } keys %{$patterns};

    my ( $numdiffs, $length ) = get_diffs( \@cMs, $patterns, $verbose );
    my @newcMs = ("0.000");
    for my $numdiff ( @{$numdiffs} ) {
        my $r         = $numdiff / $length;
        my $kosambi   = log( ( 1 + 2 * $r ) / ( 1 - 2 * $r ) ) / 4;
        my $newcMdist = $kosambi * 100;
        push @newcMs, sprintf "%.3f", $newcMs[-1] + $newcMdist;
    }
    if ($reverse) {
        my $maxcm = max @newcMs;
        @newcMs = map { sprintf "%.3f", $maxcm - $_ } @newcMs;
    }
    my %newpatterns;
    for my $i ( 0 .. $#cMs ) {
        $newpatterns{ $newcMs[$i] } = $patterns->{ $cMs[$i] };
    }
    %newpatterns;
}

sub get_diffs {
    my ( $cMs, $patterns, $verbose ) = @_;
    my @diffs;
    my @numdiffs;
    my $length = 0;
    for my $i ( 1 .. $#{$cMs} - 1 ) {
        $length = length $patterns->{ $cMs->[$i] }{pattern};
        my @last = split //, $patterns->{ $cMs->[ $i - 1 ] }{pattern};
        my @this = split //, $patterns->{ $cMs->[$i] }{pattern};
        my @next = split //, $patterns->{ $cMs->[ $i + 1 ] }{pattern};
        my $last_to_this = @diffs ? $diffs[-1] : join '', map { $last[$_] eq $this[$_] ? ' ' : 'X' } 0 .. $#this;
        if ( not @diffs ) {
            push @diffs,    $last_to_this;
            push @numdiffs, $last_to_this =~ tr/X//;
            print_marker( $cMs->[ $i - 1 ], $patterns ) if $verbose;
        }

        my $this_to_next = join '', map { $this[$_] eq $next[$_] ? ' ' : 'X' } 0 .. $#this;
        push @diffs,    $this_to_next;
        push @numdiffs, $this_to_next =~ tr/X//;

        my @last_diffs = split //, $diffs[-2];
        my @this_diffs = split //, $diffs[-1];
        for my $j ( 0 .. $#last_diffs ) {
            if ( $last_diffs[$j] ne ' ' and $this_diffs[$j] ne ' ' ) {
                $last_diffs[$j] = '*';
                $this_diffs[$j] = '*';
            }
        }
        $diffs[-2] = join '', @last_diffs;
        $diffs[-1] = join '', @this_diffs;

        print "\t$diffs[-2]\t$numdiffs[-2]\n" if $verbose;
        print_marker( $cMs->[$i], $patterns ) if $verbose;
    }
    print "\t$diffs[-1]\t$numdiffs[-1]\n" if $verbose;
    print_marker( $cMs->[-1], $patterns ) if $verbose;
    \@numdiffs, $length;
}

sub print_marker {
    my ( $cM, $patterns ) = @_;
    print "$cM\t$patterns->{$cM}{pattern}\t$patterns->{$cM}{count}\t" . get_range( $patterns->{$cM}{pos} ) . "\n";
}

sub get_range {
    my ($poslist) = @_;
    my %pos;
    for my $pos ( @{$poslist} ) {
        my ( $scaffold, $position ) = split /\t/, $pos;
        $pos{$scaffold}{$position}++;
    }
    my @pos;
    for my $scaffold ( sort { substr( $a, 5 ) <=> substr( $b, 5 ) } keys %pos ) {
        my @scfpos = sort { $a <=> $b } keys %{ $pos{$scaffold} };
        push @pos, sprintf "%s:%8d-%8d", $scaffold, $scfpos[0], $scfpos[-1];
    }
    join "\t", @pos;
}

sub write_maternals {
    my ( $maternals, $markeroutput ) = @_;
    for my $chr ( sort { $a <=> $b } keys %{$maternals} ) {
        my $maternal =
          ( sort { $maternals->{$chr}{$b}{count} <=> $maternals->{$chr}{$a}{count} } keys %{ $maternals->{$chr} } )[0];
        $maternal = $maternal;
        print $markeroutput "$chr\t$chr\t-1\t2\t$maternal\n";
    }
}

sub write_paternals {
    my ( $paternals, $markeroutput, $marker_id ) = @_;

    for my $chr ( sort { $a <=> $b } keys %{$paternals} ) {
        for my $cM ( sort { $a <=> $b } keys %{ $paternals->{$chr} } ) {
            print $markeroutput "$marker_id\t$chr\t$cM\t1\t$paternals->{$chr}{$cM}\n";
            $marker_id++;
        }
    }
}

sub hamming {
    return ( $_[0] ^ $_[1] ) =~ tr/\001-\255//;
}

sub phase {
    my ($marker) = @_;
    $marker = $marker;
    $marker =~ tr/01/10/;
    return $marker;
}
