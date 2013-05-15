#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

use Term::ExtendedColor qw/:all/;

$OUTPUT_AUTOFLUSH = 1;

my %args;

$args{block_file} = "";

my $options_okay = GetOptions( 'block=s' => \$args{block_file}, );
croak "No block file! Please specify -b $OS_ERROR\n"
  if ( $args{block_file} eq "" );

open my $blockfile, '<', $args{block_file}
  or croak "Can't open $args{block_file}! $OS_ERROR\n";

my $prev_scf;
my $prev_type;
my @prev_call;
my @blocks;
my %blocklength;

my $scf_count = 0;
while ( my $marker = <$blockfile> ) {
    next if ( $marker =~ /^(-+)$/ );
    chomp $marker;
    my ( $scf, $snp, $type, $orig, $edge, $cons ) = split /\t/,
      uncolor($marker);
    $prev_scf  = $scf  if !defined $prev_scf;
    $prev_type = $type if !defined $prev_type;

    if ( $scf ne $prev_scf or $type ne $prev_type ) {
        
        if ($scf ne $prev_scf) {
            $scf_count++;
            print "." if ($scf_count % 10 == 0);
            print "$scf_count\n" if ($scf_count % 100 == 0);
        }
        map { emit_block($blocks[$_], $scf) } 0 .. $#blocks;
        $prev_scf  = $scf;
        $prev_type = $type;
        @prev_call = ();
    }
    else {
        my @calls = split //, $cons;
        foreach my $sample ( 0 .. $#calls ) {
            if ( defined $prev_call[$sample]
                and $calls[$sample] ne $prev_call[$sample] )
            {
                emit_block($blocks[$sample], $scf);
                @prev_call = ();
            }
            push @{ $blocks[$sample] }, $calls[$sample];
            $prev_call[$sample] = $calls[$sample];
        }
    }
}

close $blockfile;

my %bl;
foreach my $scf (keys %blocklength) {
    shift @{$blocklength{$scf}};
    pop @{$blocklength{$scf}};
    map {$bl{$_}++;} @{$blocklength{$scf}};
}

foreach my $len (sort {$a<=>$b} keys %bl) {
    print "$len\t$bl{$len}\n";
}

sub emit_block {
    my ($blockref, $scf) = @_;
    my $blen = scalar @{$blockref};
    push @{$blocklength{$scf}}, $blen;
    @{$blockref} = ();
    return;
}
