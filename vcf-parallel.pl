#!/usr/bin/env perl

# vcf-parallel.pl
#
# Author: John Davey john.davey@ed.ac.uk
# Begun 01/02/13 based on convert_vcf_to_joinmap.pl

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

use Parallel::ForkManager;

# Autoflush output so reporting on progress works
$| = 1;

my $vcf_filename     = "";
my $vcf_subs         = "";
my $threads          = 1;
my $max_scf          = 0;
my $output_prefix    = "vcf-parallel-output";
my $parent_string    = "";
my $lengths_filename = "";
my $by_scaffold      = 0;

my $options_okay = GetOptions(
    'vcf=s'      => \$vcf_filename,
    'script=s'   => \$vcf_subs,
    'threads=i'  => \$threads,
    'max_scf=i'  => \$max_scf,
    'output=s'   => \$output_prefix,
    'parents=s'  => \$parent_string,
    'lengths=s'  => \$lengths_filename,
    'byscaffold' => \$by_scaffold,
);
croak "No VCF file! Please specify -v $OS_ERROR\n" if ( $vcf_filename eq "" );
croak
"No processing script! Please specify -s with process, merge and output subroutines $OS_ERROR\n"
  if ( $vcf_subs eq "" );

require $vcf_subs;

my %scflen;
my @scf;
my $genome_length = 0;

open my $vcf_file, '<', $vcf_filename
  or croak "Can't open $vcf_filename $OS_ERROR!\n";

# Parse VCF header
my $vcf_line = <$vcf_file>;
while ( $vcf_line !~ /^#CHROM/ ) {
    if ( $vcf_line =~ /^##contig=<ID=(.+),length=(\d+)>$/ ) {
        push @scf, $1;
        $scflen{$1} = $2;
        $genome_length += $2;
    }
    $vcf_line = <$vcf_file>;
}

if ( $genome_length == 0 ) {

    if ( $lengths_filename ne "" ) {
        open my $lengths_file, "<", $lengths_filename
          or croak
"No scaffold lengths in VCF header and can't open scaffold lengths file $lengths_filename! $OS_ERROR\n";
        while ( my $scflen = <$lengths_file> ) {
            chomp $scflen;
            my ( $scf, $len ) = split /\t/, $scflen;
            push @scf, $scf;
            $scflen{$scf} = $len;
            $genome_length += $len;
        }
        close $lengths_file;
        if ( $genome_length == 0 ) {
            croak
"No scaffold length information! Add scaffold lengths to VCF header or provide a scaffold lengths file with -l\n";
        }
    }
    else {
        croak
"No scaffold lengths in VCF header and no scaffold lengths file given. Please specify -l\n";
    }
}

# Last header line is the field names, including sample names
chomp $vcf_line;
my @sample_names = split /\t/, $vcf_line;
for my $i ( 0 .. 8 ) { shift @sample_names; }

close $vcf_file;

my %scf_partition;
my $partition      = 1;
my $threshold      = $genome_length / $threads;
my $partition_size = 0;
map {
    $scf_partition{$partition}{$_}++;
    $partition_size += $scflen{$_};
    if ( $partition_size > $threshold ) {
        $partition_size = 0;
        $partition++;
    }
} @scf;

my $genome_scf   = 0;
my $gen_part_len = 0;
foreach my $partition ( sort { $a <=> $b } keys %scf_partition ) {
    my $partlen = 0;
    my $numscf  = keys %{ $scf_partition{$partition} };
    foreach my $scf ( keys %{ $scf_partition{$partition} } ) {
        $partlen += $scflen{$scf};
    }
    print STDERR "$partition\t$numscf\t$partlen\n";
    $genome_scf   += $numscf;
    $gen_part_len += $partlen;
}

print STDERR "Genome\t$genome_scf\t$gen_part_len\n";

my %final_data;

my $partition_pm = new Parallel::ForkManager($threads);
$partition_pm->set_max_procs($threads);

$partition_pm->run_on_finish(
    sub {
        my ( $pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref ) =
          @_;
        merge( $data_ref, \%final_data );
    }
);

foreach my $partition ( 1 .. $threads ) {
    $partition_pm->start($partition) and next;

    open my $vcf_file, '<', $vcf_filename
      or croak "Can't open $vcf_filename $OS_ERROR!\n";

    my $cur_scf    = "";
    my $scf_count  = 0;
    my $part_start = 0;

    my %data;
    my @scf_vcf;

    while ( my $vcf_line = <$vcf_file> ) {
        next if ( $vcf_line =~ /^#/ );
        last if ( ( $max_scf > 0 ) && ( $scf_count >= $max_scf ) );
        my $scf = "";
        if ( $vcf_line =~ /^(.+?)\t/ ) {
            $scf = $1;
        }
        if ( defined $scf_partition{$partition}{$scf} ) {
            $part_start = 1;
            if ( $cur_scf ne $scf ) {
                $cur_scf = $scf;
                $scf_count++;
                if ( $scf_count % 10 == 0 ) {
                    my $scfs_remaining =
                      keys( %{ $scf_partition{$partition} } ) - $scf_count;
                    printf STDERR
                      "%3d:%4d scaffolds processed, %4d remaining\n",
                      $partition, $scf_count, $scfs_remaining;
                }
                if (($by_scaffold) && (@scf_vcf > 0)) {
                    process( \@scf_vcf, \@sample_names, \%data,
                        $parent_string );
                        @scf_vcf = ();
                }
            }
            if ( $by_scaffold ) {
                push @scf_vcf, $vcf_line;
            }
            else {
                process( $vcf_line, \@sample_names, \%data, $parent_string );
            }
        }
        else {
            last if ($part_start);
        }
    }
    if ($by_scaffold) {
        process (\@scf_vcf, \@sample_names, \%data, $parent_string);
    }
    close $vcf_file;
    print STDERR "$partition:Done, processed $scf_count scaffolds of " .
      keys( %{ $scf_partition{$partition} } ) . "\n";
    $partition_pm->finish( 0, \%data );
}

$partition_pm->wait_all_children;

output( \%final_data, \%scflen, \@sample_names,
    $genome_length, $output_prefix );
