#!/usr/bin/env perl

# convert_vcf_to_joinmap.pl
#
# Input: one VCF file containing a whole mapping cross
# Output: sex-specific markers in JoinMap format
# Author: John Davey john.davey@ed.ac.uk
# Begun 07/07/11

#############################################################################
###                                                                       ###
###                                 CODE                                  ###
###                                                                       ###
#############################################################################

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

# Autoflush output so reporting on progress works
$| = 1;

my $help               = 0;
my $vcf_filename       = "";
my $qthreshold         = 10;
my $mother_name        = "";
my $father_name        = "";
my $prefix             = "";
my $offspring_filename = "";
my $uncalled_genotypes = 0;
my $reject_singletons  = 0;

my $options_okay = GetOptions(
    'help|h'            => \$help,
    'vcf=s'             => \$vcf_filename,
    'quality=i'         => \$qthreshold,
    'mother_name=s'     => \$mother_name,
    'father_name=s'     => \$father_name,
    'prefix=s'          => \$prefix,
    'offspring=s'       => \$offspring_filename,
    'uncalled=i'        => \$uncalled_genotypes,
    'reject_singletons' => \$reject_singletons
) or pod2usage( -exitval => 3, -verbose => 0 );

pod2usage( -exitval => 5, -verbose => 1 ) if $help;

croak "No VCF file! Please specify -v $OS_ERROR\n" if ( $vcf_filename eq "" );
croak "No offspring file! Please specify -o $OS_ERROR\n"
  if ( $offspring_filename eq "" );
croak "No mother name! Please specify -m $OS_ERROR\n" if ( $mother_name eq "" );
croak "No father name! Please specify -f $OS_ERROR\n" if ( $father_name eq "" );
croak "No output prefix! Please specify -p $OS_ERROR\n" if ( $prefix eq "" );

open my $offspring_file, '<', $offspring_filename
  or croak "Can't open $offspring_filename $OS_ERROR!\n";

my @offspring;
my @females;
my @males;
my $found_mother_name = 0;
my $found_father_name = 0;
while ( my $offspring = <$offspring_file> ) {
    chomp $offspring;
    next if ( $offspring eq "" );
    my ( $sample_name, $sex ) = split /\s+/, $offspring;
    next if ( ( $sample_name eq "" ) || ( $sex eq "" ) );
    if ( $sample_name eq $mother_name ) { $found_mother_name++; next; }
    if ( $sample_name eq $father_name ) { $found_father_name++; next; }
    push @offspring, $sample_name;
    if ( $sex =~ /^[Ff]/ ) {
        push @females, $sample_name;
    }
    if ( $sex =~ /^[Mm]/ ) {
        push @males, $sample_name;
    }
}
close $offspring_file;

croak "Could not find offspring in offspring file!\n" if ( @offspring == 0 );
croak "Could not find mother $mother_name in offspring file!\n"
  if ( !$found_mother_name );
croak "Could not find father $father_name in offspring file!\n"
  if ( !$found_father_name );

open my $vcf_file, '<', $vcf_filename
  or croak "Can't open $vcf_filename $OS_ERROR!\n";

my @sample_names;

my %marker_count;

my $base_count;
my %joinmap;
my $bases_w_complete_calls_above_qthres = 0;
while ( my $vcf_line = <$vcf_file> ) {
    chomp $vcf_line;
    if ( $vcf_line =~ /^#CHROM/ ) {
        @sample_names = split /\t/, $vcf_line;
        for my $i ( 0 .. 8 ) { shift @sample_names; }
        if (@sample_names != @offspring+2) {
            print "Found " . @sample_names . " offspring in VCF file but " . (@offspring+2) . " in sample file, exiting\n";
            exit;
        }
    }
    next if ( $vcf_line =~ /^#/ );

    $base_count++;
    if ( $base_count % 10000 == 0 )   { print STDERR "." }
    if ( $base_count % 100000 == 0 )  { printf STDERR "%8d", $base_count; }
    if ( $base_count % 1000000 == 0 ) { print STDERR "\n"; }

    my $missing = 0;
    $missing++ while ( $vcf_line =~ /\.\/\./g );
    next if ( $missing > $uncalled_genotypes );

    my %base;
    my @fields = split /\t/, $vcf_line;
    next if ( $fields[4] eq "." );    # If no alternate call
    my $skip_line = 0;

    my @format_fields = split /:/, $fields[8];
    my $gt_field_num  = -1;
    my $gq_field_num  = -1;
    my $dp_field_num  = -1;
    for my $i ( 0 .. ( @format_fields - 1 ) ) {
        if ( $format_fields[$i] eq "GT" ) { $gt_field_num = $i; }
        if ( $format_fields[$i] eq "GQ" ) { $gq_field_num = $i; }
        if ( $format_fields[$i] eq "DP" ) { $dp_field_num = $i; }
    }

    for my $i ( 9 .. ( @sample_names + 8 ) ) {

        if ( $fields[$i] eq "./." ) {
            $base{ $sample_names[ $i - 9 ] }{gt} = "-/-";
            $base{ $sample_names[ $i - 9 ] }{dp} = 0;
            $base{ $sample_names[ $i - 9 ] }{gq} = 0;
            next;
        }
        my @sample_fields = split /:/, $fields[$i];
        if (($gq_field_num > -1) and 
            ( not defined $sample_fields[$gq_field_num] or 
              $sample_fields[$gq_field_num] eq '.' or 
              $sample_fields[$gq_field_num] < $qthreshold )
            ) {
            if (($sample_names[$i-9] ne $mother_name) && ($sample_names[$i-9] ne $father_name)) {$missing++;}
            if ( $missing > $uncalled_genotypes ) { $skip_line++; last; }

            $base{ $sample_names[ $i - 9 ] }{gt} = "-/-";
        }
        else {
            $base{ $sample_names[ $i - 9 ] }{gt} =
              $sample_fields[$gt_field_num];
        }

        $base{ $sample_names[ $i - 9 ] }{dp} = $sample_fields[$dp_field_num] // 0;
        $base{ $sample_names[ $i - 9 ] }{dp} = 0 if $base{ $sample_names[ $i - 9 ] }{dp} eq '.';

        if ($gq_field_num > -1) {
            $base{ $sample_names[ $i - 9 ] }{gq} = $sample_fields[$gq_field_num] // 0;
            $base{ $sample_names[ $i - 9 ] }{gq} = 0 if $base{ $sample_names[ $i - 9 ] }{gq} eq '.';
        }
        
    }

    next if ($skip_line);
    $bases_w_complete_calls_above_qthres++;
    my $f1pattern = "$base{$father_name}{gt} $base{$mother_name}{gt}";
    my $marker_type;
    my %f2_genotypes;
    $f2_genotypes{"-/-"} = "--";
    if ( $f1pattern eq "0/0 0/1" ) {
        $marker_type         = "<lmxll>";
        $f2_genotypes{"0/0"} = "ll";
        $f2_genotypes{"0/1"} = "lm";
    }
    elsif ( $f1pattern eq "1/1 0/1" ) {
        $marker_type = "<lmxll>";
        $f2_genotypes{"0/1"} = "lm";
        $f2_genotypes{"1/1"} = "ll";
    }
    elsif ( $f1pattern eq "0/1 0/0" ) {
        $marker_type         = "<nnxnp>";
        $f2_genotypes{"0/0"} = "nn";
        $f2_genotypes{"0/1"} = "np";
    }
    elsif ( $f1pattern eq "0/1 1/1" ) {
        $marker_type = "<nnxnp>";
        $f2_genotypes{"0/1"} = "np";
        $f2_genotypes{"1/1"} = "nn";
    }
    elsif ( $f1pattern eq "0/1 0/1" ) {
        $marker_type         = "<hkxhk>";
        $f2_genotypes{"0/0"} = "hh";
        $f2_genotypes{"0/1"} = "hk";
        $f2_genotypes{"1/1"} = "kk";
    }
    my $scaffold       = $fields[0];
    my $joinmap_marker = "";
    my $average_depth  = 0;
    my $average_qual   = 0;

    foreach my $offspring (@offspring) {
        if (   ( defined $base{$offspring}{gt} )
            && ( defined $f2_genotypes{ $base{$offspring}{gt} } ) )
        {
            $joinmap_marker .= "$f2_genotypes{$base{$offspring}{gt}} ";
            $average_depth += $base{$offspring}{dp} // 0;
            $average_qual  += $base{$offspring}{gq} // 0;
        }
        else {
            $skip_line++;
            last;
        }
    }

    next if ($skip_line);

    $average_depth /= @offspring;
    $average_qual  /= @offspring;
  
    my $male_depth = 0;
    my $male_qual  = 0;

    if (@males) {
        foreach my $male (@males) {
            $male_depth += $base{$male}{dp} // 0;
            $male_qual  += $base{$male}{gq} // 0;
        }
        $male_depth    /= @males;
        $male_qual     /= @males;
    }
    
    my $female_depth = 0;
    my $female_qual  = 0;
    
    if (@females) {
        foreach my $female (@females) {
            $female_depth += $base{$female}{dp} // 0;
            $female_qual  += $base{$female}{gq} // 0;
        }
        $female_depth  /= @females;
        $female_qual   /= @females;
    }
    
    my $pos_dp = 0;
    my $pos_mq = 0;
    if ( $fields[7] =~ /DP=(\d+?);/ ) {
        $pos_dp = $1;
    }
    if ( $fields[7] =~ /MQ=(.+?);/ ) {
        $pos_mq = $1;
    }
    chop $joinmap_marker;
    $joinmap{$marker_type}{$joinmap_marker}{"$scaffold:$fields[1]"}{dp} =
      int $average_depth;
    $joinmap{$marker_type}{$joinmap_marker}{"$scaffold:$fields[1]"}{gq} =
      int $average_qual;
    $joinmap{$marker_type}{$joinmap_marker}{"$scaffold:$fields[1]"}{pos_dp} =
      $pos_dp;
    $joinmap{$marker_type}{$joinmap_marker}{"$scaffold:$fields[1]"}{pos_mq} =
      $pos_mq;
    $joinmap{$marker_type}{$joinmap_marker}{"$scaffold:$fields[1]"}{male_dp} =
      int $male_depth;
    $joinmap{$marker_type}{$joinmap_marker}{"$scaffold:$fields[1]"}{male_gq} =
      int $male_qual;
    $joinmap{$marker_type}{$joinmap_marker}{"$scaffold:$fields[1]"}{female_dp} =
      int $female_depth;
    $joinmap{$marker_type}{$joinmap_marker}{"$scaffold:$fields[1]"}{female_gq} =
      int $female_qual;
}

close $vcf_file;

open my $joinmap_file, '>', "$prefix.joinmap"
  or croak "Can't open Joinmap output file $prefix.joinmap: $OS_ERROR\n";
open my $pattern_file, '>', "$prefix.patterns"
  or croak "Can't open patterns output file $prefix.patterns: $OS_ERROR\n";

print $pattern_file
"Marker\tScaffold:Position\tType\tSegregation Pattern\tAverage depth\tAverage genotype quality\tTotal depth\tMapping quality\tMale depth\tMale genotype quality\tFemale depth\tFemale genotype quality\n";
my $marker_num = 0;
foreach my $type ( sort keys %joinmap ) {
    my $scf_marker = 0;
    foreach my $marker ( sort keys %{ $joinmap{$type} } ) {

        # Output markers appearing on more than one scaffold
        # or appearing more than once on one scaffold
        my $check_marker_num = 0;
        my $num_bases        = scalar keys %{ $joinmap{$type}{$marker} };
        next if $reject_singletons and $num_bases <= 1;
        $marker_num++;
        my $output_type   = $type;
        print $joinmap_file
          "$marker_num:$num_bases\t$output_type\t$marker\n";
        foreach my $base ( sort keys %{ $joinmap{$type}{$marker} } ) {
            print $pattern_file
"$marker_num\t$base\t$output_type\t$marker\t$joinmap{$type}{$marker}{$base}{dp}\t$joinmap{$type}{$marker}{$base}{gq}\t$joinmap{$type}{$marker}{$base}{pos_dp}\t$joinmap{$type}{$marker}{$base}{pos_mq}\t$joinmap{$type}{$marker}{$base}{male_dp}\t$joinmap{$type}{$marker}{$base}{male_gq}\t$joinmap{$type}{$marker}{$base}{female_dp}\t$joinmap{$type}{$marker}{$base}{female_gq}\n";
        }
    }
}

close $pattern_file;
close $joinmap_file;

=head1 NAME

Convert VCF to JoinMap

=head1 SYNOPSIS

=over 8

=item convert_vcf_to_joinmap.pl --vcf <VCF_FILE> --offspring <OFFSPRING_FILE> --mother MOTHER_NAME --father FATHER_NAME --prefix OUTPUT_PREFIX [-q QUALITY_THRESHOLD (default 10)]

=item convert_vcf_to_joinmap.pl --help

=back

=cut
