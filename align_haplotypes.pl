#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Parallel::ForkManager;

my $config_file  = "";
my $threads      = 1;
my $dryrun       = 0;
my $stages       = "1,2,3,4";
my $options_okay = GetOptions(
    'config=s'  => \$config_file,
    'threads=i' => \$threads,
    'dryrun'    => \$dryrun,
    'stages=s'  => \$stages,
);

my %stage;
map { $stage{$_}++ } split /,/, $stages;

croak "Please specify a config file with -c" if ( $config_file eq "" );

open my $config, "<", $config_file
  or croak "Can't open $config_file:$OS_ERROR!\n";

my %samples;

while ( my $sample_line = <$config> ) {
    chomp $sample_line;
    my (
        $sample_id, $sample,    $date,   $library, $pu,
        $reference, $bamheader, $fastq1, $fastq2
    ) = split /\s+/, $sample_line;

    $samples{$sample_id}{$sample}{$date}{library}   = $library;
    $samples{$sample_id}{$sample}{$date}{pu}        = $pu;
    $samples{$sample_id}{$sample}{$date}{reference} = $reference;
    $samples{$sample_id}{$sample}{$date}{bamheader} = $bamheader;
    $samples{$sample_id}{$sample}{$date}{fastq1}    = $fastq1;
    $samples{$sample_id}{$sample}{$date}{fastq2}    = $fastq2;
}

close $config;
if ( defined $stage{1} ) {
    for my $sample_id ( sort keys %samples ) {
        for my $sample ( sort keys %{ $samples{$sample_id} } ) {
            for my $date ( sort { $a <=> $b }
                keys %{ $samples{$sample_id}{$sample} } )
            {
                my $sref = $samples{$sample_id}{$sample}{$date};

                my $refstub = $sref->{reference};
                if ( $refstub =~ /(.+)\.fa/ ) { $refstub = $1; }
                my $align_script_name = "align_files/$sample_id\.hap.align.sh";
                open my $align_script, ">", $align_script_name
                  or croak
                  "Can't write align script for $sample_id: $OS_ERROR!\n";

                print $align_script "#!/bin/bash\n";
                print $align_script
"smalt map -f bam -n $threads -o align_files/$sample_id\.hap.align.bam $refstub $sref->{fastq1} $sref->{fastq2} 2> align_files/$sample_id\.hap.align.log\n";
                close $align_script;
                chmod 0755, $align_script_name;

                system("nohup nice ./$align_script_name") if !$dryrun;
            }
        }
    }
}

my $pm = new Parallel::ForkManager($threads);
$pm->set_max_procs($threads);

if ( defined $stage{2} ) {
    for my $sample_id ( sort keys %samples ) {
        for my $sample ( sort keys %{ $samples{$sample_id} } ) {
            for my $date ( sort { $a <=> $b }
                keys %{ $samples{$sample_id}{$sample} } )
            {
                $pm->start and next;

                my $sref = $samples{$sample_id}{$sample}{$date};

                my $refstub = $sref->{reference};
                if ( $refstub =~ /(.+)\.fa/ ) { $refstub = $1; }
                my $markdup_script_name =
                  "rg.sorted.markdupes_files/$sample_id\.hap.markdup.sh";
                open my $markdup_script, ">", $markdup_script_name
                  or croak
                  "Can't write mark dupes script for $sample_id: $OS_ERROR!\n";

                print $markdup_script "#!/bin/bash\n";

                print $markdup_script
"java -Xmx20g -jar /whale-data/jd626/bin/picard/current/AddOrReplaceReadGroups.jar INPUT=align_files/$sample_id\.hap.align.bam  OUTPUT=rg.sorted_files/$sample_id\.hap.rg.sorted.bam  SORT_ORDER=coordinate ID=$sample_id PL=ILLUMINA SM=$sample DT=$date LB=$sref->{library} PU=$sref->{pu} VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true &> rg.sorted_files/$sample_id\.hap.rg.sorted.log\n";
                print $markdup_script
"java -Xmx20g -jar /whale-data/jd626/bin/picard/current/MarkDuplicates.jar INPUT=rg.sorted_files/$sample_id\.hap.rg.sorted.bam OUTPUT=rg.sorted.markdupes_files/$sample_id\.hap.rg.sorted.markdupes.bam METRICS_FILE=rg.sorted.markdupes_files/$sample_id\.hap.rg.sorted.markdupes.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &> rg.sorted.markdupes_files/$sample_id\.hap.rg.sorted.markdupes.log\n";
                print $markdup_script
"samtools index rg.sorted.markdupes_files/$sample_id\.hap.rg.sorted.markdupes.bam\n";
                close $markdup_script;

                chmod 0755, $markdup_script_name;
                system("nohup nice ./$markdup_script_name") if !$dryrun;

                $pm->finish;
            }
        }
    }
    $pm->wait_all_children;
}

if ( defined $stage{3} ) {
    for my $sample_id ( sort keys %samples ) {
        for my $sample ( sort keys %{ $samples{$sample_id} } ) {
            for my $date ( sort { $a <=> $b }
                keys %{ $samples{$sample_id}{$sample} } )
            {
                my $sref = $samples{$sample_id}{$sample}{$date};

                my $refstub = $sref->{reference};
                if ( $refstub =~ /(.+)\.fa/ ) { $refstub = $1; }
                my $rtc_script_name =
                  "rg.sorted.markdupes.realign_files/$sample_id\.hap.rtc.sh";
                open my $rtc_script, ">", $rtc_script_name
                  or croak
                  "Can't write RTC script for $sample_id: $OS_ERROR!\n";

                print $rtc_script "#!/bin/bash\n";
                print $rtc_script
"java -Xmx20g -jar /biosoft/src/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt $threads -R $sref->{reference} -I rg.sorted.markdupes_files/$sample_id\.hap.rg.sorted.markdupes.bam -o rg.sorted.markdupes.realign_files/$sample_id\.hap.intervals &> rg.sorted.markdupes.realign_files/$sample_id\.hap.RTC.log\n";
                close $rtc_script;

                chmod 0755, $rtc_script_name;
                system("nohup nice ./$rtc_script_name") if !$dryrun;
            }
        }
    }
}

if ( defined $stage{4} ) {

    $pm = new Parallel::ForkManager($threads);
    $pm->set_max_procs($threads);

    for my $sample_id ( sort keys %samples ) {
        for my $sample ( sort keys %{ $samples{$sample_id} } ) {
            for my $date ( sort { $a <=> $b }
                keys %{ $samples{$sample_id}{$sample} } )
            {
                $pm->start and next;

                my $sref = $samples{$sample_id}{$sample}{$date};

                my $refstub = $sref->{reference};
                if ( $refstub =~ /(.+)\.fa/ ) { $refstub = $1; }
                my $realign_script_name =
"rg.sorted.markdupes.realign_files/$sample_id\.hap.realign.sh";
                open my $realign_script, ">", $realign_script_name
                  or croak
                  "Can't write realign script for $sample_id: $OS_ERROR!\n";

                print $realign_script "#!/bin/bash\n";
                print $realign_script
"java -Xmx4g -jar /biosoft/src/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -T IndelRealigner -R $sref->{reference} -targetIntervals rg.sorted.markdupes.realign_files/$sample_id\.hap.intervals -I rg.sorted.markdupes.realign_files/$sample_id\.hap.rg.sorted.markdupes.bam -o rg.sorted.markdupes.realign_files/$sample_id\.hap.rg.sorted.markdupes.realign.bam &> rg.sorted.markdupes.realign_files/$sample_id\.hap.realign.log\n";
                print $realign_script
"samtools index rg.sorted.markdupes.realign_files/$sample_id\.hap.rg.sorted.markdupes.realign.bam\n";

                close $realign_script;

                chmod 0755, $realign_script_name;
                system("nohup nice ./$realign_script_name") if !$dryrun;

                $pm->finish;
            }
        }
    }
    $pm->wait_all_children;

}
