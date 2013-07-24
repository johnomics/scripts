#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Parallel::ForkManager;

my $progeny_config_file = "";
my $threads             = 1;
my $dryrun              = 0;
my $options_okay        = GetOptions(
    'config=s'  => \$progeny_config_file,
    'threads=i' => \$threads,
    'dryrun'    => \$dryrun,
);

croak "Please specify a config file with -c" if ( $progeny_config_file eq "" );

open my $progeny_config, "<", $progeny_config_file
  or croak "Can't open $progeny_config_file:$OS_ERROR!\n";

my %progeny;

while ( my $progeny_line = <$progeny_config> ) {
    chomp $progeny_line;
    my (
        $sample_id, $sample,    $date,   $library, $pu,
        $reference, $bamheader, $fastq1, $fastq2
    ) = split /\s+/, $progeny_line;

    $progeny{$sample_id}{$sample}{$date}{library}   = $library;
    $progeny{$sample_id}{$sample}{$date}{pu}        = $pu;
    $progeny{$sample_id}{$sample}{$date}{reference} = $reference;
    $progeny{$sample_id}{$sample}{$date}{bamheader} = $bamheader;
    $progeny{$sample_id}{$sample}{$date}{fastq1}    = $fastq1;
    $progeny{$sample_id}{$sample}{$date}{fastq2}    = $fastq2;
}

close $progeny_config;

for my $sample_id ( sort keys %progeny ) {
    for my $sample ( sort keys %{ $progeny{$sample_id} } ) {
        for
          my $date ( sort { $a <=> $b } keys %{ $progeny{$sample_id}{$sample} } )
        {
            my $sref = $progeny{$sample_id}{$sample}{$date};

            my $refstub = $sref->{reference};
            if ( $refstub =~ /(.+)\.fa/ ) { $refstub = $1; }
            my $progeny_align_script_name = "align_files/$sample_id\.hap.align.sh";
            open my $progeny_align_script, ">", $progeny_align_script_name
              or croak "Can't write align script for $sample_id: $OS_ERROR!\n";

            print $progeny_align_script "#!/bin/bash\n";
            print $progeny_align_script
"cat $sref->{bamheader} <(nice -n 10 smalt map -f sam -n $threads $refstub <(gunzip -c $sref->{fastq1}) <(gunzip -c $sref->{fastq2}) 2> align_files/$sample_id\.hap.align.log) | samtools view -Sb - > align_files/$sample_id\.hap.align.bam 2>> align_files/$sample_id\.hap.align.log\n";
            close $progeny_align_script;
            chmod 0755, $progeny_align_script_name;

            system("nohup ./$progeny_align_script_name") if !$dryrun;
        }
    }
}

my $pm = new Parallel::ForkManager($threads);
$pm->set_max_procs($threads);

for my $sample_id ( sort keys %progeny ) {
    for my $sample ( sort keys %{ $progeny{$sample_id} } ) {
        for
          my $date ( sort { $a <=> $b } keys %{ $progeny{$sample_id}{$sample} } )
        {
            $pm->start and next;

            my $sref = $progeny{$sample_id}{$sample}{$date};

            my $refstub = $sref->{reference};
            if ( $refstub =~ /(.+)\.fa/ ) { $refstub = $1; }
            my $progeny_markdup_script_name = "rg.sorted.markdupes_files/$sample_id\.hap.markdup.sh";
            open my $progeny_markdup_script, ">", $progeny_markdup_script_name
              or croak "Can't write mark dupes script for $sample_id: $OS_ERROR!\n";

            print $progeny_markdup_script "#!/bin/bash\n";

            print $progeny_markdup_script
"nice -n 10 java -Xmx20g -jar /whale-data/jd626/bin/picard/current/AddOrReplaceReadGroups.jar INPUT=align_files/$sample_id\.hap.align.bam  OUTPUT=rg.sorted_files/$sample_id\.hap.rg.sorted.bam  SORT_ORDER=coordinate ID=$sample_id PL=ILLUMINA SM=$sample DT=$date LB=$sref->{library} PU=$sref->{pu} VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true &> rg.sorted_files/$sample_id\.hap.rg.sorted.log\n";
            print $progeny_markdup_script
"java -Xmx20g -jar /whale-data/jd626/bin/picard/current/MarkDuplicates.jar INPUT=rg.sorted_files/$sample_id\.hap.rg.sorted.bam OUTPUT=rg.sorted.markdupes_files/$sample_id\.hap.rg.sorted.markdupes.bam METRICS_FILE=rg.sorted.markdupes_files/$sample_id\.hap.rg.sorted.markdupes.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &> rg.sorted.markdupes_files/$sample_id\.hap.rg.sorted.markdupes.log\n";
            print $progeny_markdup_script
              "samtools index rg.sorted.markdupes_files/$sample_id\.hap.rg.sorted.markdupes.bam\n";
            close $progeny_markdup_script;
            
            chmod 0755, $progeny_markdup_script_name;
            system("nohup ./$progeny_markdup_script_name") if !$dryrun;

            $pm->finish;
        }
    }
}
$pm->wait_all_children;

for my $sample_id ( sort keys %progeny ) {
    for my $sample ( sort keys %{ $progeny{$sample_id} } ) {
        for
          my $date ( sort { $a <=> $b } keys %{ $progeny{$sample_id}{$sample} } )
        {
            my $sref = $progeny{$sample_id}{$sample}{$date};

            my $refstub = $sref->{reference};
            if ( $refstub =~ /(.+)\.fa/ ) { $refstub = $1; }
            my $progeny_rtc_script_name = "rg.sorted.markdupes.realign_files/$sample_id\.hap.rtc.sh";
            open my $progeny_rtc_script, ">", $progeny_rtc_script_name
              or croak "Can't write RTC script for $sample_id: $OS_ERROR!\n";

            print $progeny_rtc_script "#!/bin/bash\n";
            print $progeny_rtc_script
"nice -n 10 java -Xmx20g -jar /biosoft/src/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt $threads -R $sref->{reference} -I rg.sorted.markdupes_files/$sample_id\.hap.rg.sorted.markdupes.bam -o rg.sorted.markdupes.realign_files/$sample_id\.hap.intervals &> rg.sorted.markdupes.realign_files/$sample_id\.hap.RTC.log\n";
            close $progeny_rtc_script;

            chmod 0755, $progeny_rtc_script_name;
            system("nohup ./$progeny_rtc_script_name") if !$dryrun;
        }
    }
}


$pm = new Parallel::ForkManager($threads);
$pm->set_max_procs($threads);

for my $sample_id ( sort keys %progeny ) {
    for my $sample ( sort keys %{ $progeny{$sample_id} } ) {
        for
          my $date ( sort { $a <=> $b } keys %{ $progeny{$sample_id}{$sample} } )
        {
            $pm->start and next;

            my $sref = $progeny{$sample_id}{$sample}{$date};

            my $refstub = $sref->{reference};
            if ( $refstub =~ /(.+)\.fa/ ) { $refstub = $1; }
            my $progeny_realign_script_name = "rg.sorted.markdupes.realign_files/$sample_id\.hap.realign.sh";
            open my $progeny_realign_script, ">", $progeny_realign_script_name
              or croak "Can't write realign script for $sample_id: $OS_ERROR!\n";

            print $progeny_realign_script "#!/bin/bash\n";
            print $progeny_realign_script
            "nice -n 10 java -Xmx4g -jar /biosoft/src/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -T IndelRealigner -R $sref->{reference} -targetIntervals rg.sorted.markdupes.realign_files/$sample_id\.hap.intervals -I rg.sorted.markdupes.realign_files/$sample_id\.hap.rg.sorted.markdupes.bam -o rg.sorted.markdupes.realign_files/$sample_id\.hap.rg.sorted.markdupes.realign.bam &> rg.sorted.markdupes.realign_files/$sample_id\.hap.realign.log\n";
            print $progeny_realign_script
              "samtools index rg.sorted.markdupes.realign_files/$sample_id\.hap.rg.sorted.markdupes.realign.bam\n";

            close $progeny_realign_script;
            
            chmod 0755, $progeny_realign_script_name;
            system("nohup ./$progeny_realign_script_name") if !$dryrun;

            $pm->finish;
        }
    }
}
$pm->wait_all_children;