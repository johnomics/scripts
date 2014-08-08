#!/usr/bin/env perl
use Data::Dumper;

my %bamfiles = (
    "offspring-1-mutated.bam"   => "1.C115_7",
    "offspring-2-mutated.bam"   => "2.C115_7",
    "offspring-31-mutated.bam"  => "31.C115_6",
    "offspring-103-mutated.bam" => "103.C115_7",
    "offspring-4-mutated.bam"   => "4.C115_7",
    "offspring-33-mutated.bam"  => "33.C115_6",
    "offspring-110-mutated.bam" => "110.C115_6",
    "offspring-37-mutated.bam"  => "37.C115_6",
    "offspring-111-mutated.bam" => "111.C115_6",
    "offspring-114-mutated.bam" => "114.C115_6",
    "offspring-74-mutated.bam"  => "74.C115_7",
    "offspring-118-mutated.bam" => "118.C115_6",
    "offspring-120-mutated.bam" => "120.C115_7",
);

my %lanes = (
    "HWI-ST1233:132" => "121103",
    "ILLUMINA:366"   => "121017",
    "ILLUMINA:381"   => "121214",
);

foreach my $bamfile ( keys %bamfiles ) {
    open my $bam, "-|", "bamtools convert -in $bamfile -format fastq"
      or die "Can't open or convert bamfile $bamfile! $OS_ERROR\n";
    my $sample = $bamfiles{$bamfile};

    my $outputfile = "";
    my %outputfiles;
    my $lines = 0;
    while ( my $fastqline = <$bam> ) {
        if ( $fastqline =~ /^@(.+?):(.+?):(.+)\/([12])$/ ) {
            my $lane_id = "$1:$2";
            my $read    = $4;
            die "Lane id $lane_id not known!" if !defined $lanes{$lane_id};

            $outputfile = "$sample.$lanes{$lane_id}.R$read.fastq";
            if ( !defined $outputfiles{$outputfile} ) {
                open $outputfiles{$outputfile}, ">", $outputfile
                  or die "Can't open output FASTQ file $outputfile! $OS_ERROR\n";
            }
        }
        print {$outputfiles{$outputfile}} $fastqline;
    }
    close $bam;
    for my $outputfile ( keys %outputfiles ) {
        close $outputfiles{$outputfile};
    }
}
