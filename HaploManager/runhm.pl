#!/usr/bin/env perl

# runhm.pl
# Run HaploMerger

# -f FASTA file of genome, zipped or unzipped
# -d Directory containing HaploMerger scripts and config files:
#    - hm.batchA-G
#    - scoreMatrix.q
#    - all_lastz.ctl
# -o Name of output directory
# -p final scaffold prefix

# Process
# 1. Create output directory
# 2. Copy config directory files to output directory
# 3. Copy FASTA file to directory, rename as genome.fa
# 4. Run hm.batchA
# 5. Run hm.batchB
# 6. Run hm.batchC
# 7. Create merged assembly of optimized scaffolds and unpaired scaffolds:
#    - If -p is provided, add prefix to optimized scaffolds
#    - Rename unpaired scaffolds to original names

# John Davey
# johnomics@gmail.com
# Begun Friday 20 March 2015

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Data::Dumper;
use File::Basename 'fileparse';
use File::Copy 'copy';
use File::Copy::Recursive 'dircopy';
use File::Path 'rmtree';
use IO::Uncompress::Gunzip qw($GunzipError);
use Cwd;

# Autoflush output so reporting on progress works
$| = 1;

my $input           = "";
my $configdir       = "";
my $outputdir       = "";
my $prefix          = "";
my $scaffold_prefix = "";
my $misassemblies   = "";
my $badmerges       = "";
my $annotation      = "";
my $g;

my $options_okay = GetOptions(
    'input=s'           => \$input,
    'configdir=s'       => \$configdir,
    'outputdir=s'       => \$outputdir,
    'prefix=s'          => \$prefix,
    'scaffold_prefix=s' => \$scaffold_prefix,
    'misassemblies=s'   => \$misassemblies,
    'badmerges=s'       => \$badmerges,
    'annotation=s'      => \$annotation,
    'g'                 => \$g,
);

croak "No output directory! Please specify -o $OS_ERROR\n" if $outputdir eq "";

if ( !-d $outputdir ) {
    croak "No config directory! Please specify -c $OS_ERROR\n" if $configdir eq "";
    croak "Config directory $configdir does not exist!\n" if !-d $configdir;
    dircopy $configdir, $outputdir;
}

if ( !-e "$outputdir/genome.fa" ) {
    croak "No FASTA file! Please specify -f $OS_ERROR\n" if $input eq "";
    croak "FASTA file $input does not exist!\n" if !-e $input;
    copy $input, "$outputdir/genome.fa";
}

my %misassemblies;
if ( $misassemblies and -e $misassemblies ) {
    open my $misfh, '<', $misassemblies or croak "Can't open misassemblies file $misassemblies $OS_ERROR\n";
    while ( my $misline = <$misfh> ) {
        chomp $misline;
        my ( $m, $scaffold, $start, $end, $type, $details ) = split /\t/, $misline;
        next if $scaffold =~ /x/;    # Ignore non-merged scaffolds for now
        if ( $scaffold =~ /.+Sc(\d+)/ ) {
            $scaffold = int($1);
        }
        push @{ $misassemblies{$scaffold} },
          { scaffold => $scaffold, start => $start, end => $end, type => $type, details => $details };
    }
    close $misfh;
}

my %badmerges;
if ( $badmerges and -e $badmerges ) {
    open my $bmfh, '<', $badmerges or croak "Can't open bad merges file $badmerges $OS_ERROR\n";
    while ( my $badmergeline = <$bmfh> ) {
        chomp $badmergeline;
        my ( $a, $b ) = split "\t", $badmergeline;
        $badmerges{$a}{$b} = 1;
        $badmerges{$b}{$a} = 1;
    }
    close $bmfh;
}

my %genes;
if ($annotation and -e $annotation) {
    open my $gff, '<', $annotation or croak "Can't open annotation file $annotation $OS_ERROR\n";
    while (my $featureline = <$gff>) {
        chomp $featureline;
        my ($scaffold, $source, $featuretype, $start, $end, $score, $strand, $phase, $attributes) = split "\t", $featureline;
        if ($featuretype eq 'gene') {
            $genes{$scaffold}{$start}{end} = $end;
            $genes{$scaffold}{$start}{attributes} = $attributes;            
        }
    }
    close $gff;
}

my ( $filename, $dirs, $suffix ) = fileparse( $input, qr/\.[^.]*/ );
chdir $outputdir;

my $log;
if ( !-e "runhm.log" ) {
    open $log, '>', "runhm.log";
}
else {
    open $log, '>>', "runhm.log";
}

print $log "HaploMerger run log\n";
print $log "Genome: $input\n";
print $log "Config directory: $configdir\n";
print $log "Output directory: $outputdir\n";
print $log "Prefix: '$prefix'\n";
print $log "Scaffold prefix: '$scaffold_prefix'\n";
print $log "Misassemblies: $misassemblies\n";
print $log "Run G: ";
print $log ( $g ? "Yes" : "No" ), "\n";

my $start_time = localtime;
my $start      = time;
print $log "Start time: $start_time\n";

my $optimized_file = "optiNewScaffolds.fa.gz";
my $unpaired_file  = "unpaired.fa.gz";

( $optimized_file, $unpaired_file, $log, my $next_start ) = run_abc( $log, $start );

( $optimized_file, $unpaired_file, $log, $next_start ) = run_e( $log, $next_start, \%badmerges )
  if $badmerges ne "";

( $optimized_file, $unpaired_file, $log, $next_start ) =
  run_f( $log, $next_start, $scaffold_prefix, \%misassemblies, \%genes )
  if $scaffold_prefix ne "" or %genes;

output_final_genome( $outputdir, $optimized_file, $unpaired_file, $prefix );

if ($g) {
    ( $optimized_file, $unpaired_file, $log, $next_start ) = run_g( $log, $next_start );
    output_final_genome( $outputdir, $optimized_file, $unpaired_file, $prefix, "refined" );
}

open $log, '>>', "runhm.log";

if ( -d 'genome.genomex.result/raw.axt' ) {
    printf $log "Removing raw.axt folder\n";
    rmtree ['genome.genomex.result/raw.axt'];
}

my $end_time = localtime;
print $log "End time: $end_time\n";

output_duration( "Total", $log, $start );
print $log "Done\n";
close $log;

sub run_abc {
    my ( $log, $start ) = @_;
    my $c_start = $start;
    if ( !-e "genome.genomex.result/mafFiltered.net.maf.tar.gz" ) {
        system "./hm.batchA.initiation_and_all_lastz genome > runhm.out 2>&1";
        my $b_start = output_duration( "A", $log, $start );
        system "./hm.batchB.chainNet_and_netToMaf genome >> runhm.out 2>&1";
        $c_start = output_duration( "B", $log, $b_start );
    }

    my $next_start = $c_start;
    if ( !-e "genome.genomex.result/optiNewScaffolds.fa.gz" ) {
        system "./hm.batchC.haplomerger genome >> runhm.out 2>&1";
        $next_start = output_duration( "C", $log, $c_start );
    }

    ( "optiNewScaffolds.fa.gz", "unpaired.fa.gz", $log, $next_start );
}

sub run_e {
    my ($log, $start, $badmerges ) = @_;
    my $next_start = $start;
    
    if (! -e "genome.genomex.result/hm.nodes_edited" ) {
        edit_nodes($badmerges);
        chdir "genome.genomex.result";
        system "ln -s hm.sc_portions hm.sc_portions_edited" if ! -l "hm.sc_portions_edited";
        system "ln -s hm.assembly_errors hm.assembly_errors_edited" if ! -l "hm.assembly_errors_edited";
        chdir "..";
    
        if ( -e "genome.genomex.result/hm.nodes_edited" ) {
            system "./hm.batchE.refine_haplomerger_updating_nodes_and_portions genome >> runhm.out 2>&1";
            $next_start = output_duration( "E", $log, $next_start );
        }
    }
    ( "optiNewScaffolds.fa.gz", "unpaired.fa.gz", $log, $next_start)
}

sub run_f {
    my ( $log, $start, $scaffold_prefix, $misassemblies, $genes ) = @_;
    my $next_start = $start;

    if (! -e "genome.genomex.result/hm.new_scaffolds_edited") {
        edit_new_scaffolds( $scaffold_prefix, $misassemblies, $genes );        
        if ( -e "genome.genomex.result/hm.new_scaffolds_edited" ) {
            system "./hm.batchF.refine_haplomerger_connections_and_Ngap_fillings genome >> runhm.out 2>&1";
            $next_start = output_duration( "F", $log, $next_start );
        }
    }

    ( "optiNewScaffolds.fa.gz", "unpaired.fa.gz", $log, $next_start );
}

sub edit_nodes {
    my ($badmerges) = @_;
    
    my $nodes_file = "genome.genomex.result/hm.nodes";
    open my $nodes, '<', $nodes_file or croak "Can't open nodes file!\n";
    open my $edited, '>', "genome.genomex.result/hm.nodes_edited" or croak "Can't open edited nodes file!\n";

    while (my $portion = <$nodes>) {
        if ($portion =~ /^#/ or $portion =~ /^$/) {
            print $edited $portion;
            next;
        }
        chomp $portion;
        my @f = split "\t", $portion;
        if ((defined $badmerges->{$f[0]} and defined $badmerges->{$f[0]}{$f[1]}) or (defined $badmerges->{$f[1]} and defined $badmerges->{$f[1]}{$f[0]})) {
            $f[-1] = 1;
        }
        my $edit = join "\t", @f;
        $edit .= "\n";
        
        print $edited $edit;
    }
    close $nodes;
    close $edited;
}

sub edit_new_scaffolds {
    my ( $scaffold_prefix, $misassemblies, $genes ) = @_;

    my $new_scaffolds_file = "genome.genomex.result/hm.new_scaffolds";
    if ( -e "genome.genomex.result/hm.new_scaffolds_updated" ) {    # If batchE has been run
        $new_scaffolds_file = "genome.genomex.result/hm.new_scaffolds_updated";
    }

    open my $new_scaffolds, '<', $new_scaffolds_file
      or croak "Can't open new scaffolds file!\n";
    open my $edited, '>', "genome.genomex.result/hm.new_scaffolds_edited"
      or croak "Can't open edited new scaffolds file!\n";

    my $scfnum    = 0;
    my $curid     = -1;
    my $scflength = 0;

    my @portions;
    while (my $portion = <$new_scaffolds>) {        
        my @f = ($portion =~ /^#/ or $portion =~ /^$/) ? ($portion) : split "\t", $portion;
        push @portions, \@f;
    }

    for my $p (0..@portions-1) {
        my @f = @{$portions[$p]};
        if (@f == 1) {
            print $edited $f[0];
            next;
        }
        
        $curid = $f[0] if $curid == -1;
        if ( $f[0] != $curid ) {
            $scfnum++;
            $curid     = $f[0];
            $scflength = 0;
        }

        my $scaffold1 = $f[5];
        my $scaffold2 = $f[12];

        my $active_portion = 0;
        if ( $scaffold_prefix ne "" and $scaffold1 ne '0' and $scaffold2 ne '0' ) {
            if ( $scaffold1 =~ /^$scaffold_prefix/ and $scaffold2 !~ /^$scaffold_prefix/ ) {
                $active_portion = 2;
            }
            elsif ( $scaffold1 !~ /^$scaffold_prefix/ and $scaffold2 =~ /^$scaffold_prefix/ ) {
                $active_portion = 1;
            }
        }

        $f[-2] = $active_portion if $active_portion;

        if ($genes) {
            find_broken_genes($p, \@portions, $genes);
        }

        $active_portion = $f[-2] ? $f[-2] : $f[4];


        my $partlength = $active_portion == 1 ? $f[11] : $f[18];

        my $partend = $scflength + $partlength;
        if ( defined $misassemblies{$scfnum} ) {
            foreach my $mis ( @{ $misassemblies{$scfnum} } ) {
                if ( $mis->{start} < $scflength + $partlength and $mis->{end} > $scflength + 1 ) {
                    $f[-1] = "-1\n";
                }
            }
        }
        $scflength += $partlength;

        my $edit = join "\t", @f;

        print $edited $edit;
    }
    close $new_scaffolds;
    close $edited;
}

sub find_broken_genes {
    my ($p, $portions, $genes) = @_;
    my $this = get_portion_region($p, $portions);
    return if !$this;
    my $scfgenes = $genes->{$this->{scaffold}};
    for $start (keys %{$scfgenes}) {
        my $end = $scfgenes->{$start}{end};
        my $neighbour = $p;
        if ($start < $this->{start} and $this->{start} <= $end) {
            $neighbour = $this->{dir} == 1 ? $p-1 : $p+1;
        }
        if ($start <= $this->{end} and $this->{end} < $end) {
            $neighbour = $this->{dir} == 1 ? $p+1 : $p-1;
        }
        next if $neighbour == $p;
        my $np = get_portion_region($neighbour, $portions);
        if ($this->{scaffold} ne $np->{scaffold}) {
            $portions->[$neighbour][-2] = $np->{active} == 1 ? 2 : 1;
        }
    }
}

sub get_portion_region {
    my ($p, $portions) = @_;
    return if !defined $portions->[$p];
    my @f = @{$portions->[$p]};
    return if (@f == 1);
    my $active_portion = $f[-2] ? $f[-2] : $f[4];
    my %region;
    $region{active} = $active_portion;
    if ($active_portion == 1) {
        $region{scaffold} = $f[5];
        $region{start}    = $f[8] + 1;
        $region{end}      = $f[9];
        $region{dir}      = $f[10];
    }
    elsif ($active_portion == 2) {
        $region{scaffold} = $f[12];
        $region{start}    = $f[15] + 1;
        $region{end}      = $f[16];
        $region{dir}      = $f[17];
    }
    return \%region;
}
sub run_g {
    my ( $log, $start ) = @_;
    my $next_start = $start;

    if ( !-e "genome.genomex.result/unpaired_refined.fa.gz" ) {
        system "./hm.batchG.refine_unpaired_sequences genome >> runhm.out 2>&1";
        $next_start = output_duration( "G", $log, $start );
    }
    ( "optiNewScaffolds.fa.gz", "unpaired_refined.fa.gz", $log, $next_start );
}

sub output_final_genome {
    my ( $outputdir, $optimized_file, $unpaired_file, $prefix, $suffix ) = @_;

    chdir "genome.genomex.result";

    my $outname = $suffix ? "$outputdir\_$suffix" : "$outputdir";

    if ( !-e "$outname.fa" ) {
        open my $finalgenome, '>', "$outname.fa";

        output_optimized_genome( $optimized_file, $finalgenome, $prefix, $outname );

        output_unpaired_genome( $unpaired_file, $finalgenome, $prefix, $outname );

        close $finalgenome;

        system "summarizeAssembly.py $outname.fa > $outname.summary";

        system "agp_from_fasta.py -f $outname.fa";
    }
    chdir "..";
}

sub output_optimized_genome {
    my ( $optimized_file, $finalgenome, $prefix, $outname ) = @_;

    if ( !-e $optimized_file ) {
        print $log "$optimized_file does not exist, abandon writing to final genome\n";
        return;
    }

    print $log "Writing $optimized_file to $outname.fa\n";
    my $opti = IO::Uncompress::Gunzip->new($optimized_file)
      or die "IO::Uncompress::Gunzip failed to open $optimized_file: $GunzipError\n";

    while ( my $fastaline = <$opti> ) {
        if ( $fastaline =~ /^>(.+) (.+)$/ ) {
            print $finalgenome ">$prefix$1 $2\n";
        }
        else {
            print $finalgenome $fastaline;
        }
    }

    close $opti;

}

sub output_unpaired_genome {
    my ( $unpaired_file, $finalgenome, $prefix, $outname ) = @_;

    if ( !-e $unpaired_file ) {
        print $log "$unpaired_file does not exist, abandon writing to final genome\n";
        return;
    }

    print $log "Writing $unpaired_file to $outname.fa\n";
    my $unpaired = IO::Uncompress::Gunzip->new($unpaired_file)
      or die "IO::Uncompress::Gunzip failed to open $unpaired_file: $GunzipError\n";

    while ( my $fastaline = <$unpaired> ) {
        if ( $fastaline =~ />(.+) old_(.+?);(.+)/ ) {
            print $finalgenome ">$prefix$1 old_$2;$3\n";
        }
        else {
            print $finalgenome $fastaline;
        }
    }

    close $unpaired;

}

sub output_duration {
    my ( $stage, $file, $start ) = @_;
    my $end      = time;
    my $duration = $end - $start;

    printf $file "$stage run time: %02d:%02d:%02d\n", int( $duration / 3600 ), int( ( $duration % 3600 ) / 60 ),
      int( $duration % 60 );

    $end;
}
