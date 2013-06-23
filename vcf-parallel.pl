#!/usr/bin/env perl

# vcf-parallel.pl
# John Davey johnomics@gmail.com

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

use Parallel::ForkManager;

$OUTPUT_AUTOFLUSH = 1;

my %args;

$args{vcf_filename}     = "";
$args{vcf_subs}         = "";
$args{threads}          = 1;
$args{output_prefix}    = "vcf-parallel-output";
$args{genetics}         = "";
$args{lengths_filename} = "";
$args{rms_filename}     = "";
$args{byscf}            = 0;
$args{maxscf}           = 0;
$args{uniquescf}        = "";

my $options_okay = GetOptions(
    'vcf=s'       => \$args{vcf_filename},
    'script=s'    => \$args{vcf_subs},
    'threads=i'   => \$args{threads},
    'output=s'    => \$args{output_prefix},
    'genetics=s'  => \$args{genetics},
    'lengths=s'   => \$args{lengths_filename},
    'rmsfile=s'   => \$args{rms_filename},
    'byscaffold'  => \$args{byscf},
    'maxscf=i'    => \$args{maxscf},
    'uniquescf=s' => \$args{uniquescf},
);
croak "No VCF file! Please specify -v $OS_ERROR\n"
  if ( $args{vcf_filename} eq "" );

croak
"No genetics file! Please specify -g to define parents, poor quality individuals and marker types\n"
  if ( $args{genetics} eq "" );

if ( $args{vcf_subs} eq "" ) {
    print STDERR "No script file given with -s, so will count lines in file\n";

    sub process {
        my ( $scf, $line, $samples, $data, $genetics ) = @_;
        $data->{lines}++;
    }

    sub merge {
        my ( $part, $all ) = @_;
        return if !defined $part;
        $all->{lines} += $part->{lines};
    }

    sub output {
        my ( $data, $genome, $outfix ) = @_;

        print "$data->{lines} lines in VCF file\n";
    }
}
else {
    require( $args{vcf_subs} );
}

my $genetics = load_genetics( $args{genetics} );

my $genome = load_genome( \%args );

my $samples = get_samples( $args{vcf_filename}, $genetics );

$genome->{scfp} = get_partitions( $genome, $args{threads} );
$args{threads} = keys %{ $genome->{scfp} };    # last thread may not be used;
     # each part is always larger than minimum size, so scaffolds for final
     # thread may be spread around the other threads

print_partitions( $genome->{scfp}, $genome->{scfl} );

output( parse_vcf( $genome, $samples, $genetics, \%args ),
    $samples, $genome, $args{output_prefix} );

sub load_genetics {
    my $geneticsfilename = shift;

    my %genetics;

    open my $geneticsfile, "<", $geneticsfilename
      or croak "Can't open marker type file $geneticsfilename: $OS_ERROR\n";

    my $infoline;
    while ( $infoline = <$geneticsfile> ) {
        chomp $infoline;
        last if ( $infoline =~ /Type/ );

        if ( $infoline =~ /^Ignore/ ) {
            my ( $ignore, $ind ) = split /\t/, $infoline;
            $genetics{ignore}{$ind} = 0;
        }
        elsif ( $infoline =~ /^Parents/ ) {
            my ( $header, $parents ) = split /\t/, $infoline;
            my @parents = split /,/, $parents;
            $genetics{parents} = \@parents;
        }
        elsif ( $infoline =~ /^Female/ or $infoline =~ /^Male/ ) {
            my ( $sex, $samples ) = split /\t/, $infoline;
            my @samples = split /,/, $samples;
            $genetics{sex}{$sex} = \@samples;
            map { $genetics{samplesex}{$_} = $sex } @samples;
        }
    }

    # Now $infoline contains type table header; ignore

    my %f2patterns;
    while ( my $marker_type = <$geneticsfile> ) {
        chomp $marker_type;
        my ( $parents, $males, $females, $type, $corrections ) = split /\t/,
          $marker_type;
        $genetics{types}{$parents}{$type}{'Male'}   = $males;
        $genetics{types}{$parents}{$type}{'Female'} = $females;

        $genetics{masks}{$type}{'Male'} = $genetics{masks}{$type}{'Male'}
          // generate_masks( $genetics{types}{$parents}{$type}{'Male'} );
        $genetics{masks}{$type}{'Female'} = $genetics{masks}{$type}{'Female'}
          // generate_masks( $genetics{types}{$parents}{$type}{'Female'} );

        if ( $corrections ne "" ) {
            my ( $orig, $fix ) = split '->', $corrections;
            $genetics{types}{$parents}{$type}{corrections}{$orig} =
              $fix;
        }

        $f2patterns{"$males:$females"}++;
    }
    close $geneticsfile;

    if ( keys %f2patterns ) {
        my $rmsfilename;
        if ( $args{rms_filename} ne "" ) {
            $rmsfilename = $args{rms_filename};
        }
        else {
            $rmsfilename = $geneticsfilename;
            $rmsfilename =~ s/txt$/rms_pval.txt/;
        }

        my $numf2patterns = keys %f2patterns;

        if ( !-e $rmsfilename ) {
            if ( !defined $args{rms_filename} ) {
                system(
"generate_rms_distributions.pl -g $geneticsfilename -t $numf2patterns -s 1000000"
                );
            }
            else {
                croak "RMS filename doesn't exist! $args{rms_filename}\n";
            }
        }

        open my $rmsfile, "<", $rmsfilename
          or croak "Can't open $rmsfilename! $OS_ERROR\n";
        my $header = <$rmsfile>;
        chomp $header;
        my @patterns = split /\t/, $header;
        shift @patterns;    # Pvalue
        while ( my $pval_line = <$rmsfile> ) {
            chomp $pval_line;
            my @vals = split /\t/, $pval_line;
            my $pval = shift @vals;
            for my $i ( 0 .. $#patterns ) {
                $genetics{rms}{ $patterns[$i] }{$pval} = $vals[$i];
            }
        }

        close $rmsfile;
    }
    \%genetics;
}

sub generate_masks {
    my ($gtstr) = @_;
    my @gts = split /,/, $gtstr;
    map { $_ = substr $_, -1; } @gts;

    my %mask;
    if ( @gts == 2 ) {
        $mask{ $gts[0] }{ $gts[1] }++;
    }
    else {
        for my $i ( 0 .. $#gts ) {
            for my $j ( 0 .. $#gts ) {
                next if $i eq $j;
                $mask{ $gts[$i] }{ $gts[$j] }++;
            }
        }
    }
    \%mask;
}

sub parse_vcf {
    my ( $genome, $samples, $genetics, $argref ) = @_;
    my %data;

    my $part_pm = new Parallel::ForkManager( $argref->{threads} );
    $part_pm->set_max_procs( $argref->{threads} );

    $part_pm->run_on_finish(
        sub {
            my $pid       = shift;
            my $exit      = shift;
            my $part      = shift;
            my $childdata = pop;
            merge( $childdata, \%data );
            print STDERR "$part:Merged\n";
        }
    );

    foreach my $part ( 1 .. $argref->{threads} ) {
        $part_pm->start($part) and next;

        my ( $partdata, $scfi ) =
          run_part( $part, $genome, $samples, $genetics, $argref );

        print STDERR "$part:Done, processed $scfi scaffolds of " .
          keys( %{ $genome->{scfp}{$part} } ) . "\n";
        $part_pm->finish( 0, $partdata );
    }
    $part_pm->wait_all_children;

    \%data;
}

sub run_part {
    my ( $part, $genome, $samples, $genetics, $argref ) = @_;

    open my $vcf_file, '<', $argref->{vcf_filename}
      or croak "Can't open $argref->{vcf_filename} $OS_ERROR!\n";

    # Skip header
    while (<$vcf_file>) {
        last if (/^#CHROM/);
    }

    my $curscf    = "";
    my $scfi      = 1;
    my $foundpart = 0;

    my %data;
    my @scf_vcf;
    my $scf;
    while ( my $vcf_line = <$vcf_file> ) {
        last if ( $argref->{maxscf} && $scfi >= $argref->{maxscf} );
        $scf = $vcf_line =~ /^(.+?)\t/ ? $1 : "";

        $curscf = $scf if $curscf eq "";    # Fill curscf at the start

        if ( defined $genome->{scfp}{$part}{$scf} ) {
            $foundpart = 1;

            if ( $curscf ne $scf ) {
                if ( $argref->{byscf} && @scf_vcf ) {
                    process( $curscf, \@scf_vcf, $samples, \%data, $genetics,
                        $genome->{scfl} )
                      if $argref->{uniquescf} eq ""
                      or $curscf eq $argref->{uniquescf};
                    @scf_vcf = ();
                    print_part_progress( $scfi, $part, $genome->{scfp} )
                      if ( $scfi % 10 == 0 );
                    $scfi++;
                }
                $curscf = $scf;
            }
            $argref->{byscf}
              ? push @scf_vcf,
              $vcf_line
              : process( $curscf, $vcf_line, $samples, \%data, $genetics,
                $genome->{scfl} );
        }
        else {
            if ($foundpart) {
                process( $curscf, \@scf_vcf, $samples, \%data, $genetics,
                    $genome->{scfl} )
                  if $argref->{byscf}
                  && ( $argref->{uniquescf} eq ""
                    or $curscf eq $argref->{uniquescf} );
                @scf_vcf = ();
                last;
            }
        }
    }
    close $vcf_file;
    process( $curscf, \@scf_vcf, $samples, \%data, $genetics, $genome->{scfl} )
      if $argref->{byscf}
      && @scf_vcf
      && ( $argref->{uniquescf} eq ""
        or $curscf eq $argref->{uniquescf} );

    return ( \%data, $scfi );
}

sub print_part_progress {
    my ( $scfi, $part, $scfpref ) = @_;

    printf STDERR "%3d:%4d scaffolds processed, %4d remaining\n",
      $part, $scfi, keys( %{ $scfpref->{$part} } ) - $scfi;

}

sub print_partitions {
    my ( $scfpref, $scflref ) = @_;

    my $genscf   = 0;
    my $genpartl = 0;
    print STDERR "Part\tScaffolds\tLength\n";
    foreach my $part ( sort { $a <=> $b } keys %{$scfpref} ) {
        my $partl  = 0;
        my $numscf = keys %{ $scfpref->{$part} };
        foreach my $scf ( keys %{ $scfpref->{$part} } ) {
            $partl += $scflref->{$scf};
        }
        print STDERR "$part\t$numscf\t$partl\n";
        $genscf   += $numscf;
        $genpartl += $partl;
    }

    print STDERR "Genome\t$genscf\t$genpartl\n";
    return;
}

sub get_partitions {
    my ( $genome, $threads ) = @_;

    my %scfp;
    my $part      = 1;
    my $threshold = $genome->{length} / $threads;
    my $part_size = 0;
    for my $scf ( @{ $genome->{scf} } ) {
        $scfp{$part}{$scf}++;
        $part_size += $genome->{scfl}{$scf};
        if ( $part_size > $threshold ) {
            $part_size = 0;
            $part++;
        }
    }

    \%scfp;
}

sub get_samples {
    my $vcf_filename = shift;
    my $genetics     = shift;

    my %parents;
    map { $parents{$_} = 0 } @{ $genetics->{parents} };

    open my $vcf_file, '<', $vcf_filename
      or croak "Can't open VCF file $vcf_filename! $OS_ERROR\n";

    my %samples;

    my $vcf_line;
    while ( $vcf_line = <$vcf_file> ) {
        if ( $vcf_line =~ /^#CHROM/ ) {
            last;
        }
    }
    close $vcf_file;

    chomp $vcf_line;
    my @sample_names = split /\t/, $vcf_line;
    map {
        if ( defined $parents{ $sample_names[$_] } ) {
            $samples{parents}{lookup}{ $sample_names[$_] } = $_;
            push @{ $samples{parents}{order} }, $sample_names[$_];
        }
        else {
            if ( !defined $genetics->{ignore}{ $sample_names[$_] } ) {
                $samples{offspring}{lookup}{ $sample_names[$_] } = $_;
                push @{ $samples{offspring}{order} }, $sample_names[$_];
            }
        }
    } 9 .. $#sample_names;

    \%samples;
}

sub load_genome {
    my $argref = shift;

    my @scflines =
         parse_vcf_header( $argref->{vcf_filename} )
      or load_lengths( $argref->{lengths_filename} )
      or croak
      "No scaffold lengths in VCF header or scaffold lengths file (-l)\n";
    return load_scflen(@scflines);
}

sub load_scflen {
    my @scflines = @_;

    my %scfl;
    my @scf;
    my $genl = 0;
    foreach my $scfline (@scflines) {
        my ( $scf, $len ) = split /\t/, $scfline;
        $scfl{$scf} = $len;
        push @scf, $scf;
        $genl += $len;
    }
    { scfl => \%scfl, scf => \@scf, length => $genl };
}

sub parse_vcf_header {
    my $vcf_filename = shift;
    open my $vcf_file, '<', $vcf_filename
      or croak "Can't open VCF file $vcf_filename! $OS_ERROR\n";

    my @scflines;
    my $vcf_line = <$vcf_file>;
    while ( $vcf_line !~ /^#CHROM/ ) {
        if ( $vcf_line =~ /^##contig=<ID=(.+),length=(\d+)>$/ ) {
            push @scflines, "$1\t$2";
        }
        $vcf_line = <$vcf_file>;
    }
    close $vcf_file;

    @scflines;
}

sub load_lengths {
    my $lengths_filename = shift;

    my @scflines;

    return @scflines if ( $lengths_filename eq "" );

    open my $lengths_file, "<", $lengths_filename
      or croak
      "Can't open scaffold lengths file $lengths_filename! $OS_ERROR\n";

    while ( my $scf_line = <$lengths_file> ) {
        if ( $scf_line =~ /^(.+)\t(.+)$/ ) {
            push @scflines, "$1\t$2";
        }
    }
    close $lengths_file;

    @scflines;
}

