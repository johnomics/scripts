use List::Util qw/sum min max/;
use POSIX qw/ceil/;
use Term::ExtendedColor qw/:all/;
use Memoize;
use Getopt::Long qw/GetOptionsFromString/;

use strict;
use warnings;

use constant MASKLEN      => 4;
use constant ERRORBLOCK   => 2000;
use constant ERRORSNPS    => 5;
use constant FS_THRESHOLD => 35;
use constant MQ_THRESHOLD => 57;
use constant GAPSIZE      => 5000;

my %swapphase = (
    "Intercross" => { "A" => "B", "B" => "A" },
    "Maternal"   => { "A" => "H", "H" => "A" },
    "Paternal"   => { "A" => "H", "H" => "A" }
);

my %callorder;

my %nullcall = ( 'GT' => './.', 'GQ' => 0, 'DP' => 0 );

my @mask;

our $setup = sub {
    my ( $vcf_filename, $extraargs ) = @_;

    my $genetics_filename = "";
    my $rms_filename      = "";
    my $options_okay      = GetOptionsFromString(
        $extraargs,
        'genetics=s'    => \$genetics_filename,
        'rmsfilename=s' => \$rms_filename,
    );
    croak "Can't process user options: $OS_ERROR\n" if !$options_okay;
    
    croak "No genetics file! Please specify -g to define parents, poor quality individuals and marker types\n"
      if ( $genetics_filename eq "" );

    my $genetics = load_genetics( $genetics_filename, $rms_filename );

    my $samples = get_samples( $vcf_filename, $genetics );

    { genetics => $genetics, samples => $samples };
};

sub load_genetics {
    my ($geneticsfilename, $rmsfilename) = @_;

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
        my ( $parents, $males, $females, $type, $corrections ) = split /\t/, $marker_type;
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
        if ( $rmsfilename eq "" ) {
            $rmsfilename = $geneticsfilename;
            $rmsfilename =~ s/txt$/rms_pval.txt/;
        }

        my $numf2patterns = keys %f2patterns;

        if ( !-e $rmsfilename ) {
            system("generate_rms_distributions.pl -g $geneticsfilename -t $numf2patterns -s 1000000");
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

our $process = sub {
    my ( $scf, $scfref, $data, $userdata, $scfl ) = @_;
    $data->{$scf} = get_markers( $scfref, $userdata->{samples}, $userdata->{genetics} );
    find_edges( $data->{$scf}, $userdata->{samples}, $userdata->{genetics} );
    collapse( $data->{$scf}, $userdata->{samples} );
    output_markers_to_file( $scf, $data->{$scf} );
    output_blocks_to_file( $scf, $data->{$scf}, $userdata->{samples}, $userdata->{genetics}, $scfl );
};

sub collapse {
    my ( $scfref, $samples ) = @_;
    foreach my $type ( keys %{$scfref} ) {
        foreach my $pos ( keys %{ $scfref->{$type} } ) {
            my $gt = "";
            my @gq;

            my @edge;
            my $cons      = "";
            my $corrected = "";
            
            foreach my $sample ( @{ $samples->{offspring}{order} } ) {
                my $ref = $scfref->{$type}{$pos}{marker}{$sample};

                $gt .= $ref->{gt} // '-';
                push @gq, $ref->{gq};

                if ( $type ne "Reject" ) {
                    push @edge, defined $ref->{edge}
                      && $ref->{edge} != 0 ? abs( $ref->{edge} ) : '.';
                    $cons      .= $ref->{cons}      // '-';
                    $corrected .= $ref->{corrected} // '-';
                }
            }
            delete $scfref->{$type}{$pos}{marker};
            $scfref->{$type}{$pos}{marker}{gt} = $gt;
            $scfref->{$type}{$pos}{marker}{gq} = \@gq;

            if ( $type ne "Reject" ) {
                $scfref->{$type}{$pos}{marker}{edge}      = \@edge;
                $scfref->{$type}{$pos}{marker}{cons}      = $cons;
                $scfref->{$type}{$pos}{marker}{corrected} = $corrected;
            }
        }
    }
}

sub get_markers {
    my ( $scfref, $samples, $genetics ) = @_;

    my %markers;
    my %prevpos;
    foreach my $snp ( @{$scfref} ) {
        chomp $snp;
        my ( $marker, $type, $parentcall, $pos, $info ) = parse_snp( $snp, $samples, $genetics );
        if ( $type ne "Reject" && $info->{'FS'} > FS_THRESHOLD ) {
            $info->{error} = "Fails FS threshold: Type $type";
            $type = "Reject";
        }
        if ( $type ne "Reject" && $info->{'MQ'} < MQ_THRESHOLD ) {
            $info->{error} = "Fails MQ threshold: Type $type";
            $type = "Reject";
        }

        check_phase( $marker, $markers{$type}{ $prevpos{$type} }{marker}, $type )
          if ( $type ne "Reject" && $prevpos{$type} );
        $markers{$type}{$pos}{marker}     = $marker;
        $markers{$type}{$pos}{parent}     = $parentcall;
        $markers{$type}{$pos}{mq}         = $info->{'MQ'};
        $markers{$type}{$pos}{fs}         = $info->{'FS'};
        $markers{$type}{$pos}{error}      = $info->{error} // "";
        $markers{$type}{$pos}{rmsobs}     = $info->{rmsobs};
        $markers{$type}{$pos}{rmspattern} = $info->{rmspattern};
        $markers{$type}{$pos}{p}          = $info->{p};
        $markers{$type}{$pos}{pgqs}       = $info->{pgqs};
        $markers{$type}{$pos}{pdps}       = $info->{pdps};
        $prevpos{$type}                   = $pos;
    }
    \%markers;
}

sub check_phase {
    my ( $cur, $prev, $type ) = @_;
    my $phasea = 0;
    my $phaseb = 0;
    $type =~ s/\-(.+)//;
    foreach my $sample ( keys %{$cur} ) {
        if ( defined $swapphase{$type}{ $cur->{$sample}{gt} } ) {
            $phasea++ if ( $cur->{$sample}{gt} eq $prev->{$sample}{gt} );
            $phaseb++
              if $swapphase{$type}{ $cur->{$sample}{gt} } eq $prev->{$sample}{gt};
        }
    }

    if ( $phaseb > $phasea ) {
        foreach my $sample ( keys %{$cur} ) {
            $cur->{$sample}{gt} = $swapphase{$type}{ $cur->{$sample}{gt} } // $cur->{$sample}{gt};
        }
    }
    return;
}

sub parse_snp {
    my ( $snp, $samples, $genetics ) = @_;

    my @f     = split /\t/, $snp;
    my $callf = $f[8];
    my $pos   = $f[1];
    my @parentcalls =
      map { get_cp( 'GT', $f[ $samples->{parents}{lookup}{$_} ], $callf ) } @{ $genetics->{parents} };

    my ( $parentcall, $vgt ) = get_parent_call( \@parentcalls );

    my @info = split /;/, $f[7];
    my %info;
    map { $info{$1} = $2 if (/^(.+)=(.+)$/); } @info;

    # Get parental GQs and DPs
    $info{pgqs} = join ':',
      map { sprintf "%2d", get_cp( 'GQ', $f[ $samples->{parents}{lookup}{$_} ], $callf ) } @{ $genetics->{parents} };

    $info{pdps} = join ':',
      map { sprintf "%3d", get_cp( 'DP', $f[ $samples->{parents}{lookup}{$_} ], $callf ) } @{ $genetics->{parents} };

    my %marker;

    if ( !defined $genetics->{types}{$parentcall} ) {
        $info{error} = "Not a valid parent call: @parentcalls";

        foreach my $sample ( @{ $samples->{offspring}{order} } ) {
            my $opos = $samples->{offspring}{lookup}{$sample};
            my $call = get_cp( 'GT', $f[$opos], $callf );
            my $a    = substr $call, 0, 1;
            my $b    = substr $call, -1;
            $marker{$sample}{gt} = $a eq '.' && $b eq '.' ? '.' : $a eq $b ? $a : 'H';
            $marker{$sample}{gq} = get_cp( 'GQ', $f[$opos], $callf );
        }
        return ( \%marker, "Reject", $parentcall, $pos, \%info );
    }

    my %gtsbysex;
    my %invalid_samples;

    foreach my $sample ( @{ $samples->{offspring}{order} } ) {
        my $opos = $samples->{offspring}{lookup}{$sample};
        my $call = get_cp( 'GT', $f[$opos], $callf );
        my $gt   = $vgt->{$call} // "X";
        $invalid_samples{$sample} = $call if $gt eq "X";
        push @{ $gtsbysex{ $genetics->{samplesex}{$sample} } }, $gt;
        $marker{$sample}{gt} = $gt;
        $marker{$sample}{gq} = get_cp( 'GQ', $f[$opos], $callf );
    }

    if ( keys %invalid_samples ) {
        $info{error} = "Invalid calls: ";
        map { $info{error} .= "$_:$invalid_samples{$_} " }
          sort keys %invalid_samples;
        return ( \%marker, "Reject", $parentcall, $pos, \%info );
    }

    my ( %types, %rms_obs, %rms_pattern );
    my @valid_types;
    for my $type ( keys %{ $genetics->{types}{$parentcall} } ) {
        ( $types{$type}, $rms_obs{$type}, $rms_pattern{$type} ) = run_rms_test(
            \@{ $gtsbysex{'Male'} },
            \@{ $gtsbysex{'Female'} },
            $genetics->{types}{$parentcall}{$type}{'Male'},
            $genetics->{types}{$parentcall}{$type}{'Female'},
            $genetics->{rms}
        );
        push @valid_types, $type if $types{$type} >= 0.05;
    }

    if ( @valid_types ne 1 ) {
        my $typestring;
        map {
            $typestring .= "$_:$types{$_}";
            $typestring .= "*" if $types{$_} >= 0.05;
            $typestring .= " "
        } sort { $types{$a} <=> $types{$b} } keys %types;
        chop $typestring;
        $info{error} =
            @valid_types == 0 ? "No valid types: "
          : @valid_types > 1  ? "Too many valid types: "
          :                     "";
        $info{error} .= $typestring;
        return ( \%marker, "Reject", $parentcall, $pos, \%info );
    }

    my $type = (@valid_types)[0];
    if ( $type eq "Reject" ) {
        $info{error} = "No valid type found";
        return ( \%marker, "Reject", $parentcall, $pos, \%info );
    }

    my @pgqs = split /:/, $info{pgqs};
    my @pdps = split /:/, $info{pdps};
    my @pcs  = split //,  $parentcall;
    my $poorqual = 0;
    my $highdp   = 0;
    map {
        my $threshold = $pcs[$_] eq 'H' ? 99 : 60;
        $poorqual = 1 if $pgqs[$_] < $threshold;
        $highdp   = 1 if $pdps[$_] >= 85;
    } 0 .. 2;

    if ($poorqual) {
        $info{error} = "Poor quality parental calls";
        return ( \%marker, "Reject", $parentcall, $pos, \%info );
    }

    if ($highdp) {
        $info{error} = "Depth>100 for at least one parent";
        return ( \%marker, "Reject", $parentcall, $pos, \%info );
    }

    $info{rmsobs}     = $rms_obs{$type};
    $info{rmspattern} = $rms_pattern{$type};
    $info{p}          = $types{$type};

    # Correct eg B/- to H
    if ( defined $genetics->{types}{$parentcall}{$type}{corrections} ) {
        foreach my $sample ( @{ $samples->{offspring}{order} } ) {
            $marker{$sample}{gt} =
              $genetics->{types}{$parentcall}{$type}{corrections}{ $marker{$sample}{gt} } // $marker{$sample}{gt};
        }
    }

    return ( \%marker, $type, $parentcall, $pos, \%info );
}

sub run_rms_test {
    my ( $malecalls, $femalecalls, $validmale, $validfemale, $rms ) = @_;

    my %exp;
    my %male_p;
    my %female_p;
    get_expected_classes( \%exp, 'M', $malecalls,   $validmale,   \%male_p );
    get_expected_classes( \%exp, 'F', $femalecalls, $validfemale, \%female_p );

    my $n = @{$malecalls} + @{$femalecalls};
    my %obs;
    map { $obs{"M$_"}++ } @{$malecalls};
    map { $obs{"F$_"}++ } @{$femalecalls};

    my $rms_obs     = get_rms( \%obs, \%exp );
    my $rms_pattern = "$validmale:$validfemale";
    my $p           = 0.99;
    while ( $p > 0 && $rms_obs >= $rms->{$rms_pattern}{$p} ) {
        $p -= 0.01;
        $p = sprintf "%.2f", $p;
    }

    return ( $p, $rms_obs, $rms->{$rms_pattern}{$p} );
}

sub get_rms {
    my ( $obs, $exp ) = @_;

    my $rms;
    for my $c ( 'MA', 'MB', 'MH', 'M.', 'FA', 'FB', 'FH', 'F.' ) {
        my $obsv = $obs->{$c} // 0;
        my $expv = $exp->{$c} // 0;
        $rms += ( $obsv - $expv )**2;
    }
    $rms = sqrt( $rms / 8 );

    return $rms;
}

sub get_expected_classes {
    my ( $exp_ref, $sex, $calls, $valid, $p ) = @_;
    my @exp = split /,/, $valid;
    my $shares = 0;
    for my $e (@exp) {
        if ( $e =~ /^(\d)?([ABH.])$/ ) {
            my $share = $1 // 1;
            $exp_ref->{"$sex$2"} = $share;
            $shares += $share;
        }
    }
    my $cum_prob = 0;
    for my $e (@exp) {
        $e = $2 if ( $e =~ /^(\d)([ABH.])$/ );
        my $prob = $exp_ref->{"$sex$e"} / $shares;
        $exp_ref->{"$sex$e"} = $prob * @{$calls};
        $p->{"$e"}{min} = $cum_prob;
        $cum_prob += $prob;
        $p->{"$e"}{max} = $cum_prob;
    }
}

sub get_parent_call {
    my ($parents) = @_;

    my @hom;
    my %hom;
    my %het;

    # Get unique homozygote alleles and heterozygote calls
    for my $pc ( @{$parents} ) {
        my ( $i, $j ) = $pc =~ /(.)\/(.)/;
        if ( $i eq $j ) {
            push @hom, $i if !defined $hom{$i};
            $hom{$i}++;
        }
        else {
            $het{$pc}++;
        }
    }

    return ( "XXX", 0 )
      if ( defined $hom{'.'} && keys %hom > 3
        or !defined $hom{'.'} && keys %hom > 2
        or keys %het > 1 );

    # At this point, we have up to two homozygous calls in @hom
    # and at most one heterozygous call in keys %het

    my @valleles = ( "A", "B" );

    # Assign 0,1,2 etc alleles to A and B symbols
    my %vallele;
    for my $hom (@hom) {
        next if ( $hom eq '.' );
        $vallele{$hom} = shift @valleles;
    }

    if ( @valleles && keys %het ) {
        for my $allele ( ( keys %het )[0] =~ /(.)\/(.)/ ) {
            $vallele{$allele} = shift @valleles
              if ( !defined $vallele{$allele} );
        }
    }
    return ( "XXX", 0 ) if keys %vallele > 2;

    my @vas = sort { $a <=> $b } keys %vallele;
    my %vgt;
    if ( defined $vas[0] ) {
        $vgt{"$vas[0]/$vas[0]"} = "$vallele{$vas[0]}";
        if ( defined $vas[1] ) {
            $vgt{"$vas[1]/$vas[1]"} = "$vallele{$vas[1]}";
            $vgt{"$vas[0]/$vas[1]"} = 'H';
        }
    }
    $vgt{"./."}   = '.';
    $vallele{"."} = ".";

    # Assign A, B and H (for heterozygote) to parental calls
    my $parentcall;
    for my $pc ( @{$parents} ) {
        my ( $i, $j ) = $pc =~ /(.)\/(.)/;
        $parentcall .= ( $i eq $j ) ? $vallele{$i} : 'H';
    }

    return ( $parentcall, \%vgt );
}

sub get_cp {
    my ( $part, $call, $order ) = @_;

    return $nullcall{$part} if $call =~ /^\.\/\./;

    my @cp = split ':', $call;
    if ( !defined $callorder{$order} ) {
        my @neworder = split ':', $order;
        map { $callorder{$order}{ $neworder[$_] } = $_ } 0 .. $#neworder;
    }
    my $callpart = $cp[ $callorder{$order}{$part} ];
    return $callpart eq '.' ? 0 : $callpart;
}

sub find_edges {
    my ( $markers, $samples, $genetics ) = @_;

    for my $type ( keys %{$markers} ) {
        next if ( $type eq "Reject" );
        my @pos = sort { $a <=> $b } keys %{ $markers->{$type} };
        for my $sample ( @{ $samples->{offspring}{order} } ) {

            for my $gt ( 'A', 'B', 'H', '.' ) {
                my %maskcall = ( 'A' => -1, 'B' => -1, 'H' => -1, '.' => -1 );
                $maskcall{$gt} = 1;
                get_sample_blocks( $markers->{$type}, $sample, \@pos, \%maskcall, MASKLEN );
            }

            call_blocks( $markers->{$type}, $sample, \@pos, MASKLEN, $type, $genetics );
            clean_blocks( $markers->{$type}, $sample, \@pos, $type );
        }
    }
}

sub call_blocks {
    my ( $marker, $sample, $pos, $masklen, $type, $genetics ) = @_;

    my $blocknum = 0;
    my @blocks   = ();
    my $validgts;
    my $prevpos = -100;
    for my $p ( @{$pos} ) {
        $validgts = $genetics->{types}{ $marker->{$p}{parent} }{$type}{ $genetics->{samplesex}{$sample} };

        # Create a new block if the edge value is greater than
        # the mask length
        $blocknum++
          if ( defined $marker->{$p}{marker}{$sample}{edge}
            && abs( $marker->{$p}{marker}{$sample}{edge} ) > $masklen );
        push @{ $blocks[$blocknum]{pos} }, $p;
        if (
            $validgts =~ /$marker->{$p}{marker}{$sample}{gt}/
            && (   $p - $prevpos > 100
                || $marker->{$prevpos}{marker}{$sample}{gt} ne $marker->{$p}{marker}{$sample}{gt} )
          )
        {
            $blocks[$blocknum]{gt}{ $marker->{$p}{marker}{$sample}{gt} }++;
            $prevpos = $p;
        }
    }

    my @blocksbysize =
      sort { @{ $blocks[$b]{pos} } <=> @{ $blocks[$a]{pos} } } 0 .. $#blocks;

    for my $i (@blocksbysize) {

        my %gts;
        my $size;

        my @sortgt =
          sort { $blocks[$i]{gt}{$b} <=> $blocks[$i]{gt}{$a} }
          keys %{ $blocks[$i]{gt} };
        my $maxgt =
            @sortgt == 0                                                      ? '~'
          : @sortgt == 1                                                      ? $sortgt[0]
          : $blocks[$i]{gt}{ $sortgt[0] } / $blocks[$i]{gt}{ $sortgt[1] } > 2 ? $sortgt[0]
          :                                                                     '~';
        map { $marker->{$_}{marker}{$sample}{cons} = $maxgt } @{ $blocks[$i]{pos} };
    }

}

sub clean_blocks {
    my ( $marker, $sample, $pos, $type ) = @_;

    my $blocknum = 0;
    my @blocks;
    my $prev_cons;
    for my $p ( @{$pos} ) {
        $blocknum++
          if ( defined $prev_cons
            && $prev_cons ne $marker->{$p}{marker}{$sample}{cons} );
        push @{ $blocks[$blocknum]{pos} }, $p;
        $blocks[$blocknum]{cons} = $marker->{$p}{marker}{$sample}{cons};
        $prev_cons = $marker->{$p}{marker}{$sample}{cons};
    }

    for my $b ( 0 .. $#blocks ) {
        if (
               @blocks >= 3
            && $b >= 1
            && $b <= $#blocks - 1
            && $blocks[ $b - 1 ]{cons} eq $blocks[ $b + 1 ]{cons}
            && $blocks[ $b - 1 ]{cons} ne $blocks[$b]{cons}
            && (   ( ( max( @{ $blocks[$b]{pos} } ) - min( @{ $blocks[$b]{pos} } ) ) < ERRORBLOCK )
                || ( @{ $blocks[$b]{pos} } <= ERRORSNPS ) )
          )
        {
            for my $p ( @{ $blocks[$b]{pos} } ) {
                $marker->{$p}{marker}{$sample}{corrected} =
                  $blocks[ $b - 1 ]{cons};
            }
            $blocks[$b]{cons} = $blocks[ $b - 1 ]{cons};
        }
        else {
            for my $p ( @{ $blocks[$b]{pos} } ) {
                $marker->{$p}{marker}{$sample}{corrected} =
                  $marker->{$p}{marker}{$sample}{cons};
            }
        }
    }
}

sub get_sample_blocks {
    my ( $marker, $sample, $pos, $maskcall, $masklen ) = @_;

    @mask = ( (-1) x $masklen, 0, (1) x $masklen ) if ( !@mask );
    my @called;
    my $prevpos = -100;
    foreach my $p ( @{$pos} ) {

        # Check if previous and current SNP are separated by a large gap
        if ( $prevpos > -100 and $p - $prevpos >= GAPSIZE ) {
            $marker->{$p}{marker}{$sample}{edge} = $masklen + 1;
        }
        next if !defined $marker->{$p}{marker}{$sample}{gt};
        next if $maskcall->{ $marker->{$p}{marker}{$sample}{gt} } == 0;

        if (   $p - $prevpos > 100
            || $marker->{$p}{marker}{$sample}{gt} ne $marker->{$prevpos}{marker}{$sample}{gt} )
        {
            push @called, $p;
            $prevpos = $p;
        }
    }

    foreach my $i ( $masklen .. $#called - $masklen ) {

        my @maskcallpos = @called[ $i - $masklen .. $i + $masklen ];

        my $edgesum = 0;

        my $edge =
          int sum map { $maskcall->{ $marker->{ $maskcallpos[$_] }{marker}{$sample}{gt} } * $mask[$_] } 0 .. $#mask;

        $marker->{ $called[$i] }{marker}{$sample}{edge} = $edge
          if !defined $marker->{ $called[$i] }{marker}{$sample}{edge}
          or $edge > $marker->{ $called[$i] }{marker}{$sample}{edge};
    }
}

our $merge = sub {
    my ( $part, $all ) = @_;
    foreach my $scf ( keys %{$part} ) {
        $all->{scf}{$scf}++;
    }
};

our $output = sub {
    my ( $data, $genome, $outfix ) = @_;
    print STDERR "Outputting...\n";

    my %outfile;
    open $outfile{'markers'}, '>', "$outfix.markers.out"
      or croak "Can't open $outfix.markers.out: $OS_ERROR\n";
    open $outfile{'blocks'}, '>', "$outfix.blocks.out"
      or croak "Can't open $outfix.blocks.out: $OS_ERROR\n";

    foreach my $scf ( sort keys %{ $data->{scf} } ) {
        foreach my $filetype ( 'markers', 'blocks' ) {
            open my $scffile, '<', "$scf.$filetype.tmp.out"
              or croak "Can't open scaffold $filetype output for $scf: $OS_ERROR\n";
            while ( my $scfline = <$scffile> ) {
                print { $outfile{$filetype} } $scfline;
            }
            close $scffile;
            system("rm $scf.$filetype.tmp.out");
        }
    }
    close $outfile{'markers'};
    close $outfile{'blocks'};
};

sub output_markers_to_file {
    my ( $scf, $data ) = @_;
    open my $markerhandle, '>', "$scf.markers.tmp.out"
      or croak "Can't open markers output file for $scf: $OS_ERROR\n";
    output_scf_markers( $scf, $data, $markerhandle );
    close $markerhandle;
}

sub output_scf_markers {
    my ( $scf, $data, $handle ) = @_;
    foreach my $type ( sort keys %{$data} ) {
        foreach my $pos ( sort { $a <=> $b } keys %{ $data->{$type} } ) {

            print $handle
"$scf\t$pos\t$type\t$data->{$type}{$pos}{parent}\t$data->{$type}{$pos}{pgqs}\t$data->{$type}{$pos}{pdps}\t$data->{$type}{$pos}{mq}\t$data->{$type}{$pos}{fs}\t";

            if ( $type ne "Reject" ) {
                printf $handle "%.2f\t%.2f\t%.2f\t",
                  $data->{$type}{$pos}{p},
                  $data->{$type}{$pos}{rmsobs},
                  $data->{$type}{$pos}{rmspattern};
            }

            my @gt = split //, $data->{$type}{$pos}{marker}{gt};
            foreach my $i ( 0 .. $#gt ) {
                my $col =
                  ceil( 255 - $data->{$type}{$pos}{marker}{gq}[$i] / 4.2 );
                print $handle fg $col, $gt[$i];
            }
            print $handle "\t";

            if ( $type eq "Reject" ) {
                print $handle "$data->{$type}{$pos}{error}\n";
                next;
            }

            my @edge      = @{ $data->{$type}{$pos}{marker}{edge} };
            my @cons      = split //, $data->{$type}{$pos}{marker}{cons};
            my @corrected = split //, $data->{$type}{$pos}{marker}{corrected};

            if ( @cons == 0 or @corrected == 0 ) {
                print $handle "\n";
                next;
            }
            my $cons      = "";
            my $corrected = "";
            foreach my $e ( 0 .. $#edge ) {
                my $col = $edge[$e] ne '.'
                  && $edge[$e] > MASKLEN ? 'red1' : 'black';
                print $handle fg $col, $edge[$e];
                $cons      .= fg $col, $cons[$e];
                $corrected .= fg $col, $corrected[$e];
            }
            print $handle "\t$cons\t$corrected\n";
        }
        print $handle '-' x 253, "\n";
    }
    print $handle '-' x 253, "\n";

}

sub output_blocks_to_file {
    my ( $scf, $data, $samples, $genetics, $scfl ) = @_;
    open my $blockhandle, '>', "$scf.blocks.tmp.out"
      or croak "Can't open blocks output file for $scf: $OS_ERROR\n";
    output_scf_blocks( $scf, $data, $blockhandle, $genetics, $samples, $scfl );
    close $blockhandle;
}

sub output_scf_blocks {
    my ( $scf, $data, $blockhandle, $genetics, $samples, $scfl ) = @_;

    my $samplenum = @{ $samples->{offspring}{order} };
    my @types     = sort keys %{ $genetics->{masks} };

    printf $blockhandle "%8s\t%8s\t%8s\t%8s", 'Scaffold', 'Start', 'End', 'Length';

    map { printf $blockhandle "\t%-${samplenum}s", $_; } @types;
    print $blockhandle "\n";

    my $empty = ' ' x $samplenum;

    my %scfpos;
    foreach my $type ( keys %{$data} ) {
        next if ( $type eq "Reject" );
        my $prevpos = 0;
        foreach my $pos ( sort { $a <=> $b } keys %{ $data->{$type} } ) {
            $scfpos{$pos}{types}{$type} =
              $data->{$type}{$pos}{marker}{corrected};
            $scfpos{$pos}{end} = $pos;

            if (   $prevpos eq 0
                or $scfpos{$prevpos}{types}{$type} ne $scfpos{$pos}{types}{$type} )
            {
                $scfpos{ $prevpos + 1 }{types}{$type} = $empty;
            }
            $prevpos = $pos;
        }
        $scfpos{ $prevpos + 1 }{types}{$type} = $empty;
    }

    my %curpat;
    map { $curpat{$_} = $empty } @types;
    my $blockpos = 1;
    my $lastpos;
    for my $p ( sort { $a <=> $b } keys %scfpos ) {
        for my $t (@types) {
            if ( defined $scfpos{$p}{types}{$t}
                and $scfpos{$p}{types}{$t} ne $curpat{$t} )
            {
                output_block( $scf, $blockpos, $p - 1, $blockhandle, \@types, \%curpat, $empty );
                $curpat{$t} = $scfpos{$p}{types}{$t};
                $blockpos = $p;
            }
        }
        $lastpos = $p;
    }
    output_block( $scf, $blockpos, $scfl->{$scf}, $blockhandle, \@types, \%curpat, $empty );
}

sub output_block {
    my ( $scf, $start, $end, $blockhandle, $types, $curpat, $empty ) = @_;
    my $length = $end - $start + 1;
    printf $blockhandle "%8s\t%8d\t%8d\t%8d", $scf, $start, $end, $length;
    for my $type ( @{$types} ) {
        print $blockhandle "\t";
        print $blockhandle $curpat->{$type} // $empty;
    }
    print $blockhandle "\n";
}

sub get_next_pat {
    my ( $i, $type, $scfpos, $pos ) = @_;
    while ( $i < @{$pos} ) {
        return $scfpos->{ $pos->[$i] }{types}{$type}
          if defined $scfpos->{ $pos->[$i] }{types}{$type};
        $i++;
    }
    return "";
}

1;
