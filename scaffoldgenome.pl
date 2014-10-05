use List::Util qw/sum min max/;
use POSIX qw/ceil/;
use Getopt::Long qw/GetOptionsFromString/;
use DBD::SQLite;

use strict;
use warnings;

use constant MASKLEN      => 4;
use constant ERRORBLOCK   => 2000;
use constant ERRORSNPS    => 5;
use constant FS_THRESHOLD => 5;
use constant MQ_THRESHOLD => 90;
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
    my $args = shift;

    my $genetics_filename = "";
    my $agp_filename      = "";
    my $rms_filename      = "";
    my $options_okay      = GetOptionsFromString(
        $args->{extraargs},
        'genetics=s'    => \$genetics_filename,
        'agp=s'         => \$agp_filename,
        'rmsfilename=s' => \$rms_filename,
    );
    croak "Can't process user options: $OS_ERROR\n" if !$options_okay;

    croak "No genetics file! Please specify -g to define parents, poor quality individuals and marker types\n"
      if ( $genetics_filename eq "" );

    croak "No AGP file! Please specify -a\n"
      if ( $agp_filename eq "" );

    my $genetics = load_genetics( $genetics_filename, $rms_filename );

    my $samples = get_samples( $args->{vcf_filename}, $genetics );

    my $agp = load_agp($agp_filename);

    my $dbfilename = create_output_database( $args->{output_prefix}, $genetics );

    { genetics => $genetics, samples => $samples, agp => $agp, dbfilename => $dbfilename };
};

sub load_genetics {
    my ( $geneticsfilename, $rmsfilename ) = @_;

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

sub load_agp {
    my $agp_filename = shift;

    open my $agp_file, "<", $agp_filename or croak "Can't open AGP file $agp_filename!\n";
    my %agp;
    while ( my $agp_line = <$agp_file> ) {
        my ( $scf, $start, $end, $part, $type, @f ) = split /\t/, $agp_line;
        $agp{$scf}{$part}{start} = $start;
        $agp{$scf}{$part}{end}   = $end;
        $agp{$scf}{$part}{type}  = $type;
    }
    close $agp_file;

    return \%agp;
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

sub create_output_database {
    my ( $output_prefix, $genetics ) = @_;

    my $dbfilename = "$output_prefix.db";
    unlink $dbfilename if ( -e $dbfilename );

    my $dbh = DBI->connect( "dbi:SQLite:dbname=$dbfilename", "", "" );

    my $sth = $dbh->prepare(
        "CREATE TABLE markers
                 (scaffold text, position integer, marker_type text,
                  parent_gt text, parent_gqs text, parent_dps text,
                  mq real, fs real,
                  p real, rms_obs real, rms_pattern real, phase integer,
                  pattern text, edge text, consensus text, corrected text,
                  error text)"
    );
    $sth->execute;

    my $statement = "CREATE TABLE blocks (scaffold text, start integer, end integer, length integer";
    map { $statement .= ", \"$_\" text" } sort keys $genetics->{masks};
    $statement .= ")";

    $sth = $dbh->prepare($statement);
    $sth->execute;

    $dbh->disconnect;

    return $dbfilename;
}

our $process = sub {
    my ( $scf, $scfref, $data, $userdata, $scfl ) = @_;
    $data->{$scf} = get_markers( $scfref, $userdata->{samples}, $userdata->{genetics} );
    find_edges( $data->{$scf}, $scf, $userdata );
    collapse( $data->{$scf}, $userdata->{samples} );
    output_to_db( $scf, $data->{$scf}, $userdata, $scfl );
    $data->{$scf} = "";
    return;
};


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

        $markers{$type}{$pos}{phase} = 0;
        $markers{$type}{$pos}{phase} = check_phase( $marker, $markers{$type}{ $prevpos{$type} }{marker}, $type )
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
    my $phase  = 0;
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
        $phase = 1;
        foreach my $sample ( keys %{$cur} ) {
            $cur->{$sample}{gt} = $swapphase{$type}{ $cur->{$sample}{gt} } // $cur->{$sample}{gt};
        }
    }
    return $phase;
}

sub parse_snp {
    my ( $snp, $samples, $genetics ) = @_;

    my @f     = split /\t/, $snp;
    my $callf = $f[8];
    my $pos   = $f[1];

    my $calls = get_calls( \@f, $samples, $genetics );

    my @parentcalls = map { $calls->{$_}{'GT'} } @{ $genetics->{parents} };

    my ( $parentcall, $vgt ) = get_parent_call( \@parentcalls );

    my @info = split /;/, $f[7];
    my %info;
    map { $info{$1} = $2 if (/^(.+)=(.+)$/); } @info;

    # Get parental GQs and DPs
    $info{pgqs} = join ':', map { sprintf "%2d", $calls->{$_}{'GQ'} } @{ $genetics->{parents} };

    $info{pdps} = join ':', map { sprintf "%3d", $calls->{$_}{'DP'} } @{ $genetics->{parents} };

    my %marker;

    if ( !defined $genetics->{types}{$parentcall} ) {
        $info{error} = "Not a valid parent call: @parentcalls";

        foreach my $sample ( @{ $samples->{offspring}{order} } ) {
            my $call = $calls->{$sample}{'GT'};
            my $a    = substr $call, 0, 1;
            my $b    = substr $call, -1;
            $marker{$sample}{gt} = $a eq '.' && $b eq '.' ? '.' : $a eq $b ? $a : 'H';
            $marker{$sample}{gq} = $calls->{$sample}{'GQ'};
        }
        return ( \%marker, "Reject", $parentcall, $pos, \%info );
    }

    my %gtsbysex;
    my %invalid_samples;

    foreach my $sample ( @{ $samples->{offspring}{order} } ) {
        my $call = $calls->{$sample}{'GT'};
        my $gt = $vgt->{$call} // "X";
        $invalid_samples{$sample} = $call if $gt eq "X";
        push @{ $gtsbysex{ $genetics->{samplesex}{$sample} } }, $gt;
        $marker{$sample}{gt} = $gt;
        $marker{$sample}{gq} = $calls->{$sample}{'GQ'};
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

sub get_calls {
    my ( $f, $samples, $genetics ) = @_;

    my %calls;
    my @parts = split ':', $f->[8];
    for my $sample ( @{ $genetics->{parents} } ) {
        my $call = $f->[ $samples->{parents}{lookup}{$sample} ];
        if ( $call =~ /^\./ ) {
            map { $calls{$sample}{$_} = $nullcall{$_} } @parts;
        }
        else {
            my @cf = split ':', $call;
            for my $partnum ( 0 .. $#parts ) {
                $calls{$sample}{ $parts[$partnum] } = $cf[$partnum] eq '.' ? 0 : $cf[$partnum];
            }
        }
    }

    for my $sample ( @{ $samples->{offspring}{order} } ) {
        my $call = $f->[ $samples->{offspring}{lookup}{$sample} ];
        if ( $call =~ /^\./ ) {
            map { $calls{$sample}{$_} = $nullcall{$_} } @parts;
        }
        else {
            my @cf = split ':', $call;
            for my $partnum ( 0 .. $#parts ) {
                $calls{$sample}{ $parts[$partnum] } = $cf[$partnum] eq '.' ? 0 : $cf[$partnum];
            }
        }
    }
    return \%calls;
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

sub find_edges {
    my ( $markers, $scf, $userdata ) = @_;

    for my $type ( keys %{$markers} ) {
        next if ( $type eq "Reject" );
        my @pos = sort { $a <=> $b } keys %{ $markers->{$type} };
        for my $sample ( @{ $userdata->{samples}{offspring}{order} } ) {

            for my $gt ( 'A', 'B', 'H', '.' ) {
                my %maskcall = ( 'A' => -1, 'B' => -1, 'H' => -1, '.' => -1 );
                $maskcall{$gt} = 1;
                get_sample_blocks( $markers->{$type}, $sample, \@pos, \%maskcall, MASKLEN, $userdata->{agp}{$scf} );
            }

            call_blocks( $markers->{$type}, $sample, \@pos, MASKLEN, $type, $userdata->{genetics} );
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

        my $mref = $marker->{$p}{marker}{$sample};

        # Create a new block if the edge value is greater than
        # the mask length
        $blocknum++
          if ( defined $mref->{edge}
            && abs( $mref->{edge} ) > $masklen );
        push @{ $blocks[$blocknum]{pos} }, $p;
        if (
            $validgts =~ /$mref->{gt}/
            && (   $p - $prevpos > 100
                || $marker->{$prevpos}{marker}{$sample}{gt} ne $mref->{gt} )
          )
        {
            my $numgts = scalar( split ',', $validgts );
            my $het_weight =
                ($numgts > 2 and $mref->{gt} eq 'H') ? 2
              : ($numgts == 2 and $mref->{gt} eq 'H' and $marker->{$p}{phase} == 0) ? 2
              : ($numgts == 2 and $mref->{gt} ne 'H' and $marker->{$p}{phase} == 1) ? 2
              :                                                                     1;
            $blocks[$blocknum]{gt}{ $mref->{gt} } += $het_weight;
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

sub get_part {
    my ( $pos, $agp ) = @_;

    map { return $_ if $pos >= $agp->{$_}{start} && $pos <= $agp->{$_}{end} } keys %{$agp};
}

sub get_sample_blocks {
    my ( $marker, $sample, $pos, $maskcall, $masklen, $scf_agp ) = @_;

    @mask = ( (-1) x $masklen, 0, (1) x $masklen ) if ( !@mask );
    my @called;
    my $prevpos = -100;
    foreach my $p ( @{$pos} ) {
        my $mref = $marker->{$p}{marker}{$sample};

        # Check if previous and current SNP are separated by a large gap
        if ( $prevpos > -100
            and ( get_part( $prevpos, $scf_agp ) != get_part( $p, $scf_agp ) or $p - $prevpos >= GAPSIZE ) )
        {
            $mref->{edge} = $masklen + 1;
        }
        next if !defined $mref->{gt} or $maskcall->{ $mref->{gt} } == 0;

        if (   $p - $prevpos > 100
            || $mref->{gt} ne $marker->{$prevpos}{marker}{$sample}{gt} )
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

our $merge = sub {
    my ( $part, $all ) = @_;
    foreach my $scf ( keys %{$part} ) {
        $all->{scf}{$scf}++;
    }
};

our $output = sub {
    my ( $data, $genome, $outfix, $userdata ) = @_;
    print STDERR "Creating indices\n";

    my $dbh = DBI->connect( "dbi:SQLite:dbname=$userdata->{dbfilename}", "", "" );
    my $sth = $dbh->prepare("create index scf_marker_type on markers (scaffold, marker_type)");
    $sth->execute;

    $sth = $dbh->prepare("create index scaffold on blocks (scaffold)");
    $sth->execute;

    $dbh->disconnect;

    print STDERR "Done\n";
};

sub output_to_db {
    my ( $scf, $data, $userdata, $scfl ) = @_;

    my $dbh = DBI->connect( "dbi:SQLite:dbname=$userdata->{dbfilename}", "", "", { AutoCommit => 0 } );

    output_markers_to_db( $scf, $data, $userdata, $dbh );
    output_blocks_to_db( $scf, $data, $userdata, $scfl, $dbh );

    $dbh->disconnect;
    return;
}

sub output_markers_to_db {
    my ( $scf, $data, $userdata, $dbh ) = @_;

    my $insert_handle = $dbh->prepare_cached(
        'INSERT INTO markers VALUES (?,?,?,?,
                                     ?,?,?,?,
                                     ?,?,?,?,
                                     ?,?,?,?,?
                                     )'
    );

    foreach my $type ( keys %{$data} ) {
        foreach my $pos ( keys %{ $data->{$type} } ) {
            my $stp = $data->{$type}{$pos};

            $stp->{edge} = defined $stp->{marker}{edge} ? join '', @{ $stp->{marker}{edge} } : "";
            $stp->{cons}      = $stp->{marker}{cons}      // "";
            $stp->{corrected} = $stp->{marker}{corrected} // "";

            if ( $type eq "Reject" ) {
                $stp->{p}          = -1;
                $stp->{rmsobs}     = -1;
                $stp->{rmspattern} = -1;
            }
        }

        foreach my $pos ( keys %{ $data->{$type} } ) {
            my $stp = $data->{$type}{$pos};

            $insert_handle->execute(
                $scf,         $pos,
                $type,        $stp->{parent},
                $stp->{pgqs}, $stp->{pdps},
                $stp->{mq},   $stp->{fs},
                sprintf( "%4.2f", $stp->{p} ),          sprintf( "%4.2f", $stp->{rmsobs} ),
                sprintf( "%4.2f", $stp->{rmspattern} ), $stp->{phase},
                $stp->{marker}{gt}, $stp->{edge},
                $stp->{cons},       $stp->{corrected},
                $stp->{error}
            );
        }
        $dbh->commit;
    }
    return;
}

sub output_blocks_to_db {
    my ( $scf, $data, $userdata, $scfl, $dbh ) = @_;

    my @types = sort keys %{ $userdata->{genetics}{masks} };

    my $insert_statement = "INSERT INTO blocks VALUES (?,?,?,?";
    map { $insert_statement .= ',?' } @types;
    $insert_statement .= ')';
    my $insert_handle = $dbh->prepare_cached($insert_statement);

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
                $scfpos{ $prevpos + 1 }{types}{$type} = "";
            }
            $prevpos = $pos;
        }
        $scfpos{ $prevpos + 1 }{types}{$type} = "";
    }

    my %curpat;
    map { $curpat{$_} = "" } @types;
    my $blockpos = 1;
    for my $p ( sort { $a <=> $b } keys %scfpos ) {
        for my $t (@types) {
            if ( defined $scfpos{$p}{types}{$t}
                and $scfpos{$p}{types}{$t} ne $curpat{$t} )
            {
                my $start = $blockpos;
                my $end   = $p - 1;
                my $len   = $end - $start + 1;

                my @out_patterns = map { $curpat{$_} } @types;
                $insert_handle->execute( $scf, $start, $end, $len, @out_patterns );
                $curpat{$t} = $scfpos{$p}{types}{$t};
                $blockpos = $p;
            }
        }
        $dbh->commit;
    }
    my @out_patterns = map { $curpat{$_} } @types;
    my $len = $scfl->{$scf} - $blockpos + 1;
    $insert_handle->execute( $scf, $blockpos, $scfl->{$scf}, $len, @out_patterns );
    $dbh->commit;
}

1;
