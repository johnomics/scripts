use List::Util qw/sum min max/;
use POSIX qw/ceil/;
use Getopt::Long qw/GetOptionsFromString/;
use DBD::SQLite;

use strict;
use warnings;

use constant FS_THRESHOLD => 5;
use constant MQ_THRESHOLD => 90;

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
        my ( $parents, $males, $females, $type, $corrections, $recombination ) = split /\t/, $marker_type;
        $genetics{types}{$parents}{$type}{'Male'}   = $males;
        $genetics{types}{$parents}{$type}{'Female'} = $females;

        croak "Recombination value for $type $parents should be Y or N, got $recombination!"
          if $recombination ne 'Y' and $recombination ne 'N';
        if ( $recombination ne "" ) {
            $genetics{types}{$parents}{$type}{recombination} = $recombination;
        }
        $genetics{hmm}{$type}{'Male'}   = $genetics{hmm}{$type}{'Male'}   // make_hmm( $males,   $recombination );
        $genetics{hmm}{$type}{'Female'} = $genetics{hmm}{$type}{'Female'} // make_hmm( $females, $recombination );

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

sub make_hmm {
    my ( $gtstr, $recombination ) = @_;

    my %hmm;

    my %gt_shares;
    my $total_shares = 0;
    my @valid_gts;

    # Get occurrence of gt (eg 'A', '2B' etc) or set share to 0
    for my $gt ( split /,/, $gtstr ) {
        if ( $gt =~ /^(\d)?([ABH.])$/ ) {
            my $share = $1 // 1;
            $gt_shares{$2} = $share;
            $total_shares += $share;
            push @valid_gts, $2;
        }
    }

    for my $gt ( 'A', 'B', 'H', '.' ) {
        $gt_shares{$gt} = 0 if !defined $gt_shares{$gt};
        $hmm{initial}{$gt} = $gt_shares{$gt} / $total_shares;

        for my $gt2 ( 'A', 'B', 'H', '.' ) {
            $hmm{transition}{$gt}{$gt2} = $gt eq $gt2 ? 1 : 0;
            $hmm{emit}{$gt}{$gt2} = 0;
        }
    }

    my $transprob = 0.0001;
    my $emitprob  = 0.0001;
    if ( $recombination eq 'Y' ) {
        for my $i ( 0 .. $#valid_gts - 1 ) {
            my $gt1 = $valid_gts[$i];
            my $gt2 = $valid_gts[ $i + 1 ];
            $hmm{transition}{$gt1}{$gt2} = $transprob;
            $hmm{transition}{$gt2}{$gt1} = $transprob;
            $hmm{transition}{$gt1}{$gt2} -= $transprob;
            $hmm{transition}{$gt2}{$gt1} -= $transprob;
        }
    }

    for my $gt (@valid_gts) {
        my $totalemit = 1;
        for my $gt2 ( 'A', 'B', 'H', '.' ) {
            if ( $gt ne $gt2 ) {
                $hmm{emit}{$gt}{$gt2} = $emitprob;
                $totalemit -= $emitprob;
            }
        }
        $hmm{emit}{$gt}{$gt} = $totalemit;
    }

    return \%hmm;
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
                  pattern text, consensus text,
                  error text)"
    );
    $sth->execute;

    my $statement = "CREATE TABLE blocks (scaffold text, start integer, end integer, length integer";
    map { $statement .= ", \"$_\" text" } sort keys $genetics->{hmm};
    $statement .= ")";

    $sth = $dbh->prepare($statement);
    $sth->execute;

    $dbh->disconnect;

    return $dbfilename;
}

our $process = sub {
    my ( $scf, $scfref, $data, $userdata, $scfl ) = @_;
    $data->{$scf} = get_markers( $scfref, $userdata->{samples}, $userdata->{genetics} );
    find_blocks( $data->{$scf}, $scf, $userdata );
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
        my ( $marker, $type, $parentcall, $pos, $pattern, $info ) = parse_snp( $snp, $samples, $genetics );
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
        $markers{$type}{$pos}{pattern}    = $pattern;
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
        my $pattern = "";
        foreach my $sample ( @{ $samples->{offspring}{order} } ) {
            my $call = $calls->{$sample}{'GT'};
            my $a    = substr $call, 0, 1;
            my $b    = substr $call, -1;
            $marker{$sample}{gt} = $a eq '.' && $b eq '.' ? '.' : $a eq $b ? $a : 'H';
            $marker{$sample}{gq} = $calls->{$sample}{'GQ'};
            $pattern .= $marker{$sample}{gt};
        }
        return ( \%marker, "Reject", $parentcall, $pos, $pattern, \%info );
    }

    my %gtsbysex;
    my %invalid_samples;

    my $pattern;
    foreach my $sample ( @{ $samples->{offspring}{order} } ) {
        my $call = $calls->{$sample}{'GT'};
        my $gt = $vgt->{$call} // "X";
        $invalid_samples{$sample} = $call if $gt eq "X";
        push @{ $gtsbysex{ $genetics->{samplesex}{$sample} } }, $gt;
        $marker{$sample}{gt} = $gt;
        $pattern .= $gt;
        $marker{$sample}{gq} = $calls->{$sample}{'GQ'};
    }

    if ( keys %invalid_samples ) {
        $info{error} = "Invalid calls: ";
        map { $info{error} .= "$_:$invalid_samples{$_} " }
          sort keys %invalid_samples;
        return ( \%marker, "Reject", $parentcall, $pos, $pattern, \%info );
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
        return ( \%marker, "Reject", $parentcall, $pos, $pattern, \%info );
    }

    my $type = (@valid_types)[0];
    if ( $type eq "Reject" ) {
        $info{error} = "No valid type found";
        return ( \%marker, "Reject", $parentcall, $pos, $pattern, \%info );
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
        return ( \%marker, "Reject", $parentcall, $pos, $pattern, \%info );
    }

    if ($highdp) {
        $info{error} = "Depth>100 for at least one parent";
        return ( \%marker, "Reject", $parentcall, $pos, $pattern, \%info );
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

    return ( \%marker, $type, $parentcall, $pos, $pattern, \%info );
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
    my @exp = split ',', $valid;
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

sub find_blocks {
    my ( $markers, $scf, $userdata ) = @_;

    for my $type ( keys %{$markers} ) {
        next if ( $type eq "Reject" );
        my @pos = sort { $a <=> $b } keys %{ $markers->{$type} };

        my @blocks = find_edges( $markers->{$type}, \@pos, $type );
        for my $block (@blocks) {
            for my $sample ( @{ $userdata->{samples}{offspring}{order} } ) {
                get_sample_consensus( $markers->{$type}, $sample, \@pos, $block, $type, $userdata->{genetics} );

#                get_sample_consensus_viterbi( $markers->{$type}, $sample, \@pos, $block, $type, $userdata->{genetics} );
            }
        }
    }
}

sub find_edges {
    my ( $markertype, $pos, $type ) = @_;

    my @blocks;
    my $block_start = $pos->[0];
    my $this        = $markertype->{ $pos->[0] };
    my $length      = length( $this->{pattern} );

    for my $i ( 1 .. $#{$pos} ) {
        my $next = $markertype->{ $pos->[$i] };
        my $distance = hamming( $this->{pattern}, $next->{pattern} );
        if ( $distance > $length * 0.25 ) {
            push @blocks, { start => $block_start, end => $pos->[ $i - 1 ] };
            $block_start = $pos->[$i];
        }
        $this = $next;
    }
    push @blocks, { start => $block_start, end => $pos->[-1] };

    @blocks;
}

sub hamming {
    my ( $a, $b ) = @_;
    my @a = split //, $a;
    my @b = split //, $b;
    my $distance = 0;
    for my $i ( 0 .. length($a) - 1 ) {
        $distance++ if $a[$i] ne '.' and $b[$i] ne '.' and $a[$i] ne $b[$i];
    }
    $distance;
}

sub get_sample_consensus {
    my ( $marker, $sample, $pos, $block, $type, $genetics ) = @_;

    my $sex = $genetics->{samplesex}{$sample};
    my $hmm = $genetics->{hmm}{$type}{$sex};

    my $seqblocks = get_sample_blocks( $marker, $sample, $pos, $block, $hmm );

    if ( @{$seqblocks} <= 1 ) {
        for my $p ( @{$pos} ) {
            next if ( $p < $block->{start} or $block->{end} < $p );
            $marker->{$p}{marker}{$sample}{cons} = $seqblocks->[0]{gt} // '.';
        }
        return;
    }

    # Get block for each position, or neighbouring blocks
    my %posblocks;
    for my $p ( @{$pos} ) {
        next if $p < $block->{start} or $block->{end} < $p;
        for my $i ( 0 .. $#{$seqblocks} ) {
            if ( $p >= $seqblocks->[$i]{start} and $p <= $seqblocks->[$i]{end} ) {
                $posblocks{$p}{block} = $i;
                delete $posblocks{$p}{left} if defined $posblocks{$p}{left};
            }
            else {
                $posblocks{$p}{left} = $i if $p > $seqblocks->[$i]{end};
                $posblocks{$p}{right} = $i
                  if !defined( $posblocks{$p}{block} )
                  and !defined( $posblocks{$p}{right} )
                  and $p < $seqblocks->[$i]{start};
            }
        }
    }

    my $maxseq = get_max_seq($seqblocks);
    my @maxgts = split //, $maxseq;

    for my $p ( keys %posblocks ) {
        $marker->{$p}{marker}{$sample}{cons} =
            defined($posblocks{$p}{block}) ? $maxgts[ $posblocks{$p}{block} ]
          : (defined($posblocks{$p}{left})
          and defined($posblocks{$p}{right})
          and $maxgts[ $posblocks{$p}{left} ] eq $maxgts[ $posblocks{$p}{right} ]) ? $maxgts[ $posblocks{$p}{left} ]
          : '.';
    }
}

sub get_max_seq {
    my $seqblocks = shift;

    my %gts;
    map { $gts{ $_->{gt} } += $_->{length} } @{$seqblocks};
    my ( $gtmax1, $gtmax2 ) = ( reverse sort { $gts{$a} <=> $gts{$b} } keys %gts )[ 0, 1 ];

    my %stateseq;
    for my $breakpoint ( 0 .. $#{$seqblocks} ) {
        my $gt1    = $gtmax1;
        my $gt2    = $gtmax2;
        my $seq1   = "";
        my $seq2   = "";
        my $score1 = 0;
        my $score2 = 0;
        for my $i ( 0 .. $#{$seqblocks} ) {
            if ( $i == $breakpoint ) {
                $gt1 = $gtmax2;
                $gt2 = $gtmax1;
            }
            $seq1 .= $gt1;
            $seq2 .= $gt2;
            $score1 += $seqblocks->[$i]{length} if $seqblocks->[$i]{gt} eq $gt1;
            $score2 += $seqblocks->[$i]{length} if $seqblocks->[$i]{gt} eq $gt2;
        }
        $stateseq{$seq1} = $score1;
        $stateseq{$seq2} = $score2;
    }

    my $maxseq = ( sort { $stateseq{$b} <=> $stateseq{$a} } keys %stateseq )[0];
    return $maxseq;
}

sub get_sample_blocks {
    my ( $marker, $sample, $pos, $block, $hmm ) = @_;

    my @seqblocks;
    my $curgt = "";
    for my $p ( @{$pos} ) {
        next if ( $p < $block->{start} or $block->{end} < $p );
        my $gt = $marker->{$p}{marker}{$sample}{gt};
        next if $hmm->{initial}{$gt} == 0;
        if ( $gt eq $curgt ) {
            $seqblocks[-1]{end} = $p;
        }
        else {
            if ( $curgt ne "" ) {
                $seqblocks[-1]{length} = $seqblocks[-1]{end} - $seqblocks[-1]{start} + 1;
                pop @seqblocks if $seqblocks[-1]{length} < 100;
            }
            $curgt = $gt;
            push @seqblocks, { gt => $curgt, start => $p, end => $p };
        }
    }
    if (@seqblocks) {
        $seqblocks[-1]{length} = $seqblocks[-1]{end} - $seqblocks[-1]{start} + 1;
        pop @seqblocks if $seqblocks[-1]{length} < 100;
    }

    return \@seqblocks if @seqblocks == 0;
    
    my @cleanblocks;
    push @cleanblocks, $seqblocks[0];
    my $i = 0;
    my $j = 1;
    my %gtlens;
    $gtlens{ $seqblocks[$i]{gt} } += $seqblocks[$i]{length};
    while ( $j <= $#seqblocks ) {
        $gtlens{ $seqblocks[$j]{gt} } += $seqblocks[$j]{length};
        if ( $seqblocks[$i]{gt} eq $seqblocks[$j]{gt} ) {
            $cleanblocks[-1]{end} = $seqblocks[$j]{end};
            $cleanblocks[-1]{length} += $seqblocks[$j]{length};
            $j++;
        }
        else {
            $i = $j;
            $j++;
            push @cleanblocks, $seqblocks[$i];
        }
    }
    \@cleanblocks;
}

sub get_sample_consensus_viterbi {
    my ( $marker, $sample, $pos, $block, $type, $genetics ) = @_;
    my $sex    = $genetics->{samplesex}{$sample};
    my $hmm    = $genetics->{hmm}{$type}{$sex};
    my @states = keys %{ $hmm->{initial} };

    my @scores;
    my @paths;

    my $initgt = $marker->{ $block->{start} }{marker}{$sample}{gt};
    my $t      = 0;
    for my $state (@states) {
        $scores[$t]{$state} = $hmm->{initial}{$state} * $hmm->{emit}{$state}{$initgt};
        $paths[$t]{$state}  = $state;
    }

    for my $p ( @{$pos} ) {
        next if ( $p <= $block->{start} or $block->{end} < $p );
        $t++;
        my $total_delta_prob = 0;                                    # For scaling
        my $gt               = $marker->{$p}{marker}{$sample}{gt};

        for my $j (@states) {
            my $delta_prob_j = 0;
            my $psi_prob_j   = 0;
            my $psi_argmax_j = '';

            for my $i (@states) {
                my $psi_prob_ij = $scores[ $t - 1 ]{$i} * $hmm->{transition}{$i}{$j};
                if ( $psi_prob_ij > $psi_prob_j ) {
                    $psi_prob_j   = $psi_prob_ij;
                    $psi_argmax_j = $i;
                }

                my $delta_prob_ij = $psi_prob_ij * $hmm->{emit}{$i}{$gt};
                if ( $delta_prob_ij > $delta_prob_j ) {
                    $delta_prob_j = $delta_prob_ij;
                }
                $total_delta_prob += $delta_prob_ij;
            }
            $scores[$t]{$j} = $delta_prob_j;
            $paths[$t]{$j}  = $psi_argmax_j;
        }

        # Scale probabilities
        map { $scores[$t]{$_} = $scores[$t]{$_} / $total_delta_prob } keys %{ $scores[$t] };
    }
    my $prob = 0;
    my $qT   = "";
    for my $i (@states) {
        if ( $scores[-1]{$i} > $prob ) {
            $prob = $scores[-1]{$i};
            $qT   = $i;
        }
    }

    $marker->{ $block->{end} }{marker}{$sample}{cons} = $qT;

    my $qtplus1 = $qT;
    $t = $#paths - 1;
    for my $p ( reverse @{$pos} ) {
        next if ( $p < $block->{start} or $block->{end} <= $p );
        $marker->{$p}{marker}{$sample}{cons} = $paths[$t]{$qtplus1};
        $qtplus1 = $paths[$t]{$qtplus1};
        $t--;
    }

    for my $p ( @{$pos} ) {
        next if ( $p < $block->{start} or $block->{end} < $p );
    }
    return;
}

sub collapse {
    my ( $scfref, $samples ) = @_;
    foreach my $type ( keys %{$scfref} ) {
        foreach my $pos ( keys %{ $scfref->{$type} } ) {
            my $gt   = "";
            my $cons = "";

            foreach my $sample ( @{ $samples->{offspring}{order} } ) {
                my $ref = $scfref->{$type}{$pos}{marker}{$sample};

                $gt .= $ref->{gt} // '-';

                if ( $type ne "Reject" ) {
                    $cons .= $ref->{cons} // '-';
                }
            }
            delete $scfref->{$type}{$pos}{marker};
            $scfref->{$type}{$pos}{marker}{gt} = $gt;

            if ( $type ne "Reject" ) {
                $scfref->{$type}{$pos}{marker}{cons} = $cons;
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
                                     ?,?,?
                                     )'
    );

    foreach my $type ( keys %{$data} ) {
        foreach my $pos ( keys %{ $data->{$type} } ) {
            my $stp = $data->{$type}{$pos};
            $stp->{cons} = $stp->{marker}{cons} // "";

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
                $stp->{marker}{gt}, $stp->{cons},
                $stp->{error}
            );
        }
        $dbh->commit;
    }
    return;
}

sub output_blocks_to_db {
    my ( $scf, $data, $userdata, $scfl, $dbh ) = @_;

    my @types = sort keys %{ $userdata->{genetics}{hmm} };

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
              $data->{$type}{$pos}{marker}{cons};
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
