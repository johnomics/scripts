use List::Util qw/sum min max/;
use POSIX qw/ceil/;
use Term::ExtendedColor qw/:all/;

use strict;
use warnings;

use constant MASKLEN      => 4;
use constant ERRORBLOCK   => 2000;
use constant ERRORSNPS    => 5;
use constant FS_THRESHOLD => 35;
use constant MQ_THRESHOLD => 57;

my %callorder;

my %chisqcrit = (
    1 => 3.84,
    2 => 5.99,
    3 => 7.82,
    4 => 9.49,
    5 => 11.07,
    6 => 12.59,
    7 => 14.07,
    8 => 15.51,
);

my %nullcall = ( 'GT' => './.', 'GQ' => 0 );

my @mask;

sub process {
    my ( $scf, $scfref, $samples, $data, $genetics ) = @_;
    $data->{$scf} = get_markers( $scfref, $samples, $genetics );
}

sub get_markers {
    my ( $scfref, $samples, $genetics ) = @_;

    my %markers;
    foreach my $snp ( @{$scfref} ) {
        chomp $snp;

        my ( $origparentcall, $inferparentcall, $type, $pos ) =
          parse_snp( $snp, $samples, $genetics );
        $markers{$inferparentcall}{$origparentcall}{$type}++;
    }

    \%markers;
}

sub parse_snp {
    my ( $snp, $samples, $genetics ) = @_;

    my %marker;

    my @f     = split /\t/, $snp;
    my $callf = $f[8];
    my $pos   = $f[1];
    my @parentcalls =
      map { get_cp( 'GT', $f[ $samples->{parents}{lookup}{$_} ], $callf ) }
      @{ $genetics->{parents} };

    my ( $parentcall, $vgt ) = get_parent_call( \@parentcalls );
    return ( "@parentcalls", $parentcall, "Reject", $pos )
      if $parentcall eq "XXX";

    my $type = get_marker_type( $parentcall, $vgt, \@f, $samples, $genetics );
#    print "Type returned: $type\n";
    return ( "@parentcalls", $parentcall, $type, $pos );
}

sub get_marker_type {
    my ( $parentcall, $vgt, $fref, $samples, $genetics ) = @_;

    my $type = "";
    
    my @femalecalls = get_f2_calls( \@{ $genetics->{sex}{"Female"} },
        $parentcall, $vgt, $fref, $samples );
    my @malecalls = get_f2_calls( \@{ $genetics->{sex}{"Male"} },
        $parentcall, $vgt, $fref, $samples );

#    print "$parentcall\nFemales:@femalecalls\nMales:@malecalls\n";
    return "Reject" if ( !defined $genetics->{types}{$parentcall} );

    my %types;
#    print "Types to test:", keys %{$genetics->{types}{$parentcall}}, "\n";
    for my $type (keys %{$genetics->{types}{$parentcall}}) {
#        print "Testing $parentcall for $type\n";
        $types{$type}++ if (get_chisq(\@malecalls, \@femalecalls, $genetics->{types}{$parentcall}{$type}{males}, $genetics->{types}{$parentcall}{$type}{females}));
    }
    return (keys %types == 1) ? (keys %types)[0] : "Reject";
}

sub get_f2_calls {
    my ( $offspring, $parentcall, $vgt, $fref, $samples ) = @_;
    my @calls;

    for my $sample ( @{$offspring} ) {
        my $opos = $samples->{offspring}{lookup}{$sample};
        my $gt = get_cp( 'GT', $fref->[$opos], $fref->[8] );
        push @calls, $vgt->{$gt} // 'X';
    }
    return @calls;
}

sub get_chisq {
    my ($malecalls, $femalecalls, $validmale, $validfemale) = @_;
    
    my $n = @{$malecalls} + @{$femalecalls};
    my %obs;
    map {$obs{"M$_"}++} @{$malecalls};
    map {$obs{"F$_"}++} @{$femalecalls};
    
    my %exp;
    get_expected_classes(\%exp, 'M', $malecalls, $validmale);
    get_expected_classes(\%exp, 'F', $femalecalls, $validfemale);
#    print "\n";
    my $chisq;
    
    for my $e (keys %exp) {
        my $obsv = $obs{$e} // 0;
        $chisq += (($obsv-$exp{$e}) ** 2) / $exp{$e};
    }
    my $df = (keys %exp) - 1;
    
#    print "Chisq=$chisq, df=$df, crit=$chisqcrit{$df}\n";
    return $chisq < $chisqcrit{$df};
}

sub get_expected_classes {
    my ($exp_ref, $sex, $calls, $valid) = @_;
    my @exp = split /,/, $valid;
    my $shares = 0;
    for my $e (@exp) {
        if ($e =~ /^(\d)?([ABH.])$/) {
            my $share = $1 // 1;
            $exp_ref->{"$sex$2"} = $share;
            $shares += $share;
#            print "$sex$2 x $share : "
        }
    }
    for my $e (@exp) {
        $e = $2 if ($e =~ /^(\d)([ABH.])$/);
        $exp_ref->{"$sex$e"} = $exp_ref->{"$sex$e"}/$shares * @{$calls};
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

    return $nullcall{$part} if $call eq './.';

    my @cp = split ':', $call;
    if ( !defined $callorder{$order} ) {
        my @neworder = split ':', $order;
        map { $callorder{$order}{ $neworder[$_] } = $_ } 0 .. $#neworder;
    }

    return $cp[ $callorder{$order}{$part} ];
}

sub merge {
    my ( $part, $all ) = @_;
    foreach my $scf ( keys %{$part} ) {
        foreach my $inferpattern ( keys %{ $part->{$scf} } ) {
            foreach my $origpattern ( keys %{ $part->{$scf}{$inferpattern} } ) {
                foreach my $type (keys %{$part->{$scf}{$inferpattern}{$origpattern}}) {
                $all->{patterns}{$inferpattern}{$origpattern}{$type} +=
                  $part->{$scf}{$inferpattern}{$origpattern}{$type};
              }
            }
        }
    }
}

sub output {
    my ( $data, $samples, $genome, $outfix ) = @_;
    print STDERR "Outputting...\n";

    open my $allout, '>', "$outfix.out"
      or croak "Can't open $outfix.out: $OS_ERROR\n";

    #    foreach my $scf ( sort keys %{ $data->{scf} } ) {
    #        open my $scffile, '<', "$scf.tmp.out"
    #          or croak "Can't open scf output for $scf: $OS_ERROR\n";
    #        while ( my $scfline = <$scffile> ) {
    #            print $allout $scfline;
    #        }
    #        close $scffile;
    #        system("rm $scf.tmp.out");
    #    }

    my %summarypattern;
    my %summarytype;
    foreach my $inferpattern ( sort keys %{ $data->{patterns} } ) {
        foreach
          my $origpattern ( sort keys %{ $data->{patterns}{$inferpattern} } )
        {
            foreach my $type (sort keys %{$data->{patterns}{$inferpattern}{$origpattern}}) {
                print $allout "$inferpattern\t$origpattern\t$type\t$data->{patterns}{$inferpattern}{$origpattern}{$type}\n";
                $summarypattern{$inferpattern} +=
                  $data->{patterns}{$inferpattern}{$origpattern}{$type};
                $summarytype{$inferpattern}{$type} += $data->{patterns}{$inferpattern}{$origpattern}{$type};
            }
        }
    }

    foreach my $pattern (sort keys %summarytype) {
        foreach my $type (sort keys %{$summarytype{$pattern}}) {
            print $allout "$pattern\t$type\t$summarytype{$pattern}{$type}\n";
        }
    }
    foreach my $pattern ( sort keys %summarypattern ) {
        print $allout "$pattern\t$summarypattern{$pattern}\n";
    }
    close $allout;
}

sub output_scf_to_file {
    my ( $scf, $data ) = @_;
    open my $scfhandle, '>', "$scf.tmp.out"
      or croak "Can't open output file for $scf: $OS_ERROR\n";
    output_scf( $scf, $data, $scfhandle );
    close $scfhandle;
}

sub output_scf {
    my ( $scf, $data, $handle ) = @_;
    foreach my $pos ( sort { $a <=> $b } keys %{$data} ) {

        print $handle
"$scf\t$pos\t$data->{$pos}{origparentcall}\t$data->{$pos}{inferparentcall}\n";
    }
    print $handle '-' x 253, "\n";

}

1;
