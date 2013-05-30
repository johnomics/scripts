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
        my ( $origparentcall, $inferparentcall ) =
          parse_snp( $snp, $samples, $genetics );
        $markers{$inferparentcall}{$origparentcall}++;
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

    return ( "@parentcalls", $parentcall, $pos );
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
    return ("XXX", 0) if keys %vallele > 2;
    $vallele{"."} = ".";

    # Assign A, B and H (for heterozygote) to parental calls
    my %vgt;
    my $parentcall;
    for my $pc ( @{$parents} ) {
        my ( $i, $j ) = $pc =~ /(.)\/(.)/;
        if ( $i eq $j ) {
            $vgt{$pc} = $vallele{$i};
            $parentcall .= $vallele{$i};
        }
        else {
            $vgt{$pc} = "H";
            $parentcall .= "H";
        }
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
        foreach my $inferpattern (keys %{$part->{$scf}}) {
            foreach my $origpattern (keys %{$part->{$scf}{$inferpattern}}) {
                $all->{patterns}{$inferpattern}{$origpattern} += $part->{$scf}{$inferpattern}{$origpattern};
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
    foreach my $inferpattern (sort keys %{$data->{patterns}}) {
        foreach my $origpattern (sort keys %{$data->{patterns}{$inferpattern}}) {
            print $allout "$inferpattern\t$origpattern\t$data->{patterns}{$inferpattern}{$origpattern}\n";
            $summarypattern{$inferpattern} += $data->{patterns}{$inferpattern}{$origpattern};
        }
    }

    foreach my $pattern (sort keys %summarypattern) {
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
