use List::Util qw/sum/;
use Storable qw/dclone/;
my %callorder;

my %nullcall = ( 'GT' => './.', 'GQ' => 0 );

my @mask;

my $masklen = 3;

sub process {
    my ( $scf, $scfref, $samples, $data, $genetics ) = @_;

    $data->{$scf} = get_markers( $scfref, $samples, $genetics );
    find_edges( $data->{$scf}, $samples );
}

sub get_markers {
    my ( $scfref, $samples, $genetics ) = @_;

    my %markers;
    foreach my $snp ( @{$scfref} ) {
        chomp $snp;
        my ( $marker, $type, $pos ) = parse_snp( $snp, $samples, $genetics );

        $markers{$type}{$pos} = chisq_bc_ok($marker) ? $marker : "Fail";
    }
    \%markers;
}

sub chisq_bc_ok {
    my ($marker) = @_;
    my %allele;
    foreach my $sample (keys %{$marker}) {
        $allele{$marker->{$sample}{gt}}++;
    }
    delete $allele{0};
    my @type = keys %allele;
    return 0 if @type != 2;
    
    my $a = $allele{$type[0]};
    my $b = $allele{$type[1]};
    return 0 if $a == 0 or $b == 0;
    my $chisq = ($a - $b) ** 2 / ($a + $b);
    
    return $chisq < 3.8 # df=1, p <0.05
}

sub parse_snp {
    my ( $snp, $samples, $genetics ) = @_;

    my %marker;

    my @f          = split /\t/, $snp;
    my $callf      = $f[8];
    my $parentcall = join ' ',
      map { get_cp( 'GT', $f[ $samples->{parents}{$_} ], $callf ) }
      @{ $genetics->{parents} };

    return ( 0, 0, 0 ) if !defined $genetics->{types}{$parentcall};

    foreach my $sample ( @{ $samples->{offspring}{order} } ) {
        my $opos = $samples->{offspring}{lookup}{$sample};
        my $gt = get_cp( 'GT', $f[$opos], $callf );
        $marker{$sample}{gt} = $genetics->{types}{$parentcall}{$gt} // 0;
        $marker{$sample}{gq} = get_cp( 'GQ', $f[$opos], $callf );
    }

    return ( \%marker, $genetics->{types}{$parentcall}{type}, $f[1] );
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

sub find_edges {
    my $markers = shift;
    my $samples = shift;

    foreach my $type ( keys %{$markers} ) {
        my @pos = sort { $a <=> $b } keys %{ $markers->{$type} };
        foreach my $sample ( @{ $samples->{offspring}{order} } ) {
            get_sample_blocks( $markers->{$type}, $sample, \@pos, $masklen );
            get_block_consensus( $markers->{$type}, $sample, \@pos, $masklen );
        }
    }
}

sub get_sample_blocks {
    my ( $marker, $sample, $pos, $masklen, $mask ) = @_;

    @mask = ( (-1) x $masklen, 0, (1) x $masklen ) if ( !@mask );
    my @called;
    foreach my $p ( @{$pos} ) {
        push @called, $p if $marker->{$p}{$sample}{gt} != 0;
    }
    foreach my $i ( $masklen .. $#called - $masklen ) {
        next if defined $marker->{$called[$i]}{$sample}{edge};
        my @maskcallpos = @called[ $i - $masklen .. $i + $masklen ];
        $marker->{ $called[$i] }{$sample}{edge} = sum map {
            $marker->{ $maskcallpos[$_] }{$sample}{gt} *
              $marker->{ $maskcallpos[$_] }{$sample}{gq} / 100 *
              $mask[$_]
        } 0 .. $#mask;
    }
}

sub get_block_consensus {
    my ( $marker, $sample, $pos, $masklen ) = @_;
    my @blockpos;
    foreach my $p (@{$pos}) {
        push @blockpos, $p;
        if (defined $marker->{$p}{$sample}{edge} && $marker->{$p}{$sample}{edge} >= $masklen) {
            calculate_consensus(\@blockpos, $marker, $sample);
            @blockpos = ();
        }
    }
    calculate_consensus(\@blockpos, $marker, $sample);
}

sub calculate_consensus {
    my ($blockpos, $marker, $sample) = @_;
    my %blockgt;
    
    map {$blockgt{$marker->{$_}{$sample}{gt}} += $marker->{$_}{$sample}{gq}} @{$blockpos};
    my $maxgt = (sort {$blockgt{$b} <=> $blockgt{$a}} keys %blockgt)[0];
    map {$marker->{$_}{$sample}{cons} = $maxgt} @{$blockpos};
    return;
}

sub merge {
    my ( $part, $all ) = @_;
    foreach my $scf ( keys %{$part} ) {
        $all->{$scf} = dclone( \%{ $part->{$scf} } );
    }
}

sub output {
    my ( $data, $samples, $genome, $outfix ) = @_;

    foreach my $scf ( sort keys %{$data} ) {
        foreach my $type ( sort keys %{ $data->{$scf} } ) {
            print "$scf\t$type\n";
            foreach
              my $pos ( sort { $a <=> $b } keys %{ $data->{$scf}{$type} } )
            {
                print "$scf\t$pos\t";
                print_pattern( $data->{$scf}{$type}{$pos},
                    "gt", $samples->{offspring}{order} );
                print "\t";
                foreach my $sample ( @{ $samples->{offspring}{order} } ) {
                    print defined $data->{$scf}{$type}{$pos}{$sample}{edge}
                      && abs( $data->{$scf}{$type}{$pos}{$sample}{edge} ) >= $masklen
                      ? '-'
                      : '.';
                }
                print "\t";
                print_pattern( $data->{$scf}{$type}{$pos},
                    "cons", $samples->{offspring}{order} );
                print "\n";
            }
        }
    }
}

sub print_pattern {
    my ( $pos, $field, $samples ) = @_;
    foreach my $sample ( @{$samples} ) {
        print $pos->{$sample}{$field} == 1 ? 'A'
          : $pos->{$sample}{$field} == -1  ? 'B'
          :                                  '-';
    }
    return;
}

1;
