use List::Util qw/sum max/;
use POSIX qw/ceil/;
use Term::ExtendedColor qw/:all/;

my %callorder;

my %nullcall = ( 'GT' => './.', 'GQ' => 0 );

my @mask;

my $masklen = 4;

sub process {
    my ( $scf, $scfref, $samples, $data, $genetics ) = @_;
    $data->{$scf} = get_markers( $scfref, $samples, $genetics );
    find_edges( $data->{$scf}, $samples );
    collapse( $data->{$scf}, $samples );
}

sub collapse {
    my ( $scfref, $samples ) = @_;
    foreach my $type ( keys %{$scfref} ) {
        foreach my $pos ( keys %{ $scfref->{$type} } ) {
            my $gt   = "";
            my $edge = "";
            my $cons = "";
            my @gq;
            foreach my $sample ( @{ $samples->{offspring}{order} } ) {
                my $ref = $scfref->{$type}{$pos}{marker}{$sample};
                $edge .= $ref->{edge} != 0 ? abs( $ref->{edge} ) : '.';
                $gt .= $ref->{gt} == 1 ? 'A' : $ref->{gt} == -1 ? 'B' : '-';
                $cons .=
                  $ref->{cons} == 1 ? 'A' : $ref->{cons} == -1 ? 'B' : '-';
                push @gq, $ref->{gq};
            }
            delete $scfref->{$type}{$pos}{marker};
            $scfref->{$type}{$pos}{marker}{edge} = $edge;
            $scfref->{$type}{$pos}{marker}{cons} = $cons;
            $scfref->{$type}{$pos}{marker}{gt}   = $gt;
            $scfref->{$type}{$pos}{marker}{gq}   = \@gq;
        }
    }
}

sub get_markers {
    my ( $scfref, $samples, $genetics ) = @_;

    my %markers;
    my %prevpos;
    foreach my $snp ( @{$scfref} ) {
        chomp $snp;
        my ( $marker, $type, $parentcall, $pos ) =
          parse_snp( $snp, $samples, $genetics );
        next if !$marker;
        next if !chisq_bc_ok($marker);
        check_phase( $marker, $markers{$type}{ $prevpos{$type} }{marker} )
          if ( $prevpos{$type} );
        $markers{$type}{$pos}{marker} = $marker;
        $markers{$type}{$pos}{parent} = $parentcall;
        $prevpos{$type}               = $pos;
    }
    \%markers;
}

sub check_phase {
    my ( $cur, $prev ) = @_;
    my $phasea = 0;
    my $phaseb = 0;
    foreach my $sample ( keys %{$cur} ) {
        $phasea++ if ( $cur->{$sample}{gt} eq $prev->{$sample}{gt} );
        $phaseb++ if ( ( $cur->{$sample}{gt} * -1 ) eq $prev->{$sample}{gt} );
    }
    if ( $phaseb > $phasea ) {
        foreach my $sample ( keys %{$cur} ) {
            $cur->{$sample}{gt} *= -1;
        }
    }
    return;
}

sub chisq_bc_ok {
    my ($marker) = @_;
    my %allele;
    foreach my $sample ( keys %{$marker} ) {
        $allele{ $marker->{$sample}{gt} }++;
    }
    delete $allele{0};
    my @type = keys %allele;
    return 0 if @type != 2;

    my $a = $allele{ $type[0] };
    my $b = $allele{ $type[1] };
    return 0 if $a == 0 or $b == 0;
    my $chisq = ( $a - $b )**2 / ( $a + $b );

    return $chisq < 3.8    # df=1, p <0.05
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
        return ( 0, 0, 0 )
          if $marker{$sample}{gt} ==
          0.5;    # Reject SNP if invalid homozygote call is found
        $marker{$sample}{gq} = get_cp( 'GQ', $f[$opos], $callf );
    }

    return ( \%marker, $genetics->{types}{$parentcall}{type}, $parentcall,
        $f[1] );
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
        push @called, $p if $marker->{$p}{marker}{$sample}{gt} != 0;
    }

    #    print "@pos\n@called\n";
    foreach my $i ( $masklen .. $#called - $masklen ) {
        next if defined $marker->{ $called[$i] }{marker}{$sample}{edge};
        my @maskcallpos = @called[ $i - $masklen .. $i + $masklen ];

        my $edgesum = 0;

        $marker->{ $called[$i] }{marker}{$sample}{edge} = int sum map {

            #            $marker->{ $maskcallpos[$_] }{$sample}{gt} * $mask[$_]
            ( $marker->{ $maskcallpos[$_] }{marker}{$sample}{gt} *
                  $marker->{ $maskcallpos[$_] }{marker}{$sample}{gq} ) / 100 *
              $mask[$_]
        } 0 .. $#mask;

#        print "$sample\t$marker->{$called[$i]}{$sample}{edge} --------------------------------------------------------\n";

    }
}

sub get_block_consensus {
    my ( $marker, $sample, $pos, $masklen ) = @_;
    my @blockpos;
    foreach my $p ( @{$pos} ) {
        push @blockpos, $p;
        if (   @blockpos > 1
            && defined $marker->{$p}{marker}{$sample}{edge}
            && abs( $marker->{$p}{marker}{$sample}{edge} ) > $masklen )
        {
            calculate_consensus( \@blockpos, $marker, $sample );
            @blockpos = ();
        }
    }
    calculate_consensus( \@blockpos, $marker, $sample );
}

sub calculate_consensus {
    my ( $blockpos, $marker, $sample ) = @_;
    my %blockgt;

    map {
        $blockgt{ $marker->{$_}{marker}{$sample}{gt} } +=
          $marker->{$_}{marker}{$sample}{gq}
    } @{$blockpos};
    my $maxgt = ( sort { $blockgt{$b} <=> $blockgt{$a} } keys %blockgt )[0];
    map { $marker->{$_}{marker}{$sample}{cons} = $maxgt } @{$blockpos};
    return;
}

sub merge {
    my ( $part, $all ) = @_;
    foreach my $scf ( keys %{$part} ) {
        foreach my $type ( keys %{ $part->{$scf} } ) {
            foreach my $pos ( keys %{ $part->{$scf}{$type} } ) {
                $all->{scf}{$scf}{$type}{$pos}{parent} =
                  $part->{$scf}{$type}{$pos}{parent};
                $all->{scf}{$scf}{$type}{$pos}{marker}{gt} =
                  $part->{$scf}{$type}{$pos}{marker}{gt};
                $all->{scf}{$scf}{$type}{$pos}{marker}{edge} =
                  $part->{$scf}{$type}{$pos}{marker}{edge};
                $all->{scf}{$scf}{$type}{$pos}{marker}{cons} =
                  $part->{$scf}{$type}{$pos}{marker}{cons};
                foreach my $gq ( @{ $part->{$scf}{$type}{$pos}{marker}{gq} } ) {
                    push @{ $all->{scf}{$scf}{$type}{$pos}{marker}{gq} }, $gq;
                }
            }
        }
    }
}

sub output {
    my ( $data, $samples, $genome, $outfix ) = @_;
    print STDERR "Outputting...\n";

    foreach my $scf ( sort keys %{$data->{scf}} ) {
        foreach my $type ( sort keys %{ $data->{scf}{$scf} } ) {
            foreach
              my $pos ( sort { $a <=> $b } keys %{ $data->{scf}{$scf}{$type} } )
            {

                print "$scf\t$pos\t$type\t";
                my @gt = split //, $data->{scf}{$scf}{$type}{$pos}{marker}{gt};
                foreach my $i ( 0 .. $#gt ) {
                    my $col = ceil( 255 -
                          $data->{scf}{$scf}{$type}{$pos}{marker}{gq}[$i] / 4.2 );
                    print fg $col, $gt[$i];
                }
                print "\t";

                my @edge = split //, $data->{scf}{$scf}{$type}{$pos}{marker}{edge};
                my @cons = split //, $data->{scf}{$scf}{$type}{$pos}{marker}{cons};
                my $cons = "";
                foreach my $e (0..$#edge) {
                    my $col = $edge[$e] > $masklen ? 'red1' : 'black';
                    print fg $col, $edge[$e];
                    $cons .= fg $col, $cons[$e];
                }
                print "\t$cons\n";
            }
            print '-' x 253, "\n";
        }
        print '-' x 253, "\n";
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
