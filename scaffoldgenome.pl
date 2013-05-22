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
    find_edges( $data->{$scf}, $samples );
    collapse( $data->{$scf}, $samples );
    output_scf_to_file( $scf, $data->{$scf} );
}

sub collapse {
    my ( $scfref, $samples ) = @_;
    foreach my $type ( keys %{$scfref} ) {
        foreach my $pos ( keys %{ $scfref->{$type} } ) {
            my $gt        = "";
            my @edge;
            my $cons      = "";
            my $corrected = "";
            my @gq;
            foreach my $sample ( @{ $samples->{offspring}{order} } ) {
                my $ref = $scfref->{$type}{$pos}{marker}{$sample};
                push @edge, defined $ref->{edge}
                  && $ref->{edge} != 0 ? abs( $ref->{edge} ) : '.';
                $gt .= $ref->{gt} == 1 ? 'A' : $ref->{gt} == -1 ? 'B' : '-';
                $cons .=
                  $ref->{cons} == 1 ? 'A' : $ref->{cons} == -1 ? 'B' : '-';
                $corrected .=
                    $ref->{corrected} == 1  ? 'A'
                  : $ref->{corrected} == -1 ? 'B'
                  :                           '-';
                push @gq, $ref->{gq};
            }
            delete $scfref->{$type}{$pos}{marker};
            $scfref->{$type}{$pos}{marker}{edge}      = \@edge;
            $scfref->{$type}{$pos}{marker}{cons}      = $cons;
            $scfref->{$type}{$pos}{marker}{corrected} = $corrected;
            $scfref->{$type}{$pos}{marker}{gt}        = $gt;
            $scfref->{$type}{$pos}{marker}{gq}        = \@gq;
        }
    }
}

sub get_markers {
    my ( $scfref, $samples, $genetics ) = @_;

    my %markers;
    my %prevpos;
    foreach my $snp ( @{$scfref} ) {
        chomp $snp;
        my ( $marker, $type, $parentcall, $pos, $info ) =
          parse_snp( $snp, $samples, $genetics );
        next if !$marker;
        next if !chisq_bc_ok($marker);
        next if $info->{'FS'} > FS_THRESHOLD;
        next if $info->{'MQ'} < MQ_THRESHOLD;
        check_phase( $marker, $markers{$type}{ $prevpos{$type} }{marker} )
          if ( $prevpos{$type} );
        $markers{$type}{$pos}{marker} = $marker;
        $markers{$type}{$pos}{parent} = $parentcall;
        $markers{$type}{$pos}{mq}     = $info->{'MQ'};
        $markers{$type}{$pos}{fs}     = $info->{'FS'};
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
    my %info;

    my @f          = split /\t/, $snp;
    my $callf      = $f[8];
    my $parentcall = join ' ',
      map { get_cp( 'GT', $f[ $samples->{parents}{lookup}{$_} ], $callf ) }
      @{ $genetics->{parents} };

    return ( 0, 0, 0 ) if !defined $genetics->{types}{$parentcall};

    my @info = split /;/, $f[7];
    map { $info{$1} = $2 if (/^(.+)=(.+)$/); } @info;

    foreach my $sample ( @{ $samples->{offspring}{order} } ) {
        my $opos = $samples->{offspring}{lookup}{$sample};
        my $gt = get_cp( 'GT', $f[$opos], $callf );
        $marker{$sample}{call} = $gt;
        $marker{$sample}{gt} = $genetics->{types}{$parentcall}{$gt} // 0;
        return ( 0, 0, 0 )
          if $marker{$sample}{gt} ==
          0.5;    # Reject SNP if invalid homozygote call is found
        $marker{$sample}{gq} = get_cp( 'GQ', $f[$opos], $callf );
    }

    return ( \%marker, $genetics->{types}{$parentcall}{type},
        $parentcall, $f[1], \%info );
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
            get_sample_blocks( $markers->{$type}, $sample, \@pos, MASKLEN );
            get_block_consensus( $markers->{$type}, $sample, \@pos, MASKLEN );
            clean_blocks( $markers->{$type}, $sample, \@pos );
        }
    }
}

sub clean_blocks {
    my ( $marker, $sample, $pos ) = @_;

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
            && ( ( max( @{ $blocks[$b]{pos} } ) - min( @{ $blocks[$b]{pos} } ) )
                < ERRORBLOCK )
            || ( @{ $blocks[$b]{pos} } <= ERRORSNPS )
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
            $marker->{ $maskcallpos[$_] }{marker}{$sample}{gt} * $mask[$_]

 #            ( $marker->{ $maskcallpos[$_] }{marker}{$sample}{gt} *
 #                  $marker->{ $maskcallpos[$_] }{marker}{$sample}{gq} ) / 100 *
 #              $mask[$_]
        } 0 .. $#mask;
    }
}

sub get_block_consensus {
    my ( $marker, $sample, $pos, $masklen ) = @_;
    my @blockpos;
    my @blocks;
    foreach my $p ( @{$pos} ) {
        if ( defined $marker->{$p}{marker}{$sample}{edge}
            && abs( $marker->{$p}{marker}{$sample}{edge} ) > $masklen )
        {
            # Don't include the edge position in the block to date;
            # calculate the consensus for the edge block on its own,
            # then move on to the following positions
            calculate_consensus( \@blockpos, $marker, $sample )
              if @blockpos >= 1;
            @blockpos = ($p);
            calculate_consensus( \@blockpos, $marker, $sample );
            @blockpos = ();
        }
        else {
            push @blockpos, $p;
        }

    }
    calculate_consensus( \@blockpos, $marker, $sample );
}

sub calculate_consensus {
    my ( $blockpos, $marker, $sample ) = @_;
    my %blockgt;

    map {
        $blockgt{ $marker->{$_}{marker}{$sample}{gt} }++;

        #        $blockgt{ $marker->{$_}{marker}{$sample}{gt} } +=
        #          $marker->{$_}{marker}{$sample}{gq}
    } @{$blockpos};

    my $maxgt = ( sort { $blockgt{$b} <=> $blockgt{$a} } keys %blockgt )[0];
    map { $marker->{$_}{marker}{$sample}{cons} = $maxgt } @{$blockpos};
}

sub merge {
    my ( $part, $all ) = @_;
    foreach my $scf ( keys %{$part} ) {
        $all->{scf}{$scf}++;

#        foreach my $type ( keys %{ $part->{$scf} } ) {
#            foreach my $pos ( keys %{ $part->{$scf}{$type} } ) {
#                $all->{scf}{$scf}{$type}{$pos}{parent} =
#                  $part->{$scf}{$type}{$pos}{parent};
#                $all->{scf}{$scf}{$type}{$pos}{marker}{gt} =
#                  $part->{$scf}{$type}{$pos}{marker}{gt};
#                $all->{scf}{$scf}{$type}{$pos}{marker}{edge} =
#                  $part->{$scf}{$type}{$pos}{marker}{edge};
#                $all->{scf}{$scf}{$type}{$pos}{marker}{cons} =
#                  $part->{$scf}{$type}{$pos}{marker}{cons};
#                foreach my $gq ( @{ $part->{$scf}{$type}{$pos}{marker}{gq} } ) {
#                    push @{ $all->{scf}{$scf}{$type}{$pos}{marker}{gq} }, $gq;
#                }
#            }
#        }
    }
}

sub output {
    my ( $data, $samples, $genome, $outfix ) = @_;
    print STDERR "Outputting...\n";

    open my $allout, '>', "$outfix.out"
      or croak "Can't open $outfix.out: $OS_ERROR\n";
    foreach my $scf ( sort keys %{ $data->{scf} } ) {
        open my $scffile, '<', "$scf.tmp.out"
          or croak "Can't open scf output for $scf: $OS_ERROR\n";
        while ( my $scfline = <$scffile> ) {
            print $allout $scfline;
        }
        close $scffile;
        system("rm $scf.tmp.out");
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
    foreach my $type ( sort keys %{$data} ) {
        foreach my $pos ( sort { $a <=> $b } keys %{ $data->{$type} } ) {

            print $handle
"$scf\t$pos\t$type\t$data->{$type}{$pos}{mq}\t$data->{$type}{$pos}{fs}\t";
            my @gt = split //, $data->{$type}{$pos}{marker}{gt};
            foreach my $i ( 0 .. $#gt ) {
                my $col =
                  ceil( 255 - $data->{$type}{$pos}{marker}{gq}[$i] / 4.2 );
                print $handle fg $col, $gt[$i];
            }
            print $handle "\t";

            my @edge      = @{$data->{$type}{$pos}{marker}{edge}};
            my @cons      = split //, $data->{$type}{$pos}{marker}{cons};
            my @corrected = split //, $data->{$type}{$pos}{marker}{corrected};
            my $cons      = "";
            my $corrected = "";
            foreach my $e ( 0 .. $#edge ) {
                my $col =
                  $edge[$e] ne '.' && $edge[$e] > MASKLEN ? 'red1' : 'black';
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

1;
