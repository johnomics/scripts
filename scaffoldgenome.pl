my %callorder;

my %nullcall = ('GT' => './.', 'GQ' => 0);

sub process {
    my ( $scfref, $samples, $data, $genetics ) = @_;

    my $markers = get_markers( $scfref, $samples, $genetics );
    find_edges($markers);
    print_edges($markers);

}

sub get_markers {
    my ( $scfref, $samples, $genetics ) = @_;

    my %markers;
    foreach my $snp ( @{$scfref} ) {
        chomp $snp;
        my ( $marker, $type, $pos ) = parse_snp( $snp, $samples, $genetics );
        next if $marker == 0;
        $markers{$type}{$pos} = $marker;
    }
    \%markers;
}

sub parse_snp {
    my ( $snp, $samples, $genetics ) = @_;

    my %marker;

    my @f = split /\t/, $snp;
    my $parentcall = join ' ',
      map { get_cp( 'GT', $f[ $samples->{parents}{$_} ], $f[8] ) }
      @{ $genetics->{parents} };
    
    return ( 0, 0, 0 ) if !defined $genetics->{types}{$parentcall};

    foreach my $sample (@{$samples->{offspring}{order}}) {
        my $opos = $samples->{offspring}{lookup}{$sample};
        my $gt = get_cp('GT', $f[$opos], $f[8]);
        $marker{$sample} = $genetics->{types}{$parentcall}{$gt} // 0;
    }

    return (\%marker, $genetics->{types}{$parentcall}{type}, $f[1]);
}

sub get_cp {
    my ( $part, $call, $order ) = @_;

    return $nullcall{$part} if $call eq './.';

    my @cp = split ':', $call;
    if ( !defined $callorder{$order} ) {
        my @neworder = split ':', $order;
        map { $callorder{$order}{ $neworder[$_] } = $_ } 0 .. $#neworder;
    }

    return $cp[$callorder{$order}{$part}];
}

sub find_edges {
    my $markers = shift;
}

sub print_edges {
    my $markers = shift;
}

sub merge {
    my ( $part, $all ) = @_;
}

sub output {
    my ( $data, $genome, $outfix ) = @_;
}

1;
