my %callorder;
my %nullcall = ( 'GT' => './.', 'GQ' => 0, 'DP' => 0 );

our $setup = sub {
    my ( $vcf_filename, $extraargs ) = @_;

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
        $samples{lookup}{ $sample_names[$_] } = $_;
        push @{ $samples{order} }, $sample_names[$_];
    } 9 .. $#sample_names;

    \%samples;
};

our $process = sub {
    my ( $scf, $vcf_line, $data_ref, $samples ) = @_;
    chomp $vcf_line;

    my @f = split /\t/, $vcf_line;
    my $callf = $f[8];

    foreach my $sample ( @{ $samples->{order} } ) {
        my $opos = $samples->{lookup}{$sample};
        my $gt   = get_cp( 'GT', $f[$opos], $callf );
        my $gq   = get_cp( 'GQ', $f[$opos], $callf );
        my $dp   = get_cp( 'DP', $f[$opos], $callf );
        $data_ref->{$gq}{$dp}{$gt}++;
    }

    return;
};

sub get_cp {
    my ( $part, $call, $order ) = @_;

    return $nullcall{$part} if $call =~ '^\./\.';

    my @cp = split ':', $call;
    if ( !defined $callorder{$order} ) {
        my @neworder = split ':', $order;
        map { $callorder{$order}{ $neworder[$_] } = $_ } 0 .. $#neworder;
    }

    my $callpart = $cp[ $callorder{$order}{$part} ];
    return $callpart eq '.' ? 0 : $callpart;
}

our $merge = sub {
    my ( $data_ref, $final_ref ) = @_;

    if ( defined($data_ref) ) {
        foreach my $gq ( keys %{$data_ref} ) {
            foreach my $dp ( keys %{ $data_ref->{$gq} } ) {
                foreach my $gt ( keys %{ $data_ref->{$gq}{$dp} } ) {
                    $final_ref->{$gq}{$dp}{$gt} += $data_ref->{$gq}{$dp}{$gt};
                }
            }
        }
    }
    return;
};

our $output = sub {
    my ( $data, $genome, $output_prefix ) = @_;

    open my $stat_file, ">", "$output_prefix.txt"
      or croak "Can't open $output_prefix.txt\n";

    print $stat_file "GQ\tDP\tGT\tCount\n";
    foreach my $gq ( sort { $a <=> $b } keys %{$data} ) {
        foreach my $dp ( sort { $a <=> $b } keys %{ $data->{$gq} } ) {
            foreach my $gt ( sort keys %{ $data->{$gq}{$dp} } ) {
                print $stat_file "$gq\t$dp\t$gt\t$data->{$gq}{$dp}{$gt}\n";
            }
        }
    }
    close $stat_file;
    return;
};

1;    # Needed to return a true value when imported by vcf-parallel
