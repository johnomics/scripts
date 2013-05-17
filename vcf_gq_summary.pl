my %callorder;
my %nullcall = ( 'GT' => './.', 'GQ' => 0, 'DP' => 0 );

sub process {
    my ( $scf, $vcf_line, $samples, $data_ref, $genetics ) = @_;
    chomp $vcf_line;

    my @f          = split /\t/, $vcf_line;
    my $callf      = $f[8];

    foreach my $sample (@{$samples->{parents}{order}}) {
        my $opos = $samples->{parents}{lookup}{$sample};
        my $gt = get_cp( 'GT', $f[$opos], $callf );
        my $gq = get_cp( 'GQ', $f[$opos], $callf );
        my $dp = get_cp( 'DP', $f[$opos], $callf );
        $data_ref->{$gq}{$dp}{$gt}++;
    }
    
    return;
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
}

sub output {
    my (
        $data, $samples, $genome, $output_prefix
    ) = @_;

    open my $stat_file, ">", "$output_prefix.txt"
      or croak "Can't open $output_prefix.txt\n";

    print $stat_file "GQ\tDP\tGT\tCount\n";
    foreach my $gq ( sort {$a<=>$b} keys %{$data} ) {
        foreach my $dp ( sort {$a<=>$b} keys %{ $data->{$gq} } ) {
            foreach my $gt ( sort keys %{ $data->{$gq}{$dp} } ) {
                print $stat_file "$gq\t$dp\t$gt\t$data->{$gq}{$dp}{$gt}\n";
            }
        }
    }
    close $stat_file;
    return;
}

1;    # Needed to return a true value when imported by vcf-parallel
