sub process {
    my ( $scf, $vcf_line, $samples, $data, $genetics ) = @_;
    chomp $vcf_line;

    my @fields = split /\t/, $vcf_line;

    my @info = split /;/, $fields[7];
    map { $data->{$1}{$2}++ if (/^(.+)=(.+)$/); } @info;

    return;
}

sub merge {
    my ( $part, $all ) = @_;

    return if !defined $part;
    
    foreach my $stat (keys %{$part}) {
        foreach my $val (keys %{$part->{$stat}}) {
            $all->{$stat}{$val} += $part->{$stat}{$val};
        }
    }
}

sub output {
    my ( $data, $samples, $genome, $outfix ) = @_;

    foreach my $stat ( keys %{$data} ) {
        open my $stat_file, ">", "$outfix.$stat.txt"
          or croak "Can't open $outfix.$stat.txt\n";

        foreach my $val (
            sort { $a <=> $b }
            keys %{ $data->{$stat} }
          )
        {
            print $stat_file "$val\t$data->{$stat}{$val}\n";
        }
        close $stat_file;
    }
    return;
}

1;    # Needed to return a true value when imported by vcf-parallel
