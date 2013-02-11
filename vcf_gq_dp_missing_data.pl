sub process {
    my ( $vcf_line, $sample_names_ref, $data_ref ) = @_;
    chomp $vcf_line;

    my $dp_max = 99;
    my $gq_max = 99;
    my @fields = split /\t/, $vcf_line;

    my %position;
    my @format_fields = split /:/, $fields[8];
    my $gt_field_num  = -1;
    my $gq_field_num  = -1;
    my $dp_field_num  = -1;
    for my $i ( 0 .. ( @format_fields - 1 ) ) {
        if ( $format_fields[$i] eq "GT" ) { $gt_field_num = $i; }
        if ( $format_fields[$i] eq "GQ" ) { $gq_field_num = $i; }
        if ( $format_fields[$i] eq "DP" ) { $dp_field_num = $i; }
    }

    my $missing = 0;
    my %pos_dp_filter;
    my %pos_gq_filter;
    for my $i ( 9 .. ( @{$sample_names_ref} + 8 ) ) {
        if ( $fields[$i] eq "./." ) {

            $data_ref->{stat}{dp}{ $sample_names_ref->[ $i - 9 ] }{0}++;
            $data_ref->{stat}{gq}{ $sample_names_ref->[ $i - 9 ] }{0}++;
        }
        else {
            my @sample_f = split /:/, $fields[$i];

            # At threshold $_, increment the number of individuals present for this position;
            # individual will be filtered at every threshold above its own depth/quality
            map { $pos_dp_filter{$_}++ }
              ( 0 .. ($sample_f[$dp_field_num] < $dp_max ? $sample_f[$dp_field_num] - 1 : $dp_max) );
            map { $pos_gq_filter{$_}++ }
              ( 0 .. ($sample_f[$gq_field_num] < $gq_max ? $sample_f[$gq_field_num] - 1 : $gq_max) );

            $data_ref->{stat}{dp}{ $sample_names_ref->[ $i - 9 ] }
              { $sample_f[$dp_field_num] }++;
            $data_ref->{stat}{gq}{ $sample_names_ref->[ $i - 9 ] }
              { $sample_f[$gq_field_num] }++;
            $data_ref->{range}{dp}{ $sample_f[$dp_field_num] }++;
            $data_ref->{range}{gq}{ $sample_f[$gq_field_num] }++;
        }
    }

    # At threshold $_, $pos_xx_filter{$_} individuals would be present at this position, which is on scaffold $fields[0]
    map {
        $data_ref->{present}{dp}{$_}{ $pos_dp_filter{$_} }{ $fields[0] }++
    } keys %pos_dp_filter;
    map {
        $data_ref->{present}{gq}{$_}{ $pos_gq_filter{$_} }{ $fields[0] }++
    } keys %pos_gq_filter;
    return;
}

sub merge {
    my ($data_ref, $final_ref) = @_;

    if ( defined($data_ref) ) {
        foreach my $stat ( keys %{ $data_ref->{present} } ) {
            foreach my $threshold ( keys %{ $data_ref->{present}{$stat} } )
            {
                foreach my $missing_ind (
                    keys %{ $data_ref->{present}{$stat}{$threshold} } )
                {

                    foreach my $scf (
                        keys %{
                            $data_ref->{present}{$stat}{$threshold}
                              {$missing_ind}
                        }
                      )
                    {
                        $final_ref->{present}{$stat}{$threshold}
                          {$missing_ind}{$scf} +=
                          $data_ref->{present}{$stat}{$threshold}
                          {$missing_ind}{$scf};
                    }
                }
            }
        }
        foreach my $stat ( keys %{ $data_ref->{stat} } ) {
            foreach my $sample_name ( keys %{ $data_ref->{stat}{$stat} } ) {
                foreach my $val (
                    keys %{ $data_ref->{stat}{$stat}{$sample_name} } )
                {
                    $final_ref->{stat}{$stat}{$sample_name}{$val} +=
                      $data_ref->{stat}{$stat}{$sample_name}{$val};
                }
            }
        }
        foreach my $stat ( keys %{ $data_ref->{range} } ) {
            foreach my $val ( keys %{ $data_ref->{range}{$stat} } ) {
                $final_ref->{range}{$stat}{$val} +=
                  $data_ref->{range}{$stat}{$val};
            }
        }
    }
}


sub output {
    my ( $final_data_ref, $scflen_ref, $sample_names_ref, $genome_length, $output_prefix ) = @_;
    foreach my $stat ( "dp", "gq" ) {
        open my $stat_file, ">", "$output_prefix.sample.$stat.txt" or croak "Can't open $output_prefix.sample.$stat.txt\n";
        
        print $stat_file "Ind\\$stat";
        foreach my $stat_val (
            sort { $a <=> $b }
            keys %{ $final_data_ref->{range}{$stat} }
          )
        {
            print $stat_file "\t$stat_val";
        }
        print $stat_file "\n";

        foreach my $sample (
            sort { $a <=> $b }
            keys %{ $final_data_ref->{stat}{$stat} }
          )
        {
            print $stat_file "$sample";
            foreach my $stat_val (
                sort { $a <=> $b }
                keys %{ $final_data_ref->{range}{$stat} }
              )
            {
                if (
                    defined $final_data_ref->{stat}{$stat}{$sample}{$stat_val} )
                {
                    print $stat_file
                      "\t$final_data_ref->{stat}{$stat}{$sample}{$stat_val}";
                }
                else {
                    print $stat_file "\t0";
                }
            }
            print $stat_file "\n";
        }
        close $stat_file;
    }

    open my $present_file, ">", "$output_prefix.present.txt" or croak "Can't open $output_prefix.present.txt\n";
    print $present_file "Stat\tThreshold\tSamples\tBases\tCumulative.Bases\tCumulative.Scaffolds\tCumulative.Scaffold.pc\tGenome.Coverage\tGenome.Coverage.pc\n";
    foreach my $stat ( sort keys %{ $final_data_ref->{present} } ) {

        foreach my $threshold (
            sort { $a <=> $b }
            keys %{ $final_data_ref->{present}{$stat} }
          )
        {
            my $cumul_bases = 0;
            my %cumul_scf;

            foreach my $present_ind ( reverse 1 .. @{$sample_names_ref} ) {

                my $bases_w_ind_present = 0;
                if (
                    defined $final_data_ref->{present}{$stat}{$threshold}
                    {$present_ind} )
                {
                    foreach my $scf (
                        keys %{
                            $final_data_ref->{present}{$stat}{$threshold}
                              {$present_ind}
                        }
                      )
                    {
                        $bases_w_ind_present +=
                          $final_data_ref->{present}{$stat}{$threshold}
                          {$present_ind}{$scf};
                        $cumul_bases +=
                          $final_data_ref->{present}{$stat}{$threshold}
                          {$present_ind}{$scf};
                        $cumul_scf{$scf}++;
                    }
                }

                my $base_coverage = 0;
                foreach my $scf ( keys %cumul_scf ) {
                    $base_coverage += $scflen_ref->{$scf};
                }
                my $pc_base_coverage = sprintf "%5.2f",
                  $base_coverage / $genome_length * 100;
                my $cumul_scf_num = keys %cumul_scf;
                my $pc_scfs       = sprintf "%5.2f",
                  $cumul_scf_num / keys( %{$scflen_ref} ) * 100;
                print $present_file
"$stat\t$threshold\t$present_ind\t$bases_w_ind_present\t$cumul_bases\t$cumul_scf_num\t$pc_scfs\t$base_coverage\t$pc_base_coverage\n";
            }
        }
    }
    close $present_file;
    return;
}

1; # Needed to return a true value when imported by vcf-parallel