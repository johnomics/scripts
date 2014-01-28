use Memoize;

memoize 'mat_inter_match';

my %mstheader = (
    population_type              => "DH",
    population_name              => "HeliconiusWGS",
    distance_function            => "kosambi",
    cut_off_p_value              => "0.000001",
    no_map_dist                  => "0",
    no_map_size                  => "0",
    missing_threshold            => "1",
    estimation_before_clustering => "yes",
    detect_bad_data              => "yes",
    objective_function           => "ML",
);

my @mstheader = (
    "population_type",   "population_name",
    "distance_function", "cut_off_p_value",
    "no_map_dist",       "no_map_size",
    "missing_threshold", "estimation_before_clustering",
    "detect_bad_data",   "objective_function",
    "number_of_loci",    "number_of_individual",
);

my %ignore = (
    "90"  => 0,
    "1N"  => 0,
    "108" => 0,
);

my %calls = (
    "0/1 0/1 0/0" =>
      { type => "maternal", "0/0" => "A", "0/1" => "B", "1/1" => "-" },
    "0/0 0/1 0/0" =>
      { type => "maternal", "0/0" => "B", "0/1" => "A", "1/1" => "-" },
    "0/1 0/1 1/1" =>
      { type => "maternal", "0/0" => "-", "0/1" => "B", "1/1" => "A" },
    "1/1 0/1 1/1" =>
      { type => "maternal", "0/0" => "-", "0/1" => "A", "1/1" => "B" },
    "0/1 0/0 0/1" =>
      { type => "paternal", "0/0" => "A", "0/1" => "B", "1/1" => "-" },
    "0/0 0/0 0/1" =>
      { type => "paternal", "0/0" => "B", "0/1" => "A", "1/1" => "-" },
    "0/1 1/1 0/1" =>
      { type => "paternal", "0/0" => "-", "0/1" => "B", "1/1" => "A" },
    "1/1 1/1 0/1" =>
      { type => "paternal", "0/0" => "-", "0/1" => "A", "1/1" => "B" },
    "0/1 0/1 0/1" =>
      { type => "intercross", "0/0" => "A", "0/1" => "X", "1/1" => "B" },
    "0/0 0/1 0/1" =>
      { type => "intercross", "0/0" => "A", "0/1" => "X", "1/1" => "B" },
    "1/1 0/1 0/1" =>
      { type => "intercross", "0/0" => "A", "0/1" => "X", "1/1" => "B" },
);

my %types = ();
map { $types{ $calls{$_}{type} }++ } keys %calls;

sub process {
    my ( $scf_vcf_ref, $sample_names_ref, $data_ref, $parent_string ) = @_;

    my $bin_size = 5000;
    my $bin_max  = $bin_size;
    my %bin;
    my $scaffold;
    foreach my $vcf_line ( @{$scf_vcf_ref} ) {
        chomp $vcf_line;

        my @fields = split /\t/, $vcf_line;
        $scaffold = $fields[0];
        my $pos = $fields[1];
        if ( $pos < $bin_max ) {
            my @format_fields = split /:/, $fields[8];
            my %field_num;
            for my $i ( 0 .. ( @format_fields - 1 ) ) {
                $field_num{ $format_fields[$i] } = $i;
            }

            my %parents;
            my @parents;
            for my $i (
                ( @{$sample_names_ref} + 6 ) .. ( @{$sample_names_ref} + 8 ) )
            {
                push @parents, $sample_names_ref->[ $i - 9 ];
                if ( $fields[$i] eq "./." ) {
                    $parents{ $sample_names_ref->[ $i - 9 ] } = "-";
                }
                else {
                    my @sample_f = split /:/, $fields[$i];
                    $parents{ $sample_names_ref->[ $i - 9 ] } =
                      $sample_f[ $field_num{GT} ];
                }
            }
            if ($parent_string) { @parents = split / /, $parent_string }
            my $parent_call = "";
            foreach my $parent (@parents) {
                $parent_call .= "$parents{$parent} ";
            }
            chop $parent_call;

            next if ( !defined $calls{$parent_call} );

            for my $i ( 9 .. ( @{$sample_names_ref} + 5 ) ) {
                my $offspring = $sample_names_ref->[ $i - 9 ];
                next if ( defined $ignore{$offspring} );

                if ( $fields[$i] eq "./." ) {
                    $data_ref->{snps}{ $calls{$parent_call}{type} }{$scaffold}
                      {$pos}{pattern} .= "-";
                    $bin{ $calls{$parent_call}{type} }
                      { $sample_names_ref->[ $i - 9 ] }{"-"}{$pos} = 0;
                }
                else {
                    my @sample_f = split /:/, $fields[$i];
                    $data_ref->{snps}{ $calls{$parent_call}{type} }{$scaffold}
                      {$pos}{pattern} .=
                      $calls{$parent_call}{ $sample_f[ $field_num{GT} ] };
                    $bin{ $calls{$parent_call}{type} }
                      { $sample_names_ref->[ $i - 9 ] }
                      { $calls{$parent_call}{ $sample_f[ $field_num{GT} ] } }
                      {$pos} = $sample_f[ $field_num{GQ} ];
                }
            }
        }
        else {
            process_bin(
                $scaffold,         $bin_max - $bin_size, \%bin,
                $sample_names_ref, $data_ref,            $parent_string
            );
            %bin = ();
            $bin_max += $bin_size;
        }
    }

# Process end of scaffold; this may not be a good idea, because the region covered
# may be very small and only a few SNPs are present
#    if ( keys %bin > 0 ) {
#        process_bin(
#            $scaffold,         $bin_max - $bin_size, \%bin,
#            $sample_names_ref, $data_ref,            $parent_string
#        );
#    }

    return;
}

sub process_bin {
    my (
        $scaffold,         $bin_pos,  $bin_ref,
        $sample_names_ref, $data_ref, $parent_string
    ) = @_;

    foreach my $type ( keys %types ) {
        if ( !defined $bin_ref->{$type} ) {
            for my $i ( 0 .. @{$sample_names_ref} - 4 ) {
                next if ( defined $ignore{ $sample_names_ref->[$i] } );
                $data_ref->{bins}{$type}{$scaffold}{$bin_pos}{pattern} .= "-";
                $data_ref->{bins}{$type}{$scaffold}{$bin_pos}{quals}   .= "\t0";
            }
            $data_ref->{bins}{$type}{$scaffold}{$bin_pos}{count} = 0;
            next;
        }
        my $snps = 0;
        for my $i ( 0 .. @{$sample_names_ref} - 4 ) {
            $snps = 0;
            my $sample = $sample_names_ref->[$i];
            next if ( defined $ignore{$sample} );
            my %sum_gt;
            foreach my $gt ( "A", "B", "X", "-" ) {
                foreach my $pos ( keys %{ $bin_ref->{$type}{$sample}{$gt} } ) {
                    $sum_gt{$gt} += $bin_ref->{$type}{$sample}{$gt}{$pos};
                    $snps++;
                }
            }
            my @gts = reverse sort { $sum_gt{$a} <=> $sum_gt{$b} } keys %sum_gt;
            my $max_gt  = $gts[0];
            my $max_gtq = $sum_gt{ $gts[0] } - $sum_gt{ $gts[1] };
            $data_ref->{bins}{$type}{$scaffold}{$bin_pos}{pattern} .= $max_gt;
            $data_ref->{bins}{$type}{$scaffold}{$bin_pos}{quals} .=
              "\t$max_gtq";
        }
        $data_ref->{bins}{$type}{$scaffold}{$bin_pos}{count} = $snps;

        my @gqs = split /\t/,
          $data_ref->{bins}{$type}{$scaffold}{$bin_pos}{quals};
        shift @gqs;
        my $bad_pos = 0;
        map { $bad_pos = 1 if ( $_ < 99 ); } @gqs;
        if ($bad_pos) {
            $data_ref->{bins}{$type}{$scaffold}{$bin_pos}{bad}++;
            next;
        }

        $data_ref->{patterns}{$type}
          { $data_ref->{bins}{$type}{$scaffold}{$bin_pos}{pattern} }{count}++;
        $data_ref->{patterns}{$type}
          { $data_ref->{bins}{$type}{$scaffold}{$bin_pos}{pattern} }{bins}
          {$scaffold}{$bin_pos}++;
    }
    return;
}

sub merge {
    my ( $data_ref, $final_ref ) = @_;

    foreach my $type ( keys %{ $data_ref->{snps} } ) {
        foreach my $scf ( keys %{ $data_ref->{snps}{$type} } ) {
            foreach my $pos ( keys %{ $data_ref->{snps}{$type}{$scf} } ) {
                $final_ref->{snps}{$type}{$scf}{$pos}{pattern} =
                  $data_ref->{snps}{$type}{$scf}{$pos}{pattern};
            }
        }
    }

    foreach my $type ( keys %{ $data_ref->{bins} } ) {
        foreach my $scf ( keys %{ $data_ref->{bins}{$type} } ) {
            foreach my $pos ( keys %{ $data_ref->{bins}{$type}{$scf} } ) {
                $final_ref->{bins}{$type}{$scf}{$pos}{pattern} =
                  $data_ref->{bins}{$type}{$scf}{$pos}{pattern};
                $final_ref->{bins}{$type}{$scf}{$pos}{quals} =
                  $data_ref->{bins}{$type}{$scf}{$pos}{quals};
                $final_ref->{bins}{$type}{$scf}{$pos}{count} =
                  $data_ref->{bins}{$type}{$scf}{$pos}{count};
                if ( defined $data_ref->{bins}{$type}{$scf}{$pos}{bad} ) {
                    $final_ref->{bins}{$type}{$scf}{$pos}{bad}++;
                }
            }
        }
    }
    foreach my $type ( keys %{ $data_ref->{patterns} } ) {
        foreach my $pattern ( keys %{ $data_ref->{patterns}{$type} } ) {
            $final_ref->{patterns}{$type}{$pattern}{count} +=
              $data_ref->{patterns}{$type}{$pattern}{count};
            foreach my $scaffold (
                keys %{ $data_ref->{patterns}{$type}{$pattern}{bins} } )
            {
                foreach my $bin_pos (
                    keys
                    %{ $data_ref->{patterns}{$type}{$pattern}{bins}{$scaffold} }
                  )
                {
                    $final_ref->{patterns}{$type}{$pattern}{bins}{$scaffold}
                      {$bin_pos} +=
                      $data_ref->{patterns}{$type}{$pattern}{bins}{$scaffold}
                      {$bin_pos};
                }
            }
        }
    }
    return;
}

sub output {
    my (
        $final_ref,     $scflen_ref, $sample_names_ref,
        $genome_length, $output_prefix
    ) = @_;

    foreach my $type ( keys %{ $final_ref->{bins} } ) {
        open my $snps, ">", "$output_prefix.$type.snps.txt"
          or croak
          "Can't open SNPs file $output_prefix.$type.snps.txt! $OS_ERROR\n";
        foreach my $scf ( sort keys %{ $final_ref->{snps}{$type} } ) {
            foreach my $pos (
                sort { $a <=> $b }
                keys %{ $final_ref->{snps}{$type}{$scf} }
              )
            {
                print $snps
"$scf\t$pos\t$final_ref->{snps}{$type}{$scf}{$pos}{pattern}\n";
            }
        }
        close $snps;
        open my $output, ">", "$output_prefix.$type.bins.txt"
          or croak
          "Can't open output file $output_prefix.$type.bins.txt! $OS_ERROR\n";
        foreach my $scf ( sort keys %{ $final_ref->{bins}{$type} } ) {
          POS:
            foreach my $pos (
                sort { $a <=> $b }
                keys %{ $final_ref->{bins}{$type}{$scf} }
              )
            {
                my $pattern = $final_ref->{bins}{$type}{$scf}{$pos}{pattern};
                my $valid =
                    ( $pattern =~ "-" ) ? "Missing"
                  : ( $final_ref->{patterns}{$type}{$pattern}{count} == 1 )
                  ? "Single"
                  : defined $final_ref->{bins}{$type}{$scf}{$pos}{bad} ? "LowGQ"
                  :   "Valid";

                print $output
"$scf\t$pos\t$valid\t$final_ref->{bins}{$type}{$scf}{$pos}{count}\t$final_ref->{bins}{$type}{$scf}{$pos}{pattern}$final_ref->{bins}{$type}{$scf}{$pos}{quals}\n";
            }
        }
    }
    close $output;

    open my $markerpos, ">", "$output_prefix.markerpos"
      or croak "Can't open $output_prefix.markerpos: $OS_ERROR\n";

    my $marker_out;
    $marker_out .= "locus_name";
    my $num_samples = 0;
    map {
        if ( !defined $ignore{ $sample_names_ref->[$_] } ) {
            $marker_out .= "\t$sample_names_ref->[$_]";
            $num_samples++;
        }
    } 0 .. $#{$sample_names_ref} - 3;
    $marker_out .= "\n";

    # Write out paternal markers
    my $marker_num = 0;
    foreach my $pattern (
        sort {
            $final_ref->{patterns}{"paternal"}{$a}{count}
              <=> $final_ref->{patterns}{"paternal"}{$b}{count}
        } keys %{ $final_ref->{patterns}{"paternal"} }
      )
    {
        next
          if ( $pattern =~ "-" );
        next if ( $final_ref->{patterns}{"paternal"}{$pattern}{count} <= 1 );

        ( $marker_num, $marker_out, $markerpos ) =
          output_marker( $marker_num, $marker_out, $pattern, $markerpos,
            $final_ref->{patterns}{paternal}{$pattern} );
    }

    # Convert and output intercross markers

    my %mat_int_pat = (
        "A" => { "A" => "A", "X" => "B" },
        "B" => { "B" => "B", "X" => "A" }
    );

# Convert intercross markers to paternal markers where they match maternal markers
    my $intercross_count = keys %{ $final_ref->{patterns}{"intercross"} };
    print STDERR "$intercross_count intercross markers to process...\n";
    my $count        = 0;
    my $notsingleton = 0;
    my $nomatch      = 0;
    my $multiplemat  = 0;
    my $output       = 0;
    foreach my $intercross ( keys %{ $final_ref->{patterns}{"intercross"} } ) {
        $count++;
        next
          if ( $final_ref->{patterns}{"intercross"}{$intercross}{count} <= 1 );
        $notsingleton++;

        if ( $count % 100 == 0 ) {
            print STDERR
"$count intercross markers processed, $notsingleton not singletons, $multiplemat have multiple maternal patterns, $nomatch have no match, $output converted\n";
        }
        my @maternal = ();
        foreach my $maternal ( keys %{ $final_ref->{patterns}{"maternal"} } ) {
            next
              if ( $final_ref->{patterns}{"maternal"}{$maternal}{count} <= 1 );
            if ( mat_inter_match( $intercross, $maternal ) ) {
                push @maternal, $maternal;
            }
        }
        if ( @maternal != 1 ) {
            if ( @maternal == 0 ) {
                $nomatch++;
                my $intercross_missing = $intercross =~ s/X/-/g;
                ( $marker_num, $marker_out, $markerpos ) = output_marker(
                    $marker_num,
                    $marker_out,
                    $intercross_missing,
                    $markerpos,
                    $final_ref->{patterns}{intercross}{$intercross},
                    $intercross
                );
            }
            else {
                $multiplemat++;
            }
            next;
        }
        my $paternal = "";
        my @matgt    = split //, $maternal[0];
        my @intgt    = split //, $intercross;
        for my $i ( 0 .. $#matgt ) {
            $paternal .= $mat_int_pat{ $matgt[$i] }{ $intgt[$i] };
        }
        $output++;
        ( $marker_num, $marker_out, $markerpos ) =
          output_marker( $marker_num, $marker_out, $paternal, $markerpos,
            $final_ref->{patterns}{intercross}{$intercross}, $intercross );
    }

    close $markerpos;

    my %mstheader = (
        population_type              => "DH",
        population_name              => "HeliconiusWGS",
        distance_function            => "kosambi",
        cut_off_p_value              => "0.000001",
        no_map_dist                  => "15",
        no_map_size                  => "2",
        missing_threshold            => "0.25",
        estimation_before_clustering => "yes",
        detect_bad_data              => "yes",
        objective_function           => "ML",
    );

#    print
#"distance_function\tcut_off_p_value\tno_map_dist\tno_map_size\tobjective_function\tTotalLGs\tLGs\tScaffolds\tLength\n";

#    foreach my $distance_function ("kosambi","haldane") {
#        foreach my $pvalue ("0.001", "0.0001", "0.00001", "0.000001", "0.0000001", "0.00000001") {
#            foreach my $map_dist (5,10,15,20,25) {
#                foreach my $map_size (0,1,2,3,4) {
#                    foreach my $objective_function ("COUNT","ML") {
#                        $mstheader{distance_function} = $distance_function;
#                        $mstheader{cut_off_p_value} = $pvalue;
#                        $mstheader{no_map_dist} = $map_dist;
#                        $mstheader{no_map_size} = $map_size;
#                        $mstheader{objective_function} = $objective_function;
    open my $markers, ">", "$output_prefix.markers"
      or croak "Can't open $output_prefix.markers: $OS_ERROR\n";

    $mstheader{"number_of_loci"}       = $marker_num;
    $mstheader{"number_of_individual"} = $num_samples;
    map { print $markers "$_ $mstheader{$_}\n"; } @mstheader;

    print $markers $marker_out;
    close $markers;

    system(
"MSTMap.exe $output_prefix.markers $output_prefix.mstmap.map > $output_prefix.mstmap.log"
    );

    system(
"assess_map_genome_coverage.pl -a ../../Hmel1-1_Release_20120601/AGP/Hmel1-1_primaryScaffolds.agp -s $output_prefix.markerpos -m $output_prefix.mstmap.map -c ../../Hmel1-1_Release_20120601/AGP/Hmel1-1_chromosome.agp -v $output_prefix"
    );

    #                    }
    #                }
    #            }
    #        }
    #    }

    return;
}

sub mat_inter_match {
    my ( $intercross, $maternal ) = @_;
    my @intercross = split //, $intercross;
    my @maternal   = split //, $maternal;

    return 0 if ( @intercross != @maternal );
    foreach my $i ( 0 .. $#maternal ) {
        next if ( $intercross[$i] eq "X" );
        return 0 if ( $maternal[$i] ne $intercross[$i] );
    }
    return 1;
}

sub output_marker {
    my ( $marker_num, $marker_out, $pattern, $markerpos, $data_ref,
        $intercross ) = @_;
    $marker_num++;
    $marker_out .= "P$marker_num\_$data_ref->{count}";
    map { $marker_out .= "\t$_"; } split //, $pattern;
    $marker_out .= "\n";

    foreach my $scaffold ( keys %{ $data_ref->{bins} } ) {
        foreach my $bin_pos (
            sort { $a <=> $b }
            keys %{ $data_ref->{bins}{$scaffold} }
          )
        {
            print $markerpos
              "P$marker_num\_$data_ref->{count}\t$scaffold\t$bin_pos\t$pattern";
            print $markerpos "\t$intercross" if defined $intercross;
            print $markerpos "\n";
        }
    }
    return ( $marker_num, $marker_out, $markerpos );
}

1;    # Needed to return a true value when imported by vcf-parallel
