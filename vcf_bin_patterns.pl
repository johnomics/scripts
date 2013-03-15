
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

);

sub process {
    my ( $scf_vcf_ref, $sample_names_ref, $data_ref, $parent_string ) = @_;

    my $bin_size = 10000;
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
                    $bin{ $calls{$parent_call}{type} }
                      { $sample_names_ref->[ $i - 9 ] }{"-"}{$pos} = 0;
                }
                else {
                    my @sample_f = split /:/, $fields[$i];
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

    foreach my $type ( keys %{$bin_ref} ) {
        next if (keys %{$bin_ref->{$type}} == 0);
        for my $i ( 0 .. @{$sample_names_ref} - 4 ) {
            my $sample = $sample_names_ref->[$i];
            next if ( defined $ignore{$sample} );
            my %sum_gt;
            foreach my $gt ( "A", "B", "-" ) {
                foreach my $pos ( keys %{ $bin_ref->{$type}{$sample}{$gt} } ) {
                    $sum_gt{$gt} += $bin_ref->{$type}{$sample}{$gt}{$pos};
                }
            }
            my @gts = reverse sort { $sum_gt{$a} <=> $sum_gt{$b} } keys %sum_gt;
            my $max_gt  = $gts[0];
            my $max_gtq = abs( $sum_gt{ $gts[0] } - $sum_gt{ $gts[1] } );
            $data_ref->{bins}{$type}{$scaffold}{$bin_pos}{pattern} .= $max_gt;
            $data_ref->{bins}{$type}{$scaffold}{$bin_pos}{quals} .=
              "\t$max_gtq";
        }
        $data_ref->{patterns}{$type}
          { $data_ref->{bins}{$type}{$scaffold}{$bin_pos}{pattern} }++;
    }
    return;
}

sub merge {
    my ( $data_ref, $final_ref ) = @_;

    foreach my $type ( keys %{ $data_ref->{bins} } ) {
        foreach my $scf ( keys %{ $data_ref->{bins}{$type} } ) {
            foreach my $pos ( keys %{ $data_ref->{bins}{$type}{$scf} } ) {
                $final_ref->{bins}{$type}{$scf}{$pos}{pattern} =
                  $data_ref->{bins}{$type}{$scf}{$pos}{pattern};
                $final_ref->{bins}{$type}{$scf}{$pos}{quals} =
                  $data_ref->{bins}{$type}{$scf}{$pos}{quals};
            }
        }
    }
    foreach my $type ( keys %{ $data_ref->{patterns} } ) {
        foreach my $pattern ( keys %{ $data_ref->{patterns}{$type} } ) {
            $final_ref->{patterns}{$type}{$pattern} +=
              $data_ref->{patterns}{$type}{$pattern};
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
        open my $output, ">", "$output_prefix.$type.bins.txt"
          or croak
          "Can't open output file $output_prefix.$type.bins.txt! $OS_ERROR\n";
        foreach my $scf ( sort keys %{ $final_ref->{bins}{$type} } ) {
          POS:
            foreach my $pos ( sort { $a <=> $b }
                keys %{ $final_ref->{bins}{$type}{$scf} } )
            {
                next
                  if ( $final_ref->{bins}{$type}{$scf}{$pos}{pattern} =~ "-" );
                next if ( $final_ref->{patterns}{$type}{$final_ref->{bins}{$type}{$scf}{$pos}{pattern}} <= 1 );
                my @gqs = split /\t/,
                  $final_ref->{bins}{$type}{$scf}{$pos}{quals};
                shift @gqs;
                my $pos_ok = 1;
                map { $pos_ok = 0 if ( $_ < 99 ); } @gqs;
                if ($pos_ok) {
                    print $output
"$scf\t$pos\t$final_ref->{bins}{$type}{$scf}{$pos}{pattern}$final_ref->{bins}{$type}{$scf}{$pos}{quals}\n";
                }
            }
        }
    }
    close $output;

    open my $markers, ">", "$output_prefix.markers"
      or croak "Can't open $output_prefix.markers: $OS_ERROR\n";
    print $markers "locus_name";
    map {
        if ( !defined $ignore{ $sample_names_ref->[$_] } ) {
            print $markers "\t$sample_names_ref->[$_]";
        }
    } 0 .. $#{$sample_names_ref} - 3;
    print $markers "\n";

    my $pnum = 0;
    foreach my $pattern (
        sort {
            $final_ref->{patterns}{"paternal"}{$a}
              <=> $final_ref->{patterns}{"paternal"}{$b}
        } keys %{ $final_ref->{patterns}{"paternal"} }
      )
    {
        next
          if ( $pattern =~ "-" );
        next if ( $final_ref->{patterns}{"paternal"}{$pattern} <= 1 );
        $pnum++;
        print $markers "P$pnum\_$final_ref->{patterns}{paternal}{$pattern}";
        map { print $markers "\t$_"; } split //, $pattern;
        print $markers "\n";
    }
    close $markers;
    return;
}

1;    # Needed to return a true value when imported by vcf-parallel
