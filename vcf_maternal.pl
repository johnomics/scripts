my %ignore = (
    "90"  => 0,
    "1N"  => 0,
    "108" => 0,
);

my %maternal = (
    "0/1 0/1 0/0" => { "0/0" => "A", "0/1" => "B", "1/1" => "-" },
    "0/0 0/1 0/0" => { "0/0" => "B", "0/1" => "A", "1/1" => "-" },
    "0/1 0/1 1/1" => { "0/0" => "-", "0/1" => "B", "1/1" => "A" },
    "1/1 0/1 1/1" => { "0/0" => "-", "0/1" => "A", "1/1" => "B" },
);

sub process {
    my ( $vcf_line, $sample_names_ref, $data_ref, $parent_string ) = @_;
    chomp $vcf_line;

    #    next if ( $vcf_line =~ /\.\/\./ );
    my @fields = split /\t/, $vcf_line;

    my @format_fields = split /:/, $fields[8];
    my %field_num;
    for my $i ( 0 .. ( @format_fields - 1 ) ) {
        $field_num{ $format_fields[$i] } = $i;
    }

    my %parents;
    my @parents;
    for my $i ( ( @{$sample_names_ref} + 6 ) .. ( @{$sample_names_ref} + 8 ) ) {
        push @parents, $sample_names_ref->[ $i - 9 ];
        if ( $fields[$i] eq "./." ) {
            $parents{ $sample_names_ref->[ $i - 9 ] } = "./.";
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

    next if ( !defined $maternal{$parent_call} );

    my $offspring_call;
    my @offspring_gq;
    my @offspring_dp;
    my @offspring_ada;
    my @offspring_adb;
    for my $i ( 9 .. ( @{$sample_names_ref} + 5 ) ) {
        next if ( defined $ignore{ $sample_names_ref->[ $i - 9 ] } );
        if ( $fields[$i] eq "./." ) {
            $offspring_call .= ".";
            push @offspring_gq, 0;
            push @offspring_dp, 0;
        }
        else {
            my @sample_f = split /:/, $fields[$i];
            $offspring_call .=
              $maternal{$parent_call}{ $sample_f[ $field_num{GT} ] };
            push @offspring_gq, $sample_f[ $field_num{GQ} ];
            push @offspring_dp, $sample_f[ $field_num{DP} ];
            my ( $ada, $adb ) = split /,/, $sample_f[ $field_num{AD} ];
            push @offspring_ada, $ada;
            push @offspring_adb, $adb;
        }
    }
    $data_ref->{pattern_scf}{ $fields[0] }{$offspring_call}++;
    $data_ref->{scaffold}{ $fields[0] }{ $fields[1] }{offspring} =
      $offspring_call;
    $data_ref->{scaffold}{ $fields[0] }{ $fields[1] }{parent} = $parent_call;
    $data_ref->{maternal}{$offspring_call}{count}++;

    my @attr = split /;/, $fields[7];
    map { $data_ref->{maternal}{$offspring_call}{mq} += $1 if (/^MQ=(.+)$/); }
      @attr;
    $data_ref->{maternal}{$offspring_call}{snpq} += $fields[5];

    map {
        $data_ref->{maternal}{$offspring_call}{gq}[$_]  += $offspring_gq[$_];
        $data_ref->{maternal}{$offspring_call}{dp}[$_]  += $offspring_dp[$_];
        $data_ref->{maternal}{$offspring_call}{ada}[$_] += $offspring_ada[$_];
        $data_ref->{maternal}{$offspring_call}{adb}[$_] += $offspring_adb[$_];
    } 0 .. $#offspring_gq;

    return;
}

sub merge {
    my ( $data_ref, $final_ref ) = @_;

    if ( defined($data_ref) ) {
        foreach my $scf ( keys %{ $data_ref->{pattern_scf} } ) {
            foreach my $pattern ( keys %{ $data_ref->{pattern_scf}{$scf} } ) {
                $final_ref->{pattern_scf}{$scf}{$pattern} +=
                  $data_ref->{pattern_scf}{$scf}{$pattern};
            }
        }
        foreach my $scf ( keys %{ $data_ref->{scaffold} } ) {
            foreach my $pos ( keys %{ $data_ref->{scaffold}{$scf} } ) {
                $final_ref->{scaffold}{$scf}{$pos}{offspring} =
                  $data_ref->{scaffold}{$scf}{$pos}{offspring};
                $final_ref->{scaffold}{$scf}{$pos}{parent} =
                  $data_ref->{scaffold}{$scf}{$pos}{parent};

            }
        }
        foreach my $val ( keys %{ $data_ref->{maternal} } ) {
            $final_ref->{maternal}{$val}{count} +=
              $data_ref->{maternal}{$val}{count};
            $final_ref->{maternal}{$val}{mq} += $data_ref->{maternal}{$val}{mq};
            $final_ref->{maternal}{$val}{snpq} +=
              $data_ref->{maternal}{$val}{snpq};
            map {
                $final_ref->{maternal}{$val}{gq}[$_] +=
                  $data_ref->{maternal}{$val}{gq}[$_];
                $final_ref->{maternal}{$val}{dp}[$_] +=
                  $data_ref->{maternal}{$val}{dp}[$_];
                $final_ref->{maternal}{$val}{ada}[$_] +=
                  $data_ref->{maternal}{$val}{ada}[$_];
                $final_ref->{maternal}{$val}{adb}[$_] +=
                  $data_ref->{maternal}{$val}{adb}[$_];

            } 0 .. $#{ $data_ref->{maternal}{$val}{gq} };
        }
    }

    return;
}

sub output {
    my (
        $final_data_ref, $scflen_ref, $sample_names_ref,
        $genome_length,  $output_prefix
    ) = @_;

    my %scf_chrs;
    open my $chromosome_agp, "<",
      "/whale-data/jd626/Hmel1-1_Release_20120601/AGP/Hmel1-1_chromosome.agp"
      or croak "Can't open chromosome AGP file! $OS_ERROR\n";
    while ( my $agp_line = <$chromosome_agp> ) {
        if ( $agp_line =~ /^chr/ && $agp_line =~ /\tD\t/ ) {
            chomp $agp_line;
            my (
                $chr, $start,     $end,     $part, $d,
                $scf, $scf_start, $scf_end, $dir,  $note
            ) = split /\t/, $agp_line;
            if ( $chr =~ /(.+)_unmapped/ ) {
                $chr = $1;
            }
            $scf_chrs{$scf}{chr}   = $chr;
            $scf_chrs{$scf}{start} = $start;
            $scf_chrs{$scf}{end}   = $end;
        }
    }

    my %pattern_scf;
    my %pattern_chr;
    open my $pattern_scf_file, ">", "$output_prefix.pattern_scf.txt"
      or croak "Can't open $output_prefix.pattern_scf.txt\n";
    open my $pattern_scf_max_file, ">", "$output_prefix.pattern_scf_max.txt"
      or croak "Can't open $output_prefix.pattern_scf_max.txt\n";

    foreach my $scf ( sort keys %{ $final_data_ref->{pattern_scf} } ) {
        my $max = 1;
        foreach my $pattern (
            reverse sort {
                $final_data_ref->{pattern_scf}{$scf}{$a}
                  <=> $final_data_ref->{pattern_scf}{$scf}{$b}
            } keys %{ $final_data_ref->{pattern_scf}{$scf} }
          )
        {
            $pattern_scf{$pattern}{scf}++;
            $pattern_scf{$pattern}{snps} +=
              $final_data_ref->{pattern_scf}{$scf}{$pattern};
            print $pattern_scf_file
"$scf\t$pattern\t$final_data_ref->{pattern_scf}{$scf}{$pattern}\n";
            if ($max) {
                print $pattern_scf_max_file
"$scf\t$pattern\t$final_data_ref->{pattern_scf}{$scf}{$pattern}\n";
                $max = 0;
            }
            $pattern_chr{ $scf_chrs{$scf}{chr} }{$pattern}{scf}++;
            $pattern_chr{ $scf_chrs{$scf}{chr} }{$pattern}{snps} +=
              $final_data_ref->{pattern_scf}{$scf}{$pattern};
        }
    }
    close $pattern_scf_file;
    close $pattern_scf_max_file;

    open my $pattern_chr_file, ">", "$output_prefix.pattern_chr.txt"
      or croak "Can't open $output_prefix.pattern_chr.txt\n";
    foreach my $chr ( sort keys %pattern_chr ) {
        foreach my $pattern (
            reverse sort {
                $pattern_chr{$chr}{$a}{snps} <=> $pattern_chr{$chr}{$b}{snps}
            } keys %{ $pattern_chr{$chr} }
          )
        {
            next if $pattern_chr{$chr}{$pattern}{snps} == 1;
            print $pattern_chr_file
"$chr\t$pattern\t$pattern_chr{$chr}{$pattern}{scf}\t$pattern_chr{$chr}{$pattern}{snps}\n";
        }
    }
    close $pattern_chr_file;

    open my $scf_file, ">", "$output_prefix.scaffolds.txt"
      or croak "Can't open $output_prefix.scaffolds.txt\n";
    foreach my $scf ( keys %{ $final_data_ref->{scaffold} } ) {
        foreach my $pos (
            sort { $a <=> $b }
            keys %{ $final_data_ref->{scaffold}{$scf} }
          )
        {
            print $scf_file
"$scf\t$pos\t$final_data_ref->{scaffold}{$scf}{$pos}{offspring}\t$final_data_ref->{scaffold}{$scf}{$pos}{parent}\n";
        }
    }
    close $scf_file;

    open my $stat_file, ">", "$output_prefix.maternal.txt"
      or croak "Can't open $output_prefix.maternal.txt\n";

    my @stats = ( "ada", "adb", "dp", "gq" );
    foreach my $val ( keys %{ $final_data_ref->{maternal} } ) {
        print $stat_file
          "$val\t$pattern_scf{$val}{scf}\t$pattern_scf{$val}{snps}\t";
        my @calls = split //, $val;
        my %av;
        my %min_type;
        my %av_type;
        my %calls;

        foreach my $type ( "A", "B" ) {
            foreach my $stat (@stats) {
                $min_type{$type}{$stat} = 250;
            }
        }

        map {
            foreach my $stat (@stats) {
                my $avstat =
                  $final_data_ref->{maternal}{$val}{count} > 0
                  ? int( $final_data_ref->{maternal}{$val}{$stat}[$_] /
                      $final_data_ref->{maternal}{$val}{count} )
                  : 0;
                push @{ $av{$stat} }, $avstat;
                $av_type{ $calls[$_] }{$stat} += $avstat;
                if ( $avstat < $min_type{ $calls[$_] }{$stat} ) {
                    $min_type{ $calls[$_] }{$stat} = $avstat;
                }
            }
            $calls{ $calls[$_] }++;
        } 0 .. $#calls;

        foreach my $type ( "A", "B" ) {
            foreach my $stat (@stats) {
                if ( $min_type{$type}{$stat} == 250 ) {
                    $min_type{$type}{$stat} = 0;
                }
            }
        }

        print $stat_file "$final_data_ref->{maternal}{$val}{count}\t";

        my $avmq =
          $final_data_ref->{maternal}{$val}{count} > 0
          ? int( $final_data_ref->{maternal}{$val}{mq} /
              $final_data_ref->{maternal}{$val}{count} )
          : 0;
        my $avsnpq =
          $final_data_ref->{maternal}{$val}{count} > 0
          ? int( $final_data_ref->{maternal}{$val}{snpq} /
              $final_data_ref->{maternal}{$val}{count} )
          : 0;
        print $stat_file "$avmq\t$avsnpq";

        foreach my $type ( "A", "B" ) {
            foreach my $stat (@stats) {
                my $av =
                  $calls{$type} > 0
                  ? int( $av_type{$type}{$stat} / $calls{$type} )
                  : 0;
                print $stat_file "\t$av\t$min_type{$type}{$stat}";
            }
        }
        map {
            printf $stat_file " %2d,%2d,%2d,%2d", $av{ada}[$_], $av{adb}[$_],
              $av{dp}[$_], $av{gq}[$_]
        } 0 .. $#calls;
        print $stat_file "\n";

    }
    close $stat_file;

    return;
}

1;    # Needed to return a true value when imported by vcf-parallel
