#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use English;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

use Term::ExtendedColor qw/:all/;
use Number::Format qw/:subs :vars/;

$OUTPUT_AUTOFLUSH = 1;

my %chrommarker = (
    "BBAABBAABAAABBBAAAAAAAABABBBAABBBABBBBBBBAABBBABBBBABBAAAABABABABBAAA" =>
      "7",
    "AABBAABBABBBAAABBBBBBBBABAAABBAAABAAAAAAABBAAABAAAABAABBBBABABABAABBB" =>
      "7",
    "ABBAABBAABBBAABBBBBAAABBBBAABBBAABAABBBABABBAAABAABBAABBAAAAAABBAABAA" =>
      "19",
    "BAABBAABBAAABBAAAAABBBAAAABBAAABBABBAAABABAABBBABBAABBAABBBBBBAABBABB" =>
      "19",
    "ABAAAABAAAABAABAAAAABAABBAABABBAABAABBAABBBAAAABBAABABBBABBABBBAAABBB" =>
      "18",
    "BABBBBABBBBABBABBBBBABBAABBABAABBABBAABBAAABBBBAABBABAAABAABAAABBBAAA" =>
      "18",
    "BBBAABBAAAABAAABABABAABAABBBBBABBBAABABBAAABABBABBBABBBAAAAAABBAAAABB" =>
      "12",
    "AAABBAABBBBABBBABABABBABBAAAAABAAABBABAABBBABAABAAABAAABBBBBBAABBBBAA" =>
      "12",
    "ABBBBBBABABAABAAAABBBBBAABABBABABAABAAABABBBABBABBABBAABBABBAABBAABBB" =>
      "1",
    "BAAAAAABABABBABBBBAAAAABBABAABABABBABBBABAAABAABAABAABBAABAABBAABBAAA" =>
      "1",
    "BABBAAAABBBBABABBAAABBBBBBAAAAAAAABBBBBABBBABABAAABABABBABBBAABBAAABB" =>
      "13",
    "ABAABBBBAAAABABAABBBAAAAAABBBBBBBBAAAAABAAABABABBBABABAABAAABBAABBBAA" =>
      "13",
    "BBABAAAAABAABBABABBABBABABABBBBBABBBBABBAAAABBABBAAAABBBABABBABBBBABA" =>
      "6",
    "AABABBBBBABBAABABAABAABABABAAAAABAAAABAABBBBAABAABBBBAAABABAABAAAABAB" =>
      "6",
    "BAAAAABAABBBBAAAAAABAABAAABBBABBAABABBBBBBBBAABABBABAABAABBBBBABBAABB" =>
      "15",
    "ABBBBBABBAAAABBBBBBABBABBBAAABAABBABAAAAAAAABBABAABABBABBAAAAABAABBAA" =>
      "15",
    "ABABABAABAABBBBBBABBBAABBBBAABAABAABAABAAAAABAABBAAAAAAAABBBBABBABBBA" =>
      "3",
    "BABABABBABBAAAAAABAAABBAAAABBABBABBABBABBBBBABBAABBBBBBBBAAAABAABAAAB" =>
      "3",
    "ABBAABBABBBABBBBBBAAAABABBAAABBBBAABABAABBAAABAAABBBBBBBBBABBBAAABABB" =>
      "16",
    "BAABBAABAAABAAAAAABBBBABAABBBAAAABBABABBAABBBABBBAAAAAAAAABAAABBBABAA" =>
      "16",
    "BABBAAABAABBABBAABBBBBBABAABBABBBBBAABBBBBBBABAABABBBAABABBBBABABAABA" =>
      "11",
    "ABAABBBABBAABAABBAAAAAABABBAABAAAAABBAAAAAAABABBABAAABBABAAAABABABBAB" =>
      "11",
    "AABBABBAABBBAAABAABABAAAAABBABBAAAABBABAABABABBBAABAABAAABABAAAAAAAAB" =>
      "8",
    "BBAABAABBAAABBBABBABABBBBBAABAABBBBAABABBABABAAABBABBABBBABABBBBBBBBA" =>
      "8",
    "BAAABBABBAABABBABABABAABBABBAABAABBAABBBBABBBAAABBBBABBAAABBBAAABBABA" =>
      "17",
    "ABBBAABAABBABAABABABABBAABAABBABBAABBAAAABAAABBBAAAABAABBBAAABBBAABAB" =>
      "17",
    "BBABAABAAABABABAAAABAAAAABBAABABABBABAABABBBBBAAABBBAABBAAAABAAAAABAB" =>
      "5",
    "AABABBABBBABABABBBBABBBBBAABBABABAABABBABAAAAABBBAAABBAABBBBABBBBBABA" =>
      "5",
    "BBAAABAABABBAABAABAAABBABAAABAABABBBBBBABABABBBBBAAAAABBAABAABBABAABB" =>
      "9",
    "AABBBABBABAABBABBABBBAABABBBABBABAAAAAABABABAAAAABBBBBAABBABBAABABBAA" =>
      "9",
    "BAABAAAAABAABABBBBBABBAAAAAABBBBBABAABBABABABBABAABBABAABBBAAABBAABAA" =>
      "4",
    "BAABABAAAABBAAAAAAABBBBBABABBBAAAABABABAAABAABABAABBBBBBABAABBAABBABA" =>
      "Z",
    "ABBABABBBBAABBBBBBBAAAAABABAAABBBBABABABBBABBABABBAAAAAABABBAABBAABAB" =>
      "Z",
    "ABBBABABBAABBBBBBABBAABBBBBBABBBABBABABAABBAAABBABBABBAAAAAAAABBAABAB" =>
      "2",
    "BAAABABAABBAAAAAABAABBAAAAAABAAABAABABABBAABBBAABAABAABBBBBBBBAABBABA" =>
      "2",
    "BAAABABBBAAAABABBAAAABBAABAABAABAABAAAAAAAAABBBAABBBBAAAAAAABAAAABBBB" =>
      "20",
    "ABBBABAAABBBBABAABBBBAABBABBABBABBABBBBBBBBBAAABBAAAABBBBBBBABBBBAAAA" =>
      "20"
);

my %args;

$args{markerfile} = "";
$args{lengthfile} = "";

my $options_okay = GetOptions(
    'markers=s' => \$args{markerfile},
    'lengths=s' => \$args{lengthfile}
);
croak "No marker file! Please specify -m $OS_ERROR\n"
  if ( $args{markerfile} eq "" );

croak "No lengths file! Please specify -l $OS_ERROR\n"
  if ( $args{lengthfile} eq "" );

open my $lengthfile, '<', $args{lengthfile}
  or croak "Can't open $args{lengthfile}! $OS_ERROR\n";

my %scflen;
my $genomelen = 0;
while ( my $lengthline = <$lengthfile> ) {
    chomp $lengthline;
    my ( $scf, $len ) = split /\t/, $lengthline;
    $scflen{$scf} = $len;
    $genomelen += $len;
}
close $lengthfile;

my $scfnum = keys %scflen;
print "Genome has $scfnum scaffolds and is $genomelen bp long\n";

open my $markerfile, '<', $args{markerfile}
  or croak "Can't open $args{markerfile}! $OS_ERROR\n";

my @scflines;
my $prev_scf;
my %chrom;
my $scf_count = 0;
while ( my $marker = <$markerfile> ) {
    next
      if ( $marker =~ /^(-+)$/ );
    chomp $marker;
    my ( $scf, $pos, $type, $mq, $fs, $orig, $edge, $cons, $corrected ) =
      split /\t/,
      uncolor($marker);
    next if ( $type ne "maternal" );
    $prev_scf = $scf if !defined $prev_scf;

    if ( $scf ne $prev_scf ) {
        $scf_count++;
        print "."            if ( $scf_count % 10 == 0 );
        print "$scf_count\n" if ( $scf_count % 100 == 0 );

        process_scf( $prev_scf, \@scflines );
        @scflines = ();
    }
    $prev_scf = $scf;
    push @scflines, "$pos:$corrected";
}
process_scf( $prev_scf, \@scflines );
close $markerfile;

print "\n";

my $genfound = 0;
my %scffound;
foreach my $chr ( sort {$a<=>$b} keys %chrom ) {
    foreach my $scf ( keys %{ $chrom{$chr}{scf} } ) {
        $chrom{$chr}{scflen} += $scflen{$scf};
        $scffound{$scf}++;
    }
    my $scfnum = keys %{ $chrom{$chr}{scf} };
    print "$chr\t";
    printf "%4d\t", $scfnum;
    print format_picture( $chrom{$chr}{markerlen}, '###,###,###' ) . "\t"
      . format_picture( $chrom{$chr}{scflen}, '###,###,###' ) . "\n";
    $genfound += $chrom{$chr}{markerlen};
}

my $scffound = keys %scffound;
printf "TOTAL\t%4d", $scffound;
print "\t" . format_picture( $genfound, "###,###,###" ) . "\n";

sub process_scf {
    my ( $scf, $scflines ) = @_;

    my $prev_pattern;
    my $start;
    my $prev_pos;

    my $snps = @{$scflines};
    foreach my $marker ( @{$scflines} ) {
        my ( $pos, $pattern ) = split /:/, $marker;
        $prev_pattern = $pattern if !defined $prev_pattern;
        $start        = $pos     if !defined $start;
        if ( $pattern ne $prev_pattern ) {
            process_marker( $scf, $start, $prev_pos, $prev_pattern );
            $start = $pos;
        }
        $prev_pattern = $pattern;
        $prev_pos     = $pos;
    }
    process_marker( $scf, $start, $prev_pos, $prev_pattern );
}

sub process_marker {
    my ( $scf, $start, $prev_pos, $prev_pattern ) = @_;
    my $length = $prev_pos - $start + 1;
    my $chr = $chrommarker{$prev_pattern} // "0";
    $chrom{$chr}{scf}{$scf}++;
    $chrom{$chr}{markerlen} += $length;
}
