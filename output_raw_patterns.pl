use List::Util qw/sum min max/;
use POSIX qw/ceil/;
use Term::ExtendedColor qw/:all/;

use strict;
use warnings;

my %callorder;

my %nullcall = ( 'GT' => './.', 'GQ' => 0 );

my @mask;

sub process {
    my ( $scf, $scfref, $samples, $data, $genetics ) = @_;
    $data->{$scf} = get_markers( $scfref, $samples, $genetics );
    output_scf_to_file( $scf, $data->{$scf}, $samples );
}

sub get_markers {
    my ( $scfref, $samples, $genetics ) = @_;

    my %markers;
    foreach my $snp ( @{$scfref} ) {
        chomp $snp;
        my ( $marker, $type, $pos ) =
          parse_snp( $snp, $samples, $genetics );
        $markers{$type}{$pos}{marker} = $marker;
        $markers{$type}{$pos}{parent} = $type;
    }
    \%markers;
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

    foreach my $sample ( @{ $samples->{offspring}{order} } ) {
        my $opos = $samples->{offspring}{lookup}{$sample};
        $marker{$sample}{gt} = get_cp( 'GT', $f[$opos], $callf );
        $marker{$sample}{gq} = get_cp( 'GQ', $f[$opos], $callf );
    }

    return ( \%marker, $parentcall, $f[1] );
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
    my ( $part, $all ) = @_;
    foreach my $scf ( keys %{$part} ) {
        $all->{scf}{$scf}++;
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
    my ( $scf, $data, $samples ) = @_;
    open my $scfhandle, '>', "$scf.tmp.out"
      or croak "Can't open output file for $scf: $OS_ERROR\n";
    output_scf( $scf, $data, $scfhandle, $samples );
    close $scfhandle;
}

sub output_scf {
    my ( $scf, $data, $handle, $samples ) = @_;
    foreach my $type ( sort keys %{$data} ) {
        foreach my $pos ( sort { $a <=> $b } keys %{ $data->{$type} } ) {

            print $handle "$scf\t$pos\t$type\t";
            foreach my $sample ( @{ $samples->{offspring}{order} } ) {
                my $col =
                  ceil( 255 - $data->{$type}{$pos}{marker}{$sample}{gq} / 4.2 );
                print $handle fg $col, $data->{$type}{$pos}{marker}{$sample}{gt};
                print $handle " ";
            }
            print $handle "\n";
        }
        print $handle '-' x 253, "\n";
    }
    print $handle '-' x 253, "\n";
}

1;
