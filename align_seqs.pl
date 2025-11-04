#!/usr/bin/perl -w
use strict;

# align_seqs.pl reference seqs.fa

my $ref = shift(@ARGV);
my $seq = shift(@ARGV);

# Cheap fasta parsing to get length!
my $ref_len = `samtools view $ref|awk '{print length(\$10)}'`;
$ref_len =~ s/\n//;

open(FH, "minimap2 -a $ref $seq 2>/dev/null |") || die;

my %uses_ref = qw/S 0 H 0 I 0 M 1 = 1 X 1 D 1/;
my %uses_seq = qw/S 1 H 0 I 1 M 1 = 1 X 1 D 0/;

while(<FH>) {
    next if (/^@/);
    chomp($_);
    my @F = split("\t", $_);
    next if $F[1] & 0x900; # supplementary or secondary

    # Leading rseq is "N" to POS
    my $rpos = $F[3]-1;
    my $spos = 0;
    my $seq = $F[9];

    # Compute reference-aligned seq
    my $aseq = "N" x $rpos;
    while ($F[5] =~ m/(\d+)(\S)/g) {
	if ($uses_ref{$2}) {
	    $rpos += $1;
	    if ($uses_seq{$2}) {
		$aseq .= substr($seq, $spos, $1);
	    } else {
		#$aseq .= "-" x $1;
		$aseq .= "-" x $1;
	    }
	}

	$spos += $1 * $uses_seq{$2};
    }

    $aseq .= "N" x ($ref_len - $rpos);

    #my $subseq = $aseq;
    #$subseq =~ tr/-N//dc;
    #next if length($subseq) > 500;

    $aseq =~ tr/ACGT/NNNN/c;

    # FASTA
    print ">$F[0]\n$aseq\n";

    # Minimally aligned SAM
    #print "$F[0]\t0\tref\t1\t0\t${ref_len}M\t*\t0\t0\t$aseq\t*\n";
}
