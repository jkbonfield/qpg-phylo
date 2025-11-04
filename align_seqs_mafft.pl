#!/usr/bin/perl -w
use strict;

# align_seqs_mafft.pl seqs.fa

my $seq = shift(@ARGV);

open(FH, "mafft --auto --thread 6 $seq 2> mafft.out |") || die;
my $first = 1;
while (<FH>) {
    if (/^>/) {
	print "\n" unless $first;
	print;
	next;
    }

    chomp();
    tr/acgt-/ACGTN/;
    tr/ACGTN/N/c;
    print "$_";
    $first = 0;
}
close(FH);

print "\n";
