#!/usr/bin/perl -w
use strict;

# An experiment.  Encode bases as 4-nucleotide strings.
#     A = AAAC
#     C = AACA
#     G = ACAA
#     T = CAAA
#     N = NNNN

# The test is to see whether the tree changes shape.  I believe it
# should not as the information on variation is still valid, although
# the distances may change as an G to T mutation in a column becomes
# a pair of mutations (AC in both cases).

while(<>) {
    if (/^>/) {
	print;
	next;
    }

    chomp($_);
    $_=~tr/ACGTN/acgtn/;
    $_=~s/a/AAAC/g;
    $_=~s/c/AACA/g;
    $_=~s/g/ACAA/g;
    $_=~s/t/CAAA/g;
    $_=~s/n/NNNN/g;
    print "$_\n";
}
