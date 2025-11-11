#!/usr/bin/perl -w
use strict;

# An experiment.  Encode bases as 4-nucleotide strings.
#     A = AAAAC
#     C = AAACA
#     G = AACAA
#     T = ACAAA
#     N = NAAAA

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
    $_=~s/a/AAAAC/g;
    $_=~s/c/AAACA/g;
    $_=~s/g/AACAA/g;
    $_=~s/t/ACAAA/g;
    $_=~s/n/NAAAA/g;
    print "$_\n";
}
