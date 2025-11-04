#!/usr/bin/perl -w
use strict;

# A noddy tool for generating uncorrelated random mutations in a
# set of sequences such that no column contains Ns and no column
# contains more than two base pairs.
#
# This is not a tree!

my $seq_len=1000;
my $snp_prob=0.1;
my $nseqs = 100;

# Create a base sequence to copy
my @base=qw/A C G T/;
my $seq_base = "";
for (my $i=0;$i<$seq_len;$i++) {
    $seq_base .= $base[int(rand(4))];
}

# Create an array of copies of the base sequence
my @seqs;
for (my $n=0;$n<$nseqs;$n++) {
    $seqs[$n]=$seq_base;
}

# Mutate specific columns in the seq array
for (my $i=0;$i<$seq_len;$i++) {
    next if rand()>$snp_prob;
    my $b1 = $base[int(rand(4))];
    my $b2 = $b1;
    while ($b1 eq $b2) {
	$b2 = $base[int(rand(4))];
    }
    for (my $n=0; $n<$nseqs; $n++) {
	substr($seqs[$n], $i, 1) = rand()<0.5 ? $b1 : $b2;
    }
}

for (my $n=0;$n<$nseqs;$n++) {
    print ">seq_$n\n$seqs[$n]\n";
}
