#!/usr/bin/perl -w
use strict;

# Adds random noise to sequencing.  Error rate in consensus basically.

# Usage add_noise fraction seed count file.fa
# Produces file.fa.<n> for count consecutive <n>s (0 upwards)

my $fract = shift(@ARGV);
srand(shift(@ARGV));
my $count = shift(@ARGV);
my $outfn = $ARGV[0];

my @bases = qw/A C G T/;

# Load seqs
my @name;
my @seq;
while (<>) {
    chomp($_);
    if (/^>/) {
	push(@name, $_);
    } else {
	push(@seq, $_);
    }
}


for (my $n = 0; $n < $count; $n++) {
    open(FH, ">", $outfn . "." . $n) || die;

    for (my $snum = 0; $snum < scalar(@name); $snum++) {
	my $name = $name[$snum];
	my $seq  = $seq[$snum];

	my $len = length($seq);
	my $nedits = int($len * $fract);
	for (my $i=0; $i<$nedits; $i++) {
	    my $p = int(rand($len));
	    my $b = substr($seq, $p, 1);
	    my $s = $b;
	    while ($b eq $s) {
		$s = $bases[int(rand(4))];
	    }
	    substr($seq, $p, 1) = $s;
	}
	print FH "$name\n$seq\n";
    }
    close(FH);
}
