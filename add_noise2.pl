#!/usr/bin/perl -w
use strict;

# Adds random noise to columns, simulating changes that are non-inherited.
# Eg homoplasy - convergent evolution rather than inherited.

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

# Pick columns to add, sorted in descending order to make adding easy
my $len = length($seq[0]);
my $nedits = int($len * $fract);
my @pos;   # position of bases
my @base1; # base 1 at pos
my @base2; # base 2 at pos
my @freq;  # probability of base1 (vs base2)

for (my $i = 0; $i < $nedits; $i++) {
    $pos[$i] = int(rand($len));
    $base1[$i] = $bases[int(rand(4))];
    do {
	$base2[$i] = $bases[int(rand(4))];
    } while $base2[$i] eq $base1[$i];
    $freq[$i] = rand();
}
@pos = sort {$b <=> $a} @pos;

for (my $n = 0; $n < $count; $n++) {
    open(FH, ">", $outfn . "." . $n) || die;

    for (my $snum = 0; $snum < scalar(@name); $snum++) {
	my $name = $name[$snum];
	my $seq  = $seq[$snum];

	for (my $i = 0; $i < $nedits; $i++) {
	    my $b = rand() < $freq[$i] ? $base1[$i] : $base2[$i];
	    substr($seq, $pos[$i], 0) = $b;
	}
	print FH "$name\n$seq\n";
    }
    close(FH);
}
