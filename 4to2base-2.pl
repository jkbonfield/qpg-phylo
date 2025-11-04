#!/usr/bin/perl -w
use strict;

# Like 4to2base.pl but we only do the translation for columns with more than
# 2 bp present (excluding Ns).

srand(0);

# Load sequences
my %seqs; # indexed on name, value is seq
my $nseq=0;

sub add_seq {
    my ($name,$seq) = @_;
    #print "Add $name ",length($seq),"\n";
    $seqs{$name}=$seq;
    $nseq++;
}

my $name = "";
my $seq = "";
while (<>) {
    chomp();
    if (/^>/) {
	add_seq($name,$seq) if ($name);
	$name=$_;
	$seq ="";
    } else {
	$seq .= $_;
    }
}

add_seq($name,$seq) if ($name);


# Scan columns and edit
sub base_freq {
    my ($pos) = @_;
    my %bases = qw/A 0 C 0 G 0 T 0 N 0/;
    foreach (keys %seqs) {
	my $base = substr($seqs{$_}, $pos, 1);
	$bases{$base}++;
    }
    return %bases;
}

my %sub = qw/A CAAA C ACAA G AACA T AAAC N NNNN/;

$"="\t";
my $seq_len = length($seq);
for (my $i = $seq_len-1; $i>=0; $i--) {
    my %b = base_freq($i);
    my $na = $b{A};
    my $nc = $b{C};
    my $ng = $b{G};
    my $nt = $b{T};
    my $nn = $b{N};
    my $n = $na + $nc + $ng + $nt + $nn;
    my $nv = 0;
    $nv++ if $na;
    $nv++ if $nc;
    $nv++ if $ng;
    $nv++ if $nt;

    if ($nv > 2) {
	# Replace with a 4bp recognition sequence instead that guarantees
	# only 2 bases per column
	foreach (keys %seqs) {
	    my $base = substr($seqs{$_},$i,1);
	    substr($seqs{$_},$i,1) = $sub{$base};
	}
    }
}

# Print up the new seqs
foreach (keys %seqs) {
    print "$_\n";
    my $s = $seqs{$_};
    $s=~s/(.{70})/$1\n/g; #line wrap
    print "$s\n";
}
