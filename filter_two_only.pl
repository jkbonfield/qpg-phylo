#!/usr/bin/perl -w
use strict;

srand(0);

# Reads a FASTQ file and edits the sequence to only contain two symbols.
# When Ns appear, we also replace those too


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

$"="\t";
my $seq_len = length($seq);
for (my $i = 0; $i<$seq_len; $i++) {
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
    $nv++ if $nn;

    # Exclude N from our chosen pair
    $b{N} = 0;
    my @a = sort {$b{$b} <=> $b{$a}} keys %b;
    
    if ($nv > 2) {
	# Guess!  We just assign unknown bases to the called one
	# in proportion to their observations.
	# THIS IS NOT CORRECT, but it gives us consistent data sets.
	my $freq1 = $b{$a[0]} / ($b{$a[0]}+$b{$a[1]});
	foreach (keys %seqs) {
	    my $base = substr($seqs{$_},$i,1);
	    next if ($base eq $a[0] || $base eq $a[1]);
	    my $edit = rand()<$freq1 ? $a[0] : $a[1];
	    substr($seqs{$_},$i,1) = $edit;
	    #print "Edit $base to $edit: $a[0]/$a[1] $freq1\n";
	}
    }
#    print "$i $nv\t$na\t$nc\t$ng\t$nt\t$nn\t$a[0]$a[1]\n";
}

# Print up the new seqs
foreach (keys %seqs) {
    print "$_\n";
    my $s = $seqs{$_};
    $s=~s/(.{70})/$1\n/g; #line wrap
    print "$s\n";
}
