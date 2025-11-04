#!/usr/bin/perl -w

use strict;

# Reads a FASTA file with COVID lineages on the 2nd word per line.

# Loads lineages
my %lineage2seq;
my %lineage;
my $NR = 0;
while (<>) {
    next unless /^>/;
    chomp($_);
    tr/>//d;
    my @F = split("\t", $_);
    $lineage{$F[1]}=1;
    push(@{$lineage2seq{$F[1]}}, $F[0]);
    #push(@{$lineage2seq{$F[1]}}, sprintf("%04d#%s", ++$NR, $F[1]));
}

# Iterates through lineages in sorted order.
# The naming means while not perfectly sorted (1, 10, 11...19, 2, 20...),
# the dots mean we still have things in tree order with parents then children.
my $depth = 0;
my $last_depth = 0;
my $parent = "";
my $last_node = "";

my %tree;

sub populate_parent {
    ($_) = @_;
    return if $_ eq "root";
    $parent = $_;
    $parent =~ s/\.[^.]*$//;
    $parent = "root" if ($parent eq $_);

    # Don't add dups.  Change to a hash
    foreach my $n (@{$tree{$parent}}) {
	return if ($n eq $_);
    }

    #print "push(@\$tree{$parent}, $_)\n";
    push(@{$tree{$parent}},$_);
    populate_parent($parent) unless $parent eq "root";
}

foreach (sort keys(%lineage)) {
    #$depth = $_;
    #$depth=~tr/.//dc;
    #$depth=length($depth);
    $parent = $_;
    $parent =~ s/\.[^.]*$//;
    $parent = "root" if $parent eq $_;
    #print $depth," "," "x$depth, $_, " "x(80-$depth-length($_)),"$parent\n";

    push(@{$tree{$parent}},$_);
    populate_parent($parent);
}

sub print_tree {
    local $"=", ";

    my @members = ();
    foreach (sort @_) {
	if (exists($lineage2seq{$_})) {
	    push(@members, @{$lineage2seq{$_}});
	}
	if (exists($tree{$_})) {
	    push(@members, @{$tree{$_}});
	}
    }
    my @ele = ();
    foreach (@members) {
	my $t = print_tree($_);
	if ($t) {
	    push(@ele, "($t)$_");
	} else {
	    push(@ele, $_);
	}
    }

    return "@ele";
}

my $tree = print_tree("root");
print "($tree)root;\n";
