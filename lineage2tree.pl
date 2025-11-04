#!/usr/bin/perl -w

use strict;

# Reads COVID lineages (named [AB].<num>.<num>...) and creates a Newick tree.

# Loads lineages
my %lineage = ();
while (<>) {
    chomp($_);
    $lineage{$_}=1;
}

# Iterates through lineages in sorted order.
# The naming means while not perfectly sorted (1, 10, 11...19, 2, 20...),
# the dots mean we still have things in tree order with parents then children.
my $depth = 0;
my $last_depth = 0;
my $tree = "";
my $parent = "";
my $last_node = "";

my %tree;

sub populate_parent {
    ($_) = @_;
    $parent = $_;
    $parent =~ s/\.[^.]*$//;
    $parent = "root" if ($parent eq $_);

    if (!exists($tree{$parent})) {
	#print "push(@\$tree{$parent}, $_)\n";
	push(@{$tree{$parent}},$_);
	populate_parent($parent) unless $parent eq "root";
    }
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
    my $tree = "";
    my $first = 1;
    #print "print_tree @_\n";
    foreach (@_) {
	if (exists($tree{$_})) {
	    $tree .= ", " if ($first == 0);
	    $tree .= "( " . print_tree(@{$tree{$_}}) . " )" . $_;
	    $first = 0;
	} else {
	    if ($first == 0) {
		$tree .= ", ";
	    } else {
		$first = 0;
	    }
	    $tree .= "$_";
	}
    }
    return $tree;
}

$tree = print_tree("root");
print "$tree;\n";
