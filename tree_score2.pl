#!/usr/bin/perl -w
use strict;

# WIP: not useful

# From https://en.wikipedia.org/wiki/Newick_format
#
# Tree -> Subtree ";"
# Subtree -> Leaf | Internal
# Leaf -> Name
# Internal -> "(" BranchSet ")" Name
# BranchSet -> Branch | Branch "," BranchSet
# Branch -> Subtree Length
# Name -> empty | string
# Length -> empty | ":" number
sub parse_tree {
    my ($tree) = @_;
    my %tree = ();    # internal node name -> children
    my %parent = ();  # child name -> parent node name
    my %leaf = ();    # set if leaf node

    my $node_num=0;

    # Remove trailing semicolon
    $tree =~ s/;$//;

    my @pos;
    do {
	#print "=== $tree ===\n";
	@pos = ();

	$_ = $tree;
	# Cache positions of a leaf node branchSet
	while (/\(([^()]+)\)([^,:()]*)((?::[0-9.e]*)?)/gc) {
	    #print "@- // @+\n";
	    #print "First: ",substr($tree,$-[1], $+[1]-$-[1]), "\n";
	    #print "All:   ",substr($tree,$-[0], $+[0]-$-[0]), "\n";
	    my $set    = substr($tree,$-[1], $+[1]-$-[1]);
	    my $name   = substr($tree,$-[2], $+[2]-$-[2]);
	    my $length = substr($tree,$-[3], $+[3]-$-[3]);
	    push(@pos, [$-[0], $+[0]-$-[0], $set, $name, $length]);
	}

	# Modify the tree
	foreach (reverse @pos) {
	    my ($st, $len, $list, $node, $length)=@$_;
	    $length =~ s/^://;
	    $list =~ s/\s+//g;

	    $node = ":$node_num" unless $node ne "";
	    $node = "#$node";
	    $node_num++;

	    #print "tree{$node} = \"$list\"\n";
	    $tree{$node}=$list;
	    foreach (split(",", $list)) {
		$parent{$_} = $node;
		$leaf{$_} = 1 if (!/^#/);
	    }

	    # Replace branchSet with an internal node
	    substr($tree,$st,$len) = $node;
	    #print "Tree is now <$tree>\n";
	}
    } while (scalar(@pos));

#    $tree{"root"} = $tree;
#    foreach (split(",", $tree)) {
#	print "parent{$_} = root\n";
#	$parent{$_} = "root";
#    }
    
    return ($tree, \%tree, \%parent, \%leaf);
}

my $depth=0;
sub print_tree {
    my ($tree, $node) = @_;
    foreach (sort $tree->{$node}) {
	print "  "x$depth,"$node=($_)\n";
	my @nodes = split(",",$_);
	foreach (@nodes) {
	    if (defined(%{$tree}{$_})) {
		$depth++;
		print_tree($tree,$_);
		$depth--;
	    }
	}
    }
}

sub print_traceback {
    my ($parent, $node) = @_;
    print "$node\n";
    print_traceback($parent, $parent->{$node}) if (defined $parent->{$node});
}

# NB: Change lineage regexp below to /^(.)/
#my $str = "((A1,A2),(A3,A4));";
#my $str = "((A1,A2,(A3,A4)),((B1,B2),(B3,B4)));";
#my $str = "((A1,B2,(A3,A4)),((B1,A2),(B3,B4)));";
#my $str = "(A,B,(C,D)E)F;";
#my $str = "(A,B,(C,D));";

my $str = shift(@ARGV);

my ($root, $tree, $parent, $leaf) = parse_tree($str);

#print_tree($tree, $root);
#print "\n\n\n";
#print_traceback($parent, "B3");

# Now bubble up to the root from every leaf node, accumulating lineage stats
# per node.
my %lineages = ();
foreach (sort keys(%{$leaf})) {
    #my ($lineage) = $_ =~ /_.*_([^:]*)/;
    #my ($lineage) = $_ =~ /^(.)/; # Ax, Bx, Cx etc.
    my ($lineage) = $_ =~ /^([A-Z.]*)/; # A.B:s1 A.B:s2 A:s3 
    #print "$_ $lineage $parent->{$_}\n";
    my $node = $parent->{$_};
    my $x = 0;
    while ($node) {
	$lineages{$node}->{$lineage}++;
	#print "  "x++$x, "$node\n";
	$node = $parent->{$node};
    }
}

sub min {
    my ($a, $b) = @_;
    return $a < $b ? $a : $b;
}

# We want one path for each lineage from child to root, excluding
# where the forks are all the same lineage.
# Eg (A,(A,(A,A))) is fine as internal nodes are still all "A".
# But ((A,B),(A,B)) is not as siblings both have a mixture.

my %root_counts;
foreach (sort keys(%{$leaf})) {
    #my ($lineage) = $_ =~ /_.*_([^:]*)/;
    my ($lineage) = $_ =~ /^(.)/; # Ax, Bx, Cx etc.
    my $node = $parent->{$_};
    my $x = 0;
    print "$_ $lineage $parent->{$_}\n";
    while ($node) {
	my @plin = sort(keys(%{$lineages{$node}}));
	print "  "x++$x, "$node // @plin\n";
	$node = $parent->{$node};
    }
}
