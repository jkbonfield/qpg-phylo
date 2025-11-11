#!/usr/bin/perl -w
use strict;

my $max_dist = 5e-5;
if ($ARGV[0] eq "-d") {
    shift(@ARGV);
    $max_dist = shift(@ARGV);
}

# NB doesn't (yet) flatten the tree, so it's misnamed.
# We assign lineages as before based on what's present in the children,
# but we then walk the tree, sorting child ordering by lineage.
# We then count the number of times we change lineage (ie run-length encode
# and count lines) and also count the number of times we switch lineage
# between leaf-children of the same parent.  These two factors produce a
# score on the clustering potential, although to doesn't factor in
# known covid lineage ancestry.

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
    my %length = ();

    my $node_num=0;

    # Remove trailing semicolon
    $tree =~ s/;$//;

    my @pos;
    do {
	#print "=== $tree ===\n";
	@pos = ();

	$_ = $tree;
	# Cache positions of a leaf node branchSet
	while (/\(([^()]+)\)([^,:()]*)((:[0-9.e+-]*)?)/gc) {
	    #print substr($_,$-[0], $+[0]-$-[0]),"\n";
	    #print "@- // @+\n";
	    #print "First: ",substr($tree,$-[1], $+[1]-$-[1]), "\n";
	    #print "All:   ",substr($tree,$-[0], $+[0]-$-[0]), "\n";
	    my $set    = substr($tree,$-[1], $+[1]-$-[1]);
	    my $name   = substr($tree,$-[2], $+[2]-$-[2]);
	    my $size   = substr($tree,$-[3], $+[3]-$-[3]);
	    my $length = $4;
	    push(@pos, [$-[0], $+[0]-$-[0], $set, $name, $size, $length]);
	}

	# Modify the tree string
	foreach (reverse @pos) {
	    my ($st, $nlen, $list, $node, $size, $length)=@$_;
	    $size =~ s/^://;
	    $list =~ s/\s+//g;
	    $length = 0 unless defined($length);
	    $length =~ s/^://;

	    $node = "/$node_num" unless $node ne "";
	    $node = "#$node";
	    $node_num++;

	    $length{$node} = $length;

	    #print "tree{$node} = \"$list\"\n";
	    foreach (split(",", $list)) {
		m/(^[^:]*)(?::(.*))?/;
		$_=$1 if (defined($1));
		#print "LLL $node $size $_ $1 // $2 //\n";
		$parent{$_} = $node;
		$leaf{$_} = 1 if (!/^#/);
		if (defined($2)) {
		    $length{$_} = $2;
		} else {
		    $length{$_} = 0 unless defined $length{$_};
		}
	    }
	    $list =~ s/:[0-9.e+-]*//g;
	    $tree{$node}=$list;

	    # Replace branchSet with an internal node
	    substr($tree,$st,$nlen) = $node;
	    #print "Tree is now <$tree>\n";
	}
    } while (scalar(@pos));

#    $tree{"root"} = $tree;
#    foreach (split(",", $tree)) {
#	print "parent{$_} = root\n";
#	$parent{$_} = "root";
#    }
    
    return ($tree, \%tree, \%parent, \%leaf, \%length);
}

my $depth=0;
sub print_tree {
    my ($tree, $length, $node) = @_;
    foreach (sort $tree->{$node}) {
	print "  "x$depth,"$node=($_) DIST: $length->{$node}\n";
	my @nodes = split(",",$_);
	foreach (@nodes) {
	    if (defined(%{$tree}{$_})) {
		$depth++;
		print_tree($tree, $length, $_);
		$depth--;
	    } else {
		print "  "x($depth+1), "CHILD: $_\tDIST: $length->{$_}\n";
	    }
	}
    }
}

sub tree_to_newick2 {
    my ($tree, $length, $node) = @_;
    local $"=",";
    my @ele = ();

    foreach my $n (split(",",$tree->{$node})) {
	if (defined(%{$tree}{$n})) {
	    my $name = $n;
	    if ($name =~ /^#:/) {
		$name = "";
	    } else {
		substr($name, 0, 1) = "";
	    }
	    my $len = "";
	    if (defined($length->{$n})) {
		my $l = $length->{$n};
		$len = ":$l";
	    }
	    push(@ele, "(" . tree_to_newick2($tree, $length, $n) . ")$name$len");
	} else {
	    my $len = "";
	    if (defined($length->{$n})) {
		my $l = $length->{$n};
		$len = ":$l";
	    }
	    push(@ele, $n . $len);
	}
    }

    return "@ele";
}

sub tree_to_newick {
    my ($tree, $length, $node) = @_;
    return "(" . tree_to_newick2($tree, $length, $node) . ");";
}

#print_tree($tree, $length, $root);

# Finds all nodes within a specific distance.
# Returns a list of close nodes,  list of far nodes, and a hash of node-to-distance.
# (All 3 as references)
sub nodes_within_dist {
    my @near;
    my @far;
    my %dist;
    my ($tree, $length, $node, $total_dist, $max_dist) = @_;
    foreach (sort $tree->{$node}) {
	#print "  "x$depth,"$node=($_) DIST: $length->{$node}\n";
	my @nodes = sort {$length->{$a} <=> $length->{$b}} split(",",$_);
	foreach (@nodes) {
	    my $dist = $total_dist + $length->{$_};
	    $dist{$_} = $dist;
	    if ($dist > $max_dist) {
		#print "FAR:  $_ $dist\n";
		push(@far, $_);
		next;
	    }
	    if (defined(%{$tree}{$_})) {
		$depth++;
		my ($nearr, $farr, $distr) = 
		    nodes_within_dist($tree, $length, $_, $total_dist + $length->{$_}, $max_dist);
		push(@near, @{$nearr});
		push(@far,  @{$farr});
		foreach (keys(%{$distr})) {
		    $dist{$_} = $distr->{$_};
		}
		$depth--;
	    } else {
		#print "  "x($depth+1), "CHILD: $_\tDIST: $length->{$_}  \t",
		#    $length->{$_} + $total_dist, "\n";
		#print "NEAR: $_ $dist <= $max_dist\n";
		push(@near, $_);
	    }
	}
    }
    return (\@near, \@far, \%dist);
}

my $merge_num = 0;
sub merge_nodes {
    my ($tree, $parent, $length, $node, $near, $far, $dist) = @_;

    print "MERGE $node <- (@$near) NOT (@$far)\n";

    # Check if this is a NOP already
    my @children = sort split(",", $tree->{$node});
    my @mergeable = sort @$near;
    return if ("@children" eq "@mergeable");
#    print "1: @children\n";
#    print "2: @mergeable\n";

    # FIXME: also check all @children vs all @children's direct parent and that
    # this parent is a direct child of node?

    # Find the minimum distance of the children in the near / far lists
    my $min_near_dist = 9e9;
    my $min_far_dist = 9e9;
    foreach (@$near) { $min_near_dist = $dist->{$_} if $min_near_dist > $dist->{$_}; }
    foreach (@$far)  { $min_far_dist  = $dist->{$_} if $min_far_dist  > $dist->{$_}; }
    $min_near_dist /= 2; # so tree view is still visible
    $min_far_dist /= 2;
    print "Min near = $min_near_dist,  min far = $min_far_dist\n";

    # Adjust the children distances
    if (scalar(@$near)) {
	foreach (@$near) { $length->{$_} = $dist->{$_} - $min_near_dist + 1e-9; }
    }
    if (scalar(@$far)) {
	foreach (@$far)  { $length->{$_} = $dist->{$_} - $min_far_dist + 1e-9; }
    }

    # FIXME: if near or far is size <= 1, don't make a new parent node.
    @children = ();
    if (scalar(@$near) > 0) {
	my $mp = "N/$merge_num";
	$length->{$mp} = $min_near_dist;
	$merge_num++;
	my @c = sort {$dist->{$a} <=> $dist->{$b}} @$near;
	foreach (@c) {
	    $parent->{$_} = $mp;
	}
	$tree->{$mp} = join(",", @c);
	$parent->{$mp} = $node;
	print "NEW NODE: tree{$mp} = $tree->{$mp}\n";
	push(@children, $mp);
    } else {
	# @$near nodes are directly children of $node rather than via their own parent
	# FIXME: make sure parents of @$near are $node
	push(@children, @$near);
    }


    # FIXME:
    # Instead of below find the best old node for the fork of FAR.  What was it closest to?
    # Make it a child of that instead?  Or a fork and child?
    # Ie distance nodes are still anchored in the same place they were before, relative to
    # what they were paired with.

    if (scalar(@$far) && 0) {
	my $mp = "F/$merge_num";
	$merge_num++;
	$length->{$mp} = $min_far_dist;
	my @c = sort {$dist->{$a} <=> $dist->{$b}} @$far;
	foreach (@c) {
	    $parent->{$_} = $mp;
	}
	$tree->{$mp} = join(",", @c);
	$parent->{$mp} = $node;
	print "NEW NODE: tree{$mp} = $tree->{$mp}\n";
	push(@children, $mp);
    } else {
	# Do this first.
	# Bubble up to find sibling that's in NEAR.
	# Then take this out of NEAR set and replace with the parent of NEAR and this FAR.
	# That parent then gets put into the new @childen, and we have no need to do
	# anything with FAR then as it's included already.
	foreach(@$far) {
	    print "FAR: $_ ",$parent->{$_},"\n";
	}
	push(@children, @$far);
    }


    $tree->{$node} = join(",", @children);

    foreach (@children) {
	$parent->{$_} = $node;
    }
}

sub print_tree_ordered {
    my ($tree, $length, $node, $total_dist) = @_;

    foreach (sort $tree->{$node}) {
	print "  "x$depth,"$node=($_) DIST: $length->{$node}\n";
	my @nodes = sort {$length->{$a} <=> $length->{$b}} split(",",$_);
	foreach (@nodes) {
	    if (defined(%{$tree}{$_})) {
		$depth++;
		print_tree_ordered($tree, $length, $_, $total_dist + $length->{$_});
		$depth--;
	    } else {
		print "  "x($depth+1), "CHILD: $_\tDIST: $length->{$_}  \t",
		    $length->{$_} + $total_dist, "\n";
	    }
	}
    }
}

# Flatten tree by tree walk.
#
# 1. Walk from any given root and accumulate total distance from that
#    root of child nodes as we walk.
# 2. Don't walk into child nodes > max_dist from root.
# 3. Aggregate all the nodes we visited in that walk as a single new
#    child.
# 4. Then continue the walk into children which are > max_dist away
#    so every node is visited.
sub flatten_tree_ordered {
    my ($tree, $parent, $length, $node, $total_dist) = @_;
    my ($near, $far, $dist) = nodes_within_dist($tree, $length, $node, 0, $max_dist);
    merge_nodes($tree, $parent, $length, $node, $near, $far, $dist) if (scalar(@{$near}));
    foreach (sort $tree->{$node}) {
	print "  "x$depth,"$node=($_) DIST: $length->{$node}\n";
	my @nodes = sort {$length->{$a} <=> $length->{$b}} split(",",$_);
	foreach (@nodes) {
	    if (defined(%{$tree}{$_})) {
		$depth++;
		flatten_tree_ordered($tree, $parent, $length, $_, $total_dist + $length->{$_});
		$depth--;
	    } else {
		print "  "x($depth+1), "CHILD: $_\tDIST: $length->{$_}  \t",
		    $length->{$_} + $total_dist, "\n";
	    }
	}
    }
}

#my $str = "((A1,A2),(A3,A4));";
#my $str = "((A1,A2,(A3,A4)),((B1,B2),(B3,B4)));";
#my $str = "(A,B,(C,D)E)F;";
#my $str = "(A,B,(C,D));";

my $str = "";
if ($ARGV[0] eq "-t") {
    $str = $ARGV[1];
} else {
    while (<>) {
	$str .= $_;
    }
}
$str =~ tr/\012//d;

my ($root, $tree, $parent, $leaf, $length) = parse_tree($str);

print "\n\n";
print_tree_ordered($tree, $length, $root, 0);
print "\n\n";
flatten_tree_ordered($tree, $parent, $length, $root, 0);
print "\n\n";
print_tree_ordered($tree, $length, $root, 0);

print "===\n";
print tree_to_newick($tree, $length, $root), "\n";
