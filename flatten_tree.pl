#!/usr/bin/perl -w
use strict;

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

sub tree_to_newick2 {
    my ($tree, $node) = @_;
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
	    push(@ele, "(" . tree_to_newick2($tree, $n) . ")$name");
	} else {
	    push(@ele, $n);
	}
    }

    return "@ele";
}

sub tree_to_newick {
    my ($tree, $node) = @_;
    return "(" . tree_to_newick2($tree, $node) . ");";
}

#my $str = "((A1,A2),(A3,A4));";
#my $str = "((A1,A2,(A3,A4)),((B1,B2),(B3,B4)));";
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
    my ($lineage) = $_ =~ /^(.)/;
    #my ($lineage) = $_ =~ /_.*_([^:]*)/;
    #print "$_ $lineage $parent->{$_}\n";
    my $node = $parent->{$_};
    my $x = 0;
    while ($node) {
	$lineages{$node}->{$lineage}++;
	#print "  "x++$x, "$node\n";
	$node = $parent->{$node};
    }
}

# Flatten tree descending from node $node into the parent
sub flatten {
    my ($node) = @_;

    return $node unless (defined($tree->{$node}));

    my @list = ();
    foreach my $n (split(",", $tree->{$node})) {
	if (defined($leaf->{$n})) {
	    push(@list, $n);
	    next;
	}

	foreach my $c (split(",", $tree->{$n})) {
	    push(@list, flatten($c));
	}
    }

    #print "Flatten $node [$tree->{$node}] => @list\n";
    return @list;
}

for (my $i=0;$i<1; $i++) {

print STDERR "### LOOP $i ###\n";

#foreach my $node (sort {$b cmp $a} keys(%{$parent})) {
foreach my $node (sort keys(%{$parent})) {
    next unless defined $lineages{$node};

    # Find lineage list for this node
    my $lref = $lineages{$node};
    my $nlin = scalar(keys(%$lref));
    print STDERR "\n$node $nlin\n";
    foreach my $lin (keys(%$lref)) {
	print STDERR "    $lref->{$lin} $lin\n";
    }

    # If it's a single lineage, look at all the siblings to check they
    # also have the same lineage and only one.
    my $p;
    if ($nlin == 1) {
	#print "$tree->{$node}\n";
	$p = $parent->{$node};
	#print "$p = ($tree->{$p})\n";

	my $ok = 1;
	foreach my $n (split(",", $tree->{$p})) {
	    next unless defined($lineages{$n});
	    my $cref = $lineages{$n};
	    my @kl = keys(%$lref);
	    my @kc = keys(%$cref);
	    print STDERR "  $n @kl @kc\n";
	    $ok = 0 if (@kl ne @kc);
	}
	next unless $ok;

	print STDERR "MERGE $p [$tree->{$p}]\n";
	my @l = flatten($p);
	print STDERR "New list = @l";

	# Do the merge and cull old
	$tree->{$p} = join(",",@l);
	foreach my $n (@l) {
	    undef($tree->{$n});
	    $parent->{$n} = $p;
	}
    }
    print STDERR "\n";
}
}; # 10 times

print tree_to_newick($tree, $root), "\n";
exit;

__END__

TODO: also need to collapse higher up nodes.

Eg a child E,F may have a mixture, but if other siblings in a cascade
A->(B,(C,(D,(E,F)))) are all the same lineage (A,B,C,D) then this could still
be collapsed.

This implies merging when the number of lineages doesn't change rather than is
equal and 1.  (Merge all bar the offending E,F which is still forked off?)


@ seq22-head2[quantum/phylo]; ./treedist.r _tmp/out.500.20/lineage.nwk _tmp/out.500.20/aligned.compat.
nwk
[1] 136
[1] 67.15732
[1] 6.564895
[1] 0.4919742


@ seq22-head2[quantum/phylo]130; ./treedist.r _tmp/out.500.20/lineage.nwk _tmp/out.500.20/aligned.fa.treefile
[1] 494
[1] 387.7573
[1] 60.93719
[1] 0.8096969


./flatten_tree.pl "`cat _tmp/out.500.20/aligned.fa.treefile`" 2>/dev/null > _tmp/out.500.20/aligned.fa.treefile2

@ seq22-head2[quantum/phylo]; ./treedist.r _tmp/out.500.20/lineage.nwk _tmp/out.500.20/aligned.fa.treefile2
[1] 372
[1] 222.6094
[1] 51.13566
[1] 0.7897409


So a big reduction, but still not close to compat.
It also proves the point: tree comparison isn't the right tactic.
