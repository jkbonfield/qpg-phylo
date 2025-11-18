#!/usr/bin/perl -w
use strict;

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

    my $node_num=0;

    # Remove trailing semicolon
    $tree =~ s/;$//;

    my @pos;
    do {
	#print "=== $tree ===\n";
	@pos = ();

	$_ = $tree;
	# Cache positions of a leaf node branchSet
	while (/\(([^()]+)\)([^,:()]*)((?::[0-9.e+-]*)?)/gc) {
	    #print substr($_,$-[0], $+[0]-$-[0]),"\n";
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

my $str = "";
if ($ARGV[0] eq "-t") {
    $str = $ARGV[1];
} else {
    while (<>) {
	$str .= $_;
    }
}
$str =~ tr/\012//d;
$str =~ s/\)[^(),:;]*/)/g;

my ($root, $tree, $parent, $leaf) = parse_tree($str);

#print_tree($tree, $root);
#print "\n\n\n";
#print_traceback($parent, "B3");

my %lineages = ();
sub lineage {
    ($_) = @_;
    #/^(.)/;

    # Leaf nodes first, and then sub-trees
    /_.*_([^:]*)/;
    return $1 if $1;

    if (defined($lineages{$_})) {
	my @l = sort keys(%{$lineages{$_}});
	#return "@l";
	return "~$l[0]";
    }
    return $_;
}

# Now bubble up to the root from every leaf node, accumulating lineage stats
# per node.
my %leaf_lineage;
foreach (sort keys(%{$leaf})) {
    my ($lineage) = lineage($_);
    $leaf_lineage{$lineage}++;
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
my $nlineages = scalar(keys %leaf_lineage);

# Find leftmost to start
my $node = $root;
while (exists($tree->{$node})) {
    print "node=$node\n";
    $_ = $tree->{$node};
    #my @nodes = sort(split(",", $_));
    my @nodes = sort {lineage($a) cmp lineage($b)} split(",", $_);
    print "children = @nodes\n";
    $node = $nodes[0];
}
print "=== $node",lineage($node),"\n";

# Walk tree.
# TODO: repeat walk by first / last node (left / right)?
#
# Just counting runs of matching nodes works as a score, but it gives an
# optimal score when we don't have a tree at all (all in one flat level).
# We also need to penalise for siblings differing.

my $last = "";
my $lineage_switch = 0;
my $switch_sibling = 0;
do {
    my $sibling = 0;
    for (;;) {
	my $p = $parent->{$node};
	if (!$p) {
	    $node = undef;
	    last;
	}
	print "p=$p\n";
	#my @nodes = sort (split(",",$tree->{$p}));
	my @nodes = sort {lineage($a) cmp lineage($b)} (split(",",$tree->{$p}));
	my $i = 0;
	for (; $i <= $#nodes; $i++) {
	    last if $nodes[$i] eq $node;
	}
	print ">>>Node $node, parent $p, $i of <$#nodes, lineages";
	foreach my $l (@nodes) {
	    print " ",lineage($l);
	}
	print "\n";

	# Next right sibling
	if ($i < $#nodes) {
	    $node = $nodes[++$i];
	    print "Sibling $node\n";
	    # If sibling is a tree, then recurse via first child
	    $sibling = defined($leaf->{$node}) ? 1 : 0;
	    while (!defined($leaf->{$node})) {
		#@nodes = sort(split(",",$tree->{$node}));
		@nodes = sort {lineage($a) cmp lineage($b)} split(",", $tree->{$node});
		$node = $nodes[0];
	    }
	    last;
	} else {
	    # Or else up one level
	    $node = $p;
	}
    }

    # NB just counting the number of lineage-runs is a good score.
    # ./flatten_tree2.pl "`cat file.nwk`" |awk '/===/ {print $3}'|uniq|wc -l
    print "=== $node ",lineage($node)," $sibling\n" if $node;
    if ($node && lineage($node) ne $last) {
	$lineage_switch++;
	$last = lineage($node);
	$switch_sibling++ if ($sibling);
    }
} while ($node);

my $score = $lineage_switch - ($nlineages-1) + $switch_sibling;
print "### SCORE $score\t$lineage_switch switches of $nlineages, $switch_sibling in siblings\n";

__END__

@ seq22-head2[quantum/phylo]; ./flatten_tree2.pl "`cat _tmp/out-xz.500.1/lineage.nwk`"|grep '###'
### SCORE 0	116 switches of 117, 0 in siblings


@ seq22-head2[quantum/phylo]; ./flatten_tree2.pl "`cat _tmp/out-xz.500.1/maple_tree.tree`"|grep '###'
### SCORE 52	150 switches of 117, 18 in siblings
@ seq22-head2[quantum/phylo]130; ./treedist.r _tmp/out-xz.500.1/lineage.nwk _tmp/out-xz.500.1/maple_tree.tree
[1] 504
[1] 396.3003


@ seq22-head2[quantum/phylo]; ./flatten_tree2.pl "`cat _tmp/out-xz.500.1/aligned.fa.treefile`"|grep '###'
### SCORE 57	148 switches of 117, 25 in siblings
@ seq22-head2[quantum/phylo]134; ./treedist.r _tmp/out-xz.500.1/lineage.nwk _tmp/out-xz.500.1/aligned.fa.treefile
[1] 500
[1] 389.724


@ seq22-head2[quantum/phylo]; ./flatten_tree2.pl "`cat _tmp/out-xz.500.1/aligned.compat.nwk`"|grep '###'
### SCORE 79	137 switches of 117, 58 in siblings
@ seq22-head2[quantum/phylo]134; ./treedist.r _tmp/out-xz.500.1/lineage.nwk _tmp/out-xz.500.1/aligned.compat.nwk
[1] 144
[1] 60.18528
[1] 8.283076
[1] 0.2546317


Ie costs:
	   RF    JRF    FT2(this)
compat	   144	  60	79
iqtree3	   500   389	57
maple	   504	 396	52
