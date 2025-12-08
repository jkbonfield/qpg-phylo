#!/usr/bin/perl -w
use strict;

#TODO:
# Score (A,A,A,B,B,C) as 3 for C and B,B vs A,A,A rather than 2.
# This means flattened costs more than non-flat.
# See MG8_eg1.nwk vs MG8_eg2.nwk for example of why this would help.
# We may even want an extra multiplier on in-node switches vs walk-switches
# to penalise not attempting to form a tree more.

# As per flatten_tree2.pl but sort nodes when recursing by branch length
# instead of by lineage name.

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

	    $node = "/$node_num" unless $node ne "";
	    $node = "#$node";
	    $node_num++;

	    #print "tree{$node} = \"$list\"\n";
	    foreach (split(",", $list)) {
		m/(^[^:]*)(?::(.*))?/;
		$_=$1 if (defined($1));
		#print "LLL $node $length $_ $1 // $2 //\n";
		$parent{$_} = $node;
		$leaf{$_} = 1 if (!/^#/);
		if (defined($2)) {
		    $length{$_} = $2;
		} else {
		    $length{$_} = 0;
		}
	    }
	    $list =~ s/:[0-9.e+-]*//g;
	    $tree{$node}=$list;

	    # Replace branchSet with an internal node
	    substr($tree,$st,$len) = $node;
	    print "Tree is now <$tree>\n";
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

my ($root, $tree, $parent, $leaf, $length) = parse_tree($str);

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
	# Alphabetical is best
	my @l = sort keys(%{$lineages{$_}});

	# However by frequency naively feels like it should be better.
#	my %l = %{$lineages{$_}};
#	my @l = sort {$l{$b} <=> $l{$a}} sort {$a cmp $b} keys %l;

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

sub hash {
    ($_)=@_;
    my $h=0;
    foreach (split("",$_)) {
	$h=(($h*31 + ord($_))) & 0xffff;
    }
    return $h;
}

sub linlength {
    ($_) = @_;

    # Cluster by length and resolve ties by lineage.
    my $lineage = lineage($_);
    $lineage = 1-(length($lineage) + hash($lineage) / 0xffff)*0.1;

    if (defined($length->{$_})) {
	return $length->{$_} + $lineage;
    } else {
	return $lineage;
    }
}

# Find leftmost to start
my $node = $root;
while (exists($tree->{$node})) {
    print "node=$node\n";
    $_ = $tree->{$node};
    #my @nodes = sort(split(",", $_));
    my @nodes = sort {linlength($a) <=> linlength($b)} split(",", $_);
    print "children = @nodes\n";
    $node = $nodes[0];
}
print "=== $node ",lineage($node),"\n";

# Walk tree.
# TODO: repeat walk by first / last node (left / right)?
#
# Just counting runs of matching nodes works as a score, but it gives an
# optimal score when we don't have a tree at all (all in one flat level).
# We also need to penalise for siblings differing.

my $last = "";
my $lineage_switch = 0;
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
	my @nodes = sort {linlength($a) <=> linlength($b)} (split(",",$tree->{$p}));
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
		@nodes = sort {linlength($a) <=> linlength($b)} split(",", $tree->{$node});
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
	print "SWITCH sibling=$sibling: $last ",lineage($node),"\n";
	$lineage_switch++ if $last;
	$last = lineage($node);
    }
} while ($node);

$depth = 0;
sub sibling_switches {
    my ($tree, $node) = @_;
    my $switches = 0;
    foreach (sort $tree->{$node}) {
	print "  "x$depth,"$node=($_)\n";
	my @nodes = split(",",$_);
	my %linfreq;
	foreach (@nodes) {
	    $linfreq{lineage($_)}++ if lineage($_) !~ /^~/;
	}
	my @linkeys = sort {$linfreq{$b} <=> $linfreq{$a}} keys(%linfreq);
	foreach (@linkeys) {
	    print "LINKEY $_ = $linfreq{$_}\n";
	}
	foreach (sort {lineage($a) cmp lineage($b)} @nodes) {
	    if (defined(%{$tree}{$_})) {
		$depth++;
		$switches += sibling_switches($tree,$_);
		$depth--;
	    } else {
		my $lin = lineage($_);
		print "SIBLING $node // $lin $linkeys[0]\n"
		    if $lin ne $linkeys[0] && $lin !~ /~/;
		$switches++ if $lin ne $linkeys[0] && $lin !~ /~/;
		print "   "x($depth+1),"CHILD $lin\n";
	    }
	}
    }

    return $switches;
}

print "ZZZ\n"; # For ... |sed -n '/ZZZ/,$p'
my $switches = sibling_switches($tree, $root);
print "Switches = $switches\n";

my $score = $lineage_switch - ($nlineages-1) + $switches;
print "### SCORE $score\t$lineage_switch switches of $nlineages, $switches in siblings\n";
