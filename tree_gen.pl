#!/usr/bin/perl -w
use strict;

# Usage: tree_gen.pl width pchild seed

# Width is the average number of sequences per depth
# Pchild is the probability of node producing a child

my $width = shift(@ARGV);
my $pchild = shift(@ARGV);
my $seed = shift(@ARGV);

my $max_depth = 100;
my $max_seq = 500;

srand($seed);

# Newick format tree is:
# Tree      -> Subtree ";"
# Subtree   -> Leaf | Internal
# Leaf      -> Name
# Internal  -> "(" BranchSet ")" Name
# BranchSet -> Branch | Branch "," BranchSet
# Branch    -> Subtree Length
# Name      -> empty | string
# Length    -> empty | ":" number
#
# eg ((A:1,B:1):1,(C:1,D:2):2,E:0);
#
#     .---A
# .---|
# |   `---B
# E
# |       .---C
# `-------|
#         `-------D
# 0...1...2...3...4

# sequence number for naming
my $count = 0;

sub graph {
    my ($nseq,$depth) = @_;
    my $str = "(";
    for (my $i=0; $i<$nseq; $i++) {
	$str .= "," if $i>0;
	if (rand()<$pchild && $depth < $max_depth && $count < $max_seq) {
	    # subgraph instead of leaf
	    $str .= graph(int(rand(2*$width)+2), $depth+1);
	} else {
	    # leaf node
	    $str .= "S".$count++;
	}
	$str .= ":".rand()/1000;
    }
    $str .= ")";
    return $str;
}

print graph(int(rand($width)+2),0), ";\n";
print STDERR "Produced $count sequences\n";
