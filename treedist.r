#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

#install.packages('TreeDist')
library('TreeDist')
t1 <- ape::read.tree(file = args[1])
t2 <- ape::read.tree(file = args[2])
#ExpectedVariation(t1, t2)
RobinsonFoulds(t1,t2)
JaccardRobinsonFoulds(t1,t2)
ClusteringInfoDistance(t1,t2)
TreeDistance(t1,t2)

