#!/usr/bin/python2.7

import ete2
import sys

def fudge_newick(newick_string):
	fudged_tree = ete2.Tree(newick_string)
	for node in fudged_tree.get_descendants():
		if node.is_leaf():
			taxon_number = int(node.name[1:])
			node.name = "s%02d" % (taxon_number)
	fudged_string = fudged_tree.write()
	return fudged_string

fudge_paths = sys.argv[1:]

for path in fudge_paths:
	newick_file = open(path)
	newick = newick_file.read().strip()
	newick_file.close()
	fudged = fudge_newick(newick)
	newick_file = open(path, "w")
	newick_file.write(fudged + "\n")
	newick_file.close()