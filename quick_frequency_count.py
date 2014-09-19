#!/usr/bin/python2.7

import sys

newick_path = sys.argv[1]
newick_file = open(newick_path)
newick_text = newick_file.read().strip()

newick_file.close()

newick_strings = newick_text.split("\n")
unique_strings = set(newick_strings)

for ns in unique_strings:
	print newick_strings.count(ns)