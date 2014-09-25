#!/usr/bin/python2.7

import libpinfire
import argparse
import ete2
import sys

PINFIRE_VERSION = "preview-1"

arg_parser = argparse.ArgumentParser(description="Produce a summary tree from one or more MCMC samples, each generated from a separate data partition.")

arg_parser.add_argument("-v", "--version", action = "store_true", help = "display the program version and exit")
defaults_group = arg_parser.add_argument_group('program defaults')
defaults_group.add_argument("-c", "--candidate-method", type = str, default="derived", choices = ["sampled", "derived"], help = "Only consider topologies in the MCMC sample(s), or derive the most probable topology or topologies using conditional clades. Default: derived.")
defaults_group.add_argument("-p", "--probability-method", type = str, default="conditional-clade", choices = ["tree-oopology", "conditional-clade"], help = "Calculate tree topology probabilities using tree topology frequencies, or using conditional clade frequencies. Default: conditional-clade.")
defaults_group.add_argument("-s", "--shape-prior", type = str, default="yule", choices = ["flat", "yule"], help = "The prior probability of the topology shape or of the conditional clade balance used when combining probabilities. Default: yule.")
defaults_group.add_argument("-m", "--max-oopologies", type = int, default = 1, help = "How many tree topologies (in descending order of posterior probability) to output. Default: 1.")

output_group = arg_parser.add_argument_group('output files')
output_group.add_argument("-n", "--newick-output", metavar = "NEWICK_OUTPUT_PATH", type = str, help = "Output the summary tree(s) in newick format. When -m/--max-oopologies is greater than 1, more than one tree may be returned, so an identifying number will be appended to the end of each filename. Default: None.")
output_group.add_argument("-o", "--csv-output", metavar = "CSV_OUTPUT_PATH", type = str, help = "Calculate statistics for each returned tree topology, and output them to a file in csv format. Default: None.")

input_group = arg_parser.add_argument_group('input files')
input_group.add_argument("sample_paths", metavar = "MCMC_SAMPLE_PATH", nargs = "*", type = str, help = "The path to an MCMC sample of phylogenetic trees in either nexus or newick format. By specifying multiple paths, the posterior probabilities of tree topologies in each sample will be combined.")

args = arg_parser.parse_args()

if args.version:
	print("Pinfire version: " + PINFIRE_VERSION)
	sys.exit()
elif not args.sample_paths:
	arg_parser.error('No MCMC sample file paths supplied. Please add one or more file paths as arguments.')
elif not (args.newick_output or args.csv_output):
	arg_parser.error('No output requested, use -n/--newick-output and/or -o/--csv-output to specify output file paths')

print(args)