#!/usr/bin/python2.7

PROGRAM_VERSION = "summary_trees.py, part of 'pinfire' preview 1"

import libpinfire
import argparse

arg_parser = argparse.ArgumentParser(description = "Produce a summary tree from one or more MCMC samples, each generated from a separate data partition.")
arg_parser.add_argument("-v", "--version", action = "version", version = PROGRAM_VERSION)

defaults_group = arg_parser.add_argument_group("program defaults")
defaults_group.add_argument("-c", "--candidate-method", type = str, default="derived", choices = ["sampled", "derived"], help = "Only consider topologies in the MCMC sample(s), or derive the most probable topology or topologies using conditional clades. Default: derived.")
defaults_group.add_argument("-p", "--probability-method", type = str, default="conditional-clade", choices = ["tree-topology", "conditional-clade"], help = "Calculate clade and tree topology probabilities using tree topology frequencies, or using conditional clade frequencies. Default: conditional-clade.")
defaults_group.add_argument("-s", "--shape-prior", type = str, default="yule", choices = ["flat", "yule"], help = "The prior probability on topology/conditional clade shape used when combining probabilities. Default: yule.")
defaults_group.add_argument("-u", "--support-values", type = str, default="both", choices = ["clade", "conditional-clade", "both"], help = "The prior probability on the topology shape or of the conditional clade balance used when combining probabilities. Default: both.")

output_group = arg_parser.add_argument_group('output files')
output_group.add_argument("-n", "--newick-output", metavar = "NEWICK_OUTPUT_PATH", type = str, help = "Output the summary tree(s) in newick format. When --max-topologies is greater than 1, more than one tree may be returned, so an identifying number will be appended to the end of each filename. This output and/or a CSV output must be supplied.")
output_group.add_argument("-o", "--csv-output", metavar = "CSV_OUTPUT_PATH", type = str, help = "Calculate statistics for each returned tree topology, and output them to a file in CSV format. This output and/or a newick output must be supplied.")

limits_group = arg_parser.add_argument_group('output limits')
limits_group.add_argument("-m", "--max-hpd", type = float, default = 1.0, help = "The size of the Highest Posterior Density (HPD) interval to return tree topologies from. The number of topologies returned will still be limited by -t/--max-topologies. Default: 1.0")
limits_group.add_argument("-t", "--max-topologies", type = int, default = 1, help = "How many tree topologies (in descending order of posterior probability) to output. The number of topologies returned will still be limited by -m/--max-hpd. Default: 1.")

input_group = arg_parser.add_argument_group('input files')
input_group.add_argument("-b", "--burn-in", type = str, default = "0", help = "A single integer or a comma-separated list of integers, specifying the number of trees to discard from the beginning of each MCMC sample. Default: 0.")
input_group.add_argument("sample_paths", metavar = "MCMC_SAMPLE_PATH", nargs = "+", type = str, help = "The path to an MCMC sample of phylogenetic trees in either nexus or newick format. By specifying multiple paths, the posterior probabilities of tree topologies in each sample will be combined.")

args = arg_parser.parse_args()

if not (args.newick_output or args.csv_output):
	arg_parser.error("specify output file paths using -n/--newick-output and -o/--csv-output")
elif args.max_topologies <= 0:
	arg_parser.error("argument -t/--max-topologies: must be equal to or greater than 1")
elif args.max_hpd <= 0.0 or args.max_hpd > 1.0:
	arg_parser.error("argument -m/--max-hpd: must be greater than 0.0 and less than 1.0")

n_samples = len(args.sample_paths)
n_burn_in = []
if "," in args.burn_in:
	burn_in = args.burn_in.split(",")
	if len(burn_in) != n_samples:
		arg_parser.error("argument -b/--burn-in: number of integers != number of input files")
else:
	burn_in = [args.burn_in] * n_samples

for i in range(n_samples):
	if burn_in[i].isdigit():
		burn_in_int = int(burn_in[i])
		n_burn_in.append(burn_in_int)
	else:
		arg_parser.error("argument -b/--burn-in: not a non-negative integer")

print(args)
print(n_burn_in)