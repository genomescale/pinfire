PROGRAM_VERSION = "scculs.py, part of SCCULS preview-1"

import libscculs
import argparse
import os

arg_parser = argparse.ArgumentParser(description = "SCCULS: Scalable Conditional-Clade Ultrametric Summary trees. Produces an annotate summary tree from an MCMC sample, using the highest tree probability calculated from conditional clade frequencies, or from topology frequencies.")
arg_parser.add_argument("-v", "--version", action = "version", version = PROGRAM_VERSION)

defaults_group = arg_parser.add_argument_group("program defaults")
defaults_group.add_argument("-c", "--candidate-method", type = str, default = "derived", choices = ["sampled", "derived"], help = "Only consider topologies in the MCMC sample, or derive the most probable topology or topologies using conditional clades. Default: derived.")
defaults_group.add_argument("-e", "--node-heights", type = str, choices = ["median", "mean"], help = "Specify the method used to calculate node heights. Without this option, node heights will not be calculated, and trees of equal branch lengths will be returned.")
defaults_group.add_argument("-p", "--probability-method", type = str, choices = ["tree-topology", "conditional-clade"], help = "When candidate method is 'sampled', infer clade and tree topology probabilities using either tree topology frequencies or conditional clade frequencies. Default: conditional-clade.")
defaults_group.add_argument("-u", "--support-values", type = str, default = "clade", choices = ["clade", "conditional-clade", "both"], help = "The type of support values added to summary tree nodes: clade monophyly probabilities, or conditional clade probabilities, or both. Default: both.")

output_group = arg_parser.add_argument_group('output files')
output_group.add_argument("-n", "--newick-output", metavar = "NEWICK_OUTPUT_PATH", type = str, help = "Output the summary tree(s) in newick format. When --max-topologies is greater than 1, more than one tree may be returned, so an identifying number will be appended to the end of each filename. This output and/or a CSV output must be supplied.")
output_group.add_argument("-o", "--csv-output", metavar = "CSV_OUTPUT_PATH", type = str, help = "Calculate statistics for each returned tree topology, and output them to a file in CSV format. This output and/or a newick output must be supplied.")

limits_group = arg_parser.add_argument_group('output limits')
limits_group.add_argument("-m", "--max-probability", type = float, default = 1.0, help = "The size of the credible set in total posterior probability to output. The number of topologies returned will still be limited by -t/--max-topologies. Default: 1.0")
limits_group.add_argument("-t", "--max-topologies", type = int, default = 1, help = "The size of the credible set in the number of unique topologies to output. The number of topologies returned will still be limited by -m/--max-probability. Default: 1.")

input_group = arg_parser.add_argument_group('program input')
input_group.add_argument("-b", "--burn-in", type = int, default = 0, help = "The number of trees to discard from the beginning of the MCMC sample. Default: 0.")
input_group.add_argument("-d", "--calibration-date", type = float, default = 0.0, help = "If any tip dates are not contemporary (including tip date sampling), set a fixed date for the calibration taxon so that the tree height is correctly calculated. Negative numbers are used for past dates, positive numbers for future dates. Default: 0.0.")
input_group.add_argument("-x", "--calibration-taxon", type = str, default = "", help = "If any tip dates are not contemporary (including tip date sampling), set the calibration taxon so that the tree height is correctly calculated.")
input_group.add_argument("sample_path", metavar = "MCMC_SAMPLE_PATH", type = str, help = "The path to an MCMC sample of phylogenetic trees in either nexus or newick format.")

args = arg_parser.parse_args()

# raise errors in response to incomplete or nonsensical user-supplied arguments
if not (args.newick_output or args.csv_output):
	arg_parser.error("specify output file paths using -n/--newick-output and -o/--csv-output")
elif args.probability_method is None and args.candidate_method == "sampled":
	arg_parser.error("argument -p/--probability-method is required when -c/--candidate-method is 'sampled'")
elif args.probability_method is not None and args.candidate_method == "derived":
	arg_parser.error("argument -p/--probability-method must not be used when -c/--candidate-method is 'derived'")
elif args.max_topologies <= 0:
	arg_parser.error("argument -t/--max-topologies: must be equal to or greater than 1")
elif args.max_probability <= 0.0 or args.max_probability > 1.0:
	arg_parser.error("argument -m/--max-probability: must be greater than 0.0 and less than 1.0")
elif not os.path.isfile(args.sample_path):
	arg_parser.error("argument MCMC_SAMPLE_PATH: not a file path")

print "Reading MCMC sample..."
calibration_taxon = args.calibration_taxon
calibration_date = args.calibration_date
sample_path = args.sample_path
sample_burn_in = args.burn_in

mcmc_sample = libscculs.trees_from_path(sample_path)
mcmc_post = mcmc_sample[sample_burn_in:] # discard burn-in
ultrametric_sample = libscculs.UltrametricSample(mcmc_post, calibration_taxon, calibration_date)
taxon_order = ultrametric_sample.taxon_order

print("Counting topologies and conditional clades...")
topology_set, topology_counts, cc_sets, cc_counts = libscculs.calculate_topology_probabilities(ultrametric_sample)

max_tree_topologies = args.max_topologies
max_probability = args.max_probability

if args.probability_method == "conditional-clade" or args.candidate_method == "derived":
	print("Calculating conditional clade probabilities...")
	for parent_hash, split_counts in cc_counts.items():
		cc_sets[parent_hash].probabilities_from_counts(split_counts)

if args.candidate_method == "sampled":
	if args.probability_method == "tree-topology":
		print("Calculating topology probabilities...")
		topology_set.probabilities_from_counts(topology_counts)
	else:
		print("Calculating topology probabilities from conditional clade probabilities...")
		topology_set.probabilities_from_ccs(cc_sets)

	topology_set.cull_probabilities(max_tree_topologies, max_probability)
	output_topology_set = topology_set

#else:
#	print("Deriving probable topologies from conditional clades...")
#	output_topology_set = libscculs.derive_best_topologies(cc_sets, max_tree_topologies, max_probability)

if args.newick_output is not None:
	for i in range(output_topology_set.n_features):
		ns = output_topology_set.data_array[i]
		prob = output_topology_set.probabilities_array[i]
		print("%g %s" % (prob, ns))