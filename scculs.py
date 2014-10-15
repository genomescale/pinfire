PROGRAM_VERSION = "scculs.py, part of SCCULS preview-1"

import libscculs
import argparse
import os
import ete2

arg_parser = argparse.ArgumentParser(description = "Produce a summary tree from an MCMC sample of rooted, ultrametric trees.")
arg_parser.add_argument("-v", "--version", action = "version", version = PROGRAM_VERSION)

defaults_group = arg_parser.add_argument_group("program defaults")
defaults_group.add_argument("-c", "--candidate-method", type = str, default="derived", choices = ["sampled", "derived"], help = "Only consider topologies in the MCMC sample, or derive the most probable topology or topologies using conditional clades. Default: derived.")
defaults_group.add_argument("-e", "--node-heights", type = str, choices = ["median", "mean"], help = "Specify the method used to calculate node heights. If none is specified, node heights will not be calculated, and trees of equal branch lengths will be returned.")
defaults_group.add_argument("-p", "--probability-method", type = str, default="conditional-clade", choices = ["tree-topology", "conditional-clade"], help = "Calculate clade and tree topology probabilities using tree topology frequencies, or using conditional clade frequencies. Default: conditional-clade.")
defaults_group.add_argument("-u", "--support-values", type = str, default="both", choices = ["clade", "conditional-clade", "both"], help = "The type of support values added to summary tree nodes: clade monophyly probabilities, or conditional clade probabilities, or both. Default: both.")

output_group = arg_parser.add_argument_group('output files')
output_group.add_argument("-n", "--newick-output", metavar = "NEWICK_OUTPUT_PATH", type = str, help = "Output the summary tree(s) in newick format. When --max-topologies is greater than 1, more than one tree may be returned, so an identifying number will be appended to the end of each filename. This output and/or a CSV output must be supplied.")
output_group.add_argument("-o", "--csv-output", metavar = "CSV_OUTPUT_PATH", type = str, help = "Calculate statistics for each returned tree topology, and output them to a file in CSV format. This output and/or a newick output must be supplied.")

limits_group = arg_parser.add_argument_group('output limits')
limits_group.add_argument("-m", "--max-hpd", type = float, default = 1.0, help = "The size of the Highest Posterior Density (HPD) interval to return tree topologies from. The number of topologies returned will still be limited by -t/--max-topologies. Default: 1.0")
limits_group.add_argument("-t", "--max-topologies", type = int, default = 1, help = "How many tree topologies (in descending order of posterior probability) to output. The number of topologies returned will still be limited by -m/--max-hpd. Default: 1.")

input_group = arg_parser.add_argument_group('program input')
input_group.add_argument("-b", "--burn-in", type = int, default = 0, help = "The number of trees to discard from the beginning of the MCMC sample. Default: 0.")
input_group.add_argument("-d", "--calibration-date", type = float, default = 0.0, help = "If any tip dates are not contemporary (including tip date sampling), set a fixed date for the calibration taxon so that the tree height is correctly calculated. Default: 0.0.")
input_group.add_argument("-x", "--calibration-taxon", type = str, default = "", help = "If any tip dates are not contemporary (including tip date sampling), set the calibration taxon with a fixed date so that the tree height is correctly calculated. Negative numbers are used for past dates, positive numbers for future dates.")
input_group.add_argument("sample_path", metavar = "MCMC_SAMPLE_PATH", type = str, help = "The path to an MCMC sample of phylogenetic trees in either nexus or newick format.")

args = arg_parser.parse_args()

if not (args.newick_output or args.csv_output):
	arg_parser.error("specify output file paths using -n/--newick-output and -o/--csv-output")
elif args.max_topologies <= 0:
	arg_parser.error("argument -t/--max-topologies: must be equal to or greater than 1")
elif args.max_hpd <= 0.0 or args.max_hpd > 1.0:
	arg_parser.error("argument -m/--max-hpd: must be greater than 0.0 and less than 1.0")

sample_path = args.sample_paths
if not os.path.isfile(sample_path):
	arg_parser.error("argument MCMC_SAMPLE_PATH: not a file path")

derive_trees_thresh = args.max_topologies
derive_post_thresh = args.max_hpd
sample_burn_in = args.burn_in
calibration_taxon = args.calibration_taxon
calibration_date = args.calibration_date

print "Reading MCMC sample %s..." % (sample_path)
sample_newick = libscculs.read_newick(sample_path)[sample_burn_in:]
sample_ultrametric = libscculs.UltrametricSample(sample_newick, calibration_taxon, calibration_date)

print("Counting topologies and conditional clades...")
tpf, tpd, ccf, ccd = libscculs.count_topologies(sample_ultrametric)
txo = ultrametric_samples[0].taxon_order
print("Calculating topology probabilities...")
tpp = libscculs.calculate_discrete_probabilities(tpf)

print("Calculating conditional clade probabilities...")
ccp = {}
for cid in ccf:
	conditional_frequencies = ccf[cid]
	conditional_data = ccd[cid]
	ccp[cid] = libscculs.calculate_discrete_probabilities(conditional_frequencies, libscculs.cc_yule_log_prior, conditional_data)

print("Calculating topology probabilities from conditional clade probabilities...")
stp, std = libscculs.defined_topology_probabilities(ccp, ultrametric_samples)

print("Deriving probable topologies from conditional clades...")
dtp, dtd = libscculs.best_topology_probabilities(ccp, txo, derive_trees_thresh, derive_post_thresh)

print("Counting number of derived topologies")
print(libscculs.nonzero_derived_topologies(ccp, txo))

if args.newick_output:
	annotated_topologies = libscculs.add_derived_probabilities(dtd, txo, ccp)
	for t_hash, ns in annotated_topologies.items():
		t_prob = dtp[t_hash]
		print(ns)
