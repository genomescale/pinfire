PROGRAM_VERSION = "scculs.py, part of SCCULS preview-1"

import libscculs
import argparse
import os
import csv

def safe_open(file_path, overwrite):
	if (overwrite == False) and os.path.exists(file_path):
		if os.path.isfile:
			raise Exception("This file already exists: " + file_path)
		elif os.path.isfolder:
			raise Exception("Please specify full paths including file names. This path is a folder: " + file_path)
		else:
			raise Exception()
	
	safe_open_file = open(file_path, "w")

	return safe_open_file

arg_parser = argparse.ArgumentParser(description = "SCCULS: Scalable Conditional-Clade Ultrametric Summary trees. Distills a summary tree from an MCMC sample, using the highest tree probability calculated from conditional clade frequencies, or from topology frequencies.")
arg_parser.add_argument("-v", "--version", action = "version", version = PROGRAM_VERSION)

defaults_group = arg_parser.add_argument_group("program defaults")
defaults_group.add_argument("-c", "--candidate-method", type = str, default = "derived", choices = ["derived", "sampled"], help = "Only consider topologies in the MCMC sample, or derive the most probable topology or topologies using conditional clades. Default: derived.")
defaults_group.add_argument("-g", "--node-heights", type = str, choices = ["median", "mean"], help = "Specify the method used to calculate node heights. Without this option, node heights will not be calculated, and trees of equal branch lengths will be returned.")
defaults_group.add_argument("-p", "--probability-method", type = str, choices = ["conditional-clade", "tree-topology"], help = "Infer tree topology probabilities using either tree topology probabilities or conditional clade probabilities. When -c/--candidate-method is 'derived', default is conditional-clade. When -c/--candidate-method is 'sampled', default is tree-topology.")
defaults_group.add_argument("-s", "--no-support-values", action = "store_true", help = "Do not calculate or add clade monophyly support values to the summary tree(s).")

output_group = arg_parser.add_argument_group('output files')
output_group.add_argument("-i", "--info-output", metavar = "INFO_OUTPUT_PATH", type = str, help = "Calculate whole-sample statistics and output them to a text format file.")
output_group.add_argument("-n", "--newick-output", metavar = "NEWICK_OUTPUT_PATH", type = str, help = "Output the summary tree(s) to newick format file(s). When -l/--max-topologies is greater than 1, more than one tree may be returned, so an identifying number will be appended to the end of each filename.")
output_group.add_argument("-o", "--csv-output", metavar = "CSV_OUTPUT_PATH", type = str, help = "Calculate statistics for each returned tree topology, and output them to CSV format file.")
output_group.add_argument("-w", "--overwrite", action = "store_true", help = "If output file paths point to existing files, overwrite the existing files.")

limits_group = arg_parser.add_argument_group('output limits')
limits_group.add_argument("-l", "--max-topologies", type = int, default = 1, help = "The size of the credible set in the number of unique topologies to output. The number of topologies returned will still be limited by -m/--max-probability. Default: 1.")
limits_group.add_argument("-m", "--max-probability", type = float, default = 1.0, help = "The size of the credible set in total posterior probability to output. The number of topologies returned will still be limited by -l/--max-topologies. Default: 1.0")

input_group = arg_parser.add_argument_group('program input')
input_group.add_argument("-b", "--burn-in", type = int, default = 0, help = "The number of trees to discard from the beginning of the MCMC sample. Default: 0.")
input_group.add_argument("-d", "--calibration-date", type = float, default = 0.0, help = "If any tip dates are not contemporary (including tip date sampling), set a fixed date for the calibration taxon so that the tree height is correctly calculated. Negative numbers are used for past dates, positive numbers for future dates. Default: 0.0.")
input_group.add_argument("-t", "--calibration-taxon", type = str, default = "", help = "If any tip dates are not contemporary (including tip date sampling), set the calibration taxon so that the tree height is correctly calculated.")
input_group.add_argument("sample_path", metavar = "MCMC_SAMPLE_PATH", type = str, help = "The path to an MCMC sample of phylogenetic trees in either nexus or newick format.")

args = arg_parser.parse_args()

# raise errors in response to incomplete or nonsensical user-supplied arguments
if args.max_topologies <= 0:
	arg_parser.error("argument -l/--max-topologies: must be equal to or greater than 1")
elif args.max_probability <= 0.0 or args.max_probability > 1.0:
	arg_parser.error("argument -m/--max-probability: must be greater than 0.0 and less than 1.0")
elif not os.path.isfile(args.sample_path):
	arg_parser.error("argument MCMC_SAMPLE_PATH: not a file path")

# set probability method
if args.probability_method is None: # defaults if not supplied
	if args.candidate_method == "derived":
		probability_method = "conditional-clade"
	else:
		probability_method = "tree-topology"
else: # user-supplied method
	probability_method = args.probability_method

calibration_taxon = args.calibration_taxon
calibration_date = args.calibration_date
sample_path = args.sample_path
sample_burn_in = args.burn_in
max_tree_topologies = args.max_topologies
max_probability = args.max_probability
overwrite = args.overwrite

print "Reading MCMC sample..."
mcmc_sample = libscculs.trees_from_path(sample_path)
mcmc_post = mcmc_sample[sample_burn_in:] # discard burn-in
ultrametric_sample = libscculs.UltrametricSample(mcmc_post, calibration_taxon, calibration_date)
taxon_order = ultrametric_sample.taxon_order
n_taxa = len(taxon_order)

print("Counting topologies and conditional clades...")
topology_set, topology_counts, cc_sets, cc_counts, clade_set = libscculs.calculate_topology_probabilities(ultrametric_sample)
n_unique_topologies = topology_set.n_features

# all circumstances where conditional clade probabilities are required
# don't bother to calculate if not needed
if (args.candidate_method == "derived") or (probability_method == "conditional-clade") or (args.support_values == "conditional-clade"):
	print("Calculating conditional clade probabilities...")
	for parent_hash, split_counts in cc_counts.items():
		cc_sets[parent_hash].probabilities_from_counts(split_counts)

if args.candidate_method == "derived": # derive credible topologies from conditional clades
	print("Deriving probable topologies from conditional clades...")
	output_topology_set = libscculs.derive_best_topologies(cc_sets, taxon_order, max_tree_topologies, max_probability)
else: # base credible topologies on frequency in MCMC sample
	output_topology_set = topology_set

if probability_method == "conditional-clade":
	print("Calculating topology probabilities from conditional clade probabilities...")
	output_topology_set.probabilities_from_ccs(cc_sets)
	if not args.no_support_values:
		print("Calculating clade probabilities from conditional clade probabilities...")
		clade_set.derive_clade_probabilities(cc_sets, n_taxa)
else:
	print("Calculating topology probabilities...")
	output_topology_set.probabilities_from_counts(topology_counts)
	if not args.no_support_values:
		print("Calculating topology and clade probabilities from MCMC sample...")
		topology_set.probabilities_from_counts(topology_counts)
		clade_set.melt_clade_probabilities(topology_set, n_taxa)

# once probabilities have been calculated for each topology in the sampled set
# then topologies that exceed maximum topology/probability limits can be removed
if args.candidate_method == "sampled":
	print("Limiting output topologies to credible set...")
	output_topology_set.cull_probabilities(max_tree_topologies, max_probability)

if not args.no_support_values:
	print("Adding clade support values to tree topologies...")
	output_topology_set.add_clade_support(clade_set, taxon_order)

if args.info_output is not None:
	print("Writing MCMC sample statistics file...")
	info_output_path = args.info_output
	info_output_file = safe_open(info_output_path, overwrite)
	info_output_file.write("Number of taxa in each tree: %i\n" % (n_taxa))
	info_output_file.write("Number of unique tree topologies in MCMC sample: %i\n" % (n_unique_topologies))

	if args.candidate_method == "derived": # calculate summary statistics for topologies
		n_nonzero_topologies = libscculs.n_derived_topologies(cc_sets, n_taxa)
		info_output_file.write("Number of topologies derived from conditional clades: %i\n" % (n_nonzero_topologies))

		#n_derived_topologies = libscculs.n_derived_topologies(cc_sets, n_taxa, include_zero_probability = True)
		#n_nonzero_topologies = libscculs.n_derived_topologies(cc_sets, n_taxa)
		#info_output_file.write("Number of topologies derived from conditional clades: %i\n" % (n_derived_topologies))
		#info_output_file.write("Number of topologies derived from conditional clades (with non-zero probabilities): %i\n" % (n_nonzero_topologies))

	info_output_file.close()

if args.newick_output is not None:
	print("Writing tree topology files...")
	newick_path_prefix = args.newick_output
	for i in range(output_topology_set.n_features):
		newick_string = output_topology_set.data_array[i]
		newick_output_path = newick_path_prefix + "." + str(i)
		newick_output_file = safe_open(newick_output_path, overwrite)
		newick_output_file.write(newick_string + "\n")
		newick_output_file.close()

if args.csv_output is not None:
	print("Writing tree statistics file...")
	csv_output_path = args.csv_output
	csv_output_file = safe_open(csv_output_path, overwrite)
	csv_writer = csv.writer(csv_output_file)

	header_row = ["topology", "probability"]
	csv_writer.writerow(header_row)

	for i in range(output_topology_set.n_features):
		topology_probability = output_topology_set.probabilities_array[i]
		output_row = [i, topology_probability]
		csv_writer.writerow(output_row)

	csv_output_file.close()
