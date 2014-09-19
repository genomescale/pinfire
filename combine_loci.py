#!/usr/bin/python2.7

import libpinfire
import os
import csv
import sys
import ete2
import subprocess
import numpy
import socket

nan = float("nan")

def accuracy(truth_path, sample_paths, derive_threshold, bhv_resolution, bhv_length, burn_in):
	print("Reading true tree...")
	truth_newick = libpinfire.read_newick(truth_path)[0]
	true_topology_newick = reset_tree_branch_lengths(truth_newick)

	ultrametric_samples = []
	for sample_path in sample_paths:
		print "Reading ultrametric sample %s..." % (sample_path)
		sample_newick = libpinfire.read_newick(sample_path)[burn_in:]
		sample_ultrametric = libpinfire.UltrametricSample(sample_newick)
		ultrametric_samples.append(sample_ultrametric)

	print("Counting topologies and conditional clades...")
	tpf, tpd, ccf, ccd = libpinfire.count_topologies(ultrametric_samples)
	txo = ultrametric_samples[0].taxon_order
	print("Calculating topology probabilities...")
	tpp = libpinfire.calculate_discrete_probabilities(tpf, libpinfire.yule_log_prior, tpd)

	print("Calculating conditional clade probabilities...")
	ccp = {}
	for cid in ccf:
		conditional_frequencies = ccf[cid]
		conditional_data = ccd[cid]
		ccp[cid] = libpinfire.calculate_discrete_probabilities(conditional_frequencies, libpinfire.cc_yule_log_prior, conditional_data)

	print("Calculating topology probabilities from conditional clade probabilities...")
	stp, std = libpinfire.sampled_topology_probabilities(ccp, ultrametric_samples)

	print("Deriving probable topologies from conditional clades...")
	dtp, dtd = libpinfire.derive_topology_probabilities(ccp, txo, derive_threshold)

	total_tpp = sum(tpp.values())
	total_stp = sum(stp.values())
	total_dtp = sum(dtp.values())

	n_tpp = len(numpy.flatnonzero(numpy.array(tpp.values())))
	n_stp = len(numpy.flatnonzero(numpy.array(stp.values())))
	n_dtp = len(numpy.flatnonzero(numpy.array(dtp.values())))

	if n_tpp > 0:
		print("Calculating precision and bias of sampled, directly inferred topology probabilities...")
		tp_sd, tp_mean_distance = precision_and_bias(tpp, tpd, true_topology_newick, bhv_resolution, bhv_length)
		print("Calculating minimum HPD interval & clade credibility of true tree using sampled, directly inferred topology probabilities...")
		t_prob_rank, t_prob_value, t_cc_rank, t_cc_value = truth_prob_and_cred(tpp, tpd, truth_newick, txo)
	else:
		tp_sd, tp_mean_distance, t_prob_rank, t_prob_value, t_cc_rank, t_cc_value = (nan, nan, nan, nan, nan, nan)

	if n_stp > 0:
		print("Calculating precision and bias of sampled, conditional-clade derived topology probabilities...")
		st_sd, st_mean_distance = precision_and_bias(stp, std, true_topology_newick, bhv_resolution, bhv_length)
		print("Calculating HPD interval & clade credibility of true tree using sampled, conditional-clade derived topology probabilities...")
		s_prob_rank, s_prob_value, s_cc_rank, s_cc_value = truth_prob_and_cred(stp, std, truth_newick, txo)
	else:
		st_sd, st_mean_distance, s_prob_rank, s_prob_value, s_cc_rank, s_cc_value = (nan, nan, nan, nan, nan, nan)

	if n_dtp > 0:
		print("Calculating precision and bias of conditional-clade derived topology probabilities...")
		dt_sd, dt_mean_distance = precision_and_bias(dtp, dtd, true_topology_newick, bhv_resolution, bhv_length)
		print("Calculating HPD interval & clade credibility of true tree using conditional-clade derived topology probabilities...")
		d_prob_rank, d_prob_value, d_cc_rank, d_cc_value = truth_prob_and_cred(dtp, dtd, truth_newick, txo)
	else:
		dt_sd, dt_mean_distance, d_prob_rank, d_prob_value, d_cc_rank, d_cc_value = (nan, nan, nan, nan, nan, nan)

	results  = [n_tpp, total_tpp, tp_sd, tp_mean_distance, t_prob_rank, t_prob_value, t_cc_rank, t_cc_value]
	results += [n_stp, total_stp, st_sd, st_mean_distance, s_prob_rank, s_prob_value, s_cc_rank, s_cc_value]
	results += [n_dtp, total_dtp, dt_sd, dt_mean_distance, d_prob_rank, d_prob_value, d_cc_rank, d_cc_value]

	return results

def precision_and_bias(topology_probabilities, topology_strings, true_topology_newick, bhv_resolution, bhv_length):
	equal_branch_length_trees = reset_branch_lengths(topology_strings)
	topology_random_sample = randomly_sample_topologies(topology_probabilities, equal_branch_length_trees, bhv_resolution)
	n_unique_trees = len(set(topology_random_sample))

	if n_unique_trees == 1:
		sample_mean_newick = topology_random_sample[0]
		sample_sd = 0.0
	else:
		n_iterations, sample_mean_newick, sample_sd = bhv_central_tendancy(topology_random_sample, bhv_length)
		print("Converged on mean tree after %d iterations" % (n_iterations))

	sample_mean_distance = bhv_distance(true_topology_newick, sample_mean_newick)

	return sample_sd, sample_mean_distance

def truth_prob_and_cred(topology_probabilities, topology_strings, truth_newick, taxon_order):
	prob_ranks = libpinfire.rank_discrete(topology_probabilities)

	truth_array = libpinfire.generate_tree_array(truth_newick, taxon_order)
	truth_hash = truth_array["f0"].tostring()

	if truth_hash in prob_ranks:
		clade_credibilities = libpinfire.calculate_clade_credibility(topology_probabilities, topology_strings, topology_strings, taxon_order)

		cred_ranks = libpinfire.rank_discrete(clade_credibilities)

		truth_prob_rank  = prob_ranks[truth_hash]
		truth_prob_value = topology_probabilities[truth_hash]

		truth_cc_rank  = cred_ranks[truth_hash]
		truth_cc_value = clade_credibilities[truth_hash]
	else:
		truth_strings = {truth_hash: truth_newick}
		clade_credibilities = libpinfire.calculate_clade_credibility(topology_probabilities, topology_strings, truth_strings, taxon_order)

		truth_prob_rank  = nan
		truth_prob_value = 0.0

		truth_cc_rank  = nan
		truth_cc_value = clade_credibilities[truth_hash]

	return truth_prob_rank, truth_prob_value, truth_cc_rank, truth_cc_value

def reset_branch_lengths(all_trees, new_branch_length = 1.0):
	reset_trees = {}
	for tree_id, newick_string in all_trees.items():
		reset_tree_string = reset_tree_branch_lengths(newick_string, new_branch_length)
		reset_trees[tree_id] = reset_tree_string

	return reset_trees

def reset_tree_branch_lengths(newick_string, new_branch_length = 1.0):
	tree = ete2.Tree(newick_string)
	for node in tree.get_descendants():
		node.dist = new_branch_length

	reset_tree_string = tree.write(format = 5)

	return reset_tree_string

def randomly_sample_topologies(tree_frequencies, tree_strings, resolution):
	ordered_hashes = sorted(tree_frequencies.keys())

	ordered_frequencies = []
	ordered_trees = []
	for tree_hash in ordered_hashes:
		frequency = tree_frequencies[tree_hash]
		newick_string = tree_strings[tree_hash]

		ordered_frequencies.append(frequency)
		ordered_trees.append(newick_string)

	scaling_factor = sum(ordered_frequencies)
	scaled_frequencies = [f / scaling_factor for f in ordered_frequencies]

	sampled_trees = numpy.random.choice(ordered_trees, size = resolution, replace = True, p = scaled_frequencies)

	return sampled_trees

def bhv_central_tendancy(tree_strings, cauchy_length):
	sm_jar_path = "SturmMean/SturmMean.jar"

	sample_suffix = unique_file_suffix()
	sample_filename = "weighted-sample%s.newick" % (sample_suffix)
	sample_path = os.path.join(temp_folder, sample_filename)
	sample_file = open(sample_path, "w")

	for ts in tree_strings:
		sample_file.write(ts + "\n")

	sample_file.flush()
	sample_file.close()

	sm_cmd = ["java", "-Xmx4096m", "-Xms512m", "-jar", sm_jar_path, "-c", str(cauchy_length), sample_path]
	try:
		sm_output = subprocess.check_output(sm_cmd)
	except subprocess.CalledProcessError as error:
		print("SturmMean error:")
		print(error.output)
		sys.exit(1)

	sm_lines = sm_output.split("\n")
	convergence_line = sm_lines[9]
	n_iterations = int(convergence_line.strip().split()[2])

	sd_line = sm_lines[10]
	sd = float(sd_line[sd_line.find(":") + 2:sd_line.find("  ")])
	mean_tree = ete2.Tree(sm_lines[11])
	mean_newick = mean_tree.write(format = 5)

	return n_iterations, mean_newick, sd

def bhv_distance(first_tree, second_tree):
	gtp_jar_path = "GTP/gtp.jar"

	distance_suffix = unique_file_suffix()
	distance_filename = "bhv-pair%s.newick" % (distance_suffix)
	distance_path = os.path.join(temp_folder, distance_filename)
	distance_file = open(distance_path, "w")

	distance_file.write(first_tree + "\n")
	distance_file.write(second_tree + "\n")

	distance_file.flush()
	distance_file.close()

	gtp_cmd = ["java", "-Xmx4096m", "-Xms512m", "-jar", gtp_jar_path, "-v", distance_path]
	try:
		gtp_output = subprocess.check_output(gtp_cmd)
	except subprocess.CalledProcessError as error:
		print("GTP error:")
		print(error.output)
		sys.exit(1)

	gtp_lines = gtp_output.split("\n")
	gd_line = gtp_lines[-3]
	geodesic_distance = float(gd_line.strip().split()[-1])

	return geodesic_distance

def unique_file_suffix():
	process_id = os.getpid()

	host_domain_name = socket.gethostname()
	if "." in host_domain_name:
		host_name = host_domain_name[:host_domain_name.find(".")]
	else:
		host_name = host_domain_name

	suffix = ".%s.%d" % (host_name, process_id)

	return suffix

derive_threshold = 10.0 ** -6
bhv_resolution = 1000
bhv_length = 25
burn_in = 7223

simulation_folder = "newick"
temp_folder = "."

output_path = sys.argv[1]
output_file = open(output_path, "w")
output_writer = csv.writer(output_file)

output_header  = ["replicate", "n_species", "individuals_per_species", "n_partitions", "loci_per_partition"]
output_header += [ "direct_n_topologies",  "direct_total_probability",  "direct_spread",  "direct_bias",  "direct_truth_prob_rank",  "direct_truth_probability",  "direct_truth_cred_rank",  "direct_truth_credibility"]
output_header += ["sampled_n_topologies", "sampled_total_probability", "sampled_spread", "sampled_bias", "sampled_truth_prob_rank", "sampled_truth_probability", "sampled_truth_cred_rank", "sampled_truth_credibility"]
output_header += ["derived_n_topologies", "derived_total_probability", "derived_spread", "derived_bias", "derived_truth_prob_rank", "derived_truth_probability", "derived_truth_cred_rank", "derived_truth_credibility"]

output_writer.writerow(output_header)

first_rep = int(sys.argv[2])
if len(sys.argv) > 3:
	last_rep = int(sys.argv[3])
else:
	last_rep = first_rep

true_replicates = range(first_rep, last_rep + 1)

loci_exponents = range(1, 8)
partitions_exponents = range(0, 5)
max_loci = 128

max_tips_exp = 6
species_exp_range = range(3, 6)
individuals_exp_range = range(1, 4)

for r in true_replicates:
	for s_exp in species_exp_range:
		for i_exp in individuals_exp_range:
			if (s_exp + i_exp) <= max_tips_exp:
				s = 2**s_exp
				i = 2**i_exp

				simulation = "%02dr_%02ds_%02di" % (r, s, i)
				truth_filename = simulation + "_species.newick"
				truth_path = os.path.join(simulation_folder, truth_filename)

				for l_e in loci_exponents:
					for k_e in partitions_exponents:
						l = 2**l_e
						k = 2**k_e

						if (l * k) <= max_loci:
							sample_paths = []
							for p in range(k):
								sample_filename = "species_%s_%02dk_%03dl_%02dp.newick" % (simulation, k, l, p)
								sample_path = os.path.join(simulation_folder, sample_filename)
								sample_paths.append(sample_path)

							parameters = [r, s, i, k, l]
							output = accuracy(truth_path, sample_paths, derive_threshold, bhv_resolution, bhv_length, burn_in)
							output_writer.writerow(parameters + output)
							output_file.flush()

output_file.close()