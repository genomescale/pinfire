#!/usr/bin/python2.7

import libsapphire
import os
import csv
import sys

derive_threshold = 10.0 ** -5
burn_in = 7223

simulation_folder = "newick"

output_path = sys.argv[1]
output_file = open(output_path, "w")
output_writer = csv.writer(output_file)

#output_header  = ["replicate", "n_species", "individuals_per_species", "n_partitions", "loci_per_partition"]
#output_header += ["n_direct_topologies",  "total_direct_probability",  "direct_spread",  "direct_bias"]
#output_header += ["n_sampled_topologies", "total_sampled_probability", "sampled_spread", "sampled_bias"]
#output_header += ["n_derived_topologies", "total_derived_probability", "derived_spread", "derived_bias"]

#output_writer.writerow(output_header)

first_rep = int(sys.argv[2])
if len(sys.argv) > 3:
	last_rep = int(sys.argv[3])
else:
	last_rep = first_rep

true_replicates = range(first_rep, last_rep + 1)

loci_exponents = range(1, 8)
max_le = max(loci_exponents)

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
					partitions_exponents = range(0, 1 + max_le - l_e)
					for k_e in partitions_exponents:
						l = 2**l_e
						k = 2**k_e

						if k <= 16:
							sample_paths = []
							for p in range(k):
								sample_filename = "species_%s_%02dk_%03dl_%02dp.newick" % (simulation, k, l, p)
								sample_path = os.path.join(simulation_folder, sample_filename)
								sample_paths.append(sample_path)

							parameters = [r, s, i, k, l]
							output = libsapphire.hpd_sets(truth_path, sample_paths, derive_threshold, burn_in)
							#output_writer.writerow(parameters + output)
							#output_file.flush()

output_file.close()