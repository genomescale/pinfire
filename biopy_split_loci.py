#!/usr/bin/python2.7
import os
import subprocess
import starbeast
import sys
import copy

simulation_folder = "split_loci"

def divide_labels(labels, segments, members):
	n_labels     = len(labels)
	label_range  = range(n_labels)
	segment_size = n_labels / segments
	select_every = segment_size / members

	if n_labels % segments != 0:
		"invalid number of segments"
		sys.exit()

	if segment_size % members != 0:
		"invalid number of members"
		sys.exit()

	partitions = []
	for p in range(segments):
		p_start = segment_size * (p)
		p_end   = segment_size * (p + 1)

		partition = []
		for n in label_range:
			if (n % select_every == 0) and (n >= p_start) and (n < p_end):
				partition.append(labels[n])

		partitions.append(partition)

	return partitions

def read_nexus(nexus_path):
	nexus_file = open(nexus_path)

	reading = False
	sequences = {}
	for r in nexus_file.readlines():
		r_strip = r.strip()
		if r_strip == ";":
			reading = False

		if reading:
			r_split = r_strip.split()
			label = r_split[0]
			sequence = r_split[1]

			label_split = label.split("_")
			species_n = int(label_split[0][1:])
			tip_n     = int(label_split[1][3:])
			species = "s%02d" % (species_n)
			tip     = "t%02d" % (tip_n)

			if species in sequences:
				sequences[species][tip] = sequence
			else:
				sequences[species] = {tip: sequence}

		if r_strip == "Matrix":
			reading = True

	return sequences

true_replicates = range(0, 10)

loci_nt = 200

loci_exponents = range(1, 8)
max_le = max(loci_exponents)
max_loci = 2**max_le

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
				sim_folder = os.path.join(simulation_folder, simulation)
				#os.makedirs(sim_folder)

				print("Simulating " + simulation)
				st_cmd = "sample_species_tree -n 1 -b 2 -d 0.4 -p l,1,0.3 %d -o %s_species.trees" % (s, simulation)
				gt_cmd = "genetree_in_sptree_sim -n %d -t %d --nexus=%s_genes.trees %s_species.trees" % (max_loci, i, simulation, simulation)
				sq_cmd = "simulate_sequences -m HKY,0.01,4 -n %d %s_genes.trees" % (loci_nt, simulation)

				#subprocess.call(st_cmd, cwd=sim_folder, shell=True)
				#subprocess.call(gt_cmd, cwd=sim_folder, shell=True)
				#subprocess.call(sq_cmd, cwd=sim_folder, shell=True)

				print("Reading    " + simulation)
				sequences = {}
				for locus_n in range(max_loci):
					locus = "locus%03d" % (locus_n)

					matrix_filename = "tree_0_%d.nex" % (locus_n)
					matrix_path = os.path.join(sim_folder, matrix_filename)

					sequences[locus] = read_nexus(matrix_path)

				locus_labels   = sorted(sequences.keys())
				species_labels = sorted(sequences[locus_labels[0]].keys())

				print("Writing    " + simulation)
				for l_e in loci_exponents:
					partitions_exponents = range(0, 1 + max_le - l_e)
					for k_e in partitions_exponents:
						l = 2**l_e
						k = 2**k_e
						if l > 16 or k > 8:
							locus_partitions = divide_labels(locus_labels, k, l)

							for p in range(k):
								loci = locus_partitions[p]

								simulation_sequences = {}
								for locus in loci:
									simulation_sequences[locus] = {}
									for species in species_labels:
										simulation_sequences[locus][species] = copy.copy(sequences[locus][species])

								beast_run = "%s_%02dk_%03dl_%02dp" % (simulation, k, l, p)
								print("Writing    " + beast_run)
								output_filename = beast_run + ".xml"
								output_path = os.path.join(sim_folder, output_filename)

								output_file = open(output_path, "w")
								starbeast_xml = starbeast.make_xml(beast_run, simulation_sequences)
								starbeast_xml.write(output_file, encoding="UTF-8", xml_declaration=True)
								output_file.close()