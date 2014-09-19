#!/usr/bin/python2.7

import dendropy
import os

root_folder = "split_loci"
simulation_folders = sorted(os.listdir(root_folder))
n_loci = 128

for simulation in simulation_folders:
	concatenation_matrix = {}
	for locus_n in range(n_loci):
		locus_filename = "tree_0_%d.nex" % (locus_n)
		locus_path = os.path.join(root_folder, simulation, locus_filename)
		locus_matrix = dendropy.DnaCharacterMatrix.get_from_path(locus_path, "nexus")
		locus_taxa = sorted(locus_matrix.taxon_set)
		for taxon in locus_taxa:
			dna_sequence = str(locus_matrix[taxon])

			species_name, individual = str(taxon).split()
			fixed_name = "s%02d" % (int(species_name[1:]))

			if individual == "tip0":
				if locus_n == 0:
					concatenation_matrix[fixed_name] = dna_sequence
				else:
					concatenation_matrix[fixed_name] += dna_sequence

	fasta_filename = "concatenation_%s.fasta" % (simulation)
	print("Writing " + fasta_filename)
	fasta_file = open(fasta_filename, "w")
	for taxon, dna_sequence in concatenation_matrix.items():
		fasta_file.write(">%s\n%s\n" % (taxon, dna_sequence))
	fasta_file.close()