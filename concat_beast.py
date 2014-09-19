#!/usr/bin/python2.7

import standard_beast
import os

def read_fasta(fasta_path):
	fasta_file = open(fasta_path)

	label = ""
	sequence = ""
	sequences = {}
	l = fasta_file.readline()
	while l != "":
		if l[0] == ">":
			if label != "":
				sequences[label] = sequence.upper()

			label = l[1:].rstrip()
			sequence = ""
		else:
			sequence += l.rstrip()

		l = fasta_file.readline()

	sequences[label] = sequence.upper()

	fasta_file.close()

	return sequences

output_folder = "concat_xml"
concat_folder = "concat_seq"
concat_filenames = os.listdir(concat_folder)
for filename in concat_filenames:
	simulation = filename[filename.find("_") + 1:filename.rfind(".")]
	concat_path = os.path.join(concat_folder, filename)
	concat_sequences = read_fasta(concat_path)

	beast_xml = standard_beast.make_xml(simulation, concat_sequences)

	output_filename = simulation + ".xml"
	output_path = os.path.join(output_folder, output_filename)
	output_file = open(output_path, "w")
	beast_xml.write(output_file, encoding="UTF-8", xml_declaration=True)
	output_file.close()