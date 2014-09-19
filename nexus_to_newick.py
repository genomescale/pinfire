#!/usr/bin/python2.7

import os
import sys
import subprocess

paup_binary = sys.argv[1]
nexus_paths = sys.argv[2:]

for nexus_path in nexus_paths:
	folder_path, file_name = os.path.split(nexus_path)
	base, ext = os.path.splitext(file_name)

	newick_file_name = base + ".newick"
	newick_file_path = os.path.join(folder_path, newick_file_name)

	paup_file_name = base + ".paup_input"
	paup_file_path = os.path.join(folder_path, paup_file_name)

	paup_block = """#NEXUS

begin paup;
  Set Maxtrees=100000;
  GetTrees rooted storeBrLens=yes file=%s;
  SaveTrees file=%s format=Newick brLens=user root=yes;
end;

quit;
""" % (file_name, newick_file_name)

	paup_file = open(paup_file_path, "w")
	paup_file.write(paup_block)
	paup_file.close()

	paup_cmd = [paup_binary, "-n", paup_file_path]

	try:
		paup_result = subprocess.check_output(paup_cmd)
		print(0)
		print(paup_result)
	except CalledProcessError as paup_error:
		print(paup_error.return_code)
		print(paup_error.output)
