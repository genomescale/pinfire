#!/usr/bin/python2.7

import csv
import os
import sys

def read_time(time_output):
	time_output = time_output.replace("+", " ")
	time_output = time_output.replace("(", " ")
	time_output = time_output.replace(")", "")
	time_results = time_output.strip().split()

	result_structure = {}
	for result in time_results:
		result_size = len(result)
		key_i = 0
		for i in range(result_size):
			c = result[i]
			if c.isdigit():
				key_i = i + 1

		result_key = result[key_i:]
		result_value = result[:key_i]

		if result_value.isdigit():
			result_numeric = int(result_value)
		elif ":" in result_value:
			result_numeric = result_value
		else:
			result_numeric = float(result_value)

		result_structure[result_key] = result_numeric

	return result_structure

sim_results_paths = sys.argv[1:]
time_folder = "time"
ess_folder = "ess"

header = []
results = []
for results_path in sim_results_paths:
	sim_results_file = open(results_path)
	sim_results_reader = csv.reader(sim_results_file)

	n_columns = 0
	row_i = 0
	for row in sim_results_reader:
		row_i += 1
		if row_i == 1:
			header = row
			n_columns = len(row)
		else:
			row_structure = {}
			for column_i in range(n_columns):
				key = header[column_i]
				value = row[column_i]
				if value.isdigit():
					numeric = int(value)
				else:
					numeric = float(value)
				row_structure[key] = numeric

			results.append(row_structure)

new_header = header + ["ess_per_hour", "seconds_per_ess"]

combined_filename = "concat_combined.csv"
combined_file = open(combined_filename, "w")
combined_writer = csv.writer(combined_file)
combined_writer.writerow(new_header)

for result in results:
	i = result['individuals_per_species']
	r = result['replicate']
	s = result['n_species']

	base = "%02dr_%02ds_%02di" % (r, s, i)
	time_filename = base + ".xml.time"
	time_path = os.path.join(time_folder, time_filename)
	time_file = open(time_path)
	time_output = time_file.read()
	time_file.close()

	time_results = read_time(time_output)
	if "user" in time_results:
		cpu_time = time_results['user'] + time_results['system']

		ess_filename = "ess_" + base + ".csv"
		ess_path = os.path.join(ess_folder, ess_filename)
		ess_file = open(ess_path)
		ess_file.readline()
		ess_reader = csv.reader(ess_file, dialect = csv.excel_tab)

		treeheight_ess = 0.0
		for row in ess_reader:
			if len(row) > 0:
				statistic = row[0]
				if statistic == "TreeHeight":
					treeheight_ess = float(row[6])
					break

		ess_file.close()

		ess_per_hour = (treeheight_ess / cpu_time) * 3600.0
		seconds_per_ess = cpu_time / treeheight_ess

		new_row = []
		for h in header:
			new_row.append(result[h])

		new_row.append(ess_per_hour)
		new_row.append(seconds_per_ess)

		combined_writer.writerow(new_row)

combined_file.close()