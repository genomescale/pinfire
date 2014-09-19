#!/usr/bin/python2.7
# coding=utf8

import csv
import os
import numpy
import matplotlib.pyplot as plt
import math
import sys

graphs_folder = "graphs"

nan = float("nan")

pantone = {
"red": 0xED2939,
"orange": 0xFF5800,
"yellow": 0xFEDF00,
"green": 0x00AD83,
"cyan": 0x00A6CB,
"blue": 0x0018A8,
"magenta": 0xD0417E}

marker_styles = {1: "o", 2: "s", 4: "d", 8: "^", 16: "v", 32: "<", 64: ">"}
plt_color_order = {
2: "#%06x" % (pantone["red"]),
4: "#%06x" % (pantone["orange"]),
8: "#%06x" % (pantone["yellow"]),
16: "#%06x" % (pantone["green"]),
32: "#%06x" % (pantone["cyan"]),
64: "#%06x" % (pantone["blue"]),
128: "#%06x" % (pantone["magenta"])}

result_fn = sys.argv[1]
result_file = open(result_fn)
result_reader = csv.reader(result_file)

header = []
results = {}
n_columns = 0
row_i = 0
x_variables = []
y_variables = []
for row in result_reader:
	if row_i == 0:
		header = row
		x_variables = row[-2:]
		y_variables = row[5:-2]
		n_columns = len(row)
	else:
		row_integers = {}
		row_floats = {}
		for col_i in range(n_columns):
			key = header[col_i]
			value = row[col_i]
			if col_i < 5:
				row_integers[key] = int(value)
			else:
				row_floats[key] = float(value)

		graph_id = "%02dx%02d" % (row_integers["n_species"], row_integers["individuals_per_species"])
		n_partitions = row_integers["n_partitions"]
		n_loci = n_partitions * row_integers["loci_per_partition"]

		if graph_id not in results:
			results[graph_id] = {}

		if n_loci not in results[graph_id]:
			results[graph_id][n_loci] = {}

		if n_partitions not in results[graph_id][n_loci]:
			results[graph_id][n_loci][n_partitions] = {}
			for key in row_floats:
				results[graph_id][n_loci][n_partitions][key] = []

		for key, value in row_floats.items():
			if not math.isnan(value):
				results[graph_id][n_loci][n_partitions][key].append(value)

	row_i += 1

def labelize(var_name):
	with_spaces = var_name.replace("_", " ")
	with_capital = with_spaces.capitalize()
	with_acronym = with_capital.replace("ess", "ESS")
	with_hash = with_acronym.replace(" n ", " #")

	return with_hash

for graph_id in results:
	for x_var in x_variables:
		for y_var in y_variables:
			fig, ax = plt.subplots()
			max_y = 0.0

			for n_loci in results[graph_id]:
				partition_order = sorted(results[graph_id][n_loci].keys())
				x_means  = [numpy.mean(results[graph_id][n_loci][np][x_var]) for np in partition_order]
				y_means  = [numpy.mean(results[graph_id][n_loci][np][y_var]) for np in partition_order]
				max_y    = max([max_y] + y_means)

				x_sterr  = [numpy.std(results[graph_id][n_loci][np][x_var]) / numpy.sqrt(len(results[graph_id][n_loci][np][x_var])) for np in partition_order]
				y_sterr  = [numpy.std(results[graph_id][n_loci][np][y_var]) / numpy.sqrt(len(results[graph_id][n_loci][np][y_var])) for np in partition_order]

				ax.errorbar(x_means, y_means, c = plt_color_order[n_loci], xerr = x_sterr, yerr = y_sterr, label = str(n_loci))

				for i in range(len(partition_order)):
					n_partitions = partition_order[i]
					loci_per_partition = n_loci / n_partitions
					ax.scatter(x_means[i], y_means[i], s = 50.0, c = plt_color_order[n_loci], marker = marker_styles[n_partitions], linewidth = '0')
					ax.scatter(x_means[i], y_means[i], s = 300.0, c = "k", marker = u"$%d×%d$" % (n_partitions, loci_per_partition), linewidth = '0', zorder = 100)

			ax.set_xscale('log')
			if max_y > 10.0:
				ax.set_yscale('log')

			plt.xlabel(labelize(x_var), fontsize=16)
			plt.ylabel(labelize(y_var), fontsize=16)

			figure_filename = "%s_%s_over_%s.png" % (graph_id, y_var, x_var)
			figure_path = os.path.join(graphs_folder, figure_filename)
			fig.savefig(figure_path, dpi = 300)

			plt.close()
