import numpy
import ete2
import math

order_differs_error = Exception("The set of taxa differs between tree samples, so they may not be directly compared")

class UltrametricSample():
	def __init__(self, newick_strings):
		self.taxon_order = []
		self.newick_strings = []
		self.tree_arrays = []
		self.feature_counts = {}

		n_trees = len(newick_strings)

		self.newick_strings = newick_strings
		self.root_heights = numpy.zeros(n_trees, dtype = numpy.float64)

		for i in range(n_trees):
			ns = self.newick_strings[i]
			tree = ete2.Tree(ns)
			if i == 0:
				taxa = tree.get_leaf_names()
				self.taxon_order = sorted(taxa)

			tree_array = generate_tree_array(ns, self.taxon_order)

			self.tree_arrays.append(tree_array)

class DiscreteFrequencies():
	feature_struct_format = ""
	feature_frequencies = numpy.array(0)

	def __init__(self, frequencies):
		frequency_labels = frequencies.keys()
		frequency_labels.sort()
		n_features = len(frequency_labels)

		first_feature = frequency_labels[0]
		n_partitions = 0
		max_label_size = 0

		frequencies_table = []
		for i in range(n_features):
			feature = frequency_labels[i]
			counts = frequencies[feature]
			feature_frequencies = tuple([feature] + counts)
			frequencies_table.append(feature_frequencies)

			if i == 0:
				n_partitions = len(counts)

			if len(feature) > max_label_size:
				max_label_size = len(feature)

		column_dtypes = ["a" + str(max_label_size)] + ["u8"] * n_partitions
		self.feature_struct_format = ",".join(column_dtypes)
		self.feature_frequencies = numpy.array(frequencies_table, self.feature_struct_format)

def sampled_topology_probabilities(cc_probabilities, tree_samples):
	n_partitions = len(tree_samples)
	topology_probabilities = {}
	topology_strings = {}

	first_taxon_order = []
	for i in range(n_partitions):
		ts = tree_samples[i]
		if i == 0:
			first_taxon_order = ts.taxon_order[:]
		elif first_taxon_order != ts.taxon_order:
			raise order_differs_error

		n_trees = len(ts.tree_arrays)
		for j in range(n_trees):
			tree_array = ts.tree_arrays[j]
			newick_string = ts.newick_strings[j]
			topology_hash = tree_array["f0"].tostring()

			if topology_hash not in topology_probabilities:
				split_probabilities = []

				topology_tree = ete2.Tree(newick_string)
				topology_string = topology_tree.write(format = 9) # strip branch lengths
				topology_strings[topology_hash] = topology_string

				for node in tree_array:
					parent_id = node[0]
					split_id = node[1]

					n_node_taxa = clade_size(parent_id)
					if n_node_taxa > 1:
						split_probability = cc_probabilities[parent_id][split_id]
						split_probabilities.append(split_probability)

				topology_probabilities[topology_hash] = numpy.product(split_probabilities)

	return topology_probabilities, topology_strings

def elucidate_cc_split(parent_id, split_id):
	parent_id_bytes = numpy.array(tuple(parent_id)).view(dtype = numpy.uint8)
	split_id_bytes = numpy.array(tuple(split_id)).view(dtype = numpy.uint8)

	parent_id_bits = numpy.unpackbits(parent_id_bytes)
	split_id_bits = numpy.unpackbits(split_id_bytes)

	n_parent_bits = len(parent_id_bits)
	n_split_bits = len(split_id_bits)

	child1_bits = numpy.zeros(n_parent_bits, dtype = numpy.uint8)
	child2_bits = numpy.zeros(n_parent_bits, dtype = numpy.uint8)

	j = 0
	for i in range(n_parent_bits):
		if parent_id_bits[i] == 1:
			if j < n_split_bits:
				if split_id_bits[j] == 1:
					child1_bits[i] = 1
				else:
					child2_bits[i] = 1
			else:
				child2_bits[i] = 1

			j += 1

	child1_bytes = numpy.packbits(child1_bits)
	child2_bytes = numpy.packbits(child2_bits)

	child1_id = child1_bytes.tostring().rstrip("\x00") # urgh C (null terminated strings)
	child2_id = child2_bytes.tostring().rstrip("\x00") # vs Python (not null terminated) strings

	return child1_id, child2_id

def derive_topology_probabilities(cc_probabilities, taxon_order, probability_threshold = 0.0):
	n_taxa = len(taxon_order)
	n_bytes, remainder = divmod(n_taxa, 8)
	if remainder > 0:
		n_bytes += 1

	derived_struct_format = "a%d,a%d,u1,f8" % (n_bytes, n_bytes)

	root_hash_bits = numpy.ones(n_taxa, dtype = numpy.uint8)
	root_hash_bytes = numpy.packbits(root_hash_bits)
	root_hash = root_hash_bytes.tostring()

	star_tree = numpy.array([(root_hash, "", 1, 1.0)], dtype=derived_struct_format)
	derived_topologies = [star_tree]

	all_resolved = False
	first_step = True
	while not all_resolved:
		n_derived_topologies = len(derived_topologies)
		if first_step:
			first_step = False
		else:
			print "Derived %d probable topologies..." % (n_derived_topologies)

		all_resolved = True
		for topology_i in reversed(range(n_derived_topologies)):
			topology = derived_topologies[topology_i]
			unresolved_nodes = numpy.flatnonzero(topology["f2"])
			topology_probability = numpy.prod(topology["f3"])
			if topology_probability > probability_threshold:
				if len(unresolved_nodes) > 0:
					all_resolved = False

					unresolved_index = unresolved_nodes[0]
					unresolved_hash = topology[unresolved_index]["f0"].tostring()
					split_probabilities = cc_probabilities[unresolved_hash]

					for split_hash in split_probabilities:
						split_probability = split_probabilities[split_hash]
						if split_probability > 0.0:
							child1_hash, child2_hash = elucidate_cc_split(unresolved_hash, split_hash)

							if clade_size(child1_hash) > 1: # unresolved (more than one taxon)
								child1_row = numpy.array([(child1_hash, "", 1, 1.0)], dtype=derived_struct_format)
							else: # resolved (leaf/tip node)
								child1_row = numpy.array([(child1_hash, "", 0, 1.0)], dtype=derived_struct_format)

							if clade_size(child2_hash) > 1: # unresolved (more than one taxon)
								child2_row = numpy.array([(child2_hash, "", 1, 1.0)], dtype=derived_struct_format)
							else: # resolved (leaf/tip node)
								child2_row = numpy.array([(child2_hash, "", 0, 1.0)], dtype=derived_struct_format)

							new_topology = numpy.concatenate((topology, child1_row, child2_row))

							# record information to resolve node
							new_topology[unresolved_index]["f1"] = split_hash
							new_topology[unresolved_index]["f2"] = 0
							new_topology[unresolved_index]["f3"] = split_probability

							derived_topologies.append(new_topology)

					derived_topologies.pop(topology_i)

			else:
				derived_topologies.pop(topology_i)

	print("Computing newick strings and tree probabilities for derived topologies...")
	derived_topology_probabilities = {}
	derived_topology_newick = {}
	for dt in derived_topologies:
		dt_hash = numpy.sort(dt["f0"]).tostring()
		dt_probability = numpy.prod(dt["f3"])

		splits = {}
		for node in dt:
			parent_id = node["f0"]
			split_id = node["f1"]
			splits[parent_id] = split_id

		tree_model = ete2.Tree()
		derive_tree_from_splits(tree_model, root_hash, taxon_order, splits)
		dt_newick = tree_model.write(format = 9)

		derived_topology_probabilities[dt_hash] = dt_probability
		derived_topology_newick[dt_hash] = dt_newick

	return derived_topology_probabilities, derived_topology_newick

def derive_tree_from_splits(current_node, parent_hash, taxon_order, splits):
	split_hash = splits[parent_hash]
	child1_hash, child2_hash = elucidate_cc_split(parent_hash, split_hash)

	child1_node = ete2.Tree()
	child2_node = ete2.Tree()

	current_node.add_child(child1_node)
	current_node.add_child(child2_node)

	if clade_size(child1_hash) == 1:
		child1_node.name = clade_taxon_names(child1_hash, taxon_order)[0]
	else:
		derive_tree_from_splits(child1_node, child1_hash, taxon_order, splits)

	if clade_size(child2_hash) == 1:
		child2_node.name = clade_taxon_names(child2_hash, taxon_order)[0]
	else:
		derive_tree_from_splits(child2_node, child2_hash, taxon_order, splits)

def clade_size(clade_hash):
	clade_node_bytes = numpy.array(tuple(clade_hash)).view(dtype = numpy.uint8)
	clade_node_bits = numpy.unpackbits(clade_node_bytes)
	n_clade_taxa = sum(clade_node_bits)

	return n_clade_taxa

def clade_taxon_names(clade_hash, taxon_order):
	taxon_names = []

	clade_node_bytes = numpy.array(tuple(clade_hash)).view(dtype = numpy.uint8)
	clade_node_bits = numpy.unpackbits(clade_node_bytes)
	n_bits = len(clade_node_bits)

	for i in range(n_bits):
		if clade_node_bits[i] == 1:
			taxon_names.append(taxon_order[i])

	return taxon_names

def count_topologies(tree_samples):
	n_partitions = len(tree_samples)
	topology_counts = {}
	topology_strings = {}
	cc_counts = {}
	cc_id_pairs = {}

	first_taxon_order = []
	for i in range(n_partitions):
		ts = tree_samples[i]
		if i == 0:
			first_taxon_order = ts.taxon_order[:]
		elif first_taxon_order != ts.taxon_order:
			raise order_differs_error

		n_trees = len(ts.tree_arrays)
		for j in range(n_trees):
			tree_array = ts.tree_arrays[j]
			newick_string = ts.newick_strings[j]
			topology_hash = tree_array["f0"].tostring()

			if topology_hash not in topology_counts:
				topology_tree = ete2.Tree(newick_string)
				topology_string = topology_tree.write(format = 9) # strip branch lengths

				topology_strings[topology_hash] = topology_string
				topology_counts[topology_hash] = [0] * n_partitions

			topology_counts[topology_hash][i] += 1

			for node in tree_array:
				parent_id = node[0]
				split_id = node[1]

				n_node_taxa = clade_size(parent_id)

				if n_node_taxa > 1:
					if parent_id not in cc_counts:
						cc_id_pairs[parent_id] = {split_id: (parent_id, split_id)}
						cc_counts[parent_id] = {split_id: [0] * n_partitions}
					elif split_id not in cc_counts[parent_id]:
						cc_id_pairs[parent_id][split_id] = (parent_id, split_id)
						cc_counts[parent_id][split_id] = [0] * n_partitions

					cc_counts[parent_id][split_id][i] += 1

	topology_frequencies = DiscreteFrequencies(topology_counts)
	cc_frequencies = {}

	for parent_id in cc_counts:
		conditional_count = cc_counts[parent_id]
		cc_frequencies[parent_id] = DiscreteFrequencies(conditional_count)

	return topology_frequencies, topology_strings, cc_frequencies, cc_id_pairs

def count_vertices(this_node, all_nodes):
	children = this_node.get_children()
	left_child = children[0]
	right_child = children[1]
	n = 1
	if not left_child.is_leaf():
		n += count_vertices(left_child, all_nodes)
	if not right_child.is_leaf():
		n += count_vertices(right_child, all_nodes)
	all_nodes.append(n)
	return n

def yule_log_prior(newick_string):
	tree = ete2.Tree(newick_string)
	f_vertices = []
	n_vertices = count_vertices(tree, f_vertices)

	log_prior = n_vertices * math.log(2) - math.lgamma(n_vertices + 2)
	for f_vertex in f_vertices:
		log_prior -= math.log(f_vertex)

	return log_prior

def log_labelled_histories(n_taxa):
	llh  = 0.0
	llh += math.lgamma(n_taxa + 1)
	llh += math.lgamma(n_taxa)
	llh -= (n_taxa - 1) * math.log(2)
	return llh

def log_subclade_rankings(n_taxa_1, n_taxa_2):
	lsr  = 0.0
	lsr += math.lgamma(n_taxa_1 + n_taxa_2 - 1)
	lsr -= math.lgamma(n_taxa_1)
	lsr -= math.lgamma(n_taxa_2)
	return lsr

def cc_yule_log_prior(cc_id_pair):
	parent_id, split_id = cc_id_pair
	child1_id, child2_id = elucidate_cc_split(parent_id, split_id)

	n_child1_taxa = clade_size(child1_id)
	n_child2_taxa = clade_size(child2_id)

	n_parent_taxa = n_child1_taxa + n_child2_taxa

	log_prior  = 0.0
	log_prior += log_labelled_histories(n_child1_taxa)
	log_prior += log_labelled_histories(n_child2_taxa)
	log_prior += log_subclade_rankings(n_child1_taxa, n_child2_taxa)
	log_prior -= log_labelled_histories(n_parent_taxa)

	return log_prior

def calculate_discrete_probabilities(discrete_frequencies, log_prior_function, log_prior_data, *log_prior_args):
	raw_log_probabilities = {}
	normalized_probabilities = {}

	df_array = discrete_frequencies.feature_frequencies
	n_features = len(df_array)
	for i in range(n_features):
		feature_hash   = df_array[i][0]
		feature_counts = list(df_array[i])[1:]
		feature_data   = log_prior_data[feature_hash]

		n_partitions = len(feature_counts)

		if 0 in feature_counts:
			normalized_probabilities[feature_hash] = 0.0
		else:
			log_probability = sum([math.log(p) for p in feature_counts]) - (n_partitions - 1) * log_prior_function(feature_data, *log_prior_args)
			raw_log_probabilities[feature_hash] = log_probability

	if len(raw_log_probabilities) > 0:
		log_sum_of_probabilities = numpy.logaddexp.reduce(raw_log_probabilities.values())

		for feature_hash in raw_log_probabilities:
			normalized_probability = math.exp(raw_log_probabilities[feature_hash] - log_sum_of_probabilities)
			normalized_probabilities[feature_hash] = normalized_probability

	return normalized_probabilities

def read_newick(newick_path):
	newick_file = open(newick_path)
	newick_strings = []
	l = newick_file.readline()
	while l != "":
		tree_string = l.strip()
		newick_strings.append(tree_string)
		l = newick_file.readline()
	return newick_strings

def generate_tree_array(newick_string, taxon_order, calibration_taxon = "", calibration_time = 0.0):
	if calibration_taxon == "":
		calibration_taxon = taxon_order[0]

	tree = ete2.Tree(newick_string)
	calibration_node = tree.get_leaves_by_name(calibration_taxon)[0]
	root_height = tree.get_distance(calibration_node) + calibration_time

	n_taxa = len(taxon_order)
	id_bytes, id_remainder = divmod(n_taxa, 8)
	if id_remainder == 0:
		id_size = id_bytes
	else:
		id_size = id_bytes + 1

	node_struct_format = "a%d,a%d,f8" % (id_size, id_size)

	tree_values = []
	recurse_node_properties(tree, taxon_order, root_height, tree_values)

	tree_array = numpy.array(tree_values, node_struct_format)
	tree_array.sort()

	return tree_array

def recurse_node_properties(node, taxon_order, node_height, tree_values):
	if node.is_leaf():
		child1_clade = set([node.name])
		child2_clade = set()
	else:
		child1, child2 = node.get_children() # assumes strictly bifurcating tree
		child1_clade = set(child1.get_leaf_names())
		child2_clade = set(child2.get_leaf_names())

		child1_height = node_height - child1.dist
		child2_height = node_height - child2.dist

		recurse_node_properties(child1, taxon_order, child1_height, tree_values)
		recurse_node_properties(child2, taxon_order, child2_height, tree_values)

	parent_id, split_id = calculate_node_ids(child1_clade, child2_clade, taxon_order)

	tree_values.append((parent_id, split_id, node_height))

def calculate_node_ids(children_a, children_b, taxon_order):
	n_taxa = len(taxon_order)
	children = set.union(children_a, children_b)

	parent_boolean = numpy.zeros(n_taxa, dtype=numpy.uint8)
	split_boolean = numpy.zeros(n_taxa, dtype=numpy.uint8)

	i = 0
	for j in range(n_taxa):
		t = taxon_order[j]
		if t in children:
			parent_boolean[j] = 1

			if i == 0:
				if t in children_a:
					a_first = True
				else:
					a_first = False

			if (t in children_b) ^ a_first: # first child always "True"
				split_boolean[i] = 1

			i += 1

	parent_packed = numpy.packbits(parent_boolean)
	split_packed = numpy.packbits(split_boolean)

	parent_id = parent_packed.tostring()
	split_id = split_packed.tostring()

	return parent_id, split_id

def calculate_clade_credibility(sample_probabilities, sample_strings, prediction_strings, taxon_order):
	flat_topology_hashes = []
	flat_topology_probabilities = []
	flat_topology_strings = []

	clade_probabilities = {}
	for topology_hash, topology_newick in sample_strings.items():
		topology_probability = sample_probabilities[topology_hash]
		topology_array = generate_tree_array(topology_newick, taxon_order)
		topology_clades = topology_array["f0"]

		for clade_id in topology_clades:
			if clade_id in clade_probabilities:
				clade_probabilities[clade_id] += topology_probability
			else:
				clade_probabilities[clade_id] = topology_probability

	prediction_credibilities = {}
	for prediction_hash, prediction_newick in prediction_strings.items():
		prediction_array = generate_tree_array(prediction_newick, taxon_order)
		prediction_clades = prediction_array["f0"]

		prediction_clade_probs = []
		for clade_id in prediction_clades:
			if clade_id in clade_probabilities:
				clade_probability = clade_probabilities[clade_id]
			else:
				clade_probability = 0.0

			prediction_clade_probs.append(clade_probability)

		prediction_clade_credibility = numpy.product(prediction_clade_probs)
		prediction_credibilities[prediction_hash] = prediction_clade_credibility

	return prediction_credibilities

def rank_discrete(discrete_probabilities):
	flat_discrete_hashes = []
	flat_discrete_probs = []
	for discrete_hash, discrete_probability in discrete_probabilities.items():
		if discrete_probability > 0.0:
			flat_discrete_hashes.append(discrete_hash)
			flat_discrete_probs.append(discrete_probability)

	ascending_prob_order = numpy.argsort(flat_discrete_probs)
	declining_prob_order = ascending_prob_order[::-1]

	discrete_rank = 0

	discrete_ranking = {}
	for i in declining_prob_order:
		discrete_rank += 1
		discrete_hash = flat_discrete_hashes[i]
		discrete_ranking[discrete_hash] = discrete_rank

	return discrete_ranking
