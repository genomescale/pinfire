import numpy
import ete2
import math

order_differs_error = Exception("The set of taxa differs between tree samples, so they may not be directly compared")

class UltrametricSample():
	def __init__(self, newick_strings, calibration_taxon, calibration_date):
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
				if calibration_taxon == "":
					calibration_taxon = self.taxon_order[0]

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
		max_label_size = 0

		frequencies_table = []
		for i in range(n_features):
			feature = frequency_labels[i]
			count = frequencies[feature]
			feature_frequencies = tuple(feature, count)
			frequencies_table.append(feature_frequencies)

			if len(feature) > max_label_size:
				max_label_size = len(feature)

		self.feature_struct_format = "a" + str(max_label_size) + "u8"
		self.feature_frequencies = numpy.array(frequencies_table, self.feature_struct_format)

def defined_topology_probabilities(cc_probabilities, ts):
	topology_probabilities = {}
	topology_strings = {}

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
				if n_node_taxa > 2:
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

def integrate_probability(original_probs, original_data, new_probs, new_data):
	for i in range(len(new_probs)):
		probability = new_probs[i]
		observation = new_data[i]

		rank = numpy.searchsorted(original_probs, probability)

		original_probs.insert(rank, probability)
		original_data.insert(rank, observation)

def best_topology_probabilities(cc_probabilities, taxon_order, trees_threshold = 1, post_threshold = 1.0):
	n_taxa = len(taxon_order)
	n_bytes, remainder = divmod(n_taxa, 8)
	if remainder > 0:
		n_bytes += 1

	derived_struct_format = "a%d,a%d,u1,f8" % (n_bytes, n_bytes)

	cherry_hash = "\x80"
	root_hash = calculate_root_hash(n_taxa)
	star_tree = numpy.array([(root_hash, "", 1, 1.0)], dtype=derived_struct_format)
	unresolved_topologies = [star_tree]
	unresolved_inv_probs = [0.0]

	resolved_topologies = []
	resolved_inv_probs = []

	best_topologies = []
	best_probabilities = []
	best_posterior = 0.0
	while ((len(unresolved_topologies) > 0) or (len(resolved_topologies) > 0)) and (len(best_topologies) < trees_threshold) and (best_posterior < post_threshold):
		print(len(unresolved_topologies), len(resolved_topologies), len(best_topologies), sum([1.0 - p for p in unresolved_inv_probs]), sum([1.0 - p for p in resolved_inv_probs]), best_posterior)

		if len(unresolved_topologies) == 0:
			max_unresolved_prob = 0.0
		else:
			max_unresolved_prob = 1.0 - unresolved_inv_probs[0]

		if len(resolved_topologies) == 0:
			max_resolved_prob = 0.0
		else:
			max_resolved_prob = 1.0 - resolved_inv_probs[0]

		if max_resolved_prob >= max_unresolved_prob:
			resolved_topology = resolved_topologies.pop(0)
			resolved_inv_probs = resolved_inv_probs[1:]
			best_topologies.append(resolved_topology)
			best_probabilities.append(max_resolved_prob)
			best_posterior += max_resolved_prob
		else:
			unresolved_topology = unresolved_topologies.pop(0)
			unresolved_inv_probs = unresolved_inv_probs[1:]

			unresolved_nodes = numpy.flatnonzero(unresolved_topology["f2"])
			unresolved_index = unresolved_nodes[0]
			unresolved_hash = unresolved_topology[unresolved_index]["f0"].tostring()
			split_probabilities = cc_probabilities[unresolved_hash]

			new_unresolved_topologies = []
			new_unresolved_inv_probs = []
			new_resolved_topologies = []
			new_resolved_inv_probs = []
			for split_hash in split_probabilities:
				split_probability = split_probabilities[split_hash]
				if split_probability > 0.0:
					child1_hash, child2_hash = elucidate_cc_split(unresolved_hash, split_hash)
					child1_size = clade_size(child1_hash)
					child2_size = clade_size(child2_hash)

					new_topology_rows = []
					if child1_size > 1:
						if child1_size == 2: # resolved (cherry)
							child1_row = numpy.array([(child1_hash, cherry_hash, 0, 1.0)], dtype=derived_struct_format)
						else: # unresolved (more than two taxa)
							child1_row = numpy.array([(child1_hash, "", 1, 1.0)], dtype=derived_struct_format)
						new_topology_rows.append(child1_row)

					if child2_size > 1:
						if child2_size == 2: # resolved (cherry)
							child2_row = numpy.array([(child2_hash, cherry_hash, 0, 1.0)], dtype=derived_struct_format)
						else: # unresolved (more than two taxa)
							child2_row = numpy.array([(child2_hash, "", 1, 1.0)], dtype=derived_struct_format)
						new_topology_rows.append(child2_row)

					new_topology = numpy.concatenate([unresolved_topology] + new_topology_rows)
					new_topology[unresolved_index]["f1"] = split_hash
					new_topology[unresolved_index]["f2"] = 0
					new_topology[unresolved_index]["f3"] = split_probability

					n_unresolved_nodes = numpy.sum(new_topology["f2"])
					new_topology_inv_probability = 1.0 - numpy.prod(new_topology["f3"])
					if n_unresolved_nodes == 0:
						new_resolved_topologies.append(new_topology)
						new_resolved_inv_probs.append(new_topology_inv_probability)
					else:
						new_unresolved_topologies.append(new_topology)
						new_unresolved_inv_probs.append(new_topology_inv_probability)

			integrate_probability(resolved_inv_probs, resolved_topologies, new_resolved_inv_probs, new_resolved_topologies)
			integrate_probability(unresolved_inv_probs, unresolved_topologies, new_unresolved_inv_probs, new_unresolved_topologies)

	print(len(unresolved_topologies), len(resolved_topologies), len(best_topologies), sum([1.0 - p for p in unresolved_inv_probs]), sum([1.0 - p for p in resolved_inv_probs]), best_posterior)

	derived_topology_probabilities = {}
	derived_topology_newick = {}
	for i in range(len(best_topologies)):
		topology = best_topologies[i]
		probability = best_probabilities[i]
		topology_hash = numpy.sort(topology["f0"]).tostring()

		splits = {}
		for node in topology:
			parent_id = node["f0"]
			split_id = node["f1"]
			splits[parent_id] = split_id

		tree_model = ete2.Tree()
		derive_tree_from_splits(tree_model, root_hash, taxon_order, splits)
		newick = tree_model.write(format = 9)

		derived_topology_probabilities[topology_hash] = probability
		derived_topology_newick[topology_hash] = newick

	return derived_topology_probabilities, derived_topology_newick

def derive_tree_from_splits(current_node, parent_hash, taxon_order, splits):
	split_hash = splits[parent_hash]
	child1_hash, child2_hash = elucidate_cc_split(parent_hash, split_hash)

	child1_node = ete2.Tree()
	child2_node = ete2.Tree()

	current_node.add_child(child1_node)
	current_node.add_child(child2_node)

	child1_size = clade_size(child1_hash)
	child2_size = clade_size(child2_hash)

	if child1_size == 1:
		child1_node.name = clade_taxon_names(child1_hash, taxon_order)[0]
	else:
		derive_tree_from_splits(child1_node, child1_hash, taxon_order, splits)

	if child2_size == 1:
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

def count_topologies(ts):
	topology_counts = {}
	topology_strings = {}
	cc_counts = {}
	cc_id_pairs = {}

	n_trees = len(ts.tree_arrays)
	for j in range(n_trees):
		tree_array = ts.tree_arrays[j]
		newick_string = ts.newick_strings[j]
		topology_hash = tree_array["f0"].tostring()

		if topology_hash not in topology_counts:
			topology_tree = ete2.Tree(newick_string)
			topology_string = topology_tree.write(format = 9) # strip branch lengths

			topology_strings[topology_hash] = topology_string
			topology_counts[topology_hash] = 1
		else:
			topology_counts[topology_hash] += 1

		for node in tree_array:
			parent_id = node[0]
			split_id = node[1]

			n_node_taxa = clade_size(parent_id)

			if n_node_taxa > 2:
				if parent_id not in cc_counts:
					cc_id_pairs[parent_id] = {split_id: (parent_id, split_id)}
					cc_counts[parent_id] = {split_id: 1}
				elif split_id not in cc_counts[parent_id]:
					cc_id_pairs[parent_id][split_id] = (parent_id, split_id)
					cc_counts[parent_id][split_id] = 1
				else:
					cc_counts[parent_id][split_id][i] += 1

	topology_frequencies = DiscreteFrequencies(topology_counts)
	cc_frequencies = {}

	for parent_id in cc_counts:
		conditional_count = cc_counts[parent_id]
		cc_frequencies[parent_id] = DiscreteFrequencies(conditional_count)

	return topology_frequencies, topology_strings, cc_frequencies, cc_id_pairs

def calculate_discrete_probabilities(discrete_frequencies):
	raw_log_probabilities = {}
	normalized_probabilities = {}

	df_array = discrete_frequencies.feature_frequencies
	n_features = len(df_array)
	for i in range(n_features):
		feature_hash  = df_array[i][0]
		feature_count = df_array[i][1]

		log_probability = sum([math.log(p) for p in feature_counts])
		raw_log_probabilities[feature_hash] = log_probability

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

def generate_tree_array(newick_string, taxon_order, calibration_taxon, calibration_date):
	tree = ete2.Tree(newick_string)
	calibration_node = tree.get_leaves_by_name(calibration_taxon)[0]
	root_height = tree.get_distance(calibration_node) + calibration_date

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
	if not node.is_leaf():
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

def reverse_cc_probabilities(cc_probabilities):
	reverse_ccp = {}
	for parent_id in cc_probabilities:
		for split_id in cc_probabilities[parent_id]:
			cc_probability = cc_probabilities[parent_id][split_id]
			child1_hash, child2_hash = elucidate_cc_split(parent_id, split_id)

			if child1_hash in reverse_ccp:
				reverse_ccp[child1_hash][parent_id] = cc_probability
			else:
				reverse_ccp[child1_hash] = {parent_id: cc_probability}

			if child2_hash in reverse_ccp:
				reverse_ccp[child2_hash][parent_id] = cc_probability
			else:
				reverse_ccp[child2_hash] = {parent_id: cc_probability}

	return reverse_ccp

def derive_clade_probabilities(clade_id, n_taxa, reverse_ccp):
	n_bytes, remainder = divmod(n_taxa, 8)
	if remainder > 0:
		n_bytes += 1

	derived_struct_format = "a%d,f8" % (n_bytes)

	root_hash = calculate_root_hash(n_taxa)

	target_clade = numpy.array([(clade_id, 1.0)], dtype=derived_struct_format)
	paths_to_root = [target_clade]

	target_clade_probability = 0.0
	while len(paths_to_root) > 0:
		incomplete_path = paths_to_root.pop(0)
		incomplete_path_head = incomplete_path[0]
		incomplete_path_hash = incomplete_path_head["f0"].tostring()
		possible_parents = reverse_ccp[incomplete_path_hash]

		for parent_hash in possible_parents:
			split_probability = possible_parents[parent_hash]
			if split_probability > 0.0:
				parent_row = numpy.array([(parent_hash, split_probability)], dtype=derived_struct_format)
				new_path = numpy.concatenate((parent_row, incomplete_path))

				if parent_hash == root_hash:
					path_probability = numpy.prod(new_path["f1"])
					target_clade_probability += path_probability
				else:
					paths_to_root.append(new_path)

	return target_clade_probability

def add_derived_probabilities(newick_strings, taxon_order, ccp):
	n_taxa = len(taxon_order)
	reverse_ccp = reverse_cc_probabilities(ccp)

	annotated_topologies = {}
	for topology_hash, topology_newick in newick_strings.items():
		root_node = ete2.Tree(topology_newick)
		for node in root_node.get_descendants():
			if not node.is_leaf():
				child1, child2 = node.get_children()
				child1_clade = set(child1.get_leaf_names())
				child2_clade = set(child2.get_leaf_names())

				parent_id, split_id = calculate_node_ids(child1_clade, child2_clade, taxon_order)
				clade_probability = derive_clade_probabilities(parent_id, n_taxa, reverse_ccp)
				node.support = clade_probability

		annotated_newick = root_node.write(format = 2)
		annotated_topologies[topology_hash] = annotated_newick

	return annotated_topologies

def calculate_root_hash(n_taxa):
	root_hash_bits = numpy.ones(n_taxa, dtype = numpy.uint8)
	root_hash_bytes = numpy.packbits(root_hash_bits)
	root_hash = root_hash_bytes.tostring()

	return root_hash

def nonzero_derived_topologies(ccp, taxon_order):
	n_taxa = len(taxon_order)
	reverse_ccp = reverse_cc_probabilities(ccp)
	clades_by_size = []
	n_subtrees = {}

	for i in range(n_taxa - 2):
		clades_by_size.append(set())

	for parent_id in ccp:
		parent_size = clade_size(parent_id)
		if parent_size > 2:
			clades_by_size[parent_size - 3].add(parent_id)

	for i in range(n_taxa - 2):
		for parent_id in clades_by_size[i]:
			n_parent_subtrees = 0
			for split_id in ccp[parent_id]:
				if ccp[parent_id][split_id] > 0.0:
					child1_id, child2_id = elucidate_cc_split(parent_id, split_id)

					n_split_subtrees = 1
					child1_size = clade_size(child1_id)
					if child1_size > 2:
						n_split_subtrees = n_split_subtrees * n_subtrees[child1_id]

					child2_size = clade_size(child2_id)
					if child2_size > 2:
						n_split_subtrees = n_split_subtrees * n_subtrees[child2_id]

					n_parent_subtrees += n_split_subtrees

			n_subtrees[parent_id] = n_parent_subtrees

	root_id = calculate_root_hash(n_taxa)
	n_root_topologies = n_subtrees[root_id]

	return n_root_topologies