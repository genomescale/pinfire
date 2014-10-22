import numpy
import ete2
import dendropy
import math

class TopologySample():
	def __init__(self, newick_strings):
		self.taxon_order = []
		self.newick_strings = []
		self.topology_arrays = []

		self.newick_strings = newick_strings
		self.n_topologies = len(self.newick_strings)

		for i in range(self.n_topologies):
			ns = self.newick_strings[i]
			tree = ete2.Tree(ns)
			if i == 0:
				taxa = tree.get_leaf_names()
				self.taxon_order = sorted(taxa)

			self.generate_topology_array(ns)

	def generate_topology_array(self, newick_string):
		topology_root = ete2.Tree(newick_string)

		n_taxa = len(self.taxon_order)
		id_bytes, id_remainder = divmod(n_taxa, 8)
		if id_remainder == 0:
			id_size = id_bytes
		else:
			id_size = id_bytes + 1

		node_struct_format = "a%d,a%d" % (id_size, id_size)

		topology_values = []
		self.recurse_node_properties(topology_root, topology_values)

		topology_array = numpy.array(topology_values, node_struct_format)
		topology_array.sort() # clades should be sorted to that topology hashes are consistent
		self.topology_arrays.append(topology_array)

	def recurse_node_properties(self, node, topology_values):
		if not node.is_leaf():
			child1, child2 = node.get_children() # assumes strictly bifurcating tree
			child1_clade = set(child1.get_leaf_names())
			child2_clade = set(child2.get_leaf_names())

			self.recurse_node_properties(child1, topology_values)
			self.recurse_node_properties(child2, topology_values)

			parent_id, split_id = calculate_node_ids(child1_clade, child2_clade, self.taxon_order)
			topology_values.append((parent_id, split_id))

class UltrametricSample(TopologySample):
	def __init__(self, newick_strings, calibration_taxon, calibration_date):
		self.taxon_order = []
		self.newick_strings = []
		self.tree_arrays = []

		self.newick_strings = newick_strings
		self.n_trees = len(self.newick_strings)

		for i in range(self.n_trees):
			ns = self.newick_strings[i]
			tree = ete2.Tree(ns)
			if i == 0:
				taxa = tree.get_leaf_names()
				self.taxon_order = sorted(taxa)
				if calibration_taxon == "":
					calibration_taxon = self.taxon_order[0]

			self.generate_tree_array(ns, calibration_taxon, calibration_date)

	def generate_tree_array(self, newick_string, calibration_taxon, calibration_date):
		tree_root = ete2.Tree(newick_string)
		calibration_node = tree_root.get_leaves_by_name(calibration_taxon)[0]
		root_height = tree_root.get_distance(calibration_node) + calibration_date

		n_taxa = len(self.taxon_order)
		id_bytes, id_remainder = divmod(n_taxa, 8)
		if id_remainder == 0:
			id_size = id_bytes
		else:
			id_size = id_bytes + 1

		node_struct_format = "a%d,a%d,f8" % (id_size, id_size)

		tree_values = []
		self.recurse_node_properties(tree_root, root_height, tree_values)

		tree_array = numpy.array(tree_values, node_struct_format)
		tree_array.sort() # clades should be sorted to that topology hashes are consistent
		self.tree_arrays.append(tree_array)

	def recurse_node_properties(self, node, node_height, tree_values):
		if not node.is_leaf():
			child1, child2 = node.get_children() # assumes strictly bifurcating tree
			child1_clade = set(child1.get_leaf_names())
			child2_clade = set(child2.get_leaf_names())

			child1_height = node_height - child1.dist
			child2_height = node_height - child2.dist

			self.recurse_node_properties(child1, child1_height, tree_values)
			self.recurse_node_properties(child2, child2_height, tree_values)

			parent_id, split_id = calculate_node_ids(child1_clade, child2_clade, self.taxon_order)
			tree_values.append((parent_id, split_id, node_height))

class DiscreteProbabilities():
	def __init__(self, data):
		sorted_hashes = sorted(data.keys())
		self.n_features = len(sorted_hashes)
		self.hashes_array = numpy.array(sorted_hashes)

		sorted_data = [data[feature_hash] for feature_hash in self.hashes_array]
		self.data_array = numpy.array(sorted_data)
		self.probabilities = {}
		self.probabilities_array = numpy.zeros(self.n_features, dtype = numpy.float64)

	def probabilities_from_counts(self, counts):
		sorted_counts = [counts[feature_hash] for feature_hash in self.hashes_array]
		counts_array = numpy.array(sorted_counts, dtype = numpy.uint64)

		log_counts = {}
		for i in range(self.n_features):
			feature_hash  = self.hashes_array[i]
			feature_count = counts_array[i]

			log_counts[feature_hash] = math.log(feature_count)

		log_sum_of_counts = numpy.logaddexp.reduce(log_counts.values())
		for feature_hash in log_counts:
			normalized_probability = math.exp(log_counts[feature_hash] - log_sum_of_counts)
			self.probabilities[feature_hash] = normalized_probability

		self.convert_probabilities()

	def convert_probabilities(self):
		sorted_probabilities = [self.probabilities[feature_hash] for feature_hash in self.hashes_array]
		self.probabilities_array = numpy.array(sorted_probabilities, dtype = numpy.float64)

	def cull_probabilities(self, max_features, max_probability):
		topology_ascending_order = numpy.argsort(self.probabilities_array)
		topology_descending_order = topology_ascending_order[::-1]

		posterior_features = 0
		posterior_probability = 0.0

		cull_indices = []
		for i in topology_descending_order:
			if (posterior_features >= max_features) or (posterior_probability >= max_probability):
				cull_hash = self.hashes_array[i]
				self.probabilities.pop(cull_hash)
				cull_indices.append(i)

			posterior_features += 1
			posterior_probability += self.probabilities_array[i]

		self.probabilities_array = numpy.delete(self.probabilities_array, cull_indices)
		self.hashes_array = numpy.delete(self.hashes_array, cull_indices)
		self.data_array = numpy.delete(self.data_array, cull_indices)

		self.n_features = len(self.probabilities_array)

class TopologyProbabilities(DiscreteProbabilities):
	def probabilities_from_ccs(self, cc_probabilities):
		topology_sample = TopologySample(self.data_array)

		for i in range(self.n_features):
			topology_array = topology_sample.topology_arrays[i]
			topology_hash = self.hashes_array[i]
			node_probabilities = []
			for node in topology_array:
				parent_hash = node[0].tostring() # the hash for the clade
				split_hash = node[1].tostring() # the hash for the bifurcation

				n_node_taxa = clade_size(parent_hash)
				if n_node_taxa >= 3: # conditional clade
					split_probability = cc_probabilities[parent_hash].probabilities[split_hash]
					node_probabilities.append(split_probability)

			topology_probability = numpy.prod(node_probabilities)
			self.probabilities[topology_hash] = topology_probability

		self.convert_probabilities()

def calculate_topology_probabilities(ts):
	topology_counts = {}
	topology_data = {}
	cc_counts = {}
	cc_data = {}

	for i in range(ts.n_trees):
		tree_array = ts.tree_arrays[i]
		topology_hash = tree_array["f0"].tostring() # topology hash is concatenated, sorted clade hashes

		if topology_hash not in topology_counts: # record topology
			tree_newick = ts.newick_strings[i]
			tree_root = ete2.Tree(tree_newick)
			topology_newick = tree_root.write(format = 9) # strip branch lengths
			topology_data[topology_hash] = topology_newick
			topology_counts[topology_hash] = 1
		else:
			topology_counts[topology_hash] += 1

		topology_array = tree_array[["f0", "f1"]] # we are only interested in clade & split hashes, not node heights
		for node in topology_array:
			parent_hash = node[0].tostring() # the hash for the clade
			split_hash = node[1].tostring() # the hash for the bifurcation

			n_node_taxa = clade_size(parent_hash)
			if n_node_taxa >= 3: # record conditional clade
				if parent_hash not in cc_counts:
					cc_data[parent_hash] = {split_hash: node}
					cc_counts[parent_hash] = {split_hash: 1}
				elif split_hash not in cc_counts[parent_hash]:
					cc_data[parent_hash][split_hash] = node
					cc_counts[parent_hash][split_hash] = 1
				else:
					cc_counts[parent_hash][split_hash] += 1

	topology_set = TopologyProbabilities(topology_data)

	cc_sets = {}
	for parent_hash, splits_data in cc_data.items():
		cc_sets[parent_hash] = DiscreteProbabilities(splits_data)

	return topology_set, topology_counts, cc_sets, cc_counts

def derive_best_topologies(cc_probabilities, taxon_order, trees_threshold = 1, post_threshold = 1.0):
	cherry_hash = "\x80"

	n_taxa = len(taxon_order)
	root_hash = calculate_root_hash(n_taxa)
	n_bytes = len(root_hash) # maximum number of bytes required to store parent clade hashes or split hashes

	derived_struct_format = "a%d,a%d,u1,f8" % (n_bytes, n_bytes)

	star_tree = numpy.array([(root_hash, "", 1, 1.0)], dtype=derived_struct_format)
	candidate_topologies = [star_tree]
	candidate_inv_probs = [0.0]

	best_topologies = []
	best_probabilities = []
	best_posterior = 0.0
	while (len(candidate_topologies) > 0) and (len(best_topologies) < trees_threshold) and (best_posterior < post_threshold):
		candidate_topology = candidate_topologies.pop(0)
		candidate_inv_prob = candidate_inv_probs.pop(0)
		candidate_nodes = numpy.flatnonzero(candidate_topology["f2"])

		if len(candidate_nodes) == 0: # candidate topology is fully resolved
			candidate_probability = 1.0 - candidate_inv_prob
			best_topologies.append(candidate_topology)
			best_probabilities.append(candidate_probability)
			best_posterior += candidate_probability
		else: # candidate topology is not fully resolved
			unresolved_node_index = candidate_nodes[0]
			unresolved_node_hash = candidate_topology[unresolved_node_index]["f0"].tostring()
			split_probabilities = cc_probabilities[unresolved_node_hash]

			new_candidate_topologies = []
			new_candidate_inv_probs = []
			for split_hash in split_probabilities:
				split_probability = split_probabilities[split_hash]
				if split_probability > 0.0:
					child1_hash, child2_hash = elucidate_cc_split(unresolved_node_hash, split_hash)
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

					new_topology = numpy.concatenate([candidate_topology] + new_topology_rows)
					new_topology[unresolved_node_index]["f1"] = split_hash
					new_topology[unresolved_node_index]["f2"] = 0
					new_topology[unresolved_node_index]["f3"] = split_probability
					new_candidate_topologies.append(new_topology)

					new_topology_inv_probability = 1.0 - numpy.prod(new_topology["f3"])
					new_candidate_inv_probs.append(new_topology_inv_probability)

			integrate_probability(candidate_inv_probs, candidate_topologies, new_candidate_inv_probs, new_candidate_topologies)
		print(len(candidate_topologies), len(best_topologies), sum([1.0 - p for p in candidate_inv_probs]), best_posterior) # number of candidate and best topologies, total posterior of candidate and best topologies

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

# read a nexus or newick format file containing phylogenetic trees
# if the file does not begin with a nexus header, assumes it is a newick file
# returns a list of newick strings, in the same order as the input file
def trees_from_path(trees_filepath):
	nexus_header = "#NEXUS"

	trees_file = open(trees_filepath)
	first_line = trees_file.readline().strip().upper()
	trees_file.seek(0)

	if first_line == nexus_header: # looks like a nexus file, convert to newick
		trees_list = dendropy.TreeList.get_from_stream(trees_file, schema = "nexus")
		trees_file.close()

		newick_blob = trees_list.as_string("newick", suppress_rooting = True)
	else: # assume file is already in newick format
		newick_blob = trees_file.read()
		trees_file.close()

	newick_strings = newick_blob.strip().split("\n")
	return newick_strings

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