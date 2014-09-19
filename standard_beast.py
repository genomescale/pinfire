import xml.etree.ElementTree

Element = xml.etree.ElementTree.Element

chain_length = str(22222 * 5000)
store_every = str(5000)
log_every = str(5000)

def indent(elem, level=0):
	i = "\n" + level*"  "
	if len(elem):
		if not elem.text or not elem.text.strip():
			elem.text = i + "  "
		if not elem.tail or not elem.tail.strip():
			elem.tail = i
		for elem in elem:
			indent(elem, level+1)
		if not elem.tail or not elem.tail.strip():
			elem.tail = i
	else:
		if level and (not elem.tail or not elem.tail.strip()):
			elem.tail = i

def make_xml(tree_name, taxa_alignment):
	beast_element = Element('beast', {'beautitemplate': 'Standard', 'beautistatus': '', 'version': '2.0', 'namespace': 'beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood'})
	beast_xml = xml.etree.ElementTree.ElementTree(beast_element)

	taxa_alignment_element = Element('data', {'id': '' + tree_name, 'name': 'alignment'})
	for taxon_name, taxon_sequence in taxa_alignment.items():
		sequence_element = Element('sequence', {'totalcount': '4', 'id': 'seq_' + taxon_name, 'value': taxon_sequence, 'taxon': taxon_name})
		taxa_alignment_element.append(sequence_element)

	map_beta_element = Element('map', {'name': 'Beta'})
	map_beta_element.text = 'beast.math.distributions.Beta'
	map_exponential_element = Element('map', {'name': 'Exponential'})
	map_exponential_element.text = 'beast.math.distributions.Exponential'
	map_inversegamma_element = Element('map', {'name': 'InverseGamma'})
	map_inversegamma_element.text = 'beast.math.distributions.InverseGamma'
	map_lognormal_element = Element('map', {'name': 'LogNormal'})
	map_lognormal_element.text = 'beast.math.distributions.LogNormalDistributionModel'
	map_gamma_element = Element('map', {'name': 'Gamma'})
	map_gamma_element.text = 'beast.math.distributions.Gamma'
	map_uniform_element = Element('map', {'name': 'Uniform'})
	map_uniform_element.text = 'beast.math.distributions.Uniform'
	map_prior_element = Element('map', {'name': 'prior'})
	map_prior_element.text = 'beast.math.distributions.Prior'
	map_laplacedistribution_element = Element('map', {'name': 'LaplaceDistribution'})
	map_laplacedistribution_element.text = 'beast.math.distributions.LaplaceDistribution'
	map_oneonx_element = Element('map', {'name': 'OneOnX'})
	map_oneonx_element.text = 'beast.math.distributions.OneOnX'
	map_normal_element = Element('map', {'name': 'Normal'})
	map_normal_element.text = 'beast.math.distributions.Normal'

	mcmc_element = Element('run', {'chainLength': chain_length, 'storeEvery': store_every, 'id': 'mcmc', 'spec': 'MCMC'})
	state_element = Element('state', {'storeEvery': store_every, 'id': 'state'})
	tree_element = Element('tree', {'id': 'Tree.t:' + tree_name, 'name': 'stateNode'})
	taxon_set_element = Element('taxonset', {'id': 'TaxonSet.' + tree_name, 'spec': 'TaxonSet'})
	data_alignment_element = Element('data', {'idref': '' + tree_name, 'name': 'alignment'})
	taxon_set_element.append(data_alignment_element)
	tree_element.append(taxon_set_element)

	birth_rate_element = Element('parameter', {'upper': '10000.0', 'id': 'birthRate.t:' + tree_name, 'name': 'stateNode'})
	birth_rate_element.text = '213.2021'
	kappa_element = Element('parameter', {'lower': '0.0', 'id': 'kappa.s:' + tree_name, 'name': 'stateNode'})
	kappa_element.text = '2.0'
	state_element.append(tree_element)
	state_element.append(birth_rate_element)
	state_element.append(kappa_element)

	random_tree_element = Element('init', {'taxa': '@' + tree_name, 'estimate': 'false', 'initial': '@Tree.t:' + tree_name, 'id': 'RandomTree.t:' + tree_name, 'spec': 'beast.evolution.tree.RandomTree'})
	constant_population0_element = Element('populationModel', {'id': 'ConstantPopulation0.t:' + tree_name, 'spec': 'ConstantPopulation'})
	random_pop_size_element = Element('parameter', {'id': 'randomPopSize.t:' + tree_name, 'name': 'popSize'})
	random_pop_size_element.text = '1.0'
	constant_population0_element.append(random_pop_size_element)
	random_tree_element.append(constant_population0_element)

	posterior_element = Element('distribution', {'id': 'posterior', 'spec': 'util.CompoundDistribution'})
	prior_element = Element('distribution', {'id': 'prior', 'spec': 'util.CompoundDistribution'})
	yule_model_element = Element('distribution', {'birthDiffRate': '@birthRate.t:' + tree_name, 'tree': '@Tree.t:' + tree_name, 'id': 'YuleModel.t:' + tree_name, 'spec': 'beast.evolution.speciation.YuleModel'})
	yule_birth_rate_prior_element = Element('prior', {'x': '@birthRate.t:' + tree_name, 'id': 'YuleBirthRatePrior.t:' + tree_name, 'name': 'distribution'})
	uniform_element = Element('Uniform', {'upper': 'Infinity', 'id': 'Uniform.0', 'name': 'distr'})
	yule_birth_rate_prior_element.append(uniform_element)

	kappa_prior_element = Element('prior', {'x': '@kappa.s:' + tree_name, 'id': 'KappaPrior.s:' + tree_name, 'name': 'distribution'})
	log_normal_distribution_model_element = Element('LogNormal', {'id': 'LogNormalDistributionModel.0', 'name': 'distr'})
	m_element = Element('parameter', {'estimate': 'false', 'id': 'RealParameter.0', 'name': 'M'})
	m_element.text = '1.0'
	s_element = Element('parameter', {'estimate': 'false', 'id': 'RealParameter.01', 'name': 'S'})
	s_element.text = '1.25'
	log_normal_distribution_model_element.append(m_element)
	log_normal_distribution_model_element.append(s_element)
	kappa_prior_element.append(log_normal_distribution_model_element)
	prior_element.append(yule_model_element)
	prior_element.append(yule_birth_rate_prior_element)
	prior_element.append(kappa_prior_element)

	likelihood_element = Element('distribution', {'id': 'likelihood', 'spec': 'util.CompoundDistribution'})
	tree_likelihood_element = Element('distribution', {'data': '@' + tree_name, 'tree': '@Tree.t:' + tree_name, 'id': 'treeLikelihood.' + tree_name, 'spec': 'TreeLikelihood'})
	site_model_element = Element('siteModel', {'id': 'SiteModel.s:' + tree_name, 'spec': 'SiteModel'})
	mutation_rate_element = Element('parameter', {'estimate': 'false', 'id': 'mutationRate.s:' + tree_name, 'name': 'mutationRate'})
	mutation_rate_element.text = '1.0'
	gamma_shape_element = Element('parameter', {'estimate': 'false', 'id': 'gammaShape.s:' + tree_name, 'name': 'shape'})
	gamma_shape_element.text = '1.0'
	proportion_invariant_element = Element('parameter', {'upper': '1.0', 'lower': '0.0', 'estimate': 'false', 'name': 'proportionInvariant', 'id': 'proportionInvariant.s:' + tree_name})
	proportion_invariant_element.text = '0.0'
	hky_element = Element('substModel', {'kappa': '@kappa.s:' + tree_name, 'id': 'hky.s:' + tree_name, 'spec': 'HKY'})
	empirical_freqs_element = Element('frequencies', {'data': '@' + tree_name, 'id': 'empiricalFreqs.s:' + tree_name, 'spec': 'Frequencies'})
	hky_element.append(empirical_freqs_element)
	site_model_element.append(mutation_rate_element)
	site_model_element.append(gamma_shape_element)
	site_model_element.append(proportion_invariant_element)
	site_model_element.append(hky_element)

	strict_clock_element = Element('branchRateModel', {'id': 'StrictClock.c:' + tree_name, 'spec': 'beast.evolution.branchratemodel.StrictClockModel'})
	clock_rate_element = Element('parameter', {'estimate': 'false', 'id': 'clockRate.c:' + tree_name, 'name': 'clock.rate'})
	clock_rate_element.text = '1.0'
	strict_clock_element.append(clock_rate_element)
	tree_likelihood_element.append(site_model_element)
	tree_likelihood_element.append(strict_clock_element)
	likelihood_element.append(tree_likelihood_element)
	posterior_element.append(prior_element)
	posterior_element.append(likelihood_element)

	yule_birth_rate_scaler_element = Element('operator', {'parameter': '@birthRate.t:' + tree_name, 'scaleFactor': '0.75', 'id': 'YuleBirthRateScaler.t:' + tree_name, 'weight': '3.0', 'spec': 'ScaleOperator'})
	tree_scaler_element = Element('operator', {'scaleFactor': '0.5', 'tree': '@Tree.t:' + tree_name, 'id': 'treeScaler.t:' + tree_name, 'weight': '3.0', 'spec': 'ScaleOperator'})
	tree_root_scaler_element = Element('operator', {'weight': '3.0', 'rootOnly': 'true', 'scaleFactor': '0.5', 'tree': '@Tree.t:' + tree_name, 'spec': 'ScaleOperator', 'id': 'treeRootScaler.t:' + tree_name})
	uniform_operator_element = Element('operator', {'tree': '@Tree.t:' + tree_name, 'id': 'UniformOperator.t:' + tree_name, 'weight': '30.0', 'spec': 'Uniform'})
	subtree_slide_element = Element('operator', {'tree': '@Tree.t:' + tree_name, 'id': 'SubtreeSlide.t:' + tree_name, 'weight': '15.0', 'spec': 'SubtreeSlide'})
	narrow_element = Element('operator', {'tree': '@Tree.t:' + tree_name, 'id': 'narrow.t:' + tree_name, 'weight': '15.0', 'spec': 'Exchange'})
	wide_element = Element('operator', {'isNarrow': 'false', 'tree': '@Tree.t:' + tree_name, 'spec': 'Exchange', 'weight': '3.0', 'id': 'wide.t:' + tree_name})
	wilson_balding_element = Element('operator', {'tree': '@Tree.t:' + tree_name, 'id': 'WilsonBalding.t:' + tree_name, 'weight': '3.0', 'spec': 'WilsonBalding'})
	kappa_scaler_element = Element('operator', {'parameter': '@kappa.s:' + tree_name, 'scaleFactor': '0.5', 'id': 'KappaScaler.s:' + tree_name, 'weight': '0.1', 'spec': 'ScaleOperator'})

	tracelog_element = Element('logger', {'sort': 'smart', 'logEvery': log_every, 'sanitiseHeaders': 'true', 'model': '@posterior', 'id': 'tracelog', 'fileName': 'beast_' + tree_name + '.log'})
	tracelog_posterior_element = Element('log', {'idref': 'posterior'})
	tracelog_likelihood_element = Element('log', {'idref': 'likelihood'})
	tracelog_prior_element = Element('log', {'idref': 'prior'})
	tracelog_treelikelihood_element = Element('log', {'idref': 'treeLikelihood.' + tree_name})
	tracelog_treeheight_element = Element('log', {'tree': '@Tree.t:' + tree_name, 'id': 'TreeHeight.t:' + tree_name, 'spec': 'beast.evolution.tree.TreeHeightLogger'})
	tracelog_yulemodel_element = Element('log', {'idref': 'YuleModel.t:' + tree_name})
	tracelog_birthrate_element = Element('parameter', {'idref': 'birthRate.t:' + tree_name, 'name': 'log'})
	tracelog_kappa_element = Element('parameter', {'idref': 'kappa.s:' + tree_name, 'name': 'log'})
	tracelog_element.append(tracelog_posterior_element)
	tracelog_element.append(tracelog_likelihood_element)
	tracelog_element.append(tracelog_prior_element)
	tracelog_element.append(tracelog_treelikelihood_element)
	tracelog_element.append(tracelog_treeheight_element)
	tracelog_element.append(tracelog_yulemodel_element)
	tracelog_element.append(tracelog_birthrate_element)
	tracelog_element.append(tracelog_kappa_element)

	screenlog_element = Element('logger', {'logEvery': log_every, 'id': 'screenlog'})
	screenlog_posterior_element = Element('log', {'idref': 'posterior'})
	screenlog_ess_element = Element('log', {'arg': '@posterior', 'spec': 'util.ESS', 'id': 'ESS.0'})
	screenlog_likelihood_element = Element('log', {'idref': 'likelihood'})
	screenlog_prior_element = Element('log', {'idref': 'prior'})
	screenlog_element.append(screenlog_posterior_element)
	screenlog_element.append(screenlog_ess_element)
	screenlog_element.append(screenlog_likelihood_element)
	screenlog_element.append(screenlog_prior_element)

	treelog_element = Element('logger', {'fileName': 'species_' + tree_name + '.trees', 'logEvery': log_every, 'mode': 'tree', 'id': 'treelog.t:' + tree_name})
	treelog_treewithmetadatalogger_element = Element('log', {'tree': '@Tree.t:' + tree_name, 'id': 'TreeWithMetaDataLogger.t:' + tree_name, 'spec': 'beast.evolution.tree.TreeWithMetaDataLogger'})
	treelog_element.append(treelog_treewithmetadatalogger_element)

	mcmc_element.append(state_element)
	mcmc_element.append(random_tree_element)
	mcmc_element.append(posterior_element)
	mcmc_element.append(yule_birth_rate_scaler_element)
	mcmc_element.append(tree_scaler_element)
	mcmc_element.append(tree_root_scaler_element)
	mcmc_element.append(uniform_operator_element)
	mcmc_element.append(subtree_slide_element)
	mcmc_element.append(narrow_element)
	mcmc_element.append(wide_element)
	mcmc_element.append(wilson_balding_element)
	mcmc_element.append(kappa_scaler_element)
	mcmc_element.append(tracelog_element)
	mcmc_element.append(screenlog_element)
	mcmc_element.append(treelog_element)
	beast_element.append(taxa_alignment_element)
	beast_element.append(map_beta_element)
	beast_element.append(map_exponential_element)
	beast_element.append(map_inversegamma_element)
	beast_element.append(map_lognormal_element)
	beast_element.append(map_gamma_element)
	beast_element.append(map_uniform_element)
	beast_element.append(map_prior_element)
	beast_element.append(map_laplacedistribution_element)
	beast_element.append(map_oneonx_element)
	beast_element.append(map_normal_element)
	beast_element.append(mcmc_element)

	indent(beast_element)
	return beast_xml