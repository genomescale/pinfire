import xml.etree.ElementTree
import math

Element = xml.etree.ElementTree.Element

#chain_length = 11111 * 10000
#store_every = 10000
#log_every = 10000

chain_length = 22222 * 5000
store_every = 5000
log_every = 5000

beast_namespace = "beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood"
beast_maps = {"Beta": "beast.math.distributions.Beta", "Exponential": "beast.math.distributions.Exponential", "InverseGamma": "beast.math.distributions.InverseGamma", "LogNormal": "beast.math.distributions.LogNormalDistributionModel", "Gamma": "beast.math.distributions.Gamma", "Uniform": "beast.math.distributions.Uniform", "prior": "beast.math.distributions.Prior", "LaplaceDistribution": "beast.math.distributions.LaplaceDistribution", "OneOnX": "beast.math.distributions.OneOnX", "Normal": "beast.math.distributions.Normal"}

all_locus_weights  = {"kappa": 0.1, "narrow": 15.0, "subtreeslide": 15.0, "tree": 3.0, "treeroot": 3.0, "uniformoperator": 30.0, "wide": 3.0, "wilsonbalding": 3.0}
additional_weights = {"strictclockrate": 3.0, "strictclockupdownoperator": 3.0, "updown": 3.0}
gt_weights = {}
gt_weights.update(all_locus_weights)
gt_weights.update(additional_weights)

st_weights = {"popmean": 3.0, "popsizebottom": 5.0, "popsizetop": 5.0, "reheight": 94.0, "updown": 20.0, "yulebirthrate": 3.0}

def get_weight(element, multiplier, gene = True):
    if gene:
        return str(gt_weights[element] * multiplier)
    else:
        return str(st_weights[element] * multiplier)

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

def make_xml(tree_name, species_tree):
    n_loci = len(species_tree)

    sum_st_weights = sum(st_weights.values())
    sum_gt_weights = sum(all_locus_weights.values()) * n_loci + sum(additional_weights.values()) * (n_loci - 1)
    old_ratio = sum_st_weights / sum_gt_weights
    new_ratio = 0.25
    st_gt_multiplier = new_ratio / old_ratio

    root_element = Element("beast", {"beautitemplate": "StarBeast", "beautistatus": "", "namespace": beast_namespace, "version": "2.0"})
    beast_xml = xml.etree.ElementTree.ElementTree(root_element)

    for m in beast_maps:
        map_element = Element("map", {"name": m})
        map_element.text = beast_maps[m]
        root_element.append(map_element)

    run_element = Element("run", {"chainLength": str(chain_length), "id": "mcmc", "spec": "MCMC", "storeEvery": str(store_every)})
    root_element.append(run_element)

    state_element = Element("state", {"id": "state", "storeEvery": str(store_every)})
    init_element = Element("init", {"birthRate": "@birthRate.t:Species", "estimate": "false", "id": "SBI", "popMean": "@popMean", "spec": "beast.evolution.speciation.StarBeastStartState", "speciesTree": "@Tree.t:Species"})
    distribution_element = Element("distribution", {"id": "posterior", "spec": "util.CompoundDistribution"})
    beast_logger_element = Element("logger", {"fileName": "beast_" + tree_name + ".log", "id": "tracelog", "logEvery": str(log_every), "model": "@posterior", "sort": "smart"})
    species_tree_logger_element = Element("logger", {"fileName": "species_" + tree_name + ".trees", "id": "speciesTreeLogger", "logEvery": str(log_every), "mode": "tree"})
    screen_logger_element = Element("logger", {"id": "screenlog", "logEvery": str(log_every), "model": "@posterior"})

    run_element.append(state_element)
    run_element.append(init_element)
    run_element.append(distribution_element)
    run_element.append(beast_logger_element)
    run_element.append(species_tree_logger_element)
    run_element.append(screen_logger_element)

    species_tree_element = Element("tree", {"id": "Tree.t:Species", "name": "stateNode"})
    species_tree_birthrate_element = Element("parameter", {"id": "birthRate.t:Species", "name": "stateNode", "lower": "0.0", "upper": "10000.0"})
    species_tree_birthrate_element.text = "213.2021"

    taxonset_element = Element("taxonset", {"id": "taxonsuperset", "spec": "TaxonSet"})
    species_tree_element.append(taxonset_element)
    state_element.append(species_tree_element)
    state_element.append(species_tree_birthrate_element)

    popmean_element = Element("parameter", {"id": "popMean", "name": "stateNode"})
    popsize_element = Element("parameter", {"id": "popSize", "name": "stateNode"})
    popsizetop_element = Element("parameter", {"id": "popSizeTop", "name": "stateNode"})

    popmean_element.text = "1.0"
    popsize_element.text = "1.0"
    popsizetop_element.text = "1.0"

    state_element.append(popmean_element)
    state_element.append(popsize_element)
    state_element.append(popsizetop_element)

    species_tree_popsize_element = Element("speciesTreePrior", {"bottomPopSize": "@popSize", "gammaParameter": "@popMean", "id": "SpeciesTreePopSize.Species", "popFunction": "linear_with_constant_root", "spec": "beast.evolution.speciation.SpeciesTreePrior", "taxonset": "@taxonsuperset", "topPopSize": "@popSizeTop", "tree": "@Tree.t:Species"})
    init_element.append(species_tree_popsize_element)

    species_coalescent_element = Element("distribution", {"id": "speciescoalescent", "spec": "util.CompoundDistribution"})
    prior_element = Element("distribution", {"id": "prior", "spec": "util.CompoundDistribution"})
    likelihood_element = Element("distribution", {"id": "likelihood", "spec": "util.CompoundDistribution"})

    distribution_element.append(species_coalescent_element)
    distribution_element.append(prior_element)
    distribution_element.append(likelihood_element)

    species_tree_popsize_element = Element("distribution", {"idref": "SpeciesTreePopSize.Species"})
    species_coalescent_element.append(species_tree_popsize_element)

    yule_birthrate_element = Element("distribution", {"birthDiffRate": "@birthRate.t:Species", "id": "YuleModel.t:Species", "spec": "beast.evolution.speciation.YuleModel", "tree": "@Tree.t:Species"})
    yule_birthrate_prior_element = Element("prior", {"id": "YuleBirthRatePrior.t:Species", "name": "distribution", "x": "@birthRate.t:Species"})
    yule_birthrate_oneonx_element = Element("OneOnX", {"id": "OneOnX.0", "name": "distr"})

    prior_element.append(yule_birthrate_element)
    prior_element.append(yule_birthrate_prior_element)
    yule_birthrate_prior_element.append(yule_birthrate_oneonx_element)

    popmean_prior_element = Element("prior", {"id": "popMean.prior", "name": "distribution", "x": "@popMean"})
    popmean_oneonx_element = Element("OneOnX", {"id": "OneOnX.01", "name": "distr"})

    prior_element.append(popmean_prior_element)
    popmean_prior_element.append(popmean_oneonx_element)

    reheight_operator_element = Element("operator", {"id": "Reheight.t:Species", "spec": "NodeReheight", "taxonset": "@taxonsuperset", "tree": "@Tree.t:Species", "weight": get_weight("reheight", st_gt_multiplier, False)})
    popsizetop_operator_element = Element("operator", {"degreesOfFreedom": "1", "id": "popSizeTopScaler.t:Species", "parameter": "@popSizeTop", "scaleFactor": "0.5", "spec": "ScaleOperator", "weight": get_weight("popsizetop", st_gt_multiplier, False)})
    popsizebottom_operator_element = Element("operator", {"degreesOfFreedom": "1", "id": "popSizeBottomScaler.t:Species", "parameter": "@popSize", "scaleFactor": "0.5", "spec": "ScaleOperator", "weight": get_weight("popsizebottom", st_gt_multiplier, False)})
    all_updown_operator_element = Element("operator", {"id": "updown.all.Species", "scaleFactor": "0.75", "spec": "UpDownOperator", "weight": get_weight("updown", st_gt_multiplier, False)})
    yule_operator_element = Element("operator", {"id": "YuleBirthRateScaler.t:Species", "parameter": "@birthRate.t:Species", "scaleFactor": "0.75", "spec": "ScaleOperator", "weight": get_weight("yulebirthrate", st_gt_multiplier, False)})
    popmean_operator_element = Element("operator", {"id": "popMeanScale.t:Species", "parameter": "@popMean", "scaleFactor": "0.75", "spec": "ScaleOperator", "weight": get_weight("popmean", st_gt_multiplier, False)})

    run_element.append(reheight_operator_element)
    run_element.append(popsizetop_operator_element)
    run_element.append(popsizebottom_operator_element)
    run_element.append(all_updown_operator_element)
    run_element.append(yule_operator_element)
    run_element.append(popmean_operator_element)

    all_updown_operator_element.append(Element("parameter", {"idref": "birthRate.t:Species", "name": "up"}))
    all_updown_operator_element.append(Element("parameter", {"idref": "popMean", "name": "down"}))
    all_updown_operator_element.append(Element("parameter", {"idref": "popSize", "name": "down"}))
    all_updown_operator_element.append(Element("parameter", {"idref": "popSizeTop", "name": "down"}))
    all_updown_operator_element.append(Element("tree", {"idref": "Tree.t:Species", "name": "down"}))

    beast_logger_element.append(Element("log", {"idref": "posterior"}))
    beast_logger_element.append(Element("log", {"idref": "likelihood"}))
    beast_logger_element.append(Element("log", {"idref": "prior"}))
    beast_logger_element.append(Element("log", {"idref": "speciescoalescent"}))
    beast_logger_element.append(Element("parameter", {"idref": "birthRate.t:Species", "name": "log"}))
    beast_logger_element.append(Element("log", {"idref": "YuleModel.t:Species"}))
    beast_logger_element.append(Element("log", {"id": "TreeHeight.Species", "spec": "beast.evolution.tree.TreeHeightLogger", "tree": "@Tree.t:Species"}))

    log_element = Element("log", {"id": "SpeciesTreeLoggerX", "popSize": "@popSize", "popSizeTop": "@popSizeTop", "spec": "beast.evolution.speciation.SpeciesTreeLogger", "speciesTreePrior": "@SpeciesTreePopSize.Species", "tree": "@Tree.t:Species"})
    treetop_element = Element("treetop", {"id": "treeTopFinder", "spec": "beast.evolution.speciation.TreeTopFinder"})

    species_tree_logger_element.append(log_element)
    log_element.append(treetop_element)

    screen_logger_element.append(Element("log", {"idref": "posterior"}))
    screen_logger_element.append(Element("log", {"arg": "@posterior", "id": "ESS.0", "spec": "util.ESS"}))
    screen_logger_element.append(Element("log", {"idref": "likelihood"}))
    screen_logger_element.append(Element("log", {"idref": "prior"}))

    n = 0
    for locus in sorted(species_tree):
        n += 1
        locus_element = Element("data", {"id": locus, "name": "alignment"})
        root_element.append(locus_element)

        gene_tree_element = Element("tree", {"id": "Tree.t:" + locus, "name": "stateNode"})
        gene_tree_taxonset_element = Element("taxonset", {"id": "TaxonSet." + locus, "spec": "TaxonSet"})
        gene_tree_data_element = Element("data", {"idref": locus, "name": "alignment"})
        gene_tree_kappa_element = Element("parameter", {"id": "kappa.s:" + locus, "lower": "0.0", "name": "stateNode"})
        gene_tree_kappa_element.text = "2.0"

        state_element.append(gene_tree_element)
        state_element.append(gene_tree_kappa_element)
        gene_tree_element.append(gene_tree_taxonset_element)
        gene_tree_taxonset_element.append(gene_tree_data_element)

        gene_element = Element("tree", {"idref": "Tree.t:" + locus, "name": "gene"})
        init_element.append(gene_element)

        gene_tree_prior_element = Element("distribution", {"id": "treePrior.t:" + locus, "spec": "beast.evolution.speciation.GeneTreeForSpeciesTreeDistribution", "speciesTree": "@Tree.t:Species", "speciesTreePrior": "@SpeciesTreePopSize.Species", "tree": "@Tree.t:" + locus})
        species_coalescent_element.append(gene_tree_prior_element)

        tree_likelihood_element = Element("distribution", {"data": "@" + locus, "id": "treeLikelihood." + locus, "spec": "TreeLikelihood", "tree": "@Tree.t:" + locus})
        site_model_element = Element("siteModel", {"id": "SiteModel.s:" + locus, "spec": "SiteModel"})
        mutation_rate_element = Element("parameter", {"estimate": "false", "id": "mutationRate.s:" + locus, "name": "mutationRate"})
        gamma_shape_element = Element("parameter", {"estimate": "false", "id": "gammaShape.s:" + locus, "name": "shape"})
        proportion_invariant_element = Element("parameter", {"estimate": "false", "id": "proportionInvariant.s:" + locus, "lower": "0.0", "name": "proportionInvariant", "upper": "1.0"})
        hky_element = Element("substModel", {"id": "hky.s:" + locus, "kappa": "@kappa.s:" + locus, "spec": "HKY"})
        empirical_freqs_element = Element("frequencies", {"data": "@" + locus, "id": "empiricalFreqs.s:" + locus, "spec": "Frequencies"})

        mutation_rate_element.text = "1.0"
        gamma_shape_element.text = "1.0"
        proportion_invariant_element.text = "0.0"

        likelihood_element.append(tree_likelihood_element)
        tree_likelihood_element.append(site_model_element)
        site_model_element.append(mutation_rate_element)
        site_model_element.append(gamma_shape_element)
        site_model_element.append(proportion_invariant_element)
        site_model_element.append(hky_element)
        hky_element.append(empirical_freqs_element)

        reheight_operator_element.append(Element("tree", {"idref": "Tree.t:" + locus, "name": "genetree"}))
        all_updown_operator_element.append(Element("tree", {"idref": "Tree.t:" + locus, "name": "down"}))

        treescaler_operator_element = Element("operator", {"id": "treeScaler.t:" + locus, "scaleFactor": "0.5", "spec": "ScaleOperator", "tree": "@Tree.t:" + locus, "weight": get_weight("tree", 1.0)})
        treerootscaler_operator_element = Element("operator", {"id": "treeRootScaler.t:" + locus, "rootOnly": "true", "scaleFactor": "0.5", "spec": "ScaleOperator", "tree": "@Tree.t:" + locus, "weight": get_weight("treeroot", 1.0)})
        uniform_operator_element = Element("operator", {"id": "UniformOperator.t:" + locus, "spec": "Uniform", "tree": "@Tree.t:" + locus, "weight": get_weight("uniformoperator", 1.0)})
        subtreeslide_operator_element = Element("operator", {"id": "SubtreeSlide.t:" + locus, "spec": "SubtreeSlide", "tree": "@Tree.t:" + locus, "weight": get_weight("subtreeslide", 1.0)})
        narrow_operator_element = Element("operator", {"id": "narrow.t:" + locus, "spec": "Exchange", "tree": "@Tree.t:" + locus, "weight": get_weight("narrow", 1.0)})
        wide_operator_element = Element("operator", {"id": "wide.t:" + locus, "isNarrow": "false", "spec": "Exchange", "tree": "@Tree.t:" + locus, "weight": get_weight("wide", 1.0)})
        wilsonbalding_operator_element = Element("operator", {"id": "WilsonBalding.t:" + locus, "spec": "WilsonBalding", "tree": "@Tree.t:" + locus, "weight": get_weight("wilsonbalding", 1.0)})
        kappa_operator_element = Element("operator", {"id": "KappaScaler.s:" + locus, "parameter": "@kappa.s:" + locus, "scaleFactor": "0.5", "spec": "ScaleOperator", "weight": get_weight("kappa", 1.0)})

        run_element.append(treescaler_operator_element)
        run_element.append(treerootscaler_operator_element)
        run_element.append(uniform_operator_element)
        run_element.append(subtreeslide_operator_element)
        run_element.append(narrow_operator_element)
        run_element.append(wide_operator_element)
        run_element.append(wilsonbalding_operator_element)
        run_element.append(kappa_operator_element)

        beast_logger_element.append(Element("log", {"idref": "treeLikelihood." + locus}))
        beast_logger_element.append(Element("log", {"idref": "treePrior.t:" + locus}))
        beast_logger_element.append(Element("log", {"id": "TreeHeight.t:" + locus, "spec": "beast.evolution.tree.TreeHeightLogger", "tree": "@Tree.t:" + locus}))
        beast_logger_element.append(Element("parameter", {"idref": "kappa.s:" + locus, "name": "log"}))

        treetop_element.append(Element("tree", {"idref": "Tree.t:" + locus}))

        gene_tree_logger_element = Element("logger", {"fileName": "$(tree)." + tree_name + ".trees", "id": "treelog.t:" + locus, "logEvery": str(log_every), "mode": "tree"})
        gene_tree_logger_element.append(Element("log", {"id": "TreeWithMetaDataLogger.t:" + locus, "spec": "beast.evolution.tree.TreeWithMetaDataLogger",  "tree": "@Tree.t:" + locus}))
        # Don't log gene trees... wastes way too much space
        # run_element.append(gene_tree_logger_element)

        if n == 1:
            kappa_prior_element = Element("prior", {"id": "KappaPrior.s:" + locus, "name": "distribution", "x": "@kappa.s:" + locus})
            kappa_lognormal_element = Element("LogNormal", {"id": "LogNormalDistributionModel.0", "name": "distr"})
            kappa_m_element = Element("parameter", {"estimate": "false", "id":"RealParameter.0", "name": "M"})
            kappa_s_element = Element("parameter", {"estimate": "false", "id":"RealParameter.01", "name": "S"})
            branch_rate_element = Element("branchRateModel", {"id": "StrictClock.c:" + locus, "spec": "beast.evolution.branchratemodel.StrictClockModel"})
            clock_rate_element = Element("parameter", {"estimate": "false", "id": "clockRate.c:" + locus, "name": "clock.rate"})

            kappa_m_element.text = "1.0"
            kappa_s_element.text = "1.25"
            clock_rate_element.text = "1.0"

            prior_element.append(kappa_prior_element)
            kappa_prior_element.append(kappa_lognormal_element)
            kappa_lognormal_element.append(kappa_m_element)
            kappa_lognormal_element.append(kappa_s_element)
            tree_likelihood_element.append(branch_rate_element)
            branch_rate_element.append(clock_rate_element)

        else:
            gene_tree_clockrate_element = Element("parameter", {"id": "clockRate.c:" + locus, "name": "stateNode"})
            kappa_prior_element = Element("prior", {"distr": "@LogNormalDistributionModel.0", "id": "KappaPrior.s:" + locus, "name": "distribution", "x": "@kappa.s:" + locus})
            clockrate_prior_element = Element("prior", {"id": "ClockPrior.c:" + locus, "name": "distribution", "x": "@clockRate.c:" + locus})
            clockrate_uniform_element = Element("Uniform", {"id": "Uniform.%02d" % (n - 2), "name": "distr", "upper": "Infinity"})
            gt_clockrate_operator_element = Element("operator", {"id": "StrictClockRateScaler.c:" + locus, "parameter": "@clockRate.c:" + locus, "scaleFactor": "0.75", "spec": "ScaleOperator", "weight": get_weight("strictclockrate", 1.0)})
            gt_updown_operator_element = Element("operator", {"id": "updown." + locus, "scaleFactor": "0.75", "spec": "UpDownOperator", "weight": get_weight("updown", 1.0)})
            clock_updown_operator_element = Element("operator", {"id": "strictClockUpDownOperator.c:" + locus, "scaleFactor": "0.75", "spec": "UpDownOperator", "weight": get_weight("strictclockupdownoperator", 1.0)})
            gene_tree_clockrate_element.text = "1.0"

            state_element.append(gene_tree_clockrate_element)
            prior_element.append(kappa_prior_element)
            prior_element.append(clockrate_prior_element)
            clockrate_prior_element.append(clockrate_uniform_element)
            run_element.append(gt_clockrate_operator_element)
            run_element.append(gt_updown_operator_element)
            run_element.append(clock_updown_operator_element)

            tree_likelihood_element.append(Element("branchRateModel", {"clock.rate": "@clockRate.c:" + locus, "id": "StrictClock.c:" + locus, "spec": "beast.evolution.branchratemodel.StrictClockModel"}))
            all_updown_operator_element.append(Element("parameter", {"idref": "clockRate.c:" + locus, "name": "up"}))
            gt_updown_operator_element.append(Element("parameter", {"idref": "clockRate.c:" + locus, "name": "up"}))
            gt_updown_operator_element.append(Element("tree", {"idref": "Tree.t:" + locus, "name": "down"}))
            clock_updown_operator_element.append(Element("parameter", {"idref": "clockRate.c:" + locus, "name": "up"}))
            clock_updown_operator_element.append(Element("tree", {"idref": "Tree.t:" + locus, "name": "down"}))

            beast_logger_element.append(Element("parameter", {"idref": "clockRate.c:" + locus, "name": "log"}))

        for species in sorted(species_tree[locus]):
            if n == 1:
                taxon_element = Element("taxon", {"id": species, "spec": "TaxonSet"})
                taxonset_element.append(taxon_element)

            for individual in species_tree[locus][species]:
                sample_id = species + individual

                if n == 1:
                    taxon_element.append(Element("taxon", {"id": sample_id, "spec": "Taxon"}))

                sequence_id = sample_id + locus
                sequence = species_tree[locus][species][individual]
                sequence_element = Element("sequence", {"id": sequence_id, "taxon": sample_id, "totalcount": "4", "value": sequence})
                locus_element.append(sequence_element)

    indent(root_element)
    return beast_xml
