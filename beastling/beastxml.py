# -*- encoding: utf-8 -*-
import datetime
import itertools
from math import log, exp
import sys
import xml.etree.ElementTree as ET

from six import BytesIO, PY3

from clldutils.path import Path

from beastling import __version__
import beastling.beast_maps as beast_maps

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

class BeastXml(object):

    def __init__(self, config):
        self.config = config
        if not self.config.processed:
            self.config.process()
        self._covarion_userdatatype_created = False
        # Tell everybody about ourselves
        for model in self.config.all_models:
            model.beastxml = self
        for clock in self.config.clocks:
            clock.beastxml = self
        self._taxon_sets = {}
        self.build_xml()

    def build_xml(self):
        """
        Creates a complete BEAST XML configuration file as an ElementTree,
        descending from the self.beast element.
        """
        attribs = {}
        attribs["beautitemplate"] = "Standard"
        attribs["beautistatus"] = ""
        attribs["namespace"] = "beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood"
        attribs["version"] ="2.0"
        self.beast = ET.Element("beast", attrib=attribs)
        self.add_beastling_comment()
        self.embed_data()
        self.add_maps()
        for model in self.config.models:
            model.add_master_data(self.beast)
            model.add_misc(self.beast)
        for clock in self.config.clocks:
            clock.add_branchrate_model(self.beast)
        self.add_run()

    def add_beastling_comment(self):
        """
        Add a comment at the root level of the XML document indicating the
        BEASTling version used to create the file, the time and date of
        generation and the original configuration file text.
        """
        comment_lines = []
        comment_lines.append("Generated by BEASTling %s on %s.\n" % (__version__,datetime.datetime.now().strftime("%A, %d %b %Y %I:%M %p")))
        if self.config.configfile:
            comment_lines.append("Original config file:\n")
            comment_lines.append(self.config.configfile.write_string())
            comment_lines.append("Please DO NOT manually edit this file without removing this message or editing")
            comment_lines.append("it to describe the changes made.  Otherwise attempts to replicate your")
            comment_lines.append("analysis using BEASTling and the above configuration may not be valid.\n")
        else:
            comment_lines.append("Configuration built programmatically.")
            comment_lines.append("No config file to include.")
        self.beastling_comment = ET.Comment("\n".join(comment_lines))
        self.beast.append(self.beastling_comment)

    def embed_data(self):
        """
        Embed a copy of each data file in a comment at the top of the XML
        document.
        """
        if not self.config.embed_data:
            return
        for filename in self.config.files_to_embed:
            self.beast.append(self.format_data_file(filename))
        for model in self.config.models:
            self.beast.append(self.format_data_file(model.data_filename))

    def format_data_file(self, filename):
        """
        Return an ElementTree node corresponding to a comment containing
        the text of the specified data file.
        """
        header = "BEASTling embedded data file: %s" % filename
        fp = Path(filename).open("r")
        data_block = "\n".join([header, fp.read()])
        fp.close()
        return ET.Comment(data_block)

    def add_maps(self):
        """
        Add <map> elements aliasing common BEAST classes.
        """
        for a, b in beast_maps.maps:
            mapp = ET.SubElement(self.beast, "map", attrib={"name":a})
            mapp.text = b

    def add_run(self):
        """
        Add the <run> element and all its descendants, which is most of the
        analysis.
        """
        if self.config.path_sampling:
            self.add_path_sampling_run()
        else:
            self.add_standard_sampling_run()
        self.estimate_tree_height()
        self.add_state()
        self.add_init()
        self.add_distributions()
        self.add_operators()
        self.add_loggers()

    def add_standard_sampling_run(self):
        """
        Add the <run> element (only) for a standard analysis, i.e. without
        path sampling.  The <state>, <init> etc. are added to whatever this
        method names self.run.
        """
        attribs = {}
        attribs["id"] = "mcmc"
        attribs["spec"] = "MCMC"
        attribs["chainLength"] = str(self.config.chainlength)
        attribs["numInitializationAttempts"] = "1000"
        if self.config.sample_from_prior:
            attribs["sampleFromPrior"] = "true"
        self.run = ET.SubElement(self.beast, "run", attrib=attribs)

    def add_path_sampling_run(self):
        """
        Add the <run> element (only) for a path sampling analysis.  We call
        this self.ps_run and assign the nested <mcmc> element to self.run,
        so that <state>, <init> etc. will be correctly added there.
        """
        attribs = {
            "id": "ps",
            "spec": "beast.inference.PathSampler",
            "chainLength": str(self.config.chainlength),
            "nrOfSteps": str(self.config.steps),
            "alpha": str(self.config.alpha),
            "rootdir": self.config.basename+"_path_sampling",
            "preBurnin": str(int((self.config.preburnin/100)*self.config.chainlength)),
            "burnInPercentage": str(self.config.log_burnin),
            "deleteOldLogs": "true",
            }
        if self.config.do_not_run:
            attribs["doNotRun"] = "true"
        self.ps_run = ET.SubElement(self.beast, "run", attribs)
        self.ps_run.text = """cd $(dir)
java -cp $(java.class.path) beast.app.beastapp.BeastMain $(resume/overwrite) -java -seed $(seed) beast.xml"""

        attribs = {}
        attribs["id"] = "mcmc"
        attribs["spec"] = "MCMC"
        attribs["chainLength"] = str(self.config.chainlength)
        self.run = ET.SubElement(self.ps_run, "mcmc", attrib=attribs)

    def add_state(self):
        """
        Add the <state> element and all its descendants.
        """
        self.state = ET.SubElement(self.run, "state", {"id":"state","storeEvery":"5000"})
        self.add_tree_state()
        for clock in self.config.clocks:
            clock.add_state(self.state)
        for model in self.config.all_models:
            model.add_state(self.state)

    def add_tip_heights(self):
        string_bits = []
        for cal in self.config.tip_calibrations.values():
            initial_height = cal.mean()
            string_bits.append("{:s} = {:}".format(next(cal.langs.__iter__()), initial_height))
        trait_string = ",\n".join(string_bits)

        datetrait = ET.SubElement(self.tree, "trait",
                      {"id": "datetrait",
                       "spec": "beast.evolution.tree.TraitSet",
                       "taxa": "@taxa",
                       "traitname": "date-backward"})
        datetrait.text = trait_string

    def add_tree_state(self):
        """
        Add tree-related <state> sub-elements.
        """
        self.tree = ET.SubElement(self.state, "tree", {"id":"Tree.t:beastlingTree", "name":"stateNode"})
        self.add_taxon_set(self.tree, "taxa", self.config.languages, define_taxa=True)
        if self.config.tree_prior in ["yule", "birthdeath"]:
            param = ET.SubElement(self.state, "parameter", {"id":"birthRate.t:beastlingTree","name":"stateNode"})
            if self.birthrate_estimate is not None:
                param.text=str(self.birthrate_estimate)
            else:
                param.text="1.0"
            if self.config.tree_prior in ["birthdeath"]:
                ET.SubElement(self.state, "parameter",
                              {"id": "deathRate.t:beastlingTree",
                               "name": "stateNode"}).text = "0.5"
                ET.SubElement(self.state, "parameter",
                              {"id": "sampling.t:beastlingTree",
                               "name": "stateNode"}).text = "0.2"

        elif self.config.tree_prior == "coalescent":
            param = ET.SubElement(self.state, "parameter", {"id":"popSize.t:beastlingTree","name":"stateNode"})
            param.text="1.0"
        if self.config.tip_calibrations:
            self.add_tip_heights()

    def add_init(self):
        """
        Add the <init> element and all its descendants.
        """

        # If a starting tree is specified, use it...
        if self.config.starting_tree:
            self.init = ET.SubElement(self.run, "init", {"estimate":"false", "id":"startingTree", "initial":"@Tree.t:beastlingTree", "spec":"beast.util.TreeParser","IsLabelledNewick":"true", "newick":self.config.starting_tree})
        # ...if not, use the simplest random tree initialiser possible
        else:
            # If we have non-trivial monophyly constraints, use ConstrainedRandomTree
            if self.config.monophyly and len(self.config.languages) > 2:
                self.add_constrainedrandomtree_init()
            # If we have hard-bound calibrations, use SimpleRandomTree
            elif any([c.dist == "uniform" for c in self.config.calibrations.values()]):
                self.add_simplerandomtree_init()
            # Otherwise, just use RandomTree
            else:
                self.add_randomtree_init()

    def estimate_tree_height(self):
        """
        Make a rough estimate of what the starting height of the tree should
        be so we can initialise somewhere decent.
        """
        birthrate_estimates = []
        for cal in self.config.calibrations.values():
            if len(cal.langs) == 1 or cal.dist not in ("normal", "lognormal"):
                continue
            # Find the midpoint of this cal
            mid = cal.mean()
            # Find the Yule birthrate which results in an expected height for
            # a tree of this many taxa which equals the midpoint of the
            # calibration.
            # The expected height of a Yule tree with n taxa and
            # birthrate λ is 1/λ * (Hn - 1), where Hn is the nth
            # harmonic number.  Hn can be asymptotically approximated
            # by Hn = log(n) + 0.5772156649. So λ = (Hn - 1) / h.
            birthrate = (log(len(cal.langs)) + 0.5772156649 - 1) / mid
            birthrate_estimates.append(birthrate)
        # If there were no calibrations that could be used, return a non-esitmate
        if not birthrate_estimates:
            self.birthrate_estimate = None
            self.treeheight_estimate = None
            return
        # Find the mean birthrate estimate
        self.birthrate_estimate = round(sum(birthrate_estimates) / len(birthrate_estimates), 4)
        # Find the expected height of a tree with this birthrate
        self.treeheight_estimate = round((1.0/self.birthrate_estimate)*(log(len(self.config.languages)) + 0.5772156649 - 1), 4)

    def add_randomtree_init(self):
        attribs = {"estimate":"false", "id":"startingTree", "initial":"@Tree.t:beastlingTree", "taxonset":"@taxa", "spec":"beast.evolution.tree.RandomTree"}
        if self.birthrate_estimate is not None:
            attribs["rootHeight"] = str(self.treeheight_estimate)
        self.init = ET.SubElement(self.run, "init", attribs)
        popmod = ET.SubElement(self.init, "populationModel", {"spec":"ConstantPopulation"})
        ET.SubElement(popmod, "popSize", {"spec":"parameter.RealParameter","value":"1"})

    def add_simplerandomtree_init(self):
        attribs = {"estimate":"false", "id":"startingTree", "initial":"@Tree.t:beastlingTree", "taxonset":"@taxa", "spec":"beast.evolution.tree.SimpleRandomTree"}
        if self.birthrate_estimate is not None:
            attribs["rootHeight"] = str(self.treeheight_estimate)
        self.init = ET.SubElement(self.run, "init", attribs)

    def add_constrainedrandomtree_init(self):
        attribs = {"estimate":"false", "id":"startingTree", "initial":"@Tree.t:beastlingTree", "taxonset":"@taxa", "spec":"beast.evolution.tree.ConstrainedRandomTree", "constraints":"@constraints"}
        if self.birthrate_estimate is not None:
            attribs["rootHeight"] = str(self.treeheight_estimate)
        self.init = ET.SubElement(self.run, "init", attribs)
        popmod = ET.SubElement(self.init, "populationModel", {"spec":"ConstantPopulation"})
        ET.SubElement(popmod, "popSize", {"spec":"parameter.RealParameter","value":"1"})

    def add_distributions(self):
        """
        Add all probability distributions under the <run> element.
        """
        self.posterior = ET.SubElement(self.run,"distribution",{"id":"posterior","spec":"util.CompoundDistribution"})
        self.add_prior()
        self.add_likelihood()

    def add_prior(self):
        """
        Add all prior distribution elements.
        """
        self.prior = ET.SubElement(self.posterior,"distribution",{"id":"prior","spec":"util.CompoundDistribution"})
        self.add_monophyly_constraints()
        self.add_calibrations()
        self.add_tree_prior()
        for clock in self.config.clocks:
            clock.add_prior(self.prior)
        for model in self.config.all_models:
            model.add_prior(self.prior)

    def add_monophyly_constraints(self):
        """
        Add monophyly constraints to prior distribution.
        """
        if self.config.monophyly:
            attribs = {}
            attribs["id"] = "constraints"
            attribs["spec"] = "beast.math.distributions.MultiMonophyleticConstraint"
            attribs["tree"] = "@Tree.t:beastlingTree"
            attribs["newick"] = self.config.monophyly_newick
            ET.SubElement(self.prior, "distribution", attribs)

    def add_calibrations(self):
        """
        Add timing calibrations to prior distribution.
        """
        # This itertools.cchain is a bit ugly, I wonder if we can get away with sticking them all in one list...
        for clade, cal in sorted(itertools.chain(self.config.calibrations.items(), self.config.tip_calibrations.items())):
            # Don't add an MRCA cal for point calibrations, those only exist to
            # cause the initial tip height to be set
            if cal.dist == "point":
                continue
            # BEAST's logcombiner chokes on spaces...
            clade = clade.replace(" ","_")
            # Create MRCAPrior node
            attribs = {}
            attribs["id"] = clade + "MRCA"
            attribs["monophyletic"] = "true"
            attribs["spec"] = "beast.math.distributions.MRCAPrior"
            attribs["tree"] = "@Tree.t:beastlingTree"
            if cal.originate:
                attribs["useOriginate"] = "true"
            elif len(cal.langs) == 1:   # If there's only 1 lang and it's not an originate cal, it must be a tip cal
                attribs["tipsonly"] = "true"

            cal_prior = ET.SubElement(self.prior, "distribution", attribs)

            # Create "taxonset" param for MRCAPrior
            taxonsetname = clade[:-len("_originate")] if clade.endswith("_originate") else clade
            self.add_taxon_set(cal_prior, taxonsetname, cal.langs)

            cal.generate_xml_element(cal_prior)
            
    def add_taxon_set(self, parent, label, langs, define_taxa=False):
        """
        Add a TaxonSet element with the specified set of languages.

        If a TaxonSet previously defined by this method contains exactly the
        same set of taxa, a reference to that TaxonSet will be added instead.
        By default, each TaxonSet will contain references to the taxa,
        assuming that they have been defined previously (most probably in the
        definition of the tree).  If this is not the case, passing
        define_taxa=True will define, rather than refer to, the taxa.
        """
        # Kill duplicates
        langs = sorted(list(set(langs)))

        # If we've been asked to build an emtpy TaxonSet, something is very wrong,
        # so better to die loud and early
        assert(langs)
        # Refer to any previous TaxonSet with the same languages
        for idref, taxa in self._taxon_sets.items():
            if langs == taxa:
                ET.SubElement(parent, "taxonset", {"idref" : idref})
                return
        if len(langs) == 1 and label == langs[0]:
            # Single taxa are IDs already. They cannot also be taxon set ids.
            label = "tx_{:}".format(label)
        # Otherwise, create and register a new TaxonSet
        taxonset = ET.SubElement(parent, "taxonset", {"id" : label, "spec":"TaxonSet"})
        ## If the taxonset is more than 3 languages in size, use plate notation to minimise XML filesize
        if len(langs) > 3:
            plate = ET.SubElement(taxonset, "plate", {
                "var":"language",
                "range":",".join(langs)})
            ET.SubElement(plate, "taxon", {"id" if define_taxa else "idref" :"$(language)"})
        ## Otherwise go for the more readable notation...
        else:
            for lang in langs:
                ET.SubElement(taxonset, "taxon", {"id" if define_taxa else "idref" : lang})
        self._taxon_sets[label] = langs

    def add_tree_prior(self):
        if self.config.tree_prior.lower() == "yule":
            self.add_yule_tree_prior()
        elif self.config.tree_prior.lower() == "birthdeath":
            self.add_birthdeath_tree_prior()
        elif self.config.tree_prior.lower() == "coalescent":
            self.add_coalescent_tree_prior()
        elif self.config.tree_prior.lower() == "uniform":
            pass
        else:
            raise ValueError("Tree prior {:} is unknown.".format(
                self.config.tree_prior.lower()))

    def add_birthdeath_tree_prior(self):
        """Add a (calibrated) birth-death tree prior."""
        # Tree prior

        attribs = {}
        attribs["id"] = "BirthDeathModel.t:beastlingTree"
        attribs["tree"] = "@Tree.t:beastlingTree"
        attribs["spec"] = "beast.evolution.speciation.BirthDeathGernhard08Model"
        attribs["birthRate"] = "@birthRate.t:beastlingTree"
        attribs["relativeDeathRate"] = "@deathRate.t:beastlingTree"
        attribs["sampleProbability"] = "@sampling.t:beastlingTree"
        attribs["type"] = "restricted"

        # Birth rate prior
        attribs = {}
        attribs["id"] = "BirthRatePrior.t:beastlingTree"
        attribs["name"] = "distribution"
        attribs["x"] = "@birthRate.t:beastlingTree"
        sub_prior = ET.SubElement(self.prior, "prior", attribs)
        uniform = ET.SubElement(sub_prior, "Uniform",
                                {"id": "Uniform.0",
                                 "name": "distr",
                                 "upper": "Infinity"})

        # Relative death rate prior
        attribs = {}
        attribs["id"] = "relativeDeathRatePrior.t:beastlingTree"
        attribs["name"] = "distribution"
        attribs["x"] = "@deathRate.t:beastlingTree"
        sub_prior = ET.SubElement(self.prior, "prior", attribs)
        uniform = ET.SubElement(sub_prior, "Uniform",
                                {"id": "Uniform.1",
                                 "name": "distr",
                                 "upper": "Infinity"})

        # Sample probability prior
        attribs = {}
        attribs["id"] = "samplingPrior.t:beastlingTree"
        attribs["name"] = "distribution"
        attribs["x"] = "@sampling.t:beastlingTree"
        sub_prior = ET.SubElement(self.prior, "prior", attribs)
        uniform = ET.SubElement(sub_prior, "Uniform",
                                {"id": "Uniform.3",
                                 "name": "distr",
                                 "lower": "0",
                                 "upper": "1"})


    def add_yule_tree_prior(self):
        """
        Add Yule birth-process tree prior.
        """
        # Tree prior
        ## Decide whether to use the standard Yule or the fancy calibrated one
        if len(self.config.calibrations) == 1:
            yule = "calibrated"
        elif len(self.config.calibrations) == 2:
            # Two calibrations can be handled by the calibrated Yule if they
            # are nested
            langs1, langs2 = [c.langs for c in self.config.calibrations.values()]
            if len(set(langs1) & set(langs2)) in (len(langs1), len(langs2)):
                yule = "calibrated"
            else:
                yule = "standard"
        else:
            yule = "standard"

        attribs = {}
        attribs["id"] = "YuleModel.t:beastlingTree"
        attribs["tree"] = "@Tree.t:beastlingTree"
        if yule == "standard":
            attribs["spec"] = "beast.evolution.speciation.YuleModel"
            attribs["birthDiffRate"] = "@birthRate.t:beastlingTree"
            if "root" in self.config.calibrations:
                attribs["conditionalOnRoot"] = "true"
        elif yule == "calibrated":
            attribs["spec"] = "beast.evolution.speciation.CalibratedYuleModel"
            attribs["birthRate"] = "@birthRate.t:beastlingTree"
        ET.SubElement(self.prior, "distribution", attribs)

        # Birth rate prior
        attribs = {}
        attribs["id"] = "YuleBirthRatePrior.t:beastlingTree"
        attribs["name"] = "distribution"
        attribs["x"] = "@birthRate.t:beastlingTree"
        sub_prior = ET.SubElement(self.prior, "prior", attribs)
        uniform = ET.SubElement(sub_prior, "Uniform", {"id":"Uniform.0","name":"distr","upper":"Infinity"})

    def add_coalescent_tree_prior(self):

        coalescent = ET.SubElement(self.prior, "distribution", {
            "id": "Coalescent.t:beastlingTree",
            "spec": "Coalescent",
            })
        popmod = ET.SubElement(coalescent, "populationModel", {
            "id": "ConstantPopulation:beastlingTree",
            "spec": "ConstantPopulation",
            })
        ET.SubElement(popmod, "parameter", {
            "idref": "popSize.t:beastlingTree",
            "name": "popSize",
            })
        ET.SubElement(coalescent, "treeIntervals", {
            "id": "TreeIntervals",
            "spec": "TreeIntervals",
            "tree": "@Tree.t:beastlingTree",
            })

    def add_likelihood(self):
        """
        Add all likelihood distribution elements.
        """
        self.likelihood = ET.SubElement(self.posterior,"distribution",{"id":"likelihood","spec":"util.CompoundDistribution"})
        for model in self.config.all_models:
            model.add_likelihood(self.likelihood)

    def add_operators(self):
        """
        Add all <operator> elements.
        """
        self.add_tree_operators()
        for clock in self.config.clocks:
            clock.add_operators(self.run)
        for model in self.config.all_models:
            model.add_operators(self.run)
        # Add one DeltaExchangeOperator for feature rates per clock
        for clock in self.config.clocks:
            clock_models = [m for m in self.config.models if m.rate_variation and m.clock == clock]
            if not clock_models:
                continue
            # Add one big DeltaExchangeOperator which operates on all
            # feature clock rates from all models
            delta = ET.SubElement(self.run, "operator", {"id":"featureClockRateDeltaExchanger:%s" % clock.name, "spec":"DeltaExchangeOperator", "weight":"3.0"})
            for model in clock_models:
                plate = ET.SubElement(delta, "plate", {
                    "var":"rate",
                    "range":",".join(model.all_rates)})
                ET.SubElement(plate, "parameter", {"idref":"featureClockRate:%s:$(rate)" % model.name})
            # Add weight vector if there has been any binarisation
            if any([w != 1 for w in itertools.chain(*[m.weights for m in clock_models])]):
                weightvector = ET.SubElement(delta, "weightvector", {
                    "id":"featureClockRateWeightParameter:%s" % clock.name,
                    "spec":"parameter.IntegerParameter",
                    "dimension":str(sum([len(m.weights) for m in clock_models])),
                    "estimate":"false"
                })
                weightvector.text = " ".join(itertools.chain(*[map(str, m.weights) for m in clock_models]))


    def add_tree_operators(self):
        """
        Add all <operator>s which act on the tree topology and branch lengths.
        """
        # Tree operators
        # Operators which affect the tree must respect the sample_topology and
        # sample_branch_length options.
        if self.config.sample_topology:
            ## Tree topology operators
            ET.SubElement(self.run, "operator", {"id":"SubtreeSlide.t:beastlingTree","spec":"SubtreeSlide","tree":"@Tree.t:beastlingTree","markclades":"true", "weight":"15.0"})
            ET.SubElement(self.run, "operator", {"id":"narrow.t:beastlingTree","spec":"Exchange","tree":"@Tree.t:beastlingTree","markclades":"true", "weight":"15.0"})
            ET.SubElement(self.run, "operator", {"id":"wide.t:beastlingTree","isNarrow":"false","spec":"Exchange","tree":"@Tree.t:beastlingTree","markclades":"true", "weight":"3.0"})
            ET.SubElement(self.run, "operator", {"id":"WilsonBalding.t:beastlingTree","spec":"WilsonBalding","tree":"@Tree.t:beastlingTree","markclades":"true","weight":"3.0"})
        if self.config.sample_branch_lengths:
            ## Branch length operators
            ET.SubElement(self.run, "operator", {"id":"UniformOperator.t:beastlingTree","spec":"Uniform","tree":"@Tree.t:beastlingTree","weight":"30.0"})
            ET.SubElement(self.run, "operator", {"id":"treeScaler.t:beastlingTree","scaleFactor":"0.5","spec":"ScaleOperator","tree":"@Tree.t:beastlingTree","weight":"3.0"})
            ET.SubElement(self.run, "operator", {"id":"treeRootScaler.t:beastlingTree","scaleFactor":"0.5","spec":"ScaleOperator","tree":"@Tree.t:beastlingTree","rootOnly":"true","weight":"3.0"})
            ## Up/down operator which scales tree height
            if self.config.tree_prior in ["yule", "birthdeath"]:
                updown = ET.SubElement(self.run, "operator", {"id":"UpDown","spec":"UpDownOperator","scaleFactor":"0.5", "weight":"3.0"})
                ET.SubElement(updown, "tree", {"idref":"Tree.t:beastlingTree", "name":"up"})
                ET.SubElement(updown, "parameter", {"idref":"birthRate.t:beastlingTree", "name":"down"})
                ### Include clock rates in up/down only if calibrations are given
                if self.config.calibrations:
                    for clock in self.config.clocks:
                        if clock.estimate_rate:
                            ET.SubElement(updown, "parameter", {"idref":clock.mean_rate_id, "name":"down"})

        if self.config.tree_prior in ["yule", "birthdeath"]:
            # Birth rate scaler
            # Birth rate is *always* scaled.
            ET.SubElement(self.run, "operator", {"id":"YuleBirthRateScaler.t:beastlingTree","spec":"ScaleOperator","parameter":"@birthRate.t:beastlingTree", "scaleFactor":"0.5", "weight":"3.0"})
        elif self.config.tree_prior == "coalescent":
            ET.SubElement(self.run, "operator", {"id":"PopulationSizeScaler.t:beastlingTree","spec":"ScaleOperator","parameter":"@popSize.t:beastlingTree", "scaleFactor":"0.5", "weight":"3.0"})

        if self.config.tree_prior in ["birthdeath"]:
            ET.SubElement(self.run, "operator",
                          {"id": "SamplingScaler.t:beastlingTree",
                           "spec": "ScaleOperator",
                           "parameter": "@sampling.t:beastlingTree",
                           "scaleFactor": "0.8",
                           "weight": "1.0"})
            ET.SubElement(self.run, "operator",
                          {"id": "DeathRateScaler.t:beastlingTree",
                           "spec": "ScaleOperator",
                           "parameter": "@deathRate.t:beastlingTree",
                           "scaleFactor": "0.5",
                           "weight": "3.0"})
 
        # Add a Tip Date scaling operator if required
        if self.config.tip_calibrations and self.config.sample_branch_lengths:
            # Get a list of taxa with non-point tip cals
            tip_taxa = [next(cal.langs.__iter__()) for cal in self.config.tip_calibrations.values() if cal.dist != "point"]
            for taxon in tip_taxa:
                tiprandomwalker = ET.SubElement(self.run, "operator",
                    {"id": "TipDatesandomWalker:%s" % taxon,
                     "spec": "TipDatesRandomWalker",
                     "windowSize": "1",
                     "tree": "@Tree.t:beastlingTree",
                     "weight": "3.0",
                     })
                self.add_taxon_set(tiprandomwalker, taxon, (taxon,))

    def add_loggers(self):
        """
        Add all <logger> elements.
        """
        self.add_screen_logger()
        self.add_tracer_logger()
        self.add_tree_loggers()

        # Log individual reconstructed traits (and possibly other per-generation metadata)
        if any([model.metadata for model in self.config.models]):
            self.add_trait_logger("_reconstructed")

    def add_screen_logger(self):
        """
        Add the screen logger, if configured to do so.
        """
        if not self.config.screenlog:
            return
        screen_logger = ET.SubElement(self.run, "logger", attrib={"id":"screenlog", "logEvery":str(self.config.log_every)})
        log = ET.SubElement(screen_logger, "log", attrib={"arg":"@posterior", "id":"ESS.0", "spec":"util.ESS"})
        log = ET.SubElement(screen_logger, "log", attrib={"idref":"prior"})
        log = ET.SubElement(screen_logger, "log", attrib={"idref":"likelihood"})
        log = ET.SubElement(screen_logger, "log", attrib={"idref":"posterior"})

    def add_tracer_logger(self):
        """
        Add file logger, if configured to do so.
        """
        if not(self.config.log_probabilities or self.config.log_params):
            return
        tracer_logger = ET.SubElement(self.run,"logger",{"id":"tracelog","fileName":self.config.basename+".log","logEvery":str(self.config.log_every),"sort":"smart"})
        # Log prior, likelihood and posterior
        if self.config.log_probabilities:
            ET.SubElement(tracer_logger,"log",{"idref":"prior"})
            ET.SubElement(tracer_logger,"log",{"idref":"likelihood"})
            ET.SubElement(tracer_logger,"log",{"idref":"posterior"})
        # Log Yule birth rate
        if self.config.log_params:
            if self.config.tree_prior in ["yule", "birthdeath"]:
                ET.SubElement(tracer_logger,"log",{"idref":"birthRate.t:beastlingTree"})
                if self.config.tree_prior in ["birthdeath"]:
                    ET.SubElement(tracer_logger, "log",
                                  {"idref": "deathRate.t:beastlingTree"})
                    ET.SubElement(tracer_logger, "log",
                                  {"idref": "sampling.t:beastlingTree"})
            elif self.config.tree_prior == "coalescent":
                ET.SubElement(tracer_logger,"log",{"idref":"popSize.t:beastlingTree"})
            for clock in self.config.clocks:
                clock.add_param_logs(tracer_logger)
            for model in self.config.all_models:
                model.add_param_logs(tracer_logger)

        # Log tree height
        if not self.config.tree_logging_pointless:
            ET.SubElement(tracer_logger,"log",{
                "id":"treeStats",
                "spec":"beast.evolution.tree.TreeStatLogger",
                "tree":"@Tree.t:beastlingTree"})

        # Log calibration clade heights
        for clade, cal in sorted(itertools.chain(self.config.calibrations.items(), self.config.tip_calibrations.items())):
            # Don't log unchanging tip heights
            if cal.dist == "point":
                continue
            clade = clade.replace(" ","_")
            ET.SubElement(tracer_logger,"log",{"idref":"%sMRCA" % clade})

        # Fine-grained logging
        if self.config.log_fine_probs:
            ET.SubElement(tracer_logger,"log",{"idref":"YuleModel.t:beastlingTree"})
            ET.SubElement(tracer_logger,"log",{"idref":"YuleBirthRatePrior.t:beastlingTree"})

    def add_tree_loggers(self):
        """
        Add tree logger, if configured to do so.
        """
        if not self.config.log_trees or self.config.tree_logging_pointless:
            return

        pure_tree_done = False
        non_strict_clocks = set([m.clock for m in self.config.models if not m.clock.is_strict])
        if not non_strict_clocks:
            # All clocks are strict, so we just do one pure log file
            self.add_tree_logger()
            pure_tree_done = True
        else:
            # There are non-strict clocks, so we do one log file each with branch rates
            for clock in non_strict_clocks:
                if len(non_strict_clocks) == 1:
                    self.add_tree_logger("", clock.branchrate_model_id)
                else:
                    self.add_tree_logger("_%s_rates" % clock.name, clock.branchrate_model_id)

        # If asked, do a topology-only tree log (i.e. no branch rates)
        if self.config.log_pure_tree and not pure_tree_done:
            self.add_tree_logger("_pure")

        # Log reconstructed traits (and possibly other per-node metadata)
        if any([model.treedata for model in self.config.models]):
            self.add_trait_tree_logger("_reconstructed")

        # Created a dedicated geographic tree log if asked to log locations,
        # or if the geo model's clock is non-strict
        if not self.config.geo_config:
            return
        if self.config.geo_config["log_locations"] or not self.config.geo_model.clock.is_strict:
            self.add_tree_logger("_geography", self.config.geo_model.clock.branchrate_model_id, True)

    def add_tree_logger(self, suffix="", branchrate_model_id=None, locations=False):
        tree_logger = ET.SubElement(self.run, "logger", {"mode":"tree", "fileName":self.config.basename + suffix + ".nex", "logEvery":str(self.config.log_every),"id":"treeLogger" + suffix})
        log = ET.SubElement(tree_logger, "log", attrib={"id":"TreeLoggerWithMetaData"+suffix,"spec":"beast.evolution.tree.TreeWithMetaDataLogger","tree":"@Tree.t:beastlingTree", "dp":str(self.config.log_dp)})
        if branchrate_model_id:
            ET.SubElement(log, "branchratemodel", {"idref":branchrate_model_id})
        if locations:
            ET.SubElement(log, "metadata", {
                "id":"location",
                "spec":"sphericalGeo.TraitFunction",
                "likelihood":"@sphericalGeographyLikelihood"}).text = "0.0"

    def add_trait_tree_logger(self, suffix=""):
        tree_logger = ET.SubElement(self.run, "logger", {"mode":"tree", "fileName":self.config.basename + suffix + ".nex", "logEvery":str(self.config.log_every),"id":"treeLogger" + suffix})
        log = ET.SubElement(tree_logger, "log", attrib={"id":"ReconstructedStateTreeLogger","spec":"beast.evolution.tree.TreeWithTraitLogger","tree":"@Tree.t:beastlingTree"})
        for model in self.config.models:
            for md in model.treedata:
                ET.SubElement(log, "metadata", {"idref": md})

    def add_trait_logger(self, suffix=""):
        """Add a logger referencing all AncestralStateLogger likelihoods in the tree."""
        trait_logger = ET.SubElement(self.run, "logger",
                                     {"fileName": self.config.basename + suffix + ".log",
                                      "logEvery": str(self.config.log_every),
                                      "id":"traitLogger" + suffix})
        for model in self.config.models:
            for reference in model.metadata:
                ET.SubElement(trait_logger, "log", {
                    "idref": reference})

    def tostring(self):
        """
        Return a string representation of the entire XML document.
        """
        out = BytesIO()
        self.write(out)
        out.seek(0)
        return out.read()

    def write(self, stream):
        indent(self.beast)
        tree = ET.ElementTree(self.beast)
        tree.write(stream, encoding='UTF-8', xml_declaration=True)

    def write_file(self, filename=None):
        """
        Write the XML document to a file.
        """
        if filename in ("stdout", "-"):
            # See https://docs.python.org/3/library/sys.html#sys.stdout
            self.write(getattr(sys.stdout, 'buffer', sys.stdout) if PY3 else sys.stdout)
        else:
            with open(filename or self.config.basename + ".xml", "wb") as stream:
                self.write(stream)
