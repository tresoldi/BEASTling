import datetime
import sys
import xml.etree.ElementTree as ET

from six import BytesIO, PY3

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
        for clock in self.config.clocks:
            clock.add_branchrate_model(self.beast)
        for model in self.config.all_models:
            model.add_misc(self.beast)
        self.add_run()

    def add_beastling_comment(self):
        """
        Add a comment at the root level of the XML document indicating the
        BEASTling version used to create the file, the time and date of
        generation and the original configuration file text.
        """
        comment_lines = []
        comment_lines.append("Generated by BEASTling %s on %s." % (__version__,datetime.datetime.now().strftime("%A, %d %b %Y %l:%M %p")))
        if self.config.configfile:
            comment_lines.append("Original config file:")
            comment_lines.append(self.config.configfile.write_string())
        else:
            comment_lines.append("Configuration built programmatically.")
            comment_lines.append("No config file to include.")
        self.beast.append(ET.Comment("\n".join(comment_lines)))

    def embed_data(self):
        """
        Embed a copy of each data file in a comment at the top of the XML
        document.
        """
        if not self.config.embed_data:
            return
        for filename in self.config.files_to_embed:
            self.beast.append(self.format_data_file(model.data_filename))
        for model in self.config.models:
            self.beast.append(self.format_data_file(model.data_filename))

    def format_data_file(self, filename):
        """
        Return an ElementTree node corresponding to a comment containing
        the text of the specified data file.
        """
        header = "BEASTling embedded data file: %s" % filename
        fp = open(filename, "r")
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
        attribs = {}
        attribs["id"] = "mcmc"
        attribs["spec"] = "MCMC"
        attribs["chainLength"] = str(self.config.chainlength)
        if self.config.sample_from_prior:
            attribs["sampleFromPrior"] = "true"
        self.run = ET.SubElement(self.beast, "run", attrib=attribs)
        self.add_state()
        self.add_init()
        self.add_distributions()
        self.add_operators()
        self.add_loggers()

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

    def add_tree_state(self):
        """
        Add tree-related <state> sub-elements.
        """
        tree = ET.SubElement(self.state, "tree", {"id":"Tree.t:beastlingTree", "name":"stateNode"})
        taxonset = ET.SubElement(tree, "taxonset", {"id":"taxa"})
        for lang in self.config.languages:
            ET.SubElement(taxonset, "taxon", {"id":lang,})

        param = ET.SubElement(self.state, "parameter", {"id":"birthRate.t:beastlingTree","name":"stateNode"})
        param.text="1.0"

    def add_init(self):
        """
        Add the <init> element and all its descendants.
        """
        # If a starting tree is specified, use it...
        if self.config.starting_tree:
            init = ET.SubElement(self.run, "init", {"estimate":"false", "id":"startingTree", "initial":"@Tree.t:beastlingTree", "spec":"beast.util.TreeParser","IsLabelledNewick":"true", "newick":self.config.starting_tree})
        # ...if not, use a random tree
        else:
            # But the random tree must respect any constraints!
            if self.config.monophyly:
                init = ET.SubElement(self.run, "init", {"estimate":"false", "id":"startingTree", "initial":"@Tree.t:beastlingTree", "taxonset":"@taxa", "spec":"beast.evolution.tree.ConstrainedRandomTree", "constraints":"@constraints"})
            else:
                init = ET.SubElement(self.run, "init", {"estimate":"false", "id":"startingTree", "initial":"@Tree.t:beastlingTree", "taxonset":"@taxa", "spec":"beast.evolution.tree.RandomTree"})
            popmod = ET.SubElement(init, "populationModel", {"spec":"ConstantPopulation"})
            ET.SubElement(popmod, "popSize", {"spec":"parameter.RealParameter","value":"1"})

    def add_distributions(self):
        """
        Add all probability distributions under the <run> element.
        """
        self.master_distribution = ET.SubElement(self.run,"distribution",{"id":"posterior","spec":"util.CompoundDistribution"})
        self.add_prior()
        self.add_likelihood()

    def add_prior(self):
        """
        Add all prior distribution elements.
        """
        self.prior = ET.SubElement(self.master_distribution,"distribution",{"id":"prior","spec":"util.CompoundDistribution"})
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
        glotto_iso_langs = [l for l in self.config.languages if l.lower() in self.config.classifications]
        if self.config.monophyly and glotto_iso_langs:
            attribs = {}
            attribs["id"] = "constraints"
            attribs["spec"] = "beast.math.distributions.MultiMonophyleticConstraint"
            attribs["tree"] = "@Tree.t:beastlingTree"
            attribs["newick"] = self.make_monophyly_newick(glotto_iso_langs)
            ET.SubElement(self.prior, "distribution", attribs)

    def add_calibrations(self):
        """
        Add timing calibrations to prior distribution.
        """
        if not self.config.calibrations:
            return
        for n, clade in enumerate(self.config.calibrations):
            if clade == "root":
                langs = self.config.languages
            else:
                langs = []
                for l in self.config.languages:
                    for name, glottocode in self.config.classifications.get(l.lower(),""):
                        if clade == name.lower() or clade == glottocode:
                            langs.append(l)
                            break
            if not langs:
                continue
            lower, upper = self.config.calibrations[clade]
            mean = (upper + lower) / 2.0
            stddev = (upper - mean) / 2.0
            attribs = {}
            attribs["id"] = clade + "MRCA"
            attribs["monophyletic"] = "true"
            attribs["spec"] = "beast.math.distributions.MRCAPrior"
            attribs["tree"] = "@Tree.t:beastlingTree"
            cal_prior = ET.SubElement(self.prior, "distribution", attribs)

            taxonset = ET.SubElement(cal_prior, "taxonset", {"id" : clade, "spec":"TaxonSet"})
            for lang in langs:
                ET.SubElement(taxonset, "taxon", {"idref":lang})
            normal = ET.SubElement(cal_prior, "Normal", {"id":"CalibrationNormal.%d" % n, "name":"distr", "offset":str(mean)})
            ET.SubElement(normal, "parameter", {"id":"parameter.hyperNormal-mean-%s.prior" % clade, "name":"mean", "estimate":"false"}).text = "0.0"
            ET.SubElement(normal, "parameter", {"id":"parameter.hyperNormal-sigma-%s.prior" % clade, "name":"sigma", "estimate":"false"}).text = str(stddev)

    def add_tree_prior(self):
        """
        Add Yule birth-process tree prior.
        """
        # Tree prior
        attribs = {}
        attribs["birthDiffRate"] = "@birthRate.t:beastlingTree"
        attribs["id"] = "YuleModel.t:beastlingTree"
        attribs["spec"] = "beast.evolution.speciation.YuleModel"
        attribs["tree"] = "@Tree.t:beastlingTree"
        ET.SubElement(self.prior, "distribution", attribs)

        # Birth rate
        attribs = {}
        attribs["id"] = "YuleBirthRatePrior.t:beastlingTree"
        attribs["name"] = "distribution"
        attribs["x"] = "@birthRate.t:beastlingTree"
        sub_prior = ET.SubElement(self.prior, "prior", attribs)
        uniform = ET.SubElement(sub_prior, "Uniform", {"id":"Uniform.0","name":"distr","upper":"Infinity"})

    def add_likelihood(self):
        """
        Add all likelihood distribution elements.
        """
        self.likelihood = ET.SubElement(self.master_distribution,"distribution",{"id":"likelihood","spec":"util.CompoundDistribution"})
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
            updown = ET.SubElement(self.run, "operator", {"id":"UpDown","spec":"UpDownOperator","scaleFactor":"0.5", "weight":"3.0"})
            ET.SubElement(updown, "tree", {"idref":"Tree.t:beastlingTree", "name":"up"})
            ET.SubElement(updown, "parameter", {"idref":"birthRate.t:beastlingTree", "name":"down"})
            ### Include clock rates in up/down only if calibrations are given
            if self.config.calibrations:
                for clock in self.config.clocks:
                    ET.SubElement(updown, "parameter", {"idref":clock.mean_rate_id, "name":"down"})

        # Birth rate scaler
        # Birth rate is *always* scaled.
        ET.SubElement(self.run, "operator", {"id":"YuleBirthRateScaler.t:beastlingTree","spec":"ScaleOperator","parameter":"@birthRate.t:beastlingTree", "scaleFactor":"0.5", "weight":"3.0"})

    def add_loggers(self):
        """
        Add all <logger> elements.
        """
        self.add_screen_logger()
        self.add_tracer_logger()
        self.add_tree_logger()

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
        if not(self.config.log_probabilities or self.config.log_params or self.config.log_all):
            return
        tracer_logger = ET.SubElement(self.run,"logger",{"id":"tracelog","fileName":self.config.basename+".log","logEvery":str(self.config.log_every),"model":"@posterior","sanitiseHeaders":"true","sort":"smart"})
        # Log prior, likelihood and posterior
        if self.config.log_probabilities or self.config.log_all:
            ET.SubElement(tracer_logger,"log",{"idref":"prior"})
            ET.SubElement(tracer_logger,"log",{"idref":"likelihood"})
            ET.SubElement(tracer_logger,"log",{"idref":"posterior"})
        # Log Yule birth rate
        if self.config.log_params or self.config.log_all:
            ET.SubElement(tracer_logger,"log",{"idref":"birthRate.t:beastlingTree"})
            for clock in self.config.clocks:
                clock.add_param_logs(tracer_logger)
            for model in self.config.all_models:
                    model.add_param_logs(tracer_logger)

        # Log tree height
        if not self.config.tree_logging_pointless:
            ET.SubElement(tracer_logger,"log",{
                "id":"treeHeight",
                "spec":"beast.evolution.tree.TreeHeightLogger",
                "tree":"@Tree.t:beastlingTree"})

        # Log calibration clade heights
        for clade in self.config.calibrations:
            ET.SubElement(tracer_logger,"log",{"idref":"%sMRCA" % clade})

    def add_tree_logger(self):
        """
        Add tree logger, if configured to do so.
        """
        if ((self.config.log_trees or self.config.log_all) and not
            self.config.tree_logging_pointless):
            tree_logger = ET.SubElement(self.run, "logger", {"mode":"tree", "fileName":self.config.basename+".nex","logEvery":str(self.config.log_every),"id":"treeWithMetaDataLogger"})
            log = ET.SubElement(tree_logger, "log", attrib={"id":"TreeLogger","spec":"beast.evolution.tree.TreeWithMetaDataLogger","tree":"@Tree.t:beastlingTree"})
            if self.config.geo_config.get("log_locations",False):
                ET.SubElement(log, "metadata", {
                    "id":"location",
                    "spec":"sphericalGeo.TraitFunction",
                    "likelihood":"@sphericalGeographyLikelihood"}).text = "0.0"

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

    def make_monophyly_newick(self, langs):
        """
        Transforms a list of list of language identifiers (ISO or Glottocodes)
        into a Newick tree representing the corresponding Glottolog
        family structure, suitable as use for a BEAST monophyly constraint.
        """
        # First we build a "monophyly structure".  This can be done in either
        # a "top-down" or "bottom-up" way.
        if self.config.monophyly_end_depth is not None:
            # A power user has explicitly provided start and end depths
            start = self.config.monophyly_start_depth
            end = self.config.monophyly_end_depth
        elif self.config.monophyly_direction == "top_down":
            # Compute start and end in a top-down fashion
            start = self.config.monophyly_start_depth
            end = start + self.config.monophyly_levels
        elif self.config.monophyly_direction == "bottom_up":
            # Compute start and end in a bottom-up fashion
            classifications = [self.config.classifications[name.lower()] for name in langs]
            end = max([len(c) for c in classifications]) - self.config.monophyly_start_depth
            start = max(0, end - self.config.monophyly_levels)
        struct = self.make_monophyly_structure(langs, depth=start, maxdepth=end)
        # Now we serialise the "monophyly structure" into a Newick tree.
        return self.make_monophyly_string(struct)

    def make_monophyly_structure(self, langs, depth, maxdepth):
        """
        Recursively partition a list of languages (ISO or Glottocodes) into
        lists corresponding to their Glottolog classification.  The process
        may be halted part-way down the Glottolog tree.
        """
        if depth > maxdepth:
            # We're done, so terminate recursion
            return langs

        def subgroup(name, depth):
            ancestors = self.config.classifications[name.lower()]
            return ancestors[depth][0] if depth < len(ancestors) else ''

        def sortkey(i):
            """
            Callable to pass into `sorted` to port sorting behaviour from py2 to py3.

            :param i: Either a string or a list (of lists, ...) of strings.
            :return: Pair (nesting level, first string)
            """
            d = 0
            while isinstance(i, list):
                d -= 1
                i = i[0] if i else ''
            return d, i

        # Find the ancestor of all the given languages at at particular depth 
        # (i.e. look `depth` nodes below the root of the Glottolog tree)
        levels = list(set([subgroup(l, depth) for l in langs]))
        if len(levels) == 1:
            # If all languages belong to the same classificatio at this depth,
            # there are two possibilities
            if levels[0] == "":
                # If the common classification is an empty string, then we know
                # that there is no further refinement possible, so stop
                # the recursion here.
                langs.sort()
                return langs
            else:
                # If the common classification is non-empty, we need to
                # descend further, since some languages might get
                # separated later
                return self.make_monophyly_structure(langs, depth+1, maxdepth)
        else:
            # If the languages belong to multiple classifications, split them
            # up accordingly and then break down each classification
            # individually.

            partition = [[l for l in langs if subgroup(l, depth) == level] for level in levels]
            partition = [part for part in partition if part]
            return sorted(
                [self.make_monophyly_structure(group, depth+1, maxdepth)
                 for group in partition],
                key=sortkey)

    def make_monophyly_string(self, struct, depth=0):
        """
        Converts a structure of nested lists into Newick string.
        """
        if not type([]) in [type(x) for x in struct]:
            return "(%s)" % ",".join(struct)
        else:
            return "(%s)" % ",".join([self.make_monophyly_string(substruct) for substruct in struct])
