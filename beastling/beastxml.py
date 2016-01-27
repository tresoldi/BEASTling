import codecs
import datetime
import os
import pdb
import sys
import xml.etree.ElementTree as ET

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

class BeastXml:

    def __init__(self, config):
        self.config = config
        if not self.config.processed:
            self.config.process()
        self.build_xml()

    def build_xml(self):

        # Root "beast" node
        attribs = {}
        attribs["beautitemplate"] = "Standard"
        attribs["beautistatus"] = ""
        attribs["namespace"] = "beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood"
        attribs["version"] ="2.0"
        self.beast = ET.Element("beast", attrib=attribs)

        # "Generated by..." comment
        comment_lines = []
        comment_lines.append("Generated by BEASTling %s on %s." % (__version__,datetime.datetime.now().strftime("%A, %d %b %Y %l:%M %p")))
        if self.config.configfile_text:
            comment_lines.append("Original config file:")
            comment_lines.append(self.config.configfile_text)
        else:
            comment_lines.append("Configuration built programmatically.")
            comment_lines.append("No config file to include.")
        self.beast.append(ET.Comment("\n".join(comment_lines)))

        # Maps
        for a, b in beast_maps.maps:
            mapp = ET.SubElement(self.beast, "map", attrib={"name":a})
            mapp.text = b

        # Taxon set
        taxonset = ET.SubElement(self.beast, "taxonset", {"id":"taxa"})
        for lang in self.config.languages:
            ET.SubElement(taxonset, "taxon", {"id":lang,})

        tree = ET.SubElement(self.beast, "tree", {"id":"Tree.t:beastlingTree"})
        ET.SubElement(tree, "taxonset", {"idref":"taxa"})

        # Monophyly constraints
        if self.config.monophyly:
            tree = ET.SubElement(self.beast, "tree", {"initial":"@Tree.t:beastlingTree", "taxonset":"@taxa", "spec":"beast.evolution.tree.ConstrainedRandomTree", "id":"coalescentSimulator","constraints":"@constraints"})
        else:
            tree = ET.SubElement(self.beast, "tree", {"initial":"@Tree.t:beastlingTree", "taxonset":"@taxa", "spec":"beast.evolution.tree.RandomTree", "id":"coalescentSimulator"})
        popmod = ET.SubElement(tree, "populationModel", {"spec":"ConstantPopulation"})
        ET.SubElement(popmod, "popSize", {"spec":"parameter.RealParameter","value":"1"})

        # Run
        attribs = {}
        attribs["chainLength"] = str(self.config.chainlength)
        attribs["id"] = "mcmc"
        attribs["spec"] = "MCMC"
        self.run = ET.SubElement(self.beast, "run", attrib=attribs)

        # Init
        init = ET.SubElement(self.run, "init", {"idref":"coalescentSimulator"})

        # State
        state = ET.SubElement(self.run, "state", {"id":"state","storeEvery":"5000"})
        ET.SubElement(state, "stateNode", {"idref":"Tree.t:beastlingTree"})
        param = ET.SubElement(state, "parameter", {"id":"birthRate.t:beastlingTree","name":"stateNode"})
        param.text="1.0"

        # Common clocks
        self.common_clocks = []
        for model in self.config.models:
            if not model.rate_variation:
                if model.clock not in self.common_clocks:
                    self.common_clocks.append(model.clock)
        for clock in self.common_clocks:
            attribs = {}
            attribs["id"] = "clockRate_%s.c" % clock
            attribs["name"] = "stateNode"
            parameter = ET.SubElement(state, "parameter", attribs)
            parameter.text="1.0"

        for model in self.config.models:
            model.add_state(state)
           
        for model in self.config.models:
            model.add_misc(self.beast)

        # Distributions
        self.master_distribution = ET.SubElement(self.run,"distribution",{"id":"posterior","spec":"util.CompoundDistribution"})

        ## Prior
        self.add_prior()
        for model in self.config.models:
            model.add_prior(self.prior)

        ## Likelihood
        self.likelihood = ET.SubElement(self.master_distribution,"distribution",{"id":"likelihood","spec":"util.CompoundDistribution"})
        for model in self.config.models:
            model.add_likelihood(self.likelihood)

        # Operators
        self.add_operators()

        # Logging
        self.add_loggers()

    def add_prior(self):

        """Add "master prior" features, independent of any data
        or models.  E.g. monophyly constraints, clade calibration
        dates, tree priors, etc."""

        self.prior = ET.SubElement(self.master_distribution,"distribution",{"id":"prior","spec":"util.CompoundDistribution"})

        # Monophyly
        glotto_iso_langs = [l for l in self.config.languages if l.lower() in self.config.classifications]
        if self.config.monophyly and glotto_iso_langs:
            attribs = {}
            attribs["id"] = "constraints"
            attribs["spec"] = "beast.math.distributions.MultiMonophyleticConstraint"
            attribs["tree"] = "@Tree.t:beastlingTree"
            attribs["newick"] = self.make_monophyly_newick(glotto_iso_langs)
            ET.SubElement(self.prior, "distribution", attribs)

        # Calibration
        if self.config.calibrations:
            for n, clade in enumerate(self.config.calibrations):
                if clade == "root":
                    langs = self.config.languages
                else:
                    langs = [l for l in self.config.languages if any([c==clade for c in [x.lower() for x in self.config.classifications[l.lower()].split(",")]])]
                if not langs:
                    continue
                lower, upper = self.config.calibrations[clade]
                mean = (upper + lower) / 2.0
                stddev = (upper - mean) / 2.0
                attribs = {}
                attribs["id"] = clade + "-calibration.prior"
                attribs["monophyletic"] = "true"
                attribs["spec"] = "beast.math.distributions.MRCAPrior"
                attribs["tree"] = "@Tree.t:beastlingTree"
                cal_prior = ET.SubElement(self.prior, "distribution", attribs)

                taxonset = ET.SubElement(cal_prior, "taxonset", {"id" : clade, "spec":"TaxonSet"})
                for lang in langs:
                    ET.SubElement(taxonset, "taxon", {"idref":lang, "spec":"Taxon"})
                normal = ET.SubElement(cal_prior, "Normal", {"id":"CalibrationNormal.%d" % n, "name":"distr", "offset":str(mean)})
                ET.SubElement(normal, "parameter", {"id":"parameter.hyperNormal-mean-%s.prior" % clade, "name":"mean", "estimate":"false"}).text = "0.0"
                ET.SubElement(normal, "parameter", {"id":"parameter.hyperNormal-sigma-%s.prior" % clade, "name":"sigma", "estimate":"false"}).text = str(stddev)

        # Common clocks
        for clock in self.common_clocks:
            sub_prior = ET.SubElement(self.prior, "prior", {"id":"clockPrior_%s.s" % clock, "name":"distribution","x":"@clockRate_%s.c" % clock})
            uniform = ET.SubElement(sub_prior, "Uniform", {"id":"UniformClockPrior_%s" % clock, "name":"distr", "upper":"Infinity"})

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


    def add_operators(self):

        # Tree topology operators
        ET.SubElement(self.run, "operator", {"id":"UniformOperator.t:beastlingTree","spec":"Uniform","tree":"@Tree.t:beastlingTree","weight":"600.0"})
        ET.SubElement(self.run, "operator", {"id":"SubtreeSlide.t:beastlingTree","spec":"SubtreeSlide","tree":"@Tree.t:beastlingTree","markclades":"true", "weight":"600.0"})
        ET.SubElement(self.run, "operator", {"id":"narrow.t:beastlingTree","spec":"Exchange","tree":"@Tree.t:beastlingTree","markclades":"true", "weight":"600.0"})
        ET.SubElement(self.run, "operator", {"id":"wide.t:beastlingTree","isNarrow":"false","spec":"Exchange","tree":"@Tree.t:beastlingTree","markclades":"true", "weight":"600.0"})
        ET.SubElement(self.run, "operator", {"id":"WilsonBalding.t:beastlingTree","spec":"WilsonBalding","tree":"@Tree.t:beastlingTree","markclades":"true","weight":"600.0"})

        # Scalers for numeric parameters
        ET.SubElement(self.run, "operator", {"id":"treeScaler.t:beastlingTree","scaleFactor":"1.0","spec":"ScaleOperator","tree":"@Tree.t:beastlingTree","weight":"10.0"})

        # Birth rate scaler
        ET.SubElement(self.run, "operator", {"id":"YuleBirthRateScaler.t:beastlingTree","spec":"ScaleOperator","parameter":"@birthRate.t:beastlingTree", "scaleFactor":"1.0", "weight":"3.0"})

        # Clock scalers (only for calibrated analyses)
        if self.config.calibrations:
            for clock in self.common_clocks:
                ET.SubElement(self.run, "operator", {"id":"geoMuScaler.c:clockRate_%s" % clock, "spec":"ScaleOperator","parameter":"@clockRate_%s.c" % clock, "scaleFactor":"1.0","weight":"10.0"})

        # Up/down
        updown = ET.SubElement(self.run, "operator", {"id":"UpDownCommon","spec":"UpDownOperator","scaleFactor":"1.0", "weight":"30.0"})
        ET.SubElement(updown, "tree", {"idref":"Tree.t:beastlingTree", "name":"up"})
        ET.SubElement(updown, "parameter", {"idref":"birthRate.t:beastlingTree", "name":"down"})
        if self.config.calibrations:
            # Estimate clocks if calibrations given
            for clock in self.common_clocks:
                ET.SubElement(updown, "parameter", {"idref":"clockRate_%s.c" % clock, "name":"down"})

        # Model specific operators
        for model in self.config.models:
            model.add_operators(self.run)

    def add_loggers(self):

        # Screen logger
        if self.config.screenlog:
            screen_logger = ET.SubElement(self.run, "logger", attrib={"id":"screenlog", "logEvery":"10000"})
            log = ET.SubElement(screen_logger, "log", attrib={"arg":"@posterior", "id":"ESS.0", "spec":"util.ESS"})
            log = ET.SubElement(screen_logger, "log", attrib={"idref":"prior"})
            log = ET.SubElement(screen_logger, "log", attrib={"idref":"likelihood"})
            log = ET.SubElement(screen_logger, "log", attrib={"idref":"posterior"})

        # Tracer log
        if self.config.log_probabilities or self.config.log_params:
            tracer_logger = ET.SubElement(self.run,"logger",{"id":"tracelog","fileName":self.config.basename+".log","logEvery":"10000","model":"@posterior","sanitiseHeaders":"true","sort":"smart"})
            if self.config.log_probabilities:
                ET.SubElement(tracer_logger,"log",{"idref":"prior"})
                ET.SubElement(tracer_logger,"log",{"idref":"likelihood"})
                ET.SubElement(tracer_logger,"log",{"idref":"posterior"})
            if self.config.log_params:
                ET.SubElement(tracer_logger,"log",{"idref":"birthRate.t:beastlingTree"})
                for clock in self.common_clocks:
                    ET.SubElement(tracer_logger,"log",{"idref":"clockRate_%s.c" % clock})
                for model in self.config.models:
                    model.add_param_logs(tracer_logger)
                
        # Tree log
        if self.config.log_trees:
            tree_logger = ET.SubElement(self.run, "logger", {"mode":"tree", "fileName":self.config.basename+".nex","logEvery":"10000","id":"treeWithMetaDataLogger"})
            log = ET.SubElement(tree_logger, "log", attrib={"id":"TreeLogger","spec":"beast.evolution.tree.TreeWithMetaDataLogger","tree":"@Tree.t:beastlingTree"})

    def write_file(self, filename=None):
        indent(self.beast)
        xml_string = ET.tostring(self.beast, encoding="UTF-8")
        if not filename:
            filename = self.config.basename+".xml"
        if filename in ("stdout", "-"):
            sys.stdout.write(unicode(xml_string, "utf-8"))
        else:
            fp = codecs.open(filename, "w", "UTF-8")
            fp.write(unicode(xml_string, "utf-8"))
            fp.close()

    def make_tight_monophyly_structure(self, langs, depth=0, maxdepth=sys.maxint):
        if depth > maxdepth:
            return langs
        levels = list(set([self.config.classifications[l.lower()].split(",")[depth] for l in langs]))
        if len(levels) == 1:
            if levels[0] == "":
                langs.sort()
                return langs
            else:
                return self.make_tight_monophyly_structure(langs, depth+1, maxdepth)
        else:
            partition = [[l for l in langs if self.config.classifications[l.lower()].split(",")[depth] == level] for level in levels]
            partition = [part for part in partition if part]
            return sorted([self.make_tight_monophyly_structure(group, depth+1, maxdepth) for group in partition])

    def make_loose_monophyly_structure(self, langs):
        points = self.config.monophyly
        return [[l for l in langs if point in self.config.classifications[l.lower()] ] for point in points]

    def make_monophyly_string(self, struct, depth=0):
        if not type([]) in [type(x) for x in struct]:
            return "(%s)" % ",".join(struct)
        else:
            return "(%s)" % ",".join([self.make_monophyly_string(substruct) for substruct in struct])

    def make_monophyly_newick(self, langs):
        if self.config.monophyly_grip == "tight":
            struct = self.make_tight_monophyly_structure(langs, depth=self.config.monophyly_start_depth, maxdepth=self.config.monophyly_end_depth)
        elif self.config.monophyly_grip == "loose":
            struct = self.make_loose_monophyly_structure(langs)
        return self.make_monophyly_string(struct)

