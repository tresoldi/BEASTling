import codecs
import os
import itertools
import math
import xml.etree.ElementTree as ET

import scipy.stats

from .basemodel import BaseModel

class IncompatibleSettingOptions(Warning): pass

# We make some use of a product function.
from operator import mul
try:
    from functools import reduce
except ImportError:
    pass
def product(iterable):
    return reduce(mul, iterable, 1)
    

class CorrelatedModel(BaseModel):

    def __init__(self, model_config, global_config):
        if model_config.get("rate_variation", False):
            raise IncompatibleSettingOptions("Please do not specify rate_variation.")
        BaseModel.__init__(self, model_config, global_config)

    def preprocess(self):
        # Remove features which are in the config but not the
        # data file
        if self.config.languages == ["*"]:
            self.languages = self.data.keys()
        else:
            self.languages = set(self.config.languages) & set(self.data.keys())
        traits = [t for t in self.features if
                any([t in self.data[lang] for lang in self.languages]
                    )]

        # Remove features which are in the config but are entirely
        # question marks for the specified languages.
        traits = [t for t in traits if
                not all([self.data[lang][t] == "?" for lang in self.languages]
                    )]

        # To start with, this trait has one value ("") and zero rates
        # in between.
        self.valuecounts = {}
        self.dimensions = {}
        self.counts = {}
        self.codemaps = {}
        bad_traits = []
        for trait in self.features:
            all_values = [self.data[l][trait] for l in self.languages]
            all_values = [v for v in all_values if v != "?" if v != "-"]
            uniq = list(set(all_values))
            counts = {}
            for v in all_values:
                counts[v] = all_values.count(v)
            uniq = list(set(all_values))
            # Sort uniq carefully.
            # Possibly all feature values are numeric strings, e.g. "1", "2", "3".
            # If we sort these as strings then we get weird things like "10" < "2".
            # This can actually matter for things like ordinal models.
            # So sort by their int values...
            if all([v.isdigit() for v in uniq]):
                uniq.sort(key=int)
            # ...otherwise, just sort normally
            else:
                uniq.sort()
            if len(uniq) == 0 or (len(uniq) == 1 and self.remove_constant_features):
                bad_traits.append(trait)
                continue
            N = len(uniq)

            self.valuecounts[trait] = N
            self.dimensions[trait] = N*(N-1)/2
            self.codemaps[trait] = self.build_codemap(uniq)
            self.counts[trait] = counts
        self.features = [t for t in self.features if t not in bad_traits]
        self.features.sort()
        self.max_n_rates = (product(self.valuecounts.values()) * sum([i-1 for i in self.valuecounts.values()]))

    def add_state(self, state):
        self.add_clock_state(state)

        trait = "compound"
        traitname = "%s:%s" % (self.name, trait)
        
        parameter = ET.SubElement(
            state, "parameter",
            {"dimension": str(self.max_n_rates),
             "id": "rawRates.s:%s" % traitname,
             "lower": "0.0",
             "name": "stateNode"})
        parameter.text="1.0"

        groupings = ET.SubElement(
            state, "parameter",
            {"spec": "beast.core.parameter.IntegerParameter",
             "dimension" : str(self.max_n_rates),
             "id" : "rateGroupings.s:%s" % traitname,
             "name" : "stateNode"})
        groupings.text=" ".join(map(str, range(self.max_n_rates)))

        sizes = ET.SubElement(
            state, "parameter",
            {"spec": "beast.core.parameter.IntegerParameter",
             "dimension" : str(self.max_n_rates),
             "id" : "rateGroupingSizes.s:%s" % traitname,
             "name" : "stateNode"})
        sizes.text="1"

    def add_prior(self, prior):
        self.add_clock_prior(prior)

        trait = "compound"
        n = 0
        traitname = "%s:%s" % (self.name, trait)

        # Relative rate
        sub_prior = ET.SubElement(prior, "prior", {"id":"ratesPrior.s:%s" % traitname, "name":"distribution","x":"@rawRates.s:%s"% traitname})
        dirichlet  = ET.SubElement(
            sub_prior, "distr",
            {"id":"Dirichlet:%s.%d.0" % (traitname, n),
             "sizes" : "@rateGroupingSizes.s:%s" % traitname,
             "spec": "parameterclone.helpers.RescaledDirichlet"})

    def add_data(self, distribution, trait, traitname):
        data = ET.SubElement(distribution,"data",{"id":traitname, "spec":"Alignment"})
        # How many characters do we need to reserve?
        for lang in self.languages:
            valuestring = ";;".join("%s" % 
                self.data[lang][trait]
                for trait in self.features)
            seq = ET.SubElement(data, "sequence", {"id":"seq_%s_%s" % (lang, traitname), "taxon":lang, "value":valuestring})
        userdatatype = ET.SubElement(
            data, "userDataType",
            {"id": "traitDataType.%s"%traitname,
             "split": ";;",
             # The componentSizes are bigger than the valuecounts because of the way ambiguities ("?" and "-", so 2) are encoded.
             # NOTE: Once we start to support ambiguities, this line needs adapting! 
             "componentSizesIncludingAmbiguities": " ".join(str(self.valuecounts[trait]+2) for trait in self.features),
             "spec":
              "correlatedcharacters.polycharacter.CompoundDataType"})
        for trait in self.features:
            subdatatype = ET.SubElement(
                userdatatype, "components",
                {"id":"traitDataType.%s:%s"%(traitname,trait),
                 "spec":"beast.evolution.datatype.UserDataType",
                 "codeMap":self.codemaps[trait],
                 "states":str(self.valuecounts[trait]),
                 "characterName": trait,
                })

    def add_likelihood(self, likelihood):
        trait = "compound"
        traitname = "%s:%s" % (self.name, trait)
        distribution = ET.SubElement(likelihood, "distribution",{"id":"traitedtreeLikelihood.%s" % traitname,"spec":"TreeLikelihood","useAmbiguities":"true"})

        tree = ET.SubElement(distribution, "tree", {"idref":"Tree.t:beastlingTree"})

        # Sitemodel
        self.add_sitemodel(distribution, trait, traitname)

        # Branchrate
        branchrate = ET.SubElement(distribution, "branchRateModel", {"id":"StrictClockModel.c:%s"%traitname,"spec":"beast.evolution.branchratemodel.StrictClockModel","clock.rate":"@clockRate.c:%s" % self.name})

        # Data
        self.add_data(distribution, trait, traitname)
        
    def add_sitemodel(self, distribution, trait, traitname):
        if trait != "compound":
            raise ValueError("There should only be a 'compound' trait")

        # Sitemodel
        sitemodel = ET.SubElement(distribution, "siteModel", {"id":"SiteModel.%s"%traitname,"spec":"SiteModel", "shape":"1","proportionInvariant":"0"})

        substmodel = ET.SubElement(
            sitemodel, "substModel",
            {"id":"substitutionmodel.s:%s"%traitname,
             "spec":"correlatedcharacters.polycharacter.CorrelatedSubstitutionModel"})
        actualrates = ET.SubElement(
            substmodel, "rates",
            {"id": "actualRates.s:%s"%traitname,
             "spec": "parameterclone.selector.Selector",
             "entry": " ".join(str(i) for i in range(self.max_n_rates)),
             "parameters": "@rawRates.s:%s" % traitname,
             "groupings": "@rateGroupings.s:%s" % traitname,
            })
        shape = ET.SubElement(substmodel, "datatype",
                              {"idref":"traitDataType.%s"%traitname})
        freq = ET.SubElement(substmodel,"frequencies",{"id":"traitfreqs.s:%s"%traitname,"spec":"Frequencies"})
        states = product(self.valuecounts.values())
        freq_string = str(1./states)
        ET.SubElement(freq,"parameter",{
            "dimension":str(states),
            "id":"traitfrequencies.s:%s"%traitname,
            "name":"frequencies"}).text=freq_string

    def add_operators(self, run):

        BaseModel.add_operators(self, run)

        trait = "compound"
        traitname = "%s:%s" % (self.name, trait)

        ET.SubElement(
            run, "operator",
            {"id": "merger:%s" % traitname,
             "spec": "parameterclone.splitandmerge.MergeOperator",
             "parameters": "@rawRates.s:%s" % traitname,
             "groupings": "@rateGroupings.s:%s" % traitname,
             "sizes": "@rateGroupingSizes.s:%s" % traitname,
             "weight": "10.0"
             })

        ET.SubElement(
            run, "operator",
            {"id": "splitter:%s" % traitname,
             "spec": "parameterclone.splitandmerge.SplitOperator",
             "parameters": "@rawRates.s:%s" % traitname,
             "groupings": "@rateGroupings.s:%s" % traitname,
             "sizes": "@rateGroupingSizes.s:%s" % traitname,
             "weight": "10.0"
             })

        ET.SubElement(
            run, "operator",
            {"id": "exchanger:%s" % traitname,
             "spec": "beast.evolution.operators.DeltaExchangeOperator",
             "parameter": "@rawRates.s:%s" % traitname,
             "delta": "1.0",
             "weightvector": "@rateGroupingSizes.s:%s" % traitname,
             "weight": "10.0"
             })

    def add_param_logs(self, logger):
        BaseModel.add_param_logs(self, logger)
        trait = "compound"
        traitname = "%s:%s" % (self.name, trait)
        ET.SubElement(logger,"log",{"idref":"actualRates.s:%s" % traitname})
        ET.SubElement(
            logger, "log",
            {"id": "independencyLogger.s:%s" % traitname,
             "spec": "correlatedcharacters.polycharacter.IndependencyLogger",
             "model": "@substitutionmodel.s:%s" % traitname})
