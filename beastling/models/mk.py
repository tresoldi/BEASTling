import xml.etree.ElementTree as ET

from .basemodel import BaseModel


class MKModel(BaseModel):

    package_notice = """[DEPENDENCY]: The Lewis Mk substitution model is implemented in the BEAST package "morph-models"."""

    def __init__(self, model_config, global_config):

        BaseModel.__init__(self, model_config, global_config)

    def add_sitemodel(self, distribution, feature, fname):

            # Sitemodel
            mr = self.get_mutation_rate(feature, fname)
            sitemodel = ET.SubElement(distribution, "siteModel", {"id":"SiteModel.%s"%fname,"spec":"SiteModel", "mutationRate":mr,"shape":"1","proportionInvariant":"0"})

            substmodel = ET.SubElement(sitemodel, "substModel",{"id":"mk.s:%s"%fname,"spec":"LewisMK","datatype":"@featureDataType.%s" % fname})
            # Do empirical frequencies
            # We don't need to do anything for uniform freqs
            # as the implementation of LewisMK handles it
            if self.frequencies == "empirical":
                freq = ET.SubElement(substmodel,"frequencies",{"id":"feature_freqs.s:%s"%fname,"spec":"Frequencies", "data":"@data_%s"%fname})
            elif self.frequencies == "estimate":
                freq = ET.SubElement(substmodel,"frequencies",{"id":"feature_freqs.s:%s"%fname,"spec":"Frequencies", "frequencies":"@feature_freqs_param.s:%s"%fname})
