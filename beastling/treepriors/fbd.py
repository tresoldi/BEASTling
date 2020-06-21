import typing as t

from beastling.treepriors.base import TreePrior
from beastling.util import xml


class FossilizedBirthDeathTree (TreePrior):
    """A generalization of the Fossilized Birth Death (FBD) tree prior.

    The prior is implemented in the [SA
    package](https://github.com/CompEvol/sampled-ancestors) for Beast2 and
    based on

    > Alexandra Gavryushkina, David Welch, Tanja Stadler, Alexei J. Drummond
    > (2014) Bayesian Inference of Sampled Ancestor Trees for Epidemiology and
    > Fossil Calibration. PLoS Computational Biology
    > https://doi.org/10.1371/journal.pcbi.1003919

    """
    __type__ = "fbd"
    package_notice = ("fbd", "Sampled Ancestors")

    def add_parameters(self, state):
        """ Add tree-related <state> sub-elements.

        Add real parameters for the following FBD parameters:
         - diversification rate
         - origin
         - sampling proportion
         - turnover
         - removal probability
         - rho
        """
        param = xml.parameter(
            state,
            id="diversificationRateFBD.t:beastlingTree", lower="0.0", name="stateNode",
            text="0.001")
        if self.birthrate_estimate is not None:
            # We have already estimated a better starting point than the arbitrary 1.0
            param.text = str(self.birthrate_estimate)
        xml.parameter(
            state,
            id="originFBD.t:beastlingTree", lower="0.0", name="stateNode",
            text="10000.0")
        if self.treeheight_estimate is not None:
            # We have already estimated a better starting point than the arbitrary 10000.0
            param.text = str(self.treehight_estimate)
        xml.parameter(
            state,
            id="samplingProportionFBD.t:beastlingTree", lower="0.0", name="stateNode", upper="1.0",
            text="0.5")
        xml.parameter(
            state,
            id="turnoverFBD.t:beastlingTree", lower="0.0", name="stateNode", upper="1.0",
            text="0.5")
        xml.parameter(
            state,
            id="rFBD.t:beastlingTree", lower="0.0", name="stateNode", upper="1.0",
            text="0.0")
        xml.parameter(
            state,
            # TODO: Why estimate=False? This is adapted from the bears example from SA.
            # The parameter describes "Probability of an individual to be sampled at present"
            id="rhoFBD.t:beastlingTree", lower="0.0", name="stateNode", upper="1.0",
            estimate="false", text="1.0")

    def add_prior(self, beastxml):
        """Add a calibrated fossilized birth-death tree prior.

        In addition to the tree prior, add uniform priors on the parameters,
        some of them are improper priors.

        """
        xml.distribution(
            beastxml.prior,
            id="FBD.t:beastlingTree", spec="beast.evolution.speciation.SABirthDeathModel",
            conditionOnRhoSampling="true",
            diversificationRate="@diversificationRateFBD.t:beastlingTree",
            origin="@originFBD.t:beastlingTree",
            samplingProportion="@samplingProportionFBD.t:beastlingTree",
            tree="@Tree.t:beastlingTree",
            turnover="@turnoverFBD.t:beastlingTree",
            removalProbability="@rFBD.t:beastlingTree",
            rho="@rhoFBD.t:beastlingTree")
        sub_prior = xml.prior(
            beastxml.prior,
            id="diversificationRatePriorFBD.t:beastlingTree", name="distribution",
            x="@diversificationRateFBD.t:beastlingTree")
        xml.Uniform(sub_prior, name="distr", upper="Infinity")
        sub_prior = xml.prior(
            beastxml.prior,
            id="originPriorFBD.t:beastlingTree", name="distribution",
            x="@originFBD.t:beastlingTree")
        xml.Uniform(sub_prior, name="distr", upper="Infinity")
        sub_prior = xml.prior(
            beastxml.prior,
            id="samplingProportionPriorFBD.t:beastlingTree", name="distribution",
            x="@samplingProportionFBD.t:beastlingTree")
        xml.Beta(sub_prior, name="distr", alpha="20.0", beta="0.5")
        sub_prior = xml.prior(
            beastxml.prior,
            id="turnoverPriorFBD.t:beastlingTree", name="distribution",
            x="@turnoverFBD.t:beastlingTree")
        xml.Uniform(sub_prior, name="distr")

    def add_fine_logging(self, tracer_logger):
        """Log the parameters of the tree prior.

        """
        xml.log(tracer_logger, idref="diversificationRateFBD.t:beastlingTree")
        xml.log(tracer_logger, idref="originFBD.t:beastlingTree")
        xml.log(tracer_logger, idref="samplingProportionFBD.t:beastlingTree")
        xml.log(tracer_logger, idref="turnoverFBD.t:beastlingTree")
        xml.log(tracer_logger, idref="rFBD.t:beastlingTree")
        xml.log(tracer_logger, idref="rhoFBD.t:beastlingTree")

    def add_operators(self, beastxml):
        super().add_operators(beastxml)
        if beastxml.config.languages.sample_branch_lengths:
            updown = xml.operator(
                beastxml.run,
                attrib={
                    "id": "UpDown",
                    "spec": "UpDownOperator",
                    "scaleFactor": "0.5",
                    "weight": "3.0"})

            xml.tree(updown, idref="Tree.t:beastlingTree", name="up")
            xml.parameter(updown, idref="originFBD.t:beastlingTree", name="up")
            xml.parameter(updown, idref="diversificationRateFBD.t:beastlingTree", name="down")
            # Include clock rates in up/down only if calibrations are given
            if beastxml.config.calibrations:
                for clock in beastxml.config.clocks:
                    if clock.estimate_rate:
                        xml.parameter(updown, idref=clock.mean_rate_id, name="down")

        # Birth rate scaler
        # Birth rate is *always* scaled.
        xml.operator(
            beastxml.run,
            attrib={
                "id": "YuleBirthRateScaler.t:beastlingTree",
                "spec": "ScaleOperator",
                "parameter": "@diversificationRateFBD.t:beastlingTree",
                "scaleFactor": "0.5",
                "weight": "3.0"})
        xml.operator(
            beastxml.run,
            attrib={
                "id": "SamplingScaler.t:beastlingTree",
                "spec": "ScaleOperator",
                "parameter": "@samplingProportionFBD.t:beastlingTree",
                "scaleFactor": "0.8",
                "weight": "1.0"})
        xml.operator(
            beastxml.run,
            attrib={
                "id": "DeathRateScaler.t:beastlingTree",
                "spec": "ScaleOperator",
                "parameter": "@turnoverFBD.t:beastlingTree",
                "scaleFactor": "0.5",
                "weight": "3.0"})
        xml.operator(
            beastxml.run,
            attrib={
                "id": "RemovalScaler.t:beastlingTree",
                "spec": "ScaleOperator",
                "parameter": "@rFBD.t:beastlingTree",
                "scaleFactor": "0.5",
                "weight": "3.0"})

        # Continue to assume complete sampling for contemporary leaves
        #         "parameter": "@rhoFBD.t:beastlingTree",
