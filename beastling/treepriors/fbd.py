import typing as t

from beastling.treepriors.base import TreePrior
from beastling.util import xml

# The definition in the BEAST package, for reference
BEAST_DEF = """
    //'direct' parameters
    public Input<RealParameter> originInput =
            new Input<RealParameter>("origin", "The time when the process started", (RealParameter)null);
    public Input<RealParameter> birthRateInput =
            new Input<RealParameter>("birthRate", "Birth rate", Input.Validate.REQUIRED);
    public Input<Function> deathRateInput =
            new Input<Function>("deathRate", "Death rate", Input.Validate.REQUIRED);
    public Input<RealParameter> samplingRateInput =
            new Input<RealParameter>("samplingRate", "Sampling rate per individual", Input.Validate.REQUIRED);

    //transformed parameters:
    public Input<RealParameter> expectedNInput =
            new Input<RealParameter>("expectedN", "The expected-N-at-present parameterisation of T",(RealParameter)null);
    public Input<RealParameter> diversificationRateInput =
            new Input<RealParameter>("diversificationRate", "Net diversification rate. Birth rate - death rate", Input.Validate.XOR, birthRateInput);
    public Input<Function> turnoverInput =
            new Input<Function>("turnover", "Turnover. Death rate/birth rate", Input.Validate.XOR, deathRateInput);
    public Input<RealParameter> samplingProportionInput =
            new Input<RealParameter>("samplingProportion", "The probability of sampling prior to death. Sampling rate/(sampling rate + death rate)", Input.Validate.XOR, samplingRateInput);


    // r parameter
    public Input<RealParameter> removalProbability =
            new Input<RealParameter>("removalProbability", "The probability that an individual is removed from the process after the sampling", Input.Validate.REQUIRED);

    public Input<RealParameter> rhoProbability =
            new Input<RealParameter>("rho", "Probability of an individual to be sampled at present", (RealParameter)null);

    // if the tree likelihood is condition on sampling at least one individual then set to true one of the inputs:
    public Input<Boolean> conditionOnSamplingInput = new Input<Boolean>("conditionOnSampling", "the tree " +
            "likelihood is conditioned on sampling at least one individual if condition on origin or at least one individual on both sides of the root if condition on root", false);
    public Input<Boolean> conditionOnRhoSamplingInput = new Input<Boolean>("conditionOnRhoSampling", "the tree " +
            "likelihood is conditioned on sampling at least one individual in present if condition on origin or at lease one extant individual on both sides of the root if condition on root", false);

    public Input<Boolean> conditionOnRootInput = new Input<Boolean>("conditionOnRoot", "the tree " +
            "likelihood is conditioned on the root height otherwise on the time of origin", false);

    public Input<Taxon> taxonInput = new Input<Taxon>("taxon", "a name of the taxon for which to calculate the prior probability of" +
            "being sampled ancestor under the model", (Taxon) null);

    public final Input<IntegerParameter> SATaxonInput = new Input<IntegerParameter>("SAtaxon", "A binary parameter is equal to zero " +
            "if the taxon is not a sampled ancestor (that is, it does not have sampled descendants) and to one " +
            "if it is a sampled ancestor (that is, it has sampled descendants)", (IntegerParameter)null);
"""


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
    package_notice = ("fbd", "Sampled Ancestors (SA)")

    def add_parameters(self, state):
        """ Add tree-related <state> sub-elements.

        Add real parameters for the following FBD parameters:
         - diversification rate
         - sampling proportion
         - turnover
         - removal probability
         - probability of an extant languages to be in the sample (ρ)
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
        # Probability of removal of a fossile after sampling – Fixed to 1 for
        # the FBD tree prior. The SA prior may modify this.
        xml.parameter(
            state,
            id="rFBD.t:beastlingTree", lower="0.0", name="stateNode", upper="1.0",
            text="1.0")
        xml.parameter(
            state,
            # TODO: Why estimate=False? This is adapted from the bears example
            # from SA. The parameter describes "Probability of an individual to
            # be sampled at present". So if we have no reason to assume we
            # sampled all languages, we should reflect that here.
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
            tree="@Tree.t:beastlingTree",
            # We have no knowledge when this diversification process started or
            # how many extant languages we should have been looking for, this
            # dataset is all we have.
            conditionOnRoot="true",
            # If there were less than 2 extant languages, we would not look at this tree.
            conditionOnRhoSampling="true",
            diversificationRate="@diversificationRateFBD.t:beastlingTree",
            samplingProportion="@samplingProportionFBD.t:beastlingTree",
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
        xml.log(tracer_logger, idref="samplingProportionFBD.t:beastlingTree")
        xml.log(tracer_logger, idref="turnoverFBD.t:beastlingTree")
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
        # Continue to assume complete sampling for contemporary leaves
        #         "parameter": "@rhoFBD.t:beastlingTree",


class SampledAncestorsTree (FossilizedBirthDeathTree):
    """A sampled ancestors tree is a generalization of the FBD tree.

    In a FBD tree, a thing either fossilizes or leaves descendants. In a SA
    tree, a proportion 0≤r≤1 of fossils also has later descendants.

    """
    __type__ = "sampled-ancestors"
    package_notice = ("Sampled Ancestors tree prior", "Sampled Ancestors (SA)")

    def add_operators(self, beastxml):
        super().add_operators(beastxml)
        xml.operator(
            beastxml.run,
            attrib={
                "id": "RemovalScaler.t:beastlingTree",
                "spec": "ScaleOperator",
                "parameter": "@rFBD.t:beastlingTree",
                "scaleFactor": "0.5",
                "weight": "3.0"})

    def add_fine_logging(self, tracer_logger):
        super().add_fine_logging(tracer_logger)
        xml.log(tracer_logger, idref="rFBD.t:beastlingTree")
