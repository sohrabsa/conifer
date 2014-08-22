package conifer;

import java.io.File;

import bayonet.distributions.Exponential;
import bayonet.distributions.Exponential.RateParameterization;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.mcmc.Move;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;
import conifer.moves.SPRMove;
import conifer.moves.SingleNNI;

public class DNAModelFactory {

	public UnrootedTreeLikelihoodOptions options = null;

	public DNAModelFactory(UnrootedTreeLikelihoodOptions options) {
		this.options = options;
	}
	
	/**
	 * Fixed topology and fixed latent branch length
	 * @author Sohrab Salehi (sohrab.salehi@gmail.com)
	 *
	 */
	public class BasePhyloModel implements SimplePhyloModelContainer 
	{
		@DefineFactor(onObservations = true)
		public UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood = 
				UnrootedTreeLikelihood.fromOptions(options);
		
		@DefineFactor
		public final IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior =
		IIDRealVectorGenerativeFactor
		.iidNormalOn(likelihood.evolutionaryModel.rateMatrixMixture.parameters);
		
		public UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> getLikelihood() {
			return this.likelihood;
		}
	}
	
	public class FullLatentModel extends BasePhyloModel 
	{
		@DefineFactor
		NonClockTreePrior<RateParameterization> treePrior = 
		NonClockTreePrior
		.on(likelihood.tree);

		@DefineFactor
		Exponential<Exponential.MeanParameterization> branchLengthHyperPrior = 
		Exponential
		.on(treePrior.branchDistributionParameters.rate)
		.withMean(0.1);
	}
	
	SimplePhyloModelContainer model;

	public static SimplePhyloModelContainer getModelWithOptions(UnrootedTreeLikelihoodOptions options) 
	{
		DNAModelFactory mf = new DNAModelFactory(options);
		
		if (mf.options.topologyLatencyMode == TopologyLatencyModes.FIXED_TOPOLOGY_FIXED_BRANCH_LENGHTS) {
			mf.model = mf.new BasePhyloModel();
		} else {
			mf.model = mf.new FullLatentModel();
		}
			
		return mf.model; 
	}
	
	public static void main(String[] args) 
	{
		SimplePhyloModelContainer model = DNAModelFactory.getModelWithOptions(
				//UnrootedTreeLikelihoodOptions.makeWithOptions(TopologyUtils.makeLeaves(10, "s_"), 5, ExpFamMixture.dnaGTR(), null, null, TopologyLatencyModes.DEFAULT));
				UnrootedTreeLikelihoodOptions.makeForForwardSimulation(20, 100, ExpFamMixture.kimura1980(), TopologyLatencyModes.FIXED_TOPOLOGY_FIXED_BRANCH_LENGHTS)
				);
		
		final MCMCFactory factory = new MCMCFactory();
		MCMCAlgorithm mcmc = factory.build(model, false);
		System.out.println(mcmc.model);
		System.out.println();
		System.out.println(model.getLikelihood().observations.nSites());
	}
}
