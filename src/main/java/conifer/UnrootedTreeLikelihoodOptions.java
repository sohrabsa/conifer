package conifer;

import java.util.List;

import conifer.ctmc.expfam.ExpFamMixture;


/**
 * Collection of options to create an UnrootedTreeLikelihood object
 * 
 *  
 * 
 * @author Sohrab Salehi (sohrab.salehi@gmail.com)
 *
 */
public class UnrootedTreeLikelihoodOptions {

	public List<TreeNode> leaves = null;
	public int nSites = 0;
	public ExpFamMixture mixtureModel = null;
	public String alignmentFilePath = null;
	public String initialTreeFilePath = null;
	public TopologyLatencyModes topologyLatencyMode;

	public UnrootedTreeLikelihoodOptions(List<TreeNode> leaves, int nSites, ExpFamMixture mixtureModel, String alignmentFilePath, String initialTreeFilePath, TopologyLatencyModes topologyLatencyMode) 
	{
		this.leaves = leaves;
		this.nSites = nSites;
		this.mixtureModel = mixtureModel;
		this.alignmentFilePath = alignmentFilePath;
		this.initialTreeFilePath = initialTreeFilePath;
		this.topologyLatencyMode = topologyLatencyMode;
	}

	public static UnrootedTreeLikelihoodOptions makeWithOptions(List<TreeNode> leaves, int nSites, ExpFamMixture mixtureModel, String alignmentFilePath, String initialTreeFilePath, TopologyLatencyModes topologyLatencyMode) {
		return new UnrootedTreeLikelihoodOptions(leaves, nSites, mixtureModel, alignmentFilePath, initialTreeFilePath, topologyLatencyMode);
	}

	public static UnrootedTreeLikelihoodOptions makeForPosteriorSimulation(StandardEvolutionaryModels evolutionaryModel, String alignmentFilePath, String initialTreeFilePath, 
			boolean fixedTopology, boolean fixedBranchLength) {
		
		// set the mixture model
		ExpFamMixture mixtureModel = evolutionaryModel.getExpFamily();
		
		// set the topology latency mode
		TopologyLatencyModes topologyLatencyMode = TopologyLatencyModes.initFromBooleans(fixedTopology, fixedBranchLength);
		
		return new UnrootedTreeLikelihoodOptions(null, 0, mixtureModel, alignmentFilePath, initialTreeFilePath, topologyLatencyMode);
	}
	
	public static UnrootedTreeLikelihoodOptions makeForForwardSimulation(List<TreeNode> leaves, int nSites, ExpFamMixture mixtureModel, TopologyLatencyModes topologyLatencyMode) {
		return new UnrootedTreeLikelihoodOptions(leaves, nSites, mixtureModel, null, null, topologyLatencyMode);
	}
	
	public static UnrootedTreeLikelihoodOptions makeForForwardSimulation(int nTaxa, int nSites,  ExpFamMixture mixtureModel, TopologyLatencyModes topologyLatencyMode) {
		return UnrootedTreeLikelihoodOptions.makeForForwardSimulation(TopologyUtils.makeLeaves(nTaxa, "t_"), nSites, mixtureModel, topologyLatencyMode);
	}	
}
