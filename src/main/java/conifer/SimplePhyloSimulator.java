package conifer;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;

import briefj.run.Results;

/**
 * This class will simulate data from TestPhyloModel
 * @author Sohrab Salehi (sohrab.salehi@gmail.com)
 *
 */
public class SimplePhyloSimulator {
	
		public SimplePhyloModelContainer model;

	public static void main(String [] args) throws ClassNotFoundException, IOException
	{
		int[] numberOfTaxa = {5, 10, 50};
		int nSites = 10;
		StandardEvolutionaryModels evolutionaryModel = StandardEvolutionaryModels.DNAGTR;
		
		TestPhyloModel runner;
		
		for (int ntaxa : numberOfTaxa) {
			for (TopologyLatencyModes mode :  TopologyLatencyModes.values()) {
				runner = new TestPhyloModel();
				runner.makeSyntheticData(UnrootedTreeLikelihoodOptions.makeForForwardSimulation(ntaxa, nSites, evolutionaryModel.getExpFamily(), mode));
				
				// copy the results to another folder 
				copyDirectory(ntaxa, mode);
			}
		}
	}

	/**
	 * copy the contents of the current directory to a new one to prevent overriding of the new experiments.
	 * @throws IOException 
	 */
	public static void copyDirectory(int ntaxa, TopologyLatencyModes mode) throws IOException {
		File newDirectory = new File(
				Results.getResultFolder().getParent() + "/experiment." + 
		Results.getResultFolder().getName() + "." + 
		ntaxa + "_" + mode.toString() + "_" + System.currentTimeMillis());
		newDirectory.mkdir();
		FileUtils.copyDirectory(Results.getResultFolder(), newDirectory);
		System.out.println(newDirectory.getAbsolutePath());
	}

}
