package conifer;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import bayonet.distributions.Exponential;
import bayonet.distributions.Exponential.RateParameterization;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.ForwardSampler;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.mcmc.Move;
import blang.processing.Processor;
import blang.processing.ProcessorContext;
import briefj.BriefIO;
import briefj.opt.Option;
import briefj.opt.OptionSet;
import briefj.run.Mains;
import briefj.run.Results;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.io.FastaUtils;
import conifer.models.MultiCategorySubstitutionModel;
import conifer.moves.SPRMove;
import conifer.moves.SingleNNI;


public class TestPhyloModel implements Runnable, Processor
{
	@Option
	public String initialTreeFilePath;
	//= new File("/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES_4.fasta");

	@Option
	public String alignmentFilePath = "";

	@Option
	public int nMCMCSweeps = 100;

	@Option
	public int burnIn = (int) Math.round(.1 * nMCMCSweeps);

	@Option
	public int thinningPeriod = 10;

	@Option
	public StandardEvolutionaryModels evolutionaryModel = StandardEvolutionaryModels.DNAGTR;

	@OptionSet(name = "factory")
	public final MCMCFactory factory = new MCMCFactory();

	@Option
	public boolean fixedTopology = false;

	@Option
	public boolean fixedBranchLength = false;

	@Option
	public boolean plotCoda = false;
	

	// we can build different models here and instantiate them according to the input parameter
	// in the run block.

	// Note: only instantiate this in run to avoid problems with command line argument parsing
	public SimplePhyloModelContainer model;

	private PrintWriter treeWriter;
	private PrintWriter detailWriter;

	@SuppressWarnings("unchecked")
	public MCMCAlgorithm getMCMCAlgorithm(UnrootedTreeLikelihoodOptions options) {
		
		// set the evolutionary model
		model = DNAModelFactory.getModelWithOptions(options);		
		return factory.build(model, false);		
	}
	
	public void setupFactory() {
		factory.addProcessor(this);
		factory.mcmcOptions.nMCMCSweeps = nMCMCSweeps;
		factory.mcmcOptions.burnIn = burnIn;
		factory.mcmcOptions.CODA = plotCoda;
		factory.mcmcOptions.thinningPeriod = thinningPeriod;

		// TODO: move to where the model selection is being taken care of
		if (fixedTopology && !fixedBranchLength) {
			factory.excludeNodeMove((Class<? extends Move>) (Object)SingleNNI.class);
			factory.excludeNodeMove(SPRMove.class);
		}
	}
	
	//TODO: move all this to a factory/template class
	@SuppressWarnings("unchecked")
	@Override
	public void run()
	{
		long startTime = System.currentTimeMillis();

		treeWriter = BriefIO.output(Results.getFileInResultFolder("FES.trees.newick"));
		detailWriter = BriefIO.output(Results.getFileInResultFolder("experiment.details.txt"));

		setupFactory();
		// setup options for posterior sampling
		UnrootedTreeLikelihoodOptions options = 
				UnrootedTreeLikelihoodOptions.makeForPosteriorSimulation(evolutionaryModel, alignmentFilePath, initialTreeFilePath, fixedTopology, fixedBranchLength);
		MCMCAlgorithm mcmc = getMCMCAlgorithm(options);
		
		// log experiment information
		logModel(mcmc);
		System.out.println(mcmc.model);

		mcmc.run();

		// log running time
		logToFile("Total time in minutes: " + ((System.currentTimeMillis() - startTime) / 60000.0));

		// compute the consensus tree
		MajorityRuleTree.buildAndWriteConsensusTree(
				Results.getFileInResultFolder("FES.trees.newick"),
				Results.getFileInResultFolder("FESConsensusTree.Nexus"));
	}
	
	public static void main(String [] args) throws ClassNotFoundException, IOException
	{
		System.out.println("Running the newest version...");
		
		args = new String[8];
		args[0] = "-initialTreeFilePath";
		//		args[1] = "/Users/sohrab/project/conifer/src/main/resources/conifer/sampleInput/FES.ape.4.nwk";
		args[1] = "/Users/sohrab/project/conifer/simulated.data/5_DEFAULT/SimulatedDataTree.newick";
		args[2] = "-alignmentFilePath";
		//		args[3] = "/Users/sohrab/project/conifer/src/main/resources/conifer/sampleInput/FES_4.fasta";
		args[3] = "/Users/sohrab/project/conifer/simulated.data/5_DEFAULT/SimulatedData.fasta";
		args[4] = "-fixedTopology";
		args[5] = "true";
		args[6] = "-evolutionaryModel";
		args[7] = "KIMURA1980";
		 

//		args = new String[1];
//		args[0] = "-help";

		// TODO: remove this
		for (int i = 0; i < args.length; i++) {
			System.out.println(args[i]);	
		}

		System.out.println(args.length);

		Mains.instrumentedRun(args, new TestPhyloModel());

		// TODO: test if all topologies are the same
		// TODO: test if all branch lengths are the same/different 
	}


	public void writeTree(UnrootedTree tree) {
		PrintWriter theWriter = BriefIO.output(Results.getFileInResultFolder("SimulatedDataTree.newick"));
		theWriter.println(tree.toNewick());
		theWriter.flush();
	}


	@Override
	public void process(ProcessorContext context)
	{
		System.out.println("Writing the tree...");
		treeWriter.println(model.getLikelihood().tree.toNewick());
		treeWriter.flush();
	}

	public void logToFile(String someline) {
		this.detailWriter.println(someline);
		this.detailWriter.flush();
	}

	public String getAlignmentFile() {
		return alignmentFilePath;
	}

	void logModel(MCMCAlgorithm mcmc) {
		// TODO: move this to a model logger method/class
		logToFile("Over AlignmentFile:" + getAlignmentFile());
		logToFile("Current Date:" + RunFacility.getCurrentDateString());
		logToFile("Over AlignmentFile:" + getAlignmentFile());
		logToFile("");
		logToFile("MCMC burnIn:\n" + factory.mcmcOptions.burnIn);
		logToFile("MCMC nMCMCSweeps:\n" + factory.mcmcOptions.nMCMCSweeps);
		logToFile("MCMC thinningPeriod:\n" + factory.mcmcOptions.thinningPeriod);
		logToFile("Fixed Topology:" + this.fixedTopology + "\n");
		logToFile("Fixed BranchLength:" + this.fixedBranchLength + "\n");
		logToFile("");
		logToFile("Model:\n" + mcmc.model);
	}

	public void makeSyntheticData(UnrootedTreeLikelihoodOptions options) throws IOException {

		detailWriter = BriefIO.output(Results.getFileInResultFolder("experiment.details.txt"));

		MCMCAlgorithm mcmc = getMCMCAlgorithm(options);
		ForwardSampler f = new ForwardSampler(mcmc.model);

		Map<Object,List<Double>> results = Maps.newHashMap();
		f.simulate(mcmc.options.random);
		writeTree(model.getLikelihood().tree);
			// write the FASTA file corresponding to the simulation
			FastaUtils.writeFasta(model.getLikelihood().observations, Results.getFileInResultFolder("SimulatedData.fasta"));
		}
		
		detailWriter.write("Total BranchLength: " +
		UnrootedTreeUtils.totalTreeLength((runner.model.getLikelihood().tree)) + "\n");

		// TODO: get the rest of the simulation parameters (total_branch_length?)
		for (Object var: algo.model.getLatentVariables()) {
			detailWriter.write(mcmc.model.getName(var) + " : " + var.toString() + "\n");
			System.out.println(mcmc.model.getName(var) + " : " + var.toString());
		}
		
		runner.detailWriter.flush();
	}
	
}
