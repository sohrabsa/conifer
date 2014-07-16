package conifer;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;


import org.apache.commons.io.FileUtils;

import bayonet.distributions.Exponential;
import bayonet.distributions.Exponential.RateParameterization;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.MCMCAlgorithm;
import blang.MCMCFactory;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
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
import conifer.models.MultiCategorySubstitutionModel;




public class TestPhyloModel implements Runnable, Processor
{
	@Option
	public String initialTreeFilePath;
	//= new File("/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES_4.fasta");

	@Option
	public String alignmentFilePath = "";
	
	@OptionSet(name = "factory")
	public final MCMCFactory factory = new MCMCFactory();
	
	public class Model
	{
		@DefineFactor(onObservations = true)
		public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood = 
		UnrootedTreeLikelihood
		.fromFastaFile(new File(alignmentFilePath))
		.withExpFamMixture(ExpFamMixture.kimura1980())
		//.withTree(new File("/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES.ape.4.nwk"));
		.withTree(new File(initialTreeFilePath));
		
		@DefineFactor
		NonClockTreePrior<RateParameterization> treePrior = 
		NonClockTreePrior
		.on(likelihood.tree);

		@DefineFactor
		Exponential<Exponential.MeanParameterization> branchLengthHyperPrior = 
		Exponential
		.on(treePrior.branchDistributionParameters.rate)
		.withMean(10.0);

		@DefineFactor
		public final IIDRealVectorGenerativeFactor<MeanVarianceParameterization> prior =
		IIDRealVectorGenerativeFactor
		.iidNormalOn(likelihood.evolutionaryModel.rateMatrixMixture.parameters);		
	}
	
	// Note: only instantiate this in run to avoid problems with command line argument parsing
	public Model model;
	
	private final PrintWriter treeWriter = BriefIO.output(Results.getFileInResultFolder("FES.trees.newick"));

	private final PrintWriter detailWriter = BriefIO.output(Results.getFileInResultFolder("experiment.details.txt"));

	@Override
	 public void run()
	 {
		logToFile("Over AlignmentFile:" + getAlignmentFile());
		
		factory.addProcessor(this);
		factory.mcmcOptions.nMCMCSweeps = 10000;
		factory.mcmcOptions.burnIn = (int) Math.round(.1 * factory.mcmcOptions.nMCMCSweeps);
		
	    model = new Model();
	    MCMCAlgorithm mcmc = factory.build(model, false);
	    System.out.println(mcmc.model);
	    
	    long startTime = System.currentTimeMillis();
		String excluding = "Simple";

		// log experiment information
		logToFile("Experiment Title:" + excluding);
		logToFile("Current Date:" + RunFacility.getCurrentDateString());
		logToFile("Over AlignmentFile:" + getAlignmentFile());
		logToFile("");
		logToFile("MCMC burnIn:\n" + factory.mcmcOptions.burnIn);
		logToFile("MCMC nMCMCSweeps:\n" + factory.mcmcOptions.nMCMCSweeps);
		logToFile("MCMC thinningPeriod:\n" + factory.mcmcOptions.thinningPeriod);
		logToFile("");
		logToFile("Model:\n" + mcmc.model);
	    
	    mcmc.run();
	    
	    // log running time
	    logToFile("Total time in minutes: " + ((System.currentTimeMillis() - startTime)/60000.0));

		// compute the tree
		MajorityRuleTree.buildAndWriteConsensusTree(
				Results.getFileInResultFolder("FES.trees.newick"),
				Results.getFileInResultFolder("FESConsensusTree.Nexus"));

		// copy the results to another folder 
		File newDirectory = new File(Results.getResultFolder().getParent() + "/experiment." + Results.getResultFolder().getName() + "." + 1 + "." + System.currentTimeMillis());
		newDirectory.mkdir();
		try {
			FileUtils.copyDirectory(Results.getResultFolder(), newDirectory);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	 }
	
	
	
	public static void main(String [] args) throws ClassNotFoundException, IOException
	{
		System.out.println("Running the new version...");
		
//		args = new String[4];
//		args[0] = "-initialTreeFilePath";
//		args[1] = "/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES.ape.4.nwk";
//		args[2] = "-alignmentFile";
//		args[3] = "/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES_4.fasta";
//		
		//args = new String[1];
		//args[0] = "-help";
		
		// TODO: remove this
		for (int i = 0; i < args.length; i++) {
			System.out.println(args[i]);	
		}
		
		System.out.println(args.length);
		
		Mains.instrumentedRun(args, new TestPhyloModel());
	}
	
	@Override
	public void process(ProcessorContext context)
	{
		treeWriter.println(model.likelihood.tree.toNewick());
		treeWriter.flush();
	}

	public void logToFile(String someline) {
		this.detailWriter.println(someline);
		this.detailWriter.flush();
	}

	public String getAlignmentFile() {
		return alignmentFilePath;
	}
}
