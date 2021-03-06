package conifer;

import java.io.File;
import java.io.PrintWriter;

import bayonet.distributions.Exponential;
import bayonet.distributions.Exponential.RateParameterization;
import bayonet.distributions.Normal.MeanVarianceParameterization;
import blang.MCMCRunner;
import blang.annotations.DefineFactor;
import blang.factors.IIDRealVectorGenerativeFactor;
import blang.processing.ProcessorContext;
import briefj.BriefIO;
import briefj.run.Results;
import conifer.ctmc.expfam.ExpFamMixture;
import conifer.factors.NonClockTreePrior;
import conifer.factors.UnrootedTreeLikelihood;
import conifer.models.MultiCategorySubstitutionModel;



public class SimpleProteinModel extends MCMCRunner
{
  File inputFile 
    = new File("/Users/crystal/Dropbox/protein/javaProteinTry/alignment.fasta");
  
  @DefineFactor(onObservations = true)
  public final UnrootedTreeLikelihood<MultiCategorySubstitutionModel<ExpFamMixture>> likelihood = 
    UnrootedTreeLikelihood
    .fromFastaProteinFile(inputFile)
    .withExpFamMixture(ExpFamMixture.proteinGTR())
    .withTree(new File("/Users/crystal/Dropbox/protein/javaProteinTry/tree.nwk"));
  
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
  
  private final PrintWriter treeWriter = BriefIO.output(Results.getFileInResultFolder("/Users/crystal/Dropbox/protein/javaProteinTry"));
  
  public static void main(String [] args)
  {
    SimpleProteinModel runner = new SimpleProteinModel();
    runner.factory.mcmcOptions.nMCMCSweeps = 100;
    runner.run();
    
  }

  protected void process(ProcessorContext context)
  {
    treeWriter.println(likelihood.tree.toNewick());
    treeWriter.flush();
  }

}
