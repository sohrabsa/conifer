package conifer.processors;

import conifer.UnrootedTree;
import blang.processing.NodeProcessor;
import blang.processing.ProcessorContext;
import blang.variables.RealValued;


/**
 * Computes the sum of the branch lengths of a tree.
 * Used for testing.
 * 
 * @author Alexandre Bouchard (alexandre.bouchard@gmail.com)
 *
 */
public class TotalTreeLengthProcessor implements NodeProcessor<UnrootedTree>, RealValued
{
  private UnrootedTree tree;
  private double currentValue;

  @Override
  public void process(ProcessorContext context)
  {
    double sum = 0.0;
    for (double branchLength : tree.getBranchLengths().values())
      sum += branchLength;
    this.currentValue = sum;
  }

  @Override
  public void setReference(UnrootedTree variable)
  {
    this.tree = variable;
  }

  @Override
  public double getValue()
  {
    return currentValue;
  }

}
