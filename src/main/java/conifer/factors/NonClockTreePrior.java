package conifer.factors;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import com.google.common.collect.Lists;


import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.UnrootedTree;

import bayonet.distributions.Exponential;
import bayonet.distributions.Exponential.RateParameterization;
import bayonet.distributions.UnivariateRealDistribution;
import blang.annotations.FactorArgument;
import blang.annotations.FactorComponent;
import blang.factors.GenerativeFactor;
import blang.variables.RealVariable;



public class NonClockTreePrior<P> implements GenerativeFactor
{
  @FactorArgument(makeStochastic = true)
  public final UnrootedTree tree;
  
  @FactorComponent
  public final P branchDistributionParameters;
  
  private final UnivariateRealDistribution branchDistribution;

  
  public NonClockTreePrior(
      UnrootedTree tree,
      P branchDistributionParameters,
      UnivariateRealDistribution branchDistribution)
  {
    this.tree = tree;
    this.branchDistribution = branchDistribution;
    this.branchDistributionParameters = branchDistributionParameters;
  }

  /**
   * TODO: missing a normalization factor with respect to the 
   * number of topologies, but since this parameter is generally not changed
   * in application, can leave out for now.
   */
  @Override
  public double logDensity()
  {
    double result = 0.0;
    for (double edgeLength : tree.getBranchLengths().values())
      result += branchLengthLogDensity(edgeLength);
    return result;
  }
  
  private double branchLengthLogDensity(double length)
  {
    branchDistribution.getRealization().setValue(length);
    return branchDistribution.logDensity();
  }
  
  public static UnrootedTree generate(
      Random random, 
      UnivariateRealDistribution branchDistribution,
      Collection<TreeNode> leaves)
  {
    UnrootedTree result = new UnrootedTree();
    for (TreeNode leaf : leaves)
      result.addNode(leaf);
    List<TreeNode> roots = Lists.newArrayList(leaves);
    loop:while (roots.size() > 1)
    {
      
      if (roots.size() == 2)
      {
        result.addEdge(roots.get(0), roots.get(1), sample(branchDistribution, random));
        break loop;
      }
      else
      {
        Pair<TreeNode,TreeNode> pair = popRandomPair(random, roots);
        TreeNode newInternal = TreeNode.nextUnlabelled();
        result.addNode(newInternal);
        result.addEdge(pair.getLeft(), newInternal, sample(branchDistribution, random));
        result.addEdge(pair.getRight(), newInternal, sample(branchDistribution, random));
        roots.add(newInternal);
      }
    }
    return result;
  }
  
  private static double sample(UnivariateRealDistribution distribution, Random rand)
  {
    ((GenerativeFactor) distribution).generate(rand);
    return distribution.getRealization().getValue();
  }
  
  public static <T> Pair<T,T> popRandomPair(Random rand, List<T> items)
  {
    if (items.size() < 2) 
      throw new RuntimeException();
    List<Integer> indices = sampleWithoutReplacement(rand, items.size(), 2);
    Collections.sort(indices);
    Pair<T,T> result = Pair.of(items.get(indices.get(0)), items.get(indices.get(1)));
    items.remove((int)indices.get(1)); // remove second before first to avoid shifts
    items.remove((int)indices.get(0));
    return result;
  }
  
  /**
   * Returns n samples from {0, ..., s-1}, without replacement
   * 
   * TODO: could be more efficient
   * TODO: move to Bayonet
   * 
   * @param n
   * @param s
   * @return
   */
  public static List<Integer> sampleWithoutReplacement(Random rand, int s, int n)
  {
    if (n > s || s < 0 || n < 0)
      throw new RuntimeException();
    List<Integer> list = new ArrayList<Integer>(s),
                  result=new ArrayList<Integer>(n);
    for (int i = 0; i < s; i++)
      list.add(i);
    Collections.shuffle(list,rand);
    for (int i = 0; i < n; i++)
      result.add(list.get(i));
    return result;
  }
  
  public static NonClockTreePrior<RateParameterization> on(UnrootedTree tree)
  {
    Exponential<RateParameterization> branchDist = Exponential.on(RealVariable.real());
    return new NonClockTreePrior<RateParameterization>(tree,branchDist.parameters, branchDist);
  }

  /**
   * Creates a new tree in place according to the prior.
   * The labels at the leaves of the current tree will be the same as
   * those before.
   */
  @Override
  public void generate(Random random)
  {
    List<TreeNode> leaves = TopologyUtils.leaves(tree.getTopology());
    UnrootedTree sampled = generate(random, branchDistribution, leaves);
    this.tree.setTo(sampled);
  }
  
}
