package conifer.io;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;

import briefj.BriefCollections;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import conifer.TreeNode;



public class FixedTreeObservations implements TreeObservations
{
  private final LinkedHashMap<TreeNode, double[][]> data = Maps.newLinkedHashMap();
  private final int nSites;
  
  public FixedTreeObservations(int nSites)
  {
    this.nSites = nSites;
  }

  @Override
  public List<TreeNode> getObservedTreeNodes()
  {
    List<TreeNode> result = Lists.newArrayList();
    for (TreeNode key : data.keySet())
      result.add(key);
    return result;
  }

  @Override
  public Object get(TreeNode leaf)
  {
    return data.get(leaf);
  }

  @Override
  public void set(TreeNode leaf, Object point)
  {
    double[][] cast = (double[][])point;
    if (cast.length != nSites)
      throw new RuntimeException();
    data.put(leaf, cast);
  }

  @Override
  public void clear()
  {
    data.clear();
  }

  @Override
  public String toString()
  {
    StringBuilder result = new StringBuilder();
    for (TreeNode node : data.keySet())
      result.append(node.toString() + " : " + Arrays.deepToString(data.get(node)) + " ");
    return result.toString();
  }

  @Override
  public int nSites()
  {
    return nSites;
  }
}
