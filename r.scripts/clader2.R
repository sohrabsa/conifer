library(phytools)
library(ape)
library(coda)

# write these valuse as pecial node comments to be visualized by figtree (special node comments, ref:https://pythonhosted.org/DendroPy/scripts/sumtrees.html)
# what is majority-rule clade consensus tree 


get.sub.trees <- function(trees) {
  all.sub.trees <- unlist(lapply( trees, function(x) subtrees(x)), recursive=F)
  class(all.sub.trees) <- "multiPhylo"
  all.sub.trees
}

get.count.array <- function(unique.sub.trees, all.sub.trees) {
  unlist(lapply( unique.sub.trees, function(st) sum(unlist(lapply(all.sub.trees, function(x) all.equal(st, x, use.edge.length = F))))   ))
}

index.of.clade <- function(clade, unique.sub.trees) {
  for (i in seq(length(unique.sub.trees))) {
    if (all.equal(unique.sub.trees[[i]], clade, use.edge.length=F)) {
      return(i)
    }
  }
}

count.of.clade <- function(clade, unique.sub.trees, counts) {
  return(counts[index.of.clade(clade, unique.sub.trees)])
}


reload.context <- function() {
  all.sub.trees <<- get.sub.trees(trees)
  unique.sub.trees <<- unique(all.sub.trees)
  counts <<- get.count.array(unique.sub.trees, all.sub.trees)
}

all.sub.trees <- NA
unique.sub.trees <- NA
counts <- NA
supports <- function(tree, trees, reload.context=F) {
  if(reload.context) reload.context()
  unlist(lapply(seq(length(subtrees(tree))), function(j) count.of.clade(subtrees(tree)[[j]], unique.sub.trees, counts)))
}



get.count.for.tree <- function(tree, all.sub.trees) {
  counts <- list()
  for (clade in subtrees(tree)) {
    index <- length(counts) + 1
    counts[[index]] <- 0
    
    for (subtree in all.sub.trees) {
      if (all.equal(clade, subtree, use.edge.length = F)) {
        counts[[index]] <- counts[[index]] + 1
      }
    }
  }
  
  unlist(counts)
}
set.node.labels <- function(tree, counts) {
  s <- unlist(lapply(subtrees(tree), function(t)  min(t$node.label)))
  nodelabels(text=counts, node=s)
}

plot.side.by.side <- function(tree1, tree2) {
  par(mfrow=c(1,2))
  plot(tree1)
  nodelabels(tree1$node.label)
  #set.node.labels(tree1, get.count.for.tree(tree1, all.sub.trees))
  plot(tree2, direction = "leftwards")
  nodelabels(tree2$node.label)
  #set.node.labels(tree2, get.count.for.tree(tree1, all.sub.trees))
}
