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
  s <- unlist(lapply(subtrees(tree), function(t) min(t$node.label)))
  nodelabels(text=counts, node=s)
}

plot.side.by.side <- function(tree1, tree2, names, plotPath) {
  jpeg(plotPath, width=1250, height=460)
  par(mfrow=c(1,2))
  plot(tree1, main=names[1])
  nodelabels(tree1$node.label)
  #set.node.labels(tree1, get.count.for.tree(tree1, all.sub.trees))
  plot(tree2, direction = "leftwards", main=names[2])
  nodelabels(tree2$node.label)
  #set.node.labels(tree2, get.count.for.tree(tree1, all.sub.trees))
  
  dev.off()
}


# assumes unique tip.labels
random.tree.from.fasta <- function(fastaPath) {
  p <- readLines(fastaPath)

  # choose lines that start with >
  t <- p[grep(">.*", p)]
  species.names <- gsub(">", "", t)
  species.names
  
  tree <- rtree(length(species.names), rooted=F, tip.label=species.names)
  outputFile <- paste0(file_path_sans_ext(fastaPath), ".nwk")
  write.tree(tree, outputFile)
  
  outputFile
}


# keep only sequences for unique species
trim.fasta.file <- function(fastaPath) {
  p <- readLines(fastaPath)
  p <- gsub(">[^|]+\\|[^|]+\\|[^|]+\\|[^|]+\\|\ (.*)\ isolate\ .*", ">\\1", p)
  
  ## remove duplicates
  # species' names
  t <- p[grep(">.*", p)]
  s.i <- grep(">.*", p)
  s.e <- c((s.i - 1)[-c(1)], length(s.i))
  p[s.i[1]:s.e[1]]
  p[s.i[2]:s.e[2]]
  
  duplicated <- unlist(lapply(which(duplicated(t)), function(x)  s.i[x]:s.e[x]))
  
  p <- p[-(duplicated)]
  
  
  writeLines(p, paste0(file_path_sans_ext(fastaPath), "_trimmed.fasta"))
}


#fastaPath <- "/home/sohrab/Downloads/FES_full.fasta"
#fastaPath <- "/home/sohrab/Downloads/UTY_full.fasta"

#trim.fasta.file(fastaPath)

#fastaPath <- "/home/sohrab/Downloads/FES_full_trimmed.fasta"
#fastaPath <- "/home/sohrab/Downloads/UTY_full_trimmed.fasta"
#random.tree.from.fasta(fastaPath)

#/home/sohrab/conifer_fork/src/main/resources/conifer/sampleInput/UTY_full_trimmed.nwk