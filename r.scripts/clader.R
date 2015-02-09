library(phytools)
library(ape)
library(coda)
library(distory)
# install.packages('distory')

# write these valuse as pecial node comments to be visualized by figtree (special node comments, ref:https://pythonhosted.org/DendroPy/scripts/sumtrees.html)
# what is majority-rule clade consensus tree 

trees <- pbtree(n = 6, nsim = 10, scale = 1)
N <- length(trees)

write.tree(trees, "~/Downloads/test.txt")



all.sub.trees <- get.sub.trees(trees)

clade.N <- length(all.sub.trees)
clade.N

length(unique(all.sub.trees))

# give me the confidence fore the clade


unique.sub.trees <- unique(all.sub.trees)
counts <- get.count.array(unique.sub.trees, all.sub.trees)

normalized.counts <- counts/clade.N
counts <- counts / Ntree
counts
q.strings <- unlist(lapply(unique(all.sub.trees), write.tree))
d <- data.frame(clades=q.strings, counts=counts)
write.csv(d, "clade.counts.csv", row.names=F)

# report the .95 percent confidence
o <- order(counts, decreasing = T)


t <- 0.0
treshold <- .95
for (i in seq(length(counts))) {
  cat("at (", i, ")", "sum is ", t + normalized.counts[o[i]], "\n")
  if (t + normalized.counts[o[i]] > treshold) break
  t <- t + normalized.counts[o[i]]
  cat("t is ", t, "\n")
}
t

hist(normalized.counts)


c <- (consensus(trees))
str(c)
print(c)

c


# given a majority tree, return the clade support for each node:


plot(read.tree("~/Downloads/majority.txt"))
plot(read.nexus("~/Downloads/majority.txt"))
str(read.nexus("~/Downloads/majority.txt"))
cophyloplot(trees[[1]], trees[[2]])
plot(consensus(trees), show.node.label=T)


c <- consensus(trees)
c1 <- subtrees(c)[[1]]
prop.clades(c, trees, rooted =T)
prop.clades(trees[[1]], trees, rooted =T)
subtrees(trees[[1]])

plot(c)
subtrees(c)

s <- subtrees(c)[[2]]
for (i in seq(length(all.sub.trees))) {
  if (all.equal(all.sub.trees[[i]], s, use.edge.length = F)) {
    print(i)
    break
  }
}



q <- unique(all.sub.trees)




system.time( {index.of.clade(s)})
system.time( {count.of.clade(s)})

system.time( {prop.clades(c, trees, rooted =T)})

index.of.clade(subtrees(c)[[1]])


q <- unique(all.sub.trees)
q <- lapply(q, function(x) {x$edge.length <- NULL; return(x)} )
class(q) <- "multiPhylo"
q
q.strings.clade <- unlist(lapply(q, function(x) write.tree(x, )  ))

q1 <- q.strings.clade[[5]]
q1.tree <- read.tree(text=q1)
plot(q1.tree)

q1


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


t <- trees[[2]] 




unlist(lapply(seq(length(subtrees(t))), function(j) count.of.clade(subtrees(t)[[j]])))

prop.part(t)
reload.context <- function() {
  all.sub.trees <<- get.sub.trees(trees)
  unique.sub.trees <<- unique(all.sub.trees)
  counts <<- get.count.array(unique.sub.trees, all.sub.trees)
}

find.all.clade.s

?prop.clades




all.sub.trees <- NA
unique.sub.trees <- NA
counts <- NA
upports <- function(tree, trees, reload.context=F) {
  if(reload.context) reload.context()
  unlist(lapply(seq(length(subtrees(tree))), function(j) count.of.clade(subtrees(tree)[[j]], unique.sub.trees, counts)))
}

trees <- pbtree(n = 6, nsim = 11, scale = 1, rooted=T)
tree <- trees[[5]]

find.all.clade.supports(tree, trees, reload.context=F)
prop.clades(tree, trees, rooted=T)

for (tree in trees) {
  print(find.all.clade.supports(tree, trees, reload.context=F))
  print(prop.part(tree, trees))
  print(prop.clades(tree, trees, rooted=T))
  print("---------")
}

trees


print(prop.part(tree, trees[c(1,2,3)]))

# for all possible subtrees s
  # for all subtrees of a tree t
    # if (s == t) count.t ++

all.sub.trees <- get.sub.trees(trees)


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

counts <- unlist(counts)

trees <- pbtree(n = 6, nsim = 10, scale = 1)

t1 <- trees[[1]]
t2 <- trees[[2]]
# visualize, given two trees, write the clade support on them and then create plot
par(mfrow=c(1,2))
plot(t2)
nodelabels()
plot(t1, direction = "leftwards")
ls()

plot(t2)
nodelabels(text=c("S"), node=c(7))
nodelabels()

# given counts, return the text and node arguments
plot(tre
# for each internal node
s <- subtrees(t2)
class(s) <- "multiPhylo"

str(s[[1]])
min(s[[1]]$node.label)
plot(s[[1]])
nodelabee)
set.node.labels <- function(tree, counts) {
  s <- unlist(lapply(subtrees(tree), function(t)  min(t$node.label)))
  nodelabels(text=counts, node=s)
}

plot(t2)
set.node.labels(tree, counts)



trees <- pbtree(n = 6, nsim = 11, scale = 1, rooted=T)
c <- consensus(trees)

write.tree(c)

all.sub.trees <<- get.sub.trees(trees)

tree1 <- trees[[1]]
tree2 <- trees[[2]]

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

add.support.to.tree <- function(tree, )


plot.side.by.side(tree1, tree2, all.sub.trees)


### experiment on distance between two phylo trees
## what model are we using for simulation? we used GTR to simulate, it had kimura

## original simulated tree
# /Users/sohrab/project/conifer/simulated.data/simulation.4_DEFAULT_DNAGTR/SimulatedDataTree.newick
original.tree <- read.tree('/Users/sohrab/project/conifer/simulated.data/simulation.4_DEFAULT_DNAGTR/SimulatedDataTree.newick')

other.conifer.tree <- read.tree('/Users/sohrab/project/conifer/results/all/2015-02-02-04-58-10-zDsuYnpo.exec/FESConsensusTree.Nexus')

# westrun trees
# '~/Downloads/configRuns.newick'
westrun.trees <- read.tree('~/Downloads/configRuns.newick')
westrun.tree <- westrun.trees[[1]]

original.tree$edge.length
westrun.tree$edge.length

sometree <- original.tree
sometree$edge.length <- original.tree$edge.length / max(original.tree$edge.length)
plot(sometree)
dist.topo(original.tree, sometree, method = 'score')

normalized.tree <- westrun.tree
normalized.tree$edge.length <- normalized.tree$edge.length / max(normalized.tree$edge.length)

plot(normalized.tree)


plottrees(sometree, normalized.tree, names=c('original', 'westrun'))
dist.topo(sometree, normalized.tree)

# normalize trees
west.trees.norm <- westrun.trees
l <- length(westrun.trees)
for (i in 1:l) {
  west.trees.norm[[i]]$edge.length <- west.trees.norm[[i]]$edge.length/ max(west.trees.norm[[i]]$edge.length)
}

plottrees(west.trees.norm[[45]], west.trees.norm[[33]], c('45', '33'))

# 33 is non trivial
## this actually works
dist.topo(sometree, west.trees.norm[[33]], 'score')
plottrees(sometree, west.trees.norm[[33]], c('original', '33'))
dist.topo(original.tree, westrun.trees[[33]], 'score')

dist.topo(original.tree, original.tree, method = 'score')
dist.topo(original.tree, westrun.tree, method = 'score')


# distroy
dist.multiPhylo(list(original.tree, westrun.trees[[33]]))

plottrees(westrun.tree, original.tree)



dist.topo(original.tree, westrun.tree, method = 'score')
dist.topo(westrun.trees[[2]], westrun.tree, method = 'score')
dist.topo(westrun.trees[[2]], westrun.tree)

dist.topo(westrun.trees[[2]], westrun.trees[[10]], 'score')

ll <- length(westrun.trees)
for (i in 1:ll) {
  for (j in 1:ll) {
    d <- (dist.topo(westrun.trees[[i]], westrun.trees[[j]]))
    if (d >= 1) print(paste(i, j, d))
  }
}
  

non.trivials.ids <- c()
non.trivials.dist <- c()

for (i in 1:ll) {
    d <- (dist.topo(westrun.trees[[i]], westrun.trees[[1]]))
    if (d >= 1) {
      print(paste(i, 1, d))
      d2 <- 1;
      tryCatch(d2 <<- dist.topo(westrun.trees[[i]], original.tree, 'score'), 
                  error = function(err) {
                 # error handler picks up where error was generated
                 print(paste("MY_ERROR:  ",err))
                 return(1)
               })
      
      non.trivials.ids <- c(non.trivials.ids, i) 
      non.trivials.dist <- c(non.trivials.dist, d2)
      print(paste(i, 1, d2))
    }
}
plottrees(sometree, west.trees.norm[[2]], c('original', '2'))

result <- data.frame(id=non.trivials.ids, dist=non.trivials.dist)
write.csv(result, '~/Desktop/west.run.dist.csv')

plottrees <- function(tree1, tree2, names) {
  par(mfrow=c(1,2))
  plot(tree1, main=names[1])
  nodelabels(tree1$node.label)
  #set.node.labels(tree1, get.count.for.tree(tree1, all.sub.trees))
  plot(tree2, direction = "leftwards", main=names[2])
  nodelabels(tree2$node.label)
  #set.node.labels(tree2, get.count.for.tree(tree1, all.sub.trees))
}

### show the trees
for (i in non.trivials.ids) {
  plottrees(westrun.trees[[i]], original.tree, c(i, 'original'))
}


