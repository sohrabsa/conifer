library(phytools)
library(ape)
library(coda)

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

?prop.clades




all.sub.trees <- NA
unique.sub.trees <- NA
counts <- NA

reload.context <- function() {
  all.sub.trees <<- get.sub.trees(trees)
  unique.sub.trees <<- unique(all.sub.trees)
  counts <<- get.count.array(unique.sub.trees, all.sub.trees)
}

find.all.clade.supports <- function(tree, trees, reload.context=F) {
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


counts <- list()
for (clade in subtrees(tree)) {
  index <- length(counts) + 1
  counts[index] <- 0
  
  for (subtree in all.sub.trees) {
    if (all.equal(clade, subtree, use.edge.length = F)) {
      counts[index] <- counts[index] + 1
    }
  }
}

counts

