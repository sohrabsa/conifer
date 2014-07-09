library(tools)
library(coda)
source("/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/conifer_fork/confier_fork/r.scripts/MrBayesBatch.R")

makeDataFileForMrBayes <- function(alignmentFile, treeFile) {
  if (file_ext(alignmentFile) != "nex") {
    # convert the alignment file from fasta to 
    alignment.nex.file <- paste0(file_path_sans_ext(alignmentFile), ".nex")
    commandString <- paste("seqmagick convert --output-format nexus --alphabet dna", alignmentFile, alignment.nex.file)
    system(commandString)
    alignmentFile <- alignment.nex.file
  }
  
  # append the tree
  theTree <- paste(readLines(treeFile), collapse = "")
  treeBlock <- unlist(list("", "begin trees;", paste0("tree mm =", theTree), "end;"))
  writeLines(treeBlock, "treeBlock.nex")
  commandString <- paste("cat", alignmentFile, "treeBlock.nex", ">", "datafile.nex")
  
  system(commandString)
}


compileBatchScriptForMrBayes <- function(data.file, bath.script.file.name) {
  mrbayes <- new("mrbayesbatch")
  setTags(mrbayes) <- list(command="set", list=list(autoclose="yes", nowarn="yes"))
  addMultiPartCommand(mrbayes) <- (list("execute", data.file))
  setTags(mrbayes) <- list(command="startvals", list=list(tau="mm", V="mm"))
  setTags(mrbayes) <- list(command="lset", list=list(nst=6, rates="Equal"))
  setTags(mrbayes) <- list(command="mcmc", list=list(ngen=10000, samplefreq="10"))
  setTags(mrbayes) <- list(command="sump", list=list(burnin=1000))
  setTags(mrbayes) <- list(command="sumt", list=list(burnin=1000))
  writeToDisk(mrbayes, bath.script.file.name)  
}

mrbayes.analysis <- function(batch.script) {
  system(paste0("mb ", batch.script), ignore.stdout = F)
  #system(paste0("mb ", batch.script), ignore.stdout = T)
}


mrbayes.calculate.ESS <- function(elapsed.time) {
  file.names <- list.files(".", pattern = "*\\.nex\\.run1\\.p")
  mrbayes.posterior <- read.table(file.names[1], skip=1, header=T, check.names=F)
  mrbayes.posterior <- mrbayes.posterior[, c(3:9)]
  
  ESS <- effectiveSize(mrbayes.posterior)
  ESSPS <- effectiveSize(mrbayes.posterior)/elapsed.time
  
  # save the values
  write.table(data.frame(ESS, check.names=F), "ess.txt")
  write.table(data.frame(ESSPS, check.names=F), "ess_per_second.txt")
}

mrbayes.driver.function <- function(tree.file.name, alignment.file.name) {
  # set dir where the mrbayes files should be placed  
  batch.dir <- "/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/conifer/extras/mrbayes/jul.analysis"
  setwd(batch.dir)
  
  batch.file.name <- "july.compile.batch"
  //alignment.file.name <- "FES_4.nex"
  //tree.file.name <- "FES.ape.4.nwk"
  makeDataFileForMrBayes(alignmentFile=alignment.file.name, treeFile=tree.file.name)
  compileBatchScriptForMrBayes(data.file="datafile.nex", bath.script.file.name=batch.file.name)
  
  # get the elapsed time
  the.run.time <- system.time(mrbayes.analysis(batch.file.name))
  elapsed.time <- the.run.time[3]
  
  # calculate ESS and ESS per second
  mrbayes.calculate.ESS(elapsed.time)
  
  # calculate the concensus tree
  
  
  # keep record of the analysis time
  writeLines(c(paste("Analysis for ", batch.file.naem, ", took", elapsed.time, "seconds.")), "experiment.details.txt")
}


consensus <- read.nexus("/Users/sohrab/Me/Apply/Canada\ Apply/Courses/Third\ Semester/conifer/extras/mrbayes/FES_8_batch.GTR.two/FES_8_batch.GTR.two.nex.tree1.con.tre")
all.trees <- read.nexus("/Users/sohrab/Me/Apply/Canada\ Apply/Courses/Third\ Semester/conifer/extras/mrbayes/FES_8_batch.GTR.two/FES_8_batch.GTR.two.nex.tree1.run1.t")
# get rid of the burn.in period
burn.in <- 1000
thining <- 10
burn.in <-  as.integer(burn.in / thining)
all.trees <- all.trees[-c(1:burn.in)]
N <- length(all.trees)

# get all sub.trees
all.sub.trees <- get.sub.trees(all.trees)

# get the unique ones
unique.sub.trees <- unique(all.sub.trees)

# get the counts
//counts <- get.count.array(unique.sub.trees, all.sub.trees)
counts <- get.count.array(all.sub.trees, all.sub.trees)

# calculate support for all clades of the consensus tree
get.count.array <- function(unique.sub.trees, all.sub.trees) {
  unlist(lapply( unique.sub.trees, function(st) sum(unlist(lapply(all.sub.trees, function(x) all.equal(st, x, use.edge.length = F))))   ))
}


# count the supports
counts <- list()
for (clade in subtrees(consensus)) {
  index <- length(counts) + 1
  counts[[index]] <- 0
  
  for (subtree in all.sub.trees) {
    
    if (all.equal(clade, subtree, use.edge.length = F)) {
      counts[[index]] <- counts[[index]] + 1
    }
  }
}

counts



