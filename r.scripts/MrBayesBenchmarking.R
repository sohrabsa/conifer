library(tools)
library(coda)
source(file.path(mainDIR, "MrBayesBatch.R"))
source(file.path(mainDIR, "clader2.R"))

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
  file.names <- list.files(".", pattern = "*.run1.p")
  mrbayes.posterior <- read.table(file.names[1], skip=1, header=T, check.names=F)
  mrbayes.posterior <- mrbayes.posterior[, c(3:9)]
  
  ESS <- effectiveSize(mrbayes.posterior)
  
  # standardize column names
  ESS <- mrbayes.standardizeColumns.ESS(ESS)
  
  ESSPS <- effectiveSize(mrbayes.posterior)/elapsed.time
  
  # save the values
  write.table(data.frame(ESS, check.names=F), "ess.txt")
  write.table(data.frame(ESSPS, check.names=F), "ess_per_second.txt")
}

# map column names from mrbayes's output to those from conifer
mrbayes.standardizeColumns.ESS <- function(ESS) {
  
  b <- names(ESS)
  # change from r(A<->C) to q(A(0),C(0))
  names(ESS) <- gsub("r\\(", "q\\(", names(ESS))
  names(ESS) <- gsub("<->", ",", names(ESS))
  names(ESS) <- gsub("([ACTG])", "\1(0)", names(ESS), perl = T)
  names(ESS) <- gsub("([ACTG])", "\\1(0)", names(ESS), perl = T)
  ESS
}


# TODO save support as node.lables
mrbayes.calculate.consensus.tree <- function(burn.in, thinning) {
  file.names <- list.files(".", pattern = "*.run1.t")
  all.trees <- read.nexus(file.names[1])
  
  # get rid of the burn.in period
  burn.in <-as.integer(burn.in / thinning)
  all.trees <- all.trees[-c(1:burn.in)]
  N <- length(all.trees)
  
  # calculate the consensus tree
  consensus.tree <- consensus(all.trees)

  # calculate the clade support for the consensus tree
  all.sub.trees <- get.sub.trees(all.trees)
  counts <- get.count.for.tree(consensus.tree, all.sub.trees)
  
  # write clade support to the disk
  support.strings <- unlist(lapply(subtrees(consensus.tree), write.tree))
  d <- data.frame(clades=support.strings, counts=counts)
  write.csv(d, "counsensus.clade.counts.csv", row.names=F)
  
  # write the consensus tree to the disk
  # this will write clade support as [internal] node labels 
  consensus.tree$node.label <- counts
  write.tree(consensus.tree, "consensus.tree")
}


mrbayes.driver.function <- function(tree.file.name, alignment.file.name) {
  # set dir where the mrbayes files should be placed  
  batch.dir <- "/home/sohrab/conifer_fork/mrbayes"
  setwd(batch.dir)
  
  batch.file.name <- "july.compile.batch"
  alignment.file.name <- "FES_4.fasta"
  tree.file.name <- "FES.ape.4.nwk"
  
  # make symbolik links to the data.files
  system(paste0("ln -s ", file.path(dataDir, alignment.file.name), " ", batch.dir))
  system(paste0("ln -s ", file.path(dataDir, tree.file.name), " ", batch.dir))
  
  makeDataFileForMrBayes(alignmentFile=alignment.file.name, treeFile=tree.file.name)
  compileBatchScriptForMrBayes(data.file="datafile.nex", bath.script.file.name=batch.file.name)
  
  # get the elapsed time
  the.run.time <- system.time(mrbayes.analysis(batch.file.name))
  elapsed.time <- the.run.time[3]
  
  # calculate ESS and ESS per second
  mrbayes.calculate.ESS(elapsed.time)
  
  # calculate the concensus tree
  mrbayes.calculate.consensus.tree(100, 10)
  
  # keep record of the analysis time
  writeLines(c(paste("Analysis for ", batch.file.name, ", took", elapsed.time, "seconds.")), "experiment.details.txt")
}


mrbayes.driver.function("FES.ape.4.nwk", "FES_4.fasta")

#consensus <- read.nexus("/Users/sohrab/Me/Apply/Canada\ Apply/Courses/Third\ Semester/conifer/extras/mrbayes/FES_8_batch.GTR.two/FES_8_batch.GTR.two.nex.tree1.con.tre")