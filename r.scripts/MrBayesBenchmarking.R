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


driver.function <- function() {
  # set dir where the mrbayes files should be placed  
  batch.dir <- "/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/conifer/extras/mrbayes/jul.analysis"
  setwd(batch.dir)
  
  batch.file.name <- "july.compile.batch"
  data.file.name <- "FES_4.nex"
  tree.file.name <- "FES.ape.4.nwk"
  makeDataFileForMrBayes(alignmentFile=data.file.name, treeFile=tree.file.name)
  compileBatchScriptForMrBayes(data.file="datafile.nex", bath.script.file.name=batch.file.name)
  
  # get the elapsed time
  the.run.time <- system.time(mrbayes.analysis(batch.file.name))
  elapsed.time <- the.run.time[3]
  
  # calculate ESS and ESS per second
  file.names <- list.files(".", pattern = "*\\.nex\\.run1\\.p")
  mrbayes.posterior <- read.table(file.names[1], skip=1, header=T, check.names = F)
  mrbayes.posterior <- mrbayes.posterior[, c(3:9)]
  
  ESS <- effectiveSize(mrbayes.posterior)
  ESSPS <- effectiveSize(mrbayes.posterior)/elapsed.time
  
  # save the values
  write.table(data.frame(ESS, check.names = F), "ess.txt")
  write.table(data.frame(ESSPS, check.names = F), "ess_per_second.txt")
  
  writeLines(c(paste("Analysis took", elapsed.time, "seconds.")), "experiment.details.txt")
}
