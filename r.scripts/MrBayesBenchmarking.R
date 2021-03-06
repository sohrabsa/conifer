library(tools)
library(coda)
source(file.path(mainDIR, "MrBayesBatch.R"))
source(file.path(mainDIR, "clader2.R"))

MRBAYES_CONSENSUS_TREE_PATH <- ""
MRBAYES_ESS_PATH <- ""
MRBAYES_ESSPERSECOND_PATH <- ""
MRBAYES_EXPERIMENT_PATH <- ""

mrbayes.load.consensus.tree <- function() {
  read.tree(MRBAYES_CONSENSUS_TREE_PATH)
}

mrbayes.load.ess <- function() {
  read.table(MRBAYES_ESS_PATH)
}

mrbayes.load.esspersecond <- function() {
  read.table(MRBAYES_ESSPERSECOND_PATH)
}

mrbayes.make.symlink <- function(target.dir) {
  system(paste0("ln -s ", MRBAYES_EXPERIMENT_PATH, " ", target.dir))
}


makeDataFileForMrBayes <- function(alignmentFile, treeFile) {
  if (file_ext(alignmentFile) != "nex") {
    # convert the alignment file from fasta to 
    alignment.nex.file <- paste0(file_path_sans_ext(alignmentFile), ".nex")
    commandString <- paste("seqmagick convert --output-format nexus --alphabet dna", alignmentFile, " ", alignment.nex.file)
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


# source: e.g. look at:
# http://mrbayes.sourceforge.net/wiki/index.php/Evolutionary_Models_Implemented_in_MrBayes_3
# http://bodegaphylo.wikispot.org/MrBayes_Tutorial_(Brown)
# Assumption: ALL SITE RATES ARE EQUAL
# create different phylogenetic models: (rate matrix madules)
# according to Probabilistic Graphical Model Representation in Phylogenetics Hohna et al., 2013
# 1. Jukes-Cantor (JC69): all exchangability rates and base frequencie are equal. -> free parameter \mu
# 2. Flenstein (F81): equal exchangability rates but different base frequencies drawn from a dirichlet -> free parameters \mu and \pi_1 ... \pi_4  
# 3. Hasegawa (HKY85): different transition transversion valuse and base frequencies drawn from a dirichlet -> freee parameters: \kappa
# 4. Kimura (K80): base different transition tranversion value but equal base frequencies -> free parameter -> \kappa
# 5. GTR (Tavare 86): base (stationary) drawn from a dirichlet as well as the exchangability rates -> free parapeters -> \alpha and \beta for the dirichlet

# Kappa for mrbayes and branch lengths: http://hydrodictyon.eeb.uconn.edu/eebedia/index.php/Phylogenetics:_MrBayes_Lab
# TL is sum of all branch lengths
# LnL is the log likelihood of the cold chain

# 1. fixed topo and fixed branch length
# 2. fixed topo and latent branch length


# exponential in mrbayes is parametrized by rate, so exponential(10) has an expected value of .1 
# ref: http://mrbayes.sourceforge.net/wiki/index.php/Evolutionary_Models_Implemented_in_MrBayes_3

# TODO: FIX rate to mean in mrbayes

mrbayes.possible.models <- function() c("JC69", "F81", "HKY85", "K80", "GTR")
# ref: http://mrbayes.sourceforge.net/wiki/index.php/FAQ-3.2#How_do_I_fix_the_tree_topology_during_an_analysis.3F
# fixed topoly:
# MrBayes> showmoves
# ...
# 5 -- Move     = eTBR(Tau)
# Type        = Extending TBR
# Parameters  = Tau [param. 5] (Topology)
# V [param. 6] (Branch lengths)
# Tuningparam = p_ext (Extension probability)
# lambda (Multiplier tuning parameter)
# p_ext = 0.800
# lambda = 0.940
# Rel. prob.  = 15.0
# MrBayes> propset eTBR(Tau)$prob=0
# MrBayes> showmoves allavailable=yes
# ...
# 7 -- Move        = Nslider(V)
# Type        = Node slider (uniform on possible positions)
# Parameter   = V [param. 6] (Branch lengths)
# Tuningparam = lambda (Multiplier tuning parameter)
# lambda = 0.191
# Rel. prob.  = 0.0
# MrBayes> propset nslider(V)$prob=5


# @model: could be in c("JC69", "F81", "HKY85", "K80", "GTR)
compileBatchScriptForMrBayes <- function(data.file, 
                                         bath.script.file.name, 
                                         model="JC69", 
                                         thinning = 10, 
                                         burn.in = 1000, 
                                         ngenerations = 10000, 
                                         fixed.topology = FALSE, 
                                         fixed.branch.length = FALSE) {
  mrbayes <- new("mrbayesbatch")
  setTags(mrbayes) <- list(command="set", list=list(autoclose="yes", nowarn="yes"))
  addMultiPartCommand(mrbayes) <- (list("execute", data.file))
  
  # setup the model
  switch(model, 
         JC69={ ##
           # equal exchangaiblity rates, no among-site rate variation
           setTags(mrbayes) <- list(command="lset", list=list(nst=1, rates="Equal"))
           
           # equal & fixed base (state) frequencies
           setTags(mrbayes) <- list(command="prset", list=list(statefreqpr="fixed(equal)"))
         }, 
         F81={ ##
           setTags(mrbayes) <- list(command="lset", list=list(nst=1, rates="Equal"))
         },
         HKY85={ ##
           setTags(mrbayes) <- list(command="lset", list=list(nst=2, rates="Equal"))
         },
         K80={ ##
           setTags(mrbayes) <- list(command="lset", list=list(nst=2, rates="Equal"))
           setTags(mrbayes) <- list(command="prset", list=list(statefreqpr="fixed(equal)"))
         },
         GTR={            
           setTags(mrbayes) <- list(command="lset", list=list(nst=6, rates="Equal"))
         },
{
  cat("Model ", model, " isn't supported yet.")
}
  )

# set the start values to be read from the provided initial tree
setTags(mrbayes) <- list(command="startvals", list=list(tau="mm", V="mm"))

# fixe topology?
if (fixed.topology == T) {
  # topology
  setTags(mrbayes) <- list(command="propset", list=list("ParsSPR(Tau,V)$prob"=0))
  setTags(mrbayes) <- list(command="propset", list=list("NNI(Tau,V)$prob"=0))
}

if (fixed.branch.length == T) {
  # branch_length
  setTags(mrbayes) <- list(command="propset", list=list("Nodeslider(V)$prob"=0))
  setTags(mrbayes) <- list(command="propset", list=list("Multiplier(V)$prob"=0))
  setTags(mrbayes) <- list(command="propset", list=list("TLMultiplier(V)$prob"=0))
}

setTags(mrbayes) <- list(command="mcmc", list=list(ngen=ngenerations, samplefreq=thinning))
setTags(mrbayes) <- list(command="sump", list=list(burnin=burn.in))
setTags(mrbayes) <- list(command="sumt", list=list(burnin=burn.in))
writeToDisk(mrbayes, bath.script.file.name)  
}

mrbayes.analysis <- function(batch.script) {
  #system(paste0("mb ", batch.script), ignore.stdout = F, )
  system(paste0("mbs ", batch.script), ignore.stdout = F)
}


mrbayes.calculate.ESS <- function(elapsed.time) {
  file.names <- list.files(".", pattern = "*.run1.p")

  print(file.names);
  
  mrbayes.posterior <- read.table(file.names[1], skip=1, header=T, check.names=F)
  
  mrbayes.posterior <- mrbayes.posterior[, c(3:ncol(mrbayes.posterior))]
  
  ESS <- effectiveSize(mrbayes.posterior)
  
  # standardize column names
  ESS <- mrbayes.standardizeColumns.ESS(ESS)
  
  ESSPS <- effectiveSize(mrbayes.posterior)/elapsed.time
  names(ESSPS) <- names(ESS)
  # save the values
  MRBAYES_ESS_PATH <<- file.path(MRBAYES_EXPERIMENT_PATH, "ess.txt")
  MRBAYES_ESSPERSECOND_PATH <<- file.path(MRBAYES_EXPERIMENT_PATH, "ess_per_second.txt") 
  write.table(data.frame(ESS, check.names=F), MRBAYES_ESS_PATH)
  write.table(data.frame(ESSPS, check.names=F), MRBAYES_ESSPERSECOND_PATH)
}

# http://hydrodictyon.eeb.uconn.edu/eebedia/index.php/Phylogenetics:_MrBayes_Lab
# map column names from mrbayes's output to those from conifer
mrbayes.standardizeColumns.ESS <- function(ESS) {
  
  n <- names(ESS)
  # 1. change from r(A<->C) to q(A(0),C(0))
  indexes <- grep("r\\([ACTG]<->[ACTG]\\)", n)
  if (length(indexes) > 0) {
    n[indexes] <- gsub("r\\(", "q\\(", n[indexes])
    n[indexes] <- gsub("<->", ",", n[indexes])
    n[indexes] <- gsub("([ACTG])", "\\1(0)", n[indexes], perl=T)
  }
  # 2. change from pi(A) to stationary(A(0))
  indexes <- grep("pi\\([ACTG]\\)", n)
  if (length(indexes) > 0) {
    n[indexes] <- gsub("pi", "stationary", n[indexes])
    n[indexes] <- gsub("([ACTG])", "\\1(0)", n[indexes], perl=T)
  }
  
  names(ESS) <- n
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
  consensus.tree <- consensus(all.trees, p=.5)
  
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
  MRBAYES_CONSENSUS_TREE_PATH <<- file.path(MRBAYES_EXPERIMENT_PATH, "consensus.tree")
  write.tree(consensus.tree, MRBAYES_CONSENSUS_TREE_PATH)
}


mrbayes.driver.function <- function(treeFilePath, 
                                    alignmentFilePath, 
                                    batch.dir, 
                                    model, 
                                    thinning, 
                                    burn.in, 
                                    numofgen, 
                                    fixed.topology,
                                    fixed.branch.length) {
  print(deparse(match.call()))
  currentWD <- getwd()
  
  if (is.null(treeFilePath)) treeFilePath <- random.tree.from.fasta(alignmentFilePath)
  
  # set dir where the mrbayes files should be placed
  MRBAYES_EXPERIMENT_PATH <<- file.path(batch.dir, format(Sys.time(), "mrbayes-%y-%m-%d-%H-%M-%S")) 
  dir.create(MRBAYES_EXPERIMENT_PATH)
  
  setwd(MRBAYES_EXPERIMENT_PATH)
  batch.file.name <- "mbbatch.txt"
  
  print(getwd())
  
  # make symbolik links to the data.files
  system(paste0("ln -s ", alignmentFilePath, " ", MRBAYES_EXPERIMENT_PATH))
  system(paste0("ln -s ", treeFilePath, " ", MRBAYES_EXPERIMENT_PATH))
  
  makeDataFileForMrBayes(alignmentFile=alignmentFilePath, treeFile=treeFilePath)
  compileBatchScriptForMrBayes(data.file="datafile.nex", 
                               bath.script.file.name=batch.file.name, 
                               model=model, 
                               thinning=thinning, 
                               burn.in=burn.in, 
                               ngenerations=numofgen,
                               fixed.topology=fixed.topology,
                               fixed.branch.length=fixed.branch.length)
  
  # get the elapsed time
  the.run.time <- system.time(mrbayes.analysis(batch.file.name))
  elapsed.time <- the.run.time[3]
  
  # calculate ESS and ESS per second
  mrbayes.calculate.ESS(elapsed.time)
  
  # calculate the concensus tree
  mrbayes.calculate.consensus.tree(burn.in, thinning)
  
  # keep record of the analysis time
  writeLines(c(paste("Analysis for ", batch.file.name, ", took", elapsed.time, "seconds.")), "experiment.details.txt")
  
  # restore to the current working directory
  setwd(currentWD)
  
  MRBAYES_EXPERIMENT_PATH
}


# /Users/sohrab/project/conifer_fork/src/main/resources/conifer/sampleInput/FES_4.fasta
# /Users/sohrab/project/conifer_fork/src/main/resources/conifer/sampleInput/FES.ape.4.nwk
# /Users/sohrab/project/conifer_fork/mrbayes

# mrbayes.driver.function(treeFilePath = "/Users/sohrab/project/conifer_fork/src/main/resources/conifer/sampleInput/FES.ape.4.nwk", alignmentFilePath = "/Users/sohrab/project/conifer_fork/src/main/resources/conifer/sampleInput/FES_4.fasta", model = "K80", batch.dir = "/Users/sohrab/project/conifer_fork/mrbayes", fixed.topology = F, fixed.branch.length = F, numofgen = 10000, burn.in = 1000, thinning = 10)

# mrbayes.driver.function(treeFilePath = NULL, alignmentFilePath = "/Users/sohrab/project/conifer/simulated.data/simulation.5_FIXED_TOPOLOGY_FIXED_BRANCH_LENGHTS_DNAGTR/SimulatedData.fasta", model = "GTR", batch.dir = "/Users/sohrab/project/conifer_fork/mrbayes", fixed.topology = T, fixed.branch.length = T, numofgen = 10000, burn.in = 1000, thinning = 10)

#mrbayes.driver.function("FES.ape.4.nwk", "FES_4.fasta")

#consensus <- read.nexus("/Users/sohrab/Me/Apply/Canada\ Apply/Courses/Third\ Semester/conifer/extras/mrbayes/FES_8_batch.GTR.two/FES_8_batch.GTR.two.nex.tree1.con.tre")

#system("mb")

# sample GTR batch for mrbayes (with data.file including a tree called mm)
# fixed topology
# Execute data.nex
# lset nst=6 rates=Equal
# startvals tau=mm V=mm;
# 
# showmoves
# 
# # branch_length
# propset Nodeslider(V)$prob=0
# 
# # topology
# propset ParsSPR(Tau,V)$prob=0
# 
# # topology
# propset NNI(Tau,V)$prob=0
# 
# # branch_length
# propset Multiplier(V)$prob=0
# 
# # branch_length
# propset TLMultiplier(V)$prob=0
# 
# mcmc ngen=10000 samplefreq=10
# sump burnin=1000
# sumt burnin=1000

# /Users/sohrab/project/conifer/simulated.data/simulation.4_DEFAULT_DNAGTR/SimulatedData.fasta
# mrbayes.driver.function(NULL, "/Users/sohrab/project/conifer/simulated.data/simulation.4_DEFAULT_DNAGTR/SimulatedData.fasta", ".", "GTR", 10, 1000, 10000, T, T)
