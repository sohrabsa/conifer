# Benchmark scaffold
library(methods)
#source("/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/conifer_fork/confier_fork/r.scripts/MrBayesBatch.R")
mainDIR <- "/Users/sohrab/project/conifer_fork/r.scripts/"
mrbayes.possible.models <- function() c("JC69", "F81", "HKY85", "K80", "GTR")

# sample run
# Rscript driver.R -alignmentPath ~/conifer/src/main/resources/conifer/sampleInput/FES_4.fasta -initialTreePath ~/conifer/src/main/resources/conifer/sampleInput/FES.ape.4.nwk


dataDir <- "/home/sohrab/conifer/src/main/resources/conifer/sampleInput"

# sample data files
"/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES_4.fasta"
"/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES.ape.4.nwk"  

# compute and write the ess for each processable css file in the folder lists

#setwd("/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/conifer/results/all/")


# point of entry to the script
# given the inputs, runs all the analysis and conncets differnt parts.
# 1. parse the args for the input files
#   1.1. read from stdin
#   1.2. read from a config text file
# 2. run mrbayes and conifer with the given inputs
#   2.1. mrbayes
#     2.1.1. run mrbayes with the inputs
#     2.1.2. parse the outputs of mrbayes and produce ESS, ESSperSec, consensus tree with clade support, and clade support csv
#     2.1.3. create symlinks in the output folder
#   2.2 conifer
#     2.2.1. run conifer with the inputs
#     2.2.2. parse the outputs of conifer and produce ESS, ESSperSec, consensus tree with clade support, and clade support csv
#     2.2.3. create symlinks in the output folder
# 3. comparison
#     3.1. ESS
#       3.1.1. combined mrbayes and conifer ess data.frame 
#       3.1.2. barchart of mrbayes and conifer
#     3.2. ESSperSec
#     3.3. consensus tree
#       3.3.1. head to head graphs
#     3.4. clade support
#       3.4.1. table with clades with high posterior
#       3.4.2. table with common clade, probability in both mrbayes and conifer
#       3.4.3. over the same tips, head to head clades with highest probability
#

driver <- function() {
  # 1. parse the args for the input files
  #   1.1. read from stdin
  # read the input and output fasta full paths from stdin
  
  supportedArguments <- names(formals(runner))
  
  arg.desc <- list(model=paste0("Phylogenetic Model to use, could be one of the following: ",
                                paste0(mrbayes.possible.models(), collapse = ", "), 
                                ". Default value: ", formals(runner)$model),
                   numofgen=paste0("Number of MCMC sweeps. ",
                                   "Default value: ", formals(runner)$numofgen),
                   thinning=paste0("How often keep the samples. ",
                                   "Default value: ", formals(runner)$thinning),
                   burnin=paste0("Number of samples of initial samples to discard. ",
                                 "Default value: ", "numofgen * 0.1"),
                   initialTreePath=paste0("Path to the initla tree. Should be in Newick format.",
                                          "Default value: Creates a random tree if not set.", formals(runner)$initialTreePath),
                   alignmentPath=paste0("Path to the sequence alignment file in fasta format. ",
                                        "Default value: ", formals(runner)$alignmentPath),
                   coniferProjectDir=paste0("Root directory of your conifer repository. Will try to expand the given path. ",
                                            "Default value: ", eval(formals(runner)$coniferProjectDir)),
                   mrbayesOutputDir=paste0("Where to create subdirectories with mrbayes outputs. ", 
                                           "Default value: ", eval(formals(runner)$mrbayesOutputDir))
  )
  
  # read std input
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    stop("Please provide the path to the initial tree, and alignment file.")
  } else if (args[1] == "-h") {
    cat("Input arguments are as follows:\n")
    #cat(paste0(names(arg.desc), ":\t\t", arg.desc, collapse="\n "), "\n")
    formatting <- unlist(lapply(names(arg.desc), function(x) paste0(rep("+", 60 - nchar(x)), collapse="") ))
    t <- paste0(formatting, names(arg.desc), " : ", arg.desc, collapse="\n")
    #cat(paste0(formatting, names(arg.desc), " : ", arg.desc, collapse="\n"), "\n")
    t <- strsplit(t, "\n")[[1]]
    #cat(paste0(t, collapse="\n"), "\n")
    t <- gsub("\\+", " ", paste0(strwrap(t, exdent = 63, width=130), collapse="\n"))
    cat(t, "\n")
  } else {
    # parse command-line switches

    switches <- gsub("(-|--)", "", args[seq(1, length(args), 2)]) 
    
    if (!all(switches %in% supportedArguments)) stop("Unrecognized argument!")
    # TODO: further check the input values
    arg.indexes <- seq(2,length(args), 2)
    arg.list <- as.list(args[arg.indexes])
    names(arg.list) <- switches
    
    if (length(grep("alignmentPath", switches)) == 0) stop("Empty alignmentPath. Alignment file has to be provided.")
    if (!is.null(arg.list$model) && ! arg.list$model %in% mrbayes.possible.models()) stop(paste0("Don't support provided model ", arg.list$model), " yet.")
    if (!file.exists(arg.list$alignmentPath)) stop("Provided alignment file doesn't exists!")
    
    source(file.path(mainDIR, "MrBayesBatch.R"))
    source(file.path(mainDIR, "clader2.R"))
    source(file.path(mainDIR, "essGeneric.R"))
    source(file.path(mainDIR, "ConiferBenchmarking.R"))
    source(file.path(mainDIR, "MrBayesBenchmarking.R"))
    
    print(arg.list)
    do.call(runner, arg.list)
  }
}


runner <- function(model = "GTR", 
                   numofgen=10000, 
                   thinning=10, 
                   burnin=as.integer(numofgen*.1), 
                   initialTreePath=NULL, 
                   alignmentPath=NULL, 
                   coniferProjectDir=file.path("~/project", "conifer"),
                   mrbayesOutputDir=file.path("~", "mrbayes"), 
                   fixed.topology=FALSE,
                   fixed.branch.length=FALSE,
                   jumpmrbayes=FALSE, 
                   jumpconifer=FALSE) {

#  "/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES_4.fasta"
#  "/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES.ape.4.nwk"  
  
  coniferProjectDir <- path.expand(coniferProjectDir)
  print(coniferProjectDir)
  
  mrbayesOutputDir <- path.expand(mrbayesOutputDir)
  dir.create(mrbayesOutputDir, showWarnings=F)
  
  numofgen <- as.integer(numofgen)
  burnin <- as.integer(burnin)
  thinning <- as.integer(thinning)
  jumpmrbayes <- as.logical(as.integer(jumpmrbayes))
  jumpconifer <- as.logical(as.integer(jumpconifer))
  fixed.topology <- as.logical(as.integer(fixed.topology))
  fixed.branch.length <- as.logical(as.integer(fixed.branch.length))
  
  
  cat("numofgen=", numofgen, " burnin=", burnin)
  
  # TODO: validate inputs
  # handle a NULL tree
  if (is.null(initialTreePath)) {
    initialTreePath <- random.tree.from.fasta(alignmentPath)
  }
# Rscript driver.R -alignmentPath ~/conifer/src/main/resources/conifer/sampleInput/FES_4.fasta --numofgen 100 -jumpmarbayes 1
  

  print(initialTreePath)
  
  # 2. run mrbayes and conifer with the given inputs
  #   2.1. mrbayes
  #     2.1.1. run mrbayes with the inputs
  if (jumpmrbayes == FALSE) {  
    mrbayes.driver.function(alignmentFilePath=alignmentPath, 
                            treeFilePath=initialTreePath, 
                            batch.dir=mrbayesOutputDir, 
                            thinning=thinning, 
                            numofgen=numofgen, 
                            burn.in=burnin, 
                            model=model,
                            fixed.topology=fixed.topology,
                            fixed.branch.length=fixed.branch.length)
  }
  
  # TODO: check if mrbayes.driver.function worked
  
  #     2.1.2. parse the outputs of mrbayes and produce ESS, ESSperSec, consensus tree with clade support, and clade support csv
  #     2.1.3. create symlinks in the output folder
    
  #   2.2 conifer
  #     2.2.1. run conifer with the inputs
  
  if (!all(file.exists(coniferProjectDir, 
                       file.path(coniferProjectDir, "settings.gradle"), 
                       file.path(coniferProjectDir, "build.gradle")))) {
    stop(paste0("Couldn't find conifer at the provided path ", coniferProjectDir))
  }

  conifer.output.dir <- NULL
  
  if (jumpconifer == FALSE) {
    # TODO: fixed branch length
    conifer.output.dir <- conifer.driver.function(alignmentPath, 
                                                  initialTreePath, 
                                                  coniferProjectDir, 
                                                  thinning=thinning, 
                                                  burn.in=burnin, 
                                                  numofgen=numofgen, 
                                                  model=model,
                                                  fixed.topology=fixed.topology,
                                                  fixed.branch.length=fixed.branch.length)
  }
    #     2.2.2. parse the outputs of conifer and produce ESS, ESSperSec, consensus tree with clade support, and clade support csv
    #     2.2.3. create symlinks in the output folder
    
    # create a symbolic link in conifer's directory
    mrbayes.make.symlink(conifer.output.dir)

  if (!jumpconifer && !jumpmrbayes) {
    # 3. comparison
    #     3.1. ESS
    #       3.1.1. combined mrbayes and conifer ess data.frame
    
    mrbayes.ess <- mrbayes.load.ess()
    conifer.ess <- conifer.load.ess()
    ess <- combine.ess(mrbayes.ess, conifer.ess)
    #       3.1.2. barchart of mrbayes and conifer
    print("-=--")
    print(conifer.output.dir)
    make.ess.barchart(ess, numberOfRuns = numofgen, conifer.output.dir)
    
    #     3.2. ESSperSec
    mrbayes.ess.persecond <- mrbayes.load.esspersecond()
    conifer.ess.persecond <- conifer.load.esspersecond()
    essPerSecond <- combine.ess(mrbayes.ess.persecond, conifer.ess.persecond)
    make.ess.per.second.barchart(essPerSecond, numofgen, conifer.output.dir)
    
    
    #     3.3. consensus tree
    #       3.3.1. head to head graphs
    
    # load consensus trees
    conifer.consensus.tree <- conifer.load.consensus.tree()
    mrbayes.consensus.tree <- mrbayes.load.consensus.tree()
    
    # are the consensus trees the same?
    message <- ifelse(all.equal(conifer.consensus.tree, mrbayes.consensus.tree, use.edge.length=F), 
                      "Consensus trees are the same.", "Consensus trees do not match!") 
    print(message)
    writeLines(message, file.path(conifer.output.dir, "correctness.tests.txt"))
  }

  
    
#   if (length(conifer.consensus.tree$tip.label) > 20) {
#     mrbayes.consensus.tree <- subtrees(mrbayes.consensus.tree)[[which.max(mrbayes.consensus.tree$node.label)]]
#     conifer.consensus.tree <- subtrees(conifer.consensus.tree)[[which.max(conifer.consensus.tree$node.label)]] 
#   } 
#   
# 
#   plot.side.by.side(conifer.consensus.tree, mrbayes.consensus.tree, 
#                     c("Conifer", "MrBayes"), 
#                     file.path(conifer.output.dir, "consensus.sidebyside.jpg")) 
#   
  #     3.4. clade support
  #       3.4.1. table with clades with high posterior
  #       3.4.2. table with common clade, probability in both mrbayes and conifer
  #       3.4.3. over the same tips, head to head clades with highest probability
  #
  # CLARIFY
}
#
  
# TODO: use full data set
#/home/sohrab/conifer_fork/src/main/resources/conifer/sampleInput/UTY_full_trimmed.nwk
# a <- "/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES_4.fasta"
# t <- "/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES.ape.4.nwk"
#runner(a, t)
#conifer.output.dir <- conifer.driver.function(a, t, "/home/sohrab/conifer")


driver()


# Rscript Driver.R -alignmentPath /Users/sohrab/Desktop/extreme.fasta  -jumpconifer 1


#  Driver.R -alignmentPath /Users/sohrab/project/conifer/simulated.data/simulation.4_DEFAULT_DNAGTR/SimulatedData.fasta -jumpconifer 1
# Rscript Driver.R -alignmentPath /Users/sohrab/project/conifer/simulated.data/simulation.4_DEFAULT_DNAGTR/SimulatedData.fasta -initialTreePath /Users/sohrab/project/conifer/simulated.data/simulation.4_DEFAULT_DNAGTR/SimulatedDataTree.newick