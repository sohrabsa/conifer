# Benchmark scaffold
library(methods)
#source("/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/conifer_fork/confier_fork/r.scripts/MrBayesBatch.R")
mainDIR <- "/home/sohrab/conifer_fork/r.scripts"
source(file.path(mainDIR, "MrBayesBatch.R"))
source(file.path(mainDIR, "clader2.R"))
source(file.path(mainDIR, "essGeneric.R"))
source(file.path(mainDIR, "ConiferBenchmarking.R"))
source(file.path(mainDIR, "MrBayesBenchmarking.R"))

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
                                            "Default value: ", eval(formals(runner)$coniferProjectDir))
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

    switches <- gr(ep)("") args[seq(1, length(args), 2)]
    
    print(switches)
    if (!all(switches %in% supportedArguments)) stop("Unrecognized argument!")
    # TODO: further check the input values
    arg.list <- list(args[seq(2,length(args), 2)])
    names(arg.list) <- switches
    
    print(arg.list$model)
    if (length(grep("alignmentPath", switches)) == 0) stop("Empty alignmentPath. Alignment file has to be provided.")
    if (! arg.list$model %in% mrbayes.possible.models()) stop(paste0("Don't support provided model ", arg.list$model), "yet.")
    
    
    #do.call(runner, arg.list)
  }
}


runner <- function(model = "GTR", 
                   numofgen=10000, 
                   thinning=10, 
                   burnin=as.integer(numofgen*.1), 
                   initialTreePath=NULL, 
                   alignmentPath=NULL, 
                   coniferProjectDir=file.path("~", system('whoami', intern=TRUE), "conifer")) {

#  "/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES_4.fasta"
#  "/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES.ape.4.nwk"  
  
  
  # TODO: validate inputs
  # handle a NULL tree
  if (is.null(initialTreePath)) {
    initialTreePath <- random.tree.from.fasta(alignmentPath)
  }
  
  
  # 2. run mrbayes and conifer with the given inputs
  #   2.1. mrbayes
  #     2.1.1. run mrbayes with the inputs
  mrbayes.driver.function(alignmentFilePath = alignment.path, treeFilePath = initial.tree.path, "/home/sohrab/conifer_fork/mrbayes")
  
  # TODO: check if mrbayes.driver.function worked
  
  #     2.1.2. parse the outputs of mrbayes and produce ESS, ESSperSec, consensus tree with clade support, and clade support csv
  #     2.1.3. create symlinks in the output folder
    
  #   2.2 conifer
  #     2.2.1. run conifer with the inputs
  conifer.output.dir <- conifer.driver.function(alignment.path, initial.tree.path, "/home/sohrab/conifer")
  #     2.2.2. parse the outputs of conifer and produce ESS, ESSperSec, consensus tree with clade support, and clade support csv
  #     2.2.3. create symlinks in the output folder
  
  # create a symbolic link in conifer's directory
  mrbayes.make.symlink(conifer.output.dir)
  
  
  # 3. comparison
  #     3.1. ESS
  #       3.1.1. combined mrbayes and conifer ess data.frame
  mrbayes.ess <- mrbayes.load.ess()
  conifer.ess <- conifer.load.ess()
  ess <- combine.ess(mrbayes.ess, conifer.ess)
  #       3.1.2. barchart of mrbayes and conifer
  make.ess.barchart(ess, numberOfRuns = 1, conifer.output.dir)
  
  

  
  #     3.2. ESSperSec
  mrbayes.ess.persecond <- mrbayes.load.esspersecond()
  conifer.ess.persecond <- conifer.load.esspersecond()
  essPerSecond <- combine.ess(mrbayes.ess.persecond, conifer.ess.persecond)
  make.ess.per.second.barchart(essPerSecond, 1, conifer.output.dir)
  
  
  
  #     3.3. consensus tree
  #       3.3.1. head to head graphs
  
  # load consensus trees
  conifer.consensus.tree <- conifer.load.consensus.tree()
  mrbayes.consensus.tree <- mrbayes.load.consensus.tree()
  plot.side.by.side(conifer.consensus.tree, mrbayes.consensus.tree, c("Conifer", "MrBayes"), file.path(conifer.output.dir, "consensus.sidebyside.jpg"))
  
  
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
a <- "/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES_4.fasta"
t <- "/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES.ape.4.nwk"
  
#runner(a, t)

#conifer.output.dir <- conifer.driver.function(a, t, "/home/sohrab/conifer")


driver()