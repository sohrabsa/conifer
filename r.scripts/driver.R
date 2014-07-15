# Benchmark scaffold

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
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    stop("Please provide the path to the initial tree, and alignment file.")
  } else if (args[1] == "-h") {
    print("")
  } else {
    initial.tree.path <- args[1]
    alignment.path <- args[2]
  
    runner(alignment.path, initial.tree.path)
  }
}


runner <- function(alignment.path, initial.tree.path) {

#  "/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES_4.fasta"
#  "/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES.ape.4.nwk"  
  
  
  # TODO: validate inputs
  
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
  plot.side.by.side(conifer.consensus.tree, mrbayes.consensus.tree, c("Conifer", "MrBayes"), file.path(conifer.output.dir, "consensus.sidebyside.jpg") )
  
  
  #     3.4. clade support
  #       3.4.1. table with clades with high posterior
  #       3.4.2. table with common clade, probability in both mrbayes and conifer
  #       3.4.3. over the same tips, head to head clades with highest probability
  #
  
  # CLARIFY
}


#
a <- "/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES_4.fasta"
t <- "/home/sohrab/conifer/src/main/resources/conifer/sampleInput/FES.ape.4.nwk"
  
runner(a, t)

#conifer.output.dir <- conifer.driver.function(a, t, "/home/sohrab/conifer")


#driver()