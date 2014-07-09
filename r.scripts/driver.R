# Benchmark scaffold

source("/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/conifer_fork/confier_fork/r.scripts/MrBayesBatch.R")



# read the input and output fasta full paths from stdin
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Please provide the path to the experiment's directory.")
} else if (args[1] == "-h") {
  print("For every sub-dir with '-csv' in their name, will calculate ESS for every file with .csv extension and gather the results in an ess.txt file.")
} else {
  input.path <- args[1]
  
  calculateESS(input.path)
  calculateESSperSecond(input.path)
}

# compute and write the ess for each processable css file in the folder lists
setwd("/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/conifer/results/all/")



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
  
  # TODO: validate inputs
  
  # 2. run mrbayes and conifer with the given inputs
  #   2.1. mrbayes
  #     2.1.1. run mrbayes with the inputs
  mrbayes.driver.function(alignment.path, initial.tree.path)
  
  # TODO: check if mrbayes.driver.function worked
  #     2.1.2. parse the outputs of mrbayes and produce ESS, ESSperSec, consensus tree with clade support, and clade support csv
  #     2.1.3. create symlinks in the output folder
  
    
  }
  
  
  
  
  
}





driver()