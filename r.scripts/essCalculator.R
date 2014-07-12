# Calculate ESS and ESS per second for csv files in a folder
require(tools)
if (!require(coda)) {
  # set a CRAN mirror
  local({r <- getOption("repos")
         r["CRAN"] <- "http://cran.us.r-project.org"
         options(repos=r)})
  install.packages("coda", dep=TRUE)
}
library(coda)
library(lattice)



"/home/sohrab/conifer/results/all/2014-07-12-11-36-46-GeyVzxr6.exec"


calculateESS <- function(input.dir) {
  experiment.files <- c(input.dir)
  
  for (experiment in experiment.files) {
    sub.folders <- list.files(experiment)
    sub.folders <- sub.folders[grepl("-csv", sub.folders)]
    
    for (parameter.folder in sub.folders) {
      # for any csv file
      csv.files <- list.files(file.path(experiment, parameter.folder))
      csv.files <-csv.files[grepl(".csv", csv.files)]
      rowNamesList <- list()
      l <- data.frame()
      for (csv.file in csv.files) {
        if (csv.file == "ess.csv") next
        
        the.file <- file.path(experiment, parameter.folder, csv.file)
        p <- read.table(the.file, row.names = 1, header=T, sep=",")
        l <- rbind(l, effectiveSize(p))
        rowNamesList[length(rowNamesList) + 1] <- csv.file
      }
      
      colnames(l) <- c("ESS")
      rownames(l) <- unlist(rowNamesList)
      write.csv(l, file.path(experiment, parameter.folder, "ess.txt"))
    }
  }
}

calculateESSperSecond <- function(input.dir) {
  # calculate ess perseconds
  experiment.files <- c(input.dir)
  for (experiment in experiment.files) {
    sub.folders <- list.files(experiment)
    sub.folders <- sub.folders[grepl("-csv", sub.folders)]
    
    for (parameter.folder in sub.folders) {
      # for any csv file
      csv.files <- list.files(file.path(experiment, parameter.folder))
      csv.files <-csv.files[grepl(".csv", csv.files)]
      rowNamesList <- list()
      l <- data.frame()
      for (csv.file in csv.files) {
        
        if (csv.file == "ess.csv") next
        
        the.file <- file.path(experiment, parameter.folder, csv.file)
        p <- read.table(the.file, row.names = 1, header=T, sep=",")
        l <- rbind(l, effectiveSize(p))
        rowNamesList[length(rowNamesList) + 1] <- csv.file      
      }
      
      colnames(l) <- c("ESS")
      rownames(l) <- unlist(rowNamesList)
      k2 <- readLines(file.path(experiment, "experiment.details.txt"))
      k2 <- k2[length(k2)]
      experiment.length.in.seconds <- as.numeric(gsub("Total time in minutes: ", "", k2)) * 60 
      l[, 1] <- l[, 1] / experiment.length.in.seconds
      write.csv(l, file.path(experiment, parameter.folder, "ess_per_second.txt"))
    }
  } 
}



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
#setwd("/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/conifer/results/all/")




