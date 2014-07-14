# conifer benchmarking
require(XML)


conifer.class.path.string <- function(class.path.dir) {
  # find class.path 
  data <- xmlParse(file.path(class.path.dir, ".classpath"))
  ldata <- xmlToList(data)
  paths <- lapply(ldata, function(x) x['path'])
  classpaths <- paste0(paths, collapse = ":")
  
  # TODO: fix paths
  gsub("/home/sohrab/conifer/build/libs/conifer.jar", "/home/sohrab/conifer/build/classes/main", classpaths)
}

# make conifer
conifer.run <- function(conifer.project.dir, ) {

  classpath <- conifer.class.path.string(conifer.project.dir)
  
  # only compile the changed class for now)
  system(paste0("javac -classpath ", classpaths, " TestPhyloModel.java"))
  system(paste0("mv ", file.path(conifer.project.dir, "home/sohrab/conifer/src/main/java/conifer/TestPhyloModel*.class"), " ",  file.path(conifer.project.dir, "/build/classes/main/conifer/"))
  
  setwd("/home/sohrab/conifer/")
  
  # set the input values
  alignmentFilePath   <- file.path(conifer.project.dir, "/src/main/resources/conifer/sampleInput/FES_4.fasta")
  initialTreeFilePath <- file.path(conifer.project.dir, "/src/main/resources/conifer/sampleInput/FES.ape.4.nwk")
  
  classpaths <- gsub("/home/sohrab/conifer/build/libs/conifer.jar:", "", classpaths)
  commandString <- paste0("java -classpath ", classpaths, " conifer.TestPhyloModel", " --initialTreeFilePath '", initialTreeFilePath, "' --alignmentFilePath '", alignmentFilePath, "'")
  commandString
  f <- system(commandString, intern = T)
  outputfolder <- gsub("outputFolder : ", "",  tail(f, n = 1))
  outputfolder
}

conifer.calculate.ESS <- function(input.dir) {
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
      rownames(l) <- gsub(".csv", "", unlist(rowNamesList))
      k2 <- readLines(file.path(experiment, "experiment.details.txt"))
      k2 <- k2[length(k2)]
      experiment.length.in.seconds <- as.numeric(gsub("Total time in minutes: ", "", k2)) * 60 
      write.csv(l, file.path(experiment, parameter.folder, "ess.txt"))
      l[, 1] <- l[, 1] / experiment.length.in.seconds
      write.csv(l, file.path(experiment, parameter.folder, "ess_per_second.txt"))
    }
  } 
}


conifer.calculate.consensus.tree <- function(input.dir) {
  file.names <- list.files(input.dir, pattern = "*.newick")
  all.trees <- read.nexus(file.names[1])
  
  # no need to discard any trees, it's been already taken care of.
  
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


conifer.driver.function <- function(tree.file.name, alignment.file.name) {
  # calculate ESS and ESS per second
  batch.dir <- "/home/sohrab/conifer_fork/mrbayes"
  setwd(batch.dir)
  
  batch.file.name <- "july.compile.batch"
  alignment.file.name <- "FES_4.fasta"
  tree.file.name <- "FES.ape.4.nwk"
  
  # run conifer and return output directory
  output.folder <- conifer.run()
  
  # calculate ESS and ESS per second
  conifer.calculate.ESS(output.folder)
  
  # calculate the concensus tree
  conifer.calculate.consensus.tree(output.folder)  
}
