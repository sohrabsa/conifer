# conifer benchmarking
require(XML)

# SENSUS_TREE_PATH <- ""
CONIFER_ESS_PATH <- ""
CONIFER_ESSPERSECOND_PATH <- ""
CONIFER_PROJECT_DIR <- ""

conifer.load.concensus.tree <- function() {
  read.tree(CONIFER_CONSENSUS_TREE_PATH)
}

conifer.load.ess <- function() {
  read.table(CONIFER_ESS_PATH)
}

conifer.load.esspersecond <- function() {
  read.table(CONIFER_ESSPERSECOND_PATH)
}


conifer.class.path.string <- function() {
  # find class.path 
  data <- xmlParse(file.path(CONIFER_PROJECT_DIR, ".classpath"))
  ldata <- xmlToList(data)
  paths <- lapply(ldata, function(x) x['path'])
  classpaths <- paste0(paths, collapse = ":")
  
  # TODO: fix paths
  gsub(file.path(CONIFER_PROJECT_DIR, "/build/libs/conifer.jar"), file.path(CONIFER_PROJECT_DIR, "/build/classes/main"), classpaths)
}

# compile and run conifer
conifer.run <- function(alignmentFilePath, treeFilePath) {
  classpath <- conifer.class.path.string()
  
  # only compile the changed class for now)
  system(paste0("javac -classpath ", classpaths, " ", file.path(CONIFER_PROJECT_DIR, "src/main/java/conifer/TestPhyloModel.java")))
  system(paste0("mv ", file.path(conifer.project.dir, "src/main/java/conifer/TestPhyloModel*.class"), " ",  file.path(conifer.project.dir, "/build/classes/main/conifer/"))
    
  classpaths <- gsub(file.path(CONIFER_PROJECT_DIR, "/build/libs/conifer.jar:"), "", classpaths)
  commandString <- paste0("java -classpath ", classpaths, " conifer.TestPhyloModel", " --initialTreeFilePath '", treeFilePath, "' --alignmentFilePath '", alignmentFilePath, "'")
  commandString
  f <- system(commandString, intern = T)
  outputfolder <- gsub("outputFolder : ", "",  tail(f, n = 1))
  outputfolder
}

conifer.calculate.ESS <- function(input.dir) {
  # calculate ess perseconds
  experiment.files <- c(input.dir)
  master.ess < - data.frame()
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
      master.ess <- rbind(master.ess, l)
      l[, 1] <- l[, 1] / experiment.length.in.seconds
      write.csv(l, file.path(experiment, parameter.folder, "ess_per_second.txt"))
    }
  } 
  
  #TODO
  CONIFER_ESS_PATH <<- file.path(input.dir, "master.ess.txt")
  CONIFER_ESSPERSECOND_PATH <<- file.path(input.dir, "master.ess.per.second.txt")
  write.csv(master.ess, CONIFER_ESS_PATH)
  master.ess[, 1] <- master.ess[, 1] / experiment.length.in.seconds
  write.csv(master.ess, CONIFER_ESSPERSECOND_PATH)
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
  CONIFER_CONSENSUS_TREE_PATH <<- file.path(input.dir, "consensus.tree")
  write.tree(consensus.tree, CONIFER_CONSENSUS_TREE_PATH)
}

# sivhf300lfmvnGOGO

conifer.driver.function <- function(treeFilePath, alignmentFilePath, conifer.project.dir) {
  # calculate ESS and ESS per second
  CONIFER_PROJECT_DIR <<- conifer.project.dir
  
  # run conifer and return output directory
  output.folder <- conifer.run(treeFilePath, alignmentFilePath)
  
  # calculate ESS and ESS per second
  conifer.calculate.ESS(output.folder)
  
  # calculate the concensus tree
  conifer.calculate.consensus.tree(output.folder)  
}
