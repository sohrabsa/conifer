# conifer benchmarking
require(XML)

# SENSUS_TREE_PATH <- ""
CONIFER_ESS_PATH <- ""
CONIFER_ESSPERSECOND_PATH <- ""
CONIFER_PROJECT_DIR <- ""
CONIFER_EXPERIMENT_PATH <- ""


conifer.load.consensus.tree <- function() {
  read.tree(CONIFER_CONSENSUS_TREE_PATH)
}

conifer.load.ess <- function() {
  read.csv(CONIFER_ESS_PATH, header = T, row.names=1)
}

conifer.load.esspersecond <- function() {
  read.csv(CONIFER_ESSPERSECOND_PATH, header = T, row.names=1)
}


conifer.class.path.string <- function() {
  # find class.path 
  classPath <- file.path(CONIFER_PROJECT_DIR, ".classpath")
  if (!file.exists(classPath)) {
    # try building the project
    c <- getwd()
    setwd(CONIFER_PROJECT_DIR)
    system("gradle build")
    system("gradle eclipse")
    setwd(c)
  }
  
  data <- xmlParse(classPath)
  ldata <- xmlToList(data)
  paths <- lapply(ldata, function(x) x['path'])
  classpaths <- paste0(paths, collapse = ":")
  
  # TODO: fix paths
  gsub(file.path(CONIFER_PROJECT_DIR, "/build/libs/conifer.jar"), file.path(CONIFER_PROJECT_DIR, "/build/classes/main"), classpaths)
}

# compile and run conifer
conifer.run <- function(alignmentFilePath, treeFilePath, thinning, burn.in, numofgen, model, fixed.topology, fixed.branch.length) {
  classpaths <- conifer.class.path.string()
  
  t <- getwd()
  setwd(file.path(CONIFER_PROJECT_DIR, "src/main/java/conifer"))
  print(file.path(CONIFER_PROJECT_DIR, "src/main/java/conifer"))
  cat("Warning! Ignoring model for the time being! NOT IMPLEMENTED IN CONIFER SIDE!\n")
  
  # only compile the changed class for now)
  system(paste0("javac ", " -classpath ", classpaths, " ", "-source 1.7 TestPhyloModel.java"))
  system(paste0("mv ", file.path(CONIFER_PROJECT_DIR, "src/main/java/conifer/TestPhyloModel*.class"), " ",  file.path(CONIFER_PROJECT_DIR, "/build/classes/main/conifer/")))
    
  classpaths <- gsub(file.path(CONIFER_PROJECT_DIR, "/build/libs/conifer.jar:"), "", classpaths)
  commandString <- paste0("java", " -classpath ", classpaths, " conifer.TestPhyloModel",
                          " --initialTreeFilePath '", treeFilePath,
                          "' --alignmentFilePath '", alignmentFilePath, "'", 
                          " --nMCMCSweeps ", numofgen,
                          " --burnIn ", burn.in, 
                          " --thinningPeriod ", thinning, 
                          " --fixedTopology", fixed.topology,
                          " --fixedBranchLength", fixed.branch.length)
  commandString
  setwd(CONIFER_PROJECT_DIR)
  f <- system(commandString, intern = T)
  setwd(t)
  outputfolder <- gsub("outputFolder : ", "",  tail(f, n = 1))
  CONIFER_EXPERIMENT_PATH <<- outputfolder
  cat("Results may be accessed in ", outputfolder, "\n")
  outputfolder
}

conifer.calculate.ESS <- function(input.dir) {
  # calculate ess perseconds
  input.dir <- "/Users/sohrab/project/conifer/results/all/2014-08-19-09-54-13-j0Yos85y.exec"
  experiment.files <- c(input.dir)
  master.ess <- data.frame()
  //experiment <- experiment.files[[1]]
  for (experiment in experiment.files) {
    sub.folders <- list.files(experiment)
    sub.folders <- sub.folders[grepl("-csv", sub.folders)]
    
    //parameter.folder <- sub.folders[[1]]
    
    for (parameter.folder in sub.folders) {
      # for any csv file
      csv.files <- list.files(file.path(experiment, parameter.folder))
      csv.files <-csv.files[grepl(".csv", csv.files)]
      rowNamesList <- list()
      l <- data.frame()
      
      //csv.file <- csv.files[[1]]
      
      for (csv.file in csv.files) {
        if (csv.file == "ess.csv" || grep("BivariateIdentity", csv.file) == 1) next
        
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
  
  # combine different ESS in a master one
  CONIFER_ESS_PATH <<- file.path(input.dir, "master.ess.txt")
  CONIFER_ESSPERSECOND_PATH <<- file.path(input.dir, "master.ess.per.second.txt")
  write.csv(master.ess, CONIFER_ESS_PATH)
  master.ess[, 1] <- master.ess[, 1] / experiment.length.in.seconds
  #names(master.ess) <- gsub("ESS", "ESSPS", names(master.ess))
  write.csv(master.ess, CONIFER_ESSPERSECOND_PATH)
}


conifer.calculate.consensus.tree <- function(input.dir) {
  file.names <- list.files(input.dir, pattern = "*.newick", full.names = T)
  #print(file.names[1])
  all.trees <- read.tree(file.names[1])
  
  # no need to discard any trees, it's been already taken care of.
  
  # calculate the consensus tree
  consensus.tree <- consensus(all.trees, p=.5)
  
  # calculate the clade support for the consensus tree
  #all.sub.trees <- get.sub.trees(all.trees)
  #counts <- get.count.for.tree(consensus.tree, all.sub.trees)
  
  # write clade support to the disk
  #support.strings <- unlist(lapply(subtrees(consensus.tree), write.tree))
  #d <- data.frame(clades=support.strings, counts=counts)
  #write.csv(d, "counsensus.clade.counts.csv", row.names=F)
  
  # write the consensus tree to the disk
  # this will write clade support as [internal] node labels 
  #consensus.tree$node.label <- counts
  CONIFER_CONSENSUS_TREE_PATH <<- file.path(input.dir, "consensus.tree")
  write.tree(consensus.tree, CONIFER_CONSENSUS_TREE_PATH)
}

# sivhf300lfmvnGOGO

conifer.driver.function <- function(treeFilePath, 
                                    alignmentFilePath, 
                                    conifer.project.dir, 
                                    thinning, 
                                    burn.in, 
                                    numofgen, 
                                    model,
                                    fixed.topology,
                                    fixed.branch.length) {
  # calculate ESS and ESS per second
  CONIFER_PROJECT_DIR <<- conifer.project.dir
  
  # run conifer and return output directory
  output.folder <- conifer.run(treeFilePath, alignmentFilePath, thinning = thinning, burn.in=burn.in, numofgen=numofgen, model=model, fixed.topology=fixed.topology, fixed.branch.length=fixed.branch.length)

  # calculate ESS and ESS per second
  conifer.calculate.ESS(output.folder)
  
  # calculate the concensus tree
  conifer.calculate.consensus.tree(output.folder)
  
  CONIFER_EXPERIMENT_PATH
}
