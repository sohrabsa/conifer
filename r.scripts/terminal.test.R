#source("clader2.R")
#source("MrBayesBenchmarking.R")
mrbayes.possible.models <- function() c("JC69", "F81", "HKY85", "K80", "GTR")
runner <- function(model = "GTR", 
                   numofgen=10000, 
                   thinning=10, 
                   burnin=as.integer(numofgen*.1), 
                   initialTreePath=NULL, 
                   alignmentPath=NULL, 
                   coniferProjectDir=file.path("~", system('whoami', intern=TRUE), "conifer")) {
  
  #print(as.list(match.call())[-c(1)] )
  
  #print(burnin)
  #print(coniferProjectDir)
}


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
    
  if (length(arches %in% supportedArguments)) stop("Unrecognized argument!")
  # TODO: further check the input values
  arg.list <- list(args[seq(2,length(args), 2)])
  names(arg.list) <- switches
  
  if (length(grep("alignmentPath", switches)) == 0) stop("Empty alignmentPath. Alignment file has to be provided.")
  
  do.call(runner, arg.list)
}





#do.call(bbc, list(b=1,a=2))



