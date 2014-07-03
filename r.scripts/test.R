

h <- "Hello Darling 

You're everything to me

"

p = "0.05253213485445601
0.17296522759112995
0.623745672969554
0.5115519496702345
0.7330154035058204
0.4254152910077026
0.15240305490905195
0.3691203294984145
0.5400585825889358
0.6553210155123839
0.630784618016299
0.6626519550956006
0.04448081243340539
0.9837535051809084
0.26367600467647234
0.5120081384824237
0.4897951278719058
0.61551437713703
0.6514010601503654
0.5623398713804768
0.4285071619661588
0.4932510169869517
0.32967576165308193
0.5046821233010397
0.4329560197971871
0.3357758363545899
0.889724875074086
0.3247479306742006
0.33817814958739173
0.36761197617261593
0.5972343295307629
0.31521652912873904
0.609013770671869
0.3886614640280912
0.4831810739934551
0.8390333428075964
0.7896206446018039
0.6766591284720659
0.5422894396573362
0.1425380141610554
0.5432194969932758
0.5587543353536408
0.5003556894910499
0.7537346448325046
0.6377727446878806
0.17328065258716044
0.3756340257106375
0.6393674409184215
0.6955578509087829
0.4235049718573491
0.457759659178936
0.43100850339427865
0.7395561319098098
0.7277310149619456
0.7767949359985122
0.8669535110557018
0.8850284003095333
0.3182472595721377
0.6509916239098651
0.33376783518448827
0.09317918621539878
0.2035022710442091
0.2958417916202911
0.23152281305879682
0.43046634238927456
0.14397415773277203
0.5176390993851431
0.8813718618522608
0.6576154213578678
0.37632112359214587
0.2851758627643598
0.7532364042788996
0.9747699610704141
0.6081448876276693
0.2692055838513698
0.3978350102044119
0.5306115414839787
0.1386514204573123
0.292773842745565
0.4317183401011175
0.7959803403233782
0.5210452397453798
0.8215131371007323
0.263639748022271
0.7367534536793166
0.3370755984290361
0.4584781494205735
0.4899780455384934
0.4412320604102836
0.2168207280810465
0.31398164104207804
0.6929825288313386
0.3078402147930963
0.31112356295594923
0.2340523469762367
0.5439460104251186
0.8526225076457143
0.295086272460221
0.5998547840150817
0.31707041289935983"


setwd("/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/bayonet")
#f <- read.table("testbeta.txt")
f <- read.table("testchisquared.txt")
hist(f[, 1], freq=F)
mean(f[, 1])
z <- as.numeric(format(f[, 1], digits = 1))

f1 <- data.frame(d=z, count=1)
f1.all <- aggregate(count~z, f1, FUN=sum)
f1$d[which.max(f1.all$count)]




library(coda)
#setwd("/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/bayonet/results/all/bayonet.distributions.TestLaplace-7BKMAiUa.exec/realization-csv")
setwd("/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/conifer/results/all/conifer.BlangSmallExample-ksPEmWww.exec/max-csv")
p <- read.table("CODAchain1.txt", row.names = 1)
(ess <- effectiveSize(p))
unlist(unname(ess["V2"]))/nrow(p)

library(LaplacesDemon)
ESS(p)


# a function, given start and end, 
# make the inbetween entries
# download them all
# write the result in fasta file

# DENND5A  HM759172:HM759362
# AXIN1	HM764222:HM764387
# FES	HM761675:HM761833
# TYR	HM757455:HM757629
# UTY	HM757170:HM757281
library(RNCBI)

setwd("/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/3rd Rotation/NIPS/")

retrieve.ncbi.sequences <- function(seqStart, seqEnd, outputFile) {
  
  # 0. initate all lines list
  all.lines = list()
  
  # 1. make the range
  startPos <- gsub("HM", "", seqStart)
  stopPos <- gsub("HM", "", seqEnd)
  
  (seqNames <- unlist(lapply(startPos:stopPos, function(x) { paste0("HM", x)  })))
  
  # 2. for each sequence, fetch the results
  for(seq in seqNames) {
    ncbi <- NCBI()
    rq <- EFetch(ncbi, "sequence")
    rq <- setRequestParameter(rq, "db", "nuccore")
    rq <- setRequestParameter(rq, "rettype", "fasta")
    rq <- setRequestParameter(rq, "id", c(seq))
    rq <- requestEFetch(rq)
    res <- getResults(rq)
    
    # add the header
    # sample header:
    #>gi|302454543|gb|HM759180.1| Cercopithecus cephus isolate CCE-1 DENND5A gene, partial sequence
    s <- res$tseqset$tseqset$tseqsetsequence$tseqsetsequence$tseq
    all.lines[length(all.lines) + 1] <- 
      paste0(">gi|", s$tseq_gi,
             "|gb|", s$tseq_accver,
             "| ", s$tseq_defline)
    # add the sequence
    all.lines[length(all.lines) + 1] <- s$tseq_sequence
    all.lines[length(all.lines) + 1] <- ""
  }

  # write to file
  fileConn<-file(outputFile)
  writeLines(unlist(all.lines), fileConn)
  close(fileConn)
}


# DENND5A  HM759172:HM759362
# AXIN1  HM764222:HM764387
# FES	HM761675:HM761833
# TYR	HM757455:HM757629
# UTY	HM757170:HM757281

geneSequences <- list("AXIN1"=c("HM764222", "HM764387"), 
                      "FES"=c("HM761675","HM761833"),
                      "TYR"=c("HM757455","HM757629"),
                      "UTY"=c("HM757170","HM757281"))                      

for (geneSequence in names(geneSequences)) {
  print(geneSequence)
  print(geneSequences[[geneSequence]][1])
  print(geneSequences[[geneSequence]][2])
  #print(geneSequence[1], geneSequence[2])
  retrieve.ncbi.sequences(geneSequences[[geneSequence]][1], geneSequences[[geneSequence]][2],
                          geneSequence)
}

retrieve.ncbi.sequences("HM759172", "HM759362")


# make phylogenetic tree with 112 leaves
library(ape)

setwd("/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/3rd Rotation/NIPS/data_files/sequences/the_four")

# UTY
gene.labels <- list("UYT"=c("Homo_sapiens",
                            "Pan_troglodytes_chimp", 
                            "Gorilla_gorilla", 
                            "Pongo_pygmaeus_orangutan"), 
                    "FES"=c("Homo_sapiens", 
                            "Pan_troglodytes_PTR-60_chimp", 
                            "Gorilla_gorilla", 
                            "Pongo_pygmaeus_PPY-155_orangutan"))
for (gene.label in names(gene.labels)) {
  print(gene.label) 
  print(gene.labels[[gene.label]])
  
  the.tree <- rtree(4, rooted=F, tip.label = gene.labels[[gene.label]])
  plot(the.tree)
  write.tree(the.tree, paste0(gene.label, ".ape.4.nwk"), digits=3)
}

FES8=c("Homo_sapiens", 
      "Pan_troglodytes_chimp", 
      "Gorilla_gorilla", 
      "Pongo_pygmaeus_orangutan",
      "Tupaia_glis_treeshrew", 
      "Cacajao_melanocephalus_uakari",
      "Loris_tardigradus_red_slender_loris",
      "Presbytis_melalophos_sumatran_surili"
)

the.tree <- rtree(8, rooted=F, tip.label = FES8)
plot(the.tree)
write.tree(the.tree, "FES.ape.8.nwk", digits=3)


for (gene.label in gene.labels) {
  
  
}



the.tree <- rtree(3, br=NULL)

# ref for nexus
# http://mrbayes.sourceforge.net/wiki/index.php/Tutorial_3.2
fastaToNexus <- function() {
  # use seqmagick
  # seqmagick convert --output-format nexus --alphabet dna UYT_4.fasta UYT_4.nex
}



readLines("hyper.prior.txt")
k <- readLines("hyper.prior.txt")
k3 <- k[1:length(k) %% 3 == 0]
k3 <- gsub("real", "", k3)
k3 = gsub("\\(", "", k3)
k3 = gsub("\\)", "", k3)
mean(as.numeric(k3))


library(PearsonDS)
e = exp(1)



test.plot <- function(the.alpha, the.beta) {
  plot.new()
  #beta = 2.8
  #alpha =2.333  
  alpha = the.alpha
  #beta = 2.358 
  beta = the.beta
  #beta = 2.58
  #alpha = e - .06
  #alpha = 2.58
  #scale = .79
  scale = 1
  the.seq <- seq(0.00,10,by=0.001)
  pVIpars <- list(a=alpha, b=beta, location=0, scale=scale)
  betaprime <- dpearsonVI(the.seq,params=pVIpars)
  lines(betaprime, x=the.seq, col="blue")
  lines(dlnorm(the.seq), x=the.seq, col="red")
}

# from limiting distribution
alpha = 1
beta = .5
scale = .5
the.seq <- seq(0.001,6,by=0.01)
pVIpars <- list(a=alpha, b=beta, location=0, scale=scale)
betaprime <- dpearsonVI(the.seq,params=pVIpars)
plot(betaprime, x=the.seq)

lines(dlnorm(the.seq), x=the.seq)

# http://robinlovelace.net/2013/10/23/nls-demonstation.html

x <- dlnorm(the.seq)

fm <- nls(y ~ cbind(1, ), start = c(B = 1), alg = "plinear"); fm

ss <- function(aalpha, bbeta) {
  dpearsonVI(the.seq, a= aalpha, b=bbeta, location=0, scale=1)
}
x <- the.seq
y <- dlnorm(the.seq)
df <- data.frame(x, y)

# fit a model
fit <- nls(y ~ ss(aalpha, bbeta), data = df, start = list(aalpha = 2.3, bbeta=2.3), algorithm="default")
test.plot((coef(summary(fit)))['aalpha', 1], (coef(summary(fit)))['bbeta', 1])
summary(fit)

# seqmagick mogrify --squeeze-threshold .1 FES_8.g.fasta


# compute and write the ess for each processable css file in the folder lists
setwd("/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/conifer/results/all/")

experiment.files <- list("conifer.SimplePhyloModel-Yp5ke8X9.exec", 
                         "experiment.conifer.SimplePhyloModel-Yp5ke8X9.exec.1.1401997544310", 
                         "experiment.conifer.SimplePhyloModel-Yp5ke8X9.exec.2.1401997893576",
                         "experiment.conifer.SimplePhyloModel-Yp5ke8X9.exec.3.1401998135461", 
                         "experiment.conifer.SimplePhyloModel-Yp5ke8X9.exec.4.1401998402376", 
                         "experiment.conifer.SimplePhyloModel-Yp5ke8X9.exec.5.1401998655232", 
                         "experiment.conifer.SimplePhyloModel-Yp5ke8X9.exec.6.1401999028796", 
                         "experiment.conifer.SimplePhyloModel-Yp5ke8X9.exec.7.1401999241764", 
                         "experiment.conifer.SimplePhyloModel-Yp5ke8X9.exec.8.1401999650578",
                         "conifer.TestPhyloModel-0hKG30ua.exec", "conifer.TestPhylo2-RMkQi7vJ.exec",
                         "conifer.TestPhyloExtraTaxa-NNA0Y5In.exec",
                         "conifer.HierarchicalPhyloModel-t6N3LgKL.exec")

experiment.files <- c("conifer.TestPhyloModel-4Q5hg4Sh.exec")
  
for (experiment in experiment.files) {
  sub.folders <- list.files(experiment)
  sub.folders <- sub.folders[grepl("-csv", sub.folders)]
  
  
  #sub.folders <- sub.folders[2]
  
  for (parameter.folder in sub.folders) {
    # for any csv file
    csv.files <- list.files(file.path(experiment, parameter.folder))
    #print(csv.files)
    csv.files <-csv.files[grepl(".csv", csv.files)]
    #the.file <- file.path(experiment, parameter.folder, "CODAchain1.txt");
    rowNamesList <- list()
    l <- data.frame()
    for (csv.file in csv.files) {
      # print(the.file)
      #print(readLines(the.file))
      
      if (csv.file == "ess.csv") next
      
      the.file <- file.path(experiment, parameter.folder, csv.file)
      print("the.file is")
      print(the.file)
      p <- read.table(the.file, row.names = 1, header=T, sep=",")
      print(head(p))
      l <- rbind(l, effectiveSize(p))
      print(l)
      #print(colnames(p))
      # rowNamesList[length(rowNamesList) + 1] <- colnames(p)
      rowNamesList[length(rowNamesList) + 1] <- csv.file
      
      #print((ess <- effectiveSize(p)))
      #print(unlist(unname(ess["V2"]))/nrow(p))
    }
    
    print("rowNamesList is")
    print(rowNamesList)
    print("l is")
    print(l)
    
    colnames(l) <- c("ESS")
    rownames(l) <- unlist(rowNamesList)
    write.csv(l, file.path(experiment, parameter.folder, "ess.txt"))
   
  }
}

library(coda)

# calculate ess perseconds
for (experiment in experiment.files) {
  sub.folders <- list.files(experiment)
  sub.folders <- sub.folders[grepl("-csv", sub.folders)]
  
  
  #sub.folders <- sub.folders[2]
  
  for (parameter.folder in sub.folders) {
    # for any csv file
    csv.files <- list.files(file.path(experiment, parameter.folder))
    #print(csv.files)
    csv.files <-csv.files[grepl(".csv", csv.files)]
    #the.file <- file.path(experiment, parameter.folder, "CODAchain1.txt");
    rowNamesList <- list()
    l <- data.frame()
    for (csv.file in csv.files) {
      # print(the.file)
      #print(readLines(the.file))
      
      if (csv.file == "ess.csv") next
      
      the.file <- file.path(experiment, parameter.folder, csv.file)
      print("the.file is")
      print(the.file)
      p <- read.table(the.file, row.names = 1, header=T, sep=",")
      print(head(p))
      l <- rbind(l, effectiveSize(p))
      print(l)
      #print(colnames(p))
      # rowNamesList[length(rowNamesList) + 1] <- colnames(p)
      rowNamesList[length(rowNamesList) + 1] <- csv.file      
    }
    
    print("rowNamesList is")
    print(rowNamesList)
    print("l is")
    print(l)
    
    colnames(l) <- c("ESS")
    rownames(l) <- unlist(rowNamesList)
    k2 <- readLines(file.path(experiment, "experiment.details.txt"))
    k2 <- k2[length(k2)]
    experiment.length.in.seconds <- as.numeric(gsub("Total time in minutes: ", "", k2)) * 60 
    l[, 1] <- l[, 1] / experiment.length.in.seconds
    write.csv(l, file.path(experiment, parameter.folder, "ess_per_second.txt"))
  }
}







list.files("experiment.conifer.SimplePhyloModel-Yp5ke8X9.exec.8.1401999650578/parameters-csv")

grepl("-csv", list.files("experiment.conifer.SimplePhyloModel-Yp5ke8X9.exec.8.1401999650578"))
# if file name has csv

p <- read.table("CODAchain1.txt", row.names = 1)
(ess <- effectiveSize(p))
unlist(unname(ess["V2"]))/nrow(p)



ff <- "/Users/sohrab/Me/Apply/Canada\ Apply/Courses/Third\ Semester/conifer/results/all/conifer.SimplePhyloModel-Yp5ke8X9.exec/parameters-csv/ess.csv"

k2 <- "Total time in minutes: 6.79305"

as.numeric(gsub("Total time in minutes: ", "", k2)) * 60 



# read one q matrix from conifer
# calculate l2 norm
# calculate spit it
# spit
# copy this vector
# shift to right
# spit the difference
# plot?

//f <- "conifer.TestPhyloModel-0hKG30ua.exec/parameters-csv"
setwd("/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/conifer/results/all/")
for (experiment in experiment.files) {
  f <- file.path(experiment, "parameters-csv")
  print(f)
  p <- list.files(f)
  p <- p[grepl("q\\(", p)]
  
  ll <- read.csv(file.path(f, p[1]), row.names = 1, header = T)
  for (q.file in p[-c(1)]) {
    ff <- read.csv(file.path(f, q.file), row.names = 1, header = T)
    #print(head(ff))
    ll <- cbind(ll, ff)
  }
  
  d <- as.matrix(ll)
  dd <- dist(d)
  length(dd)
  m <- as.matrix(dd)
  
  # calculate euclidean distance
  dist.list <- list()
  for (i in 2:nrow(ll)) {
    dist.list[length(dist.list) + 1] <- (m[i, i-1])
  }
  
  png(file.path(experiment, "parameters-csv", "q.euclidean.dist.png"))
  plot(unlist(dist.list), main="Euclidean distance between successive iterations for rate matrix Q", 
       xlab="MCMC iteration (thinned)", 
       ylab="Euclidean Distance")
  dev.off()
  
  
  # calculate the l2-norm-distance
  dist.list2 <- list()
  for (i in 2:nrow(ll)) {
    dist.list2[length(dist.list2) + 1] <- abs(sqrt(sum(ll[i, ]^2)) - sqrt(sum(ll[i-1, ]^2)))
  }
  
  png(file.path(experiment, "parameters-csv", "q.norm.dist.png"))
  plot(unlist(dist.list2), main="Distance between l2-norms of successive iterations for rate matrix Q", 
       xlab="MCMC iteration (thinned)", 
       ylab="L2 norm Distance")
  dev.off()
  
  # combine the two
  q.dist <- data.frame ("euclidean.dist" = unlist(dist.list), "norm.diff"= unlist(dist.list2))
  
  write.csv(q.dist, file.path(experiment, "parameters-csv", "q.dist.txt"))
  
  colnames(q.dist)
}


## for mrbayes

#setwd("/Users/sohrab/Me/Apply/Canada\ Apply/Courses/Third\ Semester/conifer/extras/mrbayes/FES_4_batch.GTR.simple")
setwd("/Users/sohrab/Me/Apply/Canada\ Apply/Courses/Third\ Semester/conifer/extras/mrbayes/FES_8_batch.GTR.two")
#ff  <- "FES_4_batch.GTR.simple.nex.run1.p"
ff  <- "FES_8_batch.GTR.two.nex.run1.p"
kk <- read.table(ff, skip=1, header=T)
head(kk)
kk <- kk[, c(4:9)]
head(kk)


d <- as.matrix(kk)
dd <- dist(d)
length(dd)
m <- as.matrix(dd)

# calculate euclidean distance
dist.list <- list()
for (i in 2:nrow(ll)) {
  dist.list[length(dist.list) + 1] <- (m[i, i-1])
}

png(file.path("q.euclidean.dist.png"))
plot(unlist(dist.list), main="Euclidean distance between successive iterations for rate matrix Q", 
     xlab="MCMC iteration (thinned)", 
     ylab="Euclidean Distance")
dev.off()

dist.list2 <- list()
for (i in 2:nrow(ll)) {
  dist.list2[length(dist.list2) + 1] <- abs(sqrt(sum(ll[i, ]^2)) - sqrt(sum(ll[i-1, ]^2)))
}

png(file.path("q.norm.dist.png"))
plot(unlist(dist.list2), main="Distance between l2-norms of successive iterations for rate matrix Q", 
     xlab="MCMC iteration (thinned)", 
     ylab="L2 norm Distance")
dev.off()

q.dist <- data.frame ("euclidean.dist" = unlist(dist.list), "norm.diff"= unlist(dist.list2))
write.csv(q.dist, file.path("q.dist.txt"))


### make all the fig.trees
setwd("/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/conifer/results/all/")
# add fig tree blocks to all fig tree files
t <- "/Users/sohrab/Me/Apply/Canada\ Apply/Courses/Third\ Semester/conifer/extras/figtree/lib/figtree.block.nexus"
for (exp in experiment.files) {
  f <- file.path(exp, "FESConsensusTree.Nexus")
  r <- file.path(exp, "FIG.FESConsensusTree.Nexus")
  (command.string <- paste0("cat \"", f, "\"",  " \"", t, "\" > \"", r, "\""))
  system(command.string)
  
  # for fig.tree augmented files, create the pdf version
  command.string <- "java -jar \"/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/conifer/extras/figtree/lib/figtree.jar\" -graphic PDF "
  command.string <- paste0(command.string, r, " ", file.path(exp, "FES.consensus.tree.pdf"))
  print(command.string)
  system(command.string)
}

# for uty consensus tree
for (exp in experiment.files) {
  f <- file.path(exp, "UTYConsensusTree.Nexus")
  r <- file.path(exp, "FIG.UTYConsensusTree.Nexus")
  (command.string <- paste0("cat \"", f, "\"",  " \"", t, "\" > \"", r, "\""))
  system(command.string)
  
  # for fig.tree augmented files, create the pdf version
  command.string <- "java -jar \"/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/conifer/extras/figtree/lib/figtree.jar\" -graphic PDF "
  command.string <- paste0(command.string, r, " ", file.path(exp, "UTY.consensus.tree.pdf"))
  print(command.string)
  system(command.string)
}
  
# for the mrbayes ones


  
system("cat ")



experiment.files


### exponential simulator
x = seq(-5,5, .01)
plot(x =x , y=exp(x))

# the 
J <- function(x) (x-.5)^2
plot(J, main="J")

T <- function(t) t
plot(t, main="t")

t = 1
i = .1
j = .2

P <- function(i, j, t) (exp(-(J(j) - J(i)) / T(t)))
plot(P())

which.min(rnorm(1:10))

######### TKF91
s <- "ACTGACTGAA"
chars <- c("A", "C", "G", "T")
n <- nchar(s)
r.ins <- 1
r.del <- 1.1
r.sub <- 1.2
types <- c("ins", "del", "sub")

s.ins <- rexp(1:(n+1), rate = r.ins)
s.del <- rexp(1:n, rate = r.del)
s.sub <- rexp(1:n, rate = r.sub)
deltaT <- min(c(s.ins, s.del, s.sub))
index <- which.min(c(s.ins, s.del, s.sub))
type <- types[ceiling(index / n)]

i <- index %% n

if (type == "ins") {
  char <- chars[as.vector(rmultinom(1, 1, c(.25, .25, .25, .25)))]
  s <- ins.char(s, char, i - 1)
} else if (type == "del") {
  s <- del.char(s, i)
} else {
  char <- chars[as.vector(rmultinom(1, 1, c(.25, .25, .25, .25)))]
  s <- sub.char(s, char, i)
}

del.char.at <- function(some.string, index)  {
  head.str <- substring(some.string, 1, index - 1)
  tail.str <- substring(some.string, index + 1)
  paste0(head.str, tail.str)
}

del.last.char <- function(some.string)  {
  index <- nchar(some.string)
  return(del.char.at(some.string, index))
}

list.str <- function(some.string) {
  as.list(strsplit(some.string, "")[[1]])
}



del.char <- function(some.string, index) {
  s <- list.str(some.string)
  s[index] <- NULL
  paste0(s, collapse="")
}

ins.char <- function(some.string, char, index) {
  s <- list.str(some.string)
  s <- append(s, char, after=index)
  paste0(s, collapse="")
}

sub.char <- function(some.string, char, index) {
  s <- list.str(some.string)
  s[[index]] <- char
  paste0(s, collapse="")
}

# paste0(append(list(99,66,88,11), 4, after=0), collapse=",")




require(tools)
if (!require(seqinr)) {
  install.packages("seqinr", dep=TRUE)
}

# read the input and output fasta full paths from stdin
args <- commandArgs(trailingOnly = TRUE)

input.path <- read.fasta(args[1])
output.path <- ifelse(!is.na(args[2]), 
                      args[2], 
                      paste0(file_path_sans_ext(args[1]), 
                                         "_purified.", file_ext(args[1])))
# read the fasta file
x <- read.fasta(input.path)

# replace anything chars that are not in the standard alphabet, namely (a, c, g, t) with a dash
x.replaced <- gsub("[^actg-]", "-", tolower(unlist(getSequence(x, as.string=T))), fixed = F)

# write the changed fasta file
write.fasta(x.replaced, names=names(x),file=output.path)


# given a character, if not the reference alphabet, will replace it with a dash
replace <- function(x) {
  if (!tolower(x) %in% c("a", "c", "g", "t")) {
    x = "-"
  } 
  x
}

# foreach species, check 
new=lapply(seq(length(x)), function(i) {
  unlist((lapply(getSequence(x)[[i]], function(c) replace(c))))  
})




# Pseudo code for conifer/Mr.Bayes comparison
# For each model on the same initial values
#   run the test on Mr. Bayes
#   run the test on Conifer
#   report ESS and final trees for both in a table
#   clade posterior probabilities)



# the model GTR




### --- a chunk of code to append a tree block at the end of the alignment file
## output the resulting trees in pdf
# 1. 



# run from command line
# conifer , input tree, input alignment file
# input parameters, and go, read parameters from an xml file
# first run mrbayes
# then run conifer
# probably read the parameters in JSON formats


# compute maximum clade tree
# and ESS for the branch length and rate matrix parameters
  









