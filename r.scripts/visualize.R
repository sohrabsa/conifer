## S series figure
library(ggplot2)
install.packages("ggplot2")
library(reshape2)
library(coda)
library(ggmcmc)
require(gridExtra)

setwd("/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/conifer/results/all/")
experimentDir <- "conifer.TestPhyloModel-0hKG30ua.exec"

setwd(experimentDir)

dat <- read.table("../../../data/s1/data.csv", sep = ",", colClasses="numeric")
dat <- data.frame(dat)
colnames(dat) <- c("x", "y")
ds <- read.csv("predictive.csv")


## Traceplots
#alpha0 <- read.coda("alpha0-csv/CODAchain1.txt", "alpha0-csv/CODAindex.txt")
#discount <- read.coda("discount-csv/CODAchain1.txt", "discount-csv/CODAindex.txt")
#kappa <- read.coda("kappa-csv/CODAchain1.txt", "kappa-csv/CODAindex.txt")
#nu <- read.coda("nu-csv/CODAchain1.txt", "nu-csv/CODAindex.txt")
rate <- read.coda("rate-csv/CODAchain1.txt", "rate-csv/CODAindex.txt")

#ch <- mcmc(cbind(alpha0, discount, kappa, nu))
ch <- mcmc(cbind(rate))

ch <- ggs(ch)

p <- ggs_traceplot(ch) + theme_minimal() + ggtitle("a")
qq <- ggs_density(ch, rug = T) + theme_minimal() + ggtitle("b")

#traces <- grid.arrange(q, arrangeGrob(p, qq, ncol = 1, nrow =1),  ncol = 1)
traces <- grid.arrange(p, p)


# have traces for the complex model
experiment.dir <- "conifer.TestPhyloModel-0hKG30ua.exec"
setwd(experiment.dir)

# the two trees
rate.FES <- read.coda("treePriorFES.rate-csv/CODAchain1.txt", "treePriorFES.rate-csv/CODAindex.txt")
rate.UTY <- read.coda("treePriorUTY.rate-csv/CODAchain1.txt", "treePriorUTY.rate-csv/CODAindex.txt")

#ch <- mcmc(cbind(alpha0, discount, kappa, nu))
ch <- mcmc(cbind(rate.FES, rate.UTY))

ch <- ggs(ch)
p <- ggs_traceplot(ch) + theme_minimal() + ggtitle("a")
qq <- ggs_density(ch, rug = T) + theme_minimal() + ggtitle("b")

#traces <- grid.arrange(q, arrangeGrob(p, qq, ncol = 1, nrow =1),  ncol = 1)
traces <- grid.arrange(p, qq, ncol=2, nrow=1)

# have traces for the complex model
experiment.dir <- "conifer.HierarchicalPhyloModel-t6N3LgKL.exec"
setwd(experiment.dir)

# mrbayes has only these
# r(A<->C)  r(A<->G)	r(A<->T)	r(C<->G)	r(C<->T)	r(G<->T)

# 
rate.simple <- read.coda("parameters-csv/CODAchain1.txt", "parameters-csv/CODAindex.txt")

rate.simple <- rate.simple[, c(2,3,4,6, 7, 10, 14, 15, 16, 17)]
#ch <- mcmc(cbind(alpha0, discount, kappa, nu))
ch <- mcmc(cbind(rate.simple))

ch <- ggs(ch)
p <- ggs_traceplot(ch) + theme_minimal() + ggtitle("a")
qq <- ggs_density(ch, rug = T) + theme_minimal() + ggtitle("b")

#traces <- grid.arrange(q, arrangeGrob(p, qq, ncol = 1, nrow =1),  ncol = 1)
traces <- grid.arrange(p, qq, ncol=2, nrow=1)

# not enough, read the actual ones

setwd("/Users/sohrab/Me/Apply/Canada Apply/Courses/Third Semester/conifer/results/all/")
f <- "conifer.TestPhyloModel-0hKG30ua.exec/parameters-csv"
print(f)
p <- list.files(f)
p <- p[grepl("q\\(", p)]

ll <- read.csv(file.path(f, p[1]), row.names = 1, header = T)
for (q.file in p[-c(1)]) {
  ff <- read.csv(file.path(f, q.file), row.names = 1, header = T)
  #print(head(ff))
  ll <- cbind(ll, ff)
}

# remove unnecessary from ll
head(ll)
ll <- ll[, c(1, 2, 3, 5, 6, 9)]

# add stationary
p <- list.files(f)
p <- p[grepl("statio", p)]

lll <- read.csv(file.path(f, p[1]), row.names = 1, header = T)
for (q.file in p[-c(1)]) {
  ff <- read.csv(file.path(f, q.file), row.names = 1, header = T)
  #print(head(ff))
  lll <- cbind(lll, ff)
}
head(lll)

# combine stationary and qs
######################################################################
llll <- cbind(ll, lll)
colnames(llll) <- colnames(kk)
colnames(ll)
colnames(kk)
######################################################################
ch <- mcmc(cbind(llll))
ch <- ggs(ch)
p <- ggs_traceplot(ch) + theme_minimal() + ggtitle("a")
qq <- ggs_density(ch, rug = T) + theme_minimal() + ggtitle("b")

#traces <- grid.arrange(q, arrangeGrob(p, qq, ncol = 1, nrow =1),  ncol = 1)
traces <- grid.arrange(p, qq, ncol=2, nrow=1)
######################################################################




# mrbayes simple GTR model
setwd("/Users/sohrab/Me/Apply/Canada\ Apply/Courses/Third\ Semester/conifer/extras/mrbayes/FES_4_batch.GTR.simple")
ff  <- "FES_4_batch.GTR.simple.nex.run1.p"
#ff  <- "FES_8_batch.GTR.two.nex.run1.p"
kk <- read.table(ff, skip=1, header=T)
head(kk)
kk <- kk[, -c(1:3)]
head(kk)

colnames(kk) <- colnames(rate.simple)
head(kk)
kk <- kk[-c(1:100), ]

ch <- mcmc(cbind(kk))

ch <- ggs(ch)
p <- ggs_traceplot(ch) + theme_minimal() + ggtitle("a")
qq <- ggs_density(ch, rug = T) + theme_minimal() + ggtitle("b")

#traces <- grid.arrange(q, arrangeGrob(p, qq, ncol = 1, nrow =1),  ncol = 1)
traces <- grid.arrange(p, qq, ncol=2, nrow=1)
