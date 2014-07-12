# ess generic
library(lattice)

# input is a combined, melted data.frame that contains
# 1. parameter, ess-value, source
make.ess.barchart <- function(ess, numberOfRuns) {

  # TODO: give the 
  ess.barchart <- barchart(ess_value~parameter, ess, groups=source, auto.key=T, 
           main=paste0("MCMC with n = ", numberOfRuns), xlab="Parameters", ylab="ESS")
  
  trellis.device(device="png", filename="essHead2Head.png")
  print(ess.barchart)
  dev.off()
}

combine.ess <- function(ess1, ess2) {
  d1 <- rbind(ess1, ess2)
  d1$source <- rep(c("mrbayes", "conifer"), times=c(nrow(ess1), nrow(ess2)))
  d1
}


fake.data <- function() {
  data.frame(parameter= rep(paste0("c", seq(1:10)), times=2) , ess_value=sample(20), source=rep(c("mrbayes", "conifer"), each=10))  
}

test <- function() {
  data <- fake.data()
  make.ess.barchart(data, 1000)
}

test()