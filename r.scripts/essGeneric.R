# ess generic
library(lattice)

# input is a combined, melted data.frame that contains
# 1. parameter, ess-value, source

make.bar.chart <- function(data, numberOfRuns, plotPath, ylab) {
  # TODO: remove this  
  #data <- data[data$source == "conifer", ]
  
  barchart <- barchart(ESS~parameter, data, groups=source, auto.key=T, 
                           main=paste0("MCMC with n = ", numberOfRuns), xlab="Parameters", ylab=ylab, scales=list(x=list(rot=90)))
  
  print(plotPath)
  trellis.device(device="png", filename=plotPath)
  print(barchart)
  dev.off()
}

make.ess.barchart <- function(ess, numberOfRuns, plotPath) {
  make.bar.chart(ess, numberOfRuns, file.path(plotPath, "essHead2Head.png"), ylab="ESS")
}

make.ess.per.second.barchart <- function(essPerSecond, numberOfRuns, plotPath) {
  colnames(essPerSecond) <- gsub("ESSPS", "ESS", colnames(essPerSecond))
  make.bar.chart(essPerSecond, numberOfRuns, file.path(plotPath, "essPerSecondHead2Head.png"), ylab="ESS per second")
}




combine.ess <- function(ess1, ess2) {
  
  colnames(ess1) <- c("ESS")
  colnames(ess2) <- c("ESS")
  
  ess1$parameter <- rownames(ess1)
  ess2$parameter <- rownames(ess2)
  
  d1 <- rbind(ess1, ess2)
  d1$source <- rep(c("mrbayes", "conifer"), times=c(nrow(ess1), nrow(ess2)))
  
  rownames(d1) <- NULL
  d1
}


fake.data <- function() {
  data.frame(parameter= rep(paste0("c", seq(1:10)), times=2) , ESS=sample(20), source=rep(c("mrbayes", "conifer"), each=10))  
}

test <- function() {
  data <- fake.data()
  make.ess.barchart(data, 1000)
}

#test()