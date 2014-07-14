# conifer benchmarking
require(XML)






conifer.driver.function <- function(tree.file.name, alignment.file.name) {

  # run conifer with given input
  data <- xmlParse("~/conifer/.classpath")
  ldata <- xmlToList(data)
  paths <- lapply(ldata, function(x) x['path'])
  paths <- append(".", paths, after=0)
  classpaths <- paste0(paths, collapse = ":")
  
  
  setwd("/home/sohrab/conifer/build/classes/main")
  system(paste0("java -classpath ", classpaths, " ", "conifer.TestPhyloModel"))
  
  # calculate ESS and ESS per second
  
  
  
}