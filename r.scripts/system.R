#! /usr/bin/Rscript
# set a CRAN mirror
 local({r <- getOption("repos")
       r["CRAN"] <- "http://cran.rstudio.com/"
       options(repos=r)})

library(tools)
library(base)
library(utils)
library(stats)
library(methods)
library(grDevices)
library(graphics)
library(datasets)
library(coda)
library(rstudio)
# system("mb -with-beagle=no")
#system("ls")


print("yes...")
#
system("mb" )

#system("./p.sh")
writeLines("hello")

 
