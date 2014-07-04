# define functions
setClass("mrbayesbatch", representation(lines="list"), 
         prototype(lines=list("begin mrbayes;", "end;")))
setGeneric("getLines", function(object) standardGeneric("getLines"))
setMethod("getLines", "mrbayesbatch", 
          function(object) {            
            print(object@lines)
          })

setGeneric("writeToDisk", function(x, filePath) standardGeneric("writeToDisk"))
setMethod("writeToDisk", "mrbayesbatch", 
          function(x, filePath) {            
            writeLines(unlist(x@lines), filePath)
          })

setGeneric("MCMC<-", function(x, value) standardGeneric("MCMC<-"))
setReplaceMethod("MCMC", "mrbayesbatch", function(x, value) {
  x@lines <- append(x@lines, value, after=(length(x@lines) - 1))
  x
})

# command shouldn't end with ";". It will be added automatically.
setGeneric("addCommand<-", function(x, value) standardGeneric("addCommand<-"))
setReplaceMethod("addCommand", "mrbayesbatch", function(x, value) {
  x@lines <- append(x@lines, paste0(value, ";"), after=(length(x@lines) - 1))
  x
})

setGeneric("setKeyValue<-", function(x, value) standardGeneric("setKeyValue<-"))
setReplaceMethod("setKeyValue", "mrbayesbatch", function(x, value) {
  addCommand(x) <- paste0(value$key, " = ", value$value)
  x
})

# value should be a list
setGeneric("addMultiPartCommand<-", function(x, value) standardGeneric("addMultiPartCommand<-"))
setReplaceMethod("addMultiPartCommand", "mrbayesbatch", function(x, value) {
  addCommand(x) <- paste(value, collapse = " ")
  x
})

keyValueStringFromList <- function(some.list) {
  l <- lapply(seq(length(some.list)), function(i)  paste0(names(some.list)[i], "=", some.list[i]))
  paste0(l, collapse = " ")
}

# value should be list(commandString, list(key1=value1, key2=value2, ...))
setGeneric("setTags<-", function(x, value) standardGeneric("setTags<-"))
setReplaceMethod("setTags", "mrbayesbatch", function(x, value) {
  addCommand(x) <- paste(value$command, keyValueStringFromList(value$list))
  x
})

setGeneric("dataFile<-", function(x, value) standardGeneric("dataFile<-"))
setReplaceMethod("dataFile", "mrbayesbatch", function(x, value) {
  x@lines <- append(x@lines, paste0("execute", " ", value, ";"), after=(length(x@lines) - 1))
  x
})

