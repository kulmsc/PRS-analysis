library(parallel)


includeEnet <- TRUE
includeStack <- FALSE
goNet <- FALSE #when making an ensemble to use enet over logistic regression

source("singleFunc.R")
load("readData.phenos.RData")
load("readData.trainScores.RData")
load("readData.testScores.RData")
load("readData.supports.RData")




if(!includeStack){
  fileDefs <- fileDefs[grep("stackCT", fileDefs[,1], invert=T),]
}

authorList <- read.table("allAuthors",stringsAsFactors=F)
for(goAuthor in authorList[,1]){
  print(goAuthor)
  runLikeTheWind(goAuthor, includeEnet, goNet)
}
