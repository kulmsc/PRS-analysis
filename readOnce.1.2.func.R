library(stringr)
library(vroom)

#args <- commandArgs(trailingOnly=TRUE)
#shortRead <- args[1]
#nameApp <- args[2]

shortRead = FALSE
nameApp = "scott"


#the definitions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
phenoDefs=read.table("phenoDefs",stringsAsFactors=F,header=T)
fileDefs=read.table("fileDefs",stringsAsFactors = F)
covarDefs=read.table("covarDefs",stringsAsFactors=F,header=T)
colnames(fileDefs)=c("fName","methods",paste0("p",1:(ncol(fileDefs)-2)))
if(shortRead){phenoDefs <- phenoDefs[,1:5]}


#phenotypes +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
trainPhenos=vector("list", length = (ncol(phenoDefs) - 3))
testPhenos=vector("list", length = (ncol(phenoDefs) - 3))
names(trainPhenos)=colnames(phenoDefs)[4:ncol(phenoDefs)]
names(testPhenos)=colnames(phenoDefs)[4:ncol(phenoDefs)]
col=1
for(path in colnames(phenoDefs[4:ncol(phenoDefs)])){
  ans <- readRDS(paste0(nameApp,"/",nameApp,"Pheno.1.",path,".events.RDS"))
  ans <- ans[,colSums(ans) > 0]
  trainPhenos[[col]] <- ans
  
  ans <- readRDS(paste0(nameApp,"/",nameApp,"Pheno.2.",path,".events.RDS"))
  ans <- ans[,colSums(ans) > 0]
  testPhenos[[col]] <- ans
  rm(ans)
  col=col+1
}

save.image("readData.phenos.RData")
rm(trainPhenos)
rm(testPhenos)

#the covariates ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
index1=read.table(paste0(nameApp,"/",nameApp,".phenIndex.1"))
index2=read.table(paste0(nameApp,"/",nameApp,".phenIndex.2"))
covar1=readRDS("covars.1.RDS")
covar2=readRDS("covars.2.RDS")
covar1=covar1[index1[,1],]
covar2=covar2[index2[,1],]

#the colors +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
allFiles <- list.files("diseaseFiles/","test")
allMethods <- str_split(gsub(".","-",allFiles,fixed = T),"-",simplify=T)[,2]
cp <- gg_color_hue(length(unique(allMethods)))
names(cp)<-unique(allMethods)

#the scores ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
diseaseTrainList=list()
diseaseTestList=list()
i=1
nameFiles <- list.files("diseaseFiles/","train")
if(shortRead){nameFiles <- grep("clump|ldpred|grabld|winnersCurse2d", nameFiles, value = T)}
for(f in nameFiles){
  print(f)
  #diseaseTrainList[[i]] <- read.table(paste0("diseaseFiles/",f),stringsAsFactors=F,header=T)
  diseaseTrainList[[i]] <- vroom(paste0("diseaseFiles/",f))
  diseaseTrainList[[i]] <- diseaseTrainList[[i]][index1[,1], , drop=F]
  toChange=grepl(".",colnames(diseaseTrainList[[i]]),fixed=T)
  colnames(diseaseTrainList[[i]])[toChange]=gsub(".","-",colnames(diseaseTrainList[[i]])[toChange],fixed=T)
  i=i+1
}
names(diseaseTrainList) <- nameFiles
save.image("readData.trainScores.RData")
rm(diseaseTrainList)

i=1
for(f in nameFiles){ 
  f <- gsub("train","test",f)
  print(f)
  #diseaseTestList[[i]] <- read.table(paste0("diseaseFiles/",f),stringsAsFactors=F,header=T)
  diseaseTestList[[i]] <- vroom(paste0("diseaseFiles/",f))
  diseaseTestList[[i]] <- diseaseTestList[[i]][index2[,1], ,drop=F]
  print(dim(diseaseTestList[[i]]))
  toChange=grepl(".",colnames(diseaseTestList[[i]]),fixed=T)
  colnames(diseaseTestList[[i]])[toChange]=gsub(".","-",colnames(diseaseTestList[[i]])[toChange],fixed=T)
  i=i+1
}
names(diseaseTestList) <- gsub("train","test",nameFiles)
  

  
save.image("readyData.scores.RData")
