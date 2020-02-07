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
phenNames <- colnames(phenoDefs)[4:ncol(phenoDefs)]
if(!"noncancer" %in% phenNames){
  phenNames <- c(phenNames, "noncancer", "cancer")
}
trainPhenos=vector("list", length = length(phenNames))
testPhenos=vector("list", length = length(phenNames))
names(trainPhenos)=phenNames
names(testPhenos)=phenNames
col=1
for(path in phenNames){
  ans <- readRDS(paste0(nameApp,"/",nameApp,"Pheno.1.",path,".score.events.RDS"))
  ans <- ans[,colSums(ans) > 0]
  trainPhenos[[col]] <- ans
  
  ans <- readRDS(paste0(nameApp,"/",nameApp,"Pheno.2.",path,".score.events.RDS"))
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
allMethods <- str_split(allFiles,fixed("."),simplify=T)[,2]
allMethods <- c(allMethods,"enet")
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
  diseaseTestList[[i]] <- vroom(paste0("diseaseFiles/",f))
  diseaseTestList[[i]] <- diseaseTestList[[i]][index2[,1], ,drop=F]
  print(dim(diseaseTestList[[i]]))
  toChange=grepl(".",colnames(diseaseTestList[[i]]),fixed=T)
  colnames(diseaseTestList[[i]])[toChange]=gsub(".","-",colnames(diseaseTestList[[i]])[toChange],fixed=T)
  i=i+1
}
names(diseaseTestList) <- gsub("train","test",nameFiles)
save.image("readData.testScores.RData")
rm(diseaseTestList)

#support files ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
supportTrainList=list()
i=1
for(f in list.files("supportFiles/","train")){
  #supportTrainList[[i]] <- read.table(paste0("supportFiles/",f),stringsAsFactors=F,header=T)
  supportTrainList[[i]] <- vroom(paste0("supportFiles/",f))
  supportTrainList[[i]] <- supportTrainList[[i]][index1[,1], ,drop=F]
  toChange=grepl(".",colnames(supportTrainList[[i]]),fixed=T)
  colnames(supportTrainList[[i]])[toChange]=gsub(".","-",colnames(supportTrainList[[i]])[toChange],fixed=T)
  i=i+1
}
names(supportTrainList) <- list.files("supportFiles/","train")

supportTestList=list()
i=1
for(f in list.files("supportFiles/","test")){
  #supportTestList[[i]] <- read.table(paste0("supportFiles/",f),stringsAsFactors=F,header=T)
  supportTestList[[i]] <- vroom(paste0("supportFiles/",f))
  supportTestList[[i]] <- supportTestList[[i]][index2[,1], ,drop=F]
  toChange=grepl(".",colnames(supportTestList[[i]]),fixed=T)
  colnames(supportTestList[[i]])[toChange]=gsub(".","-",colnames(supportTestList[[i]])[toChange],fixed=T)
  i=i+1
}
names(supportTestList) <- list.files("supportFiles/","test")


corrTrainRefs <- list()
corrTestRefs <- list()
i <- 1; j <- 1
for(f in list.files("scott/","supp")){
  if(strsplit(f,".",fixed=T)[[1]][2] == "1"){
    corrTrainRefs[[i]] <- readRDS(paste0("scott/",f))
    names(corrTrainRefs)[i] <- strsplit(f,".",fixed=T)[[1]][3]
    i <- i + 1
  } else if(strsplit(f,".",fixed=T)[[1]][2] == "2"){
    corrTestRefs[[j]] <- readRDS(paste0("scott/",f))
    names(corrTestRefs)[j] <- strsplit(f,".",fixed=T)[[1]][3]
    j <- j + 1
  }
}

save.image("readData.supports.RData")
