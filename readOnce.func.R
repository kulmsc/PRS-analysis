library(stringr)

readOnce <- function(){
  #For compare
  decoder=read.table("fileDecoder",stringsAsFactors=F,header=T)
  pheno1=read.table("pheno.phase1",stringsAsFactors = F)
  key=read.table("key",stringsAsFactors = F)
  colnames(key)=c("fName","methods",paste0("p",1:(ncol(key)-3)),"type")
  
  
  diseaseTrainList=list()
  i=1
  for(f in list.files("diseaseFiles/","train")){
    print(f)
    diseaseTrainList[[i]] <- read.table(paste0("diseaseFiles/",f),stringsAsFactors=F,header=T)
    toChange=grepl(".",colnames(diseaseTrainList[[i]]),fixed=T)
    colnames(diseaseTrainList[[i]])[toChange]=gsub(".","-",colnames(diseaseTrainList[[i]])[toChange],fixed=T)
    i=i+1
  }
  
  diseaseTestList=list()
  i=1
  for(f in list.files("diseaseFiles/","test")){
    print(f)
    diseaseTestList[[i]] <- read.table(paste0("diseaseFiles/",f),stringsAsFactors=F,header=T)
    toChange=grepl(".",colnames(diseaseTestList[[i]]),fixed=T)
    colnames(diseaseTestList[[i]])[toChange]=gsub(".","-",colnames(diseaseTestList[[i]])[toChange],fixed=T)
    i=i+1
  }
  
  supportTrainList=list()
  i=1
  for(f in list.files("supportFiles/","train")){
    supportTrainList[[i]] <- read.table(paste0("supportFiles/",f),stringsAsFactors=F,header=T)
    toChange=grepl(".",colnames(supportTrainList[[i]]),fixed=T)
    colnames(supportTrainList[[i]])[toChange]=gsub(".","-",colnames(supportTrainList[[i]])[toChange],fixed=T)
    i=i+1
  }
  
  supportTestList=list()
  i=1
  for(f in list.files("supportFiles/","test")){
    supportTestList[[i]] <- read.table(paste0("supportFiles/",f),stringsAsFactors=F,header=T)
    toChange=grepl(".",colnames(supportTestList[[i]]),fixed=T)
    colnames(supportTestList[[i]])[toChange]=gsub(".","-",colnames(supportTestList[[i]])[toChange],fixed=T)
    i=i+1
  }
  
  
  #for pheno
  pheno2=read.table("pheno.phase2",stringsAsFactors = F)
  
  
  
  #for boost
  allSupportDecode=read.table("supportDecoder",stringsAsFactors = F,header=T)
  oldTrait=allSupportDecode$trait
  dups=unique(allSupportDecode[duplicated(allSupportDecode[,2]),2])
  for(dup in dups){
    offending=allSupportDecode[allSupportDecode[,2]==dup,2]
    allSupportDecode[allSupportDecode[,2]==dup,2] <- paste0(offending,as.character(1:length(offending)))
  }
  allSupportDecode$oldTrait=oldTrait
  
  allDiseaseDecoder=read.table("fileDecoder",stringsAsFactors=F,header=T)
  oldTrait=allDiseaseDecoder$trait
  dups=unique(allDiseaseDecoder[duplicated(allDiseaseDecoder[,2]),2])
  for(dup in dups){
    offending=allDiseaseDecoder[allDiseaseDecoder[,2]==dup,2]
    allDiseaseDecoder[allDiseaseDecoder[,2]==dup,2] <- paste0(offending,as.character(1:length(offending)))
  }
  allDiseaseDecoder$oldTrait=oldTrait
  
  
  #for corr
  allScoresTrain=read.table("allScoresCompiled.train.gz",header=T)
  allScoresTest=read.table("allScoresCompiled.test.gz",header=T)
  allSupportsTrain=read.table("subsetSupportsCompilted.train.gz",header=T)
  allSupportsTest=read.table("subsetSupportsCompilted.test.gz",header=T)
  allSupportsTest=allSupportsTest[,colnames(allSupportsTest) %in% colnames(allSupportsTrain)]
  
  allDiseaseFiles=list.files("diseaseFiles/",pattern="train")
  commonNames=c()
  for(f in allDiseaseFiles){
    fOpen = file(paste0("diseaseFiles/",f),'r')
    firstLine=readLines(fOpen,n=1)
    close(fOpen)
    firstLine=strsplit(firstLine,split='\t')[[1]]
    commonNames = union(commonNames, firstLine)
  }
  
  #misc
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  allFiles <- list.files("diseaseFiles/","test")
  allMethods <- str_split(gsub(".","-",allFiles,fixed = T),"-",simplify=T)[,2]
  cp <- gg_color_hue(length(unique(allMethods)))
  names(cp)<-unique(allMethods)
  
  return(list(diseaseTrainList,diseaseTestList,supportTrainList,supportTestList,
              decoder,pheno1,key,pheno2,allSupportDecode,allDiseaseDecoder,allScoresTrain,allScoresTest,
              allSupportsTrain,allSupportsTest,commonNames,cp))
}

