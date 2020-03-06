library(ggplot2)
library(forcats)
library(pROC)
library(reshape2)
library(cowplot)
library(glmnet)
library(caret)
library(fmsb)
theme_set(theme_cowplot())


compareEval<-function(authorComp,finalCrit,plotAll,phenoDefs,fileDefs,trainPhenos,diseaseTrainList,cp,
			shortRead,covarDefs,includeEnet,goNet){

source("compareEval.1.2.source.R")
set.seed(rand_addon)

pdf(paste(authorComp,"compare","pdf",sep='.'))
options(warn=1)

possibleCrits=c("pval","totalOR","groupOR","splitOR","Prev","AUC")
stopifnot(finalCrit %in% possibleCrits)
critIndex=which(possibleCrits %in% finalCrit)

#check all of the current files to see if the authorComp is within
authorFiles=c()
for(f in names(diseaseTrainList)){
  fOpen = file(paste0("diseaseFiles/",f),'r')
  firstLine=readLines(fOpen,n=1)
  close(fOpen)
  firstLine=strsplit(firstLine,split='\t')[[1]]
  if(authorComp %in% firstLine){
    authorFiles=c(authorFiles,f)
  }
}

#make a matrix to contain all of the scores from the files detected above
allScores=matrix(nrow=nrow(diseaseTrainList[[1]]),ncol=length(authorFiles))
colnames(allScores)=authorFiles
for(f in authorFiles){
  goodDiseaseFile=diseaseTrainList[[which(names(diseaseTrainList)==f)]]
  if(any(is.na(goodDiseaseFile[nrow(goodDiseaseFile),]))){
    goodDiseaseFile <- rbind(apply(goodDiseaseFile[complete.cases(goodDiseaseFile),],2,mean), goodDiseaseFile[1:(nrow(goodDiseaseFile)-1),])
  }
  allScores[,which(colnames(allScores)==f)]=goodDiseaseFile[,which(colnames(goodDiseaseFile)==authorComp),drop=T]
}

#check if any score files do not have any nonzero values
checkZero=apply(allScores,2,sum)
print(sum(checkZero==0))
if(sum(checkZero==0)>0){
  badFile=names(checkZero)[checkZero==0]
  allScores=allScores[,-which(checkZero==0)]
  authorFiles=authorFiles[-which(authorFiles %in% badFile)]
}

#normalize the scores
allScores=apply(allScores,2,function(x) (x-mean(x))/sd(x)) #maybe don't do this 

#in the decoder pull out the line corresponding to the author
saveRDS(phenoDefs, "phenoDefs.RDS")
saveRDS(authorComp, "authorComp.RDS")
decoder=phenoDefs[phenoDefs[,1]==authorComp,]
decoder=unname(unlist(decoder))
print("still good")
names(decoder)=colnames(phenoDefs)
print("did good")

#do the same for covariates
covarDecoder=covarDefs[covarDefs[,1]==authorComp,]
covarDecoder=unname(unlist(covarDecoder))
names(covarDecoder)=colnames(covarDefs)

#pull out info from fileDefs
key=fileDefs[fileDefs[,1] %in% authorFiles,]

#get the easy name of the trait
currentTrait=decoder[2]

#count number of time in all reporting the disease code comes up
predLength=nrow(allScores)
if(is.null(predLength)){predLength=length(allScores)}
response=rep(0,predLength)


#get the case/control response for the pheno
for(phenIndex in 1:length(trainPhenos)){
  if(names(trainPhenos)[phenIndex] %in% names(decoder)[!is.na(decoder)]){
    fullKeyword = decoder[names(decoder) == names(trainPhenos)[phenIndex]]
    includeKeyword = strsplit(fullKeyword,";",fixed=TRUE)[[1]][1]
    excludeKeyword = strsplit(fullKeyword,";",fixed=TRUE)[[1]][2]
    if(includeKeyword != "X"){
      for(keyword in strsplit(includeKeyword,"|",fixed=TRUE)[[1]]){
        if( any(grepl(keyword,colnames(trainPhenos[[phenIndex]]))) ){
          response = response + rowSums(trainPhenos[[phenIndex]][,grepl(keyword,colnames(trainPhenos[[phenIndex]])),drop=F])
        }
      }
    }
    if(!is.na(excludeKeyword)){
      for(keyword in strsplit(excludeKeyword,"|",fixed=TRUE)[[1]]){
        if( any(grepl(keyword,colnames(trainPhenos[[phenIndex]]))) ){
          response = response - rowSums(trainPhenos[[phenIndex]][,grepl(keyword,colnames(trainPhenos[[phenIndex]])),drop=F])
        }
      }
    }
  }
}

#THIS IS WHERE THINGS CHANGE!!! on NOV 12 2019 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#response[response<1]=0      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
response[response>0]=1      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
response=factor(response)   #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print("STARTING!")
print(covarDecoder)

covarResp=c()
itemsRemove=rep(0,length(response))
#get the covariate vector, data
if(length(covarDecoder) > 3){
  #for when a medical condiation is a covariate
  for(phenIndex in 1:length(trainPhenos)){
    if(names(trainPhenos)[phenIndex] %in% names(covarDecoder)){
      fullKeyword = covarDecoder[names(covarDecoder) == names(trainPhenos)[phenIndex]]
      for(keyword in strsplit(fullKeyword,"|",fixed=TRUE)[[1]]){
        if(substr(keyword,1,1)=="-"){
          realKeyword <- substr(keyword,2,nchar(keyword))
          if( any(grepl(realKeyword,colnames(trainPhenos[[phenIndex]]))) ){
            intermedResp =  rowSums(trainPhenos[[phenIndex]][,grep(realKeyword,colnames(trainPhenos[[phenIndex]]),value=T),drop=F])
            itemsRemove = itemsRemove + intermedResp
            itemsRemove[itemsRemove>1] <- 1
          } else {
            print("keyword not in")
            print(keyword)
          }
        } else {
          if( any(grepl(keyword,colnames(trainPhenos[[phenIndex]]))) ){
            intermedResp =  rowSums(trainPhenos[[phenIndex]][,grep(keyword,colnames(trainPhenos[[phenIndex]]),value=T),drop=F])
            intermedResp[intermedResp>1] = 1
            #covarResp = covarResp + intermedResp
            covarResp <- cbind(covarResp,factor(intermedResp))  #not sure which covar system is best
          } else {
            print("keyword not in")
            print(keyword)
          }
        }
      }
    }
  }
  if(!is.null(ncol(covarResp))){
    covarResp <- covarResp[,apply(covarResp,2,function(x) length(unique(x))>1),drop=F]
    colnames(covarResp)=paste("covarResp",1:ncol(covarResp),sep="")
    dfMain=data.frame(covar1,covarResp,response)
  } else {
    dfMain=data.frame(covar1,response)
  }

  if(any(itemsRemove==1)){
    dfMain <- dfMain[itemsRemove==0,]
    allScores <- allScores[itemsRemove==0,]
  }
} else {
  #when you just have age, sex, PCs
  dfMain=data.frame(covar1,response)
}


#subset the analysis by sex if nescessary
if(decoder[3]!="A"){
  covar1 <- covar1[dfMain$sex==decoder[3],]
  allScores <- allScores[dfMain$sex==decoder[3],]
  dfMain <- dfMain[dfMain$sex==decoder[3],]
  dfMain <- dfMain[,-3]
}

totalResMat=c()
totalHiMat=c()
totalLoMat=c()
totalRocDf=c("tpr"=0,"fpr"=0,"group"=factor(0))
constantPlotAll=plotAll
#iterate over each method (clump,ldpred,etc.) (NOTE: will have to change later to compare both methods at end)
print("BEGINNING THE METHODS FOR LOOP")
for(method in unique(key$methods)){
  print("METHODS")
  print(method)
  newKey=key[key$method==method,]
  plotAll=constantPlotAll  

  if(method=="clump" | method=="oldClump"){
    paramNames=c("pval","rho","none")
  } else if (method=="ldpred"){
    paramNames=c("rho","nothing","none")
  } else if (method=="winnersCurseLike"){
    paramNames=c("nothing","nothing","none")
  } else if (method=="winnersCurseLasso"){
    paramNames=c("lambda","nothing","none")
  } else if (method=="winnersCurse2d"){
    paramNames=c("pval","rho","none")
  } else if (method=="tweedy"){
    paramNames=c("method","nothing","none")
  } else if (method=="sumer"){
    paramNames=c("method","nothing","none")
  } else if (method=="prscs"){
    paramNames=c("phi","nothing","none")
  } else if (method=="grabld"){
    paramNames=c("method","nothing","none")
  } else if (method=="lassosum"){
    paramNames=c("method","nothing","none")
  } else if (method=="annoPred"){
    paramNames=c("m1","m2","none")
  } else if (method=="LDPredFunct"){
    paramNames=c("method","nothing","none")
  } else if (method=="sblup"){
    paramNames=c("m1","m2","none")
  } else if (method=="sbayesr"){
    paramNames=c("m1","m2","none")
  } else if (method=="oldLdpred"){
    paramNames=c("m1","none","none")
  }
  
  #iterating over each parameter
  filesChecked=c()
  for(paramCol in 3:(ncol(key))){
    print("paramcol")
    print(paramCol)


    if(length(unique(newKey[,paramCol]))>1 | paramCol==3 | length(unique(filesChecked)) < nrow(newKey) ){
      print("past the unique key if")
      print("files checked")
      print(filesChecked)
      print("firstCond")
      print(length(unique(newKey[,paramCol]))>1)
      print("secondCond")
      print(paramCol==3)
      print("thirdCond")
      print(length(filesChecked)<nrow(newKey))
      #subset the key to only include the parameter that is currently changing
      #if we are at the last parameter, just include all the files that have not been checked yet

      if(paramCol==ncol(key)){
        if(length(unique(filesChecked)) < nrow(key)){
          workKey <- newKey[!(newKey$fName %in% filesChecked),]
          plotAll=FALSE
        }
      } else if (paramCol < ncol(key)) {
        #for the paramter extract files from newKey so that only one thing is changing
        otherCol <- ifelse(paramCol==4,3,4)
        otherTable=table(newKey[,otherCol])
        mostCommon=names(which(max(otherTable)==otherTable)[1]) #changed this from numeric to character
        workKey=newKey[newKey[,otherCol]==mostCommon,] #changeed this from numeric to character
        filesChecked=c(filesChecked,workKey[,1])
      } else {
        print("breaking")
        break
      }

      print("work key!!!")
      print(workKey)
      

      predDf=allScores[,colnames(allScores) %in% workKey$fName, drop=F]
      predDf=predDf[,order(colnames(predDf))[rank(workKey[,1])], drop=F]

      for(i in 1:ncol(predDf)){
        if(mean(predDf[dfMain$response==1,i]) < mean(predDf[dfMain$response==0,i])){
          print("swap in comp")
          predDf[,i]=predDf[,i]*-1
        }
      }
      
      resThings <- makeResMat(predDf, dfMain, workKey, paramCol, paramNames, currentTrait, method, plotAll) 

      totalResMat=rbind(totalResMat,cbind(workKey,resThings[[1]]))
      totalHiMat=rbind(totalHiMat,cbind(workKey,resThings[[2]]))
      totalLoMat=rbind(totalLoMat,cbind(workKey,resThings[[3]]))
      totalRocDf=rbind(totalRocDf,resThings[[4]])
    }   
  }
}

print("OUT OF METHODS FLOR LOOP")


#######################################################################################################################
# Try enet regression with all scores 
print("DOING ENET REGRESSION")
dfList <- validation(dfMain, "cross", 3)
print("288")
aucCovar <- matrix(0,nrow=3,ncol=3)
tprCovar <- matrix(0,nrow=1000,ncol=3)
fprCovar <- matrix(0,nrow=1000,ncol=3)
for(iFold in 1:3){
  dfTrain <- dfList[[1]][[iFold]]
  dfTest <- dfList[[2]][[iFold]]
  print(head(dfTrain))
  print(head(dfTest))
  covarLogistic <- glm(response ~ . ,data=dfTrain,family="binomial")
  print("did covarLogistic")
  dfTest$prob=predict(covarLogistic,dfTest,type=c("response"))
  g=pROC::roc(response ~ prob,data=dfTest)
  tprCovar[,iFold]=g$sensitivities[round(seq(1,length(g$sensitivities),length.out = 1000))]
  fprCovar[,iFold]=(1-g$specificities[round(seq(1,length(g$specificities),length.out = 1000))])
  aucCovar[iFold,] <- as.numeric(ci.auc(g))
}
aucCovar <- apply(aucCovar,2,mean)
fprCovar <- apply(fprCovar,1,mean)
tprCovar <- apply(tprCovar,1,mean)

print("308")
#get index of the best parameter set for each method
totalResMat=totalResMat[!duplicated(totalResMat$fName),]
totalResMat <- totalResMat[!is.na(totalResMat$AUC),]
bestIndices <- rep(0,length(unique(totalResMat$methods)))
i=1
for(umethod in unique(totalResMat$methods)){
  print(umethod)
  possAuc <- totalResMat[totalResMat$methods == umethod, 12]
  print(possAuc)
  if(any(is.finite(possAuc))){
  umethodAuc <- max(possAuc[is.finite(possAuc)])
  #umethodAuc <- max(totalResMat[totalResMat$methods == umethod, 12])
  print(umethodAuc)
  print(totalResMat$AUC)
  bestIndices[i] <- which(colnames(allScores) == totalResMat$fName[which(totalResMat$AUC==umethodAuc)[1]])
  } else {
  bestIndices[i] <- which(colnames(allScores) == colnames(allScores)[grepl(umethod, colnames(allScores))])
  }
  i=i+1
}

#prepare dfs for the enet
dfMain$array <- as.numeric(dfMain$array)
if(decoder[3]=="A"){
  dfMain$sex <- as.numeric(as.factor(dfMain$sex))
}

if("covarResp" %in% colnames(dfMain)){
  for(i in grep("covarResp",colnames(dfMain))){
    print(colnames(dfMain)[i])
    dfMain[,i] <- as.numeric(as.factor(dfMain[,i]))
  }
}

#fullDfInd <- as.matrix(cbind(dfMain[,-ncol(dfMain)],allScores))
#bestDfInd <- as.matrix(cbind(dfMain[,-ncol(dfMain)],allScores[,bestIndices]))

print("335")
#ENET BEST ################################
#dep <- as.matrix(dfMain$response)
print(head(allScores))
print(head(dfMain))
print("341")
if(sum(totalResMat$AUC > 0.01)>2){
  aucCut <- 0.01
} else if(sum(totalResMat$AUC > 0)>2){
  aucCut <- 0
} else {
  aucCut <- -1
}

fullDfInd <- as.matrix(cbind(dfMain[,-ncol(dfMain)],allScores))
print(dim(fullDfInd))
fullDfInd <- fullDfInd[,colnames(fullDfInd) %in% totalResMat$fName[totalResMat$AUC > aucCut]]
print("fullDfInd")
print(head(fullDfInd))
print("totalResMat")
print(totalResMat$AUC)
print(totalResMat)
print("368")
cvAns <- cv.glmnet(y=dfMain$response, x=fullDfInd, family="binomial",nfolds=3,maxit=5000)
bestLamb <- cvAns$lambda.min

print("enet - best")
bestEnet <- enetProc(fullDfInd, dfMain, bestLamb, goNet)
print("done enet best") #this is the last spot parallele made it to

#ENET FULL ################################
fullDfInd <- as.matrix(cbind(dfMain[,-ncol(dfMain)],allScores))
fullDfInd <- fullDfInd[,colnames(fullDfInd) %in% totalResMat$fName[order(totalResMat$AUC, decreasing=T)][1:5]]
cvAnsFull <- cv.glmnet(y=dfMain$response, x=fullDfInd, family="binomial",nfolds=3,maxit=5000)
fullLamb <- cvAnsFull$lambda.min

print("enet - full")
allEnet <- enetProc(fullDfInd, dfMain, fullLamb, goNet)





enetPredDf <- cbind(bestEnet[[3]], allEnet[[3]])
enetWorkKey <- data.frame(fName=c("bestEnet","allEnet"),methods="enet",p1=c(0,1),p2=c(0,0),p3=c(0,0))
enetResThings <- makeResMat(enetPredDf,dfMain,enetWorkKey, 3, c("realSplit","none","none"), currentTrait, "enet", plotAll, "enet")
if(includeEnet){
  totalResMat <- rbind(totalResMat,cbind(enetWorkKey,enetResThings[[1]]))
  totalHiMat <- rbind(totalHiMat,cbind(enetWorkKey,enetResThings[[2]]))
  totalLoMat <- rbind(totalLoMat,cbind(enetWorkKey,enetResThings[[3]]))
}

plotDf <- data.frame(rbind(bestEnet[[1]],allEnet[[1]],cbind(tpr=tprCovar,fpr=fprCovar)))
plotDf$type <- c(rep("best",nrow(bestEnet[[1]])),rep("all",nrow(allEnet[[1]])),rep("covar",length(tprCovar)))
thePlot <- ggplot(plotDf,aes(fpr,tpr))+geom_line(aes(color=type))+geom_abline(slope=1)+
  labs(x="False Positive Rate",y="True Positive Rate",title="ROC for Enet of All Methods",
       caption=paste("AUC Imp: all =",round(allEnet[[2]]-aucCovar[2],3)," ~ best =",round(bestEnet[[2]]-aucCovar[2],3)))
plot(thePlot)

##########################################################################################################################
#compare a specified statistic at the end

print("final comparison")
totalRocDf=totalRocDf[-1,]
toTake=which(totalResMat$AUC==max(totalResMat$AUC))[1]
totalRocDf=totalRocDf[(((toTake-1)*1000)+1):(toTake*1000),]
totalAuc=totalResMat$AUC[toTake]

colnames(totalResMat)[6+critIndex]="toComp"
print("ONCE MORE")
print(totalResMat)
totalResMat=totalResMat[order(totalResMat$toComp),]
totalResMat$parameter=paste(totalResMat$method, as.character(totalResMat$p1), as.character(totalResMat$p2), as.character(totalResMat$p3),sep='-')
totalResMat$parameter=factor(totalResMat$parameter,levels = totalResMat$parameter[order(totalResMat$toComp)])

print("366")
totalLoMat=totalLoMat[!(duplicated(totalLoMat$fName)),]
totalHiMat=totalHiMat[!(duplicated(totalHiMat$fName)),]
colnames(totalLoMat)[6+critIndex]="toComp"
colnames(totalHiMat)[6+critIndex]="toComp"
totalResMat=totalResMat[order(totalResMat$fName),]
totalLoMat=totalLoMat[order(totalLoMat$fName),]
totalHiMat=totalHiMat[order(totalHiMat$fName),]

print("376")
if(finalCrit=="AUC"){finalCrit="AUC Improvement"}
thePlot=ggplot(totalResMat,aes(parameter,toComp))+geom_point(aes(color=methods),size=4)+
  geom_errorbar(aes(ymin=totalLoMat$toComp,ymax=totalHiMat$toComp,width=0.2), color="grey50")+
  labs(x="Parameters",y=finalCrit,title=paste("Total Analysis of",currentTrait), caption=paste("best:",round(max(totalResMat$toComp),3)))+
  guides(color=FALSE)+theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5), plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
  scale_color_manual(values=cp)
plot(thePlot)


thePlot=ggplot(totalResMat,aes(parameter,toComp))+geom_point(aes(color=methods),size=4)+
  geom_errorbar(aes(ymin=totalLoMat$toComp,ymax=totalHiMat$toComp,width=0.2), color="grey50")+
  labs(x="Parameters",y=finalCrit,title=paste("Total Analysis of",currentTrait), caption=paste("best:",round(max(totalResMat$toComp),3)))+
  theme(axis.text.x=element_blank(), legend.position="bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 8))+
  scale_color_manual(values=cp)
plot(thePlot)



print("385")
bestName=totalResMat$fName[totalResMat$toComp==max(totalResMat$toComp)][1]
if(bestName=="allEnet"){
  bestScore <- enetPredDf[,2]
  enetBetas <- bestEnet[[4]]
  enetBetas <- enetBetas[enetBetas!=0]
} else if(bestName=="bestEnet"){
  bestScore <- enetPredDf[,1]
  enetBetas <- allEnet[[4]]
  enetBetas <- enetBetas[enetBetas!=0]
} else {
  bestScore=allScores[,colnames(allScores)==bestName]
  enetBetas <- NULL
}

print("398")
goBackDf=cbind(dfMain,bestScore)
colnames(goBackDf)[ncol(goBackDf)]="predictor"
goBackDf=goBackDf[,c(1:(ncol(dfMain)-1),(ncol(dfMain)+1),ncol(dfMain))]
if(mean(goBackDf$predictor[goBackDf$response==1]) < mean(goBackDf$predictor[goBackDf$response==0])){
  goBackDf$predictor=goBackDf$predictor*-1
  revScore=TRUE
} else {
  revScore=FALSE
}

print("409")
goBackDf <- cbind(covar1,predictor=goBackDf$predictor,response=goBackDf$response)
if(decoder[3]!="A"){
  goBackDf <- goBackDf[,-3]
}
totalRocDf=totalRocDf[,1:2]
totalRocDf$group="train"
dev.off()

return(list(currentTrait,bestName,goBackDf,decoder,totalRocDf,totalAuc[1],revScore,totalResMat, enetBetas))

}


