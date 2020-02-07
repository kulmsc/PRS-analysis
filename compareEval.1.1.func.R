library(ggplot2)
library(forcats)
library(pROC)
library(reshape2)
library(cowplot)
library(glmnet)
library(caret)


compareEval<-function(authorComp,finalCrit,plotAll,phenoDefs,fileDefs,trainPhenos,diseaseTrainList,cp,shortRead,covarDefs){
  
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
  print("MUST GO BACK TO CHECK IF THIS IS RIGHT!!!")
  goodDiseaseFile=diseaseTrainList[[which(names(diseaseTrainList)==f)]]
  if(any(is.na(goodDiseaseFile[nrow(goodDiseaseFile),]))){
    goodDiseaseFile[nrow(goodDiseaseFile),] <- apply(goodDiseaseFile[complete.cases(goodDiseaseFile),],2,mean)
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
decoder=phenoDefs[phenoDefs[,1]==authorComp,]
decoder=unname(unlist(decoder))
names(decoder)=colnames(phenoDefs)

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
response[response>1]=1
response[response<0]=0
response=factor(response)

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

  if(method=="clump"){
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
  }
  
  #iterating over each parameter
  filesChecked=c()
  for(paramCol in 3:(ncol(key))){
    print("paramcol")
    print(paramCol)
    if(length(unique(newKey[,paramCol]))>1 | paramCol==3 | length(filesChecked)<nrow(newKey) ){
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
        if(length(filesChecked) < nrow(key)){
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

      for(i in 1:ncol(predDf)){
        if(mean(predDf[dfMain$response==1,i]) < mean(predDf[dfMain$response==0,i])){
          print("swap in comp")
          print(i)
          print(table(dfMain$response))
          print(mean(predDf[dfMain$response==1,i]))
          print(mean(predDf[dfMain$response==0,i]))
          predDf[,i]=predDf[,i]*-1
        }
      }

      #for the all important df
      df=cbind(dfMain,predDf)
      print("THIS IS DF")
      print(head(df))
      
      #collect the pval and OR by logistic regression
      folds <- 3
      resMat=rep(list(data.frame(matrix(0,nrow=ncol(predDf),ncol=5,dimnames=list(NULL,c("index","pval","or","or2","or3"))))),folds)
      #resMat=lapply(resMat,function(x){ x[,1]<-factor(workKey[,paramCol]); x})
      sdLowMat=rep(list(data.frame(resMat[[1]])),folds)
      sdHiMat=rep(list(data.frame(resMat[[1]])),folds)
      aucs=rep(list(matrix(0,nrow=ncol(predDf),ncol=3)),folds)
      rocDf=rep(list(matrix(0,nrow=1000*ncol(predDf),ncol=3)),folds)
      aucCovar=rep(list(0),folds)

      #change from the current code to cross validation
      validation <- function(theType, folds, caseInBoth = 2000, controlMultiplier = 5){
        caseInds <- sample(which(df$response==1))
        controlInds <- sample(which(df$response==0))

        caseSplits <- ceiling(seq_along(caseInds) / ceiling(length(caseInds)/folds))
        controlSplits <- ceiling(seq_along(controlInds) / ceiling(length(controlInds)/folds))
        if(caseInBoth < length(caseInds)){caseInBoth <- length(caseInds)}
        if(length(controlInds) < controlMultiplier*length(caseInds)){controlMultiplier <- length(controlInds)/length(caseInds)}

        dfTrainList <- rep(list(0),folds)
        dfTestList <- rep(list(0),folds)

        for(i in 1:folds){
          if(theType=="cross"){
            dfTrainList[[i]] <- df[c(caseInds[caseSplits==i],controlInds[controlSplits==i]),]
            dfTestList[[i]] <- df[c(caseInds[caseSplits!=i],controlInds[controlSplits!=i]),]
      	  } else if(theType=="bag"){
            trainCaseSample <- sample(caseInds,round(caseInBoth/2))
            trainControlSample <- sample(controlInds,round((caseInBoth/2)*5))
            testCaseSample <- sample(setdiff(caseInds,trainCaseSample),round(caseInBoth/2))
            testControlSample <- sample(setdiff(controlInds,trainControlInds),round(caseInBoth/2)*5)
            dfTrainList[[i]] <- df[c(trainCaseSample, trainControlSample),]
            dfTestList[[i]] <- df[c(testCaseSample,testControlSample),]
          } 
        }
        if(all(c(sapply(dfTrainList, function(x) length(unique(x[,1]))),sapply(dfTestList, function(x) length(unique(x[,1]))))>1)){
          return(list(dfTrainList,dfTestList))
        } else {
	  return(validation(theType, folds, caseInBoth, controlMultiplier))
        }
      }

      dfList <- validation("cross", folds)

      print("STARTING FOLDS")
      for(iFold in 1:folds){
        print("fold")
        for(i in 1:ncol(predDf)){
          #normal logistic regression of score against response
          dfTrain <- dfList[[1]][[iFold]]
          dfTest <- dfList[[2]][[iFold]]
          dfMainTrain <- dfTrain[,grepl(paste(c("array","age","sex","PC1","PC2","PC3","PC4","response","covarResp"),collapse="|"),colnames(dfTrain)),drop=F]
          dfMainTest <- dfTest[,grepl(paste(c("array","age","sex","PC1","PC2","PC3","PC4","response","covarResp"),collapse="|"),colnames(dfTest)),drop=F]
          predDfTrain <- dfTrain[,!grepl(paste(c("array","age","sex","PC1","PC2","PC3","PC4","response","covarResp"),collapse="|"),colnames(dfTrain)),drop=F]
          predDfTest <- dfTest[,!grepl(paste(c("array","age","sex","PC1","PC2","PC3","PC4","response","covarResp"),collapse="|"),colnames(dfTest)),drop=F]

          completeLogistic <- glm(response ~ . ,data=cbind(dfMainTrain,pred=predDfTrain[,i]),family="binomial")
          compLogitCI=confint(completeLogistic,level=0.95)
          theCoefs <- coef(summary(completeLogistic))
          resMat[[iFold]][i,2]=-log10(theCoefs[nrow(theCoefs),4]) #pval
          resMat[[iFold]][i,3]=exp(theCoefs[nrow(theCoefs),1]) #odds ratio
          sdLowMat[[iFold]][i,3]=exp(compLogitCI[nrow(theCoefs),1])
          sdHiMat[[iFold]][i,3]=exp(compLogitCI[nrow(theCoefs),2])
        
          #odds ratio equation of top 5% to bottom 95%
          dfTrain=dfTrain[order(dfTrain[,ncol(dfMainTrain)+i],decreasing = T),]
          exposeGroup=dfTrain[1:round(nrow(dfTrain)*0.05),ncol(dfMainTrain)]
          safeGroup=dfTrain[(round(nrow(dfTrain)*0.05)+1):nrow(dfTrain),ncol(dfMainTrain)]
          se=sqrt((1/sum(exposeGroup==1))+(1/sum(exposeGroup==0))+(1/sum(safeGroup==1))+(1/sum(safeGroup==0)))
          resMat[[iFold]][i,4]=(sum(exposeGroup==1)*sum(safeGroup==0))/(sum(exposeGroup==0)*sum(safeGroup==1))
          sdLowMat[[iFold]][i,4]=exp(log(resMat[[iFold]][i,4])+(1.96*se))
          sdHiMat[[iFold]][i,4]=exp(log(resMat[[iFold]][i,4])-(1.96*se))
        
          #odds ratio equation of top 50% to bottom 50%
          exposeGroup=dfTrain[1:round(nrow(dfTrain)*0.5),ncol(dfMainTrain)]
          safeGroup=dfTrain[(round(nrow(dfTrain)*0.5)+1):nrow(dfTrain),ncol(dfMainTrain)]
          se=sqrt((1/sum(exposeGroup==1))+(1/sum(exposeGroup==0))+(1/sum(safeGroup==1))+(1/sum(safeGroup==0)))
          resMat[[iFold]][i,5]=(sum(exposeGroup==1)*sum(safeGroup==0))/(sum(exposeGroup==0)*sum(safeGroup==1))
          sdLowMat[[iFold]][i,5]=exp(log(resMat[[iFold]][i,5])+(1.96*se))
          sdHiMat[[iFold]][i,5]=exp(log(resMat[[iFold]][i,5])-(1.96*se))

          #roc/auc stuff
          #ROC AND AUC PLOTS!!!
          covarLogistic <- glm(response ~ . ,data=dfMainTrain,family="binomial")
          dfMainTest$prob=predict(covarLogistic,dfMainTest,type=c("response"))
          g=pROC::roc(response ~ prob,data=dfMainTest)
          aucCovar[[iFold]]=as.numeric(ci.auc(g))

          dfTest$prob=predict(completeLogistic, cbind(dfMainTest,pred=predDfTest[,i]), type=c("response"))
          g=pROC::roc(response ~ prob,data=dfTest)
          aucs[[iFold]][i,]=as.numeric(ci.auc(g$auc))
          rocDfTemp=cbind(g$sensitivities,1-g$specificities,rep(workKey[i,paramCol],length(g$specificities)))
          rocDf[[iFold]][(i*1000-999):(i*1000),]=rocDfTemp[seq(1,nrow(rocDfTemp),length.out = 1000),]
          #rocDf[[iFold]]$group=factor(rocDf[[iFold]]$group)
        }
      }  


      resMat=Reduce("+", resMat) / folds
      resMat[,1]=factor(workKey[,paramCol])
      sdLowMat=Reduce("+", sdLowMat) / folds
      sdLowMat[,1]=factor(workKey[,paramCol])
      sdHiMat=Reduce("+", sdHiMat) / folds
      sdHiMat[,1]=factor(workKey[,paramCol])
      aucs=Reduce("+", aucs) / folds
      rocDf=data.frame(Reduce("+", rocDf) / folds)
      colnames(rocDf) <- c("tpr","fpr","group")
      rocDf$group=factor(rocDf$group)
      aucCovar=Reduce("+", aucCovar) / folds
      print("WE HAVE AUCCOVAR!")
      print(aucCovar)
      print("AND THIS IS AUCS")
      print(aucs)

      print("PLOTTING")
      print(plotAll)
      print("that was status of plotAll")
      print(resMat)
      #PVAL PLOT!
      thePlot=ggplot(resMat,aes(index,pval))+geom_col(aes(fill=index))+
        labs(x=paramNames[paramCol-2],y="-log10(Pval)",title=paste("Analysis of",currentTrait),
             subtitle=paste("Using method",method,"- comparing",paramNames[paramCol-2]))+
        scale_fill_discrete(guide=FALSE)
      if(plotAll){plot(thePlot)}
      

      #OR PLOT!
      orMat=melt(resMat[,-2],id.vars = "index")
      orMat$sdHi=melt(sdHiMat[,-2],id.vars="index")[,3]
      orMat$sdLo=melt(sdLowMat[,-2],id.vars="index")[,3]
      thePlot=ggplot(orMat,aes(variable,value,color=index))+geom_point(position=position_dodge(width=0.3), size=4)+
        geom_errorbar(aes(ymin=sdLo,ymax=sdHi),width=0.2,position = position_dodge(0.3))+
        labs(x="OR comparison",y="odds ratio per SD",title=paste("Analysis of",currentTrait),
          subtitle=paste("Using method",method,"- comparing",paramNames[paramCol-2]))+
        scale_x_discrete(labels=c("by total score","top 5%, bottom 95%","top 5%, bottom 50%"))+
        guides(fill=guide_legend(title=paramNames[paramCol-2]))
      if(plotAll){plot(thePlot)}
  

      #UNADJUSTED DISESE PREVALANCE!
      groupDf=data.frame(matrix(0,ncol=ncol(predDf),nrow=100))
      groupRanges=c(round(seq(1,nrow(df),by = nrow(df)/100)),nrow(df))
      repDescribe=c()
      for(i in 1:ncol(predDf)){
        subDf=df[order(df[,ncol(dfMain)+i]),c(ncol(dfMain),ncol(dfMain)+i)]
        subDf$group=rep(0,nrow(subDf))
        repDescribe=c(repDescribe,rep(workKey[i,paramCol],100))
        for(j in 1:100){
          subDf$group[groupRanges[j]:groupRanges[j+1]]=j
          groupDf[j,i]=sum(subDf[subDf$group==j,1]==1)
        }
      }
      
      groupDf$index=1:100 
      sampSize=sum(subDf$group==1)
      prevPlotDf=melt(groupDf,id.vars="index")
      prevPlotDf$value=prevPlotDf$value/sampSize
      prevPlotDf$variable=factor(repDescribe)
      thePlot=ggplot(prevPlotDf,aes(index,value,color=variable))+geom_smooth(fill="grey80")+
        labs(x="PRS percentile",y="disease prevalance",title=paste("Prevalance of",currentTrait),
             subtitle=paste("Using method",method,"- comparing",paramNames[paramCol-2]))
      if(plotAll){plot(thePlot)}
      
      topGroup=apply(groupDf[90:100,],2,sum)
      bottomGroup=apply(groupDf[1:10,],2,sum)
      quotPrev=topGroup/bottomGroup
      sePrev=sqrt((1/topGroup)+(1/bottomGroup)+(2/(10*sampSize)))
      resMat$prev=quotPrev[-length(quotPrev)]
      sdHiMat$prev=(quotPrev+exp(1.96*sePrev))[-length(quotPrev)]
      sdLowMat$prev=(quotPrev-exp(1.96*sePrev))[-length(quotPrev)]
      

      #ROC AND AUC PLOTS!!!
      totalRocDf=rbind(totalRocDf,rocDf)
      thePlot=ggplot(rocDf,aes(y=tpr,x=fpr,color=group))+geom_line()+
        geom_abline(intercept=c(0,0),slope=1)+
        labs(x="false positive rate",y="true positive rate",
             title=paste("ROC Curves for",currentTrait),
             subtitle=paste("Using method",method,"- comparing",paramNames[paramCol-2]))
      if(plotAll){plot(thePlot)}


      resMat$AUC=aucs[,2]-aucCovar[2]
      sdLowMat$AUC=aucs[,1]-aucCovar[1]
      sdHiMat$AUC=aucs[,3]-aucCovar[3]
      thePlot=ggplot(resMat,aes(index,AUC))+geom_point(aes(color=index),size=4)+
        geom_errorbar(aes(ymin=sdLowMat$AUC,ymax=sdHiMat$AUC,width=0.2),color="grey50")+
        labs(x=paramNames[paramCol-2],y="AUC improvement",title=paste("Analysis of",currentTrait),
             subtitle=paste("Using method",method,"- comparing",paramNames[paramCol-2]))+
        scale_fill_discrete(guide=FALSE)
      if(plotAll){plot(thePlot)}

      print("MADE IT TO BOTTOM")
      print(resMat) 
      totalResMat=rbind(totalResMat,cbind(workKey,resMat))
      totalHiMat=rbind(totalHiMat,cbind(workKey,sdHiMat))
      totalLoMat=rbind(totalLoMat,cbind(workKey,sdLowMat))
    }   
  }
}

print("OUT OF METHODS FLOR LOOP")


#######################################################################################################################
# Try enet regression with all scores 
print("DOING ENET REGRESSION")
print("auccovar")
print(aucCovar)

#get index of the best parameter set for each method
print("totalResMat 1")
print(totalResMat)
totalResMat=totalResMat[!duplicated(totalResMat$fName),]
print("totalResMat 2")
print(totalResMat)
bestIndices <- rep(0,length(unique(totalResMat$methods)))
i=1
for(umethod in unique(totalResMat$methods)){
  umethodAuc <- max(totalResMat[totalResMat$methods==umethod,12])
  bestIndices[i] <- which(totalResMat$AUC==umethodAuc)[1]
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

fullDfInd <- as.matrix(cbind(dfMain[,-ncol(dfMain)],allScores))
bestDfInd <- as.matrix(cbind(dfMain[,-ncol(dfMain)],allScores[,bestIndices]))
print("BESTDFIND")
print(head(bestDfInd))
print("FULLDFIND")
print(head(fullDfInd))

dep <- as.matrix(dfMain$response)
cvAns <- cv.glmnet(y=dfMain$response, x=bestDfInd, family="binomial",nfolds=3,maxit=5000)
bestLamb <- cvAns$lambda.min

#may have to implement recursion, similar to that of boostEval
enetProc <- function(df){
  theFolds=createFolds(dfMain$response,k=3)
  tpr=c(); fpr=c(); auc=c()
  for(i in 1:3){
    trainInd=df[-theFolds[[i]],]; trainDep=dfMain$response[-theFolds[[i]]]
    testInd=df[theFolds[[i]],]; testDep=dfMain$response[theFolds[[i]]]
    innerFunc <- function(trainDep, trainInd, bestLamb){
      netLog = glmnet(y=trainDep, x=trainInd,family="binomial", lambda = bestLamb, maxit = 5000)
      if(netLog[1]==0){
        toRemove <-sample(which(trainDep==0),round(length(trainDep)/4))
        trainDep <- trainDep[-toRemove]
        trainInd <- trainInd[-toRemove,]
        return(innerFunc(trainDep,trainInd,bestLamb))
      }  else {
        return(netLog)
      }
    }
    netLog <- innerFunc(trainDep, trainInd, bestLamb)
    bestPred=predict(netLog, newx = testInd)
    g=pROC::roc(testDep ~ bestPred[,1])
    
    auc=c(auc,g$auc)
    tpr=cbind(tpr,g$sensitivities[round(seq(1,length(g$sensitivities),length.out = 1000))])
    fpr=cbind(fpr,1-g$specificities[round(seq(1,length(g$specificities),length.out = 1000))])
  }
  return(list(data.frame(tpr=apply(tpr,1,mean),fpr=apply(fpr,1,mean)), mean(auc)))
}

print("enet - best")
bestEnet <- enetProc(bestDfInd)
print("enet - full")
allEnet <- enetProc(fullDfInd)


plotDf <- data.frame(rbind(bestEnet[[1]],allEnet[[1]]))
plotDf$type <- c(rep("best",nrow(bestEnet[[1]])),rep("all",nrow(allEnet[[1]])))
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

totalLoMat=totalLoMat[!(duplicated(totalLoMat$fName)),]
totalHiMat=totalHiMat[!(duplicated(totalHiMat$fName)),]
colnames(totalLoMat)[6+critIndex]="toComp"
colnames(totalHiMat)[6+critIndex]="toComp"
totalResMat=totalResMat[order(totalResMat$fName),]
totalLoMat=totalLoMat[order(totalLoMat$fName),]
totalHiMat=totalHiMat[order(totalHiMat$fName),]

if(finalCrit=="AUC"){finalCrit="AUC Improvement"}
thePlot=ggplot(totalResMat,aes(parameter,toComp))+geom_point(aes(color=methods),size=4)+
  geom_errorbar(aes(ymin=totalLoMat$toComp,ymax=totalHiMat$toComp,width=0.2), color="grey50")+
  labs(x="Parameters",y=finalCrit,title=paste("Total Analysis of",currentTrait), caption=paste("best:",round(max(totalResMat$toComp),3)))+
  guides(color=FALSE)+theme(axis.text.x = element_text(angle = 45, hjust = 1, size=5), plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
  scale_color_manual(values=cp)
plot(thePlot)

bestName=totalResMat$fName[totalResMat$toComp==max(totalResMat$toComp)][1]
bestScore=allScores[,colnames(allScores)==bestName]
goBackDf=cbind(df[,1:ncol(dfMain)],bestScore)
colnames(goBackDf)[ncol(goBackDf)]="predictor"
goBackDf=goBackDf[,c(1:(ncol(dfMain)-1),(ncol(dfMain)+1),ncol(dfMain))]
if(mean(goBackDf$predictor[goBackDf$response==1]) < mean(goBackDf$predictor[goBackDf$response==0])){
  goBackDf$predictor=goBackDf$predictor*-1
  revScore=TRUE
} else {
  revScore=FALSE
}

totalRocDf=totalRocDf[,1:2]
totalRocDf$group="train"
dev.off()

return(list(currentTrait,bestName,goBackDf,decoder,totalRocDf,totalAuc[1],revScore,totalResMat))

}


