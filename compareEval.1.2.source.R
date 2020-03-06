validation <- function(df, theType, folds, caseInBoth = 800, controlMultiplier = 5){
  responseNum <- which(colnames(df)=="response")
  caseInds <- sample(which(df[,responseNum]==1))
  controlInds <- sample(which(df[,responseNum]==0))
  
  caseSplits <- ceiling(seq_along(caseInds) / ceiling(length(caseInds)/folds))
  controlSplits <- ceiling(seq_along(controlInds) / ceiling(length(controlInds)/folds))
  if(caseInBoth > length(caseInds)){caseInBoth <- length(caseInds)}
  if(length(controlInds) < controlMultiplier*(caseInBoth/2)){controlMultiplier <- length(controlInds)/(caseInBoth/2)}
  print("in validation") 
  print(caseInBoth)
  print(controlMultiplier)
 
  dfTrainList <- rep(list(0),folds)
  dfTestList <- rep(list(0),folds)
  
  for(i in 1:folds){
    if(theType=="cross"){
      dfTrainList[[i]] <- df[c(caseInds[caseSplits==i],controlInds[controlSplits==i]),]
      dfTestList[[i]] <- df[c(caseInds[caseSplits!=i],controlInds[controlSplits!=i]),]
    } else if(theType=="bag"){
      trainCaseSample <- sample(caseInds,round(caseInBoth/2))
      trainControlSample <- sample(controlInds,round((caseInBoth/2)*controlMultiplier))
      testCaseSample <- sample(setdiff(caseInds,trainCaseSample),round(caseInBoth/2))
      testControlSample <- sample(setdiff(controlInds,trainControlSample),round(caseInBoth/2)*controlMultiplier)
      dfTrainList[[i]] <- df[c(trainCaseSample, trainControlSample),]
      dfTestList[[i]] <- df[c(testCaseSample,testControlSample),]
    } 
  }
  if(all(c(sapply(dfTrainList, function(x) length(unique(x[,1]))),sapply(dfTestList, function(x) length(unique(x[,1]))))>1)){
    return(list(dfTrainList,dfTestList))
  } else {
    print("bad return")
    return(validation(df, theType, folds, caseInBoth, controlMultiplier))
  }
}

confintChecker <- function(theMod){

  out <- tryCatch(
     {
       compLogitCI <- confint(theMod, level=0.95)
       message("this is try")
       return(compLogitCI)
     },
     error=function(cond) {
       message("this is the error")
       theCoef <- coef(summary(theMod))
       theCoef[,2] <- theCoef[,1] * 1.5
       theCoef[,1] <- theCoef[,1] * 0.5
       theCoef <- theCoef[,1:2]
       return(theCoef)
     },
     #warning=function(cond) {
     #  message("this is warning")
     #  return(NULL)
     #},
     finally={
       message("final message")
     }
  )
  return(out)

}


makeResMat <- function(predDf, dfMain, workKey, paramCol, paramNames, currentTrait, method, plotAll, runType = "normal"){
  #collect the pval and OR by logistic regression
  folds <- 3
  df <- data.frame(dfMain,predDf)
  resMat=rep(list(data.frame(matrix(0,nrow=ncol(predDf),ncol=5,dimnames=list(NULL,c("index","pval","or","or2","or3"))))),folds)
  #resMat=lapply(resMat,function(x){ x[,1]<-factor(workKey[,paramCol]); x})
  sdLowMat=rep(list(data.frame(resMat[[1]])),folds)
  sdHiMat=rep(list(data.frame(resMat[[1]])),folds)
  aucs=rep(list(matrix(0,nrow=ncol(predDf),ncol=3)),folds)
  rocDf=rep(list(matrix(0,nrow=1000*ncol(predDf),ncol=3)),folds)
  aucCovar=rep(list(0),folds)
  
  #change from the current code to cross validation

  
  dfList <- validation(df, "cross", folds)
  
  print("STARTING FOLDS")
  for(iFold in 1:folds){
    for(i in 1:ncol(predDf)){

      #normal logistic regression of score against response
      dfTrain <- dfList[[1]][[iFold]]
      dfTest <- dfList[[2]][[iFold]]
      dfMainTrain <- dfTrain[,grepl(paste(c("array","age","sex","PC1","PC2","PC3","PC4","response","covarResp"),collapse="|"),colnames(dfTrain)),drop=F]
      dfMainTest <- dfTest[,grepl(paste(c("array","age","sex","PC1","PC2","PC3","PC4","response","covarResp"),collapse="|"),colnames(dfTest)),drop=F]
      predDfTrain <- dfTrain[,!grepl(paste(c("array","age","sex","PC1","PC2","PC3","PC4","response","covarResp"),collapse="|"),colnames(dfTrain)),drop=F]
      predDfTest <- dfTest[,!grepl(paste(c("array","age","sex","PC1","PC2","PC3","PC4","response","covarResp"),collapse="|"),colnames(dfTest)),drop=F]

      
      completeLogistic <- glm(response ~ . ,data=cbind(dfMainTrain,pred=predDfTrain[,i]),family="binomial")
      #compLogitCI=confint(completeLogistic,level=0.95)
      compLogitCI <- confintChecker(completeLogistic)
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
      dfTest$prob=predict(completeLogistic, cbind(dfMainTest,pred=predDfTest[,i]), type=c("response"))
      g=pROC::roc(response ~ prob,data=dfTest, quiet=TRUE)
      aucs[[iFold]][i,]=as.numeric(ci.auc(g$auc, progress = "none"))
      rocDfTemp=cbind(g$sensitivities,1-g$specificities,rep(workKey[i,paramCol],length(g$specificities)))
      rocDf[[iFold]][(i*1000-999):(i*1000),]=rocDfTemp[seq(1,nrow(rocDfTemp),length.out = 1000),]


      if(covarLogistic$converged){
        dfMainTest$prob=predict(covarLogistic,dfMainTest,type=c("response"))
        g=pROC::roc(response ~ prob,data=dfMainTest, quiet = TRUE)
        aucCovar[[iFold]]=as.numeric(ci.auc(g, progress = "none"))
      } else {
        #aucCovar will be the completeLogistic AUC, not good, although will make AUC Improvement 0
        print("AUC Covar did not converge")
        break
        aucCovar[[iFold]]=g$auc
      }
      #rocDf[[iFold]]$group=factor(rocDf[[iFold]]$group)
    }
  }  
  
  resMat=Reduce("+", resMat) / folds
  sdLowMat=Reduce("+", sdLowMat) / folds
  sdHiMat=Reduce("+", sdHiMat) / folds
  
  if(runType=="normal"){
    resMat[,1]=factor(workKey[,paramCol])
    sdLowMat[,1]=factor(workKey[,paramCol])
    sdHiMat[,1]=factor(workKey[,paramCol])
  } else {
    resMat[,1] <- factor(c(0,1))
    sdLowMat[,1] <- factor(c(0,1))
    sdHiMat[,1] <- factor(c(0,1))
  }
  
  aucs=Reduce("+", aucs) / folds
  rocDf=data.frame(Reduce("+", rocDf) / folds)
  colnames(rocDf) <- c("tpr","fpr","group")
  rocDf$group=factor(rocDf$group)
  aucCovar=Reduce("+", aucCovar) / folds
  
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
  #write.table(repDescribe,"222.print")
  prevPlotDf$variable=factor(repDescribe)
  thePlot=ggplot(prevPlotDf,aes(index,value,color=variable))+geom_smooth(fill="grey80")+
    labs(x="PRS percentile",y="disease prevalance",title=paste("Prevalance of",currentTrait),
         subtitle=paste("Using method",method,"- comparing",paramNames[paramCol-2]))
  if(plotAll){plot(thePlot)}
  
  #topGroup=apply(groupDf[90:100,],2,sum)
  #bottomGroup=apply(groupDf[1:10,],2,sum)
  #quotPrev=topGroup/bottomGroup
  #sePrev=sqrt((1/topGroup)+(1/bottomGroup)+(2/(10*sampSize)))
  dfMainTest <- df[,grepl(paste(c("array","age","sex","PC1","PC2","PC3","PC4","response","covarResp"),collapse="|"),colnames(df)),drop=F]
  predDfTest <- df[,!grepl(paste(c("array","age","sex","PC1","PC2","PC3","PC4","response","covarResp"),collapse="|"),colnames(df)),drop=F]
  for(i in 1:ncol(predDfTest)){
    simpleLogistic <- glm(response ~ pred ,data=cbind(dfMainTest,pred=predDfTest[,i]),family="binomial")
    resMat$prev[i] <- NagelkerkeR2(simpleLogistic)$R2
    sdHiMat$prev[i] <- NagelkerkeR2(simpleLogistic)$R2*1.1
    sdLowMat$prev[i] <- NagelkerkeR2(simpleLogistic)$R2*0.9
  }
  #resMat$prev=quotPrev[-length(quotPrev)]
  #sdHiMat$prev=(quotPrev+exp(1.96*sePrev))[-length(quotPrev)]
  #sdLowMat$prev=(quotPrev-exp(1.96*sePrev))[-length(quotPrev)]
  

  #write.table("246.print","246")  
  #ROC AND AUC PLOTS!!!
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
  #write.table("265.print","265")
  return(list(resMat,sdHiMat,sdLowMat, rocDf))
}




enetProc <- function(df, dfMain, bestLamb, goNet){
  numFolds <- 3
  df <- cbind(df,response=as.numeric(dfMain$response)-1)
  dfList <- validation(df,"cross",numFolds)
  print("done validation")

  #theFolds=createFolds(dfMain$response,k=numFolds)
  compiledPrs <- matrix(0,nrow=nrow(df),ncol=numFolds)
  tpr=c(); fpr=c(); auc=c()
  allScoreBeta <- matrix(0,nrow=numFolds,ncol=length(grep("disease",colnames(df))))
  for(i in 1:numFolds){
    print(i)
    #trainInd=df[-theFolds[[i]],]; trainDep=dfMain$response[-theFolds[[i]]]
    #testInd=df[theFolds[[i]],]; testDep=dfMain$response[theFolds[[i]]]
    trainInd <- dfList[[1]][[i]][,-ncol(df)]
    testInd <- dfList[[2]][[i]][,-ncol(df)]
    trainDep <- dfList[[1]][[i]][,ncol(df)]
    testDep <- dfList[[2]][[i]][,ncol(df)]

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

    if(goNet){
      netLog <- innerFunc(trainDep, trainInd, bestLamb)
      scoreBetas <- as.matrix(netLog$beta)[grep("disease",colnames(df)),]
    } else {
      netLog = glm(trainDep ~ ., data=data.frame(trainDep,trainInd),family="binomial")
      scoreBetas <- netLog$coefficients[2:length(netLog$coefficients)]
      scoreBetas <- as.matrix(scoreBetas)[grep("disease",names(scoreBetas)),]
    }

    compiledPrs[,i] <- rowSums(t(t(df[,grep("disease",colnames(df))])*scoreBetas))
    allScoreBeta[i,] <- scoreBetas
    trainScore <- rowSums(t(t(trainInd[,grep("disease",colnames(trainInd))])*scoreBetas))
    testScore <- rowSums(t(t(testInd[,grep("disease",colnames(testInd))])*scoreBetas))
    trainScore <- (trainScore-mean(trainScore))/sd(trainScore)
    testScore <- (testScore-mean(testScore))/sd(testScore)
    trainInd <- data.frame(trainInd[,grep("disease",colnames(trainInd),invert=T)], score=trainScore, resp=as.factor(trainDep))
    testInd <- data.frame(testInd[,grep("disease",colnames(testInd),invert=T)], score=testScore, resp=as.factor(testDep))
    
    theModel <- glm(resp ~ ., data=trainInd, family="binomial")
    bestPred=predict(theModel, testInd, type="response")
    g=pROC::roc(testInd$resp ~ bestPred, quiet = TRUE)
    
    auc=c(auc,g$auc)
    tpr=cbind(tpr,g$sensitivities[round(seq(1,length(g$sensitivities),length.out = 1000))])
    fpr=cbind(fpr,1-g$specificities[round(seq(1,length(g$specificities),length.out = 1000))])
  }

  print("done folds")
  returnScoreBetas <- apply(allScoreBeta, 2, function(x) weighted.mean(x,auc))
  names(returnScoreBetas) <- names(scoreBetas)
  finalScore <- rowSums(t(t(df[,grep("disease",colnames(df))])*returnScoreBetas))
  finalScore <- (finalScore-mean(finalScore))/sd(finalScore)
  print("returning")
  return(list(data.frame(tpr=apply(tpr,1,mean),fpr=apply(fpr,1,mean)), mean(auc), finalScore, returnScoreBetas))
}
