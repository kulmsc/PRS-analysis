library(fmsb)
library(PRROC)
library(ggplot2)
library(pROC)
library(reshape2)
library(glmnet)
library(Boruta)
library(randomForest)
library(caret)
library(e1071)
library(ROCR)
library(parallel)

boostEval <- function(authorComp,bestFile,compStat,doRF,doSVM,trainDf,testDf,allSupportDecode,allDiseaseDecoder,
                      trainROC,trainAUC,testROC,testAUC){

pdf(paste(authorComp,"boost","pdf",sep='.'))
  

#determine all the unique files names in the support files
allSuppFiles=list.files("supportFiles/",pattern = "train")
allDiseaseFiles=list.files("diseaseFiles/",pattern="train")
allFiles=c(allSuppFiles,allDiseaseFiles)
names(allFiles)=c(rep("supportFiles",length(allSuppFiles)),rep("diseaseFiles",length(allDiseaseFiles)))

commonNames=c()
for(i in 1:length(allFiles)){
  f=allFiles[i]
  fOpen = file(paste0(names(f),"/",f),'r')
  firstLine=readLines(fOpen,n=1)
  close(fOpen)
  firstLine=strsplit(firstLine,split='\t')[[1]]
  commonNames = union(commonNames, firstLine)
}
commonNames=commonNames[-which(commonNames==authorComp)]

#create a holder data frame for all the support files and possible scores
bestTracker=data.frame(matrix(0,ncol=length(allFiles),nrow=length(commonNames)))
colnames(bestTracker)=allFiles
row.names(bestTracker)=sort(commonNames)

#for each support file check out the author score for the stat specificed, recording performance in the holder made above
funky <- function(pullCol,trainDf,supportTrain,compStat,roc){
  evalDf=cbind(trainDf,supp=scale(supportTrain[,pullCol]))
  completeLogistic <- glm(response ~ age + sex + array + PC1 + PC2 + PC3 + PC4 + predictor * supp,data=evalDf,family="binomial")
  evalDf$pred=predict(completeLogistic,evalDf,type="response")
  evalDf=evalDf[order(evalDf$pred,decreasing = T),]
  
  if(compStat=="OR"){
    exposeGroup=evalDf[1:round(nrow(evalDf)*0.05),9]
    safeGroup=evalDf[(round(nrow(evalDf)*0.05)+1):nrow(evalDf),9]
    theStat=(sum(exposeGroup==1)/sum(exposeGroup==0))/(sum(safeGroup==1)/sum(safeGroup==0))
  } else if (compStat=="Prevalance") {
    theStat = sum(head(evalDf,round(nrow(evalDf)/10))$response==1)/sum(tail(evalDf,round(nrow(evalDf)/10))$response==1)
  } else {
    g=roc(response ~ pred,data=evalDf)
    theStat=as.numeric(g$auc)
  }
  return(theStat) 
} 

cl <- makeCluster(2)
for(j in 1:length(allFiles)){
  theFile=allFiles[j]
  if(names(theFile)=="supportFiles"){
    supportTrain=supportTrainList[[which(allSuppFiles==theFile)]]
  } else {
    supportTrain=diseaseTrainList[[which(allDiseaseFiles==theFile)]]
  }
  if(authorComp %in% colnames(supportTrain)){
    supportTrain=supportTrain[,-which(colnames(supportTrain)==authorComp)]
  }
  if(sum(colSums(supportTrain)==0)>0){
    supportTrain=supportTrain[,-which(colSums(supportTrain)==0)]
  }
  toChange=grepl(".",colnames(supportTrain),fixed=T)
  colnames(supportTrain)[toChange]=gsub(".","-",colnames(supportTrain)[toChange],fixed=T)
  theStats=parLapply(cl,1:ncol(supportTrain),funky,trainDf,supportTrain,compStat,roc)
  theStats=unlist(theStats)
  bestTracker[which(row.names(bestTracker) %in% colnames(supportTrain)),j]=theStats 
}

#take the best file for each author and made a data frame that contains the scores
bestCols=apply(bestTracker,1,function(x) which(x==max(x)))
bestTracker$bestFile=colnames(bestTracker)[bestCols]

#This should remove the predictors that obviously are not predictive
scaledBest=scale(apply(bestTracker[,-ncol(bestTracker)],1,max))[,1]
if(any(scaledBest< -1)){
  bestTracker=bestTracker[-which(scaledBest< -1),]
}
if(length(scaledBest)>20){
  bestTracker=bestTracker[-which(scaledBest < sort(scaledBest,decreasing = T)[20]),]
}
  
suppTrain=data.frame(matrix(0,nrow=nrow(supportTrain),ncol=nrow(bestTracker)))
suppTest=data.frame(matrix(0,nrow=278685,ncol=nrow(bestTracker)))
colnames(suppTrain)=row.names(bestTracker)
colnames(suppTest)=row.names(bestTracker)


#reading in the test support scores and distributing out the train scores which are already read in
for(suppFile in unique(bestTracker$bestFile)){
  if(str_split(suppFile,pattern=fixed("."))[[1]][1]=="support"){
    supportTrain=supportTrainList[[which(allSuppFiles==suppFile)]]
    supportTest=supportTestList[[which(allSuppFiles==suppFile)]]
  } else {
    supportTrain=diseaseTrainList[[which(allDiseaseFiles==suppFile)]]
    supportTest=diseaseTestList[[which(allDiseaseFiles==suppFile)]]
  }
  toChange=grepl(".",colnames(supportTrain),fixed=T)
  colnames(supportTrain)[toChange]=gsub(".","-",colnames(supportTrain)[toChange],fixed=T)
  toChange=grepl(".",colnames(supportTest),fixed=T)
  colnames(supportTest)[toChange]=gsub(".","-",colnames(supportTest)[toChange],fixed=T)
  
  suppTrain[,bestTracker$bestFile == suppFile] = supportTrain[,which(colnames(supportTrain) %in% rownames(bestTracker)[bestTracker$bestFile==suppFile])]
  suppTest[,bestTracker$bestFile == suppFile] = supportTest[,which(colnames(supportTest) %in% rownames(bestTracker)[bestTracker$bestFile==suppFile])]
}
suppTrain=apply(suppTrain,2,scale)
suppTest=apply(suppTest,2,scale)

#do the decoder stuff
supportDecode=allSupportDecode[allSupportDecode$author %in% colnames(suppTrain),]
temp=allDiseaseDecoder[allDiseaseDecoder$author %in% colnames(suppTrain),1:2]
colnames(temp)=colnames(supportDecode)[1:2]
supportDecode=rbind(supportDecode[,1:2],temp)
supportDecode=supportDecode[order(supportDecode$author),]

#get the baseline stat for plotting the relative change of adding the covariates
evalDf=data.frame(trainDf)
completeLogistic <- glm(response ~ age + sex + array + PC1 + PC2 + PC3 + PC4 + predictor,data=evalDf,family="binomial")
evalDf$pred=predict(completeLogistic,evalDf,type="response")
evalDf=evalDf[order(evalDf$pred,decreasing = T),]

if(compStat=="OR"){
  exposeGroup=evalDf[1:round(nrow(evalDf)*0.05),9]
  safeGroup=evalDf[(round(nrow(evalDf)*0.05)+1):nrow(evalDf),9]
  theStat=(sum(exposeGroup==1)/sum(exposeGroup==0))/(sum(safeGroup==1)/sum(safeGroup==0))
} else if (compStat=="Prevalance") {
  theStat = sum(head(evalDf,round(nrow(evalDf)/10))$response==1)/sum(tail(evalDf,round(nrow(evalDf)/10))$response==1)
} else {
  g=roc(response ~ pred,data=evalDf)
  theStat=as.numeric(g$auc)
}

#now make the plot of all supports from the best files by the stat it was analyzed by
plotDf=data.frame(stat=apply(bestTracker[,-ncol(bestTracker)],1,max),trait=factor(supportDecode$trait))
plotDf$trait=factor(plotDf$trait,levels = plotDf$trait[order(plotDf$stat)])
plotDf$stat=plotDf$stat-theStat
thePlot=ggplot(plotDf,aes(trait,stat))+geom_col()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
  labs(x="Trait",y=paste("Difference in",compStat),title=paste("Support Improvement on",decoder[2],"Prediction"))
plot(thePlot)  


#Now we have all the best support Features according to either or, prev, or auc - can now do feature selection
evalDf=cbind(trainDf,suppTrain)
evalDf$array=as.numeric(as.factor(evalDf$array))
evalDf$sex=as.numeric(as.factor(evalDf$sex))
holdOutEvalDf=evalDf

###################################################################################################################
#Feature selection through logistic lasso regression ##############################################################
#the model making
netCv = cv.glmnet(y=evalDf$response,x=as.matrix(evalDf[,-9]),family="binomial",nfolds = 3)
theFolds=createFolds(evalDf[,9],k=3)
tpr=c(); fpr=c(); auc=c()
for(i in 1:3){
  train=evalDf[-theFolds[[i]],]
  test=evalDf[theFolds[[i]],]
  netLog = glmnet(y=train$response,x=as.matrix(train[,-9]),family="binomial",lambda = netCv$lambda.min)
  bestPred=predict(netLog, newx = as.matrix(test[,-9]))
  g=roc(test$response ~ bestPred[,1])
  tpr=cbind(tpr,g$sensitivities[seq(1,nrow(test),length.out = 1000)])
  fpr=cbind(fpr,1-g$specificities[seq(1,nrow(test),length.out = 1000)])
  auc=c(auc,g$auc)
}

#the plotting
#The betas of each feature
forPlot=data.frame(traits=as.character(rownames(netLog$beta)),beta=as.numeric(netLog$beta))
forPlot$traits=as.character(forPlot$traits)
forPlot=forPlot[forPlot$beta!=0,]
forPlot[forPlot$traits %in% supportDecode$author,1]=supportDecode[supportDecode$author %in% forPlot$traits,2]
forPlot$traits=factor(forPlot$traits,levels=forPlot$traits[order(forPlot$beta)])
thePlot=ggplot(forPlot,aes(traits,beta))+geom_col()+
  labs(title="Elastic Net Logistic Regreesion Feature Analysis",subtitle=paste("For Trait",currentTrait))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))
plot(thePlot)

#The ROC plot
forPlot=data.frame(fpr=apply(fpr,1,mean),tpr=apply(tpr,1,mean))
trainROC[,2]=rev(trainROC[,2])
forPlot=rbind(cbind(forPlot,type="enet"),cbind(trainROC[,c(2,1)],type="logreg"))
forPlot$type=factor(forPlot$type)
thePlot=ggplot(forPlot,aes(fpr,tpr,color=type))+geom_point()+
  labs(title="Regularized Logistic Regression Performance - Training",
       subtitle=paste("For Trait",currentTrait),
       caption=paste("AUC of",as.character(round(mean(auc),4)),"an improvement of",as.character(round(auc-trainAUC,4))),
       y="true positive rate",x="false positive rate")+
  geom_abline(intercept = c(0,0),slope = 1)
plot(thePlot)
enetAuc=mean(auc)


##########################################################################################################################
#Feature selection with boruta ############################################################################################
#the Feature Selection
if(doRF){
casesInds=which(evalDf$response==1)
controlInds=which(evalDf$response==0)
controlInds=controlInds[sort(sample(1:length(controlInds),length(casesInds)*5))]
goDf=evalDf[sort(c(casesInds,controlInds)),]
borutaFeats <- Boruta(response ~ ., data=goDf,doTrace=2,ntree=300)
borutaHist <- borutaFeats$ImpHistory[,1:(ncol(trainDf)+ncol(suppTrain)-1)]
impRating=apply(borutaHist,2,function(x) max(which(x!=-Inf)))     #Want more than the features they give me
impRating=sort(impRating,decreasing = T)                          #So will expand to nearly approved features
finalImp=names(impRating[impRating>(impRating[1]/2)])
treeDf=evalDf[,which(colnames(evalDf) %in% c(finalImp,"response"))]
respCol=which(colnames(treeDf)=="response")

#the model making
tuneRes=tune(randomForest,response~.,data=treeDf,ranges=list(mtry=2:3),ntree=100,folds=3)
auc=c()
perfList=list()
for(i in 1:3){
  train=treeDf[-theFolds[[i]],]
  test=treeDf[theFolds[[i]],]
  rfRes=randomForest(response ~ .,data=train,ntree=1000,mtry=unlist(tuneRes$best.parameters))
  rfPred=predict(rfRes,test,type="prob")[,2]
  rfPredon = prediction(rfPred, test$response)
  rfPerf = performance(rfPredon,"tpr","fpr")
  rfAuc = performance(rfPredon,"auc")
  auc=c(auc,rfAuc@y.values[[1]])
  perfList[[i]]=list(rfPerf@x.values[[1]],rfPerf@y.values[[1]])
}

allTpr=c(perfList[[1]][[1]],perfList[[2]][[1]],perfList[[3]][[1]])
allFpr=c(perfList[[1]][[2]],perfList[[2]][[2]],perfList[[3]][[2]])
newFpr=loess(allFpr ~ allTpr)$fitted

#the plotting
impBorutaHist=borutaFeats$ImpHistory[,colnames(borutaFeats$ImpHistory) %in% finalImp]
rfFeatRanks=apply(impBorutaHist,2,function(x) x[max(which(x!=-Inf))])
forPlot=data.frame(traits=names(rfFeatRanks),beta=rfFeatRanks)
forPlot$traits=as.character(forPlot$traits)
forPlot[forPlot$traits %in% allSupportDecode$author,1]=allSupportDecode[allSupportDecode$author %in% forPlot$traits,2]
forPlot$traits=factor(forPlot$traits,levels=forPlot$traits[order(forPlot$beta)])
ggplot(forPlot,aes(traits,beta))+geom_col()+
  labs(title="Random Forest Feature Analysis",subtitle=paste("For Trait",currentTrait))

forPlot=data.frame(tpr=newFpr,fpr=allTpr)
ggplot(forPlot,aes(fpr,tpr))+geom_point()+
  geom_abline(intercept=c(0,0),slope=1)+
  labs(title="Random Forest Performance",
       subtitle=paste("For Trait",currentTrait),
       caption=paste("AUC of",as.character(mean(auc))),
       x="true positive rate",y="false positive rate")
rfAuc=mean(auc)

} else {
rfAuc=0
}

##############################################################################################################################
#Feature selection with rfe ##################################################################################################
#The functions

if(doSVM){
myCV <- function(df,guessGamma=0.4,guessCost=2,folds=3,returnROC=FALSE){
  auc=c(); tpr=c(); fpr=c()
  theFolds=createFolds(df[,1],k=folds)
  for(i in 1:folds){
    train=df[-theFolds[[i]],]
    test=df[theFolds[[i]],]
    svmRes=svm(response~.,train,kernel="radial",probability=T,gamma=guessGamma,cost=guessCost)
    svmPred=predict(svmRes,test,probability = T)
    svmPred=attr(svmPred,"prob")[,2]
    g=roc(test$response ~ svmPred)
    auc=c(auc,g$auc)
    if(returnROC){
      tpr=cbind(tpr,g$sensitivities[round(seq(1,length(g$sensitivities),length.out = 1000))])
      fpr=cbind(fpr,1-g$specificities[round(seq(1,length(g$specificities),length.out = 1000))])
    }
  }
  if(returnROC){
    return(list(auc,tpr,fpr))
  } else {
    return(mean(auc))
  }
}

myRFE <- function(df,varToTest){
  prof=c()
  for(var in varToTest){
    auc=myCV(df[,-var])
    names(auc)=colnames(df)[var]
    prof=c(prof,auc)
  }
  worstVar=varToTest[prof==max(prof)]
  nextDf=df[,-worstVar]
  if(ncol(nextDf)>=3){
    if(length(varToTest)==ncol(nextDf)){
      newProf=prof[-which(prof==max(prof))]
      nextNames=names(newProf)[newProf %in% sort(newProf,decreasing = T)[1:4]]
      nextVars=which(colnames(nextDf) %in% nextNames)
    } else {
      nextVars=which(colnames(nextDf)!="response")
    }
    ans=myRFE(nextDf,nextVars)
  } else {
    ans=NULL
  }
  return(c(ans,prof[prof==max(prof)]))
}

bayesWrapper <- function(gamma,cost){
  auc=myCV(df=svmDf,guessGamma=gamma,guessCost=cost)
  return(list(Score=auc,Pred=0))
}
  

#feature selection
rfeRes <- myRFE(goDf,which(colnames(goDf) != "response"))
survFeat=colnames(goDf[!(colnames(goDf) %in% c(names(rfeRes),"response"))])
otherFeats=names(rfeRes)[1:which(rfeRes==max(rfeRes[1:6]))]
svmDf=goDf[,colnames(goDf) %in% c(survFeat,otherFeats,"response")]
svmDf=evalDf[,colnames(goDf) %in% c(survFeat,otherFeats,"response")]

#casesInds=which(tempDf$response==1)
#controlInds=which(tempDf$response==0)
#controlInds=controlInds[sort(sample(1:length(controlInds),length(casesInds)*12))]
#svmDf3=tempDf[sort(c(casesInds,controlInds)),]
#BayesianOptimization(bayesWrapper, bounds = list(gamma = c(0, 1),cost=c(1,3)), n_iter = 1, init_points = 3)

#making the model
bestAuc=0
for(gam in seq(0.1,0.9,0.2)^2){
  for(co in seq(0.5,3,0.75)){
    auc=myCV(svmDf,gam,co)
    if(auc>bestAuc){
      bestAuc=auc; bestGam=gam; bestCo=co
    }
  }
}
allROC=myCV(svmDf,bestGam,bestCo,returnROC = TRUE)

#plotting
forPlot=data.frame(fpr=apply(allROC[[3]],1,mean),tpr=apply(allROC[[2]],1,mean))
ggplot(forPlot,aes(tpr,fpr))+geom_point()+
  labs(title="SVM Performance",
       subtitle=paste("For Trait",currentTrait),
       caption=paste("AUC of",as.character(mean(allROC[[1]]))),
       x="true positive rate",y="false positive rate")

finalFeatAuc=myCV(goDf[,which(colnames(goDf) %in% c("response",survFeat))])
baselineAuc=myCV(goDf[,which(colnames(goDf) %in% c("response","predictor"))])
rfeRes=c(survFeat=finalFeatAuc,rfeRes)
names(rfeRes)[1]=survFeat
forPlot=data.frame(traits=names(rfeRes),rfeRes)
forPlot$traits=as.character(forPlot$traits)
forPlot[forPlot$traits %in% allSupportDecode$author,1]=allSupportDecode[allSupportDecode$author %in% forPlot$traits,2]
forPlot$rfeRes <- forPlot$rfeRes-baselineAuc
forPlot$traits=factor(forPlot$traits,levels=forPlot$traits)
ggplot(forPlot,aes(traits,rfeRes))+geom_col()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x="Trait",y="AUC change on Predictor Only Model",title="SVM RFE Ranking",
       subtitle="full model on right with each trait cumulatively removed going to the left")
svmAuc=mean(allROC[[1]])
} else {
svmAuc=0
}

######################################################################################################################
######################################################################################################################
######################################################################################################################
testEvalDf=cbind(testDf,suppTest)
testEvalDf$array=as.numeric(as.factor(testEvalDf$array))
testEvalDf$sex=as.numeric(as.factor(testEvalDf$sex))
evalDf=holdOutEvalDf

allAuc=c(enetAuc,rfAuc,svmAuc)
if(max(allAuc)==enetAuc){
  #LASSO Logistic Regression
  netLog = glmnet(y=evalDf$response,x=as.matrix(evalDf[,-9]),family="binomial",lambda = netCv$lambda.min)
  bestPred=predict(netLog, newx = as.matrix(testEvalDf[,-9]))
  g=roc(testEvalDf$response ~ bestPred[,1])
  auc=c(auc,g$auc)
  
  boostPlot=cbind(tpr=g$sensitivities,fpr=1-g$specificities)
  boostPlot=as.data.frame(boostPlot[seq(1,nrow(boostPlot),length.out = 1001),])
  testROC=testROC[testROC$group=="PRS+covars",1:2]
  testROC[,2]=rev(testROC[,2])
  forPlot=rbind(cbind(boostPlot,type="enet"),cbind(testROC,type="logreg"))
  
  thePlot=ggplot(forPlot,aes(fpr,tpr,color=type))+geom_point()+
    labs(title="Regularized Logistic Regression Performance - Testing",
         subtitle=paste("For Trait",currentTrait),
         caption=paste("AUC of",as.character(round(g$auc,4)),"an improvement of",as.character(round(g$auc-testAUC,4))),
         x="true positive rate",y="false positive rate")+
    geom_abline(intercept=c(0,0),slope=1)
  plot(thePlot)
  
  
  # betas=as.matrix(netLog$beta)
  # featsToRemove=unname(which(betas[,1]!=0))
  # featRemoveMat=matrix(0,nrow=length(bestPred),ncol=length(featsToRemove)+1)
  # featRemoveMat[,1]=bestPred
  # i=2
  # for(feat in featsToRemove){
  #   netLog = glmnet(y=evalDf$response,x=as.matrix(evalDf[,-c(feat,9)]),family="binomial",lambda = netCv$lambda.min)
  #   bestPred=predict(netLog, newx = as.matrix(testEvalDf[,-c(feat,9)]))
  #   featRemoveMat[,i]=bestPred
  #   i=i+1
  # }
  # newNames=rownames(betas)[featsToRemove]
  # newNames[newNames %in% supportDecode$author]=supportDecode[supportDecode$author %in% newNames,2]
  # colnames(featRemoveMat)=c("bestPred",newNames)
  # featRemoveMat=as.data.frame(featRemoveMat)
  # featRemoveMat$response=testEvalDf$response
  # 
  # featRemoveMat[,ncol(featRemoveMat)]=as.numeric(featRemoveMat[,ncol(featRemoveMat)])
  # prevHolder=matrix(0,nrow=100,ncol=ncol(featRemoveMat))
  # colnames(prevHolder)=colnames(featRemoveMat)
  # prevHolder[,ncol(prevHolder)]=1:100
  # for(i in 1:(ncol(prevHolder)-1)){
  #   perc=quantile(featRemoveMat[,i],seq(0,100,1)/100)
  #   for(j in 1:100){
  #     prevHolder[j,i]=sum(featRemoveMat[featRemoveMat[,i]>perc[j] & featRemoveMat[,i]<perc[j+1],ncol(featRemoveMat)])
  #   }
  # }
  # prevHolder=as.data.frame(prevHolder)
  # topSum=apply(prevHolder[90:100,],2,sum)
  # keepNames=names(sort(topSum)[c(1:4,(ncol(prevHolder)-3):ncol(prevHolder))])
  # prevHolderPlot=prevHolder[,which(colnames(prevHolder) %in% keepNames)]
  # toPlot=melt(prevHolderPlot,id.vars = "response")
  # thePlot=ggplot(toPlot,aes(response,value,color=variable))+geom_smooth(se=F)+
  #   labs(x="Score Percentile",y="Total Cases",title="Effect of Removing the Score From the Model")+
  #   scale_color_discrete(name = "Score Removed")
  # plot(thePlot)
  
  
  
  
} else if(max(allAuc)==rfAuc){
  rfRes=randomForest(response ~ .,data=evalDf,ntree=300,mtry=unlist(tuneRes$best.parameters))
  rfPred=predict(rfRes,test,type="prob")[,2]
  
  
  
  
} else if(max(allAuc)==svmAuc){
  svmRes=svm(response~., evalDf,kernel="radial",probability=T,gamma=bestGam,cost=bestCo)
  svmPred=predict(svmRes,testEvalDf,probability = T)
  svmPred=attr(svmPred,"prob")[,2]
}

dev.off()

return(list(bestTracker,suppTrain,suppTest,boostPlot,g$auc))
}
