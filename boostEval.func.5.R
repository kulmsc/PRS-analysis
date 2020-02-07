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


boostEval <- function(authorComp,bestFile,compStat,doRF,doSVM){

pdf(paste(authorComp,"boost","pdf",sep='.'))
  
#read in the decoder and make sure duplicate disease names are changed
decoder=read.table("fileDecoder",stringsAsFactors=F,header=T)
decoder=decoder[decoder[,1]==authorComp,]
decoder=unname(unlist(decoder))
currentTrait=decoder[2]

#work done to make the train and test dfs of the primary score
prepWork=function(dataSet){
  if(dataSet=="train"){
    scores=read.table(paste0("diseaseFiles/",bestFile),stringsAsFactors = F,header = T)
    pheno=read.table("pheno.phase1",stringsAsFactors = F)
  } else {
    scores=read.table(paste0("diseaseFiles/",gsub('train','test',bestFile)),stringsAsFactors = F,header = T)
    pheno=read.table("pheno.phase2",stringsAsFactors = F)
  }
  
  scores=as.data.frame(scores[!is.na(pheno[,6]),])
  pheno=pheno[!is.na(pheno[,6]),]
  
  #split up pheno into easily pulled from groups
  phenoCovar=pheno[,1:6]
  phenoAges=pheno[,7:9]
  phenoAges[is.na(phenoAges[,2]),2]=phenoAges[is.na(phenoAges[,2]),1]+mean(na.omit(phenoAges[,2]-phenoAges[,1]))
  phenoAges[is.na(phenoAges[,3]),3]=phenoAges[is.na(phenoAges[,3]),1]+mean(na.omit(phenoAges[,3]-phenoAges[,1]))
  
  splitPhenos=list()
  splitLabel=c("20004","20002","20001")
  splitPhenos[[1]]=list(pheno[,10:41],pheno[,106:133],pheno[,193:198])
  splitPhenos[[2]]=list(pheno[,42:73],pheno[,134:163],pheno[,199:204])
  splitPhenos[[3]]=list(pheno[,74:105],pheno[,164:192],pheno[,205:210])
  
  
  #for each trait/score, get the uid and code that will pull out the correct myPheno and response vector
  code=decoder[3]
  currentTrait=decoder[2]
  codeMat=sapply(strsplit(code,',')[[1]],function(x) strsplit(x,"-")[[1]])
  
  #count number of time in all reporting the disease code comes up
  predictor=scores[,which(authorComp==colnames(scores))]
  response=rep(0,length(predictor))
  takenFrom=rep(1,length(predictor))
  for(crit in 1:ncol(codeMat)){
    phenoIndex=which(splitLabel==codeMat[1,crit])
    scoreCode=as.numeric(codeMat[2,crit])
    for(g in 1:length(splitPhenos)){
      compPheno=splitPhenos[[g]][[phenoIndex]]
      instResponse=apply(compPheno,1,function(x) ifelse(scoreCode %in% x,1,0))
      takenFrom[instResponse==1]=g
      response=response+instResponse
    }
  }
  response[response>1]=1
  response=factor(response)
  ages=apply(cbind(phenoAges,takenFrom),1,function(x) x[x[4]])
  
  df=data.frame(phenoCovar,age=ages,predictor=as.numeric(predictor),response)
  #df=df[order(df$predictor,decreasing=T),]
  df$predictor=scale(df$predictor)
  return(df)
}

trainDf=prepWork("train")
testDf=prepWork("test")



#determine all the unique files names in the support files
allSuppFiles=list.files("supportFiles/",pattern = "train")
commonNames=c()
for(f in allSuppFiles){
  fOpen = file(paste0("supportFiles/",f),'r')
  firstLine=readLines(fOpen,n=1)
  close(fOpen)
  firstLine=strsplit(firstLine,split='\t')[[1]]
  commonNames = union(commonNames, firstLine)
}

#create a holder data frame for all the support files and possible scores
bestTracker=data.frame(matrix(0,ncol=length(allSuppFiles),nrow=length(commonNames)))
colnames(bestTracker)=allSuppFiles
row.names(bestTracker)=commonNames

#for each support file check out the author score for the stat specificed, recording performance in the holder made above
j=1
for(suppFile in allSuppFiles){
  supportTrain=read.table(paste0("supportFiles/",suppFile),stringsAsFactors = F,header = T)
  for(i in 1:ncol(supportTrain)){
    evalDf=cbind(trainDf,supp=scale(supportTrain[,i]))
    completeLogistic <- glm(response ~ age + sex + array + PC1 + PC2 + PC3 + PC4 + predictor + supp,data=evalDf,family="binomial")
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
    bestTracker[which(row.names(bestTracker)==colnames(supportTrain)[i]),j]=theStat 
  }  
  j=j+1
}

#take the best file for each author and made a data frame that contains the scores
bestCols=apply(bestTracker,1,function(x) which(x==max(x)))
bestTracker$bestFile=colnames(bestTracker)[bestCols]
suppTrain=data.frame(matrix(0,nrow=nrow(supportTrain),ncol=nrow(bestTracker)))
suppTest=data.frame(matrix(0,nrow=278685,ncol=nrow(bestTracker)))
colnames(suppTrain)=row.names(bestTracker)
colnames(suppTest)=row.names(bestTracker)

for(suppFile in unique(bestTracker$bestFile)){
  supportTrain=read.table(paste0("supportFiles/",suppFile),stringsAsFactors = F,header = T)
  supportTest=read.table(paste0("supportFiles/",gsub('train','test',suppFile)),stringsAsFactors = F,header = T)
  
  suppTrain[,bestTracker$bestFile == suppFile] = supportTrain[,which(colnames(supportTrain) %in% rownames(bestTracker)[bestTracker$bestFile==suppFile])]
  suppTest[,bestTracker$bestFile == suppFile] = supportTest[,which(colnames(supportTrain) %in% rownames(bestTracker)[bestTracker$bestFile==suppFile])]
}
suppTrain=apply(suppTrain,2,scale)
suppTest=apply(suppTest,2,scale)

#do the decoder stuff
allSupportDecode=read.table("supportDecoder",stringsAsFactors = F,header=T)
supportDecode=allSupportDecode[allSupportDecode$author %in% colnames(suppTrain),]


#get the baseline stat for plotting the relative change of adding the covariates
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
plotDf$trait=factor(plotDf$trait,levels = plotDf$trait[order(plotDf$stat)] )
plotDf$stat=plotDf$stat-theStat
thePlot=ggplot(plotDf,aes(trait,stat))+geom_col()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x="Trait",y=paste("Difference in",compStat),title=paste("Support Improvement on",decoder[2],"Prediction"))
plot(thePlot)  

#Now we have all the best support Features according to either or, prev, or auc - can now do feature selection
evalDf=cbind(trainDf,suppTrain)
evalDf$array=as.numeric(as.factor(evalDf$array))
evalDf$sex=as.numeric(as.factor(evalDf$sex))


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
forPlot=data.frame(traits=as.character(rownames(netLog$beta)),beta=as.numeric(netLog$beta))
forPlot$traits=as.character(forPlot$traits)
forPlot=forPlot[forPlot$beta!=0,]
forPlot[forPlot$traits %in% allSupportDecode$author,1]=allSupportDecode[allSupportDecode$author %in% forPlot$traits,2]
forPlot$traits=factor(forPlot$traits,levels=forPlot$traits[order(forPlot$beta)])
thePlot=ggplot(forPlot,aes(traits,beta))+geom_col()+
  labs(title="Elastic Net Logistic Regreesion Feature Analysis",subtitle=paste("For Trait",currentTrait))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot(thePlot)

forPlot=data.frame(fpr=apply(fpr,1,mean),tpr=apply(tpr,1,mean))
thePlot=ggplot(forPlot,aes(fpr,tpr))+geom_point()+
  labs(title="Regularized Logistic Regression Performance - Training",
       subtitle=paste("For Trait",currentTrait),
       caption=paste("AUC of",as.character(round(mean(auc),4))),
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
    print(i)
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
  print(ncol(df))
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
    print(c(ans,ncol(nextDf)))
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
testEvalDf=cbind(testDf,supportTest)
testEvalDf$array=as.numeric(as.factor(testEvalDf$array))
testEvalDf$sex=as.numeric(as.factor(testEvalDf$sex))

allAuc=c(enetAuc,rfAuc,svmAuc)
if(max(allAuc)==enetAuc){
  #LASSO Logistic Regression
  netLog = glmnet(y=evalDf$response,x=as.matrix(evalDf[,-9]),family="binomial",lambda = netCv$lambda.min)
  bestPred=predict(netLog, newx = as.matrix(testEvalDf[,-9]))
  g=roc(testEvalDf$response ~ bestPred[,1])
  auc=c(auc,g$auc)
  
  forPlot=cbind(tpr=g$sensitivities,fpr=1-g$specificities)
  forPlot=as.data.frame(forPlot[seq(1,nrow(forPlot),100),])
  thePlot=ggplot(forPlot,aes(fpr,tpr))+geom_point()+
    labs(title="Regularized Logistic Regression Performance - Testing",
         subtitle=paste("For Trait",currentTrait),
         caption=paste("AUC of",as.character(round(g$auc,4))),
         x="true positive rate",y="false positive rate")+
    geom_abline(intercept=c(0,0),slope=1)
  plot(thePlot)
  
  betas=as.matrix(netLog$beta)
  featsToRemove=unname(which(betas[,1]!=0))
  featRemoveMat=matrix(0,nrow=length(bestPred),ncol=length(featsToRemove)+1)
  featRemoveMat[,1]=bestPred
  i=2
  for(feat in featsToRemove){
    netLog = glmnet(y=evalDf$response,x=as.matrix(evalDf[,-c(feat,9)]),family="binomial",lambda = netCv$lambda.min)
    bestPred=predict(netLog, newx = as.matrix(testEvalDf[,-c(feat,9)]))
    featRemoveMat[,i]=bestPred
    i=i+1
  }
  newNames=rownames(betas)[featsToRemove]
  newNames[newNames %in% supportDecode$author]=supportDecode[supportDecode$author %in% newNames,2]
  colnames(featRemoveMat)=c("bestPred",newNames)
  featRemoveMat=featRemoveMat[seq(1,nrow(featRemoveMat),100),]
  featRemoveMat=as.data.frame(featRemoveMat)
  forPlot=melt(featRemoveMat,id.vars = "bestPred")
  thePlot=ggplot(forPlot,aes(bestPred,value))+geom_point(aes(color=variable),alpha=0.1)+
    labs(x="full model score",y="one feature removed score",title="Effect of Features on Individual Scores")+
    guides(color = guide_legend(override.aes = list(alpha=1)))
  plot(thePlot)
  
  
  
  
  
} else if(max(allAuc)==rfAuc){
  rfRes=randomForest(response ~ .,data=evalDf,ntree=300,mtry=unlist(tuneRes$best.parameters))
  rfPred=predict(rfRes,test,type="prob")[,2]
  
  
  
  
} else if(max(allAuc)==svmAuc){
  svmRes=svm(response~., evalDf,kernel="radial",probability=T,gamma=bestGam,cost=bestCo)
  svmPred=predict(svmRes,testEvalDf,probability = T)
  svmPred=attr(svmPred,"prob")[,2]
}

dev.off()

return(bestTracker)
}
