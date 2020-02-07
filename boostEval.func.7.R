library(fmsb)
library(PRROC)
library(ggplot2)
library(pROC)
library(reshape2)
library(glmnet)
library(caret)
library(ROCR)
library(parallel)
library(ggrepel)

boostEval <- function(authorComp,bestFile,compStat,trainDf,testDf,allSupportDecode,allDiseaseDecode,
                      trainROC,trainAUC,testROC,testAUC,maxBest){

pdf(paste(authorComp,"boost","pdf",sep='.'))
  

#determine all the unique files names in the support files
allSuppFiles=list.files("supportFiles/",pattern = "train")
allDiseaseFiles=list.files("diseaseFiles/",pattern="train")
allFiles=c(allSuppFiles,allDiseaseFiles)
names(allFiles)=c(rep("supportFiles",length(allSuppFiles)),rep("diseaseFiles",length(allDiseaseFiles)))

#Get all of the authors
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

#check if the author has data in both a test and train file
checkStatusTrain=rep(0,length(commonNames))
for(checkDf in c(supportTrainList,diseaseTrainList)){
    checkStatusTrain[commonNames %in% names(checkDf)[apply(checkDf,2,sum)>0]]=1
}
checkStatusTest=rep(0,length(commonNames))
for(checkDf in c(supportTestList,diseaseTestList)){
  checkStatusTest[commonNames %in% names(checkDf)[apply(checkDf,2,sum)>0]]=1
}
badNames=commonNames[checkStatusTrain==0 | checkStatusTest==0]
if(sum(checkStatusTrain==0 | checkStatusTrain==0)>0){
  commonNames=commonNames[-which(commonNames %in% badNames)]
}

##########################################################################################################################
# Make best tracker ######################################################################################################

#create a holder data frame for all the support files and possible scores
bestTracker=data.frame(matrix(0,ncol=length(allFiles),nrow=length(commonNames)))
colnames(bestTracker)=allFiles
row.names(bestTracker)=sort(commonNames)

#for each support file check out the author score for the stat specificed, recording performance in the holder made above
parallelBT <- function(pullCol,trainDf,supportTrain,compStat,roc){
  evalDf=cbind(trainDf,supp=scale(supportTrain[,pullCol]))
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
  return(theStat) 
} 

cl <- makeCluster(1)
for(j in 1:length(allFiles)){
  theFile=allFiles[j]
  if(names(theFile)=="supportFiles"){
    supportTrain=supportTrainList[[which(allSuppFiles==theFile)]]
    supportTest=supportTestList[[which(allSuppFiles==theFile)]]
  } else {
    supportTrain=diseaseTrainList[[which(allDiseaseFiles==theFile)]]
    supportTest=diseaseTestList[[which(allDiseaseFiles==theFile)]]
  }
  supportTrain=supportTrain[,colnames(supportTrain) %in% colnames(supportTest)]
  supportTest=supportTest[,colnames(supportTest) %in% colnames(supportTest)]
  if(sum(colSums(supportTrain)==0)>0){
    supportTrain=supportTrain[,-which(colSums(supportTrain)==0 | colSums(supportTest)==0)]
  }
  if(authorComp %in% colnames(supportTrain)){
    supportTrain=supportTrain[,-which(colnames(supportTrain)==authorComp)]
  }
  if(sum(colnames(supportTrain) %in% badNames)>0){
    supportTrain=supportTrain[,-which(colnames(supportTrain) %in% badNames)]
  }
  
  toChange=grepl(".",colnames(supportTrain),fixed=T)
  colnames(supportTrain)[toChange]=gsub(".","-",colnames(supportTrain)[toChange],fixed=T)
  theStats=parLapply(cl,1:ncol(supportTrain),parallelBT,trainDf,supportTrain,compStat,roc)
  theStats=unlist(theStats)
  bestTracker[which(row.names(bestTracker) %in% colnames(supportTrain)),j]=theStats 
}



#take the best file for each author and made a data frame that contains the scores
bestCols=apply(bestTracker,1,function(x) which(x==max(x))[1])
bestTracker$bestFile=colnames(bestTracker)[bestCols]
sendBackBestTracker=data.frame(bestTracker)

#This should remove the predictors that obviously are not predictive
scaledBest=apply(bestTracker[,-ncol(bestTracker)],1,max)
if(length(scaledBest)>maxBest){
  bestTracker=bestTracker[-which(scaledBest < sort(scaledBest,decreasing = T)[maxBest]),]
}

print("a")
####################################################################################################################
#construct the larger evalDf #######################################################################################

suppTrain=data.frame(matrix(0,nrow=129568,ncol=nrow(bestTracker)))
suppTest=data.frame(matrix(0,nrow=278685,ncol=nrow(bestTracker)))
colnames(suppTrain)=row.names(bestTracker)
colnames(suppTest)=row.names(bestTracker)

print("b")
#reading in the test support scores and distributing out the train scores which are already read in
for(suppFile in unique(bestTracker$bestFile)){
  print(suppFile)
  if(str_split(suppFile,pattern=fixed("."))[[1]][1]=="support"){
    supportTrain=supportTrainList[[which(allSuppFiles==suppFile)]]
    supportTest=supportTestList[[which(allSuppFiles==suppFile)]]
    print("q")
  } else {
    supportTrain=diseaseTrainList[[which(allDiseaseFiles==suppFile)]]
    supportTest=diseaseTestList[[which(allDiseaseFiles==suppFile)]]
    print("w")
  }
  toChange=grepl(".",colnames(supportTrain),fixed=T)
  colnames(supportTrain)[toChange]=gsub(".","-",colnames(supportTrain)[toChange],fixed=T)
  print("e")
  toChange=grepl(".",colnames(supportTest),fixed=T)
  colnames(supportTest)[toChange]=gsub(".","-",colnames(supportTest)[toChange],fixed=T)
  
  print("r")
  suppTrain[,bestTracker$bestFile == suppFile] = supportTrain[,which(colnames(supportTrain) %in% rownames(bestTracker)[bestTracker$bestFile==suppFile])]
  print("t")
  print(dim(suppTest))
  suppTest[,bestTracker$bestFile == suppFile] = supportTest[,which(colnames(supportTest) %in% rownames(bestTracker)[bestTracker$bestFile==suppFile])]
  print(dim(suppTest))
  print("y")
}
suppTrain=apply(suppTrain,2,scale)
suppTest=apply(suppTest,2,scale)

print("c")
#do the decoder stuff - get matrix going from author to disease name
supportDecode=allSupportDecode[allSupportDecode$author %in% colnames(suppTrain),]
temp=allDiseaseDecode[allDiseaseDecode$author %in% colnames(suppTrain),1:2]
colnames(temp)=colnames(supportDecode)[1:2]
supportDecode=rbind(supportDecode[,1:2],temp)
supportDecode=supportDecode[order(supportDecode$author),]

#####################################################################################################################
#make the plot of the individual logreg features ####################################################################

#get the baseline stat for plotting the relative change of adding the covariates
evalDf=data.frame(trainDf)
completeLogistic <- glm(response ~ age + sex + array + PC1 + PC2 + PC3 + PC4 + predictor,data=evalDf,family="binomial")
evalDf$pred=predict(completeLogistic,evalDf,type="response")
evalDf=evalDf[order(evalDf$pred,decreasing = T),]


print("d")
if(compStat=="OR"){
  exposeGroup=evalDf[1:round(nrow(evalDf)*0.05),9]
  safeGroup=evalDf[(round(nrow(evalDf)*0.05)+1):nrow(evalDf),9]
  theStat=(sum(exposeGroup==1)/sum(exposeGroup==0))/(sum(safeGroup==1)/sum(safeGroup==0))
} else if (compStat=="Prevalance") {
  theStat = sum(head(evalDf,round(nrow(evalDf)/10))$response==1)/sum(tail(evalDf,round(nrow(evalDf)/10))$response==1)
} else {
  g=roc(response ~ pred,data=evalDf)
  theStat=as.numeric(g$auc)
  trainROC=data.frame(tpr=g$sensitivities,fpr=1-g$specificities,group="train")
  trainROC=trainROC[seq(1,nrow(trainROC),length.out = 1000),]
}

save(bestTracker,file="bestTracker")
save(supportDecode,file="supportDecode")
save(allSupportDecode,file="allSupportDecode")
save(allDiseaseDecode,file="allDiseaseDecode")
save(suppTrain,file="suppTrain")

#now make the plot of all supports from the best files by the stat it was analyzed by
plotDf=data.frame(stat=apply(bestTracker[,-ncol(bestTracker)],1,max),trait=factor(supportDecode$trait))
plotDf$trait=factor(plotDf$trait,levels = plotDf$trait[order(plotDf$stat)])
plotDf$stat=plotDf$stat-theStat
plotDf=plotDf[1:10,]
thePlot=ggplot(plotDf,aes(trait,stat))+geom_col()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
  labs(x="Unrelated Score",y=paste("Difference in",compStat),
       title="Improvement by adding unrelated score to model",
       subtitle=paste(decoder[2],"Model"))
plot(thePlot)  

print("e")
################################################################################################################################
# prove in depth analysis of each score's interaction with disease #############################################################

parallelA<-function(pullCol,allTrain,allTest,trainDf,testDf,roc){
  #pval
  print("a")
  trainDf$predictor=allTrain[,pullCol]
  testDf$predictor=allTest[,pullCol]
  completeLogistic <- glm(response ~ . ,data=trainDf,family="binomial")
  covarLogistic <- glm(response ~ age + sex + array + PC1 + PC2 + PC3 + PC4,data=trainDf,family="binomial")
  funcPval= -log10(coef(summary(completeLogistic))[9,4])
  
  #AUC
  print("b")
  testDf$prob=predict(completeLogistic,testDf,type="response")
  testDf$base=predict(covarLogistic,testDf,type="response")
  gScore=roc(response ~ prob,data=testDf)
  gCovar=roc(response ~ base,data=testDf)
  funcAuc=gScore$auc-gCovar$auc
  
  #OR
  print("c")
  testDf=testDf[order(testDf$predictor,decreasing = T),]
  allCutOffs=c(0.5,0.2,0.1,0.05,0.01,0.005)
  i=1
  funcOR=rep(0,length(allCutOffs))
  for(cutOff in allCutOffs){
    splitPoint=round(nrow(testDf)*cutOff)
    exposeGroup=testDf[1:splitPoint,9]
    safeGroup=testDf[(splitPoint+1):nrow(testDf),9]
    funcOR[i]=(sum(exposeGroup==1)/sum(exposeGroup==0))/(sum(safeGroup==1)/sum(safeGroup==0))
    i=i+1
  }
  
  #Prev
  print("d")
  funcPrev=rep(0,100)
  testDf$group=rep(0,nrow(testDf))
  groupRanges=c(round(seq(1,nrow(testDf),by = nrow(testDf)/100)),nrow(testDf))
  for(i in 1:100){
    testDf$group[groupRanges[i]:groupRanges[i+1]]=i
  }
  testDf$group=rev(testDf$group)
  for(i in 1:100){
    withDis=sum(testDf[testDf$group==i,9]==1)
    totalGroup=sum(testDf$group==i)
    funcPrev[i]=withDis/totalGroup
  }
  
  return(list(funcPval,funcAuc,funcOR,funcPrev))
}

makePlotsColor <- function(allRes,trainScores,translateLabels,totalTraits,numberColor){
  colnames(trainScores)=translateLabels[colnames(trainScores) %in% translateLabels[,1],2]
  
  pvals=unlist(lapply(allRes, `[[`, 1))
  aucImps=unlist(lapply(allRes, `[[`, 2))
  ors=t(matrix(unlist(lapply(allRes, `[[`, 3)),nrow=6))
  prevs=t(matrix(unlist(lapply(allRes, `[[`, 4)),nrow=100))
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  specColors=c(gg_color_hue(numberColor),rep("#808080",totalTraits-numberColor))
  
  #AUC
  fullAucImpsDf=data.frame(aucImps,traits=colnames(trainScores),stringsAsFactors = F)
  aucImpsDf=fullAucImpsDf[order(fullAucImpsDf[,1],decreasing = T),]
  aucImpsDf=aucImpsDf[1:totalTraits,]
  names(specColors)=aucImpsDf$traits
  aucImpsDf[,2]=factor(aucImpsDf[,2],levels = aucImpsDf[order(aucImpsDf[,1]),2])
  thePlot=ggplot(aucImpsDf,aes(traits,aucImps,fill=traits))+geom_col()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
    labs(x="Scores",y="AUC Improvement",title="AUC Improvement",caption=paste("Models exclude",currentTrait,"score"))+
    scale_fill_manual(values=specColors)+guides(fill=FALSE)
  plot(thePlot)
  
  #PVALS
  fullPvalsDf=data.frame(pvals,traits=colnames(trainScores),stringsAsFactors = F)
  pvalsDf=fullPvalsDf[fullPvalsDf$traits %in% as.character(aucImpsDf$traits),]
  pvalsDf=pvalsDf[1:totalTraits,]
  pvalsDf[,2]=factor(pvalsDf[,2],levels = pvalsDf[order(pvalsDf[,1]),2])
  thePlot=ggplot(pvalsDf,aes(traits,pvals,fill=traits))+geom_col()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
    labs(x="Scores",y="logistic regression PVal",title=paste("Logistic Regression PVals of PRS"),
         caption=paste("Models exclude",currentTrait,"score"))+
    scale_fill_manual(values=specColors)+guides(fill=FALSE)
  plot(thePlot)
  
  #OR vs PVALS
  fullOrPvalDf=data.frame(pvals,ors=ors[,4],traits=colnames(trainScores),stringsAsFactors = F)
  orPvalDf=data.frame(fullOrPvalDf)
  orPvalDf$word=orPvalDf$traits
  orPvalDf$word[orPvalDf$pvals < sort(orPvalDf[,1],decreasing = T)[3] & orPvalDf$ors < sort(orPvalDf[,2],decreasing = T)[3]]=""
  thePlot=ggplot(orPvalDf,aes(pvals,ors,label=word))+
    geom_point(data=orPvalDf[orPvalDf$traits %in% names(specColors[1:3]),],aes(color=traits),size=4)+
    geom_point(data=orPvalDf[!(orPvalDf$traits %in% names(specColors[1:3])),])+
    geom_label_repel(size=4,force=3)+
    labs(x="PValue",y="Odds Ratio",title="Effect and Certainty",caption=paste("Models exclude",currentTrait,"score"))+
    theme(legend.position="bottom",legend.direction="vertical")+
    scale_color_manual(values=specColors[1:3])
  plot(thePlot)
  
  #OR
  colnames(ors)=paste("cut0ff",c(0.5,0.2,0.1,0.05,0.01,0.005),sep=".")
  fullOrsDf=data.frame(ors,traits=colnames(trainScores),stringsAsFactors = F)
  orsDf=fullOrsDf[fullOrsDf$traits %in% as.character(aucImpsDf$traits),]
  orsDf=melt(orsDf,id.vars="traits")
  orsDf[,2]=sort(rep(c(0.5,0.2,0.1,0.05,0.01,0.005),totalTraits),decreasing = T)
  orsDf[,2]=factor(orsDf[,2],levels=c(0.5,0.2,0.1,0.05,0.01,0.005))
  thePlot=ggplot(orsDf,aes(variable,value,color=traits))+geom_line(aes(group=traits))+
    labs(x="Percentile Cut-Off",y="Odds Ratio",title="Odds ratios",caption=paste("Models exclude",currentTrait,"score"))+
    scale_color_manual(values=specColors)+guides(color=FALSE)
  plot(thePlot)
  
  #PREV
  colnames(prevs)=1:100/100
  fullPrevsDf=data.frame(prevs,traits=colnames(trainScores),stringsAsFactors = F)
  prevsDf=fullPrevsDf[fullPrevsDf$traits %in% as.character(aucImpsDf$traits),]
  prevsDf=melt(prevsDf,id.vars = "traits")
  prevsDf[,2]=sort(rep((1:100/100),totalTraits))
  thePlot=ggplot(prevsDf,aes(variable,value,color=traits))+geom_smooth(se=F)+
    labs(x="Score Percentile",y="Trait Prevalence",title="Disease Prevalences",
         caption=paste("Models exclude",currentTrait,"score"))+
    scale_color_manual(values=specColors)+guides(color=FALSE)
  plot(thePlot)
}

makePlotsBest <- function(allRes,trainScores,translateLabels){
  colnames(trainScores)=translateLabels[colnames(trainScores) %in% translateLabels[,1],2]
  
  pvals=unlist(lapply(allRes, `[[`, 1))
  aucImps=unlist(lapply(allRes, `[[`, 2))
  ors=t(matrix(unlist(lapply(allRes, `[[`, 3)),nrow=6))
  prevs=t(matrix(unlist(lapply(allRes, `[[`, 4)),nrow=100))
  
  #PVALS
  fullPvalsDf=data.frame(pvals,traits=colnames(trainScores),stringsAsFactors = F)
  pvalsDf=fullPvalsDf[order(fullPvalsDf[,1],decreasing = T),]
  pvalsDf=pvalsDf[1:8,]
  pvalsDf[,2]=factor(pvalsDf[,2],levels = pvalsDf[order(pvalsDf[,1]),2])
  thePlot=ggplot(pvalsDf,aes(traits,pvals))+geom_col(aes(fill=traits))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
    labs(x="traits",y="logistic regression PVal",title=paste("Pval of isolated scores on",currentTrait))+
    guides(fill=FALSE)
  plot(thePlot)
  
  #OR vs PVALS
  fullOrPvalDf=data.frame(pvals,ors=ors[,4],traits=colnames(trainScores),stringsAsFactors = F)
  orPvalDf=data.frame(fullOrPvalDf)
  orPvalDf$traits[orPvalDf$pvals < sort(orPvalDf[,1],decreasing = T)[3] & orPvalDf$ors < sort(orPvalDf[,2],decreasing = T)[3]]=""
  thePlot=ggplot(orPvalDf,aes(pvals,ors,label=traits))+geom_point()+geom_label_repel(size=3,force=3)+
    labs(x="PValue",y="Odds Ratio",title=paste("Effect and Certainty for isolated scores on",currentTrait))
  plot(thePlot)
  
  #AUC
  fullAucImpsDf=data.frame(aucImps,traits=colnames(trainScores),stringsAsFactors = F)
  aucImpsDf=fullAucImpsDf[order(fullAucImpsDf[,1],decreasing = T),]
  aucImpsDf=aucImpsDf[1:8,]
  aucImpsDf[,2]=factor(aucImpsDf[,2],levels = aucImpsDf[order(aucImpsDf[,1]),2])
  thePlot=ggplot(aucImpsDf,aes(traits,aucImps))+geom_col(aes(fill=traits))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
    labs(x="Traits",y="AUC Improvement",title=paste("AUC Improvement of isolated scores on",currentTrait))+
    guides(fill=FALSE)
  plot(thePlot)
  
  #OR
  colnames(ors)=paste("cut0ff",c(0.5,0.2,0.1,0.05,0.01,0.005),sep=".")
  fullOrsDf=data.frame(ors,traits=colnames(trainScores),stringsAsFactors = F)
  orsDf=fullOrsDf[order(rowSums(fullOrsDf[,1:6]),decreasing = T)[1:8],]
  orsDf=melt(orsDf,id.vars="traits")
  orsDf[,2]=sort(rep(c(0.5,0.2,0.1,0.05,0.01,0.005),8),decreasing = T)
  orsDf[,2]=factor(orsDf[,2],levels=c(0.5,0.2,0.1,0.05,0.01,0.005))
  thePlot=ggplot(orsDf,aes(variable,value,color=traits))+geom_line(aes(group=traits))+
    labs(x="Percentile Cut-Off",y="Odds Ratio",title=paste("Odds ratios of isolated scores on",currentTrait))
  plot(thePlot)
  
  #PREV
  colnames(prevs)=1:100/100
  fullPrevsDf=data.frame(prevs,traits=colnames(trainScores),stringsAsFactors = F)
  prevsDf=fullPrevsDf[order(rowSums(fullPrevsDf[,91:100])/rowSums(fullPrevsDf[,1:8]),decreasing = T)[1:8],]
  prevsDf=melt(prevsDf,id.vars = "traits")
  prevsDf[,2]=sort(rep((1:100/100),8))
  thePlot=ggplot(prevsDf,aes(variable,value,color=traits))+geom_smooth(se=F)+
    labs(x="Score Percentile",y="Trait Prevalence",title=paste("Disease Prevalences for isolated scores on",currentTrait))
  plot(thePlot)
}

cl <- makeCluster(2)
allResReturn=parLapply(cl,1:ncol(suppTrain),parallelA,suppTrain,suppTest,trainDf,testDf,roc)
makePlotsColor(allResReturn,suppTrain,supportDecode,8,3)

################################################################################################################################
###############################################################################################################################

print("h")
#Now we have all the best support Features according to either or, prev, or auc - can now do feature selection
evalDf=cbind(trainDf,suppTrain)
evalDf$array=as.numeric(as.factor(evalDf$array))
evalDf$sex=as.numeric(as.factor(evalDf$sex))


######################################################################################################################
# split each support feature into deciles ############################################################################

makeSplit <- function(suppDf){
  splitTrain=data.frame(matrix(0,nrow=nrow(suppDf),ncol=ncol(suppDf)*3))
  colnames(splitTrain)=sort(paste0(rep(colnames(suppDf),3),c("_q1","_q2","_q3")))
  j=1
  for(i in 1:ncol(suppDf)){
    splitTrain[,j:(j+2)]=suppDf[,i]
    splitTrain[,j][splitTrain[,j] > quantile(suppDf[,i],0.25)] = mean(suppDf[,i][suppDf[,i] > quantile(suppDf[,i],0.25)])
    splitTrain[,j+1][splitTrain[,j+1] < quantile(suppDf[,i],0.25)] = mean(suppDf[,i][suppDf[,i] < quantile(suppDf[,i],0.25)])
    splitTrain[,j+1][splitTrain[,j+1] > quantile(suppDf[,i],0.75)] = mean(suppDf[,i][suppDf[,i] > quantile(suppDf[,i],0.75)])
    splitTrain[,j+2][splitTrain[,j+2] < quantile(suppDf[,i],0.75)] = mean(suppDf[,i][suppDf[,i] < quantile(suppDf[,i],0.75)])
    j=j+3
  }
  return(splitTrain)
}
splitTrain=makeSplit(suppTrain)
splitEvalDf=cbind(trainDf,splitTrain)
splitEvalDf$array=as.numeric(as.factor(splitEvalDf$array))
splitEvalDf$sex=as.numeric(as.factor(splitEvalDf$sex))

###################################################################################################################
#Feature selection through logistic lasso regression ##############################################################
###################################################################################################################

###################################################################################################################
#enetRegression function ##########################################################################################

enetRegression <- function(evalDf,bestLamb=NULL,justAuc=F,checkMissing=F,possIndex=1,makePlot=T,addTitle=""){
  if(is.null(bestLamb)){
    netCv = try(cv.glmnet(y=evalDf$response,x=as.matrix(evalDf[,-9]),family="binomial",nfolds = 3,maxit=5000),outFile = "test")
    if(class(netCv)=="try-error"){
      if(justAuc==T){
        return(0.5)
      } else {
        return(rep(0.5,nrow(evalDf)))
      }
    } else {
      bestLamb=netCv$lambda.min 
    }
  }
  evalDf$pred=rep(0,nrow(evalDf))
  theFolds=createFolds(evalDf[,9],k=3)
  tpr=c(); fpr=c(); auc=c()
  for(i in 1:3){
    train=evalDf[-theFolds[[i]],]
    test=evalDf[theFolds[[i]],]
    if(checkMissing){
      checkedTrain=apply(train,2,function(z) length(unique(z)))
      checkedTest=apply(test,2,function(z) length(unique(z)))
      evalChecked = checkedTest==1 | checkedTrain==1
      if(any(evalChecked)){
        train=train[,-which(evalChecked)]
        test=test[,-which(evalChecked)]
      }
    }
    respCol=which(colnames(train)=="response")
    netLog = glmnet(y=train$response,x=as.matrix(train[,-respCol]),family="binomial",lambda = netCv$lambda.min,maxit = 1000)
    bestPred=predict(netLog, newx = as.matrix(test[,-respCol]))
    g=roc(test$response ~ bestPred[,1])
    print(g$auc)
    print(length(g$sensitivities))
    
    if(g$auc==0.5){
      print("was 0.5")
      g$specificities=seq(1,0,length.out = 1000)
      g$sensitivities=seq(0,1,length.out = 1000)
    } else if (length(g$sensitivities)<1000){
      print("was short")
      g$sensitivities=approx(x=1:length(g$sensitivities),y=g$sensitivities,n=1000)$y
      g$specificities=approx(x=1:length(g$specificities),y=g$specificities,n=1000)$y
    } 
    tpr=cbind(tpr,g$sensitivities[round(seq(1,length(g$sensitivities),length.out = 1000))])
    fpr=cbind(fpr,1-g$specificities[round(seq(1,length(g$specificities),length.out = 1000))])
    
    auc=c(auc,g$auc)
  }
  print("will return")
  print(mean(auc))
  if(justAuc){return(mean(auc))}
  rocPlot=data.frame(fpr=apply(fpr,1,mean),tpr=apply(tpr,1,mean))
  
  #the plotting of betas of each feature
  if(makePlot){
    forPlot=data.frame(traits=as.character(rownames(netLog$beta)),beta=as.numeric(netLog$beta))
    forPlot$traits=as.character(forPlot$traits)
    forPlot[forPlot$traits %in% supportDecode$author,1]=supportDecode[supportDecode$author %in% forPlot$traits,2]
    dupTraits=forPlot$traits[duplicated(forPlot$traits)]
    forPlot$traits[duplicated(forPlot$traits)]=paste0(dupTraits,1:length(dupTraits))
    forPlot$traits=factor(forPlot$traits,levels=forPlot$traits[order(forPlot$beta)])
    cutOff=sort(abs(forPlot$beta),decreasing = T)[15]
    forPlot=forPlot[abs(forPlot$beta) > cutOff,]
    
    print(forPlot)
    thePlot=ggplot(forPlot,aes(traits,beta))+geom_col()+
      labs(title=paste("Elastic Net Logistic Regression -",addTitle),subtitle=paste("For Trait",currentTrait))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))
    plot(thePlot)
  }
    
  return(list(as.numeric(g$auc),rocPlot))
}

######################################################################################################################
#analyze the supports and split supports #############################################################################

#The ROC plot
evalReturn=enetRegression(evalDf,addTitle = "Full Scores")
splitEvalReturn=enetRegression(splitEvalDf,addTitle = "Split Scores")

#produce roc plot to compare splitting the supports
evalReturn[[2]]$group="noSplit"
splitEvalReturn[[2]]$group="split"
allAuc=c("noSplit"=evalReturn[[1]],"split"=splitEvalReturn[[1]])

forPlot=rbind(trainROC,evalReturn[[2]],splitEvalReturn[[2]])
forPlot$group=factor(forPlot$group)
thePlot=ggplot(forPlot,aes(fpr,tpr,color=group))+geom_point()+
  labs(title="Regularized Logistic Regression Performance - Training",
       subtitle=paste("For Trait",currentTrait),
       caption=paste("AUCs: norm =",round(theStat,3),"~ scores = ",round(allAuc[1],3),"~ split = ",round(allAuc[2],3)),
       y="true positive rate",x="false positive rate")+
  geom_abline(intercept = c(0,0),slope = 1)
plot(thePlot)


######################################################################################################################
### apply better of split vs no split to the testing data ############################################################ 
if(names(allAuc)[allAuc==max(allAuc)]=="noSplit"){
  testEvalDf=cbind(testDf,suppTest)
} else {
  evalDf=splitEvalDf
  splitTest=makeSplit(suppTest)
  testEvalDf=cbind(testDf,splitTest)
}
testEvalDf$array=as.numeric(as.factor(testEvalDf$array))
testEvalDf$sex=as.numeric(as.factor(testEvalDf$sex))


#LASSO Logistic Regression
netCv = cv.glmnet(y=evalDf$response,x=as.matrix(evalDf[,-9]),family="binomial",nfolds = 3)
netLog = glmnet(y=evalDf$response,x=as.matrix(evalDf[,-9]),family="binomial",lambda = netCv$lambda.min)
bestPred=predict(netLog, newx = as.matrix(testEvalDf[,-9]))
g=roc(testEvalDf$response ~ bestPred[,1])
auc=c(auc,g$auc)

firstBoostAuc=g$auc
boostPlot=cbind(tpr=g$sensitivities,fpr=1-g$specificities)
boostPlot=as.data.frame(boostPlot[seq(1,nrow(boostPlot),length.out = 1001),])
testROCSmall=testROC[testROC$group=="PRS+covars",1:2]
forPlot=rbind(cbind(boostPlot,type="enet"),cbind(testROCSmall,type="logreg"))

thePlot=ggplot(forPlot,aes(fpr,tpr,color=type))+geom_point()+
  labs(title="Regularized Logistic Regression Performance - Testing",
       subtitle=paste("For Trait",currentTrait),
       caption=paste("AUCs: enet =",round(g$auc,3)," ~ logreg =",round(testAUC,4)),
       x="true positive rate",y="false positive rate")+
  geom_abline(intercept=c(0,0),slope=1)
plot(thePlot)
bestScoreAUC=g$auc
  

######################################################################################################################
######################################################################################################################
#Disease status analysis #############################################################################################

#####################################################################################################################
#create the data frame and analyze each feature separately ##########################################################
evalDf=cbind(trainDf,suppTrain)
evalDf$array=as.numeric(as.factor(evalDf$array))
evalDf$sex=as.numeric(as.factor(evalDf$sex))

sumScoresTrain=apply(allScoresTrain,2,sum)
allScoresTrain=allScoresTrain[,sumScoresTrain>500]
allScoresTest=allScoresTest[,sumScoresTrain>500]

sumSupportsTrain=apply(allSupportsTrain,2,sum)
allSupportsTrain=allSupportsTrain[,sumSupportsTrain>500]
allSupportsTest=allSupportsTest[,sumSupportsTrain>500]

totalPhenoTrain=cbind(allScoresTrain,allSupportsTrain)
totalPhenoTest=cbind(allScoresTest,allSupportsTest)

parallelD <- function(pullCol,trainDf,supportTrain,compStat,roc){
  evalDf=cbind(trainDf,supp=factor(supportTrain[,pullCol]))
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
  return(theStat) 
} 

cl <- makeCluster(1)
scoresStats=parLapply(cl,1:ncol(totalPhenoTrain),parallelD,trainDf,totalPhenoTrain,compStat,roc)
scoresStats=unlist(scoresStats)
scoreStatsReturn=as.vector(scoresStats)

if(length(scoresStats)>maxBest){
  totalPhenoTrain=totalPhenoTrain[,scoresStats>sort(scoresStats,decreasing = T)[maxBest] & scoresStats!=1]
  totalPhenoTest=totalPhenoTest[,scoresStats>sort(scoresStats,decreasing = T)[maxBest] & scoresStats!=1]
}

scoresStatsKeep=scoresStats[scoresStats>sort(scoresStats,decreasing = T)[maxBest] & scoresStats!=1]
fullEvalDf=cbind(evalDf,totalPhenoTrain)
fullEvalDfTest=cbind(testDf,suppTest,totalPhenoTest)
fullEvalDfTest$array=as.numeric(as.factor(fullEvalDfTest$array))
fullEvalDfTest$sex=as.numeric(as.factor(fullEvalDfTest$sex))


forPlot=data.frame(cbind(traits=colnames(totalPhenoTrain),scores=scoresStatsKeep),stringsAsFactors = F)
forPlot[,2]=as.numeric(forPlot[,2])-theStat
forPlot=forPlot[forPlot$scores>sort(forPlot$scores,decreasing = T)[15],]
forPlot[,1]=factor(forPlot[,1],levels=forPlot[order(forPlot[,2]),1])
thePlot=ggplot(forPlot,aes(traits,scores))+geom_col()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
  labs(x="Unrelated Trait",y="AUC Improvement",title=paste("Improvement by adding unrelated trait to model"),
       subtitle=paste("For trait:",currentTrait))
plot(thePlot)

####################################################################################################################
#apply decision stump on each disease feature ######################################################################
print("391")
#the full model
fullRes=enetRegression(fullEvalDf,addTitle = "With All Traits")

#the stump model
decisionStump <- function(splitByCol,x,enetRegression){
  print(splitByCol)
  library(glmnet)
  library(caret)
  library(pROC)
  splitBy=x[,splitByCol]
  x=x[,-splitByCol]
  x1=x[splitBy==0,]
  x2=x[splitBy==1,]
  
  auc1=enetRegression(x1,justAuc = T,checkMissing = T,makePlot = F)
  print("auc1")
  print(auc1)
  auc2=enetRegression(x2, justAuc = T,checkMissing = T,makePlot = F)
  print("auc2")
  print(auc2)
  return((auc1*nrow(x1)+auc2*nrow(x2))/nrow(x))
}

stumpStats=parLapply(cl,(ncol(evalDf)+1):ncol(fullEvalDf),decisionStump,fullEvalDf,enetRegression)
stumpStats=unlist(stumpStats)
compNames1=colnames(totalPhenoTrain)
splitByCol1=which(colnames(fullEvalDf)==compNames1[which(stumpStats==max(stumpStats))])
splitter1=fullEvalDf[,splitByCol1]
nextEvalDf=fullEvalDf[splitter1==0,-splitByCol1]
compNames2=colnames(nextEvalDf)[(ncol(evalDf)+1):ncol(nextEvalDf)]
print("416")

stumpStats2=parLapply(cl,(ncol(evalDf)+1):ncol(nextEvalDf),decisionStump,nextEvalDf,enetRegression)
stumpStats2=unlist(stumpStats2)
splitByCol2=which(colnames(nextEvalDf)==compNames2[which(stumpStats2==max(stumpStats2))])
splitter2=nextEvalDf[,splitByCol2]
lastEvalDf=nextEvalDf[splitter2==0,-splitByCol2]
compNames3=colnames(lastEvalDf)[(ncol(evalDf)+1):ncol(lastEvalDf)]

print("422")
stumpStats3=parLapply(cl,(ncol(evalDf)+1):ncol(lastEvalDf),decisionStump,lastEvalDf,enetRegression)
stumpStats3=unlist(stumpStats3)
splitByCol3=which(colnames(lastEvalDf)==compNames3[which(stumpStats3==max(stumpStats3))])
splitter3=lastEvalDf[,splitByCol3]


print("427")
fullEvalDf$assign=0 
fullEvalDf$assign[splitter1==1]=1
fullEvalDf$assign[splitter1==0][splitter2==1]=2
fullEvalDf$assign[splitter1==0][splitter2==0][splitter3==1]=3
fullEvalDf$assign[splitter1==0][splitter2==0][splitter3==0]=4

partRes=list()
aucRes=c()
for(i in 1:4){
  enetRes = enetRegression(fullEvalDf[fullEvalDf$assign==i,],checkMissing = T,makePlot = F)
  aucRes = c(aucRes,enetRes[[1]])
  partRes[[i]]=enetRes[[2]]
}
print("441")
splitAuc=(aucRes[1]*sum(fullEvalDf$assign==1) + aucRes[2]*sum(fullEvalDf$assign==2) +
            aucRes[3]*sum(fullEvalDf$assign==3) + aucRes[4]*sum(fullEvalDf$assign==4))/nrow(fullEvalDf)
fullRoc=partRes[[1]]*sum(fullEvalDf$assign==1) + partRes[[2]]*sum(fullEvalDf$assign==2) + 
  partRes[[3]]*sum(fullEvalDf$assign==3) + partRes[[4]]*sum(fullEvalDf$assign==4)
fullRoc = fullRoc/nrow(fullEvalDf)
print("447")

forPlot=rbind(cbind(fullRes[[2]],type="fullApply"),cbind(fullRoc,type="stumpApply"))
print("450")
thePlot=ggplot(forPlot,aes(fpr,tpr,color=type))+geom_point()+geom_abline(intercept = c(0,0),slope = 1)
  labs(x="false positive rate",y="true positive rate",
       caption=paste("AUCs: split =",round(splitAuc,4),"~ Full =",round(fullRes[[1]],4)),
       title="Regularized Logistic Regression Performance - Training",subtitle=paste("For Trait",currentTrait))
plot(thePlot)



################################################################################################################################
# now apply best to testing data ###############################################################################################
#fullEvalDfTest$array=as.numeric(as.factor(fullEvalDfTest$array))
#fullEvalDfTest$sex=as.numeric(as.factor(fullEvalDfTest$sex))
print("460")
fullEvalDfTest=cbind(testDf,suppTest,totalPhenoTest)
fullEvalDfTest$array=as.numeric(as.factor(fullEvalDfTest$array))
fullEvalDfTest$sex=as.numeric(as.factor(fullEvalDfTest$sex))

if(fullRes[[1]] > 0){ #MADE A CHANGE HERE!
  #Now stumping works best
  fullEvalDf=fullEvalDf[,-ncol(fullEvalDf)]
  netCv = cv.glmnet(y=fullEvalDf$response,x=as.matrix(fullEvalDf[,-9]),family="binomial",nfolds = 3)
  netLog = glmnet(y=fullEvalDf$response,x=as.matrix(fullEvalDf[,-9]),family="binomial",lambda = netCv$lambda.min)
  print("465")
  bestPred=predict(netLog, newx = as.matrix(fullEvalDfTest[,-9]))
  g=roc(fullEvalDfTest$response ~ bestPred[,1])
  print("468")
  finalRoc=cbind(g$sensitivities,1-g$specificities)
  finalAuc=g$auc
} else {
  print("469")
  #stumping works best
  testSplitter1=fullEvalDfTest[,colnames(fullEvalDfTest)==colnames(fullEvalDf)[splitByCol1]]
  testSplitter2=fullEvalDfTest[,colnames(fullEvalDfTest)==colnames(nextEvalDf)[splitByCol2]]
  testSplitter3=fullEvalDfTest[,colnames(fullEvalDfTest)==colnames(lastEvalDf)[splitByCol3]]
  print("477")
  fullEvalDfTest$assign=0 
  fullEvalDfTest$assign[testSplitter1==1]=1
  fullEvalDfTest$assign[testSplitter1==0 &testSplitter2==1]=2
  fullEvalDfTest$assign[testSplitter1==0 & testSplitter2==0 & testSplitter3==1]=3
  fullEvalDfTest$assign[testSplitter1==0 & testSplitter2==0 & testSplitter3==0]=4
  print("483")
  partRes=list()
  aucRes=c()
  for(i in 1:4){
    trainGroup=fullEvalDf[fullEvalDf$assign==i,]
    testGroup=fullEvalDfTest[fullEvalDfTest$assign==i,]
    
    checkedTrain=apply(trainGroup,2,function(z) length(unique(z)))
    checkedTest=apply(testGroup,2,function(z) length(unique(z)))
    evalChecked = checkedTest==1 | checkedTrain==1
    if(any(evalChecked)){
      trainGroup=trainGroup[,-which(evalChecked)]
      testGroup=testGroup[,-which(evalChecked)]
    }
    respCol=which(colnames(trainGroup)=="response")
    
    netCv = cv.glmnet(y=trainGroup$response,x=as.matrix(trainGroup[,-respCol]),family="binomial",nfolds = 3)
    netLog = glmnet(y=trainGroup$response,x=as.matrix(trainGroup[,-respCol]),family="binomial",lambda = netCv$lambda.min)
    print("491")
    bestPred=predict(netLog, newx = as.matrix(testGroup[,-respCol]))
    g=roc(testGroup$response ~ bestPred[,1])
    finalAuc=g$auc
    finalRoc=cbind(g$sensitivities,1-g$specificities)
    print("496")
    if(g$auc==0.5){
      g$specificities=seq(1,0,length.out = 1000)
      g$sensitivities=seq(0,1,length.out = 1000)
    } else if (length(g$sensitivities)<1000){
      g$sensitivities=approx(x=1:length(g$sensitivities),y=g$sensitivities,n=1000)$y
      g$specificities=approx(x=1:length(g$specificities),y=g$specificities,n=1000)$y
    }
    print("504")
    partRes[[i]]=cbind(g$sensitivities[round(seq(1,length(g$sensitivities),length.out = 1000))],1-g$specificities[round(seq(1,length(g$specificities),length.out = 1000))])
    aucRes=c(aucRes,g$auc)
  }
  print("503")
  finalAuc=(aucRes[1]*sum(fullEvalDfTest$assign==1) + aucRes[2]*sum(fullEvalDfTest$assign==2) +
              aucRes[3]*sum(fullEvalDfTest$assign==3) + aucRes[4]*sum(fullEvalDfTest$assign==4))/nrow(fullEvalDfTest)
  finalRoc=partRes[[1]]*sum(fullEvalDfTest$assign==1) + partRes[[2]]*sum(fullEvalDfTest$assign==2) + 
    partRes[[3]]*sum(fullEvalDfTest$assign==3) + partRes[[4]]*sum(fullEvalDfTest$assign==4)
  finalRoc = finalRoc/nrow(fullEvalDfTest)
}

#ROC Curve of the testing data set
addROC=testROC[testROC$group=="PRS+covars",]
scoreOnly=cbind(boostPlot,group="all_Scores")
forPlot=data.frame(finalRoc)
colnames(forPlot)=c("tpr","fpr")
forPlot$group="all_Info"
forPlot=rbind(forPlot,addROC,scoreOnly)
thePlot=ggplot(forPlot,aes(fpr,tpr,color=group))+geom_line()+geom_abline(intercept = c(0,0),slope = 1)+
  labs(x="False Positive Rate",y="True Positive Rate",title=paste("Testing Using All Info for",currentTrait),
       caption=paste0("AUCs: PRS+covars=",round(testAUC,3)," ~ allScores=",round(bestScoreAUC,3)," ~ allInfo=",round(finalAuc,3)))
       
plot(thePlot)

#final comparison of the best eval test against best full test

dev.off()

return(list(sendBackBestTracker,suppTrain,suppTest,forPlot,finalAuc,bestScoreAUC,allResReturn,scoreStatsReturn,suppTrain,totalPhenoTrain))
}
