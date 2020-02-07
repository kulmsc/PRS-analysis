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
library(cowplot)
theme_set(theme_cowplot())

boostEval <- function(authorComp,bestFile,compStat,trainDf,testDf,allSupportDecode,allDiseaseDecode,
                      trainROC,trainAUC,testROC,testAUC,maxBest, removeSimilar,decoder,currentTrait,covarAUC){

pdf(paste(authorComp,"boost","pdf",sep='.'))
 
#MUST REMOVE THIS LATER!!!
for(i in 1:length(diseaseTrainList)){
if(any(is.na(diseaseTrainList[[i]][nrow(diseaseTrainList[[i]]),]))){
    diseaseTrainList[[i]][nrow(diseaseTrainList[[i]]),] <- apply(diseaseTrainList[[i]][complete.cases(diseaseTrainList[[i]]),],2,mean)
}
}
 

#determine all the unique files names in the support files
allFiles=c(names(supportTrainList),names(diseaseTrainList))
names(allFiles)=c(rep("supportFiles",length(supportTrainList)),rep("diseaseFiles",length(diseaseTrainList)))

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
commonNames=gsub(".","-",commonNames,fixed = T)
commonNames=commonNames[-which(commonNames==authorComp)]
commonNames=unique(commonNames)

#check if the author has data in both a test and train file, add these to badfiles
checkStatusTrain=rep(0,length(commonNames))
for(checkDf in c(supportTrainList,diseaseTrainList)){
    checkStatusTrain[commonNames %in% names(checkDf)[apply(checkDf,2,function(x) any(x!=0))]]=1
}
checkStatusTest=rep(0,length(commonNames))
for(checkDf in c(supportTestList,diseaseTestList)){
  checkStatusTest[commonNames %in% names(checkDf)[apply(checkDf,2,function(x) any(x!=0))]]=1
}
badNames=commonNames[checkStatusTrain==0 | checkStatusTest==0]

#check which other authors relate to the phenotype I am analyzing
similarBadNames <- c()
if(removeSimilar){
  for(decodeCol in 4:ncol(phenoDefs)){
    for(matchTrait in strsplit(decoder[decodeCol],"|",fixed=TRUE)[[1]]){
      if(!is.na(matchTrait)){
      matchGrep <- grep(matchTrait,phenoDefs[,decodeCol])
      if(length(matchGrep) > 0){
        similarBadNames <- c(similarBadNames,phenoDefs[matchGrep,1])
      }
      }
    }
  }
}
badNames=unique(c(badNames,similarBadNames))

#remove the badNames if there are any
if(any(commonNames %in% badNames)){
  commonNames=commonNames[-which(commonNames %in% badNames)]
}


##########################################################################################################################
# Make best tracker ######################################################################################################

#create a holder data frame for all the support files and possible scores
bestTracker=data.frame(matrix(0,ncol=length(allFiles),nrow=length(commonNames)))
colnames(bestTracker)=allFiles
row.names(bestTracker)=sort(commonNames)

#for each support file check out the author score for the stat specificed, recording performance in the holder made above
parallelBT <- function(pullCol,trainDf,supportTrain,compStat,roc,caseIndices, controlIndices, caseSplits, controlSplits){
  normNew <- (supportTrain[,pullCol,drop=T] - mean(supportTrain[,pullCol,drop=T]))/sd(supportTrain[,pullCol,drop=T])
  evalDf=cbind(trainDf,supp=normNew)
  theStat <- rep(0,3) #am only building this for auc use
  
  for(i in 1:3){
    trainEvalDf <- evalDf[c(caseIndices[caseSplits==i], controlIndices[controlSplits==i]),]
    testEvalDf <- evalDf[c(caseIndices[caseSplits!=i], controlIndices[controlSplits!=i]),]
    
    completeLogistic <- glm(response ~ .,data=trainEvalDf,family="binomial")
    testEvalDf$pred=predict(completeLogistic,testEvalDf,type="response")
    testEvalDf=testEvalDf[order(testEvalDf$pred,decreasing = T),]
    
    if(compStat=="OR"){
      exposeGroup=evalDf[1:round(nrow(evalDf)*0.05),ncol(trainDf)]
      safeGroup=evalDf[(round(nrow(evalDf)*0.05)+1):nrow(evalDf),ncol(trainDf)]
      theStat=(sum(exposeGroup==1)/sum(exposeGroup==0))/(sum(safeGroup==1)/sum(safeGroup==0))
    } else if (compStat=="Prevalance") {
      theStat = sum(head(evalDf,round(nrow(evalDf)/10))$response==1)/sum(tail(evalDf,round(nrow(evalDf)/10))$response==1)
    } else {
      g=pROC::roc(response ~ pred,data=testEvalDf)
      theStat[i]=as.numeric(g$auc)
    }
  }
  
  return(mean(theStat)) 
} 

for(j in 1:length(allFiles)){
  print("start")
  theFile=allFiles[j]
  
  #get the matched train and test file
  if(names(theFile)=="supportFiles"){
    supportTrain=supportTrainList[[which(names(supportTrainList)==theFile)]]
    supportTest=supportTestList[[which(names(supportTestList)==gsub("train","test",theFile))]]
  } else {
    supportTrain=diseaseTrainList[[which(names(diseaseTrainList)==theFile)]]
    supportTest=diseaseTestList[[which(names(diseaseTestList)==gsub("train","test",theFile))]]
  }
  
  supportTrain <- supportTrain[,colnames(supportTrain) %in% colnames(supportTest)]
  supportTest <- supportTest[,colnames(supportTest) %in% colnames(supportTrain)]
  
  if(decoder[3]!="A"){
    supportTrain <- supportTrain[covar1$sex==decoder[3],]
    supportTest <- supportTest[covar2$sex==decoder[3],]
    supportTrain <- supportTrain[,-3]
    supportTest <- supportTest[,-3]
  }
  
  if(ncol(supportTrain) != ncol(supportTest)){print(theFile)}
  #remove any columns that have zero results
  if(sum(colSums(supportTrain)==0)>0){
    supportTrain=supportTrain[,-which(colSums(supportTrain)==0 | colSums(supportTest)==0)]
  }

  #make sure columns match
  supportTrain=supportTrain[,colnames(supportTrain) %in% colnames(supportTest)]
  supportTest=supportTest[,colnames(supportTest) %in% colnames(supportTest)]

  #remove the authorComp and the badnames
  if(authorComp %in% colnames(supportTrain)){
    supportTrain=supportTrain[,-which(colnames(supportTrain)==authorComp)]
  }
  if(sum(colnames(supportTrain) %in% badNames)>0){
    supportTrain=supportTrain[,-which(colnames(supportTrain) %in% badNames)]
  }

  #fix the column names
  toChange=grepl(".",colnames(supportTrain),fixed=T)
  colnames(supportTrain)[toChange]=gsub(".","-",colnames(supportTrain)[toChange],fixed=T)

  #get the case and controls
  caseIndices <- which(trainDf$response==1)
  if(sum(trainDf$response==1)*25 > sum(trainDf$response==0)){
    controlIndices <- which(trainDf$response==0)
  } else {
    controlIndices <- sample(which(trainDf$response==0),length(caseIndices)*25)
  }
  caseSplits <- sample(ceiling(seq_along(caseIndices)/round(length(caseIndices)/3)))
  controlSplits <- sample(ceiling(seq_along(controlIndices)/round(length(controlIndices)/3)))
  
  #run the analysis
  #cl <- makeCluster(2,outfile="")
  theStats=lapply(1:ncol(supportTrain),parallelBT,trainDf,supportTrain,compStat,roc,
                  caseIndices, controlIndices, caseSplits, controlSplits)
  #stopCluster(cl)
  theStats=unlist(theStats)
  bestTracker[which(row.names(bestTracker) %in% colnames(supportTrain)),j]=theStats
}


#take the best file for each author and made/ a data frame that contains the scores
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

suppTrain=data.frame(matrix(0,nrow=nrow(trainDf),ncol=nrow(sendBackBestTracker)))
suppTest=data.frame(matrix(0,nrow=nrow(testDf),ncol=nrow(sendBackBestTracker)))
colnames(suppTrain)=row.names(sendBackBestTracker)
colnames(suppTest)=row.names(sendBackBestTracker)

print("b")
#reading in the test support scores and distributing out the train scores which are already read in
for(suppFile in unique(sendBackBestTracker$bestFile)){
  print(suppFile)
  print(dim(suppTrain))
  print(dim(suppTest))

  if(str_split(suppFile,pattern=fixed("."))[[1]][1]=="support"){
    supportTrain=supportTrainList[[which(names(supportTrainList)==suppFile)]]
    supportTest=supportTestList[[which(names(supportTestList)==gsub("train","test",suppFile))]]
  } else {
    supportTrain=diseaseTrainList[[which(names(diseaseTrainList)==suppFile)]]
    supportTest=diseaseTestList[[which(names(diseaseTestList)==gsub("train","test",suppFile))]]
  }
  toChange=grepl(".",colnames(supportTrain),fixed=T)
  colnames(supportTrain)[toChange]=gsub(".","-",colnames(supportTrain)[toChange],fixed=T)
  toChange=grepl(".",colnames(supportTest),fixed=T)
  colnames(supportTest)[toChange]=gsub(".","-",colnames(supportTest)[toChange],fixed=T)

  if(decoder[3]!="A"){
    supportTrain <- supportTrain[covar1$sex==decoder[3],]
    supportTest <- supportTest[covar2$sex==decoder[3],]
  }
  
  suppTrain[,sendBackBestTracker$bestFile == suppFile] <- 
    as.data.frame(supportTrain[,which(colnames(supportTrain) %in% rownames(sendBackBestTracker)[sendBackBestTracker$bestFile==suppFile])])
  suppTest[,sendBackBestTracker$bestFile == suppFile] <- 
    as.data.frame(supportTest[,which(colnames(supportTest) %in% rownames(sendBackBestTracker)[sendBackBestTracker$bestFile==suppFile])])
}
suppTrain=apply(suppTrain,2, function(x) (x-mean(x))/sd(x))
suppTest=apply(suppTest,2, function(x) (x-mean(x))/sd(x))

print("c")
#do the decoder stuff - get matrix going from author to disease name
supportDefs <- read.table("supportDefs",header=T,stringsAsFactors=F)
supportDecode=rbind(phenoDefs[,1:2],supportDefs)
supportDecode=supportDecode[order(supportDecode$author),]
colnames(supportDecode)[2] <- "trait"
supportDecode <- supportDecode[supportDecode$author %in% rownames(sendBackBestTracker),]

#####################################################################################################################
#make the plot of the individual logreg features ####################################################################

#get the baseline stat for plotting the relative change of adding the covariates
evalDf=data.frame(trainDf)

caseIndices <- which(trainDf$response==1)
if(length(caseIndices)*50 > length(controlIndices)){
  controlIndices <- sample(which(trainDf$response==0))
} else {
  controlIndices <- sample(which(trainDf$response==0),length(caseIndices)*50)
}
caseSplits <- sample(ceiling(seq_along(caseIndices)/round(length(caseIndices)/3)))
controlSplits <- sample(ceiling(seq_along(controlIndices)/round(length(controlIndices)/3)))

print("d")
theStat <- rep(0,3)
for(i in 1:3){
  trainEvalDf <- evalDf[c(caseIndices[caseSplits==i], controlIndices[controlSplits==i]),]
  testEvalDf <- evalDf[c(caseIndices[caseSplits!=i], controlIndices[controlSplits!=i]),]
  completeLogistic <- glm(response ~ .,data=trainEvalDf,family="binomial")
  testEvalDf$pred=predict(completeLogistic,testEvalDf,type="response")
  testEvalDf=testEvalDf[order(testEvalDf$pred,decreasing = T),]
  if(compStat=="OR"){
    exposeGroup=evalDf[1:round(nrow(evalDf)*0.05),9]
    safeGroup=evalDf[(round(nrow(evalDf)*0.05)+1):nrow(evalDf),9]
    theStat=(sum(exposeGroup==1)/sum(exposeGroup==0))/(sum(safeGroup==1)/sum(safeGroup==0))
  } else if (compStat=="Prevalance") {
    theStat = sum(head(evalDf,round(nrow(evalDf)/10))$response==1)/sum(tail(evalDf,round(nrow(evalDf)/10))$response==1)
  } else {
    g=pROC::roc(response ~ pred,data=testEvalDf)
    theStat[i]=as.numeric(g$auc)
  }
}
theStat <- mean(theStat)
print("made theStat")
print(rownames(bestTracker))
print(supportDecode)
print(length(apply(bestTracker[,-ncol(bestTracker)],1,max)))
print(length(factor(supportDecode$trait[supportDecode$author %in% rownames(bestTracker) ])))

#now make the plot of all supports from the best files by the stat it was analyzed by
plotDf=data.frame(stat=apply(bestTracker[,-ncol(bestTracker)],1,max),
                  trait=factor(supportDecode$trait[supportDecode$author %in% rownames(bestTracker) ]))
print(plotDf)
plotDf$trait=factor(plotDf$trait,levels = unique(plotDf$trait[order(plotDf$stat)]))
plotDf$stat=plotDf$stat-theStat
plotDf=plotDf[1:10,]
thePlot=ggplot(plotDf,aes(trait,stat))+geom_col()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
  labs(x="Unrelated Score",y=paste("Difference in",compStat),
       title="Scores - Training Application",
       subtitle=paste(decoder[2],"Model"),caption=paste("Models include",currentTrait,"score"))
plot(thePlot)

print("e")
################################################################################################################################
# prove in depth analysis of each score's interaction with disease #############################################################

parallelA<-function(pullCol,allTrain,allTest,trainDf,testDf,roc){
  library(pROC)
  #pval
  trainDf$predictorNew=allTrain[,pullCol]
  testDf$predictorNew=allTest[,pullCol]
  completeLogistic <- glm(response ~ . ,data=trainDf,family="binomial")
  covarLogistic <- glm(response ~ . - predictorNew,data=trainDf,family="binomial")
  funcPval= -log10(coef(summary(completeLogistic))[ncol(trainDf),4])
  
  #AUC
  testDf$prob=predict(completeLogistic,testDf,type="response")
  testDf$base=predict(covarLogistic,testDf,type="response")
  gScore=pROC::roc(response ~ prob,data=testDf)
  gCovar=pROC::roc(response ~ base,data=testDf)
  funcAuc=as.numeric(ci.auc(gScore$auc))-as.numeric(ci.auc(gCovar$auc))
  
  #OR
  print("at this 0.01")
  testDf=testDf[order(testDf$predictorNew,decreasing = T),]
  allCutOffs=c(0.5,0.2,0.1,0.05,0.01,0.005)
  i=1
  funcOR=rep(0,length(allCutOffs))
  for(cutOff in allCutOffs){
    splitPoint=round(nrow(testDf)*cutOff)
    exposeGroup=testDf[1:splitPoint,(ncol(trainDf)-1)]
    safeGroup=testDf[(splitPoint+1):nrow(testDf),(ncol(trainDf)-1)]
    print(length(exposeGroup))
    funcOR[i]=(sum(exposeGroup==1)/sum(exposeGroup==0))/(sum(safeGroup==1)/sum(safeGroup==0))
    i=i+1
  }
  
  #Prev
  funcPrev=rep(0,100)
  testDf$group=rep(0,nrow(testDf))
  groupRanges=c(round(seq(1,nrow(testDf),by = nrow(testDf)/100)),nrow(testDf))
  for(i in 1:100){
    testDf$group[groupRanges[i]:groupRanges[i+1]]=i
  }
  testDf$group=rev(testDf$group)
  for(i in 1:100){
    withDis=sum(testDf[testDf$group==i,(ncol(trainDf)-1)]==1)
    totalGroup=sum(testDf$group==i)
    funcPrev[i]=withDis/totalGroup
  }
  
  return(list(funcPval,funcAuc,funcOR,funcPrev))
}

makePlots <- function(allRes,trainScores,translateLabels){
  print("makingPlots")
  colnames(trainScores)=translateLabels[colnames(trainScores) %in% translateLabels[,1],2]
  
  pvals=unlist(lapply(allRes, `[[`, 1))
  aucImps=t(matrix(unlist(lapply(allRes, `[[`, 2)),nrow=3))
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
  aucImpsDf=fullAucImpsDf[order(fullAucImpsDf[,2],decreasing = T),]
  aucImpsDf=aucImpsDf[1:8,]
  aucImpsDf[,4]=factor(aucImpsDf[,4],levels = aucImpsDf[order(aucImpsDf[,1]),4])
  colnames(aucImpsDf)[1:3]=c("lo","aucImps","hi")
  thePlot=ggplot(aucImpsDf,aes(traits,aucImps))+geom_col(aes(fill=traits))+geom_errorbar(aes(ymin=lo,ymax=hi))+
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

makePlotsColor <- function(allRes,trainScores,translateLabels,totalTraits,numberColor){
  print("make plots color")
  #colnames(trainScores)=translateLabels[colnames(trainScores) %in% translateLabels[,1],2]
  for(t in 1:ncol(trainScores)){
    colnames(trainScores)[t] <- translateLabels[translateLabels[,1]==colnames(trainScores)[t],2][1]
  }
  
  pvals=unlist(lapply(allRes, `[[`, 1))
  aucImps=t(matrix(unlist(lapply(allRes, `[[`, 2)),nrow=3))
  ors=t(matrix(unlist(lapply(allRes, `[[`, 3)),nrow=6))
  prevs=t(matrix(unlist(lapply(allRes, `[[`, 4)),nrow=100))
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  print(numberColor)
  print(totalTraits)
  specColors=c(gg_color_hue(numberColor),rep("#808080",totalTraits-numberColor))
  print(specColors)

  print("AUC")  
  #AUC
  fullAucImpsDf=data.frame(aucImps,traits=colnames(trainScores),stringsAsFactors = F)
  aucImpsDf=fullAucImpsDf[order(fullAucImpsDf[,2],decreasing = T),]
  aucImpsDf=aucImpsDf[1:totalTraits,]
  names(specColors)=aucImpsDf$traits
  aucImpsDf[,4]=factor(aucImpsDf[,4],levels = unique(aucImpsDf[order(aucImpsDf[,1]),4]))
  colnames(aucImpsDf)[1:3]=c("lo","aucImps","hi")
  #thePlot=ggplot(aucImpsDf,aes(traits,aucImps))+geom_point(aes(color=traits),size=4)+geom_errorbar(aes(ymin=lo,ymax=hi),width=0.2)+
  #  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
  #  labs(x="Scores",y="AUC Improvement",title="Testing Application: AUC Improvement",caption=paste("Models include",currentTrait,"score"))+
  #  scale_fill_manual(values=specColors)+guides(color=FALSE)
  #plot(thePlot)
  
  print("pvals")
  #PVALS
  fullPvalsDf=data.frame(pvals,traits=colnames(trainScores),stringsAsFactors = F)
  pvalsDf=fullPvalsDf[fullPvalsDf$traits %in% as.character(aucImpsDf$traits),]
  pvalsDf=pvalsDf[1:totalTraits,]
  pvalsDf[,2]=factor(pvalsDf[,2],levels = unique(pvalsDf[order(pvalsDf[,1]),2]))
  thePlot=ggplot(pvalsDf,aes(traits,pvals,fill=traits))+geom_col()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
    labs(x="Scores",y="logistic regression PVal",title=paste("Testing Application: Logistic Regression PVals"),
         caption=paste("Models include",currentTrait,"score"))+
    scale_fill_manual(values=specColors)+guides(fill=FALSE)
  plot(thePlot)
  
  print("vs")
  #OR vs PVALS
  fullOrPvalDf=data.frame(pvals,ors=ors[,4],traits=colnames(trainScores),stringsAsFactors = F)
  orPvalDf=data.frame(fullOrPvalDf)
  orPvalDf$word=orPvalDf$traits
  orPvalDf$word[orPvalDf$pvals < sort(orPvalDf[,1],decreasing = T)[3] & orPvalDf$ors < sort(orPvalDf[,2],decreasing = T)[3]]=""
  thePlot=ggplot(orPvalDf,aes(pvals,ors,label=word))+
    geom_point(data=orPvalDf[orPvalDf$traits %in% names(specColors[1:3]),],aes(color=traits),size=4)+
    geom_point(data=orPvalDf[!(orPvalDf$traits %in% names(specColors[1:3])),])+
    geom_label_repel(size=4,force=3)+
    labs(x="PValue",y="Odds Ratio",title="Testing Application: Effect and Certainty",caption=paste("Models include",currentTrait,"score"))+
    theme(legend.position="bottom",legend.direction="vertical")+guides(color=F)+
    scale_color_manual(values=specColors[1:3])
  plot(thePlot)
  
  print("or")
  #OR
  print("beforeOR")
  colnames(ors)=paste("cut0ff",c(0.5,0.2,0.1,0.05,0.01,0.005),sep=".")
  fullOrsDf=data.frame(ors,traits=colnames(trainScores),stringsAsFactors = F)
  orsDf=fullOrsDf[fullOrsDf$traits %in% as.character(aucImpsDf$traits),]
  orsDf=melt(orsDf,id.vars="traits")
  print(nrow(orsDf))
  print((nrow(orsDf)/6))
  orsDf[,2]=sort(rep(c(0.5,0.2,0.1,0.05,0.01,0.005),(nrow(orsDf)/6)),decreasing = T)
  #orsDf[,2]=sort(rep(c(0.5,0.2,0.1,0.05,0.01,0.005),totalTraits),decreasing = T)
  orsDf[,2]=factor(orsDf[,2],levels=c(0.5,0.2,0.1,0.05,0.01,0.005))
  thePlot=ggplot(orsDf,aes(variable,value,color=traits))+geom_line(aes(group=traits))+
    labs(x="Percentile Cut-Off",y="Odds Ratio",title="Testing Application: Odds ratios",caption=paste("Models include",currentTrait,"score"))+
    scale_color_manual(values=specColors)+guides(color=FALSE)
  plot(thePlot)
  print("afterOR")  

  print("prev")
  #PREV
  colnames(prevs)=1:100/100
  fullPrevsDf=data.frame(prevs,traits=colnames(trainScores),stringsAsFactors = F)
  prevsDf=fullPrevsDf[fullPrevsDf$traits %in% as.character(aucImpsDf$traits),]
  prevsDf=melt(prevsDf,id.vars = "traits")
  prevsDf[,2]=sort(rep((1:100/100),(nrow(prevsDf)/100)))
  thePlot=ggplot(prevsDf,aes(variable,value,color=traits))+geom_smooth(se=F)+
    labs(x="Score Percentile",y="Trait Prevalence",title="Testing Application: Disease Prevalences",
         caption=paste("Models include",currentTrait,"score"))+
    scale_color_manual(values=specColors)+guides(color=FALSE)
  plot(thePlot)
}

#cl <- makeCluster(2,outfile="")
allResReturn=lapply(1:ncol(suppTrain),parallelA,suppTrain,suppTest,trainDf,testDf,roc)
print("490")
print(length(allResReturn))
print(dim(allResReturn))
print(dim(suppTrain))
names(allResReturn) <- colnames(suppTrain)
#stopCluster(cl)
makePlotsColor(allResReturn,suppTrain,supportDecode,8,3)

################################################################################################################################
###############################################################################################################################

print("h")
#Now we have all the best support Features according to either or, prev, or auc - can now do feature selection
evalDf=cbind(trainDf,suppTrain)
evalDf$array=as.numeric(as.factor(evalDf$array))
if(decoder[3]=="A"){
  evalDf$sex=as.numeric(as.factor(evalDf$sex))
}

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
if(decoder[3]=="A"){
  splitEvalDf$sex=as.numeric(as.factor(splitEvalDf$sex))
}

###################################################################################################################
#Feature selection through logistic lasso regression ##############################################################
###################################################################################################################

###################################################################################################################
#enetRegression function ##########################################################################################
enetWrapper <- function(evalDf,bestLamb=NULL,justAuc=F,checkMissing=F,possIndex=1,makePlot=T,
                        addTitle="",isSplit=F, caseLimit=50, iter=0, enetRegression, enetWrapper){
  print("start enetWrapper")
  print(table(evalDf$response))
  if(sum(evalDf$response==1)<10 | sum(evalDf$response==0)<10){
    print("wrap checkWrong")
    if(justAuc){
      return(list(rep(0.5,3)))
    } else {
      return(list(rep(0.5,3),data.frame(seq(0,1,length.out = 1000),seq(1,0,length.out = 1000))))
    }
  }
  
  if(sum(evalDf$response==0) > 40*(sum(evalDf$response==1))){
    print("first reduce")
    evalDf <- evalDf[unique(c(which(evalDf$response==1),sample(1:nrow(evalDf),(20*sum(evalDf$response==1))))),]
    print(table(evalDf$response))
  }
  
  ans <- enetRegression(evalDf,bestLamb,justAuc,checkMissing,possIndex,makePlot,addTitle,isSplit,caseLimit)
  
  if(ans[[3]]==FALSE){
    print("RETRY")
    iter=iter+1
    evalDf <- evalDf[unique(c(which(evalDf$response==1),sample(1:nrow(evalDf),(exp(-iter/2)*nrow(evalDf))))),]
    if(sum(evalDf$response==1) > sum(evalDf$response==0)){
      ans <- enetWrapper(evalDf,bestLamb,justAuc,checkMissing,possIndex,makePlot,addTitle,isSplit,nrow(evalDf),iter,enetRegression,enetWrapper)
    } else {
      ans <- enetWrapper(evalDf,bestLamb,justAuc,checkMissing,possIndex,makePlot,addTitle,isSplit,caseLimit,iter,enetRegression,enetWrapper)
    }
  } else {
    return(ans)
  }
}

#the issplit argument will be used to split up the evalDf into the testing and training
enetTest <- function(evalDf,bestLamb=NULL,justAuc=F,checkMissing=F,possIndex=1,makePlot=T,addTitle="",isSplit=F, caseLimit=50){
  start("enetTest")
  train <- data.frame(evalDf)
  test <- data.frame(isSplit)
  respCol <- which(colnames(evalDf)=="response")
  print("colnames")
  print(colnames(test))
  print(colnames(train))
  
  wentEnet=TRUE
  noFails=TRUE
  if(sum(train$response==1) >= caseLimit){
    print("finding the best lambda value")
    if(is.null(bestLamb)){
      netCv = try(cv.glmnet(y=evalDf$response,x=as.matrix(evalDf[,-respCol]),family="binomial",nfolds = 3,maxit=5000),outFile = "test")
      if(class(netCv)=="try-error"){
        print("UTTER FAILURE")
        #if(justAuc==T){
        #  return(rep(0.5,3))
        #} else {
          return(list(rep(0.5,3),data.frame(seq(0,1,length.out = 1000),seq(1,0,length.out = 1000)),FALSE))
        #}
      } else {
        bestLamb=netCv$lambda.min 
      }
    }
  }
  
  print("Doing test - no cross folding")
  if(sum(test$response==1)==0){
    print("VERY WRONG")
    test$response[1]=1
  }
  if(checkMissing){
    checkedTrain=apply(train,2,function(z) length(unique(z)))
    checkedTest=apply(test,2,function(z) length(unique(z)))
    checkedTrain[respCol]=FALSE; checkedTest[respCol]=FALSE
    evalChecked = checkedTest==1 | checkedTrain==1
    if(any(evalChecked)){
      print("was missing")
      train=train[,-which(evalChecked)]
      test=test[,-which(evalChecked)]
      print(colnames(train)[1:10])
    }
  }
  respCol=which(colnames(train)=="response")
  if(sum(evalDf$response==1) < caseLimit){
    print("logistic regression")
    print(table(train$response))
    print(table(test$response))
    wentEnet=FALSE
    train=train[,1:respCol]
    netLog <- glm(response ~ .,data=train,family="binomial")
    bestPred=predict(netLog, test, type="response")
    g=pROC::roc(test$response ~ bestPred)
  } else {
    print("elastic net regression")
    netLog = glmnet(y=train$response,x=as.matrix(train[,-respCol]),family="binomial",lambda = netCv$lambda.min,maxit = 1000)
    bestPred=predict(netLog, newx = as.matrix(test[,-respCol]))
    g=pROC::roc(test$response ~ bestPred[,1])
  }
  
  if(g$auc==0.5){
    print("was 0.5")
    noFails=FALSE
    g$specificities=seq(0,1,length.out = 1000)
    g$sensitivities=seq(1,0,length.out = 1000)
  } else if (length(g$sensitivities)<1000){
    print("was short")
    print(g$auc)
    g$sensitivities=approx(x=1:length(g$sensitivities),y=g$sensitivities,n=1000)$y
    g$specificities=approx(x=1:length(g$specificities),y=g$specificities,n=1000)$y
  } 
  
  if(justAuc){
    return(list(as.numeric(ci.auc(g)),noFails,noFails))
  } else {
    tpr=g$sensitivities[round(seq(1,length(g$sensitivities),length.out = 1000))]
    fpr=1-g$specificities[round(seq(1,length(g$specificities),length.out = 1000))]
    rocPlot <- data.frame(tpr,fpr)
    return(list(as.numeric(ci.auc(g)),rocPlot,noFails))
  }
    
}


enetRegression <- function(evalDf,bestLamb=NULL,justAuc=F,checkMissing=F,possIndex=1,makePlot=T,addTitle="",isSplit=F, caseLimit=50){
  start("enetRegression")
  wentEnet=TRUE
  noFails=TRUE
  respCol <- which(colnames(evalDf)=="response")
  if(sum(evalDf$response==1) >= caseLimit){
    print("finding the best lambda value")
    if(is.null(bestLamb)){
      netCv = try(cv.glmnet(y=evalDf$response,x=as.matrix(evalDf[,respCol]),family="binomial",nfolds = 3,maxit=5000),outFile = "test")
      if(class(netCv)=="try-error"){
        print("UTTER FAILURE")
        #if(justAuc==T){
        #  return(list(rep(0.5,3)))
        #} else {
          return(list(rep(0.5,3),data.frame(seq(0,1,length.out = 1000),seq(1,0,length.out = 1000)),FALSE))
        #}
      } else {
        bestLamb=netCv$lambda.min 
      }
    }
  }
  
  print("cross folding validation")
  theFolds=createFolds(evalDf[,respCol],k=3)
  tpr=c(); fpr=c(); auc=c()
  for(i in 1:3){
    train=evalDf[-theFolds[[i]],]
    test=evalDf[theFolds[[i]],]
    print("the tables")
    print(table(test$response))
    print(table(train$response))
    if(sum(as.numeric(as.character(test$response))==1)==0 | sum(as.numeric(as.character(train$response))==0)==0){
      print("VERY WRONG")
      if(justAuc==T){
          return(list(rep(0.5,3),TRUE,TRUE))
        } else {
          return(list(rep(0.5,3),data.frame(seq(0,1,length.out = 1000),seq(1,0,length.out = 1000)),FALSE))
      }
    }
    if(checkMissing){
      checkedTrain=apply(train,2,function(z) length(unique(z)))
      checkedTest=apply(test,2,function(z) length(unique(z)))
      checkedTrain[respCol]=FALSE; checkedTest[respCol]=FALSE
      evalChecked = checkedTest==1 | checkedTrain==1
      if(any(evalChecked)){
        print("was missing")
        train=train[,-which(evalChecked)]
        test=test[,-which(evalChecked)]
        print(colnames(train)[1:10])
      }
    }
    respCol=which(colnames(train)=="response")
    if(sum(evalDf$response==1) < caseLimit){
      print("logistic regression")
      print(table(train$response))
      print(table(test$response))
      wentEnet=FALSE
      train=train[,1:respCol]
      netLog <- glm(response ~ .,data=train,family="binomial")
      bestPred=predict(netLog, test, type="response")
      g=pROC::roc(test$response ~ bestPred)
    } else {
      print("elastic net regression")
      netLog = glmnet(y=train$response,x=as.matrix(train[,-respCol]),family="binomial",lambda = netCv$lambda.min,maxit = 1000)
      bestPred=predict(netLog, newx = as.matrix(test[,-respCol]))
      g=pROC::roc(test$response ~ bestPred[,1])
    }
    
    print("got the g")
    if(g$auc==0.5){
      print("was 0.5")
      noFails=FALSE
      g$specificities=seq(0,1,length.out = 1000)
      g$sensitivities=seq(1,0,length.out = 1000)
    } else if (length(g$sensitivities)<1000){
      print("was short")
      print(g$auc)
      g$sensitivities=approx(x=1:length(g$sensitivities),y=g$sensitivities,n=1000)$y
      g$specificities=approx(x=1:length(g$specificities),y=g$specificities,n=1000)$y
    } 
    print("new test")
    print(head(g$sensitivities))
    tpr=cbind(tpr,g$sensitivities[round(seq(1,length(g$sensitivities),length.out = 1000))])
    fpr=cbind(fpr,1-g$specificities[round(seq(1,length(g$specificities),length.out = 1000))])
    
    print("returning")
    auc=rbind(auc,as.numeric(ci.auc(g)))
  }
  
  if(justAuc){
    print("go justAuc")
    toReturn <- list(apply(auc,2,mean),noFails,noFails)
    print("made return")
    return(toReturn)
  }
  print("make rocPlot")
  rocPlot=data.frame(fpr=apply(fpr,1,mean),tpr=apply(tpr,1,mean))
  print("made rocPlot") 
 
  #the plotting of betas of each feature
  if(makePlot & wentEnet & noFails){
    print("making plot")
    forPlot=data.frame(traits=as.character(rownames(netLog$beta)),beta=as.numeric(netLog$beta))
    forPlot$traits=as.character(forPlot$traits)
    if(isSplit){
      splitTraits=strsplit(as.character(forPlot$traits),"_")
      ri=1
      for(st in splitTraits){
        if(length(st)==3){
          forPlot$traits[ri]=paste(supportDecode[supportDecode$author==paste(st[1],st[2],sep="_"),2],st[3],sep="_")
        }
        ri=ri+1
      }
    } else {
      forPlot[forPlot$traits %in% supportDecode$author,1]=supportDecode[supportDecode$author %in% forPlot$traits,2]
    }
    dupTraits=forPlot$traits[duplicated(forPlot$traits)]
    forPlot$traits[duplicated(forPlot$traits)]=paste0(dupTraits,1:length(dupTraits))
    forPlot$traits=factor(forPlot$traits,levels=forPlot$traits[order(forPlot$beta)])
    cutOff=sort(abs(forPlot$beta),decreasing = T)[12]
    forPlot=forPlot[abs(forPlot$beta) > cutOff,]
    
    print(forPlot)
    thePlot=ggplot(forPlot,aes(traits,beta))+geom_col(aes(fill=traits))+
      labs(title=paste("Training Application -",addTitle),subtitle=paste("For Trait",currentTrait))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+guides(fill=FALSE)
    plot(thePlot)
  }
    
  return(list(apply(auc,2,mean),rocPlot,noFails))
}

######################################################################################################################
#analyze the supports and split supports #############################################################################

#The ROC plot
print("748a")
evalReturn=enetWrapper(evalDf,bestLamb=NULL,justAuc=F,checkMissing=F,possIndex=1,makePlot=T,
  addTitle = "Full Scores",isSplit=F, caseLimit=50, iter=0, enetRegression, enetWrapper)

print("748b")
splitEvalReturn=enetWrapper(splitEvalDf,bestLamb=NULL,justAuc=F,checkMissing=F,possIndex=1,makePlot=T,
                       addTitle = "Split Scores",isSplit=T, caseLimit=50, iter=0, enetRegression, enetWrapper)
print("748c")
#produce roc plot to compare splitting the supports
evalReturn[[2]]$group="noSplit"
splitEvalReturn[[2]]$group="split"
allAuc=round(rbind(evalReturn[[1]],splitEvalReturn[[1]]),3)

forPlot=rbind(trainROC,evalReturn[[2]],splitEvalReturn[[2]])
forPlot$group=factor(forPlot$group)
thePlot=ggplot(forPlot,aes(fpr,tpr,color=group))+geom_line()+
  labs(title="Regularized LogReg with Scores - Training",
       subtitle=paste("For Trait",currentTrait),
       caption=paste("AUCs: noSplit =",allAuc[1,2],"(",allAuc[1,1],"-",allAuc[1,3],")",
                     "~ split =",allAuc[2,2],"(",allAuc[2,1],"-",allAuc[2,3],")"),
       y="true positive rate",x="false positive rate")+
  geom_abline(intercept = c(0,0),slope = 1)
plot(thePlot)

splitVnosplitAucs <- allAuc

######################################################################################################################
### apply better of split vs no split to the testing data ############################################################ 
if(allAuc[1,2]>allAuc[2,2]){
  testEvalDf=cbind(testDf,suppTest)
} else {
  evalDf=splitEvalDf
  splitTest=makeSplit(suppTest)
  testEvalDf=cbind(testDf,splitTest)
}
testEvalDf$array=as.numeric(as.factor(testEvalDf$array))
if(decoder[3]=="A"){
  testEvalDf$sex=as.numeric(as.factor(testEvalDf$sex))
}

print("790")
#LASSO Logistic Regression
netCv = cv.glmnet(y=evalDf$response,x=as.matrix(evalDf[,-ncol(trainDf)]),family="binomial",nfolds = 3)
netLog = glmnet(y=evalDf$response,x=as.matrix(evalDf[,-ncol(trainDf)]),family="binomial",lambda = netCv$lambda.min)
bestPred=predict(netLog, newx = as.matrix(testEvalDf[,-ncol(trainDf)]))
g=pROC::roc(testEvalDf$response ~ bestPred[,1])
#auc=c(auc,g$auc)

firstBoostAuc=round(as.numeric(ci.auc(g$auc)),3)
boostPlot=cbind(tpr=g$sensitivities,fpr=1-g$specificities)
boostPlot=as.data.frame(boostPlot[seq(1,nrow(boostPlot),length.out = 1001),])
testROCSmall=testROC[testROC$group=="PRS+covars",1:2]
forPlot=rbind(cbind(boostPlot,type="enet"),cbind(testROCSmall,type="logreg"))

thePlot=ggplot(forPlot,aes(fpr,tpr,color=type))+geom_line()+
  labs(title="Regularized LogReg with Scores- Testing",
       subtitle=paste("For Trait",currentTrait),
       caption=paste("AUCs: logreg =",testAUC[2],"(",testAUC[1],"-",testAUC[3],")","~",
                     "enet =",firstBoostAuc[2],"(",firstBoostAuc[1],"-",firstBoostAuc[3],")"),
       x="true positive rate",y="false positive rate")+
  geom_abline(intercept=c(0,0),slope=1)
plot(thePlot)
bestScoreAUC=firstBoostAuc
  

######################################################################################################################
######################################################################################################################
#Disease status analysis #############################################################################################

#####################################################################################################################
#create the data frame and analyze each feature separately ##########################################################
print("disease status analysis")
evalDf=cbind(trainDf,suppTrain)
evalDf$array=as.numeric(as.factor(evalDf$array))
if(decoder[3]=="A"){
  evalDf$sex=as.numeric(as.factor(evalDf$sex))
}

allScoresTest <- testPhenos[[2]]
allScoresTrain <- trainPhenos[[2]]
allSupportsTrain <- trainPhenos[[1]]
allSupportsTest <- testPhenos[[1]]

if(decoder[3]!="A"){
  allScoresTest <- allScoresTest[covar2$sex==decoder[3], ]
  allScoresTrain <- allScoresTrain[covar1$sex==decoder[3], ]
  allSupportsTest <- allSupportsTest[covar2$sex==decoder[3], ]
  allSupportsTrain <- allSupportsTrain[covar1$sex==decoder[3], ]
}

allScoresTest <- allScoresTest[colnames(allScoresTest) %in% colnames(allScoresTrain)]
allScoresTrain <- allScoresTrain[colnames(allScoresTrain) %in% colnames(allScoresTest)]
sumScoresTrain=apply(allScoresTrain,2,sum)
allScoresTrain=allScoresTrain[,sumScoresTrain>500]
allScoresTest=allScoresTest[,sumScoresTrain>500]

allSupportsTest <- allSupportsTest[,colnames(allSupportsTest) %in% colnames(allSupportsTrain)]
allSupportsTrain <- allSupportsTrain[,colnames(allSupportsTrain) %in% colnames(allSupportsTest)]
sumSupportsTrain=apply(allSupportsTrain,2,sum)
allSupportsTrain=allSupportsTrain[,sumSupportsTrain>500]
allSupportsTest=allSupportsTest[,sumSupportsTrain>500]

totalPhenoTrain=cbind(allScoresTrain,allSupportsTrain)
totalPhenoTest=cbind(allScoresTest,allSupportsTest)




if(removeSimilar){
  checkTraits <- unlist(strsplit(decoder[4:5],"|",fixed=T))
  checkTraits <- checkTraits[!(is.na(checkTraits))]
  if(any(colnames(totalPhenoTrain) %in% checkTraits)){
    totalPhenoTrain=totalPhenoTrain[,-which(colnames(totalPhenoTrain) %in% checkTraits)]
    totalPhenoTest=totalPhenoTest[,-which(colnames(totalPhenoTest) %in% checkTraits)]
  }
}

print("865")
parallelD <- function(pullCol,trainDf,supportTrain,compStat,roc, caseIndices, controlIndices, caseSplits, controlSplits){
  print("in parallelD")
  print(dim(trainDf))
  print(dim(supportTrain))
  print(pullCol)
  addCol <- supportTrain[,pullCol]
  addCol[addCol > 1] <- 1
  evalDf=cbind(trainDf,supp=factor(addCol))
  theStat <- rep(0,3)
  
  addColIndices <- which(addCol==1)
  if(length(addColIndices)>(length(caseIndices)*5)){addColIndices <- sample(addColIndices,(length(caseIndices)*5))}
  addColSplits <- ceiling(seq_along(addColIndices)/round(length(addColIndices)/3))
  
  for(i in 1:3){
    trainEvalDf <- evalDf[c(caseIndices[caseSplits==i], controlIndices[controlSplits==i], addColIndices[addColSplits==i]), ]
    testEvalDf <- evalDf[c(caseIndices[caseSplits!=i], controlIndices[controlSplits!=i], addColIndices[addColSplits!=i]),]
    
    completeLogistic <- glm(response ~ . ,data=trainEvalDf,family="binomial")
    testEvalDf$pred=predict(completeLogistic,testEvalDf,type="response")
    testEvalDf=testEvalDf[order(testEvalDf$pred,decreasing = T),]
    
    if(compStat=="OR"){
      exposeGroup=evalDf[1:round(nrow(evalDf)*0.05),ncol(trainDf)]
      safeGroup=evalDf[(round(nrow(evalDf)*0.05)+1):nrow(evalDf),ncol(trainDf)]
      theStat=(sum(exposeGroup==1)/sum(exposeGroup==0))/(sum(safeGroup==1)/sum(safeGroup==0))
    } else if (compStat=="Prevalance") {
      theStat = sum(head(evalDf,round(nrow(evalDf)/10))$response==1)/sum(tail(evalDf,round(nrow(evalDf)/10))$response==1)
    } else {
      g=pROC::roc(response ~ pred,data=testEvalDf)
      theStat[i]=as.numeric(g$auc)
    }
  }
  
  return(mean(theStat)) 
} 



#cl <- makeCluster(2,outfile="")
scoresStats=lapply(1:ncol(totalPhenoTrain),parallelD,trainDf,totalPhenoTrain,compStat,roc,
                   caseIndices, controlIndices, caseSplits, controlSplits)
#stopCluster(cl)
print("890")
scoresStats=unlist(scoresStats)
scoreStatsReturn=as.vector(scoresStats)
names(scoreStatsReturn) <- colnames(totalPhenoTrain)
if(any(scoreStatsReturn==1)){maxBest <- maxBest+1}

#if(length(scoreStatsReturn)>maxBest){
#  totalPhenoTrain=totalPhenoTrain[,scoresStats>=sort(scoresStats,decreasing = T)[maxBest] & scoreStatsReturn!=1]
#  totalPhenoTest=totalPhenoTest[,scoresStats>=sort(scoresStats,decreasing = T)[maxBest] & scoreStatsReturn!=1]
#}

#cl <- makeCluster(2,outfile="")
allResReturnDisease=lapply(1:ncol(totalPhenoTrain),parallelA,totalPhenoTrain,totalPhenoTest,trainDf,testDf,roc)
names(allResReturnDisease) <- colnames(totalPhenoTrain)
#stopCluster(cl)
print("go for color")
print(dim(totalPhenoTrain))
print(length(allResReturnDisease))
print(dim(suppTrain))
coding3 <- read.table("coding3.tsv",header=T,stringsAsFactors=F)
coding3[,2] <- strtrim(coding3[,2],25)
coding6 <- read.table("coding6.tsv",header=T,stringsAsFactors=F)
coding6[,2] <- strtrim(coding6[,2],25)
bothCoding <- rbind(coding6,coding3)
makePlotsColor(allResReturnDisease,totalPhenoTrain,bothCoding,8,3)
print("done color")

print("900")
scoresStatsKeep=scoresStats[scoresStats>=sort(scoresStats,decreasing = T)[maxBest] & scoresStats!=1]
print(head(scoresStatsKeep))
if(length(scoreStatsReturn)>maxBest){
  totalPhenoTrain=totalPhenoTrain[,scoresStats>=sort(scoresStats,decreasing = T)[maxBest] & scoreStatsReturn!=1]
  totalPhenoTest=totalPhenoTest[,scoresStats>=sort(scoresStats,decreasing = T)[maxBest] & scoreStatsReturn!=1]
}

print(length(scoresStatsKeep))
print(dim(totalPhenoTrain))
fullEvalDf=cbind(evalDf,totalPhenoTrain)
fullEvalDfTest=cbind(testDf,suppTest,totalPhenoTest)
print("904")
fullEvalDfTest$array=as.numeric(as.factor(fullEvalDfTest$array))
if(decoder[3]=="A"){
  fullEvalDfTest$sex=as.numeric(as.factor(fullEvalDfTest$sex))
}
print("909")
#coding3 <- read.table("coding3.tsv",header=T,stringsAsFactors=F)
#coding3[,2] <- strtrim(coding3[,2],25)
#coding6 <- read.table("coding6.tsv",header=T,stringsAsFactors=F)
#coding6[,2] <- strtrim(coding6[,2],25)
forPlot=data.frame(cbind(traits=colnames(totalPhenoTrain),scores=scoresStatsKeep),stringsAsFactors = F)

forPlot[forPlot$traits %in% coding6$coding,1]<-coding6[coding6$coding %in% forPlot$traits,2]
forPlot[forPlot$traits %in% coding3$coding,1]<-coding3[coding3$coding %in% forPlot$traits,2]
forPlot[,2]=as.numeric(forPlot[,2])-theStat

#forPlot=forPlot[forPlot$scores>sort(forPlot$scores,decreasing = T)[15],]
forPlot[,1]=factor(forPlot[,1],levels=forPlot[order(forPlot[,2]),1])
thePlot=ggplot(forPlot,aes(traits,scores))+geom_col()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
  labs(x="Added Trait",y="AUC Improvement",title="Traits - Training Application",subtitle=paste("For Trait",currentTrait))+
  guides(fill=FALSE)
plot(thePlot)

####################################################################################################################
#apply decision stump on each disease feature ######################################################################
print("391")
#the full model
fullRes=enetWrapper(fullEvalDf, bestLamb=NULL, justAuc = F,checkMissing = T, possIndex=1, makePlot = F,
                    addTitle = "With All Traits",isSplit=F,caseLimit=50,iter=0,enetRegression, enetWrapper)

#the stump model
decisionStump <- function(splitByCol, x, enetWrapper, enetRegression){
  print("startStump")
  print(splitByCol)
  library(glmnet)
  library(caret)
  library(pROC)
  
  splitBy=x[,splitByCol]
  x=x[,-splitByCol]
  x1=x[splitBy==0,]
  x2=x[splitBy==1,]
  
  write.table("b","test")
  auc1=enetWrapper(x1, bestLamb=NULL, justAuc = T,checkMissing = T, possIndex=1, makePlot = F,
                   addTitle="",isSplit=F,caseLimit=50,iter=0,enetRegression, enetWrapper)
  print("auc1")
  print(auc1[[1]])
  auc2=enetWrapper(x2, bestLamb=NULL, justAuc = T,checkMissing = T, possIndex=1, makePlot = F,
                   addTitle="",isSplit=F,caseLimit=50,iter=0,enetRegression, enetWrapper)
  print("auc2")
  print(auc2[[1]])
  returnVal=((auc1[[1]]*nrow(x1))+(auc2[[1]]*nrow(x2)))/nrow(x)
  print("returning")
  #print(returnVal)
  return(returnVal)
}

#get the disease columns that best split up the samples
#cl <- makeCluster(2,outfile="")
fullStumpStats=lapply((ncol(evalDf)+1):ncol(fullEvalDf), decisionStump, fullEvalDf, enetWrapper, enetRegression)
#stopCluster(cl)
fullStumpStats=t(matrix(unlist(fullStumpStats),nrow=3))
fullStumpStats[is.na(fullStumpStats)]=0.5
stumpStats=fullStumpStats[,2]
compNames1=colnames(totalPhenoTrain)
splitByCol1=which(colnames(fullEvalDf)==compNames1[which(stumpStats==max(stumpStats))])
splitter1=fullEvalDf[,splitByCol1]
nextEvalDf=fullEvalDf[splitter1==0,-splitByCol1]
compNames2=colnames(nextEvalDf)[(ncol(evalDf)+1):ncol(nextEvalDf)]
print("416")

#cl <- makeCluster(2,outfile="")
fullStumpStats2=lapply((ncol(evalDf)+1):ncol(nextEvalDf),decisionStump,nextEvalDf,enetWrapper, enetRegression)
#stopCluster(cl)
fullStumpStats2=t(matrix(unlist(fullStumpStats2),nrow=3))
fullStumpStats2[is.na(fullStumpStats2)]=0.5
stumpStats2=fullStumpStats2[,2]
splitByCol2=which(colnames(nextEvalDf)==compNames2[which(stumpStats2==max(stumpStats2))])
splitter2=nextEvalDf[,splitByCol2]
lastEvalDf=nextEvalDf[splitter2==0,-splitByCol2]
compNames3=colnames(lastEvalDf)[(ncol(evalDf)+1):ncol(lastEvalDf)]

print("422")
#cl <- makeCluster(2,outfile="")
fullStumpStats3=lapply((ncol(evalDf)+1):ncol(lastEvalDf),decisionStump,lastEvalDf,enetWrapper, enetRegression)
#stopCluster(cl)
fullStumpStats3=t(matrix(unlist(fullStumpStats3),nrow=3))
fullStumpStats3[is.na(fullStumpStats3)]=0.5
stumpStats3=fullStumpStats3[,2]
splitByCol3=which(colnames(lastEvalDf)==compNames3[which(stumpStats3==max(stumpStats3))])
splitter3=lastEvalDf[,splitByCol3]

#now report set these columns as fixed and evaluate each group
print("427")
fullEvalDf$assign=0 
fullEvalDf$assign[splitter1==1]=1
fullEvalDf$assign[splitter1==0][splitter2==1]=2
fullEvalDf$assign[splitter1==0][splitter2==0][splitter3==1]=3
fullEvalDf$assign[splitter1==0][splitter2==0][splitter3==0]=4

partRes=list()
aucRes=c()
for(i in 1:4){
  intoEnet = fullEvalDf[fullEvalDf$assign==i,]
  enetRes = enetWrapper(intoEnet[,-ncol(intoEnet)],bestLamb=NULL,justAuc=F,checkMissing=T,possIndex=1,makePlot=F,
                        addTitle="",isSplit=F, caseLimit=50, iter=0, enetRegression, enetWrapper)
  aucRes = rbind(aucRes,enetRes[[1]])
  partRes[[i]]=enetRes[[2]]
}
print("441")
splitAuc=(aucRes[1,]*sum(fullEvalDf$assign==1) + aucRes[2,]*sum(fullEvalDf$assign==2) +
            aucRes[3,]*sum(fullEvalDf$assign==3) + aucRes[4,]*sum(fullEvalDf$assign==4))/nrow(fullEvalDf)
fullRoc=partRes[[1]]*sum(fullEvalDf$assign==1) + partRes[[2]]*sum(fullEvalDf$assign==2) + 
  partRes[[3]]*sum(fullEvalDf$assign==3) + partRes[[4]]*sum(fullEvalDf$assign==4)
fullRoc = fullRoc/nrow(fullEvalDf)
print("447")

#do the actual plotting
colnames(fullRoc)=c("fpr","tpr")
print(colnames(fullRes[[2]]))
print(colnames(fullRoc))
forPlot=rbind(cbind(fullRes[[2]],type="fullApply"),cbind(fullRoc,type="stumpApply"))
print("450")
psAuc=round(splitAuc,3)
pfAuc=round(fullRes[[1]],3)
thePlot=ggplot(forPlot,aes(fpr,tpr,color=type))+geom_line()+
  labs(x="false positive rate",y="true positive rate",
       caption=paste("AUCs: split=",psAuc[2],"(",psAuc[1],"-",psAuc[3],")","~",
                     "full=",pfAuc[2],"(",pfAuc[1],"-",pfAuc[3],")"),
       title="Regularized LogReg with Traits - Training",subtitle=paste("For Score",currentTrait))+
  geom_abline(intercept = c(0,0),slope = 1)
plot(thePlot)

stumpVregAucs <- c(pfAuc,psAuc)


################################################################################################################################
# now apply best to testing data ###############################################################################################
#fullEvalDfTest$array=as.numeric(as.factor(fullEvalDfTest$array))
#fullEvalDfTest$sex=as.numeric(as.factor(fullEvalDfTest$sex))
print("460")
fullEvalDfTest=cbind(testDf,suppTest,totalPhenoTest)
fullEvalDfTest$array=as.numeric(as.factor(fullEvalDfTest$array))
if(decoder[3]=="A"){
  fullEvalDfTest$sex=as.numeric(as.factor(fullEvalDfTest$sex))
}
fullEvalDfTest <- fullEvalDfTest[,colnames(fullEvalDfTest) %in% colnames(fullEvalDf)]

#should change the enterWrapper and REgression so it can also take in set test files
if(fullRes[[1]][2] > 0){ #MADE A CHANGE HERE!
  #Now full works best
  fullEvalDf <- fullEvalDf[,colnames(fullEvalDf) %in% colnames(fullEvalDfTest)]  
  enetRes <- enetWrapper(fullEvalDf, bestLamb=NULL, justAuc = F,checkMissing = T, possIndex=1, makePlot = F,
              addTitle="",isSplit=fullEvalDfTest,caseLimit=50,iter=0,enetTest, enetWrapper)
  finalAuc <-enetRes[[1]]
  finalRoc <- enetRes[[2]]
  
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
    
    intoEnet = fullEvalDf[fullEvalDf$assign==i,]
    enetRes = enetWrapper(trainGroup[,-ncol(trainGroup)],bestLamb=NULL,justAuc=F,checkMissing=T,possIndex=1,makePlot=F,
                            addTitle="",isSplit=testGroup[,-ncol(testGroup)], caseLimit=50, iter=0, enetTest, enetWrapper)
    aucRes = rbind(aucRes,enetRes[[1]])
    partRes[[i]]=enetRes[[2]]
  }
    
  finalAuc=(aucRes[1,]*sum(fullEvalDf$assign==1) + aucRes[2,]*sum(fullEvalDf$assign==2) +
                aucRes[3,]*sum(fullEvalDf$assign==3) + aucRes[4,]*sum(fullEvalDf$assign==4))/nrow(fullEvalDf)
  finalRoc=partRes[[1]]*sum(fullEvalDf$assign==1) + partRes[[2]]*sum(fullEvalDf$assign==2) + 
    partRes[[3]]*sum(fullEvalDf$assign==3) + partRes[[4]]*sum(fullEvalDf$assign==4)
  finalRoc = fullRoc/nrow(fullEvalDf)
}

#ROC Curve of the testing data set
fpAuc=round(finalAuc,3)
addROC=testROC[testROC$group=="PRS+covars",]
scoreOnly=cbind(boostPlot,group="all_Scores")
forPlot=data.frame(finalRoc)
colnames(forPlot)=c("tpr","fpr")
forPlot$group="all_Info"
forPlot=forPlot[seq(1,nrow(forPlot),length.out = 1000),]
forPlot=rbind(forPlot,addROC,scoreOnly)
thePlot=ggplot(forPlot,aes(fpr,tpr,color=group))+geom_line()+geom_abline(intercept = c(0,0),slope = 1)+
  labs(x="False Positive Rate",y="True Positive Rate",title="Regularized LogReg with Scores and Traits - Testing",subtitle = paste("For Trait",currentTrait),
       caption=paste0("AUCs: all_Info=",fpAuc[2],"(",fpAuc[1],"-",fpAuc[3],")","\n",
                      " ~ PRS+covars=",testAUC[2],"(",testAUC[1],"-",testAUC[3],")","\n",
                      " ~ all_Scores=",firstBoostAuc[2],"(",firstBoostAuc[1],"-",firstBoostAuc[3],")"))
plot(thePlot)

#final comparison of the best eval test against best full test

dev.off()

return(list(sendBackBestTracker,suppTrain,suppTest,forPlot,finalAuc,bestScoreAUC,allResReturn,scoreStatsReturn,suppTrain,totalPhenoTrain,splitVnosplitAucs,stumpVregAucs,allResReturnDisease))
}
