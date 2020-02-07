library(gplots)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(Rtsne)
library(parallel)
library(stringr)
library(pROC)
library(DescTools)
library(cowplot)
theme_set(theme_cowplot())
#library(MendelianRandomization)



corrEval <- function(authorComp,bestFile, trainDf, testDf, currentTrait, decoder,
                     allSupportDecode, allDiseaseDecode,allScoresTrain, 
                     allScoresTest, allSupportsTrain, allSupportsTest, testROC,testAUC,enetBetas){
 
print("start corr") 
pdf(paste(authorComp,"corr","pdf",sep='.'))
supportDefs <- read.table("supportDefs",header=T,stringsAsFactors=F)

#need to split all of the files into train/test for other disease/support files ########################################################
#######################################################################################################################################
#maxBest=10
#scaledBest=apply(bestTracker[,-ncol(bestTracker)],1,max)
#if(length(scaledBest)>maxBest){
#  bestTracker=bestTracker[-which(scaledBest < sort(scaledBest,decreasing = T)[maxBest]),]
#}

print("going grep")
if(grepl("Enet",bestFile)){
  takeMethod <- names(enetBetas)[enetBetas==max(enetBetas)]
} else {
  takeMethod <- bestFile
}

#if(grepl("2d",takeMethod)){
#  takeMethod <- sub("X",".",sub(".","-",sub(".","X",takeMethod,fixed=T),fixed=T))
#}

print(takeMethod)
print("making supp")
print(takeMethod)
suppTrain <- cbind(diseaseTrainList[[which(names(diseaseTrainList)==takeMethod)]],supportTrainList[[11]])
suppTest <- cbind(diseaseTestList[[which(names(diseaseTestList)==gsub("train","test",takeMethod))]],supportTestList[[11]])
print("made supp")


  
#splitBestTracker=str_split(bestTracker$bestFile,fixed("."))
#typeOfFile=as.character(lapply(splitBestTracker,"[[",1))
typeOfFile <- rep("support",ncol(suppTrain))
typeOfFile[colnames(suppTrain) %in% phenoDefs[,1]] <- "disease"

supportTrain=suppTrain[,typeOfFile=="support",drop=F]
supportTest=suppTest[,typeOfFile=="support",drop=F]

diseaseTrain=suppTrain[,typeOfFile=="disease",drop=F]
diseaseTest=suppTest[,typeOfFile=="disease",drop=F]
print("45")


if(decoder[3]!="A"){
  supportTrain <- supportTrain[covar1$sex==decoder[3],]
  supportTest <- supportTest[covar2$sex==decoder[3],]
  #supportTrain <- supportTrain[,-3]
  #supportTest <- supportTest[,-3]
  diseaseTrain <- diseaseTrain[covar1$sex==decoder[3],]
  diseaseTest <- diseaseTest[covar2$sex==decoder[3],]
  #diseaseTrain <- diseaseTrain[,-3]
  #diseaseTest <- diseaseTest[,-3]
}


if(is.na(diseaseTrain[nrow(diseaseTrain),1])){
  diseaseTrain <- rbind(apply(diseaseTrain[complete.cases(diseaseTrain),],2,mean), diseaseTrain[1:(nrow(diseaseTrain)-1),])
}


#plot correlation between main score and each disease file ###############################################################
##############################################################################################################################
makeCorrMat <- function(score,compMat){
  totalCor=apply(compMat,2,function(x) cor(x,score,method="spearman"))
  bottomCor=apply(compMat,2,function(x) cor(x[score<summary(score)[2]],score[score<summary(score)[2]],method="spearman"))
  topCor=apply(compMat,2,function(x) cor(x[score>summary(score)[2]],score[score>summary(score)[2]],method="spearman"))
  totalPval=apply(compMat,2,function(x) cor.test(x,score,method="spearman")$p.value)
  bottomPval=apply(compMat,2,function(x) cor.test(x[score<summary(score)[2]],score[score<summary(score)[2]],method="spearman")$p.value)
  topPval=apply(compMat,2,function(x) cor.test(x[score>summary(score)[2]],score[score>summary(score)[2]],method="spearman")$p.value)
  lowCi=apply(compMat,2,function(x) SpearmanRho(x,score,conf.level =0.95)[2])
  hiCi=apply(compMat,2,function(x) SpearmanRho(x,score,conf.level =0.95)[3])
  return(rbind(totalCor,bottomCor,topCor,totalPval,bottomPval,topPval,lowCi,hiCi))
}

print("55")
print(head(trainDf))
print(head(supportTrain))
print(dim(supportTrain))
print(dim(trainDf))
#the support files
if(length(supportTrain)>0){
supportCorMat <- as.data.frame(t(makeCorrMat(trainDf$predictor,supportTrain)))
print("59")
print(dim(supportCorMat))
print(dim(supportTrain))
allSupportDecode=supportDefs[order(supportDefs$author),1:2]
print("61")
print(dim(allSupportDecode))
print(allSupportDecode)
print(rownames(supportCorMat))
colnames(allSupportDecode)[2] <- "trait"
supportCorMat$trait = allSupportDecode$trait[allSupportDecode$author %in% rownames(supportCorMat)]
print("62")
print(supportCorMat$trait)
print(length(supportCorMat$trait))
print(length(unique(supportCorMat$trait)))
print(length(supportCorMat$totalCor))
supportCorMat$trait = factor(supportCorMat$trait,levels=supportCorMat$trait[order(supportCorMat$totalCor)])
print("63")

thePlot=ggplot(supportCorMat,aes(trait,totalCor))+geom_point()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
  labs(title="Correlation Among Non-Disease Traits",x="Traits",y="Correlation")+
  geom_errorbar(aes(ymin=lowCi,ymax=hiCi),width=0.4)
plot(thePlot)

supportCorMat$trait = factor(supportCorMat$trait,levels=supportCorMat$trait[order(supportCorMat$totalPval)])
thePlot=ggplot(supportCorMat)+
  geom_point(aes(trait,-log10(totalPval),group=1,color="blue"))+
  geom_point(aes(trait,-log10(bottomPval),group=1,color="green"))+
  geom_point(aes(trait,-log10(topPval),group=1,color="red"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
  labs(title="Correlation Among Non-Disease Traits",x="Traits",y="Correlation -log10(PVal)",color="quartiles")+
  scale_color_manual(labels = c("total","bottom","top"), values = c("blue","red","green"))
plot(thePlot)
} else {
  supportCorMat <- NULL
}

print("84")
#the disease files
if(length(diseaseTrain)>0){
diseaseCorMat <- as.data.frame(t(makeCorrMat(trainDf$predictor,diseaseTrain)))
diseaseCorMat <- diseaseCorMat[order(rownames(diseaseCorMat)),]
allDiseaseDecode=phenoDefs[order(phenoDefs$author),1:2]
colnames(allDiseaseDecode)[2] <- "disease"
diseaseCorMat$trait = allDiseaseDecode$disease[allDiseaseDecode$author %in% rownames(diseaseCorMat)]
diseaseCorMat$trait = factor(diseaseCorMat$trait,levels=unique(diseaseCorMat$trait[order(diseaseCorMat$totalCor)]))

thePlot=ggplot(diseaseCorMat,aes(trait,totalCor))+geom_point()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
  labs(title="Correlation Among Disease Traits",x="Traits",y="Correlation")+
  geom_errorbar(aes(ymin=lowCi,ymax=hiCi),width=0.4)
plot(thePlot)

diseaseCorMat$trait = factor(diseaseCorMat$trait,levels=unique(diseaseCorMat$trait[order(diseaseCorMat$totalPval)]))
thePlot=ggplot(diseaseCorMat)+
  geom_point(aes(trait,-log10(totalPval),group=1,color="blue"))+
  geom_point(aes(trait,-log10(bottomPval),group=1,color="green"))+
  geom_point(aes(trait,-log10(topPval),group=1,color="red"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
  labs(title="Correlation Among Disease Traits",x="Traits",y="Correlation -log10(PVal)",color="quartiles")+
  scale_color_manual(labels = c("total","bottom","top"), values = c("blue","red","green"))
plot(thePlot)
} else {
diseaseCorMat <- NULL
}

print("113")

####################################################################################################################################
#Logistic Regression for the current score against multiple diseases################################################################

coding6 <- read.table("coding6.tsv",header=T,stringsAsFactors=F)
coding6[,2] <- strtrim(coding6[,2],25)
allScoresTrain <- corrTrainRefs[["noncancer"]] #this will change the definition of compName
allScoresTest <- corrTestRefs[["noncancer"]]

for(i in 1:ncol(allScoresTrain)){
  allScoresTrain[,i][allScoresTrain[,i] > 1] <- 1
}
for(i in 1:ncol(allScoresTest)){
  allScoresTest[,i][allScoresTest[,i] > 1] <- 1
}
allScoresTrain <- allScoresTrain[,colnames(allScoresTrain) %in% colnames(allScoresTest)]
allScoresTest <- allScoresTest[,colnames(allScoresTest) %in% colnames(allScoresTrain)]
totalCasesTrain=apply(allScoresTrain,2,sum)
totalCasesTest=apply(allScoresTest,2,sum)
compName=str_split(decoder[5],",")[[1]]

allScoresTest=allScoresTest[,totalCasesTrain>500 & totalCasesTest>500 & !(colnames(allScoresTrain) %in% compName)]
allScoresTrain=allScoresTrain[,totalCasesTrain>500 & totalCasesTest>500 & !(colnames(allScoresTrain) %in% compName)]
colnames(allScoresTrain) <- coding6[coding6[,1] %in% colnames(allScoresTrain),2]
colnames(allScoresTest) <- coding6[coding6[,1] %in% colnames(allScoresTest),2]

if(decoder[3] != "A"){
  allScoresTrain <- allScoresTrain[covar1[,3]==decoder[3],]
  allScoresTest <- allScoresTest[covar2[,3]==decoder[3],]
  takeOut <- unique(c(which(apply(allScoresTrain,2,function(y) length(unique(y)))==1),
               which(apply(allScoresTest,2,function(y) length(unique(y)))==1)))
  allScoresTrain <- allScoresTrain[,-takeOut]
  allScoresTest <- allScoresTest[,-takeOut]
}


print("144")
coding19 <- read.table("coding19.tsv",header=T,stringsAsFactors=F)
coding19[,2] <- strtrim(coding19[,2],25)
allSupportsTrain <- corrTrainRefs[["icd10"]]
allSupportsTest <- corrTestRefs[["icd10"]]
for(i in 1:ncol(allSupportsTrain)){
  allSupportsTrain[,i][allSupportsTrain[,i] > 1] <- 1
}
for(i in 1:ncol(allSupportsTest)){
  allSupportsTest[,i][allSupportsTest[,i] > 1] <- 1
}
allSupportsTrain <- allSupportsTrain[,colnames(allSupportsTrain) %in% colnames(allSupportsTest)]
allSupportsTest <- allSupportsTest[,colnames(allSupportsTest) %in% colnames(allSupportsTrain)]
totalCasesTrain=apply(allSupportsTrain,2,sum)
totalCasesTest=apply(allSupportsTest,2,sum)
compName=str_split(decoder[7],",")[[1]]

print("165")
print(head(coding19))
print(head(colnames(allSupportsTrain)))

totalCases=apply(allSupportsTrain,2,sum)
allSupportsTrain=allSupportsTrain[,totalCasesTrain>1500 & totalCasesTest>1500 & !(colnames(allSupportsTrain) %in% compName)]
allSupportsTest=allSupportsTest[,totalCasesTest>1500 & totalCasesTest>1500 & !(colnames(allSupportsTest) %in% compName)]
print("228")
colnames(allSupportsTrain) <- coding19[coding19[,1] %in% colnames(allSupportsTrain),2]
colnames(allSupportsTest) <- coding19[coding19[,1] %in% colnames(allSupportsTest),2]

if(decoder[3] != "A"){
  allSupportsTrain <- allSupportsTrain[covar1[,3]==decoder[3],]
  allSupportsTest <- allSupportsTest[covar2[,3]==decoder[3],]
  takeOut <- unique(c(which(apply(allSupportsTrain,2,function(y) length(unique(y)))==1),
               which(apply(allSupportsTest,2,function(y) length(unique(y)))==1)))
  allSupportsTrain <- allSupportsTrain[,-takeOut]
  allSupportsTest <- allSupportsTest[,-takeOut]
}


#can make all calculations unadjusted
parallelA<-function(pullCol,allTrain,allTest,trainDf,testDf,roc){
  print("start")
  print(pullCol)
  library(pROC)

  print("pval")
  #pval
  trainDf$newResponse=allTrain[,pullCol]
  testDf$newResponse=allTest[,pullCol]
  print(table(allTrain[,pullCol]))
  print(table(allTest[,pullCol]))
  completeLogistic <- glm(newResponse ~ . - response,data=trainDf,family="binomial")
  if("sex" %in% colnames(trainDf)){
    covarLogistic <- glm(newResponse ~ age + sex + array + PC1 + PC2 + PC3 + PC4,data=trainDf,family="binomial")
  } else {
    covarLogistic <- glm(newResponse ~ age + array + PC1 + PC2 + PC3 + PC4,data=trainDf,family="binomial")
  }
  sumCoef <- coef(summary(completeLogistic))
  funcPval= -log10(sumCoef[nrow(sumCoef),4])
  print(funcPval)
  
  print("auc")
  #AUC
  testDf$prob=predict(completeLogistic,testDf,type="response")
  testDf$base=predict(covarLogistic,testDf,type="response")
  print("halfAUC")
  gScore=pROC::roc(newResponse ~ prob,data=testDf)
  gCovar=pROC::roc(newResponse ~ base,data=testDf)
  funcAuc=as.numeric(ci.auc(gScore))-as.numeric(ci.auc(gCovar))
  print(funcAuc)  

  print("or")
  #OR
  testDf=testDf[order(testDf$predictor,decreasing = T),]
  allCutOffs=c(0.5,0.2,0.1,0.05,0.01,0.005)
  i=1
  funcOR=rep(0,length(allCutOffs))
  for(cutOff in allCutOffs){
    print("cutOff")
    splitPoint=round(nrow(testDf)*cutOff)
    exposeGroup=testDf[1:splitPoint,which(colnames(testDf)=="newResponse")]
    safeGroup=testDf[(splitPoint+1):nrow(testDf),which(colnames(testDf)=="newResponse")]
    print(sum(exposeGroup==0))
    print(sum(exposeGroup==1))
    print(sum(safeGroup==0))
    print(sum(safeGroup==1))
    funcOR[i]=(sum(exposeGroup==1)/sum(exposeGroup==0))/(sum(safeGroup==1)/sum(safeGroup==0))
    i=i+1
  }
  print(funcOR)
  
  print("prev")
  #Prev
  funcPrev=rep(0,100)
  testDf$group=rep(0,nrow(testDf))
  groupRanges=c(round(seq(1,nrow(testDf),by = nrow(testDf)/100)),nrow(testDf))
  for(i in 1:100){
    testDf$group[groupRanges[i]:groupRanges[i+1]]=i
  }
  testDf$group=rev(testDf$group)
  for(i in 1:100){
    withDis=sum(testDf[testDf$group==i,10]==1)
    totalGroup=sum(testDf$group==i)
    funcPrev[i]=withDis/totalGroup
  }
  print(head(funcPrev))
  
  return(list(funcPval,funcAuc,funcOR,funcPrev))
}



makePlots <- function(allRes,trainScores){
  pvals=unlist(lapply(allRes, `[[`, 1))
  aucImps=unlist(lapply(allRes, `[[`, 2))
  aucImps=t(matrix(aucImps,nrow=3))
  ors=t(matrix(unlist(lapply(allRes, `[[`, 3)),nrow=6))
  prevs=t(matrix(unlist(lapply(allRes, `[[`, 4)),nrow=100))
  
  #PVALS
  fullPvalsDf=data.frame(pvals,traits=colnames(trainScores),stringsAsFactors = F)
  if(compName %in% fullPvalsDf[,2]){
    pvalsDf=fullPvalsDf[-which(fullPvalsDf[,2]==compName),]
  } else {
    pvalsDf=fullPvalsDf
  }
  pvalsDf=pvalsDf[order(pvalsDf[,1],decreasing = F),]
  pvalsDf=pvalsDf[1:8,]
  pvalsDf[,2]=factor(pvalsDf[,2],levels = pvalsDf[order(pvalsDf[,1]),2])
  thePlot=ggplot(pvalsDf,aes(traits,-log10(pvals)))+geom_col(aes(fill=traits))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
    labs(x="traits",y="-log10(logistic regression PVal)",title=paste("Pval of",currentTrait,"Score"))+
    guides(fill=FALSE)
  plot(thePlot)
  
  #OR vs PVALS
  fullOrPvalDf=data.frame(pvals,ors=ors[,4],traits=colnames(trainScores),stringsAsFactors = F)
  orPvalDf=data.frame(fullOrPvalDf)
  orPvalDf$traits[orPvalDf$pvals > sort(orPvalDf[,1],decreasing = F)[3] & orPvalDf$ors < sort(orPvalDf[,2],decreasing = T)[3]]=""
  thePlot=ggplot(orPvalDf,aes(-log10(pvals),ors,label=traits))+geom_point()+geom_label_repel(size=3,force=3)+
    labs(x="-log10(PValue)",y="Odds Ratio",title=paste("Effect and Certainty for",currentTrait,"Score"))
  plot(thePlot)
  
  #AUC
  fullAucImpsDf=data.frame(aucImps,traits=colnames(trainScores),stringsAsFactors = F)
  refAuc=fullAucImpsDf[which(fullAucImpsDf[,4]==compName),2]
  if(compName %in% fullAucImpsDf[,4]){
    aucImpsDf=fullAucImpsDf[-which(fullAucImpsDf[,4]==compName),]
  } else {
    aucImpsDf=fullAucImpsDf
  }
  aucImpsDf=aucImpsDf[order(aucImpsDf[,2],decreasing = T),]
  aucImpsDf=aucImpsDf[1:8,]
  aucImpsDf[,4]=factor(aucImpsDf[,4],levels = aucImpsDf[order(aucImpsDf[,2]),4])
  colnames(aucImpsDf)=c("lo","aucImps","hi","traits")
  thePlot=ggplot(aucImpsDf,aes(traits,aucImps))+geom_col(aes(fill=traits))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
    labs(x="Traits",y="AUC Improvement",title=paste("AUC Improvement of",currentTrait,"Score"))+
    guides(fill=FALSE)+geom_errorbar(aes(ymin=lo,ymax=hi),width=0.2)
  plot(thePlot)
  
  #OR
  colnames(ors)=paste("cut0ff",c(0.5,0.2,0.1,0.05,0.01,0.005),sep=".")
  fullOrsDf=data.frame(ors,traits=colnames(trainScores),stringsAsFactors = F)
  orsDf=fullOrsDf[order(rowSums(fullOrsDf[,1:6]),decreasing = T)[1:8],]
  orsDf=melt(orsDf,id.vars="traits")
  orsDf[,2]=sort(rep(c(0.5,0.2,0.1,0.05,0.01,0.005),8),decreasing = T)
  orsDf[,2]=factor(orsDf[,2],levels=c(0.5,0.2,0.1,0.05,0.01,0.005))
  thePlot=ggplot(orsDf,aes(variable,value,color=traits))+geom_line(aes(group=traits))+
    labs(x="Percentile Cut-Off",y="Odds Ratio",title=paste("Odds ratios of",currentTrait,"score"))
  plot(thePlot)
  
  #PREV
  colnames(prevs)=1:100/100
  fullPrevsDf=data.frame(prevs,traits=colnames(trainScores),stringsAsFactors = F)
  prevsDf=fullPrevsDf[order(rowSums(fullPrevsDf[,91:100])/rowSums(fullPrevsDf[,1:8]),decreasing = T)[1:8],]
  prevsDf=melt(prevsDf,id.vars = "traits")
  prevsDf[,2]=sort(rep((1:100/100),8))
  thePlot=ggplot(prevsDf,aes(variable,value,color=traits))+geom_smooth(se=F)+
    labs(x="Score Percentile",y="Trait Prevalence",title=paste("Prevalences at",currentTrait,"Score Percentiles"))+
    guides(color=FALSE)
  plot(thePlot)
}


makePlotsColor <- function(allRes,trainScores,totalTraits,numberColor){
  pvals=unlist(lapply(allRes, `[[`, 1))
  aucImps=unlist(lapply(allRes, `[[`, 2))
  aucImps=t(matrix(aucImps,nrow=3))
  ors=t(matrix(unlist(lapply(allRes, `[[`, 3)),nrow=6))
  prevs=t(matrix(unlist(lapply(allRes, `[[`, 4)),nrow=100))
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  specColors=c(gg_color_hue(numberColor),rep("#808080",totalTraits-numberColor))
  
  #AUC
  print("auc")
  fullAucImpsDf=data.frame(aucImps,traits=colnames(trainScores),stringsAsFactors = F)
  refAuc=fullAucImpsDf[which(fullAucImpsDf[,4]==compName),2]
  if(any(compName %in% fullAucImpsDf[,4])){
    aucImpsDf=fullAucImpsDf[-which(fullAucImpsDf[,4] %in% compName),]
  } else {
    aucImpsDf=fullAucImpsDf
  }
  aucImpsDf=aucImpsDf[order(aucImpsDf[,2],decreasing = T),]
  aucImpsDf=aucImpsDf[1:totalTraits,]
  names(specColors)=aucImpsDf$traits
  print(aucImpsDf)
  aucImpsDf[,4]=factor(aucImpsDf[,4],levels = unique(aucImpsDf[order(aucImpsDf[,2]),4]))
  colnames(aucImpsDf)=c("lo","aucImps","hi","traits")
  thePlot=ggplot(aucImpsDf,aes(traits,aucImps))+geom_point(aes(color=traits), size=3)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
    labs(x="Traits",y="AUC Improvement",title=paste("AUC Improvement of",currentTrait,"Score"))+
    guides(fill=FALSE)+geom_errorbar(aes(ymin=lo,ymax=hi),width=0.2)+
    scale_fill_manual(values=specColors)+guides(color=FALSE)
  plot(thePlot)
  
  #PVALS
  print("pvals")
  fullPvalsDf=data.frame(pvals,traits=colnames(trainScores),stringsAsFactors = F)
  pvalsDf=fullPvalsDf[fullPvalsDf$traits %in% as.character(aucImpsDf$traits),]
  print(pvalsDf)
  pvalsDf[,2]=factor(pvalsDf[,2],levels = unique(pvalsDf[order(pvalsDf[,1]),2]))
  thePlot=ggplot(pvalsDf,aes(traits,-log10(pvals),fill=traits))+geom_col()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
    labs(x="traits",y="-log10(logistic regression PVal)",title=paste("Pval of",currentTrait,"Score"))+
    scale_fill_manual(values=specColors)+guides(fill=FALSE)
  plot(thePlot)
  
  #OR vs PVALS
  print("vs")
  fullOrPvalDf=data.frame(pvals,ors=ors[,4],traits=colnames(trainScores),stringsAsFactors = F)
  orPvalDf=data.frame(fullOrPvalDf)
  orPvalDf$words=orPvalDf$traits
  orPvalDf$words[orPvalDf$pvals > sort(orPvalDf[,1],decreasing = F)[3] & orPvalDf$ors < sort(orPvalDf[,2],decreasing = T)[3]]=""
  thePlot=ggplot(orPvalDf,aes(-log10(pvals),ors,label=words))+
    geom_point(data=orPvalDf[orPvalDf$traits %in% names(specColors[1:3]),],aes(color=traits),size=4)+
    geom_point(data=orPvalDf[!(orPvalDf$traits %in% names(specColors[1:3])),])+
    geom_label_repel(size=3,force=3)+
    labs(x="-log10(PValue)",y="Odds Ratio",title=paste("Effect and Certainty for",currentTrait,"Score"))+
    theme(legend.position="bottom",legend.direction="vertical")+
    scale_color_manual(values=specColors[1:3])
  plot(thePlot)
  
  
  #OR
  print("or")
  colnames(ors)=paste("cut0ff",c(0.5,0.2,0.1,0.05,0.01,0.005),sep=".")
  fullOrsDf=data.frame(ors,traits=colnames(trainScores),stringsAsFactors = F)
  print(fullOrsDf)
  orsDf=fullOrsDf[fullOrsDf$traits %in% as.character(aucImpsDf$traits),]
  orsDf=melt(orsDf,id.vars="traits")
  orsDf[,2]=sort(rep(c(0.5,0.2,0.1,0.05,0.01,0.005),(nrow(orsDf)/6),decreasing = T))
  orsDf[,2]=factor(orsDf[,2],levels=c(0.5,0.2,0.1,0.05,0.01,0.005))
  thePlot=ggplot(orsDf,aes(variable,value,color=traits))+geom_line(aes(group=traits))+
    labs(x="Percentile Cut-Off",y="Odds Ratio",title=paste("Odds ratios of",currentTrait,"score"))+
    scale_color_manual(values=specColors)+guides(color=FALSE)
  plot(thePlot)
  
  #PREV
  print("prev")
  colnames(prevs)=1:100/100
  fullPrevsDf=data.frame(prevs,traits=colnames(trainScores),stringsAsFactors = F)
  prevsDf=fullPrevsDf[fullPrevsDf$traits %in% as.character(aucImpsDf$traits),]
  prevsDf=melt(prevsDf,id.vars = "traits")
  prevsDf[,2]=sort(rep((1:100/100),nrow(prevsDf)/100))
  thePlot=ggplot(prevsDf,aes(variable,value,color=traits))+geom_smooth(se=F)+
    labs(x="Score Percentile",y="Trait Prevalence",title=paste("Prevalences at",currentTrait,"Score Percentiles"))+
    guides(color=FALSE)+scale_color_manual(values=specColors)
  plot(thePlot)
}


#get the results
print("make scores")
#cl <- makeCluster(2,outfile="")
allResScores=lapply(1:ncol(allScoresTrain),parallelA,allScoresTrain,allScoresTest,trainDf,testDf,roc)
names(allResScores) <- colnames(allScoresTrain)
#stopCluster(cl)
#makePlots(allResScores,allScoresTrain)
print("makeplotscolor1")
makePlotsColor(allResScores,allScoresTrain,8,3)
allResScores[[length(allResScores)+1]] <- colnames(allScoresTrain)

print("make supports")
#cl <- makeCluster(2,outfile="")
allResSupports=lapply(1:ncol(allSupportsTrain),parallelA,allSupportsTrain,allSupportsTest,trainDf,testDf,roc)
names(allResSupports) <- colnames(allSupportsTrain)
#stopCluster(cl)
#makePlots(allResSupports,allSupportsTrain)
print("makeplotscolor2")
makePlotsColor(allResSupports,allSupportsTrain,8,3)
allResSupports[[length(allResSupports)+1]] <- colnames(allSupportsTrain)

#MR
#mrDf=data.frame(trainDf)
#scores=cbind(mrScores,base=trainDf$predictor)
#traits=cbind(mrTraits,base=trainDf$response)

#betasDf=data.frame(matrix(0,nrow=ncol(scores),ncol=ncol(traits)))
#seDf=data.frame(matrix(0,nrow=ncol(scores),ncol=ncol(traits)))
#colnames(betasDf)=colnames(traits)
#colnames(seDf)=colnames(traits)
#rownames(betasDf)=colnames(scores)
#rownames(betasDf)=colnames(scores)
#for(i in 1:ncol(scores)){
#  mrDf$newScore=scores[,i]
#  for(j in 1:ncol(traits)){
#    mrDf$newTrait=traits[,j]
#    logit=glm(newTrait ~ age + sex + array + PC1 + PC2 + PC3 + PC4 + newScore,data=mrDf,family="binomial")
#    betasDf[i,j]=coef(summary(logit))[9,1]
#    seDf[i,j]=coef(summary(logit))[9,2]
#  }
#}

#corrMat=matrix(0,nrow=ncol(scores),ncol=ncol(scores))
#for(i in 1:ncol(scores)){
#  corrMat[i,]=apply(scores,2,function(x) cor(x,scores[,i]))
#}

#baseIndex=which(colnames(betasDf)=="base")
#mrIn=mr_mvinput(bx=as.matrix(betasDf[,-baseIndex]),bxse = as.matrix(seDf[,-baseIndex]),
#                by = betasDf[,baseIndex],byse = seDf[,baseIndex],
#                correlation = corrMat)
#mrOut=mr_mvivw(mrIn)
#mrPlot=data.frame(mrOut@Estimate,mrOut@CILower,mrOut@CIUpper,colnames(traits)[-baseIndex],stringsAsFactors = F)
#colnames(mrPlot)=c("est","lo","hi","trait")
#mrPlot$trait=factor(mrPlot$trait,levels=mrPlot$trait[order(mrPlot$est)])

#ggplot(mrPlot,aes(trait,est))+geom_point()+geom_errorbar(aes(ymin=lo,ymax=hi))+
#  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
#  labs(title="Mendelian Randomization Test",xlab="Trait",ylab="Estimate",
#       caption=paste("Pval of Pleiotropy",round(mrOut@Heter.Stat[2],4)))

dev.off()
return(list(supportCorMat,diseaseCorMat,allResScores,allResSupports))


}
