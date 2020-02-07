library(gplots)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(Rtsne)
library(parallel)


corrEval <- function(authorComp,bestFile,bestTracker, trainDf, testDf, currentTrait, decoder,
                     suppTrain, suppTest, allSupportDecode, allDiseaseDecode,allScoresTrain, 
                     allScoresTest,testROC,testAUC,boostROC,boostAUC){
  
pdf(paste(authorComp,"corr","pdf",sep='.'))


#need to split all of the files into train/test for other disease/support files ########################################################
#######################################################################################################################################
splitBestTracker=str_split(bestTracker$bestFile,fixed("."))
typeOfFile=as.character(lapply(splitBestTracker,"[[",1))

supportTrain=suppTrain[,typeOfFile=="support"]
supportTest=suppTest[,typeOfFile=="support"]

diseaseTrain=suppTrain[,typeOfFile=="disease"]
diseaseTest=suppTest[,typeOfFile=="disease"]


#plot correlation between main score and each disease file ###############################################################
##############################################################################################################################
makeCorrMat <- function(score,compMat){
  totalCor=apply(compMat,2,function(x) cor(x,score,method="spearman"))
  bottomCor=apply(compMat,2,function(x) cor(x[score<summary(score)[2]],score[score<summary(score)[2]],method="spearman"))
  topCor=apply(compMat,2,function(x) cor(x[score>summary(score)[2]],score[score>summary(score)[2]],method="spearman"))
  return(rbind(totalCor,bottomCor,topCor))
}
#normalize the correlation by the correlation in the betas
betaCorr=read.table("betaCorr",header=T,row.names=1,stringsAsFactors = F)
possNames=rownames(betaCorr)
bestParam=as.character(lapply(splitBestTracker,"[[",4))
bestMethod=as.character(lapply(splitBestTracker,"[[",2))
bestName=rownames(bestTracker)
bestExtract=paste(bestName,bestParam,"final","set",bestMethod,"gz",sep=".")
authorExtract=paste(authorComp,str_split(bestFile,fixed("."))[[1]][4],"final","set",str_split(bestFile,fixed("."))[[1]][2],"gz",sep=".")

comboNames=sort(c(colnames(supportTrain),colnames(diseaseTrain)))
comboBetaCorr=rep(0,length(comboNames))
if(authorExtract %in% colnames(betaCorr)){
  betaCorr=betaCorr[rownames(betaCorr) %in% bestExtract, colnames(betaCorr) %in% authorExtract]
  betaCorr=betaCorr[order(possNames[possNames %in% bestExtract])]
  comboBetaCorr[bestExtract %in% possNames]=betaCorr
  isGood="Available"
} else {
  isGood="Not Available"
}

#the support files
supportCorMat <- as.data.frame(t(makeCorrMat(trainDf$predictor,supportTrain)))
allSupportDecode=allSupportDecode[order(allSupportDecode$author),]
supportCorMat$trait = allSupportDecode$trait[allSupportDecode$author %in% rownames(supportCorMat)]
supportCorMat$trait = factor(supportCorMat$trait,levels=supportCorMat$trait[order(supportCorMat$totalCor)])
supportCorMat$correct = supportCorMat$totalCor - comboBetaCorr[comboNames %in% rownames(supportCorMat)]

thePlot=ggplot(supportCorMat)+geom_col(aes(trait,totalCor))+
  geom_line(aes(trait,bottomCor,group=1,color="blue"))+
  geom_line(aes(trait,topCor,group=1,color="red"))+
  geom_line(aes(trait,correct,group=1,color="green"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
  labs(title="Correlation Among Non-Disease Traits",x="Traits",y="Correlation",color="quartiles",caption=paste("Correlation",isGood))+
  scale_color_manual(labels = c("bottom", "top","correct"), values = c("red", "blue","green"))
plot(thePlot)  

#the disease files
diseaseCorMat <- as.data.frame(t(makeCorrMat(trainDf$predictor,diseaseTrain)))
diseaseCorMat <- diseaseCorMat[order(rownames(diseaseCorMat)),]
allDiseaseDecode=allDiseaseDecode[order(allDiseaseDecode$author),]
diseaseCorMat$trait = allDiseaseDecode$disease[allDiseaseDecode$author %in% rownames(diseaseCorMat)]
diseaseCorMat$trait = factor(diseaseCorMat$trait,levels=diseaseCorMat$trait[order(diseaseCorMat$totalCor)])
diseaseCorMat$correct = diseaseCorMat$totalCor - comboBetaCorr[comboNames %in% rownames(diseaseCorMat)]

thePlot=ggplot(diseaseCorMat)+geom_col(aes(trait,totalCor))+
  geom_line(aes(trait,bottomCor,group=1,color="blue"))+
  geom_line(aes(trait,topCor,group=1,color="red"))+
  geom_line(aes(trait,correct,group=1,color="green"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
  labs(title="Correlation Among Disease Traits",x="Traits",y="Correlation",color="quartiles")+
  scale_color_manual(labels = c("bottom", "top","correct"), values = c("red", "blue","green"))
plot(thePlot)


####################################################################################################################################
#Logistic Regression for the current score against multiple diseases################################################################

totalCases=apply(allScoresTrain,2,sum)
nameNumberDecoder=read.table("nameNumberDecoder",stringsAsFactors = F)
diseaseNumber=as.numeric(strsplit(decoder[3],"-")[[1]][2])
compName=nameNumberDecoder[str_split(nameNumberDecoder[,2],"-",simplify = T)[,2]==diseaseNumber,1]
compName=gsub("/",".",compName)

allScoresTest=allScoresTest[,totalCases>500 | colnames(allScoresTrain)==compName]
allScoresTrain=allScoresTrain[,totalCases>500 | colnames(allScoresTrain)==compName]
diseasePval=rep(0,ncol(allScoresTrain)); diseaseAuc=rep(0,ncol(allScoresTrain))
names(diseasePval)=colnames(allScoresTrain); names(diseaseAuc)=colnames(allScoresTrain)

funky<-function(pullCol,allScoresTrain,allScoresTest,trainDf,testDf,roc){
  trainDf$newResponse=allScoresTrain[,pullCol]
  testDf$newResponse=allScoresTest[,pullCol]
  completeLogistic <- glm(newResponse ~ . - response,data=trainDf,family="binomial")
  covarLogistic <- glm(newResponse ~ age + sex + array + PC1 + PC2 + PC3 + PC4,data=trainDf,family="binomial")
  funcPval=coef(summary(completeLogistic))[9,4]
  testDf$prob=predict(completeLogistic,testDf,type="response")
  testDf$base=predict(covarLogistic,testDf,type="response")
  gScore=roc(newResponse ~ prob,data=testDf)
  gCovar=roc(newResponse ~ base,data=testDf)
  funcAuc=gScore$auc-gCovar$auc
  return(as.matrix(c(funcPval,funcAuc)))
}

cl <- makeCluster(2)
allRes=parLapply(cl,1:ncol(allScoresTrain),funky,allScoresTrain,allScoresTest,trainDf,testDf,roc)
allRes=unlist(allRes)
bothRes=matrix(allRes,ncol=ncol(allScoresTrain),nrow=2)
colnames(bothRes)=colnames(allScoresTrain)
diseasePval=bothRes[1,]
diseaseAuc=bothRes[2,]

#the compName is not within names(diseasePval)

refPval=diseasePval[names(diseasePval)==compName]
refAuc=diseaseAuc[names(diseaseAuc)==compName]
diseasePval=diseasePval[-which(names(diseasePval)==compName)]
diseaseAuc=diseaseAuc[-which(names(diseaseAuc)==compName)]
sigPval=sort(diseasePval)[1:5]
sigAuc=sort(diseaseAuc,decreasing = T)[1:5]

forPlot=data.frame(sigPval,disease=factor(names(sigPval),levels=names(sigPval)[order(sigPval)]))
thePlot=ggplot(forPlot,aes(disease,sigPval))+geom_col()+
  labs(title=paste("Model Fit of the",decoder[2],"Score Against Other Traits"),
       x="Disease Compared",y="P-Value",subtitle=paste("The",decoder[2],"PValue is",as.character(round(refPval,4))))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0,0,0,2),"cm"))
plot(thePlot)

forPlot=data.frame(sigAuc,disease=factor(names(sigAuc),levels=names(sigAuc)[order(sigAuc)]))
thePlot=ggplot(forPlot,aes(disease,sigAuc))+geom_col()+
  labs(x="Disease Compared",y="AUC",title=paste("Predictive Power of",decoder[2],"Score Against Other Traits"),
       subtitle=paste("The",decoder[2],"AUC is",as.character(round(refAuc,4))))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0,0,0,2),"cm"))
plot(thePlot)



################################################################################################################################
# How much does the presenece of each disease improve the prediction of the score under analysis ###############################
completeLogistic <- glm(response ~ array+PC1+PC2+PC3+PC4+sex+age+predictor,data=trainDf,family="binomial")
testDf$prob=predict(completeLogistic,testDf,type="response")
g=roc(response ~ prob,data=testDf)
refAuc=g$auc

funky<-function(pullCol,allScoresTrain,allScoresTest,trainDf,testDf,roc){
  trainDf$newResponse=allScoresTrain[,pullCol]
  testDf$newResponse=allScoresTest[,pullCol]
  scoreLogistic <- glm(response ~ array+PC1+PC2+PC3+PC4+sex+age+predictor*newResponse,data=trainDf,family="binomial")
  testDf$prob=predict(scoreLogistic,testDf,type="response")
  g=roc(response ~ prob,data=testDf)
  return(list(testDf$prob,g$auc))
}

allRes=parLapply(cl,1:ncol(allScoresTrain),funky,allScoresTrain,allScoresTest,trainDf,testDf,roc)
storeDiseaseProbs=data.frame(matrix(0,nrow=nrow(testDf),ncol=ncol(allScoresTrain)))
diseaseAuc=rep(0,ncol(allScoresTrain))
for(i in 1:length(allRes)){
  storeDiseaseProbs[,i]=allRes[[i]][[1]]
  diseaseAuc[i]=allRes[[i]][[2]]
}


forPlot=data.frame(disease=colnames(allScoresTrain),auc=diseaseAuc)
forPlot=forPlot[-which(forPlot$auc==1),]
forPlot=forPlot[forPlot$auc>sort(diseaseAuc,decreasing = T)[6],]
forPlot$auc=forPlot$auc-refAuc
forPlot$disease=factor(forPlot$disease,levels=forPlot$disease[order(forPlot$auc)])

thePlot=ggplot(forPlot,aes(disease,auc))+geom_col()+
  labs(title=paste("Improvement in Predicting",decoder[2]," by Adding Other Diseases to the Model"),x="Other Disease",y="AUC Improvement")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0,0,0,2),"cm"))
plot(thePlot)

colnames(storeDiseaseProbs)=colnames(allScoresTrain)
if(currentTrait %in% colnames(allScoresTrain)){
  storeDiseaseProbs=storeDiseaseProbs[,-which(colnames(storeDiseaseProbs) %in% currentTrait)]
}

prevHolder=matrix(0,nrow=100,ncol=ncol(storeDiseaseProbs))
colnames(prevHolder)=colnames(storeDiseaseProbs)
for(i in 1:(ncol(prevHolder))){
  perc=quantile(storeDiseaseProbs[,i],seq(0,100,1)/100)
  for(j in 1:100){
    prevHolder[j,i]=sum(storeDiseaseProbs[storeDiseaseProbs[,i]>perc[j] & storeDiseaseProbs[,i]<perc[j+1],ncol(storeDiseaseProbs)])
  }
}
prevHolder=as.data.frame(prevHolder)
topSum=apply(prevHolder[90:100,],2,sum)
keepNames=names(sort(topSum)[c(1:3,(ncol(prevHolder)-2):ncol(prevHolder))])
prevHolderPlot=prevHolder[,which(colnames(prevHolder) %in% keepNames)]
prevHolderPlot$index=1:100
toPlot=melt(prevHolderPlot,id.vars = "index")
thePlot=ggplot(toPlot,aes(index,value,color=variable))+geom_smooth(se=F)+
  labs(x="Score Percentile",y="Total Cases",title="Effect of Adding the Score to the Model")+
  scale_color_discrete(name = "Score Added")
plot(thePlot)


#############################################################################################################################
###############################################################################################################################
# analyze with the testing set the impact of each score on the currently analyzed score

funky<-function(pullCol,suppTrain,suppTest,trainDf,testDf){
  trainDf$newScore=suppTrain[,pullCol]
  testDf$newScore=suppTest[,pullCol]
  scoreLogistic <- glm(response ~ array+PC1+PC2+PC3+PC4+sex+age+predictor*newScore,data=trainDf,family="binomial")
  testDf$prob=predict(scoreLogistic,testDf,type="response")
  return(testDf$prob)
}

allRes=parLapply(cl,1:ncol(suppTrain),funky,suppTrain,suppTest,trainDf,testDf)
storeScoreProbs=data.frame(matrix(0,nrow=nrow(testDf),ncol=ncol(suppTrain)))
colnames(storeScoreProbs)=colnames(suppTrain)
for(i in 1:length(allRes)){
  storeScoreProbs[,i]=allRes[[i]]
}


prevHolder=matrix(0,nrow=100,ncol=ncol(storeScoreProbs))
temp=allDiseaseDecode[,1:2]
colnames(temp)=colnames(allSupportDecode)[1:2]
comboDecoder=rbind(allSupportDecode[,1:2],temp)
colnames(prevHolder)=comboDecoder$trait[comboDecoder$author %in% colnames(storeScoreProbs)]
for(i in 1:(ncol(prevHolder))){
  perc=quantile(storeScoreProbs[,i],seq(0,100,1)/100)
  for(j in 1:100){
    prevHolder[j,i]=sum(storeScoreProbs[storeScoreProbs[,i]>perc[j] & storeScoreProbs[,i]<perc[j+1],ncol(storeScoreProbs)])
  }
}
prevHolder=as.data.frame(prevHolder)
topSum=apply(prevHolder[90:100,],2,sum)
keepNames=names(sort(topSum)[c(1:3,(ncol(prevHolder)-2):ncol(prevHolder))])
prevHolderPlot=prevHolder[,which(colnames(prevHolder) %in% keepNames)]
prevHolderPlot$index=1:100
toPlot=melt(prevHolderPlot,id.vars = "index")
thePlot=ggplot(toPlot,aes(index,value,color=variable))+geom_smooth(se=F)+
  labs(x="Score Percentile",y="Total Cases",title="Effect of Adding the Disease to the Model")+
  scale_color_discrete(name = "Score Added")
plot(thePlot)

##########################################################################################################################
#do a very big enet logistic regression
totalTrain=cbind(trainDf[,1:9],suppTrain,allScoresTrain[,-which(colnames(allScoresTrain)==currentTrait)])
totalTrain$array=as.numeric(as.factor(totalTrain$array))
totalTrain$sex=as.numeric(as.factor(totalTrain$sex))

totalTest=cbind(testDf[,1:9],suppTest,allScoresTest[,-which(colnames(allScoresTest)==currentTrait)])
totalTest$array=as.numeric(as.factor(totalTest$array))
totalTest$sex=as.numeric(as.factor(totalTest$sex))

netCv = cv.glmnet(y=totalTrain$response,x=as.matrix(totalTrain[,-9]),family="binomial",nfolds = 3)
netLog = glmnet(y=totalTrain$response,x=as.matrix(totalTrain[,-9]),family="binomial",lambda = netCv$lambda.min)
bestPred=predict(netLog, newx = as.matrix(totalTest[,-9]))
g=roc(totalTest$response ~ bestPred[,1])
corrRoc=as.data.frame(cbind(g$sensitivities[seq(1,nrow(totalTest),length.out = 1000)],1-g$specificities[seq(1,nrow(totalTest),length.out = 1000)]))
colnames(corrRoc)=c("tpr","fpr")
testROC=testROC[testROC$group=="PRS+covars",1:2]
testROC[,2]=rev(testROC[,2])

print(head(corrRoc))
print(dim(corrRoc))
print(head(testROC))
print(dim(testROC))
print(head(boostROC))
print(dim(boostROC))
forPlot=rbind(cbind(corrRoc,type="enet"),cbind(testROC,type="logreg"),cbind(boostROC,type="boost"))

thePlot=ggplot(forPlot,aes(fpr,tpr,color=type))+geom_point()+
  labs(title="Regularized Logistic Regression With Scores and Disease Status",
       subtitle=paste("For Trait",currentTrait),
       caption=paste("AUC of",as.character(round(g$auc,4)),", boosted AUC of",round(boostAUC,4),", normal AUC of",round(testAUC,4)),
       y="true positive rate",x="false positive rate")+
  geom_abline(intercept = c(0,0),slope = 1)
plot(thePlot)


dev.off()
return(list(corrRoc,g$auc,bothRes,diseaseAuc))
#should I be normalizing by some other variable

}


###############################################################################################################################
###############################################################################################################################
# How does each disease effect the residuals of the primary analysis ##########################################################

# completeLogistic <- glm(response ~ array+PC1+PC2+PC3+PC4+sex+age+predictor,data=trainDf,family="binomial")
# trainDf$resids=completeLogistic$residuals
# quants=quantile(trainDf$resids,c(0.1,0.9))
# diseasePvalTotal=rep(0,ncol(allScoresTrain))
# diseasePvalXtreme=rep(0,ncol(allScoresTrain))
# for(i in 1:ncol(allScores)){
#   trainDf$newResponse=allScoresTrain[,i]
#   scoreLogistic <- glm(newResponse ~ resids,data=trainDf,family="binomial")
#   diseasePvalTotal[i]=coef(summary(scoreLogistic))[2,4]
#   scoreLogistic <- glm(response ~ resids,data=trainDf[trainDf$resids<quants[1] | trainDf$resids>quants[2],],family="binomial")
#   diseasePvalXtreme[i]=coef(summary(scoreLogistic))[2,4]
# }
# 
# goodPvals=which(diseasePvalTotal %in% sort(diseasePvalTotal)[1:10])
# forPlot=data.frame(total=diseasePvalTotal[goodPvals],lesser=diseasePvalXtreme[goodPvals],disease=colnames(allScoresTrain)[goodPvals])
