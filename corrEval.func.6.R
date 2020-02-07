library(gplots)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(Rtsne)
library(parallel)
library(stringr)
library(pROC)
library(cowplot)
#library(MendelianRandomization)



corrEval <- function(authorComp,bestFile,bestTracker, trainDf, testDf, currentTrait, decoder,
                     suppTrain, suppTest, allSupportDecode, allDiseaseDecode,allScoresTrain, 
                     allScoresTest, allSupportsTrain, allSupportsTest, testROC,testAUC,boostROC,boostAUC,mrScores,mrTraits){
  
pdf(paste(authorComp,"corr","pdf",sep='.'))


#need to split all of the files into train/test for other disease/support files ########################################################
#######################################################################################################################################
maxBest=25
scaledBest=apply(bestTracker[,-ncol(bestTracker)],1,max)
if(length(scaledBest)>maxBest){
  bestTracker=bestTracker[-which(scaledBest < sort(scaledBest,decreasing = T)[maxBest]),]
}

  
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
#supportCorMat$correct = supportCorMat$totalCor - comboBetaCorr[comboNames %in% rownames(supportCorMat)]

thePlot=ggplot(supportCorMat)+geom_col(aes(trait,totalCor))+
  geom_point(aes(trait,bottomCor,group=1,color="blue"))+
  geom_point(aes(trait,topCor,group=1,colour="red"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
  labs(title="Correlation Among Non-Disease Traits",x="Traits",y="Correlation",color="quartiles")+
  scale_color_manual(labels = c("bottom", "correct"), values = c("blue", "red"))
plot(thePlot)  

#the disease files
diseaseCorMat <- as.data.frame(t(makeCorrMat(trainDf$predictor,diseaseTrain)))
diseaseCorMat <- diseaseCorMat[order(rownames(diseaseCorMat)),]
allDiseaseDecode=allDiseaseDecode[order(allDiseaseDecode$author),]
diseaseCorMat$trait = allDiseaseDecode$disease[allDiseaseDecode$author %in% rownames(diseaseCorMat)]
diseaseCorMat$trait = factor(diseaseCorMat$trait,levels=diseaseCorMat$trait[order(diseaseCorMat$totalCor)])
#diseaseCorMat$correct = diseaseCorMat$totalCor - comboBetaCorr[comboNames %in% rownames(diseaseCorMat)]

thePlot=ggplot(diseaseCorMat)+geom_col(aes(trait,totalCor))+
  geom_point(aes(trait,bottomCor,group=1,color="blue"))+
  geom_point(aes(trait,topCor,group=1,color="red"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
  labs(title="Correlation Among Disease Traits",x="Traits",y="Correlation",color="quartiles")+
  scale_color_manual(labels = c("bottom", "correct"), values = c("red", "blue"))
plot(thePlot)


####################################################################################################################################
#Logistic Regression for the current score against multiple diseases################################################################

totalCases=apply(allScoresTrain,2,sum)
nameNumberDecoder=read.table("nameNumberDecoder",stringsAsFactors = F)

if (grepl(",",decoder[3])){
  decSplit=strsplit(decoder[3],",")[[1]][1]
  diseaseNumber=as.numeric(strsplit(decSplit,"-")[[1]][2])
} else {
  diseaseNumber=as.numeric(strsplit(decoder[3],"-")[[1]][2])
}

compName=nameNumberDecoder[str_split(nameNumberDecoder[,2],"-",simplify = T)[,2]==diseaseNumber,1]
compName=gsub("/",".",compName)

save(allScoresTrain,file="allScoresTrain")
save(allScoresTest,file="allScoresTest")
save(compName,file="compName")
allScoresTest=allScoresTest[,totalCases>500 | colnames(allScoresTrain)==compName]
allScoresTrain=allScoresTrain[,totalCases>500 | colnames(allScoresTrain)==compName]

totalCases=apply(allSupportsTrain,2,sum)
allSupportsTrain=allSupportsTrain[,totalCases>500]
allSupportsTest=allSupportsTest[,totalCases>500]


#can make all calculations unadjusted

parallelA<-function(pullCol,allTrain,allTest,trainDf,testDf,roc){
  library(pROC)
  #pval
  trainDf$newResponse=allTrain[,pullCol]
  testDf$newResponse=allTest[,pullCol]
  completeLogistic <- glm(newResponse ~ . - response,data=trainDf,family="binomial")
  covarLogistic <- glm(newResponse ~ age + sex + array + PC1 + PC2 + PC3 + PC4,data=trainDf,family="binomial")
  funcPval= -log10(coef(summary(completeLogistic))[9,4])
  
  #AUC
  testDf$prob=predict(completeLogistic,testDf,type="response")
  testDf$base=predict(covarLogistic,testDf,type="response")
  gScore=roc(newResponse ~ prob,data=testDf)
  gCovar=roc(newResponse ~ base,data=testDf)
  funcAuc=as.numeric(ci.auc(gScore))-as.numeric(ci.auc(gCovar))
  
  #OR
  testDf=testDf[order(testDf$predictor,decreasing = T),]
  allCutOffs=c(0.5,0.2,0.1,0.05,0.01,0.005)
  i=1
  funcOR=rep(0,length(allCutOffs))
  for(cutOff in allCutOffs){
    splitPoint=round(nrow(testDf)*cutOff)
    exposeGroup=testDf[1:splitPoint,10]
    safeGroup=testDf[(splitPoint+1):nrow(testDf),10]
    print(length(exposeGroup))
    funcOR[i]=(sum(exposeGroup==1)/sum(exposeGroup==0))/(sum(safeGroup==1)/sum(safeGroup==0))
    i=i+1
  }
  print(funcOR)
  
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
    labs(x="Score Percentile",y="Trait Prevalence",title=paste("Prevalences at",currentTrait,"Score Percentiles"))
  plot(thePlot)
}

#get the results
cl <- makeCluster(1)
allResScores=parLapply(cl,1:ncol(allScoresTrain),parallelA,allScoresTrain,allScoresTest,trainDf,testDf,roc)
makePlots(allResScores,allScoresTrain)

allResSupports=parLapply(cl,1:ncol(allSupportsTrain),parallelA,allSupportsTrain,allSupportsTest,trainDf,testDf,roc)
makePlots(allResSupports,allSupportsTrain)

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
