library(fmsb)
library(PRROC)
library(ggplot2)
library(forcats)
library(pROC)
library(reshape2)
library(stringr)
library(unbalanced)


phenoEval <- function(authorComp,bestFile,trainDf,currentTrait,decoder,pheno1,pheno2,diseaseTestList,toRev){
  
pdf(paste(authorComp,"pheno","pdf",sep='.'))

prepWork=function(scores,pheno){
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
  df$predictor=as.numeric(scale(df$predictor))
  return(df)
}

nameDiseaseFiles=list.files("diseaseFiles/","train")
sentInScores=diseaseTestList[[which(nameDiseaseFiles==bestFile)]]
sentOutTestDf=prepWork(sentInScores,pheno2)
if(mean(sentOutTestDf$predictor[sentOutTestDf$response==1]) < (mean(sentOutTestDf$predictor[sentOutTestDf$response==0]))){
  sentOutTestDf$predictor=sentOutTestDf$predictor*-1
}
testDf=sentOutTestDf[order(sentOutTestDf$predictor,decreasing=T),]

#Density plots
thePlot=ggplot(testDf,aes(x=predictor,fill=response))+geom_density(alpha=0.3)+
  labs(y="density",x="raw score",title=paste("Raw Score Densities for",currentTrait))+
  scale_fill_discrete(name="Disease Status",breaks=c(0,1),labels=c("No Disease","Disease"))
plot(thePlot)

#Assess best sampling method
newTrainDf=data.frame(trainDf)
newTrainDf[,1]=as.numeric(as.factor(newTrainDf[,1]))
newTrainDf[,6]=as.numeric(as.factor(newTrainDf[,6]))
aucVec=matrix(0,nrow=3,ncol=6)
j=1
for(ubMeth in c("ubOver","ubUnder","ubSMOTE","ubNCL","ubTomek","nothing")){
  print(ubMeth)
  for(i in 1:3){
    sampleAmount=round(nrow(trainDf)/2)
    sampleSet=sort(sample(1:nrow(trainDf),sampleAmount))
    fullSet=1:nrow(trainDf)
    otherSet=fullSet[!(fullSet %in% sampleSet)]
    train=newTrainDf[sampleSet,]
    test=newTrainDf[otherSet,]
  
    if(ubMeth=="nothing"){
      data=list(X=train[,-9],Y=train[,9])
    } else {
      data <- ubBalance(X=train[,-9], Y=train[,9], type=ubMeth) 
    }
    outputDf=data.frame(data$X, data$Y)
    colnames(outputDf)=colnames(trainDf)
    covarLogistic <- glm(response ~ age + sex + array + PC1 + PC2 + PC3 + PC4,data=outputDf,family="binomial")
    testProb=predict(covarLogistic,test,type="response")
    test$prob=testProb
    g=roc(response ~ prob,data=test)
    print(g$auc)
    test=test[,-ncol(test)]
    aucVec[i,j]=g$auc
  }
  j=j+1
}
samplingPvals=apply(aucVec[,-6],2,function(x) wilcox.test(x,aucVec[,6])$p.value)
if(any(samplingPvals<0.1)){
  bestMethod=c("ubOver","ubUnder","ubSMOTE","ubNCL","ubTomek")[samplingPvals==min(samplingPvals)]
} else {
  bestMethod="nothing"
}

#Violin Plot
completeLogistic <- glm(response ~ age + sex + array + PC1 + PC2 + PC3 + PC4 + predictor,data=trainDf,family="binomial")
covarLogistic <- glm(response ~ age + sex + array + PC1 + PC2 + PC3 + PC4,data=trainDf,family="binomial")
simpleLogistic <- glm(response ~ predictor,data=trainDf,family="binomial")
predConfint=confint(simpleLogistic,level=0.9)
r2 <- NagelkerkeR2(completeLogistic)$R2
expVar <- round(NagelkerkeR2(completeLogistic)$R2-NagelkerkeR2(covarLogistic)$R2,4)
if(expVar==0){expVar="< 0.001"}
pval <- -log10(coef(summary(completeLogistic))[9,4])
testProb=predict(completeLogistic,testDf,type="response")
testDf$regress=testProb
testDf=testDf[order(testDf$regress,decreasing = T),]
forReturn=cbind(testDf$response,testDf$regress)
thePlot=ggplot(testDf,aes(response,regress))+geom_violin(aes(fill=response))+
  geom_boxplot(width=0.1,outlier.alpha = 0)+
  labs(y="disease prediction",title=paste("Comparing Full Model Predictions for",currentTrait),
       subtitle=paste("Explained Variance:",as.character(expVar)),
       caption=paste("R-Squared:",round(r2,4),
                     "~ Best sampling method:",bestMethod,
                     "~ PRS Beta Confint:",round(predConfint[2,1],3),"-",round(predConfint[2,2],3)))+
  scale_x_discrete(labels=c("no disease","yes disease"))+
  scale_fill_discrete(guide=FALSE)+
  annotate("text",x=1:2,y=max(testDf$regress)*1.05,label=paste(table(testDf$response),"subjects"))
plot(thePlot)




#Odds Ratios and Relative Risks
allCutOffs=c(0.5,0.2,0.1,0.05,0.01,0.005)
resHolder=data.frame(matrix(0,nrow=length(allCutOffs),ncol=4))
colnames(resHolder)=c("fullOR","fullRR","bottomOR","bottomRR")
bottomHalf=testDf[(round(nrow(testDf)*0.5)+1):nrow(testDf),9]

i=1
for(cutOff in allCutOffs){
  exposeGroup=testDf[1:round(nrow(testDf)*cutOff),9]
  safeGroup=testDf[(round(nrow(testDf)*cutOff)+1):nrow(testDf),9]
  
  resHolder[i,1]=(sum(exposeGroup==1)/sum(exposeGroup==0))/(sum(safeGroup==1)/sum(safeGroup==0))
  resHolder[i,2]=(sum(exposeGroup==1)/length(exposeGroup))/(sum(safeGroup==1)/length(safeGroup))
  resHolder[i,3]=(sum(exposeGroup==1)/sum(exposeGroup==0))/(sum(bottomHalf==1)/sum(bottomHalf==0))
  resHolder[i,4]=(sum(exposeGroup==1)/length(exposeGroup))/(sum(bottomHalf==1)/length(bottomHalf))
  i=i+1
}

resHolder$cutOff=allCutOffs
resPlot=melt(resHolder,id.vars="cutOff")
resPlot$cutOff=factor(resPlot$cutOff,levels=unique(resPlot$cutOff[order(resPlot$cutOff,decreasing = T)]))
thePlot=ggplot(resPlot,aes(cutOff,value))+geom_point(aes(color=variable))+
  labs(title=paste("Risk/Odds Ratios for",currentTrait),
       x="Top Score Fraction Cutoff")
plot(thePlot)




#unadjusted disease prevalance
diseasePrev=rep(0,100)
diseaseOrHi=rep(0,100)
diseaseOrLo=rep(0,100)
testDf$group=rep(0,nrow(testDf))
groupRanges=c(round(seq(1,nrow(testDf),by = nrow(testDf)/100)),nrow(testDf))
for(i in 1:100){
  testDf$group[groupRanges[i]:groupRanges[i+1]]=i
}
testDf$group=rev(testDf$group)
for(i in 1:100){
  withDis=sum(testDf[testDf$group==i,9]==1)
  totalGroup=sum(testDf$group==i)
  diseasePrev[i]=withDis/totalGroup
  
  diseaseOrHi[i]=(sum(testDf[testDf$group>=i,9]==1)*sum(testDf[testDf$group<i,9]==0))/
    (sum(testDf[testDf$group>=i,9]==0)*sum(testDf[testDf$group<i,9]==1))
  diseaseOrLo[i]=(sum(testDf[testDf$group>i,9]==0)*sum(testDf[testDf$group<=i,9]==1))/
    (sum(testDf[testDf$group>i,9]==1)*sum(testDf[testDf$group<=i,9]==0))
}
diseaseOr=c((diseaseOrLo*(1/diseaseOrLo[50]))[1:50],(diseaseOrHi*(1/diseaseOrHi[50]))[51:100])
prevDf=data.frame(percentile=1:100,diseasePrev,diseaseOr)
diffTopBottom=mean(prevDf$diseasePrev[90:100])/mean(prevDf$diseasePrev[1:10])
corVal=round(cor(prevDf$diseasePrev,prevDf$percentile),5)
thePlot=ggplot(prevDf,aes(percentile,diseasePrev,color=diseaseOr))+geom_point()+
  labs(title=paste("Prevalence for",currentTrait),
       y="unadjusted disease prevalence",
       x="PRS percentile",
       subtitle=paste("multiple between top/bottom deciles:",round(diffTopBottom,4)),
       caption=paste("correlation:",as.character(corVal)))+
  guides(color=guide_legend(title="Std. OR"))
plot(thePlot)

thePlot=ggplot(testDf,aes(response,group))+
  geom_violin(aes(fill=response))+geom_boxplot(width=0.1,outlier.alpha = 0)+
  scale_x_discrete(labels=c("no disease","yes disease"))+
  scale_fill_discrete(guide=FALSE)+
  labs(title="Prevalence of Disease Split by Cases/Controls",y="Percentile")
plot(thePlot)



#ROC stuff
simpleLogistic <- glm(response ~ predictor,data=trainDf,family="binomial")
testDf$prob=predict(simpleLogistic,testDf,type="response")
g=roc(response ~ prob,data=testDf)
aucSimple=as.numeric(g$auc)
rocDf=data.frame(tpr=g$sensitivities,fpr=1-g$specificities,group=rep("PRS only",length(g$specificities)))

testDf$prob=predict(covarLogistic,testDf,type="response")
g=roc(response ~ prob,data=testDf)
aucCovar=as.numeric(g$auc)
rocDf=rbind(rocDf,data.frame(tpr=g$sensitivities,fpr=1-g$specificities,group=rep("covars",length(g$specificities))))

g=roc(response ~ regress,data=testDf)
aucComp=as.numeric(g$auc)
rocDf=rbind(rocDf,data.frame(tpr=g$sensitivities,fpr=1-g$specificities,group=rep("PRS+covars",length(g$specificities))))


if(aucSimple>aucComp){ compModel=simpleLogistic } else { compModel=completeLogistic }
anovaModel=anova(compModel,covarLogistic,test="LRT")
anovaPval=as.matrix(anovaModel)[2,5]
rocDf=rocDf[seq(1,nrow(rocDf),length.out = 3000),]
rocDf$group=factor(rocDf$group)
plotRocDf=data.frame(rocDf)
plotRocDf$fpr[plotRocDf$group=="covars"]=plotRocDf$fpr[plotRocDf$group=="covars"]+0.02
thePlot=ggplot(plotRocDf,aes(fpr,tpr,color=group))+geom_point()+
  geom_abline(intercept=c(0,0),slope=1)+
  labs(x="false positive rate",
       y="true positive rate",
       title=paste("ROC Curve for",currentTrait),
       subtitle=paste("LRT Pval between best PRS and Covar:",round(anovaPval,4)),
       caption=paste0("AUCs: simple=",round(aucSimple,4)," ~ covar=",round(aucCovar,4)," ~ complete=",round(aucComp,4)))
plot(thePlot)

dev.off()

return(list(sentOutTestDf,rocDf,aucComp,aucCovar,aucSimple,prevDf,resHolder,pval,forReturn))



}

	





