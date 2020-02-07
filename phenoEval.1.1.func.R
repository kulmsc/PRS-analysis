library(fmsb)
library(PRROC)
library(ggplot2)
library(forcats)
library(pROC)
library(reshape2)
library(stringr)
library(unbalanced)
library(cowplot)
theme_set(theme_cowplot())


phenoEval <- function(authorComp,bestFile,trainDf,currentTrait,phenoDefs,
                      trainPhenos,testPhenos,diseaseTestList,toRev,shortRead,covarDefs,enetBetas,
                      decoder){

options(warn=1)  
pdf(paste(authorComp,"pheno","pdf",sep='.'))

prepWork=function(scores,phenos,covar){
  decoder=phenoDefs[phenoDefs[,1]==authorComp,]
  decoder=unname(unlist(decoder))
  names(decoder)=colnames(phenoDefs)

  covarDecoder=covarDefs[covarDefs[,1]==authorComp,]
  covarDecoder=unname(unlist(covarDecoder))
  names(covarDecoder)=colnames(covarDefs) 
 
  if(!is.null(dim(scores))){ 
    predictor=scores[,which(authorComp==colnames(scores))]
  } else {
    predictor=scores
  } 
  print("34")
  print(head(predictor))

  predLength=nrow(scores)
  if(is.null(predLength)){predLength=length(scores)}
  response=rep(0,predLength)

  #get the case/control response for the pheno
  print("getting case control")
  for(phenIndex in 1:length(testPhenos)){
    if(names(testPhenos)[phenIndex] %in% names(decoder)[!is.na(decoder)]){
      fullKeyword = decoder[names(decoder) == names(testPhenos)[phenIndex]]
      includeKeyword = strsplit(fullKeyword,";",fixed=TRUE)[[1]][1]
      excludeKeyword = strsplit(fullKeyword,";",fixed=TRUE)[[1]][2]
      if(includeKeyword != "X"){
        for(keyword in strsplit(includeKeyword,"|",fixed=TRUE)[[1]]){
          print(keyword)
          if( any(grepl(keyword,colnames(phenos[[phenIndex]]))) ){
            temp <- rowSums(phenos[[phenIndex]][,grepl(keyword,colnames(phenos[[phenIndex]])),drop=F])
            print("check lengths")
	    print(length(temp))
            print(length(response))
            response = response + temp
          }
        }
      }
      if(!is.na(excludeKeyword)){
        for(keyword in strsplit(excludeKeyword,"|",fixed=TRUE)[[1]]){
          if( any(grepl(keyword,colnames(phenos[[phenIndex]]))) ){
            response = response - rowSums(phenos[[phenIndex]][,grepl(keyword,colnames(phenos[[phenIndex]])),drop=F])
          }
        }
      }
    }
  }
  #response[response>1]=1
  #response=factor(response)
  #THIS IS WHERE THINGS CHANGE!!! on DEC 26 2019 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print("70")
  response[response>0]=1
  response=factor(response)
  #### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  #get the covariate vector
  print("getting covariates")
  covarResp=c()
  itemsRemove=rep(0,length(response))
  if(length(covarDecoder) > 3){
    print("over three")
    covarResp=rep(0,predLength)
    for(phenIndex in 1:length(testPhenos)){
      if(names(testPhenos)[phenIndex] %in% names(covarDecoder)){
        fullKeyword = covarDecoder[names(covarDecoder) == names(testPhenos)[phenIndex]]
        for(keyword in strsplit(fullKeyword,"|",fixed=TRUE)[[1]]){
          if(substr(keyword,1,1)=="-"){
            realKeyword <- substr(keyword,2,nchar(keyword))
            if( any(grepl(realKeyword,colnames(testPhenos[[phenIndex]]))) ){
              intermedResp =  rowSums(testPhenos[[phenIndex]][,grep(realKeyword,colnames(testPhenos[[phenIndex]]),value=T),drop=F])
              itemsRemove = itemsRemove + intermedResp
              itemsRemove[itemsRemove>1] <- 1
            } else {
              print("keyword not in")
              print(keyword)
            }
          } else {
            if( any(grepl(keyword,colnames(phenos[[phenIndex]]))) ){
              intermedResp =  rowSums(testPhenos[[phenIndex]][,grep(keyword,colnames(testPhenos[[phenIndex]]),value=T),drop=F])
              intermedResp[intermedResp>1] = 1
              #covarResp = covarResp + intermedResp #not sure which covar system is best
              covarResp <- cbind(covarResp,factor(intermedResp))
            }
          }
        }
      }
    }
    if(!is.null(ncol(covarResp))){
      covarResp <- covarResp[,apply(covarResp,2,function(x) length(unique(x))>1),drop=F]
      colnames(covarResp)=paste("covarResp",1:ncol(covarResp),sep="")
      df=data.frame(covar,covarResp,predictor,response)
    } else {
      df=data.frame(covar,predictor,response)
    }

    if(any(itemsRemove==1)){
      df <- df[itemsRemove==0,]
    }
  } else {
    df=data.frame(covar,predictor,response)
  }
  print("return df at 120")
  return(df)
}

print("AT25")
if(is.null(enetBetas)){
  print("not going ENET")
  sentInScores=diseaseTestList[[which(names(diseaseTestList)==gsub("train","test",bestFile))]]
  if(is.na(sentInScores[nrow(sentInScores),1])){
    sentInScores <- rbind(apply(sentInScores[complete.cases(sentInScores),],2,mean), sentInScores[1:(nrow(sentInScores)-1),])
  }
} else {
  print("going ENET")
  print(enetBetas)
  allSentInScores=matrix(0,nrow=nrow(diseaseTestList[[1]]),ncol=length(enetBetas))
  for(i in 1:length(enetBetas)){
    print("139")
    print(names(enetBetas)[i])
    print(head(allSentInScores))
    pulledFile <- diseaseTestList[[which(names(diseaseTestList)==gsub("train","test",names(enetBetas)[i]))]]
    print(head(pulledFile))

    if(is.na(pulledFile[nrow(pulledFile),1])){
      print("145")
      keepNames <- colnames(pulledFile)
      pulledFile <- rbind(apply(pulledFile[complete.cases(pulledFile),],2,mean), pulledFile[1:(nrow(pulledFile)-1),])
      colnames(pulledFile) <- keepNames
    }
    if(sum(colnames(pulledFile)==authorComp) == 0){
      print("repping 0")
      allSentInScores[,i] <- rep(0, nrow(allSentInScores))
    } else {
      print("no repping")
      allSentInScores[,i] <- pulledFile[,which(colnames(pulledFile)==authorComp),drop=T]
    }
    print("before normo")
    print(length(unique(allSentInScores[,i])))
    if(length(unique(allSentInScores[,i])) > 1){
      print("doing normalization")
      allSentInScores[,i] <- (allSentInScores[,i]-mean(allSentInScores[,i]))/sd(allSentInScores[,i])
    }

  }
  print("new 150")
  print(head(allSentInScores))
  sentInScores <- rowSums(t(t(allSentInScores)*enetBetas))
  sentInScores <- (sentInScores-mean(sentInScores))/sd(sentInScores)
}

print("done at 150")
sentOutTestDf=prepWork(sentInScores, testPhenos, covar2)
colnames(sentOutTestDf)[colnames(sentOutTestDf)==authorComp]<-"predictor"
print("154")
print(head(sentOutTestDf))
print(table(sentOutTestDf$predictor))
sentOutTestDf$predictor <-  (sentOutTestDf$predictor-mean(sentOutTestDf$predictor))/sd(sentOutTestDf$predictor)


if(mean(sentOutTestDf$predictor[sentOutTestDf$response==1]) < (mean(sentOutTestDf$predictor[sentOutTestDf$response==0]))){
  print("Swap")
  sentOutTestDf$predictor=sentOutTestDf$predictor*-1
}

if(decoder[3]!="A"){
  sentOutTestDf <- sentOutTestDf[sentOutTestDf$sex==decoder[3],]
  sentOutTestDf <- sentOutTestDf[,-3]
}

testDf=sentOutTestDf[order(sentOutTestDf$predictor,decreasing=T),]



#Density plots
print("Density")
thePlot=ggplot(testDf,aes(x=predictor,fill=response))+geom_density(alpha=0.3)+
  labs(y="Density",x="raw polygenic risk score",title=paste("Raw Score Densities for",currentTrait),
  caption=paste("no. case:",sum(testDf$response==1)," ~ no. control:",sum(testDf$response==0)))+
  scale_fill_discrete(name="Disease Status",breaks=c(0,1),labels=c("No Disease","Disease"))
plot(thePlot)

#Assess best sampling method
newTrainDf=data.frame(trainDf)
newTrainDf[,1]=as.numeric(as.factor(newTrainDf[,1]))
if(decoder[3]=="A"){
  newTrainDf[,3]=as.numeric(as.factor(newTrainDf[,3]))
}

print("aucVec")
aucVec=matrix(0,nrow=3,ncol=6)
j=1
for(ubMeth in c("ubOver","ubUnder","ubSMOTE","nothing")){
  for(i in 1:3){
    sampleAmount=round(nrow(trainDf)/2)
    sampleSet=sort(sample(1:nrow(trainDf),sampleAmount))
    fullSet=1:nrow(trainDf)
    otherSet=fullSet[!(fullSet %in% sampleSet)]
    train=newTrainDf[sampleSet,]
    test=newTrainDf[otherSet,]
  
    if(ubMeth=="nothing"){
      data=list(X=train[,-ncol(trainDf)],Y=train[,ncol(trainDf)])
    } else {
      data <- ubBalance(X=train[,-ncol(trainDf)], Y=train[,ncol(trainDf)], type=ubMeth) 
    }
    outputDf=data.frame(data$X, data$Y)
    colnames(outputDf)=colnames(trainDf)
    covarLogistic <- glm(response ~ .,data=outputDf[,-which(colnames(outputDf)=="predictor")],family="binomial")
    testProb=predict(covarLogistic,test,type="response")
    test$prob=testProb
    g=pROC::roc(response ~ prob,data=test)
    test=test[,-ncol(test)]
    aucVec[i,j]=g$auc
  }
  j=j+1
}
print("samplingPvals")
samplingPvals=apply(aucVec[,-4],2,function(x) wilcox.test(x,aucVec[,4])$p.value)
if(any(samplingPvals<0.1)){
  bestMethod=c("ubOver","ubUnder","ubSMOTE")[samplingPvals==min(samplingPvals)]
} else {
  bestMethod="nothing"
}

print("Violin")
#Violin Plot
completeLogistic <- glm(response ~ . ,data=trainDf,family="binomial")
covarLogistic <- glm(response ~ . ,data=trainDf[,-which(colnames(trainDf)=="predictor")],family="binomial")
simpleLogistic <- glm(response ~ predictor,data=trainDf,family="binomial")

predConfint=confint(simpleLogistic,level=0.9)
r2 <- NagelkerkeR2(completeLogistic)$R2
expVar <- round(NagelkerkeR2(completeLogistic)$R2-NagelkerkeR2(covarLogistic)$R2,4)
if(expVar==0){expVar="< 0.001"}
pval <- -log10(coef(summary(completeLogistic))[ncol(trainDf),4])
testProb=predict(completeLogistic,testDf,type="response")
testDf$regress=testProb

print("geom split")
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

print("forReturn")
forReturn=cbind(testDf$response,testDf$regress)
thePlot=ggplot(testDf,aes(response,regress))+geom_violin(aes(fill=response))+
  geom_boxplot(width=0.1,outlier.alpha = 0)+
  labs(y="disease prediction",title=paste("Comparing Full Model Predictions for",currentTrait),
       subtitle=paste("Explained Variance:",as.character(expVar)),
       caption=paste("R2:",round(r2,4),
                     "~ Best sampling method:",bestMethod,
                     "~ PRS Beta Confint:",round(predConfint[2,1],3),"-",round(predConfint[2,2],3)))+
  scale_x_discrete(labels=c("no disease","yes disease"))+
  scale_fill_discrete(guide=FALSE)+
  annotate("text",x=1:2,y=max(testDf$regress)*1.05,label=paste(table(testDf$response),"subjects"))

print("BIG STUFF")
print(head(testDf$regress))
print(head(testDf$regress[order(testDf$regress,decreasing=T)]))
copyTest <- data.frame(testDf)
copyTest <- copyTest[copyTest$regress < testDf$regress[order(testDf$regress,decreasing=T)][100] &
             copyTest$regress > testDf$regress[order(testDf$regress,decreasing=F)][100],]
thePlot <- ggplot(copyTest,aes(x="base",y=regress, fill=response)) +geom_split_violin()+
  geom_boxplot(width = 0.25, outlier.shape = NA, coef=0)+
  scale_fill_discrete(name="status",labels=c("Healthy","Disease"))+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  labs(y="logisctic regression prediction",title=paste("Comparing Full Model Predictions for",currentTrait),
       subtitle=paste("Explained Variance:",as.character(expVar)),
       caption=paste("PRS Beta Confint:",round(predConfint[2,1],3),"-",round(predConfint[2,2],3)))
plot(thePlot)



print("Odds")
#Odds Ratios and Relative Risks
testDf=testDf[order(testDf$regress,decreasing = T),]
allCutOffs=c(0.5,0.2,0.1,0.05,0.01,0.005)
resHolder=data.frame(matrix(0,nrow=length(allCutOffs),ncol=8))
colnames(resHolder)=c("fullOR","fullRR","bottomOR","bottomRR","fullLo","fullHi","bottomLo","bottomHi")
resHolderHi=data.frame(resHolder)
resHolderLo=data.frame(resHolder)
bottomHalf=testDf[(round(nrow(testDf)*0.5)+1):nrow(testDf),ncol(trainDf)]

i=1
for(cutOff in allCutOffs){
  exposeGroup=testDf[1:round(nrow(testDf)*cutOff),ncol(trainDf)]
  safeGroup=testDf[(round(nrow(testDf)*cutOff)+1):nrow(testDf),ncol(trainDf)]
  
  resHolder[i,1]=(sum(exposeGroup==1)/sum(exposeGroup==0))/(sum(safeGroup==1)/sum(safeGroup==0))
  se=sqrt((1/sum(exposeGroup==1))+(1/sum(exposeGroup==0))+(1/sum(safeGroup==1))+(1/sum(safeGroup==0)))
  resHolder[i,5]=exp(log(resHolder[i,1])-(1.96*se))
  resHolder[i,6]=exp(log(resHolder[i,1])+(1.96*se))
  resHolder[i,2]=(sum(exposeGroup==1)/length(exposeGroup))/(sum(safeGroup==1)/length(safeGroup))
  resHolder[i,3]=(sum(exposeGroup==1)/sum(exposeGroup==0))/(sum(bottomHalf==1)/sum(bottomHalf==0))
  se=sqrt((1/sum(exposeGroup==1))+(1/sum(exposeGroup==0))+(1/sum(bottomHalf==1))+(1/sum(bottomHalf==0)))
  resHolder[i,7]=exp(log(resHolder[i,3])-(1.96*se))
  resHolder[i,8]=exp(log(resHolder[i,3])+(1.96*se))
  resHolder[i,4]=(sum(exposeGroup==1)/length(exposeGroup))/(sum(bottomHalf==1)/length(bottomHalf))
  i=i+1
}

resHolder$cutOff=allCutOffs
resPlot=melt(resHolder,id.vars="cutOff")
resPlot$cutOff=factor(resPlot$cutOff,levels=unique(resPlot$cutOff[order(resPlot$cutOff,decreasing = T)]))
rp1=resPlot[1:24,]
rp2=resPlot[25:48,]
thePlot=ggplot(resPlot[1:24,],aes(cutOff,value,color=variable))+geom_point(position=position_dodge(width=0.3))+
  geom_errorbar(aes(ymin=c(resPlot[25:30,3],rep(NA,6),resPlot[37:42,3],rep(NA,6)),
                    ymax=c(resPlot[31:36,3],rep(NA,6),resPlot[43:48,3],rep(NA,6))),
                width=0.2,position=position_dodge(width=0.3))+
  labs(title="Risk/Odds Ratios of Being in Score Groups", x="Top Score Fraction Cutoff",y="")
plot(thePlot)


print("Prev")
print("UNADJ PREV")
#unadjusted disease prevalance
testDf=testDf[order(testDf$predictor,decreasing = T),]
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
  withDis=sum(testDf[testDf$group==i,ncol(trainDf)]==1)
  totalGroup=sum(testDf$group==i)
  diseasePrev[i]=withDis/totalGroup
  
  #diseaseOrHi[i]=(sum(testDf[testDf$group>=i,9]==1)*sum(testDf[testDf$group<i,9]==0))/
  #  (sum(testDf[testDf$group>=i,9]==0)*sum(testDf[testDf$group<i,9]==1))
  diseaseOrLo[i]=(sum(testDf[testDf$group>i,ncol(trainDf)]==0)*sum(testDf[testDf$group<=i,ncol(trainDf)]==1))/
    (sum(testDf[testDf$group>i,ncol(trainDf)]==1)*sum(testDf[testDf$group<=i,ncol(trainDf)]==0))
  diseaseOrHi[i]=(sum(testDf[testDf$group>=i,ncol(trainDf)]==0)*sum(testDf[testDf$group<i,ncol(trainDf)]==1))/
    (sum(testDf[testDf$group>=i,ncol(trainDf)]==1)*sum(testDf[testDf$group<i,ncol(trainDf)]==0))
}
diseaseOr=c(diseaseOrLo[1],((diseaseOrHi[2:99]+diseaseOrLo[2:99])/2),diseaseOrHi[100])
diseaseOr=diseaseOr/diseaseOr[50]
#diseaseOr=c((diseaseOrLo*(1/diseaseOrLo[50]))[1:50],(diseaseOrHi*(1/diseaseOrHi[50]))[51:100])
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


print("ROC STUFF")
#ROC stuff
testDf$prob=predict(simpleLogistic,testDf,type="response")
gSimple=pROC::roc(response ~ prob,data=testDf)
aucSimple=as.numeric(ci.auc(gSimple))
print(length(gSimple$sensitivities))
print(length(gSimple$specificities))
rocDf=data.frame(tpr=gSimple$sensitivities,fpr=1-gSimple$specificities,group=rep("PRS only",length(gSimple$specificities)))
print("359")

testDf$prob=predict(covarLogistic,testDf,type="response")
gCovar=pROC::roc(response ~ prob,data=testDf)
aucCovar=as.numeric(ci.auc(gCovar))
rocDf=rbind(rocDf,data.frame(tpr=gCovar$sensitivities,fpr=1-gCovar$specificities,group=rep("covars",length(gCovar$specificities))))

testDf$prob=predict(completeLogistic,testDf,type="response")
gComp=pROC::roc(response ~ prob,data=testDf)
aucComp=as.numeric(ci.auc(gComp))
rocDf=rbind(rocDf,data.frame(tpr=gComp$sensitivities,fpr=1-gComp$specificities,group=rep("PRS+covars",length(gComp$specificities))))


print("DOING ANOVA")
#if(aucSimple[2]>aucComp[2]){ compModel=simpleLogistic } else { compModel=completeLogistic }
anovaModel=anova(covarLogistic,completeLogistic,test="LRT")
anovaPval=as.matrix(anovaModel)[2,5]
if(anovaPval<0.0001){anovaPval="<0.0001"} else {anovaPval=round(anovaPval,4)}

aucSimple=round(aucSimple,3)
aucComp=round(aucComp,3)
aucCovar=round(aucCovar,3)
rocDf=rocDf[seq(1,nrow(rocDf),length.out = 3000),]
rocDf$group=factor(rocDf$group)
plotRocDf=data.frame(rocDf)
plotRocDf$fpr[plotRocDf$group=="covars"]=plotRocDf$fpr[plotRocDf$group=="covars"]
thePlot=ggplot(plotRocDf,aes(fpr,tpr,color=group))+geom_line()+
  geom_abline(intercept=c(0,0),slope=1)+
  labs(x="false positive rate",
       y="true positive rate",
       title=paste("ROC Curve for",currentTrait),
       subtitle=paste("LRT Pval between best and covar model:",anovaPval),
       caption=paste0("AUCs: simple=",aucSimple[2],"(",aucSimple[1],"-",aucSimple[3],")","\n",
                      " ~ covar=",aucCovar[2],"(",aucCovar[1],"-",aucCovar[3],")","\n",
                      " ~ complete=",aucComp[2],"(",aucComp[1],"-",aucComp[3],")"))
plot(thePlot)

dev.off()

return(list(sentOutTestDf,rocDf,aucComp,aucCovar,aucSimple,prevDf,resHolder,pval,forReturn,list(gSimple,gCovar,gComp)))


}

	





