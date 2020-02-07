library(fmsb)
library(PRROC)
library(ggplot2)
library(forcats)
library(pROC)
library(reshape2)



compareEval<-function(authorComp,finalCrit,plotAll,decoder,key,pheno,diseaseTrainList){
  
pdf(paste(authorComp,"compare","pdf",sep='.'))

possibleCrits=c("pval","totalOR","groupOR","splitOR","Prev","AUC")
stopifnot(finalCrit %in% possibleCrits)
critIndex=which(possibleCrits %in% finalCrit)

#check all of the current files to see if the authorComp is within
authorFiles=c()
nameDiseaseFiles=list.files("diseaseFiles/","train")
for(f in nameDiseaseFiles){
  fOpen = file(paste0("diseaseFiles/",f),'r')
  firstLine=readLines(fOpen,n=1)
  close(fOpen)
  firstLine=strsplit(firstLine,split='\t')[[1]]
  if(authorComp %in% firstLine){
    authorFiles=c(authorFiles,f)
  }
}

#make a matrix to contain all of the scores from the files detected above
allScores=matrix(nrow=129568,ncol=length(authorFiles))
colnames(allScores)=authorFiles
for(f in authorFiles){
  goodDiseaseFile=diseaseTrainList[[which(nameDiseaseFiles==f)]]
  allScores[,which(colnames(allScores)==f)]=goodDiseaseFile[,which(colnames(goodDiseaseFile)==authorComp)]
}


checkZero=apply(allScores,2,sum)
print(sum(checkZero==0))
if(sum(checkZero==0)>0){
  badFile=names(checkZero)[checkZero==0]
  allScores=allScores[,-which(checkZero==0)]
  authorFiles=authorFiles[-which(authorFiles %in% badFile)]
}
allScores=apply(allScores,2,scale)

#in the decoder pull out the line corresponding to the author
decoder=decoder[decoder[,1]==authorComp,]
decoder=unname(unlist(decoder))

#fix the pheno file
allScores=allScores[!is.na(pheno[,6]),]
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

#pull out info from the key file
key=key[key$fName %in% authorFiles,]


#for each trait/score, get the uid and code that will pull out the correct myPheno and response vector
code=decoder[3]
currentTrait=decoder[2]
codeMat=sapply(strsplit(code,',')[[1]],function(x) strsplit(x,"-")[[1]])


#count number of time in all reporting the disease code comes up
predLength=nrow(allScores)
response=rep(0,predLength)
takenFrom=rep(1,predLength)
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

dfMain=data.frame(phenoCovar,age=ages,response)


totalResMat=c()
totalRocDf=c("tpr"=0,"fpr"=0,"group"=0)
#iterate over each method (clump,ldpred,etc.) (NOTE: will have to change later to compare both methods at end)
for(method in unique(key$method)){
  print("method")
  print(method)
  newKey=key[key$method==method,]
  print("newKey")
  print(newKey)
  
  if(method=="clump"){
    paramNames=c("pval","rho")
  } else if (method=="ldpred"){
    paramNames=c("rho","nothing")
  } else if (method=="winnersCurseLike"){
    paramNames=c("nothing","nothing")
  } else if (method=="winnersCurseLasso"){
    paramNames=c("lambda","nothing")
  } else if (method=="winnersCurse2d"){
    paramNames=c("pval","rho")
  } else if (method=="tweedy"){
    paramNames=c("method","nothing")
  } else if (method=="sumer"){
    paramNames=c("method","nothing")
  } else if (method=="prscs"){
    paramNames=c("phi","nothing")
  } else if (method=="grabld"){
    paramNames=c("method","nothing")
  } else if (method=="lassosum"){
    paramNames=c("method","nothing")
  } else if (method=="annoPred"){
    paramNames=c("m1","m2")
  } else if (method=="LDPredFunct"){
    paramNames=c("method","nothing")
  }
  #for the two types of parameters allowed (NOTE: will have to change later to allow more)
  for(paramCol in 3:(ncol(key)-1)){
    print(paramCol)
    if(length(unique(newKey[,paramCol]))>1 | paramCol==3){
      #subset the key to only include the parameter that is currently changing
      if(paramCol==4){
        otherCol=3
      } else {
        otherCol=4
      }
      otherTable=table(newKey[,otherCol])
      mostCommon=names(which(max(otherTable)==otherTable)[1]) #changed this from numeric to character
      workKey=newKey[newKey[,otherCol]==mostCommon,] #changeed this from numeric to character
      print("workKey")
      print(workKey)
      predDf=allScores[,colnames(allScores) %in% workKey$fName]
      if(length(predDf)==nrow(dfMain)){
        predDf=as.matrix(predDf)
      }
      df=cbind(dfMain,predDf)
      
      #collect the pval and OR by logistic regression
      resMat=data.frame(matrix(0,nrow=ncol(predDf),ncol=5))
      resMat[,1]=factor(workKey[,paramCol])
      
      colnames(resMat)=c("index","pval","or","or2","or3")
      logisticList=list()
      for(i in 1:ncol(predDf)){
        #normal logistic regression of score against response
        completeLogistic <- glm(response ~ age + sex + array + PC1 + PC2 + PC3 + PC4 + predDf[,i],data=dfMain,family="binomial")
        resMat[i,2]=-log10(coef(summary(completeLogistic))[9,4]) #pval
        resMat[i,3]=exp(coef(summary(completeLogistic))[9,1]) #odds ratio
        logisticList[[i]]=completeLogistic
        
        #logistic regression of top 5% against bottom 95%
        dfNew=dfMain[order(predDf[,i],decreasing = T),]
        dfNew$group1=rep(0,nrow(dfNew))
        dfNew$group1[1:round(nrow(df)*0.05)]=1
        dfNew$group1=factor(dfNew$group1)
        completeLogistic <- glm(response ~ age + sex + array + PC1 + PC2 + PC3 + PC4 + group1,data=dfNew,family="binomial")
        resMat[i,4]=exp(coef(summary(completeLogistic))[9,1])
        
        #logistic regression of top 5% against bottom 50%
        dfNew$group1=as.numeric(dfNew$group1)-1
        dfNew$group1[round(nrow(dfNew)/2):nrow(dfNew)] = -1
        dfNew=dfNew[dfNew$group1 %in% c(-1,1),]
        dfNew$group1=factor(dfNew$group1)
        completeLogistic <- glm(response ~ age + sex + array + PC1 + PC2 + PC3 + PC4 + group1,data=dfNew,family="binomial")
        resMat[i,5]=exp(coef(summary(completeLogistic))[9,1])
      }
      
      #PVAL PLOT!
      thePlot=ggplot(resMat,aes(index,pval))+geom_col(aes(fill=index))+
        labs(x=paramNames[paramCol-2],y="-log10(Pval)",title=paste("Analysis of",currentTrait),
             subtitle=paste("Using method",method,"- comparing",paramNames[paramCol-2]))+
        scale_fill_discrete(guide=FALSE)
      if(plotAll){plot(thePlot)}
      
      #OR PLOT!
      orMat=resMat[,-2]
      orMat=melt(orMat,id.vars = "index")
      thePlot=ggplot(orMat,aes(variable,value))+geom_col(aes(fill=index),position="dodge")+
        labs(x="OR comparison",y="odds ratio per SD",title=paste("Analysis of",currentTrait),
          subtitle=paste("Using method",method,"- comparing",paramNames[paramCol-2]))+
        scale_x_discrete(labels=c("by total score","top 5%, bottom 95%","top 5%, bottom 50%"))+
        guides(fill=guide_legend(title=paramNames[paramCol-2]))
      if(plotAll){plot(thePlot)}
  
      
      #graph the unadjusted disease prevalance
      groupDf=data.frame(matrix(0,ncol=ncol(predDf),nrow=100))
      groupRanges=c(round(seq(1,nrow(df),by = nrow(df)/100)),nrow(df))
      repDescribe=c()
      for(i in 1:ncol(predDf)){
        subDf=df[order(df[,8+i]),c(8,8+i)]
        subDf$group=rep(0,nrow(subDf))
        repDescribe=c(repDescribe,rep(workKey[i,paramCol],100))
        for(j in 1:100){
          subDf$group[groupRanges[j]:groupRanges[j+1]]=j
          withDis=sum(subDf[subDf$group==j,1]==1)
          totalGroup=sum(subDf$group==j)
          groupDf[j,i]=withDis/totalGroup
        }
      }
      groupDf$index=1:100
      plotDf=melt(groupDf,id.vars="index")
      plotDf$variable=factor(repDescribe)
      thePlot=ggplot(plotDf,aes(index,value,color=variable))+geom_smooth(fill="grey80")+
        labs(x="PRS percentile",y="disease prevalance",title=paste("Prevalance of",currentTrait),
             subtitle=paste("Using method",method,"- comparing",paramNames[paramCol-2]))
      if(plotAll){plot(thePlot)}
      topPrev=apply(groupDf[90:100,],2,mean)/apply(groupDf[1:10,],2,mean)
      resMat$prev=topPrev[-length(topPrev)]
      
      #calculate ROC stuff for models including and not including covariates
      covarLogistic <- glm(response ~ age + sex + array + PC1 + PC2 + PC3 + PC4,data=df,family="binomial")
      df$prob=predict(covarLogistic,type=c("response"))
      g=roc(response ~ prob,data=df)
      aucCovar=as.numeric(g$auc)
      
      rocDf=c()
      aucs=c()
      i=1
      for(logit in logisticList){
        df$prob=predict(logit,type=c("response"))
        g=roc(response ~ prob,data=df)
        aucs=c(aucs,as.numeric(g$auc))
        rocDf=rbind(rocDf,data.frame(tpr=g$sensitivities,fpr=1-g$specificities,group=rep(workKey[i,paramCol],length(g$specificities))))
        i=i+1
      }
      rocDf=rocDf[seq(1,nrow(rocDf),length.out = 3000),]
      totalRocDf=rbind(totalRocDf,rocDf)
      rocDf$group=factor(rocDf$group)
      
      #ROC PLOT!
      thePlot=ggplot(rocDf,aes(fpr,tpr,color=group))+geom_point()+
        geom_abline(intercept=c(0,0),slope=1)+
        labs(x="false positive rate",y="true positive rate",
             title=paste("ROC Curves for",currentTrait),
             subtitle=paste("Using method",method,"- comparing",paramNames[paramCol-2]))
      
      if(plotAll){plot(thePlot)}
      
      resMat$AUC=aucs-aucCovar
      thePlot=ggplot(resMat,aes(index,AUC))+geom_col(aes(fill=index))+
        labs(x=paramNames[paramCol-2],y="AUC improvement",title=paste("Analysis of",currentTrait),
             subtitle=paste("Using method",method,"- comparing",paramNames[paramCol-2]))+
        scale_fill_discrete(guide=FALSE)
      if(plotAll){plot(thePlot)}
  
      totalResMat=rbind(totalResMat,cbind(workKey,resMat))
    }
  }
}


#compare a specified statistic at the end
totalRocDf=totalRocDf[-1,]
toTake=which(totalResMat$AUC==max(totalResMat$AUC))
totalRocDf=totalRocDf[(((toTake-1)*1000)+1):(toTake*1000),]
totalAuc=totalResMat$AUC[toTake]

totalResMat=totalResMat[!duplicated(totalResMat$fName),]
colnames(totalResMat)[7+critIndex]="toComp"
totalResMat=totalResMat[order(totalResMat$toComp),]
totalResMat$parameter=paste(totalResMat$method,as.character(totalResMat$p1),as.character(totalResMat$p2),sep='-')
print(totalResMat$parameter)
totalResMat$parameter=factor(totalResMat$parameter,levels = totalResMat$parameter[order(totalResMat$toComp)])

thePlot=ggplot(totalResMat,aes(parameter,toComp))+geom_col(aes(fill=parameter))+
  labs(x="Parameters",y=finalCrit,title=paste("Total Analysis of",currentTrait))+
  guides(fill=FALSE)+theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))
plot(thePlot)

bestName=totalResMat$fName[nrow(totalResMat)]
bestScore=allScores[colnames(allScores)==bestName]
goBackDf=cbind(df[,1:8],bestScore)
colnames(goBackDf)[9]="predictor"
goBackDf=goBackDf[,c(1:7,9,8)]
if(mean(goBackDf$predictor[goBackDf$response==1]) < mean(goBackDf$predictor[goBackDf$response==0])){
  goBackDf$predictor=goBackDf$predictor*-1
  revScore=TRUE
} else {
  revScore=FALSE
}

totalRocDf=totalRocDf[,1:2]
totalRocDf$group="train"
dev.off()

return(list(currentTrait,bestName,goBackDf,decoder,totalRocDf,totalAuc[1],revScore,totalResMat))

}


