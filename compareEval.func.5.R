library(fmsb)
library(PRROC)
library(ggplot2)
library(forcats)
library(pROC)
library(reshape2)
library(cowplot)



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
totalHiMat=c()
totalLoMat=c()
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
      for(i in 1:ncol(predDf)){
        print("swap")
        if(mean(predDf[dfMain$response==0,i])>mean(predDf[dfMain$response==1,i])){
          predDf[,i]=predDf[,i]*-1
        }
      }
      df=cbind(dfMain,predDf)
      
      #collect the pval and OR by logistic regression
      resMat=data.frame(matrix(0,nrow=ncol(predDf),ncol=5))
      colnames(resMat)=c("index","pval","or","or2","or3")
      resMat[,1]=factor(workKey[,paramCol])
      sdLowMat=data.frame(resMat)
      sdHiMat=data.frame(resMat)
      
      logisticList=list()
      for(i in 1:ncol(predDf)){
        #normal logistic regression of score against response
        df=cbind(dfMain,predDf)
        completeLogistic <- glm(response ~ age + sex + array + PC1 + PC2 + PC3 + PC4 + predDf[,i],data=dfMain,family="binomial")
        compLogitCI=confint(completeLogistic,level=0.95)
        resMat[i,2]=-log10(coef(summary(completeLogistic))[9,4]) #pval
        resMat[i,3]=exp(coef(summary(completeLogistic))[9,1]) #odds ratio
        sdLowMat[i,3]=exp(compLogitCI[9,1])
        sdHiMat[i,3]=exp(compLogitCI[9,2])
        logisticList[[i]]=completeLogistic
        
        #logistic regression of top 5% against bottom 95%
        df=df[order(df[,8+i],decreasing = T),]
        #df=df[order(predDf[,i],decreasing = T),]
        exposeGroup=df[1:round(nrow(df)*0.05),8]
        safeGroup=df[(round(nrow(df)*0.05)+1):nrow(df),8]
        se=sqrt((1/sum(exposeGroup==1))+(1/sum(exposeGroup==0))+(1/sum(safeGroup==1))+(1/sum(safeGroup==0)))
        resMat[i,4]=(sum(exposeGroup==1)*sum(safeGroup==0))/(sum(exposeGroup==0)*sum(safeGroup==1))
        sdLowMat[i,4]=exp(log(resMat[i,4])+(1.96*se))
        sdHiMat[i,4]=exp(log(resMat[i,4])-(1.96*se))
        print(sum(safeGroup==1))
        print(sum(exposeGroup==0))
        #logistic regression of top 5% against bottom 50%
        #df=df[order(predDf[,i],decreasing = T),]
        exposeGroup=df[1:round(nrow(df)*0.05),8]
        safeGroup=df[(round(nrow(df)*0.5)+1):nrow(df),8]
        se=sqrt((1/sum(exposeGroup==1))+(1/sum(exposeGroup==0))+(1/sum(safeGroup==1))+(1/sum(safeGroup==0)))
        resMat[i,5]=(sum(exposeGroup==1)*sum(safeGroup==0))/(sum(exposeGroup==0)*sum(safeGroup==1))
        sdLowMat[i,5]=exp(log(resMat[i,5])+(1.96*se))
        sdHiMat[i,5]=exp(log(resMat[i,5])-(1.96*se))
        print(sum(safeGroup==1))
        print(sum(exposeGroup==0))
      }
      
      #PVAL PLOT!
      thePlot=ggplot(resMat,aes(index,pval))+geom_col(aes(fill=index))+
        labs(x=paramNames[paramCol-2],y="-log10(Pval)",title=paste("Analysis of",currentTrait),
             subtitle=paste("Using method",method,"- comparing",paramNames[paramCol-2]))+
        scale_fill_discrete(guide=FALSE)
      if(plotAll){plot(thePlot)}
      
      #OR PLOT!
      orMat=melt(resMat[,-2],id.vars = "index")
      orMat$sdHi=melt(sdHiMat[,-2],id.vars="index")[,3]
      orMat$sdLo=melt(sdLowMat[,-2],id.vars="index")[,3]
      thePlot=ggplot(orMat,aes(variable,value,color=index))+geom_point(position=position_dodge(width=0.3), size=4)+
        geom_errorbar(aes(ymin=sdLo,ymax=sdHi),width=0.2,position = position_dodge(0.3))+
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
          groupDf[j,i]=sum(subDf[subDf$group==j,1]==1)
        }
      }
      
      groupDf$index=1:100 
      sampSize=sum(subDf$group==1)
      plotDf=melt(groupDf,id.vars="index")
      plotDf$value=plotDf$value/sampSize
      plotDf$variable=factor(repDescribe)
      thePlot=ggplot(plotDf,aes(index,value,color=variable))+geom_smooth(fill="grey80")+
        labs(x="PRS percentile",y="disease prevalance",title=paste("Prevalance of",currentTrait),
             subtitle=paste("Using method",method,"- comparing",paramNames[paramCol-2]))
      if(plotAll){plot(thePlot)}
      
      topGroup=apply(groupDf[90:100,],2,sum)
      bottomGroup=apply(groupDf[1:10,],2,sum)
      quotPrev=topGroup/bottomGroup
      sePrev=sqrt((1/topGroup)+(1/bottomGroup)+(2/(10*sampSize)))
      resMat$prev=quotPrev[-length(quotPrev)]
      sdHiMat$prev=(quotPrev+exp(1.96*sePrev))[-length(quotPrev)]
      sdLowMat$prev=(quotPrev-exp(1.96*sePrev))[-length(quotPrev)]
      
      #calculate ROC stuff for models including and not including covariates
      covarLogistic <- glm(response ~ age + sex + array + PC1 + PC2 + PC3 + PC4,data=df,family="binomial")
      df$prob=predict(covarLogistic,type=c("response"))
      g=roc(response ~ prob,data=df)
      aucCovar=as.numeric(ci.auc(g))
      
      rocDf=c()
      aucs=matrix(0,nrow=length(logisticList),ncol=3)
      df=cbind(dfMain,predDf)
      i=1
      for(logit in logisticList){
        df$prob=predict(logit,type=c("response"))
        g=roc(response ~ prob,data=df)
        aucs[i,]=as.numeric(ci.auc(g$auc))
        #seRoc=as.data.frame(ci.se(g,specificities = seq(0, 1, length.out = 100),boot.n = 200))
        #seRoc$new=rownames(seRoc)
        #seRoc$group=rep(workKey[i,paramCol],100)
        #colnames(seRoc)=c("seLo","se","seHi","sp","group")
        #rocDf=rbind(rocDf,seRoc)
        rocDf=rbind(rocDf,data.frame(tpr=g$sensitivities,fpr=1-g$specificities,group=rep(workKey[i,paramCol],length(g$specificities))))
        i=i+1
      }
      rocDf=rocDf[seq(1,nrow(rocDf),length.out = 3000),]
      #rocDf$sp=as.numeric(rocDf$sp)
      totalRocDf=rbind(totalRocDf,rocDf)
      rocDf$group=factor(rocDf$group)
      
      thePlot=ggplot(rocDf,aes(y=tpr,x=fpr,color=group))+geom_point()+
        geom_abline(intercept=c(0,0),slope=1)+
        labs(x="false positive rate",y="true positive rate",
             title=paste("ROC Curves for",currentTrait),
             subtitle=paste("Using method",method,"- comparing",paramNames[paramCol-2]))
      if(plotAll){plot(thePlot)}
      
      resMat$AUC=aucs[,2]-aucCovar[2]
      sdLowMat$AUC=aucs[,1]-aucCovar[2]
      sdHiMat$AUC=aucs[,3]-aucCovar[2]
      thePlot=ggplot(resMat,aes(index,AUC))+geom_col(aes(fill=index))+
        geom_errorbar(aes(ymin=sdLowMat$AUC,ymax=sdHiMat$AUC,width=0.2))+
        labs(x=paramNames[paramCol-2],y="AUC improvement",title=paste("Analysis of",currentTrait),
             subtitle=paste("Using method",method,"- comparing",paramNames[paramCol-2]))+
        scale_fill_discrete(guide=FALSE)
      if(plotAll){plot(thePlot)}
  
      totalResMat=rbind(totalResMat,cbind(workKey,resMat))
      totalHiMat=rbind(totalHiMat,cbind(workKey,sdHiMat))
      totalLoMat=rbind(totalLoMat,cbind(workKey,sdLowMat))
    }
  }
}


#compare a specified statistic at the end
totalRocDf=totalRocDf[-1,]
toTake=which(totalResMat$AUC==max(totalResMat$AUC))[1]
totalRocDf=totalRocDf[(((toTake-1)*100)+1):(toTake*100),]
totalAuc=totalResMat$AUC[toTake]

totalResMat=totalResMat[!duplicated(totalResMat$fName),]
colnames(totalResMat)[7+critIndex]="toComp"
totalResMat=totalResMat[order(totalResMat$toComp),]
totalResMat$parameter=paste(totalResMat$method,as.character(totalResMat$p1),as.character(totalResMat$p2),sep='-')
totalResMat$parameter=factor(totalResMat$parameter,levels = totalResMat$parameter[order(totalResMat$toComp)])

totalLoMat=totalLoMat[!(duplicated(totalLoMat$fName)),]
totalHiMat=totalHiMat[!(duplicated(totalHiMat$fName)),]
colnames(totalLoMat)[7+critIndex]="toComp"
colnames(totalHiMat)[7+critIndex]="toComp"
totalResMat=totalResMat[order(totalResMat$fName),]
totalLoMat=totalLoMat[order(totalLoMat$fName),]
totalHiMat=totalHiMat[order(totalHiMat$fName),]

thePlot=ggplot(totalResMat,aes(parameter,toComp))+geom_col(aes(fill=methods))+
  geom_errorbar(aes(ymin=totalLoMat$toComp,ymax=totalHiMat$toComp,width=0.2))+
  labs(x="Parameters",y=finalCrit,title=paste("Total Analysis of",currentTrait))+
  guides(fill=FALSE)+theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))
plot(thePlot)

bestName=totalResMat$fName[totalResMat$toComp==max(totalResMat$toComp)][1]
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


