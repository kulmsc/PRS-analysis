library(ggplot2)
library(gplots)
library(gridExtra)
library(grid)
library(reshape2)
library(ggrepel)
library(Rtsne)
library(cowplot)
library(viridis)
library(ggridges)
library(pROC)
library(stringr)
library(vroom)
library(scales)
library(umap)
library(fmsb)

rdsHome <- "done_plots_1/"


source(paste0(rdsHome,"note"))
excludeStack <- !includeStack
excludeEnet <- !includeEnet
stopAfterPheno <- FALSE
theme_set(theme_cowplot())

# PREPARE PREPARE PREPARE
#######################################################################################################################
truCompHolder=c()
truPhenoHolder=c()
truBoostHolder=c()
truCorrHolder=c()

for(rFile in list.files(path=rdsHome,pattern="done.RDS")){
  print(rFile)
  readIn <- readRDS(paste0(rdsHome,rFile))
  truCompHolder=c(truCompHolder,readIn[[1]])
  truPhenoHolder=c(truPhenoHolder,readIn[[2]])
  #truBoostHolder=c(truBoostHolder,readIn[[3]])
  truCorrHolder=c(truCorrHolder,readIn[[3]])
  print(length(readIn[[2]]))
}

compLength <- length(readIn[[1]])
phenoLength <- length(readIn[[2]])
#boostLength <- length(readIn[[3]])
corrLength <- length(readIn[[3]])
rm(readIn)

savePlots <- list()

pdf(paste0(rdsHome,"/metaPlots.comp.pdf"))



#COMP COMP COMP COMP COMP COMP COMP COMP COMP COMP COMP COMP COMP COMP COMP COMP COMP COMP COMP COMP COMP COMP COMP COMP
###########################################################################################################################
#########################################################################################################################

#Table describing methods
#method=read.table("methodsTable",stringsAsFactors = F,header=T,sep=';')
#tt1=ttheme_minimal()
#grid.table(method,theme = tt1,rows=NULL)

#############################################################################################################################
#distribution of betas


# allBeta <- list()
# j <- 1
# for(i in seq(1,length(smallCompHolder, compLength), 2)){
#   currAuthor <- authors[j]
#   fileName <- smallCompHolder[[i]]
#   if(grepl("Enet",fileName)){
#     enetBeta <- smallCompHolder[[i+7]]
#     enetBeta <- enetBeta[enetBeta!=0]
#     toCompileBeta <- list()
#     toCompileRsids <- c()
#     
#     k <- 1
#     for(enetFile in names(enetBeta)){
#       enetFile <- str_split(enetFile, fixed("."))[[1]]
#       setFile <- paste(currAuthor,enetFile[4],"final.set",enetFile[2],"gz",sep=".")
#       #readFile <- vroom(paste0("/Volumes/Macintosh_HD_2/finalSets/",setFile))
#       toCompileBeta[[k]] <- vroom(paste0("/Volumes/Macintosh_HD_2/finalSets/", setFile),
#                                   col_names = F, col_select = c("X2","X8"))
#       toCompileBeta[[k]] <- toCompileBeta[[k]][toCompileBeta[[k]][,2]!=0,]
#       toCompileRsids <- unique(c(toCompileRsids, toCompileBeta[[k]][,1,drop=T]))
#       k <- k + 1
#     }
#     
#     realBeta <- rep(0,length(toCompileRsids))
#     k <- 1
#     for(compileBeta in toCompileBeta){
#       compileBeta <- compileBeta[!duplicated(compileBeta[,1]),]
#       print(length(realBeta[toCompileRsids %in% compileBeta[,1,drop=T]]))
#       print(length(compileBeta[,2,drop=T]*enetBeta[k]))
#       
#       realBeta[toCompileRsids %in% compileBeta[,1,drop=T]] <- 
#         realBeta[toCompileRsids %in% compileBeta[,1,drop=T]] + 
#         compileBeta[,2,drop=T]*enetBeta[k]
#       k <- k + 1
#     }
#     
#   }
# }

#############################################################################################################################
# Histogram

authors=c()
for(i in seq(1,length(truCompHolder), compLength)){
  authors=c(authors,truCompHolder[[i+4]])
}

traitNames=c()
for(i in seq(1,length(truCompHolder), compLength)){
  traitNames <- c(traitNames,truCompHolder[[i+1]])
}

for(trait in traitNames){
  if(sum(traitNames==trait)>1){
    badTrait=traitNames[traitNames==trait]
    traitNames[traitNames==trait]=paste(badTrait,1:length(badTrait),sep="_")
  }
}

for(excludeEnet in c(FALSE, TRUE)){
  if(excludeEnet){theCaption <- "excluding Enet"} else {theCaption <- "with Enet"}
  bestMethods=data.frame(matrix(0,nrow=length(unique(truCompHolder[[3]]$methods)),ncol=6))
  rownames(bestMethods)=unique(truCompHolder[[3]]$methods)
  colnames(bestMethods)=colnames(truCompHolder[[3]])[8:13]
  colnames(bestMethods)[6]="AUC"
  
  for(i in seq(1,length(truCompHolder), compLength)){
    scoreTrack=truCompHolder[[i+2]]
    if(excludeEnet){
      scoreTrack <- scoreTrack[-which(scoreTrack$methods=="enet"),]
    }
    for(j in 1:6){
       singBest=scoreTrack$methods[scoreTrack[,j+6]==max(scoreTrack[,j+6])][1]
       bestMethods[which(rownames(bestMethods)==singBest),j] = bestMethods[which(rownames(bestMethods)==singBest),j] + 1
    }
  }
  
  #DOT PLOT
  bestMethods=bestMethods[,c(1,2,5,6)]
  colnames(bestMethods)[1:4]=c("PVal","OR","R2","AUC")
  bestMethods$method=rownames(bestMethods)
  if(excludeEnet){bestMethods <- bestMethods[-which(bestMethods$method=="enet"),]}
  plotDf=melt(bestMethods,id.vars = "method")
  thePlot=ggplot(plotDf,aes(variable,value/sum(plotDf$value[plotDf$variable=="AUC"]),color=method))+geom_point(position=position_dodge(width=0.2))+
    labs(x="Assessment Statistic",y="Fraction Best",caption=theCaption)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot(thePlot)
  
  
  
  #COLUMN PLOT
  bestMethods$method=factor(bestMethods$method,levels = bestMethods$method[order(bestMethods$AUC)])
  thePlotAuc=ggplot(bestMethods,aes(method, AUC))+geom_col()+
    labs(y="",x="",subtitle="AUC")+theme(axis.text.x = element_blank())
  
  bestMethods$method=factor(bestMethods$method,levels = bestMethods$method[order(bestMethods$AUC)])
  thePlotPval=ggplot(bestMethods,aes(method, PVal))+geom_col()+
    labs(y="",x="",subtitle="PVal")+theme(axis.text.x = element_blank())
  
  bestMethods$method=factor(bestMethods$method,levels = bestMethods$method[order(bestMethods$AUC)])
  thePlotOR=ggplot(bestMethods,aes(method, OR))+geom_col()+
    labs(y="",x="",subtitle="OR")+theme(axis.text.x = element_blank())
  
  bestMethods$method=factor(bestMethods$method,levels = bestMethods$method[order(bestMethods$AUC)])
  thePlotR2=ggplot(bestMethods,aes(method, R2))+geom_col()+
    labs(y="",x="",subtitle="R2",caption=theCaption)+theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  grid.arrange(thePlotAuc, thePlotPval, thePlotOR, thePlotR2, ncol=1, heights=c(5,5,5,8))
}

#############################################################################################################################
#dot plot of  the average score for each method
excludeEnet <- !includeEnet
excludeStack <- !includeStack
fileDefs <- read.table("fileDefs",stringsAsFactors = F)
if(excludeStack){
  fileDefs <- fileDefs[fileDefs[,2]!="stackCT",]
}
if(excludeEnet){
  allMethods=unique(fileDefs[,2])
} else {
  allMethods=c("enet",unique(fileDefs[,2]))
}
splitMethods=data.frame(matrix(0,nrow=((length(truCompHolder)/3)*length(allMethods)),ncol=(ncol(truCompHolder[[3]])+1)))
colnames(splitMethods)=c(colnames(truCompHolder[[3]]),"trait")

#for each trait get the best comp results for each method
j=1
k=1
for(i in seq(1,length(truCompHolder), compLength)){
  scoreTrack=truCompHolder[[i+2]]
  for(meth in allMethods){
    if(meth %in% scoreTrack$method){
      subScoreTrack=scoreTrack[scoreTrack$method==meth,]
      toAdd=subScoreTrack[subScoreTrack$toComp==max(subScoreTrack$toComp),]
      toAdd$parameter=as.character(toAdd$parameter)
      splitMethods[j,1:ncol(toAdd)]=toAdd[1,]
      splitMethods[j,ncol(splitMethods)]=traitNames[k]
      j=j+1
    }
  }
  k=k+1
}

splitMethods=splitMethods[1:(j-1),]
splitMethods[grep("report",splitMethods[,1]),2] <- "report" #REMOVE THIS LATER!!!


#Repeat process but for the absolute best
bestTraitInds <- rep(0, length(unique(splitMethods$trait)))
j <- 1
for(trait in unique(splitMethods$trait)){
  subDf <- splitMethods[splitMethods$trait==trait,]
  bestTraitInds[j] <- which(splitMethods[,1] == subDf[subDf$toComp == max(subDf$toComp),1] & splitMethods$trait == trait)
  j <- j + 1
}

bestMethodInds <- rep(0, length(unique(splitMethods$methods)))
j <- 1
for(meth in unique(splitMethods$methods)){
  subDf <- splitMethods[splitMethods$methods==meth,]
  bestMethodInds[j] <- which(splitMethods[,1] == subDf[subDf$toComp == max(subDf$toComp),1] &
                               splitMethods$methods == meth &
                               splitMethods[,14] == subDf[subDf$toComp == max(subDf$toComp),14])
  j <- j + 1
}

#get the mean of each method in a similar formatting style
summMethods <- data.frame(matrix(NA,nrow=length(unique(splitMethods$methods)),ncol=ncol(splitMethods)))
summTraits <- data.frame(matrix(NA,nrow=length(unique(splitMethods$trait)),ncol=ncol(splitMethods)))
colnames(summMethods)=colnames(splitMethods)
colnames(summTraits)=colnames(splitMethods)
summMethodsBest <- data.frame(summMethods)
summTraitsBest <- data.frame(summTraits)
j=1
for(meth in unique(splitMethods$methods)){
  summMethods[j,c(1:6,14)]=c("compiled",meth,0,0,"disease",0,"mean")
  summMethods[j,7:12]=apply(splitMethods[splitMethods$method==meth,7:12],2,mean)
  summMethodsBest[j,c(1:6,14)]=c("compiled",meth,0,0,"disease",0,"mean")
  summMethodsBest[j,7:12]=apply(splitMethods[splitMethods$method==meth & 1:nrow(splitMethods) %in% bestMethodInds,7:12],2,mean)
  j=j+1
}
j=1
for(trait in as.character(unique(splitMethods$trait))){
  summTraits[j,c(1:6,14)]=c("compiled","mean",0,0,"disease",0,trait)
  summTraits[j,7:12]=apply(splitMethods[splitMethods$trait==trait,7:12],2,mean)
  summTraitsBest[j,c(1:6,14)]=c("compiled","mean",0,0,"disease",0,trait)
  summTraitsBest[j,7:12]=apply(splitMethods[splitMethods$trait==trait & 1:nrow(splitMethods) %in% bestTraitInds,7:12],2,mean)
  j=j+1
}



#Iterating over the 4 different statistics to create 4 different sets of 4 plots
yLabels <- c("-log10(Pval)","Odds Ratio","Nagelkerke R2","AUC Improvement")
for(j in c(7,10:12)){
  theY <- yLabels[which(j == c(7,10:12))]
  
  #Method - All Param
  forPlot=rbind(splitMethods,summMethods)
  colnames(forPlot)[12] <- "AUC"
  colnames(forPlot)[j] <- "toComp"
  forPlot$methods <- factor(forPlot$methods,levels=forPlot$methods[forPlot$trait=="mean"][order(forPlot[forPlot$trait=="mean",j])])
  thePlot=ggplot(forPlot,aes(methods,toComp))+geom_point()+geom_point(data=forPlot[forPlot$trait=="mean",],aes(color="red"))+
    theme(legend.position="none")+theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
    labs(title="Comparison of each Method - All",y=theY,x="Score Method")
  plot(thePlot)
  savePlots[[1]] <- forPlot
  
  #Method - Best Param
  forPlot=splitMethods[bestMethodInds,]
  colnames(forPlot)[12] <- "AUC"
  colnames(forPlot)[j] <- "toComp"
  forPlot$methods <- factor(forPlot$methods,levels=forPlot$methods[order(forPlot$toComp)])
  thePlot=ggplot(forPlot,aes(methods,toComp))+geom_point()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
    labs(title="Comparison of each Method - Best",y=theY,x="Score Method")
  plot(thePlot)
  savePlots[[2]] <- forPlot
  
  #Trait - All Param
  forPlot=rbind(splitMethods,summTraits)
  colnames(forPlot)[12] <- "AUC"
  colnames(forPlot)[j] <- "toComp"
  forPlot$trait <- factor(forPlot$trait,levels=forPlot$trait[forPlot$methods=="mean"][order(forPlot[forPlot$methods=="mean",j])])
  thePlot=ggplot(forPlot,aes(trait,toComp))+geom_point()+geom_point(data=forPlot[forPlot$methods=="mean",],aes(color="red"))+
    theme(legend.position="none")+theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
    labs(title="Comparison of each Trait - All",y=theY,x="Trait")
  plot(thePlot)
  savePlots[[3]] <- forPlot
  
  #Trait - Best Param
  forPlot=splitMethods[bestTraitInds,]
  colnames(forPlot)[12] <- "AUC"
  colnames(forPlot)[j] <- "toComp"
  forPlot$trait <- factor(forPlot$trait,levels=forPlot$trait[order(forPlot$toComp)])
  thePlot=ggplot(forPlot,aes(trait,toComp))+geom_point()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
    labs(title="Comparison of each Trait - Best",y=theY,x="Trait")
  plot(thePlot)
  savePlots[[4]] <- forPlot
}

forPlot=rbind(splitMethods,summTraits)
simpleMeth <- c("clump","tweedy","winnersCurseLasso","winnersCurse-2d","winnersCurseLike")
compMeth <- c("grabld","lassosum","prsCS","annoPred","sblup","sbayesr","ldpred")
forPlot$type <- NA
forPlot$type[forPlot$methods %in% simpleMeth] <- "Simple"
forPlot$type[forPlot$methods %in% compMeth] <- "Complex"
graphicalAbstract <- ggplot(forPlot[!is.na(forPlot$type),],aes(type,toComp))+geom_col()
plot(graphicalAbstract)

#################################################################################################
#REGRESS AGAINST META INFO ######################################################################
ukbbSampSize <- readRDS("caseControl.RDS")
metaData=read.table("metaData",stringsAsFactors = F)
allH2=read.table("allH2",stringsAsFactors = F)
allMeta <- cbind(ukbbSampSize[1,], ukbbSampSize[1,]/ukbbSampSize[2,], metaData[,-1], allH2[,-1])
assocHolder <- data.frame(matrix(0, nrow=length(unique(splitMethods$methods)), ncol=ncol(allMeta)))
colnames(assocHolder) <- c("ukbbCases", "ukbbRatio", "sampleSize","snps","cases","herit")
i <- 1
for(meth in unique(splitMethods$methods)){
  subDf <- splitMethods[splitMethods$methods == meth, ]
  subMeta <- allMeta[rownames(allMeta) %in% subDf$trait,]
  for(j in 1:ncol(subMeta)){
   regRes <- lm(subDf$toComp ~ subMeta[,j])
   assocHolder[i,j] <- summary(regRes)$coefficients[2,4]
  }
  i <- i + 1
}

assocHolder$methods <- unique(splitMethods$methods)
plotDf <- melt(assocHolder, id.vars = "methods")
thePlot <- ggplot(plotDf, aes(value, methods, color=variable)) + geom_point() +
  labs(x="PVal", y="Methods", color = "Meta Info")
plot(thePlot)
thePlot <- ggplot(plotDf, aes(-log10(value), methods, color=variable)) + geom_point() +
  labs(x="-log10(PVal)", y="Methods", color = "Meta Info")
plot(thePlot)

#get best files
bestFiles <- rep("0", length(unique(splitMethods$trait)))
i <- 1
for(tra in unique(splitMethods$trait)){
  subDf <- splitMethods[splitMethods$trait == tra & splitMethods$fName != "bestEnet" & splitMethods$fName != "allEnet",]
  bestFiles[i] <- subDf$fName[which.max(subDf$toComp)]
  i <- i + 1
}
write.table(cbind(traitNames, authors, bestFiles), "bestFiles", row.names=F, col.names=F, quote=F, sep='\t')

#SERIES OF DOUBLE BAR PLOTS ########################################################################
#now just looking at the methods

i <- 1
meanOverTraits <- splitMethods[rep(1,length(unique(splitMethods$methods))),]
for(meth in unique(splitMethods$methods)){
  meanOverTraits[i,c(1,2,13,14)] <- splitMethods[splitMethods$methods==meth,][1,c(1,2,13,14)]
  meanOverTraits[i,3:12] <- apply(splitMethods[splitMethods$methods==meth,3:12],2,mean)
  i <- i + 1
}
forPlot <- melt(meanOverTraits[,c(2,7,9,11,12)],id.vars = "methods")
allPlots <- list()
i <- 1
xaxs <- c("Pval","OR","R2","AUC Imp")
for(splitBy in unique(forPlot$variable)){
  allPlots[[i]] <- ggplot(forPlot[forPlot$variable==splitBy,],aes(methods,value))+
    geom_col()+guides(fill=F)+labs(y=xaxs[i])
  if(splitBy!=unique(forPlot$variable)[4]){
    allPlots[[i]] <- allPlots[[i]] + theme(axis.title.x=element_blank(), axis.text.x=element_blank())
  } else {
    allPlots[[i]] <- allPlots[[i]] + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  i <- i + 1
}
grid.arrange(allPlots[[1]],allPlots[[2]],allPlots[[3]],allPlots[[4]],ncol=1,heights=c(5,5,5,8))

meanOverTraits[,3:12] <- apply(meanOverTraits[,3:12], 2, function(x) (x-min(x))/(max(x)-min(x)))
forPlot <- melt(meanOverTraits[,c(2,7,9,11,12)],id.vars = "methods")
thePlot <- ggplot(forPlot,aes(methods,value,fill=variable))+geom_bar(stat="identity", position="dodge",width=0.8)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x="Methods",y="Scaled Value",fill="Statistic")+scale_fill_discrete(labels=c("Pval","OR","R2","AUC Imp"))
plot(thePlot)

thePlot <- ggplot(forPlot,aes(variable,value,fill=methods))+geom_bar(stat="identity", position="dodge",width=0.8)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x="Statistic",y="Scaled Value")+scale_x_discrete(labels=c("Pval","OR","R2","AUC Imp"))
plot(thePlot)


#########################################################################################################################
#Compare enet to the next best

compAucs <- data.frame(matrix(0,nrow=2,ncol=length(traitNames),dimnames = list(c("enet","notEnet"),traitNames)))

j <- 1
for(i in seq(1,length(truCompHolder), compLength)){
  bestDf <- truCompHolder[[i+2]]
  bestFiles <- bestDf$fName[bestDf$toComp %in% sort(bestDf$toComp,decreasing = T)[1:3]]
  if(grepl("Enet",bestFiles[2])){
    bestFiles <- bestFiles[-2]
  } else {
    bestFiles <- bestFiles[-3]
  }
  compAucs[,j] <- bestDf$toComp[bestDf$fName %in% bestFiles]
  j <- j + 1
}

compAucs <- compAucs[,order(compAucs[1,])]
plotDf <- melt(t(compAucs))
thePlot <- ggplot(plotDf, aes(value,Var1,color=Var2))+geom_point()+
  labs(x="AUC Imp",y="Trait",color="Method")
plot(thePlot)
savePlots[[5]] <- plotDf

#############################################################################################################################
#enet weights
allEnetWeights <- data.frame(matrix(0,nrow = length(authors), ncol = length(allMethods)))
rownames(allEnetWeights) <- traitNames
colnames(allEnetWeights) <- allMethods

k=1
for(i in seq(1,length(truCompHolder), compLength)){
  pulledBetas <- truCompHolder[[i+6]]
  if(is.null(pulledBetas)){
    methodName <- str_split(truCompHolder[[i]],fixed("."))[[1]][[2]]
    allEnetWeights[k,which(colnames(allEnetWeights)==methodName)] <- 1
  } else {
    #pulledBetas <- pulledBetas[pulledBetas > 0]
    pulledBetas <- abs(pulledBetas)/sum(abs(pulledBetas))
    pulledNames <- str_split(names(pulledBetas),fixed("."),simplify = T)[,2]
    for(j in 1:length(pulledBetas)){
      allEnetWeights[k,which(colnames(allEnetWeights)==pulledNames[j])] <- pulledBetas[j] +
        allEnetWeights[k,which(colnames(allEnetWeights)==pulledNames[j])]
    }
  }
  k <- k + 1
}

allEnetWeights <- allEnetWeights[,order(colSums(allEnetWeights))]
allEnetWeights$trait <- rownames(allEnetWeights)
plotDf <- melt(allEnetWeights)
thePlot <- ggplot(plotDf, aes(variable, value)) + geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
  labs(x="Method",y="Enet Weight")
plot(thePlot)

plotDf <- plotDf[-which(plotDf$value == 0),]
thePlot <- ggplot(plotDf, aes(trait, value, color=variable)) + geom_point()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.margin = unit(c(0.5,0.5,0.5,2),"cm"))+
  labs(x="Trait",y="Enet Weight")+scale_color_discrete(name="Method")
plot(thePlot)


#############################################################################################################################
#get the mean and sd for each method
meanAndSd=data.frame(matrix(0,nrow=nrow(summMethods),ncol=13))
colnames(meanAndSd)=c(paste0(colnames(splitMethods)[7:12],"Mean"),paste0(colnames(splitMethods)[7:12],"Sd"),"method")
j=1
for(meth in summMethods$methods){
  pullMethod <- splitMethods[splitMethods$method==meth,7:12]
  pullMethod <- pullMethod[apply(pullMethod, 1,function(y) all(!is.infinite(as.numeric(y)))), ]
  meanAndSd[j,1:6]=apply(pullMethod,2,mean)
  meanAndSd[j,7:12]=apply(pullMethod,2,sd)
  meanAndSd[j,13]=meth
  j=j+1
}

toCompare <- list(c("or2","pval"),c("or2","toComp"),c("pval","toComp"),c("prev","toComp"))
toCompareNames <- list(c("Odds Ratio","-log10(Pval)"),c("Odds Ratio","AUC Imp"),
                       c("-log10(Pval)","AUC Imp"),c("Prev","AUC Imp"))
for(i in 1:length(toCompare)){
  currComp <- toCompare[[i]]
  currNames <- toCompareNames[[i]]
  plotDf <- meanAndSd[,c(grep(paste(currComp,collapse="|"),colnames(meanAndSd)),ncol(meanAndSd))]
  colnames(plotDf) <- c("mean1","mean2","sd1","sd2","method")
  
  thePlot=ggplot(plotDf,aes(mean1,mean2,color=method))+geom_point()+
    geom_errorbar(aes(ymin=mean2-sd2,ymax=mean2+sd2,width=0.1))+
    geom_errorbarh(aes(xmin=mean1-sd1,xmax=mean1+sd1))+
    labs(x=currNames[1],y=currNames[2])
  plot(thePlot)
}

#### now comparing the ranks #######
allCrit <- c("toCompMean","pvalMean","or2Mean","prevMean")
allNames <- c("AUC Imp", "Pval","OR","R2")

for(i in 1:length(allCrit)){
  currCrit <- allCrit[i]
  currName <- allNames[i]
  yInds <- which(colnames(meanAndSd) %in% allCrit & colnames(meanAndSd) != currCrit)
  xInds <- which(colnames(meanAndSd) == currCrit)
  j <- 1
  plotList <- list()
  for(y in yInds){
    plotDf <- meanAndSd[,c(xInds,y)]
    plotDf[,1] <- rank(plotDf[,1])
    plotDf[,2] <- rank(plotDf[,2])
    colnames(plotDf) <- c("x1","x2")
    
    plotList[[j]] <- ggplot(plotDf, aes(x1,x2)) + geom_point()+
      labs(y=currName,x=allNames[-which(allNames==currName)][j])
    j <- j + 1
  }
  
  plotDf <- data.frame(apply(meanAndSd[c(xInds,yInds)],2,rank))
  colnames(plotDf) <- c("main",allNames[-i])
  plotDf <- melt(plotDf,id.vars = "main")
  thePlot <- ggplot(plotDf,aes(main,value,color=variable))+geom_point()+
    labs(x=paste(allNames[i],"Rank"),y="Other Rank",color="Other Statistic")
  plot(thePlot)
  
  grid.arrange(plotList[[1]], plotList[[2]], plotList[[3]])
}

#######################################################################################################################
#Compare to meta information by method
#columns go: sample size, total snps, number of cases

ukbbSampSize <- readRDS("caseControl.RDS")
#ukbbSampSize <- ukbbSampSize[,-ncol(ukbbSampSize)]
ukbbSampSize <- ukbbSampSize[,colnames(ukbbSampSize) %in% traitNames]
metaData=read.table("metaData",stringsAsFactors = F)
allH2=read.table("allH2",stringsAsFactors = F)

allMethod <- unique(truCompHolder[[3]]$methods)
methTraitDf <- data.frame(matrix(0,nrow=length(authors),ncol=length(allMethod)))
colnames(methTraitDf) <- allMethod
rownames(methTraitDf) <- authors
k <- 1
for(meth in allMethod){
  j <- 1
  for(i in seq(1,length(truCompHolder), compLength)){
    if(sum(truCompHolder[[i+2]]$methods==meth) > 0){
      methTraitDf[j,k] <- max(truCompHolder[[i+2]][truCompHolder[[i+2]]$methods==meth,12])
    } else {
      methTraitDf[j,k] <- NA
    }
    j <- j + 1
  }
  k <- k + 1
}

i <- 1
regressRes <- rep(list(data.frame(matrix(0,nrow=ncol(methTraitDf),ncol=4, dimnames = list(NULL,c("coef","se","tstat","pval"))),allMethod)), 7)
for(meth in allMethod){
  print(meth)
  forRegress <- data.frame(methTraitDf[,colnames(methTraitDf)==meth], 
                           metaData[metaData[,1] %in% rownames(methTraitDf),],
                           allH2[allH2[,1] %in% rownames(methTraitDf),2], 
                           metaData[metaData[,1] %in% rownames(methTraitDf),2]/metaData[metaData[,1] %in% rownames(methTraitDf),4],
                           ukbbSampSize[1,], ukbbSampSize[1,]/ukbbSampSize[2,])
  colnames(forRegress) <- c("perf", "author", "sampleSize","snps","cases","herit","caseRatio","ukbbCases","ukbbRatio")
  forRegress <- forRegress[!is.na(forRegress$perf) & !is.infinite(forRegress$perf),] 
  for(j in 3:9){
    regressRes[[j-2]][i,1:4] <- summary(lm(forRegress$perf ~ forRegress[,j]))$coefficients[2,]
    regressRes[[j-2]]$stat <- colnames(forRegress)[j]
  }
  i <- i + 1
}

plotDf <- do.call("rbind",regressRes)
plotDf$allMethod <- factor(as.character(plotDf$allMethod), plotDf$allMethod[plotDf$stat=="herit"][order(-log10(plotDf$pval[plotDf$stat=="herit"]))])
thePlot <- ggplot(plotDf, aes(-log10(pval),allMethod,color=stat))+geom_point()+
  labs(y="Method",x="-log10(Pval)")
plot(thePlot)
plotDf <- plotDf[plotDf$stat %in% c("sampleSize","herit","snps"),]
thePlot <- ggplot(plotDf, aes(-log10(pval),allMethod,color=stat))+geom_point()+
  labs(y="Method",x="-log10(Pval)")
plot(thePlot)

plotDf <- do.call("rbind",regressRes)
plotDf$allMethod <- factor(as.character(plotDf$allMethod), plotDf$allMethod[plotDf$stat=="herit"][order(plotDf$coef[plotDf$stat=="herit"])])
levels(plotDf$allMethod) <- levels(plotDf$allMethod)[order(regressRes[[3]]$pval)]
thePlot <- ggplot(plotDf, aes(coef,allMethod,color=stat))+geom_point()+
  labs(y="Method",x="Coef")
plot(thePlot)


#####PCA Plot
for(i in 1:ncol(methTraitDf)){
  methTraitDf[is.infinite(methTraitDf[,i]) | is.na(methTraitDf[,i]),i] <- 
    mean(methTraitDf[!is.infinite(methTraitDf[,i]) & !is.na(methTraitDf[,i]),i])
}
plotDf <- data.frame(prcomp(t(methTraitDf))$x)
plotDf$method <- rownames(plotDf)
thePlot <- ggplot(plotDf, aes(PC1, PC2))+geom_point()+geom_label_repel(aes(label=method),size=3)
plot(thePlot)


#############################################################################################################################
# PCA of all the scores (now moved to corr)

#set up allTrainScores for later use
covarDefs <- read.table("covarDefs", stringsAsFactors = F, header = T)
covar1 <- readRDS("covars.1.RDS")
eidIndex <- read.table("scott/scott.phenIndex.1")
covar1 <- covar1[eidIndex[,1],]

allTrainScores=data.frame(matrix(0,nrow=nrow(truCompHolder[[4]]),ncol=length(traitNames)))
colnames(allTrainScores)=traitNames
j=1
for(i in seq(1,length(truCompHolder), compLength)){
  if(nrow(truCompHolder[[i+3]]) < nrow(allTrainScores)){
    chosenSex <- covarDefs[covarDefs[,1]==authors[j],3]
    allTrainScores[covar1[,3]==chosenSex,j] <- truCompHolder[[i+3]][,8]
  } else {
    allTrainScores[,j] <- truCompHolder[[i+3]][,8]
  }
  j=j+1
}




dev.off()






# PHENO PHENO PHENO PHENO PHENO PHENO PHENO PHENO PHENO PHENO PHENO PHENO PHENO PHENO PHENO PHENO PHENO PHENO PHENO PHENO PHENO
###########################################################################################################################
###########################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pdf(paste0(rdsHome,"/metaPlots.pheno.pdf"))

# Table Describing UK Biobank

#############################################################################################################################
#Distributions of Case and Control
covar2 <- readRDS("covars.2.RDS")
eidIndex <- read.table("scott/scott.phenIndex.2")
covar2 <- covar2[eidIndex[,1],]

testScores=data.frame(matrix(0,nrow=(nrow(truPhenoHolder[[9]])*length(traitNames)),ncol=3))
allTestScores=data.frame(matrix(0,nrow=nrow(truPhenoHolder[[9]]),ncol=ncol(allTrainScores)))
colnames(testScores)=c("score","status","trait")
colnames(allTestScores) <- traitNames
j=1
k=1
countingCases <- rep(0, length(traitNames))
names(countingCases) <- traitNames
cc <- 1
for(i in seq(1,length(truPhenoHolder), phenoLength)){
  testDf=truPhenoHolder[[i+8]]
  countingCases[cc] <- sum(testDf$response != 0)
  cc <- cc + 1
  if(nrow(testDf) < nrow(allTestScores)){
    chosenSex <- covarDefs[covarDefs[,1]==authors[k],3]
    allTestScores[covar2[,3]==chosenSex,k] <- testDf$predictor
  } else {
    allTestScores[,k] <- testDf$predictor
  }
  
  
  testScores[j:(j+nrow(testDf)-1),1:2]=testDf[,(ncol(testDf)-1):ncol(testDf)]
  testScores[j:(j+nrow(testDf)-1),3]=rep(traitNames[k],nrow(testDf))
  j=j+nrow(testDf)
  k=k+1
}

savePlots[[6]] <- countingCases

testScores <- testScores[testScores$trait!=0,]
testScores$status=as.factor(testScores$status)
testScores$trait=as.factor(testScores$trait)
thePlot=ggplot(testScores,aes(x=score,y=trait,fill=status,alpha=0.4))+geom_density_ridges()+
  guides(fill=F)+guides(alpha=F)
plot(thePlot)

graphAbstractPlot2 <- ggplot(testScores,aes(score,fill=status))+geom_density(alpha=0.3)+
  labs(y="Density",x="Polygenic Risk Score")+
  scale_fill_discrete(name="Disease Status",breaks=c(1,2),labels=c("No Disease","Disease"))+
  guides(color=F)+theme(legend.position = c(0.7, 0.7))
plot(graphAbstractPlot2)



#############################################################################################################################
#Odds Ratio

allOR=data.frame(matrix(0,nrow=6,ncol=(length(traitNames)+1)))
allOR[,1]=factor(c(0.5,0.2,0.1,0.05,0.01,0.005),levels=c(0.5,0.2,0.1,0.05,0.01,0.005))
colnames(allOR)=c("cutOffs",traitNames)
j=2
for(i in seq(1,length(truPhenoHolder),phenoLength)){
  allOR[,j]=truPhenoHolder[[i+5]][,3]
  j=j+1
}

keepOR=allOR[4,]
avgOr=apply(allOR[,2:ncol(allOR)],2,mean)
keepNames=names(sort(avgOr,decreasing = T))[1:10]
allOR=melt(allOR,id.vars = "cutOffs")
allOR$variable=as.character(allOR$variable)
allOR$variable[!(allOR$variable %in% keepNames)]="Other"
thePlot=ggplot(allOR,aes(cutOffs,value))+
  geom_point(data=allOR[allOR$variable=="Other",],aes(cutOffs,value),color="grey80")+
  geom_point(data=allOR[allOR$variable!="Other",],aes(cutOffs,value,color=variable))+
  labs(x="Cut Off",y="Odds Ratio")
plot(thePlot)


#############################################################################################################################
#Prevalance
allPrev=data.frame(matrix(0,nrow=100,ncol=(length(traitNames)+1)))
allPrev[,1]=1:100
colnames(allPrev)=c("percentiles",traitNames)
j=2
for(i in seq(1,length(truPhenoHolder), phenoLength)){
  allPrev[,j]=truPhenoHolder[[i+4]][,2]
  allPrev[,j]=allPrev[,j]/sum(allPrev[,j])
  j=j+1
}
mult=apply(allPrev[91:100,],2,mean)/apply(allPrev[1:10,],2,mean)
mult=mult[2:length(mult)]
keepNames=names(sort(mult,decreasing = T))[1:8]
prevSSR=rep(0,length(traitNames))
for(i in 2:ncol(allPrev)){
  #xVal=log(allPrev[,i])
  #if(any(!is.finite(xVal))){
  #  startVal=max(which(!is.finite(xVal)))+1
  #} else {
  #  startVal=1
  #}
  prevModel <- lm(log(allPrev[!is.infinite(log(allPrev[,i])),i]) ~ allPrev$percentiles[!is.infinite(log(allPrev[,i]))])
  prevSSR[i-1] = summary(prevModel)$adj.r.squared
}

allPrev=melt(allPrev,id.vars = "percentiles")
colnames(allPrev)[2]="Trait"
thePlot=ggplot(allPrev,aes(percentiles,value,group=Trait))+
  geom_smooth(data=allPrev[!(allPrev$Trait %in% keepNames),],aes(percentiles,value),se=F,color="grey80")+
  geom_smooth(data=allPrev[allPrev$Trait %in% keepNames,],aes(percentiles,value,color=Trait),se=F)
plot(thePlot)


#Just looking at R2
j <- 1
nagelTest <- data.frame(base=rep(0,length(traitNames)), full=rep(0,length(traitNames)))
for(i in seq(1,length(truPhenoHolder), phenoLength)){
  testDf <- truPhenoHolder[[i + 8]]
  testDf <- testDf[,-which(colnames(testDf)=="array")]
  baseModel <- glm(response ~ ., testDf[,-which(colnames(testDf)=="predictor")], family="binomial")
  fullModel <- glm(response ~ ., testDf, family = "binomial")
  nagelTest[j,1] <- NagelkerkeR2(baseModel)$R2
  nagelTest[j,2] <- NagelkerkeR2(fullModel)$R2
  j <- j + 1
}
nagelTest$diff <- nagelTest$full - nagelTest$base

nagelTest$trait <- factor(traitNames, levels=traitNames[order(nagelTest$diff)])
thePlot <- ggplot(nagelTest,aes(trait,diff))+geom_point()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8))+
  labs(x="Trait",y="Nagelgerke R2 Improvement")
plot(thePlot)

nagelTest$trait <- factor(traitNames, levels=traitNames[order(nagelTest$full)])
thePlot <- ggplot(nagelTest,aes(trait,full))+geom_point()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8))+
  labs(x="Trait",y="Nagelgerke R2")
plot(thePlot)

#############################################################################################################################
#AUC Curve Statistical Difference -  DO NOT RUN UNTIL FINAL - IT TAKES A LONG TIME - BUT IT DOES WORK
# aucPval=data.frame(matrix(0,nrow=length(traitNames),ncol=2))
# colnames(aucPval) <- c("pval","trait")
# aucPval$trait <- traitNames
# 
# j=1
# for(i in seq(1,length(truPhenoHolder), phenoLength)){
#   aucPval[j,1] <- roc.test(truPhenoHolder[[i+9]][[2]], truPhenoHolder[[i+9]][[3]])$p.value
#   j <- j +1
# }
# 
# aucPval$trait=factor(aucPval$trait,levels = aucPval$trait[order(aucPval$pval)])
# thePlot=ggplot(aucPval,aes(trait,-log10(pval)))+geom_point()+
#   coord_flip() + labs(x="Trait",y="-log10(Pval of ROC Curve Difference)")
# plot(thePlot)

#############################################################################################################################
#AUC Improvement
allAUC=data.frame(matrix(0,nrow=length(traitNames),ncol=3))
allAUC$trait=traitNames
colnames(allAUC)=c("lo","val","hi","trait")
allAUCImp=data.frame(allAUC)
j=1
for(i in seq(1,length(truPhenoHolder), phenoLength)){
  allAUCImp[j,1:3]=truPhenoHolder[[i+1]]-truPhenoHolder[[i+2]]
  allAUC[j,1:3]=truPhenoHolder[[i+1]]
  j=j+1
}


allAUC$trait=factor(allAUC$trait,levels = allAUC$trait[order(allAUC$val)])
thePlot=ggplot(allAUC,aes(trait,val))+geom_point()+
  coord_flip()+
  geom_errorbar(aes(ymin=lo,ymax=hi,width=0.2))+
  labs(x="Trait",y="Total AUC")
plot(thePlot)
savePlots[[7]] <- allAUC

allAUCImp$trait=factor(allAUCImp$trait,levels = allAUCImp$trait[order(allAUCImp$val)])
thePlot=ggplot(allAUCImp,aes(trait,val))+geom_point()+
  coord_flip()+
  geom_errorbar(aes(ymin=lo,ymax=hi,width=0.2))+
  labs(x="Trait",y="AUC Improvement")
plot(thePlot)
savePlots[[8]] <- allAUCImp

graphAbstractPlot3 <- ggplot(allAUCImp,aes(val))+geom_histogram(bins=5)+
  labs(x="AUC Improvement",y="Number of Traits")
plot(graphAbstractPlot3)

#############################################################################################################################
#ROC
allROC=data.frame(matrix(0,nrow=(1000*length(traitNames)),ncol=3))
colnames(allROC)=c("trait","tpr","fpr")
jStart=1
jTrait=1
for(i in seq(1,length(truPhenoHolder),phenoLength)){
  inputROC=truPhenoHolder[[i]][truPhenoHolder[[i]][,3]=="PRS only",1:2]
  allROC[jStart:(jStart+nrow(inputROC)-1),2:3]=inputROC
  allROC[jStart:(jStart+nrow(inputROC)-1),1]=traitNames[jTrait]
  jStart=jStart+nrow(inputROC)+1
  jTrait=jTrait+1
}
allROC=allROC[allROC$trait!=0,]
keepNames=as.character(allAUC$trait[order(allAUC$val,decreasing = T)][1:5])


thePlot=ggplot(allROC,aes(fpr,tpr))+
  geom_line(data=allROC[allROC$trait %in% keepNames,],aes(fpr,tpr,color=trait),size=2)+
  geom_line(data=allROC[!(allROC$trait %in% keepNames),], aes(fpr,tpr,group=trait), color="grey80")+
  labs(x="False Positive Rate",y="True Positive Rate")+geom_abline(slope = 1)
plot(thePlot)


###########################################################################################################################
#Meta data plots showing sample size against AUC

metaData=read.table("metaData",stringsAsFactors = F)
authorNames=character(length=length(traitNames))
allSampSize=numeric(length=length(traitNames))
j=1
for(i in 1:length(authors)){
  allSampSize[j] <- metaData[metaData[,1]==authors[j],2]
  j=j+1
}

sampSizeCor <- round(cor(allAUCImp$val,allSampSize),3)
allAUCImp$sampSize=allSampSize
thePlot=ggplot(allAUCImp,aes(sampSize,val))+geom_point()+
  labs(y="AUC Improvement",x="Sample Size",caption=paste("Correlation =",sampSizeCor))
plot(thePlot)

###########################################################################################################################
#Meta data plots showing sample size against AUC

allH2=read.table("allH2",stringsAsFactors = F)
allHerit=numeric(length=length(traitNames))
j=1
for(i in 1:length(authors)){
  print(i)
  allHerit[j] <- allH2[allH2[,1]==authors[j],2]
  j=j+1
}

heritCor <- round(cor(allAUCImp$val,allHerit),3)
allAUCImp$Herit=allHerit
thePlot=ggplot(allAUCImp,aes(Herit,val))+geom_point()+
  labs(y="AUC Improvement",x="Heritability",caption=paste("Correlation =",heritCor))
plot(thePlot)

###########################################################################################################################
#Consistency of statistics

pvalOnly=c()
for(i in seq(1,length(truPhenoHolder), phenoLength)){
  pvalOnly=c(pvalOnly, truPhenoHolder[[i+6]])  
}

allStats=data.frame(traitNames,pvalOnly,nagelTest$diff,allAUCImp$val,as.numeric(unname(keepOR[2:length(keepOR)])))
for(i in 2:5){
  allStats[,i]=rank(allStats[,i],ties.method = "random")
}
allStats$avg=apply(allStats[,2:5],1,mean)
allStats$lo=apply(allStats[,2:5],1,min)
allStats$hi=apply(allStats[,2:5],1,max)
allStats$traitNames=factor(allStats$traitNames, level=allStats$traitNames[order(allStats$avg)])
thePlot=ggplot(allStats,aes(traitNames,avg))+geom_point()+coord_flip()+geom_errorbar(aes(ymin=lo,ymax=hi))+
  labs(y="Rank of Trait by All Stats",x="Trait Names")
plot(thePlot)


allCaptions <- c("col normalized","row and col normalized")
for(i in 1:2){
  allStats=data.frame(pvalOnly, prevSSR, allAUCImp$val, as.numeric(unname(keepOR[2:length(keepOR)])))
  allStats[is.infinite(allStats[,1]),1] <- 100 
  allStats <- data.frame(apply(allStats, 2, function(x) (x-min(x))/(max(x)-min(x))))
  if(i==2){
    allStats <- data.frame(t(apply(allStats,1,function(x) (x-min(x))/(max(x)-min(x)))))
    for(j in !complete.cases(allStats)){
      allStats[j,] <- rep(0,4)
    }
  }
  colnames(allStats) <- c("Pval","R2 Imp","AUC Imp","OR")
  rownames(allStats) <- traitNames
  
  allStats <- allStats[hclust(dist(allStats))$order,hclust(dist(t(allStats)))$order]
  
  baseLayer=expand.grid(x=rownames(allStats),y=colnames(allStats))
  forGrid=cbind(baseLayer,unlist(allStats))
  colnames(forGrid)=c("addTrait","evalTrait","Pval")
  thePlot=ggplot(forGrid,aes(evalTrait,addTrait,fill=Pval))+geom_raster()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.text.y = element_text(size=8))+
    scale_fill_viridis()+labs(fill="Scaled\nStat",x="Metric",y="Trait",captions=allCaptions[i])
  plot(thePlot)
  
  thePlot <- ggplot(allStats,aes(prev,or))+geom_point()
}



###########################################################################################################################
#Testing and training comparison
trainTestAucs <- data.frame(traits = traitNames,train = rep(0,length(traitNames)), test=rep(0,length(traitNames)))
trainTestAucs$test <- allAUCImp$val
for(testTrait in traitNames){
  trainTestAucs[trainTestAucs$traits==testTrait,2] <- max(splitMethods$toComp[splitMethods$trait==testTrait])
}
trainTestAucs$traits[abs(trainTestAucs[,3]-trainTestAucs[,2])<0.05] <- NA

thePlot=ggplot(trainTestAucs,aes(train,test))+geom_point()+geom_abline(slope = 1)+
  geom_label_repel(aes(label=traits),size=3)+labs(x="Train AUC Imp",y="Test AUC Imp",title="Comparison Test and Train")
plot(thePlot)
savePlots[[9]] <- trainTestAucs


dev.off()




if(stopAfterPheno){stop()}



pdf(paste0(rdsHome,"/metaPlots.corr.pdf"))

#CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
#CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR CORR

corrHeatmapFunc <- function(plotDf,theTitle,xname="Polygenic Risk Score",
                            yname="Polygenic Risk Score",fillname="AUC Imp",theCaption=NULL){
  hierClust1=hclust(dist(plotDf))
  hierClust2=hclust(dist(t(plotDf)))
  plotDf <- plotDf[hierClust1$order, hierClust2$order]
  baseLayer=expand.grid(x=rownames(plotDf),y=colnames(plotDf))
  forGrid=cbind(baseLayer,unlist(plotDf))
  colnames(forGrid)=c("addTrait","evalTrait","AUCImp")
  
  #theScale <- ((-50:50)/100)^3 + 0.1*((-50:50)/100)
  #theScale <- rescale(theScale,0:1)
  thePlot=ggplot(forGrid,aes(evalTrait,addTrait,fill=AUCImp))+geom_raster()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10), axis.text.y = element_text(size=10))+
    labs(x=xname,y=yname,title=theTitle,fill=fillname,caption=theCaption)+scale_fill_viridis_c()
  plot(thePlot)
}

#Correlation Calculations ########################################################################################
#Between the best Scores
prsCorr=data.frame(matrix(0,nrow=ncol(allTrainScores),ncol=ncol(allTrainScores)))
colnames(prsCorr)=colnames(allTrainScores)
rownames(prsCorr)=colnames(allTrainScores)
j=1
for(i in 1:ncol(allTrainScores)){
  prsCorr[i,]=apply(allTrainScores,2,function(x) cor(x,allTrainScores[,i]))
  prsCorr[i,j]=0
  j=j+1
}

hierClust=hclust(dist(prsCorr))
prsCorr <- prsCorr[hierClust$order, hierClust$order]

baseLayer=expand.grid(x=rownames(prsCorr),y=colnames(prsCorr))
forGrid=cbind(baseLayer,unlist(prsCorr))
colnames(forGrid)=c("addTrait","evalTrait","AUCImp")

theScale <- ((-50:50)/100)^3 + 0.1*((-50:50)/100)
theScale <- rescale(theScale,0:1)
#plot(1:length(theScale), theScale)

thePlot=ggplot(forGrid,aes(evalTrait,addTrait,fill=AUCImp))+geom_raster()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10), axis.text.y = element_text(size=10))+
  labs(x="Polygenic Risk Score",y="Polygenic Risk Score",title="Disease PRS Correlations",caption="All best scores")+
  scale_fill_viridis_c(values=theScale)
plot(thePlot)



#Matched method along the row
prsCorr=data.frame(matrix(0,nrow=length(traitNames),ncol=nrow(truCorrHolder[[2]])))
rownames(prsCorr)=traitNames
colnames(prsCorr)=truCorrHolder[[2]]$trait
j=1
for(i in seq(1,length(truCorrHolder), corrLength)){
  prsCorr[j,] <- truCorrHolder[[i+1]][,1]
  j <- j + 1
}

corrHeatmapFunc(prsCorr,"Disease PRS Correlations",theCaption="Row matched method score")



#Between the supplamentary traits
prsCorr=data.frame(matrix(0,nrow=length(traitNames),ncol=nrow(truCorrHolder[[1]])))
rownames(prsCorr)=traitNames
colnames(prsCorr)=truCorrHolder[[1]]$trait
j=1
for(i in seq(1,length(truCorrHolder), corrLength)){
  prsCorr[j,] <- truCorrHolder[[i]][,1]
  j <- j + 1
}

corrHeatmapFunc(prsCorr,"Supp PRS Correlations")




# Unrelated Disease Prediction ##############################################################################################
#Self-Assesed Diseases
allNames <- names(truCorrHolder[[3]])
for(i in seq(1,length(truCorrHolder), corrLength)){
  allNames <- intersect(allNames,names(truCorrHolder[[i+2]]))
}
allNames <- allNames[-which(allNames %in% c("unclassifiable",""))]
statHolder <- data.frame(matrix(0,nrow=length(traitNames),ncol=length(allNames),dimnames = list(traitNames,allNames)))
statList <- rep(list(statHolder),4)
k <- 1
for(i in seq(1,length(truCorrHolder), corrLength)){
  pullList <- truCorrHolder[[i+2]]
  l <- 1
  for(thisName in allNames){
    statList[[1]][k,l] <- pullList[[thisName]][[1]]
    statList[[2]][k,l] <- pullList[[thisName]][[2]][2]
    statList[[3]][k,l] <- pullList[[thisName]][[3]][4]
    l <- l + 1
  }
  k <- k + 1
}

critNames <- c("-log10(pval)","AUC Imp","OR")
for(i in 1:3){
  if(length(which(apply(statList[[i]],2,function(x) sum(is.infinite(x)))>0))){
    statList[[i]] <- statList[[i]][,-which(apply(statList[[i]],2,function(x) sum(is.infinite(x)))>0)]
  }
  corrHeatmapFunc(statList[[i]],"","Evaluated Trait","Polygenic Risk Score",critNames[i],theCaption="Self-reported")
  
  
  subPlot <- statList[[i]][,apply(statList[[i]], 2, function(x) sum(x %in% sort(as.vector(as.matrix((statList[[i]]))),decreasing = T)[1:20])) > 0]
  corrHeatmapFunc(subPlot,"","Evaluated Trait","Polygenic Risk Score",critNames[i],theCaption="Self-reported")
}

pvalSigCount <- apply(statList[[1]],1,function(x) sum(x>2))


#ICD Diseases
allNames <- names(truCorrHolder[[4]])
for(i in seq(1,length(truCorrHolder), corrLength)){
  allNames <- intersect(allNames,names(truCorrHolder[[i+3]]))
}
allNames <- allNames[-which(allNames %in% c("unclassifiable",""))]
statHolder <- data.frame(matrix(0,nrow=length(traitNames),ncol=length(allNames),dimnames = list(traitNames,allNames)))
statList <- rep(list(statHolder),4)
k <- 1
for(i in seq(1,length(truCorrHolder), corrLength)){
  pullList <- truCorrHolder[[i+3]]
  l <- 1
  for(thisName in allNames){
    statList[[1]][k,l] <- pullList[[thisName]][[1]] #the order is pval,auc imp, or, prev
    statList[[2]][k,l] <- pullList[[thisName]][[2]][2]
    statList[[3]][k,l] <- pullList[[thisName]][[3]][4]
    l <- l + 1
  }
  k <- k + 1
}

critNames <- c("-log10(pval)","AUC Imp","OR")
for(i in 1:3){
  if(length(which(apply(statList[[i]],2,function(x) sum(is.infinite(x)))>0))){
    statList[[i]] <- statList[[i]][,-which(apply(statList[[i]],2,function(x) sum(is.infinite(x)))>0)]
  }
  corrHeatmapFunc(statList[[i]],"","Evaluated Trait","Polygenic Risk Score",critNames[i],theCaption="ICD Diagnosis")
  

  subPlot <- statList[[i]][,apply(statList[[i]], 2, function(x) sum(x %in% sort(as.vector(as.matrix((statList[[i]]))),decreasing = T)[1:20])) > 0]
  corrHeatmapFunc(subPlot,"","Evaluated Trait","Polygenic Risk Score",critNames[i],theCaption = "ICD Diagnosis")
}

#Number Significant traits associated with score
pvalSigCount <- apply(statList[[1]],1,function(x) sum(x > 2)) + pvalSigCount
plotDf <- data.frame(pvalSigCount,trait=as.character(names(pvalSigCount)))
plotDf$trait <- factor(plotDf$trait, levels=plotDf$trait[order(plotDf$pvalSigCount)])
thePlot <- ggplot(plotDf,aes(trait,pvalSigCount))+geom_col()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y="Number Traits with \nAssociation Pval < 0.01",x="PRS")
plot(thePlot)

# #2d correlation plots##########################################################################################################
comboTrainTest <- rbind(allTrainScores,allTestScores)

#pca
pcaOut <- prcomp(t(comboTrainTest))
pca=data.frame(pcaOut$x)
pca$trait=colnames(allTrainScores)
pca$trait[pca[,2] > mean(pca[,2])-sd(pca[,2]) & pca[,2] < mean(pca[,2])+sd(pca[,2]) &
            pca[,1] > mean(pca[,1])-sd(pca[,1]) & pca[,1] < mean(pca[,1])+sd(pca[,1])] <- NA
thePlot=ggplot(pca,aes(PC1,PC2))+geom_point()+
  geom_label_repel(aes(label=trait),size=3)+
  labs(x=paste0("PC1 (",summary(pcaOut)$importance[2,1],")"),
       y=paste0("PC2 (",summary(pcaOut)$importance[2,2],")"))
plot(thePlot)


#Rtsne
outTsne <- Rtsne(t(comboTrainTest), perplexity = 2.5)
tsne <- data.frame(outTsne$Y)
tsne$trait <- colnames(allTrainScores)
thePlot=ggplot(tsne,aes(X1,X2))+geom_point()+
  geom_label_repel(aes(label=trait),size=3)+
  labs(x="tSNE_1",y="tSNE_2")
plot(thePlot)


#umap
outUmap <- umap(t(comboTrainTest), metric="cosine")
umapDf <- data.frame(outUmap$layout)
umapDf$trait <- colnames(allTrainScores)
thePlot=ggplot(tsne,aes(X1,X2))+geom_point()+
  geom_label_repel(aes(label=trait),size=3)+
  labs(x="umap_1",y="umap_2")
plot(thePlot)


################### FOR POSTER ##################
set.seed(3)
df <- data.frame(x=c(rnorm(1000,mean=0),rnorm(1000,mean=0.2)),y=factor(c(rep(0,1000),rep(1,1000))))
df$x[df$y==1 & df$x > 0] <- df$x[df$y==1 & df$x > 0]*1.2
ggplot(df,aes(x,fill=y,color=y))+geom_density(alpha=0.3)+
  labs(y="Density",x="Polygenic Risk Score")+
  scale_fill_discrete(name="Disease Status",breaks=c(0,1),labels=c("No Disease","Disease"))+
  guides(color=F)+theme(legend.position = c(0.7, 0.7))



dev.off()
saveRDS(savePlots, paste0(rdsHome,"/plotDfs.RDS"))
