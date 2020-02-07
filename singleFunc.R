runLikeTheWind <- function(authorComp,includeEnet,goNet){
  goShort=F
  sink(paste0(authorComp,".para.txt"), append=TRUE)
  
 
  #COMPARE
  finalCrit="AUC"
  plotAll=TRUE
  print(ls())
  source("compareEval.1.2.func.R")
  compReturn=compareEval(authorComp,finalCrit,plotAll,phenoDefs,fileDefs,trainPhenos,diseaseTrainList,
			cp,goShort,covarDefs,includeEnet,goNet)
  currentTrait=compReturn[[1]]
  bestFile=compReturn[[2]]
  trainDf=compReturn[[3]]
  decoder=compReturn[[4]]
  trainROC=compReturn[[5]]
  trainAUC=compReturn[[6]]
  toRev=compReturn[[7]]
  compResMat=compReturn[[8]]
  enetBetas=compReturn[[9]]
  compHolder=list(bestFile,currentTrait,compResMat,trainDf,authorComp,trainAUC,enetBetas)
  print("done compare")
  #rm(diseaseTrainList)
  
  #PHENOG
  #load("readData.testScores.RData")
  source("phenoEval.1.1.func.R")
  compReturn=phenoEval(authorComp,bestFile,trainDf,currentTrait,phenoDefs,
                       trainPhenos,testPhenos,diseaseTestList,toRev,goShort,covarDefs,enetBetas,decoder)
  testDf=compReturn[[1]]
  testROC=compReturn[[2]]
  testAUC=compReturn[[3]]
  covarAUC=compReturn[[4]]
  scoreOnlyAUC=compReturn[[5]]
  prevDf=compReturn[[6]]
  orDf=compReturn[[7]]
  pvalDf=compReturn[[8]]
  violinDf=compReturn[[9]]
  allGs=compReturn[[10]]
  phenoHolder=list(testROC,testAUC,covarAUC,scoreOnlyAUC,prevDf,orDf,pvalDf,violinDf,testDf,allGs)
  print("done pheno")

  #BOOST
  #compStat="AUC"
  #maxBest=10
  #removeSimilar=TRUE
  #load("readData.supports.RData")
  #source("boostEval.1.1.func.R")
  #compReturn=boostEval(authorComp,bestFile,compStat,trainDf,testDf, allSupportDecode, allDiseaseDecode,
  #                    trainROC,trainAUC,testROC,testAUC,maxBest,removeSimilar,decoder,currentTrait,covarAUC)
  #bestTracker=compReturn[[1]]
  #suppTrain=compReturn[[2]]
  #suppTest=compReturn[[3]]
  #boostRoc=compReturn[[4]]
  #boostAUC1=compReturn[[5]]
  #boostAUC2=compReturn[[6]]
  #allScoreStats=compReturn[[7]] #allResReturn - test - add score
  #allDiseaseStats=compReturn[[8]] #scoreStatsReturn - train - add trait
  #addedScoresDf=compReturn[[9]]
  #addedTraitsDf=compReturn[[10]]
  #scoreSplits=compReturn[[11]]
  #traitsSplits=compReturn[[12]]
  #testDiseaseStat=compReturn[[13]]
  #boostHolder=list(bestTracker,boostRoc,boostAUC1,boostAUC2,allScoreStats,allDiseaseStats,scoreSplits,traitsSplits,testDiseaseStat)
  #print("done boost")
 
  #CORR
  source("corrEval.1.2.func.R")
  compReturn=corrEval(authorComp,bestFile, trainDf, testDf, currentTrait, decoder,
           allSupportDecode, allDiseaseDecode, allScoresTrain, allScoresTest, allSupportsTrain, allSupportsTest,
  	   testROC,testAUC, enetBetas)
  corrSupport=compReturn[[1]]
  corrDisease=compReturn[[2]]
  diseaseComp=compReturn[[3]]
  supportComp=compReturn[[4]]
  corrHolder=list(corrSupport,corrDisease,diseaseComp,supportComp)
  print("done corr")


  rm(currentTrait,bestFile,trainDf,decoder,trainROC,trainAUC,testDf,testROC,testAUC,improveAUC,bestTracker,suppTrain,covar1,covar2,
     suppTest,boostRoc,boostAUC,corrSupport,corrDisease,diseaseComp,supportComp,diseaseTrainList,diseaseTestList,trainPhenos,testPhenos)
  outName <- paste(authorComp,"done","RDS",sep=".")
  saveRDS(list(compHolder,phenoHolder,corrHolder),outName)
  sink()
}
