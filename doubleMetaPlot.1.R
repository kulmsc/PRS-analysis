icdDfs <- readRDS("doneIcdPhenos/plotDfs.RDS")
diseaseDfs <- readRDS("doneDiseasePhenos/plotDfs.RDS")
medDfs <- readRDS("doneMedPhenos/plotDfs.RDS")
doubleDfs <- readRDS("doneDoublePhenos/plotDfs.RDS")
allDfs <- list(icdDfs, diseaseDfs, medDfs, doubleDfs)
allNames <- c("ICD", "Disease", "Meds", "Double")



pdf("phenMethodComp.pdf")

################################################################################
# Comparison of All
strucDfs <- list()
for(i in 1:length(allNames)){
  pullDf <- allDfs[[i]][[1]]
  pullDf <- pullDf[,colnames(pullDf) %in% c("methods", "toComp", "trait")]
  pullDf$type <- allNames[i]
  strucDfs[[i]] <- pullDf
}
plotDf <- do.call("rbind", strucDfs)
thePlot <- ggplot(plotDf, aes(type, toComp, fill=type)) + geom_boxplot() +
  labs(x="Phenotyping Method", y = "AUC Improvement") + guides(fill=FALSE)
plot(thePlot)

#############################################################################
#Average Across Traits
avgOverTraits <- data.frame(plotDf)
avgOverTraits$methods <- as.character(avgOverTraits$methods)
j <- 1
for(i in 1:length(unique(plotDf$trait))){
  for(phenMeth in allNames){
    avgOverTraits[j,2] <- mean(plotDf$toComp[plotDf$trait == unique(plotDf$trait)[i] & plotDf$type == phenMeth])
    avgOverTraits[j,1] <- "mean"
    avgOverTraits[j,3] <- unique(plotDf$trait)[i]
    avgOverTraits[j,4] <- phenMeth
    j <- j + 1
  }
}

avgOverTraits <- avgOverTraits[1:(j-1),]
thePlot <- ggplot(avgOverTraits, aes(type, toComp, color = trait)) + geom_point() + geom_line(aes(group=trait))+
  labs(x="Phenotyping Method", y = "AUC Improvement", title = "Average of Traits")
plot(thePlot)

#############################################################################
#Average Across Traits
avgOverMethods<- data.frame(plotDf)
avgOverMethods$methods <- as.character(avgOverMethods$methods)
j <- 1
for(i in 1:length(unique(as.character(plotDf$methods)))){
  for(phenMeth in allNames){
    avgOverMethods[j,2] <- mean(plotDf$toComp[plotDf$methods == unique(plotDf$methods)[i] & plotDf$type == phenMeth])
    avgOverMethods[j,3] <- "mean"
    avgOverMethods[j,1] <- unique(as.character(plotDf$methods))[i]
    avgOverMethods[j,4] <- phenMeth
    j <- j + 1
  }
}

avgOverMethods <- avgOverMethods[1:(j-1),]
thePlot <- ggplot(avgOverMethods, aes(type, toComp, color = methods)) + geom_point() + geom_line(aes(group=methods)) +
  labs(x="Phenotyping Method", y = "AUC Improvement", title = "Average of Methods")
plot(thePlot)

#################################################################################
#Comparison of Best Method
strucDfs <- list()
for(i in 1:length(allNames)){
  pullDf <- allDfs[[i]][[2]]
  pullDf <- pullDf[,colnames(pullDf) %in% c("methods", "toComp", "trait")]
  pullDf$type <- allNames[i]
  strucDfs[[i]] <- pullDf
}
plotDf <- do.call("rbind", strucDfs)
thePlot <- ggplot(plotDf, aes(type, toComp, fill=type)) + geom_boxplot() + geom_point() +
  labs(x="Phenotyping Method", y = "AUC Improvement", title = "Best Traits for each Method") + guides(fill=FALSE)
plot(thePlot)
thePlot <- ggplot(plotDf, aes(type, toComp)) + geom_line(aes(color = methods, group = methods)) +
  labs(x="Phenotyping Method", y = "AUC Improvement", title = "Best Traits for each Method")
plot(thePlot)


#################################################################################
#Comparison of each Trait - All Param
strucDfs <- list()
for(i in 1:length(allNames)){
  pullDf <- allDfs[[i]][[4]]
  pullDf <- pullDf[,colnames(pullDf) %in% c("methods", "toComp", "trait")]
  pullDf$type <- allNames[i]
  strucDfs[[i]] <- pullDf
}
plotDf <- do.call("rbind", strucDfs)
thePlot <- ggplot(plotDf, aes(type, toComp, fill=type)) + geom_boxplot() + geom_point() +
  labs(x="Phenotyping Method", y = "AUC Improvement", title = "Best Methods for each Trait") + guides(fill=FALSE)
plot(thePlot)
thePlot <- ggplot(plotDf, aes(type, toComp)) + geom_line(aes(color = trait, group = trait)) +
  labs(x="Phenotyping Method", y = "AUC Improvement", title = "Best Methods for each Trait")
plot(thePlot)


#########################################################################################
##########################################################################################
#Enet vs. Not Enet
strucDfs <- list()
for(i in 1:length(allNames)){
  pullDf <- allDfs[[i]][[5]]
  newCol <- pullDf$value[pullDf$Var2 == "enet"] - pullDf$value[pullDf$Var2 == "notEnet"]
  pullDf <- data.frame(enetDiff = newCol, trait = pullDf$Var1[1:length(newCol)])
  pullDf$type <- allNames[i]
  strucDfs[[i]] <- pullDf
}
plotDf <- do.call("rbind", strucDfs)

thePlot <- ggplot(plotDf, aes(type, enetDiff, fill = type)) + geom_boxplot() + guides(fill = FALSE)+
  labs(y = "Enet AUC - Best Single Method AUC", x="Phenotyping Method")
plot(thePlot)
thePlot <- ggplot(plotDf, aes(type, enetDiff, color = trait, group = trait)) + geom_point() + geom_line()+
  labs(y = "Enet AUC - Best Single Method AUC", x="Phenotyping Method")
plot(thePlot)
#Total AUC


#########################################################################################
##########################################################################################
#Case Control
strucDfs <- data.frame(matrix(0, nrow = length(names(allDfs[[2]][[6]])), ncol = 5))
colnames(strucDfs) <- c("traits", "icd", "disease", "meds", "double")
strucDfs$traits <- sort(names(allDfs[[2]][[6]]))
for(i in 1:length(allNames)){
  pullDf <- allDfs[[i]][[6]]
  pullDf <- pullDf[order(names(pullDf))]
  strucDfs[strucDfs$traits %in% names(pullDf),i+1] <- pullDf
}
plotDf <- melt(strucDfs, id.vars = "traits")
thePlot <- ggplot(plotDf, aes(variable, value)) + geom_boxplot()+
  labs(x = "Phenotyping Method", y = "Number of Cases")
plot(thePlot)
thePlot <- ggplot(plotDf, aes(traits, value, color=variable)) + geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_discrete(name = "Phenotyping\nMethod")+
  labs(x = "Trait", y = "Number of Cases")
plot(thePlot)


#########################################################################################
##########################################################################################
#Total AUC
strucDfs <- list()
for(i in 1:length(allNames)){
  pullDf <- allDfs[[i]][[7]]
  pullDf$type <- allNames[i]
  strucDfs[[i]] <- pullDf
}
plotDf <- do.call("rbind", strucDfs)
thePlot <- ggplot(plotDf, aes(type, val, fill = type)) + geom_boxplot() + guides(fill = FALSE)+
  labs(x="Phenotyping Method", y="Total AUC")
thePlot <- ggplot(plotDf, aes(type, val, color=trait)) + geom_point() + geom_line(aes(group=trait))+
  labs(x="Phenotyping Method", y="Total AUC")
plot(thePlot)

##########################################################################################
#Total AUC
strucDfs <- list()
for(i in 1:length(allNames)){
  pullDf <- allDfs[[i]][[8]]
  pullDf$type <- allNames[i]
  strucDfs[[i]] <- pullDf
}
plotDf <- do.call("rbind", strucDfs)
thePlot <- ggplot(plotDf, aes(type, val, fill = type)) + geom_boxplot() + guides(fill = FALSE)+
  labs(x="Phenotyping Method", y="AUC Improvment")
plot(thePlot)
thePlot <- ggplot(plotDf, aes(type, val, color=trait)) + geom_point() + geom_line(aes(group=trait))+
  labs(x="Phenotyping Method", y="AUC Improvment")
plot(thePlot)


dev.off()
