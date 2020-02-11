library(ggplot2)
library(gridExtra)
library(cowplot)
library(viridis)
theme_set(theme_cowplot())

#the columns for each element in bigSetHolder are
#chr, pos, effect, pval, index, trait

bigSetHolder <- readRDS("bigSetHolder.RDS")

numSnps <- rep(0, length(bigSetHolder))
allTraits <- rep("0", length(bigSetHolder))
for(i in 1:length(bigSetHolder)){
  bigSetHolder[[i]] <- bigSetHolder[[i]][!is.na(bigSetHolder[[i]]$V8) & is.finite(bigSetHolder[[i]]$V8),]
  bigSetHolder[[i]]$ind2 <- (bigSetHolder[[i]]$index/max(bigSetHolder[[i]]$index)) * 100
  bigSetHolder[[i]]$V82 <- (bigSetHolder[[i]]$V8/max(bigSetHolder[[i]]$V8)) * 100
  numSnps[i] <- nrow(bigSetHolder[[i]])
  allTraits[i] <- bigSetHolder[[i]]$trait[1]
}
df <- do.call("rbind",bigSetHolder)

#grid of lines
allLinePlots <- list()
for(i in 1:length(bigSetHolder)){
  allLinePlots[[i]] <- ggplot(bigSetHolder[[i]], aes(index, V8)) + geom_line() +
    labs(title=allTraits[i]) +
    theme(axis.title.x=element_blank(), axis.title.y = element_blank()) +
    scale_y_continuous(breaks = round(seq(0, max(bigSetHolder[[i]]$V8)*1.1, length.out = 4), 2))
}
do.call("grid.arrange", c(allLinePlots, ncol=5))

#grid of hex
allHexPlots <- list()
for(i in 1:length(bigSetHolder)){
  allHexPlots[[i]] <- ggplot(bigSetHolder[[i]], aes(-log10(V10), V8)) + geom_hex() + scale_fill_viridis() +
    theme(axis.title.x=element_blank(), axis.title.y = element_blank(), legend.position = "none") +
    labs(title=allTraits[i]) + 
    scale_y_continuous(breaks = round(seq(0, max(bigSetHolder[[i]]$V8)*1.1, length.out = 4), 2))
}
do.call("grid.arrange", c(allHexPlots, ncol=5))


#split into 4 line plots
maxBreaks <- quantile(numSnps, c(0.25, 0.5, 0.75, 1))
minBreaks <- c(0, maxBreaks[1:3])
quartLinePlots <- list()
for(i in 1:4){
  pullTraits <- allTraits[numSnps > minBreaks[i] & numSnps <= maxBreaks[i]]
  subDf <- df[df$trait %in% pullTraits,]
  quartLinePlots[[i]] <- ggplot(subDf, aes(index, V8, color = trait)) + geom_line() +
    theme(axis.title.x=element_blank(), axis.title.y = element_blank()) +
    scale_y_continuous(breaks = round(seq(0, max(subDf$V8)*1.1, length.out = 5), 2))
}
do.call("grid.arrange", c(quartLinePlots, ncol=2))

#############################################################################################
#############################################################################################
#density of each trait
maxVal <- max(df$V8)
minVal <- min(df$V8)
allDens <- matrix(0, nrow = 1000, ncol = length(bigSetHolder))
allX <- matrix(0, nrow = 1000, ncol = length(bigSetHolder))
for(i in 1:length(bigSetHolder)){
  theDen <- density(bigSetHolder[[i]]$V8, n = 1000, from = minVal, to = maxVal)
  allDens[,i] <- theDen$y
  allDens[,i] <- (allDens[,i] - min(allDens[,i]))/(max(allDens[,i]) - min(allDens[,i]))
  allX[,i] <- theDen$x
  allX[allDens[,i] < 1e-4, i] <- NA
}
colnames(allDens) <- allTraits
absOrder <- order(apply(allX, 2, function(x) sum(!is.na(x))))

allDens <- melt(allDens)
allX <- melt(allX)
allDens$x <- allX$value
allDens$Var2 <- factor(allDens$Var2, levels=unique(allDens$Var2)[absOrder])
thePlot <- ggplot(allDens, aes(x, Var2, color = value)) + geom_point() + scale_colour_viridis_c()+
  labs(x="Absolute Effect", y="Trait", color="Trait-Normalized\nDensity")
plot(thePlot)



#############################################################################################
#############################################################################################
#density of each trait normalized
maxVal <- max(df$V82)
minVal <- min(df$V82)
allDens <- matrix(0, nrow = 1000, ncol = length(bigSetHolder))
allX <- matrix(0, nrow = 1000, ncol = length(bigSetHolder))
for(i in 1:length(bigSetHolder)){
  theDen <- density(bigSetHolder[[i]]$V82, n = 1000, from = minVal, to = maxVal)
  allDens[,i] <- theDen$y
  allDens[,i] <- (allDens[,i] - min(allDens[,i]))/(max(allDens[,i]) - min(allDens[,i]))
  allX[,i] <- theDen$x
  allX[allDens[,i] < 1e-4, i] <- NA
}
colnames(allDens) <- allTraits
absOrder <- order(apply(allX, 2, function(x) sum(!is.na(x))))

allDens <- melt(allDens)
allX <- melt(allX)
allDens$x <- allX$value
allDens$Var2 <- factor(allDens$Var2, levels=unique(allDens$Var2)[absOrder])
thePlot <- ggplot(allDens, aes(x, Var2, color = value)) + geom_point() + scale_colour_viridis_c()+
  labs(x="Absolute Normalized Effect", y="Trait", color="Trait-Normalized\nDensity")
plot(thePlot)




##########################################################################################################
##########################################################################################################
#decile of each trait
subDf <- df[df$trait == "sle",]
allQuants <- matrix(0, nrow = 5, ncol = length(bigSetHolder))
allCumQuants <- matrix(0, nrow = 5, ncol = length(bigSetHolder))
for(i in 1:length(bigSetHolder)){
  allQuants[,i] <- quantile(bigSetHolder[[i]]$V82, 1:5/5)
  cumVal <- cumsum(sort(bigSetHolder[[i]]$V82*0.01, decreasing = F))
  splitVals <- sum(bigSetHolder[[i]]$V82*0.01)*(1:5)/5
  for(j in 1:5){
    allCumQuants[j,i] <- which.min(abs(cumVal-splitVals[j]))/nrow(bigSetHolder[[i]])
  }
}

colnames(allQuants) <- allTraits
colnames(allCumQuants) <- allTraits
quantOrder <- order(apply(allQuants, 2, mean))
cumOrder <- order(apply(allCumQuants, 2, mean))
allQuants <- melt(allQuants)
allCumQuants <- melt(allCumQuants)

allQuants$Var2 <- factor(allQuants$Var2, levels = unique(allQuants$Var2)[quantOrder])
thePlot <- ggplot(allQuants, aes(value, Var2, color = as.factor(Var1))) + geom_point() + 
  scale_color_viridis(discrete = T) + labs(x = "Scaled Effect", y = "Trait", color = "Quartile")
plot(thePlot)

allCumQuants$Var2 <- factor(allCumQuants$Var2, levels = unique(allQuants$Var2)[cumOrder])
thePlot <- ggplot(allCumQuants, aes(value, Var2, color = as.factor(Var1))) + geom_point() + 
  scale_color_viridis(discrete = T) + labs(x = "Scaled Effect", y = "Trait", color = "Equal\nSplits")
plot(thePlot)

#bar plot of numSnps
plotDf <- data.frame(numSnps, allTraits)
plotDf$allTraits <- factor(allTraits, levels = allTraits[order(numSnps)])
thePlot <- ggplot(plotDf, aes(allTraits, log10(numSnps))) + geom_col() + coord_flip()
plot(thePlot)


  
