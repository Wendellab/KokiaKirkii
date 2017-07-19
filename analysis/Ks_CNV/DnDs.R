
library("ggplot2")
library("dplyr")
library("reshape2")

######
#Read Data
######

DF <- read.csv("final_dNdS.csv", header = TRUE)
Malvaceae_DF <- read.csv("Cotton_Chocolate_final_dNdS.csv", header = TRUE)


######
#Pruning
######

clusters <- dim(DF)[1] / 3
DF$X <- rep(1:clusters, each=3)
DF$`Species Comparison` <- as.factor(rep(c("G. raimondii vs. G. kirkii", "G. raimondii vs. K. drynarioides", "G. kirkii vs. K. drynarioides"), clusters))

filtered <- !(DF$X %in% DF[DF$dS > 0.6,]$X)
DF <- DF[filtered,]

Malvaceae_DF <- Malvaceae_DF[Malvaceae_DF$dS<3,]


######
#Dating
######

groups <- DF %>% group_by(`Species Comparison`)
dS.summary.table <- summarize(groups, Mean = mean(dS), StDev = sd(dS), Median = median(dS))
dN.summary.table <- summarize(groups, Mean = mean(dN), StDev = sd(dN), Median = median(dN))

malv.time.median <- median(Malvaceae_DF$dS) / (2*60000000)

Gossypioides.Kokia.time = c(dS.summary.table[[1,4]] / (2*malv.time.median) /1000000)
Gossypium.Gossypioides.time = c(dS.summary.table[[2,4]] / (2*malv.time.median) /1000000)
Gossypium.Kokia.time = c(dS.summary.table[[3,4]] / (2*malv.time.median) /1000000)
Divergence.table <- rbind(Gossypioides.Kokia.time, Gossypium.Gossypioides.time, Gossypium.Kokia.time)
colnames(Divergence.table) <- c("Malvaceae Estimate")


Df.melted <- melt(DF, id.vars = "Species Comparison", measure.vars = c("dS", "dN"))
colnames(Df.melted)[3] <- c("Synonymous Substitution Rate")

newDf.melted <- Df.melted[Df.melted$variable == "dS",]

######
#Plotting
######

png("Figure_dNdSboxplot.png", 5000, 5000, pointsize=12, res=600)
ggplot(Df.melted, aes(x=`Species Comparison`, y = `Synonymous Substitution Rate`, fill=variable)) + geom_boxplot(position = position_dodge(1), outlier.shape = NA) + 
  theme(legend.position = c(0.95,0.9), legend.title = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(color = "black", size = 13)) + 
  scale_y_continuous(limits = c(0,0.175)) +
  stat_summary(fun.y=median, geom="point", size=1, color="white", position = position_dodge(1)) + scale_fill_manual(values = c("red3", "black"))
dev.off()




png("Figure_dSDistribution.png", 5000, 5000, pointsize=12, res=600)
ggplot(newDf.melted, aes(newDf.melted$`Synonymous Substitution Rate`, colour = newDf.melted$`Species Comparison`)) + 
  geom_density(adjust = 1/5) + 
  theme(legend.position = c(0.8, 0.9), axis.title.y = element_blank(), legend.title = element_blank()) + 
  xlab("Synonymous Substitution Rate") 
dev.off()




png("Figure_Malvaceae_Divergence.png", 5000, 5000, pointsize=12, res=600)
ggplot(Malvaceae_DF, aes(Malvaceae_DF$dS)) + 
  geom_density(adjust = 1/5) + geom_vline(xintercept = median(Malvaceae_DF$dS)) +
  theme(axis.title.y = element_blank(), legend.title = element_blank()) + 
  xlab("Synonymous Substitution Rate") 
dev.off()

