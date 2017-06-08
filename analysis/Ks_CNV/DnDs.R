
library("ggplot2")
library("dplyr")
library("reshape2")

######
#Read Data
######

DF <- read.csv("final_dNdS.csv", header = TRUE)


######
#Pruning
######

clusters <- dim(DF)[1] / 3
DF$X <- rep(1:clusters, each=3)
DF$`Species Comparison` <- as.factor(rep(c("G. raimondii vs. G. kirkii", "G. raimondii vs. K. drynarioides", "G. kirkii vs. K. drynarioides"), clusters))

filtered <- !(DF$X %in% DF[DF$dS > 0.6,]$X)
DF <- DF[filtered,]




######
#Dating
######

groups <- DF %>% group_by(`Species Comparison`)
summary.table <- summarize(groups, Mean = mean(dS), StDev = sd(dS), Median = median(dS))

palms.time <- 0.0000000026
brassica.time <- 0.000000015

Gossypioides.Kokia.time = c(summary.table[[1,4]] / (2*brassica.time) /1000000, summary.table[[1,4]] / (2*palms.time) /1000000)
Gossypium.Gossypioides.time = c(summary.table[[2,4]] / (2*brassica.time) /1000000, summary.table[[2,4]] / (2*palms.time) /1000000)
Gossypium.Kokia.time = c(summary.table[[3,4]] / (2*brassica.time) /1000000, summary.table[[3,4]] / (2*palms.time) /1000000)
Divergence.table <- rbind(Gossypioides.Kokia.time, Gossypium.Gossypioides.time, Gossypium.Kokia.time)
colnames(Divergence.table) <- c("Brassica Estimate", "Palms Estimate")


######
#Plotting
######
png("Figure_dNdS.png", 5000, 5000, pointsize=12, res=600)
Df.melted <- melt(DF, id.vars = "Species Comparison", measure.vars = c("dS", "dN"))
colnames(Df.melted)[3] <- c("Substitution Rate")


p <- ggplot(Df.melted, aes(x=`Species Comparison`, y = log(`Substitution Rate`), fill=variable)) + geom_violin(position = position_dodge(1)) + 
  theme(legend.position = c(0.95,0.9), legend.title = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(color = "black", size = 13)) 
p + stat_summary(fun.y=median, geom="point", size=1, color="white", position = position_dodge(1)) + scale_fill_manual(values = c("red3", "black"))

dev.off()
