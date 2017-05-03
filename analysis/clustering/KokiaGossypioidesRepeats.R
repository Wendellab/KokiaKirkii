library(ggplot2)
library(scales)
library(ggrepel)
library(factoextra)
library(reshape2)
library(gridExtra)
library(geomorph)
library(testit)

#############################################
sessionInfo()
# R version 3.3.3 (2017-03-06)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 14393)
# 
# locale:
# [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252 LC_NUMERIC=C                           LC_TIME=English_United States.1252    
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] reshape2_1.4.2       gridExtra_2.2.1      factoextra_1.0.4     ggrepel_0.6.5        plotrix_3.6-4        BiocInstaller_1.22.3 scales_0.4.1         ggplot2_2.2.1       
# 
# loaded via a namespace (and not attached):
#  [1] Rcpp_0.12.10     digest_0.6.12    ggpubr_0.1.2     grid_3.3.3       plyr_1.8.4       gtable_0.2.0     magrittr_1.5     stringi_1.1.5    lazyeval_0.2.0   labeling_0.3     tools_3.3.3     
# [12] stringr_1.2.0    munsell_0.4.3    colorspace_1.3-2 tibble_1.3.0    
# 
##############################################

# import file with all counts
# evaluate if 0.01% cutoff is reasonable (here, cluster 274)

data <- read.table("all.counts", header = T, row.names=1, sep="\t")
data$size <- as.numeric(rowSums(data[,-1]))
data$percent <- cumsum(data$size)/sum(data$size)

png("cotton.cutoff.png", 5000, 5000, pointsize=12, res=600)
ggplot(data, aes(x=cluster, y=percent)) + geom_line(size=1) + geom_vline(xintercept=274, color='yellow3', size=1) + scale_x_log10(labels=comma) + scale_y_log10() + geom_vline(xintercept=0, color="grey")+ geom_hline(yintercept=0, color="grey")
dev.off()


################### ordination ###################
### PCoA of scaled data ###

annot_clust <- read.table("annotated.counts", header = T, row.names=1, sep="\t")

ord_table <- annot_clust[,(3:17)]
Mb_table <- ord_table*0.0095 
# 0.0095 multiplier represents # reads (x) * 95nt/read * 1 kb/1000nt * 1Mb/1000kb * 100% = # reads * 0.0095 = # Mb in entire genome for that cluster 

# convert to percent of genome to balance the numbers for ordination, asks are the relative proportions the same
# A1 (G. herbaceum) = 1667 Mb; sampling = 175474
# A2 (G. arboreum) = 1689 Mb; sampling = 180063
# D5 (G. raimondii) = 880 Mb; sampling = 92631
# G. kirkii = K. drynarioides = 587 Mb; sampling = 61,789 reads
perc_table <- data.frame( lapply(Mb_table[1:3], function(x) x/1667), lapply(Mb_table[4:8], function(x) x/1698), lapply(Mb_table[9:13], function(x) x/880), lapply(Mb_table[14:15], function(x) x/587) )

mydata <- t(perc_table)
scaledata <- scale(mydata, center = T, scale = T)
d <- dist(scaledata, method = "euclidean")

cmdfit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
x <- cmdfit$points[,1]
y <- cmdfit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS",	type="n")
text(x, y, labels = row.names(mydata), cex=.7)
points(x, y, pch=19) 

cmdpoints <- as.data.frame(cmdfit$points)
cmdpoints$genome <- sub("_.*", "", row.names(cmdpoints))
cmdpoints$species <- rownames(cmdpoints)

png("cotton.ordination.png", 5000, 5000, pointsize=12, res=600)
ggplot(cmdpoints, aes(x=V1, y=V2, color=genome, shape=genome)) + geom_point(size=2) + xlab("PCoA component 1") + ylab("PCoA component 2") + scale_color_manual(values=c("A1"="orchid", "A2"="orchid4", "D5"="slategrey","kirkii"="blue3", "kokia"="green3"))+ geom_text_repel(aes(label=species))
dev.off()

### PCA to get the variation ###

cluster.pca <- prcomp(mydata, scale = TRUE)
cluster.eig <- get_eigenvalue(cluster.pca)

ind.coord <- cluster.pca$x

subsections <- as.data.frame(mydata[,1:2])
subsections$sub <-c("A1", "A1", "A1", "A2", "A2", "A2", "A2", "A2", "D5", "D5", "D5", "D5", "D5", "kirkii", "kokia")

subfac <- as.factor(subsections[,3])

png("cotton_outgroup.PCA.direct.annot.png", 5000, 5000, pointsize=12, res=600)
fviz_pca_ind(cluster.pca, habillage=subfac, pointsize =2, invisible="quali", repel=TRUE, labelsize=3) + theme_minimal() + labs(title = "PCA of raw counts") + theme(axis.text = element_text(size = rel(1.5)), plot.margin=margin(2,2,2,2,"cm"), plot.title=element_text(face="bold", hjust=0.5), axis.title.x = element_text(face="bold", hjust=0.5), axis.title.y = element_text(face="bold", vjust=0.5), legend.position="none") +theme_set(theme_grey(base_size=12))+ scale_color_manual(breaks=c("A1", "A2", "D5","kirkii", "kokia_"), values=c("orchid", "orchid4", "slategrey","blue3", "green3"))
dev.off()

### PCoA of log-transformed data ###
# evaluates how close the overall repetitive profiles are to one another

ord_table[ord_table==0]=0.00000001 # no -Inf values in PCA, so make them really small here
logdata <- t(log(ord_table))

dlog <- dist(logdata, method = "euclidean")

logcmd <- cmdscale(dlog,eig=TRUE, k=2) # k is the number of dim
logx <- logcmd$points[,1]
logy <- logcmd$points[,2]
plot(logx, logy, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS",	type="n")
text(logx, logy, labels = row.names(logdata), cex=.7)
points(logx, logy, pch=19) 

logcmdpoints <- as.data.frame(logcmd$points)
logcmdpoints$genome <- sub("_.*", "", row.names(logcmdpoints))
logcmdpoints$species <- rownames(logcmdpoints)

png("cotton.ordination.log.png", 5000, 5000, pointsize=12, res=600)
ggplot(logcmdpoints, aes(x=V1, y=V2, color=genome, shape=genome)) + geom_point(size=2) + xlab("PCoA component 1") + ylab("PCoA component 2") + scale_color_manual(values=c("A1"="orchid", "A2"="orchid4", "D5"="slategrey","kirkii"="blue3", "kokia"="green3"))+ geom_text_repel(aes(label=species))
dev.off()


### Procrustes ANOVA and pairwise tests using complex linear models
### which species differ from one another in repeats
### https://groups.google.com/forum/#!topic/geomorph-r-package/8_B2thhxU_o

# Adams, D.C., and E. Otarola-Castillo. 2013. geomorph: an R package for the collection and analysis of geometric morphometric shape data. Methods in Ecology and Evolution. 4:393-399.
# Adams, D.C., M. L. Collyer, A. Kaliontzopoulou, and E. Sherratt. 2016 geomorph: Software for geometric morphometric analyses. R package version 3.0.2 http://cran.r-project.org/web/packages/geomorph/index.html.

# set factors for manova tests, same order is the table rows
facGroup <- as.factor(c("A", "A", "A", "A", "A", "A", "A", "A", "D", "D", "D", "D", "D", "OG", "OG"))
facSpecies <- as.factor(c("A1", "A1", "A1", "A2", "A2", "A2", "A2", "A2", "D5", "D5", "D5", "D5", "D5", "kirkii", "kokia"))

# procD.lm with data as percent genome size
species.diff <- advanced.procD.lm(mydata~facSpecies, ~1, groups=~facSpecies)

#Effect sizes (Z)
#             A1       A2        D5    kirkii     kokia
#A1     0.000000 0.575608 2.1022526 1.6452576 1.6624025
#A2     0.575608 0.000000 2.8185505 1.8842304 1.9055143
#D5     2.102253 2.818551 0.0000000 0.9200665 0.9360889
#kirkii 1.645258 1.884230 0.9200665 0.0000000 0.3545530
#kokia  1.662403 1.905514 0.9360889 0.3545530 0.0000000

#P-values
#          A1    A2    D5 kirkii  kokia
#A1     1.000 0.758 0.015 0.0260 0.0410
#A2     0.758 1.000 0.003 0.0010 0.0010
#D5     0.015 0.003 1.000 0.5520 0.5210
#kirkii 0.026 0.001 0.552 1.0000 0.7945
#kokia  0.041 0.001 0.521 0.7945 1.0000


# procD.lm with data as raw count data
mycountdata <- as.matrix(t(ord_table))
count.diff <- advanced.procD.lm(mycountdata~facSpecies, ~1, groups=~facSpecies)

#Effect sizes (Z)
#              A1        A2        D5    kirkii     kokia
#A1     0.0000000 0.6049583 2.2238744 1.6073137 1.6092942
#A2     0.6049583 0.0000000 2.9822639 1.8773175 1.8776222
#D5     2.2238744 2.9822639 0.0000000 0.4478466 0.4516000
#kirkii 1.6073137 1.8773175 0.4478466 0.0000000 0.1240738
#kokia  1.6092942 1.8776222 0.4516000 0.1240738 0.0000000

#P-values
#          A1    A2    D5 kirkii kokia
#A1     1.000 0.646 0.012  0.049 0.063
#A2     0.646 1.000 0.003  0.002 0.002
#D5     0.012 0.003 1.000  0.937 0.934
#kirkii 0.049 0.002 0.937  1.000 0.972
#kokia  0.063 0.002 0.934  0.972 1.000



########### characterize composition ###########

annot_clust$cluster <- NULL

# 9.5 multiplier represents # reads (x) * 95nt/read * 1 kb/1000nt * 100% = # reads * 9.5 = # Kb in entire genome for that class 
Kbamount <- data.frame(annot_clust[1], apply(annot_clust[2:16], 2, function (x) x*9.5))
KBsum <- aggregate(. ~Lineage, data=Kbamount, FUN=sum)

KBsum$A1 <- rowMeans(KBsum[,2:4])
KBsum$A2 <- rowMeans(KBsum[,5:9])
KBsum$D5 <- rowMeans(KBsum[,10:14])

KBsum$A1min <- apply(KBsum[,2:4], 1, min)
KBsum$A2min <- apply(KBsum[,5:9], 1, min)
KBsum$D5min <- apply(KBsum[,10:14], 1, min)
KBsum$kirkiimin <- KBsum$kirkii
KBsum$kokia_min <- KBsum$kokia_

KBsum$A1max <- apply(KBsum[,2:4], 1, max)
KBsum$A2max <- apply(KBsum[,5:9], 1, max)
KBsum$D5max <- apply(KBsum[,10:14], 1, max)
KBsum$kirkiimax <- KBsum$kirkii
KBsum$kokia_max <- KBsum$kokia_

KBsum <- KBsum[,-(2:14)]
KBm <- melt(KBsum[,-(7:16)])
min <- c(KBsum$kirkiimin, KBsum$kokia_min, KBsum$A1min, KBsum$A2min, KBsum$D5min)
max <- c(KBsum$kirkiimax, KBsum$kokia_max, KBsum$A1max, KBsum$A2max, KBsum$D5max)
KBm$min <- min
KBm$max <- max

limits <- aes(ymax=KBm$max, ymin=KBm$min)

dodge <- position_dodge(width=0.9)

png("Figure_TE.amounts.png", 7500, 5000, pointsize=12, res=600)

ggplot(KBm, aes(x=Lineage, y=value, fill = variable)) + geom_bar(stat = "identity",position = dodge) + scale_y_log10(labels=comma) + scale_fill_manual(breaks=c("kirkii", "kokia_", "A1", "A2", "D5"), values=c("blue3", "green3", "orchid", "orchid4", "slategrey")) + geom_errorbar(limits, position = dodge) + labs(title = "Aggregate amounts in each species", x="Broad element category", y="Aggregate amount (mean) in kilobases") + theme(axis.text = element_text(size = rel(1.5)), plot.margin=margin(2,2,2,2,"cm"), plot.title=element_text( face="bold", hjust=0.5), axis.title.x = element_text(face="bold", hjust=0.5), axis.title.y = element_text(face="bold", vjust=0.5))+theme_set(theme_grey(base_size=12))
dev.off()

sum(KBsum$kirkii)/1000 # convert to Mb
# [1] 110.3615
sum(KBsum$kokia_)/1000 # convert to Mb
# [1] 109.4685


########### compare kirkii vs kokia ###########

### plot a 1:1 line and display cluster comparisons relative to 1:1 expection
# remember, these have equivalent genome sizes

KdGk <- annot_clust[,c(1,15:16)]
KdGk$signK <- {ifelse((KdGk$kirkii > KdGk$kokia_), "positive", "negative")} 

png("Figure_TE.comparisons.png", 5000, 5000, pointsize=12, res=600)

p1 <- ggplot(KdGk, aes(x=kokia_, y=kirkii, shape=signK, color=signK)) + geom_point(size=2) + geom_abline(intercept=0, slope=1) + scale_color_manual(values=c("positive"="blue3", "negative"="green3")) +  scale_x_continuous(expand = c(0, 0), limits=c(0,750)) + scale_y_continuous(expand = c(0, 0), limits=c(0,750)) + theme(legend.position="none", axis.title.x = element_text(face="italic", hjust=0.5), axis.title.y = element_text(face="italic", vjust=0.5)) +labs(title = "All clusters", x="Kokia drynarioides", y="Gossypioides kirkii")
p2 <- ggplot(KdGk[KdGk$Lineage == "LTR", ], aes(x=kokia_, y=kirkii, shape= signK, color=signK)) + geom_point(size=2) + geom_abline(intercept=0, slope=1)+ scale_color_manual(values=c("positive"="blue3", "negative"="green3"))+ scale_x_continuous(expand = c(0, 0), limits=c(0,400)) + scale_y_continuous(expand = c(0, 0), limits=c(0,400))+ theme(legend.position="none", axis.title.x = element_text(face="italic", hjust=0.5), axis.title.y = element_text(face="italic", vjust=0.5)) +labs(title = "unclassified LTR elements", x="Kokia drynarioides", y="Gossypioides kirkii")
p3 <- ggplot(KdGk[KdGk$Lineage == "LTR/Gypsy", ], aes(x=kokia_, y=kirkii, shape= signK, color=signK)) + geom_point(size=2) + geom_abline(intercept=0, slope=1)+ scale_color_manual(values=c("positive"="blue3", "negative"="green3"))+ scale_x_continuous(expand = c(0, 0), limits=c(0,750)) + scale_y_continuous(expand = c(0, 0), limits=c(0,750))+ theme(legend.position="none", axis.title.x = element_text(face="italic", hjust=0.5), axis.title.y = element_text(face="italic", vjust=0.5)) +labs(title = "LTR/Gypsy", x="Kokia drynarioides", y="Gossypioides kirkii")
p4 <- ggplot(KdGk[KdGk$Lineage == "LTR/Copia", ], aes(x=kokia_, y=kirkii, shape= signK, color=signK)) + geom_point(size=2) + geom_abline(intercept=0, slope=1)+ scale_color_manual(values=c("positive"="blue3", "negative"="green3"))+ scale_x_continuous(expand = c(0, 0), limits=c(0,150)) + scale_y_continuous(expand = c(0, 0), limits=c(0,150))+ theme(legend.position="none", axis.title.x = element_text(face="italic", hjust=0.5), axis.title.y = element_text(face="italic", vjust=0.5)) +labs(title = "LTR/Copia", x="Kokia drynarioides", y="Gossypioides kirkii")

grid.arrange(p1,p2,p3,p4, ncol=2)

dev.off()


### chi2 test to determine clusters that are significantly different ###
  
chi_table <- annot_clust[,c(1,15:16)]
chi_table <- chi_table[!(rowSums(chi_table[,c(2:3)])==0),]

chi_table$p.value <- apply(chi_table[,c(2:3)], 1, function(x) chisq.test(x, simulate.p.value=TRUE, B=50000)$p.value)
chi_table$statistic <- apply(chi_table[,c(2:3)], 1, function(x) chisq.test(x)$statistic)
# note if you DO NOT use the p-value simulation
# there are 35 clusters with the warning "In chisq.test(x) : Chi-squared approximation may be incorrect"
# The "simulate.p.value=T" option (default value is FALSE) does the Monte Carlo simulation; from https://ww2.coastal.edu/kingw/statistics/R-tutorials/independ.html

chi_table$p.adjust <- p.adjust(chi_table$p.value, method="BH")

# > row.names(chi_table[(chi_table$p.adjust<0.001),])
# [1] "CL0002" "CL0005" "CL0050" "CL0084" "CL0096" "CL0097" "CL0101" "CL0105" "CL0107" "CL0110" "CL0116" "CL0117" "CL0119"
# [14] "CL0121" "CL0126" "CL0129" "CL0141" "CL0162" "CL0164" "CL0175" "CL0177" "CL0187" "CL0188"
# > row.names(chi_table[(chi_table$p.adjust<0.005),])
# [1] "CL0002" "CL0005" "CL0017" "CL0050" "CL0066" "CL0082" "CL0084" "CL0085" "CL0096" "CL0097" "CL0101" "CL0105" "CL0107"
# [14] "CL0110" "CL0116" "CL0117" "CL0119" "CL0121" "CL0126" "CL0129" "CL0136" "CL0141" "CL0158" "CL0162" "CL0164" "CL0175"
# [27] "CL0177" "CL0187" "CL0188" "CL0190" "CL0191" "CL0238" "CL0253" "CL0274"
# > row.names(chi_table[(chi_table$p.adjust<0.01),])
# [1] "CL0002" "CL0005" "CL0017" "CL0050" "CL0066" "CL0069" "CL0082" "CL0084" "CL0085" "CL0096" "CL0097" "CL0098" "CL0100"
# [14] "CL0101" "CL0105" "CL0107" "CL0110" "CL0116" "CL0117" "CL0118" "CL0119" "CL0121" "CL0126" "CL0129" "CL0136" "CL0141"
# [27] "CL0149" "CL0158" "CL0162" "CL0164" "CL0175" "CL0177" "CL0187" "CL0188" "CL0190" "CL0191" "CL0203" "CL0238" "CL0242"
# [40] "CL0253" "CL0274"p05 <- row.names(chi_table[(chi_table$p.adjust<0.05),])

p05 <- row.names(chi_table[(chi_table$p.adjust<0.05),])
# [1] "CL0002" "CL0005" "CL0010" "CL0017" "CL0050" "CL0060" "CL0066" "CL0069" "CL0082" "CL0084" "CL0085" "CL0096" "CL0097"
# [14] "CL0098" "CL0100" "CL0101" "CL0105" "CL0107" "CL0110" "CL0116" "CL0117" "CL0118" "CL0119" "CL0121" "CL0126" "CL0128"
# [27] "CL0129" "CL0131" "CL0136" "CL0141" "CL0147" "CL0149" "CL0150" "CL0153" "CL0158" "CL0161" "CL0162" "CL0164" "CL0168"
# [40] "CL0170" "CL0175" "CL0177" "CL0187" "CL0188" "CL0190" "CL0191" "CL0202" "CL0203" "CL0219" "CL0238" "CL0242" "CL0243"
# [53] "CL0253" "CL0271" "CL0274"

#  clusters with warnings if no p-value simulation
warnME <- NULL
for (i in c(1:188)) { assign("last.warning", NULL, envir = baseenv()); ifelse(has_warning(chisq.test(chi_table[i,c(2:3)])), warnME <- c(warnME,(row.names(chi_table)[i])),"") }
#  [1] "CL0001" "CL0004" "CL0013" "CL0018" "CL0019" "CL0024" "CL0031" "CL0034" "CL0037" "CL0044" "CL0054" "CL0055" "CL0064"
# [14] "CL0067" "CL0080" "CL0092" "CL0132" "CL0134" "CL0155" "CL0179" "CL0181" "CL0195" "CL0209" "CL0214" "CL0215" "CL0223"
# [27] "CL0233" "CL0234" "CL0243" "CL0246" "CL0252" "CL0259" "CL0264" "CL0267" "CL0269"

intersect(warnME, p05)
#  only CL243 overlaps between the warnings and sig diff clusters at p<0.05
                             
# make a table of only significant clusters and then learn about it 
sigtable <- chi_table[p05,c(1:3)]
sigtable$greater <- ifelse(sigtable$kokia_ > sigtable$kirkii, "kokia", ifelse(sigtable$kirkii > sigtable$kokia_, "kirkii", "same"))
table(sigtable$greater)
#kirkii  kokia 
#    21     34 

table(sigtable$Lineage)
#  *    LTR     LTR/Copia      LTR/Gypsy       non-LTR_retroposon        Unspecified 
#  2     11         6              34                  1                      1 

sum(sigtable$kirkii)*0.0095 # 0.0095 is the same factor used above, with the conversion to Mb
# [1] 70.4235
sum(sigtable$kokia_)*0.0095
# [1] 68.894
sum(sigtable$kirkii)*0.0095-sum(sigtable$kokia_)*0.0095
# [1] 1.5295
                            
                             
                             
########### relative aging of transposable elements ###########

# see TE_dating_histogram.pl for table generation

ages <- read.table("cluster.ages", header=T, row.names=1, sep="\t")
KdGkages <- read.table("KokKirk.cluster.ages", header=T, row.names=1, sep="\t")
KdGkages$age <- KdGkages$Category

#Category	age
#1	young
#3	old
#4	young
#5	old
#4*	old

KdGkages$age <- sub("3|5|4\\*","old",KdGkages$age)
KdGkages$age <- sub("1|4","young",KdGkages$age)

table(ages$age)
#  old young 
#  202    72 

table(KdGkages$age)

# old young 
# 127    49 


sigages <- merge(sigtable,ages,by="row.names")
table(sigages$age)
#  old young 
#   30    25 

sigGkKd <- merge(sigtable, KdGkages, by="row.names")
table(sigGkKd$age)

# old young 
#  31    24 

changed <- NULL
for (i in c(1:55)) { ifelse(identical(as.character(sigages$age[i]), as.character(sigGkKd$age[i])), "", changed <- c(changed,sigages$Row.names[i])) }
changed
# CL0097 young -> old

