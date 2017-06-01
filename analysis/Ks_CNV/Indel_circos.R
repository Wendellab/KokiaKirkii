setwd("~/Box Sync/Research/Singleton_analysis")
library("SOFIA")
library("plyr")

circos.location <- '/Users/jconover/software/circos/circos-0.69-5'

# Basic location tables
gff <- read.table("G.raimondii.bed.SOFIA",sep="\t", col.names = c("chr", "pos", "locus"))
chromosomes <- read.table("G.raimondii.chr.SOFIA",sep="\t", col.names = c("chr", "pos", "locus"))
chromosomes$map <- rep(1,13)
indel.table <- as.data.frame(read.table("indels.circos.table", sep = "\t", col.names = c("chr", "pos", "size", "V4")))
GKI <- indel.table[c(indel.table$V4 == "GkI"), c(1:3)]
colnames(GKI)[3] <- "G. kirkii Insertion"
GKD <- indel.table[c(indel.table$V4 == "GkD"), c(1:3)]
colnames(GKD)[3] <- "G. kirkii Deletion"
KDI <- indel.table[c(indel.table$V4 == "KdI"), c(1:3)]
colnames(KDI)[3] <- "K. drynarioides Insertion"
KDD <- indel.table[c(indel.table$V4 == "KdD"), c(1:3)]
colnames(KDD)[3] <- "K. Drynarioides Deletion"

#Gene loss and gains
KDloss <- read.delim("Kokia_loss.single.txt", sep = "\n", col.names = "locus")
KDloss$`K. Drynarioides Gene Loss` <- KDloss$locus
KDgain <- read.delim("Kokia_gain.single.txt", sep = "\n", col.names = "locus")
KDgain$`K. Drynarioides Gene Gain` <- KDgain$locus
GKloss <- read.delim("Kirkii_loss.single.SOFIA.txt", sep = "\n", col.names = "locus")
GKloss$`G. Kirkii Gene Loss` <- GKloss$locus
GKgain <- read.delim("Kirkii_gain.single.txt", sep = "\n", col.names = "locus")
GKgain$`G. Kirkii Gene Gain` <- GKgain$locus



GKI <- GKI[c(1:50),]
GKD <- GKD[c(1:50),]

# Making Dataframes for SOFIA function
base <- merge(chromosomes, gff, all = TRUE)
Kirk.indel <- merge(GKI, GKD, all = TRUE)
Kirk.indel$pos <- as.integer(Kirk.indel$pos)
Kirk.CNV <- merge(GKloss, GKgain, all = TRUE)
Kirk.CNV.base <- merge(base, Kirk.CNV, all = TRUE)
Kirk <- merge(Kirk.indel, Kirk.CNV.base, all = TRUE)

#fix formatting and data types
Kirk$map <- as.integer(c(rep(1, dim(Kirk)[1])))
Kirk$locus <- factor(Kirk$locus, levels = levels(addNA(Kirk$locus)), labels = c(levels(Kirk$locus), "a"), exclude = NULL)
Kirk$chr <- as.integer(gsub("Chr", "", Kirk$chr))
Kirk$pos <- as.numeric(Kirk$pos)
Kirk$`G. kirkii Insertion` <- as.numeric(Kirk$`G. kirkii Insertion`)
Kirk$`G. kirkii Deletion` <- as.numeric(Kirk$`G. kirkii Deletion`)
DF <- data.frame(map = Kirk$map, chr = Kirk$chr, pos = Kirk$pos, locus = Kirk$locus, Insertion = Kirk$`G. kirkii Insertion`, Deletion = Kirk$`G. kirkii Deletion`, loss = Kirk$`G. Kirkii Gene Loss`, gain = Kirk$`G. Kirkii Gene Gain`)
DF <- DF[c(1:20000),]

# defining the location of each plot (5 plots total)
plotLocation<-data.frame(r0=c(0.90,.88,.78,.68),r1=c(.99,.9,.87,.77))
# in this case, the first plot is going to be in the position 0.90-0.99, the second in the
# position 0.88-0.90, and so on.

# defining the background for the 5 plots. In this case, the first two plots are not going to have
# background. All plot are going to have horizontal axes spaced at 5 units.
plotBackground<-data.frame(backgroundShow=c(FALSE,FALSE,FALSE, FALSE),
                           backgroundColor=rep('vvlgrey',4),axisShow=rep(TRUE,4),axisSep=rep(4,4))

# defining the overall configuration for the figure. In this case, there are two maps,
# the chromosomes for map 1 are going to be ordered from 1 to 12, however map 2 is going
# to have its chromosomes in order 9 to 1. In addition, the argument rev is going to
# reverse chromosomes for map 1. Finally, map 1 is going to have two different
# kind of blues (lblue_a3 and vvdblue_a3, which mieans light blue for the first one,
# and very very dark blue for the second; a3 is a degree of transparency) for their
# chromosomes, while the map two is going to be all dark green.
Kirkii.chromoConfiguration<-data.frame(order=c(1:13),map=c(rep(1,13)),
      rev=c(rep(FALSE,13)),color=c(rep('blue3',13)),radius=rep(1,13))

Kokia.chromoConfiguration<-data.frame(order=c(1:13),map=c(rep(1,13)),
      rev=c(rep(FALSE,13)),color=c(rep('green3',13)),radius=rep(1,13))

# defining plot types
plotType<-c("scatter", "scatter", "heatmap", "heatmap")

# defining marker colors. For the first plot, the letters are going to be black,
# for the second plot, the presence of numberic data is going to be colored with
# black (for traditional heatmaps with multiple colors, a color palette has to be defined),
# for the third plot, the color 'chr' tells SOFIA to use the color of the chromosome where
# the marker is located, for the forth plot, the markers are going to be colored using the
# palette piyg-11-div (higher and lower values are going to have different colors), and
# the fifth line plot is going to be colored with dark red.
plotColor<-c('PiYG', 'RdBu','black','black')

# definifn marker size. For the text plot, the marker size defines the font size.
# For scatter, it defines the circle size while for line plots, it defines the line
# thickness. For heatmaps, any random number can be included since it is not used at all.
markerSize<-c(8,10,16,16)

# please change the argument circosLocation
SOFIA(data=DF,chromoConfiguration=Kirkii.chromoConfiguration,plotBackground=plotBackground, plotLocation=plotLocation,plotType=plotType,plotColor=plotColor,markerSize=markerSize, circosLocation=circos.location,tickSuffix='cM',returnConf=TRUE,circosDisplay=FALSE)

#SOFIA(data=Kokia,linkColor='chr',linkGeometry=c(.001,.1),linkRadius=c(.57,.57),
      #linksFlag=TRUE,chromoConfiguration=Kokia.chromoConfiguration,plotBackground=plotBackground,
      #plotLocation=plotLocation,plotType=plotType,plotColor=plotColor,markerSize=markerSize,
      #circosLocation=circos.location,tickSuffix='cM',returnConf=TRUE,circosDisplay=FALSE)
