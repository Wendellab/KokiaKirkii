
library(tidyr)
library(plyr)

df <- read.table("indels.table", header=T)
df$Gk.NULL <- gsub("\\/.*","",df$Gk.NULL)
df$Kokia.NULL <- gsub("\\/.*","",df$Kokia.NULL)
df <- df[df$Gk.NULL!=df$Kokia.NULL,]
df <- separate(df, 'ALT', paste("ALT", 1:4, sep="."), sep =",", extra="drop")
df[ df == "*" ] <- NA
df[is.na(df)] <-""
df$refLen<- nchar(as.vector(df$REF))
df$altLen.1<- nchar(as.vector(df$ALT.1))
df$altLen.2<- nchar(as.vector(df$ALT.2))
df$altLen.3<- nchar(as.vector(df$ALT.3))
df$altLen.4<- nchar(as.vector(df$ALT.4))


df$Gk.len <- { ifelse (df$Gk.NULL==0, df$refLen, 
	ifelse(df$Gk.NULL==1, df$altLen.1, 
		ifelse(df$Gk.NULL==2, df$altLen.2,
			ifelse(df$Gk.NULL==3, df$altLen.3, 
				ifelse(df$Gk.NULL==4, df$altLen.4, "err")))))}
				
df$Kokia.len <- { ifelse (df$Kokia.NULL==0, df$refLen, 
	ifelse(df$Kokia.NULL==1, df$altLen.1, 
		ifelse(df$Kokia.NULL==2, df$altLen.2,
			ifelse(df$Kokia.NULL==3, df$altLen.3, 
				ifelse(df$Kokia.NULL==4, df$altLen.4, "err")))))}
				
df$diffK.G <- as.numeric(df$Kokia.len) - as.numeric(df$Gk.len)
df <- df[df$Gk.len!=df$Kokia.len,]

###
histGk <- hist(as.numeric(df$Gk.len), plot=F, breaks=77)
histKd <- hist(as.numeric(df$Kokia.len), plot=F, breaks=133)

#png("Figure_indels.png", 5000, 5000, pointsize=12, res=600)
par(mfrow=c(2,2))

plot(histGk, col=rgb(0,0,1,2/4), main="A", xlab="Indel Size", xaxt='n')
plot(histKd, col=rgb(0,1,0,1/4), add=T)
axis(1, at = seq(0,400, by=20), las=2)

plot(histGk, col=rgb(0,0,1,2/4), main="B", xlab="Indel Size", xlim=c(10,100), ylim=c(10,50000), xaxt='n')
plot(histKd, col=rgb(0,1,0,1/4), add=T)
axis(1, at = seq(0,400, by=20), las=2)

plot(histGk, col=rgb(0,0,1,2/4), main="C", xlab="Indel Size", xlim=c(50,250), ylim=c(10,1750), xaxt='n')
plot(histKd, col=rgb(0,1,0,1/4), add=T)
axis(1, at = seq(0,400, by=20), las=2)

plot(histGk, col=rgb(0,0,1,2/4), main="D", xlab="Indel Size", xlim=c(100,400), ylim=c(10,200), xaxt='n')
plot(histKd, col=rgb(0,1,0,1/4), add=T)
axis(1, at = seq(0,400, by=20), las=2)

#dev.off()

####


smalldf <- df[(as.numeric(df$Kokia.len) < 51),]
smalldf <- smalldf[(as.numeric(smalldf$Gk.len) < 51),]

ShistGk <- hist(as.numeric(smalldf$Gk.len), plot=F, breaks=50)
ShistKd <- hist(as.numeric(smalldf$Kokia.len), plot=F, breaks=50)
plot(ShistGk, col=rgb(0,0,1,1/4), main="Indels <50 nt in K. drynarioides and G. kirkii", xlab="Indel Size")
plot(ShistKd, col=rgb(1,0,0,1/4), add=T)

####

tinydf <- df[(as.numeric(df$Kokia.len) < 11),]
tinydf <- tinydf[(as.numeric(tinydf$Gk.len) < 11),]

ThistGk <- hist(as.numeric(tinydf$Gk.len), plot=F, breaks=10)
ThistKd <- hist(as.numeric(tinydf$Kokia.len), plot=F, breaks=10)
plot(ThistGk, col=rgb(0,0,1,1/4), main="Indels <10 nt in K. drynarioides and G. kirkii", xlab="Indel Size")
plot(ThistKd, col=rgb(1,0,0,1/4), add=T)

#####

png("Figure_largeindels.png", 5000, 5000, pointsize=12, res=600)

largedf <- df[(as.numeric(df$Kokia.len) > 100 | as.numeric(df$Gk.len) > 100),]

LhistGk <- hist(as.numeric(largedf$Gk.len), plot=F, breaks=77)
LhistKd <- hist(as.numeric(largedf$Kokia.len), plot=F, breaks=133)
plot(LhistGk, col=rgb(0,0,1,2/4), main="Indels >100 nt in K. drynarioides and G. kirkii", xlab="Indel Size")
plot(LhistKd, col=rgb(0,1,0,1/4), add=T)

dev.off()

> sum(df$diffK.G)
[1] 197562
> sum(smalldf$diffK.G)
[1] -2169
> sum(tinydf$diffK.G)
[1] -5615

###

# peg as insertions or deletions
# first only take rows that are the same as the reference for either Kokia or Gk
indels <- df[(as.numeric(df$Kokia.NULL) == 0 | as.numeric(df$Gk.NULL) == 0),]

> nrow(indels[as.numeric(indels$Kokia.len) > as.numeric(indels$refLen),])
[1] 130177
> nrow(indels[as.numeric(indels$Kokia.len) < as.numeric(indels$refLen),])
[1] 159222
> nrow(indels[as.numeric(indels$Gk.len) > as.numeric(indels$refLen),])
[1] 87951
> nrow(indels[as.numeric(indels$Gk.len) < as.numeric(indels$refLen),])
[1] 113241

Kdins <- (indels[as.numeric(indels$Kokia.len) > as.numeric(indels$refLen),])
Kddel <- (indels[as.numeric(indels$Kokia.len) < as.numeric(indels$refLen),])
Gkins <- (indels[as.numeric(indels$Gk.len) > as.numeric(indels$refLen),])
Gkdel <- (indels[as.numeric(indels$Gk.len) < as.numeric(indels$refLen),])


> sum(Kdins$diffK.G)
[1] 838870
> sum(Kddel$diffK.G)
[1] -770267
> sum(Gkins$diffK.G)
[1] -424764
> sum(Gkdel$diffK.G)
[1] 537999

> abs(sum(Kdins$diffK.G))-abs(sum(Kddel$diffK.G))
[1] 68603
> abs(sum(Gkins$diffK.G))-abs(sum(Gkdel$diffK.G))
[1] -113235

histKdI <- hist(as.numeric(Kdins$Kokia.len), plot=F, breaks=532)
histGkI <- hist(as.numeric(Gkins$Gk.len), plot=F, breaks=308)

histKdD <- hist(abs(as.numeric(Kddel$diffK.G)), plot=F, breaks=221)
histGkD <- hist(abs(as.numeric(Gkdel$diffK.G)), plot=F, breaks=242)

png("Figure_indels.png", 5000, 5000, pointsize=12, res=600)
split.screen( figs = c( 2, 1 ) )

screen(1)
plot(histGkI$count, col=rgb(0,0,1,2/4), main="Insertion sizes", xlab="", xaxt='n', log="y", type='h', lwd=2, lend=2, ylab='Frequency') #blue
par(new=TRUE)
plot(histKdI$count, col=rgb(0,1,0,2/4), log="y", type='h', lwd=2, lend=2, xaxt='n', yaxt='n', ylab="", xlab="") #green
axis(1, at = seq(0,540, by=20), las=2)

screen(2)
plot(histGkD$count, col=rgb(0,0,1,2/4), main="Deletion sizes", xlab="", xaxt='n', log="y", type='h', lwd=2, lend=2, ylab='Frequency') #blue
par(new=TRUE)
plot(histKdD$count, col=rgb(0,1,0,2/4), log="y", type='h', lwd=2, lend=2, xaxt='n', yaxt='n', ylab="", xlab="") #green
axis(1, at = seq(0,400, by=10), las=2)


close.screen( all = TRUE )

dev.off()

### table for circos plot ###



KdI <- Kdins[,c("CHROM", "POS", "Kokia.len")]
KdI$type <- "KdI"
colnames(KdI)[3] <- "length"

GkI <- Gkins[,c("CHROM", "POS", "Gk.len")]
GkI$type <- "GkI"
colnames(GkI)[3] <- "length"

KdD <- Kddel[,c("CHROM", "POS", "diffK.G")]
KdD$type <- "KdD"
colnames(KdD)[3] <- "length"
KdD$length <- abs(KdD$length)

GkD <- Gkdel[,c("CHROM", "POS", "diffK.G")]
GkD$type <- "GkD"
colnames(GkD)[3] <- "length"

indel_table <- rbind(KdI, GkI, KdD, GkD)
write.table(indel_table, file="indels.circos.table", quote=FALSE, sep = "\t", row.names=FALSE)


