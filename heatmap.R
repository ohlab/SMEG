library(gplots)
data <- read.csv("merged_table.txt", header=T, sep="\t")
df = setNames(data.frame(t(data[,-1])), data[,1]) #transform the data. Here, the rowname is reomved
df <- cbind(Genome = rownames(df), df) #duplicate the 1st row and add a row name called "Sample"
rownames(df) <- NULL #remove the duplicate
row.names(df) <- df$Genome
df <- df[, -1]
df[is.na(df)] <- 0.99
scaleyellowred <- colorRampPalette(c("black","lightyellow","red"), space = "rgb")(n=299)
my.breaks <- c(seq(0.99, 0.999, length=100),seq(0.99901, max(1,((max(df)+1)/2)), length=100),seq((max(1,((max(df)+1)/2))+0.0001), max(df), length=100))
lmat = rbind(c(0,3),c(2,1),c(0,4))
lwid = c(0.5,9)
lhei = c(4,9,2.7)

pdf("merged.pdf")
heatmap.2(as.matrix(df),breaks = my.breaks,srtCol=90,denscol=NA,key.ylab=NA,dendrogram="column",key.title=NA,key.xlab="SMEG",cexCol=1,lmat = lmat, lwid = lwid, lhei = lhei,trace="none",col = scaleyellowred,labRow = "",key.ytickfun = NA)
dev.off()
