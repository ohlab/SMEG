library(gplots)
data <- read.csv("merged_table.txt", header=T, sep="\t")
df = setNames(data.frame(t(data[,-1])), data[,1]) #transform the data. Here, the rowname is removed
df <- cbind(Genome = rownames(df), df) #duplicate the 1st row and add a row name called "Genome"
rownames(df) <- NULL #remove the duplicate
row.names(df) <- df$Genome
df <- df[, -1]
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)
lmat = rbind(c(0,3),c(2,1),c(0,4))
lwid = c(0.5,9)
lhei = c(4,9,2.7)

pdf("merged.pdf")
heatmap.2(as.matrix(df),margin=c(23,0.4),denscol=NA,key.ylab=NA,dendrogram="column",key.title=NA,key.xlab="SMEG",cexCol=1,lmat = lmat, lwid = lwid, lhei = lhei,trace="none",col = scaleyellowred,labRow = "",key.ytickfun = NA)
dev.off()
