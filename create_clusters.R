library(WGCNA)
library(getopt)
library(ape)
library(dynamicTreeCut)
library(data.table)

spec = matrix(c('ifile','i',1,"character"),byrow=TRUE,ncol=4)
opt=getopt(spec)
input <- opt$ifile

tree <- read.tree(input)
data<-cophenetic(tree)
outlier <- as.data.frame(data[,1])
setDT(outlier, keep.rownames = TRUE)[]
remove <- (median(outlier$`data[, 1]`)*30)
outlier <- outlier[which(outlier[,2]>remove),]
outlier2 <- as.character(outlier[,1])

if(identical(outlier2, "character(0)")){
df = setNames(data.frame(t(data[,-1])), data[,1]) #transform the data. Here, the rowname is removed
} else {
outlier <- as.character(t(outlier[,1]))
data <- as.data.frame(data)
data <- data[ , -which(names(data) %in% c(outlier))]
data <- data[!rownames(data) %in% outlier, ]
df = setNames(data.frame(t(data[,-1])), data[,1]) #transform the data. Here, the rowname is removed
}
df <- df[complete.cases(df), ]
hc <- hclust(dist(df), "mcquitty")
dend <- as.dendrogram(hc)

########################
clusters <- cutreeDynamic(hc, distM = as.matrix(dist(df)), method = "hybrid",minClusterSize = 2,deepSplit = 4)
colorTree = labels2colors(cutreeDynamic(hc, distM = as.matrix(dist(df)), method = "hybrid",minClusterSize = 2,deepSplit = 4))
a <- as.data.frame(clusters)
b <- as.data.frame(row.names(df))
out <- cbind(b,a)

pdf("clusters_deepSplit4.pdf")
plotDendroAndColors(hc,colors=colorTree, dendroLabels = FALSE, marAll = c(1, 8, 3, 1),main = "")
dev.off()

write.table(out, file = "clusters_deepSplit4.txt",quote = F, col.names = F,row.names = F, sep='\t')
###########################
clusters <- cutreeDynamic(hc, distM = as.matrix(dist(df)), method = "hybrid",minClusterSize = 2,deepSplit = 0)
colorTree = labels2colors(cutreeDynamic(hc, distM = as.matrix(dist(df)), method = "hybrid",minClusterSize = 2,deepSplit = 0))
a <- as.data.frame(clusters)
b <- as.data.frame(row.names(df))
out <- cbind(b,a)

pdf("clusters_deepSplit0.pdf")
plotDendroAndColors(hc,colors=colorTree, dendroLabels = FALSE, marAll = c(1, 8, 3, 1),main = "")
dev.off()

write.table(out, file = "clusters_deepSplit0.txt",quote = F, col.names = F,row.names = F, sep='\t')

