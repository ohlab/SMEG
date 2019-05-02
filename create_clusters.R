library(getopt)
library(ape)
library(dynamicTreeCut)
library(data.table)

spec = matrix(c('ifile','i',1,"character",'cfile','c',1,"double"),byrow=TRUE,ncol=4)
opt=getopt(spec)
input <- opt$ifile
condition <- opt$cfile

removeStrains <- function(data,outlier,outlier2) {
  if(identical(outlier2, "character(0)")){
    df <- data.frame(t(data))
  } else {
    outlier <- as.character(t(outlier[,1]))
    data <- as.data.frame(data)
    data <- data[ , -which(names(data) %in% c(outlier))]
    data <- data[!rownames(data) %in% outlier, ]
    df <- data.frame(t(data))
  }
  df <- df[complete.cases(df), ]
  return(df)
}

empty <- (file.info("breakCluster.txt")$size )
empty[is.na(empty)] <- 0

tree <- read.tree("tree.newick")
data<-cophenetic(tree)
outlier <- as.data.frame(data[,1])
setDT(outlier, keep.rownames = TRUE)[]
if(condition == 1){
  remove <- 100000000
}else{
remove <- (median(outlier$`data[, 1]`)*30)
}
outlier <- outlier[which(outlier[,2]>remove),]
outlier2 <- as.character(outlier[,1])

df <- removeStrains(data,outlier,outlier2)
hc <- hclust(dist(df), "mcquitty")
dend <- as.dendrogram(hc)

########################
clusters <- cutreeDynamic(hc, distM = as.matrix(dist(df)), method = "hybrid",minClusterSize = 2,deepSplit = 4)

a <- as.data.frame(clusters)
b <- as.data.frame(row.names(df))
out <- cbind(b,a)
colnames(out)[1] <- "strains"

if(empty > 0)
{
  subTree <- read.csv("breakCluster.txt", header = F, sep = '\t')
  out <- read.csv(input, header = F, sep = '\t')
  colname <- c("strains","clusters")
  colnames(out) <- colname
  for (v in 1:nrow(subTree))
  {
    r <- subTree[v,]
    clus_num <- paste("^",r,"$",sep = '')
    otherClusters <- out[- grep(clus_num, out$clusters),]
    otherClusters2 <- as.character(out[- grep(clus_num, out$clusters),1]) 
    df <- removeStrains(data,otherClusters,otherClusters2)

    if(nrow(df) >= 2){

    hc <- hclust(dist(df), "mcquitty")
    dend <- as.dendrogram(hc)
    clusters <- cutreeDynamic(hc, distM = as.matrix(dist(df)), method = "hybrid",minClusterSize = 2,deepSplit = 0)
    a <- as.data.frame(clusters)
    b <- as.data.frame(row.names(df))
    results <- cbind(b,a) 
    tmp <- paste(r,".","%d",sep = '')
    results <- transform(results, clusters = sprintf(tmp, clusters))
    colnames(results)[1] <- "strains"
    colnames(otherClusters)[1] <- "strains"
    out <- rbind(otherClusters,results)
    }else{
      out <- out
    }

  }
  out <- setDT(out)[, temp := .GRP, by = clusters][,c(1,3)]
  colnames(out)[2] <- "clusters"
}

write.table(out, file = input,quote = F, col.names = F,row.names = F, sep='\t')
