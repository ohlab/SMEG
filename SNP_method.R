library(getopt)
library(dplyr)
library(gsubfn)
library(ggplot2)
library(data.table)

spec = matrix(c('ifile','i',1,"character",'cov','c',1,"double",'dnaa','d',1,"double",'min_snp','s',1,"double"),byrow=TRUE,ncol=4)
opt=getopt(spec)
outputfile <- opt$ifile
coverage <- opt$cov
dnaa <- opt$dnaa
min_snp <- opt$min_snp

theta <- c(0.02,0.04,0.06)

slidingwindowplot <- function(windowsize, inputseq)
{
  starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
  n <- length(starts)
  chunkbps <- numeric(n)
  chunkstats<- numeric(n)
  for (i in 1:n) {
    chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)] 
    chunkmean <- mean(chunk)
    chunkstdv<-sd(chunk)
    chunkbps[i] <- chunkmean
    chunkstats[i]<-chunkstdv
  }
  return (list(starts,chunkbps,chunkstats))
}
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

column.types <- c("character", "numeric", "numeric")

myFiles <- list.files(pattern="*.cov.txt")
for (i in myFiles){
  
for (y in theta){

  all.data <- read.csv(i, header=FALSE, sep="\t",colClasses=column.types)
strain_cov <- mean(all.data$V3)
cut <- nrow(all.data[which(all.data[,3]==0),])/nrow(all.data)

all.data <- all.data[which(all.data[,3]>0),] 
no_of_SNPs <- nrow(all.data)

if(no_of_SNPs < 200){
y <- y/2
}

  if(strain_cov < coverage | no_of_SNPs < min_snp) {
    SMEG <- 1 
  } else {
    
    binSize <- round(nrow(all.data)*y)
	if (binSize == 0){
	binSize <- 1
	}
    myvector_all<-as.vector(as.matrix(all.data[3]))
    windowAll<-slidingwindowplot(binSize,myvector_all)
    df<-data.frame(windowAll[[1]],windowAll[[2]],windowAll[[3]])
    colname<-c("x","mean","sd")
    colnames(df)<-colname
   
    df[,-1] <- remove_outliers(df$mean)
    df[, -1] <- lapply( df[, -1], function(x){ (x/sum(x, na.rm=TRUE))*100} )
    df[, -1] <- log2(df$mean)
    df <- df[is.finite(rowSums(df)),]
    
    P <- ggplot(data = df, aes(x = x, y = mean))
    G <- P + stat_smooth(aes(outfit=fit<<-..y..),method ="loess",span=1,method.args = list(family="symmetric")) + geom_jitter()+theme_bw()+xlab("Genome Location (bp)")+ylab("Log2(% Coverage)")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
    V<-ggplot_build(G)$data[[1]]
    V[V==""] <- NA
    V <- na.omit(V)
   
if(dnaa < 0.25 | dnaa > 0.75){

lower <- as.numeric(quantile(V$x)[2])
higher <- as.numeric(quantile(V$x)[4])

merge1<-V[which(V[,2]< lower),]
merge2<-V[which(V[,2]> higher),]
all<- rbind(merge1,merge2)

trough<- min(filter(V, x>= (quantile(V$x)[2]), x<=(quantile(V$x)[4])) %>% select (y))
peak<-filter(all, y == max(y)) %>% select(y)

} else {
   trough<-filter(V, y == min(y)) %>% select(y)
   peak<-filter(V, y == max(y)) %>% select(y)
 }
   

SMEG <- 2^peak/2^trough

   
if(SMEG < 1){
SMEG<-1.00
} else {
SMEG<-SMEG
}
  }
var <- gsub(".cov.txt","",i)

out<-paste(var, SMEG, strain_cov, no_of_SNPs, sep = '\t')
write.table(out, file = outputfile,quote = F, col.names = F,row.names = F, sep='\t',append = T)
}
}

data <- read.csv(outputfile, header = F, sep = '\t')
maxSMEG <- aggregate(V2 ~ V1, data = data, max)
minSMEG <- aggregate(V2 ~ V1, data = data, min)
data <- aggregate(.~V1,data,median)
data <- cbind(data, minSMEG[,2], maxSMEG[,2])
data <- transform(data, newcol=paste(data$`minSMEG[, 2]`, data$`maxSMEG[, 2]`, sep=" - "))[,-c(5,6)]
write.table(data, file = outputfile,quote = F, col.names = F,row.names = F, sep='\t') 
