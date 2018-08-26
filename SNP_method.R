library(getopt)
library(dplyr)
library(gsubfn)
library(ggplot2)
library(data.table)

spec = matrix(c('ifile','i',1,"character",'cov','c',1,"double",'thet','t',1,"double"),byrow=TRUE,ncol=4)
opt=getopt(spec)
outputfile <- opt$ifile
coverage <- opt$cov
theta <- opt$thet

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

myFiles <- list.files(pattern="*.final.txt")
for (i in myFiles){
  
  all.data <- read.csv(i, header=FALSE, sep="\t",colClasses=column.types)
strain_cov <- mean(all.data$V3) 
 
  if(strain_cov < coverage) {
    SMEG <- 1 
  } else {
    
    binSize <- nrow(all.data)*theta
    myvector_all<-as.vector(as.matrix(all.data[3]))
    windowAll<-slidingwindowplot(binSize,myvector_all)
    df<-data.frame(windowAll[[1]],windowAll[[2]],windowAll[[3]])
    colname<-c("x","mean","sd")
    colnames(df)<-colname
   
cut <- nrow(all.data[which(all.data[,3]==0),])/nrow(all.data)
if( cut < 0.65){
 
    df[,-1] <- remove_outliers(df$mean)
    df[, -1] <- lapply( df[, -1], function(x){ (x/sum(x, na.rm=TRUE))*100} )
    df[, -1] <- log2(df$mean)
    
    P <- ggplot(data = df, aes(x = x, y = mean))
    G <- P + stat_smooth(aes(outfit=fit<<-..y..),method ="loess",span=1,method.args = list(family="symmetric")) + geom_jitter()+theme_bw()+xlab("Genome Location (bp)")+ylab("Log2(% Coverage)")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
    V<-ggplot_build(G)$data[[1]]
    V[V==""] <- NA
    V <- na.omit(V)
    
    trough<-filter(V, y == min(y)) %>% select(y)
    peak<-filter(V, y == max(y)) %>% select(y)

SMEG <- 2^peak/2^trough

if(SMEG < 1){
SMEG<-1.00
} else {
SMEG<-SMEG
}
}else {
SMEG<-1.00
}
  }
var <- gsub(".final.txt","",i)

out<-paste(var, SMEG, strain_cov, sep = '\t')
write.table(out, file = outputfile,quote = F, col.names = F,row.names = F, sep='\t',append = T)
}
