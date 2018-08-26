library(getopt)
library(dplyr)
library(gsubfn)
library(ggplot2)
library(data.table)

spec = matrix(c('ifile','i',1,"character",'cov','c',1,"double"),byrow=TRUE,ncol=4)
opt=getopt(spec)
outputfile <- opt$ifile
coverage <- opt$cov

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
myFiles <- list.files(pattern="*.geneDosage.txt$")
for (i in myFiles){
all.data <- read.csv(i, header=FALSE, sep="\t",colClasses=column.types)
covfile <- paste(i,"cov",sep = '.')

covfile <- read.csv(covfile, header=FALSE, sep="\t")
strain_cov <- covfile[1,]

if(strain_cov < coverage) {
    SMEG <- 1
  } else {

binSize <- 10
myvector_all<-as.vector(as.matrix(all.data[3]))
windowAll<-slidingwindowplot(binSize,myvector_all)
df<-data.frame(windowAll[[1]],windowAll[[2]],windowAll[[3]])
colname<-c("x","mean","sd")
colnames(df)<-colname

df[,-1] <- remove_outliers(df$mean)
df[, -1] <- log2(df$mean)

P <- ggplot(data = df, aes(x = x, y = mean))
G <- P + stat_smooth(aes(outfit=fit<<-..y..),method ="loess",span=1,method.args = list(family="symmetric")) + geom_jitter()+theme_bw()+xlab("Genome Location (bp)")+ylab("Log2(% Coverage)")+scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
V<-ggplot_build(G)$data[[1]]

V[V==""] <- NA
V <- na.omit(V)

trough<- as.numeric(filter(V, y== min(y)) %>% select (ymax))
peak<-as.numeric(filter(V, y == max(y)) %>% select(ymin))

SMEG <- 2^peak/2^trough

if(SMEG < 1){
SMEG<-1.00
} else {
SMEG<-SMEG
}
}

var <- gsub(".geneDosage.txt","",i)
var1<-strapplyc(outputfile, "(.*)...", simplify = TRUE)
var <- gsub(var1,"",var)

out<-paste(var, SMEG, strain_cov, sep = '\t')
write.table(out, file = outputfile,quote = F, col.names = F,row.names = F, sep='\t',append = T)
}
