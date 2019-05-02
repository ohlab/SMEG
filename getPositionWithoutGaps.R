library(seqinr)
library(getopt)
library(data.table)
spec = matrix(c('ifile','i',1,"character",'ifile2','x',1,"character",'mfile','m',1,"double"),byrow=TRUE,ncol=4)
opt=getopt(spec)

inputfile1 <- opt$ifile
inputfile2 <- opt$ifile2
method <- opt$mfile

data <- as.data.table(read.fasta(file = inputfile2))
data <- as.data.frame(data)
colNam <- colnames(data)
data$position <- seq.int(nrow(data))
delete <- as.numeric(which(data[1] == "-"))

if(length(delete) != 0){
  new <- as.data.frame(data[-(delete),])
  new$pos_no_gap <- seq.int(nrow(new))
}else{
  new <- as.data.frame(data)
  new$pos_no_gap <- data$position
}

### CLUSTER SNP METHOD ###
if(method == 0){
sub <- read.csv(inputfile1, sep = '\t', header = T)
colnames(sub)[1] <- "position"

newtable <- merge(sub,new, by  = "position") 
newtable <- newtable[,c(3,6,2)]
colname <- c("strain","position","basetype")
colnames(newtable) <- colname
newtable <- transform(newtable, strain = sprintf('cluster%d', strain))
}else{
### REFERENCE METHOD #####
sub <- read.csv(inputfile1, sep = '\t', header = T)

newtable <- merge(sub,new, by  = "position")

myvars <- c("strain", "pos_no_gap","basetype")
newtable <- newtable[myvars]
colname <- c("strain","position","basetype")
colnames(newtable) <- colname

}

write.table(newtable, file = "modified_uniq_cluster_SNPs.txt",quote = F, col.names = T,row.names = F, sep='\t')
