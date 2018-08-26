library(getopt)
library(data.table)

spec = matrix(c('ifile','i',1,"character",'pfile','p',1,"character",'ofile','o',1,"character"),byrow=TRUE,ncol=4)

opt=getopt(spec)
input1 <- opt$ifile
input2 <- opt$pfile
outputfile <- opt$ofile

df1 <- read.csv(input1,sep = '\t', header = F)
df2 <- read.csv(input2,sep = '\t', header = F)

colnames(df2)[2] <- "uniq"
colnames(df1)[1] <- "uniq"
unique_positions <- setDT(df2)[uniq %chin% df1$uniq]

df1 <- df1[,c(1,3)]
merged <- merge(unique_positions, df1, by="uniq")

write.table(merged, file = outputfile,quote = F, col.names = F,row.names = F, sep = '\t')
