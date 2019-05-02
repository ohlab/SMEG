library(getopt)
library(data.table)
spec = matrix(c('ifile','i',1,"character",'ifile2','x',1,"character",'ifile3','y',1,"character"),byrow=TRUE,ncol=4)
opt=getopt(spec)

inputfile1 <- opt$ifile
inputfile2 <- opt$ifile2
inputfile3 <- opt$ifile3

uniqSNP <- as.data.table(read.csv(inputfile1, header = T, sep = '\t',stringsAsFactors=FALSE))
ref <- read.csv(inputfile2, header = F, sep = '\t',stringsAsFactors=FALSE)
des <- as.data.table(read.csv(inputfile3, header = F, sep = '\t',stringsAsFactors=FALSE))
des$V1 <- des$V1 + 1

uniqSNP$ID <- "SMEG"
des$ID <- "SMEG"

setkey(des, ID,V1, V2)
uniqSNP[, c("V1", "V2") := position]

uniqSNP <- foverlaps(uniqSNP, des)[, c("strain", "position", "basetype","V1","V2","V3"), with = FALSE]
colName <- c("strain", "position", "basetype","start_des","stop_des","gene")
colnames(uniqSNP) <- colName
################
colName <- c("reference", "start_ref", "stop_ref","gene")
colnames(ref) <- colName
uniqSNP <- merge(uniqSNP, ref, by="gene")[, c("strain", "position", "basetype","start_des","stop_des","start_ref", "stop_ref","gene"), with = FALSE]
uniqSNP$distance_form_start <- uniqSNP$position - uniqSNP$start_des
uniqSNP$check_strand <- uniqSNP$stop_ref - uniqSNP$start_ref

positive_Strand <- uniqSNP[which(uniqSNP[,10]>0),]
negative_Strand <- uniqSNP[which(uniqSNP[,10]<0),]

positive_Strand$new_position <- positive_Strand$start_ref + positive_Strand$distance_form_start
positive_Strand <- positive_Strand[, c("strain", "new_position", "basetype"), with = FALSE]

negative_Strand$new_position <- negative_Strand$start_ref - negative_Strand$distance_form_start
negative_Strand <- negative_Strand[, c("strain", "new_position", "basetype"), with = FALSE]
#########
negative_Strand$basetype <- gsub('A', 'w', negative_Strand$basetype) 
negative_Strand$basetype <- gsub('a', 'w', negative_Strand$basetype) 
negative_Strand$basetype <- gsub('T', 'x', negative_Strand$basetype) 
negative_Strand$basetype <- gsub('t', 'x', negative_Strand$basetype) 
negative_Strand$basetype <- gsub('C', 'y', negative_Strand$basetype) 
negative_Strand$basetype <- gsub('c', 'y', negative_Strand$basetype) 
negative_Strand$basetype <- gsub('G', 'z', negative_Strand$basetype) 
negative_Strand$basetype <- gsub('g', 'z', negative_Strand$basetype) 
###################
negative_Strand$basetype <- gsub('w', 'T', negative_Strand$basetype) 
negative_Strand$basetype <- gsub('x', 'A', negative_Strand$basetype) 
negative_Strand$basetype <- gsub('y', 'G', negative_Strand$basetype) 
negative_Strand$basetype <- gsub('z', 'C', negative_Strand$basetype)
###################
out <- rbind(positive_Strand,negative_Strand)
out <- out[order(out$new_position),]
write.table(out, file = "newcoordinates.txt",quote = F, col.names = F,row.names = F, sep='\t')
####################

