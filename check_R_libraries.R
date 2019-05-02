list.of.packages <- c("dplyr", "getopt", "ggplot2", "gsubfn", "gplots", "ape", "dynamicTreeCut", "seqinr", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
write.table(new.packages,file = ".checkedR",quote = F,col.names = F,row.names = F)
