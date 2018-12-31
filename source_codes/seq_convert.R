library(seqinr)
library(stringr)
args=commandArgs(T)

#input.seq <- args[1]
input.seq <- args[1]
output.seq <- args[2]
cat(input.seq)
sequence <- read.fasta(input.seq,as.string = T)

sequence <- lapply(sequence,toupper)
sequence_table <- paste0(names(sequence),',',unlist(sequence))

write.table(x = sequence_table,file = output.seq,row.names = F,quote = F,col.names = F)
