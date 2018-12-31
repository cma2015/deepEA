
library(PEA)
args=commandArgs(T)
 
 input.bam <- args[1]
 CMR <- args[2]
 refGenome <- args[3]
 outname <- args[4]
 paired <- args[5]
 level <- args[6]
 #pc_name <- arg[5]
cpus <- args[7]

paired <-  as.logical(paired)
 level <- as.numeric(level)
 cpus <- as.numeric(cpus)


cat("input bam is ",input.bam,'\n')  ## input bam
cat("Out file is ",outname,'\n')  ## input bam

cat("refGenome is ",refGenome,'\n') ## genome





results <- pseudoURatio(refGenome = refGenome,inputBAM = input.bam, cpus = cpus)


  
cmrMat <- CMRCalling(CMR = CMR, inputBAM = input.bam, 
                       mappedInput = mappedInput, refGenome = refGenome,paired=paired,level=level,cpus=cpus)


cmrMat <- cmrMat[which(cmrMat[,3]!="0"),]
#require(stringr)
#pc_name <- str_split(aaa,'\\/',simplify = T)[length(str_split(aaa,'\\/',simplify = T))]
#print(pc_name)
write.table(cmrMat, file = outname, sep = "\t", quote = F, row.names = F, col.names = F)  

