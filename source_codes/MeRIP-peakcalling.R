
library(PEA)
args=commandArgs(T)
 
 input.bam <- args[1]
 RIP.bam <- args[2]
 method <- args[3]
 GTF <- args[4]
 outname <- args[5]
 paired <- args[6]
 level <- args[7]
 ratio <- args[8]
 #pc_name <- arg[5]

paired <-  as.logical(paired)
 level <- as.numeric(level)
 ratio <-as.numeric(ratio)

cat("RIP bam is ",RIP.bam,'\n')  ## ip bam
cat("input bam is ",input.bam,'\n')  ## input bam
cat("Out file is ",outname,'\n')  ## input bam

if (method=="SlidingWindow") {
  refGenome <- args[4]
  cat("refGenome is ",refGenome,'\n') ## genome
  
  cat("Calculating mapped reads number...\n ")
  mappedInput <- system(paste0("samtools flagstat ",args[1],"  | grep 'mapped (' | cut -f 1,2,3 -d ' '"),intern = TRUE)
  mappedRIP <- system(paste0("samtools flagstat ",args[2],"  | grep 'mapped (' | cut -f 1,2,3 -d ' '"),intern = TRUE)
  
  eval(parse(text = paste0('mappedInput <- ',mappedInput)))
  eval(parse(text = paste0('mappedRIP <- ',mappedRIP)))
  
  
  cat("mappedInput is ",mappedInput,'\n') ## mappedInput
  cat("mappedRIP is ",mappedRIP,'\n') ## mappedRIP
  #system(paste0("samtools faidx ", refGenome))
  #system(paste0("cut -f1,2 ", refGenome, ".fai > genome.size"))
  
  cmrMat <- CMRCalling(CMR = "m6A", IPBAM = RIP.bam, inputBAM = input.bam, method = "SlidingWindow", 
                       mappedInput = mappedInput, mappedRIP = mappedRIP, refGenome = refGenome,paired=paired,level=level,ratio=ratio)
}

if (method=="exomePeak") {
  cmrMat <- CMRCalling(CMR = "m6A", method = "exomePeak", IPBAM = RIP.bam, inputBAM = input.bam, GTF = GTF,paired=paired,level=level,ratio=ratio)  
}
if (method=="MetPeak") {
  cmrMat <- CMRCalling(CMR = "m6A", method = "MetPeak", IPBAM = RIP.bam,  inputBAM = input.bam, GTF = GTF,paired=paired,level=level,ratio=ratio)  
}
if (method=="BayesPeak") {
  cmrMat <- CMRCalling(CMR = "m6A", method = "BayesPeak", IPBAM = RIP.bam, inputBAM = input.bam, GTF = GTF,paired=paired,level=level,ratio=ratio)  
}
if (method=="MACS2") {
  cmrMat <- CMRCalling(CMR = "m6A", method = "MACS2", IPBAM = RIP.bam, inputBAM = input.bam, GTF = GTF, ...="--nomodel",paired=paired,level=level,ratio=ratio)  
}
#require(stringr)
#pc_name <- str_split(aaa,'\\/',simplify = T)[length(str_split(aaa,'\\/',simplify = T))]
#print(pc_name)
write.table(cmrMat, file = outname, sep = "\t", quote = F, row.names = F, col.names = F)  

