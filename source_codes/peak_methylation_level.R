 # getwd()
 # setwd("/home/malab5/docker/PEAC/data/")

args=commandArgs(T)
library(DiffBind)
library(exomePeak)
#bed <- import.bed("/home/malab5/docker/PEAC/data/peak_checked_merge.bed")
input.bam <- args[1]
ip.bam <- args[2]
bed <- args[3]
minimal_counts <- args[4]
out_path <- args[5]

minimal_counts <- as.numeric(minimal_counts)

sample_path <- "/home/songjie/PEAC/galaxy/tools/Methylation/sample.csv"
# ip.bam <- "wt_ip.bam"
# input.bam <- "wt_input.bam"
# bed <- "peak_checked_merge.bed"
# minimal_counts <- 10
# out_path <- 'bed_count.txt'
# Rscript peak_methylation_level.R "wt_ip.bam" "wt_input.bam" "peak_checked_merge.bed" 10 'bed_count2.txt'
### cut bed file 
bed_whole <- read.table(bed)
bed_whole <- bed_whole[,1:3]
write.table(bed_whole,file = paste0(bed,13),quote = F,sep = '\t',row.names = F,col.names = F)
bed <- paste0(bed,13)

### change the sample list file 
sample_id <- read.csv(sample_path)
sample_id$bamReads <- ip.bam
sample_id$bamControl <- input.bam
sample_id$Peaks <- bed
write.csv(sample_id,file = sample_path,row.names = F)

### load sample 
tamoxifen <- dba(sampleSheet=sample_path)
 #plot(tamoxifen)
 #bed_count <- dba.count(DBA = tamoxifen,peaks = bed)
#bed_count <- dba.count(DBA = tamoxifen,peaks = bed,minOverlap=0.01,score = DBA_SCORE_RPKM_FOLD,bRemoveDuplicates=T)
bed_count <- dba.count(DBA = tamoxifen,minOverlap=0.01,score = DBA_SCORE_RPKM_FOLD,bRemoveDuplicates=T)
#dim(bed_count$peaks[[1]])
# bed_count$peaks[[1]][1:10,]
names(bed_count$peaks[[1]])[c(4,5,6,7,8)] <- c("MFPKM_FC","MFPKM_ip","Reads_ip","MFPKM_input",  "Reads_input")
out_peak_count <- bed_count[[1]][[1]][which(bed_count[[1]][[1]]$Reads_ip>=10),]
  
#tamoxifen <- dba.count(tamoxifen, summits=250)
TOTAL_IP <- system(command = paste0("samtools view -c -F 260 ",ip.bam ),intern = T)
TOTAL_INPUT <- system(command = paste0("samtools view -c -F 260 ",input.bam ),intern = T)

##  binomial distribution based tes
ctest_res <- ctest(IP=out_peak_count$Reads_ip,
                   INPUT=out_peak_count$Reads_input,
                   TOTAL_IP=as.numeric(TOTAL_IP), 
                   TOTAL_INPUT=as.numeric(TOTAL_INPUT),
                   FOLD = 1, minimal_counts_in_fdr = 10)
out_peak_count$Reads_FC <- exp(ctest_res$log.fc)
out_peak_count$log10.p <- log10(exp(ctest_res$log.p))
out_peak_count$log10.fdr <- log10(exp(ctest_res$log.fdr))

out_peak_count[,c(4,5,7,9,10,11)] <- round(out_peak_count[,c(4,5,7,9,10,11)],3)


write.table(out_peak_count,file = out_path,row.names = F,quote = F,sep = '\t')

#test <- read.csv(("/home/malab5/R/x86_64-pc-linux-gnu-library/3.4/DiffBind/extra/tamoxifen.csv"))
