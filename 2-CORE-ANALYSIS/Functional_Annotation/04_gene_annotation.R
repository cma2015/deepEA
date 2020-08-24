library(data.table)
library(rtracklayer)
library(argparse)

options(scipen = 20, stringsAsFactors = F)


parser <- ArgumentParser()
parser$add_argument("-cmr", default = NULL, dest = "peak",
                    help = "The directory of CMR regions.")
parser$add_argument("-gtf", default = NULL, dest = "GTF",
                    help = "The directory of GTF.")
parser$add_argument("-overlap", default = 50, dest = "overlap")
parser$add_argument("-outGene", default = NULL, dest = "outGene",
                    help = "")
parser$add_argument("-outBED", default = NULL, dest = "outBED",
                    help = "")


args <- parser$parse_args()

peakMat <- fread(file = args$peak, sep = '\t', header = F, quote = '', stringsAsFactors = F)
peakMat <- data.frame(peakMat[,1:3])
colnames(peakMat) <- c("chr", "start", "end" )
peakGRange <- GenomicRanges::GRanges(peakMat)

GTF <- import(args$GTF)
geneGTF <- subset(GTF, GTF$type == "gene")
rm(GTF)

interRes <- findOverlaps(query = peakGRange, subject = geneGTF, minoverlap = as.numeric(args$overlap))
uniqueHits <- unique(queryHits(interRes))

resDF <- data.frame()
for (i in 1:length(interRes)) {
  curPeak <- as.data.frame(peakGRange[queryHits(interRes)[i]])
  curGene <- geneGTF[subjectHits(interRes)[i]]$gene_id
  curStrand <- as.character(strand(geneGTF[subjectHits(interRes)[i]]))
  res <- c(curPeak$seqnames, curPeak$start, curPeak$end, curGene, ".", curStrand)
  resDF <- rbind(resDF, res)
}

colnames(resDF) <- c("chr", "start", "end", "name", "N", "strand")
geneID <- unique(resDF$name)
write.table(geneID, file = args$outGene, sep = "\t", quote = F, row.names = F, col.names = F)
write.table(resDF, file = args$outBED, sep = "\t", quote = F, row.names = F, col.names = F)




