#library(RCAS, quietly = T,warn.conflicts = F)
library(GenomicFeatures, quietly = T,warn.conflicts = F)
library(stringr, quietly = T,warn.conflicts = F)
require(rtracklayer, quietly = T,warn.conflicts = F)
#library(GenomicRanges,lib.loc = "~/R/x86_64-pc-linux-gnu-library/3.5.1/")


###########  test data ##########
# peak1 <- "./Galaxy99-[m6A_peaks].encodepeak"
# 
# gff <- "./Araport11_GFF3_genes_transposons.201606.gff"
# Rscript annotation.R "./Galaxy99-[m6A_peaks].encodepeak"  "./Araport11_GFF3_genes_transposons.201606.gff" 
cat('Loading function...\n')
getTxdbFeaturesFromGRanges <- function (gffData,txdb) {
  transcripts <- GenomicFeatures::transcripts(txdb)
  # if (!is.null(gffData$transcript_id)) {
  #   m <- match(transcripts$tx_name, gffData$transcript_id)
  # }else{
  #   m <- match(transcripts$tx_name, gffData$ID)
  # }
  # m <- na.omit(m)
  # transcripts$gene_name <- gffData[m]$gene_name
  transcripts$gene_name <- str_split(transcripts@elementMetadata@listData$tx_name,"\\.",simplify = T)[,1]
  tmp <- GenomicFeatures::exonsBy(x = txdb, by = "tx", use.names = TRUE)
  exons <- BiocGenerics::unlist(tmp)
  exons$tx_name <- names(exons)
  #  exons$gene_name <- gffData[match(names(exons), gffData$transcript_id)]$gene_name
  exons$gene_name <- str_split(exons@elementMetadata@listData$tx_name,"\\.",simplify = T)[,1]
  
  tmp <- GenomicFeatures::intronsByTranscript(x = txdb, use.names = TRUE)
  introns <- BiocGenerics::unlist(tmp)
  introns$tx_name <- names(introns)
  # m <- match(names(introns), gffData$transcript_id)
  # introns$gene_name <- gffData[m]$gene_name
  introns$gene_name <- str_split(introns@elementMetadata@listData$tx_name,"\\.",simplify = T)[,1]
  
  promoters <- GenomicFeatures::promoters(txdb)
  # m <- match(promoters$tx_name, gffData$transcript_id)
  # promoters$gene_name <- gffData[m]$gene_name
  promoters$gene_name <- str_split(promoters@elementMetadata@listData$tx_name,"\\.",simplify = T)[,1]
  
  tmp <- range(GenomicFeatures::fiveUTRsByTranscript(x = txdb, 
                                                     use.names = TRUE))
  fiveUTRs <- BiocGenerics::unlist(tmp)
  fiveUTRs$tx_name <- names(fiveUTRs)
  # m <- match(names(fiveUTRs), gffData$transcript_id)
  # fiveUTRs$gene_name <- gffData[m]$gene_name
  fiveUTRs$gene_name <- str_split(fiveUTRs@elementMetadata@listData$tx_name,"\\.",simplify = T)[,1]
  
  tmp <- range(GenomicFeatures::threeUTRsByTranscript(x = txdb, 
                                                      use.names = TRUE))
  threeUTRs <- BiocGenerics::unlist(tmp)
  threeUTRs$tx_name <- names(threeUTRs)
  # m <- match(names(threeUTRs), gffData$transcript_id)
  # threeUTRs$gene_name <- gffData[m]$gene_name
  threeUTRs$gene_name <- str_split(threeUTRs@elementMetadata@listData$tx_name,"\\.",simplify = T)[,1]
  
  tmp <- GenomicFeatures::cdsBy(x = txdb, by = "tx", use.names = TRUE)
  cds <- BiocGenerics::unlist(tmp)
  cds$tx_name <- names(cds)
  # cds$gene_name <- gffData[match(names(cds), gffData$transcript_id)]$gene_name
  cds$gene_name <- str_split(cds@elementMetadata@listData$tx_name,"\\.",simplify = T)[,1]
  
  txdbFeatures = list(transcripts = transcripts, exons = exons, 
                      promoters = promoters, fiveUTRs = fiveUTRs, introns = introns, 
                      cds = cds, threeUTRs = threeUTRs)
  return(txdbFeatures)
}

args=commandArgs(T)

peak1 <- args[1]
gff <- args[2]
type <- args[3]
outfile_anno <- args[4]
outfile_feature <- args[5]


# peak1 <-  "./Galaxy99-[m6A_peaks].encodepeak"  
# gff <- "./Araport11_GFF3_genes_transposons.201606.gff"
# type <- "peak_to_gene"
# outfile_anno <-  "peak_to_gene.txt"
# outfile_feature <-  "peak_to_gene.transcript.txt"


# peak1 <- "/home/songjie/PEAC/galaxy/database/files/000/dataset_37.dat"
# gff <- '/home/songjie/PEAC/galaxy/database/files/000/dataset_59.dat.gff'
#gene_to_type /home/songjie/PEAC/galaxy/database/files/000/dataset_155.dat /home/songjie/PEAC/galaxy/database/files/000/dataset_156.dat ;


#setwd(workpath)
#print(2)
#region <- Guitar::narrowPeaktoGRanges(peak1) 
cat("Loading data ...\n")
testregion <- read.table(peak1,header=F)
testregion <- testregion[,1:3]
colnames(testregion) <- c('chrom','start','end')
region <- as(testregion, "GRanges")
#txdb_ara <- makeTxDbFromGFF(gff)

gff_region <- import(gff)
#try(txdb_ara <- makeTxDbFromGFF(gff))
############## gene peak 
cat("txdb_ara complete",'\n')
#gene_gff <- GenomicFeatures::genes(txdb_ara)
gene_gff <- gff_region[which(!is.na(gff_region@elementMetadata@listData$gene_id))]
if(length(gene_gff)==0){
gene_gff <- gff_region[which(gff_region@elementMetadata@listData$type=="gene")]
gene_gff@elementMetadata@listData$gene_id <- gene_gff@elementMetadata@listData$ID
}

#print(3)
if(type=="peak_to_gene"){
  peak_to_gene <- findOverlaps(region, gene_gff)
  region_data <- as.data.frame(region)
  region_data$strand[peak_to_gene@from] <- as.data.frame(gene_gff)$strand[peak_to_gene@to]
  region_data$gene <- NA
  region_data$gene[peak_to_gene@from] <- (as.data.frame(gene_gff))$gene_id[peak_to_gene@to]
#  print(4)
  write.table(as.matrix(region_data),file = outfile_anno,row.names = F,quote = F,sep = '\t')
}
if(type=="gene_to_type"){
  gene_to_peak <- findOverlaps(gene_gff,region)
  region_data <- gene_gff[gene_to_peak@from]@elementMetadata@listData$gene_id
  peak_name <- paste0(testregion[,1],":",testregion[,2],"-",testregion[,3])
  region_data <- cbind(unique(region_data),
                       sapply(unique(region_data), function(x) paste0(peak_name[gene_to_peak@to][which(region_data==x)],collapse = ';')))
  colnames(region_data) <- c('gene','peaks')
  write.table((region_data),file = outfile_anno,row.names = F,quote = F,sep = '\t')
}




############## transcrpt feature #####
#txdbFeatures <- getTxdbFeaturesFromGRanges(gffData = gff_region,txdb = txdb_ara)
#dt <- getTargetedGenesTable(queryRegions = region, 
#                            txdbFeatures = txdbFeatures)
#dt <- dt[order(transcripts, decreasing = TRUE)]
#dt <- as.data.frame(dt)
#
#write.table(as.matrix(dt),file = outfile_feature,row.names = F,quote = F,sep = '\t')
