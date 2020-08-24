library(m6ALogisticModel)
library(GenomicRanges)
library(pipeR)
library(GenomicFeatures)
library(Biostrings)
library(rtracklayer)
mainDic <- "/home/galaxy/tools/3-ADVANCED-ANALYSIS/ML-based_Modelling_Analysis/"
source(paste0(mainDic, '00_feature_depend.R'))
source(paste0(mainDic, '00_extract_motif.R'))
library(argparse)

MHmakeRandomString <- function(n=1, lenght=12)
{
  randomString <- c(1:n)                  # initialize vector
  for (i in 1:n)
  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                    lenght, replace=TRUE),
                             collapse="")
  }
  return(randomString)
}


parser <- ArgumentParser()
parser$add_argument("-genome", default = NULL, dest = "genome", help = "The genome sequence in FASTA format.")
parser$add_argument("-gtf", default = NULL, dest = "gtf", help = "The annotation file in GTF format.")
parser$add_argument("-cmr", default = NULL, dest = "peak", help = "The directory of CMR regions.")
parser$add_argument("-motif", default = "RRACH", dest = "motif", help = "The motif of CMR, optional, only required for peak regions.")
parser$add_argument("-pos", default = "3", dest = "pos", help = "The relative position of CMR on the motif.")
parser$add_argument("-ratio", default = "1", dest = "ratio", help = "The ration between positive and negative samples.")
parser$add_argument("-posOut", default = NULL, dest = "posOut", help = "The output directory of positve samples in BED format.")
parser$add_argument("-negOut", default = NULL, dest = "negOut", help = "The output directory of negative samples in BED format.")


args <- parser$parse_args()
genomePath <- args$genome
if(file.exists(paste0(genomePath, ".2bit"))){
    genome <- TwoBitFile(paste0(genomePath, ".2bit"))
}else{
    dna <- readDNAStringSet(genomePath)
    dna <- replaceAmbiguities(dna, new = "N")
    export(dna, paste0(genomePath, ".2bit"))
    rm(dna)
    genome <- TwoBitFile(paste0(genomePath, ".2bit"))
}

GTF <- import(args$gtf)
exonGTF <- subset(GTF, GTF$type == "exon" & GTF$gene_biotype == "protein_coding")

# read CMR
rawPeak <- read.table(file = args$peak, sep = "\t", header = F, quote = "", stringsAsFactors = F)
colnames(rawPeak)[1:3] <- c("chr", "start", "end" )
peaks <- GRanges(rawPeak[,1:3])
strand(peaks) <- rawPeak$V6

motif <- args$motif
pos <- as.numeric(args$pos)
positive_sites <- extract_motif_list(seqRegion = peaks,
                                genome = genome,
                                motif = motif)
# retain sites located in exons
interRes <- findOverlaps(query = positive_sites, subject = exonGTF)
positive_sites <- positive_sites[unique(queryHits(interRes))]

##########################Extract negative sites####################################
interRes <- findOverlaps(query = positive_sites, subject = exonGTF)
posGeneID <- unique(exonGTF$gene_id[subjectHits(interRes)])

library(snowfall)
ratio <- as.numeric(args$ratio)
.oneNeg <- function(i, posGeneID, exonGTF, positive_sites, genome, motif, ratio){
  curGene <- posGeneID[i]
  curGTF <- exonGTF[which(exonGTF$gene_id == curGene)]
  curPosOverlap <- findOverlaps(query = positive_sites, subject = curGTF)
  curPosSites <- positive_sites[unique(queryHits(curPosOverlap))]
  curPosSitesLen <- length(curPosSites)*ratio
  
  curNegSites <- extract_motif_list(seqRegion = curGTF,
                                    genome = genome,
                                    motif = motif)
  interRes <- findOverlaps(query = curNegSites, subject = curPosSites)
  idx <- unique(queryHits(interRes))
  curNegSites <- curNegSites[setdiff(1:length(curNegSites), idx)]
  
  if(length(curNegSites) >= curPosSitesLen){
    curNegSites <- curNegSites[sample(1:length(curNegSites), curPosSitesLen)]
  }
  curList <- list(posSites = curPosSites, negSites = curNegSites)
  # cat(i, ":", curPosSitesLen, ": ", length(curNegSites), "\n")
  curList
}

Neg <- function(geneNum, posGeneID, exonGTF, positive_sites, genome, motif, CPUs = 2, ratio){
  resList <- list()
  sfInit(parallel = TRUE, cpus = CPUs)
  sfLibrary( "dplyr", character.only = TRUE)
  sfLibrary( "m6ALogisticModel", character.only = TRUE)
  sfLibrary( "GenomicRanges", character.only = TRUE)
  sfLibrary( "pipeR", character.only = TRUE)
  sfLibrary( "GenomicFeatures", character.only = TRUE)
  sfLibrary( "Biostrings", character.only = TRUE)
  sfLibrary( "rtracklayer", character.only = TRUE)
  sfSource(file = paste0(mainDic, '00_feature_depend.R'))
  sfSource(file = paste0(mainDic, '00_extract_motif.R'))
  cvRes <- sfApply( matrix(1:geneNum, ncol = 1), 1,  .oneNeg,
                    posGeneID = posGeneID, 
                    exonGTF = exonGTF, 
                    positive_sites = positive_sites,
                    genome = genome, 
                    motif = motif,
		    ratio = ratio)
  sfStop()
  cvRes
}

res <- Neg(geneNum = length(posGeneID), posGeneID = posGeneID,
           exonGTF = exonGTF,
           positive_sites = positive_sites,
           genome = genome, 
           motif = motif,
           CPUs = 2,
	   ratio = ratio)

names(res) <- posGeneID

resPos <- NULL
resNeg <- NULL
for(i in 1:length(res)){
  resPos <- c(resPos, res[[i]]$posSites)
  resNeg <- c(resNeg, res[[i]]$negSites)
}

resPos <- do.call(what = c, args = resPos)
resPos <- unique(resPos)

resNeg <- do.call(what = c, args = resNeg)
resNeg <- unique(resNeg)
export(resPos - (pos-1), con = args$posOut, format = "BED")
export(resNeg - (pos-1), con = args$negOut, format = "BED")

