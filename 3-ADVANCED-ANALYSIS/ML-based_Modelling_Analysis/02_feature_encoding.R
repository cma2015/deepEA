suppressPackageStartupMessages(library(magrittr)) # need to run every time you start R and want to use %>%
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(m6ALogisticModel))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(pipeR))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(rtracklayer))
library(argparse)
mainDic <- "/home/galaxy/tools/3-ADVANCED-ANALYSIS/ML-based_Modelling_Analysis/"
source(paste0(mainDic, '00_feature_depend.R'))
source(paste0(mainDic, '00_extract_motif.R'))
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
parser$add_argument("-cmrType", default = NULL, dest = "cmrType", help = "The CMR type: peaks or sites.")
parser$add_argument("-cmr", default = NULL, dest = "peak", help = "The directory of CMR regions.")
parser$add_argument("-type", default = "genomic", dest = "featureType", help = "The feature type: genomic, sequence or both.")
parser$add_argument("-motif", default = "RRACH", dest = "motif", help = "The motif of CMR, optional, only required for peak regions.")
parser$add_argument("-gtf", default = NULL, dest = "gtf", help = "The annotation file in GTF format.")
parser$add_argument("-genome", default = NULL, dest = "genome", help = "The genome sequence in FASTA format.")
parser$add_argument("-pos", default = "3", dest = "pos", help = "The relative position of CMR on the motif.")
parser$add_argument("-seqFeatures", default = NULL, dest = "seqFeatures", help = "The selected sequence features.")
parser$add_argument("-genomicFeatures", default = NULL, dest = "genomicFeatures", help = "The selected genomic features.")
parser$add_argument("-out", default = NULL, dest = "out", help = "A string specifying the directory of output file.")
parser$add_argument("-flank", default = 50, type = "integer", dest = "flank", help = "The sequence length of upstream and downstream centered on CMR.")
parser$add_argument("-scale", default = "F", dest = "scale", help = "Logical, select if standardization will be performed.")

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

# read CMR
rawPeak <- read.table(file = args$peak, sep = "\t", header = F, quote = "", stringsAsFactors = F)
colnames(rawPeak)[1:3] <- c("chr", "start", "end" )
peaks <- GRanges(rawPeak[,1:3])
strand(peaks) <- rawPeak$V6
flank <- as.numeric(args$flank)
if(args$cmrType == "peaks"){
        motif <- args$motif
        pos <- as.numeric(args$pos)
        flank <- flank - round(nchar(motif)/(pos - 1))
        positive_sites <- extract_motif_list(seqRegion = peaks,
                                        genome = genome,
                                        motif = motif)
    }else{
        positive_sites <- peaks
}
if(args$featureType == "genomic" | args$featureType == "both"){
    if(is.null(args$gtf)){
        stop("Genomic features are selected, please provide annotation file (GTF format) file.")
    }
    GTF <- args$gtf
    gffTxdb <- suppressMessages(suppressWarnings(makeTxDbFromGFF(GTF)))
    matureSE <- SummarizedExperiment()
    rowRanges(matureSE) <- positive_sites

    ##### genomic feature
    genomicFeatures <-suppressWarnings(predictors_annot(se = matureSE,
                                        txdb = gffTxdb,
                                        bsgnm = genome,
                                        genes_ambiguity_method = "average",
                                        self_clustering = T,
                                        standardization = args$scale,
                                        motif_clustering = "RRACH"))
    genomicFeatures <- mcols(genomicFeatures)
    genomicFeatures <- sapply(genomicFeatures, as.numeric)

    geFeatures <- unlist(strsplit(args$genomicFeatures, ","))
    geFeatures <- geFeatures[which(geFeatures != "None")]
    genomicFeatures <- genomicFeatures[,geFeatures]

    if(args$featureType == "both"){
        Seq <- Biostrings::getSeq(genome, positive_sites + flank)
        names(Seq) <- paste0("seq_", 1:length(Seq))
        idx <- setdiff(1:length(Seq), grep("N", Seq))
        Seq <- Seq[idx]
        outSeq = MHmakeRandomString()
        writeXStringSet(Seq, filepath = paste0(mainDic, outSeq, ".fasta"))
        cmd <- paste0("sed -i 's/T/U/g' ", mainDic, outSeq, ".fasta")
        system(cmd)
        seqFeatureName <- unlist(strsplit(args$seqFeatures, ","))
        seqFeatureName <-seqFeatureName[which(seqFeatureName != "None")]
        positive_sites <- positive_sites[idx]
        if(is.element("one_hot_SCP", seqFeatureName)){
            seqFeatures <- sequence_features(positive_sites + flank, genome)
            seqFeatureName <- seqFeatureName[which(seqFeatureName != "one_hot_SCP")]
        }else{
            seqFeatures <- NULL
        }

        for(i in 1:length(seqFeatureName)){
            out = MHmakeRandomString()
            curFeatures <- runBioSeq(method = seqFeatureName[i], inputSeq = paste0(mainDic, outSeq, ".fasta"), out = paste0(mainDic, out, ".txt"))
            seqFeatures <- cbind(seqFeatures, curFeatures)
        }
        resFeatures <- cbind(genomicFeatures, seqFeatures)
    }else{
        resFeatures <- genomicFeatures
    }
}else{

    Seq <- Biostrings::getSeq(genome, positive_sites + flank)
    names(Seq) <- paste0("seq_", 1:length(Seq))
    idx <- setdiff(1:length(Seq), grep("N", Seq))
    Seq <- Seq[idx]
    outSeq = MHmakeRandomString()
    writeXStringSet(Seq, filepath = paste0(mainDic, outSeq, ".fasta"))
    cmd <- paste0("sed -i 's/T/U/g' ", mainDic, outSeq, ".fasta")
    system(cmd)
    seqFeatureName <- unlist(strsplit(args$seqFeatures, ","))
    seqFeatureName <- seqFeatureName[which(seqFeatureName != "None")]
    positive_sites <- positive_sites[idx]
    if(is.element("one_hot_SCP", seqFeatureName)){
        seqFeatures <- sequence_features(positive_sites + flank, genome)
        seqFeatureName <- seqFeatureName[which(seqFeatureName != "one_hot_SCP")]
    }else{
        seqFeatures <- NULL
    }

    for(i in 1:length(seqFeatureName)){
        out = MHmakeRandomString()
        curFeatures <- runBioSeq(method = seqFeatureName[i], inputSeq = paste0(mainDic, outSeq, ".fasta"), out = paste0(mainDic, out, ".txt"))
        seqFeatures <- cbind(seqFeatures, curFeatures)
    }
    resFeatures <- seqFeatures
}

write.table(resFeatures, file = args$out, sep = "\t", quote = F)

