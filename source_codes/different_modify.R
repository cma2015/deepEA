library('getopt')
library('DiffBind')
library('rjson')
#})
getwd()
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'infile' , 'i', 1, "character",
  'outfile' , 'o', 1, "character",
  'scorecol', 'n', 1, "integer",
  'lowerbetter', 'l', 1, "logical",
  'summits', 's', 1, "integer",
  'th', 't', 1, "double",
  'format', 'f', 1, "character",
  'plots' , 'p', 2, "character",
  'bmatrix', 'b', 0, "logical",
  "rdaOpt", "r", 0, "logical",
  'infoOpt' , 'a', 0, "logical",
  'verbose', 'v', 2, "integer",
  'help' , 'h', 0, "logical"
), byrow=TRUE, ncol=4);

opt = getopt(spec);

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}

parser <- newJSONParser()
parser$addData(opt$infile)
factorList <- parser$getObject()
filenamesIn <- unname(unlist(factorList[[1]][[2]]))



# save(parser,opt,filenamesIn,file = "diftest.RData")
# load("diftest.RData")

peaks <- filenamesIn[grepl("peaks.bed", filenamesIn)]
bams <- filenamesIn[grepl("bamreads.bam", filenamesIn)]
ctrls <- filenamesIn[grepl("bamcontrol.bam", filenamesIn)]
reps <- filenamesIn[grepl("rep", filenamesIn)]

groups <- sapply(strsplit(peaks,"-"), `[`, 1)
reps <- sapply(strsplit(peaks,"-"), `[`, 2)
samples <- sapply(strsplit(peaks,"-"), `[`, 3)

if ( length(ctrls) != 0 ) {
  sampleTable <- data.frame(SampleID=samples,
                            Condition=groups,
                            bamReads=bams,
                            bamControl=ctrls,
                            Replicate=reps,
                            Peaks=peaks,
                            Tissue=samples, # using "Tissue" column to display ids as labels in PCA plot
                            stringsAsFactors=FALSE)
} else {
  sampleTable <- data.frame(SampleID=samples,
                            Replicate=samples,
                            Condition=groups,
                            bamReads=bams,
                            Replicate=reps,
                            Peaks=peaks,
                            Tissue=samples,
                            stringsAsFactors=FALSE)
}


print(sampleTable)
write.csv(sampleTable,file="sampleTable.csv",row.names = F)

  options(warn =-1)
#sample = dba(sampleSheet=tamoxifen$samples)
sample = dba(sampleSheet="sampleTable.csv", peakFormat='bed', scoreCol=opt$scorecol, bLowerScoreBetter=opt$lowerbetter)
#if ( !is.null(opt$summits) ) {
#  sample_count = dba.count(sample, summits=opt$summits)
#} else {
  sample_count = dba.count(sample)
#}
sample_contrast = dba.contrast(sample_count, categories=DBA_CONDITION, minMembers=2)
sample_analyze = dba.analyze(sample_contrast)
diff_bind = dba.report(sample_analyze, th=opt$th)

if ( !is.null(opt$plots) ) {
  pdf(opt$plots)
  orvals = dba.plotHeatmap(sample_analyze, contrast=1, correlations=FALSE, cexCol=0.8, th=opt$th)
  dba.plotPCA(sample_analyze, contrast=1, th=opt$th, label=DBA_TISSUE, labelSize=0.3)
  dba.plotMA(sample_analyze, th=opt$th)
  dba.plotVolcano(sample_analyze, th=opt$th)
#  options(warn =-1)
  dba.plotBox(sample_analyze, th=opt$th,pvalMethod=t.test,notch=FALSE)
  dev.off()
}
resSorted <- diff_bind[order(diff_bind$FDR),]
out_file <- as.data.frame(resSorted)
write.table(out_file,file = opt$outfile,row.names = F,quote = F,sep = '\t',col.names = T)

if (opt$format == "bed") {
  out_file_bed  <- data.frame(Chrom=seqnames(resSorted),
                           Start=start(resSorted) - 1,
                           End=end(resSorted),
                           Name=paste0("Diffpeak", 1:length(resSorted)),
                           Score=-log10(out_file$FDR),
                           Strand=gsub("\\*", ".", strand(resSorted)))
  write.table(out_file_bed,file = opt$outfile,row.names = F,quote = F,sep = '\t',col.names = T)
  
} else if (opt$format == "interval") {
  # Output as interval
  df <- as.data.frame(resSorted)
  extrainfo <- NULL
  for (i in 1:nrow(df)) {
    extrainfo[i] <- paste0(c(df$width[i], df[i, 6:ncol(df)]), collapse="|")
  }
  out_file_interval  <- data.frame(Chrom=seqnames(resSorted),
                           Start=start(resSorted) - 1,
                           End=end(resSorted),
                           Name=paste0("Diffpeak", 1:length(resSorted)),
                           Score=-log10(out_file$FDR),
                           Strand=gsub("\\*", ".", strand(resSorted)),
                           Comment=extrainfo)
  write.table(out_file_interval,file = opt$outfile,row.names = F,quote = F,sep = '\t',col.names = T)
}



# Output RData file
if (!is.null(opt$rdaOpt)) {
    save.image(file = "DiffBind_analysis.RData")
}

# Output analysis info
if (!is.null(opt$infoOpt)) {
    info <- "DiffBind_analysis_info.txt"
    cat("dba.count Info\n\n", file=info, append = TRUE)
    capture.output(sample, file=info, append=TRUE)
    cat("\ndba.analyze Info\n\n", file=info, append = TRUE)
    capture.output(sample_analyze, file=info, append=TRUE)
    cat("\nSessionInfo\n\n", file=info, append = TRUE)
    capture.output(sessionInfo(), file=info, append=TRUE)
}
