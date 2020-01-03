args <- commandArgs(trailingOnly = TRUE)

peaks <- args[1]
IP <- strsplit(args[2], ",")[[1]]
input <- strsplit(args[3], ",")[[1]]
output <- args[4]
minOverlap <- as.numeric(args[5])
method <- args[6]


library(DiffBind)
library(rtracklayer)

genePeak <- rtracklayer::import(args[1])
sampleSheet <- data.frame(SampleID = paste0("sample_", 1:length(IP)),
                          Tissue = "Any",
                          Factor = "m6A",
                          Condition = "Normal",
                          Replicate = 1:length(IP),
                          bamReads = IP,
                          ControlID = paste0("control_", 1:length(IP)),
                          bamControl = input,
                          Peaks = peaks,
                          PeakCaller = "bed",
                          stringsAsFactors = FALSE)
tamoxifen <- dba(sampleSheet = sampleSheet, minOverlap = minOverlap)
tamoxifen_count <- dba.count(tamoxifen, score = method, peaks = genePeak,
                       minOverlap = minOverlap)
      
peaks <- tamoxifen_count$peaks[[1]]

write.table(peaks, file=args[4], sep = "\t", quote = F, row.names = F, col.names = T)