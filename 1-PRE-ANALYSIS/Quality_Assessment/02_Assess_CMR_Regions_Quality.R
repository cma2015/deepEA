args <- commandArgs(trailingOnly = TRUE)

peaks <- args[1]
IP <- strsplit(args[2], ",")[[1]]
input <- strsplit(args[3], ",")[[1]]
output <- args[4]
minOverlap <- as.numeric(args[5])
method <- args[6]


library(DiffBind)
library(rtracklayer)
mainDic <- "/home/galaxy/tools/1-PRE-ANALYSIS/Quality_Assessment/"
rmarkdown::render(paste0(mainDic, "02_Assess_CMR_Regions_Quality.Rmd"),
                    output_format = 'html_document',
		    output_file = output)
