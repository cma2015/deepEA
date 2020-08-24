library(argparse)
library(flexdashboard)
parser <- ArgumentParser()
parser$add_argument("-cmrDir", default = NULL, dest = "cmrDir", help = "The directory of CMR regions")
parser$add_argument("-out", default = NULL, dest = "output", help = "The directory of output")
parser$add_argument("-method", default = NULL, dest = "method", help = "Normalization method")
parser$add_argument("-k", default = NULL, dest = "k", help = "The number of clusters in kmeans")
parser$add_argument("-log", default = NULL, dest = "log", help = "Log")

args <- parser$parse_args()
mainDic <- "/home/galaxy/tools/3-ADVANCED-ANALYSIS/Multi-omics_Integrative_Analysis/"

method <- args$method
cmrDir <- args$cmrDir
k <- args$k
log <- args$log
rmarkdown::render(paste0(mainDic, "01_integrate_two_data_sets.Rmd"), output_file = args$out)
