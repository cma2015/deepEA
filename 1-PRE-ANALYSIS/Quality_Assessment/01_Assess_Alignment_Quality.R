if(!require(Trumpet)) devtools::install_github("skyhorsetomoon/Trumpet"); suppressMessages(library(Trumpet))

args <- commandArgs(trailingOnly = TRUE)

cat(args[2], "111", args[3], "222", args[4], "333", "\n")

outpath <- paste0(args[1], "/", strsplit(basename(args[4]), "\\.")[[1]][1])
trumpet_report <- Trumpet_report(IP_BAM = unlist(strsplit(args[2], ",")), 
                                 Input_BAM = unlist(strsplit(args[3], ",")), 
                                 contrast_IP_BAM =NULL, 
                                 contrast_Input_BAM = NULL,
                                 condition1 = "sample",condition2 = NULL,
                                 GENE_ANNO_GTF = args[4], OUTPUT_DIR = outpath)

system(paste0("cp ", outpath, "/Trumpet_report.html ", args[5]))
system(paste0("rm -rf ", outpath))
