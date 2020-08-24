library(argparse)
library(flexdashboard)
parser <- ArgumentParser()
parser$add_argument("-cmrDir", default = NULL, dest = "cmrDir", help = "The directory of CMR regions")
parser$add_argument("-method", default = NULL, dest = "method", help = "Normalization method")
parser$add_argument("-ifSubgenome", default = TRUE, dest = "ifSubgenome", help = "If subgenome is available")
parser$add_argument("-subgenome", default = NULL, dest = "subgenome", help = "The directory of subgenome information")
parser$add_argument("-homoeologs", default = NULL, dest = "homoeologs", help = "The directory of homoeologs")
parser$add_argument("-CDS", default = NULL, dest = "CDS", help = "The directory of CDS")
parser$add_argument("-ifGSEA", default = NULL, dest = "ifGSEA", help = "If performing functional gene sets enrichment analysis")
parser$add_argument("-enrichMethod", default = "fisher", dest = "enrichMethod", help = "Functional enrichment method")
parser$add_argument("-geneList", default = NULL, dest = "geneList", help = "The directory of funcitonal gene sets")
parser$add_argument("-bgGene", default = NULL, dest = "bgGene", help = "The directory of background gene sets")
parser$add_argument("-out", default = NULL, dest = "output", help = "The directory of output")

args <- parser$parse_args()
mainDic <- "/home/galaxy/tools/3-ADVANCED-ANALYSIS/Multi-omics_Integrative_Analysis/"
if(args$ifSubgenome == FALSE & args$ifGSEA == FALSE){
	method <- args$method
	cmrDir <- args$cmrDir
	rmarkdown::render(paste0(mainDic, "02-multi-omics-integration-no-subgenome-no-gene-sets.Rmd"), output_file = args$output)
}else if(args$ifSubgenome == FALSE & args$ifGSEA == TRUE){
	method <- args$method
	cmrDir <- args$cmrDir
	enrichMethod <- args$enrichMethod
	geneList <- args$geneList
	bgGeneDir <- args$bgGene
	rmarkdown::render(paste0(mainDic, "02-multi-omics-integration-no-subgenome-with-gene-sets.Rmd"), output_file = args$output)
}else if(args$ifSubgenome == TRUE & args$ifGSEA == FALSE){
	method <- args$method
	cmrDir <- args$cmrDir
	dup_dir <- args$homoeologs
	rmarkdown::render(paste0(mainDic, "02-multi-omics-integration-with-subgenome-no-gene-sets.Rmd"), output_file = args$output)
}else{
	method <- args$method
	cmrDir <- args$cmrDir
	enrichMethod <- args$enrichMethod
	geneList <- args$geneList
	bgGeneDir <- args$bgGene
	subgenomeDir <- args$subgenome
	dupDir <- args$homoeologs
	CDSDir <- args$CDS
	rmarkdown::render(paste0(mainDic, "02-multi-omics-integration-with-subgenome-with-gene-sets.Rmd"), output_file = args$output)
}
