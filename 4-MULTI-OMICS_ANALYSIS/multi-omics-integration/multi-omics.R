library(flexdashboard)
args <- commandArgs(T)
CMRGeneDir <- args[1]
GTFDir <- args[2]
genome <- args[3]
threshold <- as.numeric(args[4])
featureType <- args[5]

mainDic <- "/home/zhaijj/DeepEA/galaxy/tools/4-MULTI-OMICS_ANALYSIS/multi-omics-integration/"
if(featureType == "geneLength"){
  rmarkdown::render(paste0(mainDic, "00_multi-omics_geneLength.Rmd"),
  					output_file = "CMR_multi-omics_genomic_feature.html",
  					output_dir = mainDic)
}else if(featureType == "exonLength"){
  rmarkdown::render(paste0(mainDic, "01_multi-omics_exonLength.Rmd"),
  					output_file = "CMR_multi-omics_genomic_feature.html",
  					output_dir = mainDic)
}else if(featureType == "intronLength"){
  rmarkdown::render(paste0(mainDic, "02_multi-omics_intronLength.Rmd"),
  					output_file = "CMR_multi-omics_genomic_feature.html",
  					output_dir = mainDic)
}else if(featureType == "exonNumber"){
  rmarkdown::render(paste0(mainDic, "03_multi-omics_exonNumber.Rmd"),
  					output_file = "CMR_multi-omics_genomic_feature.html",
  					output_dir = mainDic)
}else{
  rmarkdown::render(paste0(mainDic, "04_multi-omics_GCContent.Rmd"),
  					output_file = "CMR_multi-omics_genomic_feature.html",
  					output_dir = mainDic)
}
