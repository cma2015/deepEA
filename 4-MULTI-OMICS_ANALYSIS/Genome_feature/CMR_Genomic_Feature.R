library(flexdashboard)
args <- commandArgs(T)
CMRGeneDir <- args[1]
GTFDir <- args[2]
genome <- args[3]
seqLength <- as.numeric(args[4])
stepSize <- as.numeric(args[5])
ratio <- as.numeric(args[6])
featureType <- as.character(args[7])

mainDic <- "/home/zhaijj/DeepEA/galaxy/tools/4-MULTI-OMICS_ANALYSIS/Genome_feature/"
if(featureType == "geneLength"){
  rmarkdown::render(paste0(mainDic, "00_geneLength.Rmd"),
                    output_file = 'CMR_Genomic_Feature.html',
                    output_dir = mainDic)
}else if(featureType == "exonLength"){
  rmarkdown::render(paste0(mainDic, "01_exonLength.Rmd"),
                    output_file = 'CMR_Genomic_Feature.html',
                    output_dir = mainDic)
}else if(featureType == "intronLength"){
  rmarkdown::render(paste0(mainDic, "02_intronLength.Rmd"),
                    output_file = 'CMR_Genomic_Feature.html',
                    output_dir = mainDic)
}else if(featureType == "exonNumber"){
    rmarkdown::render(paste0(mainDic, "03_exonNumber.Rmd"),
                      output_file = 'CMR_Genomic_Feature.html',
                      output_dir = mainDic)
}else{
	rmarkdown::render(paste0(mainDic, "04_GCContent.Rmd"),
	                  output_file = 'CMR_Genomic_Feature.html',
	                  output_dir = mainDic)
}