args=commandArgs(T)
peaks <- args[1]
gff <- args[2]
html <- args[3]
single_case <- args[4]
name <- args[5]
noBins <- as.numeric(args[6])
maximalAmbiguity <- as.numeric(args[7])
minimalComponentLength <- as.numeric(args[8])
minimalNcRNALength <- as.numeric(args[9])
output <- args[10]
mainDic <- '/home/galaxy/tools/2-CORE-ANALYSIS/Functional_Annotation/'

if(single_case){
    if(html){
        rmarkdown::render(paste0(mainDic, "00_CMR_distribution_single_case_html.Rmd"),
                    output_format = 'html_document',
                    output_file = output)
    }else{
        rmarkdown::render(paste0(mainDic, "00_CMR_distribution_single_case_pdf.Rmd"),
                    output_format = 'pdf_document',
                    output_file = output)
    }
}else{
    peaks <- unlist(strsplit(peaks, ','))
    name <- unlist(strsplit(name, ','))
    if(html){
        rmarkdown::render(paste0(mainDic, "00_CMR_distribution_group_case_html.Rmd"),
                    output_format = 'html_document',
                    output_file = output)
    }else{
        rmarkdown::render(paste0(mainDic, "00_CMR_distribution_group_case_pdf.Rmd"),
                    output_format = 'pdf_document',
                    output_file = output)
    }
}


