### motif plot
### extract sequence 
library(seqinr)
library(stringr)
library(reshape2)

args=commandArgs(T)

posseq <- args[1]
negseq <- args[2]
type <- args[3]
plottype <- args[4]
#ext <- args[5]
Figure_logo <- args[5]
workpath <- args[6]
#setwd(workpath)


source(paste0(workpath,"/seqlogo.R"))
### test data 
# posseq <- "./extract.fa"
# #negseq <- args[2]
# type <- "seqlogo"
# plottype <- "weblogo"
 ext <- "pdf"

# test code Rscript motif_plot.R "./extract.fa" NULL "seqlogo" "weblogo" "png"
#posseq <- "/home/songjie/PEAC/galaxy/database/files/000/dataset_199.dat"
# /home/songjie/PEAC/galaxy/database/files/000/dataset_199.dat 
# type <- "seqlogo" 
# plottype <- "barlogo" 
#/home/songjie/PEAC/galaxy/database/files/000/dataset_205.dat /
#Figure_logo <- "/home/songjie/PEAC/galaxy/tools/visualization"
#Figure_logo <- "./test.pdf"

if (ext=="pdf") {
  n=1
}else{
  n=50
}


if(type=='seqlogo'){
  posseq <- read.fasta(file = posseq,as.string = T)
  if (plottype=='weblogo') {
    cat(plottype,'\n')
    #eval(parse(text = paste0(ext,"('./Figure_logo.",ext,"',height = 5*n,width = 5*nchar(posseq[[1]])/20*n)")))
    pdf(Figure_logo,height = 5*n,width = 5*nchar(posseq[[1]])/20*n)
    logoplot(sequences = posseq,type = 'weblogo')
    dev.off()
  }
  if (plottype=='barlogo'){
    cat(plottype)
    pdf(Figure_logo,height = 5*n,width = 5*nchar(posseq[[1]])/20*n)
    p <- logoplot(sequences = posseq,type = 'barlogo')
    print(p)
    dev.off()
    }
}


if(type=='twosamplelogo'){
cat(1,'\n')
  posseq <- read.fasta(file = posseq,as.string = T)
  negseq <- read.fasta(file = negseq,as.string = T)
  
#  eval(parse(text = paste0(ext,"('./Figure_logo.",ext,"',height = 5*n,width = 5*nchar(posseq[[1]])/20*n)")))
  pdf(Figure_logo,height = 5*n,width = 5*nchar(posseq[[1]])/20*n)
  p <- twologoplot(pos_sequences = posseq,neg_sequences = negseq)
   print(p)
  dev.off()
  
}
