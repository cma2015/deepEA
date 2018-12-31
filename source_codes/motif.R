## Scanning the motifs in the peak regions
library(rGADEM)
library(motifRG)
library(seqinr)
library(parallel)
library(reshape2)

args=commandArgs(T)

input.seq <- args[1]
genome <- args[2]
motif.method <- args[3]
motif.number <- args[4]
motif.length <- args[5]
plot_motif <- args[6]
motif <- args[7]
Figure_motif <- args[8]
workpath <- args[9]
#briefly <- args[6]
#random.genome <- args[6]

motif.number <- as.numeric(motif.number)
motif.length <- as.numeric(motif.length)
### test data 
# input.seq <- "/home/songjie/PEAC/galaxy/database/files/000/dataset_202.dat"
# genome <- "/home/songjie/PEAC/galaxy/database/files/000/dataset_216.dat"
# motif.method <- 'motifRG'
# motif.number <- "3"
# motif.length <- "5"
# workpath <- '/home/songjie/PEAC/galaxy/tools/Annotation'
#/home/songjie/PEAC/galaxy/tools/Annotation/motif.R /home/songjie/PEAC/galaxy/database/files/000/dataset_202.dat /home/songjie/PEAC/galaxy/database/files/000/dataset_34.dat motifRG 3 5 #TRUE /home/songjie/PEAC/galaxy/database/files/000/dataset_212.dat /home/songjie/PEAC/galaxy/database/files/000/dataset_213.dat  /home/songjie/PEAC/galaxy/tools/Annotation ;


MD.peak.seq <- readDNAStringSet(input.seq)
MD.control.seq <- readDNAStringSet(genome)
MD.control.seq <- MD.control.seq[sample(1:length(MD.control.seq),length(MD.peak.seq),replace = T)]

if (motif.method=='motifRG') {
  category <- c(rep(1, length(MD.peak.seq)), rep(0, length(MD.control.seq)))
  MD.motifs <- findMotif(append(MD.peak.seq, MD.control.seq),category,
                         max.motif=motif.number,enriched=FALSE,
                         max.width=motif.length,is.parallel=T,mc.cores=3)
  results <- summaryMotif(MD.motifs$motifs, MD.motifs$category)
if(length(MD.motifs$motifs)==0){
print("None motif can be discovery,please increase the suquence and genome size")
}
  res.table <- data.frame(motif=rownames(results),results)
  
  # ?createControlRegions
  # MD.control.seq <- createControlRegions(MD.peak.seq)
  # 
  
}
  


if (motif.method=='GADEM') {
  results <- GADEM(Sequences = MD.peak.seq,genome = MD.control.seq, verbose = T,
                   numTop3mer=10,numTop4mer=20,numTop5mer=30,numEM=50,maxSpaceWidth = 5,
                   numGeneration=5,
                   nmotifs=as.numeric(motif.number),
                   extTrim=0,
                   slideWinPWM = as.numeric(motif.length),
                   fEM=0.3,
                   weightType=1
                   )
  res.table <- data.frame(
    motif= toupper(sapply(results@motifList, function(x) x@consensus)), #names(results),
    scores=sapply(results@motifList, function(y) mean(-log(sapply(y@alignList, function(x) x@pval)))),
    signs=sapply(results@motifList, function(y) median((sapply(y@alignList, function(x) x@pval))))
  )
  res.table
}

if (plot_motif=='TRUE') {
  source(paste0(workpath,"/extra_seq.R"))
  source(paste0(workpath,"/seqlogo.R"))
  all_seq <- read.fasta(input.seq,as.string = T)
  #all_seq <- lapply(all_seq, toupper)
 # eval(parse(text = paste0("pdf","('./Figure_motif.","pdf","',height = 3,width = as.numeric(motif.length) )")))
pdf(Figure_motif,height = 3,width = as.numeric(motif.length))
print(res.table$motif)
  for (i in as.character(res.table$motif)) {
    re_motif <- tolower(merger.base(i))
	    motif_pos <- lapply(all_seq, function(x) {
	      #x <- all_seq[[1]]
	      pos <- words.pos(re_motif,x)
	      res <- sapply(pos, function(y) substr(x,y,y+motif.length))
	      return(res)
	    }
	    )
    motif_pos <- motif_pos[which(lengths(motif_pos)!=0)]
    motif_pos <- as.list(unlist(motif_pos))
    logoplot(sequences = motif_pos,type = 'weblogo')
  }
  dev.off()
}



write.table(res.table,file = motif,row.names = F,quote = F,sep = '\t')
