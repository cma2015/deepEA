require(rtracklayer)
library(seqinr)
library(stringr)
## load raw data and transform
## load raw data and transform
options(scipen = 200)
args=commandArgs(T)

#input.seq <- args[1]
gtf <- args[1]
genome <- args[2]
motif <- args[3]
out <-  'out'

# gtf <- "/home/malab5/data/arabidopsis/Araport11_GFF3_genes_transposons.201606.gtf" 
# genome <- "/home/malab5/data/arabidopsis/TAIR10_Chr.all.fasta"
# motif <- "[AG][AG]AC[ACT]"
# out <- "genome.fa"

merger.base <- function(motif){
  R<- '[AG]'
  Y<- '[CT]'
  M<- '[AC]'
  K<- '[GT]'
  S<- '[GC]'
  W<- '[AT]'
  H<- '[ATC]'
  B<- '[GTC]'
  V<- '[GAC]'
  D<- '[GAT]'
  N<- '[ATCG]'
  
  for(i in c("R","Y","M","K","S","W","H","B","V","D","N")){
    motif <- gsub(pattern = i, replacement = eval(parse(text = i)), x = motif)
  }
  return(motif)
}

motif <- merger.base(motif)

gtf_file <- import.gff(gtf)
gtf_file <- intersect(gtf_file,gtf_file)

export.bed(gtf_file,con = "peaks_with_sites.bed")

# extract sequence
# bedtools getfasta -s -fi ~/data/rerio/Danio_rerio.Zv9.dna.chromosome.1.fa -bed peaks_with_sites.bed > peaks_with_sites.fa
script <- paste0("bedtools getfasta -s -fi ",genome,' -bed peaks_with_sites.bed -fo peaks_with_sites.fa') 
system(script)
# system("cat  ~/data/rerio/Danio_rerio.Zv9.dna.chromosome.1.fa  > test123.fa")
# read.table('test123.tet')


### 

peak_seq <- read.fasta('peaks_with_sites.fa',as.string = T)
seq.res <- lapply(peak_seq,toupper)
#seq.res <- lapply(seq.res,c2s)
seq.pos <- sapply(seq.res, function(x) words.pos(text = x,pattern = motif))
#seq.pos <- seq.pos[which(lengths(seq.pos)!=0)]

seq.pos.trans <- seq.pos

### + strand    start  +2
for(i in grep(pattern = "+",names(seq.pos),fixed = T)){
  seq_name <- names(seq.pos[i])
  seq_start <- as.numeric(str_split(seq_name,"[:-]",simplify = T)[2])
  seq_chr <- (str_split(seq_name,"[:-]",simplify = T)[1])
  seq_strand <- (str_split(seq_name,"[()]",simplify = T)[2])
  seq.pos.trans[[i]] <- paste(seq_chr,seq.pos[[i]]+seq_start+2,seq_strand,sep = "_")
  if (i%%100==0) {cat(i/length(seq.pos)*100,"% \n")}
}
### + strand      end - 2
for(i in grep(pattern = "(-)",names(seq.pos),fixed = T)){
  seq_name <- names(seq.pos[i])
  seq_end <- as.numeric(str_split(seq_name,"[-(]",simplify = T)[2])
  seq_chr <- (str_split(seq_name,"[:-]",simplify = T)[1])
  seq_strand <- (str_split(seq_name,"[()]",simplify = T)[2])
  seq.pos.trans[[i]] <- paste(seq_chr,seq_end-seq.pos[[i]]-1,seq_strand,sep = "_")
  if (i%%100==0) {cat(i/length(seq.pos)*100,"% \n")}
}
#### extract sequence   
seq.pos.trans.name <- names(seq.pos.trans)

seq.pos.trans.bed <- list()
for(i in 1:length(seq.pos.trans.name)){
  x <- seq.pos.trans.name[i]
  ache <- str_split(seq.pos.trans[[x]],"_",simplify = T)
  bag.site.bed <- cbind(ache[,1],as.numeric(ache[,2])-1,ache[,2],".",paste0(x,":",ache[,2],",",i,",1"),ache[,3])
  seq.pos.trans.bed[[i]] <- bag.site.bed
}
names(seq.pos.trans.bed) <- seq.pos.trans.name
# seq.pos.trans.bed <- sapply(seq.pos.trans.name, function(x){
#   ache <- str_split(seq.pos.trans[[x]],"_",simplify = T)
#   bag.site.bed <- cbind(ache[,1],as.numeric(ache[,2])-1,ache[,2],".",paste0(x,":",ache[,2]),ache[,3])
#   return(bag.site.bed)
# })
# peakTosite <- findOverlaps(peaks_with_sites, m6a_sites_from_peak)
# for (i in unique(peakTosite@from)) {
#   Msites <- peakTosite@to[which(peakTosite@from==i)]
#   Msites <- m6a_sites_from_peak@ranges@start[Msites]
#   M_bag_sites <- !is.na(match(as.numeric(seq.pos.trans.bed[[i]][,2]), Msites))
#   if (length(which(!M_bag_sites))==0) {
#     seq.pos.trans.bed[[i]] <- NA
#     next()
#   }
#   if (is.na(seq.pos.trans.bed[[i]][1,2])) {
#     seq.pos.trans.bed[[i]] <- NA
#     next()
#   }
#   seq.pos.trans.bed[[i]][M_bag_sites,5] <- paste0(seq.pos.trans.bed[[i]][M_bag_sites,5],":Me:T")
#   seq.pos.trans.bed[[i]][!M_bag_sites,5] <- paste0(seq.pos.trans.bed[[i]][!M_bag_sites,5],":nonMe:T")
#   if (i%%100==0) {cat(i/length(unique(peakTosite@from))*100,"% \n")}
# }
## delet error peak ( no T motif may be near )

seq.pos.trans.bed <- seq.pos.trans.bed[which(!is.na(seq.pos.trans.bed))]
seq.pos.trans.bed <- do.call("rbind", seq.pos.trans.bed) 
seq.pos.trans.bed <- seq.pos.trans.bed[which(!is.na(seq.pos.trans.bed[,2])),]

motif.seq.pos.trans.bed <- seq.pos.trans.bed
motif.seq.pos.trans.bed[,2] <- as.numeric(seq.pos.trans.bed[,2])-50
motif.seq.pos.trans.bed[,3] <- as.numeric(seq.pos.trans.bed[,3])+50

write.table(motif.seq.pos.trans.bed,file = "motif.seq.pos.trans.bed",quote = F,sep = '\t',row.names = F,col.names = F)
#bedtools getfasta -s -fi ~/data/rerio/Danio_rerio.Zv9.dna.chromosome.1.fa -bed motif.seq.pos.trans.bed -name > motif.seq.pos.trans.fa
script2 <- paste0("bedtools getfasta -s -fi ",genome,' -bed motif.seq.pos.trans.bed -fo motif.seq.pos.trans.fa') 
system(script2)

motif.seq.pos.trans.fa <- read.fasta("./motif.seq.pos.trans.fa",as.string = T)
names(motif.seq.pos.trans.fa) <- motif.seq.pos.trans.bed[,5]
### check and plot 
#source("~/scriptlib/seqlogo.R")
#logoplot(sequences = motif.seq.pos.trans.fa[1:500])
motif.seq.pos.trans.fa <- lapply(motif.seq.pos.trans.fa, toupper)
write.fasta(names = names(motif.seq.pos.trans.fa),sequences = motif.seq.pos.trans.fa,file.out = 'out')