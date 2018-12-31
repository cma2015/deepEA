


require(seqinr)
#cdna <- read.fasta("~/data/arabidopsis/TAIR10_chr_all.fas")
#cdnalength <- lengths(cdna)
#names(cdna) <- c(1:5,"Mt","Pt")

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


extra_seq <- function(seq,n){
  posseq <- list()
  if (ncol(seq)!=3) {
    seq <- cbind(seq,"+")
  }
  name <- paste(seq[,1],seq[,2],sep = "_")
  for (i in 1:nrow(seq)){
    pos <- seq[i,c(1,2,3)]
    transcript <- cdna[[as.character(pos[1])]]
    if(length(transcript[[1]])==0){
      cat("transcrpt  ",pos[1],"  is none, so skip\n")
      next()
    }
    transcriptlength <- cdnalength[attr(transcript,"name")]
    
    if (transcriptlength < (as.numeric(pos[2])+n)) {
      endgap = as.numeric(pos[2])+n-transcriptlength
    }else{
      endgap = 0 
    }
    if ((as.numeric(pos[2])-n) <= 0 ) {
      begingap = -(as.numeric(pos[2])-n)+1
    }else{
      begingap = 0 
    }
    sequence <- transcript[as.numeric(seq[i,2])+((-n+begingap):(n-endgap))]
    sequence <- c(rep("n",begingap),sequence,rep("n",endgap))
    if(pos[3] == "-"){
      sequence <- sequence[length(sequence):1]
      sequence <- chartr("atgc","tacg",sequence)
    }
    posseq[[name[i]]] <- sequence
    #cat(i)
  }
  posseq
}

extra_seq_bed <- function(seq,name=NULL,input_bed=T){
  posseq <- list()
  #cat('if input bed file, start position need to add 1 ! \n')
  if(input_bed){seq[,2] <- as.numeric(seq[,2])+1}
  if(is.null(name)&ncol(seq)==4){name <- paste(seq[,1],seq[,2],seq[,3],seq[,4],sep = "_")}
  if(is.null(name)){name <- paste(seq[,1],seq[,2],seq[,3],sep = "_")}
  
  if (ncol(seq)!=4) {
    seq <- cbind(seq,"+")
  }

  for (i in 1:nrow(seq)){
    pos <- seq[i,c(1,2,3,4)]
    transcript <- cdna[[as.character(pos[1])]]
    if(length(transcript[[1]])==0){
      cat("transcrpt  ",pos[1],"  is none, so skip\n")
      next()
    }
    #transcriptlength <- cdnalength[attr(transcript,"name")]
    
    # if (transcriptlength < (as.numeric(pos[3]))) {
    #   endgap = as.numeric(pos[3])-transcriptlength
    # }else{
    #   endgap = 0 
    # }
    # if ((as.numeric(pos[2])) <= 0 ) {
    #   begingap = -(as.numeric(pos[2]))+1
    # }else{
    #   begingap = 0 
    # }
    sequence <- transcript[as.numeric(seq[i,2]):as.numeric(seq[i,3])]
    #sequence <- c(rep("n",begingap),sequence,rep("n",endgap))
    if(pos[4] == "-"){
      sequence <- sequence[length(sequence):1]
      sequence <- chartr("atgc","tacg",sequence)
    }
    posseq[[name[i]]] <- sequence
    #cat(i)
  }
  posseq
}

extra_seq_reverse <- function(seq){
  name <- names(seq)
  new_seq <- lapply(seq, function(x) x[length(x):1] )
  new_seq <- lapply(new_seq, function(x) chartr("atgc","tacg",x) )
  names(new_seq) <- paste0(name,"_reverse")
  new_seq
}
