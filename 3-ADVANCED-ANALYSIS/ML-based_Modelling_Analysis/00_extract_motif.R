extract_motif_list <- function(seqRegion,genome,motif="RRACH",strictN=T,setdiff=NULL){
  negative_bag <- list()
  seq_num <- length(seqRegion)
  for(i in 1:seq_num){
    negative_bag[[i]] <- extract_motif(seqRegion =seqRegion[i] ,
                                   genome = genome,motif=motif,
                                   strictN=strictN,setdiff=setdiff)
                     
    # if(i%%500==0){cat(round(i/seq_num,2)*100,'% complete\n')}
  }
  negative_bag <- do.call('c',negative_bag[which(lengths(negative_bag)!=0)])
  return(negative_bag)
}



extract_motif <- function(seqRegion,genome,motif="RRACH",strictN=T,setdiff=NULL){
  full.negative.gene <- seqRegion
  FaFile <- genome 
  frank_motif_length <- round((nchar(motif)-1)/2)
  full.negative.gene <- full.negative.gene+frank_motif_length
  geneInformation.seq <- getSeq(FaFile, full.negative.gene)
  full.negative.motif <- vmatchPattern(motif, geneInformation.seq, fixed=F)
  if(length(full.negative.motif[[1]])==0){
    return(full.negative.motif[[1]])
  } 
  indel.num <- round(nchar(motif)/2)
  if (as.logical(full.negative.gene@strand=="+")) {
    full.negative.sites <- full.negative.gene@ranges@start+full.negative.motif[[1]]@start+indel.num-1
  }else{
    full.negative.sites <- full.negative.gene@ranges@start+
      full.negative.gene@ranges@width-
      full.negative.motif[[1]]@start-indel.num
  }
  full.negative.region <- GRanges(seqnames = full.negative.gene@seqnames,
                                  ranges = IRanges(start = full.negative.sites,width = 1),
                                  strand = full.negative.gene@strand)+indel.num

  # if(!is.null(mark_sites)){
  #   mark_num <- rep(0,length(full.negative.region))
  #   mark_num[findOverlaps(full.negative.region,mark_sites)@from] <- 2
  #   full.negative.region <- GRanges(full.negative.region,
  #                                   motif=getSeq(FaFile, full.negative.region),
  #                                   peak=as.character(full.negative.gene),
  #                                   mark_sites=mark_num)
  # }else{
    full.negative.region <- GRanges(full.negative.region,
                                    motif=getSeq(FaFile, full.negative.region),
                                    peak=as.character(full.negative.gene-frank_motif_length)
    )
  #   )
  # }
  if(strictN){
    full.negative.region <- full.negative.region[which(!grepl("N",full.negative.region@elementMetadata@listData$motif))]
  }else{
    full.negative.region <- full.negative.region[which(!full.negative.region@elementMetadata@listData$motif==paste0(rep("N",indel.num*2+1),collapse = ""))]
  }
  if(!is.null(setdiff)){
    out_num <- findOverlaps(full.negative.region,setdiff)@from
    full.negative.region@elementMetadata@listData$mark_sites[out_num] <- 1 #full.negative.region[setdiff(seq_along(full.negative.region),out_num)]
  }
  return(full.negative.region)
}
