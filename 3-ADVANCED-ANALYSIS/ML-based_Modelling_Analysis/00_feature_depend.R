getMilrProb <- function(beta, X, bag) {
  .Call('milr_getMilrProb', PACKAGE = 'milr', beta, X, bag)
}

FindMode <- function(x) {

  ux <- unique(x)

  ux[which.max(tabulate(match(x, ux)))]

}

sample_sequence <- function (Sequence, Subset_GR, BSgnm, Fixed = F, N = NULL, Replace = F)
{
  regions_reduced <- reduce(Subset_GR)
  regions_DNASS <- DNAStringSet(getSeq(BSgnm, regions_reduced))
  regions_GRL <- split(regions_reduced, paste0("EX_", 1:length(regions_reduced)))
  regions_GRL <- regions_GRL[paste0("EX_", 1:length(regions_reduced))]
  MIndx <- vmatchPattern(Sequence, regions_DNASS, fixed = Fixed)
  #print(MIndx)
  MIndx_gr <- GRanges(seqnames = rep(names(regions_GRL), elementNROWS(MIndx)),
                      ranges = unlist(MIndx))
  sequence_on_regions <- mapFromTranscripts(MIndx_gr, regions_GRL)
  mcols(sequence_on_regions) = NULL
  if (is.null(N)) {
  }
  else {
    if (Replace == F)
      N = min(N, length(sequence_on_regions))
    Indx <- sample.int(length(sequence_on_regions), N, replace = Replace)
    sequence_on_regions <- sequence_on_regions[Indx]
  }
  return(sequence_on_regions)
}
properties_map <- function(query_gr, feature_gr, feature_properties, no_map_val = NA, normalize = F) {

  fol <- findOverlaps(query_gr,feature_gr)

  return_vec <- rep(no_map_val,length(query_gr))

  features_mapped <- feature_properties[subjectHits(fol)]

  if(is.logical(feature_properties)) {
    Weighted_vec  <- tapply(features_mapped, queryHits(fol), any)
  }else{
    Weighted_vec <- tapply(features_mapped, queryHits(fol), mean)
  }

  return_vec[ as.numeric( names(Weighted_vec) )] <- Weighted_vec

  if(normalize) {
    known_indx <- !is.na(return_vec)
    return_vec[known_indx] = (return_vec[known_indx] - mean(return_vec[known_indx]))/sd(return_vec[known_indx])
  }

  return(return_vec)
}
distance_map <- function(query_gr,
                         subject_gr,
                         end = c("five","three"),
                         maximum_distance = 2000,
                         standardize = T){
  end <- match.arg(end)

  if(end == "five"){
    match_indx <- precede(query_gr,subject_gr,ignore.strand = F)
  }

  if(end == "three"){
    match_indx <- follow(query_gr,subject_gr,ignore.strand = F)
  }

  na_indx <- is.na(match_indx)

  dist_all <- vector("integer",length(query_gr))

  dist_all[!na_indx] <- SummarizedExperiment::distance(query_gr[!na_indx],
                                 subject_gr[match_indx[!na_indx]])

  dist_all[is.na(dist_all)] <- maximum_distance

  dist_all <- pmin(dist_all,maximum_distance)

  if(standardize) {
    return((dist_all-mean(dist_all))/sd(dist_all))
  } else {
    return(dist_all)
  }
}
relative_pos_map <- function(query_gr, feature_grl, no_map_val = NA, standardize = TRUE) {

  return_vec <- rep(no_map_val,length(query_gr))

  map2tx <- mapToTranscripts(query_gr,feature_grl)

  relat_pos <- start(map2tx)/ sum(width(feature_grl))[map2tx$transcriptsHits]

  relpos_idx <- tapply( relat_pos, map2tx$xHits, mean)

  return_vec[as.numeric(names(relpos_idx))] = relpos_idx

  if(standardize) {
    return((return_vec - mean(return_vec,na.rm = T))/sd(return_vec,na.rm = T))
  } else {
    return(return_vec)
  }
}


predictors_annot <-
function (se, txdb, bsgnm, fc = NULL, pc = NULL, struct_hybridize = NULL,
          feature_lst = NULL, motif = c("AAACA", "GAACA", "AGACA",
                                        "GGACA", "AAACT", "GAACT", "AGACT", "GGACT", "AAACC",
                                        "GAACC", "AGACC", "GGACC"), motif_clustering = "DRACH",
          self_clustering = FALSE, hk_genes_list = NULL, isoform_ambiguity_method = c("longest_tx",
                                                                                      "average"), genes_ambiguity_method = c("drop_overlap",
                                                                                                                             "average"), standardization = TRUE)
{
 # print(1)

  isoform_ambiguity_method <- match.arg(isoform_ambiguity_method)
  genes_ambiguity_method <- match.arg(genes_ambiguity_method)
  row_gr =  SummarizedExperiment::rowRanges(se)
  genes_txdb <- GenomicFeatures::genes(txdb)
  exbytx_txdb <- exonsBy(txdb, by = "tx")
  if (isoform_ambiguity_method == "longest_tx") {
    Longest_tx_table <- find_longest_transcript(exbytx_txdb,
                                                txdb)
    Kept_tx_indx <- Longest_tx_table$TXID[Longest_tx_table$longest]
    rm(Longest_tx_table)
  } else {
    Kept_tx_indx <- T
  }
  exbytx_txdb <- exbytx_txdb[Kept_tx_indx]
  if (genes_ambiguity_method == "drop_overlap") {
    exbytx_txdb <- exbytx_txdb[countOverlaps(exbytx_txdb,
                                             exbytx_txdb) == 1]
  }
  exbg_txdb <- exonsBy(txdb, by = "gene")
  if (genes_ambiguity_method == "drop_overlap") {
    keep_indx <- countOverlaps(exbg_txdb, exbg_txdb) == 1
    removed_exbg_txdb <- exbg_txdb[!keep_indx]
    exbg_txdb <- exbg_txdb[keep_indx]
    rm(keep_indx)
  }
  exs_txdb <- unlist(exbytx_txdb)
  txid <- names(exbytx_txdb)
  cds <- cdsBy(txdb, by = "tx")
  cds <- cds[names(cds) %in% txid]
  Stop_codons <- resize(unlist(range(cds)), 3, fix = "end")
  Start_codons <- resize(unlist(range(cds)), 3, fix = "start")
  TSS <- resize(unlist(range(exbytx_txdb)), 1, fix = "start")
  A_idx <- vcountPattern("A", DNAStringSet(getSeq(bsgnm, TSS))) >  0


  TSS <- resize(TSS, 100, fix = "start")
  UTR3 <- threeUTRsByTranscript(txdb)
  UTR3 <- UTR3[names(UTR3) %in% txid]
  UTR5 <- fiveUTRsByTranscript(txdb)
  UTR5 <- UTR5[names(UTR5) %in% txid]
  Feature_matrix = data.frame(UTR5 = row_gr %over% UTR5)
  Intron <- intronsByTranscript(txdb)
  scale_i <- function(x, stand = T) {
    if (stand)
      return((x - mean(x))/sd(x))
    else return(x)
  }
  i = 1
  Speak <- function(Fname, I) {
   # cat(paste0("feature ", i, " : ", Fname, " is generated.\n"))
    return(I + 1)
  }
  i = Speak("5'utr", i)
  Feature_matrix$UTR3 <- row_gr %over% UTR3
  i = Speak("3'utr", i)
  Feature_matrix$cds <- row_gr %over% cds
  i = Speak("cds", i)
  Feature_matrix$Stop_codons <- row_gr %over% (Stop_codons +
                                                 100)
  i = Speak("stop codons 203bp", i)
  Feature_matrix$Start_codons <- row_gr %over% (Start_codons +
                                                  100)
  i = Speak("start codons 203bp", i)
  Feature_matrix$TSS <- row_gr %over% TSS
  i = Speak("downstream 100bp of transcription start site",
            i)
  Feature_matrix$TSS_A <- row_gr %over% TSS[A_idx]
  i = Speak("downstream 100bp of transcription start site with sequence A",
            i)
  Feature_matrix$exon_stop <- row_gr %over% subsetByOverlaps(exbg_txdb,
                                                             Stop_codons)
  i = Speak("exon with stop codon", i)
  exbytx_unlist <- unlist(exbytx_txdb)
  Feature_matrix$alternative_exon <- row_gr %over% subsetByOverlaps(exs_txdb,
                                                                    unlist(Intron) + 1, type = "within", maxgap = 0L)
  i = Speak("alternatively spliced exons", i)
  Feature_matrix$constitutive_exon <- row_gr %over% subsetByOverlaps(exbytx_unlist,
                                                                     unlist(Intron) + 1, type = "within", maxgap = 0L, invert = T)
  i = Speak("constitutively spliced exons", i)
  ex_ranks <- exbytx_unlist$exon_rank
  names(ex_ranks) = 1:length(ex_ranks)
  Idx_last_exon <- tapply(ex_ranks, names(exbytx_unlist), function(x) names(x)[which.max(x)])
  Idx_last_exon <- as.numeric(Idx_last_exon[unique(names(exbytx_unlist))])
  Indx_last_exon <- vector("logical", length(ex_ranks))
  Indx_last_exon[Idx_last_exon] <- T
  Indx_first_exon <- ex_ranks == 1
  Feature_matrix$internal_exon <- row_gr %over% exbytx_unlist[!(Indx_last_exon |
                                                                  Indx_first_exon)]
  i = Speak("internal exons", i)
  long_exs_txdb <- exs_txdb[width(exs_txdb) > 400]
  Feature_matrix$long_exon <- row_gr %over% long_exs_txdb
  i = Speak("long exons (length > 400bp)", i)
  last_exons <- exbytx_unlist[Indx_last_exon]
  last_exons_400bp <- resize(last_exons[width(last_exons) >=
                                          400], 400, fix = "start")
  Feature_matrix$last_exon <- row_gr %over% last_exons
  i = Speak("the last exon", i)
  Feature_matrix$last_exon_400bp <- row_gr %over% last_exons_400bp
  i = Speak("5' start 400 bp of the last exon", i)
  Feature_matrix$last_exon_sc400 <- row_gr %over% subsetByOverlaps(last_exons_400bp,
                                                                   Stop_codons)
  i = Speak("5' start 400 bp of the last exon including stop codons",
            i)
  Feature_matrix$intron <- row_gr %over% Intron[names(Intron) %in%
                                                  txid]
  i = Speak("intron", i)
  Feature_matrix$pos_UTR5 <- relative_pos_map(row_gr, UTR5,
                                              0, standardization)
  i = Speak("relative positioning on 5'utr", i)
  Feature_matrix$pos_UTR3 <- relative_pos_map(row_gr, UTR3,
                                              0, standardization)
  i = Speak("relative positioning on 3'utr", i)
  Feature_matrix$pos_cds <- relative_pos_map(row_gr, cds, 0,
                                             standardization)
  i = Speak("relative positioning on cds", i)
  exs_grl <- GenomicRanges::split(exs_txdb, 1:length(exs_txdb))
  Feature_matrix$pos_exons <- relative_pos_map(row_gr, exs_grl,
                                               0, standardization)
  i = Speak("relative positioning on exon", i)
  Splicing_junctions <- reduce(c(resize(exbytx_unlist, 1, fix = "start"),
                                 resize(exbytx_unlist, 1, fix = "end")), min.gapwidth = 0L)
  Splicing_junctions <- subsetByOverlaps(Splicing_junctions,
                                         c(resize(TSS, 1, fix = "start"), resize(unlist(range(exbytx_txdb)),
                                                                                 1, fix = "end")), invert = T)
  Feature_matrix$dist_sj_5_p2000 <- distance_map(row_gr, Splicing_junctions,
                                                 end = "five", maximum_distance = 2000, standardize = standardization)
  i = Speak("distance to the upstream (5' end) splicing junction",
            i)
  Feature_matrix$dist_sj_3_p2000 <- distance_map(row_gr, Splicing_junctions,
                                                 end = "three", maximum_distance = 2000, standardize = standardization)
  i = Speak("distance to the downstream (3' end) splicing junction",
            i)
  Feature_matrix$length_UTR3 <- properties_map(query_gr = row_gr,
                                               feature_gr = UTR3, feature_properties = sum(width(UTR3)),
                                               no_map_val = NA, normalize = standardization)
  Feature_matrix$length_UTR3[is.na(Feature_matrix$length_UTR3)] = 0
  i = Speak("3'UTR length (z-score)", i)
  Feature_matrix$length_UTR5 <- properties_map(query_gr = row_gr,
                                               feature_gr = UTR5, feature_properties = sum(width(UTR5)),
                                               no_map_val = NA, normalize = standardization)
  Feature_matrix$length_UTR5[is.na(Feature_matrix$length_UTR5)] = 0
  i = Speak("5'UTR length (z-score)", i)
  Feature_matrix$length_cds <- properties_map(query_gr = row_gr,
                                              feature_gr = cds, feature_properties = sum(width(cds)),
                                              no_map_val = NA, normalize = standardization)
  Feature_matrix$length_cds[is.na(Feature_matrix$length_cds)] = 0
  i = Speak("cds length (z-score)", i)
  Feature_matrix$length_gene_ex <- properties_map(query_gr = row_gr,
                                                  feature_gr = exbg_txdb, feature_properties = sum(width(exbg_txdb)),
                                                  no_map_val = NA, normalize = standardization)
  Feature_matrix$length_gene_ex[is.na(Feature_matrix$length_gene_ex)] = 0
  i = Speak("gene length-exons (z-score)", i)
  Feature_matrix$length_gene_full <- properties_map(query_gr = row_gr,
                                                    feature_gr = genes_txdb, feature_properties = width(genes_txdb),
                                                    no_map_val = NA, normalize = standardization)
  Feature_matrix$length_gene_full[is.na(Feature_matrix$length_gene_full)] = 0
  i = Speak("gene length full transcript (z-score)", i)
  if (any(width(row_gr) != 1)) {
    warning("At least 1 range with width > 1, the motifs are not attached.")
  }
  else {
    is_motif <- function(motif, dss, exact = F) vcountPattern(DNAString(motif),
                                                              dss, fixed = exact) > 0
    DSS <- DNAStringSet(getSeq(bsgnm, row_gr + 2))
    for (m_i in motif) {
      Feature_matrix[[m_i]] <- is_motif(m_i, DSS, F)
      i = Speak(paste0("motif --- ", m_i), i)
    }
  }
  if (self_clustering) {
    Feature_matrix$clust_f1000 <- as.numeric(scale_i(countOverlaps(row_gr +
                                                                     1000, row_gr) - 1, standardization))
    i = Speak("clustering indicators --- number of neighboors within 1000bp flanking regions",
              i)
    Feature_matrix$clust_f100 <- as.numeric(scale_i(countOverlaps(row_gr +
                                                                    100, row_gr) - 1, standardization))
    i = Speak("clustering indicators ---  number of neighboors within 100bp flanking regions",
              i)
    dist_self <- rep(2000, length(row_gr))
    match_obj <- distanceToNearest(row_gr)
    dist_self[queryHits(match_obj)] <- mcols(match_obj)$distance
    Feature_matrix$dist_nearest_p2000 <- as.numeric(scale_i(pmin(dist_self,
                                                                 2000), standardization))
    i = Speak("clustering indicators --- distance to the nearest neigboors (peaked at 2000bp)",
              i)
    Feature_matrix$dist_nearest_p200 <- as.numeric(scale_i(pmin(dist_self,
                                                                200), standardization))
    i = Speak("clustering indicators --- distance to the nearest neigboors (peaked at 200bp)",
              i)
    rm(match_obj)
    rm(dist_self)
  }
  if (!is.null(motif_clustering)) {
    tx_reduced <- reduce(transcripts(txdb))
    row_gr_expanded <- reduce(row_gr + 2000)
    motif_regions <- GenomicRanges::intersect(row_gr_expanded, tx_reduced)
    rm(row_gr_expanded, tx_reduced)
    # test1 <<- motif_clustering
    # test2 <<- motif_regions
    # test3 <<- bsgnm
    motif_transcripts <- sample_sequence(motif_clustering,
                                         motif_regions, bsgnm)
    rm(motif_regions)
    Feature_matrix[[paste0("clust_", motif_clustering, "_f1000")]] <- as.numeric(scale_i(countOverlaps(row_gr +
                                                                                                         1000, motif_transcripts), standardization))
    i = Speak(paste0("motif clustering --- number of ", motif_clustering,
                     " motif neighboors within 1000bp flanking regions"),
              i)
    Feature_matrix[[paste0("clust_", motif_clustering, "_f100")]] <- as.numeric(scale_i(countOverlaps(row_gr +
                                                                                                        100, motif_transcripts), standardization))
    i = Speak(paste0("motif clustering ---  number of ",
                     motif_clustering, " motif neighboors within 100bp flanking regions"),
              i)
    motif_transcripts <- subsetByOverlaps(motif_transcripts,
                                          row_gr, invert = T)
    match_obj <- distanceToNearest(row_gr, motif_transcripts)
    dist_motif <- rep(2000, length(row_gr))
    dist_motif[queryHits(match_obj)] <- mcols(match_obj)$distance
    dist_self <- rep(2000, length(row_gr))
    match_obj <- distanceToNearest(row_gr)
    dist_self[queryHits(match_obj)] <- mcols(match_obj)$distance
    less_index <- dist_self < dist_motif
    dist_motif[less_index] <- dist_self[less_index]
    Feature_matrix[[paste0("dist_", motif_clustering, "_p2000")]] <- as.numeric(scale_i(pmin(dist_motif,
                                                                                             2000), standardization))
    i = Speak(paste0("motif clustering --- distance to the nearest ",
                     motif_clustering, " motif (peaked at 2000bp)"), i)
    Feature_matrix[[paste0("dist_", motif_clustering, "_p200")]] <- as.numeric(scale_i(pmin(dist_motif,
                                                                                            200), standardization))
    i = Speak(paste0("motif clustering --- distance to the nearest ",
                     motif_clustering, " motif (peaked at 200bp)"), i)
    rm(match_obj)
    rm(dist_self)
    rm(dist_motif)
    rm(less_index)
    rm(motif_transcripts)
  }
  if (!is.null(pc)) {
    Feature_matrix$PC_1bp <- score(pc, row_gr)
    Feature_matrix$PC_1bp[is.na(Feature_matrix$PC_1bp)] = mean(na.omit(Feature_matrix$PC_1bp))
    Feature_matrix$PC_1bp = as.numeric(scale_i(Feature_matrix$PC_1bp,
                                               standardization))
    i = Speak("phast cons scores 1bp", i)
    Feature_matrix$PC_101bp <- score(pc, row_gr + 50)
    Feature_matrix$PC_101bp[is.na(Feature_matrix$PC_101bp)] = mean(na.omit(Feature_matrix$PC_101bp))
    Feature_matrix$PC_101bp = as.numeric(scale_i(Feature_matrix$PC_101bp,
                                                 standardization))
    i = Speak("phast cons scores 101bp", i)
  }
  if (!is.null(fc)) {
    suppressWarnings(Feature_matrix$FC_1bp <- score(fc, row_gr))
    Feature_matrix$FC_1bp[is.na(Feature_matrix$FC_1bp)] = mean(na.omit(Feature_matrix$FC_1bp))
    Feature_matrix$FC_1bp = as.numeric(scale_i(Feature_matrix$FC_1bp,
                                               standardization))
    i = Speak("fitness consequences scores 1bp z score",
              i)
    suppressWarnings(Feature_matrix$FC_101bp <- score(fc,
                                                      row_gr + 50))
    Feature_matrix$FC_101bp[is.na(Feature_matrix$FC_101bp)] = mean(na.omit(Feature_matrix$FC_101bp))
    Feature_matrix$FC_101bp = as.numeric(scale_i(Feature_matrix$FC_101bp,
                                                 standardization))
    i = Speak("fitness consequences scores 101bp z score",
              i)
  }
  if (is.null(struct_hybridize)) {
  }
  else {
    Feature_matrix$struct_hybridize <- row_gr %over% struct_hybridize
    i = Speak("RNA structure --- predicted hybridized region",
              i)
    Feature_matrix$struct_loop <- row_gr %over% infer_loop(struct_hybridize)
    i = Speak("RNA structure --- inferred loop structures between hybridized region",
              i)
  }
  if (!is.null(feature_lst)) {
    for (i_flst in names(feature_lst)) {
      suppressWarnings(Feature_matrix[[i_flst]] <- row_gr %over%
                         feature_lst[[i_flst]])
      i = Speak(paste0("annotation feature --- ", i_flst),
                i)
    }
  }
  coding_tx <- names(cds)
  exbytx_nc <- exbytx_txdb[!names(exbytx_txdb) %in% coding_tx]
  lnc_idx <- sum(width(exbytx_nc)) > 200
  Feature_matrix$sncRNA <- row_gr %over% exbytx_nc[!lnc_idx]
  i = Speak("snc RNA (<= 200bp)", i)
  Feature_matrix$lncRNA <- row_gr %over% exbytx_nc[lnc_idx]
  i = Speak("lnc RNA (> 200bp)", i)
  txbygenes <- transcriptsBy(txdb, by = "gene")
  Feature_matrix$isoform_num <- properties_map(query_gr = row_gr,
                                               feature_gr = txbygenes, feature_properties = pmin(elementNROWS(txbygenes),
                                                                                                 20), no_map_val = 0, normalize = standardization)
  i = Speak("isoform number z score", i)
  Feature_matrix$exon_num <- properties_map(query_gr = row_gr,
                                            feature_gr = exbytx_txdb, feature_properties = elementNROWS(exbytx_txdb),
                                            no_map_val = 0, normalize = standardization)
  i = Speak("exon number z score", i)
  if (!is.null(hk_genes_list)) {
    Feature_matrix$HK_genes <- row_gr %over% exbg_txdb[names(exbg_txdb) %in%
                                                         hk_genes_list]
    i = Speak("house keeping genes", i)
  }
  exbg_gr <- unlist(exbg_txdb)
  exs_GC_freq <- letterFrequency(DNAStringSet(getSeq(bsgnm,
                                                    exbg_gr)), letters = "CG")
  Genes_GC_freq <- tapply(exs_GC_freq, names(exbg_gr), sum)
  Genes_length <- tapply(width(exbg_gr), names(exbg_gr), sum)
  Genes_GC_cont = Genes_GC_freq/Genes_length
  Feature_matrix$GC_cont_genes <- properties_map(query_gr = row_gr,
                                                 feature_gr = exbg_txdb, feature_properties = Genes_GC_cont,
                                                 no_map_val = 0.5, normalize = standardization)
  i = Speak("gene level GC content z score", i)
  Feature_matrix$GC_cont_101bp <- as.numeric(letterFrequency(DNAStringSet(getSeq(bsgnm,
                                                                                row_gr + 50)), letters = "CG", as.prob = T))
  if (standardization) {
    Feature_matrix$GC_cont_101bp <- (Feature_matrix$GC_cont_101bp -
                                       mean(Feature_matrix$GC_cont_101bp))/sd(Feature_matrix$GC_cont_101bp)
  i = Speak("101bp GC content z score", i)
  }
  if (standardization) {
    Feature_matrix$GC_cont_101bp_abs <- abs(Feature_matrix$GC_cont_101bp)
    i = Speak("absolute value of the 101bp GC content z score",
              i)
  }
  if (genes_ambiguity_method == "drop_overlap") {
    Feature_matrix[row_gr %over% removed_exbg_txdb, ] = NA
  }
  mcols(se) = Feature_matrix
  return(se)
}

sequence_features <- function (query_gr, bsgnm)
{
  stopifnot(all(width(query_gr) == width(query_gr)[1]))
  bsgnm_view <- getSeq(bsgnm, query_gr)
  sequences <- as.character(DNAStringSet(bsgnm_view))
  sequences_lst <- strsplit(sequences, "")
  chemical_properties <- function(NT = c("A", "T", "C", "G",
                                         "N")) {
    nt <- match.arg(NT)
    feature_chem = c(nt %in% c("A", "G"), nt %in% c("A",
                                                    "C"), nt %in% c("A", "T"))
    return(feature_chem)
  }
  nt_frequency <- function(NT_vec) {
    D <- vector("numeric", length(NT_vec))
    for (i in seq_along(NT_vec)) {
      D[i] <- mean(NT_vec[1:i] == NT_vec[i])
    }
    return(D)
  }
  features_lst <- lapply(sequences_lst, function(x) {
    M_features <- rbind(sapply(x, chemical_properties), nt_frequency(x))
    return(as.vector(M_features))
  })
  nt_num <- width(query_gr)[1]
  feature_M <- matrix(unlist(features_lst), ncol = nt_num *
                        4, byrow = T)
  feature_M <- as.data.frame(feature_M)
  colnames(feature_M) <- paste0(rep(c("purine", "amino", "weakHyb",
                                      "cumFreq")), "_", rep(1:nt_num, each = 4))
  return(as.data.frame(feature_M))
}


runBioSeq <- function(method, inputSeq, out){
  BioSeqDic <- "/home/galaxy/tools/DeepEA_software/BioSeq-Analysis2.0-Seq/"
  if(method == "1mer"){
	  cmd <- paste0("cd ", BioSeqDic, " && /home/miniconda2/bin/python ", BioSeqDic, "feature.py ",
	                inputSeq, " RNA -method Kmer -k 1 -out ", out)
  }else if(method == "2mer"){
	  cmd <- paste0("cd ", BioSeqDic, " && /home/miniconda2/bin/python ", BioSeqDic, "feature.py ",
	                inputSeq, " RNA -method Kmer -k 2 -out ", out)
  }else if(method == "3mer"){
	  cmd <- paste0("cd ", BioSeqDic, " && /home/miniconda2/bin/python ", BioSeqDic, "feature.py ",
	                inputSeq, " RNA -method Kmer -k 3 -out ", out)
  }else if(method == "4mer"){
	  cmd <- paste0("cd ", BioSeqDic, " && /home/miniconda2/bin/python ", BioSeqDic, "feature.py ",
	                inputSeq, " RNA -method Kmer -k 4 -out ", out)
  }else{
	  cmd <- paste0("cd ", BioSeqDic, " && /home/miniconda2/bin/python ", BioSeqDic, "feature.py ",
	                inputSeq, " RNA -method ", method, " -out ", out)
  }
  system(cmd)
  featureMat <- read.table(file = out, sep = "\t", header = F, quote = "",
                           stringsAsFactors = F)
  sequences <- seqinr::read.fasta(inputSeq, as.string = T)
  rownames(featureMat) <- names(sequences)
  colnames(featureMat) <- paste0(method, "_", 1:ncol(featureMat))
  featureMat <- as.matrix(featureMat)
  return(featureMat)
}
