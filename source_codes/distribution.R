### distribute
# test code #
# Rscript distribution.R "./Galaxy99-[m6A_peaks].encodepeak" "./Araport11_GFF3_genes_transposons.201606.gff" "m6A" "pdf" TRUE TRUE TRUE TRUE TRUE TRUE

library(Guitar)
library(RCAS)
library(stringr)
library(GenomicFeatures)

args=commandArgs(T)

peak1 <- args[1]
gff <- args[2]
name <- args[3]
ext <- args[4]

distribute1 <- args[5]
distribute2 <- args[6]
peaks_over_Chromosomes <- args[7]
distanse_between_TSS <- args[8]
RNA_type <- args[9]
RNA_region_type  <- args[10]
outfile <- args[11]
#setwd(workpath)


cat("peak1 == ",peak1,'\n')
cat("gff == ",gff,'\n')
cat("name == ",name,'\n')
cat("distribute1 == ",distribute1,'\n')
cat("distribute2 == ",distribute2,'\n')
cat("peaks_over_Chromosomes == ",peaks_over_Chromosomes,'\n')
cat("distanse_between_TSS == ",distanse_between_TSS,'\n')

figure.1 <- NULL
figure.2 <- NULL
figure.3 <- NULL
figure.4 <- NULL
figure.5 <- NULL
figure.6 <- NULL


if (ext=="pdf") {
  n=1
}else{
  n=50
}
#### variation function ###
getTxdbFeaturesFromGRanges <- function (gffData,txdb) {
  transcripts <- GenomicFeatures::transcripts(txdb)
  # if (!is.null(gffData$transcript_id)) {
  #   m <- match(transcripts$tx_name, gffData$transcript_id)
  # }else{
  #   m <- match(transcripts$tx_name, gffData$ID)
  # }
  # m <- na.omit(m)
  # transcripts$gene_name <- gffData[m]$gene_name
  transcripts$gene_name <- str_split(transcripts@elementMetadata@listData$tx_name,"\\.",simplify = T)[,1]
  tmp <- GenomicFeatures::exonsBy(x = txdb, by = "tx", use.names = TRUE)
  exons <- BiocGenerics::unlist(tmp)
  exons$tx_name <- names(exons)
#  exons$gene_name <- gffData[match(names(exons), gffData$transcript_id)]$gene_name
  exons$gene_name <- str_split(exons@elementMetadata@listData$tx_name,"\\.",simplify = T)[,1]
  
  tmp <- GenomicFeatures::intronsByTranscript(x = txdb, use.names = TRUE)
  introns <- BiocGenerics::unlist(tmp)
  introns$tx_name <- names(introns)
  # m <- match(names(introns), gffData$transcript_id)
  # introns$gene_name <- gffData[m]$gene_name
  introns$gene_name <- str_split(introns@elementMetadata@listData$tx_name,"\\.",simplify = T)[,1]
  
  promoters <- GenomicFeatures::promoters(txdb)
  # m <- match(promoters$tx_name, gffData$transcript_id)
  # promoters$gene_name <- gffData[m]$gene_name
  promoters$gene_name <- str_split(promoters@elementMetadata@listData$tx_name,"\\.",simplify = T)[,1]
  
  tmp <- range(GenomicFeatures::fiveUTRsByTranscript(x = txdb, 
                                                     use.names = TRUE))
  fiveUTRs <- BiocGenerics::unlist(tmp)
  fiveUTRs$tx_name <- names(fiveUTRs)
  # m <- match(names(fiveUTRs), gffData$transcript_id)
  # fiveUTRs$gene_name <- gffData[m]$gene_name
  fiveUTRs$gene_name <- str_split(fiveUTRs@elementMetadata@listData$tx_name,"\\.",simplify = T)[,1]
  
  tmp <- range(GenomicFeatures::threeUTRsByTranscript(x = txdb, 
                                                      use.names = TRUE))
  threeUTRs <- BiocGenerics::unlist(tmp)
  threeUTRs$tx_name <- names(threeUTRs)
  # m <- match(names(threeUTRs), gffData$transcript_id)
  # threeUTRs$gene_name <- gffData[m]$gene_name
  threeUTRs$gene_name <- str_split(threeUTRs@elementMetadata@listData$tx_name,"\\.",simplify = T)[,1]
  
  tmp <- GenomicFeatures::cdsBy(x = txdb, by = "tx", use.names = TRUE)
  cds <- BiocGenerics::unlist(tmp)
  cds$tx_name <- names(cds)
  # cds$gene_name <- gffData[match(names(cds), gffData$transcript_id)]$gene_name
  cds$gene_name <- str_split(cds@elementMetadata@listData$tx_name,"\\.",simplify = T)[,1]
  
  txdbFeatures = list(transcripts = transcripts, exons = exons, 
                      promoters = promoters, fiveUTRs = fiveUTRs, introns = introns, 
                      cds = cds, threeUTRs = threeUTRs)
  return(txdbFeatures)
}

# getFeatureBoundaryCoverage <- 
#   function (queryRegions, featureCoords, flankSize = 500, sampleN = 0) 
#   {
#     if (sampleN > 0 && sampleN < length(featureCoords)) {
#       featureCoords <- sort(featureCoords[sample(length(featureCoords), 
#                                                  sampleN)])
#     }
#     fivePrimeFlanks <- GenomicRanges::flank(x = featureCoords, 
#                                             width = flankSize, start = TRUE, both = TRUE)
#     threePrimeFlanks <- GenomicRanges::flank(x = featureCoords, 
#                                              width = flankSize, start = FALSE, both = TRUE)
#     cvgFivePrime <- genomation::ScoreMatrix(target = queryRegions, 
#                                             windows = fivePrimeFlanks, strand.aware = TRUE)
#     cvgThreePrime <- genomation::ScoreMatrix(target = queryRegions, 
#                                              windows = threePrimeFlanks, strand.aware = TRUE)
#     mdata <- data.frame(fivePrime = colSums(cvgFivePrime), threePrime = colSums(cvgThreePrime), 
#                         bases = c(-flankSize:(flankSize - 1)))
#     return(mdata)
#   }

###########  test data ##########
# peak1 <- "./Galaxy99-[m6A_peaks].encodepeak"
# 
# gff <- "./Araport11_GFF3_genes_transposons.201606.gff"
# name <- "m6A"
# ext <- "pdf"

######### prepare data ############
####

testregion <- read.table(peak1,header=F)
testregion <- testregion[,1:3]
colnames(testregion) <- c('chrom','start','end')
region <- as(testregion, "GRanges")


#region <- narrowPeaktoGRanges(peak1) 
feature_list <- list(region) 
names(feature_list) <- c(name)
txdb_ara <- makeTxDbFromGFF(gff)
gc_txdb_ara <- makeGuitarCoordsFromTxDb(txdb_ara,noBins =100)

gff_region <- import(gff)
txdbFeatures <- getTxdbFeaturesFromGRanges(gffData = gff_region,txdb = txdb_ara)
#gc_txdb_ara@seqnames@values <- as.factor(tolower(gc_txdb_ara@seqnames@values))
################################################################
if(distribute1=="TRUE"){
cat("built distribute1 figure... \n")
eval(parse(text = paste0(ext,"('./Figure1_distribute1.",ext,"',height = 5*n,width = 10*n)")))
figure.1 <- GuitarPlot(feature_list, 
                       GuitarCoordsFromTxDb = gc_txdb_ara,noBins=100)
dev.off()
}



#### peaks over Chromosomes
if(peaks_over_Chromosomes=="TRUE"){
cat("built peaks_over_Chromosomes figure... \n")
figure.2 <-  ChIPseeker::covplot(region,title = "Peaks over Chromosomes")   #specific chr
eval(parse(text = paste0(ext,"('./Figure2_peaks_over_Chromosomes.",ext,"',height = 5*n,width = 10*n)")))
print(figure.2)
dev.off()
}

#### the distanse between TSS/TES and peaks
if(distanse_between_TSS=="TRUE"){
  cat("built distanse_between_TSS figure... \n")
  
# cvg <- getFeatureBoundaryCoverage(queryRegions = region,
#                                   featureCoords = txdbFeatures$transcripts,
#                                   flankSize = 1000, 
#                                   sampleN = 10000)
# 
# yLimit <- (as.integer(max(c(cvg$fivePrime, cvg$threePrime))/10)+1)*10
# 
# figure.3 <- subplot(
#   plot_ly(data = cvg, x = ~bases, y = ~fivePrime, type = 'scatter', mode = 'lines'),
#   plot_ly(data = cvg, x = ~bases, y = ~threePrime, type = 'scatter', mode = 'lines'),
#   margin = 0.05
# ) %>% layout (xaxis = list(title = 'Distance (bp) to TSS'), 
#               xaxis2 = list(title = 'Distance (bp) to TES'), 
#               yaxis = list(title = 'coverage', range = c(0, yLimit)),
#               yaxis2 = list(title = 'coverage', range = c(0, yLimit)),
#               showlegend = FALSE) 
# 
# #eval(parse(text = paste0(ext,"('./distanse_between_TSS.",ext,"',height = 5*n,width = 10*n)")))
# plotly::export(figure.3,"Figure3_distanse_between_TSS.pdf")
# htmlwidgets::saveWidget(as_widget(figure.3), "distanse_between_TSS.html")
# webshot::webshot(url = "./distanse_between_TSS.html", file = paste0("Figure3_distanse_between_TSS.",ext), vwidth = 10, vheight = 5,zoom = 2)
# 

cvgF <- getFeatureBoundaryCoverage(queryRegions = region, 
                                   featureCoords = txdbFeatures$transcripts, 
                                   flankSize = 1000, 
                                   boundaryType = 'fiveprime', 
                                   sampleN = 10000)
cvgT <- getFeatureBoundaryCoverage(queryRegions = region, 
                                   featureCoords = txdbFeatures$transcripts, 
                                   flankSize = 1000, 
                                   boundaryType = 'threeprime', 
                                   sampleN = 10000)

cvgF$boundary <- 'fiveprime'
cvgT$boundary <- 'threeprime'

df <- rbind(cvgF, cvgT)       



figure.3 <- ggplot2::ggplot(df, aes(x = bases, y = meanCoverage)) + 
  geom_ribbon(fill = 'lightgreen', 
              aes(ymin = meanCoverage - standardError * 1.96, 
                  ymax = meanCoverage + standardError * 1.96)) + 
  geom_line(color = 'black') + 
  facet_grid(~ boundary) + theme_bw(base_size = 14)

eval(parse(text = paste0(ext,"('./Figure3_distanse_between_TSS.",ext,"',height = 5*n,width = 10*n)")))
print(figure.3)
dev.off()

}
# plotly::export(figure.3, "plot.pdf")

####  RNA type
if(RNA_type=="TRUE"){
  cat("built RNA_type figure... \n")
#   
# overlaps <- queryGff(queryRegions = region, gffData = gff_region)  
# #data.table is used to do quick summary operations
# overlaps.dt <- data.table(as.data.frame(overlaps)) 
# label.name <- colnames(overlaps.dt)[which(overlaps.dt=="mirna",arr.ind = T)[1,2]]
# biotype_col <- grep(label.name, colnames(overlaps.dt), value = T)
# df <- overlaps.dt[,length(unique(overlappingQuery)), by = biotype_col]
# colnames(df) <- c("feature", "count")
# df$percent <- round(df$count / length(region) * 100, 1)
# df <- df[order(count, decreasing = TRUE)]
# figure.4 <- plot_ly(data = df, 
#              type = "bar",
#              x = df$feature,
#              y = df$percent,
#              text = paste("count:", df$count), color=df$feature)
# figure.4 <- layout(p = figure.4, 
#        margin = list(l=100, r=100, b=150), 
#        xaxis = list(showticklabels = TRUE,  tickangle = 90), 
#        yaxis = list(title = paste("percentage of query regions,", 
#                                   "n =", length(region))))
# htmlwidgets::saveWidget(as_widget(figure.4), "RNA_type.html")
# webshot::webshot(url = "./RNA_type.html", file = paste0("Figure4_RNA_type.",ext), vwidth = 10, vheight = 5,zoom = 1)

##
overlaps <- as.data.table(queryGff(queryRegions = region, gffData = gff_region))
label.name <- colnames(overlaps)[which(overlaps=="mirna",arr.ind = T)[1,2]]
overlaps <- overlaps[which(!is.na(overlaps[[label.name]])),]
biotype_col <- grep(label.name, colnames(overlaps), value = T)
df <- overlaps[,length(unique(overlappingQuery)), by = biotype_col]
colnames(df) <- c("feature", "count")
df$percent <- round(df$count / length(queryRegions) * 100, 1)
df <- df[order(count, decreasing = TRUE)]
figure.4 <- ggplot2::ggplot(df, aes(x = reorder(feature, -percent), y = percent)) + 
  geom_bar(stat = 'identity', aes(fill = feature)) + 
  geom_label(aes(y = percent + 0.5), label = df$count) + 
  labs(x = 'Transcript Feature', y = paste0('percent overlap (n = ', length(region), ')')) + 
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 90))
eval(parse(text = paste0(ext,"('./Figure4_RNA_type.",ext,"',height = 7*n,width = 10*n)")))
print(figure.4)
dev.off()

}
### RNA region type 
if(RNA_region_type=="TRUE"){
  cat("built RNA_region_type figure... \n")
  
summary <- summarizeQueryRegions(queryRegions = region, 
                                 txdbFeatures = txdbFeatures)
df <- data.frame(summary)
df$percent <- round((df$count / length(region)), 3) * 100
df$feature <- rownames(df)
figure.5 <- ggplot2::ggplot(df, aes(x = reorder(feature, -percent), y = percent)) + 
  geom_bar(stat = 'identity', aes(fill = feature)) + 
  geom_label(aes(y = percent + 3), label = df$count) + 
  labs(x = 'transcript feature', y = paste0('percent overlap (n = ', length(queryRegions), ')')) + 
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 90))
eval(parse(text = paste0(ext,"('./Figure5_RNA_region_type.",ext,"',height = 5*n,width = 8*n)")))
print(figure.5)
dev.off()

# 
# figure.5 <- plot_ly( data = df, 
#               x = rownames(df), 
#               y = df$percent, 
#               type = 'bar',
#               text = paste("count:", df$count), 
#               color = rownames(df)
# )
# figure.5 <- layout(p = figure.5, 
#        xaxis = list(title = 'features'),
#        yaxis = list(title = paste("percentage of query regions,", 
#                                   "n =", length(region)
#        )
#        ), 
#        margin = list(b = 150, r = 50)
# )
# htmlwidgets::saveWidget(as_widget(figure.5), "RNA_region_type.html")
# webshot::webshot(url = "./RNA_region_type.html", file = paste0("Figure5_RNA_region_type.",ext), vwidth = 10, vheight = 5,zoom = 2)

}
### distribute2 
if(distribute2=="TRUE"){
  cat("built distribute2 figure... \n")
  
covList <- calculateCoverageProfileList(queryRegions = region, 
                                        targetRegionsList = txdbFeatures, 
                                        sampleN = length(region))

figure.6 <-ggplot2::ggplot(covList, aes(x = bins, y = meanCoverage)) + 
  geom_ribbon(fill = 'lightgreen', 
              aes(ymin = meanCoverage - standardError * 1.96, 
                  ymax = meanCoverage + standardError * 1.96)) + 
  geom_line(color = 'black') + theme_bw(base_size = 14) +
  facet_wrap( ~ feature, ncol = 3) 

eval(parse(text = paste0(ext,"('./Figure6_distribute2.",ext,"',height = 5*n,width = 10*n)")))
print(figure.6)
dev.off()


###
}
######### plot distribute data ############
if(distribute1=="TRUE"){
  cat("built distribute1 figure... \n")
#  eval(parse(text = paste0(ext,"('./Figure_distribute.",ext,"',height = 5*n,width = 10*n)")))
  pdf(outfile,height = 5*n,width = 10*n)
  figure.1 <- GuitarPlot(feature_list, 
                         GuitarCoordsFromTxDb = gc_txdb_ara,noBins=100)
  print(figure.2)
  print(figure.3)
  print(figure.4)
  print(figure.5)
  print(figure.6)
  dev.off()
}

