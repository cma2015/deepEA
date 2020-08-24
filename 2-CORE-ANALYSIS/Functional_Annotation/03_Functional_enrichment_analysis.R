library(grid)
library(enrichplot)
library(clusterProfiler)
library(stringr)
library(AnnotationHub)
library(grid)
library(biomaRt)

args <- commandArgs(T)
species <- args[1]
inputGene <- args[2]
local_hub <- args[3]
type <- args[4]
pvalueCutoff <- args[5]
qvalueCutoff <- args[6]
nodesNumber <- args[7]
ont <- args[8]
outFigure <- args[9]
outTable <- args[10]

output <- "/home/galaxy/tools/2-CORE-ANALYSIS/Functional_Annotation/"
genes <- read.table(file = inputGene, sep = '\t', header = F, quote = '',
                    stringsAsFactors = F)
genes <- genes$V1

if(local_hub == "local"){
  hub <- AnnotationHub::AnnotationHub(AnnotationHub::AnnotationHub(hub = AnnotationHub::getAnnotationHubOption("CACHE")),
                                      localHub=T)
}else{
  hub <- AnnotationHub()
}

# Convert Ensembl ID to Entrez ID
Ensembl2Entrez <- function(species_name = "Zea mays", geneID = NULL, drop = TRUE){
  	load("/home/galaxy/tools/2-CORE-ANALYSIS/Functional_Annotation/ensembl_species.RData")
	species_type <- species$species[grep(species_name, species$description, ignore.case = TRUE)]
	species_id <- species$dataset[grep(species_name, species$description, ignore.case = TRUE)]
	if(species_type == "ensembl"){
		mart <- useMart(biomart = 'ensembl', dataset = species_id)
	}else{
		biomart <- paste0(species_type, "_mart")
		host <- paste0(species_type, ".ensembl.org")
		mart <- useMart(biomart = biomart, host = host, dataset = species_id)
	}
	IDTable <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), mart = mart)
	idx <- match(geneID, IDTable$ensembl_gene_id)
	entrezID <- IDTable$entrezgene_id[idx]

	if(drop){
		entrezID <- na.omit(entrezID)
	}
	res <- list(entrezID = entrezID, IDTable = IDTable)
	res
}


multiplot <- function(..., plotlist = NULL, cols = 1, layout=NULL) {


	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)

	numPlots = length(plots)

	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) {
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
						ncol = cols, nrow = ceiling(numPlots/cols))
	}

	if (numPlots==1) {
		print(plots[[1]])

	} else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

		# Make each plot, in the correct location
		for (i in 1:numPlots) {
		# Get the i,j matrix positions of the regions that contain this subplot
		matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

		print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
										layout.pos.col = matchidx$col))
		}
	}
}

res <- Ensembl2Entrez(species_name = species,
                           geneID = genes, drop = TRUE)
entrezID <- res$entrezID
IDTable <- res$IDTable


if(length(grep("Oryza sativa", species, ignore.case = TRUE)) != 0){
	species <- "Oryza sativa"
}

AH <-  AnnotationHub::query(hub, c(species,'OrgDb'))[1]
OrgDb_AH <- hub[[AH$ah_id]]

if(type == "GO"){
	pdf(file = outFigure, width = 13, height = 13)
	options(digits = 3)
	ego2 <- enrichGO(gene = entrezID,
					OrgDb = OrgDb_AH,
					keyType = 'ENTREZID',
					ont = ont,
					pAdjustMethod = "fdr",
					pvalueCutoff = as.numeric(pvalueCutoff),
					qvalueCutoff = as.numeric(qvalueCutoff))
	tt <- lapply(ego2@result$geneID, function(x){
		curID <- strsplit(x, "/")[[1]]
		resID <- IDTable$ensembl_gene_id[match(curID, IDTable$entrezgene_id)]
		resID <- paste(resID, collapse = "/")
		})
	ego2@result$geneID <- unlist(tt)
	write.table(ego2@result, file = outTable,
				sep = '\t', quote = F, row.names = F, col.names = T)
	# go_data <- rbind(go_data,cbind(ego2@result[,1:7], i))
	p1 <- barplot(ego2, showCategory = as.numeric(nodesNumber))
	p2 <- dotplot(ego2)
	p3 <- emapplot(ego2)
	multiplot(plotlist = list(p1, p2, p3),cols = 1)
	dev.off()

}else{
	pdf(file = outFigure, height = 13, width = 13)
	ath <- search_kegg_organism(species, by = 'scientific_name')[1,]
	eKEGG <- enrichKEGG(gene = entrezID,
						organism = ath$kegg_code,
						pvalueCutoff = as.numeric(pvalueCutoff))
	tt <- lapply(eKEGG@result$geneID, function(x){
                curID <- strsplit(x, "/")[[1]]
                resID <- IDTable$ensembl_gene_id[match(curID, IDTable$entrezgene_id)]
                resID <- paste(resID, collapse = "/")
                })
        eKEGG@result$geneID <- unlist(tt)
	write.table(eKEGG@result, file = outTable,
				sep = '\t', quote = F, row.names = F, col.names = T)
	eKEGG@result$Description <- eKEGG@result$ID
	p1 <- barplot(eKEGG, showCategory = as.numeric(nodesNumber))
	p2 <- dotplot(eKEGG)
	p3 <- emapplot(eKEGG)
	multiplot(plotlist = list(p1, p2, p3), cols = 1)
	dev.off()
}

