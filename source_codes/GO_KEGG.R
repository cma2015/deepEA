#library(AnnotationHub)
#Rscript gotest.R Arabidopsis /home/malab5/docker/PEAC/R/R/Table_genetable.txt Annotation TAIR 3 TRUE TRUE 
# devtools::install_github("GuangchuangYu/enrichplot")
# devtools::install_github(c("GuangchuangYu/DOSE", "GuangchuangYu/clusterProfiler"))

library(enrichplot)
library(clusterProfiler)
library(stringr)
####################

args=commandArgs(T)
species <- args[1]
gene_names <- args[2]
inputType <-  args[3]
if (inputType=="Annotation") {
  gene_names <- read.delim(gene_names)
  gene_names <- unique(gene_names[,1]) #unique(str_split(as.character(gene_names$tx_name[which(gene_names$exons!=0)]),'\\.',simplify = T)[,1])
}
geneNameType <- args[4]
go_level <- args[5]
local_hub <- args[6]
kegg <- args[7]


Table_GO <- args[8]
Figure_GO <- args[9]
Figure_GO_level <- args[10]
if(kegg=="TRUE"){
Table_Kegg <-  args[11]
}
#####################

# 
  # species <- 'Arabidopsis'
  # gene_names <- read.delim('/home/malab5/docker/PEAC/R/R/Table_genetable.txt')
  # gene_names <- unique(str_split(as.character(gene_names$tx_name[which(gene_names$exons!=0)]),'\\.',simplify = T)[,1])
  # #gene_names <- as.character(unlist(read.delim('/home/malab5/CAtransGene.txt')))
  # geneNameType <- "TAIR"
  # go_level <- 3

# "ARACYC"       "ARACYCENZYME" "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"     "GO"          
# "GOALL"        "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PMID"         "REFSEQ"       "SYMBOL"       "TAIR"   
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
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

if (local_hub=="TRUE") {
  hub <- AnnotationHub::AnnotationHub(AnnotationHub::AnnotationHub(hub = AnnotationHub::getAnnotationHubOption("CACHE")),
                                      localHub=T)
}else{
  hub <- AnnotationHub::AnnotationHub()
}

AH <-  AnnotationHub::query(hub, c(species,'OrgDb'))[1]
OrgDb_AH <- hub[[AH$ah_id]]
length(keys(OrgDb_AH))
#bitr(keys(OrgDb_AH,'ENTREZID')[1], 'ENTREZID', c("REFSEQ", "GO", "ONTOLOGY"), OrgDb_AH)
#columns(OrgDb_AH)

GO_group <- bitr(gene_names, geneNameType, c("REFSEQ", "GO", "ONTOLOGY"), OrgDb_AH)
sample_genes <- bitr(gene_names, geneNameType, 'ENTREZID', OrgDb_AH)[[2]]
####
go_data <- c()
pdf(file = Figure_GO,width = 13,height = 13)
for(i in c("MF", "BP", "CC")){
  options(digits = 3)  
  ego2 <- enrichGO(gene         = sample_genes,
                   OrgDb         = OrgDb_AH,
                   keyType       = 'ENTREZID',
                   ont           = i,
                   pAdjustMethod = "fdr",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
  go_data <- rbind(go_data,cbind(ego2@result[,1:7],i))
  p1 <- barplot(ego2, showCategory=12)+ggtitle(paste0('Sub-Ontologies: ',i))
  p2 <- dotplot(ego2)
  p3 <- emapplot(ego2)
  multiplot(plotlist = list(p1,p2,p3),cols = 1)
  cat("GO plot ",i,'\n')
}
dev.off()
names(go_data)[8] <- "Type"
write.table(go_data,file = Table_GO,quote = F,row.names = F,col.names = T,sep = '\t')


####
pdf(file = Figure_GO_level,width = 8,height = 5)
for(i in c("MF", "BP", "CC")){
ggo <- groupGO(gene     = sample_genes,
               OrgDb    = OrgDb_AH,
               ont      = i,
               level    = as.numeric(go_level),
               readable = TRUE)
print(barplot(ggo, drop=TRUE, showCategory=12)+ggtitle(paste0('Sub-Ontologies: ',i)))
}
dev.off()

#cnetplot(ego2, categorySize="pvalue",showCategory)
if (kegg=='TRUE') {
  ath <- search_kegg_organism(species, by='scientific_name')[1,]
  kk <- enrichKEGG(gene         = sample_genes,
                   keyType      = 'ncbi-geneid',
                   organism     = ath$kegg_code,
                   pvalueCutoff = 0.05)
  kegg_res <- kk@result[,1:7]
  write.table(kegg_res,file = Table_Kegg,quote = F,row.names = F,col.names = T,sep = '\t')
}


# aaa <- geneannot.map(in.ids=sample_genes, in.type="ENTREZID", out.type='PATH', org='ath', unique.map=TRUE, na.rm=TRUE)#, unique.map=TRUE, na.rm=TRUE)
# bbb <- setdiff(unique(unlist(str_split(aaa[,2],"[ ;]"))),c(" ",NA,""))
#table(unlist(gene.idtype.bods))
