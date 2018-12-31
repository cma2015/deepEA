library(RandomWalkRestartMH)
library(igraph)
library(stringr)

args=commandArgs(T)

#bed <- import.bed("/home/malab5/docker/PEAC/data/peak_checked_merge.bed")
Ara_PPI_raw <- args[1]
Gene_list <- args[2]
analysis_num <- args[3]
Restart_probability <- args[4]
out_list <- args[5]
out_figure <- args[6]

analysis_num <- as.numeric(analysis_num)
Restart_probability <- as.numeric(Restart_probability)

Ara_PPI_raw <- read.delim(Ara_PPI_raw)
Genefile <- read.delim(Gene_list)
# analysis_num <- 50

Ara_PPI <- Ara_PPI_raw[,1:2]
normalize.multiplex.adjacency <- function (x) {
  if (!is(x, "dgCMatrix")) {
    stop("Not a dgCMatrix object of Matrix package")
  }
  Adj_Matrix_Norm <- BiocGenerics::t((BiocGenerics::t(x))/(pmax(Matrix::colSums(x, na.rm = FALSE, 
                                             dims = 1, sparseResult = FALSE),1)))
  return(Adj_Matrix_Norm)
}


Ara_PPI <- str_split(unique(apply(Ara_PPI, 1, function(x) paste0(sort(x)[1],"_",sort(x)[2]))),"_",simplify = T)
Ara_PPI_graph <- graph.data.frame(Ara_PPI, directed = FALSE)


PPI_MultiplexObject <- create.multiplex(Ara_PPI_graph,Layers_Name=c("PPI"))


AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)
AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)

SeedGene <- intersect(Genefile[,1],Ara_PPI[,1])
# print(SeedGene)
# print(AdjMatrix_PPI)
# print(PPI_MultiplexObject)
#  SeedGene <-sample(Ara_PPI[,1],100)
  #SeedGene <- c("AT1G01040","AT1G09700")
  ## We launch the algorithm with the default parameters (See details on manual)
  RWR_PPI_Results <- Random.Walk.Restart.Multiplex(AdjMatrix_PPI,r=Restart_probability,
                                                   PPI_MultiplexObject,SeedGene)
#  res[[i]] <- RWR_PPI_Results$RWRM_Results[1:10,]$NodeNames
#  length(unique(unlist(res)))

TopResults_PPI <- create.multiplexNetwork.topResults(RWR_PPI_Results,PPI_MultiplexObject,k=analysis_num)
#as.data.frame(TopResults_PPI)
write.table(RWR_PPI_Results[[1]][1:analysis_num,],file = out_list,row.names = F,quote = F,sep = '\t')

pdf(file = out_figure)
par(mar=c(0.1,0.1,0.1,0.1))
plot(TopResults_PPI, vertex.label.color="black",vertex.frame.color="#ffffff",
     vertex.size= 15, edge.curved=.3,
     vertex.color = ifelse(igraph::V(TopResults_PPI)$name %in% SeedGene,"yellow",
                           "#00CCFF"), edge.color="blue",edge.width=1)
dev.off()
