args <- commandArgs(T)
release_version <- args[1]
# Type <- args[1] #"Ensembl Plants"
Type <- "plants"
Species <- args[2] #"Zea mays"
DataType <- args[3]  #"Genome" #cDNA, CDS, Protein, GTF, GFF3
out <- args[4]
releaseVersion <- args[]
library(threadr)
mainDic <- "/home/DeepEA/galaxy/tools/1-PRE-ANALYSIS/Data_Preparation/"

options(stringsAsFactors = FALSE)
MHmakeRandomString <- function(n=1, lenght=12)
{
  randomString <- c(1:n)                  # initialize vector
  for (i in 1:n)
  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                    lenght, replace=TRUE),
                             collapse="")
  }
  return(randomString)
}

generateURL <- function(Type, DataType, Species){
  if(Type == "plants"){
    ftp <- paste0("ftp://ftp.ensemblgenomes.org/pub/release-", release_version, "/plants/")
    species_df <- read.delim(file = paste0(mainDic, "species_EnsemblPlants.txt"),
                             sep = '\t', header = FALSE)
  }else if(Type == "fungi"){
    ftp <- "ftp://ftp.ensemblgenomes.org/pub/release-44/fungi/"
    species_df <- read.delim(file = paste0(mainDic, "species_EnsemblFungi.txt"),
                             sep = '\t', header = FALSE)
  }else if(Type == "metazoa"){
    ftp <- "ftp://ftp.ensemblgenomes.org/pub/release-44/metazoa/"
    species_df <- read.delim(file = paste0(mainDic, "species_EnsemblMetazoa.txt"),
                             sep = '\t', header = FALSE)
  }else if(Type == "protists"){
    ftp <- "ftp://ftp.ensemblgenomes.org/pub/release-44/protists/"
    species_df <- read.delim(file = paste0(mainDic, "species_EnsemblProtists.txt"),
                             sep = '\t', header = FALSE)
  }else{
    ftp <- "ftp://ftp.ensembl.org/pub/release-97/"
    species_df <- read.delim(file = paste0(mainDic, "species_EnsemblVertebrates.txt"),
                             sep = '\t', header = FALSE)
  }

  species_value <- Species
  assembly <- paste(strsplit(species_df$V3[which(species_df$V2 == Species)], " ")[[1]], collapse = "_")
  Species <- species_df$V1[which(species_df$V2 == species_value)]

  if(DataType == "Genome"){
    curFTP <- paste0(ftp, "fasta/", species_value, "/dna/")
    filenames <-  list_files_ftp(url = curFTP,
                                 credentials = "",
                                 sleep = NA,
                                 sort = FALSE,
                                 verbose = FALSE)
    url <- filenames[grep("dna_sm.toplevel.fa.gz", filenames)]
  }else if(DataType == "cDNA"){
    curFTP <- paste0(ftp, "fasta/", species_value, "/cdna/")
    filenames <-  list_files_ftp(url = curFTP,
                                 credentials = "",
                                 sleep = NA,
                                 sort = FALSE,
                                 verbose = FALSE)
    url <- filenames[grep(".cdna.all.fa.gz", filenames)]
  }else if(DataType == "CDS"){
    curFTP <- paste0(ftp, "fasta/", species_value, "/cds/")
    filenames <-  list_files_ftp(url = curFTP,
                                 credentials = "",
                                 sleep = NA,
                                 sort = FALSE,
                                 verbose = FALSE)
    url <- filenames[grep(".cds.all.fa.gz", filenames)]
  }else if(DataType == "Protein"){
    curFTP <- paste0(ftp, "fasta/", species_value, "/pep/")
    filenames <-  list_files_ftp(url = curFTP,
                                 credentials = "",
                                 sleep = NA,
                                 sort = FALSE,
                                 verbose = FALSE)
    url <- filenames[grep(".pep.all.fa.gz", filenames)]
  }else if(DataType == "GTF"){
    curFTP <- paste0(ftp, "gtf/", species_value, "/")
    filenames <- list_files_ftp(url = curFTP, credentials = "", sleep = NA, sort = FALSE, verbose = FALSE)
    url <- filenames[grep(".44.gtf.gz", filenames)]
  }else{
    curFTP <- paste0(ftp, "gff3/", species_value, "/")
    filenames <- list_files_ftp(url = curFTP, credentials = "", sleep = NA, sort = FALSE, verbose = FALSE)
    url <- filenames[grep(".44.gff3.gz", filenames)]
  }
  url
}

url <- generateURL(Type = Type, DataType = DataType, Species = Species)

dirname <- MHmakeRandomString()
dir.create(paste0(mainDic, dirname), showWarnings = FALSE)

if(DataType == "GFF3"){
  outName <- paste0(paste(strsplit(Species, " ")[[1]], collapse = "_"),
                    "_", tolower(DataType), ".gff3.gz")
}else if(DataType == "GTF"){
  outName <- paste0(paste(strsplit(Species, " ")[[1]], collapse = "_"),
                    "_", tolower(DataType), ".gtf.gz")
}else{
  outName <- paste0(paste(strsplit(Species, " ")[[1]], collapse = "_"),
                    "_", tolower(DataType), ".fasta.gz")
}
cmd <- paste0("axel -n 20 -o ", paste0(mainDic, dirname, "/"), outName, " ", url, " > tmp.txt")
system(command = cmd)
system(command = paste0("gunzip ", paste0(mainDic,dirname, "/", outName)))

outName <- gsub(".gz", "", outName)
system(command = paste0("mv ", paste0(mainDic, dirname, "/", outName), " ", out))
system(command = paste0("rm -r ", paste0(mainDic, dirname)))

