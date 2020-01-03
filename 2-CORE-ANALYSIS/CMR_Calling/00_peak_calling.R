library(seqinr)
library(argparse)

parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add a help option
parser$add_argument("-ip" , default = NULL, dest = "IPBAM", help = "The IP BAM")
parser$add_argument("-input" , default = NULL, dest = "inputBAM", help = "The input BAM")
parser$add_argument("-gtf" , default = NULL, dest = "gtf")
parser$add_argument("-genome" , default = NULL, dest = "genome")
parser$add_argument("-m" , default = NULL, dest = "method")
parser$add_argument("-fdr" , default = 0.05, dest = "fdr")
parser$add_argument("-bin" , default = 100, dest = "bin")
parser$add_argument("-iterations" , default = 10000, dest = "iterations")
parser$add_argument("-readsCount", default = 5, dest = "readsCount")
parser$add_argument("-cpu", default = 1, dest = "cpu")
parser$add_argument("-concatenate", default = 4, dest = "concatenate")
parser$add_argument("-ratio" , default = -1, dest = "ratio")
parser$add_argument("-paired" , default = F, dest = "paired")
parser$add_argument("-out" , default = NULL, dest = "outbed")

args <- parser$parse_args()
source('/home/DeepEA/galaxy/tools/2-CORE-ANALYSIS/CMR_Calling/00_CMRCalling.R')

if(args$method=="SlidingWindow"){
  mappedInput <- system(paste0("samtools flagstat ", args$inputBAM,"  | grep 'mapped (' | cut -f 1,2,3 -d ' '"), intern = TRUE)
  mappedRIP <- system(paste0("samtools flagstat ", args$IPBAM,"  | grep 'mapped (' | cut -f 1,2,3 -d ' '"), intern = TRUE)
  eval(parse(text = paste0('mappedInput <- ',mappedInput)))
  eval(parse(text = paste0('mappedRIP <- ',mappedRIP)))
  cat("MappedRIP: ", mappedRIP, "\n!!")
  cat("MappedInput: ", mappedInput, "\n!!")
  res <- peakCalling(IPBAM = args$IPBAM, inputBAM = args$inputBAM,
                     method = "SlidingWindow",
                     refGenome = args$genome,
                     mappedInput = mappedInput,
                     mappedRIP = mappedRIP, level = args$fdr,
                     ratio = args$ratio, cpus = args$cpu,
					           readsCount = args$readsCount,
                     concatenate = args$concatenate)
}else if(args$method == "exomePeak"){
  cat("!!!", args$fdr, '\n')
  res <- peakCalling(IPBAM = args$IPBAM,
                     inputBAM = args$inputBAM,
                     GTF = args$gtf, method = "exomePeak",
                     PEAK_CUTOFF_FDR = as.numeric(args$fdr))
}else if(args$method == "MetPeak"){
  cat("!!!", args$fdr, '\n')
  res <- peakCalling(IPBAM = args$IPBAM,
                     inputBAM = args$inputBAM,
                     GTF = args$gtf, method = "MetPeak",
                     PEAK_CUTOFF_FDR = as.numeric(args$fdr))
}else if(args$method == "BayesPeak"){
  res <- peakCalling(IPBAM = args$IPBAM,
                     inputBAM = args$inputBAM,
                     GTF = args$gtf, method = "BayesPeak",
                     bin.size = as.numeric(args$bin),
                     iterations = as.numeric(args$iterations))
}

write.table(res$peaks, file = args$outbed, sep = '\t',
            quote = F, row.names = F, col.names = T)

