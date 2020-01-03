library(seqLogo)
library(DiffLogo)
library(seqinr)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-pos", default = NULL, dest = "pos_seq", help = "positive sequences")
parser$add_argument("-neg", default = NULL, dest = "neg_seq", help = "negative sequences")
parser$add_argument("-s", default = 1, dest = "start", help = "start position")
parser$add_argument("-e", default = NULL, dest = "end", help = "end position")
parser$add_argument("-o", default = NULL, dest = "out", help = "output directory")
parser$add_argument("-twosample", default = "TRUE", dest = "twosample", help = "two samples or single sample")
pfm <- function(sequences, start, end){
  sequences <- lapply(sequences, s2c)
  ##### calculate positon frequence matrix
  sequences_pfw <- list()
  sequences_row <- length(lengths(sequences))
  sequences_col <- length(sequences[[1]])
  pfw <- matrix(NA, 4, sequences_col)
  rownames(pfw) <- c("a","c","g","t")
  for (i in 1:sequences_col){
    pfw[,i] <- table(sapply(sequences, function(x) x[i]))[c("a","c","g","t")]/sequences_row
  }
  pfw[which(is.na(pfw))] <- 0
  pfw <- apply(pfw,2,function(x) x/sum(x))[,start:end]
  pfw
}

twologoplot <- function(pos_sequences, neg_sequences, start = 1,
                        end = as.numeric(nchar(pos_sequences[1])),
                        stackHeight=sumOfAbsICDifferences,...){

  PWM1 <- makePWM(pfm(sequences = pos_sequences, start = start, end = end))
  PWM2 <- makePWM(pfm(sequences = neg_sequences, start = start, end = end))

  diffLogoFromPwm(pwm1 = PWM1, pwm2 = PWM2,stackHeight=stackHeight, ...)
}


logoplot <- function(sequence, start = 1, end = as.numeric(nchar(sequence[1]))){
  pwm <- makePWM(pfm(sequence, start = start, end = end))
  seqLogo(pwm = pwm)
}


args <- parser$parse_args()
posSeq <- read.fasta(file = args$pos_seq, as.string = T)
pdf(file = args$out, height = 5, width = 10)
if(args$twosample){
  negSeq <- read.fasta(file = args$neg_seq, as.string = T)
  if(is.null(args$end)){
    twologoplot(pos_sequences = posSeq, neg_sequences = negSeq, start = args$start)
  }else{
    twologoplot(pos_sequences = posSeq, neg_sequences = negSeq, start = args$start, end = args$end)
  }
}else{
  if(is.null(args$end)){
    logoplot(sequence = posSeq)
  }else{
    logoplot(sequence = posSeq, start = args$start, end = args$end)
  }
}
dev.off()
