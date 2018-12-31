
require("seqinr")
#aaa <- read.fasta(file = "~/Desktop/test",as.string = T)
pfm <- function(sequences,RNA=F){
  require(seqLogo)
  sequences <- lapply(sequences, function(x) s2c(x))
  ##### calculate positon frequence matrix
  sequences_pfw <- list()
  sequences_row <- length(lengths(sequences))
  sequences_col <- length(sequences[[1]])
  pfw <- matrix(,4,sequences_col)
  if (RNA) {
    rownames(pfw) <- c("a","c","g","u")
    for (i in 1:sequences_col){
      pfw[,i] <- table(sapply(sequences, function(x) x[i]))[c("a","c","g","u")]/sequences_row
    }
  }else{
  rownames(pfw) <- c("a","c","g","t")
  for (i in 1:sequences_col){
    pfw[,i] <- table(sapply(sequences, function(x) x[i]))[c("a","c","g","t")]/sequences_row
  }
  }
  pfw[which(is.na(pfw))] <- 0
 # pfw <- apply(pfw,2,function(x) x/sum(x))[,start:end]
  pfw
}
logoplot <- function(sequences,type="weblogo",start = 1,end = NULL,...){
  if(is.null(end)){
    end=nchar(sequences[[1]])
  }
  
  pfm1 <- function(sequences,RNA=F){
    require(seqLogo)
    sequences <- lapply(sequences, function(x) s2c(x))
    ##### calculate positon frequence matrix
    sequences_pfw <- list()
    sequences_row <- length(lengths(sequences))
    sequences_col <- length(sequences[[1]])
    pfw <- matrix(0,4,sequences_col)
    # rownames(pfw) <- c("a","c","g","t")
    # for (i in 1:sequences_col){
    #   pfw[,i] <- table(sapply(sequences, function(x) x[i]))[c("a","c","g","t")]/sequences_row
    # }
    
    if (RNA) {
      rownames(pfw) <- c("a","c","g","u")
      for (i in 1:sequences_col){
        pfw[,i] <- table(sapply(sequences, function(x) x[i]))[c("a","c","g","u")]/sequences_row
      }
    }else{
      rownames(pfw) <- c("a","c","g","t")
      for (i in 1:sequences_col){
        pfw[,i] <- table(sapply(sequences, function(x) x[i]))[c("a","c","g","t")]/sequences_row
      }
    }
    
    pfw[which(is.na(pfw))] <- 0
    pfw <- apply(pfw,2,function(x) x/sum(x))[,start:end]
    pfw
  }
  
  weblogo <- function(pfm,...){
    ##### calculate positon weight matrix 
    require(seqLogo)
    pfw <- pfm1(sequences=sequences,...)
    pwm <- makePWM(pfw)
    #####  plot weblogo figure 
    seqLogo(pwm)
  }
  
  #####  plot weblogo bar figure
  barlogo <- function(pfm){
    library(ggplot2)
    freqs<-pfm[,start:end]
    rownames(freqs) <- c("A","C","G","T")
    freqdf <- as.data.frame(t(freqs))
    
    freqdf$pos = as.numeric(as.character(rownames(freqdf)))
    
    freqdf$height <- apply(freqdf[,c('A', 'C','G','T')], MARGIN=1,
                           FUN=function(x){2--sum(log(x^x,base=2))})
    
    logodf <- data.frame(A=freqdf$A*freqdf$height, C=freqdf$C*freqdf$height,
                         G=freqdf$G*freqdf$height, T=freqdf$T*freqdf$height, 
                         position=freqdf$pos)
    
    lmf <- melt(logodf, id.var='position')
    names(lmf)[3] <- "Sequence_Conservation"
    p <- ggplot(data=lmf, aes(x=position, y=Sequence_Conservation))  +
      geom_bar(aes(fill=variable,order=Sequence_Conservation), position='stack', 
               stat='identity', alpha=0.5) + 
      scale_fill_manual(values = c("green","red","yellow","blue"), limits = c('A','T','G','C')) + 
      theme_bw()

    return(p)
  }
  #####  export plot
  pfminfo <- pfm1(sequences,...)
  switch (type,
          weblogo = weblogo(pfm = pfminfo,...),
          barlogo = barlogo(pfm = pfminfo,...)
  )
  
  
}
twologoplot <- function(pos_sequences,neg_sequences,start = 1,end = nchar(pos_sequences[1]),stackHeight=sumOfAbsICDifferences,...){
  pfm <- function(sequences){
    require(seqLogo)
    sequences <- lapply(sequences, function(x) s2c(x))
    ##### calculate positon frequence matrix
    sequences_pfw <- list()
    sequences_row <- length(lengths(sequences))
    sequences_col <- length(sequences[[1]])
    pfw <- matrix(,4,sequences_col)
    rownames(pfw) <- c("a","c","g","t")
    for (i in 1:sequences_col){
      pfw[,i] <- table(sapply(sequences, function(x) x[i]))[c("a","c","g","t")]/sequences_row
    }
    pfw[which(is.na(pfw))] <- 0
    pfw <- apply(pfw,2,function(x) x/sum(x))[,start:end]
    pfw
  }
  
  ########
  require(seqLogo)
  PWM1 <- makePWM(pfm(sequences = pos_sequences))
  PWM2 <- makePWM(pfm(sequences = neg_sequences))
  
  require(DiffLogo)
  diffLogoFromPwm(pwm1 = PWM1, pwm2 = PWM2,stackHeight=stackHeight, ...)
}

PosNegFrequence_compare <- function(pos_seq,neg_seq){
  
  pfm <- function(sequences){
    require(seqLogo)
    require(stringr)
    require(seqinr)
    sequences <- lapply(sequences, function(x) s2c(x))
    ##### calculate positon frequence matrix
    sequences_pfw <- list()
    sequences_row <- length(lengths(sequences))
    sequences_col <- length(sequences[[1]])
    pfw <- matrix(,4,sequences_col)
    colnames(pfw) <- c(-20:20)
    rownames(pfw) <- c("a","c","g","t")
    for (i in 1:sequences_col){
      pfw[,i] <- table(sapply(sequences, function(x) x[i]))[c("a","c","g","t")]/sequences_row
    }
    pfw[which(is.na(pfw))] <- 0
    #pfw <- apply(pfw,2,function(x) x/sum(x))[,start:end]
    pfw
  }
  require(ggplot2)
  require(reshape2)
  pos_seq_pfm <- pfm(pos_seq)
  neg_seq_pfm <- pfm(neg_seq)
  
  
  data1 <- as.data.frame(pos_seq_pfm)
  data1_plot <- melt(data1)
  data1_plot[["base"]] <- rep(c("A","C","G","T"),41)
  names(data1_plot) <- c("Position","Frequency","base")
  
  data2 <- as.data.frame(-neg_seq_pfm)
  data2_plot <- melt(data2)
  data2_plot[["base"]] <- rep(c("A","C","G","T"),41)
  
  x <- 1:41
  label <- -20:20
  axis_x <- data.frame(x=x,label=label,base="A")
  
  p <- ggplot(data = data1_plot, mapping = aes(x = Position, y = Frequency, fill = base)) + 
    geom_bar(stat= 'identity', position = 'fill')  + ylim(-1,1) + 
    geom_bar(data = data2_plot, mapping = aes(x = factor(variable), y = value, fill = base),stat= 'identity', position = 'fill')+
    geom_hline(aes(yintercept = 0),color = "white",lwd=5) 
  #geom_text(data = data2_plot, mapping = aes(x = factor(variable),y = value, fill = base),
  #          label=paste(round(data2_plot$value,2)*100,'%',sep = ''),
  #          size=2.8, colour = "black",position=position_stack(.2), vjust=1.5)
  
  p + geom_text(data=axis_x,aes(x=x,y=0,label=label)) +
    theme( axis.text.x  = element_text(size=0) ,
           axis.ticks = element_blank(),
           panel.grid.minor=element_blank(),
           panel.grid.major=element_blank())
}

#logoplot(sequences = a1,type = "barlogo")
#twologoplot(pos_sequences = aaa[1:50],aaa[50:100],sparse=T)

