# functions required in multi-omics integration analysis module
# define mode of points in ternary plot
defineMod <- function(inputVec){
  cmr <- as.numeric(inputVec[1])/sum(inputVec)*100 # cmr
  trans <- as.numeric(inputVec[3])/sum(inputVec)*100 # translation
  exp <- as.numeric(inputVec[2])/sum(inputVec)*100 # expression
  if(cmr >= 20 & cmr <= 60 & trans >= 20 & trans <= 60 & exp >= 20 & exp <= 60){
    mod <- "balanced"
  }else if(exp < 20 & trans > 60 & cmr < 20){
    mod <- "Translation dominant"
  }else if(trans > 20 & trans < 80 & exp < 20 & cmr > 20 & cmr < 80){
    mod <- "Expression suppressed"
  }else if(trans > 20 & trans < 80 & exp > 20 & exp < 80 & cmr< 20){
    mod <- "m6A suppressed"
  }else if(trans < 20 & cmr > 60 & exp < 20){
    mod <- "m6A dominant"
  }else if(trans < 20 & cmr > 20 & cmr < 80 & exp > 20 & exp < 80){
    mod <- "Translation suppressed"
  }else if(trans < 20 & cmr < 20 & exp > 60){
    mod <- "Expression dominant"
  }else{
    mod <- "others"
  }
  mod 
}

# A hyperlink for GO ID to ebi
id2link <- function(ID){
  res <- paste0("<a href=", "'https://www.ebi.ac.uk/QuickGO/term/",
                ID, "/' target='_blank'>", ID,
                "</a>")
  res
}

# permutation test
permuTest <- function(B = 1000, gene_all, N, targetGene, targetNumber){
  res <- NULL
  for(i in 1:B){
    curGene <- sample(gene_all, N)
    interGene <- intersect(curGene, targetGene)
    res <- c(res, length(interGene))
  }
  p <- length(which(res >= targetNumber))/B
  p
}

label <- function(txt) {
  list(
    text = txt, 
    x = 0.1, y = 1,
    ax = 0, ay = 0,
    xref = "paper", yref = "paper", 
    align = "center",
    font = list(family = "serif", size = 15, color = "white"),
    bgcolor = "#b3b3b3", bordercolor = "black", borderwidth = 2
  )
}

# reusable function for axis formatting
axis <- function(txt) {
  list(
    title = txt, tickformat = ".0%", tickfont = list(size = 10)
  )
}

# ternary plot each type of CMR genes
plotTernary <- function(cmrMat, geneType = "Singletons", sub = TRUE, subgenome = subgenome){
	df <- subset(cmrMat, cmrMat$Type == geneType)
	df <- na.omit(df)
	if(sub){
		df$subgenome <- subgenome$Subgenome[match(df$geneID, subgenome$B73.v4.ID)]
		df$subgenome <- subgenome$Subgenome[match(df$geneID, subgenome$B73.v4.ID)]
	}

	tt <- as.matrix(df[,2:4])
	# tt <- tt/100
	res <- apply(tt, 1, defineMod)
	df$mod <- res
	p <- plot_ly(
	df, 
	a = ~cmrLevel, 
	b = ~trans.level, 
	c = ~exp.level, 
	color = ~mod, 
	type = "scatterternary",
	text = ~geneID
) %>% 
  layout(
    ternary = list(
      aaxis = axis("CMR"), 
      baxis = axis("Translation"), 
      caxis = axis("Expression")
    ),
    title = geneType
  )
  resList <- list(df = df, p = p)
  return(resList)
}

# ternary plot for functional genes
plotTernarySubset <- function(df, geneType = "Singletons",  domestication = domestication, improvement = improvement, housekeeping = housekeeping, TF = TF) {
	df.subset <- df %>% mutate(domestication = geneID %in% domestication$AGPv4,
						improvement = geneID %in% improvement$AGPv4,
						housekeeping = geneID %in% housekeeping$AGPv4,
						TF = geneID %in% TF$Gene_ID) %>% 
		filter((domestication + improvement + housekeeping + TF)>=1) %>% 
		mutate(domestication = ifelse(domestication, "domestication", NA),
			improvement = ifelse(improvement, "improvement", NA),
			housekeeping = ifelse(housekeeping, "housekeeping", NA),
			TF = ifelse(TF, "TF", NA)) 
	function_type <-  df.subset %>% rowwise() %>% do(a = paste(na.omit(c(.$domestication, .$improvement, .$housekeeping, .$TF)), collapse = ";"))
	df.subset$function_type <- unlist(function_type$a)
	df.subset$text <- paste0(df.subset$geneID, "\n", df.subset$function_type)
	# Initiating a plotly visualization 
	p <- plot_ly(
		df.subset, 
		a = ~cmrLevel, 
		b = ~trans.level, 
		c = ~exp.level, 
		color = ~function_type,
		type = "scatterternary",
		text = ~text
	) %>% 
		layout(
		ternary = list(
			aaxis = axis("CMR"), 
			baxis = axis("Translation"), 
			caxis = axis("Expression")
		),
		title = geneType
		)
	return(p)
}

GOTable <- function(resGO, ontology = "BP", geneType = "Singletons", pageLength = 5){
	dfGO <- resGO %>% filter(Ontology == ontology)
	dfMat <- matrix(1, nrow = length(unique(dfGO$GO.ID)), ncol = length(unique(dfGO$Type)))
	rownames(dfMat) <- unique(dfGO$GO.ID)
	colnames(dfMat) <- unique(dfGO$Type)
	for (i in unique(dfGO$Type)) {
		curID <- dfGO$GO.ID[which(dfGO$Type == i)]
		curP <- dfGO$elimFisher[which(dfGO$Type == i)]
		dfMat[curID, i] <- formatC(as.numeric(curP), format = "e", digits = 2)
	}
	dfMat <- data.frame(dfMat)
	dfMat$Term <- dfGO$Term[match(rownames(dfMat), dfGO$GO.ID)]
	dfMat <- dfMat[,c(ncol(dfMat),1:(ncol(dfMat)-1))]

	ID <- rownames(dfMat)
	dfMat$ID <- sapply(rownames(dfMat), id2link)
	dfMat <- dfMat[,c(ncol(dfMat),1:(ncol(dfMat)-1))]
	rownames(dfMat) <- 1:nrow(dfMat)
	new_dfmat <- dfMat

	filter_para_col <- colnames(new_dfmat)[-1]
	filter_para_col <- colnames(new_dfmat)[-1]

	set_list <- lapply(1:length(filter_para_col),function(cc){
		formatter("span",
				style = x ~ style(color = ifelse(as.numeric(x) < 0.05, "red", "green")),
				x ~ icontext(ifelse(as.numeric(x) < 0.05, "ok", "remove"),x) )
	})
	names(set_list) <- colnames(new_dfmat)[-1]

	tsis_event_ifTable <- formattable(new_dfmat,set_list)
	res <- as.datatable(tsis_event_ifTable,extensions = 'Buttons',
				options =   list(dom = 'Bfrtip',
		pageLength = pageLength,
		buttons = list(c('copy', 'csv', 'excel', 'pdf')),
		searchHighlight = TRUE,
		colReorder = TRUE,
		scrollX = TRUE,
		fixedColumns = TRUE,
		extensions = 'Responsive',
		deferRender = TRUE,
		scroller = TRUE,
		lengthChange = FALSE
		))
	resList <- list(dfMat = dfMat, Table = res, ID = ID)
	return(resList)
}

functionalGene <- function(df, domestication = domestication, improvement = improvement, housekeeping = housekeeping, TF = TF){
  stat <- df %>% select(geneID, mod) %>% group_by(mod) %>% 
  	summarise(Number = n(),
            domestication = length(intersect(geneID, domestication$AGPv4)),
            improvement = length(intersect(geneID, improvement$AGPv4)),
            housekeeping = length(intersect(geneID, housekeeping$AGPv4)),
            TF = length(intersect(geneID, TF$Gene_ID)))
	colnames(stat)[1] <- "Mode"
	# for domestication
	p.domestication <- NULL
	for (i in 1:nrow(stat)) {
		curNumber <- stat$domestication[i]
		if(curNumber == 0){
			p.value <- 1
		}else{
			p.value <- permuTest(B = 1000, gene_all = all_genes, N = stat$Number[i],
								targetGene = domestication$AGPv4, targetNumber = curNumber)
		}
		p.domestication <- c(p.domestication, p.value)
	}

	# for improvement
	p.improvement <- NULL
	for (i in 1:nrow(stat)) {
		curNumber <- stat$improvement[i]
		if(curNumber == 0){
			p.value <- 1
		}else{
			p.value <- permuTest(B = 1000, gene_all = all_genes, N = stat$Number[i],
								targetGene = improvement$AGPv4, targetNumber = curNumber)
		}
		p.improvement <- c(p.improvement, p.value)
	}

	# for housekeeping
	p.housekeeping <- NULL
	for (i in 1:nrow(stat)) {
		curNumber <- stat$housekeeping[i]
		if(curNumber == 0){
			p.value <- 1
		}else{
			p.value <- permuTest(B = 1000, gene_all = all_genes, N = stat$Number[i],
								targetGene = housekeeping$AGPv4, targetNumber = curNumber)
		}
		p.housekeeping <- c(p.housekeeping, p.value)
	}

	# for TF
	p.TF <- NULL
	for (i in 1:nrow(stat)) {
		curNumber <- stat$TF[i]
		if(curNumber == 0){
			p.value <- 1
		}else{
			p.value <- permuTest(B = 1000, gene_all = all_genes, N = stat$Number[i],
								targetGene = TF$AGPv4, targetNumber = curNumber)
		}
		p.TF <- c(p.TF, p.value)
	}

	stat$domestication.pvalue <- p.domestication
	stat$improvement.pvalue <- p.improvement
	stat$housekeeping.pvalue <- p.housekeeping
	stat$TF.pvalue <- p.TF

	stat <- stat[,c(1:3,7,4,8,5,9,6,10)]

	x <-(stat %>% mutate(domestication.pvalue = cell_spec(domestication.pvalue, "html",
													color = "white",
													background = ifelse(domestication.pvalue < 0.05, "#D7261E", "#1f77b4"),
													bold = TRUE),
					improvement.pvalue = cell_spec(improvement.pvalue, "html",
													color = "white",
													background = ifelse(improvement.pvalue < 0.05, "#D7261E", "#1f77b4"),
													bold = TRUE),
					housekeeping.pvalue = cell_spec(housekeeping.pvalue, "html",
													color = "white",
													background = ifelse(housekeeping.pvalue < 0.05, "#D7261E", "#1f77b4"),
													bold = TRUE),
					TF.pvalue = cell_spec(TF.pvalue, "html",
													color = "white",
													background = ifelse(TF.pvalue < 0.05, "#D7261E", "#1f77b4"),
													bold = TRUE)) %>% 
	knitr::kable(format = "html", escape = FALSE) %>%  
	kableExtra::kable_styling(bootstrap_options = "striped",full_width = FALSE, 
								stripe_color = 'black', latex_options = "bordered") %>% 
	add_header_above(c(" ", " ", "Domestication" = 2,
						"Improvement" = 2,
						"Housekeeping" = 2,
						"TF" = 2)))
	x <- gsub(pattern = "improvement.pvalue", "p.value", x)
	x <- gsub(pattern = "housekeeping.pvalue", "p.value", x)
	x <- gsub(pattern = "domestication.pvalue", "p.value", x)
	x <- gsub(pattern = "TF.pvalue", "p.value", x)
	x <- gsub(pattern = "domestication", "Number", x)
	x <- gsub(pattern = "housekeeping", "Number", x)
	x <- gsub(pattern = "improvement", "Number", x)
	x <- gsub(pattern = "TF", "Number", x)
	return(x)
}

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
	
