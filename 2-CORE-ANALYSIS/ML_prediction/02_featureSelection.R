options(stringsAsFactors = F, warnings = -1)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(DALEX))
suppressPackageStartupMessages(library(pheatmap))
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-posMat" , default = NULL, dest = "posMat",
                    help = "Feature matrix of positive samples.")
parser$add_argument("-negMat" , default = NULL, dest = "negMat",
                    help = "Feature matrix of negative samples.")
parser$add_argument("-cutoff", default = 0.5, dest = "cutoff",
                    help = "Pearson correlation coefficient threshold.")
parser$add_argument("-method", default = NULL, dest = "method",
                    help = "Feature selection method.")
parser$add_argument("-algorithm", default = NULL, dest = "algorithm",
                    help = "Machine learning algorithm for feature selection, only required for variable_importance")
parser$add_argument("-plot", default = NULL, dest = "plot",
                    help = "If plot")
# parser$add_argument("-topFeatures", default = 10, dest = "topFeatures",
#                     help = "")
parser$add_argument("-outFigure", default = NULL, dest = "outFigure",
                    help = "")
parser$add_argument("-outImportance", default = NULL, dest = "outImportance",
                    help = "")
parser$add_argument("-outFeatures", default = NULL, dest = "outFeatures",
                    help = "")
parser$add_argument("-featureNum", default = NULL, dest = "featureNum",
                    help = "")
parser$add_argument("-evaluate", default = "repeatedcv", dest = "evaluate",
                    help = "")
parser$add_argument("-measure", default = NULL, dest = "measure", 
                    help = "")
parser$add_argument("-parallel", default = NULL, dest = "parallel", 
                    help = "")
parser$add_argument("-standardize", default = NULL, dest = "standardize", 
                    help = "")
parser$add_argument("-repeats", default = NULL, type="integer", dest = "repeats", 
                    help = "")
parser$add_argument("-iters", default = NULL, type="integer", dest = "iters", 
                    help = "")
parser$add_argument("-ntree", default = 100, type="integer", dest = "ntree", 
                    help = "")

args <- parser$parse_args()

posMat <- read.table(file = args$posMat, sep = '\t', header = T,
                     quote = "")
negMat <- read.table(file = args$negMat, sep = '\t', header = T,
                     quote = "")

posMat$label <- 1 #assign 1 for positives
negMat$label <- 0 # assign 0 for negatives
featureMat <- rbind(posMat, negMat)

if(args$method == "remove_redundant"){
  correlationMatrix <- cor(featureMat[,1:(ncol(featureMat)-1)])
  if(args$plot){
	  pheatmap::pheatmap(correlationMatrix, filename = paste0(args$outFigure, ".pdf"), height = 6, width = 6, show_rownames = F, show_colnames = F)
  }
  highlyCorrelated <- findCorrelation(correlationMatrix,
                                      cutoff = as.numeric(args$cutoff),
                                      names = FALSE)
  resMat <- featureMat[,setdiff(1:ncol(featureMat), highlyCorrelated)]
}else if(args$method == "variable_importance"){
  Mod <- caret::train(x = featureMat[,1:(ncol(featureMat)-1)],
                      y = factor(featureMat$label),
                      method = as.character(args$algorithm))
  Imp <- varImp(Mod)
  if(args$plot == TRUE){
    pdf(args$outFigure, 5, 5)
    plot(Imp, top = 10, main = 'Variable Importance')
    dev.off()
  }
  importMat <- Imp$importance
  importMat <- importMat[order(importMat$Overall, decreasing = T), , drop = FALSE]
  write.table(importMat, file = args$outImportance, sep = '\t',
              quote = F, row.names = T, col.names = F)
  idx <- rownames(importMat)[1:as.numeric(args$featureNum)]
  resMat <- featureMat[,idx]
}else if(args$method == "RFE"){
  ctrl <- rfeControl(functions = rfFuncs,
                     method = args$evaluate,
                     verbose = FALSE)
  Profile <- rfe(x = featureMat[,1:(ncol(featureMat)-1)],
                 y = factor(featureMat$label),
                 sizes = as.numeric(args$featureNum),
                 rfeControl = ctrl)
  resMat <- featureMat[, Profile$optVariables[1:as.numeric(args$featureNum)]]
}else if(args$method == "lasso_regression"){
    if(args$parallel == "yes"){
      parallel <- TRUE
    }else{
      parallel <- FALSE
    }
    if(args$standardize == "yes"){
      standardize <- TRUE
    }else{
      standardize <- FALSE
    }
    cv.lasso <- cv.glmnet(x = as.matrix(featureMat[,1:(ncol(featureMat)-1)]),
                          y = as.double(featureMat$label),
                          family = 'binomial',
                          alpha = 1,
                          parallel = parallel,
                          standardize = standardize,
                          type.measure = args$measure)
    df_coef <- round(as.matrix(coef(cv.lasso, s=cv.lasso$lambda.min)), 2)
    write.table(df_coef, file = args$outImportance, sep = '\t',
                quote = F, row.names = T, col.names = F)
    # See all contributing variables
    resFeatures <- names(df_coef[df_coef[, 1] != 0, ])
    resMat <- featureMat[,resFeatures]
}else if(args$method == "genetic_algorithm"){
    # Define control function
    ga_ctrl <- gafsControl(functions = rfGA,  # another option is `caretGA`.
                            method = args$evaluate,
                            repeats = args$repeats)
    # Genetic Algorithm feature selection
    set.seed(100)
    ga_obj <- gafs(x = featureMat[,1:(ncol(featureMat)-1)],
                   y = factor(featureMat$label),
                   iters = args$iters,   # normally much higher (100+)
                   gafsControl = ga_ctrl)
    resMat <- featureMat[, ga_obj$optVariables]
}else if(args$method == "simulated_annealing"){
    # Define control function
    sa_ctrl <- safsControl(functions = rfSA,
                            method = args$evaluate,
                            repeats = as.numeric(args$repeats),
                            improve = as.numeric(args$iters)) # n iterations without improvement before a reset

    # Genetic Algorithm feature selection
    set.seed(100)
    sa_obj <- safs(x = featureMat[,1:(ncol(featureMat)-1)],
                   y = factor(featureMat$label),
                   safsControl = sa_ctrl)
    resMat <- featureMat[, sa_obj$optVariables]
}else{
    rf_mod <- randomForest(x = featureMat[,1:(ncol(featureMat)-1)],
                           y = factor(featureMat$label),
                           ntree = as.numeric(args$ntree))
    explained_rf <- explain(rf_mod, data = featureMat[,1:(ncol(featureMat)-1)], y = featureMat$label)
    varimps <- variable_dropout(explained_rf, type='raw')
    write.table(varimps, file = args$outImportance, sep = '\t',
                quote = F, row.names = T, col.names = F)
    varimps <- varimps[which(varimps$variable != "_full_model_" & varimps$variable != "_baseline_"), ]
    resFeatures <- varimps[order(varimps$dropout_loss, decreasing = T), 1]
    resMat <- featureMat[,resFeatures]
}

write.table(resMat, file = args$outFeatures, sep = "\t", quote = F, row.names = F, col.names = T)
