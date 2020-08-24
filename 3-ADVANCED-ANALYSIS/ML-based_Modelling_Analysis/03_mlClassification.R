library(argparse)
library(randomForest)
library(ROCR)
library(pROC)
library(ggplot2)
library(e1071)
library(rpart)
library(xgboost)
parser <- ArgumentParser()
parser$add_argument("-type", default = "train", dest = "type", help = "Analysis type")
parser$add_argument("-method", default = "randomForest", dest = "method", help = "Machine learning algorithms for building CMR predictor")
parser$add_argument("-posMat", default = NULL, dest = "posMat", help = "Feature matrix of positive samples")
parser$add_argument("-negMat", default = NULL, dest = "negMat", help = "Feature matrix of negative samples")
parser$add_argument("-candidate", default = NULL, dest = "candidate", help = "Feature matrix of candidate samples")
parser$add_argument("-clf", default = NULL, dest = "clf", help = "The model used to predict")
parser$add_argument("-k", default = 0, dest = "k", help = "The K-fold cross-validation")
parser$add_argument("-perc", default = 0.2, dest = "perc", help = "The percentage of hold-out test samples")
parser$add_argument("-cpus", default = 1, dest = "cpus", help = "The number of threads used for parallel computing")
parser$add_argument("-outScore", default = NULL, dest = "outScore", help = "The probability of being a CMR")
parser$add_argument("-outModel", default = NULL, dest = "outModel", help = "The directory of model")
parser$add_argument("-outCV", default = NULL, dest = "outCV", help = "The directory of cross validation")
parser$add_argument("-outTest", default = NULL, dest = "outTest", help = "The directory of hold-out test evaluation")

args <- parser$parse_args()
type <- args$type
method <- args$method
posDir <- args$posMat
negDir <- args$negMat
perc <- args$perc
k <- as.numeric(args$k)
cpus <- args$cpus

mainDic <- "/home/galaxy/tools/3-ADVANCED-ANALYSIS/ML-based_Modelling_Analysis/"
source(paste0(mainDic, '03_cross_validation.R'))

if(type == "train"){
    posMat <- read.table(file = posDir, sep = '\t', header = T, quote = "", stringsAsFactors = F)
    negMat <- read.table(file = negDir, sep = '\t', header = T, quote = "", stringsAsFactors = F)
    posTrain <- sample(rownames(posMat), nrow(posMat)*0.8)
    negTrain <- sample(rownames(negMat), nrow(negMat)*0.8)
    posTest <- setdiff(rownames(posMat), posTrain)
    negTest <- setdiff(rownames(negMat), negTrain)

    pos.train.mat <- posMat[posTrain, ]
    neg.train.mat <- negMat[negTrain, ]
    pos.test.mat <- posMat[posTest, ]
    neg.test.mat <- negMat[negTest, ]

    featureMat <- rbind(pos.train.mat, neg.train.mat)
    rm(pos.train.mat)
    rm(neg.train.mat)

    clf <- classifier(method = method, featureMat = featureMat, positiveSamples = posTrain, negativeSamples = negTrain)
    # predicting on test samples
    pos.test.score <- .predictor(method = method, classifier = clf, featureMat = pos.test.mat)
    neg.test.score <- .predictor(method = method, classifier = clf, featureMat = neg.test.mat)

    if(k != 0){
        cvRes <- cross_validation(seed = 1, method = method, featureMat, posTrain, negTrain, cross = k, cpus = cpus)
        posScore <- NULL
        negScore <- NULL
        for(i in 1:length(cvRes)){
            posScore <- c(posScore, cvRes[[i]]$positives.test.score) 
            negScore <- c(negScore, cvRes[[i]]$negatives.test.score)
        }
        resMat <- NULL
        for(i in 1:100){
            curThreshold <- i/100
            curRes <- evalPrediction(threshold = curThreshold, posScore = posScore,
                                    negScore = negScore)
            resMat <- rbind(resMat, curRes)
        }
        threshold <- as.numeric(which(resMat[,"F-score"] == max(resMat[,"F-score"], na.rm = T)))/100
        pdf(file = args$outCV, height = 5, width = 10)
        par(mfrow = c(1,2))
        plotROC(cvRes)
        plot((1:100)/100, resMat[,1], type = "l", col = "red",  ylab = "", lwd = 2, xlab = "threshold")
        lines((1:100)/100, resMat[,2], col = "blue", lwd = 2)
        legend("topleft", legend = c("Sn", "Sp"), col = c("red", "blue"), lwd = rep(2, 2))
        dev.off()
        
        evalRes <- evalPrediction(threshold = threshold, posScore = pos.test.score, negScore = neg.test.score)
        df <- data.frame(measures = factor(names(evalRes), levels = names(evalRes)),
                 values = evalRes)
        p <- ggplot(data = df, aes(x = measures, y = values)) +
        geom_bar(stat = "identity", fill = "steelblue")+
        theme_classic()
	ggsave(p, filename = paste0(args$outTest, ".pdf"), height = 5, width = 5)
    }else{
        threshold <- 0.5
        evalRes <- evalPrediction(threshold = 0.5, posScore = pos.test.score, negScore = neg.test.score)
        df <- data.frame(measures = factor(names(evalRes), levels = names(evalRes)),
                 values = evalRes)
        p <- ggplot(data = df, aes(x = measures, y = values)) +
        geom_bar(stat = "identity", fill = "steelblue")+
        theme_classic()
	ggsave(p, filename = paste0(args$outTest, ".pdf"), height = 5, width = 5)
    }
    save(clf, file = args$outModel)
}else if(type == "bulit_in"){
    load(paste0(mainDic, "model.RData"))
    featureMat <- read.table(file = args$candidate, sep = '\t', header = T, quote = "", stringsAsFactors = F)
    predScore <- predict(object = obj, newdata = featureMat, type = "vote")[,"1"]
    df <- data.frame(ID = names(predScore), Score = predScore)
    write.table(df, file = args$outScore, sep = "\t", quote = F, row.names = F, col.names = T)
}else{
    load(args$clf)
    featureMat <- read.table(file = args$candidate, sep = '\t', header = T, quote = "", stringsAsFactors = F)
    predScore <- .predictor(method = args$method, classifier = clf, featureMat = featureMat)
    df <- data.frame(ID = names(predScore), Score = predScore)
    write.table(df, file = args$outScore, sep = "\t", quote = F, row.names = F, col.names = T)

}
