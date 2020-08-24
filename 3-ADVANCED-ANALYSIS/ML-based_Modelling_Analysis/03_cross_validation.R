
#Function: get sample index for cv cross validation
.cvSampleIndex <- function( len, cross = 5, seed = 1 ) {
  
  cv <- cross
  sample_matrix <- matrix(0, nrow = len, ncol = cv)
  colnames(sample_matrix) <- paste("cv", c(1:cv), sep = "" )
  
  #random samples 
  set.seed(seed)
  index <- sample(1:len, len, replace = FALSE )
  step = floor( len/cv )
  
  start <- NULL
  end <- NULL
  train_lens <- rep(0, cv)
  for( i in c(1:cv) ) {
    start <- step*(i-1) + 1
    end <- start + step - 1
    if( i == cv ) 
      end <- len
    
    train <- index[-c(start:end)]
    test <- index[start:end]
    train_lens[i] <- length(train)
    
    sample_matrix[,i] <- c(train, test)
  }#end for i
  
  return( list( train_lens = train_lens, sample_matrix = sample_matrix))
}



classifier <- function( method = c("randomForest", "svm", "decision_tree", "Logistic_regression"), featureMat, positiveSamples, 
                        negativeSamples, ...) {
  call <- match.call()
  
  if( length(method) > 1){
    method <- method[1]
  } 
    
  
  if( is.null(rownames(featureMat) ) )
    stop("Error: no row names (i.e., sample IDs) were assigned for featureMat." )
  if( is.null(colnames(featureMat) ) )
    stop("Error: no colnames were defined for featureMat." )
  
  positiveSamples <- intersect( rownames(featureMat), positiveSamples )
  negativeSamples <- intersect( rownames(featureMat), negativeSamples )
  posLen <- length(positiveSamples)
  negLen <- length(negativeSamples)
  if( posLen == 0 )
    stop("Error: no positive samples included in featureMat." )
  if( negLen == 0 )
    stop("Error: no negative samples were included in featureMat." )
  
  label <- c( rep(1, posLen), rep(0, negLen) )
  fmat <- data.frame( featureMat[c(positiveSamples, negativeSamples), ] )
  tmpData <- cbind( fmat, label )
  colnames(tmpData) <- c(colnames(fmat), "Class")
  if( method == "randomForest" ) {
    obj <- randomForest(x = fmat, y = factor(label), ... )
  }else if(method == "svm"){
    obj <- svm(x = fmat, y = factor(label), ... )
  }else if(method == "decision_tree"){
    obj <- rpart(Class~., data = tmpData, method = 'class')
  }else if(method == "Logistic_regression"){
    obj <- glm(Class ~.,family = binomial(link = 'logit'), data = tmpData)
  }else{
    obj <- xgboost(data = as.matrix(featureMat[c(positiveSamples, negativeSamples), ]),
                     label = label,
                     nthread = 1, nrounds = 1000,
                     objective = "binary:logistic")
  }
  obj
}


.predictor <- function( method = c("randomForest", "svm", "decision_tree", "Logistic_regression"), classifier, featureMat ) {
  
  if(length(method) > 1){
    method <- method[1]
  }
  
  if( method == "randomForest") {
    res <- predict(classifier, data.frame(featureMat), type= "vote" )[,"1"]
  }else if (method == "svm") {
    res <- predict( classifier, data.frame(featureMat), type = "raw") 
  }else if(method == "decision_tree"){
    res <- predict( classifier, data.frame(featureMat))[,"1"] 
  }else if(method == "Logistic_regression"){
    res <- predict( classifier, data.frame(featureMat), type = "response")
  }else{
    res <- predict(classifier, as.matrix(featureMat))
  }
  names(res) <- rownames(featureMat)
  res
}

##get system time for seed and then generate random index
.randomSeed <- function() {
  curtime <- format(Sys.time(), "%H:%M:%OS4")
  XXX <- unlist(strsplit(curtime, ":"))
  curtimeidx <- (as.numeric(XXX[1])*3600 + as.numeric(XXX[2])*60 + as.numeric(XXX[3]))*10000
  curtimeidx
}



.one_cross_validation <- function( cv, method, featureMat, positives, negatives, posSample_cv, negSample_cv, balanced = TRUE, ratio = 10, ... ) {
  call <- match.call()
  j <- cv
  
  #for train samples
  train_genes_p <- positives[ (posSample_cv$sample_matrix[,j][1:posSample_cv$train_lens[j]] ) ]
  test_genes_p <- positives[ (posSample_cv$sample_matrix[,j][-c(1:posSample_cv$train_lens[j])]) ]
  
  #trained negatives randomly selected, and tested on all negatives
  train_genes_n <- negatives[(negSample_cv$sample_matrix[,j][1:negSample_cv$train_lens[j]] ) ]
  test_genes_n <- negatives[ (negSample_cv$sample_matrix[,j][-c(1:negSample_cv$train_lens[j])]) ]
  
  #select part of train_genes_n
  if( balanced == TRUE ) {
    if( length(train_genes_n) > ratio*length(train_genes_p) ) {
      train_genes_n <- train_genes_n[sample(1:length(train_genes_n), replace=FALSE)[1:(ratio*length(train_genes_p))]]
    }
  }
  
  
  
  obj <- classifier( method = method, featureMat = featureMat, positiveSamples = train_genes_p, negativeSamples = train_genes_n, ... )
  bestmodel <- obj
  
  positives.train.score <- .predictor( method = method, classifier = bestmodel, featureMat = featureMat[train_genes_p,])
  negatives.train.score <- .predictor( method = method, classifier = bestmodel, featureMat = featureMat[train_genes_n,])
  positives.test.score <- .predictor( method = method, classifier = bestmodel, featureMat = featureMat[test_genes_p,])
  negatives.test.score <- .predictor( method = method, classifier = bestmodel, featureMat = featureMat[test_genes_n,])
  
  
  
  train.AUC <- roc( c(rep(1, length(train_genes_p)), rep(0, length(train_genes_n))), 
                    c(positives.train.score, negatives.train.score) )$auc[1]
  test.AUC <- roc( c(rep(1, length(test_genes_p)), rep(0, length(test_genes_n))), 
                   c(positives.test.score, negatives.test.score) )$auc[1]
  
  res <- ( list( positves.train = train_genes_p, negatives.train = train_genes_n, 
                 positives.test = test_genes_p, negatives.test = test_genes_n, 
                 ml = method, classifier = bestmodel, 
                 positives.train.score = positives.train.score,
                 negatives.train.score = negatives.train.score,
                 positives.test.score = positives.test.score,
                 negatives.test.score = negatives.test.score,
                 train.AUC = train.AUC,
                 test.AUC = test.AUC) )
  
  res
}


cross_validation <- function( seed = 1, method = c("randomForest", "svm"), 
                              featureMat, positives, negatives, cross = 5, 
                              cpus = 1, ... ){
  
  call <- match.call()
  
  if( length(method) > 1){
    method <- method[1]
  } 
  
  #sample index for cv
  posSample_cv <- .cvSampleIndex(length(positives), cross = cross, seed = seed)
  negSample_cv <- .cvSampleIndex(length(negatives), cross = cross, seed = seed)
  
  cvRes <- list()
  if( cpus > 1 ) {
    #require(snowfall)
    sfInit(parallel = TRUE, cpus = cpus)
    sfExport("classifier")
    sfExport(".predictor")
    sfExport(".one_cross_validation")
    sfLibrary( "pROC", character.only = TRUE)
    sfLibrary( "e1071", character.only = TRUE)
    sfLibrary( "randomForest", character.only = TRUE )
    
    cvRes <- sfApply( matrix(1:cross, ncol = 1), 1,  .one_cross_validation, method = method, featureMat = featureMat, positives = positives, negatives = negatives, posSample_cv = posSample_cv, negSample_cv = negSample_cv, ...)
    sfStop()
  }else {
    for( j in 1:cross ) {
      cvRes[[j]] <- .one_cross_validation( cv = j, method = method, featureMat = featureMat, positives = positives, negatives = negatives, posSample_cv = posSample_cv, negSample_cv = negSample_cv, ... )
    }
  }
  cvRes
}

evalPrediction <- function(threshold,
                           posScore = NULL, negScore = NULL,
                           beta = 1,
                           TP, TN, FP, FN){
  
  if(is.null(posScore) & is.null(negScore)){
    Sn <- TP/(TP+FN)
    Sp <- TN/(TN+FP)
    Pr <- TP/(TP+FP)
    Acc <- (TP+TN)/(TP+TN+FP+FN)
    
    Fscore <- ((1+beta^2)*Pr*Sn)/(beta^2*Pr+Sn)
    MCC <- (TP*TN-FP*FN)/sqrt(((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
    
    res <- c(Sn, Sp, Pr, Acc, Fscore, MCC)
    names(res) <- c("Sn", "Sp", "Pr", "Acc", "F-score", "MCC")
  }else{
    TP <- as.numeric(length(which(posScore >= threshold)))
    FN <- as.numeric(length(which(posScore < threshold)))
    TN <- as.numeric(length(which(negScore < threshold)))
    FP <- as.numeric(length(which(negScore >= threshold)))
    
    Sn <- TP/(TP+FN)
    Sp <- TN/(TN+FP)
    Pr <- TP/(TP+FP)
    Acc <- (TP+TN)/(TP+TN+FP+FN)
    
    Fscore <- ((1+beta^2)*Pr*Sn)/(beta^2*Pr+Sn)
    MCC <- (TP*TN-FP*FN)/sqrt(((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
    AUC <- pROC::roc( c(rep(1, length(posScore)), rep(0, length(negScore))),
                      c(posScore, negScore) )$auc[1]
    res <- c(Sn, Sp, Pr, Acc, Fscore, MCC, AUC)
    names(res) <- c("Sn", "Sp", "Pr", "Acc", "F-score", "MCC", "AUC")
  }
  res
}

threshold_identify <- function(cvRes){
  posScore <- NULL
  negScore <- NULL
  for(i in 1:length(cvRes)){
    posScore <- c(posScore, cvRes[[i]]$positives.test.score) 
    negScore <- c(negScore, cvRes[[i]]$negatives.test.score)
  }
  res <- NULL
  for(i in 1:100){
    curThreshold <- i/100
    curRes <- evalPrediction(threshold = curThreshold, posScore = posScore,
                             negScore = negScore)["F-score"]
    res <- c(res, curRes)
  }
  threshold <- as.numeric(which(res == max(res, na.rm = T)))/100
  threshold
}



plotROC <- function(cvRes) {
  
  
  cvListPredictions <- list()
  cvListLabels <- list()
  AUCVec <- rep(0, length(cvRes) )
  for( i in 1:length(cvRes) ) {
    curCV <- cvRes[[i]]
    cvListPredictions[[i]] <- c( curCV$positives.test.score, curCV$negatives.test.score )
    cvListLabels[[i]] <- c( rep(1, length(curCV$positives.test.score)), rep(0, length(curCV$negatives.test.score) ) )
    AUCVec[i] <- curCV$test.AUC
  }
  mAUC <- format( mean(AUCVec), digits= 3)
  
  #if( !require(ROCR) ) {
  #   install.packages("ROCR")
  #   library(ROCR)
  #}
  pred <- prediction(cvListPredictions, cvListLabels)
  perf <- performance(pred,"tpr","fpr")
  
  
  par(mar=c(5,6,4,2))   
  plot(perf, col= "gray", lty=3, main = paste( "AUC = ", mAUC, sep = ""), cex.lab = 2.5, cex.axis = 2, cex.main = 3, mgp = c(4,1.8,0) )
  plot(perf, col = "black",  lwd= 3, avg="vertical", spread.estimate="none", add=TRUE)  
  
}
