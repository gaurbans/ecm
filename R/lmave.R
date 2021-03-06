#'Build multiple lm models and average them
#'
#'Builds k lm models on k partitions of the data and averages their coefficients to get create one model. 
#'Each partition excludes k/nrow(data) observations. See links in the References section for further details on
#'this methodology.
#'@param formula The formula to be passed to lm
#'@param data The data to be used
#'@param k The number of models or data partitions desired
#'@param method Whether to split data by folds ("fold"), nested folds ("nestedfold"), or bootstrapping ("boot")
#'@param seed Seed for reproducibility (only needed if method is "boot")
#'@param weights Optional vector of weights to be passed to the fitting process
#'@param ... Additional arguments to be passed to the 'lm' function
#'@return an lm object
#'@details
#'In some cases--especially in some time series modeling (see ecmave function)--rather than building one model on the entire dataset, it may be preferable to build multiple models on subsets 
#'of the data and average them. The lmave function splits the data into k partitions of size (k-1)/k*nrow(data), builds k models, and then averages the coefficients of these 
#'models to get a final model. This is similar to averaging multiple tree regression models in algorithms like random forest. 
#'
#'Unlike the 'ecm' functin, this function only works with the 'lm' linear fitter. 
#'@references 
#'Jung, Y. & Hu, J. (2016). "A K-fold Averaging Cross-validation Procedure". https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5019184/
#'
#'Cochrane, C. (2018). "Time Series Nested Cross-Validation". https://towardsdatascience.com/time-series-nested-cross-validation-76adba623eb9
#'@seealso \code{lm}
#'@examples
#'##Not run
#'
#'#Build linear models to predict Wilshire 5000 index based on corporate profits, 
#'#Federal Reserve funds rate, and unemployment rate
#'data(Wilshire)
#'
#'#Build one model on the entire dataset
#'modelall <- lm(Wilshire5000 ~ ., data = Wilshire[-1])
#'
#'#Build a five fold averaged linear model on the entire dataset
#'modelave <- lmave('Wilshire5000 ~ .', data = Wilshire[-1], k = 5) 
#'
#'@export
#'@importFrom stats lm coef
lmave <- function(formula, data, k, method = 'boot', seed = 5, weights = NULL, ...){
  lmall <- lm(formula, data, ...)
  modellist <- 1:k
  
  if (method == 'fold' | method == 'nestedfold'){
    folds <- cut(seq(1, nrow(data)), breaks=k, labels=FALSE)
    if (method == 'fold'){
      operator <- "!="
    } else if (method == 'nestedfold'){
      modellist <- modellist[-1]
      operator <- "<"
    }
    trnIdxExp <- paste0("which(folds", operator, "i)")
  
  } else if (method == 'boot'){
    set.seed(seed)
    trnIdxExp <- "sample(nrow(data), (k-1) / k * nrow(data))"
  }
   
  models <- lapply(modellist, function(i) {
    trnIdx <- eval(parse(text=trnIdxExp))
    trn <- data[trnIdx, ]
    if (!is.null(weights)){
      w <- weights[trnIdx]
      lm(as.formula(formula), data = trn, weights = w, ...)
    } else if (is.null(weights)) {
      lm(as.formula(formula), data = trn, ...)
    }
  })
  
  lmnames <- names(lmall$coefficients)
  lmall$coefficients <- rowMeans(as.data.frame(sapply(models, function(m) coef(m))))
   
  names(lmall$coefficients) <- lmnames
  lmall$fitted.values <- predict(lmall, data)
  target <- trimws(gsub('~.*$', '', formula))
  lmall$residuals <- data[, target] - lmall$fitted.values
  
  return(lmall)
}
