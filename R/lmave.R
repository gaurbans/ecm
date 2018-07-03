#'Build multiple lm models and average them
#'
#'Builds k lm models on k partitions of the data and averages their coefficients to get create one model. 
#'Each partition excludes k/nrow(data) observations. 
#'@param formula The formula to be passed to lm
#'@param data The data to be used
#'@param k The number of models or data partitions desired
#'@param method Whether to split data by folds ("fold") or by bootstrapping ("boot")
#'@param seed Seed for reproducibility (only needed if method is "boot")
#'@param weights Optional vector of weights to be passed to the fitting process
#'@param ... Additional arguments to be passed to the 'lm' function
#'@return an lm object
#'@details
#'In some cases--especially in some time series modeling (see ecmave function)--rather than building one model on the entire dataset, it may be preferable to build multiple models on subsets 
#'of the data and average them. The lmave function splits the data into k partitions, builds k models, each on the data without one of the partitions, and then averages the coefficients of these 
#'models to get a final model. This is similar to averaging multiple tree regression models in algorithms like random forest. 
#'@references Jung, Y. & Hu, J. (2016). "A K-fold Averaging Cross-validation Procedure". https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5019184/
#'@seealso \code{lm}
#'@examples
#'#Build linear models to predict Wilshire 5000 index based on corporate profits, 
#'#Federal Reserve funds rate, and unemployment rate
#'data(Wilshire)
#'
#'#Build one model on the entire dataset
#'modelall <- lm(Wilshire5000 ~ ., data = Wilshire)
#'
#'#Build a five fold averaged linear model on the entire dataset
#'modelave <- lmave('Wilshire5000 ~ .', data = Wilshire, k = 5) 
#'
#'@export
#'@importFrom stats lm coef
lmave <- function(formula, data, k, method = 'boot', seed = 5, weights = NULL, ...){
  lmall <- lm(formula, data, ...)
  modellist <- 1:k
  
  if (method == 'fold'){
    folds <- cut(seq(1, nrow(data)), breaks=k, labels=FALSE)
    models <- lapply(modellist, function(i) {
      tstIdx <- which(folds==i, arr.ind = TRUE)
      trn <- data[-tstIdx, ]
      if (!is.null(weights)){
        w <- weights[-tstIdx]
        lm(as.formula(formula), data = trn, weights = w, ...)
      } else if (is.null(weights)) {
        lm(as.formula(formula), data = trn, ...)
      }
      
    })
  } else if (method == 'boot'){
    set.seed(seed)
    models <- lapply(modellist, function(i) {
      tstIdx <- sample(nrow(data), 1 / k * nrow(data))
      trn <- data[-tstIdx, ]
      if (!is.null(weights)){
        w <- weights[-tstIdx]
        lm(as.formula(formula), data = trn, weights = w, ...)
      } else if (is.null(weights)){
        lm(as.formula(formula), data = trn, ...)
      }
    })
  }
  
  lmnames <- names(lmall$coefficients)
  lmall$coefficients <- rowMeans(as.data.frame(sapply(models, function(m) coef(m))))
  names(lmall$coefficients) <- lmnames
  lmall$fitted.values <- predict(lmall, data)
  target <- trimws(gsub('~.*$', '', formula))
  lmall$residuals <- data[, target] - lmall$fitted.values
  
  return(lmall)
}

