#'Build an averaged error correction model
#'
#'Builds multiple ECM models on subsets of the data and averages them. See the lmave function for more details
#'on the methodology and use cases for this approach. 
#'@param y The target variable
#'@param xeq The variables to be used in the equilibrium term of the error correction model
#'@param xtr The variables to be used in the transient term of the error correction model
#'@param lags The number of lags to use
#'@param includeIntercept Boolean whether the y-intercept should be included
#'@param k The number of models or data partitions desired
#'@param method Whether to split data by folds ("fold"), nested folds ("nestedfold"), or bootstrapping ("boot")
#'@param seed Seed for reproducibility (only needed if method is "boot")
#'@param weights Optional vector of weights to be passed to the fitting process
#'@param ... Additional arguments to be passed to the 'lm' function (careful in that these may need to be modified for ecm or may not be appropriate!)
#'@return an lm object representing an error correction model
#'@details
#'In some cases, instead of building an ECM on the entire dataset, it may be preferable to build k ECM models on k subsets of the data, each subset containing (k-1)/k*nrow(data)
#'observations of the full dataset, and then average their coefficients. Reasons to do this include controlling for overfitting or extending the training sample. For example, 
#'in many time series modeling exercises, the holdout test sample is often the latest few months or years worth of data. Ideally, it's desirable to include these data since 
#'they likely have more future predictive power than older observations. However, including the entire dataset in the training sample could result in overfitting, or using a 
#'different time period as the test sample may be even less representative of future performance. One potential solution is to build multiple ECM models using the entire dataset, 
#'each with a different holdout test sample, and then average them to get a final ECM. This approach is somewhat similar to the idea of random forest regression, in which 
#'multiple regression trees are built on subsets of the data and then averaged. 
#'
#'This function only works with the 'lm' linear fitter. 
#'@seealso \code{lm}
#'@examples
#'##Not run
#'
#'#Use ecm to predict Wilshire 5000 index based on corporate profits, 
#'#Federal Reserve funds rate, and unemployment rate
#'data(Wilshire)
#'
#'#Use 2015-12-01 and earlier data to build models
#'trn <- Wilshire[Wilshire$date<='2015-12-01',]
#'
#'#Build five ECM models and average them to get one model
#'xeq <- xtr <- trn[c('CorpProfits', 'FedFundsRate', 'UnempRate')]
#'model1 <- ecmave(trn$Wilshire5000, xeq, xtr, includeIntercept=TRUE, k=5)
#'
#'@export
#'@importFrom stats lm
ecmave <- function (y, xeq, xtr, lags=1, includeIntercept = TRUE, k, method = 'boot', seed = 5, weights = NULL, ...) {
  if (!is.null(xtr)){
    if (sum(grepl("^delta|Lag[0-9]$", names(xtr))) > 0) {
      warning("You have column name(s) in xtr that begin with 'delta' or end with 'Lag[0-9]'. It is strongly recommended that you change this, otherwise the function 'ecmpredict' may result in errors or incorrect predictions.")
    }
    if (!inherits(xtr, "data.frame")) {
      stop("xtr does not inherit class 'data.frame'. See details on how to input them as data frames.")
    }
  }
  
  if (!is.null(xeq)){
    if (sum(grepl("^delta", names(xeq))) > 0) {
      warning("You have column name(s) in xeq that begin with 'delta' or end with 'Lag[0-9]'. It is strongly recommended that you change this, otherwise the function 'ecmpredict' may result in errors or incorrect predictions.")
    }
    if (!inherits(xeq, "data.frame")) {
      stop("xeq does not inherit class 'data.frame'. See details on how to input them as data frames.")
    }
    if (nrow(xeq) < (lags+1)) {
      stop("Insufficient data for the lags specified.")
    }
  }
  
  if (!is.null(xeq)){
    xeqnames <- names(xeq)
    xeqnames <- paste0(xeqnames, paste0("Lag", as.character(lags)))
    xeq <- data.frame(sapply(xeq, lagpad, lags))
  }
  
  if (!is.null(xtr)){
    xtrnames <- names(xtr)
    xtrnames <- paste0("delta", xtrnames)
    xtr <- data.frame(apply(xtr, 2, diff, lags))
  }
  
  if (class(y)=='data.frame'){
    if (ncol(y) > 1){
      warning("You have more than one column in y, only the first will be used")
    }
    y <- y[,1]
  }
  yLag <- y[1:(length(y) - lags)]
  
  if (!is.null(xtr) & !is.null(xeq)){
    x <- cbind(xtr, xeq[complete.cases(xeq), ])
    xnames <- c(xtrnames, xeqnames)
  } else if (!is.null(xtr) & is.null(xeq)){
    x <- xtr
    xnames <- xtrnames
  } else if (is.null(xtr) & !is.null(xeq)){
    x <- xeq[complete.cases(xeq), ]
    xnames <- xeqnames
  }
  
  x <- cbind(x, yLag)
  names(x) <- c(xnames, paste0("yLag", as.character(lags)))
  x$dy <- diff(y, lags)
  
  if (includeIntercept){
    formula <- 'dy ~ .'
  } else {
    formula <- 'dy ~ . - 1'
  }
  ecm <- lmave(formula, data = x, k = k, method = method, seed = seed, weights = weights, ...)
  
  return(ecm)
}
