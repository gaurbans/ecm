#'Build an averaged error correction model
#'
#'Builds multiple ECM models on subsets of the data and averages them.
#'@param y The target variable
#'@param xeq The variables to be used in the equilibrium term of the error correction model
#'@param xtr The variables to be used in the transient term of the error correction model
#'@param includeIntercept Boolean whether the y-intercept should be included
#'@param k The number of models or data partitions desired
#'@param method Whether to split data by folds ("fold") or by bootstrapping ("boot")
#'@param seed Seed for reproducibility (only needed if method is "boot")
#'@param weights Optional vector of weights to be passed to the fitting process
#'@param ... Additional arguments to be passed to the 'lm' function (careful in that these may need to be modified for ecm or may not be appropriate!)
#'@return an lm object representing an error correction model
#'@details
#'In some cases, instead of building an ECM on the entire dataset, it may be preferable to build k ECM models on k subsets of the data and then average their coefficients. 
#'Reasons to do this include controlling for overfitting or extending the training sample. For example, in many time series modeling exercises, the holdout test sample is 
#'often the latest few months or years worth of data. Ideally it's desirable to include these data since they likely have more future predictive power than older observations.
#'However, including the entire dataset in the training sample could result in overfitting, or using a different time period as the test sample may be even less representative 
#'of future performance. One potential solution is to build multiple ECM models, each with a different holdout test sample, and then average them to get a final ECM. This 
#'approach is somewhat similar to the idea of random forest regression, in which multiple regression trees are built on subsets of the data and then averaged. 
#'
#'@seealso \code{lm}
#'@examples
#'#Use ecm to predict Wilshire 5000 index based on corporate profits, 
#'#Federal Reserve funds rate, and unemployment rate
#'data(Wilshire)
#'
#'#Use 2014-12-01 and earlier data to build models
#'trn <- Wilshire[Wilshire$date<='2014-12-01',]
#'
#'#Build five ECM models and average them to get one model
#'xeq <- xtr <- trn[c('CorpProfits', 'FedFundsRate', 'UnempRate')]
#'model1 <- ecmave(trn$Wilshire5000, xeq, xtr, includeIntercept=TRUE, k=5)
#'
#'@export
#'@importFrom stats lm
ecmave <- function (y, xeq, xtr, includeIntercept = TRUE, k, method = 'boot', seed = 5, weights = NULL, ...) {
  if (sum(grepl("^delta|Lag1$", names(xtr))) > 0 | sum(grepl("^delta", names(xeq))) > 0) {
    warning("You have column name(s) in xeq or xtr that begin with 'delta' or end with 'Lag1'. It is strongly recommended that you change this, otherwise the function 'ecmpredict' may result in errors or incorrect predictions.")
  }
  
  if (class(xtr) != "data.frame" | class(xeq) != "data.frame") {
    stop("xeq or xtr is not of class 'data.frame'. See details on how to input them as data frames.")
  }
  
  xeqnames <- names(xeq)
  xeqnames <- paste0(xeqnames, "Lag1")
  ifelse(ncol(xeq) > 1, xeq <- rbind(rep(NA, ncol(xeq)), xeq[1:(nrow(xeq) - 1), ]), xeq <- data.frame(c(NA, xeq[1:(nrow(xeq) - 1), ])))
  
  xtrnames <- names(xtr)
  xtrnames <- paste0("delta", xtrnames)
  xtr <- data.frame(apply(xtr, 2, diff, 1))
  
  yLag1 <- y[1:(length(y) - 1)]
  x <- cbind(xtr, xeq[complete.cases(xeq), ])
  x <- cbind(x, yLag1)
  names(x) <- c(xtrnames, xeqnames, "yLag1")
  x$dy <- diff(y, 1)
  
  if (includeIntercept){
    formula <- 'dy ~ .'
  } else {
    formula <- 'dy ~ . - 1'
  }
  ecm <- lmave(formula, data = x, k = k, method = method, seed = seed, weights = weights, ...)
  
  return(ecm)
}











