#'Build an error correction model
#'
#'Builds an lm object that represents an error correction model (ECM) by automatically differencing and
#'lagging predictor variables according to ECM methodology.
#'@param y The target variable
#'@param xeq The variables to be used in the equilibrium term of the error correction model
#'@param xtr The variables to be used in the transient term of the error correction model
#'@param includeIntercept Boolean whether the y-intercept should be included (should be set to TRUE if using 'earth' as linearFitter)
#'@param weights Optional vector of weights to be passed to the fitting process
#'@param linearFitter Whether to use 'lm' or 'earth' to fit the model
#'@param ... Additional arguments to be passed to the 'lm' or 'earth' function (careful that some arguments may not be appropriate for ecm!)
#'@return an lm object representing an error correction model
#'@details
#'The general format of an ECM is \deqn{\Delta y_{t} = \beta_{0} + \beta_{1}\Delta x_{1,t} +...+ \beta_{i}\Delta x_{i,t} + \gamma(y_{t-1} - (\alpha_{1}x_{1,t-1} +...+ \alpha_{i}x_{i,t-1})).}
#'The ecm function here modifies the equation to the following: \deqn{\Delta y = \beta_{0} + \beta_{1}\Delta x_{1,t} +...+ \beta_{i}\Delta x_{i,t} + \gamma y_{t-1} + \gamma_{1}x_{1,t-1} +...+ \gamma_{i}x_{i,t-1},}
#'\deqn{where \gamma_{i} = -\gamma \alpha_{i},} 
#'so it can be modeled as a simpler ordinary least squares (OLS) function using R's lm function.
#'
#'Ordinarily, the ECM uses lag=1 when differencing the transient term and lagging the equilibrium term, as specified in the equation above. However, the ecm function here gives the user the ability to specify a lag greater than 1. 
#'
#'Notice that an ECM models the change in the target variable (y). This means that the predictors will be lagged and differenced,
#'and the model will be built on one observation less than what the user inputs for y, xeq, and xtr. If these arguments contain vectors with too few observations (eg. one single observation),
#'the function will not work. Additionally, for the same reason, if using weights in the ecm function, the length of weights should be one less than the number of rows in xeq or xtr.
#'
#'When inputting a single variable for xeq or xtr in base R, it is important to input it in the format "xeq=df['col1']" so they inherit the class 'data.frame'. Inputting such as "xeq=df[,'col1']" or "xeq=df$col1" will result in errors in the ecm function. You can load data via other R packages that store data in other formats, as long as those formats also inherit the 'data.frame' class.
#'
#'By default, base R's 'lm' is used to fit the model. However, users can opt to use 'earth', which uses Jerome Friedman's Multivariate Adaptive Regression Splines (MARS) to build a regression model, which transforms each continuous variable into piece-wise linear hinge functions. This allows for non-linear features in both the transient and equilibrium terms.
#'
#'ECM models are used for time series data. This means the user may need to consider stationarity and/or cointegration before using the model.
#'@seealso \code{lm, earth}
#'@examples
#'##Not run
#'
#'#Use ecm to predict Wilshire 5000 index based on corporate profits, 
#'#Federal Reserve funds rate, and unemployment rate.
#'data(Wilshire)
#'
#'#Use 2015-12-01 and earlier data to build models
#'trn <- Wilshire[Wilshire$date<='2015-12-01',]
#'
#'#Assume all predictors are needed in the equilibrium and transient terms of ecm.
#'xeq <- xtr <- trn[c('CorpProfits', 'FedFundsRate', 'UnempRate')]
#'model1 <- ecm(trn$Wilshire5000, xeq, xtr, includeIntercept=TRUE)
#'
#'#Assume CorpProfits and FedFundsRate are in the equilibrium term, 
#'#UnempRate has only transient impacts.
#'xeq <- trn[c('CorpProfits', 'FedFundsRate')]
#'xtr <- trn['UnempRate']
#'model2 <- ecm(trn$Wilshire5000, xeq, xtr, includeIntercept=TRUE)
#'
#'@export
#'@importFrom stats lm
ecm <- function (y, xeq, xtr, includeIntercept = TRUE, weights = NULL, linearFitter = 'lm', ...) {
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
    if (nrow(xeq) < 2) {
      stop("Insufficient data for the lags specified.")
    }
  }
  
  if (!is.null(xeq)){
    xeqnames <- names(xeq)
    xeqnames <- paste0(xeqnames, "Lag1")
    xeq <- data.frame(sapply(xeq, lagpad))
  }
  
  if (!is.null(xtr)){
    xtrnames <- names(xtr)
    xtrnames <- paste0("delta", xtrnames)
    xtr <- data.frame(apply(xtr, 2, diff))
  }
  
  if (is.data.frame(y)){
    if (ncol(y) > 1){
      warning("You have more than one column in y, only the first will be used")
    }
    y <- y[,1]
  }
  yLag <- y[1:(length(y) - 1)]
  
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
  names(x) <- c(xnames, paste0("yLag1"))
  x$dy <- diff(y)
  
  if (linearFitter=='lm'){
    if (includeIntercept){
      ecm <- lm(dy ~ ., data = x, weights = weights, ...)
    } else {
      ecm <- lm(dy ~ . - 1, data = x, weights = weights, ...)
    }
  } else if (linearFitter=='earth'){
    if (includeIntercept){
      ecm <- earth::earth(dy ~ ., data = x, weights = weights, linpreds='yLag1', ...)
    } else {
      ecm <- earth::earth(dy ~ . - 1, data = x, weights = weights, linpreds='yLag1', ...)
    }
  }
  ecm$linearFitter <- linearFitter
  
  return(ecm)
}
