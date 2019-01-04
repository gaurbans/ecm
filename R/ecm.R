#'Build an error correction model
#'
#'Builds an lm object that represents an error correction model (ECM) by automatically differencing and
#'lagging predictor variables according to ECM methodology.
#'@param y The target variable
#'@param xeq The variables to be used in the equilibrium term of the error correction model
#'@param xtr The variables to be used in the transient term of the error correction model
#'@param includeIntercept Boolean whether the y-intercept should be included
#'@param weights Optional vector of weights to be passed to the fitting process
#'@param ... Additional arguments to be passed to the 'lm' function (careful in that these may need to be modified for ecm or may not be appropriate!)
#'@return an lm object representing an error correction model
#'@details
#'The general format of an ECM is \deqn{\Delta y = \beta_{0} + \beta_{1}\Delta x_{1,t} +...+ \beta_{i}\Delta x_{i,t} + \gamma(y_{t-1} - (\alpha_{1}x_{1,t-1} +...+ \alpha_{i}x_{i,t-1})).}
#'The ecm function here modifies the equation to the following: \deqn{\Delta y = \beta_{0} + \beta_{1}\Delta x_{1,t} +...+ \beta_{i}\Delta x_{i,t} + \gamma y_{t-1} + \gamma_{1}x_{1,t-1} +...+ \gamma_{i}x_{i,t-1},}
#'\deqn{where \gamma_{i} = -\gamma \alpha_{i},} 
#'so it can be modeled as a simpler ordinary least squares (OLS) function using R's lm function.
#'
#'Notice that an ECM models the change in the target variable (y). This means that the predictors will be lagged and differenced,
#'and the model will be built on one observation less than what the user inputs for y, xeq, and xtr. If these arguments contain vectors with too few observations (eg. one single observation),
#'the function will not work. Additionally, for the same reason, if using weights in the ecm function, the length of weights should be one less than the number of rows in xeq or xtr.
#'
#'When inputting a single variable for xeq or xtr, it is important to input it in the format "xeq=df['col1']" in order to retain the data frame class. Inputting such as "xeq=df[,'col1']" or "xeq=df$col1" will result in errors in the ecm function.
#'
#'ECM models are used for time series data. This means the user may need to consider stationarity and/or cointegration before using the model.
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
#'#Assume all predictors are needed in the equilibrium and transient terms of ecm
#'xeq <- xtr <- trn[c('CorpProfits', 'FedFundsRate', 'UnempRate')]
#'model1 <- ecm(trn$Wilshire5000, xeq, xtr, includeIntercept=TRUE)
#'
#'#Assume CorpProfits and FedFundsRate are in the equilibrium term, 
#'#UnempRate has only transient impacts
#'xeq <- trn[c('CorpProfits', 'FedFundsRate')]
#'xtr <- trn['UnempRate']
#'model2 <- ecm(trn$Wilshire5000, xeq, xtr, includeIntercept=TRUE)
#'
#'@export
#'@importFrom stats lm
ecm <- function (y, xeq, xtr, includeIntercept = TRUE, weights = NULL, ...) {
  if (sum(grepl("^delta|Lag1$", names(xtr))) > 0 | sum(grepl("^delta", names(xeq))) > 0) {
    warning("You have column name(s) in xeq or xtr that begin with 'delta' or end with 'Lag1'. It is strongly recommended that you change this, otherwise the function 'ecmpredict' may result in errors or incorrect predictions.")
  }
  
  if (class(xtr) != "data.frame" | class(xeq) != "data.frame") {
    stop("xeq or xtr is not of class 'data.frame'. See details on how to input them as data frames.")
  }
  
  xeqnames <- names(xeq)
  xeqnames <- paste0(xeqnames, "Lag1")
  if(ncol(xeq) > 1) {
    xeq <- rbind(rep(NA, ncol(xeq)), xeq[1:(nrow(xeq) - 1), ])
  } else {
    xeq <- data.frame(c(NA, xeq[1:(nrow(xeq) - 1), ]))
  }

  xtrnames <- names(xtr)
  xtrnames <- paste0("delta", xtrnames)
  xtr <- data.frame(apply(xtr, 2, diff, 1))
  
  if (class(y)=='data.frame'){
    if (ncol(y) > 1){
      warning("You have more than one column in y, only the first will be used")
    }
    y <- y[,1]
  }
  yLag1 <- y[1:(length(y) - 1)]
  
  x <- cbind(xtr, xeq[complete.cases(xeq), ])
  x <- cbind(x, yLag1)
  names(x) <- c(xtrnames, xeqnames, "yLag1")
  x$dy <- diff(y, 1)
  
  if (includeIntercept){
    ecm <- lm(dy ~ ., data = x, weights = weights, ...)
  } else {
    ecm <- lm(dy ~ . - 1, data = x, weights = weights, ...)
  }
  
  return(ecm)
}
