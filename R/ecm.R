#'Build an error correction model
#'
#'Builds an lm object that represents an error correction model (ECM) by automatically differencing and
#'lagging predictor variables according to ECM methodology.
#'@param y The target variable
#'@param xeq The variables to be used in the equilibrium term of the error correction model
#'@param xtr The variables to be used in the transient term of the error correction model
#'@param lags The number of lags to use
#'@param includeIntercept Boolean whether the y-intercept should be included
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
#'#From a strictly statistical standpoint, Wilshire data may not be stationary
#'#and hence model1 and model2 may have heteroskedasticity in the residuals.
#'#Let's check for that.
#'lmtest::bptest(model1)
#'lmtest::bptest(model2)
#'#The Breush-Pagan tests suggest we should reject homoskedasticity in the residuals for both models.
#'
#'lmtest::bgtest(model1)
#'lmtest::bgtest(model2)
#'#The Breusch-Godfrey tests suggest we should reject that there is no serial correlation 
#'#in the residuals.
#'
#'#Given the above issues, see adjusted std. errors and p-values for our models.
#'lmtest::coeftest(model1, vcov=sandwich::NeweyWest)
#'lmtest::coeftest(model2, vcov=sandwich::NeweyWest)
#'
#'@export
#'@importFrom stats lm
ecm <- function (y, xeq, xtr, lags=1, includeIntercept = TRUE, weights = NULL, linearFitter = 'lm', ...) {
  if (sum(grepl("^delta|Lag[0-9]$", names(xtr))) > 0 | sum(grepl("^delta", names(xeq))) > 0) {
    warning("You have column name(s) in xeq or xtr that begin with 'delta' or end with 'Lag[0-9]'. It is strongly recommended that you change this, otherwise the function 'ecmpredict' may result in errors or incorrect predictions.")
  }
  
  if (!inherits(xtr, "data.frame") | !inherits(xeq, "data.frame")) {
    stop("xeq or xtr does not inherit class 'data.frame'. See details on how to input them as data frames.")
  }
  
  if (nrow(xeq) < (lags+1)) {
    stop("Insufficient data for the lags specified.")
  }
  
  xeqnames <- names(xeq)
  xeqnames <- paste0(xeqnames, paste0("Lag", as.character(lags)))
  xeq <- data.frame(sapply(xeq, lagpad, lags))
  
  xtrnames <- names(xtr)
  xtrnames <- paste0("delta", xtrnames)
  xtr <- data.frame(apply(xtr, 2, diff, lags))
  
  if (class(y)=='data.frame'){
    if (ncol(y) > 1){
      warning("You have more than one column in y, only the first will be used")
    }
    y <- y[,1]
  }
  yLag <- y[1:(length(y) - lags)]
  
  x <- cbind(xtr, xeq[complete.cases(xeq), ])
  x <- cbind(x, yLag)
  names(x) <- c(xtrnames, xeqnames, paste0("yLag", as.character(lags)))
  x$dy <- diff(y, lags)
  
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
  
  
  return(ecm)
}
