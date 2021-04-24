#'Backwards selection using an averaged error correction model
#'
#'Much like the ecmback function, ecmaveback uses backwards selection to build an error correction model.
#'However, it uses the averaging method of ecmave to build models and then choose variables based on lowest AIC or BIC, or highest adjusted R-squared.
#'@param y The target variable
#'@param xeq The variables to be used in the equilibrium term of the error correction model
#'@param xtr The variables to be used in the transient term of the error correction model
#'@param lags The number of lags to use
#'@param includeIntercept Boolean whether the y-intercept should be included
#'@param k The number of models or data partitions desired
#'@param criterion Whether AIC (default), BIC, or adjustedR2 should be used to select variables
#'@param method Whether to split data by folds ("fold"), nested folds ("nestedfold"), or bootstrapping ("boot")
#'@param seed Seed for reproducibility (only needed if method is "boot")
#'@param weights Optional vector of weights to be passed to the fitting process
#'@param keep Optional character vector of variables to forcibly retain
#'@param ... Additional arguments to be passed to the 'lm' function (careful in that these may need to be modified for ecm or may not be appropriate!)
#'@return an lm object representing an error correction model using backwards selection
#'@details
#'When inputting a single variable for xeq or xtr, it is important to input it in the format "xeq=df['col1']" in order to retain the data frame class. Inputting such as "xeq=df[,'col1']" or "xeq=df$col1" will result in errors in the ecm function.
#'
#'If using weights, the length of weights should be one less than the number of rows in xeq or xtr. 
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
#'#Use backwards selection to choose which predictors are needed 
#'xeq <- xtr <- trn[c('CorpProfits', 'FedFundsRate', 'UnempRate')]
#'modelaveback <- ecmaveback(trn$Wilshire5000, xeq, xtr, k = 5)
#'print(modelaveback)
#'#Backwards selection chose CorpProfits and FedFundsRate in the equilibrium term, 
#'#CorpProfits and UnempRate in the transient term.
#'
#'modelavebackFFR <- ecmaveback(trn$Wilshire5000, xeq, xtr, k = 5, keep = 'UnempRate')
#'print(modelavebackFFR)
#'#Backwards selection was forced to retain UnempRate in both terms.
#'
#'@export
#'@importFrom stats lm complete.cases AIC as.formula drop1
ecmaveback <- function (y, xeq, xtr, lags=1, includeIntercept = T, criterion = "AIC", k, method = 'boot', seed = 5, weights = NULL, keep = NULL, ...) {
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
    formula <- "dy ~ ."
  } else {
    formula <- "dy ~ . - 1"
  }
  full <- lmave(formula, data = x, k = k, method = method, seed = seed, weights = weights, ...)
  dontdropIdx <- numeric(2)
  
  if (criterion == "AIC" | criterion == "BIC") {
    if (criterion == "AIC") {
      kIC <- 2
    } else if (criterion == "BIC") {
      kIC <- log(nrow(x))
    }
    
    fullAIC <- partialAIC <- AIC(full, k=kIC)
    while (partialAIC <= fullAIC & length(rownames(drop1(full))) > length(dontdropIdx)){
      dontdropVars <- "^<none>$|^yLag1$"
      if (!is.null(keep)) {
        for (i in 1:length(keep)){
          dontdropVars <- paste0(dontdropVars, "|^delta", keep[i], "$", "|^", keep[i], "Lag1$")
        }
      } 
      dontdropIdx <- grep(dontdropVars, rownames(drop1(full, k = kIC)))
      todrop <- rownames(drop1(full, k = kIC))[-dontdropIdx][which.min(drop1(full, k = kIC)$AIC[-dontdropIdx])]
      x <- x[-which(names(x) %in% todrop)]
      possible <- lmave(formula, data = x, k = k, method = method, seed = seed, weights = weights, ...)
      partialAIC <- AIC(possible)
      if (partialAIC < fullAIC & length(rownames(drop1(full))) > length(dontdropIdx)){
        fullAIC <- partialAIC
        full <- possible
      } else {
        ecm <- full
      }
    }
  } else if (criterion == 'adjustedR2'){
    fullAdjR2 <- partialAdjR2 <- summary(full)$adj.r.sq
    while (partialAdjR2 >= fullAdjR2 & length(full$coefficients) > length(dontdropIdx)) {
      fullAdjR2 <- summary(full)$adj.r.sq
      if (!is.null(keep)) {
        dontdropVars <- paste0("^delta", keep, "$", "|^", keep, "Lag1$")
        dontdropVars <- paste0(dontdropVars, collapse = '|')
        dontdropIdx <- grep(dontdropVars, rownames(summary(full)$coef))
        if (includeIntercept) {
          dontdropIdx <- c(1, dontdropIdx)
        }
        todrop <- which.max(summary(full)$coef[-dontdropIdx, 4])
      } else {
        if (includeIntercept) {
          todrop <- which.max(summary(full)$coef[-1, 4])
        } else {
          todrop <- which.max(summary(full)$coef[, 4])
        }
      }
      newx <- x[-todrop]
      partial <- lmave(formula, data = x, k = k, method = method, seed = seed, weights = weights, ...)
      partialAdjR2 <- summary(partial)$adj.r.sq
      if (partialAdjR2 >= fullAdjR2 & length(full$coefficients) > length(dontdropIdx)) {
        x <- newx
        full <- partial
      } else {
        ecm <- full
      }
    }
  }
  
  if (sum(grepl("^delta", names(ecm$coefficients))) == 0) {
    warning("Backwards selection has opted to leave out all transient terms from the final model. This means you have a first order differenced autoregressive model of sorts, not a full error correction model.")
  } else if (sum(grepl("Lag1$", names(ecm$coefficients))) == 0) {
    warning("Backwards selection has opted to leave out all equilibrium terms from the final model. This means you have a first order differenced autoregressive model of sorts, not a full error correction model.")
  }
  return(ecm)
}
