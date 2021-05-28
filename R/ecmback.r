#'Backwards selection to build an error correction model
#'
#'Much like the ecm function, this builds an error correction model.
#'However, it uses backwards selection to select the optimal predictors based on lowest AIC or BIC, or highest adjusted R-squared, rather than using all predictors.
#'@param y The target variable
#'@param xeq The variables to be used in the equilibrium term of the error correction model
#'@param xtr The variables to be used in the transient term of the error correction model
#'@param lags The number of lags to use
#'@param includeIntercept Boolean whether the y-intercept should be included
#'@param criterion Whether AIC (default), BIC, or adjustedR2 should be used to select variables
#'@param weights Optional vector of weights to be passed to the fitting process
#'@param keep Optional character vector of variables to forcibly retain
#'@param ... Additional arguments to be passed to the 'lm' function (careful in that these may need to be modified for ecm or may not be appropriate!)
#'@return an lm object representing an error correction model using backwards selection
#'@details
#'When inputting a single variable for xeq or xtr, it is important to input it in the format "xeq=df['col1']" in order to retain the data frame class. Inputting such as "xeq=df[,'col1']" or "xeq=df$col1" will result in errors in the ecm function.
#'
#'If using weights, the length of weights should be one less than the number of rows in xeq or xtr. 
#'
#'This function only works with the 'lm' linear fitter. The 'earth' linear fitter already does some variable selection, so one can use that via that 'ecm' function. 
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
#'modelback <- ecmback(trn$Wilshire5000, xeq, xtr)
#'print(modelback)
#'#Backwards selection chose CorpProfits and FedFundsRate in the equilibrium term, 
#'#CorpProfits and UnempRate in the transient term.
#'
#'modelbackFFR <- ecmback(trn$Wilshire5000, xeq, xtr, keep = 'UnempRate')
#'print(modelbackFFR)
#'#Backwards selection was forced to retain UnempRate in both terms.
#'
#'@export
#'@importFrom stats lm complete.cases
ecmback <- function (y, xeq, xtr, lags=1, includeIntercept = T, criterion = "AIC", weights = NULL, keep = NULL, ...) {
  if (sum(grepl("^delta|Lag1$", names(xtr))) > 0 | sum(grepl("^delta|Lag1$", names(xeq))) > 0) {
    warning("You have column name(s) in xeq or xtr that begin with 'delta' or end with 'Lag1'. It is strongly recommended that you change this, otherwise the function 'ecmpredict' will result in errors or incorrect predictions.")
  }
  
  if (class(xtr) != "data.frame" | class(xeq) != "data.frame") {
    stop("xeq or xtr is not of class 'data.frame'. See details on how to input them as data frames.")
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
  
  if (includeIntercept) {
    formula <- 'dy ~ .'
  } else {
    formula <- 'dy ~ . - 1'
  }
  full <- lm(as.formula(formula), data = x, weights = weights, ...)
  dontdropIdx <- numeric(2)
  
  if (criterion == "AIC" | criterion == "BIC") {
    if (criterion == "AIC") {
      kIC = 2
    } else if (criterion == "BIC") {
      kIC = log(nrow(x))
    }
    
    fullAIC <- partialAIC <- AIC(full, k = kIC)
    while (partialAIC <= fullAIC & length(rownames(drop1(full))) > length(dontdropIdx)) {
      dontdropVars <- "^<none>$|^yLag1$"
      if (!is.null(keep)) {
        for (i in 1:length(keep)){
          dontdropVars <- paste0(dontdropVars, "|^delta", keep[i], "$", "|^", keep[i], "Lag1$")
        }
      }
      dontdropIdx <- grep(dontdropVars, rownames(drop1(full, k = kIC)))
      todrop <- rownames(drop1(full, k = kIC))[-dontdropIdx][which.min(drop1(full, k = kIC)$AIC[-dontdropIdx])]
      x <- x[-which(names(x) %in% todrop)]
      possible <- lm(as.formula(formula), data = x, weights = weights, ...)
      partialAIC <- AIC(possible)
      if (partialAIC < fullAIC & length(rownames(drop1(full))) > length(dontdropIdx)) {
        fullAIC <- partialAIC
        full <- possible
        ecmeq <- full
      } else {
        ecmeq <- full
      }
    }
  } else if (criterion == "adjustedR2") {
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
      partial <- lm(dy ~ ., data = newx, weights = weights, ...)
      partialAdjR2 <- summary(partial)$adj.r.sq
      if (partialAdjR2 >= fullAdjR2 & length(full$coefficients) > length(dontdropIdx)) {
        x <- newx
        full <- partial
        ecmeq <- full
      } else {
        ecmeq <- full
      }
    }
  }
  
  if (sum(grepl("^delta", names(ecmeq$coefficients))) == 0) {
    warning("Backwards selection has opted to leave out all transient terms from the final model. This means you have a first order differenced autoregressive model of sorts, not a full error correction model.")
  } else if (sum(grepl("Lag1$", names(ecmeq$coefficients))) == 0) {
    warning("Backwards selection has opted to leave out all equilibrium terms from the final model. This means you have a first order differenced autoregressive model of sorts, not a full error correction model.")
  }
  
  return(ecmeq)
}
