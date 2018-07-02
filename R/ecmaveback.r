#'Backwards selection to build an error correction model
#'
#'Much like the ecm function, this builds an error correction model.
#'However, it uses backwards selection to select the optimal predictors based on lowest AIC or BIC, or highest adjusted R-squared, rather than using all predictors.
#'ecmback has the same parameters and output as ecm.
#'@param y The target variable
#'@param xeq The variables to be used in the equilibrium term of the error correction model
#'@param xtr The variables to be used in the transient term of the error correction model
#'@param includeIntercept Boolean whether the y-intercept should be included
#'@param criterion Whether AIC (default), BIC, or adjustedR2 should be used to select variables
#'@param weights Optional vector of weights to be passed to the fitting process
#'@param ... Additional arguments to be passed to the 'lm' function
#'@return an lm object representing an error correction model using backwards selection
#'@details
#'When inputting a single variable for xeq or xtr, it is important to input it in the format "xeq=df['col1']" in order to retain the data frame class. Inputting such as "xeq=df[,'col1']" or "xeq=df$col1" will result in errors in the ecm function.
#'@seealso \code{lm}
#'@examples
#'#Use ecm to predict Wilshire 5000 index based on corporate profits, 
#'#Federal Reserve funds rate, and unemployment rate
#'data(Wilshire)
#'
#'#Use 2014-12-01 and earlier data to build models
#'trn <- Wilshire[Wilshire$date<='2014-12-01',]
#'
#'#Use backwards selection to choose which predictors are needed 
#'xeq <- xtr <- trn[c('CorpProfits', 'FedFundsRate', 'UnempRate')]
#'modelback <- ecmback(trn$Wilshire5000, xeq, xtr, includeIntercept=TRUE, criterion = 'AIC')
#'print(modelback)
#'#Backwards selection chose CorpProfits in the equilibrium term, 
#'#CorpProfits and UnempRate in the transient term.
#'
#'@export
#'@importFrom stats lm complete.cases step
ecmaveback <- function (y, xeq, xtr, includeIntercept = T, criterion = "AIC", k = NULL, weights = NULL, ...) {
  if (sum(grepl("^delta|Lag1$", names(xtr))) > 0 | sum(grepl("^delta|Lag1$", names(xeq))) > 0) {
    warning("You have column name(s) in xeq or xtr that begin with 'delta' or end with 'Lag1'. It is strongly recommended that you change this, otherwise the function 'ecmpredict' will result in errors or incorrect predictions.")
  }
  
  xeqnames <- names(xeq)
  xeqnames <- paste0(xeqnames, "Lag1")
  xeq <- as.data.frame(xeq)
  ifelse(ncol(xeq) > 1, xeq <- rbind(rep(NA, ncol(xeq)), xeq[1:(nrow(xeq) - 1), ]), xeq <- data.frame(c(NA, xeq[1:(nrow(xeq) - 1), ])))
  
  xtrnames <- names(xtr)
  xtrnames <- paste0("delta", xtrnames)
  xtr <- as.data.frame(xtr)
  xtr <- data.frame(apply(xtr, 2, diff, 1))
  
  yLag1 <- y[1:(length(y) - 1)]
  x <- cbind(xtr, xeq[complete.cases(xeq), ])
  x <- cbind(x, yLag1)
  names(x) <- c(xtrnames, xeqnames, "yLag1")
  x$dy <- diff(y, 1)
  
  if (includeIntercept){
    formula <- "dy ~ ."
    full <- lmave(formula = formula, data = x, k = k, ...)
  } else {
    formula <- "dy ~ . - 1"
    full <- lmave(formula = formula, data = x, k = k, ...)
  }
  
  if (criterion == "AIC" | criterion == "BIC") {
    if (criterion == "AIC") {
      kIC <- 2
    } else if (criterion == "BIC") {
      kIC <- log(nrow(x))
    }
    
    fullAIC <- partialAIC <- AIC(full, k=kIC)
    while (partialAIC <= fullAIC){
      todrop <- rownames(drop1(full, k=kIC))[-grep('none|yLag1', rownames(drop1(full, k=kIC)))][which.min(drop1(full, k=kIC)$AIC[-grep('none|yLag1', rownames(drop1(full, k=kIC)))])]
      formula <- paste0(formula, ' - ', todrop)
      full <- lmave(formula = formula, data = x, k = k, ...)
      partialAIC <- AIC(full)
      if (partialAIC < fullAIC & length(rownames(drop1(full))) > 2){
        fullAIC <- partialAIC
      } else {
        ecm <- full
      }
    }
  } else if (criterion == 'adjustedR2'){
    fullAdjR2 <- partialAdjR2 <- summary(full)$adj.r.sq
    while (partialAdjR2 >= fullAdjR2) {
      fullAdjR2 <- summary(full)$adj.r.sq
      if (includeIntercept){
        todrop <- which.max(summary(full)$coef[-1, 4])
      } else {
        todrop <- which.max(summary(full)$coef[, 4])
      }
      newx <- x[which(!names(x) %in% names(todrop))]
      full <- lmave(formula = formula, data = newx, k = k, ...)
      partialAdjR2 <- summary(full)$adj.r.sq
      if (partialAdjR2 >= fullAdjR2) {
        x <- newx
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