#'Predict using an ecm object
#'
#'Takes an ecm object and uses it to predict based on new data. This prediction does the undifferencing required to transform the change in y back to y itself.
#'@param model ecm object used to make predictions
#'@param newdata Data frame to on which to predict
#'@param init Initial value for prediction
#'@return Numeric predictions on new data based ecm object
#'@details 
#'Since error correction models only model the change in the target variable, an initial value must be specified. Additionally, the 'newdata' parameter should have at least 3 rows of data.
#'@examples
#'##Not run
#'
#'data(Wilshire)
#'
#'#Rebuilding model1 from ecm example
#'trn <- Wilshire[Wilshire$date<='2015-12-01',]
#'xeq <- xtr <- trn[c('CorpProfits', 'FedFundsRate', 'UnempRate')]
#'model1 <- ecm(trn$Wilshire5000, xeq, xtr)
#'model2 <- ecm(trn$Wilshire5000, xeq, xtr, linearFitter='earth')
#'
#'#Use 2015-12-01 and onwards data as test data to predict
#'tst <- Wilshire[Wilshire$date>='2015-12-01',]
#'
#'#predict on tst using model1 and initial FedFundsRate
#'tst$model1Pred <- ecmpredict(model1, tst, tst$Wilshire5000[1])
#'tst$model2Pred <- ecmpredict(model2, tst, tst$Wilshire5000[1])
#'
#'@export
#'@importFrom stats predict
#'@importFrom utils tail
ecmpredict <- function (model, newdata, init) {
  if (class(model)=='earth') {
    coef_names <- model$namesx
  } else if (class(model)=='lm') {
    coef_names <- names(model$coefficients)
  }
  
  lags <- as.integer(substring(tail(coef_names, 1), nchar(tail(coef_names, 1))))
  if(nrow(newdata) < (lags+2)){
    stop(paste0("Your input for 'newdata' has fewer rows than necessary to predict on a model with ", lags, " lags."))
  }
  if (sum(grepl('^delta', coef_names)) >= 1) {
    form <- coef_names
    xtrnames <- form[grep("^delta", form)]
    xtrnames <- substr(xtrnames, 6, max(nchar(xtrnames)))
    xtr <- newdata[which(names(newdata) %in% xtrnames)]
    xtr <- data.frame(apply(xtr, 2, diff, lags))
    names(xtr) <- paste0("delta", names(xtr))
    xtrnames <- names(xtr)
  }
  
  if (sum(grepl('Lag[0-9]$', coef_names)) > 1) {
    form <- coef_names
    xeqnames <- form[grep("^(?!delta).*", form, perl = T)]
    if ('(Intercept)' %in% xeqnames){
      xeqnames <- xeqnames[-c(1, length(xeqnames))]
    }
    xeqnames <- substr(xeqnames, 1, unlist(lapply(gregexpr("Lag", xeqnames), function(x) x[length(x)])) - 1)
    xeq <- newdata[which(names(newdata) %in% xeqnames)]
    names(xeq) <- paste0(names(xeq), "Lag", lags)
    xeqnames <- names(xeq)
  }
  
  if (exists('xeq')) {
    xeq <- data.frame(sapply(xeq, lagpad, lags))
  }   
  
  if (exists('xeq') & exists('xtr')) {
    x <- cbind(xtr, xeq[complete.cases(xeq), ])
    x$yLag1 <- init
    names(x) <- c(xtrnames, xeqnames, paste0("yLag", lags))
  } else if (!exists('xeq') & exists('xtr')) {
    x <- xtr 
    x$yLag1 <- init
    names(x) <- c(xtrnames, paste0("yLag", lags))
  } else if (exists('xeq') & !exists('xtr')) {
    x <- as.data.frame(xeq[complete.cases(xeq),])
    x$yLag1 <- init
    names(x) <- c(xeqnames, paste0("yLag", lags))
  }
  
  modelpred <- predict(model, x[1, ])
  for (i in 2:nrow(x)) {
    x$yLag1[i] <- x$yLag1[i - 1] + modelpred
    modelpred <- predict(model, x[i, ])
  }
  modelpred <- predict(model, x)
  modelpred <- cumsum(c(init, modelpred))
  return(modelpred)
}

