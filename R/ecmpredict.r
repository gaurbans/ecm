#'Predict using an ecm object
#'
#'Takes an ecm object and uses it to predict based on new data. This prediction does the undifferencing required to transform the change in y back to y itself.
#'@param model ecm object used to make predictions
#'@param newdata Data frame to on which to predict
#'@param init Initial value(s) for prediction 
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
  if(nrow(newdata) < 2){
    stop("Your input for 'newdata' has insufficient data.")
  }
  
  if (model$linearFitter=='earth') {
    form <- model$namesx
  } else if (model$linearFitter=='lm') {
    form <- names(model$coefficients)
  }
  
  if (sum(grepl('^delta', form)) >= 1) {
    xtrfctnames <- form[grep("^delta", form)]
    xtrfctnames <- substr(xtrfctnames, 6, max(nchar(xtrfctnames)))
    xtrfct <- newdata[which(names(newdata) %in% xtrfctnames)]
    xtrfct <- data.frame(apply(xtrfct, 2, diff))
    names(xtrfct) <- paste0("delta", names(xtrfct))
    xtrfctnames <- names(xtrfct)
  }
  
  if (sum(grepl('Lag1$', form)) > 1) {
    xeqfctnames <- form[grep("^(?!delta).*", form, perl = T)]
    xeqfctnames <- substr(xeqfctnames, 1, unlist(lapply(gregexpr("Lag", xeqfctnames), function(x) x[length(x)])) - 1)
    xeqfct <- newdata[which(names(newdata) %in% xeqfctnames)]
    names(xeqfct) <- paste0(names(xeqfct), "Lag1")
    xeqfct <- data.frame(sapply(xeqfct, lagpad))
    xeqfctnames <- names(xeqfct)
  }
  
  if (exists('xeqfct') & exists('xtrfct')) {
    x <- cbind(xtrfct, xeqfct[complete.cases(xeqfct), ])
    x$yLag1 <- init
    names(x) <- c(xtrfctnames, xeqfctnames, "yLag1")
  } else if (!exists('xeqfct') & exists('xtrfct')) {
    x <- xtrfct 
    x$yLag1 <- init
    names(x) <- c(xtrfctnames, "yLag1")
  } else if (exists('xeqfct') & !exists('xtrfct')) {
    x <- as.data.frame(xeqfct[complete.cases(xeqfct),])
    x$yLag1 <- init
    names(x) <- c(xeqfctnames, "yLag1")
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

