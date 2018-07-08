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
#'trn <- Wilshire[Wilshire$date<='2014-12-01',]
#'xeq <- xtr <- trn[c('CorpProfits', 'FedFundsRate', 'UnempRate')]
#'model1 <- ecm(trn$Wilshire5000, xeq, xtr)
#'
#'#Use 2014-12-01 and onwards data as test data to predict
#'tst <- Wilshire[Wilshire$date>='2014-12-01',]
#'
#'#predict on tst using model1 and initial FedFundsRate
#'tst$model1Pred <- ecmpredict(model1, tst, tst$Wilshire5000[1])
#'
#'@export
#'@importFrom stats predict
ecmpredict <- function (model, newdata, init) {
  if(nrow(newdata) < 3){
    stop("Your input for 'newdata' has less that 3 rows. This insufficient to make proper ECM predictions.")
  }
  if (sum(grepl('^delta', names(model$coefficients))) >= 1) {
    form <- names(model$coefficients)
    xtrnames <- form[grep("^delta", form)]
    xtrnames <- substr(xtrnames, 6, max(nchar(xtrnames)))
    xtr <- newdata[which(names(newdata) %in% xtrnames)]
    xtr <- data.frame(apply(xtr, 2, diff, 1))
    names(xtr) <- paste0("delta", names(xtr))
    xtrnames <- names(xtr)
  }
  
  if (sum(grepl('Lag1$', names(model$coefficients))) > 1) {
    form <- names(model$coefficients)
    xeqnames <- form[grep("^(?!delta).*", form, perl = T)]
    xeqnames <- xeqnames[-c(1, length(xeqnames))]
    xeqnames <- substr(xeqnames, 1, unlist(lapply(gregexpr("Lag", xeqnames), function(x) x[length(x)])) - 1)
    xeq <- newdata[which(names(newdata) %in% xeqnames)]
    names(xeq) <- paste0(names(xeq), "Lag1")
    xeqnames <- names(xeq)
  }
  
  if (exists('xeq')) {
    if (ncol(xeq) > 1) {
      xeq <- rbind(rep(NA, ncol(xeq)), xeq[1:(nrow(xeq) - 1), ])
    } else if (ncol(xeq) == 1) {
      xeq <- data.frame(c(NA, xeq[1:(nrow(xeq) - 1), ]))
    }
  }   
  
  if (exists('xeq') & exists('xtr')) {
    x <- cbind(xtr, xeq[complete.cases(xeq), ])
    x$yLag1 <- init
    names(x) <- c(xtrnames, xeqnames, "yLag1")
  } else if (!exists('xeq') & exists('xtr')) {
    x <- xtr 
    x$yLag1 <- init
    names(x) <- c(xtrnames, "yLag1")
  } else if (exists('xeq') & !exists('xtr')) {
    x <- as.data.frame(xeq[complete.cases(xeq),])
    x$yLag1 <- init
    names(x) <- c(xeqnames, "yLag1")
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

