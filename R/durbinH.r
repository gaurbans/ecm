#'Calculate Durbin's h-statistic
#'
#'Calculates Durbin's h-statistic for autoregressive models.
#'@param model The model being assessed
#'@param ylag1var The variable in the model that represents the lag of the y-term
#'@return Numeric Durbin's h statistic 
#'@details
#'Using the Durbin-Watson (DW) test for autoregressive models is inappropriate because the DW test itself 
#'tests for first order autocorrelation. This doesn't apply to an ECM model, for which the DW test is 
#'still valid, but the durbinH function in included here in case an autoregressive model has been built.
#'If Durbin's h-statistic is greater than 1.96, it is likely that autocorrelation exists.
#'
#'@seealso \code{lm}
#'@examples
#'#Build a simple AR1 model to predict performance of the Wilshire 5000 Index
#'data(Wilshire)
#'Wilshire$Wilshire5000Lag1 <- c(NA, Wilshire$Wilshire5000[1:(nrow(Wilshire)-1)])
#'Wilshire <- Wilshire[complete.cases(Wilshire),]
#'AR1model <- lm(Wilshire5000 ~ Wilshire5000Lag1, data=Wilshire)
#'
#'#Check Durbin's h-statistic on AR1model
#'durbinH(AR1model, "Wilshire5000Lag1")
#'#The h-statistic is 4.23, which means there is likely autocorrelation in the data.
#'
#'@export
#'@importFrom car durbinWatsonTest
durbinH <- function(model, ylag1var){
  d <- car::durbinWatsonTest(model)
  n <- length(model$fitted.values) + 1
  v <- summary(model)$coef[which(row.names(summary(model)$coef)==ylag1var),2]^2
  durbinH <- (1 - 0.5 * d$dw) * sqrt(n / (1 - n*v))
  return(durbinH)
}



