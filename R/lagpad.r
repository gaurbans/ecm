#'Lag a vector
#'
#'Create a vector of the lag of a variable and fill missing values with NA's.
#'@param x A vector to be lagged
#'@param k The number of lags to output
#'@return The lagged vector with NA's in missing values
#'@export
lagpad <- function(x, k=1) {
  if (k < 1) {
    stop("k needs to be 1 or greater")
  }
  
  i<-is.vector(x)
  if(is.vector(x)) x <- matrix(x) else x <- matrix(x, nrow(x))
  x <- rbind(matrix(rep(NA, k*ncol(x)), ncol=ncol(x)), matrix(x[1:(nrow(x) - k),], ncol = ncol(x)))
  if(i) x[1:length(x)] else x
}
