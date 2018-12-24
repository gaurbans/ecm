#'Get cumulative count
#'
#'Get the cumulative count of a variable of interest
#'@param x A vector for which to get cumulative count
#'@return The cumulative count of all items in x
#'@export
cumcount <- function(x){
  cumcount <- numeric(length(x))
  names(cumcount) <- x
  
  for(i in 1:length(x)){
    cumcount[i] <- sum(x[1:i]==x[i])
  }
  
  return(cumcount)
}
