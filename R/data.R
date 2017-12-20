#'FRED data on the Wilshire 5000 index and other economic factors
#'
#'A dataset containing quarterly performance of the Wilshire 5000 index, corporate profits, Federal Reserve funds rate, and the unemployment rate.
#'@usage data(Wilshire)
#'@format A data frame with 182 rows and 6 variables:
#'\describe{
#' \item{date}{monthly date}
#' \item{Wilshire5000}{quarterly Wilshire 5000 index, in value}
#' \item{CorpProfits}{quarterly corporate profits, in value}
#' \item{FedFundsRate}{quarterly federal funds rate, in percent}
#' \item{UnempRate}{quarterly unemployment rate, in percent}
#'}
#'@source \url{https://fred.stlouisfed.org/}
"Wilshire"