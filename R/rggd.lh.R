#' @name rggd.lh
#' @aliases rggd.lh
#' @title loglikelihood for rggd
#' @param data Vector or matrix of observations. If x is a matrix, each row is taken to be a new observation.
#' @param par Location, scale, and shape parameters. Can be vectors, but the lengths must be appropriate.
#'
#' @return log-likelihood value
#' @export
#'
rggdLh <- function(data,par) {


  R <-ncol(data)
  nr<-nrow(data)

  mu <-par[1]
  sc <-par[2]
  h  <-par[3]

  ri  <- (R-seq(1:(R))) # r-i
  cr  <- (1-ri*h)     # c_r

  y <- exp(-(data - mu)/sc)
  f <- 1 - h * exp(-(data[,R] - mu)/sc)


  y <- log(sc) - log(y) - log(cr)
  y <- rowSums(y, na.rm = TRUE)

  log.den3 = (R*h - 1)/h * log(f) + y

  return(-log.den3)

}
