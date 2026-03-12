#' @name rk4d.lh
#' @aliases rk4d.lh
#' @title loglikelihood for rk4d
#' @param data Vector or matrix of observations. If x is a matrix, each row is taken to be a new observation.
#' @param par Location, scale, and shape parameters. Can be vectors, but the lengths must be appropriate.
#'
#' @return log-likelihood value
#' @export
#'
rk4d.lh <- function(data,par) {


  R <-ncol(data)
  nr<-nrow(data)

  mu <-par[1]
  sc <-par[2]
  xi <-par[3]
  h  <-par[4]

  ri  <- (R - seq(1:(R)))
  cr  <- (1 - ri * h)

  y <- 1 - xi * (data - mu) / sc
  f <- 1 - h * (1 - xi * (data[,R] - mu) / sc)^(1/xi)

  y <- log(sc) + (1 - 1/xi) * log(y) - log(cr)
  y <- rowSums(y, na.rm=TRUE)

  log.den3 <- (R * h - 1) / h * log(f) + y

  return(-log.den3)

}
