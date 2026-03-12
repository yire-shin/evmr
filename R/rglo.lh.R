#' @name rglo.lh
#' @aliases rglo.lh
#' @title loglikelihood for rglo
#' @param data Vector or matrix of observations. If x is a matrix, each row is taken to be a new observation.
#' @param par Location, scale, and shape parameters. Can be vectors, but the lengths must be appropriate.
#'
#' @return log-likelihood value
#' @export
#'
rgloLh <- function(data,par) {

  R <-ncol(data)
  nr<-nrow(data)

  mu <-par[1]
  sc <-par[2]
  xi <-par[3]

  ri  <- (R-seq(1:(R))) # r-i
  cr  <- (1+ri)    # c_r

  log.den3 = -R*log(sc) + sum(log(cr)) + (1+R) * log(1 / (1+((1-xi*(data[,R]-mu)/sc)^(1/xi))))+
    rowSums(((1/xi) -1)*log(1-xi*((data-mu)/sc)))

  return(log.den3)

}
