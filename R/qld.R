#' @name qld
#' @aliases qld
#' @title qunatile of logistic
#' @param p Vector of quantiles.
#' @param loc,scale Location, scale parameters. Can be vectors, but the lengths must be appropriate.
#'
#' @return quantile function for the 1-LD distribution
#' @export
#'
#' @examples
#' \dontrun{
#' qld(0.5, loc = 0, scale=1)
#' }
qld<-function (p,loc = 0, scale = 1)
{
  mu = loc
  sig= scale

  z = mu + sig * log( 1/(expm1(-log(p))))
  z
}
