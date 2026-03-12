#' @name qggd
#' @aliases qggd
#' @title qunatile of ggd
#' @param p Vector of quantiles.
#' @param loc,scale,shape Location, scale, and shape parameters. Can be vectors, but the lengths must be appropriate.
#'
#' @return quantile function for the 1-GGD distribution
#' @export
#'
#' @examples
#' \dontrun{
#' qggd(0.5, loc = 0, scale=1, shape = 0.1)
#' }
qggd<-function (p,loc = 0, scale = 1, shape = 0)
{
  mu = loc
  sig= scale
  h  = shape

  if(h>0){
    f = p
    f_inv = 1/f
    z = mu - sig * log( (1-exp(h*log(f)))/h )
  }else{
    f = p
    f_inv = 1/f
    z = mu + sig * log( (-h)/(expm1(h*log(f))))
  }
  z
}
