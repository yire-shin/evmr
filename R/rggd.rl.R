#' @name rggd.rl
#' @aliases rggd.rl
#' @title return level for ggd
#' @param z An object returned by \code{rggd.fit}. The object should represent a stationary model.
#' @param year The return level (i.e.\ the profile likelihood is for the value that is exceeded with probability 1/\code{m}).
#' @param show Logical; if \code{TRUE} (the default), print details of the result.
#' @return The return level and numeric vector giving the standard errors for the return level.
#' @export
#'
#' @examples
#' \dontrun{
#' data(bangkok)
#' result<-rk4d.fit(bangkok)
#' rggd.rl(result$mle,result$cov,100)
#' }
rggd.rl<-function(z, year=c(20,50,100,200), show=F){

  del <-matrix(ncol=length(year),nrow=3)

  mu = z$mle[1]
  sig= z$mle[2]
  h  = z$mle[3]

  f = 1-(1/year) # F
  f_inv = 1/f

  del[1,] <- 1

  del[2,] <- -log( (1-exp(h*log(f)))/h )

  del[3,] <- sig * (1/h + ((f^h)*log(f))/(1-(f^h)))

  del.t <-t(del)

  options(digits=8)
  if(h>0){
    z$rl = mu - sig * log( (1-exp(h*log(f)))/h )
  }else{
    z$rl = mu + sig * log( (-h)/(expm1(h*log(f))))
  }

  z$rlse <-diag(sqrt((del.t %*% z$cov %*% del)))

  names(z$rl)    <-paste0(as.character(year),"y")
  names(z$rlse) <-paste0(as.character(year),"y")

  if(show) print(z[c(12,13)])
  invisible(z)

}
