#' @name rgd.rl
#' @aliases rgd.rl
#' @title return level for gumbel
#' @param z An object returned by \code{rgd.fit}. The object should represent a stationary model.
#' @param year The return level (i.e.\ the profile likelihood is for the value that is exceeded with probability 1/\code{m}).
#' @param show Logical; if \code{TRUE} (the default), print details of the result.
#' @return The return level and numeric vector giving the standard errors for the return level.
#' @export
#'
#' @examples
#' \dontrun{
#' data(bangkok)
#' result<-rgd.fit(bangkok)
#' rgd.rl(result$mle,result$cov,100)
#' }
rgd.rl<-function(z, year=c(20, 50, 100, 200),show=F){

  del <-matrix(ncol=length(year),nrow=2)

  mu = z$mle[1]
  sig= z$mle[2]

  p   = 1/year
  yp  = -log(1-p)

  del[1,] <- 1

  del[2,] <- -log(yp)

  del.t <-t(del)

  z$rl   <-mu - sig * log(yp)
  z$rlse <-diag(sqrt((del.t %*% z$cov %*% del)))

  names(z$rl)    <-paste0(as.character(year),"y")
  names(z$rlse)  <-paste0(as.character(year),"y")

  if(show) print(z[c(12,13)])
  invisible(z)


}
