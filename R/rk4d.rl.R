#' @name rk4d.rl
#' @aliases rk4d.rl
#' @title return level for k4d
#' @param z An object returned by \code{rk4dfit}. The object should represent a stationary model.
#' @param year The return level (i.e.\ the profile likelihood is for the value that is exceeded with probability 1/\code{m}).
#' @param show Logical; if \code{TRUE} (the default), print details of the result.
#' @return The return level and numeric vector giving the standard errors for the return level.
#' @export
#'
#' @examples
#' \dontrun{
#' data(bangkok)
#' result<-rk4d.fit(bangkok)
#' rk4d.rl(result$mle,result$cov,100)
#' }
rk4d.rl<-function (z, year=c(20,50,100,200),show=F)
{

  del <-matrix(ncol=length(year), nrow=4)

  mu  = z$mle[1]
  sig = z$mle[2]
  xi  = as.vector(z$mle[3])
  h   = as.vector(z$mle[4])
  p   = 1-(1/year)
  yp  = (1-((p)^h)) / h

  # d(z)/d(mu)
  del[1,] <- 1

  # d(z)/d(sig)
  del[2,] <- (1/xi) - (1/xi)*(yp^xi)

  # d(z)/d(xi)
  del[3,] <- -(sig/xi)*log(yp)*(yp^xi)+(sig/(xi^2))*(yp^xi)-(sig/(xi^2))

  # d(z)/d(h)
  del[4,] <-  sig*(yp^(xi-1))*((log(p)*(p^h)*h + (1 - (p^h))) / h^2)

  del.t <-t(del)

  z$rl    <-mu + (sig/xi) - (sig/xi) * (yp^xi)
  z$rlse <-diag(sqrt((del.t %*% z$cov %*% del)))

  names(z$rl)    <-paste0(as.character(year),"y")
  names(z$rlse) <-paste0(as.character(year),"y")

  if(show) print(z[c(12,13)])
  invisible(z)
}
