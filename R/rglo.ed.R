#' @name rglo.ed
#' @aliases rglo.ed
#' @title Entropy difference for rGLO
#' @param data Data should be contain n rows, each a rGLO observation.
#' @param par Location, scale, and shape parameters. Can be vectors, but the lengths must be appropriate.
#'
#' @return list containing the following components. A subset of these components are printed after the fit.  \code{statistic}, \code{p.value}, \code{theta} are always printed.#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(bangkok)
#' rgloEd(bangkok)
#' }
rgloEd<-function(data,par=NULL){

  if(is.null(par)){

    y      <-rglo.fit(data,show=F)

    theta1 <-y$mle
    #theta1[3] <- -theta1[3]

  }else{

    theta1<-par
    #theta1[3]<- -theta1[3]

  }

  R <-ncol(data)
  nr<-nrow(data)

  Diff <- rgloLh(data[,1:R],theta1)- rgloLh(as.matrix(data[,1:(R-1)],ncol=R-1),theta1)

  EstVar   <-sum((Diff -mean(Diff))^2)/(nr-1)

  #FirstMom <- -log(theta1[2])+log(R)+(-theta1[3]+1)*digamma(R)+(1+R)*log(1/(1+R))-R*log(1/R)
  eta <- -log(theta1[2])+log(R)-((R+1)/(R))+theta1[3]*(digamma(1)-digamma(R))

  Diff1 <-sum(Diff)/nr
  Diff <-sqrt(nr)*(Diff1-eta)/sqrt(EstVar)
  p.value <-2*(1-stats::pnorm(abs(Diff)))

  out<-list(statistics=as.numeric(Diff),p.value=as.numeric(p.value),theta=theta1,ybar=as.numeric(Diff1))

  out

}
