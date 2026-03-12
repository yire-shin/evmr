#' @name rk4d.ed
#' @aliases rk4d.ed
#' @title Entropy difference for rK4d
#' @param data Data should be contain n rows, each a rK4D observation
#'
#' @return list containing the following components. A subset of these components are printed after the fit.  \code{statistic}, \code{p.value}, \code{theta} are always printed.#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(bangkok)
#' rk4dEd(bangkok)
#' }
rk4dEd<-function(data){

  y      <-rk4d.fit(data,show=F)

  theta1 <-y$mle

  mu<-theta1[1]
  sc<-theta1[2]
  xi<-theta1[3]
  h <-theta1[4]

  R <-ncol(data)
  nr<-nrow(data)

  Diff <- rk4d.lh(data[,1:R],theta1)- rk4d.lh(as.matrix(data[,1:(R-1)],ncol=R-1),theta1)

  EstVar   <-sum((Diff -mean(Diff))^2)/(nr-1)


  if(h>0){
    ar <- (1-(R-1)*h)/h
    ar1 <- (1-(R-2)*h)/h
    term1 <-log(sc)
    term2 <-log(1-(R-1)*h)
    term3 <-((1-R*h)/h)*(digamma(ar)-digamma(ar+R))
    term4 <-((1-(R-1)*h)/h)*(digamma(ar1)-digamma(ar1+R-1))
    term5 <-(1-xi)*((digamma(R)-digamma(ar+R))-log(h))

    eta <- -term1+term2+term3-term4+term5

  }else{
    ar <- (1-(R-1)*h)/h
    ar1 <- (1-(R-2)*h)/h
    term1 <- log(sc)
    term2 <- log(1-(R-1)*h)
    term3 <-((1-R*h)/h)*(digamma((1/-h)+R)-digamma(1/-h))
    term4 <-((1-(R-1)*h)/h)*(digamma((1/-h)+R-1)-digamma(1/-h))
    term5 <-(1-xi)*(digamma(R)-digamma(1/-h)-log(-h))

    eta <- -term1+term2+term3-term4+term5

  }


  Diff1 <-sum(Diff)/nr
  Stat <-sqrt(nr)*(Diff1-eta)/sqrt(EstVar)
  p.value <-2*(1-stats::pnorm(abs(Stat)))

  out<-list(statistics=as.numeric(Stat),p.value=as.numeric(p.value),theta=theta1,ybar=as.numeric(Diff1))

  out

}
