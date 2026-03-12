#' @name rggd.ed
#' @aliases rggd.ed
#' @title Entropy difference for rGGD
#' @param data Data should be contain n rows, each a rGGD observation
#'
#' @return list containing the following components. A subset of these components are printed after the fit.  \code{statistic}, \code{p.value}, \code{theta} are always printed.#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(bangkok)
#' rggdEd(bangkok)
#' }
rggdEd<-function(data){

  y      <-rggd.fit(data,show=F)

  theta1 <-y$mle

  mu<-theta1[1]
  sc<-theta1[2]
  h <-theta1[3]

  R <-ncol(data)
  nr<-nrow(data)

  Diff <- rggdLh(data[,1:R],theta1)- rggdLh(as.matrix(data[,1:(R-1)],ncol=R-1),theta1)
  EstVar   <-sum((Diff -mean(Diff))^2)/(nr-1)

  ri  <- (R-seq(1:(R))) # r-i
  cr  <- prod((1-ri*h))     # c_r

  ri1  <- (R-1-seq(1:(R-1)))
  cr1  <- prod((1-ri1*h))

  if(h>0){
    ar <- (1-(R-1)*h)/h
    ar1 <- (1-(R-2)*h)/h
    term1 <-log(sc)
    term2 <-log(1-(R-1)*h)
    term3 <-((1-R*h)/h)*(digamma(ar)-digamma(ar+R)) # ((1/h)^(R))*(cr/factorial(R-1))*beta(ar,R) = 1
    term4 <-((1-(R-1)*h)/h)*(digamma(ar1)-digamma(ar1+R-1)) # ((1/h)^(R-1))*cr1/factorial(R-2)*beta(ar1,R-1) = 1
    term5 <-((digamma(R)-digamma(ar+R))-log(h)) # cr/factorial(R-1)*((1/h)^(R))*beta(ar,R) =1

    eta <- -term1+term2+term3-term4+term5

  }else{
    ar <- (1-(R-1)*h)/h
    ar1 <- (1-(R-2)*h)/h
    term1 <- log(sc)
    term2 <- log(1-(R-1)*h)
    term3 <-((1-R*h)/h)*(digamma((1/-h)+R)-digamma(1/-h)) # cr/factorial(R-1)*(-1)*((1/-h)^(R))*beta(R,1/-h)=1
    term4 <-((1-(R-1)*h)/h)*(digamma((1/-h)+R-1)-digamma(1/-h)) # *cr1/factorial(R-2)*(-1)*((1/h)^(R-1))*beta(R-1,1/-h)=1
    term5 <-(digamma(R)-digamma(1/-h)-log(-h)) # cr/factorial(R-1)*((1/-h)^(R))*beta(R,1/-h) =1

    eta <- -term1+term2+term3-term4+term5

  }

  Diff1 <-sum(Diff)/nr
  Diff <-sqrt(nr)*(Diff1-eta)/sqrt(EstVar)
  p.value <-2*(1-stats::pnorm(abs(Diff)))

  out<-list(statistics=as.numeric(Diff),p.value=as.numeric(p.value),theta=theta1,ybar=as.numeric(Diff1))

  out

}

