#' @name rk4d.edtest
#' @aliases rk4d.edtest
#' @title Entropy difference test for rK4D
#' @param data Data should be contain n rows, each a K4Dr observation.
#'
#' @return {
#'Function returns a dataframe containing the test statistics, estimates, and p-value results of the sequential tests.
#' \item{r}{Value of r to be tested.}
#' \item{p.values}{Raw p-values from the individual tests at each value of r.}
#' \item{ForwardStop}{Transformed p-values according to the ForwardStop stopping rule.}
#' \item{StrongStop}{Transformed p-values according to the StrongStop stopping rule.}
#' \item{statistic}{Returned test statistics of each individual test.}
#' \item{est.loc}{Estimated location parameter for the given r.}
#' \item{est.scale}{Estimated scale parameter for the given r.}
#' \item{est.shape}{Estimated shape parameter for the given r.}
#'}
#' @export
#'
#' @examples
#' \dontrun{
#' data(bangkok)
#' rk4dEdtest(bangkok)
#' }
rk4dEdtest<-function(data){

  data <- as.matrix(data)
  R    <- ncol(data)

  result <-matrix(0, R-1,10)

  for(i in 2:R){

    result[i-1,1] <- i
    fit <-rk4dEd(data[,1:i])
    result[i-1,2]  <-fit$p.value
    result[i-1,3]  <-fit$statistics
    result[i-1,4:7]<-fit$theta
    result[i-1,8]  <-fit$ybar

  }

  result[,9]  <-rev(eva::pSeqStop(rev(result[, 2]))$ForwardStop)
  result[,10]  <-rev(eva::pSeqStop(rev(result[, 2]))$StrongStop)


  colnames(result) <-c("r","p.values","statistic","est.loc","est.scale","est.shape1","est.shape2","ybar","ForwardStop", "StrongStop")
  as.data.frame(result)
}
