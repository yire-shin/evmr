#' @name rggd.edtest
#' @aliases rggd.edtest
#' @title Entropy difference test for rGGD
#' @param data Data should be contain n rows, each a GGDr observation.
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
#' rggdEdtest(bangkok)
#' }
rggdEdtest<-function(data){

  data <- as.matrix(data)
  R    <- ncol(data)

  result <-matrix(0, R-1,9)

  for(i in 2:R){

    result[i-1,1] <- i
    fit <-rggdEd(data[,1:i])
    result[i-1,2]  <-fit$p.value
    result[i-1,3]  <-fit$statistics
    result[i-1,4:6]<-fit$theta
    result[i-1,7]  <-fit$ybar

  }

  result[,8]  <-rev(eva::pSeqStop(rev(result[, 2]))$ForwardStop)
  result[,9]  <-rev(eva::pSeqStop(rev(result[, 2]))$StrongStop)


  colnames(result) <-c("r","p.values","statistic","est.loc","est.scale","est.shape","ybar","ForwardStop", "StrongStop")
  as.data.frame(result)
}
