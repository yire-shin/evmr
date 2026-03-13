#' Sequential Entropy Difference Test for rK4D Models
#'
#' Performs the sequential entropy difference (ED) test for selecting the
#' number of order statistics in the r-largest four-parameter kappa
#' distribution (rK4D) model.
#'
#' The procedure computes ED tests sequentially for \eqn{r = 2, \dots, R} and
#' applies the ForwardStop and StrongStop stopping rules to control the
#' false discovery rate.
#'
#' @param data A numeric matrix or data frame containing the r-largest order
#'   statistics. Each row represents one observation (or block), and columns
#'   must be ordered from largest to smallest.
#'
#' @return A data frame containing:
#' \itemize{
#'   \item \code{r}: value of \eqn{r} tested
#'   \item \code{p.values}: raw p-values from the entropy difference tests
#'   \item \code{statistic}: test statistics for each value of \eqn{r}
#'   \item \code{est.loc}: estimated location parameter
#'   \item \code{est.scale}: estimated scale parameter
#'   \item \code{est.shape1}: estimated first shape parameter
#'   \item \code{est.shape2}: estimated second shape parameter
#'   \item \code{ybar}: mean entropy difference
#'   \item \code{ForwardStop}: adjusted values from the ForwardStop rule
#'   \item \code{StrongStop}: adjusted values from the StrongStop rule
#' }
#'
#' @details
#' The function sequentially applies the entropy difference test
#' (\code{\link{rk4dEd}}) for increasing values of \eqn{r}. The resulting
#' p-values are adjusted using the ForwardStop and StrongStop procedures
#' to help determine an appropriate value of \eqn{r}.
#'
#' @references
#' Bader, B., Yan, J., & Zhang, X. (2017).
#' Automated selection of \eqn{r} for the r-largest order statistics approach.
#' \emph{Statistics and Computing}.
#' \doi{10.1007/s11222-016-9697-3}
#'
#' Shin, Y., & Park, J.-S. (2023).
#' Modeling climate extremes using the four-parameter kappa distribution
#' for r-largest order statistics.
#' \emph{Weather and Climate Extremes}.
#' \doi{10.1016/j.wace.2022.100533}
#'
#' @seealso \code{\link{rk4dEd}}, \code{\link{rk4d.fit}}
#' @export
#'
#' @examples
#' \dontrun{
#' data(bangkok)
#' rk4dEdtest(bangkok)
#' }
rk4dEdtest <- function(data) {

  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a matrix or data frame.")
  }

  data <- as.matrix(data)

  if (!is.numeric(data)) {
    stop("'data' must be numeric.")
  }

  if (ncol(data) < 2) {
    stop("'data' must contain at least two columns.")
  }

  if (any(!is.finite(data), na.rm = TRUE)) {
    stop("'data' must contain only finite values, except possible NA padding.")
  }

  R <- ncol(data)

  result <- matrix(NA_real_, nrow = R - 1, ncol = 10)

  for (i in 2:R) {
    result[i - 1, 1] <- i

    fit <- rk4dEd(data[, 1:i, drop = FALSE])

    result[i - 1, 2] <- fit$p.value
    result[i - 1, 3] <- fit$statistics
    result[i - 1, 4:7] <- fit$theta
    result[i - 1, 8] <- fit$ybar
  }

  result[, 9] <- rev(eva::pSeqStop(rev(result[, 2]))$ForwardStop)
  result[, 10] <- rev(eva::pSeqStop(rev(result[, 2]))$StrongStop)

  colnames(result) <- c(
    "r", "p.values", "statistic", "est.loc", "est.scale",
    "est.shape1", "est.shape2", "ybar", "ForwardStop", "StrongStop"
  )

  as.data.frame(result)
}
