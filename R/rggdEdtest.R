#' Sequential Entropy Difference Test for rGGD Models
#'
#' Performs the sequential entropy difference (ED) test for selecting the
#' number of order statistics in the r-largest generalized Gumbel distribution
#' (rGGD) model.
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
#' \item \code{r} Value of \eqn{r} tested
#' \item \code{p.values} Raw p-values from the entropy difference tests
#' \item \code{statistic} Test statistics for each value of \eqn{r}
#' \item \code{est.loc} Estimated location parameter
#' \item \code{est.scale} Estimated scale parameter
#' \item \code{est.shape} Estimated shape parameter
#' \item \code{ybar} Mean entropy difference
#' \item \code{ForwardStop} Adjusted values from the ForwardStop rule
#' \item \code{StrongStop} Adjusted values from the StrongStop rule
#' }
#'
#' @details
#' The function sequentially applies the entropy difference test
#' (\code{\link{rggdEd}}) for increasing values of \eqn{r}.
#' The columns of \code{data} must represent decreasing order statistics
#' within each row, with the first column containing the block maximum.
#' The resulting p-values are adjusted using the ForwardStop and StrongStop
#' procedures to help determine an appropriate value of \eqn{r}.
#'
#' @references
#' Shin, Y., & Park, J.-S. (2025).
#' Generalized Gumbel model for r-largest order statistics with application
#' to peak streamflow.
#' \emph{Scientific Reports}.
#' \doi{10.1038/s41598-024-83273-y}
#'
#' Bader, B., Yan, J., & Zhang, X. (2017).
#' Automated selection of \eqn{r} for the r-largest order statistics approach.
#' \emph{Statistics and Computing}.
#' \doi{10.1007/s11222-016-9697-3}
#'
#' @seealso \code{\link{rggdEd}}, \code{\link{rggd.fit}}
#' @export
#'
#' @examples
#' \dontrun{
#' data(bangkok)
#' rggdEdtest(bangkok)
#' }
rggdEdtest <- function(data) {

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

  result <- matrix(NA_real_, R - 1, 9)

  for (i in 2:R) {
    result[i - 1, 1] <- i

    fit <- rggdEd(data[, 1:i, drop = FALSE])

    result[i - 1, 2] <- fit$p.value
    result[i - 1, 3] <- fit$statistics
    result[i - 1, 4:6] <- fit$theta
    result[i - 1, 7] <- fit$ybar
  }

  result[, 8] <- rev(eva::pSeqStop(rev(result[, 2]))$ForwardStop)
  result[, 9] <- rev(eva::pSeqStop(rev(result[, 2]))$StrongStop)

  colnames(result) <- c(
    "r", "p.values", "statistic", "est.loc", "est.scale",
    "est.shape", "ybar", "ForwardStop", "StrongStop"
  )

  as.data.frame(result)
}
