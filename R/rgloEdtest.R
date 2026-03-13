#' Sequential Entropy Difference Test for rGLO Models
#'
#' Performs the sequential entropy difference (ED) test for selecting the
#' number of order statistics in the r-largest generalized logistic
#' distribution (rGLO) model.
#'
#' The procedure computes ED tests sequentially for \eqn{r = 2, \dots, R} and
#' applies the ForwardStop and StrongStop stopping rules to control the
#' false discovery rate.
#'
#' @param data A numeric matrix or data frame containing the r-largest order
#'   statistics. Each row represents one observation (or block), and columns
#'   must be ordered from largest to smallest.
#' @param par An optional numeric vector of length 3 giving the location,
#'   scale, and shape parameters. If \code{NULL}, parameters are estimated
#'   separately at each value of \eqn{r} using \code{\link{rgloEd}}.
#'
#' @return A data frame containing:
#' \itemize{
#'   \item \code{r}: value of \eqn{r} tested
#'   \item \code{p.values}: raw p-values from the entropy difference tests
#'   \item \code{statistic}: test statistics for each value of \eqn{r}
#'   \item \code{est.loc}: estimated location parameter
#'   \item \code{est.scale}: estimated scale parameter
#'   \item \code{est.shape}: estimated shape parameter
#'   \item \code{ybar}: mean entropy difference
#'   \item \code{ForwardStop}: adjusted values from the ForwardStop rule
#'   \item \code{StrongStop}: adjusted values from the StrongStop rule
#' }
#'
#' @details
#' The function sequentially applies the entropy difference test
#' (\code{\link{rgloEd}}) for increasing values of \eqn{r}. The resulting
#' p-values are adjusted using the ForwardStop and StrongStop procedures
#' to help determine an appropriate value of \eqn{r}.
#'
#' @references
#' Bader, B., Yan, J., & Zhang, X. (2017).
#' Automated selection of \eqn{r} for the r-largest order statistics approach.
#' \emph{Statistics and Computing}.
#' \doi{10.1007/s11222-016-9697-3}
#'
#' Shin, Y., & Park, J-S. (2024).
#' Generalized logistic model for r-largest order statistics with
#' hydrological application.
#' \emph{Stochastic Environmental Research and Risk Assessment}.
#' \doi{10.1007/s00477-023-02642-7}
#'
#' @seealso \code{\link{rgloEd}}, \code{\link{rglo.fit}}
#' @export
#'
#' @examples
#' \dontrun{
#' data(bangkok)
#' rgloEdtest(bangkok)
#' }
rgloEdtest <- function(data, par = NULL) {

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

  if (!is.null(par)) {
    if (!is.numeric(par) || length(par) != 3 || any(!is.finite(par))) {
      stop("'par' must be a numeric vector of length 3.")
    }
  }

  R <- ncol(data)

  result <- matrix(NA_real_, nrow = R - 1, ncol = 9)

  for (i in 2:R) {
    result[i - 1, 1] <- i

    fit <- rgloEd(data[, 1:i, drop = FALSE], par = par)

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
