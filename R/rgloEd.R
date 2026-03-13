#' Entropy Difference Test for rGLO Models
#'
#' Performs the entropy difference (ED) test for selecting the number of
#' order statistics in the r-largest generalized logistic distribution (rGLO)
#' model.
#'
#' The test compares the entropy of models fitted with \eqn{r} and
#' \eqn{r-1} order statistics and evaluates whether the additional order
#' statistic provides significant information.
#'
#' @param data A numeric matrix or data frame containing the r-largest
#'   order statistics. Each row represents one block or observation, and
#'   columns must be ordered from largest to smallest.
#' @param par An optional numeric vector of length 3 giving the location,
#'   scale, and shape parameters. If \code{NULL}, the parameters are estimated
#'   using \code{\link{rglo.fit}}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{statistics}: the entropy difference test statistic
#'   \item \code{p.value}: the two-sided p-value
#'   \item \code{theta}: the estimated or supplied parameter vector
#'   \item \code{ybar}: the sample mean entropy difference
#' }
#'
#' @details
#' This function applies the entropy difference test to the r-largest
#' generalized logistic model. If \code{par} is not supplied, the model
#' parameters are first estimated using \code{\link{rglo.fit}}.
#'
#' @references
#' Bader, B., Yan, J., & Zhang, X. (2017).
#' Automated selection of \eqn{r} for the r-largest order statistics approach.
#' \emph{Statistics and Computing}.
#' \doi{10.1007/s11222-016-9697-3}
#'
#' Shin, Y., & Park, J-S. (2024).
#' Generalized logistic model for r-largest order statistics with hydrological
#' application.
#' \emph{Stochastic Environmental Research and Risk Assessment}.
#' \doi{10.1007/s00477-023-02642-7}
#'
#' @seealso \code{\link{rglo.fit}}, \code{\link{rgloLh}}
#' @export
#'
#' @examples
#' \dontrun{
#' data(bangkok)
#' rgloEd(bangkok)
#' }
rgloEd <- function(data, par = NULL) {

  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a matrix or data frame.")
  }

  data <- as.matrix(data)

  if (!is.numeric(data)) {
    stop("'data' must be numeric.")
  }

  if (ncol(data) < 2) {
    stop("'data' must have at least two columns for the entropy difference test.")
  }

  if (any(!is.finite(data), na.rm = TRUE)) {
    stop("'data' must contain only finite values, except possible NA padding.")
  }

  if (is.null(par)) {
    y <- rglo.fit(data, show = FALSE)
    theta1 <- y$mle
  } else {
    if (!is.numeric(par) || length(par) != 3 || any(!is.finite(par))) {
      stop("'par' must be a numeric vector of length 3.")
    }
    theta1 <- par
  }

  R <- ncol(data)
  nr <- nrow(data)

  Diff <- rgloLh(data[, 1:R, drop = FALSE], theta1) -
    rgloLh(as.matrix(data[, 1:(R - 1), drop = FALSE]), theta1)

  EstVar <- sum((Diff - mean(Diff))^2) / (nr - 1)

  if (!is.finite(EstVar) || EstVar <= 0) {
    stop("Estimated variance is not positive; the entropy difference statistic cannot be computed.")
  }

  eta <- -log(theta1[2]) + log(R) - ((R + 1) / R) +
    theta1[3] * (digamma(1) - digamma(R))

  Diff1 <- sum(Diff) / nr
  Diff <- sqrt(nr) * (Diff1 - eta) / sqrt(EstVar)
  p.value <- 2 * (1 - stats::pnorm(abs(Diff)))

  out <- list(
    statistics = as.numeric(Diff),
    p.value = as.numeric(p.value),
    theta = theta1,
    ybar = as.numeric(Diff1)
  )

  out
}
