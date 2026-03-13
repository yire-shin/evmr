#' Entropy Difference Test for rGGD Models
#'
#' Performs the entropy difference (ED) test for selecting the number of
#' order statistics in the r-largest generalized Gumbel distribution (rGGD)
#' model.
#'
#' The test compares the entropy of models fitted with \eqn{r} and
#' \eqn{r-1} order statistics and evaluates whether the additional order
#' statistic provides significant information.
#'
#' @param data A numeric matrix or data frame containing the r-largest
#'   order statistics. Each row represents one block or observation, and
#'   columns must be ordered from largest to smallest.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{statistics}: the entropy difference test statistic
#'   \item \code{p.value}: the two-sided p-value
#'   \item \code{theta}: the estimated parameter vector of the rGGD model
#'   \item \code{ybar}: the sample mean entropy difference
#' }
#'
#' @details
#' This function fits the rGGD model using \code{\link{rggd.fit}} and then
#' computes the entropy difference test statistic by comparing the fitted
#' likelihood contributions from models with \eqn{r} and \eqn{r-1} order
#' statistics.
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
#' @export
#'
#' @examples
#' \dontrun{
#' data(bangkok)
#' rggdEd(bangkok)
#' }
rggdEd <- function(data) {

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

  y <- rggd.fit(data, show = FALSE)
  theta1 <- y$mle

  mu <- theta1[1]
  sc <- theta1[2]
  h  <- theta1[3]

  R  <- ncol(data)
  nr <- nrow(data)

  Diff <- rggdLh(data[, 1:R, drop = FALSE], theta1) -
    rggdLh(as.matrix(data[, 1:(R - 1), drop = FALSE]), theta1)

  EstVar <- sum((Diff - mean(Diff))^2) / (nr - 1)

  if (!is.finite(EstVar) || EstVar <= 0) {
    stop("Estimated variance is not positive; the entropy difference statistic cannot be computed.")
  }

  ri  <- R - seq_len(R)
  cr  <- prod(1 - ri * h)

  ri1 <- (R - 1) - seq_len(R - 1)
  cr1 <- prod(1 - ri1 * h)

  if (h > 0) {
    ar  <- (1 - (R - 1) * h) / h
    ar1 <- (1 - (R - 2) * h) / h

    term1 <- log(sc)
    term2 <- log(1 - (R - 1) * h)

    term3 <- ((1 - R * h) / h) *
      (digamma(ar) - digamma(ar + R))

    term4 <- ((1 - (R - 1) * h) / h) *
      (digamma(ar1) - digamma(ar1 + R - 1))

    term5 <- digamma(R) - digamma(ar + R) - log(h)

    eta <- -term1 + term2 + term3 - term4 + term5

  } else {
    ar  <- (1 - (R - 1) * h) / h
    ar1 <- (1 - (R - 2) * h) / h

    term1 <- log(sc)
    term2 <- log(1 - (R - 1) * h)

    term3 <- ((1 - R * h) / h) *
      (digamma((1 / -h) + R) - digamma(1 / -h))

    term4 <- ((1 - (R - 1) * h) / h) *
      (digamma((1 / -h) + R - 1) - digamma(1 / -h))

    term5 <- digamma(R) - digamma(1 / -h) - log(-h)

    eta <- -term1 + term2 + term3 - term4 + term5
  }

  Diff1 <- sum(Diff) / nr
  Diff  <- sqrt(nr) * (Diff1 - eta) / sqrt(EstVar)
  p.value <- 2 * (1 - stats::pnorm(abs(Diff)))

  out <- list(
    statistics = as.numeric(Diff),
    p.value    = as.numeric(p.value),
    theta      = theta1,
    ybar       = as.numeric(Diff1)
  )

  out
}
