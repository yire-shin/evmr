#' Entropy Difference Test for rK4D Models
#'
#' Performs the entropy difference (ED) test for selecting the number of
#' order statistics in the r-largest four-parameter kappa distribution
#' (rK4D) model.
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
#'   \item \code{theta}: the estimated parameter vector of the rK4D model
#'   \item \code{ybar}: the sample mean entropy difference
#' }
#'
#' @details
#' This function fits the rK4D model using \code{\link{rk4d.fit}} and then
#' computes the entropy difference test statistic by comparing the fitted
#' likelihood contributions from models with \eqn{r} and \eqn{r-1} order
#' statistics.
#'
#' @references
#' Bader, B., Yan, J., & Zhang, X. (2017).
#' Automated selection of \eqn{r} for the r-largest order statistics approach.
#' \emph{Statistics and Computing}.
#' \doi{10.1007/s11222-016-9697-3}
#'
#' Shin, Y., Park, J.-S., and coauthors (2023).
#' Modeling climate extremes using the four-parameter kappa distribution
#' for r-largest order statistics.
#' \emph{Weather and Climate Extremes}.
#' \doi{10.1016/j.wace.2022.100533}
#'
#' @seealso \code{\link{rk4d.fit}}, \code{\link{rk4dLh}}
#' @export
#'
#' @examples
#' \dontrun{
#' data(bangkok)
#' rk4dEd(bangkok)
#' }
rk4dEd <- function(data) {

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

  y <- rk4d.fit(data, show = FALSE)
  theta1 <- y$mle

  if (length(theta1) < 4 || any(!is.finite(theta1))) {
    stop("Estimated parameter vector is invalid.")
  }

  mu <- theta1[1]
  sc <- theta1[2]
  xi <- theta1[3]
  h  <- theta1[4]

  R  <- ncol(data)
  nr <- nrow(data)

  Diff <- rk4dLh(data[, 1:R, drop = FALSE], theta1) -
    rk4dLh(as.matrix(data[, 1:(R - 1), drop = FALSE]), theta1)

  EstVar <- sum((Diff - mean(Diff))^2) / (nr - 1)

  if (!is.finite(EstVar) || EstVar <= 0) {
    stop("Estimated variance is not positive; the entropy difference statistic cannot be computed.")
  }

  if (h > 0) {
    ar  <- (1 - (R - 1) * h) / h
    ar1 <- (1 - (R - 2) * h) / h

    term1 <- log(sc)
    term2 <- log(1 - (R - 1) * h)
    term3 <- ((1 - R * h) / h) * (digamma(ar) - digamma(ar + R))
    term4 <- ((1 - (R - 1) * h) / h) * (digamma(ar1) - digamma(ar1 + R - 1))
    term5 <- (1 - xi) * ((digamma(R) - digamma(ar + R)) - log(h))

    eta <- -term1 + term2 + term3 - term4 + term5

  } else {
    ar  <- (1 - (R - 1) * h) / h
    ar1 <- (1 - (R - 2) * h) / h

    term1 <- log(sc)
    term2 <- log(1 - (R - 1) * h)
    term3 <- ((1 - R * h) / h) * (digamma((1 / -h) + R) - digamma(1 / -h))
    term4 <- ((1 - (R - 1) * h) / h) * (digamma((1 / -h) + R - 1) - digamma(1 / -h))
    term5 <- (1 - xi) * (digamma(R) - digamma(1 / -h) - log(-h))

    eta <- -term1 + term2 + term3 - term4 + term5
  }

  Diff1 <- sum(Diff) / nr
  Stat <- sqrt(nr) * (Diff1 - eta) / sqrt(EstVar)
  p.value <- 2 * (1 - stats::pnorm(abs(Stat)))

  out <- list(
    statistics = as.numeric(Stat),
    p.value = as.numeric(p.value),
    theta = theta1,
    ybar = as.numeric(Diff1)
  )

  out
}
