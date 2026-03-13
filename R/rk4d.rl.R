#' Return Levels for the Four-Parameter Kappa Distribution
#'
#' Computes return levels and their standard errors for a stationary
#' four-parameter kappa model fitted by \code{\link{rk4d.fit}}.
#'
#' @param z An object returned by \code{\link{rk4d.fit}}. The fitted model
#'   should represent a stationary model.
#' @param year A numeric vector of return periods for which return levels
#'   are to be computed.
#' @param show Logical. If \code{TRUE}, the estimated return levels and
#'   their standard errors are printed.
#'
#' @return The input object \code{z} with two additional components:
#' \itemize{
#'   \item \code{rl}: a numeric vector of estimated return levels
#'   \item \code{rlse}: a numeric vector of standard errors of the estimated return levels
#' }
#'
#' @details
#' For a return period \eqn{T}, the return level is defined as the quantile
#' exceeded with probability \eqn{1/T}. Under the four-parameter kappa
#' distribution, the return level is
#' \deqn{x_T = \mu + \frac{\sigma}{\xi} - \frac{\sigma}{\xi}
#' \left(\frac{1-(1-1/T)^h}{h}\right)^\xi,}
#' and standard errors are obtained using the delta method.
#'
#' @seealso \code{\link{rk4d.fit}}, \code{\link{rk4d.prof}}
#' @export
#'
#' @examples
#' x <- rk4dr(n = 50, r = 2, loc = 10, scale = 2, shape1 = 0.1, shape2 = 0.1)
#' fit <- rk4d.fit(x$rmat, num_inits = 5)
#' out <- rk4d.rl(fit, year = c(20, 50, 100))
#' out$rl
#' out$rlse
rk4d.rl <- function(z, year = c(20, 50, 100, 200), show = FALSE) {

  if (!inherits(z, "rk4d.fit")) {
    warning("'z' does not inherit from class 'rk4d.fit'.")
  }

  if (isTRUE(z$trans)) {
    stop("Return levels are only implemented for stationary models.")
  }

  if (!is.numeric(year) || any(!is.finite(year)) || any(year <= 1)) {
    stop("'year' must be a numeric vector with values greater than 1.")
  }

  if (is.null(z$mle) || length(z$mle) < 4) {
    stop("'z$mle' must contain location, scale, and two shape estimates.")
  }

  if (is.null(z$cov) || all(is.na(z$cov))) {
    stop("Covariance matrix is not available in 'z$cov'.")
  }

  mu  <- z$mle[1]
  sig <- z$mle[2]
  xi  <- z$mle[3]
  h   <- z$mle[4]

  if (!is.finite(sig) || sig <= 0) {
    stop("The fitted scale parameter must be positive.")
  }

  if (!is.finite(xi) || abs(xi) < .Machine$double.eps^0.5) {
    stop("The fitted first shape parameter is too close to 0 for 'rk4d.rl()'.")
  }

  if (!is.finite(h) || abs(h) < .Machine$double.eps^0.5) {
    stop("The fitted second shape parameter is too close to 0 for 'rk4d.rl()'.")
  }

  del <- matrix(nrow = 4, ncol = length(year))

  p  <- 1 - (1 / year)
  yp <- (1 - p^h) / h

  if (any(!is.finite(yp)) || any(yp <= 0)) {
    stop("Invalid return level calculation: transformed term must be positive.")
  }

  del[1, ] <- 1
  del[2, ] <- (1 / xi) - (1 / xi) * (yp^xi)
  del[3, ] <- -(sig / xi) * log(yp) * (yp^xi) + (sig / xi^2) * (yp^xi) - (sig / xi^2)
  del[4, ] <- sig * (yp^(xi - 1)) * ((log(p) * (p^h) * h + (1 - p^h)) / h^2)

  del.t <- t(del)

  z$rl <- mu + (sig / xi) - (sig / xi) * (yp^xi)
  z$rlse <- sqrt(diag(del.t %*% z$cov %*% del))

  names(z$rl) <- paste0(as.character(year), "y")
  names(z$rlse) <- paste0(as.character(year), "y")

  if (show) {
    print(z[c("rl", "rlse")])
  }

  invisible(z)
}
