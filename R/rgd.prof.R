#' Profile Likelihood for Return Levels under the rGD Model
#'
#' Computes and plots the profile log-likelihood for a return level under
#' a stationary r-largest Gumbel distribution model fitted by \code{rgd.fit()}.
#'
#' @param z An object returned by \code{\link{rgd.fit}}. The fitted model must
#'   be stationary.
#' @param m A return period greater than 1. The profile likelihood is computed
#'   for the corresponding return level exceeded with probability \eqn{1/m}.
#' @param xlow,xup The lower and upper bounds of the return level grid over
#'   which the profile likelihood is evaluated.
#' @param conf A numeric vector of confidence levels for profile likelihood
#'   confidence intervals.
#' @param nint The number of grid points used to evaluate the profile likelihood.
#'
#' @return A data frame containing the return period, estimated return level,
#'   confidence level, lower confidence limit, upper confidence limit, and
#'   interval width. A profile likelihood plot is also produced.
#'
#' @details
#' The function evaluates the profile log-likelihood over a grid of return
#' level values and plots the resulting curve. Horizontal and vertical lines
#' are added to indicate profile likelihood confidence intervals for the
#' confidence levels specified in \code{conf}.
#'
#' @seealso \code{\link{rgd.fit}}, \code{\link{rgd.rl}}
#' @export
#'
#' @examples
#' \dontrun{
#' x <- rgdr(n = 50, r = 2, loc = 10, scale = 2)
#' fit <- rgd.fit(x$rmat)
#' rgd.prof(fit, m = 100, xlow = 12, xup = 25)
#' }
rgd.prof <- function(z, m, xlow, xup, conf = 0.95, nint = 100) {

  if (m <= 1)
    stop("'m' must be greater than 1.")

  if (!is.numeric(xlow) || !is.numeric(xup) || length(xlow) != 1 || length(xup) != 1) {
    stop("'xlow' and 'xup' must be numeric scalars.")
  }

  if (xlow >= xup) {
    stop("'xlow' must be smaller than 'xup'.")
  }

  if (any(conf <= 0 | conf >= 1)) {
    stop("'conf' must contain values in (0,1).")
  }

  if (!is.numeric(nint) || length(nint) != 1 || nint < 2) {
    stop("'nint' must be an integer greater than 1.")
  }

  if (!inherits(z, "rgd.fit")) {
    warning("'z' does not inherit from class 'rgd.fit'.")
  }

  if (isTRUE(z$trans)) {
    stop("Profile likelihood is only implemented for stationary models.")
  }

  cat("If the routine fails, try changing the plotting interval.\n")

  r <- z$r
  p <- 1 - (1 / m)
  v <- numeric(nint)
  x <- seq(xlow, xup, length.out = nint)
  sol <- c(z$mle[2])
  rl <- rgd.rl(z, m)$rl

  rgd.plik <- function(a) {
    zr <- as.matrix(z$data[, r, drop = FALSE])

    mu <- xp + a[1] * log(-log(p))

    if (any(a[1] <= 0, na.rm = TRUE)) return(1e6)

    y <- (as.matrix(z$data) - mu) / a[1]
    y <- y + log(a[1])
    y <- rowSums(y, na.rm = TRUE)

    sum(exp(-(zr - mu) / a[1]) + y)
  }

  for (i in seq_len(nint)) {
    xp <- x[i]
    opt <- try(suppressWarnings(stats::optim(sol, rgd.plik, control = list(trace = 0))),
               silent = TRUE)
    if (inherits(opt, "try-error")) {
      v[i] <- 1e6
    } else {
      v[i] <- opt$value
      sol <- opt$par
    }
  }

  d <- data.frame(x = x, v = -v)

  cent <- which.min(abs(round(x) - round(rl)))

  result_list <- vector("list", length(conf))

  for (i in seq_along(conf)) {
    current_chisq <- -z$nllh - 0.5 * stats::qchisq(conf[i], 1)
    current_diff <- d$v - current_chisq

    ci1 <- d$x[which.min(abs(current_diff[1:cent]))]
    ci2 <- d$x[cent:length(current_diff)][which.min(abs(current_diff[cent:length(current_diff)]))]
    cidiff <- ci2 - ci1

    result_list[[i]] <- data.frame(
      year = m,
      returnlevel = rl,
      conf = conf[i],
      ci1 = ci1,
      ci2 = ci2,
      cidiff = cidiff,
      row.names = NULL
    )
  }

  w_df <- do.call(rbind, result_list)
  print(w_df)

  plot(
    d$x, d$v, type = "l", xlab = "Return level", xlim = c(xlow, xup), las = 1,
    ylab = "Profile log-likelihood",
    main = sprintf("Profile Likelihood for %s-Year Return Level under rGD (r = %s)", m, z$r)
  )

  col <- c("blue", "red", "darkgreen")
  threshold <- vector("list", length(conf))
  ci <- vector("list", length(conf))

  graphics::abline(v = rl, col = "gray", lty = 2)

  for (i in seq_along(conf)) {
    threshold[[i]] <- -z$nllh - 0.5 * stats::qchisq(conf[i], 1)
    ci[[i]] <- result_list[[i]]

    graphics::segments(
      x0 = ci[[i]]$ci1, x1 = ci[[i]]$ci2,
      y0 = threshold[[i]], y1 = threshold[[i]],
      col = col[i], lty = 1
    )

    graphics::abline(v = ci[[i]]$ci1, col = col[i], lty = 2)
    graphics::abline(v = ci[[i]]$ci2, col = col[i], lty = 2)
  }

  legend_labels <- c("Return level", sprintf("Conf: %.2f", conf))
  legend_colors <- c("gray", col[seq_along(conf)])

  graphics::legend(
    "topright",
    legend = legend_labels,
    col = legend_colors,
    lty = c(2, rep(1, length(conf))),
    bg = "white",
    cex = 0.8
  )

  invisible(w_df)
}
