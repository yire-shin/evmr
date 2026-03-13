#' Profile Likelihood for Return Levels under the rGLO Model
#'
#' Computes and plots the profile log-likelihood for a return level under
#' a stationary r-largest generalized logistic distribution (rGLO) model
#' fitted by \code{\link{rglo.fit}}.
#'
#' @param z An object returned by \code{\link{rglo.fit}}. The fitted model
#'   must represent a stationary model.
#' @param m A return period greater than 1. The profile likelihood is computed
#'   for the corresponding return level exceeded with probability \eqn{1/m}.
#' @param xlow,xup Lower and upper bounds of the return level grid over which
#'   the profile likelihood is evaluated.
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
#' @seealso \code{\link{rglo.fit}}, \code{\link{rglo.rl}}
#' @export
#'
#' @examples
#' \dontrun{
#' x <- rglor(n = 50, r = 2, loc = 10, scale = 2, shape = 0.1)
#' fit <- rglo.fit(x$rmat)
#' rglo.prof(fit, m = 100, xlow = 12, xup = 25)
#' }
rglo.prof <- function(z, m, xlow, xup, conf = 0.95, nint = 100) {

  tol <- .Machine$double.eps^0.5

  if (!inherits(z, "rglo.fit")) {
    warning("'z' does not inherit from class 'rglo.fit'.")
  }

  if (isTRUE(z$trans)) {
    stop("Profile likelihood is only implemented for stationary models.")
  }

  if (!is.numeric(m) || length(m) != 1 || !is.finite(m) || m <= 1) {
    stop("'m' must be a numeric scalar greater than 1.")
  }

  if (!is.numeric(xlow) || length(xlow) != 1 || !is.finite(xlow)) {
    stop("'xlow' must be a finite numeric scalar.")
  }

  if (!is.numeric(xup) || length(xup) != 1 || !is.finite(xup)) {
    stop("'xup' must be a finite numeric scalar.")
  }

  if (xlow >= xup) {
    stop("'xlow' must be smaller than 'xup'.")
  }

  if (!is.numeric(conf) || any(!is.finite(conf)) || any(conf <= 0 | conf >= 1)) {
    stop("'conf' must contain values in (0,1).")
  }

  if (!is.numeric(nint) || length(nint) != 1 || nint < 2) {
    stop("'nint' must be an integer greater than 1.")
  }
  nint <- as.integer(nint)

  if (is.null(z$mle) || length(z$mle) < 3) {
    stop("'z$mle' must contain location, scale, and shape estimates.")
  }

  if (is.null(z$nllh) || !is.finite(z$nllh)) {
    stop("'z$nllh' must be available and finite.")
  }

  cat("If the routine fails, try changing the plotting interval.\n")

  r <- z$r
  p <- 1 - (1 / m)
  v <- numeric(nint)
  x <- seq(xlow, xup, length.out = nint)
  sol <- c(z$mle[2], z$mle[3])   # initial values: sigma, xi
  rl <- rglo.rl(z, year = m)$rl[1]

  rglo.plik <- function(a) {

    sc <- a[1]
    xi <- a[2]

    if (!is.finite(sc) || !is.finite(xi) || sc <= 0 || abs(xi) < tol) {
      return(Inf)
    }

    # Solve return-level equation for mu
    # xp = mu + sc/xi - sc/xi * ((1-p)/p)^xi
    # => mu = xp - sc/xi + sc/xi * ((1-p)/p)^xi
    mu <- xp - sc / xi + (sc / xi) * (((1 - p) / p)^xi)

    if (!is.finite(mu)) {
      return(Inf)
    }

    xmat <- as.matrix(z$data)
    z1 <- drop(xmat[, 1, drop = FALSE])
    zr <- drop(xmat[, r, drop = FALSE])

    ri <- r - seq_len(r)
    cr <- 1 + ri

    if (any(cr <= 0, na.rm = TRUE)) {
      return(Inf)
    }

    y <- 1 - xi * ((xmat - mu) / sc)
    yr <- 1 - xi * ((zr - mu) / sc)

    if (any(y <= 0, na.rm = TRUE) || any(yr <= 0, na.rm = TRUE)) {
      return(Inf)
    }

    f <- 1 + yr^(1 / xi)

    if (any(!is.finite(f), na.rm = TRUE) || any(f <= 0, na.rm = TRUE)) {
      return(Inf)
    }

    # support check
    cond <- (sc / xi) + mu
    if (xi > 0) {
      if (any(cond <= z1, na.rm = TRUE)) {
        return(Inf)
      }
    } else {
      if (any(cond >= zr, na.rm = TRUE)) {
        return(Inf)
      }
    }

    yy <- log(sc) +
      (1 - 1 / xi) * log(y) -
      matrix(log(cr), nrow = nrow(xmat), ncol = ncol(xmat), byrow = TRUE)

    yy <- rowSums(yy, na.rm = TRUE)

    nll <- sum((r + 1) * log(f) + yy, na.rm = TRUE)

    if (!is.finite(nll)) {
      return(Inf)
    }

    nll
  }

  for (i in seq_len(nint)) {
    xp <- x[i]

    opt <- try(
      suppressWarnings(
        stats::optim(sol, rglo.plik, control = list(trace = 0))
      ),
      silent = TRUE
    )

    if (inherits(opt, "try-error") || !is.finite(opt$value)) {
      v[i] <- Inf
    } else {
      v[i] <- opt$value
      sol <- opt$par
    }
  }

  d <- data.frame(x = x, v = -v)

  if (all(!is.finite(d$v))) {
    stop("Profile likelihood failed over the entire grid. Try changing 'xlow' and 'xup'.")
  }

  cent <- which.min(abs(x - rl))

  result_list <- vector("list", length(conf))

  for (i in seq_along(conf)) {
    current_chisq <- -z$nllh - 0.5 * stats::qchisq(conf[i], 1)
    current_diff <- d$v - current_chisq

    left_idx <- seq_len(cent)
    right_idx <- cent:length(current_diff)

    ci1 <- d$x[left_idx][which.min(abs(current_diff[left_idx]))]
    ci2 <- d$x[right_idx][which.min(abs(current_diff[right_idx]))]

    result_list[[i]] <- data.frame(
      year = m,
      returnlevel = rl,
      conf = conf[i],
      ci1 = ci1,
      ci2 = ci2,
      cidiff = ci2 - ci1,
      row.names = NULL
    )
  }

  w_df <- do.call(rbind, result_list)
  print(w_df)

  plot(
    d$x, d$v, type = "l", xlab = "Return level", xlim = c(xlow, xup), las = 1,
    ylab = "Profile log-likelihood",
    main = sprintf("Profile Likelihood for %s-Year Return Level under rGLO (r = %s)", m, z$r)
  )

  cols <- c("blue", "red", "darkgreen", "purple", "orange")
  cols <- rep(cols, length.out = length(conf))

  graphics::abline(v = rl, col = "gray", lty = 2)

  for (i in seq_along(conf)) {
    thr <- -z$nllh - 0.5 * stats::qchisq(conf[i], 1)
    ci  <- result_list[[i]]

    graphics::segments(
      x0 = ci$ci1, x1 = ci$ci2,
      y0 = thr, y1 = thr,
      col = cols[i], lty = 1
    )

    graphics::abline(v = ci$ci1, col = cols[i], lty = 2)
    graphics::abline(v = ci$ci2, col = cols[i], lty = 2)
  }

  legend_labels <- c("Return level", sprintf("Conf: %.2f", conf))
  legend_colors <- c("gray", cols[seq_along(conf)])
  legend_lty <- c(2, rep(1, length(conf)))

  graphics::legend(
    "topright",
    legend = legend_labels,
    col = legend_colors,
    lty = legend_lty,
    bg = "white",
    cex = 0.8
  )

  invisible(w_df)
}
