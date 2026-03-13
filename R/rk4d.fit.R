#' Fit the Four-Parameter Kappa Distribution to r-Largest Order Statistics
#'
#' Fits the four-parameter kappa distribution to \eqn{r}-largest order
#' statistics using maximum likelihood estimation. Stationary and
#' non-stationary models are supported through generalized linear modelling
#' of the location, scale, and two shape parameters.
#'
#' @param xdat A numeric vector, matrix, or data frame of observations.
#'   Each row should contain decreasing order statistics for a given year
#'   or block. The first column therefore contains block maxima. Only the
#'   first \code{r} columns are used in the fitted model. If \code{r} is
#'   \code{NULL}, all available columns are used.
#' @param r The number of largest order statistics to use in the fitted model.
#'   If \code{NULL}, all columns of \code{xdat} are used.
#' @param penk Optional penalty for the first shape parameter. Supported values
#'   include \code{"CD"} and \code{"MS"}.
#' @param penh Optional penalty for the second shape parameter. Supported values
#'   include \code{"MS"} and \code{"MSa"}.
#' @param ydat A matrix or data frame of covariates for non-stationary
#'   modelling of the parameters, or \code{NULL} for a stationary model.
#'   The number of rows must match the number of rows of \code{xdat}.
#' @param mul,sigl,shl,hl Integer vectors indicating which columns of
#'   \code{ydat} are used as covariates for the location, scale, first shape,
#'   and second shape parameters, respectively.
#' @param mulink,siglink,shlink,hlink Inverse link functions for the location,
#'   scale, first shape, and second shape parameters, respectively.
#' @param num_inits The number of initial parameter sets used in the
#'   optimization.
#' @param muinit,siginit,shinit,hinit Numeric vectors giving initial values for
#'   the location, scale, first shape, and second shape parameters. If
#'   \code{NULL}, default initial values based on L-moments are used.
#' @param show Logical. If \code{TRUE}, details of the fitted model are printed.
#' @param method Optimization method passed to \code{\link{optim}} for
#'   stationary fits.
#' @param maxit Maximum number of iterations for \code{\link{optim}}.
#' @param ... Additional control arguments passed to the optimizer.
#'
#' @return A list with components including:
#' \item{trans}{Logical; \code{TRUE} if a non-stationary model is fitted.}
#' \item{model}{A list containing \code{mul}, \code{sigl}, \code{shl}, and \code{hl}.}
#' \item{link}{A character vector describing the inverse link functions.}
#' \item{conv}{The convergence code returned by the optimizer.}
#' \item{nllh}{The negative log-likelihood evaluated at the fitted parameters.}
#' \item{data}{The data used in the fit.}
#' \item{mle}{The maximum likelihood estimates.}
#' \item{cov}{The estimated covariance matrix when available.}
#' \item{se}{The estimated standard errors when available.}
#' \item{vals}{A matrix containing fitted values of the location, scale,
#' first shape, and second shape parameters at each observation.}
#' \item{r}{The number of order statistics used in the fitted model.}
#'
#' @references
#'
#' Hosking, J. R. M. (1994).
#' The four-parameter kappa distribution.
#' \emph{IBM Journal of Research and Development}, 38(3), 251–258.
#'
#' Martins, E. S., & Stedinger, J. R. (2000).
#' Generalized maximum-likelihood generalized extreme-value quantile
#' estimators for hydrologic data.
#' \emph{Water Resources Research}, 36(3), 737–744.
#' \doi{10.1029/1999WR900330}
#'
#' Coles, S., & Dixon, M. (1999).
#' Likelihood-based inference for extreme value models.
#' \emph{Extremes}, 2(1), 5–23.
#' \doi{10.1023/A:1009905222644}
#'
#' Coles, S. (2001).
#' An Introduction to Statistical Modeling of Extreme Values.
#' Springer.
#'
#' Shin, Y., & Park, J.-S. (2023).
#' Modeling climate extremes using the four-parameter kappa distribution
#' for r-largest order statistics.
#' \emph{Weather and Climate Extremes}.
#' \doi{10.1016/j.wace.2022.100533}
#'
#' @seealso \code{\link{optim}}
#' @export
#'
#' @examples
#' x <- rk4dr(n = 50, r = 2, loc = 10, scale = 2, shape1 = 0.1, shape2 = 0.1)
#' fit <- rk4d.fit(x$rmat, num_inits = 5)
#' fit$r
#' fit$mle
rk4d.fit <- function(xdat, r = NULL, penk = NULL, penh = NULL,
                     ydat = NULL, mul = NULL, sigl = NULL, shl = NULL, hl = NULL,
                     mulink = identity, siglink = identity,
                     shlink = identity, hlink = identity,
                     num_inits = 100,
                     muinit = NULL, siginit = NULL, shinit = NULL, hinit = NULL,
                     show = TRUE, method = "Nelder-Mead", maxit = 10000, ...) {

  z <- list()
  tol <- .Machine$double.eps^0.5
  BIG <- 1e6

  ## ----------------------------------
  ## Input handling
  ## ----------------------------------
  if (is.vector(xdat)) {
    xdat <- matrix(xdat, ncol = 1)
  } else {
    xdat <- as.matrix(xdat)
  }

  if (!is.numeric(xdat)) {
    stop("'xdat' must be numeric.")
  }

  if (is.null(r)) {
    r <- ncol(xdat)
  } else {
    if (!is.numeric(r) || length(r) != 1 || r < 1 || r != as.integer(r)) {
      stop("'r' must be a positive integer.")
    }
    r <- as.integer(r)
    if (r > ncol(xdat)) {
      stop("'r' cannot exceed the number of columns in 'xdat'.")
    }
    xdat <- xdat[, 1:r, drop = FALSE]
  }

  if (!is.null(ydat)) {
    ydat <- as.data.frame(ydat)
    if (nrow(ydat) != nrow(xdat)) {
      stop("'ydat' must have the same number of rows as 'xdat'.")
    }
  }

  if (!is.numeric(num_inits) || length(num_inits) != 1 || num_inits < 1) {
    stop("'num_inits' must be a positive integer.")
  }
  num_inits <- as.integer(num_inits)

  npmu <- length(mul) + 1
  npsc <- length(sigl) + 1
  npsh <- length(shl) + 1
  nph  <- length(hl) + 1

  z$trans <- FALSE

  mu_names <- if (is.null(mul)) "mu" else c("mu0", paste0("mu", seq_len(npmu - 1)))
  sigma_names <- if (is.null(sigl)) "sigma" else c("sigma0", paste0("sigma", seq_len(npsc - 1)))
  xi_names <- if (is.null(shl)) "xi" else c("xi0", paste0("xi", seq_len(npsh - 1)))
  h_names <- if (is.null(hl)) "h" else c("h0", paste0("h", seq_len(nph - 1)))

  kappar <- lmomco::parkap(lmomco::lmoms(xdat[, 1]))$para

  if (is.null(mul)) {
    mumat <- matrix(1, nrow = nrow(xdat), ncol = 1)
    if (is.null(muinit)) muinit <- kappar[1]
  } else {
    z$trans <- TRUE
    mumat <- cbind(1, as.matrix(ydat[, mul, drop = FALSE]))
    if (is.null(muinit)) muinit <- c(kappar[1], rep(0, length(mul)))
  }

  if (is.null(sigl)) {
    sigmat <- matrix(1, nrow = nrow(xdat), ncol = 1)
    if (is.null(siginit)) siginit <- kappar[2]
  } else {
    z$trans <- TRUE
    sigmat <- cbind(1, as.matrix(ydat[, sigl, drop = FALSE]))
    if (is.null(siginit)) siginit <- c(kappar[2], rep(0, length(sigl)))
  }

  if (is.null(shl)) {
    shmat <- matrix(1, nrow = nrow(xdat), ncol = 1)
    if (is.null(shinit)) shinit <- kappar[3]
  } else {
    z$trans <- TRUE
    shmat <- cbind(1, as.matrix(ydat[, shl, drop = FALSE]))
    if (is.null(shinit)) shinit <- c(kappar[3], rep(0, length(shl)))
  }

  if (is.null(hl)) {
    hmat <- matrix(1, nrow = nrow(xdat), ncol = 1)
    if (is.null(hinit)) hinit <- kappar[4]
  } else {
    z$trans <- TRUE
    hmat <- cbind(1, as.matrix(ydat[, hl, drop = FALSE]))
    if (is.null(hinit)) hinit <- c(kappar[4], rep(0, length(hl)))
  }

  z$model <- list(mul = mul, sigl = sigl, shl = shl, hl = hl)
  z$link <- c(
    mulink = deparse(substitute(mulink)),
    siglink = deparse(substitute(siglink)),
    shlink = deparse(substitute(shlink)),
    hlink = deparse(substitute(hlink))
  )

  zr <- drop(xdat[, r, drop = FALSE])

  init <- c(muinit, siginit, shinit, hinit)
  names(init) <- c(mu_names, sigma_names, xi_names, h_names)

  init_list <- vector("list", num_inits)
  init_list[[1]] <- init

  if (num_inits >= 2) {
    for (i in 2:num_inits) {
      init_list[[i]] <- init + c(
        stats::rnorm(npmu, mean = 0, sd = 0.3),
        abs(stats::rnorm(npsc, mean = 0, sd = 0.2)),
        stats::rnorm(npsh, mean = 0, sd = 0.05),
        stats::rnorm(nph, mean = 0, sd = 0.05)
      )
    }
  }

  ## ----------------------------------
  ## Penalty helpers
  ## ----------------------------------
  get_penalty_k <- function(xi, penk, r_now) {
    if (is.null(penk)) return(0)

    if (penk == "CD") {
      if (xi >= 0) {
        p_k <- 1
      } else if (xi > -1 && xi < 0) {
        p_k <- exp(-((1 / (1 + xi)) - 1))
      } else {
        return(BIG)
      }
      return(-r_now * log(p_k))
    }

    if (penk == "MS") {
      if (xi >= 0.5 || xi <= -0.5) return(BIG)
      p_k <- ((0.5 + xi)^(6 - 1)) * ((0.5 - xi)^(9 - 1)) / beta(6, 9)
      return(-r_now * log(p_k))
    }

    0
  }

  get_penalty_h <- function(h, penh, r_now) {
    if (is.null(penh)) return(0)

    if (r_now == 1) {
      if (penh == "MS") {
        if (h <= -0.5 || h >= 0.5) return(BIG)
        Bef <- function(x) ((0.5 + x)^(6 - 1)) * ((0.5 - x)^(9 - 1))
        Be <- stats::integrate(Bef, lower = -0.5, upper = 0.5)$value
        p_h <- ((0.5 + h)^(6 - 1)) * ((0.5 - h)^(9 - 1)) / Be
        return(-r_now * log(p_h))
      }
      if (penh == "MSa") {
        if (h <= -1.2 || h >= 1.2) return(BIG)
        Bef <- function(x) ((1.2 + x)^(6 - 1)) * ((1.2 - x)^(9 - 1))
        Be <- stats::integrate(Bef, lower = -1.2, upper = 1.2)$value
        p_h <- ((1.2 + h)^(6 - 1)) * ((1.2 - h)^(9 - 1)) / Be
        return(-r_now * log(p_h))
      }
      return(0)
    }

    B <- 1 / (r_now - 1)

    if (penh == "MS") {
      if (h <= -0.5 || h >= B) return(BIG)
      Bef <- function(x) ((0.5 + x)^(6 - 1)) * ((B - x)^(9 - 1))
      Be <- stats::integrate(Bef, lower = -0.5, upper = B)$value
      p_h <- ((0.5 + h)^(6 - 1)) * ((B - h)^(9 - 1)) / Be
      return(-r_now * log(p_h))
    }

    if (penh == "MSa") {
      if (h <= -1.2 || h >= B) return(BIG)
      Bef <- function(x) ((1.2 + x)^(6 - 1)) * ((B - x)^(9 - 1))
      Be <- stats::integrate(Bef, lower = -1.2, upper = B)$value
      p_h <- ((1.2 + h)^(6 - 1)) * ((B - h)^(9 - 1)) / Be
      return(-r_now * log(p_h))
    }

    0
  }

  ## ----------------------------------
  ## Likelihood for r = 1
  ## ----------------------------------
  k4d.lik <- function(a) {
    out <- tryCatch({
      mu <- drop(mulink(mumat %*% a[1:npmu]))
      sc <- drop(siglink(sigmat %*% a[seq(npmu + 1, length.out = npsc)]))
      xi <- drop(shlink(shmat %*% a[seq(npmu + npsc + 1, length.out = npsh)]))
      h  <- drop(hlink(hmat %*% a[seq(npmu + npsc + npsh + 1, length.out = nph)]))

      if (any(!is.finite(mu)) || any(!is.finite(sc)) || any(!is.finite(xi)) || any(!is.finite(h))) return(BIG)
      if (any(sc <= 0, na.rm = TRUE)) return(BIG)
      if (any(abs(xi) < tol, na.rm = TRUE) || any(abs(h) < tol, na.rm = TRUE)) return(BIG)

      x1 <- drop(xdat[, 1, drop = FALSE])

      y <- 1 - xi * ((x1 - mu) / sc)
      if (any(y <= 0, na.rm = TRUE)) return(BIG)

      f <- 1 - h * y^(1 / xi)
      if (any(!is.finite(f), na.rm = TRUE) || any(f <= 0, na.rm = TRUE)) return(BIG)
      if (any(f^(1 / h) > 1, na.rm = TRUE)) return(BIG)

      penalty <- get_penalty_k(xi[1], penk, 1) + get_penalty_h(h[1], penh, 1)
      if (!is.finite(penalty) || penalty >= BIG) return(BIG)

      nll <- sum(
        log(sc) +
          ((h - 1) / h) * log(f) +
          (1 - 1 / xi) * log(y),
        na.rm = TRUE
      ) + penalty

      if (!is.finite(nll)) BIG else nll
    }, error = function(e) BIG)

    out
  }

  ## ----------------------------------
  ## Likelihood for r > 1
  ## ----------------------------------
  rk4d.lik <- function(a) {
    out <- tryCatch({
      mu <- drop(mulink(mumat %*% a[1:npmu]))
      sc <- drop(siglink(sigmat %*% a[seq(npmu + 1, length.out = npsc)]))
      xi <- drop(shlink(shmat %*% a[seq(npmu + npsc + 1, length.out = npsh)]))
      h  <- drop(hlink(hmat %*% a[seq(npmu + npsc + npsh + 1, length.out = nph)]))

      if (any(!is.finite(mu)) || any(!is.finite(sc)) || any(!is.finite(xi)) || any(!is.finite(h))) return(BIG)
      if (any(sc <= 0, na.rm = TRUE)) return(BIG)
      if (any(abs(xi) < tol, na.rm = TRUE) || any(abs(h) < tol, na.rm = TRUE)) return(BIG)

      ri <- r - seq_len(r)
      cr <- 1 - ri * h[1]
      if (any(cr <= 0, na.rm = TRUE)) return(BIG)

      xmat <- xdat
      tmp <- sweep(xmat, 1, mu, "-")
      tmp <- sweep(tmp, 1, sc, "/")

      y <- 1 - xi * tmp
      if (any(y <= 0, na.rm = TRUE)) return(BIG)

      yr <- 1 - xi * ((zr - mu) / sc)
      if (any(yr <= 0, na.rm = TRUE)) return(BIG)

      f <- 1 - h[1] * yr^(1 / xi)
      if (any(!is.finite(f), na.rm = TRUE) || any(f <= 0, na.rm = TRUE)) return(BIG)
      if (any(f^(1 / h[1]) > 1, na.rm = TRUE)) return(BIG)

      if (max(h, na.rm = TRUE) > min(1 / yr^(1 / xi), na.rm = TRUE)) return(BIG)
      if (r >= 2 && min(h, na.rm = TRUE) > (1 / (r - 1))) return(BIG)

      yy <- log(sc) +
        (1 - 1 / xi) * log(y) -
        matrix(log(cr), nrow = nrow(xmat), ncol = ncol(xmat), byrow = TRUE)

      yy <- rowSums(yy, na.rm = TRUE)

      penalty <- get_penalty_k(xi[1], penk, r) + get_penalty_h(h[1], penh, r)
      if (!is.finite(penalty) || penalty >= BIG) return(BIG)

      nll <- sum(((r * h[1] - 1) / h[1]) * log(f) + yy, na.rm = TRUE) + penalty

      if (!is.finite(nll)) BIG else nll
    }, error = function(e) BIG)

    out
  }

  likfun <- if (r == 1) k4d.lik else rk4d.lik

  ## ----------------------------------
  ## Optimization
  ## ----------------------------------
  optim_results <- lapply(init_list, function(init_now) {
    if (!z$trans) {
      try(
        stats::optim(
          init_now, likfun, hessian = TRUE, method = method,
          control = list(maxit = maxit, trace = 0)
        ),
        silent = TRUE
      )
    } else {
      try(
        suppressWarnings(
          Rsolnp::solnp(
            pars = init_now,
            fun = likfun,
            control = list(trace = 0)
          )
        ),
        silent = TRUE
      )
    }
  })

  optim_results <- Filter(function(res) {
    !is.null(res) && !inherits(res, "try-error")
  }, optim_results)

  if (length(optim_results) == 0) {
    stop("All optimization attempts failed. Try num_inits = 1 or different initial values.")
  }

  ## ----------------------------------
  ## Score optimization results
  ## ----------------------------------
  optim_value <- data.frame(
    num = seq_along(optim_results),
    nllh = sapply(optim_results, function(res) {
      if (z$trans) min(res$values) else res$value
    }),
    grad = sapply(optim_results, function(res) {
      par_now <- if (z$trans) res$pars else res$par
      if (is.null(par_now) || any(!is.finite(par_now))) return(Inf)
      sum(abs(numDeriv::grad(likfun, par_now)))
    })
  )

  optim_value <- optim_value[is.finite(optim_value$nllh) & optim_value$nllh != BIG, , drop = FALSE]

  if (nrow(optim_value) == 0) {
    stop("All valid optimization results failed. Try num_inits = 1.")
  }

  optim_value <- optim_value[order(optim_value$grad, optim_value$nllh), , drop = FALSE]

  best_result <- optim_results[[optim_value$num[1]]]
  best_par <- if (z$trans) best_result$pars else best_result$par

  mu <- drop(mulink(mumat %*% best_par[1:npmu]))
  sc <- drop(siglink(sigmat %*% best_par[seq(npmu + 1, length.out = npsc)]))
  xi <- drop(shlink(shmat %*% best_par[seq(npmu + npsc + 1, length.out = npsh)]))
  h  <- drop(hlink(hmat %*% best_par[seq(npmu + npsc + npsh + 1, length.out = nph)]))

  z$r <- r
  z$conv <- if (z$trans) best_result$convergence else best_result$convergence
  z$nllh <- if (z$trans) min(best_result$values) else best_result$value
  z$data <- xdat
  z$mle <- best_par

  if (!z$trans && !is.null(best_result$hessian)) {
    z$cov <- tryCatch(solve(best_result$hessian), error = function(e) NA)
    z$se <- if (is.matrix(z$cov)) sqrt(diag(z$cov)) else NA
  } else {
    z$cov <- NA
    z$se <- NA
  }

  z$vals <- cbind(mu = mu, sigma = sc, xi = xi, h = h)

  if (show) {
    if (z$trans) {
      print(z[c("model", "link")])
    }
    if (!is.null(z$conv) && identical(z$conv, 0L)) {
      print(z[c("conv", "nllh", "mle", "se")])
    }
  }

  class(z) <- "rk4d.fit"
  invisible(z)
}
