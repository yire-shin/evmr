#' Fit the Generalized Logistic Distribution to r-Largest Order Statistics
#'
#' Fits the generalized logistic distribution to \eqn{r}-largest order
#' statistics using maximum likelihood estimation. Stationary and
#' non-stationary models are supported through generalized linear modelling
#' of the location, scale, and shape parameters.
#'
#' @param xdat A numeric vector, matrix, or data frame of observations.
#'   Each row should contain decreasing order statistics for a given year
#'   or block. The first column therefore contains block maxima. Only the
#'   first \code{r} columns are used in the fitted model. If \code{r} is
#'   \code{NULL}, all available columns are used. If some rows contain fewer
#'   order statistics than others, missing values should be appended at the
#'   end of the corresponding rows.
#' @param r The number of largest order statistics to use in the fitted model.
#'   If \code{NULL}, all columns of \code{xdat} are used.
#' @param ydat A matrix or data frame of covariates for non-stationary
#'   modelling of the parameters, or \code{NULL} for a stationary model.
#'   The number of rows must match the number of rows of \code{xdat}.
#' @param mul,sigl,shl Integer vectors indicating which columns of
#'   \code{ydat} are used as covariates for the location, scale, and shape
#'   parameters, respectively. Use \code{NULL} for stationary parameters.
#' @param mulink,siglink,shlink Inverse link functions for the location,
#'   scale, and shape parameters, respectively.
#' @param num_inits The number of initial parameter sets used in the
#'   optimization.
#' @param muinit,siginit,shinit Numeric vectors giving initial values for the
#'   location, scale, and shape parameters. If \code{NULL}, default initial
#'   values based on L-moments are used.
#' @param show Logical. If \code{TRUE}, details of the fitted model are printed.
#' @param method Optimization method passed to \code{\link{optim}} for
#'   stationary fits.
#' @param maxit Maximum number of iterations for \code{\link{optim}}.
#' @param ... Additional control arguments passed to the optimizer.
#'
#' @return A list with components including:
#' \item{trans}{Logical; \code{TRUE} if a non-stationary model is fitted.}
#' \item{model}{A list containing \code{mul}, \code{sigl}, and \code{shl}.}
#' \item{link}{A character vector describing the inverse link functions.}
#' \item{conv}{The convergence code returned by the optimizer.}
#' \item{nllh}{The negative log-likelihood evaluated at the fitted parameters.}
#' \item{data}{The data used in the fit.}
#' \item{mle}{The maximum likelihood estimates.}
#' \item{cov}{The estimated covariance matrix when available.}
#' \item{se}{The estimated standard errors when available.}
#' \item{vals}{A matrix containing fitted values of the location, scale,
#' and shape parameters at each observation.}
#' \item{r}{The number of order statistics used in the fitted model.}
#'
#' @references
#'
#' Coles, S. (2001).
#' An Introduction to Statistical Modeling of Extreme Values.
#' Springer.
#'
#' Shin, Y., & Park, J-S. (2024).
#' Generalized logistic model for r-largest order statistics with
#' hydrological application.
#' \emph{Stochastic Environmental Research and Risk Assessment}.
#' \doi{10.1007/s00477-023-02642-7}
#'
#' @seealso \code{\link{optim}}
#' @export
#'
#' @examples
#' x <- rglor(n = 50, r = 2, loc = 10, scale = 2, shape = 0.1)
#' fit <- rglo.fit(x$rmat, num_inits = 5)
#' fit$r
#' fit$mle
rglo.fit <- function(xdat, r = NULL, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL,
                     mulink = identity, siglink = identity, shlink = identity,
                     num_inits = 100, muinit = NULL, siginit = NULL, shinit = NULL,
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

  z$trans <- FALSE

  mu_names <- if (is.null(mul)) "mu" else c("mu0", paste0("mu", seq_len(npmu - 1)))
  sigma_names <- if (is.null(sigl)) "sigma" else c("sigma0", paste0("sigma", seq_len(npsc - 1)))
  xi_names <- if (is.null(shl)) "xi" else c("xi0", paste0("xi", seq_len(npsh - 1)))

  glopar <- lmomco::parglo(lmomco::lmoms(xdat[, 1]))$para

  if (is.null(mul)) {
    mumat <- matrix(1, nrow = nrow(xdat), ncol = 1)
    if (is.null(muinit)) muinit <- glopar[1]
  } else {
    z$trans <- TRUE
    mumat <- cbind(1, as.matrix(ydat[, mul, drop = FALSE]))
    if (is.null(muinit)) muinit <- c(glopar[1], rep(0, length(mul)))
  }

  if (is.null(sigl)) {
    sigmat <- matrix(1, nrow = nrow(xdat), ncol = 1)
    if (is.null(siginit)) siginit <- glopar[2]
  } else {
    z$trans <- TRUE
    sigmat <- cbind(1, as.matrix(ydat[, sigl, drop = FALSE]))
    if (is.null(siginit)) siginit <- c(glopar[2], rep(0, length(sigl)))
  }

  if (is.null(shl)) {
    shmat <- matrix(1, nrow = nrow(xdat), ncol = 1)
    if (is.null(shinit)) shinit <- glopar[3]
  } else {
    z$trans <- TRUE
    shmat <- cbind(1, as.matrix(ydat[, shl, drop = FALSE]))
    if (is.null(shinit)) shinit <- c(glopar[3], rep(0, length(shl)))
  }

  z$model <- list(mul = mul, sigl = sigl, shl = shl)

  z$link <- c(
    mulink = deparse(substitute(mulink)),
    siglink = deparse(substitute(siglink)),
    shlink = deparse(substitute(shlink))
  )

  zr <- drop(xdat[, r, drop = FALSE])

  init <- c(muinit, siginit, shinit)
  names(init) <- c(mu_names, sigma_names, xi_names)

  init_list <- vector("list", num_inits)
  init_list[[1]] <- init

  if (num_inits >= 2) {
    for (i in 2:num_inits) {
      init_list[[i]] <- init + c(
        stats::rnorm(npmu, 0, 0.3),
        abs(stats::rnorm(npsc, 0, 0.2)),
        stats::rnorm(npsh, 0, 0.1)
      )
    }
  }

  ## ----------------------------------
  ## Likelihood r=1
  ## ----------------------------------

  glo.lik <- function(a) {

    out <- tryCatch({

      mu <- drop(mulink(mumat %*% a[1:npmu]))
      sc <- drop(siglink(sigmat %*% a[seq(npmu + 1, length.out = npsc)]))
      xi <- drop(shlink(shmat %*% a[seq(npmu + npsc + 1, length.out = npsh)]))

      if (any(!is.finite(mu)) || any(!is.finite(sc)) || any(!is.finite(xi))) return(BIG)
      if (any(sc <= 0)) return(BIG)
      if (any(abs(xi) < tol)) return(BIG)

      x1 <- drop(xdat[, 1, drop = FALSE])

      y <- 1 - xi * (x1 - mu) / sc
      if (any(y <= 0)) return(BIG)

      f <- 1 + y^(1 / xi)

      if (any(!is.finite(f)) || any(f <= 0)) return(BIG)

      nll <- sum(log(sc) + 2 * log(f) + (1 - 1 / xi) * log(y))

      if (!is.finite(nll)) BIG else nll

    }, error = function(e) BIG)

    out
  }

  ## ----------------------------------
  ## Likelihood r>1
  ## ----------------------------------

  rglo.lik <- function(a) {

    out <- tryCatch({

      mu <- drop(mulink(mumat %*% a[1:npmu]))
      sc <- drop(siglink(sigmat %*% a[seq(npmu + 1, length.out = npsc)]))
      xi <- drop(shlink(shmat %*% a[seq(npmu + npsc + 1, length.out = npsh)]))

      if (any(!is.finite(mu)) || any(!is.finite(sc)) || any(!is.finite(xi))) return(BIG)
      if (any(sc <= 0)) return(BIG)
      if (any(abs(xi) < tol)) return(BIG)

      ri <- r - seq_len(r)
      cr <- 1 + ri

      if (any(cr <= 0)) return(BIG)

      xmat <- xdat

      tmp <- sweep(xmat, 1, mu, "-")
      tmp <- sweep(tmp, 1, sc, "/")

      y <- 1 - xi * tmp
      if (any(y <= 0)) return(BIG)

      zr_std <- (zr - mu) / sc
      f <- 1 + (1 - xi * zr_std)^(1 / xi)

      if (any(!is.finite(f)) || any(f <= 0)) return(BIG)

      yy <- log(sc) +
        (1 - 1 / xi) * log(y) -
        matrix(log(cr), nrow = nrow(xmat), ncol = ncol(xmat), byrow = TRUE)

      yy <- rowSums(yy)

      nll <- sum((r + 1) * log(f) + yy)

      if (!is.finite(nll)) BIG else nll

    }, error = function(e) BIG)

    out
  }

  likfun <- if (r == 1) glo.lik else rglo.lik

  ## ----------------------------------
  ## Optimization
  ## ----------------------------------

  optim_results <- lapply(init_list, function(init_now) {

    val0 <- try(likfun(init_now), silent = TRUE)

    if (inherits(val0, "try-error") || !is.finite(val0)) {
      return(NULL)
    }

    if (!z$trans) {

      try(
        stats::optim(
          init_now, likfun,
          hessian = TRUE,
          method = method,
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
    stop("All optimization attempts failed.")
  }

  ## ----------------------------------
  ## Select best model
  ## ----------------------------------

  optim_value <- data.frame(
    num = seq_along(optim_results),
    nllh = sapply(optim_results, function(res) {
      if (z$trans) min(res$values) else res$value
    }),
    grad = sapply(optim_results, function(res) {
      par_now <- if (z$trans) res$pars else res$par
      if (is.null(par_now) || any(!is.finite(par_now))) return(BIG)
      sum(abs(numDeriv::grad(likfun, par_now)))
    })
  )

  optim_value <- optim_value[is.finite(optim_value$nllh) & optim_value$nllh != BIG, ]

  if (nrow(optim_value) == 0) {
    stop("All valid optimization results failed.")
  }

  optim_value <- optim_value[order(optim_value$grad, optim_value$nllh), ]

  best_result <- optim_results[[optim_value$num[1]]]
  best_par <- if (z$trans) best_result$pars else best_result$par

  mu <- drop(mumat %*% best_par[1:npmu])
  sc <- drop(sigmat %*% best_par[seq(npmu + 1, length.out = npsc)])
  xi <- drop(shmat %*% best_par[seq(npmu + npsc + 1, length.out = npsh)])

  z$r <- r
  z$conv <- best_result$convergence
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

  z$vals <- cbind(mu = mu, sigma = sc, xi = xi)

  if (show) {

    if (z$trans) {
      print(z[c("model", "link")])
    }

    if (!is.null(z$conv) && identical(z$conv, 0L)) {
      print(z[c("conv", "nllh", "mle", "se")])
    }

  }

  class(z) <- "rglo.fit"
  invisible(z)
}
