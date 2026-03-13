#' Fit the Logistic Distribution to r-Largest Order Statistics
#'
#' Fits the logistic distribution to \eqn{r}-largest order statistics
#' using maximum likelihood estimation. Stationary and non-stationary models
#' are supported through generalized linear modelling of the location and
#' scale parameters.
#'
#' @param xdat A numeric vector, matrix, or data frame of observations.
#'   Each row should contain decreasing order statistics for a given year
#'   or block. The first column therefore contains block maxima. Only the
#'   first \code{r} columns are used in the fitted model. If \code{r} is
#'   \code{NULL}, all available columns are used.
#' @param r The number of largest order statistics to use in the fitted model.
#'   If \code{NULL}, all columns of \code{xdat} are used.
#' @param ydat A matrix or data frame of covariates for non-stationary
#'   modelling of the parameters, or \code{NULL} for a stationary model.
#'   The number of rows must match the number of rows of \code{xdat}.
#' @param mul,sigl Integer vectors indicating which columns of
#'   \code{ydat} are used as covariates for the location and scale
#'   parameters, respectively.
#' @param mulink,siglink Inverse link functions for the location and
#'   scale parameters, respectively.
#' @param num_inits The number of initial parameter sets used in the
#'   optimization.
#' @param muinit,siginit Numeric vectors giving initial values for the
#'   location and scale parameters. If \code{NULL}, default initial
#'   values based on L-moments are used.
#' @param show Logical. If \code{TRUE}, details of the fitted model are printed.
#' @param method Optimization method passed to \code{\link{optim}} for
#'   stationary fits.
#' @param maxit Maximum number of iterations for \code{\link{optim}}.
#' @param ... Additional control arguments passed to the optimizer.
#'
#' @return A list with components including:
#' \itemize{
#'   \item \code{trans}: logical; \code{TRUE} if a non-stationary model is fitted
#'   \item \code{model}: a list containing \code{mul} and \code{sigl}
#'   \item \code{link}: a character vector describing the inverse link functions
#'   \item \code{conv}: the convergence code returned by the optimizer
#'   \item \code{nllh}: the negative log-likelihood evaluated at the fitted parameters
#'   \item \code{data}: the data used in the fit
#'   \item \code{mle}: the maximum likelihood estimates
#'   \item \code{cov}: the estimated covariance matrix when available
#'   \item \code{se}: the estimated standard errors when available
#'   \item \code{vals}: a matrix containing fitted values of the location and scale
#'   \item \code{r}: the number of order statistics used in the fitted model
#' }
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
#' x <- rldr(n = 50, r = 2, loc = 10, scale = 2)
#' fit <- rld.fit(x$rmat, num_inits = 5)
#' fit$r
#' fit$mle
rld.fit <- function(xdat, r = NULL, ydat = NULL, mul = NULL, sigl = NULL,
                    mulink = identity, siglink = identity, num_inits = 100,
                    muinit = NULL, siginit = NULL, show = TRUE,
                    method = "Nelder-Mead", maxit = 10000, ...) {

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

  z$trans <- FALSE

  mu_names <- if (is.null(mul)) "mu" else c("mu0", paste0("mu", seq_len(npmu - 1)))
  sigma_names <- if (is.null(sigl)) "sigma" else c("sigma0", paste0("sigma", seq_len(npsc - 1)))

  ldpar <- lmomco::parglo(lmomco::lmoms(xdat[, 1]))$para[1:2]

  if (is.null(mul)) {
    mumat <- matrix(1, nrow = nrow(xdat), ncol = 1)
    if (is.null(muinit)) muinit <- ldpar[1]
  } else {
    z$trans <- TRUE
    mumat <- cbind(1, as.matrix(ydat[, mul, drop = FALSE]))
    if (is.null(muinit)) muinit <- c(ldpar[1], rep(0, length(mul)))
  }

  if (is.null(sigl)) {
    sigmat <- matrix(1, nrow = nrow(xdat), ncol = 1)
    if (is.null(siginit)) siginit <- ldpar[2]
  } else {
    z$trans <- TRUE
    sigmat <- cbind(1, as.matrix(ydat[, sigl, drop = FALSE]))
    if (is.null(siginit)) siginit <- c(ldpar[2], rep(0, length(sigl)))
  }

  z$model <- list(mul = mul, sigl = sigl)
  z$link <- c(
    mulink = deparse(substitute(mulink)),
    siglink = deparse(substitute(siglink))
  )

  zr <- drop(xdat[, r, drop = FALSE])

  init <- c(muinit, siginit)
  names(init) <- c(mu_names, sigma_names)

  init_list <- vector("list", num_inits)
  init_list[[1]] <- init

  if (num_inits >= 2) {
    for (i in 2:num_inits) {
      init_list[[i]] <- init + c(
        stats::rnorm(npmu, mean = 0, sd = 0.3),
        abs(stats::rnorm(npsc, mean = 0, sd = 0.2))
      )
    }
  }

  ## ----------------------------------
  ## Likelihood for r = 1
  ## ----------------------------------
  ld.lik <- function(a) {
    out <- tryCatch({
      mu <- drop(mulink(mumat %*% a[1:npmu]))
      sc <- drop(siglink(sigmat %*% a[seq(npmu + 1, length.out = npsc)]))

      if (any(!is.finite(mu)) || any(!is.finite(sc))) return(BIG)
      if (any(sc <= 0, na.rm = TRUE)) return(BIG)

      x1 <- drop(xdat[, 1, drop = FALSE])

      y <- exp(-(x1 - mu) / sc)
      f <- 1 + exp(-(x1 - mu) / sc)

      if (any(!is.finite(y), na.rm = TRUE) ||
          any(!is.finite(f), na.rm = TRUE) ||
          any(y <= 0, na.rm = TRUE) ||
          any(f <= 0, na.rm = TRUE)) {
        return(BIG)
      }

      nll <- sum(log(sc) - log(y) + 2 * log(f), na.rm = TRUE)

      if (!is.finite(nll)) BIG else nll
    }, error = function(e) BIG)

    out
  }

  ## ----------------------------------
  ## Likelihood for r > 1
  ## ----------------------------------
  rld.lik <- function(a) {
    out <- tryCatch({
      mu <- drop(mulink(mumat %*% a[1:npmu]))
      sc <- drop(siglink(sigmat %*% a[seq(npmu + 1, length.out = npsc)]))

      if (any(!is.finite(mu)) || any(!is.finite(sc))) return(BIG)
      if (any(sc <= 0, na.rm = TRUE)) return(BIG)

      ri <- r - seq_len(r)
      cr <- 1 + ri

      if (any(cr <= 0, na.rm = TRUE)) return(BIG)

      xmat <- xdat
      y <- exp(-(xmat - mu) / sc)
      f <- 1 + exp(-(zr - mu) / sc)

      if (any(!is.finite(y), na.rm = TRUE) ||
          any(!is.finite(f), na.rm = TRUE) ||
          any(y <= 0, na.rm = TRUE) ||
          any(f <= 0, na.rm = TRUE)) {
        return(BIG)
      }

      yy <- log(sc) - log(y) - matrix(log(cr), nrow = nrow(xmat), ncol = ncol(xmat), byrow = TRUE)
      yy <- rowSums(yy, na.rm = TRUE)

      nll <- sum((r + 1) * log(f) + yy, na.rm = TRUE)

      if (!is.finite(nll)) BIG else nll
    }, error = function(e) BIG)

    out
  }

  likfun <- if (r == 1) ld.lik else rld.lik

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

  z$vals <- cbind(mu = mu, sigma = sc)

  if (show) {
    if (z$trans) {
      print(z[c("model", "link")])
    }
    if (!is.null(z$conv) && identical(z$conv, 0L)) {
      print(z[c("conv", "nllh", "mle", "se")])
    }
  }

  class(z) <- "rld.fit"
  invisible(z)
}
