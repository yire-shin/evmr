#' Fit the Gumbel Distribution to r-Largest Order Statistics
#'
#' Fits the Gumbel distribution to \eqn{r}-largest order statistics using
#' maximum likelihood estimation. Stationary and non-stationary models are
#' supported through generalized linear modelling of the location and scale
#' parameters.
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
#' @param ydat A matrix or data frame of covariates for non-stationary modelling
#'   of the parameters, or \code{NULL} for a stationary model. The number of rows
#'   must match the number of rows of \code{xdat}.
#' @param mul,sigl Integer vectors indicating which columns of \code{ydat} are
#'   used as covariates for the location and scale parameters, respectively.
#'   Use \code{NULL} for stationary parameters.
#' @param mulink,siglink Inverse link functions for the location and scale
#'   parameters, respectively.
#' @param num_inits The number of initial parameter sets used in the optimization.
#' @param muinit,siginit Numeric vectors giving initial values for the location
#'   and scale parameters. If \code{NULL}, default initial values based on
#'   L-moments are used.
#' @param show Logical. If \code{TRUE}, details of the fitted model are printed.
#' @param method Optimization method passed to \code{\link{optim}} for stationary fits.
#' @param maxit Maximum number of iterations for \code{\link{optim}}.
#' @param ... Additional control arguments passed to the optimizer.
#'
#' @return A list with components including:
#' \item{trans}{Logical; \code{TRUE} if a non-stationary model is fitted.}
#' \item{model}{A list containing \code{mul} and \code{sigl}.}
#' \item{link}{A character string describing the inverse link functions.}
#' \item{conv}{The convergence code returned by the optimizer. A value of 0
#' indicates successful convergence for \code{optim}.}
#' \item{nllh}{The negative log-likelihood evaluated at the fitted parameters.}
#' \item{data}{The data used in the fit.}
#' \item{mle}{The maximum likelihood estimates.}
#' \item{cov}{The estimated covariance matrix.}
#' \item{se}{The estimated standard errors.}
#' \item{vals}{A matrix containing fitted values of the location and scale
#' parameters at each observation.}
#' \item{r}{The number of order statistics used in the fitted model.}
#'
#' @references
#'
#' Coles, S. (2001).
#' An Introduction to Statistical Modeling of Extreme Values.
#' Springer.
#'
#' Shin, Y., & Park, J.-S. (2025).
#' Generalized Gumbel model for r-largest order statistics with application
#' to peak streamflow.
#' \emph{Scientific Reports}.
#' \doi{10.1038/s41598-024-83273-y}

#'
#' @seealso \code{\link{optim}}
#' @export
#'
#' @examples
#' x <- rgdr(n = 50, r = 2, loc = 10, scale = 2)
#' fit <- rgd.fit(x$rmat)
#' fit$r
#' fit$mle
rgd.fit <- function(xdat, r = NULL, ydat = NULL, mul = NULL, sigl = NULL,
                    mulink = identity, siglink = identity, num_inits = 100,
                    muinit = NULL, siginit = NULL, show = TRUE,
                    method = "Nelder-Mead", maxit = 10000, ...) {

  z <- list()

  if (is.null(r)) {
    if (is.vector(xdat)) {
      xdat <- matrix(xdat, ncol = 1)
      r <- 1
    } else {
      xdat <- as.matrix(xdat)
      r <- ncol(xdat)
    }
  } else {
    if (is.vector(xdat)) {
      if (r != 1) {
        stop("If 'xdat' is a vector, 'r' must be 1.")
      }
      xdat <- matrix(xdat, ncol = 1)
    } else {
      xdat <- as.matrix(xdat)
      xdat <- xdat[, 1:r, drop = FALSE]
    }
  }

  xdat <- as.data.frame(xdat)

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

  gdpar <- lmomco::pargum(lmomco::lmoms(xdat[, 1]))$para

  if (is.null(mul)) {
    mumat <- as.matrix(rep(1, nrow(xdat)))
    if (is.null(muinit)) muinit <- gdpar[1]
  } else {
    z$trans <- TRUE
    mumat <- cbind(rep(1, nrow(xdat)), ydat[, mul, drop = FALSE])
    if (is.null(muinit)) muinit <- c(gdpar[1], rep(0, length(mul)))
  }

  if (is.null(sigl)) {
    sigmat <- as.matrix(rep(1, nrow(xdat)))
    if (is.null(siginit)) siginit <- gdpar[2]
  } else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, nrow(xdat)), ydat[, sigl, drop = FALSE])
    if (is.null(siginit)) siginit <- c(gdpar[2], rep(0, length(sigl)))
  }

  z$model <- list(mul = mul, sigl = sigl)
  z$link <- deparse(substitute(c(mulink, siglink)))

  zr <- as.matrix(xdat[, r, drop = FALSE])

  init <- c(muinit, siginit)
  names(init) <- c(mu_names, sigma_names)

  init_list <- list(init)
  for (i in 2:num_inits) {
    new_init <- init + c(
      stats::rnorm(npmu, mean = 0, sd = 1),
      abs(stats::rnorm(npsc, mean = 0, sd = 1))
    )
    init_list[[i]] <- new_init
  }

  rgd.lik <- function(a) {

    mu <- drop(mulink(mumat %*% a[1:npmu]))
    sc <- drop(siglink(sigmat %*% a[seq(npmu + 1, length.out = npsc)]))

    if (any(sc <= 0, na.rm = TRUE)) {
      return(1e6)
    }

    xmat <- as.matrix(xdat)
    zr1  <- drop(as.matrix(zr))

    y <- sweep(xmat, 1, mu, "-")
    y <- sweep(y,    1, sc, "/")
    y <- sweep(y,    1, log(sc), "+")
    y <- rowSums(y, na.rm = TRUE)

    sum(exp(-(zr1 - mu) / sc) + y, na.rm = TRUE)
  }

  optim_results <- lapply(init_list, function(init) {
    if (!z$trans) {
      stats::optim(
        init, rgd.lik, hessian = TRUE, method = method,
        control = list(maxit = maxit, trace = 0)
      )
    } else {
      suppressWarnings(Rsolnp::solnp(
        init, rgd.lik, control = list(trace = 0)
      ))
    }
  })

  optim_value <- data.frame(
    num = seq_along(optim_results),
    nllh = sapply(optim_results, function(res) {
      if (z$trans) min(res$values) else res$value
    }),
    grad = sapply(optim_results, function(res) {
      par_now <- if (z$trans) res$pars else res$par
      sum(abs(numDeriv::grad(rgd.lik, par_now)))
    })
  )

  optim_value <- optim_value[optim_value$nllh != 1e6, , drop = FALSE]
  optim_value <- optim_value[order(optim_value$grad, optim_value$nllh), , drop = FALSE]

  best_result <- optim_results[[optim_value$num[1]]]
  best_par <- if (z$trans) best_result$pars else best_result$par

  mu <- drop(mumat %*% best_par[1:npmu])
  sc <- drop(sigmat %*% best_par[seq(npmu + 1, length.out = npsc)])

  z$r <- r
  z$conv <- if (z$trans) best_result$convergence else best_result$convergence
  z$nllh <- if (z$trans) min(best_result$values) else best_result$value
  z$data <- xdat
  z$mle <- best_par

  if (!z$trans && !is.null(best_result$hessian)) {
    z$cov <- solve(best_result$hessian)
    z$se <- sqrt(diag(z$cov))
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

  class(z) <- "rgd.fit"
  invisible(z)
}
