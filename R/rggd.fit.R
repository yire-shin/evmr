#' Fit the Generalized Gumbel Distribution to r-Largest Order Statistics
#'
#' Fits the generalized Gumbel distribution to \eqn{r}-largest order statistics
#' using maximum likelihood estimation. Stationary and non-stationary models
#' are supported through generalized linear modelling of the location, scale,
#' and shape parameters.
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
#' @param mul,sigl,hl Integer vectors indicating which columns of \code{ydat}
#'   are used as covariates for the location, scale, and shape parameters,
#'   respectively. Use \code{NULL} for stationary parameters.
#' @param mulink,siglink,hlink Inverse link functions for the location, scale,
#'   and shape parameters, respectively.
#' @param num_inits The number of initial parameter sets used in the optimization.
#' @param muinit,siginit,hinit Numeric vectors giving initial values for the
#'   location, scale, and shape parameters. If \code{NULL}, default initial
#'   values based on L-moments are used.
#' @param show Logical. If \code{TRUE}, details of the fitted model are printed.
#' @param method Optimization method passed to \code{\link{optim}} for stationary fits.
#' @param maxit Maximum number of iterations for \code{\link{optim}}.
#' @param ... Additional control arguments passed to the optimizer.
#'
#' @return A list with components including:
#' \item{trans}{Logical; \code{TRUE} if a non-stationary model is fitted.}
#' \item{model}{A list containing \code{mul}, \code{sigl}, and \code{hl}.}
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
#'#' @references
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
#' x <- rggdr(n = 50, r = 2, loc = 10, scale = 2, shape = 0.1)
#' fit <- rggd.fit(x$rmat)
#' fit$r
#' fit$mle
rggd.fit <- function(xdat, r = NULL, ydat = NULL, mul = NULL, sigl = NULL, hl = NULL,
                     mulink = identity, siglink = identity, hlink = identity,
                     num_inits = 100, muinit = NULL, siginit = NULL, hinit = NULL,
                     show = TRUE, method = "Nelder-Mead", maxit = 10000, ...) {

  z <- list()
  tol <- .Machine$double.eps^0.5

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
  nph  <- length(hl) + 1

  z$trans <- FALSE

  mu_names <- if (is.null(mul)) "mu" else c("mu0", paste0("mu", seq_len(npmu - 1)))
  sigma_names <- if (is.null(sigl)) "sigma" else c("sigma0", paste0("sigma", seq_len(npsc - 1)))
  h_names <- if (is.null(hl)) "h" else c("h0", paste0("h", seq_len(nph - 1)))

  ggdpar <- lmomco::pargev(lmomco::lmoms(xdat[, 1]))$para

  if (is.null(mul)) {
    mumat <- matrix(1, nrow = nrow(xdat), ncol = 1)
    if (is.null(muinit)) muinit <- ggdpar[1]
  } else {
    z$trans <- TRUE
    mumat <- cbind(1, as.matrix(ydat[, mul, drop = FALSE]))
    if (is.null(muinit)) muinit <- c(ggdpar[1], rep(0, length(mul)))
  }

  if (is.null(sigl)) {
    sigmat <- matrix(1, nrow = nrow(xdat), ncol = 1)
    if (is.null(siginit)) siginit <- ggdpar[2]
  } else {
    z$trans <- TRUE
    sigmat <- cbind(1, as.matrix(ydat[, sigl, drop = FALSE]))
    if (is.null(siginit)) siginit <- c(ggdpar[2], rep(0, length(sigl)))
  }

  if (is.null(hl)) {
    hmat <- matrix(1, nrow = nrow(xdat), ncol = 1)
    if (is.null(hinit)) hinit <- ggdpar[3]
  } else {
    z$trans <- TRUE
    hmat <- cbind(1, as.matrix(ydat[, hl, drop = FALSE]))
    if (is.null(hinit)) hinit <- c(ggdpar[3], rep(0, length(hl)))
  }

  z$model <- list(mul = mul, sigl = sigl, hl = hl)
  z$link <- c(
    mulink = deparse(substitute(mulink)),
    siglink = deparse(substitute(siglink)),
    hlink = deparse(substitute(hlink))
  )

  zr <- drop(xdat[, r, drop = FALSE])

  init <- c(muinit, siginit, hinit)
  names(init) <- c(mu_names, sigma_names, h_names)

  init_list <- vector("list", num_inits)
  init_list[[1]] <- init
  if (num_inits >= 2) {
    for (i in 2:num_inits) {
      init_list[[i]] <- init + c(
        stats::rnorm(npmu, mean = 0, sd = 1),
        abs(stats::rnorm(npsc, mean = 0, sd = 1)),
        stats::rnorm(nph, mean = 0, sd = 0.5)
      )
    }
  }

  ## ----------------------------------
  ## Likelihood for r = 1
  ## ----------------------------------
  ggd.lik <- function(a) {
    mu <- drop(mulink(mumat %*% a[1:npmu]))
    sc <- drop(siglink(sigmat %*% a[seq(npmu + 1, length.out = npsc)]))
    h  <- drop(hlink(hmat %*% a[seq(npmu + npsc + 1, length.out = nph)]))

    if (any(!is.finite(mu)) || any(!is.finite(sc)) || any(!is.finite(h))) {
      return(1e6)
    }
    if (any(sc <= 0, na.rm = TRUE)) {
      return(1e6)
    }
    if (any(abs(h) < tol, na.rm = TRUE)) {
      return(1e6)
    }

    x1 <- drop(xdat[, 1, drop = FALSE])

    y <- exp(-(x1 - mu) / sc)
    f <- 1 - h * y

    if (any(y <= 0, na.rm = TRUE) || any(f <= 0, na.rm = TRUE)) {
      return(1e6)
    }

    sum(log(sc) - log(y) - ((1 - h) / h) * log(f), na.rm = TRUE)
  }

  ## ----------------------------------
  ## Likelihood for r > 1
  ## ----------------------------------
  rggd.lik <- function(a) {
    mu <- drop(mulink(mumat %*% a[1:npmu]))
    sc <- drop(siglink(sigmat %*% a[seq(npmu + 1, length.out = npsc)]))
    h  <- drop(hlink(hmat %*% a[seq(npmu + npsc + 1, length.out = nph)]))

    if (any(!is.finite(mu)) || any(!is.finite(sc)) || any(!is.finite(h))) {
      return(1e6)
    }
    if (any(sc <= 0, na.rm = TRUE)) {
      return(1e6)
    }
    if (any(abs(h) < tol, na.rm = TRUE)) {
      return(1e6)
    }

    if (r >= 2 && min(h, na.rm = TRUE) > (1 / (r - 1))) {
      return(1e6)
    }

    ri <- r - seq_len(r)
    cr <- 1 - ri * h[1]

    if (any(cr <= 0, na.rm = TRUE)) {
      return(1e6)
    }

    xmat <- xdat
    y <- sweep(xmat, 1, mu, "-")
    y <- sweep(y, 1, sc, "/")
    y <- exp(-y)

    f <- 1 - h[1] * exp(-(zr - mu) / sc)

    if (any(y <= 0, na.rm = TRUE) || any(f <= 0, na.rm = TRUE)) {
      return(1e6)
    }

    yy <- log(sc) - log(y) - matrix(log(cr), nrow = nrow(xmat), ncol = ncol(xmat), byrow = TRUE)
    yy <- rowSums(yy, na.rm = TRUE)

    sum(((r * h[1] - 1) / h[1]) * log(f) + yy, na.rm = TRUE)
  }

  likfun <- if (r == 1) ggd.lik else rggd.lik

  ## ----------------------------------
  ## Optimization
  ## ----------------------------------
  optim_results <- lapply(init_list, function(init_now) {
    if (!z$trans) {
      stats::optim(
        init_now, likfun, hessian = TRUE, method = method,
        control = list(maxit = maxit, trace = 0)
      )
    } else {
      suppressWarnings(
        Rsolnp::solnp(
          pars = init_now,
          fun = likfun,
          control = list(trace = 0)
        )
      )
    }
  })

  ## ----------------------------------
  ## Score results
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

  optim_value <- optim_value[is.finite(optim_value$nllh) & optim_value$nllh != 1e6, , drop = FALSE]

  if (nrow(optim_value) == 0) {
    stop("All optimization attempts failed. Try different initial values or reduce 'r'.")
  }

  optim_value <- optim_value[order(optim_value$grad, optim_value$nllh), , drop = FALSE]

  best_result <- optim_results[[optim_value$num[1]]]
  best_par <- if (z$trans) best_result$pars else best_result$par

  mu <- drop(mumat %*% best_par[1:npmu])
  sc <- drop(sigmat %*% best_par[seq(npmu + 1, length.out = npsc)])
  h  <- drop(hmat %*% best_par[seq(npmu + npsc + 1, length.out = nph)])

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

  z$vals <- cbind(mu = mu, sigma = sc, h = h)

  if (show) {
    if (z$trans) print(z[c("model", "link")])
    if (!is.null(z$conv) && identical(z$conv, 0L)) {
      print(z[c("conv", "nllh", "mle", "se")])
    }
  }

  class(z) <- "rggd.fit"
  invisible(z)
}
