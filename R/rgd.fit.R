#' @name rgd.fit
#' @aliases rgd.fit
#' @title Gumbel Distribution for $r$ Largest Order Statistics
#' @description
#' Maximum-likelihood fitting for the order statistic model,including generalized linear modelling of each parameter.
#'
#' @param xdat A numeric matrix of data to be fitted. Each row should be a vector of decreasing order, containing the largest order statistics for each year (or time period). The first column therefore contains annual (or period) maxima. Only the first \code{r} columns are used for the fitted model. By default, all columns are used.If one year (or time period) contains fewer order statistics than another, missing values can be appended to the end of the corresponding row.#'
#' @param r The largest \code{r} order statistics are used for the fitted model.
#' @param ydat A matrix of covariates for generalized linear modelling of the parameters (or \code{NULL} (the default) for stationary fitting). The number of rows should be the same as the number of rows of \code{xdat}.
#' @param mul,sigl Numeric vectors of integers, giving the columns of \code{ydat} that contain covariates for generalized linear modelling of the location, scale and shape parameters repectively (or \code{NULL} (the default) if the corresponding parameter is stationary).
#' @param mulink,siglink Inverse link functions for generalized linear modelling of the location, scale and shape parameters repectively.
#' @param num_inits Specifies the number of initial parameter sets to be generated for the optimization process.
#' @param muinit,siginit numeric of length equal to total number of parameters used to model the location, scale or shape parameter(s), resp.  See Details section for default (NULL) initial values.
#' @param show Logical; if \code{TRUE} (the default), print details of the fit.
#' @param method The optimization method (see \code{\link{optim}} for details).
#' @param maxit The maximum number of iterations.
#' @param \dots Other control parameters for the optimization. These are passed to components of the \code{control} argument of \code{optim}.
#'
#' @return {
#'  A list containing the following components. A subset of these components are printed after the fit. If \code{show} is \code{TRUE}, then assuming that successful convergence is indicated, the components \code{nllh}, \code{mle} and \code{se} are always printed.
#'
#'  \item{trans}{An logical indicator for a non-stationary fit.}\item{model}{A list with components \code{mul}, \code{sigl} and \code{shl}.} \item{link}{A character vector giving inverse link functions.} \item{conv}{The convergence code, taken from the list returned by \code{\link{optim}}. A zero indicates successful convergence.} \item{nllh}{The negative logarithm of the likelihood evaluated at the maximum likelihood estimates.}
#'  \item{model}{A list with components \code{mul}, \code{sigl}, \code{shl} and \code{hl}.}
#'  \item{link}{A character vector giving inverse link functions.}
#'  \item{conv}{The convergence code, taken from the list returned by \code{\link{optim}}. A zero indicates successful convergence.}
#'  \item{nllh}{The negative logarithm of the likelihood evaluated at the maximum likelihood estimates.}
#'  \item{data}{The data that has been fitted. For non-stationary models, the data is standardized.}
#'  \item{mle}{A vector containing the maximum likelihood estimates.}
#'  \item{cov}{The covariance matrix.}
#'  \item{se}{A vector containing the standard errors.}
#'  \item{vals}{A matrix with three columns containing the maximum likelihood estimates of the location, scale and shape parameters at each data point.}
#'  \item{r}{The number of order statistics used.}
#' }
#' @seealso \code{\link{optim}}
#' @export
#'
#' @examples
#' \dontrun{
#' data(bangkok)
#' rgd.fit(bangkok)
#' }
rgd.fit <- function(xdat, r = NULL, ydat = NULL, mul = NULL, sigl = NULL,
                    mulink = identity, siglink = identity, num_inits = 100,
                    muinit = NULL, siginit = NULL, show = TRUE, method = "Nelder-Mead", maxit = 10000, ...){

  options(digits=8)
  z <- list()

  # If r is NULL, set r to the number of columns in xdat
  if (is.null(r)) {
    if (is.vector(xdat)) {
      # If xdat is a vector, convert it to a matrix with one column and set r = 1
      xdat <- matrix(xdat, ncol = 1)
      r <- 1
    } else {
      # If xdat is a matrix, keep it as is and set r to the number of columns
      r <- dim(xdat)[2]
    }
  } else {
    # 2. When r is specified
    if (r == 1) {
      # If r = 1, ensure xdat remains a matrix with one column
      if (is.vector(xdat)) {
        xdat <- matrix(xdat,ncol=1)
      } else {
        # If xdat is already a matrix, ensure it remains a matrix with one column
        xdat <- as.matrix(xdat[, 1:r, drop = FALSE])
      }
    } else {
      # If r > 1, subset xdat to the first r columns, maintaining it as a matrix
      #xdat <- matrix(as.matrix(xdat[,1:r],ncol=r),ncol=r)
      xdat <- as.matrix(xdat[,1:r],ncol=r)
    }
  }

  # Determine the number of parameters for each component (mu, sigma)
  npmu <- length(mul) + 1
  npsc <- length(sigl) + 1

  z$trans <- FALSE

  # Generate parameter names based on the length of each list
  mu_names <- if (is.null(mul)) "mu" else c("mu0", paste0("mu", seq_len(npmu - 1)))
  sigma_names <- if (is.null(sigl)) "sigma" else c("sigma0", paste0("sigma", seq_len(npsc - 1)))

  # Set initial values based on L-moments of the data and user-specified predictors
  gdpar <- lmomco::pargum(lmomco::lmoms(xdat[, 1]))$para

  # Generate multiple sets of initial values, each with random perturbations
  init_list <- list(gdpar)


  # Set up the mu matrix and initial values if each component (mu, sigma) is provided

  if(is.null(mul)) {
    mumat <- as.matrix(rep(1, dim(xdat)[1]))
    if( is.null( muinit)) muinit <- gdpar[1]
  } else {
    z$trans <- TRUE
    mumat <- cbind(rep(1, dim(xdat)[1]), ydat[, mul])
    if( is.null( muinit)) muinit <- c(gdpar[1], rep(0, length(mul)))
  }
  if(is.null(sigl)) {
    sigmat <- as.matrix(rep(1, dim(xdat)[1]))
    if( is.null( siginit)) siginit <- gdpar[2]
  } else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, dim(xdat)[1]), ydat[, sigl])
    if( is.null( siginit)) siginit <- c(gdpar[2], rep(0, length(sigl)))
  }

  z$model <- list(mul, sigl)
  z$link <- deparse(substitute(c(mulink, siglink)))

  z1 <- as.matrix(xdat[, 1],ncol=1)
  zr <- as.matrix(xdat[, r],ncol=1)

  init <- c(muinit, siginit)
  names(init) <- c(mu_names, sigma_names)

  # Generate multiple sets of initial values, each with random perturbations
  init_list <- list(init)

  for (i in 2:num_inits) {
    # Apply random noise to each component of init to create diversity in starting points
    new_init <- init + c(
      stats::rnorm(npmu, mean = 0, sd = 1),  # Random noise for mu parameters
      abs(stats::rnorm(npsc, mean = 0, sd = 1))  # Random noise for sigma parameters
    )
    init_list[[i]] <- new_init
  }

  # Define the log-likelihood function for the case when r = 1
  rgd.lik <- function(a) {

    mu <- mulink(mumat %*% (a[1:npmu]))
    sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))

    if (any(sc <= 0))
      return(10^6)
    y <- (xdat - mu)/sc

    y <- y+log(sc)
    y <- rowSums(y, na.rm = TRUE)

    sum( exp(- (zr-mu)/sc) + y, na.rm=T)

  }


  # Apply optimization on each set of initial values and retain results
  optim_results <- lapply(init_list, function(init) {
    if (r == 1) {
      if(z$trans==F){stats::optim(init, rgd.lik, hessian=TRUE, method=method, control=list(maxit=maxit, trace = 0))
      }else{suppressWarnings(Rsolnp::solnp(init, rgd.lik, control = list(trace = 0)))}
    } else {
      if(z$trans==F){stats::optim(init, rgd.lik, hessian=TRUE, method=method, control=list(maxit=maxit, trace = 0))
      }else{suppressWarnings(Rsolnp::solnp(init, rgd.lik, control = list(trace = 0)))}
    }
  })

  # Collect optimization results and filter out invalid results
  optim_value <- data.frame(
    num = 1:length(optim_results),
    nllh = sapply(optim_results, function(res) {
                  if(z$trans) min(res$values) else res$value}),
    grad = sapply(optim_results, function(res) {
      sum(abs(if (r == 1) numDeriv::grad(rgd.lik, res$par) else numDeriv::grad(rgd.lik, res$par)))
    })
  )


  optim_value <- optim_value[optim_value$nllh != 10^6, ]
  optim_value <- optim_value[order(optim_value$grad, optim_value$nllh), ]
  best_result <- optim_results[[optim_value$num[1]]]

  # Extract and store the best-fit parameter values
  mu <- drop(mumat %*% (best_result$par[1:npmu]))
  sc <- drop(sigmat %*% (best_result$par[seq(npmu + 1, length = npsc)]))

  # Store results in the output list with dynamic parameter names
  z$r    <- r
  z$conv <- best_result$convergence
  z$nllh <- best_result$value
  z$data <- xdat
  z$mle  <- best_result$par
  z$cov  <- solve(best_result$hessian)
  z$se   <- sqrt(diag(z$cov))
  z$vals <- cbind(mu, sc)

  if(show) {
    if(z$trans)
      print(z[c(2, 3)])
    #else print(z[4])
    if(!z$conv)
      print(z[c(4, 6, 8, 10)])
  }

  class(z) <- "rgd.fit"
  invisible(z)
}
