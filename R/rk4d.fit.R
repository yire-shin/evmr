#' @name rk4d.fit
#' @aliases rk4d.fit
#' @title Four Parameter Kappa Distribution for $r$ Largest Order Statistics
#' @description
#' Maximum-likelihood fitting for the order statistic model,including generalized linear modelling of each parameter.
#'
#' @param xdat A numeric matrix of data to be fitted. Each row should be a vector of decreasing order, containing the largest order statistics for each year (or time period). The first column therefore contains annual (or period) maxima. Only the first \code{r} columns are used for the fitted model. By default, all columns are used.If one year (or time period) contains fewer order statistics than another, missing values can be appended to the end of the corresponding row.#'
#' @param r The largest \code{r} order statistics are used for the fitted model.
#' @param penk Maximum penalized likelihood estimation. Coles and Dixon's Penalty function, Martins and Stedinger's Penalty function
#' @param penh Maximum penalized likelihood estimation.  Coles and Dixon's Penalty function, Martins and Stedinger's Penalty function
#' @param ydat A matrix of covariates for generalized linear modelling of the parameters (or \code{NULL} (the default) for stationary fitting). The number of rows should be the same as the number of rows of \code{xdat}.
#' @param mul,sigl,shl,hl Numeric vectors of integers, giving the columns of \code{ydat} that contain covariates for generalized linear modelling of the location, scale and shape parameters repectively (or \code{NULL} (the default) if the corresponding parameter is stationary).
#' @param mulink,siglink,shlink,hlink Inverse link functions for generalized linear modelling of the location, scale and shape parameters repectively.
#' @param num_inits Specifies the number of initial parameter sets to be generated for the optimization process.
#' @param muinit,siginit,shinit,hinit numeric of length equal to total number of parameters used to model the location, scale or shape parameter(s), resp.  See Details section for default (NULL) initial values.
#' @param show Logical; if \code{TRUE} (the default), print details of the fit.
#' @param method The optimization method (see \code{\link{optim}} for details).
#' @param maxit The maximum number of iterations.
#' @param \dots Other control parameters for the optimization. These are passed to components of the \code{control} argument of \code{optim}.
#'
#' @return {
#'  A list containing the following components. A subset of these components are printed after the fit. If \code{show} is \code{TRUE}, then assuming that successful convergence is indicated, the components \code{nllh}, \code{mle} and \code{se} are always printed.
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
#' rk4d.fit(bangkok)
#' }
rk4d.fit <- function(xdat, r = NULL, penk = NULL, penh = NULL, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL, hl = NULL,
                     mulink = identity, siglink = identity, shlink = identity, hlink = identity, num_inits = 100,
                     muinit = NULL, siginit = NULL, shinit = NULL, hinit = NULL, show = TRUE, method = "Nelder-Mead", maxit = 10000, ...){

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


  # Determine the number of parameters for each component (mu, sigma, xi, h)
  npmu <- length(mul) + 1
  npsc <- length(sigl) + 1
  npsh <- length(shl) + 1
  nph  <- length(hl) + 1

  z$trans <- FALSE

  # Generate parameter names based on the length of each list
  mu_names <- if (is.null(mul)) "mu" else c("mu0", paste0("mu", seq_len(npmu - 1)))
  sigma_names <- if (is.null(sigl)) "sigma" else c("sigma0", paste0("sigma", seq_len(npsc - 1)))
  xi_names <- if (is.null(shl)) "xi" else c("xi0", paste0("xi", seq_len(npsh - 1)))
  h_names <- if (is.null(hl)) "h" else c("h0", paste0("h", seq_len(nph - 1)))

  # Set initial values based on L-moments of the data and user-specified predictors
  kappar <- lmomco::parkap(lmomco::lmoms(xdat[, 1]))$para

  # Generate multiple sets of initial values, each with random perturbations
  init_list <- list(kappar)


  # Set up the mu matrix and initial values if each component (mu, sigma, xi, h) is provided

  if(is.null(mul)) {
    mumat <- as.matrix(rep(1, dim(xdat)[1]))
    if( is.null( muinit)) muinit <- kappar[1]
  } else {
    z$trans <- TRUE
    mumat <- cbind(rep(1, dim(xdat)[1]), ydat[, mul])
    if( is.null( muinit)) muinit <- c(kappar[1], rep(0, length(mul)))
  }
  if(is.null(sigl)) {
    sigmat <- as.matrix(rep(1, dim(xdat)[1]))
    if( is.null( siginit)) siginit <- kappar[2]
  } else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, dim(xdat)[1]), ydat[, sigl])
    if( is.null( siginit)) siginit <- c(kappar[2], rep(0, length(sigl)))
  }
  if(is.null(shl)) {
    shmat <- as.matrix(rep(1, dim(xdat)[1]))
    if( is.null( shinit)) shinit <- kappar[3]
  }  else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, dim(xdat)[1]), ydat[, shl])
    if( is.null( shinit)) shinit <- c(kappar[3], rep(0, length(shl)))
  }
  if(is.null(hl)) {
    hmat <- as.matrix(rep(1, dim(xdat)[1]))
    if( is.null( hinit)) hinit <- kappar[4]
  }  else {
    z$trans <- TRUE
    hmat <- cbind(rep(1, dim(xdat)[1]), ydat[, hl])
    if( is.null( hinit)) hinit <- c(kappar[4], rep(0, length(hl)))
  }

  z$model <- list(mul, sigl, shl, hl)
  z$link <- deparse(substitute(c(mulink, siglink, shlink, hlink)))

  z1 <- as.matrix(xdat[, 1],ncol=1)
  zr <- as.matrix(xdat[, r],ncol=1)

  init <- c(muinit, siginit, shinit, hinit)
  names(init) <- c(mu_names, sigma_names, xi_names, h_names)

  # Generate multiple sets of initial values, each with random perturbations
  init_list <- list(init)

  for (i in 2:num_inits) {
    # Apply random noise to each component of init to create diversity in starting points
    new_init <- init + c(
      stats::rnorm(npmu, mean = 0, sd = 1),  # Random noise for mu parameters
      abs(stats::rnorm(npsc, mean = 0, sd = 1)),  # Random noise for sigma parameters
      stats::rnorm(npsh, mean = 0, sd = 0.1), # Random noise for xi parameters
      stats::rnorm(nph, mean = 0, sd = 0.1)   # Random noise for h parameters
    )
    init_list[[i]] <- new_init
  }

  # Define the log-likelihood function for the case when r = 1
  k4d.lik <- function(a) {

    mu <- drop(mumat %*% (a[1:npmu]))
    sc <- drop(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    xi <- drop(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    h  <- drop(hmat %*% (a[seq(npmu + npsc + nph + 1, length = nph)]))

    y <- (xdat - mu) / sc
    y <- 1 - xi * y
    f <- 1 - h * y^(1/xi)

    if (any(y <= 0) || any(sc <= 0) || any(f <= 0) || any(f^(1/h) > 1))
      return(10^6)

    # maximum penalized likelihood estimation

    if(is.null(penk)==F){

      if(penk=="CD"){

        if(xi[1]>=0){ p_k <- 1
        }else if( xi[1] > -1 && xi[1] <0 ){ p_k <- exp(-((1/(1+xi[1])) -1))
        }else if( xi[1] <= -1){ p_k <- 0}

      }else if(penk=="MS"){

        if( xi[1] >= 0.5 || xi[1] <= -0.5) return(10^6)
        p_k <-((0.5+xi[1])^(6-1)) * ((0.5-xi[1])^(9-1)) / beta(6,9) # MS_for k

      }}

    if(is.null(penh)==F){
      if(penh=="MS"){

        if( h[1] <= -0.5 || h[1] >= 0.5) return(10^6)
        Bef <-function(x) { ((0.5+x)^(6-1)) * ((0.5-x)^(9-1)) }
        Be  <-stats::integrate(Bef, lower=-0.5, upper= 0.5)[1]$value
        p_h <-((0.5+h[1])^(6-1)) * ((0.5-h[1])^(9-1)) / Be

      }else if(penh=="MSa"){

        if( h[1] >= 1.2 || h[1] <= -1.2) return(10^6)
        Bef <-function(x) { ((1.2+x)^(6-1)) * ((1.2-x)^(9-1)) }
        Be  <-stats::integrate(Bef, lower=-1.2, upper= 1.2 )[1]$value
        p_h <-((1.2+h[1])^(6-1)) * ((1.2-h[1])^(9-1)) / Be

      }}

    if(is.null(penk)==T & is.null(penh)==T){ penalty<-0
    }else if(is.null(penk)==F & is.null(penh)==F){ penalty <- r*log(p_k * p_h)
    }else if(is.null(penk)==F){ penalty <- r*log(p_k)
    }else if(is.null(penh)==F){ penalty <- r*log(p_h)}

    sum(log(sc)) + sum(log(1 - h * y^(1/xi)) * ((h - 1) / h)) + sum(log(y) * (1 - 1/xi)) - penalty

  }

  # Define the log-likelihood function for the case when r > 1
  rk4d.lik <- function(a) {

    mu <- drop(mumat %*% (a[1:npmu]))
    sc <- drop(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    xi <- drop(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    h  <- drop(hmat %*% (a[seq(npmu + npsc + nph + 1, length = nph)]))

    ri  <- (r - seq(1:(r)))
    cr  <- (1 - ri * h[1])

    if (any(sc <= 0) | any(cr < 0)) return(10^6)

    y <- 1 - xi * (xdat - mu) / sc

    if (any(y<0)) return(10^6)

    f <- 1 - h * (1 - xi * (zr - mu) / sc)^(1/xi)

    if (max(h) > min(1 / (1 - xi * (zr - mu) / sc)^(1/xi), na.rm=TRUE)) return(10^6)
    if (r >= 2 && min(h) > (1 / (r - 1))) return(10^6)
    if (min(y, na.rm=TRUE) <= 0 || min(f, na.rm=TRUE) <= 0 || max(f^(1/h), na.rm=TRUE) > 1) return(10^6)

    y <- log(sc) + (1 - 1/xi) * log(y) - log(cr)
    y <- rowSums(y, na.rm=TRUE)

    # maximum penalized likelihood estimation

    if(is.null(penk)==F){
      if(penk=="CD"){

        if(xi[1]>=0){ p_k <- 1
        }else if( xi[1] > -1 && xi[1] <0 ){ p_k <- exp(-((1/(1+xi[1])) -1))
        }else if( xi[1] <= -1){ p_k <- 0}

      }else if(penk=="MS"){

        if( xi[1] >= 0.5 || xi[1] <= -0.5) return(10^6)
        p_k <-((0.5+xi[1])^(6-1)) * ((0.5-xi[1])^(9-1)) / beta(6,9) # MS_for k

      }}

    if(is.null(penh)==F){

      B   <-1/(r-1)

      if(penh=="MS"){

        if( h[1] <= -0.5 || h[1] >= B) return(10^6)
        Bef <-function(x) { ((0.5+x)^(6-1)) * ((B-x)^(9-1)) }
        Be  <-stats::integrate(Bef, lower=-0.5, upper= B)[1]$value
        p_h <-((0.5+h[1])^(6-1)) * ((B-h[1])^(9-1)) / Be

      }else if(penh=="MSa"){

        if( h[1] <= -1.2 || h[1] >= B) return(10^6)

        Bef <-function(x) { ((1.2+x)^(6-1)) * ((B-x)^(9-1)) }
        Be  <-stats::integrate(Bef, lower=-1.2, upper= B )[1]$value
        p_h <-((1.2+h[1])^(6-1)) * ((B-h[1])^(9-1)) / Be

      }}

    if(is.null(penk)==T & is.null(penh)==T){ penalty<-0
    }else if(is.null(penk)==F & is.null(penh)==F){ penalty <- r*log(p_k * p_h)
    }else if(is.null(penk)==F){ penalty <- r*log(p_k)
    }else if(is.null(penh)==F){ penalty <- r*log(p_h)}

    sum((r * h - 1) / h * log(f) + y, na.rm=TRUE) - penalty

  }

  # Apply optimization on each set of initial values and retain results
  optim_results <- lapply(init_list, function(init) {
    if (r == 1) {
      if(z$trans==F){stats::optim(init, k4d.lik, hessian=TRUE, method=method, control=list(maxit=maxit))
      }else{suppressWarnings(Rsolnp::solnp(init, k4d.lik, control = list(trace = 0)))}
    } else {
      if(z$trans==F){stats::optim(init, rk4d.lik, hessian=TRUE, method=method, control=list(maxit=maxit))
      }else{suppressWarnings(Rsolnp::solnp(init, rk4d.lik, control = list(trace = 0)))}
    }
  })

  # Collect optimization results and filter out invalid results
  optim_value <- data.frame(
    num = 1:length(optim_results),
    nllh = sapply(optim_results, function(res) {
                  if(z$trans) min(res$values) else res$value}),
    grad = sapply(optim_results, function(res) {
      sum(abs(if (r == 1) numDeriv::grad(k4d.lik, res$par) else numDeriv::grad(rk4d.lik, res$par)))
    })
  )

  optim_value <- optim_value[optim_value$nllh != 10^6, ]
  optim_value <- optim_value[order(optim_value$grad, optim_value$nllh), ]
  best_result <- optim_results[[optim_value$num[1]]]

  # Extract and store the best-fit parameter values
  mu <- drop(mumat %*% (best_result$par[1:npmu]))
  sc <- drop(sigmat %*% (best_result$par[seq(npmu + 1, length = npsc)]))
  xi <- drop(shmat %*% (best_result$par[seq(npmu + npsc + 1, length = npsh)]))
  h  <- drop(hmat %*% (best_result$par[seq(npmu + npsc + nph + 1, length = nph)]))

  # Store results in the output list with dynamic parameter names
  z$r    <- r
  z$conv <- best_result$convergence
  z$nllh <- best_result$value
  z$data <- xdat
  z$mle  <- best_result$par
  z$cov  <- solve(best_result$hessian)
  z$se   <- sqrt(diag(z$cov))
  z$vals <- cbind(mu, sc, xi, h)

  if(show) {
    if(z$trans)
      print(z[c(2, 3)])
    #else print(z[4])
    if(!z$conv)
      print(z[c(4, 6, 8, 10)])
  }

  class(z) <- "rk4d.fit"
  invisible(z)
}
