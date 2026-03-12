#' @name rgd.prof
#' @aliases rgd.prof
#' @title profile likelihood for rgd
#' @param z An object returned by \code{rggdfit}. The object should represent a stationary model.
#' @param m The return level (i.e.\ the profile likelihood is for the value that is exceeded with probability 1/\code{m}).
#' @param xlow,xup The least and greatest value at which to evaluate the profile likelihood.
#' @param conf The confidence coefficient of the plotted profile confidence interval.
#' @param nint The number of points at which the profile likelihood is evaluated.
#'
#' @return plot of the profile likelihood is produced, with a horizontal line representing a profile confidence interval with confidence coefficient \code{conf}.
#' @export
#'
#' @examples
#' data(bangkok)
#' \dontrun{
#' z<-rgd.fit(bangkok)
#' rgd.Prof(z,100,xlow=180,xup=230)
#' }
rgd.Prof <- function(z, m, xlow, xup, conf = 0.95, nint = 100) {

  if (m <= 1)
    stop("`m` must be greater than one")
  cat("If routine fails, try changing plotting interval", fill = TRUE)

  r <- z$r
  p <- 1 - (1 / m)
  v <- numeric(nint)
  x <- seq(xlow, xup, length = nint)
  sol <- c(z$mle[2])
  rl  <- rgd.rl(z,m)$rl

  rgd.plik <- function(a) {

    zr <- as.matrix(z$data[, r],ncol=1)

    mu <- xp + a[1]*log(-log(p))

    if (any(a[1] <= 0,na.rm=T)) return(10^6)

    y <- (z$data - mu)/a[1]

    y <- y+log(a[1])
    y <- rowSums(y, na.rm = TRUE)

    sum( exp(-(zr-mu)/a[1]) + y)
  }

  for (i in 1:nint) {
    xp <- x[i]
    opt <- try(suppressWarnings(stats::optim(sol, rgd.plik, control = list(trace = 0))))
    if (inherits(opt, "try-error")) {
      v[i] <- 10^6
    } else {
      v[i] <- opt$value
      sol <- opt$par
    }
  }

  d <- data.frame(x = x, v = -v)

  cent <- which.min(abs(round(x) - round(rl)))

  # Initialize an empty list to store results for each conf value
  result_list <- list()

  # Loop through each conf to calculate w for each case and extract ci1, ci2
  for (i in seq_along(conf)) {
    current_chisq <- -z$nllh - 0.5 * stats::qchisq(conf[i], 1)
    current_diff <- d$v - current_chisq

    # Find ci1 and ci2
    ci1 <- d$x[which.min(abs(current_diff[1:cent]))]
    ci2 <- d$x[cent:length(current_diff)][which.min(abs(current_diff[cent:length(current_diff)]))]
    cidiff <- ci2 - ci1

    # Append results to the list
    result_list[[i]] <- data.frame(
      year = m,
      returnlevel = rl,
      conf = conf[i],
      ci1 = ci1,
      ci2 = ci2,
      cidiff = cidiff,
      row.names = ""
    )
  }

  # Combine all results into a single data frame
  w_df <- do.call(rbind, result_list)

  # Print the combined data frame
  print(w_df)

  # Plot the profile likelihood and confidence intervals
  plot(d$x, d$v, type = "l", xlab = "Return level", xlim = c(xlow, xup), las=1,
       ylab = "Profile Log-likelihood", main = sprintf("Profile Likelihood for CI of %syear Return level in the rGD (r = %s)", m, z$r))
  ma <- -z$nllh

  col <- c("blue", "red", "darkgreen")

  threshold<-list()
  ci       <-list()


  # Draw horizontal lines for each confidence level
  # Draw vertical lines for ci1 and ci2 for each confidence level

  graphics::abline(v = rl, col = "gray", lty = 2)

  for (i in seq_along(conf)) {
    threshold[[i]] <- -z$nllh - 0.5 * stats::qchisq(conf[i], 1)
    ci[[i]]        <- result_list[[i]]

    graphics::segments(x0 = ci[[i]]$ci1, x1= ci[[i]]$ci2, y0=threshold[[i]], y1 = threshold[[i]], col = col[i], lty = 1)
    graphics::segments(x0 = ci[[i]]$ci1, x1= ci[[i]]$ci2, y0=threshold[[i]], y1 = threshold[[i]], col = col[i], lty = 1)

    graphics::abline(v = ci[[i]]$ci1, col = col[i], lty = 2)
    graphics::abline(v = ci[[i]]$ci2, col = col[i], lty = 2)
  }

  # Automatically position legend in the top right corner
  legend_labels <- c("Return level", sprintf("Conf: %.2f", conf))
  legend_colors <- c("gray", col[seq_along(conf)])

  graphics::legend(
    "topright",                         # Legend position
    legend = legend_labels,             # Labels based on conf values
    col = legend_colors,                # Colors for each conf
    lty = 1,                            # Line style
    #title = "",                        # Title for legend
    bg = "white",                       # Background color
    cex = 0.8                           # Adjust text size
  )

  invisible(w_df)
}
