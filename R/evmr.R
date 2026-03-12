#' Fit and Compare r-Largest Order Statistics Models
#'
#' Fit multiple extreme value models for r-largest order statistics and
#' return a combined summary table including parameter estimates,
#' standard errors, and return levels.
#'
#' @param data A vector, matrix, or data frame containing r-largest order
#' statistics.
#' @param models Character vector specifying models to fit.
#' @param num_inits Number of random initial values used in optimization.
#'
#' @return A data frame summarizing fitted models.
#' @export
#'
#' @examples
#' \dontrun{
#' data(bangkok)
#' evmr(bangkok)
#' }
evmr <- function(data,
                 models = c("rk4d","rglo","rggd","rgd","rld"),
                 num_inits = 100) {

  allowed <- c("rk4d","rglo","rggd","rgd","rld")

  if (!all(models %in% allowed)) {
    stop("models must be one of: ", paste(allowed, collapse = ", "))
  }

  if (is.vector(data)) {
    data <- matrix(data, ncol = 1)
  }

  common_cols <- c(
    "model", "r", "nllh",
    "mu", "sigma", "xi", "h",
    "mu.se", "sigma.se", "xi.se", "h.se",
    "rl20", "rl20.se",
    "rl50", "rl50.se",
    "rl100", "rl100.se",
    "rl200", "rl200.se"
  )

  standardize_result <- function(x, model_name) {
    x <- as.data.frame(x)

    x$model <- model_name

    for (nm in common_cols) {
      if (!nm %in% names(x)) {
        x[[nm]] <- NA_real_
      }
    }

    x$model <- as.character(x$model)

    x <- x[, common_cols]

    x
  }

  out <- list()

  if ("rk4d" %in% models) {
    tmp <- rk4d.summary(data, num_inits = num_inits)
    out[[length(out) + 1]] <- standardize_result(tmp, "rk4d")
  }

  if ("rglo" %in% models) {
    tmp <- rglo.summary(data, num_inits = num_inits)
    out[[length(out) + 1]] <- standardize_result(tmp, "rglo")
  }

  if ("rggd" %in% models) {
    tmp <- rggd.summary(data, num_inits = num_inits)
    out[[length(out) + 1]] <- standardize_result(tmp, "rggd")
  }

  if ("rgd" %in% models) {
    tmp <- rgd.summary(data, num_inits = num_inits)
    out[[length(out) + 1]] <- standardize_result(tmp, "rgd")
  }

  if ("rld" %in% models) {
    tmp <- rld.summary(data, num_inits = num_inits)
    out[[length(out) + 1]] <- standardize_result(tmp, "rld")
  }

  result <- do.call(rbind, out)
  rownames(result) <- NULL
  result
}
