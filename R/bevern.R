#' Bevern Stream Flow Data
#'
#' Annual r-largest stream flow observations from the Bevern River in the UK.
#' The dataset contains the three largest daily stream flow values recorded
#' in each year.
#'
#' This dataset is commonly used for extreme value analysis based on
#' r-largest order statistics.
#'
#' @format A data frame with 52 rows and 4 columns:
#' \describe{
#'   \item{Year}{Year of observation}
#'   \item{r1}{Largest daily stream flow in the year}
#'   \item{r2}{Second largest daily stream flow}
#'   \item{r3}{Third largest daily stream flow}
#' }
#'
#' @details
#' The data represent annual r-largest daily stream flow observations from
#' the Bevern River. Each row corresponds to one year and contains the
#' three largest daily stream flow measurements recorded in that year.
#'
#' @source
#' United Kingdom hydrological records. This is the original data source
#' containing the daily stream flow observations.
#'
#' @references
#' Shin, Y. and Park, J.-S. (2024).
#' Generalized logistic model for r-largest order statistics,
#' with hydrological application.
#'
#' @examples
#' data(bevern)
#' head(bevern)
#'
"bevern"
