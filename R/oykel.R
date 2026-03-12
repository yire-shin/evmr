#' Oykel River Stream Flow Data
#'
#' Annual r-largest daily stream flow observations from the Oykel River
#' in the United Kingdom. The dataset contains the three largest daily
#' stream flow values recorded in each year.
#'
#' The data are used for extreme value analysis based on
#' r-largest order statistics models.
#'
#' @format A data frame with 42 rows and 4 variables:
#' \describe{
#'   \item{Year}{Year of observation}
#'   \item{r1}{Largest daily stream flow in the year}
#'   \item{r2}{Second largest daily stream flow}
#'   \item{r3}{Third largest daily stream flow}
#' }
#'
#' @details
#' Each row represents one year and contains the three largest
#' daily stream flow observations recorded in that year.
#' Missing observations are represented by \code{NA}.
#'
#' @source
#' United Kingdom hydrological records. This is the original data source
#' containing the daily stream flow data.
#'
#' @references
#' Shin, Y. and Park, J.-S. (2025).
#' Generalized Gumbel model for r-largest order statistics,
#' with hydrological application.
#'
#' @examples
#' data(oykel)
#' head(oykel)
#'
"oykel"
