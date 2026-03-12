#' Bangkok Rainfall Data
#'
#' Annual top five daily rainfall events recorded in Bangkok, Thailand,
#' from 1961 to 2018. The dataset contains the five largest daily rainfall
#' amounts observed each year.
#'
#' The data are commonly used for extreme value analysis based on
#' r-largest order statistics.
#'
#' @format A data frame with 58 rows and 5 columns:
#' \describe{
#'   \item{X1}{Largest daily rainfall in the year (mm)}
#'   \item{X2}{Second largest daily rainfall (mm)}
#'   \item{X3}{Third largest daily rainfall (mm)}
#'   \item{X4}{Fourth largest daily rainfall (mm)}
#'   \item{X5}{Fifth largest daily rainfall (mm)}
#' }
#'
#' @details
#' Each row corresponds to one year from 1961 to 2018 and contains the
#' five largest daily rainfall observations recorded in that year.
#'
#' @source
#' Rain gauge station records from Bangkok, Thailand.
#'
#' @references
#' Shin, Y and Park, J-S. (2023). *Modeling climate extremes using the four-parameter kappa distribution for r-largest order statistics*.
#'
#' @examples
#' data(bangkok)
#' head(bangkok)
#'
"bangkok"
