#' Data - Low dimensional dense simulation
#' 
#' A dataset containing the result of a large scale simulation of the low 
#' dimensional and dense setting.
#' 
#' @format a data frame with 8 variables:
#' \describe{
#'   \item{method}{the topological ordering method}
#'   \item{n}{number of observations}
#'   \item{p}{number of parameters}
#'   \item{sigma}{stnadard deviation of the common error}
#'   \item{Kendall}{Kendalls tau btween the true and esimated topological ordering}
#'   \item{Recall}{percent of arrows true positive arrows in estimated graph}
#'   \item{Flipped}{percent of arrows flipped compared to true graph}
#'   \item{FDR}{percent of arrows false positivesin the estimateed graph}
#' } 
#' 
"LowDense"

#set.seed(299972)
#scen <- CausDiscEqVar:::standard_senarios("LowDense")
#LowDense <- CausDiscEqVar::sim_graph_est(scen, 50)
#usethis::use_data(LowDense, LowDense, overwrite = TRUE)


#' Data - Low dimensional sparse simulation
#' 
#' A dataset containing the result of a large scale simulation of the low 
#' dimensional and sparse setting.
#' 
#' @format a data frame with 8 variables:
#' \describe{
#'   \item{method}{the topological ordering method}
#'   \item{n}{number of observations}
#'   \item{p}{number of parameters}
#'   \item{sigma}{stnadard deviation of the common error}
#'   \item{Kendall}{Kendalls tau btween the true and esimated topological ordering}
#'   \item{Recall}{percent of arrows true positive arrows in estimated graph}
#'   \item{Flipped}{percent of arrows flipped compared to true graph}
#'   \item{FDR}{percent of arrows false positivesin the estimateed graph}
#' } 
#' 
"LowSparse"

#set.seed(700386)
#scen <- CausDiscEqVar:::standard_senarios("LowSparse")
#LowSparse <- CausDiscEqVar::sim_graph_est(scen, 50)
#usethis::use_data(LowSparse, LowSparse)
