
#' Data - Article results
#' 
#' A dataset containing the reported results from the article.
#' 
#' @format a data frame with 18 variables:
#' \describe{
#'   \item{p}{number of parameters}
#'   \item{graph_setting}{Low dimensional sparse, dense or high dimensional A/B}
#'   \item{l}{lower bound for absolute ausal effect}
#'   \item{u}{upper bound for sbaolute causal effect}
#'   \item{uniquq_ordering}{TRUE / FALSE}
#'   \item{n}{number of observations}
#'   \item{sigma}{standard deviation of the common error}
#'   \item{measure}{Model selection parameter, "mse", "mae", "deviance", "AIC", "BIC"}
#'   \item{which}{Model selection parameter, "min" or "1se"}
#'   \item{method}{"TD", "BU", "HTD", "HBU"}
#'   \item{max.degree}{HTD parameter}
#'   \item{search}{"full", "B&B" or "OMP"}
#'   \item{M}{Tuning parameter in organic lasso}
#'   \item{Kendall}{Kendalls tau btween the true and esimated topological ordering}
#'   \item{Hamming}{Structual Hamming Distance between true and estimated graph}
#'   \item{Recall}{percent of arrows true positive arrows in estimated graph}
#'   \item{Flipped}{percent of arrows flipped compared to true graph}
#'   \item{FDR}{percent of arrows false positivesin the estimateed graph}
#' } 

"ArticleRes"

# ArticleRes <- data.frame(
#   p = c(rep(rep(c(5,5,5,20,20,20,40,40,40),2),4),
#         rep(c(40,60,80,120,160, 50,75,100,150,200,100,150,200,300,400),4)),
#   graph_setting = c(rep(rep(c("dense", "sparse"), each = 18),2),
#                     rep(c("A", "B"), each = 30)),
#   l = c(rep(0.3, 36), rep(0.1, 36), rep(0.7, 60)),
#   u = 1,
#   unique_ordering = TRUE,
#   n = c(rep(rep(c(100,500,1000,100,500,1000,100,500,1000),2),4),
#         rep(rep(c(80,100,200), each = 5),4)),
#   sigma = 1,
#   measure = "mse",
#   which = "min",
#   method = c(rep(rep(c("TD", "BU"), each = 9),4),
#              rep(rep(c("HTD", "HBU"), each = 15),2)),
#   max.degree = c(rep(NA, 72),rep(3,60)),
#   search = c(rep(NA, 72), rep(rep(c("B&B", NA), each = 15),2)),
#   M = c(rep(NA, 72), rep(rep(c(NA, 0.5), each = 15),2)),
#   Kendall = c(0.85, 0.98, 0.99, 0.92, 0.99, 1.00, 0.96, 0.99, 1.00,
#               0.82, 0.97, 0.98, 0.85, 0.97, 0.99, 0.91, 0.98, 0.99,
#               0.87, 0.98, 0.99, 0.77, 0.96, 0.99, 0.72, 0.96, 0.99,
#               0.84, 0.96, 0.98, 0.59, 0.88, 0.94, 0.44, 0.80, 0.91,
#               rep(NA, 36),
#               0.99, 0.98, 0.95, 0.84, 0.72,
#               1.00, 0.99, 0.97, 0.86, 0.73,
#               1.00, 1.00, 0.99, 0.87, 0.74,
#               0.89, 0.89, 0.87, 0.83, 0.73,
#               0.93, 0.92, 0.87, 0.84, 0.78,
#               0.95, 0.90, 0.79, 0.74, 0.64,
#               1.00, 0.99, 0.95, 0.77, 0.55,
#               1.00, 1.00, 0.97, 0.74, 0.63,
#               1.00, 1.00, 0.99, 0.80, 0.65,
#               0.70, 0.52, 0.39, 0.25, 0.16,
#               0.70, 0.50, 0.38, 0.26, 0.12,
#               0.77, 0.61, 0.48, 0.20, 0.13),
#   Hamming = c(rep(NA, 36),
#               1.3, 0.7, 0.5, 31, 22, 28, 170, 152, 136,
#               1.3, 0.7, 0.5, 32, 22, 28, 174, 155, 137,
#               1.6, 0.8, 0.6, 7, 3.5, 2.2, 14, 7, 5,
#               1.7, 0.9, 0.6, 7, 3.5, 2.2, 15, 7, 5,
#               rep(NA, 60)),
#   Recall = c(91, 99, 99, 85, 99, 100, 71, 96, 97,
#              89, 98, 99, 83, 98, 100, 69, 96, 97,
#              91, 98, 99, 85, 98, 100, 81, 98, 99,
#              89, 98, 99, 79, 96, 98, 72, 94, 98,
#              73, 80, 85, 73, 91, 94, 66, 93, 96,
#              73, 80, 84, 73, 91, 94, 65, 93, 95,
#              74, 85, 88, 69, 85, 90, 64, 84, 90,
#              73, 84, 88, 69, 84, 90, 63, 84, 89,
#              rep(NA, 60)),
#   Flipped = c(7, 1, 1, 3,  1, 0, 2,  0, 0,
#               9, 1, 1, 5,  1, 0, 3,  1, 0,
#               6, 1, 1, 9,  2, 0, 10, 2, 1,
#               7, 2, 1, 13, 4, 2, 16, 5, 2,
#               7, 4, 3, 4,  2, 2,  2, 2, 1,
#               7, 4, 3, 3,  3, 2,  3, 2, 1,
#               8, 3, 3, 4,  4, 3,  3, 3, 3,
#               8, 4, 4, 4,  4, 2,  4, 3, 3,
#               rep(NA, 60)),
#   FDR = c(17, 4, 3, 32, 28, 26, 41, 41, 40,
#           18, 4, 3, 35, 29, 26, 43, 42, 41,
#           16, 5, 3, 35, 19, 14, 38, 24, 17,
#           17, 5, 4, 40, 22, 16, 46, 31, 22,
#           16, 8, 5, 27, 24, 21, 36, 38, 36,
#           15, 7, 5, 28, 24, 21, 37, 39, 36,
#           18, 7, 6, 16,  9,  5, 16,  8,  6,
#           18, 7, 6, 17,  8,  5, 18,  7,  6,
#           rep(NA,60))
# )
# 
# usethis::use_data(ArticleRes, ArticleRes, overwrite = TRUE)




#' Data - Low dimensional simulation
#' 
#' A dataset containing the result of a large scale simulation of the low 
#' dimensional setting.
#' 
#' @format a data frame with 18 variables:
#' \describe{
#'   \item{p}{number of parameters}
#'   \item{graph_setting}{Low dimensional sparse, dense or high dimensional A/B}
#'   \item{l}{lower bound for absolute ausal effect}
#'   \item{u}{upper bound for sbaolute causal effect}
#'   \item{uniquq_ordering}{TRUE / FALSE}
#'   \item{n}{number of observations}
#'   \item{sigma}{standard deviation of the common error}
#'   \item{measure}{Model selection parameter, "mse", "mae", "deviance", "AIC", "BIC"}
#'   \item{which}{Model selection parameter, "min" or "1se"}
#'   \item{method}{"TD", "BU", "HTD", "HBU"}
#'   \item{max.degree}{HTD parameter}
#'   \item{search}{"full", "B&B" or "OMP"}
#'   \item{M}{Tuning parameter in organic lasso}
#'   \item{Kendall}{Kendalls tau btween the true and esimated topological ordering}
#'   \item{Hamming}{Structual Hamming Distance between true and estimated graph}
#'   \item{Recall}{percent of arrows true positive arrows in estimated graph}
#'   \item{Flipped}{percent of arrows flipped compared to true graph}
#'   \item{FDR}{percent of arrows false positivesin the estimateed graph}
#' } 
#' 
"LowSim"

#usethis::use_data(LowSim, LowSim, overwrite = TRUE)





#' Data - High dimensional simulation
#' 
#' A dataset containing the result of a large scale simulation of the high 
#' dimensional setting.
#' 
#' @format a data frame with 18 variables:
#' \describe{
#'   \item{p}{number of parameters}
#'   \item{graph_setting}{Low dimensional sparse, dense or high dimensional A/B}
#'   \item{l}{lower bound for absolute ausal effect}
#'   \item{u}{upper bound for sbaolute causal effect}
#'   \item{uniquq_ordering}{TRUE / FALSE}
#'   \item{n}{number of observations}
#'   \item{sigma}{standard deviation of the common error}
#'   \item{measure}{Model selection parameter, "mse", "mae", "deviance", "AIC", "BIC"}
#'   \item{which}{Model selection parameter, "min" or "1se"}
#'   \item{method}{"TD", "BU", "HTD", "HBU"}
#'   \item{max.degree}{HTD parameter}
#'   \item{search}{"full", "B&B" or "OMP"}
#'   \item{M}{Tuning parameter in organic lasso}
#'   \item{Kendall}{Kendalls tau btween the true and esimated topological ordering}
#'   \item{Hamming}{Structual Hamming Distance between true and estimated graph}
#'   \item{Recall}{percent of arrows true positive arrows in estimated graph}
#'   \item{Flipped}{percent of arrows flipped compared to true graph}
#'   \item{FDR}{percent of arrows false positivesin the estimateed graph}
#' } 
#' 
"HighSim"


#usethis::use_data(HighSim, HighSim, overwrite = TRUE)

