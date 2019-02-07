
#' Find the Topological Ordering of parameters
#' 
#' Find the Topological Ordering of parameters, which are assumed to come from 
#' a linear SEM (also known as SCM) with equal error variance.
#' 
#' MORE DETAILED DECSRIPTION
#' 
#' @param X a matrix or data frame containing the oberved variables.
#' @param max.degree an integer larger then 1 describing the maximum in-degree 
#'   in the underlying graph.
#'   
#' @return TOP ORDER RETURNS !!
#' 
#' @examples 
#' n <- 1000
#' p <- 5
#' sigma <- 1
#' 
#' # an example where all nodes in the graph are root nodes
#' X <- matrix(rnorm(n * p, 0, sigma), ncol = p)
#' order <- top_order(X)
#' 
#' # an example where pa(X_i) = X_j for all j < i
#' for (i in seq_len(p)[-1]) {
#'   X[ ,i] <- X[ ,i-1] + X[ ,i]
#' }
#' order <- top_order(X)
#' 
#' 
#' @export


top_order <- function(X, method = "TD", max.degree = 8L) {
  if (!is.numeric(max.degree)) {
    stop("'max.degreee' must be an integer larger then 1")
  } 
  if (length(max.degree) != 1) {
    stop("'max.degree' must have length 1.")
  } 
  if (max.degree < 1L) {
    stop("'max.degree' must be at least 1.")
  }
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  if (!is.matrix(X)) {
    stop("'X' must be a matrix or data.frame")
  }
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
  }
  
  cov <- cov(X)

  # TODO:
  # how should we specify the method ? From the article we have 
  # - Top Down
  # - Bottom Up
  # - Top Down in p>n setting with set max degree
  # - Bottom Up in p>n setting for sparse graphs 
  
  # We might also want to consider implementing
  # - Greedy search from Peters and Bühlmann 2014
  # - Ghoshal and Honorio 2018
  # - what about Po-Ling and Bühlmann 2014 ???
  # but do these method need an odering stage ????
  
  # TODO:
  # When method specification is done - how do we implwment this in a nice way

  
  
}
