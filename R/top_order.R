
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


top_order <- function(X, method = "TD", max.degree = 8L, ...) {
  if (!is.integer(max.degree)) {
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
  
  p <- ncol(X)
  cov <- cov(X)
  vars <- structure(list(X = X, cov = cov(X)), class = method)

  theta <- numeric(0)
  for (z in seq_len(p)) {
    est <- sapply(seq_len(p)[-theta], function(j) {
      est_step(vars, theta, j, max.degree, ...)
    })
    theta <- c(theta, min(est))
  }
  return(theta)
}


est_step <- function(vars, theta, j, max.degree, ...) {
  UseMethod("est_step", vars)
}

est_step.TD <- function(vars, theta, j, max.degree, ...) {
  if (length(theta) > max.degree) {
    C <- combn(theta, max.degree)
    tmp <- sapply(seq_len(ncol(C)), function(i) {
      1 / solve(vars$cov[c(C[ ,i], j) , c(C[ ,i], j)])[j,j]
    })
    return(min(tmp))
  } else {
    return(1 / solve(vars$cov[c(theta, j) , c(theta, j)])[j,j])
  }
} 

est_step.BU <- function(vars, theta, j, max.degree, ...) {
  if (max.degree) {
    if (missing(M)) {
      M <- 2 # TODO what is a good default value ?
    }
    lambda <- 2 * sqrt(2 * M * (1 / nrow(vars$X)) * log(ncol(X)))
    # use stats::nlm perhaps ?
    # the glmnet does not output sd estimates 
    # TODO
  } else {
    return(solve(vars$cov[c(theta, j) , c(theta, j)])[j,j])
  }
}




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
