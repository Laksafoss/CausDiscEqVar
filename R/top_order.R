
#' Find the Topological Ordering of parameters
#' 
#' Find the Topological Ordering of parameters, which are assumed to come from 
#' a linear SEM (also known as SCM) with equal error variance.
#' 
#' MORE DETAILED DECSRIPTION
#' 
#' @param X a matrix containing the oberved variables.
#' @param method the estimation method. Posible choises are TD, BU, HTD
#' @param ... terms passed to the method specific estimation steps. If the 
#'   method TD is specified it is posible to specify a \code{max.degree}. 
#'   
#' @return a vector of length equal to the number of parameters indicating the 
#'   estimated topological ordering
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


top_order <- function(X, method = "TD", ...) {
  p <- ncol(X)
  vars <- structure(list(X = X, cov = stats::cov(X)), class = method)
  
  index <- seq_len(p)
  theta <- numeric(0)
  for (z in seq_len(p)) {
    t <- which.min(sapply(index, function(j) {
      est_step(vars, theta, j, ...)
    }))
    theta <- c(theta, index[t])
    index <- index[-t]
  }
  if (method == "BU") { # add also HBU later
    theta <- rev(theta)
  }
  return(theta)
}



est_step <- function(vars, theta, j, ...) {
  UseMethod("est_step", vars)
}

est_step.TD <- function(vars, theta, j, ...) {
  set <- c(theta, j)
  return(1 / solve(vars$cov[set,set])[length(set),length(set)])
} 


est_step.BU <- function(vars, theta, j, ...) {
  set <- c(theta, j)
  return(solve(vars$cov[set,set])[length(set),length(set)])
}


est_step.HTD <- function(vars, theta, j, max.degree, ...) {
  if (missing(max.degree)) {
    warning("'max.degree' was not specified. Set to 2")
    max.degree <- 2L
  }
  if (length(theta) <= max.degree) {
    set <- c(theta, j)
    return(1 / solve(vars$cov[set,set])[length(set),length(set)])
  } else {
    CC <- rbind(utils::combn(theta, max.degree),j)
    tmp <- sapply(seq_len(ncol(CC)), function(i) {
      C <- CC[,i]
      1 / solve(vars$cov[C,C])[length(C),length(C)]
    })
    return(min(tmp))
  }
} 

