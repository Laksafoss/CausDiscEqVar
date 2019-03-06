
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


top_order <- function(X, method = "TD", ...) {
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
  vars <- structure(list(X = X, cov = cov(X)), class = method)
  
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

#if (!is.integer(max.degree)) {
#stop("'max.degreee' must be an integer larger then 1")
#} 
#if (length(max.degree) != 1) {
#  stop("'max.degree' must have length 1.")
#} 
#if (max.degree < 1L) {
#  stop("'max.degree' must be at least 1.")
#}


est_step <- function(vars, theta, j, ...) {
  UseMethod("est_step", vars)
}

est_step.TD <- function(vars, theta, j, ...) {
  set <- c(theta, j)
  return(1 / solve(vars$cov[set,set])[length(set),length(set)])
} 


est_step.HTD <- function(vars, theta, j, max.degree, ...) {
  if (missing(max.degree)) {
    warning("'max.degree' was not specified. Set to 8")
    max.degree <- 8L
  }
  if (length(theta) <= max.degree) {
    set <- c(theta, j)
    return(1 / solve(vars$cov[set,set])[length(set),length(set)])
  } else {
    CC <- rbind(combn(theta, max.degree),j)
    tmp <- sapply(seq_len(ncol(CC)), function(i) {
      C <- CC[,i]
      1 / solve(vars$cov[C,C])[length(C),length(C)]
    })
    return(min(tmp))
  }
} 

est_step.BU <- function(vars, theta, j, ...) {
  set <- c(theta, j)
  return(solve(vars$cov[set,set])[length(set),length(set)])
}


est_step.HBU <- function(vars, theta, j, ...) {
  if (max.degree == 1) {
    stop("This is not done yet")
    if (missing(M)) {
      M <- 2 # TODO what is a good default value ?
    }
    lambda <- 2 * sqrt(2 * M * (1 / nrow(vars$X)) * log(ncol(X)))
    # use stats::nlm perhaps ?
    # the glmnet does not output sd estimates 
    # TODO
  } else { 
    index <- c(theta, j)
    return(solve(vars$cov[index,index])[length(index),length(index)])
  }
}

graph_from_top <- function(X, top) {
  p <- ncol(X)
  tmp <- sapply(top[-1], function(i) {
    above <- top[seq_len(which(top == i)-1)]
    fit <- lars::lars(X[, above, drop = FALSE], X[ , i], intercept = FALSE)
    beta <- coef(fit, s=which.min(fit$Cp))
    B <- rep(0,p)
    B[above] <- beta
    B
  })
  tmp <- cbind(0, tmp)
  tmp[order(top),order(top)]
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
