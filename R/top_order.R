
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
    CC <- rbind(combn(theta, max.degree),j)
    tmp <- sapply(seq_len(ncol(CC)), function(i) {
      C <- CC[,i]
      1 / solve(vars$cov[C,C])[length(C),length(C)]
    })
    return(min(tmp))
  }
} 





#' Find graph from topological ordering
#' 
#' SHORT DESCRIPTION
#' 
#' LONG DESCRIPTION
#' 
#' @param  X a matrix containing the observed variables
#' @param top a topological ordering of the variables
#' 
#' @return the B graph estimated from the data
#' 
#' @examples 
#' 
#' # we create some data from the graph B
#' n <- 1000
#' B <- matrix(c(0,1,0,1,
#'               0,0,2,0,
#'               0,0,0,1,
#'               0,0,0,0), ncol = 4, nrow = 4, byrow = T)
#' X <- matrix(0, ncol = 4, nrow = n)
#' for (i in 1:4) {
#'   X[ ,i] <- X %*% B[ ,i] + rnorm(n)
#' }
#' 
#' # we then find the graph using a topological ordering
#' top <- c(1,2,3,4) # from B we know this to be the true ordering
#' grap_from_top(X, top)
#' 
#' @export

graph_from_top <- function(X, top) {
  p <- ncol(X)
  tmp <- sapply(top[-1], function(i) {
    above <- top[seq_len(which(top == i)-1)]
    if (length(above) == 1) {
      fit <- glmnet::cv.glmnet(cbind(1,X[ ,above, drop = FALSE]), X[ ,i], 
                               type.measure = "mse", alpha = 1, intercept=FALSE)
      beta <- fit$glmnet.fit$beta[-1,which(fit$lambda == fit$lambda.1se)]
    } else {
      fit <- glmnet::cv.glmnet(X[,above, drop = FALSE], X[ ,i], 
                               type.measure = "mse", alpha = 1, intercept=FALSE)
      beta <- fit$glmnet.fit$beta[,which(fit$lambda == fit$lambda.1se)]
    }
    B <- rep(0,p)
    B[above] <- beta
    B
  })
  tmp <- cbind(0, tmp)
  tmp[order(top),order(top)]
}
