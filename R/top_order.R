
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
#' order <- top_order(X, method = "TD")
#' order <- top_order(X, method = "BU")
#' order <- top_order(X, method = "HTD")
#' 
#' 
#' @export


top_order <- function(X, method = "TD", ...) {
  p <- ncol(X)
  vars <- structure(list(X = X, cov = stats::cov(X)), class = method)
  
  index <- seq_len(p)
  theta <- numeric(0)
  for (z in seq_len(p-1)) {
    t <- which.min(sapply(index, function(j) {
      est_step(vars, theta, j, ...)
    }))
    theta <- c(theta, index[t])
    index <- index[-t]
  }
  theta <- c(theta, index)
  if (method == "BU") {
    theta <- rev(theta)
  }
  return(theta)
}



est_step <- function(vars, theta, j, ...) {
  UseMethod("est_step", vars)
}

est_step.TD <- function(vars, theta, j, ...) {
  set <- c(theta, j)
  ind <- length(set)
  return(1 / solve(vars$cov[set,set])[ind,ind])
} 


est_step.BU <- function(vars, theta, j, ...) {
  set <- setdiff(seq_len(ncol(vars$cov)), theta)
  ind <- which(set == j)
  return(solve(vars$cov[set,set])[ind,ind])
}

est_step.HBU <- function(vars, theta, j, ...) {
  index <- setdiff(seq_len(ncol(vars$X)), c(theta, j))
  if (length(index) == 1) {
    vars$X <- cbind(vars$X, 1)
    index <- c(index, ncol(vars$X))
  }
  lambda <- 0.5 * sqrt(log(ncol(vars$X)) / nrow(vars$X))
  fit <- natural::olasso_path(vars$X[,index, drop = FALSE], 
                              vars$X[,j], 
                              lambda = c(lambda, lambda * 1.1),
                              intercept = FALSE)
  return(fit$sig_obj[1])
}

est_step.HTD <- function(vars, theta, j,
                         max.degree = 2L, search = "B&B", ...) {
  if (length(theta) <= max.degree) {
    set <- c(theta, j)
    return(1 / solve(vars$cov[set,set])[length(set),length(set)])
  } else {
    if (search == "full") {
      CC <- rbind(utils::combn(theta, max.degree),j)
      tmp <- sapply(seq_len(ncol(CC)), function(i) {
        C <- CC[,i]
        1 / solve(vars$cov[C,C])[length(C),length(C)]
      })
      return(min(tmp))
    } else if (search == "B&B") {
      if (length(theta) >= 50) {
        really.big <- TRUE
      } else {
        really.big <- FALSE
      }
      tmp <- leaps::regsubsets(vars$X[,theta], vars$X[,j], 
                               nvmax = max.degree, intercept = FALSE,
                               really.big = really.big)
      #index <- which.min(tmp$rss)
      #pp <- sum(summary(tmp)$which[index-1, ])
      #return(sqrt(tmp$rss[index])/(tmp$nn - pp))
      return(min(tmp$rss))
    } else if (search == "OMP") {
      index <- OMP(target = vars$X[,j], data = vars$X[,theta], max.degree)
      set <- c(theta[index], j)
      ind <- length(set)
      return(1 / solve(vars$cov[set,set])[ind,ind])
    }
  }
} 



OMP <- function(Y, X, max.degree) {
  d <- ncol(X)
  if (d <= max.degree) {
    return(seq_len(d))
  }
  
  R <- Y
  C <- numeric(0)
  for (i in 1:max.degree) {
    index <- setdiff(seq_len(d), C)
    vals <- abs(t(X[,index]) %*% R)
    C <- c(C, index[which.max(vals)])
    P <- X[,C] %*% solve(t(X[,C]) %*% X[,C]) %*% t(X[,C])
    R <- (diag(ncol(P)) - P) %*% Y
  }
  return(C)
}
