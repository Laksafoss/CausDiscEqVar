
#' Find the Topological Ordering of parameters
#' 
#' Find the Topological Ordering of parameters, which are assumed to come from 
#' a linear SEM (also known as SCM) where it is assumed that the independent 
#' noise terms have equal error variance.
#' 
#' MORE DETAILED DECSRIPTION
#' 
#' @param X A matrix containing the oberved variables.
#' @param method The estimation method. Posible choises are "TD", "BU", "HTD"
#'   and "HBU".
#' @param ... Terms passed to the method specific estimation steps.
#' 
#'   If the method "HBU" is specified it is possible to specify a tuning 
#'   parameter \code{M} to be used in the organic lasso.
#' 
#'   If the method "HTD" is specified it is posible to specify a 
#'   \code{max.degree} and a \code{search}. The \code{max.degree} specifies the 
#'   assumed maximal in-degree in the true causal graph. The parameter 
#'   \code{search} may be set to either "full", "B&B" or "OMP" indicating how 
#'   the algorithem should search for \code{max.degree} number of parameters to 
#'   produce the lowest conditional variance.
#'   
#'   The "full" search method simply looks at all possible subsets of the 
#'   current ansestreal set of size \code{max.degree}, the "B&B" method searches
#'   via a bound and branch method, and lastly "OMP" uses orthogonal matching 
#'   pursuit (fprward selection) to find \code{max.degree} variables from the 
#'   ansestreal set. 
#'   
#' @return A vector of length equal to the number of parameters indicating the 
#'   estimated topological ordering is returned.
#'   
#' @seealso The wrapper function \code{\link{graph_est}}combines \code{top_order} 
#'   and the model selection function \code{\link{graph_from_top}} into one 
#'   function that estimates the causal graph from the data. 
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
#' order <- top_order(X, method = "HTD", max.degree = 2L, search = "B&B")
#' order <- top_order(X, method = "HBU")
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

est_step.HBU <- function(vars, theta, j, M = 0.5, ...) {
  print(M)
  index <- setdiff(seq_len(ncol(vars$X)), c(theta, j))
  if (length(index) == 1) {
    vars$X <- cbind(vars$X, 1)
    index <- c(index, ncol(vars$X))
  }
  lambda <- M * sqrt(log(ncol(vars$X)) / nrow(vars$X))
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
      index <- OMP(vars$X[,j], vars$X[,theta], max.degree)
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
