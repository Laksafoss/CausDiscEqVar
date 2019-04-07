
#' Find the Topological Ordering of parameters
#' 
#' Find the Topological Ordering of parameters, which are assumed to come from 
#' a linear SEM (also known as SCM) where it is assumed that the independent 
#' noise terms have equal error variance.
#' 
#' The function \code{top_order} estimated the topological order of the causal
#' graph. The method is based on the theory and pseudo code from Chen, Drton and 
#' Wang (2018). The 4 methods "TD" (Top Down), "BU" (Bottom Up), "HTD" (High 
#' dimensional Top Down) and "HBU" (High dimensional Bottom Up) are implemented
#' in this function.
#' 
#' The "TD" and "BU" methods \emph{only} work when the number of observations 
#' \code{n} is larger then the number of parameters \code{p}. Moreover the
#' consistency results shown in the preprint only hold if \eqn{n > p^2 log(p)}.
#' 
#' The "HTD" and "HBU" both require extra assumption either on the maximal 
#' in-degree or on the maximal size of markov blankets in the graph respectivly. 
#' Hence these methods both has parameters that allows the user to tune 
#' and use different sub procedures that fit the assumptions on the graph. 
#' 
#' The "HBU" needs a tuning parameter \code{M} which is used in the organic 
#' lasso to estimate the conditional variance: 
#' \deqn{\hat{\sigma}^2_{j,\Theta} = \min((1/n)||X_j - X_{V \ (\Theta \cup \{j\})}\beta||_2^2 + 2\lambda||\beta||_1^2)}
#' where \eqn{\lambda = (2M(1/n)log(p))^{1/2}}.
#' 
#' The "HTD" needs both a \code{max.degree} equal to the assumed maximal 
#' in-degree in the graph and a \code{search} parameter. One of three different 
#' \code{search} procedures "full", "B&B" or "OMP" must be used. All three 
#' procedures seek to find a subset of variables of size \code{max.degree} that
#' minimize the conditional variance of the current variable in question given 
#' this subset. 
#' 
#' The "full" search method simply goes through all poissible subsets of size 
#' \code{max.degree}. The "B&B" searched by bound and branch methods. Lastly the
#' "OMP" procedure findes \code{max.degree} variables via forward selection. -
#' 
#' @param X A matrix containing the oberved variables.
#' @param method The estimation method. Possible choises are "TD", "BU", "HTD"
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
#'   pursuit (forward selection) to find \code{max.degree} variables from the 
#'   ansestreal set. 
#'   
#' @return A vector of length equal to the number of parameters with the 
#'   estimated topological ordering is returned.
#'   
#' @seealso The wrapper function \code{\link{graph_est}}combines \code{top_order} 
#'   and the model selection function \code{\link{graph_from_top}} into one 
#'   function that estimates the causal graph from the data. 
#' 
#' @section References:
#' Chen, W., Drton, M. & Wang, Y. S. (2018). On Causal Discovery with Equal 
#' Variance Assumption. \emph{arXiv preprint arXiv:1807.03419.}
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
#' 
#' top_order(X, method = "TD")
#' 
#' top_order(X, method = "BU")
#' 
#' top_order(X, method = "HTD", max.degree = 2L, search = "B&B")
#' 
#' top_order(X, method = "HBU")
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
