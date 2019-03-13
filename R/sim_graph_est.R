#' Simulation tool for graph_est
#' 
#' This is an internal simulation function for testing preformance of the 
#' \code{\link{graph_est}} function. 
#' 
#' LONG DESCRIPTION
#' 
#' @param senarios a data frame or matrix
#' @param m the number of simulations of each senario
#' 
#' @return WHAT IT RETURNS
#' 
#' @examples 
#' 
#' 

sim_graph_est <- function(scenarios, m) {
  
  res <- replicate(m, apply(scenarios, 1, function(s) {
    # now s["n"] gives n in the senario !
    # TODO
  })) 
  
}

sim_B <- function(p, graph_setting, l, u) {
  if (p <= 50) {
    if (graph_setting == "sparse") {
      pc <- 3/(2 * p - 2)
    } else if (graph_setting == "dense") {
      pc <- 0.3
    }
    Bsmall <- diag(1, p-1, p-1)
    index <- upper.tri(Bsmall, diag = F)
    Bsmall[index] <- rbinom(sum(index), 1, pc)
    Bsmall[Bsmall == 1] <- sample(c(1,-1), 1) * runif(sum(Bsmall), l, u)
    B <- rbind(cbind(0, Bsmall), 0)
  } else {
    if (graph_setting == "A") {
      # TODO
      B <- matrix(0, p, p)
    } else if (graph_setting == "B") {
      # TODO
      B <- matrix(0,p,p)
    }
  }
  return(B)
}


sim_X <- function(B, n, sigma, alpha = rep(1, ncol(B))) {
  p <- ncol(B)
  N <- matrix(rnorm(n * p, 0, rep(alpha * sigma, each = n)), ncol = p)
  X <- matrix(0, ncol = p, nrow = n)
  for (i in seq_len(p)) {
    X[ ,i] <- X%*%B[,i] + N[,i]
  }
  return(X)
}
