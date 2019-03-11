#' Estimate graph from observed data  
#' 
#' Estimates the graph of a linear SEM with assumed equal variance of the noise 
#' terms. 
#' 
#' LONG DESCRIPTION
#' 
#' @param X a matrix or data frame containing the observed variables.
#' @param method the estimation method. Posible choises are TD, BU, HTD
#' @param ... terms passed to the method specific estimations steps.
#' 
#' @return WHAT IS RETURNED 
#' 
#' @examples 
#' 
#' # we create some data from the graph B
#' n <- 1000
#' B <- matrix(c(0,1,0,1,
#'               0,0,2,0,
#'               0,0,0,1,
#'               0,0,0,0), ncol = 4, nrow = 4, byrow = TRUE)
#' X <- matrix(0, ncol = 4, nrow = n)
#' for (i in 1:4) {
#'   X[ ,i] <- X %*% B[ ,i] + rnorm(n)
#' }
#' 
#' # from the simulated data we etimate 
#' graph_est(X, method = "TD)
#' 
#' @export

graph_est <- function(X, method = "TD", ...) {
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  if (!is.matrix(X)) {
    stop("'X' must be a matrix or data.frame")
  }
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
  }
  
  top <- top_order(X, method = "TD", ...)
  B <- G <- graph_from_top(X, top) 
  G[G!=0] <- 1
  
  res <- structure(list(top_order = top, graph = G, B = B), class = "graph_est")
  return(res)
}






