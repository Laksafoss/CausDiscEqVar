#' Estimate graph from observational data
#' 
#' SHORT DESCRIPTION
#' 
#' LONG DESCRIPTION
#' 
#' @param X a matrix or data frame containing the observed variables.
#' @param method
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
#'               0,0,0,0), ncol = 4, nrow = 4, byrow = T)
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
  B <- graph_from_top(X, top) 
  
  res <- structure(list(top_order = top, graph = B), class = "graph_est")
  return(res)
}






