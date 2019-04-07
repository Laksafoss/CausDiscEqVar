#' Estimate graph from observed data  
#' 
#' Estimates the graph of a linear SEM with assumed equal variance of the noise 
#' terms. 
#' 
#' \code{graph_est} is a simple wrapper function which given data will fist 
#' estimate the topological ordering of the causal graph using the function 
#' \code{\link{top_order}}, and then estimate the causal graph itself using the 
#' function \code{\link{graph_from_top}}. 
#' 
#' @inheritParams graph_from_top
#' @inheritParams top_order
#' 
#' @return The \code{graph_est} function returns a vector \code{top_order} 
#'   with the estimated topological ordering of the causal graph, a matrix 
#'   \code{G} of the estimated causal graph, and lastly a matrix \code{B} of the 
#'   estimated regression matrix. 
#' 
#' @seealso \code{graph_est} is simply a wrapper functio for the two functions
#'   \code{\link{top_order}} and \code{\link{graph_from_top}} which givem data 
#'   will estimate the topological ordering of the causal graph and estimates 
#'   the causal graph itself given an ordering respectivly.
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
#' # from the simulated data we etimate in two different way to illustrate
#' graph_est(X, method = "TD", measure = "deviance", which = "1se")
#' graph_est(X, method = "HTD", 
#'           measure = "deviance", which = "1se", 
#'           max.degree = 2L, search = "full")
#' 
#' @export

graph_est <- function(X, method = "TD", 
                      measure = "deviance", which = "1se", ...) {
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
  B <- G <- graph_from_top(X, top, measure, which) 
  G[G != 0] <- 1
  
  res <- structure(list(top_order = top, graph = G, B = B), 
                   class = "graph_est")
  return(res)
}






