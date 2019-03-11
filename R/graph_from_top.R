
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
#'               0,0,0,0), ncol = 4, nrow = 4, byrow = TRUE)
#' X <- matrix(0, ncol = 4, nrow = n)
#' for (i in 1:4) {
#'   X[ ,i] <- X %*% B[ ,i] + rnorm(n)
#' }
#' 
#' # we then find the graph using a topological ordering
#' top <- c(1,2,3,4) # from B we know this to be the true ordering
#' graph_from_top(X, top)
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
