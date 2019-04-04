
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

graph_from_top <- function(X, top, measure = "deviance", which = "1se") {
  if (measure %in% c("deviance", "mse", "mae")) {
    type.measure <- measure
    cvmean <- function(fit) {fit$cvm}
  } else if (measure == c("AIC")) {
    type.measure <- "deviance"
    cvmean <- function(fit) { fit$cvm + 2 * fit$nzero }
  } else if (measure == c("BIC")) {
    type.measure <- "deviance"
    cvmean <- function(fit) {fit$cvm + 
        log(fit$glmnet.fit$nobs) * fit$nzero}
  }
  
  p <- ncol(X)
  
  tmp <- sapply(top[-1], function(i) {
    above <- A <- top[seq_len(which(top == i)-1)]
    if (length(above) == 1) {
      X <- cbind(X, 1)
      A <- c(A, ncol(X))
    }
    fit <- glmnet::cv.glmnet(X[ , A, drop = FALSE], X[ , i], 
                             type.measure = type.measure, 
                             alpha = 1, intercept=FALSE)
    #if (length(fit$cvm) != length(fit$nzero)) {
    #  fit$glmnet.fit$df <- fit$glmnet.fit$df[fit$glmnet.fit$df != 0]
    #  fit$glmnet.fit$beta <- fit$glmnet.fit$beta[ ,fit$glmnet.fit$df != 0]
    
    cvm <- cvmean(fit)
    min <- which.min(cvm)
    if (which == "min") {
      fitindex <- which(fit$glmnet.fit$lambda == fit$lambda[min])
      beta <- fit$glmnet.fit$beta[,fitindex]
    } else if (which == "1se") {
      index <- which.max(cvm[min] + fit$cvsd[min] >= cvm & 
                         cvm[min] - fit$cvsd[min] <= cvm)
      fitindex <- which(fit$glmnet.fit$lambda == fit$lambda[index])
      beta <- fit$glmnet.fit$beta[,fitindex]
    }
    Bcol <- rep(0,p)
    if (length(above) == 1) {
      Bcol[above] <- beta[-length(beta)]
    } else {
      Bcol[above] <- beta  
    }
    Bcol
  })
  tmp <- cbind(0, tmp)
  tmp[order(top),order(top)]
}
