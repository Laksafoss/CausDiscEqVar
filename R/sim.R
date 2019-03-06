

simB <- function(p, pc, l, u) {
  if (missing(pc)) {pc <- 3/(2*p - 2)}
  if (missing(l)) {l <- 0.1}
  if (missing(u)) {u <- 1}
  if(l >= u) {stop("'l' must be smaller then 'u'")}
  
  Bsmall <- diag(1, p-1, p-1)
  index <- upper.tri(Bsmall, diag = F)
  Bsmall[index] <- rbinom(sum(index), 1, pc)
  Bsmall[Bsmall == 1] <- sample(c(1,-1), 1) * runif(sum(Bsmall), l, u)
  rbind(cbind(0, Bsmall), 0)
}


simX <- function(B, n, sigma = 1, alpha = rep(1, ncol(B))) {
  p <- ncol(B)
  N <- matrix(rnorm(n * p, 0, rep(alpha * sigma, each = n)), ncol = p)
  X <- matrix(0, ncol = p, nrow = n)
  for (i in seq_len(p)) {
    X[ ,i] <- X%*%B[,i] + N[,i]
  }
  return(X)
} 


analysis <- function(p, pc, l = 0.3, u = 1, n = 100, sigma = 1, 
                     alpha = rep(1, ncol(B)), method = c("TD","BU"), ...) {
  B <- simB(p, pc, l, u)
  X <- simX(B, n, sigma, alpha)
  order <- sample(seq_len(ncol(X)))
  B <- B[order, order]
  X <- X[ ,order]
  
  res <- data.frame(
    expand.grid(method = method, stringsAsFactors = F),
    n = n,
    p = p,
    sigma = sigma,
    Kendall = NA,
    Recall = NA,
    Flipped = NA,
    FDR = NA)
  for (i in seq_len(nrow(res))) {
    top <- top_order(X, method = res[i,"method"], ...)
    Bhat <- graph_from_top(X, top)
    
    res[i, "Kendall"] <- 0 #TODO
    res[i, "Recall"] <- round(mean(which(B != 0) %in% which(Bhat != 0)) * 100)
    res[i, "Flipped"] <- round(mean(which(t(Bhat) != 0) %in% which(B != 0)) * 100)
    res[i, "FDR"] <- round(mean(which(Bhat != 0) %in% which(B == 0)) * 100) 
  }
  return(res)
}


wrapper_analysis <- function(m, n, p, pc, l, u, sigma, method, ...) {
  Z <- expand.grid(n = n, p = p, pc = pc, l = l, u = u, sigma = sigma)
  W <- replicate(m, mapply(function(n, p, pc, l, u, sigma){
    analysis(p = p, pc = pc, l = l, u = u, 
             n = n, sigma = sigma, method = method, ...)
  }, n = Z$n, p = Z$p, pc = Z$pc, l = Z$l, u = Z$u, sigma = Z$sigma,
  SIMPLIFY = FALSE))
  Reduce(rbind, W)
}
