

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
  order <- seq_len(p)
  #order <- sample(seq_len(ncol(X)))
  #B <- B[order, order]
  #X <- X[ ,order]
  
  res <- data.frame(
    expand.grid(method = method, stringsAsFactors = F),
    n = n,
    p = p,
    sigma = sigma,
    Kendall = NA,
    Recall = NA,
    Flipped = NA,
    FDR = NA)
  
  large_res <- list(simulated = list(order = order, B = B))
  
  for (i in seq_len(nrow(res))) {
    top <- top_order(X, method = res[i,"method"], ...)
    Bhat <- graph_from_top(X, top)
    
    large_res[[res[i,"method"]]] <- list(order = top, B = Bhat)
    
    tri <- upper.tri(matrix(0,p,p))
    res[i, "Kendall"] <- (2 / (p * (p - 1))) * 
      sum(sign(outer(order,order,"-")[tri] * outer(top,top,"-")[tri]))
    res[i, "Recall"] <- round(mean(which(B != 0) %in% which(Bhat != 0)) * 100)
    res[i, "Flipped"] <- round(mean(which(t(Bhat) != 0) %in% which(B != 0)) * 100)
    res[i, "FDR"] <- round(mean(which(Bhat != 0) %in% which(B == 0)) * 100) 
  }
  return(structure(list(small_res = res, large_res = large_res), 
                   class = "top_analysis"))
}


wrapper_analysis <- function(m, grid, method, ...) {
  W <- replicate(m, mapply(function(p,pc,l,u,n,sigma) {
    analysis(p = p, pc = pc, l = l, u = u, n = n, 
             sigma = sigma, method = method)$small_res
  }, p = grid$p, pc = grid$pc, l = grid$l, u = grid$u, 
  n = grid$n, sigma = grid$sigma, SIMPLIFY = FALSE), simplify = FALSE)
  Reduce(rbind, unlist(W, recursive = FALSE))
}


make_table <- function(x) {
  Tab <- x %>% 
    group_by(p, n , method) %>% 
    summarize(
      Kendall = mean(Kendall), 
      Recall = mean(Recall), 
      Flipped = mean(Flipped),
      FDR = mean(FDR))
  
  TKen <- spread(Tab[,c("p", "n", "method", "Kendall")], method, Kendall)
  TRec <- spread(Tab[,c("p", "n", "method", "Recall")], method, Recall)
  TFli <- spread(Tab[,c("p", "n", "method", "Flipped")], method, Flipped)
  TFDR <- spread(Tab[,c("p", "n", "method", "FDR")], method, FDR)
  
  A <- full_join(TKen, TRec, by = c("n","p"), suffix = c(".Kendall", ".Recall"))
  B <- full_join(TFli, TFDR, by = c("n", "p"), suffix  = c(".Flipped", ".FDR"))
  FULL <- full_join(A, B, by = c("n", "p"))
  round(FULL, digit = 2)
}
