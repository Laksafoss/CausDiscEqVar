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
#' @export


sim_graph_est <- function(scenarios, m) {
  nam <- c("p", "graph_setting", "l", "u", "n", "sigma", "TD", "BU", "HTD", "HBU")
  if(! all(nam %in% names(scenarios))) {
    if (class(scenarios) == "character") {
      scenarios <- standard_senarios(scenarios)
    }
  }
  if (any(apply(cbind(-scenarios$l,scenarios$u), 1, sum) <= 0 )) {
    stop("the lower bound 'l' must be smaller then the upper bound 'u' in scenarios")
  }
  large_res <- replicate(m, lapply(seq_len(nrow(scenarios)), function(k) {
    s <- scenarios[k,]
    B <- sim_B(s$p, s$graph_setting, s$l, s$u)
    X <- sim_X(B, s$n, s$sigma)
    order <- seq_len(s$p)
    
    res <- data.frame(method = names(s[7:10])[s[7:10]==TRUE],
                      n = s$n, p = s$p, sigma = s$sigma,
                      Kendall = NA, Recall = NA, Flipped = NA, FDR = NA,
                      row.names = NULL, stringsAsFactors = FALSE)
    
    for (i in seq_len(nrow(res))) {
      top <- top_order(X, method = res[i,"method"], ...)
      Bhat <- graph_from_top(X, top)
      
      tri <- upper.tri(matrix(0,s$p, s$p))
      res[i, "Kendall"] <- (2 / (s$p * (s$p - 1))) * 
        sum(sign(outer(order,order,"-")[tri] * outer(top,top,"-")[tri]))
      res[i, "Recall"] <- round(mean(which(B != 0) %in% which(Bhat != 0)) * 100)
      res[i, "Flipped"] <- round(mean(which(t(Bhat) != 0) %in% which(B != 0)) * 100)
      res[i, "FDR"] <- round(mean(which(Bhat != 0) %in% which(B == 0)) * 100)
    }
    return(res)
  }), simplify = FALSE) 
  
  large_res <- Reduce(rbind, unlist(large_res, recursive = FALSE))
  class(large_res) <- c(class(large_res), "sim_analysis")
  return(large_res)
}


#' @rdname sim_graph_est
#' @export

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
  } else {
    if (graph_setting == "A") {
      Bsmall <- diag(1, p-1, p-1)
      for (j in seq_len(p-1)[-c(1:2)]) {
        possible <- which(rowSums(Bsmall[seq_len(j-1),])<4)
        Bsmall[sample(possible, 2), j] <- 1
      }
    } else if (graph_setting == "B") {
      Bsmall <- diag(1, p-1, p-1)
      for (j in seq_len(p-1)[-c(1:2)]) {
        possible <- which(seq_len(p-1) <= min(j,10))
        Bsmall[sample(possible, 2), j]
      }
    }
  }
  Bsmall[Bsmall == 1] <- sample(c(1,-1), 1) * runif(sum(Bsmall), l, u)
  B <- rbind(cbind(0, Bsmall), 0)
  return(B)
}


#' @rdname sim_graph_est
#' @export

sim_X <- function(B, n, sigma, alpha = rep(1, ncol(B))) {
  p <- ncol(B)
  N <- matrix(rnorm(n * p, 0, rep(alpha * sigma, each = n)), ncol = p)
  X <- matrix(0, ncol = p, nrow = n)
  for (i in seq_len(p)) {
    X[ ,i] <- X%*%B[,i] + N[,i]
  }
  return(X)
}


#' @rdname sim_graph_est
#' @export

standard_senarios <- function(setting) {
  if (setting == "LowDense") {
    scenarios <- expand.grid(p = c(5, 20, 40), graph_setting = c("dense"),
                             l = 0.3, u = 1,
                             n = c(100, 500, 1000), sigma = c(1),
                             TD = TRUE, BU = TRUE, HTD = FALSE, HBU = FALSE)
  } else if (setting == "LowSparse") {
    scenarios <- expand.grid(p = c(5, 20, 40), graph_setting = c("dense"),
                             l = 0.3, u = 1,
                             n = c(100, 500, 1000), sigma = c(1),
                             TD = TRUE, BU = TRUE, HTD = FALSE, HBU = FALSE)
  } else if (setting == "HighA") {
    qq <- c(0.5, 0.75, 1, 1.5, 2)
    scenarios <- data.frame(p = c(qq * 80, qq * 100, qq * 200),
                            graph_setting = c("A"),
                            l = 0.7, u = 1,
                            n = rep(c(80,100,200), each = length(qq)),
                            sigma = 1,
                            TD = FALSE, BU = FALSE, HTD = TRUE, HBU = TRUE)
  } else if (setting == "HighB") {
    qq <- c(0.5, 0.75, 1, 1.5, 2)
    scenarios <- data.frame(p = c(qq * 80, qq * 100, qq * 200),
                            graph_setting = c("B"),
                            l = 0.7, u = 1,
                            n = rep(c(80,100,200),each = length(qq)),
                            sigma = 1,
                            TD = FALSE, BU = FALSE, HTD = TRUE, HBU = TRUE)
  } else {
    stop("this is not a recognized standard simulation setting")
  }
  return(scenarios)
}



#' @rdname  sim_graph_est
#' @export

plot_sim_analysis <- function(data) {
  data$n <- as.factor(data$n)
  
  bb <- ggplot2::scale_y_continuous(position = "right")
  
  PK <- ggplot2::ggplot(cbind(data, dummy = "Kendall"), 
                        ggplot2::aes(x = n, 
                                     y = as.numeric(Kendall), 
                                     color = method)) +
    ggplot2::geom_boxplot(show.legend = FALSE) +
    ggplot2::facet_grid(p ~ dummy, switch = "y")+ 
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    bb
  
  PR <- ggplot2::ggplot(cbind(data, dummy = "Recall"), 
                        ggplot2::aes(x = n, 
                                     y = as.numeric(Recall), 
                                     color = method)) +
    ggplot2::geom_boxplot(show.legend = FALSE) +
    ggplot2::facet_grid(p ~ dummy)+ 
    ggplot2::theme(strip.text.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) +
    bb
  
  PF <- ggplot2::ggplot(cbind(data, dummy = "Flipped"), 
                        ggplot2::aes(x = n, 
                                     y = as.numeric(Flipped), 
                                     color = method)) +
    ggplot2::geom_boxplot(show.legend = FALSE) +
    ggplot2::facet_grid(p ~ dummy) + 
    ggplot2::theme(strip.text.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) +
    bb
  
  PFDR <- ggplot2::ggplot(cbind(data, dummy = "FDR"), 
                          ggplot2::aes(x = n, 
                                       y = as.numeric(FDR), 
                                       color = method)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_grid(p ~ dummy) + 
    ggplot2::xlab("n") + 
    ggplot2::theme(strip.text.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank()) +
    bb
  
  gridExtra::grid.arrange(PK, PR, PF, PFDR, 
                          ncol = 4, widths = c(1.05,1,1,1.35))
}




#' @rdname sim_graph_est
#' @export

make_table <- function(data) {
  Tab <- x %>% 
    dplyr::group_by(p, n , method) %>% 
    dplyr::summarize(
      Kendall = mean(Kendall), 
      Recall = mean(Recall), 
      Flipped = mean(Flipped),
      FDR = mean(FDR))
  
  TKen <- tidyr::spread(Tab[,c("p", "n", "method", "Kendall")], method, Kendall)
  TRec <- tidyr::spread(Tab[,c("p", "n", "method", "Recall")], method, Recall)
  TFli <- tidyr::spread(Tab[,c("p", "n", "method", "Flipped")], method, Flipped)
  TFDR <- tidyr::spread(Tab[,c("p", "n", "method", "FDR")], method, FDR)
  
  A <- dplyr::full_join(TKen, TRec, by = c("n","p"), suffix = c(".Kendall", ".Recall"))
  B <- dplyr::full_join(TFli, TFDR, by = c("n", "p"), suffix  = c(".Flipped", ".FDR"))
  FULL <- dplyr::full_join(A, B, by = c("n", "p"))
  round(FULL, digit = 2)
}
