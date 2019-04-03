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

sim_graph_est <- function(scenarios, top, graph, m) {
  scnam <- c("p", "graph_setting", "l", "u", "unique_ordering", 
             "n", "sigma")
  topnam <- c("method", "max.degree", "search")
  graphnam <- c("measure", "which")

  if (! (is.data.frame(scenarios) & 
         is.data.frame(top) & is.data.frame(graph))) {
    stop("'scenarios', 'top' and 'graph' must all be data.frames")
  }
  if (! all(scnam %in% names(scenarios))) {
    stop("One or more parameters are missing from 'scenarios'")
  } 
  if (! all(topnam %in% names(top))) {
    stop("One or more parameters are missing from 'top'")
  }
  if (! all(graphnam %in% names(graph))) {
    stop("One or more parameters are missing from 'graph'")
  }

  res <- data.frame()
  
  for (k in 1:m) {
    cat("Loop ", k, ": ")
    
    for (i in 1:nrow(scenarios)) {
      cat(i)
      s <- scenarios[i,]
      B <- sim_B(s$p, s$graph_setting, s$l, s$u, s$unique_ordering)
      X <- sim_X(B, s$n, s$sigma)
      order <- seq_len(s$p)
      
      top_est <- apply(top, 1, function(t) {
        top_order(X, 
                  method = t["method"], 
                  max.degree = t["max.degree"],
                  search = t["search"])
      })
      
      toptop <- cbind(top, "Kendall" = kendall(order, top_est))
      
      graph_est <- unlist(lapply(1:nrow(graph), function(i) {
        lapply(1:ncol(top_est), function(j) {
          graph_from_top(X, top_est[,j], measure = graph[i,"measure"], which = graph[i,"which"])
        })
      }), recursive = FALSE)

      small_res <- sapply(graph_est, function(Bhat) {
        E <- which(B != 0) # True edges
        nE <- which(B == 0) # True non-edges
        Ehat <- which(Bhat != 0) # estimated edges
        nEhat <- which(Bhat == 0) # estimated non-edges
        FEhat <- which(t(Bhat) != 0) # estimated edges flipped
        
        if (length(E) == 0) { E <- 0}
        if (length(nE) == 0) { nE <- 0}
        if (length(Ehat) == 0) { Ehat <- 0}
        if (length(nEhat) == 0) { nEhat <- 0}
        if (length(FEhat) == 0) { FEhat <- 0}
        
        c(
          "Hamming" = sum(E %in% nEhat, Ehat %in% nE, FEhat %in% E),
          "Recall" = mean(Ehat %in% E),
          "Flipped" = mean(FEhat %in% E),
          "FDR" = mean(Ehat %in% nE)
        )
      })
      
      reps <- rep(1:nrow(graph), each = nrow(top))
      tmp <- data.frame(s, graph[reps, ], toptop, t(small_res), row.names = NULL)
      res <- rbind(res, tmp)
      cat(", ")
    } 
    cat("\n")
  } 
  class(res) <- c(class(res), "sim_analysis")
  return(res)
}


#' @rdname sim_graph_est
#' 
kendall <- function(true, est) {
  p <- length(true)
  tri <- upper.tri(matrix(0, p, p))
  apply(est, 2, function(e) {
    if (p != length(e)) {
      stop("The two orderings must have same length")
    }
    (2 / (p * (p - 1))) * 
      sum(sign(outer(true, true, "-")[tri] * outer(e, e, "-")[tri]) )
  })
}


#' @rdname sim_graph_est
#' @export
sim_B <- function(p, graph_setting, l, u, unique_ordering = TRUE) {
  if (graph_setting %in% c("sparse", "dense")) {
    if (graph_setting == "sparse") {
      pc <- 3/(2 * p - 2)
    } else if (graph_setting == "dense") {
      pc <- 0.3
    }
    Bsmall <- diag(1, p-1, p-1)
    if (unique_ordering) {
      index <- upper.tri(Bsmall, diag = FALSE)
    } else {
      index <- upper.tri(Bsmall, diag = TRUE)
      pc <- (2 - 2 * pc) / p + pc # should we do this ??
    }
    Bsmall[index] <- rbinom(sum(index), 1, pc)
  } else if (graph_setting %in% c("A", "B")) {
    if (unique_ordering) { 
      Bsmall <- diag(1, p-1, p-1)
      rm <- 1:2; ch <- 2
    } else {
      Bsmall <- diag(0, p-1, p-1)
      rm <- 1:3; ch <- 3
    }
    if (graph_setting == "A") {
      for (j in seq_len(p-1)[-rm]) {
        possible <- which(rowSums(Bsmall[seq_len(j-1),]) < 4)
        Bsmall[sample(possible, ch), j] <- 1
      }
    } else if (graph_setting == "B") {
      for (j in seq_len(p-1)[-rm]) {
        possible <- which(seq_len(p-1) <= min(j,10))
        Bsmall[sample(possible, ch), j] <- 1
      }
    }
  }
  
  Bsmall[Bsmall == 1] <- sample(c(1,-1), sum(Bsmall), replace = TRUE) * 
    runif(sum(Bsmall), l, u)
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
