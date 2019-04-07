#' Simulation tool for graph_est
#' 
#' This is an internal simulation function for testing preformance of the 
#' \code{\link{top_order}} and \code{\link{graph_from_top}} function. 
#' 
#' The \code{sim_graph_est} function is used to simulate and evaluate the 
#' simulated data. The \code{kendall} is a function that finds the *Kendall's* 
#' *Tau coefficient* between two topological orderings and is used to evaluate 
#' the preformance of the estimated topological orderings found by 
#' \code{\link{top_order}} in the internal of \code{graph_from_top}. This 
#' coefficient is however only meaningfull when the true ordering is unique.  
#' 
#' In the function \code{sim_graph_est} the three functions \code{kendall}, 
#' \code{\link{sim_B}} and \code{\link{sim_X}} are use. The functions 
#' \code{\link{sim_B}} and \code{\link{sim_X}} simulate the needed data and then
#' the kendall's tau coefficient between the estimated and true topological 
#' ordering is calculated by \code{kendall} along with a few other preformance
#' measures.
#' 
#' @param scenarios a data frame with columns \code{p}, \code{graph_setting}, 
#'   \code{l}, \code{u}, \code{unique_ordering}, \code{n} and \code{sigma}. 
#'   These variables are then used in the functions \code{\link{sim_B}} and 
#'   \code{\link{sim_X}} to simulate a deteset for each row in \code{scenarios}.
#' @param top a data frame with columns \code{method}, \code{max.degree}, 
#'   \code{search} and \code{M}. These variables are used in the function 
#'   \code{\link{top_order}} along with the data sets simulated based on 
#'   \code{scenarios}. 
#' @param graph a data frame with columns \code{measure} and \code{which}. These 
#'   variables are used in the function \code{\link{graph_from_top}} along with
#'   the data sets and their estimated topological orderings. 
#' @param m the number of realizations of each \code{scenario}. 
#' @param true the true topological ordering of the graph
#' @param est one or more estimated topological orderings in a matrix, where 
#'   each column is an estimated topological ordering.
#' 
#' @return 
#'   The \code{sim_graph_est} function returns a data frame with 18 variables and 
#'   \code{ncol(scenarios)} x \code{ncol(top)} x \code{ncol(graph)} 
#'   x m rows. The first 13 variables are simply the input parameters that 
#'   generated that particular output row. The remaning 5 variables describe the 
#'   preformance measures:
#'   
#' * \code{Kendall} Kendalls tau between the true and esimated topological ordering
#' * \code{Hamming} Structual Hamming Distance between true and estimated graph
#' * \code{Recall} percent of arrows true positive arrows in estimated graph
#' * \code{Flipped} percent of arrows flipped compared to true graph
#' * \code{FDR} percent of arrows false positivesin the estimateed graph
#' 
#'   The \code{Kendall} function returnes one numeric value between -1 and 1. 
#' 
#' @seealso Both \code{\link{sim_B}} and \code{\link{sim_X}} are used in the 
#' interior of \code{sim_graph_est}.
#' 
#' @md
#' @examples 
#' 
#' \dontrun{
#' 
#' ####  The two low dimensional settings from the article
#' 
#' scenarios <- expand.grid(p = c(5, 20, 40),
#'                          graph_setting = c("dense", "sparse"),
#'                          l = 0.3,
#'                          u = 1,
#'                          unique_ordering = TRUE,
#'                          n = c(100, 500, 1000),
#'                          sigma = 1,
#'                          stringsAsFactors = FALSE)
#'                         
#' top <- data.frame(method = c("TD", "BU"),
#'                   max.degree = NA,
#'                   search = NA,
#'                   M = NA,
#'                   stringsAsFactors = FALSE)
#'                   
#' graph <- data.frame(measure = "deviance",
#'                     which = "1se",
#'                     stringAsFactors = FALSE)
#' 
#' SIM <- sim_graph_est(scenarios, top, graph, m = 500)                  
#' 
#' 
#' 
#' 
#' ####  The two high dimensional settings from the article
#' 
#' scenarios <- data.frame(p = rep(c( 40,  60,  80, 120, 160,
#'                                    50,  75, 100, 150, 200,
#'                                   100, 150, 200, 300, 400), 2),
#'                          graph_setting = rep(c("A", "B"), each = 15),
#'                          l = 0.7,
#'                          u = 1,
#'                          unique_ordering = TRUE,
#'                          n = rep(rep(c(80, 100, 200), each = 5), 2),
#'                          sigma = 1,
#'                          stringsAsFactors = FALSE)
#'                         
#' top <- data.frame(method = c("HTD", "HBU"),
#'                   max.degree = c(3L, NA),
#'                   search = c("B&B", NA),
#'                   M = c(NA, 0.5),
#'                   stringsAsFactors = FALSE)
#'                   
#' graph <- data.frame(measure = NA,
#'                     which = NA,
#'                     stringAsFactors = FALSE)
#' 
#' SIM <- sim_graph_est(scenarios, top, graph, m = 50)                  
#' }
#' 
#' 
#' 
#' @export

sim_graph_est <- function(scenarios, top, 
                          graph = data.frame(measure = NA, which = NA), m) {
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
  
  find_graph <- ifelse(any(is.na(graph[1,])), FALSE, TRUE)

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
                  max.degree = as.numeric(t["max.degree"]),
                  search = t["search"],
                  M = as.numeric(t["M"]))
      })
      
      toptop <- cbind(top, "Kendall" = kendall(order, top_est))
      
      if (find_graph) {
        graph_est <- unlist(lapply(1:nrow(graph), function(g) {
          lapply(1:ncol(top_est), function(t) {
            graph_from_top(X, top_est[,t], 
                           measure = graph[g,"measure"], 
                           which = graph[g,"which"])
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
            "Recall" = mean(E %in% Ehat) * 100,
            "Flipped" = mean(FEhat %in% E) * 100,
            "FDR" = mean(Ehat %in% nE) * 100
          )
        })
      } else { # find_graph == FALSE
        small_res <- c(
          "Hamming" = NA,
          "Recall" = NA,
          "Flipped" = NA,
          "FDR" = NA
        )
      }
      
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


#' Simulates a graph
#' 
#' This is an internal simulation function. It creates a graph with regression 
#' coefficients.
#' 
#' The \code{sim_B} function simulates a causal graph with regression 
#' coefficients. The graph is stored as matrix, where the (i,j)th entry in the
#' matrix is non zero only if there is an arrow from j to i in the graph. The 
#' value in the (i,j)th entry corresponds to the direct causal effect from j
#' to i.  
#' 
#' @param p number of variables to simulate
#' @param graph_setting either "dense", "sparse", "A", or "B". These settings 
#'   correspons to the settings described in the article.
#' @param l lower bound for direct abolute causal effect
#' @param u upper bound for direct absolute causal effect
#' @param unique_ordering if \code{TRUE} the simulated graph will have a unique 
#'   topological ordering
#' 
#' @return A matrix encoding a graph.
#' 
#' @examples 
#' sim_B(5, "dense", 0.3, 1)
#' 
#' sim_B(10, "A", 0.7, 1)
#' 
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


#' Simulates data from a graph
#' 
#' This is an internal simulation function for simulating data from a graph. 
#' 
#' The \code{sim_B} function simulates a causal graph with regression 
#' coefficients. The graph must be encoded as an upper triangular matrix.
#' 
#' @param B a upper triangular matrix encoding a graph 
#' @param n number of observations from the graph
#' @param sigma the common standard error
#' @param alpha a vector of length equal to the number of parameters given by 
#'   \code{B}. These \code{alpha}s are multiplied to the common \code{sigma}.
#' 
#' @return The \code{sim_X} function returns a matrix with data from the graph 
#' \code{B} with columns equal to the number of columns in \code{B} and rows 
#' equal to \code{n}.
#' 
#' @examples 
#' (B <- matrix(c(0,0,0,1,0,0,1,1,0), ncol = 3))
#' 
#' sim_X(B, 10, 1)
#' 
#' sim_X(B, 10, 1, c(1, 5, 10))
#' 
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
