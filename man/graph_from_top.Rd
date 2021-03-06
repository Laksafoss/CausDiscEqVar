% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_from_top.R
\name{graph_from_top}
\alias{graph_from_top}
\title{Find graph from topological ordering}
\usage{
graph_from_top(X, top, measure = "deviance", which = "1se")
}
\arguments{
\item{X}{A matrix containing the observed variables}

\item{top}{A topological ordering of the variables}

\item{measure}{Either "mse", "mae"}

\item{which}{Either "min" or "1se"}
}
\value{
the coefficient graph estimated from the data. This is a graph which 
  has values different from zero iff there is an arrow in the estimated 
  causal graph. The values of the non-zero entries in this matrix are the 
  estimated (causal) effects.
}
\description{
Given data and a topological ordering of the causal graph this function will
estimate the causal graph itself.
}
\details{
To estimate the graph, the parents of each variable is fund via model 
selection. As the underlying grap is assumed to be a DAG the parents of a 
variable \eqn{X_{(j)}} must be found among the variables with a lower 
ordering \eqn{(X_{(i)})_{i<j}}.

All model selection is done via cross-validation lasso.
}
\examples{

# we create some data from the graph B
n <- 1000
B <- matrix(c(0,1,0,1,
              0,0,2,0,
              0,0,0,1,
              0,0,0,0), ncol = 4, nrow = 4, byrow = TRUE)
X <- matrix(0, ncol = 4, nrow = n)
for (i in 1:4) {
  X[ ,i] <- X \%*\% B[ ,i] + rnorm(n)
}

# we then find the graph using a topological ordering
top <- c(1,2,3,4) # from B we know this to be the true ordering
graph_from_top(X, top)

}
\seealso{
The wrapper function \code{\link{graph_est}} combines the function
  \code{\link{top_order}}, which estimates the topological ordering of the 
  causal graph from data, and \code{graph_from_top} into one function that 
  estimates the causal graph from the data.
}
