% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/top_order.R
\name{top_order}
\alias{top_order}
\title{Find the Topological Ordering of parameters}
\usage{
top_order(X, method = "TD", ...)
}
\arguments{
\item{X}{A matrix containing the oberved variables.}

\item{method}{The estimation method. Possible choises are "TD", "BU", "HTD"
and "HBU".}

\item{...}{Terms passed to the method specific estimation steps.

  If the method "HBU" is specified it is possible to specify a tuning 
  parameter \code{M} to be used in the organic lasso.

  If the method "HTD" is specified it is posible to specify a 
  \code{max.degree} and a \code{search}. The \code{max.degree} specifies the 
  assumed maximal in-degree in the true causal graph. The parameter 
  \code{search} may be set to either "full", "B&B" or "OMP" indicating how 
  the algorithem should search for \code{max.degree} number of parameters to 
  produce the lowest conditional variance.
  
  The "full" search method simply looks at all possible subsets of the 
  current ansestreal set of size \code{max.degree}, the "B&B" method searches
  via a bound and branch method, and lastly "OMP" uses orthogonal matching 
  pursuit (forward selection) to find \code{max.degree} variables from the 
  ansestreal set.}
}
\value{
A vector of length equal to the number of parameters with the 
  estimated topological ordering is returned.
}
\description{
Find the Topological Ordering of parameters, which are assumed to come from 
a linear SEM (also known as SCM) where it is assumed that the independent 
noise terms have equal error variance.
}
\details{
The function \code{top_order} estimated the topological order of the causal
graph. The method is based on the theory and pseudo code from Chen, Drton and 
Wang (2018). The 4 methods "TD" (Top Down), "BU" (Bottom Up), "HTD" (High 
dimensional Top Down) and "HBU" (High dimensional Bottom Up) are implemented
in this function.

The "TD" and "BU" methods \emph{only} work when the number of observations 
\code{n} is larger then the number of parameters \code{p}. Moreover the
consistency results shown in the preprint only hold if \eqn{n > p^2 log(p)}.

The "HTD" and "HBU" both require extra assumption either on the maximal 
in-degree or on the maximal size of markov blankets in the graph respectivly. 
Hence these methods both has parameters that allows the user to tune 
and use different sub procedures that fit the assumptions on the graph.

The "HBU" needs a tuning parameter \code{M} which is used in the organic 
lasso to estimate the conditional variance: 
\deqn{\hat{\sigma}^2_{j,\Theta} = \min((1/n)||X_j - X_{V \ (\Theta \cup \{j\})}\beta||_2^2 + 2\lambda||\beta||_1^2)}
where \eqn{\lambda = (2M(1/n)log(p))^{1/2}}.

The "HTD" needs both a \code{max.degree} equal to the assumed maximal 
in-degree in the graph and a \code{search} parameter. One of three different 
\code{search} procedures "full", "B&B" or "OMP" must be used. All three 
procedures seek to find a subset of variables of size \code{max.degree} that
minimize the conditional variance of the current variable in question given 
this subset.

The "full" search method simply goes through all poissible subsets of size 
\code{max.degree}. The "B&B" searched by bound and branch methods. Lastly the
"OMP" procedure findes \code{max.degree} variables via forward selection. -
}
\section{References}{

Chen, W., Drton, M. & Wang, Y. S. (2018). On Causal Discovery with Equal 
Variance Assumption. \emph{arXiv preprint arXiv:1807.03419.}
}

\examples{
n <- 1000
p <- 5
sigma <- 1

# an example where all nodes in the graph are root nodes
X <- matrix(rnorm(n * p, 0, sigma), ncol = p)
order <- top_order(X)

# an example where pa(X_i) = X_j for all j < i
for (i in seq_len(p)[-1]) {
  X[ ,i] <- X[ ,i-1] + X[ ,i]
}

top_order(X, method = "TD")

top_order(X, method = "BU")

top_order(X, method = "HTD", max.degree = 2L, search = "B&B")

top_order(X, method = "HBU")
}
\seealso{
The wrapper function \code{\link{graph_est}}combines \code{top_order} 
  and the model selection function \code{\link{graph_from_top}} into one 
  function that estimates the causal graph from the data.
}
