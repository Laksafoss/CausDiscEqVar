
library(CausDiscEqVar)
library(ggplot2)
library(dplyr)

B <- simB(10, 0.3, 0.3, 1)
X <- simX(B, n = 1000)

topTD <- top_order(X, method = "TD")
BhatTD <- graph_from_top(X, topTD)

topHTD <- top_order(X, method = "HTD", max.degree = 5L)
BhatHTD <- graph_from_top(X, topHTD)

topBU <- top_order(X, method = "BU")
BhatBU <- graph_from_top(X, topBU)

#topHBU <- top_order(X, method = "HBU")
#BhatHBU <- graph_from_top(X, topHBU)

A <- analysis(p = 5, 0.3, 0.3, 1, n = 1000)
A$small_res
A$large_res$simulated$order
A$large_res$TD$order
A$large_res$BU$order

##  ==  Low-dimensional sparse settings  ===================================  ##
grid <- expand.grid(p = c(5,20,40), pc = 0.3, l = 0.3, u = 1,
                    n = c(100,500,1000), sigma = 1)
Low <- wrapper_analysis(m = 3, grid = grid, method = c("TD", "BU"))
make_table(Low)

# Kendall's tau
ggplot(Low, aes(as.factor(p), Kendall)) + 
  geom_boxplot(aes(col = as.factor(n))) + 
  facet_wrap(aes(as.factor(method)))

# Recall
ggplot(Low, aes(as.factor(p), Recall)) + 
  geom_boxplot(aes(col = as.factor(n))) + 
  facet_wrap(aes(as.factor(method)))

# Flipped
ggplot(Low, aes(as.factor(p), Flipped)) + 
  geom_boxplot(aes(col = as.factor(n))) + 
  facet_wrap(aes(as.factor(method)))

# Fals Discovery Rate
ggplot(Low, aes(as.factor(p), FDR)) + 
  geom_boxplot(aes(col = as.factor(n))) + 
  facet_wrap(aes(as.factor(method)))





##  ==  High-dimensional setting with maximum in-degree q = 2  =============  ##
3/(2*80-2)
A <- analysis(p = 80, pc = 3/(2*80-2), l = 0.7, u = 1, n = 80, method = "HTD")
A$small_res

grid <- lapply(c(80,100), function(n) {
  f <- c(0.5, 1, 1.5, 2)
  lapply(f * n, function(p) {
    expand.grid(p = p, pc = 3 / (2 * p - 2), l = 0.7, u = 1, n = n, sigma = 1)
  })
})
grid <- Reduce(rbind, unlist(grid, recursive = F))
High <- wrapper_analysis(m = 2, grid = grid, method = c("HTD")) # add HBU


# Kendall's tau
ggplot(Y, aes(as.factor(p), Kendall)) + geom_boxplot(aes(col = as.factor(n)))

# Recall
ggplot(Y, aes(as.factor(p), Recall)) + geom_boxplot(aes(col = as.factor(n)))

# Flipped
ggplot(Y, aes(as.factor(p), Flipped)) + geom_boxplot(aes(col = as.factor(n)))

# Fals Discovery Rate
ggplot(Y, aes(as.factor(p), FDR)) + geom_boxplot(aes(col = as.factor(n)))

