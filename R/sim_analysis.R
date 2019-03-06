
library(CausDiscEqVar)
library(ggplot2)

B <- simB(10, 0.3, 0.3, 1)
X <- simX(B, n = 1000)

topTD <- top_order(X, method = "TD")
graph_from_top(X, topTD)

topHTD <- top_order(X, method = "HTD", max.degree = 5L)
graph_from_top(X, topHTD)

topBU <- top_order(X, method = "BU")
graph_from_top(X, topBU)

#topHBU <- top_order(X, method = "HBU")
#graph_from_top(X, topHBU)


##  ==  Low-dimensional sparse settings  ===================================  ##

Low <- wrapper_analysis(m = 10,   # strange error - not in raw code, as that works fine above
                        n = c(100,500,1000), 
                        p = c(5,20,40), pc = 0.3, l = 0.3, u = 1,
                        sigma = 1, method = c("TD", "BU"))


# Kendall's tau
ggplot(Low, aes(as.factor(p), Kendall)) + geom_boxplot(aes(col = as.factor(n)))

# Recall
ggplot(Low, aes(as.factor(p), Recall)) + geom_boxplot(aes(col = as.factor(n)))

# Flipped
ggplot(Low, aes(as.factor(p), Flipped)) + geom_boxplot(aes(col = as.factor(n)))

# Fals Discovery Rate
ggplot(Low, aes(as.factor(p), FDR)) + geom_boxplot(aes(col = as.factor(n)))





##  ==  High-dimensional setting with maximum in-degree q = 2  =============  ##

High <- wrapper_analysis(m = 5, 
                         n = c(100,500,1000), 
                         p = c(5,20,40), pc = 0.3, l = 0.3, u = 1,
                         sigma = 10, method = c("TD"),
                         max.degree = c(2L))

# Kendall's tau
ggplot(Y, aes(as.factor(p), Kendall)) + geom_boxplot(aes(col = as.factor(n)))

# Recall
ggplot(Y, aes(as.factor(p), Recall)) + geom_boxplot(aes(col = as.factor(n)))

# Flipped
ggplot(Y, aes(as.factor(p), Flipped)) + geom_boxplot(aes(col = as.factor(n)))

# Fals Discovery Rate
ggplot(Y, aes(as.factor(p), FDR)) + geom_boxplot(aes(col = as.factor(n)))

