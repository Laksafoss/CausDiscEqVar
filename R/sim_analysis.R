
# -------  ANALYSIS  -----------------------------------------------------------

library(ggplot2)


Low <- wrapper_analysis(m = 5, 
                        n = c(100,500,1000), 
                        p = c(5,20,40), 
                        pc = 0.3, 
                        sigma = 10)

# Kendall's tau
ggplot(Y, aes(as.factor(p), Kendall)) + geom_boxplot(aes(col = as.factor(n)))

# Recall
ggplot(Y, aes(as.factor(p), Recall)) + geom_boxplot(aes(col = as.factor(n)))

# Flipped
ggplot(Y, aes(as.factor(p), Flipped)) + geom_boxplot(aes(col = as.factor(n)))

# Fals Discovery Rate
ggplot(Y, aes(as.factor(p), FDR)) + geom_boxplot(aes(col = as.factor(n)))
