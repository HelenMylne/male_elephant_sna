library(LaplacesDemon)
library(float)

load('motnp_edgeweights_conditionalprior.RData')
rm(list = ls()[!ls() %in% 'edge_samples'])
x <- object.size(edge_samples) ; format(x, units = 'Mb')

edges_float <- fl(edge_samples)
x <- object.size(edges_float) ; format(x, units = 'Mb')

n <- ncol(edge_samples) #100
test_matrix <- matrix(data = invlogit(rnorm(n*n, mean = 0, sd = 2)),
                      nrow = n, ncol = n)
x <- object.size(test_matrix) ; format(x, units = 'Kb')

test_float <- fl(test_matrix)
x <- object.size(test_float) ; format(x, units = 'Kb')
