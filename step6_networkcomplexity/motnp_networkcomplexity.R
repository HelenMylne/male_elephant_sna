#### set up ####
#### set up ####
# load packages
library(tidyverse)
library(dplyr)
library(rstan)
library(igraph)
library(cmdstanr)
library(bisonR)
library(brms)
library(mclust)

#### load in network model ####
load('motnp_bisonr_edgescalculated_strongprior.RData')

#### network metrics ####
# Network metrics can be extracted using the function extract_metric(bison_model, metric_name). The metric_name variable is made up of two parts: the network feature and the network metric. For example, node_eigen will calculate the eigenvector centrality of each node.

### edge-level properties
# weight
edge_weight <- extract_metric(motnp_edge_weights_strongpriors, "edge_weight")
head(edge_weight)

### node-level properties
# strength
node_strength <- extract_metric(motnp_edge_weights_strongpriors, "node_strength")
head(node_strength)

# eigenvector centrality
node_eigen <- extract_metric(motnp_edge_weights_strongpriors, "node_eigen")
head(node_eigen)

# betweenness
node_betweenness <- extract_metric(motnp_edge_weights_strongpriors, "node_betweenness")
head(node_betweenness)

# closeness
node_closeness <- extract_metric(motnp_edge_weights_strongpriors, "node_closeness")
head(node_closeness)

### network-level properties
# weighted density
global_density <- extract_metric(motnp_edge_weights_strongpriors, "global_density")
head(global_density)

# coefficient of variation of edge weights (aka social differentation)
global_cv <- extract_metric(motnp_edge_weights_strongpriors, "global_cv")
head(global_cv)

# standard deviation of edge weights
global_std <- extract_metric(motnp_edge_weights_strongpriors, "global_std")
head(global_std)

#### network complexity ####
# Edge mixture models can be used to detect different types of social connections by looking for clustering in edge weights. There are two key outputs from a bisonR edge mixture model: a posterior probability distribution over the number of components (clusters, interpreted associal connection types, from 1 to a maximum 𝐾, over the entire network), the probability that each dyad belongs to each possible component.

# Once the model has fit well, we can use the bison_mixture() function to fit the edge mixture model:
motnp_mixture <- bison_mixture(motnp_edge_weights_strongpriors, num_components = 5, verbose = FALSE)
summary(motnp_mixture)

# To get network-level probabilities of the number of components, we can use get_network_component_probabilities() function:
get_network_component_probabilities(motnp_mixture)

# To use probabilities over component membership for a given number of components, we can use the get_edge_component_probabilities() function, where 5 is the number of components we assume to exist in the network:
get_edge_component_probabilities(motnp_mixture, 5) # forcing it to use 5 components it does then spread them out but basically no there is only 1 from the other tests

# It’s often useful to know where the mixture components are located, so we can get their posterior means using the get_component_means() function for a given number of components. Here we’ll calculate the posterior means of the components for the mixture model with 3 components:
(means <- get_component_means(motnp_mixture, 3))
plot(means$`50%`, las = 1, pch = '-', xlab = 'K', ylab = 'component distribution', ylim = c(0,0.1))
for(i in 1:3){
  lines(x = c(i,i), y = c(means$`5%`[i],means$`95%`[i]), col = 'red') # components don't overlap, but they are so close together in their overall magnitude that K = 1 anyway
}
