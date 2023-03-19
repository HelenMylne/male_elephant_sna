#### set up ####
# load packages
#library(tidyverse) ; library(bisonR)
library(tidyverse, lib.loc = '../packages/') # library(tidyverse)
library(bisonR, lib.loc = '../packages/')    # library(bisonR)

#### short windows ####
max_time_window <- 5
outputs <- data.frame(time_window = 1:max_time_window,
                      mean_edge = NA, sd_edge = NA,
                      median_edge = NA, q2.5_edge = NA, q97.5_edge = NA,
                      mean_strength = NA, sd_strength = NA,
                      median_strength = NA, q2.5_strength = NA, q97.5_strength = NA,
                      mean_eigen = NA, sd_eigen = NA,
                      median_eigen = NA, q2.5_eigen = NA, q97.5_eigen = NA,
                      mean_density = NA, sd_density = NA,
                      median_density = NA, q2.5_density = NA, q97.5_density = NA,
                      mean_cv = NA, sd_cv = NA,
                      median_cv = NA, q2.5_cv = NA, q97.5_cv = NA,
                      mean_std = NA, sd_std = NA,
                      median_std = NA, q2.5_std = NA, q97.5_std = NA,
                      num_components = NA,
                      prob_k1 = NA, prob_k2 = NA, prob_k3 = NA, prob_k4 = NA, prob_k5 = NA)
pdf('../outputs/mpnpshort_allwindows_components.pdf')

for( time_window in 1:max_time_window){
  #### load in network model ####
  load(paste0('mpnp_edgecalculations/mpnpshort',time_window,'_bisonr_edgescalculated.RData'))
  
  #### network metrics ####
  ### edge-level properties
  # weight
  edge_weight <- extract_metric(mpnp_edge_weights, "edge_weight")
  saveRDS(edge_weight, paste0('../data_processed/mpnp_network_properties/mpnpshort',
                              time_window,'_edgeweights.RDS'))
  outputs$mean_edge[time_window] <- mean(edge_weight)
  outputs$sd_edge[time_window] <- sd(edge_weight)
  outputs$median_edge[time_window] <- median(edge_weight)
  outputs$q2.5_edge[time_window] <- quantile(edge_weight, 0.025)
  outputs$q97.5_edge[time_window] <- quantile(edge_weight, 0.975)
  
  ### node-level properties
  # strength
  node_strength <- extract_metric(mpnp_edge_weights, "node_strength")
  saveRDS(node_strength, paste0('../data_processed/mpnp_network_properties/mpnpshort',
                                time_window,'_nodestrengths.RDS'))
  outputs$mean_strength[time_window] <- mean(node_strength)
  outputs$sd_strength[time_window] <- sd(node_strength)
  outputs$median_strength[time_window] <- median(node_strength)
  outputs$q2.5_strength[time_window] <- quantile(node_strength, 0.025)
  outputs$q97.5_strength[time_window] <- quantile(node_strength, 0.975)
  
  # eigenvector centrality
  node_eigen <- extract_metric(mpnp_edge_weights, "node_eigen")
  saveRDS(node_eigen, paste0('../data_processed/mpnp_network_properties/mpnpshort',
                             time_window,'_nodeeigenvectors.RDS'))
  node_strength <- extract_metric(mpnp_edge_weights, "node_strength")
  saveRDS(node_strength, paste0('../data_processed/mpnp_network_properties/mpnpshort',
                                time_window,'_nodestrengths.RDS'))
  outputs$mean_eigen[time_window] <- mean(node_eigen)
  outputs$sd_eigen[time_window] <- sd(node_eigen)
  outputs$median_eigen[time_window] <- median(node_eigen)
  outputs$q2.5_eigen[time_window] <- quantile(node_eigen, 0.025)
  outputs$q97.5_eigen[time_window] <- quantile(node_eigen, 0.975)
  
  ### network-level properties
  # weighted density
  global_density <- extract_metric(mpnp_edge_weights, "global_density")
  saveRDS(global_density, paste0('../data_processed/mpnp_network_properties/mpnpshort',
                                 time_window,'_globaldensity.RDS'))
  outputs$mean_density[time_window] <- mean(global_density)
  outputs$sd_density[time_window] <- sd(global_density)
  outputs$median_density[time_window] <- median(global_density)
  outputs$q2.5_density[time_window] <- quantile(global_density, 0.025)
  outputs$q97.5_density[time_window] <- quantile(global_density, 0.975)
  
  # coefficient of variation of edge weights (aka social differentation)
  global_cv <- extract_metric(mpnp_edge_weights, "global_cv")
  saveRDS(global_cv, paste0('../data_processed/mpnp_network_properties/mpnpshort',
                            time_window,'_globalcv.RDS'))
  outputs$mean_cv[time_window] <- mean(global_cv)
  outputs$sd_cv[time_window] <- sd(global_cv)
  outputs$median_cv[time_window] <- median(global_cv)
  outputs$q2.5_cv[time_window] <- quantile(global_cv, 0.025)
  outputs$q97.5_cv[time_window] <- quantile(global_cv, 0.975)
  
  # standard deviation of edge weights
  global_std <- extract_metric(mpnp_edge_weights, "global_std")
  saveRDS(global_std, paste0('../data_processed/mpnp_network_properties/mpnpshort',
                             time_window,'_globalstd.RDS'))
  outputs$mean_std[time_window] <- mean(global_std)
  outputs$sd_std[time_window] <- sd(global_std)
  outputs$median_std[time_window] <- median(global_std)
  outputs$q2.5_std[time_window] <- quantile(global_std, 0.025)
  outputs$q97.5_std[time_window] <- quantile(global_std, 0.975)
  
  #### network complexity ####
  # Edge mixture models can be used to detect different types of social connections by looking for clustering in edge weights. There are two key outputs from a bisonR edge mixture model: a posterior probability distribution over the number of components (clusters, interpreted associal connection types, from 1 to a maximum 𝐾, over the entire network), the probability that each dyad belongs to each possible component.
  
  # Once the model has fit well, we can use the bison_mixture() function to fit the edge mixture model:
  mpnp_mixture <- bison_mixture(mpnp_edge_weights, num_components = 5, verbose = FALSE)
  summary(mpnp_mixture)
  
  # To get network-level probabilities of the number of components, we can use get_network_component_probabilities() function:
  num_components <- get_network_component_probabilities(mpnp_mixture)
  outputs$num_components[time_window] <- num_components$num_components[which(num_components$probability == max(num_components$probability))]
  outputs$prob_k1[time_window] <- num_components$probability[1]
  outputs$prob_k2[time_window] <- num_components$probability[2]
  outputs$prob_k3[time_window] <- num_components$probability[3]
  outputs$prob_k4[time_window] <- num_components$probability[4]
  outputs$prob_k5[time_window] <- num_components$probability[5]
  
  # It’s often useful to know where the mixture components are located, so we can get their posterior means using the get_component_means() function for a given number of components. Here we’ll calculate the posterior means of the components for the mixture model with 3 components:
  num_components <- num_components$num_components[which(num_components$probability == max(num_components$probability))]
  (means <- get_component_means(mpnp_mixture, num_components))
  plot(means$`50%`, las = 1, pch = '-', xlab = 'K', ylab = 'component distribution', ylim = c(0,0.1))
  for(i in 1:num_components){
    lines(x = c(i,i), y = c(means$`5%`[i],means$`95%`[i]), col = 'red')
  }
}

dev.off()

write_csv(outputs, '../data_processed/mpnpshort_allwindows_networkcomplexity.csv')
#### long window ####
outputs <- data.frame(time_window = 1:max_time_window,
                      mean_edge = NA, sd_edge = NA,
                      median_edge = NA, q2.5_edge = NA, q97.5_edge = NA,
                      mean_strength = NA, sd_strength = NA,
                      median_strength = NA, q2.5_strength = NA, q97.5_strength = NA,
                      mean_eigen = NA, sd_eigen = NA,
                      median_eigen = NA, q2.5_eigen = NA, q97.5_eigen = NA,
                      mean_density = NA, sd_density = NA,
                      median_density = NA, q2.5_density = NA, q97.5_density = NA,
                      mean_cv = NA, sd_cv = NA,
                      median_cv = NA, q2.5_cv = NA, q97.5_cv = NA,
                      mean_std = NA, sd_std = NA,
                      median_std = NA, q2.5_std = NA, q97.5_std = NA,
                      num_components = NA,
                      prob_k1 = NA, prob_k2 = NA, prob_k3 = NA, prob_k4 = NA, prob_k5 = NA)

## load in network model
load(paste0('mpnp_edgecalculations/mpnplong',time_window,'_bisonr_edgescalculated.RData'))

#### network metrics
### edge-level properties
# weight
edge_weight <- extract_metric(mpnp_edge_weights, "edge_weight")
saveRDS(edge_weight, paste0('../data_processed/mpnp_network_properties/mpnplong',
                            time_window,'_edgeweights.RDS'))
outputs$mean_edge[time_window] <- mean(edge_weight)
outputs$sd_edge[time_window] <- sd(edge_weight)
outputs$median_edge[time_window] <- median(edge_weight)
outputs$q2.5_edge[time_window] <- quantile(edge_weight, 0.025)
outputs$q97.5_edge[time_window] <- quantile(edge_weight, 0.975)

### node-level properties
# strength
node_strength <- extract_metric(mpnp_edge_weights, "node_strength")
saveRDS(node_strength, paste0('../data_processed/mpnp_network_properties/mpnplong',
                              time_window,'_nodestrengths.RDS'))
outputs$mean_strength[time_window] <- mean(node_strength)
outputs$sd_strength[time_window] <- sd(node_strength)
outputs$median_strength[time_window] <- median(node_strength)
outputs$q2.5_strength[time_window] <- quantile(node_strength, 0.025)
outputs$q97.5_strength[time_window] <- quantile(node_strength, 0.975)

# eigenvector centrality
node_eigen <- extract_metric(mpnp_edge_weights, "node_eigen")
saveRDS(node_eigen, paste0('../data_processed/mpnp_network_properties/mpnplong',
                           time_window,'_nodeeigenvectors.RDS'))
node_strength <- extract_metric(mpnp_edge_weights, "node_strength")
saveRDS(node_strength, paste0('../data_processed/mpnp_network_properties/mpnplong',
                              time_window,'_nodestrengths.RDS'))
outputs$mean_eigen[time_window] <- mean(node_eigen)
outputs$sd_eigen[time_window] <- sd(node_eigen)
outputs$median_eigen[time_window] <- median(node_eigen)
outputs$q2.5_eigen[time_window] <- quantile(node_eigen, 0.025)
outputs$q97.5_eigen[time_window] <- quantile(node_eigen, 0.975)

### network-level properties
# weighted density
global_density <- extract_metric(mpnp_edge_weights, "global_density")
saveRDS(global_density, paste0('../data_processed/mpnp_network_properties/mpnplong',
                               time_window,'_globaldensity.RDS'))
outputs$mean_density[time_window] <- mean(global_density)
outputs$sd_density[time_window] <- sd(global_density)
outputs$median_density[time_window] <- median(global_density)
outputs$q2.5_density[time_window] <- quantile(global_density, 0.025)
outputs$q97.5_density[time_window] <- quantile(global_density, 0.975)

# coefficient of variation of edge weights (aka social differentation)
global_cv <- extract_metric(mpnp_edge_weights, "global_cv")
saveRDS(global_cv, paste0('../data_processed/mpnp_network_properties/mpnplong',
                          time_window,'_globalcv.RDS'))
outputs$mean_cv[time_window] <- mean(global_cv)
outputs$sd_cv[time_window] <- sd(global_cv)
outputs$median_cv[time_window] <- median(global_cv)
outputs$q2.5_cv[time_window] <- quantile(global_cv, 0.025)
outputs$q97.5_cv[time_window] <- quantile(global_cv, 0.975)

# standard deviation of edge weights
global_std <- extract_metric(mpnp_edge_weights, "global_std")
saveRDS(global_std, paste0('../data_processed/mpnp_network_properties/mpnplong',
                           time_window,'_globalstd.RDS'))
outputs$mean_std[time_window] <- mean(global_std)
outputs$sd_std[time_window] <- sd(global_std)
outputs$median_std[time_window] <- median(global_std)
outputs$q2.5_std[time_window] <- quantile(global_std, 0.025)
outputs$q97.5_std[time_window] <- quantile(global_std, 0.975)

#### network complexity
# Edge mixture models can be used to detect different types of social connections by looking for clustering in edge weights. There are two key outputs from a bisonR edge mixture model: a posterior probability distribution over the number of components (clusters, interpreted associal connection types, from 1 to a maximum 𝐾, over the entire network), the probability that each dyad belongs to each possible component.

# Once the model has fit well, we can use the bison_mixture() function to fit the edge mixture model:
motnp_mixture <- bison_mixture(mpnp_edge_weights, num_components = 5, verbose = FALSE)
summary(motnp_mixture)

# To get network-level probabilities of the number of components, we can use get_network_component_probabilities() function:
num_components <- get_network_component_probabilities(motnp_mixture)
outputs$num_components[time_window] <- num_components$num_components[which(num_components$probability == max(num_components$probability))]
outputs$prob_k1[time_window] <- num_components$probability[1]
outputs$prob_k2[time_window] <- num_components$probability[2]
outputs$prob_k3[time_window] <- num_components$probability[3]
outputs$prob_k4[time_window] <- num_components$probability[4]
outputs$prob_k5[time_window] <- num_components$probability[5]

# It’s often useful to know where the mixture components are located, so we can get their posterior means using the get_component_means() function for a given number of components. Here we’ll calculate the posterior means of the components for the mixture model with 3 components:
pdf('../outputs/mpnplongwindow_components.pdf')
num_components <- num_components$num_components[which(num_components$probability == max(num_components$probability))]
(means <- get_component_means(motnp_mixture, num_components))
plot(means$`50%`, las = 1, pch = '-', xlab = 'K', ylab = 'component distribution', ylim = c(0,0.1))
for(i in 1:num_components){
  lines(x = c(i,i), y = c(means$`5%`[i],means$`95%`[i]), col = 'red')
}

write_csv(outputs, '../data_processed/mpnplongwindow_networkcomplexity.csv')