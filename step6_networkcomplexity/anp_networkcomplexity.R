#### set up ####
# load packages
#library(tidyverse) ; library(dplyr) ; library(rstan) ; library(igraph) ; library(cmdstanr) ; library(bisonR) ; library(brms) ; library(mclust)
library(tidyverse, lib.loc = '../packages/') # library(tidyverse)
library(bisonR, lib.loc = '../packages/')    # library(bisonR)

#### short windows ####
max_time_window <- 36
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
pdf('../outputs/anpshort_allwindows_components.pdf')

for( time_window in 1:max_time_window){
  #### load in network model ####
  load(paste0('anp_edgecalculations/anpshort',time_window,'_bisonr_edgescalculated.RData'))
  
  #### network metrics ####
  ### edge-level properties
  # weight
  edge_weight <- extract_metric(anp_edge_weights, "edge_weight")
  saveRDS(edge_weight, paste0('../data_processed/anp_network_properties/anpshort',
                              time_window,'_edgeweights.RDS'))
  outputs$mean_edge[time_window] <- mean(edge_weight)
  outputs$sd_edge[time_window] <- sd(edge_weight)
  outputs$median_edge[time_window] <- median(edge_weight)
  outputs$q2.5_edge[time_window] <- quantile(edge_weight, 0.025)
  outputs$q97.5_edge[time_window] <- quantile(edge_weight, 0.975)
  
  ### node-level properties
  # strength
  node_strength <- extract_metric(anp_edge_weights, "node_strength")
  saveRDS(node_strength, paste0('../data_processed/anp_network_properties/anpshort',
                                time_window,'_nodestrengths.RDS'))
  outputs$mean_strength[time_window] <- mean(node_strength)
  outputs$sd_strength[time_window] <- sd(node_strength)
  outputs$median_strength[time_window] <- median(node_strength)
  outputs$q2.5_strength[time_window] <- quantile(node_strength, 0.025)
  outputs$q97.5_strength[time_window] <- quantile(node_strength, 0.975)
  
  # eigenvector centrality
  node_eigen <- extract_metric(anp_edge_weights, "node_eigen")
  saveRDS(node_eigen, paste0('../data_processed/anp_network_properties/anpshort',
                             time_window,'_nodeeigenvectors.RDS'))
  node_strength <- extract_metric(anp_edge_weights, "node_strength")
  saveRDS(node_strength, paste0('../data_processed/anp_network_properties/anpshort',
                                time_window,'_nodestrengths.RDS'))
  outputs$mean_eigen[time_window] <- mean(node_eigen)
  outputs$sd_eigen[time_window] <- sd(node_eigen)
  outputs$median_eigen[time_window] <- median(node_eigen)
  outputs$q2.5_eigen[time_window] <- quantile(node_eigen, 0.025)
  outputs$q97.5_eigen[time_window] <- quantile(node_eigen, 0.975)
  
  ### network-level properties
  # weighted density
  global_density <- extract_metric(anp_edge_weights, "global_density")
  saveRDS(global_density, paste0('../data_processed/anp_network_properties/anpshort',
                                 time_window,'_globaldensity.RDS'))
  outputs$mean_density[time_window] <- mean(global_density)
  outputs$sd_density[time_window] <- sd(global_density)
  outputs$median_density[time_window] <- median(global_density)
  outputs$q2.5_density[time_window] <- quantile(global_density, 0.025)
  outputs$q97.5_density[time_window] <- quantile(global_density, 0.975)
  
  # coefficient of variation of edge weights (aka social differentation)
  global_cv <- extract_metric(anp_edge_weights, "global_cv")
  saveRDS(global_cv, paste0('../data_processed/anp_network_properties/anpshort',
                            time_window,'_globalcv.RDS'))
  outputs$mean_cv[time_window] <- mean(global_cv)
  outputs$sd_cv[time_window] <- sd(global_cv)
  outputs$median_cv[time_window] <- median(global_cv)
  outputs$q2.5_cv[time_window] <- quantile(global_cv, 0.025)
  outputs$q97.5_cv[time_window] <- quantile(global_cv, 0.975)
  
  # standard deviation of edge weights
  global_std <- extract_metric(anp_edge_weights, "global_std")
  saveRDS(global_std, paste0('../data_processed/anp_network_properties/anpshort',
                             time_window,'_globalstd.RDS'))
  outputs$mean_std[time_window] <- mean(global_std)
  outputs$sd_std[time_window] <- sd(global_std)
  outputs$median_std[time_window] <- median(global_std)
  outputs$q2.5_std[time_window] <- quantile(global_std, 0.025)
  outputs$q97.5_std[time_window] <- quantile(global_std, 0.975)
  
  #### network complexity ####
  # Edge mixture models can be used to detect different types of social connections by looking for clustering in edge weights. There are two key outputs from a bisonR edge mixture model: a posterior probability distribution over the number of components (clusters, interpreted associal connection types, from 1 to a maximum 𝐾, over the entire network), the probability that each dyad belongs to each possible component.
  
  # Once the model has fit well, we can use the bison_mixture() function to fit the edge mixture model:
  anp_mixture <- bison_mixture(anp_edge_weights, num_components = 5, verbose = FALSE)
  summary(anp_mixture)
  
  # To get network-level probabilities of the number of components, we can use get_network_component_probabilities() function:
  num_components <- get_network_component_probabilities(anp_mixture)
  outputs$num_components[time_window] <- num_components$num_components[which(num_components$probability == max(num_components$probability))]
  outputs$prob_k1[time_window] <- num_components$probability[1]
  outputs$prob_k2[time_window] <- num_components$probability[2]
  outputs$prob_k3[time_window] <- num_components$probability[3]
  outputs$prob_k4[time_window] <- num_components$probability[4]
  outputs$prob_k5[time_window] <- num_components$probability[5]
  
  # It’s often useful to know where the mixture components are located, so we can get their posterior means using the get_component_means() function for a given number of components. Here we’ll calculate the posterior means of the components for the mixture model with 3 components:
  num_components <- num_components$num_components[which(num_components$probability == max(num_components$probability))]
  (means <- get_component_means(anp_mixture, num_components))
  plot(means$`50%`, las = 1, pch = '-', xlab = 'K', ylab = 'component distribution', ylim = c(0,0.1))
  for(i in 1:num_components){
    lines(x = c(i,i), y = c(means$`5%`[i],means$`95%`[i]), col = 'red')
  }
}
warnings()

dev.off()

write_csv(outputs, '../data_processed/anpshort_allwindows_networkcomplexity.csv')
#### long windows ####
max_time_window <- 7
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
pdf('../outputs/anplong_allwindows_components.pdf')

for( time_window in 1:max_time_window){
  #### load in network model ####
  load(paste0('anp_edgecalculations/anplong',time_window,'_bisonr_edgescalculated.RData'))
  
  #### network metrics ####
  ### edge-level properties
  # weight
  edge_weight <- extract_metric(anp_edge_weights, "edge_weight")
  saveRDS(edge_weight, paste0('../data_processed/anp_network_properties/anplong',
                              time_window,'_edgeweights.RDS'))
  outputs$mean_edge[time_window] <- mean(edge_weight)
  outputs$sd_edge[time_window] <- sd(edge_weight)
  outputs$median_edge[time_window] <- median(edge_weight)
  outputs$q2.5_edge[time_window] <- quantile(edge_weight, 0.025)
  outputs$q97.5_edge[time_window] <- quantile(edge_weight, 0.975)
  
  ### node-level properties
  # strength
  node_strength <- extract_metric(anp_edge_weights, "node_strength")
  saveRDS(node_strength, paste0('../data_processed/anp_network_properties/anplong',
                                time_window,'_nodestrengths.RDS'))
  outputs$mean_strength[time_window] <- mean(node_strength)
  outputs$sd_strength[time_window] <- sd(node_strength)
  outputs$median_strength[time_window] <- median(node_strength)
  outputs$q2.5_strength[time_window] <- quantile(node_strength, 0.025)
  outputs$q97.5_strength[time_window] <- quantile(node_strength, 0.975)
  
  # eigenvector centrality
  node_eigen <- extract_metric(anp_edge_weights, "node_eigen")
  saveRDS(node_eigen, paste0('../data_processed/anp_network_properties/anplong',
                             time_window,'_nodeeigenvectors.RDS'))
  node_strength <- extract_metric(anp_edge_weights, "node_strength")
  saveRDS(node_strength, paste0('../data_processed/anp_network_properties/anplong',
                                time_window,'_nodestrengths.RDS'))
  outputs$mean_eigen[time_window] <- mean(node_eigen)
  outputs$sd_eigen[time_window] <- sd(node_eigen)
  outputs$median_eigen[time_window] <- median(node_eigen)
  outputs$q2.5_eigen[time_window] <- quantile(node_eigen, 0.025)
  outputs$q97.5_eigen[time_window] <- quantile(node_eigen, 0.975)
  
  ### network-level properties
  # weighted density
  global_density <- extract_metric(anp_edge_weights, "global_density")
  saveRDS(global_density, paste0('../data_processed/anp_network_properties/anplong',
                                 time_window,'_globaldensity.RDS'))
  outputs$mean_density[time_window] <- mean(global_density)
  outputs$sd_density[time_window] <- sd(global_density)
  outputs$median_density[time_window] <- median(global_density)
  outputs$q2.5_density[time_window] <- quantile(global_density, 0.025)
  outputs$q97.5_density[time_window] <- quantile(global_density, 0.975)
  
  # coefficient of variation of edge weights (aka social differentation)
  global_cv <- extract_metric(anp_edge_weights, "global_cv")
  saveRDS(global_cv, paste0('../data_processed/anp_network_properties/anplong',
                            time_window,'_globalcv.RDS'))
  outputs$mean_cv[time_window] <- mean(global_cv)
  outputs$sd_cv[time_window] <- sd(global_cv)
  outputs$median_cv[time_window] <- median(global_cv)
  outputs$q2.5_cv[time_window] <- quantile(global_cv, 0.025)
  outputs$q97.5_cv[time_window] <- quantile(global_cv, 0.975)
  
  # standard deviation of edge weights
  global_std <- extract_metric(anp_edge_weights, "global_std")
  saveRDS(global_std, paste0('../data_processed/anp_network_properties/anplong',
                             time_window,'_globalstd.RDS'))
  outputs$mean_std[time_window] <- mean(global_std)
  outputs$sd_std[time_window] <- sd(global_std)
  outputs$median_std[time_window] <- median(global_std)
  outputs$q2.5_std[time_window] <- quantile(global_std, 0.025)
  outputs$q97.5_std[time_window] <- quantile(global_std, 0.975)
  
  #### network complexity ####
  # Edge mixture models can be used to detect different types of social connections by looking for clustering in edge weights. There are two key outputs from a bisonR edge mixture model: a posterior probability distribution over the number of components (clusters, interpreted associal connection types, from 1 to a maximum 𝐾, over the entire network), the probability that each dyad belongs to each possible component.
  
  # Once the model has fit well, we can use the bison_mixture() function to fit the edge mixture model:
  anp_mixture <- bison_mixture(anp_edge_weights, num_components = 5, verbose = FALSE)
  summary(anp_mixture)
  
  # To get network-level probabilities of the number of components, we can use get_network_component_probabilities() function:
  num_components <- get_network_component_probabilities(anp_mixture)
  outputs$num_components[time_window] <- num_components$num_components[which(num_components$probability == max(num_components$probability))]
  outputs$prob_k1[time_window] <- num_components$probability[1]
  outputs$prob_k2[time_window] <- num_components$probability[2]
  outputs$prob_k3[time_window] <- num_components$probability[3]
  outputs$prob_k4[time_window] <- num_components$probability[4]
  outputs$prob_k5[time_window] <- num_components$probability[5]
  
  # It’s often useful to know where the mixture components are located, so we can get their posterior means using the get_component_means() function for a given number of components. Here we’ll calculate the posterior means of the components for the mixture model with 3 components:
  num_components <- num_components$num_components[which(num_components$probability == max(num_components$probability))]
  (means <- get_component_means(anp_mixture, num_components))
  plot(means$`50%`, las = 1, pch = '-', xlab = 'K', ylab = 'component distribution', ylim = c(0,0.1))
  for(i in 1:num_components){
    lines(x = c(i,i), y = c(means$`5%`[i],means$`95%`[i]), col = 'red')
  }
}

warnings()

dev.off()

write_csv(outputs, '../data_processed/anplong_allwindows_networkcomplexity.csv')

#### check outputs ####
### short windows
short_windows <- read_csv('../data_processed/anpshort_allwindows_networkcomplexity.csv')
max_time_window <- nrow(short_windows)
hist(short_windows$mean_edge, breaks = max_time_window)
hist(short_windows$mean_strength, breaks = max_time_window)
hist(short_windows$mean_eigen, breaks = max_time_window)
hist(short_windows$mean_density, breaks = max_time_window)
hist(short_windows$mean_cv, breaks = max_time_window)
hist(short_windows$mean_std, breaks = max_time_window)

plot(plogis(short_windows$mean_edge) ~ short_windows$time_window, type = 'l',
     ylab = 'inv. logit mean edges', xlab = 'time window', las = 1, col = 'red',
     ylim = c(0,0.4))
lines(plogis(short_windows$q2.5_edge) ~ short_windows$time_window, col = 'red', lty = 2)
lines(plogis(short_windows$q97.5_edge) ~ short_windows$time_window, col = 'red', lty = 2)

plot(short_windows$mean_strength ~ short_windows$time_window, type = 'l',
     ylab = 'mean strength centrality', xlab = 'time window', las = 1, col = 'orange',
     ylim = c(0,15))
lines(short_windows$q2.5_strength ~ short_windows$time_window, col = 'orange', lty = 2)
lines(short_windows$q97.5_strength ~ short_windows$time_window, col = 'orange', lty = 2)

plot(short_windows$mean_eigen ~ short_windows$time_window, type = 'l', ylim = c(0.3,1),
     ylab = 'mean eigenvector centrality', xlab = 'time window', las = 1, col = 'green')
lines(short_windows$q2.5_eigen ~ short_windows$time_window, col = 'green', lty = 2)
lines(short_windows$q97.5_eigen ~ short_windows$time_window, col = 'green', lty = 2)

plot(short_windows$mean_density ~ short_windows$time_window, type = 'l', ylim = c(0,0.1),
     ylab = 'mean global density', xlab = 'time window', las = 1, col = 'blue')
lines(short_windows$q2.5_density ~ short_windows$time_window, col = 'blue', lty = 2)
lines(short_windows$q97.5_density ~ short_windows$time_window, col = 'blue', lty = 2)

plot(short_windows$mean_cv ~ short_windows$time_window, type = 'l', ylim = c(0.9, 1.3),
     ylab = 'mean global cv', xlab = 'time window', las = 1, col = 'purple2')
lines(short_windows$q2.5_cv ~ short_windows$time_window, col = 'purple2', lty = 2)
lines(short_windows$q97.5_cv ~ short_windows$time_window, col = 'purple2', lty = 2)

plot(short_windows$mean_std ~ short_windows$time_window, type = 'l', ylim = c(0.03,0.13),
     ylab = 'mean global std', xlab = 'time window', las = 1, col = 'magenta')
lines(short_windows$q2.5_std ~ short_windows$time_window, col = 'magenta', lty = 2)
lines(short_windows$q97.5_std ~ short_windows$time_window, col = 'magenta', lty = 2)

table(short_windows$num_components)
k_probs <- short_windows[,c('prob_k1','prob_k2','prob_k3','prob_k4','prob_k5')] %>%
  pivot_longer(everything(), values_to = 'probability', names_to = 'component')
table(k_probs$component, k_probs$probability)
k_probs$probability2 <- ifelse(k_probs$probability < 0.00000000001, 0, k_probs$probability)
table(k_probs$component, k_probs$probability2)

### long windows
long_windows <- read_csv('../data_processed/anplong_allwindows_networkcomplexity.csv')
max_time_window <- nrow(long_windows)
hist(long_windows$mean_edge, breaks = max_time_window)
hist(long_windows$mean_strength, breaks = max_time_window)
hist(long_windows$mean_eigen, breaks = max_time_window)
hist(long_windows$mean_density, breaks = max_time_window)
hist(long_windows$mean_cv, breaks = max_time_window)
hist(long_windows$mean_std, breaks = max_time_window)

plot(plogis(long_windows$mean_edge) ~ long_windows$time_window, type = 'l',
     ylab = 'inv. logit mean edges', xlab = 'time window', las = 1, col = 'red',
     ylim = c(0,0.2))
lines(plogis(long_windows$q2.5_edge) ~ long_windows$time_window, col = 'red', lty = 2)
lines(plogis(long_windows$q97.5_edge) ~ long_windows$time_window, col = 'red', lty = 2)

plot(long_windows$mean_strength ~ long_windows$time_window, type = 'l',
     ylab = 'mean strength centrality', xlab = 'time window', las = 1, col = 'orange',
     ylim = c(0,15))
lines(long_windows$q2.5_strength ~ long_windows$time_window, col = 'orange', lty = 2)
lines(long_windows$q97.5_strength ~ long_windows$time_window, col = 'orange', lty = 2)

plot(long_windows$mean_eigen ~ long_windows$time_window, type = 'l', ylim = c(0,1),
     ylab = 'mean eigenvector centrality', xlab = 'time window', las = 1, col = 'green')
lines(long_windows$q2.5_eigen ~ long_windows$time_window, col = 'green', lty = 2)
lines(long_windows$q97.5_eigen ~ long_windows$time_window, col = 'green', lty = 2)

plot(long_windows$mean_density ~ long_windows$time_window, type = 'l', ylim = c(0.02,0.04),
     ylab = 'mean global density', xlab = 'time window', las = 1, col = 'blue')
lines(long_windows$q2.5_density ~ long_windows$time_window, col = 'blue', lty = 2)
lines(long_windows$q97.5_density ~ long_windows$time_window, col = 'blue', lty = 2)

plot(long_windows$mean_cv ~ long_windows$time_window, type = 'l', ylim = c(0.9, 1.5),
     ylab = 'mean global cv', xlab = 'time window', las = 1, col = 'purple2')
lines(long_windows$q2.5_cv ~ long_windows$time_window, col = 'purple2', lty = 2)
lines(long_windows$q97.5_cv ~ long_windows$time_window, col = 'purple2', lty = 2)

plot(long_windows$mean_std ~ long_windows$time_window, type = 'l', ylim = c(0,0.06),
     ylab = 'mean global std', xlab = 'time window', las = 1, col = 'magenta')
lines(long_windows$q2.5_std ~ long_windows$time_window, col = 'magenta', lty = 2)
lines(long_windows$q97.5_std ~ long_windows$time_window, col = 'magenta', lty = 2)

table(long_windows$num_components)
k_probs <- long_windows[,c('prob_k1','prob_k2','prob_k3','prob_k4','prob_k5')] %>%
  pivot_longer(everything(), values_to = 'probability', names_to = 'component')
table(k_probs$component, k_probs$probability)
