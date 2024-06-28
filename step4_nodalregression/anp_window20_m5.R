#### set up #####
## load packages
library(LaplacesDemon) ; library(MASS) ; library(tidyverse) ; library(cmdstanr) ; library(sna)
set_cmdstan_path('../packages/.cmdstan/cmdstan-2.34.1/') # Viking ; set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0/') # Desktop

## set theme for plots
theme_set(theme_bw())

## create function to extract centrality
extract_eigen_centrality <- function(nodes_df, dyads_df, edgeweight_matrix, logit = TRUE, window){
  ## calculate data size parameters
  num_nodes <- nrow(nodes_df)
  num_dyads <- nrow(dyads_df)
  num_samples <- nrow(edgeweight_matrix)
  
  ## build adjacency tensor
  dyads_df$node_1_id <- as.integer(as.factor(dyads_df$node_1_randomised))
  dyads_df$node_2_id <- as.integer(as.factor(dyads_df$node_2_randomised))+1
  adj_tensor <- array(0, c(num_samples, num_nodes, num_nodes),
                      dimnames = list(NULL, nodes_df$node_window, nodes_df$node_window))
  
  ## fill adjacency tensor
  for (dyad_id in 1:num_dyads) {
    dyad_row <- dyads_df[dyad_id, ]
    adj_tensor[, dyad_row$node_1_id, dyad_row$node_2_id] <- edgeweight_matrix[, dyad_id]
    adj_tensor[, dyad_row$node_2_id, dyad_row$node_1_id] <- edgeweight_matrix[, dyad_id]
  }
  
  ## calculate centrality and store posterior samples in a matrix
  centrality_samples_invlogit <- matrix(0, num_samples, num_nodes,
                                        dimnames = list(NULL, nodes_df$node_window))
  for(i in 1:(num_samples)){
    centrality_samples_invlogit[i, ] <- sna::evcent(adj_tensor[i,,], gmode = 'graph')
  }
  
  ## convert to logit scale
  if(logit == TRUE) {
    centrality_samples <- logit(centrality_samples_invlogit)
    return(centrality_samples)
  } else {
    return(centrality_samples_invlogit)
  }
}

## create function to plot traceplots
traceplot <- function(parameter_df, parameters_to_plot, all_chains = TRUE, plot_chain = 1){
  parameters <- parameter_df %>%
    dplyr::select(all_of(parameters_to_plot),`.draw`,`.chain`,`.iteration`) %>%
    pivot_longer(cols = all_of(parameters_to_plot), names_to = 'parameter', values_to = 'draw') %>%
    rename(chain = .chain,
           chain_position = .iteration,
           draw_id = .draw)
  if(all_chains == FALSE){
    parameters <- parameters %>%
      filter(chain == plot_chain) # inspect individual chains -- check for wandery sections that might be hidden by other chains when all plotted together
  }
  ggplot(data = parameters,
         aes(x = chain_position, y = draw, colour = as.factor(chain)))+
    geom_line()+
    facet_wrap(. ~ parameter, scales = 'free_y')+
    theme_bw()+
    theme(legend.position = 'none')
}

## create posterior predictive check function
post_pred_check <- function(centrality_matrix, nodes_df, time_window, cent_cov, parameters){
  ## set up plot
  plot(density(centrality_matrix[1, which(nodes_df$window == time_window)]),
       las = 1, ylim = c(0,1.5),
       #main = paste0("Posterior check:\nblack = data, blue = predicted, window = ", time_window),
       main = paste0("Posterior check: window = ", time_window),
       col=rgb(0, 0, 0, 0.25))
  
  ## create data required
  eigen_data <- nodes_df %>%
    dplyr::select(node_random, age, window) %>%
    filter(window == time_window)
  num_nodes_window <- nrow(eigen_data)
  
  ## plot lines of predictions
  for (i in 1:100) {
    j <- sample(1:length(parameters$beta_age), 1)
    lines(density(centrality_matrix[j, which(nodes_df$window == time_window)]),
          col = rgb(0, 0, 0, 0.25))
    mu <- parameters$beta_age[j]*eigen_data$age + parameters$intercept[j]
    for(k in 1:length(mu)) {
      mu[k] <- mu[k] + as.numeric(rand_window[j,time_window]) + as.numeric(rand_node[j,eigen_data$node_random[k]])
    }
    sigma <- cent_cov + diag(rep(parameters$sigma[j], num_nodes_window))
    #nu <- params$nu[j]
    lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
    #lines(density(LaplacesDemon::rmvt(n = 1, mu = mu, S = sigma, df = nu)), col=rgb(0, 0, 1, 0.25))
  }
}

## create mean prediction function
get_mean_predictions <- function(predict_df, parameters, include_node = TRUE, include_window = TRUE){
  ## create empty matrix to fill with predictions
  mean_matrix <- matrix(NA, nrow = nrow(parameters), ncol = nrow(predict_df),
                        dimnames = list(NULL, predict_df$node_window))
  
  ## populate matrix = mean centrality values per node, predicting for real data
  for(i in 1:nrow(mean_matrix)){
    mean_matrix[i,] <- parameters$beta_age[i] * predict_df$age_std + parameters$intercept[i]
  }
  
  if(include_window == TRUE){
    window_effect <- parameters %>% dplyr::select(grep('window_random_effect', colnames(parameters), value=TRUE))
    for(i in 1:nrow(mean_matrix)){
      mean_matrix[i,] <- mean_matrix[i,] + as.numeric(window_effect[i,predict_df$window])
    }
  }
  
  if(include_node == TRUE){
    node_effect <- parameters %>% dplyr::select(grep('node_random_effect', colnames(parameters), value=TRUE))
    for(i in 1:nrow(mean_matrix)){
      mean_matrix[i,] <- mean_matrix[i,] + as.numeric(node_effect[i, predict_df$node_random])
    }
  }
  
  return(mean_matrix)
  
}

#### load window 20 ####
## import workspace image
load('anp_edgecalculations/anpshort20_edgeweights_conditionalprior.RData')

## add extra column for randomised node (though here it is not random!)
nodes <- nodes %>%
  mutate(node_1 = node, node_2 = node,
         node_1_randomised = node,
         node_2_randomised = node)

## randomise node IDs in dyads data frame
cdf <- cdf %>%
  filter(id_1 %in% nodes$id) %>%
  filter(id_2 %in% nodes$id) %>%
  left_join(nodes[,c('node_1','node_1_randomised')], by = 'node_1') %>%
  left_join(nodes[,c('node_2','node_2_randomised')], by = 'node_2')

## clear workspace a bit
rm(counts_df, edge_binary, edgelist, edges, make_edgelist, plot_network_threshold_anp) ; gc()

######## CHECK WITH ALL ELEPHANTS INCLUDED ####
#### check M5 data -- observed 14 times; 25 partners seen with once, 6 twice and 2 with M5 three times. Should be 33 other elephants that it uses the upper prior for and can't have a draw of zero. ####
m5 <- cdf %>% 
  filter(id_1 == 'M005' | id_2 == 'M005')
table(m5$event_count)
m5$period_count_1[1]

m5_partners <- m5 %>% 
  filter(event_count > 0)

m5_edge_ids <- which(cdf$dyad_id %in% m5$dyad_id)
m5_edges <- edge_samples[,m5_edge_ids]#as.character(m5$dyad_id)]
other_edges <- edge_samples[,! (1:ncol(edge_samples) %in% m5_edge_ids)]

m5_means <- apply(m5_edges, 2, mean)
other_means <- apply(other_edges, 2, mean)
mean(m5_means)
mean(other_means)

m5_means <- as.data.frame(m5_means) %>% 
  rename(mean_edge = m5_means)
m5_means$dyad_id <- rownames(m5_means)
m5_means$m5 <- 'm5'

other_means <- as.data.frame(other_means) %>% 
  rename(mean_edge = other_means)
other_means$dyad_id <- rownames(other_means)
other_means$m5 <- 'other'

mean_edges <- rbind(m5_means, other_means)
ggplot(mean_edges)+
  geom_jitter(aes(x = m5, y = mean_edge))+
  geom_boxplot(aes(x = m5, y = mean_edge),
               notch = T,
               fill = 'seagreen1',
               colour = 'red')

m5_summary <- summary %>% 
  filter(dyad_id %in% m5$dyad_id) %>% 
  mutate(together = ifelse(event_count == 0, 'never together', 'sometimes together'))
ggplot(m5_summary)+
  geom_boxplot(aes(x = as.factor(event_count), y = median, fill = together),
               notch = T,
               alpha = 0.4)+
  geom_jitter(aes(x = as.factor(event_count), y = median, colour = together))+
  scale_fill_viridis_d(end = 0.5)+
  scale_colour_viridis_d(end = 0.5)

summary <- summary %>% 
  mutate(m5 = ifelse(dyad_id %in% m5$dyad_id,
                    'm5','other'),
         together = ifelse(event_count == 0, 'never together', 'sometimes together'))
ggplot(summary)+
  geom_jitter(aes(x = as.factor(event_count),
                  y = median,
                  colour = m5,
                  group = m5))+
  geom_boxplot(aes(x = as.factor(event_count),
                   y = median,
                   fill = m5),
               notch = T,
               alpha = 0.4)+
  scale_fill_viridis_d(end = 0.5)+
  scale_colour_viridis_d(end = 0.5)

sightings14 <- summary %>% 
  filter(period_count_1 %in% c(13,14,15) | 
           period_count_2 %in% c(13,14,15))

ggplot(sightings14)+
  # geom_jitter(aes(x = as.factor(event_count),
  #                 y = median,
  #                 colour = m5,
  #                 group = m5))+
  geom_boxplot(aes(x = as.factor(event_count),
                   y = median,
                   fill = m5),
               notch = T,
               alpha = 0.4)+
  scale_fill_viridis_d(end = 0.5)+
  scale_colour_viridis_d(end = 0.5)

m5_dyads_to_plot <- sample(1:ncol(m5_edges), 16, replace = F)
other_dyads_to_plot <- sample(1:ncol(other_edges), 16, replace = F)
par(mfrow = c(4,4),
    mai = c(0.1,0.1,0.1,0.1))
for(i in 1:16){
  plot(m5_edges[1:1000,m5_dyads_to_plot[i]], type = 'l', col = 'red')
  lines(m5_edges[1001:2000,m5_dyads_to_plot[i]], col = 'green')
  lines(m5_edges[2001:3000,m5_dyads_to_plot[i]], col = 'blue')
  lines(m5_edges[3001:4000,m5_dyads_to_plot[i]], col = 'purple')
}
for(i in 1:16){
  plot(other_edges[1:1000,other_dyads_to_plot[i]], type = 'l', col = 'red')
  lines(other_edges[1001:2000,other_dyads_to_plot[i]], col = 'green')
  lines(other_edges[2001:3000,other_dyads_to_plot[i]], col = 'blue')
  lines(other_edges[3001:4000,other_dyads_to_plot[i]], col = 'purple')
}

#### extract centrality -- all as one function ####
cent_new <- extract_eigen_centrality(nodes_df = nodes,
                                     dyads_df = cdf,
                                     edgeweight_matrix = edge_samples,
                                     logit = TRUE)

#### extract centrality -- stepwise, as in function, to see if I can find the issue ####
## calculate data size parameters
num_nodes <- nrow(nodes)
num_dyads <- nrow(cdf)
num_samples <- nrow(edge_samples)

## build adjacency tensor
cdf$node_1_id <- as.integer(as.factor(cdf$node_1_randomised))
cdf$node_2_id <- as.integer(as.factor(cdf$node_2_randomised))+1
adj_tensor <- array(0, c(num_samples, num_nodes, num_nodes),
                    dimnames = list(NULL, nodes$node, nodes$node))

## fill adjacency tensor
for (dyad_id in 1:num_dyads) {
  dyad_row <- cdf[dyad_id, ]
  adj_tensor[, dyad_row$node_1_id, dyad_row$node_2_id] <- edge_samples[, dyad_id]
  adj_tensor[, dyad_row$node_2_id, dyad_row$node_1_id] <- edge_samples[, dyad_id]
}
adj_tensor[which(nodes$id == 'M005'),,]

## calculate centrality and store posterior samples in a matrix
centrality_samples_invlogit <- matrix(0, num_samples, num_nodes,
                                      dimnames = list(NULL, nodes$node))
for(i in 1:(num_samples)){
  centrality_samples_invlogit[i, ] <- sna::evcent(adj_tensor[i,,], gmode = 'graph')
}

## convert to logit scale
centrality_samples <- logit(centrality_samples_invlogit)

#### check ####
par(mfrow = c(1,1),
    mai = c(1,1,1,0.5))

## extract mean and covariance
steps_mu <- apply(centrality_samples, 2, mean)
steps_cov <- cov(centrality_samples)
function_mu <- apply(cent_new, 2, mean)
function_cov <- cov(cent_new)

## add mean estimate per ID to nodes data frame
steps_df <- as.data.frame(steps_mu) %>%
  rename(mean_eigen_steps = steps_mu)
steps_df$node <- as.numeric(rownames(steps_df))

function_df <- as.data.frame(function_mu) %>%
  rename(mean_eigen_function = function_mu)
function_df$node <- as.numeric(rownames(function_df))

nodes <- nodes %>%
  dplyr::select(-node_1, -node_2, -node_1_randomised, -node_2_randomised) %>%
  left_join(steps_df, by = 'node') %>%
  left_join(function_df, by = 'node')

## see where it falls -- identical
hist(nodes$mean_eigen_steps)
hist(nodes$mean_eigen_function)
plot(nodes$mean_eigen_function ~ nodes$sightings, pch = 19, col = rgb(0,0,0,0.4))
points(mean_eigen_function ~ sightings,
       data = nodes[nodes$id == 'M005',],
       pch = 19, col = 'red', cex = 2)

rm(list = ls()[! ls() %in% c('cdf','edge_samples','extract_eigen_centrality','fit_edges_anp','get_mean_predictions','n_chains','n_dyads','n_samples','n_windows','nodes','periods','post_pred_check','summary','time_window','traceplot')])

######## REMOVE ALL DYADS WHERE ONE OR OTHER IS <10 ####
nodes <- nodes %>% 
  filter(age > 9)
cdf <- cdf %>% 
  filter(age_start_1 > 9) %>% 
  filter(age_start_2 > 9)
edge_samples <- edge_samples[,as.character(cdf$dyad_id)]

#### check M5 data -- observed 14 times; 25 partners seen with once, 6 twice and 2 with M5 three times. Should be 33 other elephants that it uses the upper prior for and can't have a draw of zero. ####
m5 <- cdf %>% 
  filter(id_1 == 'M005' | id_2 == 'M005')
table(m5$event_count)
m5$period_count_1[1]

m5_partners <- m5 %>% 
  filter(event_count > 0)

m5_edge_ids <- which(cdf$dyad_id %in% m5$dyad_id)
m5_edges <- edge_samples[,m5_edge_ids]#as.character(m5$dyad_id)]
other_edges <- edge_samples[,! (1:ncol(edge_samples) %in% m5_edge_ids)]

m5_means <- apply(m5_edges, 2, mean)
other_means <- apply(other_edges, 2, mean)
mean(m5_means)
mean(other_means)

m5_means <- as.data.frame(m5_means) %>% 
  rename(mean_edge = m5_means)
m5_means$dyad_id <- rownames(m5_means)
m5_means$m5 <- 'm5'

other_means <- as.data.frame(other_means) %>% 
  rename(mean_edge = other_means)
other_means$dyad_id <- rownames(other_means)
other_means$m5 <- 'other'

mean_edges <- rbind(m5_means, other_means)
ggplot(mean_edges)+
  geom_jitter(aes(x = m5, y = mean_edge))+
  geom_boxplot(aes(x = m5, y = mean_edge),
               notch = T,
               fill = 'seagreen1',
               colour = 'red')

m5_summary <- summary %>% 
  filter(dyad_id %in% m5$dyad_id) %>% 
  mutate(together = ifelse(event_count == 0, 'never together', 'sometimes together'))
ggplot(m5_summary)+
  geom_boxplot(aes(x = as.factor(event_count), y = median, fill = together),
               notch = T,
               alpha = 0.4)+
  geom_jitter(aes(x = as.factor(event_count), y = median, colour = together))+
  scale_fill_viridis_d(end = 0.5)+
  scale_colour_viridis_d(end = 0.5)

summary <- summary %>% 
  mutate(m5 = ifelse(dyad_id %in% m5$dyad_id,
                     'm5','other'),
         together = ifelse(event_count == 0, 'never together', 'sometimes together'))
ggplot(summary)+
  geom_jitter(aes(x = as.factor(event_count),
                  y = median,
                  colour = m5,
                  group = m5))+
  geom_boxplot(aes(x = as.factor(event_count),
                   y = median,
                   fill = m5),
               notch = T,
               alpha = 0.4)+
  scale_fill_viridis_d(end = 0.5)+
  scale_colour_viridis_d(end = 0.5)

sightings14 <- summary %>% 
  filter(period_count_1 %in% c(13,14,15) | 
           period_count_2 %in% c(13,14,15))

ggplot(sightings14)+
  # geom_jitter(aes(x = as.factor(event_count),
  #                 y = median,
  #                 colour = m5,
  #                 group = m5))+
  geom_boxplot(aes(x = as.factor(event_count),
                   y = median,
                   fill = m5),
               notch = T,
               alpha = 0.4)+
  scale_fill_viridis_d(end = 0.5)+
  scale_colour_viridis_d(end = 0.5)

m5_dyads_to_plot <- sample(1:ncol(m5_edges), 16, replace = F)
other_dyads_to_plot <- sample(1:ncol(other_edges), 16, replace = F)
par(mfrow = c(4,4),
    mai = c(0.1,0.1,0.1,0.1))
for(i in 1:16){
  plot(m5_edges[1:1000,m5_dyads_to_plot[i]], type = 'l', col = 'red')
  lines(m5_edges[1001:2000,m5_dyads_to_plot[i]], col = 'green')
  lines(m5_edges[2001:3000,m5_dyads_to_plot[i]], col = 'blue')
  lines(m5_edges[3001:4000,m5_dyads_to_plot[i]], col = 'purple')
}
for(i in 1:16){
  plot(other_edges[1:1000,other_dyads_to_plot[i]], type = 'l', col = 'red')
  lines(other_edges[1001:2000,other_dyads_to_plot[i]], col = 'green')
  lines(other_edges[2001:3000,other_dyads_to_plot[i]], col = 'blue')
  lines(other_edges[3001:4000,other_dyads_to_plot[i]], col = 'purple')
}

#### extract centrality -- all as one function ####
cent_new <- extract_eigen_centrality(nodes_df = nodes,
                                     dyads_df = cdf,
                                     edgeweight_matrix = edge_samples,
                                     logit = TRUE)

#### check ####
par(mfrow = c(1,1),
    mai = c(1,1,1,0.5))

## extract mean and covariance
eigen_mu <- apply(cent_new, 2, mean)
eigen_cov <- cov(cent_new)

## add mean estimate per ID to nodes data frame
eigen_df <- as.data.frame(eigen_mu) %>%
  rename(mean_eigen_eigen = eigen_mu)
eigen_df$node <- as.numeric(rownames(eigen_df))

nodes <- nodes %>%
  dplyr::select(-node_1, -node_2, -node_1_randomised, -node_2_randomised) %>%
  left_join(eigen_df, by = 'node')

## see where it falls -- identical
hist(nodes$mean_eigen_eigen)
plot(nodes$mean_eigen_eigen ~ nodes$sightings)
points(mean_eigen_eigen ~ sightings,
       data = nodes[nodes$id == 'M005',],
       pch = 19, col = 'red', cex = 2)

###################################
#### load window 3 ####
## import workspace image
load('anp_edgecalculations/anpshort3_edgeweights_conditionalprior.RData')

## add extra column for randomised node (though here it is not random!)
nodes <- nodes %>%
  mutate(node_1 = node, node_2 = node,
         node_1_randomised = node,
         node_2_randomised = node)

## randomise node IDs in dyads data frame
cdf <- cdf %>%
  filter(id_1 %in% nodes$id) %>%
  filter(id_2 %in% nodes$id) %>%
  left_join(nodes[,c('node_1','node_1_randomised')], by = 'node_1') %>%
  left_join(nodes[,c('node_2','node_2_randomised')], by = 'node_2')

## clear workspace a bit
rm(counts_df, edge_binary, edgelist, edges, make_edgelist, plot_network_threshold_anp) ; gc()

######## CHECK WITH ALL ELEPHANTS INCLUDED ####
#### check M108 data -- observed 1 time; 5 partners seen with once ####
m108 <- cdf %>% 
  filter(id_1 == 'M108' | id_2 == 'M108')
table(m108$event_count)
m108$period_count_1[1]

m108_partners <- m108 %>% 
  filter(event_count > 0)

m108_edge_ids <- which(cdf$dyad_id %in% m108$dyad_id)
m108_edges <- edge_samples[,m108_edge_ids]#as.character(m108$dyad_id)]
other_edges <- edge_samples[,! (1:ncol(edge_samples) %in% m108_edge_ids)]

m108_means <- apply(m108_edges, 2, mean)
other_means <- apply(other_edges, 2, mean)
mean(m108_means)
mean(other_means)

m108_means <- as.data.frame(m108_means)
colnames(m108_means) <- 'mean_edge'
m108_means$dyad_id <- rownames(m108_means)
m108_means$m108 <- 'm108'

other_means <- as.data.frame(other_means) %>% 
  rename(mean_edge = other_means)
other_means$dyad_id <- rownames(other_means)
other_means$m108 <- 'other'

mean_edges <- rbind(m108_means, other_means)
ggplot(mean_edges)+
  geom_jitter(aes(x = m108, y = mean_edge))+
  geom_boxplot(aes(x = m108, y = mean_edge),
               notch = T,
               fill = 'seagreen1',
               colour = 'red')

m108_summary <- summary %>% 
  filter(dyad_id %in% m108$dyad_id) %>% 
  mutate(together = ifelse(event_count == 0, 'never together', 'sometimes together'))
ggplot(m108_summary)+
  geom_boxplot(aes(x = as.factor(event_count), y = median, fill = together),
               notch = T,
               alpha = 0.4)+
  geom_jitter(aes(x = as.factor(event_count), y = median, colour = together))+
  scale_fill_viridis_d(end = 0.5)+
  scale_colour_viridis_d(end = 0.5)

summary <- summary %>% 
  mutate(m108 = ifelse(dyad_id %in% m108$dyad_id,
                     'm108','other'),
         together = ifelse(event_count == 0, 'never together', 'sometimes together'))
ggplot(summary)+
  geom_jitter(aes(x = as.factor(event_count),
                  y = median,
                  colour = m108,
                  group = m108))+
  geom_boxplot(aes(x = as.factor(event_count),
                   y = median,
                   fill = m108),
               notch = T,
               alpha = 0.4)+
  scale_fill_viridis_d(end = 0.5)+
  scale_colour_viridis_d(end = 0.5)

summary$period_count_1[summary$id_1 == 'M108']
sightings1 <- summary %>% 
  filter(period_count_1 %in% c(1,2) | 
           period_count_2 %in% c(1,2))

ggplot(sightings1)+
  # geom_jitter(aes(x = as.factor(event_count),
  #                 y = median,
  #                 colour = m108,
  #                 group = m108))+
  geom_boxplot(aes(x = as.factor(event_count),
                   y = median,
                   fill = m108),
               notch = T,
               alpha = 0.4)+
  scale_fill_viridis_d(end = 0.5)+
  scale_colour_viridis_d(end = 0.5)

m108_dyads_to_plot <- sample(1:ncol(m108_edges), 16, replace = F)
other_dyads_to_plot <- sample(1:ncol(other_edges), 16, replace = F)
par(mfrow = c(4,4),
    mai = c(0.1,0.1,0.1,0.1))
for(i in 1:16){
  plot(m108_edges[1:1000,m108_dyads_to_plot[i]], type = 'l', col = 'red')
  lines(m108_edges[1001:2000,m108_dyads_to_plot[i]], col = 'green')
  lines(m108_edges[2001:3000,m108_dyads_to_plot[i]], col = 'blue')
  lines(m108_edges[3001:4000,m108_dyads_to_plot[i]], col = 'purple')
}
for(i in 1:16){
  plot(other_edges[1:1000,other_dyads_to_plot[i]], type = 'l', col = 'red')
  lines(other_edges[1001:2000,other_dyads_to_plot[i]], col = 'green')
  lines(other_edges[2001:3000,other_dyads_to_plot[i]], col = 'blue')
  lines(other_edges[3001:4000,other_dyads_to_plot[i]], col = 'purple')
}

#### extract centrality ####
cent_new <- extract_eigen_centrality(nodes_df = nodes,
                                     dyads_df = cdf,
                                     edgeweight_matrix = edge_samples,
                                     logit = TRUE)

#### check ####
par(mfrow = c(1,1),
    mai = c(1,1,1,0.5))

## extract mean and covariance
function_mu <- apply(cent_new, 2, mean)
function_cov <- cov(cent_new)

## add mean estimate per ID to nodes data frame
function_df <- as.data.frame(function_mu) %>%
  rename(mean_eigen_function = function_mu)
function_df$node <- as.numeric(rownames(function_df))

nodes <- nodes %>%
  dplyr::select(-node_1, -node_2, -node_1_randomised, -node_2_randomised) %>%
  left_join(function_df, by = 'node')

## see where it falls -- identical
hist(nodes$mean_eigen_function)
plot(nodes$mean_eigen_function ~ nodes$sightings)
points(mean_eigen_function ~ sightings,
       data = nodes[nodes$id == 'M108',],
       pch = 19, col = 'red', cex = 2)

rm(list = ls()[! ls() %in% c('cdf','edge_samples','extract_eigen_centrality','fit_edges_anp','get_mean_predictions','n_chains','n_dyads','n_samples','n_windows','nodes','periods','post_pred_check','summary','time_window','traceplot')])

######## REMOVE ALL DYADS WHERE ONE OR OTHER IS <10 ####
nodes <- nodes %>% 
  filter(age > 9)
cdf <- cdf %>% 
  filter(age_start_1 > 9) %>% 
  filter(age_start_2 > 9)
edge_samples <- edge_samples[,as.character(cdf$dyad_id)]

#### check m108 data ####
m108 <- cdf %>% 
  filter(id_1 == 'M108' | id_2 == 'M108')
table(m108$event_count)
m108$period_count_1[m108$id_1 == 'M108']

m108_partners <- m108 %>% 
  filter(event_count > 0)

m108_edge_ids <- which(cdf$dyad_id %in% m108$dyad_id)
m108_edges <- edge_samples[,m108_edge_ids]#as.character(m108$dyad_id)]
other_edges <- edge_samples[,! (1:ncol(edge_samples) %in% m108_edge_ids)]

m108_means <- apply(m108_edges, 2, mean)
other_means <- apply(other_edges, 2, mean)
mean(m108_means)
mean(other_means)

m108_means <- as.data.frame(m108_means) %>% 
  rename(mean_edge = m108_means)
m108_means$dyad_id <- rownames(m108_means)
m108_means$m108 <- 'm108'

other_means <- as.data.frame(other_means) %>% 
  rename(mean_edge = other_means)
other_means$dyad_id <- rownames(other_means)
other_means$m108 <- 'other'

mean_edges <- rbind(m108_means, other_means)
ggplot(mean_edges)+
  geom_jitter(aes(x = m108, y = mean_edge))+
  geom_boxplot(aes(x = m108, y = mean_edge),
               notch = T,
               fill = 'seagreen1',
               colour = 'red')

m108_summary <- summary %>% 
  filter(dyad_id %in% m108$dyad_id) %>% 
  mutate(together = ifelse(event_count == 0, 'never together', 'sometimes together'))
ggplot(m108_summary)+
  geom_boxplot(aes(x = as.factor(event_count), y = median, fill = together),
               notch = T,
               alpha = 0.4)+
  geom_jitter(aes(x = as.factor(event_count), y = median, colour = together))+
  scale_fill_viridis_d(end = 0.5)+
  scale_colour_viridis_d(end = 0.5)

summary <- summary %>% 
  mutate(m108 = ifelse(dyad_id %in% m108$dyad_id,
                     'm108','other'),
         together = ifelse(event_count == 0, 'never together', 'sometimes together'))
ggplot(summary)+
  geom_jitter(aes(x = as.factor(event_count),
                  y = median,
                  colour = m108,
                  group = m108))+
  geom_boxplot(aes(x = as.factor(event_count),
                   y = median,
                   fill = m108),
               notch = T,
               alpha = 0.4)+
  scale_fill_viridis_d(end = 0.5)+
  scale_colour_viridis_d(end = 0.5)

sightings1 <- summary %>% 
  filter(period_count_1 %in% c(1) | 
           period_count_2 %in% c(1))

ggplot(sightings1)+
  # geom_jitter(aes(x = as.factor(event_count),
  #                 y = median,
  #                 colour = m108,
  #                 group = m108))+
  geom_boxplot(aes(x = as.factor(event_count),
                   y = median,
                   fill = m108),
               notch = T,
               alpha = 0.4)+
  scale_fill_viridis_d(end = 0.5)+
  scale_colour_viridis_d(end = 0.5)

m108_dyads_to_plot <- sample(1:ncol(m108_edges), 16, replace = F)
other_dyads_to_plot <- sample(1:ncol(other_edges), 16, replace = F)
par(mfrow = c(4,4),
    mai = c(0.1,0.1,0.1,0.1))
for(i in 1:16){
  plot(m108_edges[1:1000,m108_dyads_to_plot[i]], type = 'l', col = 'red')
  lines(m108_edges[1001:2000,m108_dyads_to_plot[i]], col = 'green')
  lines(m108_edges[2001:3000,m108_dyads_to_plot[i]], col = 'blue')
  lines(m108_edges[3001:4000,m108_dyads_to_plot[i]], col = 'purple')
}
for(i in 1:16){
  plot(other_edges[1:1000,other_dyads_to_plot[i]], type = 'l', col = 'red')
  lines(other_edges[1001:2000,other_dyads_to_plot[i]], col = 'green')
  lines(other_edges[2001:3000,other_dyads_to_plot[i]], col = 'blue')
  lines(other_edges[3001:4000,other_dyads_to_plot[i]], col = 'purple')
}

#### extract centrality -- all as one function ####
cent_new <- extract_eigen_centrality(nodes_df = nodes,
                                     dyads_df = cdf,
                                     edgeweight_matrix = edge_samples,
                                     logit = TRUE)

#### check ####
par(mfrow = c(1,1),
    mai = c(1,1,1,0.5))

## extract mean and covariance
colnames(cent_new) <- nodes$id
function_mu <- apply(cent_new, 2, mean)
function_cov <- cov(cent_new)

## add mean estimate per ID to nodes data frame
function_df <- as.data.frame(function_mu) %>%
  rename(mean_eigen_function = function_mu)
function_df$id <- rownames(function_df)

nodes <- nodes %>%
  #dplyr::select(-node_1, -node_2, -node_1_randomised, -node_2_randomised) %>%
  left_join(function_df, by = 'id')
colnames(nodes)[5:6] <- c('mean_eigen_all','mean_eigen_nobabies')

## see where it falls
hist(nodes$mean_eigen_function)
plot(nodes$mean_eigen_function ~ nodes$sightings)
points(mean_eigen_function ~ sightings,
       data = nodes[nodes$id == 'M108',],
       pch = 19, col = 'red', cex = 2)

###################################
rm(list = ls()) ; gc()

load('anp_nodalregression/run_using_all_elephants/anp_short_nodal.RData')
min(nodes_all$age) # just confirming that this is indeed the version that includes all elephants, not limited to those over 10

calves_inc <- nodes_all
calves_inc %>% 
  mutate(age_cat_full = ifelse(age < 10, '<10 yrs', age_cat_full)) %>% 
  mutate(age_cat_full = factor(age_cat_full,
                               levels = c('<10 yrs','10-15 yrs','16-20 yrs',
                                          '21-25 yrs','25-40 yrs', '>40 yrs'))) %>% 
  ggplot()+
  geom_histogram(aes(x = mean_eigen, fill = age_cat_full),
                 #colour = 'black',
                 binwidth = 0.1)+
  scale_fill_viridis_d(direction = -1)

rm(list = ls()[!ls() %in% 'nodes_under10inc']) ; gc()
load('anp_nodalregression/excluding_under10s/anp_short_nodal.RData')
min(nodes_under10exc$age)

no_babies <- nodes_all
no_babies %>% 
  mutate(age_cat_full = factor(age_cat_full,
                               levels = c('10-15 yrs','16-20 yrs',
                                          '21-25 yrs','25-40 yrs', '>40 yrs'))) %>% 
  ggplot()+
  geom_histogram(aes(x = mean_eigen, fill = age_cat_full),
                 #colour = 'black',
                 binwidth = 0.1)+
    scale_fill_viridis_d(direction = -1, end = 5/6)

# doesn't change the overall distribution

rm(nodes_all) ; gc()
no_babies_outliers <- no_babies %>% 
  filter(mean_eigen < -4)
calves_inc_outliers <- calves_inc %>% 
  filter(mean_eigen < -4)

calves_inc_outliers %>% 
  mutate(age_cat_full = ifelse(age < 10, '<10 yrs', age_cat_full)) %>% 
  mutate(age_cat_full = factor(age_cat_full,
                               levels = c('<10 yrs','10-15 yrs','16-20 yrs',
                                          '21-25 yrs','25-40 yrs', '>40 yrs'))) %>% 
  ggplot()+
  geom_histogram(aes(x = mean_eigen, fill = age_cat_full),
                 #colour = 'black',
                 binwidth = 0.1)+
  scale_fill_viridis_d(direction = -1)
no_babies_outliers %>% 
  mutate(age_cat_full = factor(age_cat_full,
                               levels = c('10-15 yrs','16-20 yrs',
                                          '21-25 yrs','25-40 yrs', '>40 yrs'))) %>% 
  ggplot()+
  geom_histogram(aes(x = mean_eigen, fill = age_cat_full),
                 #colour = 'black',
                 binwidth = 0.1)+
  scale_fill_viridis_d(direction = -1, end = 5/6)

calves_exc <- calves_inc %>% 
  filter(age > 9)

compare <- no_babies %>% 
  dplyr::select(id, node, window, node_window, age, age_cat_full, sightings, mean_eigen) %>% 
  left_join(calves_exc[,c('id','window','mean_eigen')], by = c('id','window')) %>% 
  rename(excluded_before_eigen = mean_eigen.x,
         excluded_after_eigen = mean_eigen.y)

compare %>% 
  mutate(age_cat_full = factor(age_cat_full,
                             levels = c('10-15 yrs','16-20 yrs',
                                        '21-25 yrs','25-40 yrs', '>40 yrs'))) %>% 
  ggplot()+
  geom_point(aes(x = excluded_before_eigen,
                 y = excluded_after_eigen,
                 colour = age_cat_full))+
  geom_abline(slope = 1, intercept = 0)+
  scale_colour_viridis_d()
##########################
colnames(cents_all)

nodes_all_1and20 <- nodes_all %>% 
  filter(window == 1 | window == 20) %>% 
  mutate(mean_eigen_new = apply(cents_all, 2, mean))

cents_df <- cents_all %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'node_window', values_to = 'eigen') %>% 
  left_join(nodes_all_1and20, by = 'node_window')

ggplot(nodes_all_1and20)+
  geom_point(aes(y = mean_eigen_new, x = age))+
  facet_wrap(. ~ window)




