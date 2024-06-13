#### information #####
# script to assess effect of age on centrality in ANP population
# model takes all time windows simultaneously, so all parameters are calculated together for ease of comparison to other populations

#### set up #####
## install cmdstanr because it's being a baby about loading normally
#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))  # Install the cmdstanr package
#library(cmdstanr)
#install_cmdstan()  

## load other packages
# library(LaplacesDemon) ; library(MASS) ; library(tidyverse) ; library(cmdstanr) ; library(sna)
library(cmdstanr, lib.loc = '../packages/')
library(tidyverse, lib.loc = '../packages/')
library(LaplacesDemon, lib.loc = '../packages/')
library(MASS, lib.loc = '../packages/')
#library(sna, lib.loc = '../packages/')
set_cmdstan_path('../packages/.cmdstan/cmdstan-2.34.1/') # Viking ; set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0/') # Desktop

## set theme for plots
theme_set(theme_bw())

#### create functions for use in both models ####
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

print('functions written')

#################### SHORT WINDOWS ####################
## set seed for reproducibility
set.seed(12345)

## define PDF output
pdf('../outputs/step4_nodalregression/anpshort_nodalregression_modelprep.pdf')

#### import data and randomise node IDs -- ensures that the model is not reading a correlation between node ID and age and partitioning too much of the variance to node ID #####
## load time window 1 and remove additional data
load('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')
rm(list = ls()[! ls() %in% c('nodes','n_windows','extract_eigen_centrality','traceplot','post_pred_check', 'get_mean_predictions')]) ; gc()

## create full data frame with window variable
nodes_all <- nodes %>%
  mutate(window = 1)

## attach all other data frames
for(time_window in 2:n_windows){
  ## import workspace image for time window
  load(paste0('anp_edgecalculations/anpshort',time_window,'_edgeweights_conditionalprior.RData'))
  rm(list = ls()[! ls() %in% c('nodes','nodes_all','time_window','extract_eigen_centrality','traceplot','post_pred_check','get_mean_predictions')]) ; gc()

  ## add window variable
  nodes <- nodes %>%
    mutate(window = time_window)

  ## attach to full data frame
  nodes_all <- rbind(nodes_all, nodes)

  ## print progress marker
  print(time_window)
}

## plot node ID vs window -- increasing node ID with window and age
ggplot(nodes_all)+
  geom_point(aes(x = node, y = window, colour = age))

## create data frame in which node is randomised
nodes_random <- nodes_all %>%
  dplyr::select(node) %>%
  distinct()
nodes_random$node_random <- sample(1:nrow(nodes_random), replace = FALSE)

## join to model data frame
nodes_all <- nodes_all %>%
  left_join(nodes_random, by = 'node') %>%
  mutate(age_std = (age - mean(age)) / sd(age))

## combine with window ID
nodes_all$node_window <- paste0(nodes_all$node_random, '_', nodes_all$window)

## plot node ID vs window -- no pattern of node ID with window and age
ggplot(nodes_all)+
  geom_point(aes(x = node_random, y = window, colour = age))

## clean up
rm(nodes_random, nodes) ; gc()
save.image('anp_nodalregression/anp_short_nodal_dataprep.RData')
write_csv(nodes_all, '../data_processed/step4_nodalregression/anp_allnodes_short.csv')
print('node data created')

#### reimport data and calculate centrality #####
# load('anp_nodalregression/run_using_all_elephants/anp_short_nodal_dataprep.RData')
## load time window 1 and remove additional data
load('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')
rm(list = ls()[! ls() %in% c('nodes_all', 'cdf_1', 'edge_samples','n_windows','extract_eigen_centrality','traceplot','post_pred_check', 'get_mean_predictions')]) ; gc()

## select nodes data frame for window 1 -- do this from nodes_all, not the nodes data frame already in the saved workspace, so now working from randomised node ID, not original ID
nodes <- nodes_all %>%
  filter(window == 1) %>%
  mutate(node_1 = node,
         node_2 = node,
         node_1_randomised = node_random,
         node_2_randomised = node_random)
nodes_all <- nodes_all %>%
  mutate(mean_eigen = NA) %>%
  filter(window != 1)

## randomise node IDs in dyads data frame
cdf_1 <- cdf_1 %>%
  left_join(nodes[,c('node_1','node_1_randomised')], by = 'node_1') %>%
  left_join(nodes[,c('node_2','node_2_randomised')], by = 'node_2')

## extract centrality for window 1
cents_all <- extract_eigen_centrality(nodes_df = nodes,
                                      dyads_df = cdf_1,
                                      edgeweight_matrix = edge_samples,
                                      window = 1)

## extract mean and covariance
cent_mu1 <- apply(cents_all, 2, mean)
cent_cov1 <- cov(cents_all)

## add mean estimate per ID to nodes data frame
mean_df <- as.data.frame(cent_mu1) %>%
  rename(mean_eigen = cent_mu1)
mean_df$node_window <- rownames(mean_df)
nodes <- nodes %>%
  dplyr::select(-node_1, -node_2, -node_1_randomised, -node_2_randomised) %>%
  left_join(mean_df, by = 'node_window')

## prep for combining everything to a single data frame
covs_all <- list()
covs_all[[1]] <- cent_cov1

## combine node data back into nodes_all
nodes_all <- rbind(nodes_all, nodes)

## clean up
rm(cent_mu1, cent_cov1, nodes, cdf_1, mean_df) ; gc()

## for loop
for(time_window in 2:n_windows){
  ## import workspace image
  load(paste0('anp_edgecalculations/anpshort',time_window,'_edgeweights_conditionalprior.RData'))
  rm(list = ls()[! ls() %in% c('nodes', 'cdf', 'edge_samples','n_windows','nodes_all','covs_all','cents_all','time_window','extract_eigen_centrality','traceplot','post_pred_check', 'get_mean_predictions')]) ; gc()

  ## select randomised nodes data frame for window
  nodes <- nodes_all %>%
    filter(window == time_window) %>%
    mutate(node_1 = node, node_2 = node,
           node_1_randomised = node_random, node_2_randomised = node_random) %>%
    dplyr::select(-mean_eigen)
  nodes_all <- nodes_all %>%
    filter(window != time_window)

  ## randomise node IDs in dyads data frame
  cdf <- cdf %>%
    left_join(nodes[,c('node_1','node_1_randomised')], by = 'node_1') %>%
    left_join(nodes[,c('node_2','node_2_randomised')], by = 'node_2')

  ## extract centrality
  cent_new <- extract_eigen_centrality(nodes_df = nodes,
                                       dyads_df = cdf,
                                       edgeweight_matrix = edge_samples,
                                       window = time_window)

  ## extract mean and covariance
  cent_mu <- apply(cent_new, 2, mean)
  cent_cov <- cov(cent_new)

  ## add mean estimate per ID to nodes data frame
  mean_df <- as.data.frame(cent_mu) %>%
    rename(mean_eigen = cent_mu)
  mean_df$node_window <- rownames(mean_df)
  nodes <- nodes %>%
    dplyr::select(-node_1, -node_2, -node_1_randomised, -node_2_randomised) %>%
    left_join(mean_df, by = 'node_window')

  ## combine everything to a single data frame
  nodes_all <- rbind(nodes_all, nodes)
  cents_all <- cbind(cents_all, cent_new)
  covs_all[[time_window]] <- cent_cov

  ## add progress marker
  print(time_window)
}

## clean up and save workspace
rm(cdf, cent_cov, cent_new, edge_samples, nodes, cent_mu, time_window, extract_eigen_centrality) ; gc()
save.image('anp_nodalregression/anp_short_nodal.RData')
print('data imported')

#### visualise centralities #####
df_plot <- data.frame(cents_all) %>%
  dplyr::select(sample(nodes_all$node_random, 30)) %>%
  pivot_longer(cols = everything(),
               names_to = 'node_window', values_to = 'centrality') %>%
  separate(node_window, into = c('X','node_window'), remove = T, sep = 1) %>%
  dplyr::select(-X) %>%
  separate(node_window, into = c('node_random','window'), sep = '_', remove = F) %>%
  left_join(nodes_all[,c('node_window','age')], by = 'node_window')
ggplot(data = df_plot, aes(x = centrality, fill = age, group = node_window)) +
  geom_density(linewidth = 0.4, alpha = 0.8) +
  facet_grid(rows = vars(as.factor(node_random)), scales = 'free') +
  labs(x = "Eigenvector centrality (standardised)") +
  theme_void() +
  theme(strip.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
print('centralities plotted')

#### check normal approximation #####
par(mfrow = c(3,3))
for(time_window in 1:n_windows){
  ## define mean and covariance
  mu <- nodes_all$mean_eigen[nodes_all$window == time_window]
  cov <- covs_all[[time_window]]

  ## simulate from multivariate normal
  sim_cent_samples <- MASS::mvrnorm(1e5, mu, cov)

  ## identify random node of interest
  node_id_sample <- sample(which(nodes_all$window == time_window),1)
  nodes_timewindow <- nodes_all %>%
    filter(window == time_window)
  node_id_check <- which(nodes_timewindow$node == nodes_all$node[node_id_sample])

  ## plot true density curve
  plot(density(cents_all[,node_id_sample]), lwd = 2, las = 1,
       main = paste0('time window ', time_window),
       xlab = '', ylab = '')

  ## plot normal approximation
  lines(density(sim_cent_samples[, node_id_check]),
        col = rgb(0,0,1,0.5), lwd = 2)
}
par(mfrow = c(1,1))

## clean up and save workspace
rm(mean_df, df_plot, nodes_timewindow, node_id_check, node_id_sample, time_window, mu, cov, sim_cent_samples) ; gc()
save.image('anp_nodalregression/anp_short_nodal_modelprep.RData')
print('normal approximation complete')

#### create data list #####
#load('anp_nodalregression/anp_short_nodal_modelprep.RData')
n_data <- nrow(nodes_all)
n_nodes <- length(unique(nodes_all$node_random))
n_windows <- length(unique(nodes_all$window))

model_data_list <- list(
  # global data size
  num_data = n_data,
  num_nodes = n_nodes,
  num_windows = n_windows,
  # per time window data size
  num_nodes_window1 = length(unique(nodes_all$node_random[nodes_all$window == 1])),
  num_nodes_window2 = length(unique(nodes_all$node_random[nodes_all$window == 2])),
  num_nodes_window3 = length(unique(nodes_all$node_random[nodes_all$window == 3])),
  num_nodes_window4 = length(unique(nodes_all$node_random[nodes_all$window == 4])),
  num_nodes_window5 = length(unique(nodes_all$node_random[nodes_all$window == 5])),
  num_nodes_window6 = length(unique(nodes_all$node_random[nodes_all$window == 6])),
  num_nodes_window7 = length(unique(nodes_all$node_random[nodes_all$window == 7])),
  num_nodes_window8 = length(unique(nodes_all$node_random[nodes_all$window == 8])),
  num_nodes_window9 = length(unique(nodes_all$node_random[nodes_all$window == 9])),
  num_nodes_window10 = length(unique(nodes_all$node_random[nodes_all$window == 10])),
  num_nodes_window11 = length(unique(nodes_all$node_random[nodes_all$window == 11])),
  num_nodes_window12 = length(unique(nodes_all$node_random[nodes_all$window == 12])),
  num_nodes_window13 = length(unique(nodes_all$node_random[nodes_all$window == 13])),
  num_nodes_window14 = length(unique(nodes_all$node_random[nodes_all$window == 14])),
  num_nodes_window15 = length(unique(nodes_all$node_random[nodes_all$window == 15])),
  num_nodes_window16 = length(unique(nodes_all$node_random[nodes_all$window == 16])),
  num_nodes_window17 = length(unique(nodes_all$node_random[nodes_all$window == 17])),
  num_nodes_window18 = length(unique(nodes_all$node_random[nodes_all$window == 18])),
  num_nodes_window19 = length(unique(nodes_all$node_random[nodes_all$window == 19])),
  num_nodes_window20 = length(unique(nodes_all$node_random[nodes_all$window == 20])),
  num_nodes_window21 = length(unique(nodes_all$node_random[nodes_all$window == 21])),
  num_nodes_window22 = length(unique(nodes_all$node_random[nodes_all$window == 22])),
  num_nodes_window23 = length(unique(nodes_all$node_random[nodes_all$window == 23])),
  num_nodes_window24 = length(unique(nodes_all$node_random[nodes_all$window == 24])),
  num_nodes_window25 = length(unique(nodes_all$node_random[nodes_all$window == 25])),
  num_nodes_window26 = length(unique(nodes_all$node_random[nodes_all$window == 26])),
  num_nodes_window27 = length(unique(nodes_all$node_random[nodes_all$window == 27])),
  num_nodes_window28 = length(unique(nodes_all$node_random[nodes_all$window == 28])),
  num_nodes_window29 = length(unique(nodes_all$node_random[nodes_all$window == 29])),
  num_nodes_window30 = length(unique(nodes_all$node_random[nodes_all$window == 30])),
  num_nodes_window31 = length(unique(nodes_all$node_random[nodes_all$window == 31])),
  num_nodes_window32 = length(unique(nodes_all$node_random[nodes_all$window == 32])),
  num_nodes_window33 = length(unique(nodes_all$node_random[nodes_all$window == 33])),
  num_nodes_window34 = length(unique(nodes_all$node_random[nodes_all$window == 34])),
  num_nodes_window35 = length(unique(nodes_all$node_random[nodes_all$window == 35])),
  num_nodes_window36 = length(unique(nodes_all$node_random[nodes_all$window == 36])),
  # number of nodes in all preceding time windows for node age indexing
  num_nodes_prev_windows = c(0,
                             length(which(nodes_all$window < 2)),
                             length(which(nodes_all$window < 3)),
                             length(which(nodes_all$window < 4)),
                             length(which(nodes_all$window < 5)),
                             length(which(nodes_all$window < 6)),
                             length(which(nodes_all$window < 7)),
                             length(which(nodes_all$window < 8)),
                             length(which(nodes_all$window < 9)),
                             length(which(nodes_all$window < 10)),
                             length(which(nodes_all$window < 11)),
                             length(which(nodes_all$window < 12)),
                             length(which(nodes_all$window < 13)),
                             length(which(nodes_all$window < 14)),
                             length(which(nodes_all$window < 15)),
                             length(which(nodes_all$window < 16)),
                             length(which(nodes_all$window < 17)),
                             length(which(nodes_all$window < 18)),
                             length(which(nodes_all$window < 19)),
                             length(which(nodes_all$window < 20)),
                             length(which(nodes_all$window < 21)),
                             length(which(nodes_all$window < 22)),
                             length(which(nodes_all$window < 23)),
                             length(which(nodes_all$window < 24)),
                             length(which(nodes_all$window < 25)),
                             length(which(nodes_all$window < 26)),
                             length(which(nodes_all$window < 27)),
                             length(which(nodes_all$window < 28)),
                             length(which(nodes_all$window < 29)),
                             length(which(nodes_all$window < 30)),
                             length(which(nodes_all$window < 31)),
                             length(which(nodes_all$window < 32)),
                             length(which(nodes_all$window < 33)),
                             length(which(nodes_all$window < 34)),
                             length(which(nodes_all$window < 35)),
                             length(which(nodes_all$window < 36))),
  # centrality means per time window
  centrality_mu_1 = nodes_all$mean_eigen[nodes_all$window == 1],
  centrality_mu_2 = nodes_all$mean_eigen[nodes_all$window == 2],
  centrality_mu_3 = nodes_all$mean_eigen[nodes_all$window == 3],
  centrality_mu_4 = nodes_all$mean_eigen[nodes_all$window == 4],
  centrality_mu_5 = nodes_all$mean_eigen[nodes_all$window == 5],
  centrality_mu_6 = nodes_all$mean_eigen[nodes_all$window == 6],
  centrality_mu_7 = nodes_all$mean_eigen[nodes_all$window == 7],
  centrality_mu_8 = nodes_all$mean_eigen[nodes_all$window == 8],
  centrality_mu_9 = nodes_all$mean_eigen[nodes_all$window == 9],
  centrality_mu_10 = nodes_all$mean_eigen[nodes_all$window == 10],
  centrality_mu_11 = nodes_all$mean_eigen[nodes_all$window == 11],
  centrality_mu_12 = nodes_all$mean_eigen[nodes_all$window == 12],
  centrality_mu_13 = nodes_all$mean_eigen[nodes_all$window == 13],
  centrality_mu_14 = nodes_all$mean_eigen[nodes_all$window == 14],
  centrality_mu_15 = nodes_all$mean_eigen[nodes_all$window == 15],
  centrality_mu_16 = nodes_all$mean_eigen[nodes_all$window == 16],
  centrality_mu_17 = nodes_all$mean_eigen[nodes_all$window == 17],
  centrality_mu_18 = nodes_all$mean_eigen[nodes_all$window == 18],
  centrality_mu_19 = nodes_all$mean_eigen[nodes_all$window == 19],
  centrality_mu_20 = nodes_all$mean_eigen[nodes_all$window == 20],
  centrality_mu_21 = nodes_all$mean_eigen[nodes_all$window == 21],
  centrality_mu_22 = nodes_all$mean_eigen[nodes_all$window == 22],
  centrality_mu_23 = nodes_all$mean_eigen[nodes_all$window == 23],
  centrality_mu_24 = nodes_all$mean_eigen[nodes_all$window == 24],
  centrality_mu_25 = nodes_all$mean_eigen[nodes_all$window == 25],
  centrality_mu_26 = nodes_all$mean_eigen[nodes_all$window == 26],
  centrality_mu_27 = nodes_all$mean_eigen[nodes_all$window == 27],
  centrality_mu_28 = nodes_all$mean_eigen[nodes_all$window == 28],
  centrality_mu_29 = nodes_all$mean_eigen[nodes_all$window == 29],
  centrality_mu_30 = nodes_all$mean_eigen[nodes_all$window == 30],
  centrality_mu_31 = nodes_all$mean_eigen[nodes_all$window == 31],
  centrality_mu_32 = nodes_all$mean_eigen[nodes_all$window == 32],
  centrality_mu_33 = nodes_all$mean_eigen[nodes_all$window == 33],
  centrality_mu_34 = nodes_all$mean_eigen[nodes_all$window == 34],
  centrality_mu_35 = nodes_all$mean_eigen[nodes_all$window == 35],
  centrality_mu_36 = nodes_all$mean_eigen[nodes_all$window == 36],
  # covariance matrix per time window
  centrality_cov_1 = covs_all[[1]],
  centrality_cov_2 = covs_all[[2]],
  centrality_cov_3 = covs_all[[3]],
  centrality_cov_4 = covs_all[[4]],
  centrality_cov_5 = covs_all[[5]],
  centrality_cov_6 = covs_all[[6]],
  centrality_cov_7 = covs_all[[7]],
  centrality_cov_8 = covs_all[[8]],
  centrality_cov_9 = covs_all[[9]],
  centrality_cov_10 = covs_all[[10]],
  centrality_cov_11 = covs_all[[11]],
  centrality_cov_12 = covs_all[[12]],
  centrality_cov_13 = covs_all[[13]],
  centrality_cov_14 = covs_all[[14]],
  centrality_cov_15 = covs_all[[15]],
  centrality_cov_16 = covs_all[[16]],
  centrality_cov_17 = covs_all[[17]],
  centrality_cov_18 = covs_all[[18]],
  centrality_cov_19 = covs_all[[19]],
  centrality_cov_20 = covs_all[[20]],
  centrality_cov_21 = covs_all[[21]],
  centrality_cov_22 = covs_all[[22]],
  centrality_cov_23 = covs_all[[23]],
  centrality_cov_24 = covs_all[[24]],
  centrality_cov_25 = covs_all[[25]],
  centrality_cov_26 = covs_all[[26]],
  centrality_cov_27 = covs_all[[27]],
  centrality_cov_28 = covs_all[[28]],
  centrality_cov_29 = covs_all[[29]],
  centrality_cov_30 = covs_all[[30]],
  centrality_cov_31 = covs_all[[31]],
  centrality_cov_32 = covs_all[[32]],
  centrality_cov_33 = covs_all[[33]],
  centrality_cov_34 = covs_all[[34]],
  centrality_cov_35 = covs_all[[35]],
  centrality_cov_36 = covs_all[[36]],
  # node IDs for all time windows
  nodes_window1 = nodes_all$node_random[nodes_all$window == 1],
  nodes_window2 = nodes_all$node_random[nodes_all$window == 2],
  nodes_window3 = nodes_all$node_random[nodes_all$window == 3],
  nodes_window4 = nodes_all$node_random[nodes_all$window == 4],
  nodes_window5 = nodes_all$node_random[nodes_all$window == 5],
  nodes_window6 = nodes_all$node_random[nodes_all$window == 6],
  nodes_window7 = nodes_all$node_random[nodes_all$window == 7],
  nodes_window8 = nodes_all$node_random[nodes_all$window == 8],
  nodes_window9 = nodes_all$node_random[nodes_all$window == 9],
  nodes_window10 = nodes_all$node_random[nodes_all$window == 10],
  nodes_window11 = nodes_all$node_random[nodes_all$window == 11],
  nodes_window12 = nodes_all$node_random[nodes_all$window == 12],
  nodes_window13 = nodes_all$node_random[nodes_all$window == 13],
  nodes_window14 = nodes_all$node_random[nodes_all$window == 14],
  nodes_window15 = nodes_all$node_random[nodes_all$window == 15],
  nodes_window16 = nodes_all$node_random[nodes_all$window == 16],
  nodes_window17 = nodes_all$node_random[nodes_all$window == 17],
  nodes_window18 = nodes_all$node_random[nodes_all$window == 18],
  nodes_window19 = nodes_all$node_random[nodes_all$window == 19],
  nodes_window20 = nodes_all$node_random[nodes_all$window == 20],
  nodes_window21 = nodes_all$node_random[nodes_all$window == 21],
  nodes_window22 = nodes_all$node_random[nodes_all$window == 22],
  nodes_window23 = nodes_all$node_random[nodes_all$window == 23],
  nodes_window24 = nodes_all$node_random[nodes_all$window == 24],
  nodes_window25 = nodes_all$node_random[nodes_all$window == 25],
  nodes_window26 = nodes_all$node_random[nodes_all$window == 26],
  nodes_window27 = nodes_all$node_random[nodes_all$window == 27],
  nodes_window28 = nodes_all$node_random[nodes_all$window == 28],
  nodes_window29 = nodes_all$node_random[nodes_all$window == 29],
  nodes_window30 = nodes_all$node_random[nodes_all$window == 30],
  nodes_window31 = nodes_all$node_random[nodes_all$window == 31],
  nodes_window32 = nodes_all$node_random[nodes_all$window == 32],
  nodes_window33 = nodes_all$node_random[nodes_all$window == 33],
  nodes_window34 = nodes_all$node_random[nodes_all$window == 34],
  nodes_window35 = nodes_all$node_random[nodes_all$window == 35],
  nodes_window36 = nodes_all$node_random[nodes_all$window == 36],
  # exposure variable
  node_age = nodes_all$age_std)

## check inputs
colour <- rep(c('red','orange','yellow','green','blue','purple'), each = 6)
shapes <- rep(c(3,4,15,16,17,18),6)
plot(model_data_list$centrality_mu_1 ~ model_data_list$node_age[1:model_data_list$num_nodes_window1],
     pch = shapes[1], col = colour[1],
     ylim = c(-10,0), xlim = c(-2.5,5))
for(i in 2:n_windows){
  if(i < n_windows){
    points(model_data_list[[(i+n_windows+4)]] ~ model_data_list$node_age[(model_data_list$num_nodes_prev_windows[i]+1):model_data_list$num_nodes_prev_windows[i+1]],
           pch = shapes[i], col = colour[i])
  } else {
    points(model_data_list[[(i+n_windows+4)]] ~ model_data_list$node_age[(model_data_list$num_nodes_prev_windows[i]+1):length(model_data_list$node_age)],
           pch = shapes[i], col = colour[i])
  }
}

print('data created')

#### prior predictive check ####
n <- 100
beta_age <- rnorm(n, 0, 0.8)
min_raw <- min(unlist(model_data_list[grep('centrality_mu', names(model_data_list))]))
max_raw <- max(unlist(model_data_list[grep('centrality_mu', names(model_data_list))]))
intercept  <- rnorm(n, LaplacesDemon::logit(0.05), 3) # taking the intercept from the logit of results from Chiyo 2011 (doesn't state mean/median centrality so estimated from graph based on where correlation line would cross x = 0)
plot(NULL, las = 1, xlab = 'age (standardised)', ylab = 'eigenvector',
     #ylim = c(min_raw-2, max_raw+2),
     ylim = c(-15, 5),
     xlim = c(min(model_data_list$node_age),
              max(model_data_list$node_age)))
abline(h = min_raw, lty = 2) ; abline(h = max_raw, lty = 2)
for(i in 1:n){
  lines(x = seq(min(nodes_all$age_std), max(nodes_all$age_std), length.out = 2),
        y = intercept[i] + beta_age[i]*c(min(nodes_all$age_std),
                                         max(nodes_all$age_std)),
        col = rgb(0,0,1,0.4))
}
rm(n, max_raw, min_raw, beta_age, intercept) ; gc()
dev.off()
print('prior predictive completed, model ready to run')

#### run model ####
## set model parameters
n_chains <- 4
n_samples <- 1000

## load model
nodal_regression <- cmdstan_model('models/eigen_regression_anp.stan',
                                  cpp_options = list(stan_threads = TRUE))

## run model
fit_anp_nodal <- nodal_regression$sample(data = model_data_list,
                                         chains = n_chains, parallel_chains = n_chains,
                                         threads_per_chain = 4,
                                         iter_warmup = n_samples, iter_sampling = n_samples)

## save workspace
save.image('anp_nodalregression/anp_short_nodal.RData')
print('model run')

#### check outputs ####
#load('anp_nodalregression/anp_short_nodal.RData')

## define PDF output
pdf('../outputs/step4_nodalregression/anpshort_nodalregression_modelchecks.pdf')

## extract model fit -- window effects have very low ESS, but Rhat is good for all
( summary <- fit_anp_nodal$summary() )
par(mfrow = c(3,1))
hist(summary$rhat, breaks = 50)
hist(summary$ess_bulk, breaks = 50)
hist(summary$ess_tail, breaks = 50)
par(mfrow = c(1,1))

## extract posterior
#params <- rstan::extract(fit_anp_nodal)
params <- fit_anp_nodal$draws(format = 'draws_df')

## separate random effects from global parameters
rand_window <- params %>%
  dplyr::select(grep('window_random_effect', colnames(params), value=TRUE))
rand_node <- params %>%
  dplyr::select(grep('node_random_effect', colnames(params), value=TRUE))

## traceplot all parameters
#traceplot(fit_anp_nodal, pars = c('intercept','beta_age','sigma','predictor[1]','predictor[50]','predictor[100]'))
plot_params <- c('intercept','beta_age','sigma') #,'nu')
plot_rand_nodes <- colnames(params)[sample(grep('node_random_effect',colnames(params)), size = 12, replace = F)]
plot_rand_windows1 <- c('window_random_effect[1]','window_random_effect[2]','window_random_effect[3]','window_random_effect[4]','window_random_effect[5]','window_random_effect[6]','window_random_effect[7]','window_random_effect[8]','window_random_effect[9]')
plot_rand_windows2 <- c('window_random_effect[10]','window_random_effect[11]','window_random_effect[12]','window_random_effect[13]','window_random_effect[14]','window_random_effect[15]','window_random_effect[16]','window_random_effect[17]','window_random_effect[18]')
plot_rand_windows3 <- c('window_random_effect[19]','window_random_effect[20]','window_random_effect[21]','window_random_effect[22]','window_random_effect[23]','window_random_effect[24]','window_random_effect[25]','window_random_effect[26]','window_random_effect[27]')
plot_rand_windows4 <- c('window_random_effect[28]','window_random_effect[29]','window_random_effect[30]','window_random_effect[31]','window_random_effect[32]','window_random_effect[33]','window_random_effect[34]','window_random_effect[35]','window_random_effect[36]')

# params %>%
#   select(all_of(plot_params),`.draw`,`.chain`,`.iteration`) %>%
#   pivot_longer(cols = all_of(plot_params), names_to = 'parameter', values_to = 'draw') %>%
#   rename(chain = .chain,
#          chain_position = .iteration,
#          draw_id = .draw) %>%
#   #filter(chain == 1) %>% # inspect individual chains -- check for wandery sections that might be hidden by other chains when all plotted together
#   ggplot(aes(x = chain_position, y = draw, colour = as.factor(chain)))+
#   geom_line()+
#   facet_wrap(. ~ parameter, scales = 'free_y')+
#   theme_bw()+
#   theme(legend.position = 'none')
traceplot(parameter_df = params, parameters_to_plot = plot_params)         # intercept looks terrible, beta age and sigma look decent
traceplot(parameter_df = params, parameters_to_plot = plot_rand_nodes)     # generally good
traceplot(parameter_df = params, parameters_to_plot = plot_rand_windows1)  # not horrendous, but also definitely not good
traceplot(parameter_df = params, parameters_to_plot = plot_rand_windows2)  # not horrendous, but also definitely not good
traceplot(parameter_df = params, parameters_to_plot = plot_rand_windows3)  # not horrendous, but also definitely not good
traceplot(parameter_df = params, parameters_to_plot = plot_rand_windows4)  # not horrendous, but also definitely not good
rm(plot_params,plot_rand_nodes,plot_rand_windows1,plot_rand_windows2,plot_rand_windows3,plot_rand_windows4) ; gc()

## save workspace
save.image('anp_nodalregression/anp_short_nodal.RData')
print('outputs checked')

#### posterior predictive check ####
par(mfrow = c(3,3))
for(time_window in 1:n_windows){
  post_pred_check(centrality_matrix = cents_all,
                  nodes_df = nodes_all, time_window = time_window,
                  cent_cov = covs_all[[time_window]],
                  parameters = params)
}
par(mfrow = c(1,1))

## define PDF output
dev.off()
pdf('../outputs/step4_nodalregression/anpshort_nodalregression_modelpredictions.pdf')
print('posterior predictive complete')

#### predict from model -- standardised scale ####
## get mean predictions
mu_predictions <- get_mean_predictions(predict_df = nodes_all, parameters = params,
                                       include_window = TRUE, include_node = TRUE)

## add mean and CI of predicted means to input data frame for comparison
nodes_all$mu_mean <- apply(mu_predictions, 2, mean)
nodes_all$mu_upr <- apply(mu_predictions, 2, quantile, prob = 0.975)
nodes_all$mu_lwr <- apply(mu_predictions, 2, quantile, prob = 0.025)

## plot mean of model vs mean of raw data
compare_plot <- ggplot()+
    geom_point(data = nodes_all, #size = 0.5,
               mapping = aes(x = mean_eigen, y = mu_mean, colour = age))+
    scale_colour_viridis_c()+
    facet_wrap(facets = . ~ as.factor(window))+
    labs(colour = 'window',
         x = 'raw eigen mean',
         y = 'predicted mean')+
    geom_abline(slope = 1, intercept = 0) # add line showing where points would lie if model fit was perfect
compare_plot +
  scale_x_continuous(limits = c(-7,-1))+
  scale_y_continuous(limits = c(-7,-1))
compare_plot +
  facet_wrap(facets = . ~ as.factor(window), scales = 'free')

## put together sigma arrays, separated by time window
num_nodes_windows <- unlist(model_data_list[grep('num_nodes_window', names(model_data_list))])
sigma_list <- list()
for(time_window in 1:n_windows) {
  sigma_window <- array(NA, dim = c(num_nodes_windows[time_window],
                                 num_nodes_windows[time_window],
                                 nrow(params)),
                     dimnames = list(nodes_all$nodes_random[nodes_all$window == time_window],
                                     nodes_all$nodes_random[nodes_all$window == time_window],
                                     NULL))
  cent_cov <- covs_all[[time_window]]
  for(i in 1:nrow(params)){
    sigma_window[,,i] <- cent_cov + diag(rep(params$sigma[i], num_nodes_windows[time_window]))
  }
  sigma_list[[time_window]] <- sigma_window
}

## create empty matrix to take full set of predicted values per elephant
full_predictions <- matrix(NA, nrow = nrow(params), ncol = nrow(nodes_all),
                      dimnames = list(NULL, nodes_all$node_window))

## populate matrix using mean values in matrix mu_std, and sigma values based on time window
for(time_window in 1:n_windows){
  sigma_window <- sigma_list[[time_window]]
  columns <- which(nodes_all$window == time_window)
  for(i in 1:nrow(full_predictions)){
    full_predictions[i,columns] <- MASS::mvrnorm(1,
                                                 mu_predictions[i,columns],
                                                 sigma_window[,,i])
  }
}

## add CI of predicted data points to input data frame for comparison
nodes_all$predict_lwr_std <- apply(full_predictions, 2, quantile, prob = 0.025)
nodes_all$predict_upr_std <- apply(full_predictions, 2, quantile, prob = 0.975)

## plot predictions
ggplot(nodes_all)+
  geom_ribbon(aes(x = age, fill = as.factor(window),
                  ymin = predict_lwr_std, ymax = predict_upr_std),
              alpha = 0.2)+                    # 95% CI of all predictions
  geom_ribbon(aes(x = age, fill = as.factor(window),
                  ymin = mu_lwr, ymax = mu_upr),
              alpha = 0.4)+                    # 95% CI of predicted means
  geom_line(aes(x = age, y = mu_mean,
                colour = as.factor(window)))+  # line showing mean of predicted means
  geom_point(aes(x = age, y = mean_eigen),
             size = 0.5)+                      # raw data points
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  scale_x_continuous(expand = c(0,0))+
  facet_wrap(. ~ as.factor(window),
             scales = 'free_y')+                 # separate plots per window
  theme(legend.position = 'none')+
  #guides(colour = guide_legend(nrow = 6),
  #       fill = guide_legend(nrow = 6))+
  labs(#colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'age (years)')

## save workspace
save.image('anp_nodalregression/anp_short_nodal.RData')
print('predictions made')

#### extract original values from output -- The model works on the standardised scale. Convert predictions to unstandardised scale and then run contrasts to calculate the slope coefficient. ####
## get mean predictions
newdata <- nodes_all %>%
  mutate(age_std = age_std+1)
mu_predictions_newdata <- get_mean_predictions(predict_df = newdata, parameters = params,
                                               include_window = TRUE, include_node = TRUE)
nodes_all$mu_mean_plus1 <- apply(mu_predictions_newdata, 2, mean)

## full distribution of predictions using age_std + 1 stdev
full_predictions_newdata <- matrix(NA, nrow = nrow(params), ncol = nrow(newdata),
                                   dimnames = list(NULL, nodes_all$node_window))
for(time_window in 1:n_windows){
  sigma_window <- sigma_list[[time_window]]
  columns <- which(nodes_all$window == time_window)
  for(i in 1:nrow(full_predictions_newdata)){
    full_predictions_newdata[i,columns] <- MASS::mvrnorm(1,
                                                         mu_predictions_newdata[i,columns],
                                                         sigma_window[,,i])
  }
}

## contrast predictions on standardised scale -- check that this returns the marginal effect presented in the summary
contrast_std <- full_predictions_newdata - full_predictions    # contrast between predicted values for raw data and all same data but add 1 to age
head(contrast_std[,1:5])                              # check matrix looks right
mean(params$beta_age)                                 # output parameter direct from model
mean(contrast_std)                                    # should be very similar to mean of output parameter
quantile(contrast_std, prob = c(0.025, 0.975))        # very wide

## contrast predictions on output scale -- reportable values of effect of plus 1 year
contrast <- contrast_std / sd(nodes_all$age)          # convert to outcome scale
head(contrast[,1:5])                                  # check matrix looks right
mean(contrast)                                        # effect size on outcome scale
sd(contrast)                                          # size of uncertainty on outcome scale
quantile(contrast, prob = c(0.025, 0.975))            # very wide

## save workspace
save.image('anp_nodalregression/anp_short_nodal.RData')
print('contrasts complete')

#### final "clean" plots using hypothetical data ####
#load('anp_nodalregression/anp_short_nodal.RData')

## create dummy dataset
newdat <- nodes_all %>%
  dplyr::select(node_random, age, age_std, window)

## get mean predictions
fake_mean <- get_mean_predictions(predict_df = newdat, parameters = params,
                                  include_window = TRUE, include_node = FALSE)
newdat$predict_mean <- apply(fake_mean, 2, mean)
newdat$predict_mean_lwr <- apply(fake_mean, 2, quantile, prob = 0.025)
newdat$predict_mean_upr <- apply(fake_mean, 2, quantile, prob = 0.975)

## create prediction matrix
fake_pred <- matrix(NA, nrow = nrow(params), ncol = nrow(newdat))
for(time_window in 1:n_windows){
  sigma_window <- sigma_list[[time_window]]
  columns <- which(nodes_all$window == time_window)
  for(i in 1:nrow(full_predictions_newdata)){
    fake_pred[i,columns] <- MASS::mvrnorm(1,
                                          fake_mean[i,columns],
                                          sigma_window[,,i])
  }
}
print('predictions for clean plots made')

## add CI of predicted data points to input data frame for comparison
newdat$predict_pred_lwr <- apply(fake_pred, 2, quantile, prob = 0.025)
newdat$predict_pred_upr <- apply(fake_pred, 2, quantile, prob = 0.975)

## plot predictions
newdat_summary <- newdat %>%
  group_by(age, window) %>%
  mutate(predict_pred_lwr = mean(predict_pred_lwr),
         predict_pred_upr = mean(predict_pred_upr)) %>%
  dplyr::select(age, predict_pred_lwr, predict_pred_upr,window) %>%
  distinct()

## plot mean values
ggplot()+
  geom_ribbon(data = newdat_summary,
              aes(ymin = predict_pred_lwr, ymax = predict_pred_upr,
                  x = age, fill = as.factor(window)),
              alpha = 0.2)+                      # 95% CI of all predictions
  geom_ribbon(data = newdat,
              aes(ymin = predict_mean_lwr, ymax = predict_mean_upr,
                  x = age, fill = as.factor(window)),
              alpha = 0.4)+                      # 95% CI of predicted means
  geom_line(data = newdat,
            aes(x = age, y = predict_mean, colour = as.factor(window)))+  # predicted means
  geom_point(data = nodes_all,
             aes(x = age, y = mean_eigen))+      # original data points
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  facet_wrap(. ~ as.factor(window))+             # separate plots per window
  #theme(legend.position = 'bottom')+
  theme(legend.position = 'none')+
  labs(colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'age (years)')
print('mean values plotted')

## convert raw data to data frame
cents_all_df <- cents_all %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = 'node_window', values_to = 'eigen') %>%
  left_join(nodes_all[,c('node_window','age','window')], by = 'node_window')

## plot full distribution
(final_plot <- ggplot()+
  geom_ribbon(data = newdat_summary,
              aes(ymin = predict_pred_lwr, ymax = predict_pred_upr,
                  x = age, fill = as.factor(window)),
              alpha = 0.2)+                      # 95% CI of all predictions
  geom_ribbon(data = newdat,
              aes(ymin = predict_mean_lwr, ymax = predict_mean_upr,
                  x = age, fill = as.factor(window)),
              alpha = 0.4)+                      # 95% CI of predicted means
  geom_point(data = cents_all_df,
             aes(x = age, y = eigen),
             alpha = 0.01, size = 0.5)+         # original data points
  geom_line(data = newdat,
            aes(x = age, y = predict_mean, colour = as.factor(window)))+  # predicted means
  geom_point(data = nodes_all,
             aes(x = age, y = mean_eigen),
             size = 0.5, colour = 'white')+      # original mean data points
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  #facet_wrap(. ~ as.factor(window))+             # separate plots per window
  theme(legend.position = 'bottom')+
  labs(colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'age (years)'))
ggsave(plot = final_plot, device = 'svg', width = 800, height = 800, units = 'px',
       filename = '../outputs/step4_nodalregression/anpshort_nodalregression_finalplot.svg')
ggsave(plot = final_plot, device = 'png', width = 800, height = 800, units = 'px',
       filename = '../outputs/step4_nodalregression/anpshort_nodalregression_finalplot.png')

## clean up and save workspace
save.image('anp_nodalregression/anp_short_nodal.RData')
dev.off()
print('short windows complete')

#### smooth plot -- I think this should replace the code chunk above, but check before deleting ####
# load('anp_nodalregression/anp_short_nodal.RData')
rm(cent_cov, compare_plot, contrast, contrast_std, fake_mean, fake_pred, fit_anp_nodal, full_predictions, full_predictions_newdata, model_data_list, mu_predictions, mu_predictions_newdata, newdata, nodal_regression, rand_node, rand_window, summary, colour, shapes) ; gc()

## create age categorical variables to match MOTNP
nodes_all <- nodes_all %>% 
  mutate(age_cat = ifelse(age < 15, 1,
                          ifelse(age < 20, 2,
                                 ifelse(age < 25, 3,
                                        ifelse(age < 40, 4, 5)))),
         age_cat_full = ifelse(age < 15, '10-15 yrs',
                          ifelse(age < 20, '16-20 yrs',
                                 ifelse(age < 25, '21-25 yrs',
                                        ifelse(age < 40,
                                               '25-40 yrs',
                                               '>40 yrs')))))
         
## predict from raw data but exclude window and node effects (except for in sigma_window as I think that's unavoidable)
set.seed(12345) ; draws <- sample(x = 1:n_samples, size = 125, replace = F)
pred_mean_norandom <- get_mean_predictions(predict_df = nodes_all,
                                           parameters = params[c(draws,
                                                                 draws+1000,
                                                                 draws+2000,
                                                                 draws+3000),],
                                           include_window = FALSE, include_node = FALSE)

## create data frame for smooth predictions
norandom_summary <- nodes_all %>%
  dplyr::select(age, age_std, age_cat) %>%
  distinct() %>%
  mutate(predict_mean_logit = NA, predict_mean_invlogit = NA,
         predict_sd_logit = NA, predict_sd_invlogit = NA,
         predict_mean_lwr_logit = NA, predict_mean_lwr_invlogit = NA,
         predict_mean_upr_logit = NA, predict_mean_upr_invlogit = NA,
         predict_full_lwr_logit = NA, predict_full_lwr_invlogit = NA,
         predict_full_upr_logit = NA, predict_full_upr_invlogit = NA)
for(i in 1:nrow(norandom_summary)){
  ## select only columns that refer to elephants of that age
  mean_predictions_per_age <- pred_mean_norandom[,which(nodes_all$age == norandom_summary$age[i])]
  
  ## calculate values for mean line
  norandom_summary$predict_mean_logit[i] <- mean(mean_predictions_per_age)
  norandom_summary$predict_mean_invlogit[i] <- mean(invlogit(mean_predictions_per_age))
  
  ## calculate standard deviation
  norandom_summary$predict_sd_logit[i] <- sd(mean_predictions_per_age)
  norandom_summary$predict_sd_invlogit[i] <- sd(invlogit(mean_predictions_per_age))
  
  ## calculate values for shaded ribbon
  norandom_summary$predict_mean_lwr_logit[i] <- quantile(mean_predictions_per_age,
                                                         prob = 0.025)
  norandom_summary$predict_mean_lwr_invlogit[i] <- quantile(invlogit(mean_predictions_per_age),
                                                            prob = 0.025)
  norandom_summary$predict_mean_upr_logit[i] <- quantile(mean_predictions_per_age,
                                                         prob = 0.975)
  norandom_summary$predict_mean_upr_invlogit[i] <- quantile(invlogit(mean_predictions_per_age),
                                                            prob = 0.975)
}

## create empty matrix to take full set of predicted values per elephant
pred_full_norandom <- matrix(NA, nrow = nrow(pred_mean_norandom), ncol = nrow(nodes_all),
                             dimnames = list(NULL, nodes_all$node_window))

## populate matrix using mean values in matrix mu_std, and sigma values based on time window
for(time_window in 1:n_windows){
  sigma_window <- sigma_list[[time_window]]
  columns <- which(nodes_all$window == time_window)
  for(i in 1:nrow(pred_full_norandom)){
    pred_full_norandom[i,columns] <- MASS::mvrnorm(1,
                                                   pred_mean_norandom[i,columns],
                                                   sigma_window[,,i])
  }
}

## add CI of predicted data points to input data frame for comparison
for(i in 1:nrow(norandom_summary)){
  ## select only columns that refer to elephants of that age
  mean_predictions_per_age <- pred_full_norandom[,which(nodes_all$age == norandom_summary$age[i])]
  
  ## calculate values for shaded ribbon
  norandom_summary$predict_full_lwr_logit[i] <- quantile(mean_predictions_per_age, 
                                                        prob = 0.025)
  norandom_summary$predict_full_lwr_invlogit[i] <- quantile(invlogit(mean_predictions_per_age),
                                                            prob = 0.025)
  norandom_summary$predict_full_upr_logit[i] <- quantile(mean_predictions_per_age,
                                                         prob = 0.975)
  norandom_summary$predict_full_upr_invlogit[i] <- quantile(invlogit(mean_predictions_per_age),
                                                            prob = 0.975)
}

# ## add full distribution (probably unnecessary and becomes too big for computer to open)
# cents_all_df <- cents_all %>%
#   as.data.frame() %>%
#   pivot_longer(cols = everything(), names_to = 'node_window', values_to = 'eigen') %>%
#   left_join(nodes_all[,c('node_window','age','window')], by = 'node_window')

## faceted by window
(faceted <- ggplot()+
    geom_ribbon(data = newdat_summary,
                aes(ymin = invlogit(predict_pred_lwr),
                    ymax = invlogit(predict_pred_upr),
                    group = as.factor(window),#fill = as.factor(window),
                    x = age),
                alpha = 0.2)+                      # 95% CI of all predictions
    geom_ribbon(data = newdat,
                aes(ymin = invlogit(predict_mean_lwr),
                    ymax = invlogit(predict_mean_upr),
                    group = as.factor(window),#fill = as.factor(window),
                    x = age),
                alpha = 0.4)+                      # 95% CI of predicted means
    # geom_point(data = cents_all_df,
    #            aes(x = age, y = eigen, colour = as.factor(window)),
    #            alpha = 0.01, size = 0.5)+         # original data points
    geom_point(data = distinct(nodes_all[,c('age','mean_eigen','window','age_cat_full','sightings')]),
               #size = 0.5, colour = 'white',
               aes(x = age,
                   y = invlogit(mean_eigen),
                   colour = factor(age_cat_full,
                                   levels = c('10-15 yrs','16-20 yrs',
                                              '21-25 yrs', '25-40 yrs',
                                              '>40 yrs')),
                   #group = as.factor(window),
                   size = sightings
               ),
               alpha = 0.5)+      # original mean data points
    geom_line(data = newdat,
              aes(x = age,
                  y = invlogit(predict_mean),
                  group = as.factor(window)))+  # predicted means
    scale_colour_viridis_d(direction = -1)+ # begin = 0, end = 0.7
    theme(legend.position = 'bottom',
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 18))+
    facet_wrap(. ~ as.factor(window)) +
    guides(#colour = guide_legend(nrow = 5),
      colour = 'none',
      size = guide_legend(override.aes = list(fill = 'grey80')))+
    labs(colour = '  age category\n(as in MOTNP)',
       y = 'eigenvector centrality', x = 'age (years)'))
ggsave(plot = faceted, device = 'svg', width = 2400, height = 2400, units = 'px',
       filename = 'anpshort_nodalregression_faceted.svg',
       path = '../outputs/step4_nodalregression/')
ggsave(plot = faceted, device = 'png', width = 2400, height = 2400, units = 'px',
       filename = 'anpshort_nodalregression_faceted.png',
       path = '../outputs/step4_nodalregression/')

## overall result
(overall <- ggplot()+
    geom_ribbon(data = norandom_summary,
                aes(x = age,
                    ymin = predict_full_lwr_invlogit, #ymin = predict_full_lwr_logit,
                    ymax = predict_full_upr_invlogit), #ymax = predict_full_upr_logit),
                alpha = 0.2)+
    # geom_point(data = cents_all_df,
    #            aes(x = age, y = eigen, colour = as.factor(window)),
    #            alpha = 0.01, size = 0.5)+         # original data points
    geom_point(data = distinct(nodes_all[,c('age','mean_eigen','age_cat_full','sightings')]),
               aes(x = age,
                   y = invlogit(mean_eigen),
                   colour = factor(age_cat_full,
                                   levels = c('10-15 yrs','16-20 yrs',
                                              '21-25 yrs', '25-40 yrs',
                                              '>40 yrs')),
                   size = sightings),
               alpha = 0.5)+
    geom_ribbon(data = norandom_summary,
                aes(x = age,
                    ymin = predict_mean_lwr_invlogit, #ymin = predict_mean_lwr_logit,
                    ymax = predict_mean_upr_invlogit), #ymax = predict_mean_upr_logit),
                alpha = 0.4)+
    geom_line(data = norandom_summary,
              aes(x = age,
                  y = predict_mean_invlogit
                  ))+
    scale_colour_viridis_d(direction = -1)+ # begin = 0, end = 0.7
    theme(legend.position = 'bottom',
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 18))+
    guides(#colour = guide_legend(nrow = 5),
      colour = 'none',
      size = guide_legend(override.aes = list(fill = 'grey80')))+
    labs(colour = '  age category\n(as in MOTNP)',
         y = 'eigenvector centrality', x = 'age (years)'))
ggsave(plot = overall, device = 'svg', width = 2400, height = 1800, units = 'px',
       filename = 'anpshort_nodalregression_nofaceting.svg',
       path = '../outputs/step4_nodalregression/')
ggsave(plot = overall, device = 'png', width = 2400, height = 1800, units = 'px',
       filename = 'anpshort_nodalregression_nofaceting.png',
       path = '../outputs/step4_nodalregression/')

save.image('anp_nodalregression/anp_short_nodal.RData')
rm(list = ls()) ; gc()

# ## create dummy dataset -- CAN ONLY PLOT MEAN SHADING FROM THIS, NOT OVERALL PREDICTIONS 95% CI: PREDICTIONS USE MVRNORM WHICH REQUIRES THE SIGMA MATRIX TO HAVE THE SAME NUMBER OF COLUMNS AS THERE ARE ELEPHANTS IN THE TIME WINDOW, SO CAN'T MAKE ALL WINDOWS CONTAIN ALL ELEPHANTS FOR THE SAKE OF A NICE HYPOTHETICAL PLOT
# newdat <- expand_grid(unique(nodes_all$node_random),
#                       unique(nodes_all$age_std),
#                       unique(nodes_all$window)) %>% 
#   rename(node_random = `unique(nodes_all$node_random)`,
#          age_std = `unique(nodes_all$age_std)`,
#          window = `unique(nodes_all$window)`) %>% 
#   left_join(distinct(nodes_all[,c('age','age_std','age_cat')]), by = 'age_std') %>% 
#   mutate(node_window = paste0(node_random, '_', window))
# 
# ## get mean predictions
# set.seed(12345) ; draws <- sample(x = 1:n_samples, size = 125, replace = F)
# fake_mean <- get_mean_predictions(predict_df = newdat,
#                                   parameters = params[c(draws,
#                                                         draws+1000,
#                                                         draws+2000,
#                                                         draws+3000),],
#                                   include_window = FALSE, include_node = FALSE)
# 
# ## create data frame for smooth predictions
# fake_summary <- newdat %>% 
#   dplyr::select(age, age_std, age_cat) %>% 
#   distinct() %>% 
#   mutate(predict_mean = NA,
#          predict_sd = NA,
#          predict_mean_lwr = NA,
#          predict_mean_upr = NA)
# for(i in 1:nrow(fake_summary)){
#   mean_predictions_per_age <- fake_mean[,which(newdat$age == fake_summary$age[i])]
#   fake_summary$predict_mean[i] <- mean(mean_predictions_per_age)
#   fake_summary$predict_sd[i] <- sd(mean_predictions_per_age)
#   fake_summary$predict_mean_lwr[i] <- quantile(mean_predictions_per_age, prob = 0.025)
#   fake_summary$predict_mean_upr[i] <- quantile(mean_predictions_per_age, prob = 0.975)
# }
# 
# ## plot mean of means
# fake_summary$age_cat_full <- ifelse(fake_summary$age_cat == 1, '10-15 yrs',
#                                     ifelse(fake_summary$age_cat == 2, '16-20 yrs',
#                                            ifelse(fake_summary$age_cat == 3, '21-25 yrs',
#                                                   ifelse(fake_summary$age_cat == 4, '25-40 yrs',
#                                                          '>40 yrs'))))
# 
# #### just take the mean of all raw predictions -- DOESN'T WORK TO GIVE A SMOOTH LINE BECAUSE TOO MCUH VARIATION FROM INDIVIDUALS AND WINDOWS STILL
# smooth_mean <- data.frame(age = min(nodes_all$age):max(nodes_all$age),
#                           mean = NA, lwr_mu = NA, upr_mu = NA,
#                           lwr_full = NA, upr_full = NA)
# for(i in 1:nrow(smooth_mean)){
#   age_mean_pred <- mu_predictions[,which(nodes_all$age == smooth_mean$age[i])]
#   age_full_pred <- full_predictions[,which(nodes_all$age == smooth_mean$age[i])]
#   smooth_mean$mean[i] <- mean(age_mean_pred)
#   smooth_mean$lwr_mu[i] <- quantile(age_mean_pred, prob = 0.025)
#   smooth_mean$upr_mu[i] <- quantile(age_mean_pred, prob = 0.975)
#   smooth_mean$lwr_full[i] <- quantile(age_full_pred, prob = 0.025)
#   smooth_mean$upr_full[i] <- quantile(age_full_pred, prob = 0.975)
# }
# 
#################### LONG WINDOWS ####################
set.seed(12345)

## define PDF output
pdf('../outputs/step4_nodalregression/anplong_nodalregression_modelprep.pdf')

#### import data and randomise node IDs -- ensures that the model is not reading a correlation between node ID and age and partitioning too much of the variance to node ID #####
## load time window 1 and remove additional data
load('anp_edgecalculations/anplong1_edgeweights_conditionalprior.RData')
rm(list = ls()[! ls() %in% c('nodes','n_windows','extract_eigen_centrality','traceplot','post_pred_check',
                             'get_mean_predictions')]) ; gc()

## create full data frame with window variable
nodes_all <- nodes %>%
  mutate(window = 1)

## attach all other data frames
for(time_window in 2:n_windows){
  ## import workspace image for time window
  load(paste0('anp_edgecalculations/anplong',time_window,'_edgeweights_conditionalprior.RData'))
  rm(list = ls()[! ls() %in% c('nodes','nodes_all','time_window','extract_eigen_centrality','traceplot','post_pred_check', 'get_mean_predictions')]) ; gc()

  ## add window variable
  nodes <- nodes %>%
    mutate(window = time_window)

  ## attach to full data frame
  nodes_all <- rbind(nodes_all, nodes)

  ## print progress marker
  print(time_window)
}

## plot node ID vs window -- increasing node ID with window and age
ggplot(nodes_all)+
  geom_point(aes(x = node, y = window, colour = age))

## create data frame in which node is randomised
nodes_random <- nodes_all %>%
  dplyr::select(node) %>%
  distinct()
nodes_random$node_random <- sample(1:nrow(nodes_random), replace = FALSE)

## join to model data frame
nodes_all <- nodes_all %>%
  left_join(nodes_random, by = 'node') %>%
  mutate(age_std = (age - mean(age)) / sd(age))

## combine with window ID
nodes_all$node_window <- paste0(nodes_all$node_random, '_', nodes_all$window)

## plot node ID vs window -- no pattern of node ID with window and age
ggplot(nodes_all)+
  geom_point(aes(x = node_random, y = window, colour = age))

## clean up
rm(nodes_random, nodes) ; gc()
save.image('anp_nodalregression/anp_long_nodal_dataprep.RData')
write_csv(nodes_all, '../data_processed/step4_nodalregression/anp_allnodes_long.csv')
print('node data created')

#### reimport data and calculate centrality #####
# load('anp_nodalregression/anp_long_nodal_dataprep.RData')
## load time window 1 and remove additional data
load('anp_edgecalculations/anplong1_edgeweights_conditionalprior.RData')
rm(list = ls()[! ls() %in% c('nodes_all', 'cdf', 'edge_samples','n_windows','extract_eigen_centrality','traceplot','post_pred_check', 'get_mean_predictions')]) ; gc()

## select nodes data frame for window 1 -- do this from nodes_all, not the nodes data frame already in the saved workspace, so now working from randomised node ID, not original ID
nodes <- nodes_all %>%
  filter(window == 1) %>%
  mutate(node_1 = node,
         node_2 = node,
         node_1_randomised = node_random,
         node_2_randomised = node_random)
nodes_all <- nodes_all %>%
  mutate(mean_eigen = NA) %>%
  filter(window != 1)

## randomise node IDs in dyads data frame
cdf <- cdf %>%
  left_join(nodes[,c('node_1','node_1_randomised')], by = 'node_1') %>%
  left_join(nodes[,c('node_2','node_2_randomised')], by = 'node_2')

## extract centrality for window 1
cents_all <- extract_eigen_centrality(nodes_df = nodes,
                                      dyads_df = cdf,
                                      edgeweight_matrix = edge_samples,
                                      window = 1)

## extract mean and covariance
cent_mu1 <- apply(cents_all, 2, mean)
cent_cov1 <- cov(cents_all)

## add mean estimate per ID to nodes data frame
mean_df <- as.data.frame(cent_mu1) %>%
  rename(mean_eigen = cent_mu1)
mean_df$node_window <- rownames(mean_df)
nodes <- nodes %>%
  dplyr::select(-node_1, -node_2, -node_1_randomised, -node_2_randomised) %>%
  left_join(mean_df, by = 'node_window')

## prep for combining everything to a single data frame
covs_all <- list()
covs_all[[1]] <- cent_cov1

## combine node data back into nodes_all
nodes_all <- rbind(nodes_all, nodes)

## clean up
rm(cent_mu1, cent_cov1, nodes, cdf, mean_df) ; gc()

## for loop
for(time_window in 2:n_windows){
  ## import workspace image
  load(paste0('anp_edgecalculations/anplong',time_window,'_edgeweights_conditionalprior.RData'))
  rm(list = ls()[! ls() %in% c('nodes', 'cdf', 'edge_samples','n_windows','nodes_all','covs_all','cents_all','time_window','extract_eigen_centrality','traceplot','post_pred_check', 'get_mean_predictions')]) ; gc()

  ## select randomised nodes data frame for window
  nodes <- nodes_all %>%
    filter(window == time_window) %>%
    mutate(node_1 = node, node_2 = node,
           node_1_randomised = node_random, node_2_randomised = node_random) %>%
    dplyr::select(-mean_eigen)
  nodes_all <- nodes_all %>%
    filter(window != time_window)

  ## randomise node IDs in dyads data frame
  cdf <- cdf %>%
    left_join(nodes[,c('node_1','node_1_randomised')], by = 'node_1') %>%
    left_join(nodes[,c('node_2','node_2_randomised')], by = 'node_2')

  ## extract centrality
  cent_new <- extract_eigen_centrality(nodes_df = nodes,
                                       dyads_df = cdf,
                                       edgeweight_matrix = edge_samples,
                                       window = time_window)

  ## extract mean and covariance
  cent_mu <- apply(cent_new, 2, mean)
  cent_cov <- cov(cent_new)

  ## add mean estimate per ID to nodes data frame
  mean_df <- as.data.frame(cent_mu) %>%
    rename(mean_eigen = cent_mu)
  mean_df$node_window <- rownames(mean_df)
  nodes <- nodes %>%
    dplyr::select(-node_1, -node_2, -node_1_randomised, -node_2_randomised) %>%
    left_join(mean_df, by = 'node_window')

  ## combine everything to a single data frame
  nodes_all <- rbind(nodes_all, nodes)
  cents_all <- cbind(cents_all, cent_new)
  covs_all[[time_window]] <- cent_cov

  ## add progress marker
  print(time_window)
}

## clean up and save workspace
rm(cdf, cent_cov, cent_new, edge_samples, nodes, cent_mu, time_window, extract_eigen_centrality) ; gc()
save.image('anp_nodalregression/anp_long_nodal_modelprep.RData')
print('data imported')

#### visualise centralities #####
df_plot <- data.frame(cents_all) %>%
  dplyr::select(sample(nodes_all$node_random, 30)) %>%
  pivot_longer(cols = everything(),
               names_to = 'node_window', values_to = 'centrality') %>%
  separate(node_window, into = c('X','node_window'), remove = T, sep = 1) %>%
  dplyr::select(-X) %>%
  separate(node_window, into = c('node_random','window'), sep = '_', remove = F) %>%
  left_join(nodes_all[,c('node_window','age')], by = 'node_window')
ggplot(data = df_plot, aes(x = centrality, fill = age, group = node_window)) +
  geom_density(linewidth = 0.4, alpha = 0.8) +
  facet_grid(rows = vars(as.factor(node_random)), scales = 'free') +
  labs(x = "Eigenvector centrality (standardised)") +
  theme_void() +
  theme(strip.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
print('centralities plotted')

#### check normal approximation #####
par(mfrow = c(3,3))
for(time_window in 1:n_windows){
  ## define mean and covariance
  mu <- nodes_all$mean_eigen[nodes_all$window == time_window]
  cov <- covs_all[[time_window]]

  ## simulate from multivariate normal
  sim_cent_samples <- MASS::mvrnorm(1e5, mu, cov)

  ## identify random node of interest
  node_id_sample <- sample(which(nodes_all$window == time_window),1)
  nodes_timewindow <- nodes_all %>%
    filter(window == time_window)
  node_id_check <- which(nodes_timewindow$node == nodes_all$node[node_id_sample])

  ## plot true density curve
  plot(density(cents_all[,node_id_sample]), lwd = 2, las = 1,
       main = paste0('time window ', time_window),
       xlab = '', ylab = '')

  ## plot normal approximation
  lines(density(sim_cent_samples[, node_id_check]),
        col = rgb(0,0,1,0.5), lwd = 2)
}
par(mfrow = c(1,1))

## clean up and save workspace
rm(mean_df, df_plot, nodes_timewindow, node_id_check, node_id_sample, time_window, mu, cov, sim_cent_samples) ; gc()
save.image('anp_nodalregression/anp_long_nodal_modelprep.RData')
print('normal approximation complete')

#### create data list #####
#load('anp_nodalregression/anp_long_nodal_modelprep.RData')
n_data <- nrow(nodes_all)
n_nodes <- length(unique(nodes_all$node_random))
n_windows <- length(unique(nodes_all$window))

model_data_list <- list(
  # global data size
  num_data = n_data,
  num_nodes = n_nodes,
  num_windows = n_windows,
  # per time window data size
  num_nodes_window1 = length(unique(nodes_all$node_random[nodes_all$window == 1])),
  num_nodes_window2 = length(unique(nodes_all$node_random[nodes_all$window == 2])),
  num_nodes_window3 = length(unique(nodes_all$node_random[nodes_all$window == 3])),
  num_nodes_window4 = length(unique(nodes_all$node_random[nodes_all$window == 4])),
  num_nodes_window5 = length(unique(nodes_all$node_random[nodes_all$window == 5])),
  num_nodes_window6 = length(unique(nodes_all$node_random[nodes_all$window == 6])),
  num_nodes_window7 = length(unique(nodes_all$node_random[nodes_all$window == 7])),
  # number of nodes in all preceding time windows for node age indexing
  num_nodes_prev_windows = c(0,
                             length(which(nodes_all$window < 2)),
                             length(which(nodes_all$window < 3)),
                             length(which(nodes_all$window < 4)),
                             length(which(nodes_all$window < 5)),
                             length(which(nodes_all$window < 6)),
                             length(which(nodes_all$window < 7))),
  # centrality means per time window
  centrality_mu_1 = nodes_all$mean_eigen[nodes_all$window == 1],
  centrality_mu_2 = nodes_all$mean_eigen[nodes_all$window == 2],
  centrality_mu_3 = nodes_all$mean_eigen[nodes_all$window == 3],
  centrality_mu_4 = nodes_all$mean_eigen[nodes_all$window == 4],
  centrality_mu_5 = nodes_all$mean_eigen[nodes_all$window == 5],
  centrality_mu_6 = nodes_all$mean_eigen[nodes_all$window == 6],
  centrality_mu_7 = nodes_all$mean_eigen[nodes_all$window == 7],
  # covariance matrix per time window
  centrality_cov_1 = covs_all[[1]],
  centrality_cov_2 = covs_all[[2]],
  centrality_cov_3 = covs_all[[3]],
  centrality_cov_4 = covs_all[[4]],
  centrality_cov_5 = covs_all[[5]],
  centrality_cov_6 = covs_all[[6]],
  centrality_cov_7 = covs_all[[7]],
  # node IDs for all time windows
  nodes_window1 = nodes_all$node_random[nodes_all$window == 1],
  nodes_window2 = nodes_all$node_random[nodes_all$window == 2],
  nodes_window3 = nodes_all$node_random[nodes_all$window == 3],
  nodes_window4 = nodes_all$node_random[nodes_all$window == 4],
  nodes_window5 = nodes_all$node_random[nodes_all$window == 5],
  nodes_window6 = nodes_all$node_random[nodes_all$window == 6],
  nodes_window7 = nodes_all$node_random[nodes_all$window == 7],
  nodes_window8 = nodes_all$node_random[nodes_all$window == 8],
  # exposure variable
  node_age = nodes_all$age_std)

## check inputs
colour <- c('red','orange','yellow','green','blue','purple','magenta')
plot(model_data_list$centrality_mu_1 ~ model_data_list$node_age[1:model_data_list$num_nodes_window1],
     pch = 19, col = colour[1],
     ylim = c(-10,0), xlim = c(-2.5,5))
for(i in 2:n_windows){
  if(i < n_windows){
    points(model_data_list[[(i+n_windows+4)]] ~ model_data_list$node_age[(model_data_list$num_nodes_prev_windows[i]+1):model_data_list$num_nodes_prev_windows[i+1]],
           col = colour[i], pch = 19)
  } else {
    points(model_data_list[[(i+n_windows+4)]] ~ model_data_list$node_age[(model_data_list$num_nodes_prev_windows[i]+1):length(model_data_list$node_age)],
           col = colour[i], pch = 19)
  }
}

print('data created')

#### prior predictive check ####
n <- 100
beta_age <- rnorm(n, 0, 0.8)
min_raw <- min(unlist(model_data_list[grep('centrality_mu', names(model_data_list))]))
max_raw <- max(unlist(model_data_list[grep('centrality_mu', names(model_data_list))]))
intercept  <- rnorm(n, LaplacesDemon::logit(0.05), 3) # taking the intercept from the logit of results from Chiyo 2011 (doesn't state mean/median centrality so estimated from graph based on where correlation line would cross x = 0)
plot(NULL, las = 1, xlab = 'age (standardised)', ylab = 'eigenvector',
     #ylim = c(min_raw-2, max_raw+2),
     ylim = c(-15, 5),
     xlim = c(min(model_data_list$node_age),
              max(model_data_list$node_age)))
abline(h = min_raw, lty = 2) ; abline(h = max_raw, lty = 2)
for(i in 1:n){
  lines(x = seq(min(nodes_all$age_std), max(nodes_all$age_std), length.out = 2),
        y = intercept[i] + beta_age[i]*c(min(nodes_all$age_std),
                                         max(nodes_all$age_std)),
        col = rgb(0,0,1,0.4))
}
rm(n, max_raw, min_raw, beta_age, intercept) ; gc()
dev.off()
print('prior predictive completed, model ready to run')

#### run model ####
## set model parameters
n_chains <- 4
n_samples <- 1000

## load model
nodal_regression <- cmdstan_model('models/eigen_regression_anplong.stan',
                                  cpp_options = list(stan_threads = TRUE))

## run model
fit_anp_nodal <- nodal_regression$sample(data = model_data_list,
                                         chains = n_chains, parallel_chains = n_chains,
                                         threads_per_chain = 4,
                                         iter_warmup = n_samples, iter_sampling = n_samples)

## save workspace
save.image('anp_nodalregression/anp_long_nodal.RData')
print('model run')

#### check outputs ####
#load('anp_nodalregression/anp_long_nodal.RData')

## define PDF output
pdf('../outputs/step4_nodalregression/anplong_nodalregression_modelchecks.pdf')

## extract model fit -- very good!
( summary <- fit_anp_nodal$summary() )
par(mfrow = c(3,1))
hist(summary$rhat, breaks = 50)
hist(summary$ess_bulk, breaks = 50)
hist(summary$ess_tail, breaks = 50)
par(mfrow = c(1,1))

## extract posterior
#params <- rstan::extract(fit_anp_nodal)
params <- fit_anp_nodal$draws(format = 'draws_df')

## separate random effects from global parameters
rand_window <- params %>%
  dplyr::select(grep('window_random_effect', colnames(params), value=TRUE))
rand_node <- params %>%
  dplyr::select(grep('node_random_effect', colnames(params), value=TRUE))

## traceplot all parameters
#traceplot(fit_anp_nodal, pars = c('intercept','beta_age','sigma','predictor[1]','predictor[50]','predictor[100]'))
plot_params <- c('intercept','beta_age','sigma') #,'nu')
plot_rand_nodes <- colnames(params)[sample(grep('node_random_effect',colnames(params)), size = 12, replace = F)]
plot_rand_windows <- colnames(params)[grep('window_random_effect',colnames(params))]

traceplot(parameter_df = params, parameters_to_plot = plot_params)
traceplot(parameter_df = params, parameters_to_plot = plot_rand_nodes)
traceplot(parameter_df = params, parameters_to_plot = plot_rand_windows)
rm(plot_params,plot_rand_nodes,plot_rand_windows) ; gc()

## save workspace
save.image('anp_nodalregression/anp_long_nodal.RData')
print('outputs checked')

#### posterior predictive check ####
par(mfrow = c(3,3))
for(time_window in 1:n_windows){
  post_pred_check(centrality_matrix = cents_all,
                  nodes_df = nodes_all, time_window = time_window,
                  cent_cov = covs_all[[time_window]],
                  parameters = params)
}
par(mfrow = c(1,1))

## define PDF output
dev.off
pdf('../outputs/step4_nodalregression/anplong_nodalregression_modelpredictions.pdf')
print('posterior predictive complete')

#### predict from model -- standardised scale ####
## get mean predictions
mu_predictions <- get_mean_predictions(predict_df = nodes_all, parameters = params,
                                       include_window = TRUE, include_node = TRUE)

## add mean and CI of predicted means to input data frame for comparison
nodes_all$mu_mean <- apply(mu_predictions, 2, mean)
nodes_all$mu_upr <- apply(mu_predictions, 2, quantile, prob = 0.975)
nodes_all$mu_lwr <- apply(mu_predictions, 2, quantile, prob = 0.025)

## plot mean of model vs mean of raw data
compare_plot <- ggplot()+
  geom_point(data = nodes_all, #size = 0.5,
             mapping = aes(x = mean_eigen, y = mu_mean, colour = age))+
  scale_colour_viridis_c()+
  facet_wrap(facets = . ~ as.factor(window))+
  labs(colour = 'window',
       x = 'raw eigen mean',
       y = 'predicted mean')+
  geom_abline(slope = 1, intercept = 0) # add line showing where points would lie if model fit was perfect
compare_plot +
  scale_x_continuous(limits = c(-7,-1))+
  scale_y_continuous(limits = c(-7,-1))
compare_plot +
  facet_wrap(facets = . ~ as.factor(window), scales = 'free')

## put together sigma arrays, separated by time window
num_nodes_windows <- unlist(model_data_list[grep('num_nodes_window', names(model_data_list))])
sigma_list <- list()
for(time_window in 1:n_windows) {
  sigma_window <- array(NA, dim = c(num_nodes_windows[time_window],
                                    num_nodes_windows[time_window],
                                    nrow(params)),
                        dimnames = list(nodes_all$nodes_random[nodes_all$window == time_window],
                                        nodes_all$nodes_random[nodes_all$window == time_window],
                                        NULL))
  cent_cov <- covs_all[[time_window]]
  for(i in 1:nrow(params)){
    sigma_window[,,i] <- cent_cov + diag(rep(params$sigma[i], num_nodes_windows[time_window]))
  }
  sigma_list[[time_window]] <- sigma_window
}

## create empty matrix to take full set of predicted values per elephant
full_predictions <- matrix(NA, nrow = nrow(params), ncol = nrow(nodes_all),
                           dimnames = list(NULL, nodes_all$node_window))

## populate matrix using mean values in matrix mu_std, and sigma values based on time window
for(time_window in 1:n_windows){
  sigma_window <- sigma_list[[time_window]]
  columns <- which(nodes_all$window == time_window)
  for(i in 1:nrow(full_predictions)){
    full_predictions[i,columns] <- MASS::mvrnorm(1,
                                                 mu_predictions[i,columns],
                                                 sigma_window[,,i])
  }
}

## add CI of predicted data points to input data frame for comparison
nodes_all$predict_lwr_std <- apply(full_predictions, 2, quantile, prob = 0.025)
nodes_all$predict_upr_std <- apply(full_predictions, 2, quantile, prob = 0.975)

## plot predictions
ggplot(nodes_all)+
  geom_ribbon(aes(x = age, ymin = predict_lwr_std, ymax = predict_upr_std, fill = as.factor(window)),
              alpha = 0.2)+                      # background layer showing the 95% CI of all predictions
  geom_ribbon(aes(x = age, ymin = mu_lwr, ymax = mu_upr, fill = as.factor(window)),
              alpha = 0.4)+                      # mid layer showing the 95% CI of predicted means
  geom_line(aes(x = age, y = mu_mean, colour = as.factor(window)))+  # line showing mean of predicted means
  geom_point(aes(x = age, y = mean_eigen), size = 0.5)+              # original data points (standardised centrality, actual age)
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  scale_x_continuous(expand = c(0,0))+
  facet_wrap(. ~ as.factor(window),
             scales = 'free')+             # separate plots per window
  theme(legend.position = 'none')+
  #guides(colour = guide_legend(nrow = 6),
  #       fill = guide_legend(nrow = 6))+
  labs(#colour = 'time window', fill = 'time window',
    y = 'eigenvector centrality', x = 'age (years)')

## save workspace
save.image('anp_nodalregression/anp_long_nodal.RData')
print('predictions made')

#### extract original values from output -- simulated slope value originally used produces an effect on the unstandardised scale. The model works on the standardised scale. Convert predictions to unstandardised scale and then run contrasts to calculate the slope coefficient. ####
## get mean predictions
newdata <- nodes_all %>%
  mutate(age_std = age_std+1)
mu_predictions_newdata <- get_mean_predictions(predict_df = newdata, parameters = params,
                                               include_window = TRUE, include_node = TRUE)
nodes_all$mu_mean_plus1 <- apply(mu_predictions_newdata, 2, mean)

## full distribution of predictions using age_std + 1 stdev
full_predictions_newdata <- matrix(NA, nrow = nrow(params), ncol = nrow(newdata),
                                   dimnames = list(NULL, nodes_all$node_window))
for(time_window in 1:n_windows){
  sigma_window <- sigma_list[[time_window]]
  columns <- which(nodes_all$window == time_window)
  for(i in 1:nrow(full_predictions_newdata)){
    full_predictions_newdata[i,columns] <- MASS::mvrnorm(1,
                                                         mu_predictions_newdata[i,columns],
                                                         sigma_window[,,i])
  }
}

## contrast predictions on standardised scale -- check that this returns the marginal effect presented in the summary
contrast_std <- full_predictions_newdata - full_predictions    # contrast between predicted values for raw data and all same data but add 1 to age
head(contrast_std[,1:5])                              # check matrix looks right
mean(params$beta_age)                                 # output parameter direct from model
mean(contrast_std)                                    # should be very similar to mean of output parameter
quantile(contrast_std, prob = c(0.025, 0.975))        # very wide

## contrast predictions on output scale -- reportable values of effect of plus 1 year
contrast <- contrast_std / sd(nodes_all$age)          # convert to outcome scale
head(contrast[,1:5])                                  # check matrix looks right
mean(contrast)                                        # should be very similar to mean of output parameter
quantile(contrast, prob = c(0.025, 0.975))            # very wide

## save workspace
save.image('anp_nodalregression/anp_long_nodal.RData')
print('contrasts complete')

#### final "clean" plots using hypothetical data ####
#load('anp_nodalregression/anp_long_nodal.RData')

## create dummy dataset
newdat <- nodes_all %>%
  dplyr::select(node_random, age, age_std, window)

## get mean predictions
fake_mean <- get_mean_predictions(predict_df = newdat, parameters = params,
                                  include_window = TRUE, include_node = FALSE)
newdat$predict_mean <- apply(fake_mean, 2, mean)
newdat$predict_mean_lwr <- apply(fake_mean, 2, quantile, prob = 0.025)
newdat$predict_mean_upr <- apply(fake_mean, 2, quantile, prob = 0.975)

## create prediction matrix
fake_pred <- matrix(NA, nrow = nrow(params), ncol = nrow(newdat))
for(time_window in 1:n_windows){
  sigma_window <- sigma_list[[time_window]]
  columns <- which(nodes_all$window == time_window)
  for(i in 1:nrow(full_predictions_newdata)){
    fake_pred[i,columns] <- MASS::mvrnorm(1,
                                          fake_mean[i,columns],
                                          sigma_window[,,i])
  }
}
print('predictions for clean plots made')

## add CI of predicted data points to input data frame for comparison
newdat$predict_pred_lwr <- apply(fake_pred, 2, quantile, prob = 0.025)
newdat$predict_pred_upr <- apply(fake_pred, 2, quantile, prob = 0.975)

## plot predictions
newdat_summary <- newdat %>%
  group_by(age, window) %>%
  mutate(predict_pred_lwr = mean(predict_pred_lwr),
         predict_pred_upr = mean(predict_pred_upr)) %>%
  dplyr::select(age, predict_pred_lwr, predict_pred_upr,window) %>%
  distinct()

## plot mean values
ggplot()+
  geom_ribbon(data = newdat_summary,
              aes(x = age, ymin = predict_pred_lwr, ymax = predict_pred_upr, fill = as.factor(window)),
              alpha = 0.2)+                      # background layer showing the 95% CI of all predictions
  geom_ribbon(data = newdat,
              aes(x = age, ymin = predict_mean_lwr, ymax = predict_mean_upr, fill = as.factor(window)),
              alpha = 0.4)+                      # mid layer showing the 95% CI of predicted means
  geom_line(data = newdat,
            aes(x = age, y = predict_mean, colour = as.factor(window)))+  # line showing mean of predicted means
  geom_point(data = nodes_all,
             aes(x = age, y = mean_eigen))+      # original data points (standardised centrality, actual age)
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  facet_wrap(. ~ as.factor(window))+             # separate plots per window
  theme(legend.position = 'bottom')+
  labs(colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'age (years)')
print('mean values plotted')

## convert raw data to data frame
cents_all_df <- cents_all %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = 'node_window', values_to = 'eigen') %>%
  left_join(nodes_all[,c('node_window','age','window')], by = 'node_window')

## plot full distribution
ggplot()+
  geom_ribbon(data = newdat_summary,
              aes(x = age, ymin = predict_pred_lwr, ymax = predict_pred_upr, fill = as.factor(window)),
              alpha = 0.2)+                      # background layer showing the 95% CI of all predictions
  geom_ribbon(data = newdat,
              aes(x = age, ymin = predict_mean_lwr, ymax = predict_mean_upr, fill = as.factor(window)),
              alpha = 0.4)+                      # mid layer showing the 95% CI of predicted means
  geom_point(data = cents_all_df,
             aes(x = age, y = eigen),
             alpha = 0.01, size = 0.5)+         # original data points (standardised centrality, actual age)
  geom_line(data = newdat,
            aes(x = age, y = predict_mean, colour = as.factor(window)))+  # line showing mean of predicted means
  geom_point(data = nodes_all,
             aes(x = age, y = mean_eigen),
             size = 0.5, colour = 'white')+      # original mean data points (standardised centrality, actual age)
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  facet_wrap(. ~ as.factor(window))+             # separate plots per window
  theme(legend.position = 'bottom')+
  labs(colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'age (years)')

## clean up and save workspace
save.image('anp_nodalregression/anp_long_nodal.RData')
dev.off()
print('long windows complete')

#### smooth plot -- I think this should replace the code chunk above, but check before deleting ####
# load('anp_nodalregression/anp_long_nodal.RData')
rm(cent_cov, compare_plot, contrast, contrast_std, fake_mean, fake_pred, fit_anp_nodal, full_predictions, full_predictions_newdata, model_data_list, mu_predictions, mu_predictions_newdata, newdata, nodal_regression, rand_node, rand_window, summary, colour) ; gc()

## create age categorical variables to match MOTNP
nodes_all <- nodes_all %>% 
  mutate(age_cat = ifelse(age < 15, 1,
                          ifelse(age < 20, 2,
                                 ifelse(age < 25, 3,
                                        ifelse(age < 40, 4, 5)))),
         age_cat_full = ifelse(age < 15, '10-15 yrs',
                               ifelse(age < 20, '16-20 yrs',
                                      ifelse(age < 25, '21-25 yrs',
                                             ifelse(age < 40,
                                                    '25-40 yrs',
                                                    '>40 yrs')))))

## predict from raw data but exclude window and node effects (except for in sigma_window as I think that's unavoidable)
set.seed(12345) ; draws <- sample(x = 1:n_samples, size = 125, replace = F)
pred_mean_norandom <- get_mean_predictions(predict_df = nodes_all,
                                           parameters = params[c(draws,
                                                                 draws+1000,
                                                                 draws+2000,
                                                                 draws+3000),],
                                           include_window = FALSE, include_node = FALSE)

## create data frame for smooth predictions
norandom_summary <- nodes_all %>%
  dplyr::select(age, age_std, age_cat) %>%
  distinct() %>%
  mutate(predict_mean_logit = NA, predict_mean_invlogit = NA,
         predict_sd_logit = NA, predict_sd_invlogit = NA,
         predict_mean_lwr_logit = NA, predict_mean_lwr_invlogit = NA,
         predict_mean_upr_logit = NA, predict_mean_upr_invlogit = NA,
         predict_full_lwr_logit = NA, predict_full_lwr_invlogit = NA,
         predict_full_upr_logit = NA, predict_full_upr_invlogit = NA)
for(i in 1:nrow(norandom_summary)){
  ## select only columns that refer to elephants of that age
  mean_predictions_per_age <- pred_mean_norandom[,which(nodes_all$age == norandom_summary$age[i])]
  
  ## calculate values for mean line
  norandom_summary$predict_mean_logit[i] <- mean(mean_predictions_per_age)
  norandom_summary$predict_mean_invlogit[i] <- mean(invlogit(mean_predictions_per_age))
  
  ## calculate standard deviation
  norandom_summary$predict_sd_logit[i] <- sd(mean_predictions_per_age)
  norandom_summary$predict_sd_invlogit[i] <- sd(invlogit(mean_predictions_per_age))
  
  ## calculate values for shaded ribbon
  norandom_summary$predict_mean_lwr_logit[i] <- quantile(mean_predictions_per_age,
                                                         prob = 0.025)
  norandom_summary$predict_mean_lwr_invlogit[i] <- quantile(invlogit(mean_predictions_per_age),
                                                            prob = 0.025)
  norandom_summary$predict_mean_upr_logit[i] <- quantile(mean_predictions_per_age,
                                                         prob = 0.975)
  norandom_summary$predict_mean_upr_invlogit[i] <- quantile(invlogit(mean_predictions_per_age),
                                                            prob = 0.975)
}

## create empty matrix to take full set of predicted values per elephant
pred_full_norandom <- matrix(NA, nrow = nrow(pred_mean_norandom), ncol = nrow(nodes_all),
                             dimnames = list(NULL, nodes_all$node_window))

## populate matrix using mean values in matrix mu_std, and sigma values based on time window
for(time_window in 1:n_windows){
  sigma_window <- sigma_list[[time_window]]
  columns <- which(nodes_all$window == time_window)
  for(i in 1:nrow(pred_full_norandom)){
    pred_full_norandom[i,columns] <- MASS::mvrnorm(1,
                                                   pred_mean_norandom[i,columns],
                                                   sigma_window[,,i])
  }
}

## add CI of predicted data points to input data frame for comparison
for(i in 1:nrow(norandom_summary)){
  ## select only columns that refer to elephants of that age
  mean_predictions_per_age <- pred_full_norandom[,which(nodes_all$age == norandom_summary$age[i])]
  
  ## calculate values for shaded ribbon
  norandom_summary$predict_full_lwr_logit[i] <- quantile(mean_predictions_per_age, 
                                                         prob = 0.025)
  norandom_summary$predict_full_lwr_invlogit[i] <- quantile(invlogit(mean_predictions_per_age),
                                                            prob = 0.025)
  norandom_summary$predict_full_upr_logit[i] <- quantile(mean_predictions_per_age,
                                                         prob = 0.975)
  norandom_summary$predict_full_upr_invlogit[i] <- quantile(invlogit(mean_predictions_per_age),
                                                            prob = 0.975)
}

# ## add full distribution (probably unnecessary and becomes too big for computer to open)
# cents_all_df <- cents_all %>%
#   as.data.frame() %>%
#   pivot_longer(cols = everything(), names_to = 'node_window', values_to = 'eigen') %>%
#   left_join(nodes_all[,c('node_window','age','window')], by = 'node_window')

## faceted by window
(faceted <- ggplot()+
    geom_ribbon(data = newdat_summary,
                aes(ymin = invlogit(predict_pred_lwr),
                    ymax = invlogit(predict_pred_upr),
                    group = as.factor(window),#fill = as.factor(window),
                    x = age),
                alpha = 0.2)+                      # 95% CI of all predictions
    geom_ribbon(data = newdat,
                aes(ymin = invlogit(predict_mean_lwr),
                    ymax = invlogit(predict_mean_upr),
                    group = as.factor(window),#fill = as.factor(window),
                    x = age),
                alpha = 0.4)+                      # 95% CI of predicted means
    # geom_point(data = cents_all_df,
    #            aes(x = age, y = eigen, colour = as.factor(window)),
    #            alpha = 0.01, size = 0.5)+         # original data points
    geom_point(data = distinct(nodes_all[,c('age','mean_eigen','window','age_cat_full','sightings')]),
               #size = 0.5, colour = 'white',
               aes(x = age,
                   y = invlogit(mean_eigen),
                   colour = factor(age_cat_full,
                                   levels = c('10-15 yrs','16-20 yrs',
                                              '21-25 yrs', '25-40 yrs',
                                              '>40 yrs')),
                   #group = as.factor(window),
                   size = sightings
               ),
               alpha = 0.5)+      # original mean data points
    geom_line(data = newdat,
              aes(x = age,
                  y = invlogit(predict_mean),
                  group = as.factor(window)))+  # predicted means
    scale_colour_viridis_d(direction = -1)+ # begin = 0, end = 0.7
    theme(legend.position = 'bottom',
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 18))+
    facet_wrap(. ~ as.factor(window)) +
    guides(#colour = guide_legend(nrow = 5),
           colour = 'none',
           size = guide_legend(override.aes = list(fill = 'grey80')))+
    labs(colour = '  age category\n(as in MOTNP)',
         y = 'eigenvector centrality', x = 'age (years)'))
ggsave(plot = faceted, device = 'svg', width = 2400, height = 2400, units = 'px',
       filename = 'anplong_nodalregression_faceted_invlogit.svg',
       path = '../outputs/step4_nodalregression/')
ggsave(plot = faceted, device = 'png', width = 2400, height = 2400, units = 'px',
       filename = 'anplong_nodalregression_faceted_invlogit.png',
       path = '../outputs/step4_nodalregression/')

## overall result
(overall <- ggplot()+
    geom_ribbon(data = norandom_summary,
                aes(x = age,
                    ymin = predict_full_lwr_invlogit,
                    ymax = predict_full_upr_invlogit),
                alpha = 0.2)+
    # geom_point(data = cents_all_df,
    #            aes(x = age, y = eigen, colour = as.factor(window)),
    #            alpha = 0.01, size = 0.5)+         # original data points
    geom_point(data = distinct(nodes_all[,c('age','mean_eigen','age_cat_full','sightings')]),
               aes(x = age, y = invlogit(mean_eigen),
                   colour = factor(age_cat_full,
                                   levels = c('10-15 yrs','16-20 yrs',
                                              '21-25 yrs', '25-40 yrs',
                                              '>40 yrs')),
                   size = sightings
               ),
               alpha = 0.5)+
    geom_ribbon(data = norandom_summary,
                aes(x = age,
                    ymin = predict_mean_lwr_invlogit,
                    ymax = predict_mean_upr_invlogit),
                alpha = 0.4)+
    geom_line(data = norandom_summary,
              aes(x = age,
                  y = predict_mean_invlogit))+
    scale_colour_viridis_d(direction = -1)+ # begin = 0, end = 0.7
    theme(legend.position = 'bottom',
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 18))+
    guides(#colour = guide_legend(nrow = 5),
      colour = 'none',
      size = guide_legend(override.aes = list(fill = 'grey80')))+
    labs(colour = '  age category\n(as in MOTNP)',
         y = 'eigenvector centrality', x = 'age (years)'))
ggsave(plot = overall, device = 'svg', width = 2400, height = 1800, units = 'px',
       filename = 'anplong_nodalregression_nofaceting.svg',
       path = '../outputs/step4_nodalregression/')
ggsave(plot = overall, device = 'png', width = 2400, height = 1800, units = 'px',
       filename = 'anplong_nodalregression_nofaceting.png',
       path = '../outputs/step4_nodalregression/')

## save outputs
dev.off()

#################### RANDOM EFFECTS ####################
#### LONG: calculate proportion of variation explained by window random effect ####
## calculate values
mean_window_logit_long <- apply(rand_window, 2, mean)
stdv_window_logit_long <- apply(rand_window, 2, sd)

## convert to data frame
mean_window_logit_long <- mean_window_logit_long %>% 
  as.data.frame() %>% 
  mutate(eles_per_window = as.vector(unlist(model_data_list[grep(pattern = 'num_nodes_window',
                                                                 x = names(model_data_list))]))) %>% 
  mutate(window = 1:7) %>% 
  relocate(window) %>% 
  mutate(stdv_effect = as.vector(stdv_window_logit_long))
colnames(mean_window_logit_long)[2] <- 'mean_effect'

## plot
mean_window_logit_long %>% 
  ggplot()+
  geom_errorbar(aes(x = eles_per_window,
                    ymax = mean_effect + stdv_effect,
                    ymin = mean_effect - stdv_effect,
                    colour = as.factor(window)))+
  scale_colour_viridis_d()+
  geom_point(aes(x = eles_per_window,
                 #colour = as.factor(window),
                 y = mean_effect),
             pch = 21,
             colour = 'black',
             fill = 'white')+
  geom_smooth(aes(x = eles_per_window,
                  y = mean_effect),
              colour = 'black')+
  labs(x = 'number of elephants in time window',
       y = 'mean effect of window ID (logit scale)')+
  guides(colour = guide_legend(ncol = 2, nrow = 18,
                               title = 'window ID',
                               keyheight = 1))+
  theme(legend.position = 'right',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 14))
ggsave(plot = last_plot(), device = 'png',
       filename = 'anplong_window_random_effect.png',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')
ggsave(plot = last_plot(), device = 'svg',
       filename = 'anplong_window_random_effect.svg',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')

#### LONG: calculate proportion of variation explained by node random effect ####
## calculate values
mean_node_logit_long <- apply(rand_node, 2, mean)
stdv_node_logit_long <- apply(rand_node, 2, sd)

## convert to data frame
counts_long <- nodes_all %>% 
  mutate(observed = ifelse(sightings == 0, 0, 1)) %>% 
  group_by(node) %>% 
  summarise(total_sightings = sum(sightings),
            total_windows = sum(observed))
mean_node_logit_long <- mean_node_logit_long %>% 
  as.data.frame() %>% 
  mutate(node = sort(unique(nodes_all$node))) %>% 
  relocate(node) %>% 
  left_join(counts_long, by = 'node') %>% 
  mutate(stdv_effect = as.vector(stdv_node_logit_long))
colnames(mean_node_logit_long)[2] <- 'mean_effect'

## plot
mean_node_logit_long %>% 
  ggplot()+
  geom_errorbar(aes(x = total_sightings,#total_windows,
                    colour = node,
                    ymax = mean_effect + stdv_effect,
                    ymin = mean_effect - stdv_effect))+
  scale_colour_viridis_c(breaks = c(1,seq(25,700, by = 25)),
                         name = 'node ID')+
  geom_point(aes(x = total_sightings,#total_windows,
                 #colour = node,
                 y = mean_effect),
             pch = 21,
             colour = 'black',
             fill = 'white')+
  geom_smooth(aes(x = total_sightings,#total_windows,
                  y = mean_effect),
              colour = 'black')+
  labs(x = 'total number observations per individual',#'number of windows in which individual was observed',
       y = 'mean effect of node ID (logit scale)')+
  guides(colour = guide_legend(ncol = 2,
                               keyheight = 1))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))
ggsave(plot = last_plot(), device = 'png',
       filename = 'anplong_node_random_effect.png',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')
ggsave(plot = last_plot(), device = 'svg',
       filename = 'anplong_node_random_effect.svg',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')

## plot both
random_effects_long <- mean_node_logit_long %>% 
  select(-node, -total_sightings, -total_windows) %>% 
  mutate(type = 'node')
window_effect_long <- mean_window_logit_long %>% 
  select(-window, -eles_per_window) %>% 
  mutate(type = 'window')
random_effects_long <- rbind(random_effects_long,window_effect_long)
rm(window_effect_long)
ggplot()+
  geom_violin(data = random_effects_long,
              aes(x = type,
                  y = mean_effect))+
  geom_jitter(data = random_effects_long,
              aes(x = type,
                  y = mean_effect))

#### SHORT: calculate proportion of variation explained by window random effect ####
rm(list = ls()[!ls() %in% c('mean_node_logit_long','mean_window_logit_long')])
mean_node_logit_long <- mean_node_logit_long %>% 
  mutate(window_length = 'long')
mean_window_logit_long <- mean_window_logit_long %>% 
  mutate(window_length = 'long')
load('anp_nodalregression/anp_short_nodal.RData')

## calculate values
mean_window_logit_short <- apply(rand_window, 2, mean)
stdv_window_logit_short <- apply(rand_window, 2, sd)

## convert to data frame
mean_window_logit_short <- mean_window_logit_short %>% 
  as.data.frame() %>% 
  mutate(eles_per_window = as.vector(unlist(model_data_list[grep(pattern = 'num_nodes_window',
                                                                 x = names(model_data_list))]))) %>% 
  mutate(window = 1:36) %>% 
  relocate(window) %>% 
  mutate(stdv_effect = as.vector(stdv_window_logit_short))
colnames(mean_window_logit_short)[2] <- 'mean_effect'

## plot
mean_window_logit_short %>% 
  ggplot()+
  geom_errorbar(aes(x = eles_per_window,
                    ymax = mean_effect + stdv_effect,
                    ymin = mean_effect - stdv_effect,
                    colour = as.factor(window)))+
  scale_colour_viridis_d()+
  geom_point(aes(x = eles_per_window,
                 #colour = as.factor(window),
                 y = mean_effect),
             pch = 21,
             colour = 'black',
             fill = 'white')+
  geom_smooth(aes(x = eles_per_window,
                  y = mean_effect),
              colour = 'black')+
  labs(x = 'number of elephants in time window',
       y = 'mean effect of window ID (logit scale)')+
  guides(colour = guide_legend(ncol = 2, nrow = 18,
                               title = 'window ID',
                               keyheight = 1))+
  theme(legend.position = 'right',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 14))
ggsave(plot = last_plot(), device = 'png',
       filename = 'anpshort_window_random_effect.png',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')
ggsave(plot = last_plot(), device = 'svg',
       filename = 'anpshort_window_random_effect.svg',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')

#### SHORT: calculate proportion of variation explained by node random effect ####
## calculate values
mean_node_logit_short <- apply(rand_node, 2, mean)
stdv_node_logit_short <- apply(rand_node, 2, sd)

## convert to data frame
counts_short <- nodes_all %>% 
  mutate(observed = ifelse(sightings == 0, 0, 1)) %>% 
  group_by(node) %>% 
  summarise(total_sightings = sum(sightings),
            total_windows = sum(observed))
mean_node_logit_short <- mean_node_logit_short %>% 
  as.data.frame() %>% 
  mutate(node = sort(unique(nodes_all$node))) %>% 
  relocate(node) %>% 
  left_join(counts_short, by = 'node') %>% 
  mutate(stdv_effect = as.vector(stdv_node_logit_short))
colnames(mean_node_logit_short)[2] <- 'mean_effect'

## plot
mean_node_logit_short %>% 
  ggplot()+
  geom_errorbar(aes(x = total_sightings,#total_windows,
                    colour = node,
                    ymax = mean_effect + stdv_effect,
                    ymin = mean_effect - stdv_effect))+
  scale_colour_viridis_c(breaks = c(1,seq(25,700, by = 25)),
                         name = 'node ID')+
  geom_point(aes(x = total_sightings,#total_windows,
                 #colour = node,
                 y = mean_effect),
             pch = 21,
             colour = 'black',
             fill = 'white')+
  geom_smooth(aes(x = total_sightings,#total_windows,
                  y = mean_effect),
              colour = 'black')+
  labs(x = 'total number observations per individual',#'number of windows in which individual was observed',
       y = 'mean effect of node ID (logit scale)')+
  guides(colour = guide_legend(ncol = 2,
                               keyheight = 1))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))
ggsave(plot = last_plot(), device = 'png',
       filename = 'anpshort_node_random_effect.png',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')
ggsave(plot = last_plot(), device = 'svg',
       filename = 'anpshort_node_random_effect.svg',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')

## plot both
random_effects_short <- mean_node_logit_short %>% 
  select(-node, -total_sightings, -total_windows) %>% 
  mutate(type = 'node')
window_effect_short <- mean_window_logit_short %>% 
  select(-window, -eles_per_window) %>% 
  mutate(type = 'window')
random_effects_short <- rbind(random_effects_short,window_effect_short)
rm(window_effect_short)
ggplot()+
  geom_violin(data = random_effects_short,
              aes(x = type,
                  y = mean_effect))+
  geom_jitter(data = random_effects_short,
              aes(x = type,
                  y = mean_effect))

#### plot together ####
rm(list = ls()[!ls() %in% c('mean_node_logit_long','mean_window_logit_long',
                            'mean_node_logit_short','mean_window_logit_short')])
mean_node_logit_short <- mean_node_logit_short %>% 
  mutate(window_length = 'short')
mean_window_logit_short <- mean_window_logit_short %>% 
  mutate(window_length = 'short')
save.image('anp_nodalregression/plot_random_effects.RData')

## combine data frames
node <- rbind(mean_node_logit_long, mean_node_logit_short)
window <- rbind(mean_window_logit_long, mean_window_logit_short)

## plot
node %>% 
  ggplot()+
  geom_errorbar(aes(x = total_sightings,#total_windows,
                    colour = node,
                    ymax = mean_effect + stdv_effect,
                    ymin = mean_effect - stdv_effect))+
  scale_colour_viridis_c(breaks = c(1,seq(25,700, by = 25)),
                         name = 'node ID')+
  geom_point(aes(x = total_sightings,#total_windows,
                 #colour = node,
                 y = mean_effect),
             pch = 21,
             colour = 'black',
             fill = 'white')+
  geom_smooth(aes(x = total_sightings,#total_windows,
                  y = mean_effect),
              colour = 'black')+
  labs(x = 'total number observations per individual',#'number of windows in which individual was observed',
       y = 'mean effect of node ID (logit scale)')+
  guides(colour = guide_legend(ncol = 2,
                               keyheight = 1))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  facet_wrap(. ~ factor(window_length, levels = c('short','long')))
ggsave(plot = last_plot(), device = 'png',
       filename = 'anp_node_random_effect_longshort.png',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')
ggsave(plot = last_plot(), device = 'svg',
       filename = 'anp_node_random_effect_longshort.svg',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')

## plot
window %>% 
  ggplot()+
  geom_errorbar(aes(x = eles_per_window,
                    ymax = mean_effect + stdv_effect,
                    ymin = mean_effect - stdv_effect,
                    colour = as.factor(window)))+
  scale_colour_viridis_d()+
  geom_point(aes(x = eles_per_window,
                 #colour = as.factor(window),
                 y = mean_effect),
             pch = 21,
             colour = 'black',
             fill = 'white')+
  geom_smooth(aes(x = eles_per_window,
                  y = mean_effect),
              colour = 'black')+
  labs(x = 'number of elephants in time window',
       y = 'mean effect of window ID (logit scale)')+
  guides(colour = guide_legend(ncol = 2, nrow = 18,
                               title = 'window ID',
                               keyheight = 1))+
  theme(legend.position = 'right',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 14))+
  facet_wrap(. ~ factor(window_length,
                        levels = c('short','long')),
             scales = 'free_x')
ggsave(plot = last_plot(), device = 'png',
       filename = 'anp_window_random_effect_longshort.png',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')
ggsave(plot = last_plot(), device = 'svg',
       filename = 'anp_window_random_effect_longshort.svg',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')
