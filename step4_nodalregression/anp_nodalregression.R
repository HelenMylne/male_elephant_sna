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
set_cmdstan_path('../packages/.cmdstan/cmdstan-2.34.1/') # Viking

## set seed for reproducibility
set.seed(12345)

## set theme for plots
theme_set(theme_bw())

## define PDF output
#pdf('../outputs/anp_nodalregression_modelprep.pdf')
pdf('../outputs/anp_nodalregression_priorpredictive.pdf')

#### import data and randomise node IDs -- ensures that the model is not reading a correlation between node ID and age and partitioning too much of the variance to node ID #####
## load time window 1 and remove additional data
load('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')
rm(list = ls()[! ls() %in% c('nodes','n_windows')]) ; gc()

## create full data frame with window variable
nodes_all <- nodes %>%
  mutate(window = 1)

## attach all other data frames
for(time_window in 2:n_windows){
  ## import workspace image for time window
  load(paste0('anp_edgecalculations/anpshort',time_window,'_edgeweights_conditionalprior.RData'))
  rm(list = ls()[! ls() %in% c('nodes','nodes_all','time_window')]) ; gc()

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

#### reimport data and calculate centrality #####
# load('anp_nodalregression/anp_short_nodal_dataprep.RData')
## load time window 1 and remove additional data
load('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')
rm(list = ls()[! ls() %in% c('nodes_all', 'cdf_1', 'edge_samples','n_windows')]) ; gc()

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
  rm(list = ls()[! ls() %in% c('nodes', 'cdf', 'edge_samples','n_windows','nodes_all','covs_all','cents_all','time_window','extract_eigen_centrality')]) ; gc()

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

#### create data list #####
load('anp_nodalregression/anp_short_nodal_modelprep.RData')
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

#### run model ####
## set model parameters
n_chains <- 4
n_samples <- 1000

## load model
nodal_regression <- cmdstan_model('models/eigen_regression_anp.stan')

## run model
fit_anp_nodal <- nodal_regression$sample(data = model_data_list,
                                         chains = n_chains, parallel_chains = n_chains,
                                         threads_per_chain = 4,
                                         iter_warmup = n_samples, iter_sampling = n_samples)

## save workspace
save.image('anp_nodalregression/anp_short_nodal.RData')

#### check outputs ####
load('anp_nodalregression/anp_short_nodal.RData')

## define PDF output
pdf('../outputs/anp_nodalregression_modelchecks.pdf')

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
  dplyr::select(grep('rand_window', colnames(params), value=TRUE))
rand_node <- params %>% 
  dplyr::select(grep('rand_node', colnames(params), value=TRUE))

## traceplot all parameters
#traceplot(fit_anp_nodal, pars = c('intercept','beta_age','sigma','predictor[1]','predictor[50]','predictor[100]'))
plot_params <- c('intercept','beta_age','sigma')
plot_rand_nodes <- colnames(params)[sample(grep('rand_node',colnames(params)), size = 12, replace = F)]
plot_rand_windows1 <- c('rand_window[1]','rand_window[2]','rand_window[3]','rand_window[4]','rand_window[5]','rand_window[6]','rand_window[7]','rand_window[8]','rand_window[9]')
plot_rand_windows2 <- c('rand_window[10]','rand_window[11]','rand_window[12]','rand_window[13]','rand_window[14]','rand_window[15]','rand_window[16]','rand_window[17]','rand_window[18]')
plot_rand_windows3 <- c('rand_window[19]','rand_window[20]','rand_window[21]','rand_window[22]','rand_window[23]','rand_window[24]','rand_window[25]','rand_window[26]','rand_window[27]')
plot_rand_windows4 <- c('rand_window[28]','rand_window[29]','rand_window[30]','rand_window[31]','rand_window[32]','rand_window[33]','rand_window[34]','rand_window[35]','rand_window[36]')

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
traceplot(parameter_df = params, parameters_to_plot = plot_params)
traceplot(parameter_df = params, parameters_to_plot = plot_rand_nodes)
traceplot(parameter_df = params, parameters_to_plot = plot_rand_windows1)
traceplot(parameter_df = params, parameters_to_plot = plot_rand_windows2)
traceplot(parameter_df = params, parameters_to_plot = plot_rand_windows3)
traceplot(parameter_df = params, parameters_to_plot = plot_rand_windows4)
rm(plot_params,plot_rand_nodes,plot_rand_windows1,plot_rand_windows2,plot_rand_windows3,plot_rand_windows4) ; gc()

## save workspace
save.image('anp_nodalregression/anp_short_nodal.RData')

#### posterior predictive check ####
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
    lines(density(MASS::mvrnorm(1, mu, sigma)),
          col = rgb(0, 0, 1, 0.25))
  }
}

## check on standardised scale
par(mfrow = c(3,3))
# plot(density(cents_all[1, which(nodes_all$window == 1)]), las = 1, ylim = c(0,1),
#      #plot(density(cents_all_std[1, which(nodes_all$window == 1)]), las = 1, ylim = c(0,1),
#      main = "Posterior predictive check (standardised centrality):\nblack = data, blue = predicted",
#      col=rgb(0, 0, 0, 0.25))
# eigen_data1 <- list(num_nodes_window1 = length(which(nodes_all$window == 1)),
#                     centrality_mu_1 = sim_cent_mu_1,
#                     centrality_cov_1 = sim_cent_cov_1,
#                     node_age = sim$age_std[sim$window == 1],
#                     window = 1,
#                     nodes = sim$node_random[sim$window == 1],
#                     nodes_window1 = sim$node_random[sim$window == 1])
# for (i in 1:100) {
#   j <- sample(1:length(params$beta_age), 1)
#   lines(density(cents_all[j, which(sim$window == 1)]), col=rgb(0, 0, 0, 0.25))
#   #lines(density(cents_all_std[j, which(sim$window == 1)]), col=rgb(0, 0, 0, 0.25))
#   mu <- params$beta_age[j]*eigen_data1$node_age + params$intercept[j]
#   for(k in 1:length(mu)) {
#     mu[k] <- mu[k] + as.numeric(rand_window[j,eigen_data1$window]) + as.numeric(rand_node[j,eigen_data1$nodes[k]])
#   }
#   sigma <- sim_cent_cov_1 + diag(rep(params$sigma[j], eigen_list$num_nodes_window1))
#   lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
# }
for(time_window in 1:n_windows){
  post_pred_check(centrality_matrix = cents_all,
                  nodes_df = nodes_all, time_window = time_window,
                  cent_cov = covs_all[[time_window]],
                  parameters = params)
}
par(mfrow = c(1,1))

#### predict from model -- standardised scale ####
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
    window_effect <- parameters %>% dplyr::select(grep('rand_window', colnames(parameters), value=TRUE))
    for(i in 1:nrow(mean_matrix)){
      mean_matrix[i,] <- mean_matrix[i,] + as.numeric(window_effect[i,predict_df$window])
    }
  }
  
  if(include_node == TRUE){
    node_effect <- parameters %>% dplyr::select(grep('rand_node', colnames(parameters), value=TRUE))
    for(i in 1:nrow(mean_matrix)){
      mean_matrix[i,] <- mean_matrix[i,] + as.numeric(node_effect[i, predict_df$node_random])
    }
  }
  
  return(mean_matrix)
  
}

## get mean predictions
mu_predictions <- get_mean_predictions(predict_df = nodes_all, parameters = params,
                                       include_window = TRUE, include_node = TRUE)

## add mean and CI of predicted means to input data frame for comparison
nodes_all$mu_mean <- apply(mu_predictions, 2, mean)
nodes_all$mu_upr <- apply(mu_predictions, 2, quantile, prob = 0.975)
nodes_all$mu_lwr <- apply(mu_predictions, 2, quantile, prob = 0.025)

## plot mean of model vs mean of raw data
ggplot()+
  geom_point(data = nodes_all, #size = 0.5,
             mapping = aes(x = mean_eigen, y = mu_mean, colour = as.factor(window)))+
  scale_colour_viridis_d()+
  #facet_wrap(. ~ as.factor(window))+
  labs(colour = 'window',
       x = 'simulated mean (standardised)',
       y = 'predicted mean (standardised)')+
  geom_abline(slope = 1, intercept = 0) # add line showing where points would lie if model fit was perfect

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
save.image('anp_nodalregression/anp_short_nodal.RData')

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
mean(contrast_std)                                    # should be very nodes_allilar to mean of output parameter
quantile(contrast_std, prob = c(0.025, 0.975))        # very wide

## contrast predictions on output scale -- reportable values of effect of plus 1 year
contrast <- contrast_std / sd(nodes_all$age)          # convert to outcome scale
head(contrast[,1:5])                                  # check matrix looks right
mean(contrast)                                        # should be very nodes_allilar to mean of output parameter
quantile(contrast, prob = c(0.025, 0.975))            # very wide

## save workspace
save.image('anp_nodalregression/anp_short_nodal.RData')

#### final "clean" plots using hypothetical data ####
## create dummy dataset
newdat <- nodes_all %>% 
  select(node_random, age, age_std, window)

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

## add CI of predicted data points to input data frame for comparison
newdat$predict_pred_lwr <- apply(fake_pred, 2, quantile, prob = 0.025)
newdat$predict_pred_upr <- apply(fake_pred, 2, quantile, prob = 0.975)

## plot predictions
newdat_summary <- newdat %>% 
  group_by(age, window) %>% 
  mutate(predict_pred_lwr = mean(predict_pred_lwr),
         predict_pred_upr = mean(predict_pred_upr)) %>% 
  select(age, predict_pred_lwr, predict_pred_upr,window) %>% 
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
  geom_point(data = sim,
             aes(x = age, y = mu))+              # original data points (standardised centrality, actual age)
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  facet_wrap(. ~ as.factor(window))+             # separate plots per window
  theme(legend.position = 'bottom')+
  labs(colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'age (years)')

## convert raw data to data frame
cents_all_df <- cents_all %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'node_window', values_to = 'eigen') %>% 
  left_join(sim[,c('node_window','age','window')], by = 'node_window')

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
  geom_point(data = sim,
             aes(x = age, y = mu),
             size = 0.5, colour = 'white')+      # original mean data points (standardised centrality, actual age)
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  facet_wrap(. ~ as.factor(window))+             # separate plots per window
  theme(legend.position = 'bottom')+
  labs(colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'age (years)')

## clean up and save workspace
save.image('anp_nodalregression/anp_short_nodal.RData')
dev.off()

# ########################
# ### compute normal approximation ####
# ## check covariance
# plot(centrality_samples_std[, 1], centrality_samples_std[, 2],
#      xlab = 'standardised eigenvectors ID1',
#      ylab = 'standardised eigenvectors ID2',
#      las = 1, pch = 19, col = rgb(0,0,1,0.2))
# 
# ## compute normal approximation
# centrality_mu <- apply(centrality_samples_std, 2, mean)
# centrality_cov <- cov(centrality_samples_std)
# 
# centrality_samples_sim <- MASS::mvrnorm(1e5, centrality_mu, centrality_cov)
# 
# plot(density(centrality_samples_std[, 1]), lwd = 2, main = "Estimated tandardised centrality vs normal approximation", xlab = "Logit edge weight")
# lines(density(centrality_samples_sim[, 1]), col = rgb(0,0,1,0.5), lwd = 2)
# 
# ##### prior predictive check ####
# ## simulate
# age <- 10:60
# mu_mean <- 0
# mu_stdv <- 0.1
# mu <- rnorm(100, mu_mean, mu_stdv)
# mu <- sort(mu)
# mu_mtrx <- matrix(NA, nrow = length(mu), ncol = length(age))
# for(i in 1:nrow(mu_mtrx)){
#  mu_mtrx[i,] <- mu[i]*age
# }
# sigma_range <- rexp(25, 2)
# sigma_range <- sort(sigma_range)
# par(mfrow = c(5,5), mai = c(0.2,0.2,0.2,0.2))
# for(j in 1:25){
#  plot(NULL, xlim = c(10,60), ylim = c(0,1), las = 1,
#       xlab = '', ylab = '')
#  sigma <- diag(rep(sigma_range[j], length(age)))
#  for(i in 1:length(mu)){
#  predictor <- mu_mtrx[i,]
#  y <- MASS::mvrnorm(1, predictor, sigma)
#  lines(x = age, y = invlogit(y), col = rgb(0,0,1,1))
#  }
# }
# par(mfrow = c(1,1))
# dev.off()
# rm(list = ls())
# 
# ### run model ####
# ## standardise age variable
# nodes$age_std <- (nodes$age - mean(nodes$age))/sd(nodes$age)
# 
# ## create data list
# eigen_list <- list(num_nodes = n_eles,
#                    nodes = nodes$node_rank,
#                    centrality_mu = centrality_mu,
#                    centrality_cov = centrality_cov,
#                    node_age = nodes$age_std)
# 
# # load model
# nodal_regression <- stan_model('models/eigen_regression_intercept_standardiseall.stan') # nodal_regression <- cmdstan_model('models/eigen_regression_intercept_standardiseall.stan')
# 
# ## run model
# fit_anp1_eigen <- sampling(nodal_regression, data = eigen_list,
#                            cores = n_chains, chains = n_chains,
#                            iter = n_samples*2, warmup = n_samples) # fit_anp1_eigen <- nodal_regression$sample(data = eigen_list, chains = n_chains, parallel_chains = n_chains, iter_warmup = n_samples, iter_sampling = n_samples)
# 
# ## save output
# rm(edge_samples, adj_tensor, i) ; gc()
# save.image('anp_nodalregression/anpshort1_nodalregression_conditionaledge_rstan.RData')
# 
# ### posterior check ####
# # load('anp_nodalregression/anpshort1_nodalregression_conditionaledge.RData')
# ## summarise model
# #fit_anp1_eigen$summary()
# summary(fit_anp1_eigen)
# 
# ## extract posterior
# params <- rstan::extract(fit_anp1_eigen)
# 
# ## traceplot linear effect size
# traceplot(fit_anp1_eigen, pars = c('beta_age','sigma','predictor[1]','predictor[10]','predictor[20]'))
#  params %>%
#    select(intercept,beta_age,sigma,`predictor[1]`,`predictor[50]`,`predictor[100]`) %>%
#    pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'draw') %>%
#    mutate(chain_position = rep(rep(1:n_samples, each = 6), n_chains),
#           chain = rep(1:n_chains, each = 6*n_samples)) %>%
#    #filter(chain == 4) %>%
#    ggplot(aes(x = chain_position, y = draw, colour = as.factor(chain)))+
#    geom_line()+
#    facet_wrap(. ~ parameter, scales = 'free_y')+
#    theme_bw()+
#    theme(legend.position = 'none') # chains 1-3 are good, chain 4 poorly mixed after about 600 draws
# 
# ## check model fit
#  summary <- summary(fit_anp1_eigen)
# par(mfrow = c(3,1))
#  hist(summary$rhat, breaks = 50)
#  hist(summary$ess_bulk, breaks = 50)
#  hist(summary$ess_tail, breaks = 50)
#  par(mfrow = c(1,1))
# 
# ## posterior predictive check
# plot(density(centrality_samples_std[1, ]), main="Posterior predictive check (standardised centrality):\nblack = data, blue = predicted", col=rgb(0, 0, 0, 0.25))
# for (i in 1:100) {
#   j <- sample(1:length(params$beta_age), 1)
#   lines(density(centrality_samples_std[j, ]), col=rgb(0, 0, 0, 0.25))
#   mu <- params$beta_age[j]*eigen_list$node_age + params$intercept[j]
#   sigma <- centrality_cov + diag(rep(params$sigma[j], n_eles))
#   lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
# }
# 
# ## extract original values from output
# ( fit_slope <- params$beta_age * (sd(centrality_mu_ustd)/sd(nodes$age)) )
# ( fit_intcp <- mean(nodes$age) - fit_slope * mean(nodes$age) )
# plot((nodes$age*mean(fit_slope) + mean(fit_intcp)) ~ nodes$age,
#      main = 'mean predictions vs age')
# 
# ## plot density curves for each
# plot(density(fit_slope)) ; abline(v = 0, lty = 2)
# plot(density(fit_intcp)) ; abline(v = 0, lty = 2)
# 
# ## save image
# save.image('anp_nodalregression/anpshort1_nodalregression_conditionaledge_rstan.RData')
# 
# ### predict from model ####
# #load('anp_nodalregression/anpshort1_nodalregression_conditionaledge_rstan.RData')
# #(summary <- as.data.frame(round(summary(fit_anp1_eigen)$summary[1:2, c(1, 4, 8)], 3)))
# #summary$parameter <- rownames(summary)
# 
# ## calculate mean predictions for model -- THERE'S SOMETHING WRONG WITH MY LOWER AND UPPER PREDICTIVE BOUNDS HERE, AS THEY FLIP WHICH IS LARGER IN A COUPLE OF PLACES, BUT i HAVE NO IDEA WHAT SO COME BACK TO THIS
# nodes$mean_predict <- nodes$age_std * mean(params$beta_age) + mean(params$intercept)           # mean prediction = age * mean slope
# nodes$lwr_predict <- nodes$age_std * quantile(params$beta_age, probs = 0.025) + quantile(params$intercept, probs = 0.025)
# nodes$upr_predict <- nodes$age_std * quantile(params$beta_age, probs = 0.975) + quantile(params$intercept, probs = 0.975)
# 
# ## compare to raw data
# plot(mean_predict ~ age_std, data = nodes, type = 'l', las = 1,
#      ylim = c(min(c(nodes$mean_eigen_std,nodes$mean_predict)),
#               max(c(nodes$mean_eigen_std,nodes$mean_predict))))
# points(nodes$mean_eigen_std ~ nodes$age_std)
# polygon(y = c(nodes$lwr_predict, rev(nodes$upr_predict)),
#         x = c(nodes$age_std, rev(nodes$age_std)),
#         col = rgb(1,1,0,0.5), border = NA)
# 
# ## reverse standardisation
# nodes$mean_eigen_ustd
# (nodes$predict_ustd <- nodes$mean_predict * sd(nodes$mean_eigen_ustd) + mean(nodes$mean_eigen_ustd))
# plot(nodes$predict_ustd ~ nodes$mean_eigen_ustd,
#      ylim = c(min(nodes$mean_eigen_ustd), max(nodes$mean_eigen_ustd)))
# abline(a = 0, b = 1)
# 
# ## compare unstandardised to raw data -- unstandardising lwr and upr does not work like this!
# nodes$lwr_ustd <- nodes$age * quantile(fit_slope, probs = 0.025) + quantile(fit_intcp, probs = 0.025)
# nodes$upr_ustd <- nodes$age * quantile(fit_slope, probs = 0.975) + quantile(fit_intcp, probs = 0.975)
# plot(nodes$predict_ustd ~ nodes$age, type = 'l', las = 1,
#      ylim = c(min(nodes$mean_eigen_ustd), max(nodes$mean_eigen_ustd)))    # plot against age
# points(nodes$mean_eigen_ustd ~ nodes$age)                                 # add raw points
# polygon(y = c(nodes$lwr_ustd, rev(nodes$upr_ustd)), x = c(nodes$age, rev(nodes$age)),
#         col = rgb(1,1,0,0.5), border = NA) # THIS DOES NOT WORK CURRENTLY!
# 
# ## convert to invlogit scale
# nodes$mean_predict_invlogit <- invlogit(nodes$predict_ustd)
# nodes$lwr_invlogit <- invlogit(nodes$lwr_ustd)
# nodes$upr_invlogit <- invlogit(nodes$upr_ustd)
# nodes$mu_invlogit <- invlogit(nodes$mean_eigen_ustd)
# plot(nodes$mu_invlogit ~ nodes$age, ylim = c(0,0.3))
# lines(nodes$mean_predict_invlogit ~ nodes$age)
# polygon(y = c(nodes$lwr_invlogit, rev(nodes$upr_invlogit)), x = c(nodes$age, rev(nodes$age)),
#         col = rgb(1,1,0,0.5), border = NA)
# 
# ## simulate full predictions for model
# mod_mu <- data.frame(id = 1:100,
#                      age = seq(from = min(nodes$age),
#                                to = max(nodes$age),
#                                length.out = 100))
# sim <- matrix(NA, nrow = n_chains*n_samples, ncol = nrow(mod_mu),
#               dimnames = list(1:(n_chains*n_samples), mod_mu$age))
# for(i in 1:nrow(sim)){
#   for(j in 1:ncol(sim)){
#     sim[i,j] <- MASS::mvrnorm(n = 1, mu = params$beta_age[i]*mod_mu$age[j],
#                               Sigma = params$sigma[i])
#   }
# }
# 
# ## plot simulations
# sim_df <- as.data.frame(sim) %>%
#   pivot_longer(cols = everything(), names_to = 'age', values_to = 'eigen_sim')
# plot(sim_df$eigen_sim ~ sim_df$age, col = rgb(0,0,0,0.01), pch = 19, cex = 0.5)
# 
# ## summarise simulations
# sim_summary <- data.frame(age = as.numeric(unique(sim_df$age)),
#                           lwr = NA, mid = NA, upr = NA)
# for(i in 1:nrow(sim_summary)){
#   x <- sim_df %>% filter(age == sim_summary$age[i])
#   sim_summary$lwr[i] <- quantile(x$eigen_sim, 0.025)
#   sim_summary$mid[i] <- quantile(x$eigen_sim, 0.5)
#   sim_summary$upr[i] <- quantile(x$eigen_sim, 0.975)
# }
# 
# ## plot raw with model output
# df_long <- df_long %>%
#   group_by(node_rank) %>%
#   mutate(mean_eigen = mean(centrality)) %>%
#   ungroup() %>%
#   left_join(nodes[,c('node_rank','sightings')], by = 'node_rank')
# ggplot()+
#   # geom_ribbon(data = sim_summary, aes(x = age, ymin = lwr, ymax = upr),       # shade simulations
#   #             colour = 'transparent', fill = rgb(0,0,0,0.1))+
#   # geom_ribbon(data = mod_mu, aes(x = age, ymin = lwr, ymax = upr),            # shade mean distribution
#   #             colour = 'transparent', fill = rgb(33/255, 145/255, 140/255, 0.5))+
#   geom_point(data = df_long, aes(x = age, y = centrality),                    # all eigenvector draws
#              colour = rgb(253/255, 231/255, 37/255, 0.01))+
#   geom_point(data = df_long, aes(x = age, y = mean_eigen, size = sightings),  # mean eigenvector
#              colour = rgb(68/255, 1/255, 84/255))+
#   # geom_line(data = mod_mu, aes(x = age, y = mid),                             # mean line
#   #           colour = rgb(33/255, 145/255, 140/255), linewidth = 1)+
#   scale_x_continuous('age (years)')+
#   scale_y_continuous('eigenvector centrality (standardised)')+
#   theme_classic()+
#   theme(axis.text = element_text(size = 18),
#         axis.title = element_text(size = 22),
#         legend.text = element_text(size = 18),
#         legend.title = element_text(size = 22))
# ggsave(filename = '../outputs/anpshort1_nodalregression_raw.png', device = 'png',
#        plot = last_plot())
# 
# ## save output
# save.image('anp_nodalregression/anpshort1_nodalregression_conditionaledge_predictions.RData')
# dev.off()
# rm(list = ls()[!ls() %in% c('nodal_regression','n_windows')])
# 
# ##### run model with rstan -- time window 2-36 ####
# set.seed(12345)
# 
# for(time_window in 2:n_windows){
#   ### load data ####
#   load(paste0('anp_edgecalculations/anpshort',time_window,'_edgeweights_conditionalprior.RData'))
#   rm(counts_df, edges, edgelist, edge_binary, fit_edges_anp, summary, make_edgelist, plot_network_threshold_anp, periods) ; gc()
# 
#   ## set up pdf
#   pdf(paste0('../outputs/anpshort',time_window,'_nodalregression.pdf'))
# 
#   ## add progress marker
#   print(paste0('start window ',time_window,' at ',Sys.time()))
# 
#   ### extract centralities ####
#   ## build adjacency tensor
#   n_eles <- nrow(nodes)
#   cdf$node_1_id <- as.integer(as.factor(cdf$node_1))
#   cdf$node_2_id <- as.integer(as.factor(cdf$node_2))+1
#   adj_tensor <- array(0, c(n_chains*n_samples, n_eles, n_eles))
#   for (dyad_id in 1:n_dyads) {
#     dyad_row <- cdf[dyad_id, ]
#     adj_tensor[, dyad_row$node_1_id, dyad_row$node_2_id] <- edge_samples[, dyad_id]
#     adj_tensor[, dyad_row$node_2_id, dyad_row$node_1_id] <- edge_samples[, dyad_id]
#   }
#   head(adj_tensor[1,,])
# 
#   ## calculate centrality and store posterior samples in a matrix
#   centrality_samples_invlogit <- matrix(0, n_chains*n_samples, n_eles)
#   for(i in 1:(n_chains*n_samples)){
#     centrality_samples_invlogit[i, ] <- sna::evcent(adj_tensor[i,,], gmode = 'graph')
#   }
#   head(centrality_samples_invlogit) # Unstandardised eigenvector centrality, 0-1 bounded
# 
#   ## convert to logit scale
#   centrality_samples <- logit(centrality_samples_invlogit)
#   head(centrality_samples)          # unstandardised eigenvector centrality, natural scale
#   centrality_mu_ustd <- apply(centrality_samples, 2, mean)
# 
#   ## add mean estimate per ID to nodes data frame
#   df <- data.frame(node_1_id = 1:n_eles,
#                    mean_eigen_ustd = centrality_mu_ustd) %>%
#     left_join(distinct(cdf[,c('node_1_id','id_1')]), by = 'node_1_id') %>%
#     rename(node_rank = node_1_id,
#            id = id_1)
#   df$id <- ifelse(is.na(df$id) == FALSE, df$id,
#                   nodes$id[! nodes$id %in% df$id])
#   nodes <- nodes %>%
#     mutate(node_rank = as.integer(as.factor(node))) %>%
#     left_join(df, by = c('node_rank','id'))
#   rm(df); gc()
# 
#   ## standardise
#   centrality_samples_std <- centrality_samples
#   for(i in 1:nrow(centrality_samples_std)){
#     centrality_samples_std[i,] <- (centrality_samples[i,] - mean(centrality_samples[i,])) / sd(centrality_samples[i,])
#   }
#   head(centrality_samples_std)      # standardised eigenvector centrality, natural scale
# 
#   ## visualise centralities
#   df_wide <- data.frame(centrality_samples_std)
#   colnames(df_wide) <- 1:n_eles
#   df_long <- pivot_longer(df_wide, cols = everything(),
#                           names_to = "node_rank", values_to = "centrality") %>%
#     mutate(node_rank = as.integer(node_rank)) %>%
#     left_join(nodes[,c('node_rank','age')], by = 'node_rank')
#   df_long %>%
#     mutate(nodes_reordered = fct_reorder(.f = as.factor(node_rank), .x = age, .desc = T)) %>%
#     ggplot(aes(x = centrality, fill = age)) +
#     geom_density(linewidth = 0.4) +
#     facet_grid(rows = vars(as.factor(nodes_reordered)), scales = "free") +
#     labs(x = "Eigenvector centrality (standardised)") +
#     theme_void() +
#     theme(strip.text.y = element_text(size = 12),
#           axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
#           axis.title.x = element_text(size = 12),
#           plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
# 
#   ## add mean estimate per ID to nodes data frame
#   df <- df_long %>%
#     group_by(node_rank) %>%
#     mutate(mean_eigen_std = mean(centrality)) %>%
#     dplyr::select(node_rank, mean_eigen_std) %>%
#     distinct()
#   nodes <- nodes %>%
#     left_join(df, by = 'node_rank')
#   rm(df); gc()
# 
#     ## add progress marker
#   print('eigen values extracted')
# 
#   ### compute normal approximation ####
#   ## check covariance
#   plot(centrality_samples_std[, 1], centrality_samples_std[, 2],
#        xlab = 'standardised eigenvectors ID1',
#        ylab = 'standardised eigenvectors ID2',
#        las = 1, pch = 19, col = rgb(0,0,1,0.2))
# 
#   ## compute normal approximation
#   centrality_mu <- apply(centrality_samples_std, 2, mean)
#   centrality_cov <- cov(centrality_samples_std)
# 
#   centrality_samples_sim <- MASS::mvrnorm(1e5, centrality_mu, centrality_cov)
# 
#   plot(density(centrality_samples_std[, 1]), lwd = 2, main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
#   lines(density(centrality_samples_sim[, 1]), col = rgb(0,0,1,0.5), lwd = 2)
# 
#   ## add progress marker
#   print('normal approximation computed')
# 
#   ### run model ####
#   ## standardise age variable
#   nodes$age_std <- (nodes$age - mean(nodes$age))/sd(nodes$age)
# 
#   ## create data list
#   eigen_list <- list(num_nodes = n_eles,
#                      nodes = nodes$node_rank,
#                      centrality_mu = centrality_mu,
#                      centrality_cov = centrality_cov,
#                      node_age = nodes$age_std)
# 
#   ## run model
#   fit_anp_eigen <- sampling(nodal_regression, data = eigen_list,
#                              cores = n_chains, chains = n_chains,
#                              iter = n_samples*2, warmup = n_samples)
#   # fit_anp_eigen <- nodal_regression$sample(data = eigen_list, chains = n_chains, parallel_chains = n_chains, iter_warmup = n_samples, iter_sampling = n_samples)
# 
#   ## save output
#   rm(edge_samples, adj_tensor, i) ; gc()
#   save.image(paste0('anp_nodalregression/anpshort',time_window,'_nodalregression_conditionaledge_rstan.RData'))
# 
#   ## add progress marker
#   print('model run')
# 
#   ### posterior check -- new ####
#   ## summarise model
#   #fit_anp_eigen$summary()
#   summary(fit_anp_eigen)
# 
#   ## extract posterior
#   params <- rstan::extract(fit_anp_eigen)
# 
#   ## traceplot linear effect size
#   traceplot(fit_anp_eigen, pars = c('beta_age','sigma','predictor[1]','predictor[10]','predictor[20]'))
#   # params %>%
#   #   select(intercept,beta_age,sigma,`predictor[1]`,`predictor[50]`,`predictor[100]`) %>%
#   #   pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'draw') %>%
#   #   mutate(chain_position = rep(rep(1:n_samples, each = 6), n_chains),
#   #          chain = rep(1:n_chains, each = 6*n_samples)) %>%
#   #   #filter(chain == 4) %>%
#   #   ggplot(aes(x = chain_position, y = draw, colour = as.factor(chain)))+
#   #   geom_line()+
#   #   facet_wrap(. ~ parameter, scales = 'free_y')+
#   #   theme_bw()+
#   #   theme(legend.position = 'none') # chains 1-3 are good, chain 4 poorly mixed after about 600 draws
# 
#   ## check model fit
#   # summary <- summary(fit_anp_eigen)
#   # par(mfrow = c(3,1))
#   # hist(summary$rhat, breaks = 50)
#   # hist(summary$ess_bulk, breaks = 50)
#   # hist(summary$ess_tail, breaks = 50)
#   # par(mfrow = c(1,1))
# 
#   ## posterior predictive check
#   plot(density(centrality_samples_std[1, ]), main="Posterior predictive check (standardised centrality):\nblack = data, blue = predicted", col=rgb(0, 0, 0, 0.25))
#   for (i in 1:100) {
#     j <- sample(1:length(params$beta_age), 1)
#     lines(density(centrality_samples_std[j, ]), col=rgb(0, 0, 0, 0.25))
#     mu <- params$beta_age[j]*eigen_list$node_age + params$intercept[j]
#     sigma <- centrality_cov + diag(rep(params$sigma[j], n_eles))
#     lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
#   }
# 
#   ## add progress marker
#   print('posterior predictive check complete')
# 
#     ## extract original values from output
#   ( fit_slope <- params$beta_age * (sd(centrality_mu_ustd)/sd(nodes$age)) )
#   ( fit_intcp <- mean(nodes$age) - fit_slope * mean(nodes$age) )
#   plot((nodes$age*mean(fit_slope) + mean(fit_intcp)) ~ nodes$age,
#        main = 'mean predictions vs age')
# 
#   ## plot density curves for each
#   plot(density(fit_slope)) ; abline(v = 0, lty = 2)
#   plot(density(fit_intcp)) ; abline(v = 0, lty = 2)
# 
#   ## save output ####
#   save.image(paste0('anp_nodalregression/anpshort',time_window,'_nodalregression_conditionaledge_rstan.RData'))
#   dev.off()
#   rm(list = ls()[! ls() %in% c('time_window','nodal_regression')])
# 
#   ## add progress marker
#   print(paste0('time window ',time_window,' complete'))
# 
# }
# 
# ##### predict from models -- time windows 2-36 (short) ####
# 
# ##### run model with rstan -- time window 1-7 (long) ####
# set.seed(12345)
# 
# for(time_window in 1:7){
#   ### load data ####
#   load(paste0('anp_edgecalculations/anplong',time_window,'_edgeweights_conditionalprior.RData'))
#   rm(counts_df, edges, edgelist, edge_binary, fit_edges_anp, summary, make_edgelist, plot_network_threshold_anp, periods) ; gc()
#   
#   ## set up pdf
#   pdf(paste0('../outputs/anplong',time_window,'_nodalregression.pdf'))
#   
#   ## add progress marker
#   print(paste0('start window ',time_window,' at ',Sys.time()))
#   
#   ### extract centralities ####
#   ## build adjacency tensor
#   n_eles <- nrow(nodes)
#   cdf$node_1_id <- as.integer(as.factor(cdf$node_1))
#   cdf$node_2_id <- as.integer(as.factor(cdf$node_2))+1
#   adj_tensor <- array(0, c(n_chains*n_samples, n_eles, n_eles))
#   for (dyad_id in 1:n_dyads) {
#     dyad_row <- cdf[dyad_id, ]
#     adj_tensor[, dyad_row$node_1_id, dyad_row$node_2_id] <- edge_samples[, dyad_id]
#     adj_tensor[, dyad_row$node_2_id, dyad_row$node_1_id] <- edge_samples[, dyad_id]
#   }
#   head(adj_tensor[1,,])
#   
#   ## calculate centrality and store posterior samples in a matrix
#   centrality_samples_invlogit <- matrix(0, n_chains*n_samples, n_eles)
#   for(i in 1:(n_chains*n_samples)){
#     centrality_samples_invlogit[i, ] <- sna::evcent(adj_tensor[i,,], gmode = 'graph')
#   }
#   head(centrality_samples_invlogit) # Unstandardised eigenvector centrality, 0-1 bounded
#   
#   ## convert to logit scale
#   centrality_samples <- logit(centrality_samples_invlogit)
#   head(centrality_samples)          # unstandardised eigenvector centrality, natural scale
#   centrality_mu_ustd <- apply(centrality_samples, 2, mean)
#   
#   ## add mean estimate per ID to nodes data frame
#   df <- data.frame(node_1_id = 1:n_eles,
#                    mean_eigen_ustd = centrality_mu_ustd) %>% 
#     left_join(distinct(cdf[,c('node_1_id','id_1')]), by = 'node_1_id') %>% 
#     rename(node_rank = node_1_id,
#            id = id_1)
#   df$id <- ifelse(is.na(df$id) == FALSE, df$id,
#                   nodes$id[! nodes$id %in% df$id])
#   nodes <- nodes %>% 
#     mutate(node_rank = as.integer(as.factor(node))) %>% 
#     left_join(df, by = c('node_rank','id'))
#   rm(df); gc()
#   
#   ## standardise
#   centrality_samples_std <- centrality_samples
#   for(i in 1:nrow(centrality_samples_std)){
#     centrality_samples_std[i,] <- (centrality_samples[i,] - mean(centrality_samples[i,])) / sd(centrality_samples[i,])
#   }
#   head(centrality_samples_std)      # standardised eigenvector centrality, natural scale
#   
#   ## visualise centralities
#   df_wide <- data.frame(centrality_samples_std)
#   colnames(df_wide) <- 1:n_eles
#   df_long <- pivot_longer(df_wide, cols = everything(),
#                           names_to = "node_rank", values_to = "centrality") %>% 
#     mutate(node_rank = as.integer(node_rank)) %>% 
#     left_join(nodes[,c('node_rank','age')], by = 'node_rank')
#   df_long %>% 
#     mutate(nodes_reordered = fct_reorder(.f = as.factor(node_rank), .x = age, .desc = T)) %>% 
#     ggplot(aes(x = centrality, fill = age)) +
#     geom_density(linewidth = 0.4) +
#     facet_grid(rows = vars(as.factor(nodes_reordered)), scales = "free") +
#     labs(x = "Eigenvector centrality (standardised)") + 
#     theme_void() + 
#     theme(strip.text.y = element_text(size = 12),
#           axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
#           axis.title.x = element_text(size = 12),
#           plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#   
#   ## add mean estimate per ID to nodes data frame
#   df <- df_long %>% 
#     group_by(node_rank) %>% 
#     mutate(mean_eigen_std = mean(centrality)) %>% 
#     dplyr::select(node_rank, mean_eigen_std) %>% 
#     distinct()
#   nodes <- nodes %>% 
#     left_join(df, by = 'node_rank')
#   rm(df); gc()
#   
#   ## add progress marker
#   print('eigen values extracted')
#   
#   ### compute normal approximation ####
#   ## check covariance
#   plot(centrality_samples_std[, 1], centrality_samples_std[, 2],
#        xlab = 'standardised eigenvectors ID1',
#        ylab = 'standardised eigenvectors ID2',
#        las = 1, pch = 19, col = rgb(0,0,1,0.2))
#   
#   ## compute normal approximation
#   centrality_mu <- apply(centrality_samples_std, 2, mean)
#   centrality_cov <- cov(centrality_samples_std)
#   
#   centrality_samples_sim <- MASS::mvrnorm(1e5, centrality_mu, centrality_cov)
#   
#   plot(density(centrality_samples_std[, 1]), lwd = 2, main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
#   lines(density(centrality_samples_sim[, 1]), col = rgb(0,0,1,0.5), lwd = 2)
#   
#   ## add progress marker
#   print('normal approximation computed')
#   
#   ### run model ####
#   ## standardise age variable
#   nodes$age_std <- (nodes$age - mean(nodes$age))/sd(nodes$age)
#   
#   ## create data list
#   eigen_list <- list(num_nodes = n_eles,
#                      nodes = nodes$node_rank,
#                      centrality_mu = centrality_mu,
#                      centrality_cov = centrality_cov,
#                      node_age = nodes$age_std)
#   
#   ## run model
#   fit_anp_eigen <- sampling(nodal_regression, data = eigen_list,
#                             cores = n_chains, chains = n_chains,
#                             iter = n_samples*2, warmup = n_samples)
#   # fit_anp_eigen <- nodal_regression$sample(data = eigen_list, chains = n_chains, parallel_chains = n_chains, iter_warmup = n_samples, iter_sampling = n_samples)
#   
#   ## save output
#   rm(edge_samples, adj_tensor, i) ; gc()
#   save.image(paste0('anp_nodalregression/anplong',time_window,'_nodalregression_conditionaledge_rstan.RData'))
#   
#   ## add progress marker
#   print('model run')
#   
#   ### posterior check -- new ####
#   ## summarise model
#   #fit_anp_eigen$summary()
#   summary(fit_anp_eigen)
#   
#   ## extract posterior
#   params <- rstan::extract(fit_anp_eigen)
#   
#   ## traceplot linear effect size
#   traceplot(fit_anp_eigen, pars = c('beta_age','sigma','predictor[1]','predictor[10]','predictor[20]'))
#   # params %>% 
#   #   select(intercept,beta_age,sigma,`predictor[1]`,`predictor[50]`,`predictor[100]`) %>% 
#   #   pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'draw') %>% 
#   #   mutate(chain_position = rep(rep(1:n_samples, each = 6), n_chains),
#   #          chain = rep(1:n_chains, each = 6*n_samples)) %>% 
#   #   #filter(chain == 4) %>% 
#   #   ggplot(aes(x = chain_position, y = draw, colour = as.factor(chain)))+
#   #   geom_line()+
#   #   facet_wrap(. ~ parameter, scales = 'free_y')+
#   #   theme_bw()+
#   #   theme(legend.position = 'none') # chains 1-3 are good, chain 4 poorly mixed after about 600 draws
#   
#   ## check model fit
#   # summary <- summary(fit_anp_eigen)
#   # par(mfrow = c(3,1))
#   # hist(summary$rhat, breaks = 50)
#   # hist(summary$ess_bulk, breaks = 50)
#   # hist(summary$ess_tail, breaks = 50)
#   # par(mfrow = c(1,1))
#   
#   ## posterior predictive check
#   plot(density(centrality_samples_std[1, ]), main="Posterior predictive check (standardised centrality):\nblack = data, blue = predicted", col=rgb(0, 0, 0, 0.25))
#   for (i in 1:100) {
#     j <- sample(1:length(params$beta_age), 1)
#     lines(density(centrality_samples_std[j, ]), col=rgb(0, 0, 0, 0.25))
#     mu <- params$beta_age[j]*eigen_list$node_age + params$intercept[j]
#     sigma <- centrality_cov + diag(rep(params$sigma[j], n_eles))
#     lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
#   }
#   
#   ## add progress marker
#   print('posterior predictive check complete')
#   
#   ## extract original values from output
#   ( fit_slope <- params$beta_age * (sd(centrality_mu_ustd)/sd(nodes$age)) )
#   ( fit_intcp <- mean(nodes$age) - fit_slope * mean(nodes$age) )
#   plot((nodes$age*mean(fit_slope) + mean(fit_intcp)) ~ nodes$age,
#        main = 'mean predictions vs age')
#   
#   ## plot density curves for each
#   plot(density(fit_slope)) ; abline(v = 0, lty = 2)
#   plot(density(fit_intcp)) ; abline(v = 0, lty = 2)
#   
#   ## save output ####
#   save.image(paste0('anp_nodalregression/anplong',time_window,'_nodalregression_conditionaledge_rstan.RData'))
#   dev.off()
#   rm(list = ls()[! ls() %in% c('time_window','nodal_regression')])
#   
#   ## add progress marker
#   print(paste0('time window ',time_window,' complete'))
#   
# }
# 
# ##### predict from models -- time windows 1-7 (long) ####
# 
# 
# ##### run model with rstan -- all time windows together ####
# ## load data and remove additional data
# load('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')
# cdf_1 <- cdf_1 %>% dplyr::select(-dyad_rank)
# nodes_1 <- nodes %>% 
#   mutate(window = 1)
# edge_samples_1 <- edge_samples
# rm(counts_df, nodes, edges, edge_samples, edgelist, edge_binary, fit_edges_anp1, summary, make_edgelist, plot_network_threshold_anp, periods) ; gc()
# 
# edge_samples_list <- list()
# edge_samples_list[[1]] <- edge_samples_1
# 
# ### extract centralities ####
# ## build adjacency tensor
# n_eles <- nrow(nodes_1)
# cdf_1$node_1_id <- as.integer(as.factor(cdf_1$node_1))
# cdf_1$node_2_id <- as.integer(as.factor(cdf_1$node_2))+1
# adj_tensor <- array(0, c(n_chains*n_samples, n_eles, n_eles))
# for (dyad_id in 1:n_dyads) {
#   dyad_row <- cdf_1[dyad_id, ]
#   adj_tensor[, dyad_row$node_1_id, dyad_row$node_2_id] <- edge_samples_1[, dyad_id]
#   adj_tensor[, dyad_row$node_2_id, dyad_row$node_1_id] <- edge_samples_1[, dyad_id]
# }
# adj_tensor[1,,]
# 
# ## calculate centrality and store posterior samples in a matrix
# centrality_samples_invlogit_1 <- matrix(0, nrow = n_chains*n_samples, ncol = n_eles,
#                                       dimnames = list(NULL, paste0(nodes_1$id,'_',nodes_1$window)))
# for(i in 1:(n_chains*n_samples)){
#   centrality_samples_invlogit_1[i, ] <- sna::evcent(adj_tensor[i,,], gmode = 'graph')
# }
# head(centrality_samples_invlogit_1) # Unstandardised eigenvector centrality, 0-1 bounded
# 
# rm(adj_tensor, n_eles, dyad_row)
# cdf_1 <- cdf_1 %>% 
#   dplyr::select(-node_1_id, -node_2_id)
# 
# rm(edge_samples_1, dyad_id, i) ; gc()
# 
# ## combine data
# for(time_window in 2:n_windows){
#   # load data
#   load(paste0('anp_edgecalculations/anpshort',time_window,'_edgeweights_conditionalprior.RData'))
#   print(paste0('data loaded in for window ',time_window))
#   
#   # combine counts_df data frames -- should be exactly the same as original counts_df, but combine together like this just to be sure
#   cdf_1 <- rbind(cdf_1, cdf)
#   
#   # combine nodes data frames
#   nodes <- nodes %>% mutate(window = time_window)
#   nodes_1 <- rbind(nodes_1, nodes)
#   print(paste0('data frames combined for window ',time_window))
#   
#   # combine edge_samples matrices
#   edge_samples_list[[time_window]] <- edge_samples
#   print(paste0('edge matrix added to list for window ',time_window))
#   
#   # create adjacency tensor
#   n_eles <- nrow(nodes)
#   cdf$node_1_id <- as.integer(as.factor(cdf$node_1))
#   cdf$node_2_id <- as.integer(as.factor(cdf$node_2))+1
#   adj_tensor <- array(0, c(n_chains*n_samples, n_eles, n_eles))
#   for (dyad_id in 1:n_dyads) {
#     dyad_row <- cdf[dyad_id, ]
#     adj_tensor[, dyad_row$node_1_id, dyad_row$node_2_id] <- edge_samples[, dyad_id]
#     adj_tensor[, dyad_row$node_2_id, dyad_row$node_1_id] <- edge_samples[, dyad_id]
#   }
#   
#   # calculate centrality and store posterior samples in a matrix
#   centrality_samples_invlogit <- matrix(data = 0, nrow = n_chains*n_samples, ncol = n_eles,
#                                         dimnames = list(NULL, paste0(nodes$id,'_',nodes$window)))
#   for(i in 1:(n_chains*n_samples)){
#     centrality_samples_invlogit[i, ] <- sna::evcent(adj_tensor[i,,], gmode = 'graph')
#   }
#   centrality_samples_invlogit_1 <- cbind(centrality_samples_invlogit_1, centrality_samples_invlogit)
#   print(paste0('eigenvectors calculated for window ',time_window))
#   
# }
# cdf <- cdf_1 %>% rename(window = period)
# nodes <- nodes_1
# centrality_samples_invlogit <- centrality_samples_invlogit_1
# rm(counts_df, edges, edgelist, edge_binary, fit_edges_anp, summary, make_edgelist, plot_network_threshold_anp, periods) ; gc()
# rm(cdf_1, nodes_1, edge_samples, centrality_samples_invlogit_1, dyad_row, adj_tensor, dyad_id, i) ; gc()
# 
# ## set up pdf
# #pdf('../outputs/anpshort_nodalregression.pdf')
# 
# ## convert to logit scale
# centrality_samples <- logit(centrality_samples_invlogit)
# head(centrality_samples)          # unstandardised eigenvector centrality, natural scale
# centrality_mu_ustd <- apply(centrality_samples, 2, mean)
# 
# ## standardise
# centrality_samples_std <- centrality_samples
# for(i in 1:nrow(centrality_samples_std)){
#   centrality_samples_std[i,] <- (centrality_samples_std[i,] - mean(centrality_samples_std[i,])) / sd(centrality_samples_std[i,])
# }
# head(centrality_samples_std)          # unstandardised eigenvector centrality, natural scale
# centrality_mu_std <- apply(centrality_samples_std, 2, mean)
# 
# ## add mean estimate per ID to nodes data frame
# n_data <- nrow(nodes)
# df <- data.frame(node_1_id = 1:n_data,
#                  mean_eigen_ustd = centrality_mu_ustd,
#                  mean_eigen_std  = centrality_mu_std)
# df$node_window <- rownames(df)
# df <- df %>% 
#   separate(node_window, into = c('id','window'), sep = '_', remove = F) %>% 
#   mutate(window = as.numeric(window))
# nodes <- nodes %>%
#   left_join(df, by = c('id','window'))
# rm(df); gc()
# 
# ## visualise centralities
# df_wide <- data.frame(centrality_samples_std)
# colnames(df_wide) <- 1:n_data
# df_long <- pivot_longer(df_wide, cols = everything(),
#                         names_to = "node_1_id", values_to = "centrality") %>%
#   mutate(node_1_id = as.integer(node_1_id)) %>%
#   left_join(nodes[,c('node_1_id','age','id')], by = 'node_1_id')
# df_long %>%
#   filter(id %in% sample(nodes$id, 50, replace = F)) %>% 
#   mutate(nodes_reordered = fct_reorder(.f = as.factor(node_1_id), .x = age, .desc = T)) %>%
#   ggplot(aes(x = centrality, fill = age)) +
#   geom_density(linewidth = 0.4) +
#   facet_grid(rows = vars(as.factor(nodes_reordered)), scales = "free") +
#   labs(x = "Eigenvector centrality (standardised)") +
#   theme_void() +
#   theme(strip.text.y = element_text(size = 12),
#         axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
#         axis.title.x = element_text(size = 12),
#         plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
# 
# ### compute normal approximation ####
# ## check covariance
# plot(centrality_samples_std[, 1], centrality_samples_std[, 2],
#      xlab = 'standardised eigenvectors ID1',
#      ylab = 'standardised eigenvectors ID2',
#      las = 1, pch = 19, col = rgb(0,0,1,0.2))
# 
# plot(centrality_samples_std[, 1], centrality_samples_std[, nrow(nodes)],
#      xlab = 'standardised eigenvectors ID1',
#      ylab = 'standardised eigenvectors ID2',
#      las = 1, pch = 19, col = rgb(0,0,1,0.2))
# 
# ## compute normal approximation
# centrality_mu <- apply(centrality_samples_std, 2, mean)
# centrality_cov <- cov(centrality_samples_std)
# 
# centrality_samples_sim <- MASS::mvrnorm(1e5, centrality_mu, centrality_cov)
# 
# plot(density(centrality_samples_std[, 1]), lwd = 2, main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
# lines(density(centrality_samples_sim[, 1]), col = rgb(0,0,1,0.5), lwd = 2)
# 
# ### run model ####
# ## standardise age variable
# nodes$age_std <- (nodes$age - mean(nodes$age))/sd(nodes$age)
# 
# ## create data list
# eigen_list <- list(num_data = n_data,
#                    num_nodes = n_eles,
#                    num_windows = n_windows,
#                    centrality_mu = centrality_mu,
#                    centrality_cov = centrality_cov,
#                    node_age = nodes$age_std,
#                    window = nodes$window,
#                    nodes = nodes$node)
# 
# # load model
# nodal_regression <- stan_model('models/eigen_regression_combinewindows.stan') # nodal_regression <- cmdstan_model('models/eigen_regression_combinewindows.stan')
# 
# ## run model
# fit_anp1_eigen <- sampling(nodal_regression, data = eigen_list,
#                            cores = n_chains, chains = n_chains,
#                            iter = n_samples*2, warmup = n_samples) # fit_anp1_eigen <- nodal_regression$sample(data = eigen_list, chains = n_chains, parallel_chains = n_chains, iter_warmup = n_samples, iter_sampling = n_samples)
# 
# ## save output
# rm(edge_samples_list, df_long, df_wide, centrality_samples_sim, eigen_list, i) ; gc()
# save.image('anp_nodalregression/anpshortall_nodalregression.RData')
# 
# ### posterior check ####
# # load('anp_nodalregression/anpshort1_nodalregression_conditionaledge.RData')
# ## summarise model
# #fit_anp1_eigen$summary()
# summary(fit_anp1_eigen)
# 
# ## extract posterior
# params <- rstan::extract(fit_anp1_eigen)
# 
# ## traceplot linear effect size
# traceplot(fit_anp1_eigen, pars = c('beta_age','sigma','predictor[1]','predictor[10]','predictor[20]'))
# params %>%
#   select(intercept,beta_age,sigma,`predictor[1]`,`predictor[50]`,`predictor[100]`) %>%
#   pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'draw') %>%
#   mutate(chain_position = rep(rep(1:n_samples, each = 6), n_chains),
#          chain = rep(1:n_chains, each = 6*n_samples)) %>%
#   #filter(chain == 4) %>%
#   ggplot(aes(x = chain_position, y = draw, colour = as.factor(chain)))+
#   geom_line()+
#   facet_wrap(. ~ parameter, scales = 'free_y')+
#   theme_bw()+
#   theme(legend.position = 'none') # chains 1-3 are good, chain 4 poorly mixed after about 600 draws
# 
# ## check model fit
# summary <- summary(fit_anp1_eigen)
# par(mfrow = c(3,1))
# hist(summary$rhat, breaks = 50)
# hist(summary$ess_bulk, breaks = 50)
# hist(summary$ess_tail, breaks = 50)
# par(mfrow = c(1,1))
# 
# ## posterior predictive check
# plot(density(centrality_samples_std[1, ]), main="Posterior predictive check (standardised centrality):\nblack = data, blue = predicted", col=rgb(0, 0, 0, 0.25))
# for (i in 1:100) {
#   j <- sample(1:length(params$beta_age), 1)
#   lines(density(centrality_samples_std[j, ]), col=rgb(0, 0, 0, 0.25))
#   mu <- params$beta_age[j]*eigen_list$node_age + params$intercept[j]
#   sigma <- centrality_cov + diag(rep(params$sigma[j], n_eles))
#   lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
# }
# 
# ## extract original values from output
# ( fit_slope <- params$beta_age * (sd(centrality_mu_ustd)/sd(nodes$age)) )
# ( fit_intcp <- mean(nodes$age) - fit_slope * mean(nodes$age) )
# plot((nodes$age*mean(fit_slope) + mean(fit_intcp)) ~ nodes$age,
#      main = 'mean predictions vs age')
# 
# ## plot density curves for each
# plot(density(fit_slope)) ; abline(v = 0, lty = 2)
# plot(density(fit_intcp)) ; abline(v = 0, lty = 2)
# 
# ## save image
# save.image('anp_nodalregression/anpshort1_nodalregression_conditionaledge_rstan.RData')
# 
# # ### predict from model ####
# # #load('anp_nodalregression/anpshort1_nodalregression_conditionaledge_rstan.RData')
# # #(summary <- as.data.frame(round(summary(fit_anp1_eigen)$summary[1:2, c(1, 4, 8)], 3)))
# # #summary$parameter <- rownames(summary)
# #
# # ## calculate mean predictions for model -- THERE'S SOMETHING WRONG WITH MY LOWER AND UPPER PREDICTIVE BOUNDS HERE, AS THEY FLIP WHICH IS LARGER IN A COUPLE OF PLACES, BUT i HAVE NO IDEA WHAT SO COME BACK TO THIS
# # nodes$mean_predict <- nodes$age_std * mean(params$beta_age) + mean(params$intercept)           # mean prediction = age * mean slope
# # nodes$lwr_predict <- nodes$age_std * quantile(params$beta_age, probs = 0.025) + quantile(params$intercept, probs = 0.025)
# # nodes$upr_predict <- nodes$age_std * quantile(params$beta_age, probs = 0.975) + quantile(params$intercept, probs = 0.975)
# #
# # ## compare to raw data
# # plot(mean_predict ~ age_std, data = nodes, type = 'l', las = 1,
# #      ylim = c(min(c(nodes$mean_eigen_std,nodes$mean_predict)),
# #               max(c(nodes$mean_eigen_std,nodes$mean_predict))))
# # points(nodes$mean_eigen_std ~ nodes$age_std)
# # polygon(y = c(nodes$lwr_predict, rev(nodes$upr_predict)),
# #         x = c(nodes$age_std, rev(nodes$age_std)),
# #         col = rgb(1,1,0,0.5), border = NA)
# #
# # ## reverse standardisation
# # nodes$mean_eigen_ustd
# # (nodes$predict_ustd <- nodes$mean_predict * sd(nodes$mean_eigen_ustd) + mean(nodes$mean_eigen_ustd))
# # plot(nodes$predict_ustd ~ nodes$mean_eigen_ustd,
# #      ylim = c(min(nodes$mean_eigen_ustd), max(nodes$mean_eigen_ustd)))
# # abline(a = 0, b = 1)
# #
# # ## compare unstandardised to raw data -- unstandardising lwr and upr does not work like this!
# # nodes$lwr_ustd <- nodes$age * quantile(fit_slope, probs = 0.025) + quantile(fit_intcp, probs = 0.025)
# # nodes$upr_ustd <- nodes$age * quantile(fit_slope, probs = 0.975) + quantile(fit_intcp, probs = 0.975)
# # plot(nodes$predict_ustd ~ nodes$age, type = 'l', las = 1,
# #      ylim = c(min(nodes$mean_eigen_ustd), max(nodes$mean_eigen_ustd)))    # plot against age
# # points(nodes$mean_eigen_ustd ~ nodes$age)                                 # add raw points
# # polygon(y = c(nodes$lwr_ustd, rev(nodes$upr_ustd)), x = c(nodes$age, rev(nodes$age)),
# #         col = rgb(1,1,0,0.5), border = NA) # THIS DOES NOT WORK CURRENTLY!
# #
# # ## convert to invlogit scale
# # nodes$mean_predict_invlogit <- invlogit(nodes$predict_ustd)
# # nodes$lwr_invlogit <- invlogit(nodes$lwr_ustd)
# # nodes$upr_invlogit <- invlogit(nodes$upr_ustd)
# # nodes$mu_invlogit <- invlogit(nodes$mean_eigen_ustd)
# # plot(nodes$mu_invlogit ~ nodes$age, ylim = c(0,0.3))
# # lines(nodes$mean_predict_invlogit ~ nodes$age)
# # polygon(y = c(nodes$lwr_invlogit, rev(nodes$upr_invlogit)), x = c(nodes$age, rev(nodes$age)),
# #         col = rgb(1,1,0,0.5), border = NA)
# #
# # ## simulate full predictions for model
# # mod_mu <- data.frame(id = 1:100,
# #                      age = seq(from = min(nodes$age),
# #                                to = max(nodes$age),
# #                                length.out = 100))
# # sim <- matrix(NA, nrow = n_chains*n_samples, ncol = nrow(mod_mu),
# #               dimnames = list(1:(n_chains*n_samples), mod_mu$age))
# # for(i in 1:nrow(sim)){
# #   for(j in 1:ncol(sim)){
# #     sim[i,j] <- MASS::mvrnorm(n = 1, mu = params$beta_age[i]*mod_mu$age[j],
# #                               Sigma = params$sigma[i])
# #   }
# # }
# #
# # ## plot simulations
# # sim_df <- as.data.frame(sim) %>%
# #   pivot_longer(cols = everything(), names_to = 'age', values_to = 'eigen_sim')
# # plot(sim_df$eigen_sim ~ sim_df$age, col = rgb(0,0,0,0.01), pch = 19, cex = 0.5)
# #
# # ## summarise simulations
# # sim_summary <- data.frame(age = as.numeric(unique(sim_df$age)),
# #                           lwr = NA, mid = NA, upr = NA)
# # for(i in 1:nrow(sim_summary)){
# #   x <- sim_df %>% filter(age == sim_summary$age[i])
# #   sim_summary$lwr[i] <- quantile(x$eigen_sim, 0.025)
# #   sim_summary$mid[i] <- quantile(x$eigen_sim, 0.5)
# #   sim_summary$upr[i] <- quantile(x$eigen_sim, 0.975)
# # }
# #
# # ## plot raw with model output
# # df_long <- df_long %>%
# #   group_by(node_rank) %>%
# #   mutate(mean_eigen = mean(centrality)) %>%
# #   ungroup() %>%
# #   left_join(nodes[,c('node_rank','sightings')], by = 'node_rank')
# # ggplot()+
# #   # geom_ribbon(data = sim_summary, aes(x = age, ymin = lwr, ymax = upr),       # shade simulations
# #   #             colour = 'transparent', fill = rgb(0,0,0,0.1))+
# #   # geom_ribbon(data = mod_mu, aes(x = age, ymin = lwr, ymax = upr),            # shade mean distribution
# #   #             colour = 'transparent', fill = rgb(33/255, 145/255, 140/255, 0.5))+
# #   geom_point(data = df_long, aes(x = age, y = centrality),                    # all eigenvector draws
# #              colour = rgb(253/255, 231/255, 37/255, 0.01))+
# #   geom_point(data = df_long, aes(x = age, y = mean_eigen, size = sightings),  # mean eigenvector
# #              colour = rgb(68/255, 1/255, 84/255))+
# #   # geom_line(data = mod_mu, aes(x = age, y = mid),                             # mean line
# #   #           colour = rgb(33/255, 145/255, 140/255), linewidth = 1)+
# #   scale_x_continuous('age (years)')+
# #   scale_y_continuous('eigenvector centrality (standardised)')+
# #   theme_classic()+
# #   theme(axis.text = element_text(size = 18),
# #         axis.title = element_text(size = 22),
# #         legend.text = element_text(size = 18),
# #         legend.title = element_text(size = 22))
# # ggsave(filename = '../outputs/anpshort1_nodalregression_raw.png', device = 'png',
# #        plot = last_plot())
# #
# # ## save output
# # save.image('anp_nodalregression/anpshort1_nodalregression_conditionaledge_predictions.RData')
# # dev.off()
# # rm(list = ls()[!ls() %in% c('nodal_regression','n_windows')])
# #
