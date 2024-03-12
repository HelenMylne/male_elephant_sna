#### information ####
# script takes data input from edge weight estimation for ANP population (input = 1000 draws for each of 4 chains per dyad posterior distribution)
# runs through and generates the necessary covariance matrices for each time window to assess how much memory the model may require to run

#### set up ####
#library(tidyverse) ; library(LaplacesDemon)
library(tidyverse, lib.loc = '../packages/') # library(tidyverse)
library(LaplacesDemon, lib.loc = '../packages/')

# #### window 1 ####
# load('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')
# rm(list = ls()[!ls() %in% c('edge_samples','n_windows','cdf_1')]) ; gc()
# 
# ## convert to logit scale
# logit_weights <- apply(edge_samples, 2, LaplacesDemon::logit)
# 
# ## save in a format that can be added to in a loop
# logit_edge_draws_mu <- list()
# logit_edge_draws_cov <- list()
# 
# ## fit a multivariate normal dist to the edges -- quantities will be given to Stan model as data to model joint posteriors of edge weight in the regression
# logit_edge_draws_mu[[1]] <- apply(logit_weights, 2, mean)
# logit_edge_draws_cov[[1]] <- cov(logit_weights)
# 
# ## identify older and younger of dyad
# cdf_1$age_min <- NA ; cdf_1$age_max <- NA
# for(i in 1:nrow(cdf_1)){
#   x <- c(cdf_1$age_start_1[i],cdf_1$age_start_2[i])
#   cdf_1$age_min[i] <- min(x)
#   cdf_1$age_max[i] <- max(x)
# }
# 
# ## convert node IDs to values 1:52 not numbers based on casename
# all_node_IDs <- list()
# all_node_IDs[[1]] <- unique(c(cdf_1$id_1, cdf_1$id_2))
# 
# ## calculate number of nodes per time window
# n_nodes_windows <- list()
# n_nodes_windows[[1]] <- length(unique(c(cdf_1$id_1, cdf_1$id_2)))
# 
# ## calculate number of dyads per time window
# n_dyads_windows <- list()
# n_dyads_windows[[1]] <- nrow(cdf_1)
# 
# ## create data frame of ages
# dyads_all <- cdf_1 %>% 
#   rename(node_1_original = node_1,
#          node_2_original = node_2) %>% 
#   dplyr::select(-period_start, -period_end, -dyad_rank,
#                 -bmo_1, -bmo_2, -byr_1, -byr_2, -birth_1, -birth_2,
#                 -dmo_1, -dmo_2, -dyr_1, -dyr_2, -death_1, -death_2) %>% 
#   rename(window = period) %>% 
#   mutate(dyad_window_nonrandom = paste0(dyad_id,'_',window))
# 
# #### loop through remaining windows and calculate matrices ####
# for( time_window in 2:n_windows){
#   ## load window
#   load(paste0('anp_edgecalculations/anpshort',time_window,'_edgeweights_conditionalprior.RData'))
#   rm(list = ls()[!ls() %in% c('edge_samples','cdf',
#                               'n_windows','time_window',
#                               'logit_edge_draws_mu','logit_edge_draws_cov',
#                               'dyads_all','all_node_IDs',
#                               'n_nodes_windows','n_dyads_windows')]) ; gc()
#   
#   ## convert to logit scale ####
#   logit_weights <- apply(edge_samples, 2, LaplacesDemon::logit)
#   
#   ## fit a multivariate normal dist to the edges -- quantities will be given to Stan model as data to model joint posteriors of edge weight in the regression
#   logit_edge_draws_mu[[time_window]] <- apply(logit_weights, 2, mean)
#   logit_edge_draws_cov[[time_window]] <- cov(logit_weights)
#   
#   ## identify older and younger of dyad
#   cdf$age_min <- NA ; cdf$age_max <- NA
#   for(i in 1:nrow(cdf)){
#     x <- c(cdf$age_start_1[i],cdf$age_start_2[i])
#     cdf$age_min[i] <- min(x)
#     cdf$age_max[i] <- max(x)
#   }
#   
#   ## convert node IDs to values 1:52 not numbers based on casename
#   all_node_IDs[[time_window]] <- unique(c(cdf$id_1, cdf$id_2))
#   
#   ## calculate number of nodes per time window
#   n_nodes_windows[[time_window]] <- length(unique(c(cdf$id_1, cdf$id_2)))
#   
#   ## calculate number of dyads per time window
#   n_dyads_windows[[time_window]] <- nrow(cdf)
#   
#   ## create data frame of ages
#   cdf <- cdf %>% 
#     rename(node_1_original = node_1,
#            node_2_original = node_2) %>% 
#     dplyr::select(-period_start, -period_end,
#                   -bmo_1, -bmo_2, -byr_1, -byr_2, -birth_1, -birth_2,
#                   -dmo_1, -dmo_2, -dyr_1, -dyr_2, -death_1, -death_2) %>% 
#     rename(window = period) %>% 
#     mutate(dyad_window_nonrandom = paste0(dyad_id,'_',window))
#   
#   ## append to data frame of other ages
#   dyads_all <- rbind(dyads_all, cdf)
# }
# 
# ## save image ####
# save.image('anpshort_dyadicregression_datasizecheck.RData')

#### prep data for model ####
load('anpshort_dyadicregression_datasizecheck.RData')

## global data values
n_data <- nrow(dyads_all)
n_dyads <- length(unique(dyads_all$dyad_window_nonrandom))
n_nodes <- length(unique(c(dyads_all$node_1_original,
                           dyads_all$node_2_original)))

## randomise node ID
nodes_df_rand <- data.frame(node_1_original = unique(c(dyads_all$node_1_original,
                                                       dyads_all$node_2_original))) %>% 
  mutate(node_2_original = node_1_original) %>% 
  mutate(node_1_randomised = sample(1:n_nodes, n_nodes, replace = F)) %>% 
  mutate(node_2_randomised = node_1_randomised)
dyads_all <- dyads_all %>% 
  left_join(nodes_df_rand[,c('node_1_original','node_1_randomised')],
            by = 'node_1_original') %>% 
  left_join(nodes_df_rand[,c('node_2_original','node_2_randomised')],
            by = 'node_2_original') %>% 
  mutate(node_1_window = paste0(node_1_randomised,'_',window),
         node_2_window = paste0(node_2_randomised,'_',window))

## randomise dyad ID
dyads_df_rand <- data.frame(dyad_original = unique(dyads_all$dyad_id))
dyads_df_rand$dyad_randomised <- sample(dyads_df_rand$dyad_original,
                                        nrow(dyads_df_rand),
                                        replace = F)
dyads_all <- dyads_all %>% 
  left_join(dyads_df_rand[,c('dyad_randomised','dyad_original')],
            by = 'dyad_original')

## combine window with randomised node/dyad
dyads_all$dyad_window <- paste0(dyads_all$dyad_randomised,'_',dyads_all$window)
dyads_all$node_1_window <- paste0(dyads_all$node_1_randomised,'_',dyads_all$window)
dyads_all$node_2_window <- paste0(dyads_all$node_1_randomised,'_',dyads_all$window)

## create data list 
dyad_data <- list(
  ## global data size ####
  num_data = n_data,
  num_dyads = n_dyads,
  num_nodes = n_nodes,
  num_windows = n_windows,
  ## per time window data size: dyads ####
  num_dyads_window1 = n_dyads_windows[[1]],
  num_dyads_window2 = n_dyads_windows[[2]],
  num_dyads_window3 = n_dyads_windows[[3]],
  num_dyads_window4 = n_dyads_windows[[4]],
  num_dyads_window5 = n_dyads_windows[[5]],
  num_dyads_window6 = n_dyads_windows[[6]],
  num_dyads_window7 = n_dyads_windows[[7]],
  num_dyads_window8 = n_dyads_windows[[8]],
  num_dyads_window9 = n_dyads_windows[[9]],
  num_dyads_window10 = n_dyads_windows[[10]],
  num_dyads_window11 = n_dyads_windows[[11]],
  num_dyads_window12 = n_dyads_windows[[12]],
  num_dyads_window13 = n_dyads_windows[[13]],
  num_dyads_window14 = n_dyads_windows[[14]],
  num_dyads_window15 = n_dyads_windows[[15]],
  num_dyads_window16 = n_dyads_windows[[16]],
  num_dyads_window17 = n_dyads_windows[[17]],
  num_dyads_window18 = n_dyads_windows[[18]],
  num_dyads_window19 = n_dyads_windows[[19]],
  num_dyads_window20 = n_dyads_windows[[20]],
  num_dyads_window21 = n_dyads_windows[[21]],
  num_dyads_window22 = n_dyads_windows[[22]],
  num_dyads_window23 = n_dyads_windows[[23]],
  num_dyads_window24 = n_dyads_windows[[24]],
  num_dyads_window25 = n_dyads_windows[[25]],
  num_dyads_window26 = n_dyads_windows[[26]],
  num_dyads_window27 = n_dyads_windows[[27]],
  num_dyads_window28 = n_dyads_windows[[28]],
  num_dyads_window29 = n_dyads_windows[[29]],
  num_dyads_window30 = n_dyads_windows[[30]],
  num_dyads_window31 = n_dyads_windows[[31]],
  num_dyads_window32 = n_dyads_windows[[32]],
  num_dyads_window33 = n_dyads_windows[[33]],
  num_dyads_window34 = n_dyads_windows[[34]],
  num_dyads_window35 = n_dyads_windows[[35]],
  num_dyads_window36 = n_dyads_windows[[36]],
  ## per time window data size: nodes ####
  num_nodes_window1 = n_nodes_windows[[1]],
  num_nodes_window2 = n_nodes_windows[[2]],
  num_nodes_window3 = n_nodes_windows[[3]],
  num_nodes_window4 = n_nodes_windows[[4]],
  num_nodes_window5 = n_nodes_windows[[5]],
  num_nodes_window6 = n_nodes_windows[[6]],
  num_nodes_window7 = n_nodes_windows[[7]],
  num_nodes_window8 = n_nodes_windows[[8]],
  num_nodes_window9 = n_nodes_windows[[9]],
  num_nodes_window10 = n_nodes_windows[[10]],
  num_nodes_window11 = n_nodes_windows[[11]],
  num_nodes_window12 = n_nodes_windows[[12]],
  num_nodes_window13 = n_nodes_windows[[13]],
  num_nodes_window14 = n_nodes_windows[[14]],
  num_nodes_window15 = n_nodes_windows[[15]],
  num_nodes_window16 = n_nodes_windows[[16]],
  num_nodes_window17 = n_nodes_windows[[17]],
  num_nodes_window18 = n_nodes_windows[[18]],
  num_nodes_window19 = n_nodes_windows[[19]],
  num_nodes_window20 = n_nodes_windows[[20]],
  num_nodes_window21 = n_nodes_windows[[21]],
  num_nodes_window22 = n_nodes_windows[[22]],
  num_nodes_window23 = n_nodes_windows[[23]],
  num_nodes_window24 = n_nodes_windows[[24]],
  num_nodes_window25 = n_nodes_windows[[25]],
  num_nodes_window26 = n_nodes_windows[[26]],
  num_nodes_window27 = n_nodes_windows[[27]],
  num_nodes_window28 = n_nodes_windows[[28]],
  num_nodes_window29 = n_nodes_windows[[29]],
  num_nodes_window30 = n_nodes_windows[[30]],
  num_nodes_window31 = n_nodes_windows[[31]],
  num_nodes_window32 = n_nodes_windows[[32]],
  num_nodes_window33 = n_nodes_windows[[33]],
  num_nodes_window34 = n_nodes_windows[[34]],
  num_nodes_window35 = n_nodes_windows[[35]],
  num_nodes_window36 = n_nodes_windows[[36]],
  # ## number of nodes/dyads in all preceding time windows for node age indexing ####
  # num_nodes_prev_windows = c(0,
  #                            length(which(nodes_all$window < 2)),
  #                            length(which(nodes_all$window < 3)),
  #                            length(which(nodes_all$window < 4)),
  #                            length(which(nodes_all$window < 5)),
  #                            length(which(nodes_all$window < 6)),
  #                            length(which(nodes_all$window < 7)),
  #                            length(which(nodes_all$window < 8)),
  #                            length(which(nodes_all$window < 9)),
  #                            length(which(nodes_all$window < 10)),
  #                            length(which(nodes_all$window < 11)),
  #                            length(which(nodes_all$window < 12)),
  #                            length(which(nodes_all$window < 13)),
  #                            length(which(nodes_all$window < 14)),
  #                            length(which(nodes_all$window < 15)),
  #                            length(which(nodes_all$window < 16)),
  #                            length(which(nodes_all$window < 17)),
  #                            length(which(nodes_all$window < 18)),
  #                            length(which(nodes_all$window < 19)),
  #                            length(which(nodes_all$window < 20)),
  #                            length(which(nodes_all$window < 21)),
  #                            length(which(nodes_all$window < 22)),
  #                            length(which(nodes_all$window < 23)),
  #                            length(which(nodes_all$window < 24)),
  #                            length(which(nodes_all$window < 25)),
  #                            length(which(nodes_all$window < 26)),
  #                            length(which(nodes_all$window < 27)),
  #                            length(which(nodes_all$window < 28)),
  #                            length(which(nodes_all$window < 29)),
  #                            length(which(nodes_all$window < 30)),
  #                            length(which(nodes_all$window < 31)),
  #                            length(which(nodes_all$window < 32)),
  #                            length(which(nodes_all$window < 33)),
  #                            length(which(nodes_all$window < 34)),
  #                            length(which(nodes_all$window < 35)),
  #                            length(which(nodes_all$window < 36))),
  
  ## logit_edge means per time window ####
  logit_edge_mu_1 = logit_edge_draws_mu[[1]],
  logit_edge_mu_2 = logit_edge_draws_mu[[2]],
  logit_edge_mu_3 = logit_edge_draws_mu[[3]],
  logit_edge_mu_4 = logit_edge_draws_mu[[4]],
  logit_edge_mu_5 = logit_edge_draws_mu[[5]],
  logit_edge_mu_6 = logit_edge_draws_mu[[6]],
  logit_edge_mu_7 = logit_edge_draws_mu[[7]],
  logit_edge_mu_8 = logit_edge_draws_mu[[8]],
  logit_edge_mu_9 = logit_edge_draws_mu[[9]],
  logit_edge_mu_10 = logit_edge_draws_mu[[10]],
  logit_edge_mu_11 = logit_edge_draws_mu[[11]],
  logit_edge_mu_12 = logit_edge_draws_mu[[12]],
  logit_edge_mu_13 = logit_edge_draws_mu[[13]],
  logit_edge_mu_14 = logit_edge_draws_mu[[14]],
  logit_edge_mu_15 = logit_edge_draws_mu[[15]],
  logit_edge_mu_16 = logit_edge_draws_mu[[16]],
  logit_edge_mu_17 = logit_edge_draws_mu[[17]],
  logit_edge_mu_18 = logit_edge_draws_mu[[18]],
  logit_edge_mu_19 = logit_edge_draws_mu[[19]],
  logit_edge_mu_20 = logit_edge_draws_mu[[20]],
  logit_edge_mu_21 = logit_edge_draws_mu[[21]],
  logit_edge_mu_22 = logit_edge_draws_mu[[22]],
  logit_edge_mu_23 = logit_edge_draws_mu[[23]],
  logit_edge_mu_24 = logit_edge_draws_mu[[24]],
  logit_edge_mu_25 = logit_edge_draws_mu[[25]],
  logit_edge_mu_26 = logit_edge_draws_mu[[26]],
  logit_edge_mu_27 = logit_edge_draws_mu[[27]],
  logit_edge_mu_28 = logit_edge_draws_mu[[28]],
  logit_edge_mu_29 = logit_edge_draws_mu[[29]],
  logit_edge_mu_30 = logit_edge_draws_mu[[30]],
  logit_edge_mu_31 = logit_edge_draws_mu[[31]],
  logit_edge_mu_32 = logit_edge_draws_mu[[32]],
  logit_edge_mu_33 = logit_edge_draws_mu[[33]],
  logit_edge_mu_34 = logit_edge_draws_mu[[34]],
  logit_edge_mu_35 = logit_edge_draws_mu[[35]],
  logit_edge_mu_36 = logit_edge_draws_mu[[36]],
  ## logit_edge covariance matrix per time window ####
  logit_edge_cov_1 = logit_edge_draws_cov[[1]],
  logit_edge_cov_2 = logit_edge_draws_cov[[2]],
  logit_edge_cov_3 = logit_edge_draws_cov[[3]],
  logit_edge_cov_4 = logit_edge_draws_cov[[4]],
  logit_edge_cov_5 = logit_edge_draws_cov[[5]],
  logit_edge_cov_6 = logit_edge_draws_cov[[6]],
  logit_edge_cov_7 = logit_edge_draws_cov[[7]],
  logit_edge_cov_8 = logit_edge_draws_cov[[8]],
  logit_edge_cov_9 = logit_edge_draws_cov[[9]],
  logit_edge_cov_10 = logit_edge_draws_cov[[10]],
  logit_edge_cov_11 = logit_edge_draws_cov[[11]],
  logit_edge_cov_12 = logit_edge_draws_cov[[12]],
  logit_edge_cov_13 = logit_edge_draws_cov[[13]],
  logit_edge_cov_14 = logit_edge_draws_cov[[14]],
  logit_edge_cov_15 = logit_edge_draws_cov[[15]],
  logit_edge_cov_16 = logit_edge_draws_cov[[16]],
  logit_edge_cov_17 = logit_edge_draws_cov[[17]],
  logit_edge_cov_18 = logit_edge_draws_cov[[18]],
  logit_edge_cov_19 = logit_edge_draws_cov[[19]],
  logit_edge_cov_20 = logit_edge_draws_cov[[20]],
  logit_edge_cov_21 = logit_edge_draws_cov[[21]],
  logit_edge_cov_22 = logit_edge_draws_cov[[22]],
  logit_edge_cov_23 = logit_edge_draws_cov[[23]],
  logit_edge_cov_24 = logit_edge_draws_cov[[24]],
  logit_edge_cov_25 = logit_edge_draws_cov[[25]],
  logit_edge_cov_26 = logit_edge_draws_cov[[26]],
  logit_edge_cov_27 = logit_edge_draws_cov[[27]],
  logit_edge_cov_28 = logit_edge_draws_cov[[28]],
  logit_edge_cov_29 = logit_edge_draws_cov[[29]],
  logit_edge_cov_30 = logit_edge_draws_cov[[30]],
  logit_edge_cov_31 = logit_edge_draws_cov[[31]],
  logit_edge_cov_32 = logit_edge_draws_cov[[32]],
  logit_edge_cov_33 = logit_edge_draws_cov[[33]],
  logit_edge_cov_34 = logit_edge_draws_cov[[34]],
  logit_edge_cov_35 = logit_edge_draws_cov[[35]],
  logit_edge_cov_36 = logit_edge_draws_cov[[36]],
  ## exposure variables ####
  age_min = dyads_all$age_min,
  age_max = dyads_all$age_max,
  ## node IDs unique to time windows ####
  node_1_window = dyads_all$node_1_window,             # node IDs for multimembership effects
  node_2_window = dyads_all$node_2_window,              # node IDs for multimembership effects
  ## node IDs across time windows ####
  node_1_id = dyads_all$node_1_randomised,             # node IDs for multimembership effects
  node_2_id = dyads_all$node_2_randomised,             # node IDs for multimembership effects
  ## dyad IDs unique to time windows ####
  dyad_window = dyads_all$dyad_window,
  ## dyad IDs across time windows ####
  dyad_id = dyads_all$dyad_randomised,
)

## save image ####
save.image('anpshort_dyadicregression_datasizecheck.RData')
