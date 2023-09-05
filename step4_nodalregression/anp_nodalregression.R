#### Information ####
# 1) to estimate the ages of individuals in other populations, first came up with a survival curve for the Amboseli elephants which have good birth and death data.
# 2) selected the Gompertz bathtub distribution
# 3) use distribution to convert MOTNP ages in categories to probability distributions of true age (first part of this script)
# 4) extract node centrality distribution from draws of edge weight
# 5) use MOTNP age distributions to predict network centrality (third part of this script)

#### set up ####
options(future.globals.maxSize = 10000*(1024^2))   # running model with full age distributions = 5.01GB of globals, which exceeds default maximum allowed size. Set to allow up to 10 GB to allow model to run

# library(tidyverse) ; library(cmdstanr) ; library(brms) ; library(Rcpp) ; library(ggdist) ; library(posterior) ; library(bayesplot) ; library(igraph) ; library(LaplacesDemon) ; library(bisonR) ; library(janitor)
library(tidyverse, lib.loc = '../packages/')       # library(tidyverse)
library(cmdstanr, lib.loc = '../packages/')        # library(cmdstanr)
library(brms, lib.loc = '../packages/')            # library(brms)
library(Rcpp, lib.loc = '../packages/')            # library(Rcpp)
library(ggdist, lib.loc = '../packages/')          # library(ggdist)
library(posterior, lib.loc = '../packages/')       # library(posterior)
library(bayesplot, lib.loc = '../packages/')       # library(bayesplot)
library(igraph, lib.loc = '../packages/')          # library(igraph)
library(LaplacesDemon, lib.loc = '../packages/')   # library(LaplacesDemon)
library(bisonR, lib.loc = '../packages/')          # library(bisonR)
library(janitor, lib.loc = '../packages/')         # library(janitor)

## set cmdstan path
set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0/')

## load work space from calculating edge weights
load('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')

#### prior predictive check ####
## define PDF output
pdf('../outputs/anpshort1_nodalregression_conditionalprior_bisonR.pdf')

## simulate
age <- 1:60
beta_mu <- 0
beta_sigma <- 0.005
mean_age <- mean(nodes$age)
plot(NULL, xlim = c(10,60), ylim = c(0,1), las = 1,
     xlab = 'age', ylab = 'eigenvector centrality')
for(i in 1:100){
  intercept <- rbeta(1,2,2)   # is this right?? I've gone for a symmetrical one here that in itself explores most of the parameter space and allows it to see whether some of the lines are steep enough to go from top to bottom, but on the assumption that when combined, they will explore only the space relevant to their starting position
  beta <- rnorm(1, beta_mu, beta_sigma)
  lines(x = age, y = intercept + (age - mean_age)*beta, col = rgb(0,0,1,0.5)) # vast majority come out somewhere sensible, and those that don't would if they started at a different value for age 10 so fine in combo with posterior intercept
}

#### nodal regression -- bisonR, all time windows ####
# load edge weight model and data frames
for( time_window in 1:36) {
  
  ## load work space from calculating edge weights
  load(file = paste0('anp_edgecalculations/anpshort',time_window,'_bisonr_edgescalculated.RData'))
  rm(counts_df_model, counts_df_allwindows, df, edgelist, ew_chain, ew_edgelist, global_cv, anp_edges_null, draw98, i, ids, median98, sri98, plot_network_threshold, plot_network_threshold2) ; gc()
  
  ## create nodal data frame from counts_df (data frame of dyad counts of sightings together vs apart)
  df_nodal <- counts_df[,c('node_1','id_1')] %>% distinct()
  colnames(df_nodal) <- c('node_2','id_2')
  df_nodal <- rbind(df_nodal, counts_df[,c('node_2','id_2')]) %>% distinct()
  colnames(df_nodal) <- c('node','id')
  colnames(nodes)[1] <- 'node'
  df_nodal <- left_join(df_nodal, nodes, by = 'node')
  
  # define PDF output
  pdf(paste0('../outputs/anpshort',time_window,'_nodalregression_plots_meanage.pdf'))

  ## eigenvector only ####
  anp_edge_weights$model_data$prior_fixed_mu    <- beta_mu
  anp_edge_weights$model_data$prior_fixed_sigma <- beta_sigma
  
  ## run model ####
  anp_eigen <- bison_brm(
    bison(node_eigen(node)) ~ age,
    anp_edge_weights,
    df_nodal,
    chains = 4,
    iter = 10000,
    thin = 2
  )
  summary(anp_eigen) # FIT 100 IMPUTED MODELS, 4 CHAINS PER MODEL, EACH 1000 DRAWS LONG (+1000 DRAWS WARM UP). WARNING AT END OF MODEL RUN THAT CHAINS <3 DRAWS LONG AS ACTUAL CHAIN IS ONLY 1 WARMUP AND 1 SAMPLE. ONLY IMPUTED CHAINS FOLLOW THE SPECIFIED ITERATIONS AND THINNING
  
  ## save output
  save.image(paste0('anp_edgecalculations/anpshort',time_window,'_nodalregression.RData'))
  
  ## posterior check ####
  summary(anp_eigen)
  
  # plot
  eigen_values <- anp_eigen$data
  plot(eigen_values$bison_node_eigen ~ eigen_values$age, 
       las = 1, pch = 19, col = rgb(0,0,1,0.2),
       xlab = 'mean age estimate', ylab = 'eigenvector centrality',
       main = 'effect of age on eigenvector centrality')
  
  hist(anp_eigen$rhats[,2], las = 1, main = 'Rhat values for 100 imputed model runs', xlab = 'Rhat')
  
  post_eigen <- as.data.frame(as_draws_df(anp_eigen)) %>% clean_names()
  plot(data = post_eigen[post_eigen$chain == 1,], b_age ~ draw, type = 'l', xlim = c(0,10000))
  lines(data = post_eigen[post_eigen$chain == 2,], b_age ~ draw, col = 'red')
  lines(data = post_eigen[post_eigen$chain == 3,], b_age ~ draw, col = 'blue')
  lines(data = post_eigen[post_eigen$chain == 4,], b_age ~ draw, col = 'green')
  lines(data = post_eigen[post_eigen$chain == 5,], b_age ~ draw, col = 'purple')
  unique(post_eigen$chain)
  
  post_eigen$imputation <- rep(1:100, each = 10000)
  post_eigen$chain_imp  <- as.factor(rep(1:4, each = length(post_eigen$chain[post_eigen$chain == 1])))
  ggplot(post_eigen, aes(x = iteration, y = b_age, colour = chain_imp))+
    geom_line()+
    facet_wrap(imputation ~ .)   # looks good, well mixed
  
  hist(post_eigen$b_age)         # natural scale?
  #hist(plogis(post_eigen$b_age)) # logit scale?
  
  #plot(anp_eigen)
  
  plot(conditional_effects(anp_eigen), points = TRUE) # older elephants have lower network centrality than younger
  
  if(length(which(eigen_values$age != nodes$age)) == 0 ){
    data <- cbind(eigen_values, nodes[,c('node','sightings')])
    
    ggplot(data = data, aes(x = age, y = bison_node_eigen,
                            colour = age))+
      geom_point(aes(size = sightings))+
      geom_smooth()+
      scale_x_continuous(name = 'age')+
      scale_y_continuous(name = 'mean eigenvector centrality')+
      scale_colour_viridis_c()+
      theme_classic()+
      theme(legend.position = 'none', axis.title = element_text(size = 18),
            axis.text = element_text(size = 14))
  }
  
  # compare empirical distribution to posterior predictive distribution
  y <- eigen_values$bison_node_eigen                         # extract eigenvector centralities
  yrep <- posterior_predict(anp_eigen, draws = 500)          # make predictions of eigenvector centtrality
  dim(yrep)
  ppc_dens_overlay(y, yrep[1:1000, ])                        # plot 1000 predictions over empirical distribution (is it ok that this slightly extends beyond x = 1??)
  ppc_hist(y, yrep[1:55, ])                                  # compare 55 predictions to empirical distribution
  ppc_pit_ecdf(y, yrep[1:50,])                               # no idea....
  ppc_ecdf_overlay(y, yrep[1:50,]) + xaxis_text()            # also no idea...
  
  ## save output
  save.image(paste0('anp_edgecalculations/anpshort',time_window,'_nodalregression.RData'))
  dev.off()
  
}

#### nodal regression -- Stan ####
### extract eigenvector centralities -- time window 1 ####
pdf('../outputs/anpshort1_nodalregression_extracteigen.pdf')
cdf_1 <- cdf_1 %>% mutate(sri = event_count / period_count_dyad)
ele_ids <- unique(c(cdf_1$id_1, cdf_1$id_2))
n_eles <- length(ele_ids)
edges$chain <- ifelse(edges$chain == 'chain1', 1,
                      ifelse(edges$chain == 'chain2', 2,
                             ifelse(edges$chain == 'chain3', 3, 4)))
edges$draw_id <- edges$position + (edges$chain-1) * 1000
edges <- left_join(edges, cdf_1[,c('dyad_id','node_1','node_2','id_1','id_2')], by = 'dyad_id')

adj_tensor <- array(NA, c(length(unique(cdf_1$id_1))+1,
                         length(unique(cdf_1$id_2))+1,
                         nrow(edge_samples)),
                    dimnames = list(ele_ids, ele_ids, NULL))
for (i in 1:n_dyads) {
  dyad_row <- cdf_1[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- edge_samples[,i]
}
adj_tensor[,,1]

## create array for eigen values to be saved into
eigen <- array(data = NA, dim = c(n_eles, 4, n_samples*n_chains),
               dimnames = list(nodes$node,
                               c('node','age','sightings','eigenvector'),
                               1:(n_samples*n_chains)))
eigen[,1,] <- nodes$node
eigen[,2,] <- nodes$age
eigen[,3,] <- nodes$sightings

## fill array
for(draw in 1:(n_samples*n_chains)){
  network <- graph_from_adjacency_matrix(adjmatrix = adj_tensor[,,draw],
                                         diag = FALSE, mode = 'undirected', weighted = TRUE)
  eigen_values <- as.matrix(igraph::eigen_centrality(network, directed = FALSE)$vector)
  eigen[,4,draw] <- eigen_values[,1]
}

### plots to check outputs -- time window 1 ####
# check age against sightings
plot(nodes$sightings ~ nodes$age,
     las = 1, xlab = 'age', ylab = 'sightings',
     main = 'time window 1 (short):\nsightings ~ age',
     pch = 19, col = 'purple')

# check average group size against sightings
obs <- read_csv('../data_processed/step1_dataprocessing/anp_sightings_rawcombined.csv') %>% 
  filter(obs_date >= periods[1] & obs_date < periods[2])
nodes$grp_sze_avg <- NA ; for(i in 1:nrow(nodes)){
  ele_obs <- obs %>% filter(id == nodes$id[i])
  nodes$grp_sze_avg[i] <- mean(ele_obs$grp_size)
}
plot(nodes$grp_sze_avg ~ nodes$sightings,
     las = 1, xlab = 'sightings', ylab = 'average group size',
     main = 'time window 1 (short):\naverage group size ~ sightings',
     pch = 19, col = 'purple')

# check average group size against age
plot(nodes$grp_sze_avg ~ nodes$age,
     las = 1, xlab = 'average group size', ylab = 'age',
     main = 'time window 1 (short):\naverage group size ~ age',
     pch = 19, col = 'purple')

# check SRI eigenvector against sightings
adj_mat <- matrix(NA,nrow = length(ele_ids), ncol = length(ele_ids),
                    dimnames = list(ele_ids, ele_ids))
for (i in 1:n_dyads) {
  dyad_row <- cdf_1[i, ]
  adj_mat[dyad_row$id_1, dyad_row$id_2] <- dyad_row$sri[1]
}
network <- graph_from_adjacency_matrix(adjmatrix = adj_mat,
                                       diag = FALSE, mode = 'undirected', weighted = TRUE)
eigen_values <- as.data.frame(eigen_centrality(network, directed = FALSE)$vector) %>% 
  rename(eigen_sri = `eigen_centrality(network, directed = FALSE)$vector`)
eigen_values$id <- rownames(eigen_values)
nodes <- left_join(nodes, eigen_values, by = 'id')

plot(nodes$eigen_sri ~ nodes$sightings, ylim = c(0,1),
     las = 1, xlab = 'sightings', ylab = 'eigenvector',
     main = 'time window 1 (short):\nSRI eigenvector ~ sightings',
     pch = 19, col = 'purple')

# check average group size against SRI eigenvector
plot(nodes$eigen_sri ~ nodes$grp_sze_avg,
     las = 1, xlab = 'average group size', ylab = 'eigenvector',
     main = 'time window 1 (short):\nSRI eigenvector ~ average group size',
     pch = 19, col = 'purple')

# check BISoN eigenvector against sightings
plot(NULL, xlim = c(0,max(nodes$sightings)), ylim = c(0,1),
     las = 1, xlab = 'sightings', ylab = 'eigenvector',
     main = 'time window 1 (short):\nBISoN eigenvector ~ sightings')
for(i in 1:n_samples){
  points(eigen[,4,i] ~ eigen[,3,i], pch = 19, col = rgb(1,1,0,0.1))
}
for(i in 1:n_eles){
  x <- eigen[i,,]
  points(mean(x[4,]) ~ x[3,1], pch = 19, col = 'purple')
}

# check average group size against BISoN eigenvector
eigen2 <- array(data = NA, dim = c(n_eles, 5, n_samples*n_chains),
                dimnames = list(ele_ids, c(colnames(eigen[,,1]),'grp_sze_avg'), NULL))
eigen2[,1:4,] <- eigen[,1:4,] ; eigen2[,5,] <- nodes$grp_sze_avg
eigen <- eigen2 ; rm(eigen2) ; gc()
plot(NULL, xlim = c(0,max(nodes$grp_sze_avg)), ylim = c(0,1),
     las = 1, xlab = 'average group size', ylab = 'eigenvector',
     main = 'time window 1 (short):\nBISoN eigenvector ~ average group size')
for(i in 1:n_samples){
  points(eigen[,4,i] ~ eigen[,5,i], pch = 19, col = rgb(1,1,0,0.1))
}
for(i in 1:n_eles){
  x <- eigen[i,,]
  points(mean(x[4,]) ~ x[5,1], pch = 19, col = 'purple')
}

## save workspace for future
rm(dyad_row, edge_binary, edgelist, eigen_values, network, adj_tensor, draw, i) ; gc()
save.image('anp_nodalregression/anpshort1_nodalregression_conditionaledge.RData')
dev.off()

### extract eigenvector centralities -- time window 2-36 short ####
rm(list = ls()[! ls() %in% 'n_windows']) ; gc()
for(time_window in 2:n_windows) {
  ### extract data ####
  load(file = paste0('anp_edgecalculations/anpshort', time_window, '_edgeweights_conditionalprior.RData'))
  pdf(file = paste0('../outputs/anpshort',time_window,'_nodalregression_extracteigen.pdf'))
  cdf <- cdf %>% mutate(sri = event_count / period_count_dyad)
  
  # extract edge weights
  ele_ids <- unique(c(cdf$id_1, cdf$id_2))
  n_eles <- length(ele_ids)
  edges$chain <- ifelse(edges$chain == 'chain1', 1,
                        ifelse(edges$chain == 'chain2', 2,
                               ifelse(edges$chain == 'chain3', 3, 4)))
  edges$draw_id <- edges$position + (edges$chain-1) * 1000
  edges <- edges %>% 
    rename(dyad_id = dyad) %>% 
    left_join(cdf[,c('dyad_id','node_1','node_2','id_1','id_2')],
              by = 'dyad_id')
  
  # set up array to store edge weights
  adj_tensor <- array(NA, c(length(unique(cdf$id_1))+1,
                            length(unique(cdf$id_2))+1,
                            nrow(edge_samples)),
                      dimnames = list(ele_ids, ele_ids, NULL))
  
  # convert edge weight format to array
  for (i in 1:n_dyads) {
    dyad_row <- cdf[i, ]
    adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- edge_samples[,i]
  }
  adj_tensor[,,1]
  
  ## create array for eigen values to be saved into: 4000 layers of data frames, where each data frame contains node, age, sightings, and one draw from edge model
  eigen <- array(data = NA, dim = c(n_eles, 4, n_samples*n_chains),
                 dimnames = list(nodes$node,
                                 c('node','age','sightings','eigenvector'),
                                 1:(n_samples*n_chains)))
  eigen[,1,] <- nodes$node
  eigen[,2,] <- nodes$age
  eigen[,3,] <- nodes$sightings
  
  ## fill array: create network for each set of draws, extract eigen centrality per network, fill data frame 'eigenvector' column with values per network
  for(draw in 1:(n_samples*n_chains)){
    network <- graph_from_adjacency_matrix(adjmatrix = adj_tensor[,,draw],
                                           diag = FALSE,
                                           mode = 'undirected', weighted = TRUE)
    eigen_values <- as.matrix(eigen_centrality(network, directed = FALSE)$vector)
    eigen[,4,draw] <- eigen_values[,1]
  }
  
  ### plots to check outputs ####
  ## check age against sightings
  plot(nodes$sightings ~ nodes$age,
       las = 1, xlab = 'age', ylab = 'sightings',
       main = paste0('time window ',time_window,' (short):\nsightings ~ age'),
       pch = 19, col = 'purple')
  
  ## check average group size against sightings
  obs <- read_csv('../data_processed/step1_dataprocessing/anp_sightings_rawcombined.csv') %>% 
    filter(obs_date >= periods[time_window] & obs_date < periods[time_window+1])
  nodes$grp_sze_avg <- NA ; for(i in 1:nrow(nodes)){
    ele_obs <- obs %>% filter(id == nodes$id[i])
    nodes$grp_sze_avg[i] <- mean(ele_obs$grp_size)
  }
  plot(nodes$grp_sze_avg ~ nodes$sightings,
       las = 1, xlab = 'sightings', ylab = 'average group size',
       main = paste0('time window ',time_window,' (short):\naverage group size ~ sightings'),
       pch = 19, col = 'purple')
  
  ## check average group size against age
  plot(nodes$grp_sze_avg ~ nodes$age,
       las = 1, xlab = 'age', ylab = 'average group size',
       main = paste0('time window ',time_window,' (short):\naverage group size ~ age'),
       pch = 19, col = 'purple')
  
  ## check SRI eigenvector against sightings
  adj_mat <- matrix(NA,nrow = length(ele_ids), ncol = length(ele_ids),
                    dimnames = list(ele_ids, ele_ids))
  for (i in 1:n_dyads) {
    dyad_row <- cdf[i, ]
    adj_mat[dyad_row$id_1, dyad_row$id_2] <- dyad_row$sri[1]
  }
  network <- graph_from_adjacency_matrix(adjmatrix = adj_mat,
                                         diag = FALSE, mode = 'undirected', weighted = TRUE)
  eigen_values <- as.data.frame(eigen_centrality(network, directed = FALSE)$vector) %>% 
    rename(eigen_sri = `eigen_centrality(network, directed = FALSE)$vector`)
  eigen_values$id <- rownames(eigen_values)
  nodes <- left_join(nodes, eigen_values, by = 'id')
  
  plot(nodes$eigen_sri ~ nodes$sightings, ylim = c(0,1),
       las = 1, xlab = 'sightings', ylab = 'eigenvector',
       main = paste0('time window ',time_window,' (short):\nSRI eigenvector ~ sightings'),
       pch = 19, col = 'purple')
  
  ## check average group size against SRI eigenvector
  plot(nodes$eigen_sri ~ nodes$grp_sze_avg,
       las = 1, xlab = 'average group size', ylab = 'eigenvector',
       main = paste0('time window ',time_window,' (short):\nSRI eigenvector ~ average group size'),
       pch = 19, col = 'purple')
  
  ## check BISoN eigenvector against sightings
  plot(NULL, xlim = c(0,max(nodes$sightings)), ylim = c(0,1),
       las = 1, xlab = 'sightings', ylab = 'eigenvector',
       main = paste0('time window ',time_window,'  (short):\nBISoN eigenvector ~ sightings'))
  for(i in 1:n_samples){
    points(eigen[,4,i] ~ eigen[,3,i], pch = 19, col = rgb(1,1,0,0.1))
  }
  for(i in 1:n_eles){
    x <- eigen[i,,]
    points(mean(x[4,]) ~ x[3,1], pch = 19, col = 'purple')
  }
 
  ## check average group size against BISoN eigenvector
  eigen2 <- array(data = NA, dim = c(n_eles, 5, n_samples*n_chains),
                  dimnames = list(ele_ids, c(colnames(eigen[,,1]),'grp_sze_avg'), NULL))
  eigen2[,1:4,] <- eigen[,1:4,] ; eigen2[,5,] <- nodes$grp_sze_avg
  eigen <- eigen2 ; rm(eigen2) ; gc()
  plot(NULL, xlim = c(0,max(nodes$grp_sze_avg)), ylim = c(0,1),
       las = 1, xlab = 'average group size', ylab = 'eigenvector',
       main = paste0('time window ',time_window,' (short):\nBISoN eigenvector ~ average group size'))
  for(i in 1:n_samples){
    points(eigen[,4,i] ~ eigen[,5,i], pch = 19, col = rgb(1,1,0,0.1))
  }
  for(i in 1:n_eles){
    x <- eigen[i,,]
    points(mean(x[4,]) ~ x[5,1], pch = 19, col = 'purple')
  }
  
  ## save workspace for future
  rm(dyad_row, edge_binary, edgelist, eigen_values, network, adj_tensor, draw, i) ; gc()
  save.image(file = paste0('anp_nodalregression/anpshort',time_window,'_nodalregression_conditionaledge.RData'))
  dev.off()
  
  ## clear workspace and set up for next loop
  rm(list = ls()[! ls() %in% c('n_windows', 'time_window')])
  
  ## print progress marker
  print(paste0('time window ', time_window, ' finished at ', Sys.time()))
}

### extract eigenvector centralities -- time window 1-7 long ####
rm(list = ls()) ; gc()
n_windows <- read_csv('../data_processed/step1_dataprocessing/anp_bayesian_pairwiseevents_aggregated_allperiods_longperiods_impossiblepairsremoved.csv') ; n_windows <- length(unique(n_windows$period))
for(time_window in 1:n_windows) {
  ### extract data ####
  load(file = paste0('anp_edgecalculations/anplong', time_window, '_edgeweights_conditionalprior.RData'))
  pdf(file = paste0('../outputs/anplong',time_window,'_nodalregression_extracteigen.pdf'))
  cdf <- cdf %>% mutate(sri = event_count / count_period_dyad)
  
  # extract edge weights
  ele_ids <- unique(c(cdf$id_1, cdf$id_2))
  n_eles <- length(ele_ids)
  edges$chain <- ifelse(edges$chain == 'chain1', 1,
                        ifelse(edges$chain == 'chain2', 2,
                               ifelse(edges$chain == 'chain3', 3, 4)))
  edges$draw_id <- edges$position + (edges$chain-1) * 1000
  edges <- edges %>% 
    rename(dyad_id = dyad) %>% 
    left_join(cdf[,c('dyad_id','node_1','node_2','id_1','id_2')],
              by = 'dyad_id')
  
  # set up array to store edge weights
  adj_tensor <- array(NA, c(length(unique(cdf$id_1))+1,
                            length(unique(cdf$id_2))+1,
                            nrow(edge_samples)),
                      dimnames = list(ele_ids, ele_ids, NULL))
  
  # convert edge weight format to array
  for (i in 1:n_dyads) {
    dyad_row <- cdf[i, ]
    adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- edge_samples[,i]
  }
  adj_tensor[,,1]
  
  ## create array for eigen values to be saved into: 4000 layers of data frames, where each data frame contains node, age, sightings, and one draw from edge model
  eigen <- array(data = NA, dim = c(n_eles, 4, n_samples*n_chains),
                 dimnames = list(nodes$node,
                                 c('node','age','sightings','eigenvector'),
                                 1:(n_samples*n_chains)))
  eigen[,1,] <- nodes$node
  eigen[,2,] <- nodes$age
  eigen[,3,] <- nodes$sightings
  
  ## fill array: create network for each set of draws, extract eigen centrality per network, fill data frame 'eigenvector' column with values per network
  for(draw in 1:(n_samples*n_chains)){
    network <- graph_from_adjacency_matrix(adjmatrix = adj_tensor[,,draw],
                                           diag = FALSE,
                                           mode = 'undirected', weighted = TRUE)
    eigen_values <- as.matrix(eigen_centrality(network, directed = FALSE)$vector)
    eigen[,4,draw] <- eigen_values[,1]
  }
  
  ### plots to check outputs ####
  ## check age against sightings
  plot(nodes$sightings ~ nodes$age,
       las = 1, xlab = 'age', ylab = 'sightings',
       main = paste0('time window ',time_window,' (long):\nsightings ~ age'),
       pch = 19, col = 'purple')
  
  ## check average group size against sightings
  obs <- read_csv('../data_processed/step1_dataprocessing/anp_sightings_rawcombined.csv') %>% 
    filter(obs_date >= periods[time_window] & obs_date < periods[time_window+1])
  nodes$grp_sze_avg <- NA ; for(i in 1:nrow(nodes)){
    ele_obs <- obs %>% filter(id == nodes$id[i])
    nodes$grp_sze_avg[i] <- mean(ele_obs$grp_size)
  }
  plot(nodes$grp_sze_avg ~ nodes$sightings,
       las = 1, xlab = 'sightings', ylab = 'average group size',
       main = paste0('time window ',time_window,' (long):\naverage group size ~ sightings'),
       pch = 19, col = 'purple')
  
  ## check average group size against age
  plot(nodes$grp_sze_avg ~ nodes$age,
       las = 1, xlab = 'age', ylab = 'average group size',
       main = paste0('time window ',time_window,' (long):\naverage group size ~ age'),
       pch = 19, col = 'purple')
  
  ## check SRI eigenvector against sightings
  adj_mat <- matrix(NA,nrow = length(ele_ids), ncol = length(ele_ids),
                    dimnames = list(ele_ids, ele_ids))
  for (i in 1:n_dyads) {
    dyad_row <- cdf[i, ]
    adj_mat[dyad_row$id_1, dyad_row$id_2] <- dyad_row$sri[1]
  }
  network <- graph_from_adjacency_matrix(adjmatrix = adj_mat,
                                         diag = FALSE, mode = 'undirected', weighted = TRUE)
  eigen_values <- as.data.frame(eigen_centrality(network, directed = FALSE)$vector) %>% 
    rename(eigen_sri = `eigen_centrality(network, directed = FALSE)$vector`)
  eigen_values$id <- rownames(eigen_values)
  nodes <- left_join(nodes, eigen_values, by = 'id')
  
  plot(nodes$eigen_sri ~ nodes$sightings, ylim = c(0,1),
       las = 1, xlab = 'sightings', ylab = 'eigenvector',
       main = paste0('time window ',time_window,' (long):\nSRI eigenvector ~ sightings'),
       pch = 19, col = 'purple')
  
  ## check average group size against SRI eigenvector
  plot(nodes$eigen_sri ~ nodes$grp_sze_avg,
       las = 1, xlab = 'average group size', ylab = 'eigenvector',
       main = paste0('time window ',time_window,' (long):\nSRI eigenvector ~ average group size'),
       pch = 19, col = 'purple')
  
  ## check BISoN eigenvector against sightings
  plot(NULL, xlim = c(0,max(nodes$sightings)), ylim = c(0,1),
       las = 1, xlab = 'sightings', ylab = 'eigenvector',
       main = paste0('time window ',time_window,'  (long):\nBISoN eigenvector ~ sightings'))
  for(i in 1:n_samples){
    points(eigen[,4,i] ~ eigen[,3,i], pch = 19, col = rgb(1,1,0,0.1))
  }
  for(i in 1:n_eles){
    x <- eigen[i,,]
    points(mean(x[4,]) ~ x[3,1], pch = 19, col = 'purple')
  }
  
  ## check average group size against BISoN eigenvector
  eigen2 <- array(data = NA, dim = c(n_eles, 5, n_samples*n_chains),
                  dimnames = list(ele_ids, c(colnames(eigen[,,1]),'grp_sze_avg'), NULL))
  eigen2[,1:4,] <- eigen[,1:4,] ; eigen2[,5,] <- nodes$grp_sze_avg
  eigen <- eigen2 ; rm(eigen2) ; gc()
  plot(NULL, xlim = c(0,max(nodes$grp_sze_avg)), ylim = c(0,1),
       las = 1, xlab = 'average group size', ylab = 'eigenvector',
       main = paste0('time window ',time_window,' (long):\nBISoN eigenvector ~ average group size'))
  for(i in 1:n_samples){
    points(eigen[,4,i] ~ eigen[,5,i], pch = 19, col = rgb(1,1,0,0.1))
  }
  for(i in 1:n_eles){
    x <- eigen[i,,]
    points(mean(x[4,]) ~ x[5,1], pch = 19, col = 'purple')
  }
  
  ## save workspace for future
  rm(dyad_row, edge_binary, edgelist, eigen_values, network, adj_tensor, draw, i) ; gc()
  save.image(file = paste0('anp_nodalregression/anplong',time_window,'_nodalregression_conditionaledge.RData'))
  dev.off()
  
  ## clear workspace and set up for next loop
  rm(list = ls()[! ls() %in% c('n_windows', 'time_window')])
  
  ## print progress marker
  print(paste0('time window ', time_window, ' finished at ', Sys.time()))
}

for(time_window in 1:n_windows){
  load(paste0('anp_nodalregression/anplong',time_window,'_nodalregression_conditionaledge.RData'))
  
  ## convert eigen to single data frame
  eigen_long <- eigen[,,1] %>% 
    as.data.frame() %>% 
    mutate(draw = 1)
  for(i in 2:(n_samples*n_chains)){
    eigen_draw <- eigen[,,i] %>% 
      as.data.frame() %>% 
      mutate(draw = i)
    eigen_long <- rbind(eigen_long, eigen_draw)
    if(i %% 100 == 0) { print(i) }
  }
  
  ## save eigenvector estimates
  saveRDS(eigen_long, file = paste0('../data_processed/step4_nodalregression/anplong',time_window,'eigenvectorestimates.RDS'))
  rm(list = ls()[! ls() %in% c('n_windows','time_window')])
  
  ## add time marker
  print(paste0('time window ',time_window, ' complete at ', Sys.time()))
}

n_windows <- 36
for(time_window in 1:n_windows){
  load(paste0('anp_nodalregression/anpshort',time_window,'_nodalregression_conditionaledge.RData'))
  
  ## convert eigen to single data frame
  eigen_long <- eigen[,,1] %>% 
    as.data.frame() %>% 
    mutate(draw = 1)
  for(i in 2:(n_samples*n_chains)){
    eigen_draw <- eigen[,,i] %>% 
      as.data.frame() %>% 
      mutate(draw = i)
    eigen_long <- rbind(eigen_long, eigen_draw)
    if(i %% 500 == 0) { print(i) }
  }
  
  ## save eigenvector estimates
  saveRDS(eigen_long, file = paste0('../data_processed/step4_nodalregression/anpshort',time_window,'eigenvectorestimates.RDS'))
  rm(list = ls()[! ls() %in% c('n_windows','time_window')])
  
  ## add time marker
  print(paste0('time window ',time_window, ' complete at ', Sys.time()))
}
rm(list = ls()) ; gc()

## run model -- brms multiple ####
## create data list
eigen_list <- list()
for(i in 1:(n_chains*n_samples)){
  eigen_list[[i]] <- eigen[,,i]
}

# set prior
prior_slope <- set_prior(prior = 'normal(0, 0.005)')
#prior_intercept <- set_prior(prior = 'beta(2, 2)')

## run model
anp_eigen <- brm_multiple(
  formula = bf(eigenvector ~ age),
  data = eigen_list,
  chains = 4,
  iter = 10000,
  thin = 2,
  prior = prior_slope
)

summary(anp_eigen)

## save output
save.image('anpshort1_nodalregression_conditionaledge_brmsmultiple.RData')

## posterior check ####
# plot
post_eigen <- as.data.frame(as_draws_df(anp_eigen)) %>% clean_names()
plot(data = post_eigen[post_eigen$chain == 1,], b_age ~ draw, type = 'l', xlim = c(0,10000))
lines(data = post_eigen[post_eigen$chain == 2,], b_age ~ draw, col = 'red')
lines(data = post_eigen[post_eigen$chain == 3,], b_age ~ draw, col = 'blue')
lines(data = post_eigen[post_eigen$chain == 4,], b_age ~ draw, col = 'green')
lines(data = post_eigen[post_eigen$chain == 5,], b_age ~ draw, col = 'purple')
unique(post_eigen$chain)

post_eigen$imputation <- rep(1:100, each = 10000)
post_eigen$chain_imp  <- as.factor(rep(1:4, each = length(post_eigen$chain[post_eigen$chain == 1])))
ggplot(post_eigen, aes(x = iteration, y = b_age, colour = chain_imp))+
  geom_line()+
  facet_wrap(imputation ~ .)   # looks good, well mixed

hist(post_eigen$b_age)         # natural scale?
#hist(plogis(post_eigen$b_age)) # logit scale?

#plot(anp_eigen)

plot(conditional_effects(anp_eigen), points = TRUE) # older elephants have lower network centrality than younger

if(length(which(eigen_values$age != nodes$age)) == 0 ){
  data <- cbind(eigen_values, nodes[,c('node','sightings')])
  
  ggplot(data = data, aes(x = age, y = bison_node_eigen,
                          colour = age))+
    geom_point(aes(size = sightings))+
    geom_smooth()+
    scale_x_continuous(name = 'age')+
    scale_y_continuous(name = 'mean eigenvector centrality')+
    scale_colour_viridis_c()+
    theme_classic()+
    theme(legend.position = 'none', axis.title = element_text(size = 18),
          axis.text = element_text(size = 14))
}

# compare empirical distribution to posterior predictive distribution
y <- eigen_values$bison_node_eigen                         # extract eigenvector centralities
yrep <- posterior_predict(anp_eigen, draws = 500)          # make predictions of eigenvector centtrality
dim(yrep)
ppc_dens_overlay(y, yrep[1:1000, ])                        # plot 1000 predictions over empirical distribution (is it ok that this slightly extends beyond x = 1??)
ppc_hist(y, yrep[1:55, ])                                  # compare 55 predictions to empirical distribution
ppc_pit_ecdf(y, yrep[1:50,])                               # no idea....
ppc_ecdf_overlay(y, yrep[1:50,]) + xaxis_text()            # also no idea...

## save output
save.image(paste0('anp_edgecalculations/anpshort1_nodalregression.RData'))
dev.off()


## run model -- cmdstanr ####
## create data list
#load('anp_nodalregression/anpshort1_nodalregression_conditionaledge.RData')
eigen_all <- eigen[,,1]
for(i in 2:(n_samples*n_chains)){
  eigen_all <- rbind(eigen_all, eigen[,,i])
}
eigen_all <- as.data.frame(eigen_all) %>% 
  left_join(nodes, by = c('node','age','sightings'))
eigen_wide <- matrix(NA, nrow = n_samples*n_chains, ncol = n_eles,
                     dimnames = list(1:(n_samples*n_chains),
                                     ele_ids))
for(i in 1:n_eles){
  node_eigen <- eigen_all %>%
    filter(id == colnames(eigen_wide)[i])
  eigen_wide[,i] <- node_eigen$eigenvector
}
rm(node_eigen) ; gc()
mean_eigen <- apply(eigen_wide, 2, mean) %>% 
  as.data.frame()
mean_eigen$id <- rownames(mean_eigen)
colnames(mean_eigen)[1] <- 'mean_eigen'
nodes <- nodes %>% left_join(mean_eigen,  by = 'id')

stdv_eigen <- apply(eigen_wide, 2, sd) %>% 
  as.data.frame()
stdv_eigen$id <- rownames(stdv_eigen)
colnames(stdv_eigen)[1] <- 'stdv_eigen'
nodes <- nodes %>% left_join(stdv_eigen,  by = 'id')

eigen_list <- list(num_nodes = nrow(nodes),
                   nodes = nodes$node,
                   centrality_mu = nodes$mean_eigen,
                   centrality_sd = nodes$stdv_eigen,
                   node_age = nodes$age,
                   node_age2 = (nodes$age)^2)

# set prior
nodal_regression <- cmdstan_model('models/eigen_regression.stan')

## run model
fit_anp1_eigen <- nodal_regression$sample(
  data = eigen_list,
  chains = n_chains,
  parallel_chains = n_chains
)

## save output
save.image('anp_nodalregression/anpshort1_nodalregression_conditionaledge_cmdstan.RData')

## posterior check ####
# extract draws
post_eigen <- as.data.frame(as_draws_df(fit_anp1_eigen)) %>% 
  janitor::clean_names() %>% 
  select(-lp) %>% 
  pivot_longer(cols = 1:last_col(3))

# traceplot linear effect size
b_age <- post_eigen %>% filter(name == 'beta_age')
plot(data = b_age[b_age$chain == 1,], value ~ iteration, col = 'purple',
     type = 'l', xlim = c(0,1000))
lines(data = b_age[b_age$chain == 2,], value ~ iteration, col = 'red')
lines(data = b_age[b_age$chain == 3,], value ~ iteration, col = 'blue')
lines(data = b_age[b_age$chain == 4,], value ~ iteration, col = 'green')

# traceplot quadratic effect size
b_age2 <- post_eigen %>% filter(name == 'beta_age2')
plot(data = b_age2[b_age2$chain == 1,], value ~ iteration, col = 'purple',
     type = 'l', xlim = c(0,1000))
lines(data = b_age2[b_age2$chain == 2,], value ~ iteration, col = 'red')
lines(data = b_age2[b_age2$chain == 3,], value ~ iteration, col = 'blue')
lines(data = b_age2[b_age2$chain == 4,], value ~ iteration, col = 'green')

hist(b_age$value)   # natural scale? should it be: hist(plogis(b_age$value))
hist(b_age2$value)  # hist(plogis(b_age2$value))

# plot raw data
ggplot()+
  geom_point(data = eigen_all, mapping = aes(x = age, y = eigenvector),
             colour = rgb(253/255, 231/255, 37/255, 0.01))+
  geom_point(data = nodes,
             mapping = aes(x = age, y = mean_eigen, size = sightings),
             colour = rgb(68/255, 1/255, 84/255))+
  scale_x_continuous('age (years)')+
  scale_y_continuous('eigenvector centrality')+
  theme_bw()+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

# compare empirical distribution to posterior predictive distribution
nodes$node_rank <- as.integer(as.factor(nodes$node))
mean_predict <- post_eigen %>% 
  filter(name != 'beta_age') %>% filter(name != 'beta_age2') %>% 
  filter(name != 'sigma') %>% 
  group_by(name) %>% 
  mutate(mean_predict = mean(value)) %>% 
  ungroup() %>% 
  separate(name, into = c('predictor','node_rank'), remove = F, sep = '_') %>% 
  select(-predictor) %>% 
  mutate(node_rank = as.integer(node_rank)) %>% 
  rename(predicted_eigen = value) %>% 
  left_join(nodes, by = 'node_rank')
mean_predict %>% 
  select(mean_eigen, mean_predict, age, sightings) %>% 
  distinct() %>% 
  ggplot()+
  geom_point(aes(x = mean_eigen, y = mean_predict, colour = age, size = sightings))+
  scale_x_continuous('mean calculated eigenvector')+
  scale_y_continuous('mean predicted eigenvector')+
  theme_classic()
ggplot(mean_predict)+
  geom_point(aes(x = mean_eigen, y = predicted_eigen), colour = rgb(0,0,1,0.01))+
  geom_point(aes(x = mean_eigen, y = mean_predict), colour = 'white')+
  scale_x_continuous('mean calculated eigenvector')+
  scale_y_continuous('mean predicted eigenvector')+
  theme_classic()

## save output
save.image('anp_nodalregression/anpshort1_nodalregression_conditionaledge_cmdstan.RData')
dev.off()

#### check outputs and make graphs ####
load('anp_nodalregression/anpshort1_nodalregression_conditionaledge_cmdstan.RData')




















#### run model for all other windows ####