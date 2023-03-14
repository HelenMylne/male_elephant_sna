#### Bayesian analysis of EfA data ####
# Script to process association data from Makgadikgadi Pans National Park, Botswana
# Data collected: 2012-2021 by Elephants for Africa
# Data supplied by: Kate Evans

#### Set up ####
# load packages
library(tidyverse, lib.loc = '../packages/')   # library(tidyverse)
library(dplyr, lib.loc = '../packages/')       # library(dplyr)
#library(rstan, lib.loc = '../packages/')      # library(rstan)
library(cmdstanr, lib.loc = '../packages/')    # library(cmdstanr)
library(bisonR, lib.loc = '../packages/')      # library(bisonR)
#library(igraph, lib.loc = '../packages/')     # library(igraph)
library(asnipe, lib.loc = '../packages/')      # library(asnipe)
library(sna, lib.loc = '../packages/')         # library(sna)
library(raster, lib.loc = '../packages/')      # library(raster)

# information
sessionInfo()
R.Version()

# set seed
set.seed(12345)

#### create priors and functions for all ####
# set priors
priors <- get_default_priors('binary')
priors$edge <- 'normal(-2.5, 1.5)'    # slightly right skewed so that more density is towards 0 but still very weak and all values are possible
priors$fixed <- 'normal(0, 1.5)'  # slightly less wide that before, but zero centred so allows age to have a positive or negative impact on centrality and edge weight
prior_check(priors, 'binary')

### adapt bisonR plot_network function to give more flexibility over plotting options
plot_network_threshold <- function (obj, ci = 0.9, lwd = 2, threshold = 0.3,
                                    vertex.label.color1 = 'transparent', vertex.label.font1 = NA, 
                                    vertex.color1 = 'transparent',
                                    vertex.size1 = 1, edge.color1 = rgb(0, 0, 0, 1),
                                    vertex.label.color2 = 'black', vertex.label.font2 = 'Helvetica', 
                                    vertex.color2 = 'seagreen1',
                                    vertex.size2 = 8, edge.color2 = rgb(0, 0, 0, 0.3))
{
  edgelist <- get_edgelist(obj, ci = ci, transform = TRUE)
  net <- igraph::graph_from_edgelist(as.matrix(edgelist[, 1:2]), 
                                     directed = obj$directed)
  md <- edgelist[, 3]
  ub <- edgelist[, 5]
  coords <- igraph::layout_nicely(net)
  igraph::plot.igraph(net, layout = coords,
                      edge.width = ifelse(md < threshold, 0, md * lwd),
                      vertex.label.color = vertex.label.color1, vertex.color = vertex.color1, vertex.size = vertex.size1,
                      edge.color = rgb(0, 0, 0, 1), edge.arrow.size = 0)
  igraph::plot.igraph(net, layout = coords,
                      edge.width = ifelse(md < threshold, 0, ub * lwd),
                      vertex.label.color = vertex.label.color2, vertex.color = vertex.color2, vertex.size = vertex.size2,
                      edge.color = edge.color1, edge.arrow.size = 0, 
                      add = TRUE)
  if (obj$directed) {
    igraph::plot.igraph(net, edge.width = md * 0, layout = coords, 
                        vertex.label.color = vertex.label.color2, vertex.color = vertex.color2, 
                        edge.color = edge.color2, add = TRUE)
  }
}

### adapt to remove unconnected nodes
plot_network_threshold2 <- function (obj, ci = 0.9, lwd = 2, threshold = 0.3,
                                     label.colour = 'transparent', label.font = 'Helvetica', 
                                     node.size = 4, node.colour = 'seagreen1',
                                     link.colour1 = 'black', link.colour2 = rgb(0, 0, 0, 0.3))
{
  edgelist <- get_edgelist(obj, ci = ci, transform = TRUE)
  threshold_edges <- edgelist[edgelist$median >= threshold,]
  net <- igraph::graph_from_edgelist(as.matrix(threshold_edges[, 1:2]))
  
  if(is.data.frame(node.size) == TRUE ) {
    nodes_list <- data.frame(node = as.numeric(names(net[1])),
                             sightings = NA)
    for(i in 1:nrow(nodes_list)){
      nodes_list$sightings[i] <- nodes$sightings[which(nodes$node == nodes_list$node[i])]
    }
    node_sightings <- log(nodes_list$sightings)*5
  } else { node_sightings <- node.size }
  
  if(is.data.frame(node.colour) == TRUE ) {
    nodes_list <- data.frame(node = as.numeric(names(net[1])),
                             age = NA)
    for(i in 1:nrow(nodes_list)){
      nodes_list$age[i] <- nodes$age[which(nodes$node == nodes_list$node[i])]
    }
    node_age <- nodes_list$age
  } else { node_age <- node.colour }
  
  md <- threshold_edges[, 3]
  ub <- threshold_edges[, 5]
  coords <- igraph::layout_nicely(net)
  igraph::plot.igraph(net, layout = coords,
                      vertex.label.color = ifelse(is.null(label.colour) == TRUE,
                                                  ifelse(node_age < 20, 'black', 'white'),
                                                  label.colour),
                      label.family = label.font,
                      vertex.color = ifelse(node_age < 15, '#FDE725FF',
                                            ifelse(node_age < 20, '#55C667FF',
                                                   ifelse(node_age < 30, '#1F968BFF', 
                                                          ifelse(node_age < 40, '#39568CFF', '#440154FF')))), 
                      vertex.size = node_sightings,
                      frame.color = NA, frame.width = 0,
                      edge.color = NA, edge.arrow.size = 0, edge.width = 0)
  igraph::plot.igraph(net, layout = coords, add = TRUE,
                      vertex.label = NA, vertex.color = 'transparent', vertex.size = 0, 
                      frame.color = NA, frame.width = 0,
                      edge.color = link.colour1, edge.arrow.size = 0, edge.width = md * lwd)
  igraph::plot.igraph(net, layout = coords, add = TRUE,
                      vertex.label = NA, vertex.color = 'transparent', vertex.size = 0, 
                      frame.color = NA, frame.width = 0,
                      edge.color = link.colour2, edge.arrow.size = 0, edge.width = ub * lwd)
}

# alternative method to compare models
model_averaging <- function(models) {
  loos <- lapply(models, function(model) loo::loo(model$log_lik, r_eff=NA))
  
  results_matrix <- loo::loo_model_weights(loos, method="pseudobma")
  if (!is.null(names(models))) {
    names(results_matrix) <- names(models)
  }
  results_matrix
}

dev.off()

#### short time windows ####
for( time_window in 1:5 ){
  ## set up ####
  ### add time marker
  print(paste0('start window ', time_window ,' at ', Sys.time()))
  
  ### create output pdf
  pdf(file = paste0('../outputs/mpnpshort',time_window,'_bisonr_edgeweight.pdf'))
  
  ### read in data frame for edge weight model
  filename <- paste0('../data_processed/mpnp_period',time_window,'_pairwiseevents.csv')
  counts_df <- read_csv(filename)
  
  ### read in ages for filtering
  filename <- paste0('../data_processed/mpnp',time_window,'_ageestimates_mcmcoutput.rds')
  mpnp_ages <- readRDS(filename) %>%
    pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')
  
  ### filter counts data frame down to males of known age only
  length(unique(counts_df$id_1))
  length(unique(mpnp_ages$id))   # missing individuals as some had no usable age estimate
  counts_df <- counts_df %>% 
    filter(id_1 %in% unique(mpnp_ages$id)) %>% 
    filter(id_2 %in% unique(mpnp_ages$id))
  
  ### create model data
  counts_df_model <- counts_df[, c('node_1','node_2','together','count_dyad')] %>%
    distinct()
  colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')
  
  ## run model ####
  ### run edge weight model
  mpnp_edge_weights <- bison_model(
    ( event | duration ) ~ dyad(node_1_id, node_2_id), 
    data = counts_df_model, 
    model_type = "binary",
    priors = priors
  )

  ### run diagnostic plots
  plot_trace(mpnp_edge_weights, par_ids = 2)
  plot_predictions(mpnp_edge_weights, num_draws = 20, type = "density")
  plot_predictions(mpnp_edge_weights, num_draws = 20, type = "point")
  
  ### extract edge weight summaries
  edgelist <- get_edgelist(mpnp_edge_weights, ci = 0.9, transform = TRUE)
  plot(density(edgelist$median))
  summary(mpnp_edge_weights)
  
  ### add time marker
  print(paste0('model completed at ', Sys.time()))
  
  ## non-random edge weights ####
  ### run null model
  mpnp_edges_null <- bison_model(
    (event | duration) ~ 1, 
    data = counts_df_model, 
    model_type = "binary",
    priors = priors
  )
  
  ### compare null model with fitted model -- bisonR model stacking
  model_comparison(list(non_random_model = mpnp_edge_weights, random_model = mpnp_edges_null))
  
  ### compare null model with fitted model -- Jordan pseudo model averaging
  model_averaging(models = list(non_random_model = mpnp_edge_weights, random_model = mpnp_edges_null))
  
  ### add time marker
  print(paste0('random network comparison completed at ', Sys.time()))
  
  # save workspace image
  save.image(file = paste0('mpnp_edgecalculations/mpnpshort',time_window,'_bisonr_edgescalculated.RData'))
  
  ## plot network ####
  # create nodes data frame
  #nodes <- data.frame(bull = sort(unique(mpnp_ages$id)),
  #                    age = NA,
  #                    sightings = NA)
  #for(i in 1:nrow(nodes)){
  #  x <- mpnp_ages[mpnp_ages$id == nodes$bull[i],]
  #  nodes$age[i] <- mean(x$age)
  #  if(nodes$bull[i] != 'M99'){
  #    y <- counts_df[counts_df$id_1 == nodes$bull[i], c('id_1','count_period_1')]
  #    nodes$sightings[i] <- y[1,2]
  #  }
  #}
  #nodes$sightings <- as.numeric(nodes$sightings)
  #str(nodes)
  
  # plot network
  #plot_network_threshold(mpnp_edge_weights, lwd = 2, ci = 0.9, threshold = 0.15,
  #                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
  #                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
  #                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
  #plot_network_threshold2(obj = mpnp_edge_weights, threshold = 0.15,
  #                        node.colour = nodes$age, node.size = nodes$sightings)

  ### add time marker
  print(paste0('network plots completed at ', Sys.time()))
  
  ## compare edge weight distributions to simple SRI ####
  ### plot against SRI
  colnames(edgelist)[1:2] <- colnames(counts_df_model)[1:2]
  edgelist$node_1_id <- as.integer(edgelist$node_1_id) ; edgelist$node_2_id <- as.integer(edgelist$node_2_id)
  summary <- left_join(edgelist, counts_df_model, by = c('node_1_id','node_2_id'))
  counts <- counts_df[,c('node_1','node_2','count_period_1','count_period_2')]
  colnames(counts)[1:2] <- c('node_1_id','node_2_id')
  summary <- left_join(summary, counts, by = c('node_1_id','node_2_id'))
  summary$sri <- summary$event / (summary$duration)
  
  #plot(density(summary$sri), main = 'SRI vs model output: blue=all,\nred=both seen 8 times, green=both 12 times')
  #lines(density(summary$median), col = 'blue')
  #lines(density(summary$median[which(summary$count_1 >= 8 & summary$count_2 >= 8)]), col = 'red')
  #lines(density(summary$median[which(summary$count_1 >= 12 & summary$count_2 >= 12)]), col = 'green')
  
  # try plotting a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
  dyads <- counts_df[,c('dyad_id','node_1','node_2')]
  colnames(dyads) <- c('dyad_id','node_1_id','node_2_id')
  dyads <- left_join(dyads, counts_df_model, by = c('node_1_id','node_2_id'))
  length(which(is.na(dyads$duration) == TRUE))
  
  draws <- as.data.frame(mpnp_edge_weights$chain) %>%
    pivot_longer(cols = everything(), names_to = 'dyad_model', values_to = 'edge')
  draws$dyad_id <- rep(counts_df$dyad_id, 4000)
  draws$weight <- gtools::inv.logit(draws$edge)
  draws$draw <- rep(1:4000,  each = nrow(counts_df_model))
  draws <- left_join(draws, dyads, by = 'dyad_id')
  draws$sri <- draws$event / draws$duration
  
  #subset_draws <- draws[draws$sri > 0.2,]
  #subset_draws$median <- NA
  #for(i in 1:nrow(subset_draws)){
  #  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]
  #  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), subset_draws$median[i])
  #}
  #head(subset_draws)
  #which(is.na(subset_draws$median) == TRUE)[1]
  
  #subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)
  #ggplot(data = subset_draws, mapping = aes(x = weight))+
  #  geom_density(colour = 'blue')+
  #  facet_wrap(. ~ dyad_id, ncol = 15)+
  #  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  #  geom_vline(mapping = aes(xintercept = sri), colour = 'red')
  
  # write_csv(subset_draws, '../data_processed/mpnpshort',time_window,'_sampledyads_sri0.2_binary_vs_sri.csv')
  
  # clean environment
  #rm(draws, dyads, priors, subset_draws, x) ; gc()
  
  ## coefficient of variation of edge weights (aka social differentiation) ####
  # extract cv for model
  global_cv <- extract_metric(mpnp_edge_weights, "global_cv")
  head(global_cv)
  hist(global_cv)
  
  ### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
  counts_df_model <- counts_df[,c('node_1','node_2','together','count_dyad','id_1','id_2','dyad_id')]
  colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
  ew_chain <- as.data.frame(mpnp_edge_weights$chain)
  colnames(ew_chain) <- counts_df_model$dyad_id
  ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
  ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
  ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
  ew_chain$mean <- NA
  hist(ew_chain$draw)
  #plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'edge distribution', xlab = 'edge weight', ylab = 'density', las = 1)
  #for(i in sort(unique(ew_chain$dyad_id))){
  #  x <- ew_chain[ew_chain$dyad_id == i,]
  #  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  #  lines(density(x$draw), col = rgb(0,0,1,0.1))
  #}
  #lines(density(ew_chain$draw), lwd = 2)
  
  (draw98 <- quantile(ew_chain$draw, 0.98))
  
  ew_edgelist <- bisonR::get_edgelist(mpnp_edge_weights)
  (median98 <- quantile(ew_edgelist$median, 0.98))
  
  counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
  (sri98 <- quantile(counts_df_model$sri, 0.98))
  
  ## clean up ####
  ### save pdf
  dev.off()
  
  # save workspace image
  save.image(file = paste0('mpnp_edgecalculations/mpnpshort',time_window,'_bisonr_edgescalculated.RData'))
  
  ### clear environment
  rm(list = ls()[!(ls() %in% c('time_window','priors',
                               'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])
  
  ### add time marker
  print(paste0('time window ', time_window, ' completed at ', Sys.time()))
  
}

#### long time window ####
## set up ####
### add time marker
print(paste0('start long window at ', Sys.time()))

### create output pdf
pdf('../outputs/mpnplongwindow_bisonr_edgeweight.pdf')

### read in data frame for edge weight model
counts_df <- read_csv('../data_processed/mpnp_longtimewindow_pairwiseevents.csv')

### read in ages for filtering
ages <- readRDS('../data_processed/mpnplongwindow_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

### filter counts data frame down to males of known age only
length(unique(counts_df$id_1))
length(unique(mpnp_ages$id))   # missing individuals as some had no usable age estimate
counts_df <- counts_df %>% 
  filter(id_1 %in% unique(mpnp_ages$id)) %>% 
  filter(id_2 %in% unique(mpnp_ages$id))

### create model data
counts_df_model <- counts_df[, c('node_1','node_2','together','count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

## run model ####
### run edge weight model
mpnp_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### run diagnostic plots
plot_trace(mpnp_edge_weights, par_ids = 2)
plot_predictions(mpnp_edge_weights, num_draws = 20, type = "density")
plot_predictions(mpnp_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(mpnp_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(mpnp_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
mpnp_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = mpnp_edge_weights, random_model = mpnp_edges_null))

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = mpnp_edge_weights, random_model = mpnp_edges_null))

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

# save workspace image
save.image(file = paste0('mpnp_edgecalculations/mpnp_longwindow_bisonr_edgescalculated.RData'))

## plot network ####
# create nodes data frame
#nodes <- data.frame(bull = sort(unique(mpnp_ages$id)),
#                    age = NA,
#                    sightings = NA)
#for(i in 1:nrow(nodes)){
#  x <- mpnp_ages[mpnp_ages$id == nodes$bull[i],]
#  nodes$age[i] <- mean(x$age)
#  if(nodes$bull[i] != 'M99'){
#    y <- counts_df[counts_df$id_1 == nodes$bull[i], c('id_1','count_period_1')]
#    nodes$sightings[i] <- y[1,2]
#  }
#}
#nodes$sightings <- as.numeric(nodes$sightings)
#str(nodes)

# plot network
#plot_network_threshold(mpnp_edge_weights, lwd = 2, ci = 0.9, threshold = 0.15,
#                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
#                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
#                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
#plot_network_threshold2(obj = mpnp_edge_weights, threshold = 0.15,
#                        node.size = nodes, node.colour = nodes, lwd = 10)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## compare edge weight distributions to simple SRI ####
### plot against SRI
colnames(edgelist)[1:2] <- colnames(counts_df_model)[1:2]
edgelist$node_1_id <- as.integer(edgelist$node_1_id) ; edgelist$node_2_id <- as.integer(edgelist$node_2_id)
summary <- left_join(edgelist, counts_df_model, by = c('node_1_id','node_2_id'))
counts <- counts_df[,c('node_1','node_2','count_period_1','count_period_2')]
colnames(counts)[1:2] <- c('node_1_id','node_2_id')
summary <- left_join(summary, counts, by = c('node_1_id','node_2_id'))
summary$sri <- summary$event / (summary$duration)

#plot(density(summary$sri), main = 'SRI vs model output: blue=all,\nred=both seen 8 times, green=both 12 times')
#lines(density(summary$median), col = 'blue')
#lines(density(summary$median[which(summary$count_1 >= 8 & summary$count_2 >= 8)]), col = 'red')
#lines(density(summary$median[which(summary$count_1 >= 12 & summary$count_2 >= 12)]), col = 'green')

# try plotting a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
dyads <- counts_df[,c('dyad_id','node_1','node_2')]
colnames(dyads) <- c('dyad_id','node_1_id','node_2_id')
dyads <- left_join(dyads, counts_df_model, by = c('node_1_id','node_2_id'))
length(which(is.na(dyads$duration) == TRUE))

draws <- as.data.frame(mpnp_edge_weights$chain) %>%
  pivot_longer(cols = everything(), names_to = 'dyad_model', values_to = 'edge')
draws$dyad_id <- rep(counts_df$dyad_id, 4000)
draws$weight <- gtools::inv.logit(draws$edge)
draws$draw <- rep(1:4000,  each = nrow(counts_df_model))
draws <- left_join(draws, dyads, by = 'dyad_id')
draws$sri <- draws$event / draws$duration

#subset_draws <- draws[draws$sri > 0.2,]
#subset_draws$median <- NA
#for(i in 1:nrow(subset_draws)){
#  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]
#  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), subset_draws$median[i])
#}
#head(subset_draws)
#which(is.na(subset_draws$median) == TRUE)[1]

#subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)
#ggplot(data = subset_draws, mapping = aes(x = weight))+
#  geom_density(colour = 'blue')+
#  facet_wrap(. ~ dyad_id, ncol = 15)+
#  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
#  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

#write_csv(subset_draws, '../data_processed/mpnp_longwindow_sampledyads_sri0.2_binary_vs_sri.csv')

# clean environment
#rm(draws, dyads, priors, subset_draws, x) ; gc()

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(mpnp_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[,c('node_1','node_2','together','count_dyad','id_1','id_2','dyad_id')]
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(mpnp_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
#plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'edge distribution', xlab = 'edge weight', ylab = 'density', las = 1)
#for(i in sort(unique(ew_chain$dyad_id))){
#  x <- ew_chain[ew_chain$dyad_id == i,]
#  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
#  lines(density(x$draw), col = rgb(0,0,1,0.1))
#}
#lines(density(ew_chain$draw), lwd = 2)

(draw98 <- quantile(ew_chain$draw, 0.98))

ew_edgelist <- bisonR::get_edgelist(mpnp_edge_weights)
(median98 <- quantile(ew_edgelist$median, 0.98))

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
(sri98 <- quantile(counts_df_model$sri, 0.98))

## clean up ####
### save pdf
dev.off()

# save workspace image
save.image(file = paste0('mpnp_edgecalculations/mpnp_longwindow_bisonr_edgescalculated.RData'))

### clear environment
rm(list = ls()[!(ls() %in% c('time_window','priors',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

### add time marker
print(paste0('time window ', time_window, ' completed at ', Sys.time()))
