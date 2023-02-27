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

#### Period 1 ####
# create file of output graphs
pdf('../outputs/mpnp1_networkplots.pdf', width = 10, height = 10)

## create data list ####
counts_df1 <- read_csv('../data_processed/mpnp_period1_pairwiseevents.csv')
str(counts_df1)

### load in ages 
mpnp1_ages <- readRDS('../data_processed/mpnp1_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

### filter counts data frame down to males only
counts_df1 <- counts_df1 %>% 
  filter(id_1 %in% unique(mpnp1_ages$id)) %>% 
  filter(id_2 %in% unique(mpnp1_ages$id))

# reassign dyad numbers to remove gaps
counts_df1$node_1_males <- as.integer(as.factor(counts_df1$node_1))
counts_df1$node_2_males <- as.integer(as.factor(counts_df1$node_2))+1

### standardise dyad_id
counts_df1$dyad_males <- as.integer(as.factor(counts_df1$dyad_id))

## edge weights ####
### create data frame for edge weight model
counts_df1_model <- counts_df1[, c('node_1_males','node_2_males','together','count_dyad')] %>% distinct()
colnames(counts_df1_model) <- c('node_1_id','node_2_id','event','duration')

### set priors
priors <- get_default_priors('binary')
priors$edge <- 'normal(-2.5, 1.5)'    # slightly right skewed so that more density is towards 0 but still very weak and all values are possible
priors$fixed <- 'normal(0, 1.5)'  # slightly less wide that before, but zero centred so allows age to have a positive or negative impact on centrality and edge weight
prior_check(priors, 'binary')

### run edge weight model
mpnp1_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df1_model, 
  model_type = "binary",
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(mpnp1_edge_weights, par_ids = 2)
plot_predictions(mpnp1_edge_weights, num_draws = 20, type = "density")
plot_predictions(mpnp1_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(mpnp1_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(mpnp1_edge_weights)

### save workspace image for reloading at a later date that doesn't require running model again
save.image('mpnp_bisonr_edgescalculated.RData')

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### load in workspace image with edge weights already calculated
#load('mpnp_bisonr_edgescalculated.RData')

### run null model
mpnp1_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df1_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR comparison
model_comparison(list(non_random_model = mpnp1_edge_weights, random_model = mpnp1_edges_null))

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging <- function(models) {
  loos <- lapply(models, function(model) loo::loo(model$log_lik, r_eff=NA))
  
  results_matrix <- loo::loo_model_weights(loos, method="pseudobma")
  if (!is.null(names(models))) {
    names(results_matrix) <- names(models)
  }
  results_matrix
}
model_averaging(models = list(non_random_model = mpnp1_edge_weights, random_model = mpnp1_edges_null))

# save workspace image
save.image('mpnp_bisonr_randomnetwork.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
### adapt bisonr plot_network function to give more flexibility over plotting options
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

### load in ages 
#mpnp1_ages <- readRDS('../data_processed/mpnp1_ageestimates_mcmcoutput.rds') %>%
#  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

### create nodes data frame
nodes <- data.frame(bull = sort(unique(mpnp_ages$id)),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  x <- mpnp_ages[mpnp_ages$id == nodes$bull[i],]
  nodes$age[i] <- mean(x$age)
  if(nodes$bull[i] != 'M99'){
    y <- counts_df1[counts_df1$id_1 == nodes$bull[i], c('id_1','count_1')]
    nodes$sightings[i] <- y[1,2]
  }
}
nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

### plot network
plot_network_threshold(mpnp1_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')

### adapt to remove unconnected nodes
plot_network_threshold2 <- function (obj, ci = 0.9, lwd = 2, threshold = 0.3,
                                     vertex.label.color1 = 'transparent', vertex.label.font1 = NA, 
                                     vertex.color1 = 'transparent',
                                     vertex.size1 = 1, edge.color1 = rgb(0, 0, 0, 1),
                                     vertex.label.color2 = 'black', vertex.label.font2 = 'Helvetica', 
                                     vertex.color2 = 'seagreen1',
                                     vertex.size2 = 8, edge.color2 = rgb(0, 0, 0, 0.3))
{
  edgelist <- get_edgelist(obj, ci = ci, transform = TRUE)
  threshold_edges <- edgelist[edgelist$median >= threshold,]
  net <- igraph::graph_from_edgelist(as.matrix(threshold_edges[, 1:2]))
  md <- threshold_edges[, 3]
  ub <- threshold_edges[, 5]
  coords <- igraph::layout_nicely(net)
  igraph::plot.igraph(net, layout = coords,
                      vertex.label.color = vertex.label.color1, label.family = vertex.label.font1, 
                      vertex.color = vertex.color1, vertex.size = vertex.size1,
                      edge.color = rgb(0, 0, 0, 1), edge.arrow.size = 0, edge.width = md * lwd)
  igraph::plot.igraph(net, layout = coords,
                      vertex.label.color = vertex.label.color2, label.family = vertex.label.font2,
                      vertex.color = vertex.color2, vertex.size = vertex.size2,
                      edge.color = edge.color1, edge.arrow.size = 0, edge.width = ub * lwd,
                      add = TRUE)
  if (obj$directed) {
    igraph::plot.igraph(net, edge.width = md * 0, layout = coords, 
                        vertex.label.color = vertex.label.color2, vertex.color = vertex.color2, 
                        edge.color = edge.color2, add = TRUE)
  }
}

### plot network
plot_network_threshold2(obj = mpnp1_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## compare edge weight distributions to simple SRI ####
# load in workspace image with edge weights already calculated
#load('mpnp_bisonr_edgescalculated.RData')

### plot against SRI
head(edgelist)
colnames(edgelist)[1:2] <- colnames(counts_df_model)[1:2]
edgelist$node_1_id <- as.integer(edgelist$node_1_id) ; edgelist$node_2_id <- as.integer(edgelist$node_2_id)
summary <- left_join(edgelist, counts_df_model, by = c('node_1_id','node_2_id'))
counts <- counts_df[,c('node_1_males','node_2_males','count_1','count_2')]
colnames(counts)[1:2] <- c('node_1_id','node_2_id')
summary <- left_join(summary, counts, by = c('node_1_id','node_2_id'))
summary$sri <- summary$event / (summary$duration)

plot(density(summary$sri), main = 'SRI vs model output: blue=all,\nred=both seen 8 times, green=both 12 times')
lines(density(summary$median), col = 'blue')
lines(density(summary$median[which(summary$count_1 >= 8 & summary$count_2 >= 8)]), col = 'red')
lines(density(summary$median[which(summary$count_1 >= 12 & summary$count_2 >= 12)]), col = 'green')

# try plotting a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
dyads <- counts_df[,c('dyad_id','node_1_males','node_2_males')]
colnames(dyads) <- c('dyad_id','node_1_id','node_2_id')
dyads <- left_join(dyads, counts_df_model, by = c('node_1_id','node_2_id'))
length(which(is.na(dyads$duration) == TRUE))

draws <- as.data.frame(mpnp1_edge_weights$chain) %>% pivot_longer(cols = everything(), names_to = 'dyad_model', values_to = 'edge')
draws$dyad_id <- rep(counts_df$dyad_id, 4000) ## IS THIS RIGHT??
draws$weight <- gtools::inv.logit(draws$edge)
draws$draw <- rep(1:4000,  each = nrow(counts_df_model))
draws <- left_join(draws, dyads, by = 'dyad_id')
draws$sri <- draws$event / draws$duration

#set.seed(15)
#subset_draws <- draws[draws$dyad_id %in% sample(draws$dyad_id, 150, replace = F),]
#subset_draws$median <- NA
#for(i in 1:nrow(subset_draws)){
#  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]
#  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), #subset_draws$median[i])
#}
#head(subset_draws)
#which(is.na(subset_draws$median) == TRUE)[1]

#subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)
#ggplot(data = subset_draws, mapping = aes(x = weight))+
#  geom_density(colour = 'blue')+
#  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+
#  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
#  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

#write_csv(subset_draws, '../data_processed/mpnp1_sampledyads_random_binary_vs_sri.csv')

subset_draws <- draws[draws$sri > 0.2,]
subset_draws$median <- NA
for(i in 1:nrow(subset_draws)){
  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]
  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), subset_draws$median[i])
}
head(subset_draws)
which(is.na(subset_draws$median) == TRUE)[1]

subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)
ggplot(data = subset_draws, mapping = aes(x = weight))+
  geom_density(colour = 'blue')+
  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15)+
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

ggplot(data = subset_draws, mapping = aes(x = weight))+
  geom_density(colour = 'blue')+
  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

write_csv(subset_draws, '../data_processed/mpnp1_sampledyads_sri0.2_binary_vs_sri.csv')

# clean environment
rm(draws, dyads, priors, subset_draws, x) ; gc()
dev.off()

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(mpnp1_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[,c('node_1_males','node_2_males','event_count','count_dyad','id_1','id_2','dyad_males')]
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ms_chain <- as.data.frame(mpnp1_edge_weights$chain)
colnames(ms_chain) <- counts_df_model$dyad_id
ms_chain <- pivot_longer(ms_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ms_chain$chain_position <- rep(1:length(unique(ms_chain$dyad_id)), each = 4000)
ms_chain$draw <- LaplacesDemon::invlogit(ms_chain$draw)
ms_chain$mean <- NA
hist(ms_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ms_chain$dyad_id))){
  x <- ms_chain[ms_chain$dyad_id == i,]
  ms_chain$mean <- ifelse(ms_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ms_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ms_chain$draw), lwd = 2)

quantile(ms_chain$draw, 0.98)

ms_edgelist <- bisonR::get_edgelist(mpnp1_edge_weights)
quantile(ms_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### save pdf
dev.off()

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
rm(list= ls()[!(ls() %in% c('mpnp1_edge_weights','mpnp1_edges_null',
                            'counts_df','counts_df_model',
                            'gbi_males','m_mat',
                            #'random_nodes',
                            'random_networks'))])

### write over saved workspace image -- only keep the things you might need later that take a long time to run
save.image('mpnp1_bisonr_edgescalculated.RData')

### end pdf
dev.off()

#### Period 2 ####
# create file of output graphs
pdf('../outputs/mpnp2_networkplots.pdf', width = 10, height = 10)

## create data list ####
counts_df2 <- read_csv('../data_processed/mpnp_period2_pairwiseevents.csv')
str(counts_df2)

### load in ages 
mpnp2_ages <- readRDS('../data_processed/mpnp2_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

### filter counts data frame down to males only
counts_df2 <- counts_df2 %>% 
  filter(id_1 %in% unique(mpnp2_ages$id)) %>% 
  filter(id_2 %in% unique(mpnp2_ages$id))

# reassign dyad numbers to remove gaps
counts_df2$node_1_males <- as.integer(as.factor(counts_df2$node_1))
counts_df2$node_2_males <- as.integer(as.factor(counts_df2$node_2))+1

### standardise dyad_id
counts_df2$dyad_males <- as.integer(as.factor(counts_df2$dyad_id))

## edge weights ####
### create data frame for edge weight model
counts_df2_model <- counts_df2[, c('node_1_males','node_2_males','together','count_dyad')] %>% distinct()
colnames(counts_df2_model) <- c('node_1_id','node_2_id','event','duration')

### set priors
priors <- get_default_priors('binary')
priors$edge <- 'normal(-2.5, 1.5)'    # slightly right skewed so that more density is towards 0 but still very weak and all values are possible
priors$fixed <- 'normal(0, 1.5)'  # slightly less wide that before, but zero centred so allows age to have a positive or negative impact on centrality and edge weight
prior_check(priors, 'binary')

### run edge weight model
mpnp2_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df2_model, 
  model_type = "binary",
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(mpnp2_edge_weights, par_ids = 2)
plot_predictions(mpnp2_edge_weights, num_draws = 20, type = "density")
plot_predictions(mpnp2_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(mpnp2_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(mpnp2_edge_weights)

### save workspace image for reloading at a later date that doesn't require running model again
save.image('mpnp_bisonr_edgescalculated.RData')

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### load in workspace image with edge weights already calculated
#load('mpnp_bisonr_edgescalculated.RData')

### run null model
mpnp2_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df2_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR comparison
model_comparison(list(non_random_model = mpnp2_edge_weights, random_model = mpnp2_edges_null))

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging <- function(models) {
  loos <- lapply(models, function(model) loo::loo(model$log_lik, r_eff=NA))
  
  results_matrix <- loo::loo_model_weights(loos, method="pseudobma")
  if (!is.null(names(models))) {
    names(results_matrix) <- names(models)
  }
  results_matrix
}
model_averaging(models = list(non_random_model = mpnp2_edge_weights, random_model = mpnp2_edges_null))

# save workspace image
save.image('mpnp_bisonr_randomnetwork.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
### adapt bisonr plot_network function to give more flexibility over plotting options
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

### load in ages 
#mpnp2_ages <- readRDS('../data_processed/mpnp2_ageestimates_mcmcoutput.rds') %>%
#  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

### create nodes data frame
nodes <- data.frame(bull = sort(unique(mpnp_ages$id)),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  x <- mpnp_ages[mpnp_ages$id == nodes$bull[i],]
  nodes$age[i] <- mean(x$age)
  if(nodes$bull[i] != 'M99'){
    y <- counts_df2[counts_df2$id_1 == nodes$bull[i], c('id_1','count_1')]
    nodes$sightings[i] <- y[1,2]
  }
}
nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

### plot network
plot_network_threshold(mpnp2_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')

### adapt to remove unconnected nodes
plot_network_threshold2 <- function (obj, ci = 0.9, lwd = 2, threshold = 0.3,
                                     vertex.label.color1 = 'transparent', vertex.label.font1 = NA, 
                                     vertex.color1 = 'transparent',
                                     vertex.size1 = 1, edge.color1 = rgb(0, 0, 0, 1),
                                     vertex.label.color2 = 'black', vertex.label.font2 = 'Helvetica', 
                                     vertex.color2 = 'seagreen1',
                                     vertex.size2 = 8, edge.color2 = rgb(0, 0, 0, 0.3))
{
  edgelist <- get_edgelist(obj, ci = ci, transform = TRUE)
  threshold_edges <- edgelist[edgelist$median >= threshold,]
  net <- igraph::graph_from_edgelist(as.matrix(threshold_edges[, 1:2]))
  md <- threshold_edges[, 3]
  ub <- threshold_edges[, 5]
  coords <- igraph::layout_nicely(net)
  igraph::plot.igraph(net, layout = coords,
                      vertex.label.color = vertex.label.color1, label.family = vertex.label.font1, 
                      vertex.color = vertex.color1, vertex.size = vertex.size1,
                      edge.color = rgb(0, 0, 0, 1), edge.arrow.size = 0, edge.width = md * lwd)
  igraph::plot.igraph(net, layout = coords,
                      vertex.label.color = vertex.label.color2, label.family = vertex.label.font2,
                      vertex.color = vertex.color2, vertex.size = vertex.size2,
                      edge.color = edge.color1, edge.arrow.size = 0, edge.width = ub * lwd,
                      add = TRUE)
  if (obj$directed) {
    igraph::plot.igraph(net, edge.width = md * 0, layout = coords, 
                        vertex.label.color = vertex.label.color2, vertex.color = vertex.color2, 
                        edge.color = edge.color2, add = TRUE)
  }
}

### plot network
plot_network_threshold2(obj = mpnp2_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## compare edge weight distributions to simple SRI ####
# load in workspace image with edge weights already calculated
#load('mpnp_bisonr_edgescalculated.RData')

### plot against SRI
head(edgelist)
colnames(edgelist)[1:2] <- colnames(counts_df_model)[1:2]
edgelist$node_1_id <- as.integer(edgelist$node_1_id) ; edgelist$node_2_id <- as.integer(edgelist$node_2_id)
summary <- left_join(edgelist, counts_df_model, by = c('node_1_id','node_2_id'))
counts <- counts_df[,c('node_1_males','node_2_males','count_1','count_2')]
colnames(counts)[1:2] <- c('node_1_id','node_2_id')
summary <- left_join(summary, counts, by = c('node_1_id','node_2_id'))
summary$sri <- summary$event / (summary$duration)

plot(density(summary$sri), main = 'SRI vs model output: blue=all,\nred=both seen 8 times, green=both 12 times')
lines(density(summary$median), col = 'blue')
lines(density(summary$median[which(summary$count_1 >= 8 & summary$count_2 >= 8)]), col = 'red')
lines(density(summary$median[which(summary$count_1 >= 12 & summary$count_2 >= 12)]), col = 'green')

# try plotting a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
dyads <- counts_df[,c('dyad_id','node_1_males','node_2_males')]
colnames(dyads) <- c('dyad_id','node_1_id','node_2_id')
dyads <- left_join(dyads, counts_df_model, by = c('node_1_id','node_2_id'))
length(which(is.na(dyads$duration) == TRUE))

draws <- as.data.frame(mpnp2_edge_weights$chain) %>% pivot_longer(cols = everything(), names_to = 'dyad_model', values_to = 'edge')
draws$dyad_id <- rep(counts_df$dyad_id, 4000) ## IS THIS RIGHT??
draws$weight <- gtools::inv.logit(draws$edge)
draws$draw <- rep(1:4000,  each = nrow(counts_df_model))
draws <- left_join(draws, dyads, by = 'dyad_id')
draws$sri <- draws$event / draws$duration

#set.seed(15)
#subset_draws <- draws[draws$dyad_id %in% sample(draws$dyad_id, 150, replace = F),]
#subset_draws$median <- NA
#for(i in 1:nrow(subset_draws)){
#  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]
#  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), #subset_draws$median[i])
#}
#head(subset_draws)
#which(is.na(subset_draws$median) == TRUE)[1]

#subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)
#ggplot(data = subset_draws, mapping = aes(x = weight))+
#  geom_density(colour = 'blue')+
#  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+
#  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
#  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

#write_csv(subset_draws, '../data_processed/mpnp2_sampledyads_random_binary_vs_sri.csv')

subset_draws <- draws[draws$sri > 0.2,]
subset_draws$median <- NA
for(i in 1:nrow(subset_draws)){
  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]
  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), subset_draws$median[i])
}
head(subset_draws)
which(is.na(subset_draws$median) == TRUE)[1]

subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)
ggplot(data = subset_draws, mapping = aes(x = weight))+
  geom_density(colour = 'blue')+
  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15)+
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

ggplot(data = subset_draws, mapping = aes(x = weight))+
  geom_density(colour = 'blue')+
  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

write_csv(subset_draws, '../data_processed/mpnp2_sampledyads_sri0.2_binary_vs_sri.csv')

# clean environment
rm(draws, dyads, priors, subset_draws, x) ; gc()
dev.off()

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(mpnp2_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[,c('node_1_males','node_2_males','event_count','count_dyad','id_1','id_2','dyad_males')]
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ms_chain <- as.data.frame(mpnp2_edge_weights$chain)
colnames(ms_chain) <- counts_df_model$dyad_id
ms_chain <- pivot_longer(ms_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ms_chain$chain_position <- rep(1:length(unique(ms_chain$dyad_id)), each = 4000)
ms_chain$draw <- LaplacesDemon::invlogit(ms_chain$draw)
ms_chain$mean <- NA
hist(ms_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ms_chain$dyad_id))){
  x <- ms_chain[ms_chain$dyad_id == i,]
  ms_chain$mean <- ifelse(ms_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ms_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ms_chain$draw), lwd = 2)

quantile(ms_chain$draw, 0.98)

ms_edgelist <- bisonR::get_edgelist(mpnp2_edge_weights)
quantile(ms_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### save pdf
dev.off()

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
rm(list= ls()[!(ls() %in% c('mpnp2_edge_weights','mpnp2_edges_null',
                            'counts_df','counts_df_model',
                            'gbi_males','m_mat',
                            #'random_nodes',
                            'random_networks'))])

### write over saved workspace image -- only keep the things you might need later that take a long time to run
save.image('mpnp2_bisonr_edgescalculated.RData')

### end pdf
dev.off()

#### Period 3 ####
# create file of output graphs
pdf('../outputs/mpnp3_networkplots.pdf', width = 10, height = 10)

## create data list ####
counts_df3 <- read_csv('../data_processed/mpnp_period3_pairwiseevents.csv')
str(counts_df3)

### load in ages 
mpnp3_ages <- readRDS('../data_processed/mpnp3_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

### filter counts data frame down to males only
counts_df3 <- counts_df3 %>% 
  filter(id_1 %in% unique(mpnp3_ages$id)) %>% 
  filter(id_2 %in% unique(mpnp3_ages$id))

# reassign dyad numbers to remove gaps
counts_df3$node_1_males <- as.integer(as.factor(counts_df3$node_1))
counts_df3$node_2_males <- as.integer(as.factor(counts_df3$node_2))+1

### standardise dyad_id
counts_df3$dyad_males <- as.integer(as.factor(counts_df3$dyad_id))

## edge weights ####
### create data frame for edge weight model
counts_df3_model <- counts_df3[, c('node_1_males','node_2_males','together','count_dyad')] %>% distinct()
colnames(counts_df3_model) <- c('node_1_id','node_2_id','event','duration')

### set priors
priors <- get_default_priors('binary')
priors$edge <- 'normal(-2.5, 1.5)'    # slightly right skewed so that more density is towards 0 but still very weak and all values are possible
priors$fixed <- 'normal(0, 1.5)'  # slightly less wide that before, but zero centred so allows age to have a positive or negative impact on centrality and edge weight
prior_check(priors, 'binary')

### run edge weight model
mpnp3_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df3_model, 
  model_type = "binary",
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(mpnp3_edge_weights, par_ids = 2)
plot_predictions(mpnp3_edge_weights, num_draws = 20, type = "density")
plot_predictions(mpnp3_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(mpnp3_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(mpnp3_edge_weights)

### save workspace image for reloading at a later date that doesn't require running model again
save.image('mpnp_bisonr_edgescalculated.RData')

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### load in workspace image with edge weights already calculated
#load('mpnp_bisonr_edgescalculated.RData')

### run null model
mpnp3_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df3_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR comparison
model_comparison(list(non_random_model = mpnp3_edge_weights, random_model = mpnp3_edges_null))

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging <- function(models) {
  loos <- lapply(models, function(model) loo::loo(model$log_lik, r_eff=NA))
  
  results_matrix <- loo::loo_model_weights(loos, method="pseudobma")
  if (!is.null(names(models))) {
    names(results_matrix) <- names(models)
  }
  results_matrix
}
model_averaging(models = list(non_random_model = mpnp3_edge_weights, random_model = mpnp3_edges_null))

# save workspace image
save.image('mpnp_bisonr_randomnetwork.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
### adapt bisonr plot_network function to give more flexibility over plotting options
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

### load in ages 
#mpnp3_ages <- readRDS('../data_processed/mpnp3_ageestimates_mcmcoutput.rds') %>%
#  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

### create nodes data frame
nodes <- data.frame(bull = sort(unique(mpnp_ages$id)),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  x <- mpnp_ages[mpnp_ages$id == nodes$bull[i],]
  nodes$age[i] <- mean(x$age)
  if(nodes$bull[i] != 'M99'){
    y <- counts_df3[counts_df3$id_1 == nodes$bull[i], c('id_1','count_1')]
    nodes$sightings[i] <- y[1,2]
  }
}
nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

### plot network
plot_network_threshold(mpnp3_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')

### adapt to remove unconnected nodes
plot_network_threshold2 <- function (obj, ci = 0.9, lwd = 2, threshold = 0.3,
                                     vertex.label.color1 = 'transparent', vertex.label.font1 = NA, 
                                     vertex.color1 = 'transparent',
                                     vertex.size1 = 1, edge.color1 = rgb(0, 0, 0, 1),
                                     vertex.label.color2 = 'black', vertex.label.font2 = 'Helvetica', 
                                     vertex.color2 = 'seagreen1',
                                     vertex.size2 = 8, edge.color2 = rgb(0, 0, 0, 0.3))
{
  edgelist <- get_edgelist(obj, ci = ci, transform = TRUE)
  threshold_edges <- edgelist[edgelist$median >= threshold,]
  net <- igraph::graph_from_edgelist(as.matrix(threshold_edges[, 1:2]))
  md <- threshold_edges[, 3]
  ub <- threshold_edges[, 5]
  coords <- igraph::layout_nicely(net)
  igraph::plot.igraph(net, layout = coords,
                      vertex.label.color = vertex.label.color1, label.family = vertex.label.font1, 
                      vertex.color = vertex.color1, vertex.size = vertex.size1,
                      edge.color = rgb(0, 0, 0, 1), edge.arrow.size = 0, edge.width = md * lwd)
  igraph::plot.igraph(net, layout = coords,
                      vertex.label.color = vertex.label.color2, label.family = vertex.label.font2,
                      vertex.color = vertex.color2, vertex.size = vertex.size2,
                      edge.color = edge.color1, edge.arrow.size = 0, edge.width = ub * lwd,
                      add = TRUE)
  if (obj$directed) {
    igraph::plot.igraph(net, edge.width = md * 0, layout = coords, 
                        vertex.label.color = vertex.label.color2, vertex.color = vertex.color2, 
                        edge.color = edge.color2, add = TRUE)
  }
}

### plot network
plot_network_threshold2(obj = mpnp3_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## compare edge weight distributions to simple SRI ####
# load in workspace image with edge weights already calculated
#load('mpnp_bisonr_edgescalculated.RData')

### plot against SRI
head(edgelist)
colnames(edgelist)[1:2] <- colnames(counts_df_model)[1:2]
edgelist$node_1_id <- as.integer(edgelist$node_1_id) ; edgelist$node_2_id <- as.integer(edgelist$node_2_id)
summary <- left_join(edgelist, counts_df_model, by = c('node_1_id','node_2_id'))
counts <- counts_df[,c('node_1_males','node_2_males','count_1','count_2')]
colnames(counts)[1:2] <- c('node_1_id','node_2_id')
summary <- left_join(summary, counts, by = c('node_1_id','node_2_id'))
summary$sri <- summary$event / (summary$duration)

plot(density(summary$sri), main = 'SRI vs model output: blue=all,\nred=both seen 8 times, green=both 12 times')
lines(density(summary$median), col = 'blue')
lines(density(summary$median[which(summary$count_1 >= 8 & summary$count_2 >= 8)]), col = 'red')
lines(density(summary$median[which(summary$count_1 >= 12 & summary$count_2 >= 12)]), col = 'green')

# try plotting a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
dyads <- counts_df[,c('dyad_id','node_1_males','node_2_males')]
colnames(dyads) <- c('dyad_id','node_1_id','node_2_id')
dyads <- left_join(dyads, counts_df_model, by = c('node_1_id','node_2_id'))
length(which(is.na(dyads$duration) == TRUE))

draws <- as.data.frame(mpnp3_edge_weights$chain) %>% pivot_longer(cols = everything(), names_to = 'dyad_model', values_to = 'edge')
draws$dyad_id <- rep(counts_df$dyad_id, 4000) ## IS THIS RIGHT??
draws$weight <- gtools::inv.logit(draws$edge)
draws$draw <- rep(1:4000,  each = nrow(counts_df_model))
draws <- left_join(draws, dyads, by = 'dyad_id')
draws$sri <- draws$event / draws$duration

#set.seed(15)
#subset_draws <- draws[draws$dyad_id %in% sample(draws$dyad_id, 150, replace = F),]
#subset_draws$median <- NA
#for(i in 1:nrow(subset_draws)){
#  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]
#  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), #subset_draws$median[i])
#}
#head(subset_draws)
#which(is.na(subset_draws$median) == TRUE)[1]

#subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)
#ggplot(data = subset_draws, mapping = aes(x = weight))+
#  geom_density(colour = 'blue')+
#  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+
#  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
#  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

#write_csv(subset_draws, '../data_processed/mpnp3_sampledyads_random_binary_vs_sri.csv')

subset_draws <- draws[draws$sri > 0.2,]
subset_draws$median <- NA
for(i in 1:nrow(subset_draws)){
  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]
  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), subset_draws$median[i])
}
head(subset_draws)
which(is.na(subset_draws$median) == TRUE)[1]

subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)
ggplot(data = subset_draws, mapping = aes(x = weight))+
  geom_density(colour = 'blue')+
  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15)+
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

ggplot(data = subset_draws, mapping = aes(x = weight))+
  geom_density(colour = 'blue')+
  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

write_csv(subset_draws, '../data_processed/mpnp3_sampledyads_sri0.2_binary_vs_sri.csv')

# clean environment
rm(draws, dyads, priors, subset_draws, x) ; gc()
dev.off()

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(mpnp3_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[,c('node_1_males','node_2_males','event_count','count_dyad','id_1','id_2','dyad_males')]
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ms_chain <- as.data.frame(mpnp3_edge_weights$chain)
colnames(ms_chain) <- counts_df_model$dyad_id
ms_chain <- pivot_longer(ms_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ms_chain$chain_position <- rep(1:length(unique(ms_chain$dyad_id)), each = 4000)
ms_chain$draw <- LaplacesDemon::invlogit(ms_chain$draw)
ms_chain$mean <- NA
hist(ms_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ms_chain$dyad_id))){
  x <- ms_chain[ms_chain$dyad_id == i,]
  ms_chain$mean <- ifelse(ms_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ms_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ms_chain$draw), lwd = 2)

quantile(ms_chain$draw, 0.98)

ms_edgelist <- bisonR::get_edgelist(mpnp3_edge_weights)
quantile(ms_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### save pdf
dev.off()

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
rm(list= ls()[!(ls() %in% c('mpnp3_edge_weights','mpnp3_edges_null',
                            'counts_df','counts_df_model',
                            'gbi_males','m_mat',
                            #'random_nodes',
                            'random_networks'))])

### write over saved workspace image -- only keep the things you might need later that take a long time to run
save.image('mpnp3_bisonr_edgescalculated.RData')

### end pdf
dev.off()

#### Period 4 ####
# create file of output graphs
pdf('../outputs/mpnp4_networkplots.pdf', width = 10, height = 10)

## create data list ####
counts_df4 <- read_csv('../data_processed/mpnp_period4_pairwiseevents.csv')
str(counts_df4)

### load in ages 
mpnp4_ages <- readRDS('../data_processed/mpnp4_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

### filter counts data frame down to males only
counts_df4 <- counts_df4 %>% 
  filter(id_1 %in% unique(mpnp4_ages$id)) %>% 
  filter(id_2 %in% unique(mpnp4_ages$id))

# reassign dyad numbers to remove gaps
counts_df4$node_1_males <- as.integer(as.factor(counts_df4$node_1))
counts_df4$node_2_males <- as.integer(as.factor(counts_df4$node_2))+1

### standardise dyad_id
counts_df4$dyad_males <- as.integer(as.factor(counts_df4$dyad_id))

## edge weights ####
### create data frame for edge weight model
counts_df4_model <- counts_df4[, c('node_1_males','node_2_males','together','count_dyad')] %>% distinct()
colnames(counts_df4_model) <- c('node_1_id','node_2_id','event','duration')

### set priors
priors <- get_default_priors('binary')
priors$edge <- 'normal(-2.5, 1.5)'    # slightly right skewed so that more density is towards 0 but still very weak and all values are possible
priors$fixed <- 'normal(0, 1.5)'  # slightly less wide that before, but zero centred so allows age to have a positive or negative impact on centrality and edge weight
prior_check(priors, 'binary')

### run edge weight model
mpnp4_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df4_model, 
  model_type = "binary",
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(mpnp4_edge_weights, par_ids = 2)
plot_predictions(mpnp4_edge_weights, num_draws = 20, type = "density")
plot_predictions(mpnp4_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(mpnp4_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(mpnp4_edge_weights)

### save workspace image for reloading at a later date that doesn't require running model again
save.image('mpnp_bisonr_edgescalculated.RData')

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### load in workspace image with edge weights already calculated
#load('mpnp_bisonr_edgescalculated.RData')

### run null model
mpnp4_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df4_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR comparison
model_comparison(list(non_random_model = mpnp4_edge_weights, random_model = mpnp4_edges_null))

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging <- function(models) {
  loos <- lapply(models, function(model) loo::loo(model$log_lik, r_eff=NA))
  
  results_matrix <- loo::loo_model_weights(loos, method="pseudobma")
  if (!is.null(names(models))) {
    names(results_matrix) <- names(models)
  }
  results_matrix
}
model_averaging(models = list(non_random_model = mpnp4_edge_weights, random_model = mpnp4_edges_null))

# save workspace image
save.image('mpnp_bisonr_randomnetwork.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
### adapt bisonr plot_network function to give more flexibility over plotting options
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

### load in ages 
#mpnp4_ages <- readRDS('../data_processed/mpnp4_ageestimates_mcmcoutput.rds') %>%
#  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

### create nodes data frame
nodes <- data.frame(bull = sort(unique(mpnp_ages$id)),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  x <- mpnp_ages[mpnp_ages$id == nodes$bull[i],]
  nodes$age[i] <- mean(x$age)
  if(nodes$bull[i] != 'M99'){
    y <- counts_df4[counts_df4$id_1 == nodes$bull[i], c('id_1','count_1')]
    nodes$sightings[i] <- y[1,2]
  }
}
nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

### plot network
plot_network_threshold(mpnp4_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')

### adapt to remove unconnected nodes
plot_network_threshold2 <- function (obj, ci = 0.9, lwd = 2, threshold = 0.3,
                                     vertex.label.color1 = 'transparent', vertex.label.font1 = NA, 
                                     vertex.color1 = 'transparent',
                                     vertex.size1 = 1, edge.color1 = rgb(0, 0, 0, 1),
                                     vertex.label.color2 = 'black', vertex.label.font2 = 'Helvetica', 
                                     vertex.color2 = 'seagreen1',
                                     vertex.size2 = 8, edge.color2 = rgb(0, 0, 0, 0.3))
{
  edgelist <- get_edgelist(obj, ci = ci, transform = TRUE)
  threshold_edges <- edgelist[edgelist$median >= threshold,]
  net <- igraph::graph_from_edgelist(as.matrix(threshold_edges[, 1:2]))
  md <- threshold_edges[, 3]
  ub <- threshold_edges[, 5]
  coords <- igraph::layout_nicely(net)
  igraph::plot.igraph(net, layout = coords,
                      vertex.label.color = vertex.label.color1, label.family = vertex.label.font1, 
                      vertex.color = vertex.color1, vertex.size = vertex.size1,
                      edge.color = rgb(0, 0, 0, 1), edge.arrow.size = 0, edge.width = md * lwd)
  igraph::plot.igraph(net, layout = coords,
                      vertex.label.color = vertex.label.color2, label.family = vertex.label.font2,
                      vertex.color = vertex.color2, vertex.size = vertex.size2,
                      edge.color = edge.color1, edge.arrow.size = 0, edge.width = ub * lwd,
                      add = TRUE)
  if (obj$directed) {
    igraph::plot.igraph(net, edge.width = md * 0, layout = coords, 
                        vertex.label.color = vertex.label.color2, vertex.color = vertex.color2, 
                        edge.color = edge.color2, add = TRUE)
  }
}

### plot network
plot_network_threshold2(obj = mpnp4_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## compare edge weight distributions to simple SRI ####
# load in workspace image with edge weights already calculated
#load('mpnp_bisonr_edgescalculated.RData')

### plot against SRI
head(edgelist)
colnames(edgelist)[1:2] <- colnames(counts_df_model)[1:2]
edgelist$node_1_id <- as.integer(edgelist$node_1_id) ; edgelist$node_2_id <- as.integer(edgelist$node_2_id)
summary <- left_join(edgelist, counts_df_model, by = c('node_1_id','node_2_id'))
counts <- counts_df[,c('node_1_males','node_2_males','count_1','count_2')]
colnames(counts)[1:2] <- c('node_1_id','node_2_id')
summary <- left_join(summary, counts, by = c('node_1_id','node_2_id'))
summary$sri <- summary$event / (summary$duration)

plot(density(summary$sri), main = 'SRI vs model output: blue=all,\nred=both seen 8 times, green=both 12 times')
lines(density(summary$median), col = 'blue')
lines(density(summary$median[which(summary$count_1 >= 8 & summary$count_2 >= 8)]), col = 'red')
lines(density(summary$median[which(summary$count_1 >= 12 & summary$count_2 >= 12)]), col = 'green')

# try plotting a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
dyads <- counts_df[,c('dyad_id','node_1_males','node_2_males')]
colnames(dyads) <- c('dyad_id','node_1_id','node_2_id')
dyads <- left_join(dyads, counts_df_model, by = c('node_1_id','node_2_id'))
length(which(is.na(dyads$duration) == TRUE))

draws <- as.data.frame(mpnp4_edge_weights$chain) %>% pivot_longer(cols = everything(), names_to = 'dyad_model', values_to = 'edge')
draws$dyad_id <- rep(counts_df$dyad_id, 4000) ## IS THIS RIGHT??
draws$weight <- gtools::inv.logit(draws$edge)
draws$draw <- rep(1:4000,  each = nrow(counts_df_model))
draws <- left_join(draws, dyads, by = 'dyad_id')
draws$sri <- draws$event / draws$duration

#set.seed(15)
#subset_draws <- draws[draws$dyad_id %in% sample(draws$dyad_id, 150, replace = F),]
#subset_draws$median <- NA
#for(i in 1:nrow(subset_draws)){
#  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]
#  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), #subset_draws$median[i])
#}
#head(subset_draws)
#which(is.na(subset_draws$median) == TRUE)[1]

#subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)
#ggplot(data = subset_draws, mapping = aes(x = weight))+
#  geom_density(colour = 'blue')+
#  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+
#  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
#  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

#write_csv(subset_draws, '../data_processed/mpnp4_sampledyads_random_binary_vs_sri.csv')

subset_draws <- draws[draws$sri > 0.2,]
subset_draws$median <- NA
for(i in 1:nrow(subset_draws)){
  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]
  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), subset_draws$median[i])
}
head(subset_draws)
which(is.na(subset_draws$median) == TRUE)[1]

subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)
ggplot(data = subset_draws, mapping = aes(x = weight))+
  geom_density(colour = 'blue')+
  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15)+
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

ggplot(data = subset_draws, mapping = aes(x = weight))+
  geom_density(colour = 'blue')+
  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

write_csv(subset_draws, '../data_processed/mpnp4_sampledyads_sri0.2_binary_vs_sri.csv')

# clean environment
rm(draws, dyads, priors, subset_draws, x) ; gc()
dev.off()

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(mpnp4_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[,c('node_1_males','node_2_males','event_count','count_dyad','id_1','id_2','dyad_males')]
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ms_chain <- as.data.frame(mpnp4_edge_weights$chain)
colnames(ms_chain) <- counts_df_model$dyad_id
ms_chain <- pivot_longer(ms_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ms_chain$chain_position <- rep(1:length(unique(ms_chain$dyad_id)), each = 4000)
ms_chain$draw <- LaplacesDemon::invlogit(ms_chain$draw)
ms_chain$mean <- NA
hist(ms_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ms_chain$dyad_id))){
  x <- ms_chain[ms_chain$dyad_id == i,]
  ms_chain$mean <- ifelse(ms_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ms_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ms_chain$draw), lwd = 2)

quantile(ms_chain$draw, 0.98)

ms_edgelist <- bisonR::get_edgelist(mpnp4_edge_weights)
quantile(ms_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### save pdf
dev.off()

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
rm(list= ls()[!(ls() %in% c('mpnp4_edge_weights','mpnp4_edges_null',
                            'counts_df','counts_df_model',
                            'gbi_males','m_mat',
                            #'random_nodes',
                            'random_networks'))])

### write over saved workspace image -- only keep the things you might need later that take a long time to run
save.image('mpnp4_bisonr_edgescalculated.RData')

### end pdf
dev.off()

#### Period 5 ####
# create file of output graphs
pdf('../outputs/mpnp5_networkplots.pdf', width = 10, height = 10)

## create data list ####
counts_df5 <- read_csv('../data_processed/mpnp_period5_pairwiseevents.csv')
str(counts_df5)

### load in ages 
mpnp5_ages <- readRDS('../data_processed/mpnp5_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

### filter counts data frame down to males only
counts_df5 <- counts_df5 %>% 
  filter(id_1 %in% unique(mpnp5_ages$id)) %>% 
  filter(id_2 %in% unique(mpnp5_ages$id))

# reassign dyad numbers to remove gaps
counts_df5$node_1_males <- as.integer(as.factor(counts_df5$node_1))
counts_df5$node_2_males <- as.integer(as.factor(counts_df5$node_2))+1

### standardise dyad_id
counts_df5$dyad_males <- as.integer(as.factor(counts_df5$dyad_id))

## edge weights ####
### create data frame for edge weight model
counts_df5_model <- counts_df5[, c('node_1_males','node_2_males','together','count_dyad')] %>% distinct()
colnames(counts_df5_model) <- c('node_1_id','node_2_id','event','duration')

### set priors
priors <- get_default_priors('binary')
priors$edge <- 'normal(-2.5, 1.5)'    # slightly right skewed so that more density is towards 0 but still very weak and all values are possible
priors$fixed <- 'normal(0, 1.5)'  # slightly less wide that before, but zero centred so allows age to have a positive or negative impact on centrality and edge weight
prior_check(priors, 'binary')

### run edge weight model
mpnp5_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df5_model, 
  model_type = "binary",
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(mpnp5_edge_weights, par_ids = 2)
plot_predictions(mpnp5_edge_weights, num_draws = 20, type = "density")
plot_predictions(mpnp5_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(mpnp5_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(mpnp5_edge_weights)

### save workspace image for reloading at a later date that doesn't require running model again
save.image('mpnp_bisonr_edgescalculated.RData')

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### load in workspace image with edge weights already calculated
#load('mpnp_bisonr_edgescalculated.RData')

### run null model
mpnp5_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df5_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR comparison
model_comparison(list(non_random_model = mpnp5_edge_weights, random_model = mpnp5_edges_null))

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging <- function(models) {
  loos <- lapply(models, function(model) loo::loo(model$log_lik, r_eff=NA))
  
  results_matrix <- loo::loo_model_weights(loos, method="pseudobma")
  if (!is.null(names(models))) {
    names(results_matrix) <- names(models)
  }
  results_matrix
}
model_averaging(models = list(non_random_model = mpnp5_edge_weights, random_model = mpnp5_edges_null))

# save workspace image
save.image('mpnp_bisonr_randomnetwork.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
### adapt bisonr plot_network function to give more flexibility over plotting options
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

### load in ages 
#mpnp5_ages <- readRDS('../data_processed/mpnp5_ageestimates_mcmcoutput.rds') %>%
#  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

### create nodes data frame
nodes <- data.frame(bull = sort(unique(mpnp_ages$id)),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  x <- mpnp_ages[mpnp_ages$id == nodes$bull[i],]
  nodes$age[i] <- mean(x$age)
  if(nodes$bull[i] != 'M99'){
    y <- counts_df5[counts_df5$id_1 == nodes$bull[i], c('id_1','count_1')]
    nodes$sightings[i] <- y[1,2]
  }
}
nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

### plot network
plot_network_threshold(mpnp5_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')

### adapt to remove unconnected nodes
plot_network_threshold2 <- function (obj, ci = 0.9, lwd = 2, threshold = 0.3,
                                     vertex.label.color1 = 'transparent', vertex.label.font1 = NA, 
                                     vertex.color1 = 'transparent',
                                     vertex.size1 = 1, edge.color1 = rgb(0, 0, 0, 1),
                                     vertex.label.color2 = 'black', vertex.label.font2 = 'Helvetica', 
                                     vertex.color2 = 'seagreen1',
                                     vertex.size2 = 8, edge.color2 = rgb(0, 0, 0, 0.3))
{
  edgelist <- get_edgelist(obj, ci = ci, transform = TRUE)
  threshold_edges <- edgelist[edgelist$median >= threshold,]
  net <- igraph::graph_from_edgelist(as.matrix(threshold_edges[, 1:2]))
  md <- threshold_edges[, 3]
  ub <- threshold_edges[, 5]
  coords <- igraph::layout_nicely(net)
  igraph::plot.igraph(net, layout = coords,
                      vertex.label.color = vertex.label.color1, label.family = vertex.label.font1, 
                      vertex.color = vertex.color1, vertex.size = vertex.size1,
                      edge.color = rgb(0, 0, 0, 1), edge.arrow.size = 0, edge.width = md * lwd)
  igraph::plot.igraph(net, layout = coords,
                      vertex.label.color = vertex.label.color2, label.family = vertex.label.font2,
                      vertex.color = vertex.color2, vertex.size = vertex.size2,
                      edge.color = edge.color1, edge.arrow.size = 0, edge.width = ub * lwd,
                      add = TRUE)
  if (obj$directed) {
    igraph::plot.igraph(net, edge.width = md * 0, layout = coords, 
                        vertex.label.color = vertex.label.color2, vertex.color = vertex.color2, 
                        edge.color = edge.color2, add = TRUE)
  }
}

### plot network
plot_network_threshold2(obj = mpnp5_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## compare edge weight distributions to simple SRI ####
# load in workspace image with edge weights already calculated
#load('mpnp_bisonr_edgescalculated.RData')

### plot against SRI
head(edgelist)
colnames(edgelist)[1:2] <- colnames(counts_df_model)[1:2]
edgelist$node_1_id <- as.integer(edgelist$node_1_id) ; edgelist$node_2_id <- as.integer(edgelist$node_2_id)
summary <- left_join(edgelist, counts_df_model, by = c('node_1_id','node_2_id'))
counts <- counts_df[,c('node_1_males','node_2_males','count_1','count_2')]
colnames(counts)[1:2] <- c('node_1_id','node_2_id')
summary <- left_join(summary, counts, by = c('node_1_id','node_2_id'))
summary$sri <- summary$event / (summary$duration)

plot(density(summary$sri), main = 'SRI vs model output: blue=all,\nred=both seen 8 times, green=both 12 times')
lines(density(summary$median), col = 'blue')
lines(density(summary$median[which(summary$count_1 >= 8 & summary$count_2 >= 8)]), col = 'red')
lines(density(summary$median[which(summary$count_1 >= 12 & summary$count_2 >= 12)]), col = 'green')

# try plotting a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
dyads <- counts_df[,c('dyad_id','node_1_males','node_2_males')]
colnames(dyads) <- c('dyad_id','node_1_id','node_2_id')
dyads <- left_join(dyads, counts_df_model, by = c('node_1_id','node_2_id'))
length(which(is.na(dyads$duration) == TRUE))

draws <- as.data.frame(mpnp5_edge_weights$chain) %>% pivot_longer(cols = everything(), names_to = 'dyad_model', values_to = 'edge')
draws$dyad_id <- rep(counts_df$dyad_id, 4000) ## IS THIS RIGHT??
draws$weight <- gtools::inv.logit(draws$edge)
draws$draw <- rep(1:4000,  each = nrow(counts_df_model))
draws <- left_join(draws, dyads, by = 'dyad_id')
draws$sri <- draws$event / draws$duration

#set.seed(15)
#subset_draws <- draws[draws$dyad_id %in% sample(draws$dyad_id, 150, replace = F),]
#subset_draws$median <- NA
#for(i in 1:nrow(subset_draws)){
#  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]
#  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), #subset_draws$median[i])
#}
#head(subset_draws)
#which(is.na(subset_draws$median) == TRUE)[1]

#subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)
#ggplot(data = subset_draws, mapping = aes(x = weight))+
#  geom_density(colour = 'blue')+
#  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+
#  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
#  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

#write_csv(subset_draws, '../data_processed/mpnp5_sampledyads_random_binary_vs_sri.csv')

subset_draws <- draws[draws$sri > 0.2,]
subset_draws$median <- NA
for(i in 1:nrow(subset_draws)){
  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]
  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), subset_draws$median[i])
}
head(subset_draws)
which(is.na(subset_draws$median) == TRUE)[1]

subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)
ggplot(data = subset_draws, mapping = aes(x = weight))+
  geom_density(colour = 'blue')+
  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15)+
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

ggplot(data = subset_draws, mapping = aes(x = weight))+
  geom_density(colour = 'blue')+
  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

write_csv(subset_draws, '../data_processed/mpnp5_sampledyads_sri0.2_binary_vs_sri.csv')

# clean environment
rm(draws, dyads, priors, subset_draws, x) ; gc()
dev.off()

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(mpnp5_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[,c('node_1_males','node_2_males','event_count','count_dyad','id_1','id_2','dyad_males')]
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ms_chain <- as.data.frame(mpnp5_edge_weights$chain)
colnames(ms_chain) <- counts_df_model$dyad_id
ms_chain <- pivot_longer(ms_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ms_chain$chain_position <- rep(1:length(unique(ms_chain$dyad_id)), each = 4000)
ms_chain$draw <- LaplacesDemon::invlogit(ms_chain$draw)
ms_chain$mean <- NA
hist(ms_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ms_chain$dyad_id))){
  x <- ms_chain[ms_chain$dyad_id == i,]
  ms_chain$mean <- ifelse(ms_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ms_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ms_chain$draw), lwd = 2)

quantile(ms_chain$draw, 0.98)

ms_edgelist <- bisonR::get_edgelist(mpnp5_edge_weights)
quantile(ms_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### save pdf
dev.off()

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
rm(list= ls()[!(ls() %in% c('mpnp5_edge_weights','mpnp5_edges_null',
                            'counts_df','counts_df_model',
                            'gbi_males','m_mat',
                            #'random_nodes',
                            'random_networks'))])

### write over saved workspace image -- only keep the things you might need later that take a long time to run
save.image('mpnp5_bisonr_edgescalculated.RData')

### end pdf
dev.off()
