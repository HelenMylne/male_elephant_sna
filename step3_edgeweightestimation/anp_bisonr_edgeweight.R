# CONVERT EVERYTHING TO A SINGLE OR JUST A FEW FUNCTIONS IF YOU CAN!!!
#### set up ####
# load packages
library(tidyverse) #, lib.loc = '../packages/') # library(tidyverse)
library(dplyr, lib.loc = '../packages/')     # library(dplyr)
#library(rstan, lib.loc = '../packages/')     # library(rstan)
library(cmdstanr, lib.loc = '../packages/')  # library(cmdstanr)
library(bisonR, lib.loc = '../packages/')    # library(bisonR)
library(asnipe, lib.loc = '../packages/')    # library(asnipe)
library(sna, lib.loc = '../packages/')       # library(sna)
library(raster, lib.loc = '../packages/')    # library(raster)

# information
sessionInfo()
R.Version()
rstan::stan_version()

# set seed
set.seed(12345)

# set cmdstanr path
#set_cmdstan_path('H:/rlibs/4.2.1/')

pdf(file = '../outputs/anp_edgeweight_setup.pdf')

#### create data lists ####
# read in dyad data
counts_df <- read_csv('../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_shortwindows_impossiblepairsremoved.csv')

# read in sightings data
ate <- read_csv('../data_processed/anp_sightings_rawcombined.csv')
head(ate)
table(ate$num_bulls)
ate$num_bulls_recount <- NA
for(i in 1:nrow(ate)){
  x <- ate[ate$obs_id == ate$obs_id[i],]
  ate$num_bulls_recount[i] <- nrow(x)
}
mean(ate$num_bulls_recount) ; sd(ate$num_bulls_recount)
mean(ate$num_bulls) ; sd(ate$num_bulls)
rm(x, i) ; gc()

# count number of males per sighting
counts <- as.data.frame(table(ate$id)) ; colnames(counts) <- c('id','count')
summary(counts$count) ; sd(counts$count)

# narrow down to sightings data of interest
sightings <- ate[,c('obs_id','obs_date','correct_time_hms','grid_code','obs_type_old')] %>% distinct()
length(unique(sightings$obs_id))         # 24386 -- duplicates due to wrong obs_type
sightings$distinct <- NA
for(i in 1:nrow(sightings)) {
  x <- sightings[sightings$obs_id == sightings$obs_id[i], ]
  sightings$distinct[i] <- nrow(x)
  if(is.na(sightings$obs_type_old[i]) == TRUE ) {
    obs_type <- unique(x$obs_type_old[!is.na(x$obs_type_old)])
    sightings$obs_type_old[i] <- ifelse(length(obs_type) == 1, obs_type, 'U')
  }
}
sightings <- distinct(sightings) # now the right length
table(sightings$obs_type_old)

# read in individual male data
males <- readxl::read_excel('../data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()
males$id <- paste0('M',males$casename)
ids <- sort(unique(ate$id))
males <- males %>% dplyr::filter(id %in% ids)

# counts per male per time period
length(unique(counts_df$period))
males$count_all <- NA
males$p1 <- NA ; males$p2 <- NA ; males$p3 <- NA ; males$p4 <- NA ; males$p5 <- NA ; males$p6 <- NA ; males$p7 <- NA
males$p8 <- NA ; males$p9 <- NA ; males$p10 <- NA ; males$p11 <- NA ; males$p12 <- NA ; males$p13 <- NA
males$p14 <- NA ; males$p15 <- NA ; males$p16 <- NA ; males$p17 <- NA ; males$p18 <- NA ; males$p19 <- NA
males$p20 <- NA ; males$p21 <- NA ; males$p22 <- NA ; males$p23 <- NA ; males$p24 <- NA ; males$p25 <- NA
males$p26 <- NA ; males$p27 <- NA ; males$p28 <- NA ; males$p29 <- NA ; males$p30 <- NA ; males$p31 <- NA
males$p32 <- NA ; males$p33 <- NA ; males$p34 <- NA ; males$p35 <- NA ; males$p36 <- NA
windows <- c(sort(unique(counts_df$period_start)),max(counts_df$period_end))

for(i in 1:nrow(males)){
  counts <- ate[ate$id == males$id[i],]
  males$count_all[i] <- nrow(counts)
  for(j in 1:max(counts_df$period)){
    males[i,j+22] <- length(which(counts$obs_date >= windows[j] & counts$obs_date < windows[j+1]))
  }
}

# time windows
periods <- data.frame(period = 1:(length(windows)-1),
                      period_start = windows[1:(length(windows)-1)],
                      period_end = windows[2:(length(windows))])

# check out total sightings per individual/dyad per period and see if they are reasonable
table(counts_df$period_count_1) ; hist(counts_df$period_count_1)
table(counts_df$period_count_2) ; hist(counts_df$period_count_2)
table(counts_df$period_count_dyad) ; hist(counts_df$period_count_dyad, breaks = 30)
table(counts_df$event_count) ; hist(counts_df$event_count, breaks = 20, ylim = c(0,15000))
# many are not seen much, but repeat sightings of pairs do seem to be genuinely rare

rm(ate, counts, sightings, i, x, j, obs_type)

### add time marker
print(paste0('data read in at ', Sys.time()))

#### set up priors and functions for all ####
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

#### Period == 1 ####
# create output pdf
pdf(file = '../outputs/anp1_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 1, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp1_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp1_edge_weights, par_ids = 2)
plot_predictions(anp1_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp1_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp1_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp1_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp1_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp1_edge_weights, random_model = anp1_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp1_edge_weights, random_model = anp1_edges_null))

# save workspace image
save.image('anp1_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 1 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 1 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp1_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp1_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## compare edge weight distributions to simple SRI ####
# load in workspace image with edge weights already calculated
#load('anp1_bisonr_edgescalculated.RData')

### plot against SRI
#head(edgelist)
#colnames(edgelist)[1:2] <- colnames(counts_df_model)[1:2]
#edgelist$node_1_id <- as.integer(edgelist$node_1_id) ; edgelist$node_2_id <- as.integer(edgelist$node_2_id)
#summary <- left_join(edgelist, counts_df_model, by = c('node_1_id','node_2_id'))
#counts <- counts_df[,c('node_1_males','node_2_males','count_1','count_2')]
#colnames(counts)[1:2] <- c('node_1_id','node_2_id')
#summary <- left_join(summary, counts, by = c('node_1_id','node_2_id'))
#summary$sri <- summary$event / (summary$duration)

#plot(density(summary$sri), main = 'SRI vs model output: blue=all,\nred=both seen 8 times, green=both 12 times')
#lines(density(summary$median), col = 'blue')
#lines(density(summary$median[which(summary$count_1 >= 8 & summary$count_2 >= 8)]), col = 'red')
#lines(density(summary$median[which(summary$count_1 >= 12 & summary$count_2 >= 12)]), col = 'green')

# try plotting a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
#dyads <- counts_df[,c('dyad_id','node_1_males','node_2_males')]
#colnames(dyads) <- c('dyad_id','node_1_id','node_2_id')
#dyads <- left_join(dyads, counts_df_model, by = c('node_1_id','node_2_id'))
#length(which(is.na(dyads$duration) == TRUE))

#draws <- as.data.frame(anp1_edge_weights$chain) %>% pivot_longer(cols = everything(), names_to = 'dyad_model', values_to = 'edge')
#draws$dyad_id <- rep(counts_df$dyad_id, 4000) ## IS THIS RIGHT??
#draws$weight <- gtools::inv.logit(draws$edge)
#draws$draw <- rep(1:4000,  each = nrow(counts_df_model))
#draws <- left_join(draws, dyads, by = 'dyad_id')
#draws$sri <- draws$event / draws$duration

#set.seed(15)
#subset_draws <- draws[draws$dyad_id %in% sample(draws$dyad_id, 150, replace = F),]
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
#  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+
#  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
#  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

#write_csv(subset_draws, '../data_processed/anp1_sampledyads_random_binary_vs_sri.csv')

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
#  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15)+
#  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
#  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

#ggplot(data = subset_draws, mapping = aes(x = weight))+
#  geom_density(colour = 'blue')+
#  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+
#  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
#  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

#write_csv(subset_draws, '../data_processed/anp1_sampledyads_sri0.2_binary_vs_sri.csv')

# clean environment
#rm(draws, dyads, priors, subset_draws, x) ; gc()
#dev.off()

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp1_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[,c('node_1_males','node_2_males','event_count','count_dyad','id_1','id_2','dyad_males')]
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp1_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp1_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp1_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### Period == 2 ####
# create output pdf
pdf(file = '../outputs/anp2_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 2, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp2_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp2_edge_weights, par_ids = 2)
plot_predictions(anp2_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp2_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp2_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp2_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp2_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp2_edge_weights, random_model = anp2_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp2_edge_weights, random_model = anp2_edges_null))

# save workspace image
save.image('anp2_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 2 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 2 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp2_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp2_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp2_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 2, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp2_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp2_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp2_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 3 ####
# create output pdf
pdf(file = '../outputs/anp3_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 3, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp3_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp3_edge_weights, par_ids = 2)
plot_predictions(anp3_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp3_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp3_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp3_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp3_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp3_edge_weights, random_model = anp3_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp3_edge_weights, random_model = anp3_edges_null))

# save workspace image
save.image('anp3_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 3 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 3 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp3_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp3_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp3_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 3, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp3_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp3_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp3_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 4 ####
# create output pdf
pdf(file = '../outputs/anp4_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 4, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp4_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp4_edge_weights, par_ids = 2)
plot_predictions(anp4_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp4_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp4_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp4_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp4_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp4_edge_weights, random_model = anp4_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp4_edge_weights, random_model = anp4_edges_null))

# save workspace image
save.image('anp4_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 4 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 4 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp4_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp4_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp4_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 4, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp4_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp4_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp4_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 5 ####
# create output pdf
pdf(file = '../outputs/anp5_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 5, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp5_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp5_edge_weights, par_ids = 2)
plot_predictions(anp5_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp5_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp5_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp5_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp5_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp5_edge_weights, random_model = anp5_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp5_edge_weights, random_model = anp5_edges_null))

# save workspace image
save.image('anp5_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 5 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 5 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp5_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp5_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp5_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 5, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp5_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp5_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp5_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 6 ####
# create output pdf
pdf(file = '../outputs/anp6_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 6, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp6_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp6_edge_weights, par_ids = 2)
plot_predictions(anp6_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp6_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp6_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp6_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp6_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp6_edge_weights, random_model = anp6_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp6_edge_weights, random_model = anp6_edges_null))

# save workspace image
save.image('anp6_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 6 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 6 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp6_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp6_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp6_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 6, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp6_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp6_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp6_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 7 ####
# create output pdf
pdf(file = '../outputs/anp7_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 7, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp7_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp7_edge_weights, par_ids = 2)
plot_predictions(anp7_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp7_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp7_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp7_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp7_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp7_edge_weights, random_model = anp7_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp7_edge_weights, random_model = anp7_edges_null))

# save workspace image
save.image('anp7_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 7 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 7 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp7_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp7_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp7_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 7, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp7_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp7_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp7_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 8 ####
# create output pdf
pdf(file = '../outputs/anp8_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 8, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp8_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp8_edge_weights, par_ids = 2)
plot_predictions(anp8_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp8_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp8_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp8_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp8_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp8_edge_weights, random_model = anp8_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp8_edge_weights, random_model = anp8_edges_null))

# save workspace image
save.image('anp8_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 8 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 8 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp8_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp8_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp8_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 8, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp8_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp8_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp8_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 9 ####
# create output pdf
pdf(file = '../outputs/anp9_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 9, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp9_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp9_edge_weights, par_ids = 2)
plot_predictions(anp9_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp9_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp9_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp9_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp9_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp9_edge_weights, random_model = anp9_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp9_edge_weights, random_model = anp9_edges_null))

# save workspace image
save.image('anp9_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 9 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 9 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp9_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp9_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp9_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 9, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp9_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp9_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp9_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 10 ####
# create output pdf
pdf(file = '../outputs/anp10_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 10, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp10_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp10_edge_weights, par_ids = 2)
plot_predictions(anp10_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp10_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp10_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp10_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp10_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp10_edge_weights, random_model = anp10_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp10_edge_weights, random_model = anp10_edges_null))

# save workspace image
save.image('anp10_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 10 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 10 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp10_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp10_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp10_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 10, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp10_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp10_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp10_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 11 ####
# create output pdf
pdf(file = '../outputs/anp11_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 11, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp11_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp11_edge_weights, par_ids = 2)
plot_predictions(anp11_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp11_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp11_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp11_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp11_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp11_edge_weights, random_model = anp11_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp11_edge_weights, random_model = anp11_edges_null))

# save workspace image
save.image('anp11_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 11 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 11 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp11_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp11_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp11_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 11, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp11_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp11_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp11_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 12 ####
# create output pdf
pdf(file = '../outputs/anp12_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 12, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp12_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp12_edge_weights, par_ids = 2)
plot_predictions(anp12_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp12_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp12_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp12_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp12_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp12_edge_weights, random_model = anp12_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp12_edge_weights, random_model = anp12_edges_null))

# save workspace image
save.image('anp12_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 12 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 12 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp12_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp12_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp12_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 12, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp12_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp12_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp12_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 13 ####
# create output pdf
pdf(file = '../outputs/anp13_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 13, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp13_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp13_edge_weights, par_ids = 2)
plot_predictions(anp13_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp13_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp13_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp13_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp13_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp13_edge_weights, random_model = anp13_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp13_edge_weights, random_model = anp13_edges_null))

# save workspace image
save.image('anp13_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 13 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 13 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp13_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp13_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp13_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 13, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp13_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp13_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp13_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 14 ####
# create output pdf
pdf(file = '../outputs/anp14_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 14, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp14_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp14_edge_weights, par_ids = 2)
plot_predictions(anp14_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp14_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp14_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp14_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp14_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp14_edge_weights, random_model = anp14_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp14_edge_weights, random_model = anp14_edges_null))

# save workspace image
save.image('anp14_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 14 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 14 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp14_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp14_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp14_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 14, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp14_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp14_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp14_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 15 ####
# create output pdf
pdf(file = '../outputs/anp15_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 15, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp15_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp15_edge_weights, par_ids = 2)
plot_predictions(anp15_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp15_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp15_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp15_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp15_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp15_edge_weights, random_model = anp15_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp15_edge_weights, random_model = anp15_edges_null))

# save workspace image
save.image('anp15_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 15 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 15 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp15_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp15_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp15_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 15, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp15_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp15_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp15_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 16 ####
# create output pdf
pdf(file = '../outputs/anp16_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 16, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp16_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp16_edge_weights, par_ids = 2)
plot_predictions(anp16_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp16_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp16_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp16_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp16_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp16_edge_weights, random_model = anp16_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp16_edge_weights, random_model = anp16_edges_null))

# save workspace image
save.image('anp16_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 16 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 16 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp16_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp16_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp16_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 16, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp16_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp16_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp16_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 17 ####
# create output pdf
pdf(file = '../outputs/anp17_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 17, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp17_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp17_edge_weights, par_ids = 2)
plot_predictions(anp17_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp17_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp17_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp17_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp17_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp17_edge_weights, random_model = anp17_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp17_edge_weights, random_model = anp17_edges_null))

# save workspace image
save.image('anp17_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 17 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 17 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp17_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp17_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp17_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 17, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp17_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp17_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp17_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 18 ####
# create output pdf
pdf(file = '../outputs/anp18_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 18, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp18_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp18_edge_weights, par_ids = 2)
plot_predictions(anp18_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp18_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp18_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp18_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp18_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp18_edge_weights, random_model = anp18_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp18_edge_weights, random_model = anp18_edges_null))

# save workspace image
save.image('anp18_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 18 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 18 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp18_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp18_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp18_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 18, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp18_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp18_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp18_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 19 ####
# create output pdf
pdf(file = '../outputs/anp19_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 19, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp19_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp19_edge_weights, par_ids = 2)
plot_predictions(anp19_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp19_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp19_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp19_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp19_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp19_edge_weights, random_model = anp19_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp19_edge_weights, random_model = anp19_edges_null))

# save workspace image
save.image('anp19_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 19 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 19 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp19_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp19_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp19_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 19, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp19_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp19_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp19_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 20 ####
# create output pdf
pdf(file = '../outputs/anp20_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 20, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp20_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp20_edge_weights, par_ids = 2)
plot_predictions(anp20_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp20_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp20_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp20_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp20_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp20_edge_weights, random_model = anp20_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp20_edge_weights, random_model = anp20_edges_null))

# save workspace image
save.image('anp20_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 20 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 20 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp20_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp20_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp20_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 20, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp20_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp20_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp20_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 21 ####
# create output pdf
pdf(file = '../outputs/anp21_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 21, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp21_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp21_edge_weights, par_ids = 2)
plot_predictions(anp21_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp21_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp21_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp21_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp21_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp21_edge_weights, random_model = anp21_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp21_edge_weights, random_model = anp21_edges_null))

# save workspace image
save.image('anp21_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 21 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 21 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp21_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp21_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp21_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 21, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp21_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp21_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp21_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 22 ####
# create output pdf
pdf(file = '../outputs/anp22_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 22, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp22_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp22_edge_weights, par_ids = 2)
plot_predictions(anp22_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp22_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp22_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp22_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp22_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp22_edge_weights, random_model = anp22_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp22_edge_weights, random_model = anp22_edges_null))

# save workspace image
save.image('anp22_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 22 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 22 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp22_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp22_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp22_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 22, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp22_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp22_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp22_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 23 ####
# create output pdf
pdf(file = '../outputs/anp23_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 23, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp23_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp23_edge_weights, par_ids = 2)
plot_predictions(anp23_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp23_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp23_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp23_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp23_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp23_edge_weights, random_model = anp23_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp23_edge_weights, random_model = anp23_edges_null))

# save workspace image
save.image('anp23_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 23 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 23 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp23_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp23_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp23_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 23, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp23_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp23_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp23_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 24 ####
# create output pdf
pdf(file = '../outputs/anp24_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 24, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp24_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp24_edge_weights, par_ids = 2)
plot_predictions(anp24_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp24_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp24_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp24_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp24_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp24_edge_weights, random_model = anp24_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp24_edge_weights, random_model = anp24_edges_null))

# save workspace image
save.image('anp24_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 24 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 24 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp24_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp24_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp24_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 24, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp24_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp24_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp24_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 25 ####
# create output pdf
pdf(file = '../outputs/anp25_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 25, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp25_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp25_edge_weights, par_ids = 2)
plot_predictions(anp25_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp25_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp25_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp25_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp25_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp25_edge_weights, random_model = anp25_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp25_edge_weights, random_model = anp25_edges_null))

# save workspace image
save.image('anp25_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 25 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 25 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp25_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp25_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp25_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 25, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp25_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp25_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp25_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 26 ####
# create output pdf
pdf(file = '../outputs/anp26_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 26, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp26_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp26_edge_weights, par_ids = 2)
plot_predictions(anp26_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp26_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp26_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp26_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp26_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp26_edge_weights, random_model = anp26_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp26_edge_weights, random_model = anp26_edges_null))

# save workspace image
save.image('anp26_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 26 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 26 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp26_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp26_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp26_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 26, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp26_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp26_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp26_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 27 ####
# create output pdf
pdf(file = '../outputs/anp27_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 27, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp27_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp27_edge_weights, par_ids = 2)
plot_predictions(anp27_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp27_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp27_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp27_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp27_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp27_edge_weights, random_model = anp27_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp27_edge_weights, random_model = anp27_edges_null))

# save workspace image
save.image('anp27_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 27 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 27 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp27_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp27_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp27_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 27, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp27_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp27_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp27_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 28 ####
# create output pdf
pdf(file = '../outputs/anp28_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 28, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp28_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp28_edge_weights, par_ids = 2)
plot_predictions(anp28_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp28_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp28_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp28_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp28_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp28_edge_weights, random_model = anp28_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp28_edge_weights, random_model = anp28_edges_null))

# save workspace image
save.image('anp28_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 28 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 28 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp28_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp28_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp28_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 28, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp28_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp28_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp28_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 29 ####
# create output pdf
pdf(file = '../outputs/anp29_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 29, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp29_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp29_edge_weights, par_ids = 2)
plot_predictions(anp29_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp29_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp29_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp29_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp29_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp29_edge_weights, random_model = anp29_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp29_edge_weights, random_model = anp29_edges_null))

# save workspace image
save.image('anp29_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 29 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 29 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp29_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp29_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp29_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 29, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp29_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp29_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp29_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 30 ####
# create output pdf
pdf(file = '../outputs/anp30_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 30, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp30_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp30_edge_weights, par_ids = 2)
plot_predictions(anp30_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp30_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp30_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp30_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp30_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp30_edge_weights, random_model = anp30_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp30_edge_weights, random_model = anp30_edges_null))

# save workspace image
save.image('anp30_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 30 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 30 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp30_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp30_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp30_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 30, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp30_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp30_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp30_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 31 ####
# create output pdf
pdf(file = '../outputs/anp31_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 31, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp31_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp31_edge_weights, par_ids = 2)
plot_predictions(anp31_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp31_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp31_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp31_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp31_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp31_edge_weights, random_model = anp31_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp31_edge_weights, random_model = anp31_edges_null))

# save workspace image
save.image('anp31_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 31 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 31 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp31_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp31_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp31_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 31, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp31_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp31_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp31_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 32 ####
# create output pdf
pdf(file = '../outputs/anp32_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 32, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp32_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp32_edge_weights, par_ids = 2)
plot_predictions(anp32_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp32_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp32_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp32_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp32_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp32_edge_weights, random_model = anp32_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp32_edge_weights, random_model = anp32_edges_null))

# save workspace image
save.image('anp32_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 32 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 32 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp32_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp32_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp32_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 32, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp32_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp32_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp32_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 33 ####
# create output pdf
pdf(file = '../outputs/anp33_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 33, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp33_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp33_edge_weights, par_ids = 2)
plot_predictions(anp33_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp33_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp33_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp33_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp33_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp33_edge_weights, random_model = anp33_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp33_edge_weights, random_model = anp33_edges_null))

# save workspace image
save.image('anp33_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 33 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 33 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp33_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp33_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp33_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 33, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp33_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp33_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp33_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 34 ####
# create output pdf
pdf(file = '../outputs/anp34_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 34, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp34_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp34_edge_weights, par_ids = 2)
plot_predictions(anp34_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp34_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp34_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp34_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp34_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp34_edge_weights, random_model = anp34_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp34_edge_weights, random_model = anp34_edges_null))

# save workspace image
save.image('anp34_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 34 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 34 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp34_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp34_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp34_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 34, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp34_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp34_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp34_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 35 ####
# create output pdf
pdf(file = '../outputs/anp35_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 35, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp35_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp35_edge_weights, par_ids = 2)
plot_predictions(anp35_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp35_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp35_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp35_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp35_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp35_edge_weights, random_model = anp35_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp35_edge_weights, random_model = anp35_edges_null))

# save workspace image
save.image('anp35_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 35 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 35 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp35_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp35_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp35_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 35, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp35_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp35_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp35_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

#### period == 36 ####
# create output pdf
pdf(file = '../outputs/anp36_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[counts_df$period == 36, c('node_1','node_2','event_count','period_count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### run edge weight model
anp36_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp36_edge_weights, par_ids = 2)
plot_predictions(anp36_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp36_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp36_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp36_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp36_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp36_edge_weights, random_model = anp36_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging(models = list(non_random_model = anp36_edge_weights, random_model = anp36_edges_null))

# save workspace image
save.image('anp36_bisonr_edgescalculated.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  df <- counts_df[(counts_df$period == 36 & counts_df$node_1 == nodes$bull[i]),
                  c('id_1','age_start_1','period_count_1')] %>% distinct()
  if(nrow(df) > 0){
    nodes$age[i] <- df$age_start_1[1]
    nodes$sightings[i] <- df$period_count_1[1]
  }
  else {
    df <- counts_df[(counts_df$period == 36 & counts_df$node_2 == nodes$bull[i]),
                    c('id_2','age_start_2','period_count_2')] %>% distinct()
    nodes$age[i] <- df$age_start_2[1]
    nodes$sightings[i] <- df$period_count_2[1]
  }
}

nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(anp36_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = anp36_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## coefficient of variation of edge weights (aka social differentiation) ####
# extract cv for model
global_cv <- extract_metric(anp36_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[counts_df$period == 36, c('node_1','node_2','event_count','period_count_dyad','id_1','id_2','dyad_id')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ew_chain <- as.data.frame(anp36_edge_weights$chain)
colnames(ew_chain) <- counts_df_model$dyad_id
ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
ew_chain$mean <- NA
hist(ew_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ew_chain$dyad_id))){
  x <- ew_chain[ew_chain$dyad_id == i,]
  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ew_chain$draw), lwd = 2)

quantile(ew_chain$draw, 0.98)

ew_edgelist <- bisonR::get_edgelist(anp36_edge_weights)
quantile(ew_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
### save pdf
dev.off()

### save workspace image
save.image('anp36_bisonr_edgescalculated.RData')

### clear environment
rm(list = ls()[!(ls() %in% c('counts_df','priors','periods','windows','males',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

