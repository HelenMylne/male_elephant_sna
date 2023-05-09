#### set up ####
# load packages
# library(tidyverse) ; library(dplyr) ; library(cmdstanr) ; library(bisonR) ; library(asnipe) ; library(sna) ; library(raster)
library(cmdstanr, lib.loc = '../packages/')      # library(cmdstanr)
library(tidyverse, lib.loc = '../packages/')     # library(tidyverse)
library(dplyr, lib.loc = '../packages/')         # library(dplyr)
#library(rstan, lib.loc = '../packages/')        # library(rstan)
library(extraDistr, lib.loc = '../packages/')    # library(extraDistr)
library(bisonR, lib.loc = '../packages/')        # library(bisonR)
library(asnipe, lib.loc = '../packages/')        # library(asnipe)
library(sna, lib.loc = '../packages/')           # library(sna)
library(raster, lib.loc = '../packages/')        # library(raster)
library(spatsoc, lib.loc = '../packages/')       # library(spatsoc)
library(LaplacesDemon, lib.loc = '../packages/') # library(LaplacesDemon)

# information
sessionInfo()
R.Version()
rstan::stan_version()

# set seed
set.seed(12345)

#### create data lists ####
### import data for aggregated model (binomial) -- counts of positive associations and total sightings
counts_df <- read_csv('../data_processed/motnp_binomialpairwiseevents_malesonly.csv')

### load in ages 
motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

### add time marker
print(paste0('data read in at ', Sys.time()))

#### edge weights -- binomial conjugate priors, males only ####
# add pdf output file
pdf(file = '../outputs/motnp_bisonr_edgeweight_widebinomconj.pdf', height = 10, width = 14)

### create data frame for edge weight model
counts_df_model <- counts_df[, c('node_1_males','node_2_males','event_count','count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### new wide prior
prior_alpha <- 0.7
prior_beta <- 3
plot(density( rbeta(50000, prior_alpha, prior_beta) ), xlim = c(0,1), col = 'red')

### set priors for dyads seen together at any point
priors <- get_default_priors('binary_conjugate') # obtain structure for bison model priors
priors$edge <- paste0('beta(',prior_alpha,',',prior_beta,')') 
prior_check(priors, 'binary_conjugate')

### create dummy variable indicating if ever seen together or not -- use to select which prior to draw edge weight from
counts_df_model$seen_together <- ifelse(counts_df_model$event > 0, 1, 0)

### run edge weight model
motnp_fit_edges <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id),   # count of sightings together given count of total sightings as a result of the individuals contained within the dyad
  data = counts_df_model,
  model_type = "binary_conjugate",
  priors = priors
)

### save workspace
save.image('motnp_bisonr_edgescalculated_widebinomconj_v2.RData')

### run diagnostic plots
#plot_trace(motnp_fit_edges, par_ids = 2)                            # trace plot -- not needed for conjugate model??
plot_predictions(motnp_fit_edges, num_draws = 20, type = "density")  # compare predictions to raw data -- predictions are more variable than observations, both overestimating the number of pairs only seen together a few times, but also allowing for together scores higher than observed
plot_predictions(motnp_fit_edges, num_draws = 20, type = "point")    # compare predictions to raw data -- anything below the line is predicting below SRI (e.g. upper end SHOULD be massively below line because we don't want scores of 1 from pairs seen once)

### compare edge weights to prior
edges <- as.data.frame(motnp_fit_edges$chain) %>%               # extract chain of values from model
  pivot_longer(cols = everything(), names_to = 'dyad', values_to = 'draw') %>%  # convert to long format
  mutate(dyad_id = rep(counts_df$dyad_id, 4000),                                # add column to identify dyad per draw
         draw = plogis(draw))                                                   # convert draws to proportion
head(edges)
plot(NULL, xlim = c(0,1), ylim = c(0,50), las = 1,                              # prepare plot window
     main = 'edge distributions',ylab = 'density', xlab = 'edge weight')
plot_seen <- sample(counts_df$dyad_id[which(counts_df$event_count > 0)], 500, replace = F)      # sample 5000 dyads that were seen together to plot
plot_unseen <- sample(counts_df$dyad_id[which(counts_df$event_count == 0)], 500, replace = F)   # sample 5000 dyads that were NOT seen together to plot
plot_samples <- c(plot_seen,plot_unseen)                                        # sample 10000 dyads to plot
for(dyad in plot_samples) {                                                     # plot probability density for sampled dyads
  x <- edges[edges$dyad_id == dyad,]
  lines(density(x$draw),
        col = ifelse(dyad %in% plot_seen, 'red','blue'))
}
lines(density(edges$draw), lwd = 3)                                             # add probability density line for all draws from all dyads

edges$chain_position <- rep(1:4000, each = length(unique(edges$dyad)))          # add column for position in chain
edges <- edges[edges$chain_position < 1001,]                                    # only take first 1000 values to speed up next plot
edges$mean <- NA
for(dyad in unique(edges$dyad)) {                                               # calculate mean edge per dyad
  if(is.na(edges$mean[dyad]) == TRUE){
    edges$mean[edges$dyad == dyad] <- mean(edges$draw[edges$dyad == dyad])
  }
}
plot_low <- edges[edges$mean == min(edges$mean),]                               # dyad with minimum mean edge weight
plot_q25 <- edges[edges$mean == quantile(edges$mean, 0.25),]                    # dyad with 25th percentile mean edge weight
plot_mid <- edges[edges$mean == quantile(edges$mean, 0.501),]                   # dyad with median mean edge weight (0.5 exactly couldn't identify a single dyad, but 0.501 found it)
plot_q75 <- edges[edges$mean == quantile(edges$mean, 0.75),]                    # dyad with 75th percentile mean edge weight
plot_q98 <- edges[edges$mean == quantile(edges$mean, 0.98),]                    # dyad with 98th percentile mean edge weight
plot_high <- edges[edges$mean == max(edges$mean),]                              # dyad with maximum mean edge weight
plot_data <- rbind(plot_low, plot_q25, plot_mid, plot_q75, plot_q98, plot_high) # combine all to single data frame

ggplot(plot_data,                                                               # plot distributions of selected dyads
       aes(x = draw, colour = as.factor(mean), fill = as.factor(mean)))+
  geom_density(linewidth = 1)+
  theme_classic()+
  scale_fill_viridis_d(alpha = 0.2)+                                            # colour distributions by dyad (translucent)
  scale_color_viridis_d()+                                                      # colour distributions by dyad (solid)
  theme(legend.position = 'none',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_x_continuous(name = 'edge weight', limits = c(0,1))

prior_plot <- data.frame(dyad = 'prior',                                        # create prior dataframe with draws from prior distribution, formatted to match columns in plot_data
                         draw = plogis(rbeta(5000,prior_alpha,prior_beta)),
                         dyad_id = 1000000,
                         mean = NA,
                         chain_position = 1:1000)
prior_plot$mean <- mean(prior_plot$draw)                                        # calculate mean of edge prior distribution
ggplot(plot_data,                                                               # plot distributions of selected dyads
       aes(x = draw, colour = as.factor(mean), fill = as.factor(mean)),
       #fill = viridis, colour = colours
)+
  geom_density(linewidth = 1, alpha = 0.2)+
  geom_density(data = prior_plot, linewidth = 1, colour = 'black', linetype = 2, alpha = 0)+  # add probability density line for prior istribution -- no fill, dashed line not solid, different colour scale to others
  theme_classic()+
  scale_fill_viridis_d(option = 'plasma',                                       # use alternative viridis colour pallette to avoid confusion regarding ages (all other plots, viridis scale = age)
                       alpha = 0.2)+
  scale_color_viridis_d(
    option = 'plasma'
  )+
  theme(legend.position = 'none',
        axis.text = element_text(size = 14),
      axis.title = element_text(size = 18))+
  scale_x_continuous(name = 'edge weight', limits = c(0,1))

plot_chains <- rbind(plot_high, plot_mid, plot_low)                             # make data set for trace plots
plot_chains$order <- factor(plot_chains$dyad_id,
                            levels = c('73914','86520','85862'))                # reorder to colour by chain so matches plot above, and also has narrowest chain on top and most uncertain at the back
ggplot(#data = plot_chains,
  #aes(y = draw, x = chain_position, colour = order)
)+
  geom_line(data = plot_high, aes(y = draw, x = chain_position),                # chain for most uncertain
            #colour = '#fde725',                                                # yellow --  viridis scale B
            colour = '#f0f921'                                                  # yellow -- viridis scale 'plasma'
  )+
  geom_line(data = plot_mid, aes(y = draw, x = chain_position),                 # chain for median
            #colour = '#21918c',                                                # turquoise -- viridis scale B
            colour = '#e16462'                                                  # orange -- viridis scale 'plasma'
  )+
  geom_line(data = plot_low, aes(y = draw, x = chain_position),                 # chain for minimum edge weight
            #colour = '#440154',                                                # purple -- viridis scale B
            colour = '#0d0887'                                                  # red -- viridis scale 'plasma'
  )+
  #scale_color_viridis_d()+
  theme_classic()+
  theme(#legend.position = 'none',
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18))+
  scale_x_continuous(name = 'chain position')+
  scale_y_continuous(name = 'chain position',
                     limits = c(0,1))

### add time marker
print(paste0('binomial conjugate model completed at ', Sys.time()))

## non-random edge weights ####
## run null model -- RUN THIS SEPARATELY AS CURRENTLY SCREWING THINGS UP
#load('motnp_bisonr_edgescalculated_widebinomconj.RData')
#priors <- get_default_priors('binary_conjugate') # obtain structure for bison model priors
#priors$edge <- 'beta(0.7,1.5)'
#motnp_edges_null <- bison_model(
#  (event | duration) ~ 1,                       # response variable does not vary with dyad, all dyads drawn from the same distribution -- THIS RUNS BUT COMPARISONS TO FULL MODEL DON'T WORK
#  data = counts_df_model, 
#  model_type = "binary_conjugate",
#  priors = priors
#)

### compare null model with fitted model -- bisonR model stacking
#model_comparison(list(non_random_model = motnp_fit_edges, random_model = motnp_edges_null)) # compare fit for model allowed to vary by dyad vs model that draws all dyad strengths from the same distribution -- vast majority of network is best explained by random model, but 4.3% of dyads better explained by non-random. Network is therefore non-random, but with only a small proportion of dyads associating more strongly than random distribution will allow

### compare null model with fitted model -- Jordan pseudo model averaging
#model_averaging <- function(models) {
#  loos <- lapply(models, function(model) loo::loo(model$log_lik, r_eff=NA))
#  
#  results_matrix <- loo::loo_model_weights(loos, method="pseudobma")
#  if (!is.null(names(models))) {
#    names(results_matrix) <- names(models)
#  }
#  results_matrix
#}  # produce alternative method for comparing models
#model_averaging(models = list(non_random_model = motnp_fit_edges, random_model = motnp_edges_null))            # 100% confidence that random model is better

# save workspace image
#save.image('motnp_bisonr_edgescalculated_widebinomconj_v2.RData')

### add time marker
#print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
### adapt bisonR plot_network function to give more flexibility over plotting options
plot_network_threshold <- function (obj, ci = 0.9, lwd = 2, threshold = 0.3,  # threshold = median edge weight must be over X
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
                      vertex.label.color = vertex.label.color1, vertex.color = vertex.color1,
                      vertex.size = vertex.size1,
                      edge.color = rgb(0, 0, 0, 1), edge.arrow.size = 0)
  igraph::plot.igraph(net, layout = coords,
                      edge.width = ifelse(md < threshold, 0, ub * lwd),
                      vertex.label.color = vertex.label.color2, vertex.color = vertex.color2,
                      vertex.size = vertex.size2,
                      edge.color = edge.color1, edge.arrow.size = 0, 
                      add = TRUE)
  if (obj$directed) {
    igraph::plot.igraph(net, edge.width = md * 0, layout = coords, 
                        vertex.label.color = vertex.label.color2, vertex.color = vertex.color2, 
                        edge.color = edge.color2, add = TRUE)
  }
}

### load in ages 
motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

### create nodes data frame
nodes <- data.frame(id = sort(unique(c(counts_df$id_1,counts_df$id_2))),  # all unique individuals
                    node = NA, age = NA, sightings = NA)                  # data needed on each
for(i in 1:nrow(nodes)){
  # add age data
  x <- motnp_ages[motnp_ages$id == nodes$id[i],]                          # vector of ages per individual
  nodes$age[i] <- mean(x$age)                                             # for initial test purposes, age = mean only
  # add node id and count data
  if(nodes$id[i] != 'M99') { y <- counts_df[counts_df$id_1 == nodes$id[i], c('id_1','count_1','node_1_males')] }   # add data for M99 -- only male who is only present in counts_df in ID2, all others are ID1
  else { y <- counts_df[counts_df$id_2 == nodes$id[i], c('id_2','count_2','node_2_males')] }  # add data for all other males
  nodes$sightings[i] <- y[1,2]
  nodes$node[i] <- y[1,3]
}
nodes$sightings <- as.numeric(nodes$sightings)   # convert to numeric
nodes$node <- as.character(nodes$node)           # convert to character
str(nodes)                                       # check structure

### plot network
plot_network_threshold(motnp_fit_edges, lwd = 15, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')

### adapt to remove unconnected nodes and simplify input requirements
plot_network_threshold2 <- function (obj, ci = 0.95, lwd = 2, threshold = 0.3,
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

### plot network
plot_network_threshold2(obj = motnp_fit_edges, threshold = 0.15,
                        node.size = nodes, node.colour = nodes, lwd = 10)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## compare edge weight distributions to simple SRI ####
# load in workspace image with edge weights already calculated
#load('motnp_bisonr_edgescalculated_widebinomconj.RData')

### plot against SRI
edgelist <- get_edgelist(motnp_fit_edges)
head(edgelist)                                                                      # check structure of edgelist
colnames(edgelist)[1:2] <- colnames(counts_df_model)[1:2]                           # rename for join function
edgelist$node_1_id <- as.integer(edgelist$node_1_id)                                # convert to integers
edgelist$node_2_id <- as.integer(edgelist$node_2_id)                                # convert to integers
summary <- left_join(edgelist, counts_df_model, by = c('node_1_id','node_2_id'))    # combine distribution data with raw counts
counts <- counts_df[,c('node_1_males','node_2_males','count_1','count_2')]          # add individual count data
colnames(counts)[1:2] <- c('node_1_id','node_2_id')                                 # rename for join function
summary <- left_join(summary, counts, by = c('node_1_id','node_2_id'))              # combine summary data with individual counts
summary$sri <- summary$event / (summary$duration)                                   # calculate basic SRI

plot(density(summary$sri), main = 'SRI vs model output: blue=all,\nred=both seen 8 times, green=both 12 times') # plot SRI distribution
lines(density(summary$median), col = 'blue')                                                                    # distribution of median estimates
lines(density(summary$median[which(summary$count_1 >= 8 & summary$count_2 >= 8)]), col = 'red')                 # distribution of median estimates for frequently sighted males
lines(density(summary$median[which(summary$count_1 >= 12 & summary$count_2 >= 12)]), col = 'green')             # distribution of median estimates for frequently sighted males

## try plotting a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
dyads <- counts_df[,c('dyad_id','node_1_males','node_2_males')]                        # create data frame of dyads
colnames(dyads) <- c('dyad_id','node_1_id','node_2_id')                                # rename for joining
dyads <- left_join(dyads, counts_df_model, by = c('node_1_id','node_2_id'))            # join with counts data (creates counts_df_model but with dyad_id)
length(which(is.na(dyads$duration) == TRUE))

draws <- as.data.frame(motnp_fit_edges$chain) %>%                      # generate dataframe of draws per dyad
  pivot_longer(cols = everything(), names_to = 'dyad_model', values_to = 'edge')       # convert to a long format
draws$dyad_id <- rep(counts_df$dyad_id, 4000)                                          # add dyad_id data
draws$weight <- gtools::inv.logit(draws$edge)                                          # convert to 0-1 bounded values
draws$draw <- rep(1:4000,  each = nrow(counts_df_model))                               # add chain position number
draws <- left_join(draws, dyads, by = 'dyad_id')                                       # combine to add node data
draws$sri <- draws$event / draws$duration                                              # calculate basic SRI to compare with distributions

set.seed(15)                                                                           # make subsetting reproducible
subset_draws <- draws[draws$dyad_id %in% sample(draws$dyad_id, 150, replace = F),]     # sample 150 dyads at random to plot (takes a very long time)
subset_draws$median <- NA ; for(i in 1:nrow(subset_draws)){                            # set up for loop
  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]                  # generate dataset cut down to just the dyad in question
  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), subset_draws$median[i]) }   # record median draw value
head(subset_draws)                                                                     # check structure of data frame
which(is.na(subset_draws$median) == TRUE)[1]                                           # check that all are filled in

subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)           # order dataframe based on total sightings per pair
ggplot(data = subset_draws, mapping = aes(x = weight))+                                # plot dataframe
  geom_density(colour = 'blue')+                                                       # plot density plots of probability distirbtuions output per dyad
  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+                    # split by dyad, allow to vary height depending on dyad
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+            # add line showing where the median estimate is
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')                          # add line showing where the SRI value is

write_csv(subset_draws, '../data_processed/motnp_sampledyads_random_binary_vs_sri_widebinomconj_v2.csv')    # save output for future reference

subset_draws <- read_csv('../data_processed/motnp_sampledyads_random_binary_vs_sri_widebinomconj.csv')
subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)           # order dataframe based on total sightings per pair
ggplot(data = subset_draws, mapping = aes(x = weight))+                                # plot dataframe
  geom_density(colour = 'blue')+                                                       # plot density plots of probability distirbtuions output per dyad
  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+                    # split by dyad, allow to vary height depending on dyad
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+            # add line showing where the median estimate is
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')                          # add line showing where the SRI value is

subset_draws <- draws[draws$sri > 0.2,]                                                # repeat subsetting but selecting only the dyads with the highest SRI values (generally some of the least sighted individuals)
subset_draws$median <- NA ; for(i in 1:nrow(subset_draws)){                            # set up for loop
  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]                  # generate dataset cut down to just the dyad in question
  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), subset_draws$median[i]) }   # record median draw value
head(subset_draws)                                                                     # check structure of data frame
which(is.na(subset_draws$median) == TRUE)[1]                                           # check that all are filled in

subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)           # order dataframe based on total sightings per pair
ggplot(data = subset_draws, mapping = aes(x = weight))+                                # plot dataframe
  geom_density(colour = 'blue')+                                                       # plot density plots of probability distributions output per dyad
  facet_wrap(. ~ dyad_id, ncol = 15)+                                                  # split by dyad, allow to vary height depending on dyad
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+            # add line showing where the median estimate is
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')                          # add line showing where the SRI value is

ggplot(data = subset_draws, mapping = aes(x = weight))+                                # plot dataframe
  geom_density(colour = 'blue')+                                                       # plot density plots of probability distirbtuions output per dyad
  facet_wrap(. ~ dyad_id, ncol = 15, scales = 'free_y')+                    # split by dyad, allow to vary height depending on dyad
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+            # add line showing where the median estimate is
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')                          # add line showing where the SRI value is

write_csv(subset_draws, '../data_processed/motnp_sampledyads_sri0.2_binary_vs_sri_widebinomconj_v2.csv')    # save output for future reference

# clean environment
rm(draws, dyads, priors, subset_draws, x) ; gc()
dev.off()

## coefficient of variation of edge weights (aka social differentiation) ####
#load('motnp_bisonr_edgescalculated_widebinomconj.RData')

# add pdf output file
pdf(file = '../outputs/motnp_bisonr_edgeweight_widebinomconj_cv.pdf')

# extract cv for model
global_cv_widebinomconj <- extract_metric(motnp_fit_edges, "global_cv", num_draws = 10000)
head(global_cv_widebinomconj)
hist(global_cv_widebinomconj)

# calculate SRI for all dyads
counts_df$sri <- counts_df$event_count / counts_df$count_dyad
summary(counts_df$sri)

# calculate CV of SRI for all dyads
raster::cv(counts_df$sri)         # very high, but lower than gbi_matrix
#raster::cv(m$sri[m$sri > 0])     # still massive even when I remove the 0s -- zero inflation is real

### create SRI matrix
# generate matrix
N <- length(unique(c(counts_df$id_1, counts_df$id_2)))      # number of elephants
ids <- unique(c(counts_df$id_1, counts_df$id_2))            # IDs of elephants
m_mat <- diag(nrow = N)                                     # matrix of males -- NxN to fill with SRI
colnames(m_mat) <- ids ; rownames(m_mat) <- ids             # dimnames = elephant IDs

# populate matrix with SRI values
for( i in 1:N ) {                                           # rows
  for( j in 1:N ) {                                         # columns
    if( i >= j ) { m_mat[i,j] <- m_mat[j,i] }               # make symmetrical about diagonal
    else {
      id1 <- colnames(m_mat)[i]                             # identify elephant 1 of dyad
      id2 <- rownames(m_mat)[j]                             # identify elephant 2 of dyad
      m_mat[i,j] <- counts_df$sri[which(counts_df$id_1 == id1 & counts_df$id_2 == id2)]      # add SRI values -- match IDs and dyad
    }
  }
}

# calculate CV for matrix
raster::cv(m_mat)                # still way too high
which(is.nan(m_mat) == TRUE)     # check if any are blank
sd(m_mat)/mean(m_mat)            # matches raster version -- check it is doing what I think it is doing

# create gbi_matrix
eles <- read_csv(file = '../data_processed/motnp_eles_long.csv')       # read in long format of raw data
eles$location <- paste(eles$gps_s, eles$gps_e, sep = '_')              # make single variable for unique locations
eles <- eles[,c(1,16,2,3,17,4,5,14,7,8,10,13)]                         # rearrange variables
nodes_data <- read_csv(file = '../data_processed/motnp_elenodes.csv' ) # read in node data
colnames(nodes_data)
str(nodes_data)
eles_asnipe <- eles[,c(3,4,2,5)]                                       # date, time, elephant, location
eles_asnipe$Date <- as.integer(eles_asnipe$date)                       # make date numeric
eles_asnipe$Date <- 1+eles_asnipe$Date - min(eles_asnipe$Date)         # start from 1, not 1st January 1970
eles_asnipe$Time <- ifelse(eles_asnipe$time > 1, NA, eles_asnipe$time) # time = proportion of day so anything >1 has to be wrong
eles_asnipe$Time <- eles_asnipe$Time*(24*60*60)                        # convert time values to seconds through day
which(is.na(eles_asnipe$Time))                                         # 161 698 1122 1469 1770
eles_asnipe[c(161,698,1122,1469,1770),]                                # all 1 sighting of B7+M44
eles_asnipe <- eles_asnipe[,c(5,6,3,4)]                                # create data frame to produce gbi matrix from
colnames(eles_asnipe) <- c('Date','Time','ID','Location')              # rename variables for get_gbi
eles_asnipe$ID <- as.character(eles_asnipe$ID)                         # correct data type
str(eles_asnipe)
eles_asnipe$d_pad <- str_pad(eles_asnipe$Date, 3, pad = '0')  # 0-pad dates
eles_asnipe$encounter <- paste(eles_asnipe$d_pad, eles_asnipe$Time, eles_asnipe$Location, sep = '_') # unique value for each sighting
eles_asnipe$group <- as.integer(as.factor(eles_asnipe$encounter)) # unique factor for every sighting
max(eles_asnipe$group)                         # 574 -- agrees with number of different sightings for which elephants were identified
eles_dt <- eles_asnipe[,c(3,7)]                # create data table for gbi matrix
eles_dt <- data.table::setDT(eles_dt)          # create data table for gbi matrix
gbi_matrix <- spatsoc::get_gbi(DT = eles_dt, group = 'group', id = 'ID')  # create gbi matrix
gbi_males <- gbi_matrix[,colnames(gbi_matrix) %in% ids]
gbi_males <- gbi_males[rowSums(gbi_males) > 0,] # remove male-only sightings

# set up permutations
N_networks <- 10000
# rm(counts_df_model, edgelist, eles_asnipe, eles_dt, females_df, gbi_matrix, motnp_ages, motnp_fit_edges, motnp_edges_null, i, id1, id2, j, model_averaging) ; gc()

# create vector of days for each sighting
sightings <- eles[,c('elephant','encounter','date')] %>%  # set up dataframe of actual encounters
  filter(elephant %in% ids) %>%                           # remove encounters that only include females, youngsters, or dead males
  dplyr::select(-elephant) %>%                                   # remove ID column
  distinct()                                              # cut down to one row per encounter

# run network permutations
random_networks <- asnipe::network_permutation(association_data = gbi_males, # permute network -- raw data = gbi_matrix
                                               association_matrix = m_mat,   # SRI matrix
                                               permutations = N_networks,    # 10000 times
                                               days = sightings$date,        # permute within day only
                                               within_day = TRUE)            # permute within day only
print('random networks generated')

# calculate coefficient of variation per randomised network
cv_random_networks <- rep(0,N_networks) ; for (i in c(1:N_networks)) {       # set up for loop
  net_rand <- random_networks[i,,]                                           # select random network
  cv_random_networks[i] <- raster::cv(net_rand)                              # calculate CV and save into vector
  if(i %% 1000 == 0) {print(i)}                                              # track progress
}
print(paste0('network permutations for entire network completed at ', Sys.time()))

# compare permuted networks to actual network
hist(cv_random_networks, xlim = c(min(cv_random_networks)-50,max(cv_random_networks)+50),   # plot histogram of random network CVs
     #ylim = c(0,10),  # there are 3 networks out of 10000 that have cv > 290, but then all are in the same bar as cv(m_mat)
     main = 'permutations for network of all male elephants')
abline(v = raster::cv(m_mat), col = 'red')                                                  # add line for measured network
text(round(raster::cv(m_mat),3), col = 'red', x = median(cv_random_networks)+50, y = 2000)  # add text to specify value for measured network

plot_cv <- as.data.frame(cv_random_networks)                # create dataframe
cv_network <- raster::cv(m_mat)                             # save value of measured CV
cv_random_networks <- plot_cv$cv_random_networks            # vector of randomised CV values
ggplot(data = plot_cv)+                                     # plot CV values
  geom_vline(xintercept = cv_network, linewidth = 1.5,      # add line for true value
             colour = rgb(68/255, 1/255, 84/255))+          # dark purple (viridis scale)
  #geom_text(x = cv_network - 50, y = 700,
  geom_histogram(aes(cv_random_networks),                   # random cv values
                 fill = rgb(253/255, 231/255, 37/255),      # yellow (viridis scale)
                 colour = 'black', bins = 100)+             # set bin width
  theme_classic()+
  scale_x_continuous(name = 'coefficient of variation',
                     expand = c(0,0),
                     limits = c(min(cv_random_networks)-10,
                                cv_network+10))+
  scale_y_continuous(name = 'frequency',
                     expand = c(0,0))+
  #          colour = rgb(68/255, 1/255, 84/255),
  #          label = 'coefficient of/nvariation for/nmeasured network')+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))

### write out outputs for future reference
write_csv(plot_cv, '../data_processed/motnp_networkpermutations_cv_widebinomconj_v2.csv')

### combine with global_cv
plot_cv_model <- data.frame(cv = c(plot_cv$cv_random_networks, (global_cv_widebinomconj*100)),             # combine permutations and model values
                            iteration = rep(1:length(global_cv_widebinomconj), 2),
                            type = rep(c('permutation','model_draw'), each = length(global_cv_widebinomconj)))
write_csv(plot_cv_model, '../data_processed/motnp_networkpermutations_cv_widebinomconj_v2.csv')               # save for future reference
ggplot()+
  geom_vline(xintercept = cv_network, linewidth = 1.5,                                                   # line for measured network from SRI values
             colour = rgb(68/255, 1/255, 84/255))+                                                       # dark purple (viridis)
  #annotate('text', x = cv_network - 50, y = 700), colour = rgb(68/255, 1/255, 84/255),
  #         label = 'coefficient of/nvariation for/nmeasured network')+
  theme_classic()+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))+
  geom_histogram(data = plot_cv_model[plot_cv_model$type == 'model_draw',], aes(x = cv),                 # histogram of CV draws from model
                 fill = rgb(33/255, 145/255, 140/255),                                                   # yellow (viridis)
                 bins = 100, colour = 'black')+
  geom_histogram(data = plot_cv_model[plot_cv_model$type == 'permutation',], aes(x = cv),                # histogram of CV draws from permutations
                 fill = rgb(253/255, 231/255, 37/255),                                                   # turquoise  (viridis)
                 bins = 100, colour = 'black')+
  scale_x_continuous(name = 'coefficient of variation',
                     #limits = c(min(cv_random_networks)-10,
                     #cv_network+10),
                     expand = c(0,0))+
  scale_y_continuous(name = 'frequency',
                     expand = c(0,0))

# repeat plot for poster -- larger fonts, fewer draws
plot_cv_model2 <- plot_cv_model[c(1:10000, sample(10001:20000, 4000, replace = F)),]                     # select only 4000 values from dataframe
ggplot()+
  geom_vline(xintercept = cv_network, linewidth = 1.5,                                                   # line for measured network from SRI values
             colour = rgb(68/255, 1/255, 84/255))+                                                       # dark purple (viridis)
  #annotate('text', x = cv_network - 50, y = 700), colour = rgb(68/255, 1/255, 84/255),
  #         label = 'coefficient of/nvariation for/nmeasured network')+
  theme_classic()+
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 22))+
  geom_histogram(data = plot_cv_model2[plot_cv_model2$type == 'model_draw',], aes(x = cv),               # histogram of CV draws from model
                 fill = rgb(33/255, 145/255, 140/255),                                                   # yellow (viridis)
                 bins = 100, colour = 'black')+
  geom_histogram(data = plot_cv_model2[plot_cv_model2$type == 'permutation',], aes(x = cv),              # histogram of CV draws from permutations
                 fill = rgb(253/255, 231/255, 37/255),                                                   # turquoise  (viridis)
                 bins = 100, colour = 'black')+
  scale_x_continuous(name = 'coefficient of variation',
                     #limits = c(min(cv_random_networks)-10,
                     #cv_network+10),
                     breaks = round(seq(100, 400, by = 50),-1),
                     expand = c(0,0))+
  scale_y_continuous(name = 'frequency',
                     expand = c(0,0))

### run statistical test to confirm that difference is significant
# plot_cv <- read_csv('../data_processed/motnp_networkpermutations_cv_widebinomconj.csv')
hist(plot_cv$cv_random_networks)    # plot cv values from random networks -- 1 tiny bar very high value and outside of the main distribution
cv_network                          # plot line for measured value -- falls over the tiny bar

mean(plot_cv$cv_random_networks)    # mean is far lower than measured value
sd(plot_cv$cv_random_networks)      # mean + 2*SD is still far lower than measured value

length(which(plot_cv$cv_random_networks >= cv_network))/length(plot_cv$cv_random_networks) # 3 / 10000 = 3e-04 = p-value (don't need an additional test because that will further randomise the data which are already randomised, just need to know proportion that are greater than or equal to your value)

mean(global_cv_widebinomconj) ; sd(global_cv_widebinomconj)     # totally outside either the random distribution and also doesn't overlap at all with the SRI measured network
max(global_cv_widebinomconj) ; min(plot_cv$cv_random_networks)     # no overlap

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2% ####
counts_df_model <- counts_df[,c('node_1_males','node_2_males','event_count','count_dyad','id_1','id_2','dyad_males')]   # select variables of interest
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')                      # rename variables
ms_chain <- as.data.frame(motnp_fit_edges$chain)                                                        # extract chain values
colnames(ms_chain) <- counts_df_model$dyad_id                                                                           # rename with dyad IDs
ms_chain <- pivot_longer(ms_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')                       # pivot longer
ms_chain$chain_position <- rep(1:length(unique(ms_chain$dyad_id)), each = 4000)                                         # add chain positions
ms_chain$draw <- LaplacesDemon::invlogit(ms_chain$draw)                                                                 # convert values to 0-1 bounded
ms_chain$mean <- NA                                                                                                     # set up for loop for mean values
hist(ms_chain$draw)                                                                                                     # plot histogram of draw values
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using binomial conjugate prior',           # empty plot for density lines
     xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sample(unique(ms_chain$dyad_id), size = 1000, replace = F)){                                                   # set up for loop using 1000 sample dyads
  x <- ms_chain[ms_chain$dyad_id == i,]                                                                                 # select values for dyad weights
  ms_chain$mean <- ifelse(ms_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ms_chain$mean)                             # calculate mean values
  lines(density(x$draw), col = rgb(0,0,1,0.1))                                                                          # plot sample line
}
lines(density(ms_chain$draw), lwd = 2)                                                                                  # plot average line

quantile(ms_chain$draw, 0.98)                                                                                           # identify 98th percentile draw
quantile(edgelist$median, 0.98)                                                                                         # identify 98th percentile median
counts_df_model$sri <- counts_df_model$event / counts_df_model$duration                                                 # calculate SRI values
quantile(counts_df_model$sri, 0.98)                                                                                     # identify 98th percentile SRI

## clean up ####
### save pdf
dev.off()

### add time marker
print(paste0('binomial conjugate model completed at ', Sys.time()))

rm(list= ls()[!(ls() %in% c('motnp_fit_edges','motnp_edges_null',
                            'counts_df','counts_df_model',
                            'gbi_males','m_mat',
                            #'random_nodes',
                            'random_networks'))])

### write over saved workspace image -- only keep the things you might need later that take a long time to run
save.image('motnp_bisonr_edgescalculated_widebinomconj_v2.RData')
