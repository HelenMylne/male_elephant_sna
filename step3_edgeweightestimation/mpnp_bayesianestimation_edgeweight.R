#### Bayesian analysis of ATE data ####
# Script to process association data from Makgadikgadi Pans National Park, Botswana.
# Data collected: 1st November 2019 to 5th August 2021
# Collected by: Elephants for Africa (Dr Kate Evans)
# Data supplied by: EfA and Dr Kate Evans, 19th October 2021
# Data input: raw data provided by EfA processed in scripts mpnp_dataprocessing1.R and mpnp_dataprocessing2.R

#### set up ####
### install packages
#install.packages('tidyverse', lib = '../packages/')
#install.packages('dplyr', lib = '../packages/')
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages(c("StanHeaders","rstan"),type="source")
#install.packages("remotes", lib = '../packages/')
#remotes::install_github("stan-dev/cmdstanr", lib = '../packages/') ## OR USE # install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages('rethinking')
#install.packages('igraph', lib = '../packages/')
#install.packages('dagitty', lib = '../packages/')
#install.packages('janitor', lib = '../packages/')
#install.packages('lubridate', lib = '../packages/')
#install.packages('hms', lib = '../packages/')
#install.packages('readxl', lib = '../packages/')
#install.packages('LaplacesDemon', lib = '../packages/')

### load packages
# library(tidyverse) ; library(dplyr) ; library(cmdstanr) ; library(igraph) ; library(janitor) ; library(lubridate) ; library(hms) ; library(readxl)
library(cmdstanr, lib.loc = '../packages/')    # library(cmdstanr)
library(tidyverse, lib.loc = '../packages/')   # library(tidyverse)
library(dplyr, lib.loc = '../packages/')       # library(dplyr)
#library(rstan, lib.loc = '../packages/')      # library(rstan)
library(igraph, lib.loc = '../packages/')      # library(igraph)
library(janitor, lib.loc = '../packages/')     # library(janitor)
library(lubridate, lib.loc = '../packages/')   # library(lubridate)
library(readxl, lib.loc = '../packages/')      # library(readxl)

### set stan path
#set_cmdstan_path('/Users/helen/.cmdstan/cmdstan-2.32.1')
set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0')

### set seed
set.seed(12345)

### create file of output graphs
pdf('../outputs/mpnplong_edgeweights_conditionalprior.pdf', width = 20, height = 15)

### set ggplot theme
theme_set(theme_light())

### load edge weight model
edge_binary <- cmdstan_model("models/edge_binary_basic.stan")   # load model
edge_binary                                                     # check model priors etc.

### add progress marker
print(paste0('model loaded in at ', Sys.time()))

#### import data ####
### import aggregated counts of sightings together and apart, and info about individuals
counts_df <- read_csv('../data_processed/step1_dataprocessing/mpnp_longtimewindow_pairwiseevents.csv') %>%   # counts_df = data frame of all dyad information, specifically including aggregated counts of together vs apart
  distinct()
nrow(counts_df) - length(unique(counts_df$dyad_id)) # 0

### create nodes data frame for period 1
nodes <- data.frame(id = sort(unique(c(counts_df$id_1,counts_df$id_2))),    # all unique individuals
                    node = NA, sightings = NA)                              # data needed on each
for(i in 1:nrow(nodes)){
  # extract data about individual from counts_df data frame
  if(nodes$id[i] %in% counts_df$id_1) {
    x <- counts_df[counts_df$id_1 == nodes$id[i], c('id_1','node_1','count_1')] %>% distinct()
  } else { x <- counts_df[counts_df$id_2 == nodes$id[i], c('id_2','node_2','count_2')] %>% distinct() }
  colnames(x) <- c('id','node','count')
  # add individual data
  nodes$node[i] <- x$node          # node ID number
  nodes$sightings[i] <- x$count    # number of sightings within the time window
}

### import age data
eles <- readRDS('../data_processed/step2_ageestimation/mpnp_longwindow_ageestimates_mcmcoutput.rds')
mean_ages <- data.frame(id = colnames(eles),
                        age = apply(eles, 2, mean))

# combine node counts with age data
nodes <- left_join(nodes, mean_ages, by = 'id') # join age data to sightings
missing_age <- nodes[is.na(nodes$age),]         # 52 eles don't have ages -- these are the 52 individuals for whom no age estimate was available in the raw data (category 10 = unknown)

### create data list
n_chains <- 4                            # number of MCMC chains to run
n_samples <- 1000                        # number of samples per chain
n_dyads <- nrow(counts_df)               # number of dyads in time window
counts_ls <- list(                       # create data list
  n_dyads    = n_dyads,                  # total number of times one or other of the dyad was observed
  dyad_ids   = counts_df$dyad_id,        # identifier for each dyad
  together   = counts_df$together,       # count number of sightings seen together
  count_dyad = counts_df$count_dyad      # count total number of times seen
)

### add progress marker
print(paste0('data list created at ', Sys.time()))

#### run model on real standardised data ####
### Fit model
fit_edges_mpnp1 <- edge_binary$sample(
  data = counts_ls, 
  chains = n_chains, parallel_chains = n_chains,
  iter_warmup = n_samples, iter_sampling = n_samples)

### check model
fit_edges_mpnp1

### add progress marker
print(paste0('model run completed at ', Sys.time()))

### Extract posterior samples
posterior_samples <- fit_edges_mpnp1$draws()

### extract edge weights
edge_weights_matrix <- posterior_samples[,,2:(nrow(counts_df)+1)]  # remove column of intercepts

### clean up and save workspace
rm(posterior_samples, x, i) ; gc()                                 # clean workspace
saveRDS(edge_weights_matrix, '../data_processed/step3_edgeweightestimation/mpnplong_edgedistributions_conditionalprior_wideformat.RDS')
save.image('mpnp_edgecalculations/mpnplong_edgeweights_conditionalprior.RData')
#load('mpnp_edgecalculations/mpnplong_edgeweights_conditionalprior.RData')
#edge_weights_matrix <- readRDS('../data_processed/step3_edgeweightestimation/mpnplong_edgedistributions_conditionalprior_wideformat.RDS')

### convert edge samples to data frame, putting all chains into 1 column
#edges <- as.data.frame(edge_weights_matrix[,,1])               # convert draws for first dyad to data frame
#colnames(edges) <- c('chain1','chain2','chain3','chain4')      # rename columns
#edges <- pivot_longer(edges, everything(), values_to = 'edge_draw', names_to = 'chain') # tidy draws for first dyad
#edges$dyad <- counts_ls$dyad_ids[1]                            # add dyad id column
#edges$position <- rep(1:n_samples, each = n_chains)            # position within chain for traceplots
#for(i in 2:n_dyads){                                           # repeat for the rest of the dyads
#  x <- as.data.frame(edge_weights_matrix[,,i])
#  colnames(x) <- c('chain1','chain2','chain3','chain4')
#  x <- pivot_longer(x, everything(), values_to = 'edge_draw', names_to = 'chain')
#  x$dyad <- counts_ls$dyad_ids[i]
#  x$position <- rep(1:n_samples, each = n_chains)
#  edges <- rbind(edges, x)                                     # append each to the last to create a single dataframe with all dyads
#}

### save edge samples for use in dyadic/nodal regressions
#saveRDS(edges, '../data_processed/step3_edgeweightestimation/mpnplong_edgedistributions_conditionalprior.RDS')
#edges <- readRDS('../data_processed/step3_edgeweightestimation/mpnplong_edgedistributions_conditionalprior.RDS')

### add progress marker
print(paste0('edge values extracted at ', Sys.time()))

#### check outputs: edge weights ####
### Assign random set of columns to check
if(length(which(counts_df$together >= 1)) >= 200){ n_test <- 200 } else { n_test <- length(which(counts_df$together >= 1)) } # set number of dyads to sample: will be equal numbers of dyads that have and have not ever been seen together, so if there are 200 or more dyads that have been seen together at some point then cap at 200, otherwise use as many as there are
plot_dyads <- c(sample(counts_df$dyad_id[counts_df$together >= 1], size = n_test, replace = F), # randomly sample dyads that have been seen together
                sample(counts_df$dyad_id[counts_df$together == 0], size = n_test, replace = F)) # randomly sample dyads never seen together
#plot_edges <- edges[edges$dyad %in% plot_dyads,]                    # subset edge weight data using randomly drawn dyads
#plot_edges$seen_together <- NA ; for(i in 1:length(plot_dyads)){    # set up for loop
#  plot_edges$seen_together[plot_edges$dyad == plot_dyads[i]] <- ifelse(counts_df$together[counts_df$dyad_id == plot_dyads[i]] > 0, 1, 0) # the value of seen_together is 1 if the dyad has ever been seen in the same group, and 0 if they have not
#}

plot_edges <- as.data.frame(edge_weights_matrix[,,which(counts_df$dyad_id == plot_dyads[1])])        # convert draws for first dyad to data frame
colnames(plot_edges) <- c('chain1','chain2','chain3','chain4')                                       # rename columns
plot_edges <- pivot_longer(data = plot_edges, cols = everything(), values_to = 'edge_draw', names_to = 'chain')    # tidy draws for first dyad
plot_edges$dyad <- plot_dyads[1]                                                                     # add dyad id column
plot_edges$position <- rep(1:n_samples, each = n_chains)                                             # position within chain for traceplots
for(i in 2:length(plot_dyads)){                                                                      # repeat for the rest of the dyads
  x <- as.data.frame(edge_weights_matrix[,,which(counts_df$dyad_id == plot_dyads[i])])
  colnames(x) <- c('chain1','chain2','chain3','chain4')
  x <- pivot_longer(x, everything(), values_to = 'edge_draw', names_to = 'chain')
  x$dyad <- plot_dyads[i]
  x$position <- rep(1:n_samples, each = n_chains)
  plot_edges <- rbind(plot_edges, x)                                                                 # append each to the last to create a single dataframe with all dyads
}
plot_edges <- plot_edges %>% 
  rename(dyad_id = dyad) %>% 
  left_join(counts_df[,c('dyad_id','count_dyad','together','apart')]) %>% 
  mutate(seen_together = ifelse(together > 0, 1, 0))

### build traceplots
ggplot(data = plot_edges[plot_edges$seen_together == 1,],   # plot only dyads that have been seen together
       aes(y = edge_draw, x = position, colour = chain))+   # plot all chains over each other for each dyad to check mixing
  geom_line()+                                              # draw as line plot
  facet_wrap(. ~ dyad_id)+                                     # new panel per dyad as each pair has own edge weight parameter
  theme_classic()+                                          # make it look nicer
  theme(legend.position = 'none',                           # remove legend
        strip.background = element_blank(), strip.text = element_blank())    # remove facet strips
ggplot(data = plot_edges[plot_edges$seen_together == 0,],   # repeat for dyads that have never been seen together
       aes(y = edge_draw, x = position, colour = chain))+
  geom_line()+
  facet_wrap(. ~ dyad_id)+
  theme_classic() + theme(legend.position = 'none', strip.background = element_blank(), strip.text = element_blank())

### density plots
plot(NULL, xlim = c(0,1), ylim = c(0,30), las = 1, xlab = 'edge weight', ylab = 'density')         # set up plot window
for(i in 1:length(plot_dyads)){                                                                    # plot randomly sampled dyads
  x <- plot_edges[plot_edges$dyad_id == plot_dyads[i],]                                               # select data to plot
  lines(density(x$edge_draw), col = ifelse(x$seen_together == 1, rgb(0,0,1,0.1), rgb(1,0,0,0.1)))  # draw edge weight probability plot. blue = seen together at least once, red = never seen together
}

### add progress marker
print(paste0('edges checked at ', Sys.time()))

#### check outputs: plot network ####
### create custom network plotting function
plot_network_threshold_mpnp <- function (edge_samples, dyad_data, lwd = 2, threshold = 0.3,
                                     label.colour = 'transparent', label.font = 'Helvetica', 
                                     node.size = 4, node.colour = 'seagreen1',
                                     link.colour1 = 'black', link.colour2 = rgb(0, 0, 0, 0.3))
{
  dyad_name <- do.call(paste, c(dyad_data[c("node_1", "node_2")], sep=" <-> "))
  edge_lower <- apply(edge_samples, 2, function(x) quantile(x, probs=0.025))
  edge_upper <- apply(edge_samples, 2, function(x) quantile(x, probs=0.975))
  edge_median <- apply(edge_samples, 2, function(x) quantile(x, probs=0.5))
  edge_list <- cbind(
    "median"=round(edge_median, 3), 
    "2.5%"=round(edge_lower, 3), 
    "97.5%"=round(edge_upper, 3)
  )
  rownames(edge_list) <- dyad_name
  edgelist <- as.data.frame(edge_list)
  edgelist$node_1 <- as.character(dyad_data$node_1)
  edgelist$node_2 <- as.character(dyad_data$node_2)
  edgelist <- edgelist[,c(4:5,1:3)]
  #net_all <- igraph::graph_from_edgelist(as.matrix(edgelist[, 1:2]), directed = F)
  threshold_edges <- edgelist[edgelist$median >= threshold,]
  if(nrow(threshold_edges) == 0) { stop('No edges above threshold') }
  net <- igraph::graph_from_edgelist(as.matrix(threshold_edges[, 1:2]), directed = F)
  
  if(is.data.frame(node.size) == TRUE ) {
    nodes_list <- data.frame(node = rep(NA, length(unique(c(threshold_edges$node_1, threshold_edges$node_2)))), #as.numeric(names(net_all[[1]])),
                             sightings = NA)
    for(i in 1:nrow(nodes_list)){
    nodes_all <- rep(NA, 2*nrow(threshold_edges))  
      for(a in 1:2){
        for(b in 1:nrow(threshold_edges)){
          nodes_all[a + (b-1)*2] <- threshold_edges[b,a]
        }
      }
      nodes_list$node <- unique(nodes_all)
      nodes_list$sightings[i] <- nodes$sightings[which(nodes$node == nodes_list$node[i])]
   }
    node_sightings <- nodes_list$sightings
  } else { node_sightings <- node.size }
  
  if(is.data.frame(node.colour) == TRUE ) {
    nodes_list <- data.frame(node = rep(NA, length(unique(c(threshold_edges$node_1, threshold_edges$node_2)))), #as.numeric(names(net_all[[1]])),
                             age = NA)
    for(i in 1:nrow(nodes_list)){
      nodes_all <- rep(NA, 2*nrow(threshold_edges))  
      for(a in 1:2){
        for(b in 1:nrow(threshold_edges)){
          nodes_all[a + (b-1)*2] <- threshold_edges[b,a]
        }
      }
      nodes_list$node <- unique(nodes_all)
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

### create single matrix of edge samples
edge_samples <- matrix(data = NA, nrow = n_samples*n_chains, ncol = n_dyads)   # matrix for storing edge samples
for(j in 1:n_dyads){                                                           # for every dyad, fill matrix with weights (currently 4 columns per dyad as saved each chain separately)
  edge_samples[,j] <- edge_weights_matrix[,,j]
}
colnames(edge_samples) <- counts_df$dyad_id                                    # match to dyad ID numbers

edge_weights_matrix[1,1,]

### plot network across 6 different threshold values for comparison to other networks
plot_network_threshold_mpnp(edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.05,
                           node.size = nodes, node.colour = nodes, lwd = 15)
plot_network_threshold_mpnp(edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.10,
                           node.size = nodes, node.colour = nodes, lwd = 15)
plot_network_threshold_mpnp(edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.15,
                           node.size = nodes, node.colour = nodes, lwd = 15)
plot_network_threshold_mpnp(edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.20,
                           node.size = nodes, node.colour = nodes, lwd = 15)
plot_network_threshold_mpnp(edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.25,
                           node.size = nodes, node.colour = nodes, lwd = 15)

### clean workspace
rm(x, edge_weights_matrix, i, j) ; gc()

### save image
save.image('mpnp_edgecalculations/mpnplong_edgeweights_conditionalprior.RData')

### add progress marker
print(paste0('networks plotted at ', Sys.time()))

#### check outputs: compare edge weight distributions to simple SRI ####
### load in workspace image with edge weights already calculated
#load('mpnp_edgecalculations/mpnplong_edgeweights_conditionalprior.RData')

### create data frame with SRI and bison weights in it
make_edgelist <- function (edge_samples, dyad_data) # function pulled directly from BISoN for extracting edge lists
{
  dyad_name <- do.call(paste, c(dyad_data[c("node_1", "node_2")], sep=" <-> "))
  edge_lower <- apply(edge_samples, 2, function(x) quantile(x, probs=0.025))
  edge_upper <- apply(edge_samples, 2, function(x) quantile(x, probs=0.975))
  edge_median <- apply(edge_samples, 2, function(x) quantile(x, probs=0.5))
  edge_list <- cbind(
    "median"=round(edge_median, 3), 
    "2.5%"=round(edge_lower, 3), 
    "97.5%"=round(edge_upper, 3)
  )
  rownames(edge_list) <- dyad_name
  edgelist <- as.data.frame(edge_list)
  edgelist$node_1 <- as.character(dyad_data$node_1)
  edgelist$node_2 <- as.character(dyad_data$node_2)
  edgelist <- edgelist[,c(4:5,1:3)]
  return(edgelist)
}
edgelist <- make_edgelist(edge_samples = edge_samples, dyad_data = counts_df)   # obtain edge list
head(edgelist)                                                                  # check structure of edgelist
edgelist$node_1 <- as.integer(edgelist$node_1)                                  # convert to integers
edgelist$node_2 <- as.integer(edgelist$node_2)                                  # convert to integers
summary <- left_join(edgelist, counts_df, by = c('node_1','node_2'))            # combine distribution data with raw counts
head(summary)
summary$sri <- summary$together / summary$count_dyad

### plot bison against SRI
plot(density(summary$sri), main = 'red = SRI, blue = BISoN', col = 'red', lwd = 2, las = 1) # plot SRI distribution
lines(density(summary$median), col = 'blue', lwd = 2)                                       # distribution of median estimates

### plot a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
#plot_edges$dyad_id <- as.numeric(plot_edges$dyad)
draws <- left_join(plot_edges, summary, by = 'dyad_id')# %>%                     # combine to add node and median value data
#  select(-dyad.x) %>% rename(dyad = dyad.y)                                     # clean up names
which(is.na(draws) == TRUE)                                                     # check fully merged
   
#draws$dyad_id <- reorder(draws$dyad_id, draws$count_dyad)                       # order data frame based on total sightings per pair
ggplot(data = draws, mapping = aes(x = edge_draw)) +                            # plot data frame
  geom_density(colour = 'blue') +                                               # plot density plots of probability distributions output per dyad
  facet_wrap(. ~ dyad_id, nrow = 20) +                                          # split by dyad, allow to vary height depending on dyad
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3) +    # add line showing where the median estimate is
  geom_vline(mapping = aes(xintercept = sri), colour = 'red') +                 # add line showing where the SRI value is
  theme(strip.text.x = element_blank(), strip.background = element_blank())     # remove facet labels so plots can be bigger

#draws$dyad_id <- reorder(draws$dyad_id, draws$sri)                              # order dataframe based on SRI
ggplot(data = draws, mapping = aes(x = edge_draw))+                             # plot dataframe
  geom_density(colour = 'blue')+                                                # plot density plots of probability distributions output per dyad
  facet_wrap(. ~ dyad_id, nrow = 5)+                                            # split by dyad, allow to vary height depending on dyad
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+     # add line showing where the median estimate is
  geom_vline(mapping = aes(xintercept = sri), colour = 'red') +                 # add line showing where the SRI value is
  theme(strip.text.x = element_blank(), strip.background = element_blank())     # remove facet labels so plots can be bigger

write_csv(draws, '../data_processed/step3_edgeweightestimation/mpnplong_sampledyads_conditionalprior.csv') # save output for future reference

### clean environment
rm(draws, plot_edges, n_test, plot_dyads) ; gc()
save.image('mpnp_edgecalculations/mpnplong_edgeweights_conditionalprior.RData')
dev.off()
rm(counts_df, edge_samples, edgelist, edges, fit_edges_mpnp1, nodes, summary, n_dyads) ; gc()

### add progress marker
print(paste0('long window completed at ', Sys.time()))

################ generate loop to run through short windows 3 to end (1 and 2 to be run separately using rstan instead of cmdstanr) ################
### set up values for running loop
n_windows <- 5      # number of short time windows in mpnp data

for(time_window in 3:n_windows){
  #### set up ####
  ### set seed
  set.seed(12345)
  
  ### create file of output graphs
  pdf(file = paste0('../outputs/mpnpshort',time_window,'_edgeweights_conditionalprior.pdf'), width = 20, height = 15)
  
  #### import data ####
  ### subset by time window
  counts_df <- read_csv(paste0('../data_processed/step1_dataprocessing/mpnp_period',time_window,'_pairwiseevents.csv'))
  
  ### create nodes data frame
  nodes <- data.frame(id = sort(unique(c(counts_df$id_1,counts_df$id_2))),  # all unique individuals
                      node = NA, sightings = NA)                  # data needed on each
  for(i in 1:nrow(nodes)){
    # extract data about individual from counts_df data frame
    if(nodes$id[i] %in% counts_df$id_1) {
      x <- counts_df[counts_df$id_1 == nodes$id[i], c('id_1','node_1','count_period_1')] %>% distinct()
    } else { x <- counts_df[counts_df$id_2 == nodes$id[i], c('id_2','node_2','count_period_2')] %>% distinct() }
    colnames(x) <- c('id','node','count')
    # add individual data
    nodes$node[i] <- x$node          # node ID number
    nodes$sightings[i] <- x$count    # number of sightings within the time window
  }
  
  ### import age data
  eles <- readRDS(paste0('../data_processed/step2_ageestimation/mpnp',time_window,'_ageestimates_mcmcoutput.rds'))
  mean_ages <- data.frame(id = colnames(eles),
                          age = apply(eles, 2, mean))
  
  # combine node counts with age data
  nodes <- left_join(nodes, mean_ages, by = 'id') # join age data to sightings
  length(which(is.na(nodes$age) == TRUE))         # count elephants without age data
  
  ### set values for model
  n_dyads <- nrow(counts_df)
  
  ### create data list
  counts_ls <- list(
    n_dyads    = n_dyads,                  # total number of times one or other of the dyad was observed
    dyad_ids   = counts_df$dyad_id,        # identifier for each dyad
    together   = counts_df$together,       # count number of sightings seen together
    count_dyad = counts_df$count_dyad      # count total number of times seen
  )
  
  ### add progress marker
  print(paste0('data imported for time window ',time_window,' at ', Sys.time()))
  
  #### run model on real standardised data ####
  ### Fit model
  fit_edges_mpnp <- edge_binary$sample(
    data = counts_ls, 
    chains = n_chains, 
    parallel_chains = n_chains,
    iter_warmup = n_samples,
    iter_sampling = n_samples)
  
  ### save model
  save.image(paste0('mpnp_edgecalculations/mpnpshort',time_window,'_edgeweights_conditionalprior.RData'))
  
  ### check model
  fit_edges_mpnp
  
  ### add progress marker
  print(paste0('model run for time window ',time_window,' at ', Sys.time()))
  
  # view summary
  fit_edges_mpnp$summary()
  
  # Extract posterior samples
  posterior_samples <- fit_edges_mpnp$draws()
  
  # Convert the array to a matrix -- save for eigenvector centralities
  edge_weights_matrix <- posterior_samples[,,2:(nrow(counts_df)+1)]
  rm(posterior_samples, x, i) ; gc()
  
  # convert matrix to data frame -- save samples and plot outputs
  edges <- as.data.frame(edge_weights_matrix[,,1])
  colnames(edges) <- c('chain1','chain2','chain3','chain4')
  edges <- pivot_longer(edges, everything(), values_to = 'edge_draw', names_to = 'chain')
  edges$dyad <- counts_ls$dyad_ids[1]
  edges$position <- rep(1:n_samples, each = n_chains)
  for(i in 2:n_dyads){
    x <- as.data.frame(edge_weights_matrix[,,i])
    colnames(x) <- c('chain1','chain2','chain3','chain4')
    x <- pivot_longer(x, everything(), values_to = 'edge_draw', names_to = 'chain')
    x$dyad <- counts_ls$dyad_ids[i]
    x$position <- rep(1:n_samples, each = n_chains)
    edges <- rbind(edges, x)
  }
  
  ### save data 
  saveRDS(edges, paste0('../data_processed/step3_edgeweightestimation/mpnpshort',time_window,'_edgedistributions_conditionalprior.RDS'))
  
  ### add progress marker
  print(paste0('edges extracted for time window ',time_window,' at ', Sys.time()))
  
  #### check outputs ####
  # Assign random set of columns to check -- maximum 200 of each type, but with same number of 'zero' dyads as 'non-zeroes'
  if(length(which(counts_df$together >= 1)) >= 200){ n_test <- 200 } else { n_test <- length(which(counts_df$together >= 1)) } # identify number of samples to include
  plot_dyads <- c(sample(counts_df$dyad_id[counts_df$together >= 1], size = n_test, replace = F),    # sample 'non-zero' dyads
                  sample(counts_df$dyad_id[counts_df$together == 0], size = n_test, replace = F))    # sample 'zero' dyads
  plot_edges <- edges[edges$dyad %in% plot_dyads,]                                                   # extract edges for sampled dyads
  plot_edges$seen_together <- NA
  for(i in 1:length(plot_dyads)){                                      # set up loop to make 0/1 dummy variable for seen vs not seen together
    plot_edges$seen_together[plot_edges$dyad == plot_dyads[i]] <- ifelse(counts_df$together[counts_df$dyad_id == plot_dyads[i]] > 0, 1, 0)
  }
  
  ### build traceplots -- split seen and not seen to check they are different and look for anomalies
  ggplot(data = plot_edges[plot_edges$seen_together == 1,], aes(y = edge_draw, x = position, colour = chain))+
    geom_line()+
    facet_wrap(. ~ dyad_id)+
    theme_classic()+
    theme(legend.position = 'none', strip.background = element_blank(), strip.text = element_blank())
  ggplot(data = plot_edges[plot_edges$seen_together == 0,], aes(y = edge_draw, x = position, colour = chain))+
    geom_line()+
    facet_wrap(. ~ dyad_id)+
    theme_classic()+
    theme(legend.position = 'none', strip.background = element_blank(), strip.text = element_blank())
  
  ### density plots -- split seen and not seen to check they are different and look for anomalies
  plot(NULL, xlim = c(0,1), ylim = c(0,50), las = 1, xlab = 'edge weight', ylab = 'density')
  for(i in 1:length(plot_dyads)){
    x <- plot_edges[plot_edges$dyad == plot_dyads[i],]
    lines(density(x$edge_draw), col = ifelse(x$seen_together == 1, rgb(0,0,1,0.1), rgb(1,0,0,0.1)))
  }
  
  ### add progress marker
  print(paste0('edges checked for time window ',time_window,' at ', Sys.time()))
  
  #### check outputs: plot network ####
  ### create single matrix of edge samples
  edge_samples <- matrix(data = NA, nrow = n_samples*n_chains, ncol = n_dyads)   # matrix for storing edge samples
  for(j in 1:n_dyads){                                                           # for every dyad, fill matrix with weights (currently 4 columns per dyad as saved each chain separately)
    edge_samples[,j] <- edge_weights_matrix[,,j]
  }
  colnames(edge_samples) <- counts_df$dyad_id                                    # match to dyad ID numbers
  
  save.image(paste0('mpnp_edgecalculations/mpnpshort',time_window,'_edgeweights_conditionalprior.RData'))
  
  ### plot network across 6 different threshold values for comparison to other networks
  plot_network_threshold_mpnp(edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.05,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  plot_network_threshold_mpnp(edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.10,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  plot_network_threshold_mpnp(edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.15,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  plot_network_threshold_mpnp(edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.20,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  #plot_network_threshold_mpnp(edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.25,
  #                           node.size = nodes, node.colour = nodes, lwd = 15)
  
  ### clean workspace
  rm(counts_ls, x, edge_weights_matrix, i, j) ; gc()
  
  ### save image
  save.image(paste0('mpnp_edgecalculations/mpnpshort',time_window,'_edgeweights_conditionalprior.RData'))
  
  ### add progress marker
  print(paste0('network plotted for time window ',time_window,' at ', Sys.time()))
  
  #### check outputs: compare edge weight distributions to simple SRI ####
  ### create data frame with SRI and bison weights in it
  edgelist <- make_edgelist(edge_samples = edge_samples, dyad_data = counts_df) # obtain edge list
  head(edgelist)                                                                # check structure of edgelist
  edgelist$node_1 <- as.integer(edgelist$node_1)                                # convert to integers
  edgelist$node_2 <- as.integer(edgelist$node_2)                                # convert to integers
  summary <- left_join(edgelist, counts_df, by = c('node_1','node_2'))          # combine distribution data with raw counts
  head(summary)
  summary$sri <- summary$together / summary$count_dyad
  
  ### plot bison against SRI
  plot(density(summary$sri), main = 'red = SRI, blue = BISoN', col = 'red', lwd = 2, las = 1) # plot SRI distribution
  lines(density(summary$median), col = 'blue', lwd = 2)                                       # distribution of median estimates
  
  ### plot a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
  plot_edges$dyad_id <- as.numeric(plot_edges$dyad)
  draws <- left_join(plot_edges, summary, by = 'dyad_id') %>%                  # combine to add node and median value data
    select(-dyad.x) %>% rename(dyad = dyad.y)                                  # clean up names
  which(is.na(draws) == TRUE)                                                  # check fully merged
  
  draws$dyad_id <- reorder(draws$dyad_id, draws$count_dyad)                    # order data frame based on total sightings per pair
  ggplot(data = draws, mapping = aes(x = edge_draw)) +                         # plot data frame
    geom_density(colour = 'blue') +                                            # plot density plots of probability distributions output per dyad
    facet_wrap(. ~ dyad_id, nrow = 20) +                                       # split by dyad, allow to vary height depending on dyad
    geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3) + # add line showing where the median estimate is
    geom_vline(mapping = aes(xintercept = sri), colour = 'red') +              # add line showing where the SRI value is
    theme(strip.text.x = element_blank(), strip.background = element_blank())  # remove facet labels so plots can be bigger
  
  draws$dyad_id <- reorder(draws$dyad_id, draws$sri)                           # order dataframe based on SRI
  ggplot(data = draws, mapping = aes(x = edge_draw))+                          # plot dataframe
    geom_density(colour = 'blue')+                                             # plot density plots of probability distirbtuions output per dyad
    facet_wrap(. ~ dyad_id, nrow = 20)+                                        # split by dyad, allow to vary height depending on dyad
    geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+  # add line showing where the median estimate is
    geom_vline(mapping = aes(xintercept = sri), colour = 'red') +              # add line showing where the SRI value is
    theme(strip.text.x = element_blank(), strip.background = element_blank())  # remove facet labels so plots can be bigger
  ggplot(data = draws, mapping = aes(x = edge_draw))+                          # repeat but with free y axes
    geom_density(colour = 'blue')+                                             # plot density plots of probability distirbtuions output per dyad
    facet_wrap(. ~ dyad_id, nrow = 20, scales = 'free_y')+                     # split by dyad, allow to vary height depending on dyad
    geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+  # add line showing where the median estimate is
    geom_vline(mapping = aes(xintercept = sri), colour = 'red')+               # add line showing where the SRI value is
    theme(strip.text.x = element_blank(), strip.background = element_blank())  # remove facet labels so plots can be bigger
  
  write_csv(draws, paste0('../data_processed/step3_edgeweightestimation/mpnpshort',time_window,'_sampledyads_conditionalprior.csv')) # save output for future reference
  
  ### add progress marker
  print(paste0('time window ',time_window,' completed at ', Sys.time()))
  
  #### clean environment ####
  rm(draws, plot_edges, n_test, plot_dyads) ; gc()
  save.image(paste0('mpnp_edgecalculations/mpnpshort',time_window,'_edgeweights_conditionalprior.RData'))
  dev.off()
  rm(counts_df, edge_samples, edgelist, edges, fit_edges_mpnp, nodes, summary, n_dyads) ; gc()
  
}

#### plot all short outputs of 0.15 for thesis ####
n_windows <- 5
for(time_window in 1:n_windows){
  pdf(paste0('../outputs/step3_edgeweightestimation/mpnpshort',time_window,'_network_0.15.png'))
  
  load(paste0('mpnp_edgecalculations/mpnpshort',time_window,'_edgeweights_conditionalprior.RData'))
  
  plot_network_threshold_mpnp(edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.15,
                              node.size = nodes, node.colour = nodes, lwd = 15)
  
  dev.off()
  rm(list = ls()[!ls() %in% c('time_window','n_windows')]) ; gc()
  
  print(time_window)
}
