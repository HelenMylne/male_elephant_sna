#### Bayesian analysis of ATE data ####
# Script to process association data from Amboseli National Park, Kenya
# Data collected: 1972-2021 by Amboseli Trust for Elephants
# Data supplied by: Vicki Fishlock (March 2022) and Phyllis Lee (December 2021)
# Data input: raw data provided by ATE processed in scripts 22.03.20_anp_dataprocessing1.R and 22.03.22_anp_dataprocessing2.R

#### set up ####
### install packages
# install.packages('tidyverse')
# install.packages('dplyr')
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages(c("StanHeaders","rstan"),type="source")
# install.packages("remotes") ; remotes::install_github("stan-dev/cmdstanr") ## OR USE # install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages('rethinking')
# install.packages('igraph')
# install.packages('dagitty')
# install.packages('janitor')
# install.packages('lubridate')
# install.packages('hms')
# install.packages('readxl')

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
#set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0')

### set seed
set.seed(12345)

### create file of output graphs
pdf('../outputs/anpshort1_edgeweights_conditionalprior.pdf', width = 20, height = 15)

### load edge weight model
edge_binary <- cmdstan_model("models/edge_binary_basic.stan")   # load model
edge_binary                                                     # check model priors etc.

#### import data ####
### import aggregated counts of sightings together and apart, and info about individuals
counts_df <- read_csv('../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_shortwindows_impossiblepairsremoved.csv')  # counts_df = data frame of all dyad information, spceifically including aggregated counts of together vs apart

### set up values for running loop
periods <- sort(unique(c(counts_df$period_start, counts_df$period_end)))  # dates  of short time windows in ANP data
n_windows <- length(unique(counts_df$period))                             # number of short time windows in ANP data

### subset by time window
cdf_1 <- counts_df[counts_df$period == 1,]                                # select data for first time window only

### create nodes data frame for period 1
nodes <- data.frame(id = sort(unique(c(cdf_1$id_1,cdf_1$id_2))),          # all unique individuals
                    node = NA, age = NA, sightings = NA)                  # data needed on each
for(i in 1:nrow(nodes)){
  # extract data about individual from cdf_1 data frame
  if(nodes$id[i] %in% cdf_1$id_1) {
    x <- cdf_1[cdf_1$id_1 == nodes$id[i], c('id_1','node_1','period_count_1','age_start_1')] %>% distinct()
  } else { x <- cdf_1[cdf_1$id_2 == nodes$id[i], c('id_2','node_2','period_count_2','age_start_2')] %>% distinct() }
  colnames(x) <- c('id','node','period_count','age_start')
  # add individual data
  nodes$node[i] <- x$node                 # node ID number
  nodes$age[i] <- x$age_start             # estimated age in years at the start of the time window
  nodes$sightings[i] <- x$period_count    # number of sightings within the time window
}

### create data list
n_chains <- 4                            # number of MCMC chains to run
n_samples <- 1000                        # number of samples per chain
n_dyads <- nrow(cdf_1)                   # number of dyads in time window
counts_ls <- list(                       # create data list
  n_dyads    = n_dyads,                  # total number of times one or other of the dyad was observed
  dyad_ids   = cdf_1$dyad_id,            # identifier for each dyad
  together   = cdf_1$event_count,        # count number of sightings seen together
  count_dyad = cdf_1$period_count_dyad   # count total number of times seen
)

#### run model on real standardised data -- period 1 ####
### Fit model
fit_edges_anp1 <- edge_binary$sample(
  data = counts_ls, 
  chains = n_chains, parallel_chains = n_chains,
  iter_warmup = n_samples, iter_sampling = n_samples)

### check model
fit_edges_anp1

### Extract posterior samples
posterior_samples <- fit_edges_anp1$draws()

### extract edge weights
edge_weights_matrix <- posterior_samples[,,2:(nrow(cdf_1)+1)]  # remove column of intercepts
rm(posterior_samples, x, i) ; gc()                             # clean workspace

### convert edge samples to data frame, putting all chains into 1 column
edges <- as.data.frame(edge_weights_matrix[,,1])               # convert draws for first dyad to data frame
colnames(edges) <- c('chain1','chain2','chain3','chain4')      # rename columns
edges <- pivot_longer(edges, everything(), values_to = 'edge_draw', names_to = 'chain') # tidy draws for first dyad
edges$dyad <- counts_ls$dyad_ids[1]                            # add dyad id column
edges$position <- rep(1:n_samples, each = n_chains)            # position within chain for traceplots
for(i in 2:n_dyads){                                           # repeat for the rest of the dyads
  x <- as.data.frame(edge_weights_matrix[,,i])
  colnames(x) <- c('chain1','chain2','chain3','chain4')
  x <- pivot_longer(x, everything(), values_to = 'edge_draw', names_to = 'chain')
  x$dyad <- counts_ls$dyad_ids[i]
  x$position <- rep(1:n_samples, each = n_chains)
  edges <- rbind(edges, x)                                     # append each to the last to create a single dataframe with all dyads
}

### save edge samples for use in dyadic/nodal regressions 
saveRDS(edges, '../data_processed/anpshort1_edgedistributions_conditionalprior.RDS')
#edges <- readRDS('../data_processed/anpshort1_edgedistributions_conditionalprior.RDS')

#### check outputs: edge weights ####
### Assign random set of columns to check
if(length(which(cdf_1$event_count >= 1)) >= 200){ n_test <- 200 } else { n_test <- length(which(cdf_1$event_count >= 1)) } # set number of dyads to sample: will be equal numbers of dyads that have and have not ever been seen together, so if there are 200 or more dyads that have been seen together at some point then cap at 200, otherwise use as many as there are
plot_dyads <- c(sample(cdf_1$dyad_id[cdf_1$event_count >= 1], size = n_test, replace = F), # randomly sample dyads that have been seen together
                sample(cdf_1$dyad_id[cdf_1$event_count == 0], size = n_test, replace = F)) # randomly sample dyads never seen together
plot_edges <- edges[edges$dyad %in% plot_dyads,]                    # subset edge weight data using randomly drawn dyads
plot_edges$seen_together <- NA ; for(i in 1:length(plot_dyads)){    # set up for loop
  plot_edges$seen_together[plot_edges$dyad == plot_dyads[i]] <- ifelse(cdf_1$event_count[cdf_1$dyad_id == plot_dyads[i]] > 0, 1, 0) # the value of seen_together is 1 if the dyad has ever been seen in the same group, and 0 if they have not
}

### build traceplots
ggplot(data = plot_edges[plot_edges$seen_together == 1,],   # plot only dyads that have been seen together
       aes(y = edge_draw, x = position, colour = chain))+   # plot all chains over each other for each dyad to check mixing
  geom_line()+                                              # draw as line plot
  facet_wrap(. ~ dyad)+                                     # new panel per dyad as each pair has own edge weight parameter
  theme_classic()+                                          # make it look nicer
  theme(legend.position = 'none',                           # remove legend
        strip.background = element_blank(), strip.text = element_blank())    # remove facet strips
ggplot(data = plot_edges[plot_edges$seen_together == 0,],   # repeat for dyads that have never been seen together
       aes(y = edge_draw, x = position, colour = chain))+
  geom_line()+
  facet_wrap(. ~ dyad)+
  theme_classic() + theme(legend.position = 'none', strip.background = element_blank(), strip.text = element_blank())

### density plots
plot(NULL, xlim = c(0,1), ylim = c(0,30), las = 1, xlab = 'edge weight', ylab = 'density')         # set up plot window
for(i in 1:length(plot_dyads)){                                                                    # plot randomly sampled dyads
  x <- plot_edges[plot_edges$dyad == plot_dyads[i],]                                               # select data to plot
  lines(density(x$edge_draw), col = ifelse(x$seen_together == 1, rgb(0,0,1,0.1), rgb(1,0,0,0.1)))  # draw edge weight probability plot. blue = seen together at least once, red = never seen together
}

#### check outputs: plot network ####
### create plotting function
plot_network_threshold_anp <- function (edge_samples, dyad_data, lwd = 2, threshold = 0.3,
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
    node_sightings <- nodes_list$sightings*8 #log(nodes_list$sightings)*5
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
edge_samples <- matrix(data = NA, nrow = n_samples*n_chains, ncol = n_dyads)
for(j in 1:n_dyads){
  edge_samples[,j] <- edge_weights_matrix[,,j]
}
colnames(edge_samples) <- cdf_1$dyad_id

### plot network
plot_network_threshold_anp(edge_samples = edge_samples, dyad_data = cdf_1, threshold = 0.05,
                           node.size = nodes, node.colour = nodes, lwd = 15)
plot_network_threshold_anp(edge_samples = edge_samples, dyad_data = cdf_1, threshold = 0.10,
                           node.size = nodes, node.colour = nodes, lwd = 15)
plot_network_threshold_anp(edge_samples = edge_samples, dyad_data = cdf_1, threshold = 0.15,
                           node.size = nodes, node.colour = nodes, lwd = 15)
plot_network_threshold_anp(edge_samples = edge_samples, dyad_data = cdf_1, threshold = 0.20,
                           node.size = nodes, node.colour = nodes, lwd = 15)
plot_network_threshold_anp(edge_samples = edge_samples, dyad_data = cdf_1, threshold = 0.25,
                           node.size = nodes, node.colour = nodes, lwd = 15)
plot_network_threshold_anp(edge_samples = edge_samples, dyad_data = cdf_1, threshold = 0.30,
                           node.size = nodes, node.colour = nodes, lwd = 15)

### clean workspace
rm(counts_ls, x, edge_weights_matrix, i, j) ; gc()

### save image
save.image('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')
#load('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')

#### check outputs: compare edge weight distributions to simple SRI ####
### load in workspace image with edge weights already calculated
#load('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')

### create data frame with SRI and bison weights in it
make_edgelist <- function (edge_samples, dyad_data)
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
edgelist <- make_edgelist(edge_samples = edge_samples, dyad_data = cdf_1)     # obtain edge list
head(edgelist)                                                                # check structure of edgelist
edgelist$node_1 <- as.integer(edgelist$node_1)                                # convert to integers
edgelist$node_2 <- as.integer(edgelist$node_2)                                # convert to integers
summary <- left_join(edgelist, cdf_1, by = c('node_1','node_2'))    # combine distribution data with raw counts
head(summary)
summary$sri <- summary$event_count / summary$period_count_dyad

### plot bison against SRI
plot(density(summary$sri), main = 'red = SRI, blue = BISoN', col = 'red', lwd = 2, las = 1) # plot SRI distribution
lines(density(summary$median), col = 'blue', lwd = 2)                                       # distribution of median estimates

### plot a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
plot_edges$dyad_id <- as.numeric(plot_edges$dyad)
draws <- left_join(plot_edges, summary, by = 'dyad_id') %>%                  # combine to add node and median value data
  select(-dyad.x) %>% rename(dyad = dyad.y)                                  # clean up names
which(is.na(draws) == TRUE)                                                  # check fully merged

draws$dyad_id <- reorder(draws$dyad_id, draws$period_count_dyad)             # order data frame based on total sightings per pair
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

write_csv(draws, '../data_processed/anpshort1_sampledyads_conditionalprior.csv') # save output for future reference

### clean environment
rm(draws, plot_edges, n_test, plot_dyads) ; gc()
save.image('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')
dev.off()
rm(cdf_1, edge_samples, edgelist, edges, fit_edges_anp1, nodes, summary, n_dyads) ; gc()

################ generate loop to run through windows 2:end ################
for(time_window in 2:n_windows){
  #### set up ####
  ### set seed
  set.seed(12345)
  
  ### create file of output graphs
  pdf(file = paste0('../outputs/anpshort',time_window,'_edgeweights_conditionalprior_TEST.pdf'), width = 20, height = 15)
  
  #### import data ####
  ### subset by time window
  cdf <- counts_df[counts_df$period == time_window,]
  
  ### create nodes data frame
  nodes <- data.frame(id = sort(unique(c(cdf$id_1,cdf$id_2))),          # all unique individuals
                      node = NA, age = NA, sightings = NA)              # data needed on each
  for(i in 1:nrow(nodes)){
    # extract data about individual from cdf data frame
    if(nodes$id[i] %in% cdf$id_1) {
      x <- cdf[cdf$id_1 == nodes$id[i], c('id_1','node_1','period_count_1','age_start_1')] %>%
        distinct()
    } else { 
      x <- cdf[cdf$id_2 == nodes$id[i], c('id_2','node_2','period_count_2','age_start_2')] %>%
        distinct()
      }
    colnames(x) <- c('id','node','period_count','age_start')
    # add individual data
    nodes$node[i] <- x$node
    nodes$age[i] <- x$age_start
    nodes$sightings[i] <- x$period_count
  }
  
  ### set values for model
  n_dyads <- nrow(cdf)
  
  ### create data list
  counts_ls <- list(
    n_dyads    = n_dyads,                # total number of times one or other of the dyad was observed
    dyad_ids   = cdf$dyad_id,            # identifier for each dyad
    together   = cdf$event_count,        # count number of sightings seen together
    count_dyad = cdf$period_count_dyad)  # count total number of times seen
  
  #### run model on real standardised data -- period 1 ####
  ### Fit model
  fit_edges_anp <- edge_binary$sample(
    data = counts_ls, 
    chains = n_chains, 
    parallel_chains = n_chains,
    iter_warmup = n_samples,
    iter_sampling = n_samples)
  
  ### check model
  fit_edges_anp
  
  # Extract posterior samples
  posterior_samples <- fit_edges_anp$draws()
  
  # Convert the array to a matrix -- save for eigenvector centralities
  edge_weights_matrix <- posterior_samples[,,2:(nrow(cdf)+1)]
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
  saveRDS(edges, paste0('../data_processed/anpshort',time_window,'_edgedistributions_conditionalprior.RDS'))
  
  #### check outputs ####
  # Assign random set of columns to check
  if(length(which(cdf$event_count >= 1)) >= 200){ n_test <- 200 } else { n_test <- length(which(cdf$event_count >= 1)) }
  plot_dyads <- c(sample(cdf$dyad_id[cdf$event_count >= 1], size = n_test, replace = F),
                  sample(cdf$dyad_id[cdf$event_count == 0], size = n_test, replace = F))
  plot_edges <- edges[edges$dyad %in% plot_dyads,]
  plot_edges$seen_together <- NA
  for(i in 1:length(plot_dyads)){
    plot_edges$seen_together[plot_edges$dyad == plot_dyads[i]] <- ifelse(cdf$event_count[cdf$dyad_id == plot_dyads[i]] > 0, 1, 0)
  }
  
  ### build traceplots
  ggplot(data = plot_edges[plot_edges$seen_together == 1,], aes(y = edge_draw, x = position, colour = chain))+
    geom_line()+
    facet_wrap(. ~ dyad)+
    theme_classic()+
    theme(legend.position = 'none', strip.background = element_blank(), strip.text = element_blank())
  ggplot(data = plot_edges[plot_edges$seen_together == 0,], aes(y = edge_draw, x = position, colour = chain))+
    geom_line()+
    facet_wrap(. ~ dyad)+
    theme_classic()+
    theme(legend.position = 'none', strip.background = element_blank(), strip.text = element_blank())
  
  ### density plots
  plot(NULL, xlim = c(0,1), ylim = c(0,50), las = 1, xlab = 'edge weight', ylab = 'density')
  for(i in 1:length(plot_dyads)){
    x <- plot_edges[plot_edges$dyad == plot_dyads[i],]
    lines(density(x$edge_draw), col = ifelse(x$seen_together == 1, rgb(0,0,1,0.1), rgb(1,0,0,0.1)))
  }
  
  #### check outputs: plot network ####
  ### create single matrix of edge samples
  edge_samples <- matrix(data = NA, nrow = n_samples*n_chains, ncol = n_dyads)
  for(j in 1:n_dyads){
    edge_samples[,j] <- edge_weights_matrix[,,j]
  }
  colnames(edge_samples) <- cdf$dyad_id
  
  ### plot network
  plot_network_threshold_anp(edge_samples = edge_samples, dyad_data = cdf, threshold = 0.05,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  plot_network_threshold_anp(edge_samples = edge_samples, dyad_data = cdf, threshold = 0.10,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  plot_network_threshold_anp(edge_samples = edge_samples, dyad_data = cdf, threshold = 0.15,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  plot_network_threshold_anp(edge_samples = edge_samples, dyad_data = cdf, threshold = 0.20,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  plot_network_threshold_anp(edge_samples = edge_samples, dyad_data = cdf, threshold = 0.25,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  plot_network_threshold_anp(edge_samples = edge_samples, dyad_data = cdf, threshold = 0.30,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  
  ### clean workspace
  rm(counts_ls, x, edge_weights_matrix, i, j) ; gc()
  
  ### save image
  save.image(paste0('anp_edgecalculations/anpshort',time_window,'_edgeweights_conditionalprior.RData'))
  
  #### check outputs: compare edge weight distributions to simple SRI ####
  ### create data frame with SRI and bison weights in it
  edgelist <- make_edgelist(edge_samples = edge_samples, dyad_data = cdf)     # obtain edge list
  head(edgelist)                                                                # check structure of edgelist
  edgelist$node_1 <- as.integer(edgelist$node_1)                                # convert to integers
  edgelist$node_2 <- as.integer(edgelist$node_2)                                # convert to integers
  summary <- left_join(edgelist, cdf, by = c('node_1','node_2'))    # combine distribution data with raw counts
  head(summary)
  summary$sri <- summary$event_count / summary$period_count_dyad
  
  ### plot bison against SRI
  plot(density(summary$sri), main = 'red = SRI, blue = BISoN', col = 'red', lwd = 2, las = 1) # plot SRI distribution
  lines(density(summary$median), col = 'blue', lwd = 2)                                       # distribution of median estimates
  
  ### plot a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
  plot_edges$dyad_id <- as.numeric(plot_edges$dyad)
  draws <- left_join(plot_edges, summary, by = 'dyad_id') %>%                  # combine to add node and median value data
    select(-dyad.x) %>% rename(dyad = dyad.y)                                  # clean up names
  which(is.na(draws) == TRUE)                                                  # check fully merged
  
  draws$dyad_id <- reorder(draws$dyad_id, draws$period_count_dyad)             # order data frame based on total sightings per pair
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
  
  write_csv(draws, paste0('../data_processed/anpshort',time_window,'_sampledyads_conditionalprior.csv')) # save output for future reference
  
  ### clean environment ####
  rm(draws, plot_edges, n_test, plot_dyads) ; gc()
  save.image(paste0('anp_edgecalculations/anpshort',time_window,'_edgeweights_conditionalprior.RData'))
  dev.off()
  rm(cdf, edge_samples, edgelist, edges, fit_edges_anp1, nodes, summary, n_dyads) ; gc()
  
}

## clean up for long windows
rm(list = ls()[! ls() %in% c('make_edgelist','plot_network_threshold','n_samples','n_chains','edge_binary')])

################ generate loop to run through long windows ################
### import aggregated counts of sightings together and apart, and info about individuals
counts_df <- read_csv('../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_longperiods_impossiblepairsremoved.csv')  # counts_df = data frame of all dyad information, spceifically including aggregated counts of together vs apart

### set up values for running loop
periods <- sort(unique(c(counts_df$period_start, counts_df$period_end)))  # dates  of long time windows in ANP data
n_windows <- length(unique(counts_df$period))                             # number of long time windows in ANP data

### run loop
for(time_window in 1:n_windows){
  #### set up ####
  ### set seed
  set.seed(12345)
  
  ### create file of output graphs
  pdf(paste0('../outputs/anplong',time_window,'_edgeweights_conditionalprior.pdf'), width = 20, height = 15)
  
  #### import data ####
  ### subset by time window
  cdf <- counts_df[counts_df$period == time_window,]
  
  ### create nodes data frame
  nodes <- data.frame(id = sort(unique(c(cdf$id_1,cdf$id_2))),          # all unique individuals
                      node = NA, age = NA, sightings = NA)              # data needed on each
  for(i in 1:nrow(nodes)){
    # extract data about individual from cdf data frame
    if(nodes$id[i] %in% cdf$id_1) {
      x <- cdf[cdf$id_1 == nodes$id[i], c('id_1','node_1','period_count_1','age_start_1')] %>%
        distinct()
    } else { 
      x <- cdf[cdf$id_2 == nodes$id[i], c('id_2','node_2','period_count_2','age_start_2')] %>%
        distinct()
    }
    colnames(x) <- c('id','node','period_count','age_start')
    # add individual data
    nodes$node[i] <- x$node
    nodes$age[i] <- x$age_start
    nodes$sightings[i] <- x$period_count
  }
  
  ### set values for model
  n_dyads <- nrow(cdf)
  
  ### create data list
  counts_ls <- list(
    n_dyads    = n_dyads,                # total number of times one or other of the dyad was observed
    dyad_ids   = cdf$dyad_id,            # identifier for each dyad
    together   = cdf$event_count,        # count number of sightings seen together
    count_dyad = cdf$period_count_dyad)  # count total number of times seen
  
  #### run model on real standardised data -- period 1 ####
  ### Fit model
  fit_edges_anp <- edge_binary$sample(
    data = counts_ls, 
    chains = n_chains, 
    parallel_chains = n_chains,
    iter_warmup = n_samples,
    iter_sampling = n_samples)
  
  ### check model
  fit_edges_anp
  
  # Extract posterior samples
  posterior_samples <- fit_edges_anp$draws()
  
  # Convert the array to a matrix -- save for eigenvector centralities
  edge_weights_matrix <- posterior_samples[,,2:(nrow(cdf)+1)]
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
  saveRDS(edges, paste0('../data_processed/anplong',time_window,'_edgedistributions_conditionalprior.RDS'))
  
  #### check outputs ####
  # Assign random set of columns to check
  if(length(which(cdf$event_count >= 1)) >= 200){ n_test <- 200 } else { n_test <- length(which(cdf$event_count >= 1)) }
  plot_dyads <- c(sample(cdf$dyad_id[cdf$event_count >= 1], size = n_test, replace = F),
                  sample(cdf$dyad_id[cdf$event_count == 0], size = n_test, replace = F))
  plot_edges <- edges[edges$dyad %in% plot_dyads,]
  plot_edges$seen_together <- NA
  for(i in 1:length(plot_dyads)){
    plot_edges$seen_together[plot_edges$dyad == plot_dyads[i]] <- ifelse(cdf$event_count[cdf$dyad_id == plot_dyads[i]] > 0, 1, 0)
  }
  
  ### build traceplots
  ggplot(data = plot_edges[plot_edges$seen_together == 1,], aes(y = edge_draw, x = position, colour = chain))+
    geom_line()+
    facet_wrap(. ~ dyad)+
    theme_classic()+
    theme(legend.position = 'none', strip.background = element_blank(), strip.text = element_blank())
  ggplot(data = plot_edges[plot_edges$seen_together == 0,], aes(y = edge_draw, x = position, colour = chain))+
    geom_line()+
    facet_wrap(. ~ dyad)+
    theme_classic()+
    theme(legend.position = 'none', strip.background = element_blank(), strip.text = element_blank())
  
  ### density plots
  plot(NULL, xlim = c(0,1), ylim = c(0,50), las = 1, xlab = 'edge weight', ylab = 'density')
  for(i in 1:length(plot_dyads)){
    x <- plot_edges[plot_edges$dyad == plot_dyads[i],]
    lines(density(x$edge_draw), col = ifelse(x$seen_together == 1, rgb(0,0,1,0.1), rgb(1,0,0,0.1)))
  }
  
  #### check outputs: plot network ####
  ### create single matrix of edge samples
  edge_samples <- matrix(data = NA, nrow = n_samples*n_chains, ncol = n_dyads)
  for(j in 1:n_dyads){
    edge_samples[,j] <- edge_weights_matrix[,,j]
  }
  colnames(edge_samples) <- cdf$dyad_id
  
  ### plot network
  plot_network_threshold_anp(edge_samples = edge_samples, dyad_data = cdf, threshold = 0.05,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  plot_network_threshold_anp(edge_samples = edge_samples, dyad_data = cdf, threshold = 0.10,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  plot_network_threshold_anp(edge_samples = edge_samples, dyad_data = cdf, threshold = 0.15,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  plot_network_threshold_anp(edge_samples = edge_samples, dyad_data = cdf, threshold = 0.20,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  plot_network_threshold_anp(edge_samples = edge_samples, dyad_data = cdf, threshold = 0.25,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  plot_network_threshold_anp(edge_samples = edge_samples, dyad_data = cdf, threshold = 0.30,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  
  ### clean workspace
  rm(counts_ls, x, edge_weights_matrix, i, j) ; gc()
  
  ### save image
  save.image(paste0('anp_edgecalculations/anplong',time_window,'_edgeweights_conditionalprior.RData'))
  
  #### check outputs: compare edge weight distributions to simple SRI ####
  ### create data frame with SRI and bison weights in it
  edgelist <- make_edgelist(edge_samples = edge_samples, dyad_data = cdf)     # obtain edge list
  head(edgelist)                                                                # check structure of edgelist
  edgelist$node_1 <- as.integer(edgelist$node_1)                                # convert to integers
  edgelist$node_2 <- as.integer(edgelist$node_2)                                # convert to integers
  summary <- left_join(edgelist, cdf, by = c('node_1','node_2'))    # combine distribution data with raw counts
  head(summary)
  summary$sri <- summary$event_count / summary$period_count_dyad
  
  ### plot bison against SRI
  plot(density(summary$sri), main = 'red = SRI, blue = BISoN', col = 'red', lwd = 2, las = 1) # plot SRI distribution
  lines(density(summary$median), col = 'blue', lwd = 2)                                       # distribution of median estimates
  
  ### plot a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
  plot_edges$dyad_id <- as.numeric(plot_edges$dyad)
  draws <- left_join(plot_edges, summary, by = 'dyad_id') %>%                  # combine to add node and median value data
    select(-dyad.x) %>% rename(dyad = dyad.y)                                  # clean up names
  which(is.na(draws) == TRUE)                                                  # check fully merged
  
  draws$dyad_id <- reorder(draws$dyad_id, draws$period_count_dyad)             # order data frame based on total sightings per pair
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
  
  write_csv(draws, paste0('../data_processed/anplong',time_window,'_sampledyads_conditionalprior.csv')) # save output for future reference
  
  ### clean environment
  rm(draws, plot_edges, n_test, plot_dyads) ; gc()
  save.image(paste0('anp_edgecalculations/anplong',time_window,'_edgeweights_conditionalprior.RData'))
  dev.off()
  rm(cdf, edge_samples, edgelist, edges, fit_edges_anp1, nodes, summary, n_dyads) ; gc()
  
}
