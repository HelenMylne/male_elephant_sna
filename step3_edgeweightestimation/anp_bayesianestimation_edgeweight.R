#### Bayesian analysis of ATE data ####
# Script to process association data from Amboseli National Park, Kenya
# Data collected: 1972-2021 by Amboseli Trust for Elephants
# Data supplied by: Vicki Fishlock (March 2022) and Phyllis Lee (December 2021)
# Data input: raw data provided by ATE processed in scripts 22.03.20_anp_dataprocessing1.R and 22.03.22_anp_dataprocessing2.R

#### set up ####
# install packages
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

# load packages
# library(tidyverse) ; library(dplyr) ; library(cmdstanr) ; library(igraph) ; library(janitor) ; library(lubridate) ; library(hms) ; library(readxl)
library(tidyverse, lib.loc = 'packages/')   # library(tidyverse)
library(dplyr, lib.loc = 'packages/')       # library(dplyr)
#library(rstan, lib.loc = 'packages/')      # library(rstan)
library(cmdstanr, lib.loc = 'packages/')    # library(cmdstanr)
library(igraph, lib.loc = 'packages/')      # library(igraph)
library(janitor, lib.loc = 'packages/')     # library(janitor)
library(lubridate, lib.loc = 'packages/')   # library(lubridate)
library(hms, lib.loc = 'packages/')         # library(hms)
library(readxl, lib.loc = 'packages/')      # library(readxl)

# set stan path
#set_cmdstan_path('/Users/helen/.cmdstan/cmdstan-2.32.1')
#set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0')

# load edge weight model
edge_binary <- cmdstan_model("models/edge_binary_basic.stan")
edge_binary

# set seed
set.seed(12345)

# create file of output graphs
pdf('../data_processed/anp_edgeweights_period1.pdf', width = 20, height = 15)

#### create data lists ####
counts_df <- read_csv('../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_shortwindows_impossiblepairsremoved.csv')

periods <- sort(unique(c(counts_df$period_start, counts_df$period_end)))

males <- readxl::read_excel('../data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()
males$id <- paste0('M',males$casename)
ids <- sort(unique(ate$id))
males <- males %>% dplyr::filter(id %in% ids)

# softcode certain values
n_windows <- length(unique(counts_df$period))
n_males <- length(unique(males$id))

#### run model on real standardised data -- period 1 ####
### subset by time window
cdf_1 <- counts_df[counts_df$period == 1,]

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
  nodes$node[i] <- x$node
  nodes$age[i] <- x$age_start                                             # for initial test purposes, age = mean only
  nodes$sightings[i] <- x$period_count
}

### create data list
n_chains <- 4
n_samples <- 1000
n_dyads <- nrow(cdf_1)
counts_ls <- list(
  n_dyads    = n_dyads,                  # total number of times one or other of the dyad was observed
  dyad_ids   = cdf_1$dyad_id,            # identifier for each dyad
  together   = cdf_1$event_count,        # count number of sightings seen together
  count_dyad = cdf_1$period_count_dyad)  # count total number of times seen
n_eles <- length(unique(c(cdf_1$id_1, cdf_1$id_2)))

### Fit model
fit_edges_anp1 <- edge_binary$sample(
  data = counts_ls, 
  chains = n_chains, 
  parallel_chains = n_chains)

### check model
fit_edges_anp1
#variable        mean     median    sd   mad         q5        q95 rhat ess_bulk ess_tail
#lp__           -5125.17 -5124.42 29.94 29.99 -5174.70 -5077.04 1.00     1403     2395
#edge_weight[1]     0.04     0.02  0.04  0.03     0.00     0.13 1.00     4931     2188
#edge_weight[2]     0.05     0.03  0.06  0.04     0.00     0.17 1.00     4879     2130
#edge_weight[3]     0.05     0.03  0.06  0.04     0.00     0.17 1.00     4974     2342
#edge_weight[4]     0.04     0.02  0.04  0.03     0.00     0.12 1.00     4572     2201
#edge_weight[5]     0.05     0.03  0.05  0.04     0.00     0.16 1.00     4881     2182
#edge_weight[6]     0.04     0.03  0.05  0.03     0.00     0.15 1.00     4942     2154
#edge_weight[7]     0.05     0.03  0.06  0.04     0.00     0.17 1.00     4128     2029
#edge_weight[8]     0.05     0.03  0.06  0.04     0.00     0.16 1.00     4586     2144
#edge_weight[9]     0.05     0.03  0.06  0.04     0.00     0.17 1.00     4684     2108

# Extract posterior samples
posterior_samples <- fit_edges_anp1$draws()

# Convert the array to a matrix -- save for eigenvector centralities
edge_weights_matrix <- posterior_samples[,,2:1327]

# convert matrix to data frame -- save samples and plot outputs
edges <- as.data.frame(edge_weights_matrix[,,1])
colnames(edges) <- c('chain1','chain2','chain3','chain4')
edges <- pivot_longer(edges, everything(), values_to = 'edge_draw', names_to = 'chain')
edges$dyad <- counts_ls$dyad_ids[1]
edges$position <- rep(1:(n_samples/4), each = 4)
for(i in 2:n_dyads){
  x <- as.data.frame(edge_weights_matrix[,,i])
  colnames(x) <- c('chain1','chain2','chain3','chain4')
  x <- pivot_longer(x, everything(), values_to = 'edge_draw', names_to = 'chain')
  x$dyad <- counts_ls$dyad_ids[i]
  x$position <- rep(1:(n_samples/4), each = 4)
  edges <- rbind(edges, x)
}

### save data 
saveRDS(edges, '../data_processed/anp1_edgedistributions.RDS')

#### check outputs -- NETWORK PLOTS USING MY ADAPTED BISONR CODE NOT CURRENTLY WORKING, BUT COME BACK TO THAT ####
# Assign random set of columns to check
if(length(which(cdf_1$event_count >= 1)) >= 200){ n_test <- 200 } else { n_test <- length(which(cdf_1$event_count >= 1)) }
plot_dyads <- c(sample(cdf_1$dyad_id[cdf_1$event_count >= 1], size = n_test, replace = F),
                sample(cdf_1$dyad_id[cdf_1$event_count == 0], size = n_test, replace = F))
plot_edges <- edges[edges$dyad %in% plot_dyads,]
plot_edges$seen_together <- NA
for(i in 1:length(plot_dyads)){
  plot_edges$seen_together[plot_edges$dyad == plot_dyads[i]] <- ifelse(cdf_1$event_count[cdf_1$dyad_id == plot_dyads[i]] > 0, 1, 0)
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

### create plotting function
#plot_network_threshold_anp <- function (edge_samples, dyad_data, lwd = 2, threshold = 0.3,
#                                     label.colour = 'transparent', label.font = 'Helvetica', 
#                                     node.size = 4, node.colour = 'seagreen1',
#                                     link.colour1 = 'black', link.colour2 = rgb(0, 0, 0, 0.3))
#{
#  dyad_name <- do.call(paste, c(dyad_data[c("node_1", "node_2")], sep=" <-> "))
#  edge_lower <- apply(edge_samples, 2, function(x) quantile(x, probs=0.025))
#  edge_upper <- apply(edge_samples, 2, function(x) quantile(x, probs=0.975))
#  edge_median <- apply(edge_samples, 2, function(x) quantile(x, probs=0.5))
#  edge_list <- cbind(
#    "median"=round(edge_median, 3), 
#    "2.5%"=round(edge_lower, 3), 
#    "97.5%"=round(edge_upper, 3)
#  )
#  rownames(edge_list) <- dyad_name
#  edgelist <- as.data.frame(edge_list)
#  edgelist$node_1 <- dyad_data$node_1
#  edgelist$node_2 <- dyad_data$node_2
#  edgelist <- edgelist[,c(4:5,1:3)]
#  threshold_edges <- edgelist[edgelist$median >= threshold,]
#  net <- igraph::graph_from_edgelist(as.matrix(threshold_edges[, 1:2]), directed = F)
#  
#  if(is.data.frame(node.size) == TRUE ) {
#    nodes_list <- data.frame(node = as.numeric(names(net[1])),
#                             sightings = NA)
#    for(i in 1:nrow(nodes_list)){
#      nodes_list$sightings[i] <- nodes$sightings[which(nodes$node == nodes_list$node[i])]
#    }
#    node_sightings <- log(nodes_list$sightings)*5
#  } else { node_sightings <- node.size }
#  
#  if(is.data.frame(node.colour) == TRUE ) {
#    nodes_list <- data.frame(node = as.numeric(names(net[1])),
#                             age = NA)
#    for(i in 1:nrow(nodes_list)){
#      nodes_list$age[i] <- nodes$age[which(nodes$node == nodes_list$node[i])]
#    }
#    node_age <- nodes_list$age
#  } else { node_age <- node.colour }
#  
#  md <- threshold_edges[, 3]
#  ub <- threshold_edges[, 5]
#  coords <- igraph::layout_nicely(net)
#  igraph::plot.igraph(net, layout = coords,
#                      vertex.label.color = ifelse(is.null(label.colour) == TRUE,
#                                                  ifelse(node_age < 20, 'black', 'white'),
#                                                  label.colour),
#                      label.family = label.font,
#                      vertex.color = ifelse(node_age < 15, '#FDE725FF',
#                                            ifelse(node_age < 20, '#55C667FF',
#                                                   ifelse(node_age < 30, '#1F968BFF', 
#                                                          ifelse(node_age < 40, '#39568CFF', '#440154FF')))), 
#                      vertex.size = node_sightings,
#                      frame.color = NA, frame.width = 0,
#                      edge.color = NA, edge.arrow.size = 0, edge.width = 0)
#  igraph::plot.igraph(net, layout = coords, add = TRUE,
#                      vertex.label = NA, vertex.color = 'transparent', vertex.size = 0, 
#                      frame.color = NA, frame.width = 0,
#                      edge.color = link.colour1, edge.arrow.size = 0, edge.width = md * lwd)
#  igraph::plot.igraph(net, layout = coords, add = TRUE,
#                      vertex.label = NA, vertex.color = 'transparent', vertex.size = 0, 
#                      frame.color = NA, frame.width = 0,
#                      edge.color = link.colour2, edge.arrow.size = 0, edge.width = ub * lwd)
#}

### create single matrix of edge samples
#edge_samples <- matrix(data = NA, nrow = n_samples, ncol = n_dyads)
#for(j in 1:n_dyads){
#  edge_samples[,j] <- edge_weights_matrix[,,j]
#}
#colnames(edge_samples) <- cdf_1$dyad_id

### plot network
#plot_network_threshold2(edge_samples = edge_samples, dyad_data = cdf_1, threshold = 0.15,
#                        node.size = nodes, node.colour = nodes, lwd = 15)

### end pdf
dev.off()

### clear workspace
rm(cdf_1, n_dyads, counts_ls, n_eles, fit_edges_anp1, posterior_samples, edges, plot_dyads, plot_edges, nodes, edge_weights_matrix, edge_samples)

################ generate loop to run through windows 2:end ################
for(time_window in 2:n_windows){
  ### set up pdf of outputs
  pdf('../data_processed/anp_edgeweights_period',time_window,'.pdf', width = 20, height = 15)
  
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
  n_eles <- length(unique(c(cdf$id_1, cdf$id_2)))
  
  ### create data list
  counts_ls <- list(
    n_dyads    = n_dyads,                  # total number of times one or other of the dyad was observed
    dyad_ids   = cdf$dyad_id,            # identifier for each dyad
    together   = cdf$event_count,        # count number of sightings seen together
    count_dyad = cdf$period_count_dyad)  # count total number of times seen
  
  ### Fit model
  fit_edges <- edge_binary$sample(
    data = counts_ls, 
    chains = n_chains, 
    parallel_chains = n_chains,
    iter_warmup = n_samples,
    iter_sampling = n_samples)
  
  ### check model
  fit_edges
  
  # Extract posterior samples
  posterior_samples <- fit_edges$draws()
  
  # Convert the array to a matrix -- save for eigenvector centralities
  edge_weights_matrix <- posterior_samples[,,2:(length(posterior_samples)/n_chains*n_samples)]
  
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
  saveRDS(edges, paste0('../data_processed/anp',time_window,'_edgedistributions.RDS'))
  
  #### check outputs ####
  # Assign random set of columns to check
  if(length(which(cdf$event_count >= 1)) >= 200){
    n_test <- 200
  } else {
      n_test <- length(which(cdf$event_count >= 1))
      }
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
  
  ### create nodes data frame
  nodes <- data.frame(id = sort(unique(c(cdf$id_1,cdf$id_2))),          # all unique individuals
                      node = NA, age = NA, sightings = NA)                  # data needed on each
  for(i in 1:nrow(nodes)){
    # extract data about individual from cdf data frame
    if(nodes$id[i] %in% cdf$id_1) {
      x <- cdf[cdf$id_1 == nodes$id[i], c('id_1','node_1','period_count_1','age_start_1')] %>% distinct()
    } else { x <- cdf[cdf$id_2 == nodes$id[i], c('id_2','node_2','period_count_2','age_start_2')] %>% distinct() }
    colnames(x) <- c('id','node','period_count','age_start')
    # add individual data
    nodes$node[i] <- x$node
    nodes$age[i] <- x$age_start                                             # for initial test purposes, age = mean only
    nodes$sightings[i] <- x$period_count
  }
  
  ### create single matrix of edge samples
  edge_samples <- matrix(data = NA, nrow = n_samples, ncol = n_dyads)
  for(j in 1:n_dyads){
    edge_samples[,j] <- edge_weights_matrix[,,j]
  }
  colnames(edge_samples) <- cdf$dyad_id
  
  ### plot network
  #plot_network_threshold2(edge_samples = edge_samples, dyad_data = cdf, threshold = 0.15,
  #                        node.size = nodes, node.colour = nodes, lwd = 15)
  
  ### end pdf
  dev.off()
  
  ### clear workspace
  rm(cdf, n_dyads, counts_ls, n_eles, fit_edges, posterior_samples, edges, plot_dyads, plot_edges, nodes, edge_weights_matrix, edge_samples)
  
}
