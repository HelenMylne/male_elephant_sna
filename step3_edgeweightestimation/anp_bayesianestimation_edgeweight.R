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
library(tidyverse, lib.loc = '../packages/')   # library(tidyverse)
library(dplyr, lib.loc = '../packages/')       # library(dplyr)
#library(rstan, lib.loc = '../packages/')      # library(rstan)
library(cmdstanr, lib.loc = '../packages/')    # library(cmdstanr)
library(igraph, lib.loc = '../packages/')      # library(igraph)
library(janitor, lib.loc = '../packages/')     # library(janitor)
library(lubridate, lib.loc = '../packages/')   # library(lubridate)
library(hms, lib.loc = '../packages/')         # library(hms)
library(readxl, lib.loc = '../packages/')      # library(readxl)

# set stan path
#set_cmdstan_path('/Users/helen/.cmdstan/cmdstan-2.32.1')
#set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0')

# set seed
set.seed(12345)

# create file of output graphs
pdf('../outputs/anp_edgeweights_period1_conditionalprior_dyadreg.pdf', width = 20, height = 15)

#### prior predictive check ####
age_min <- 5:50
age_max <- seq(10, 50, length.out = 12)
par(mfrow = c(4,3))
for(age in age_max){
  plot(NULL, xlim = c(5,50), ylim = c(0,1), las = 1,
     xlab = 'age of younger', ylab = 'edge weight', main = paste0('oldest = ',round(age)))
  for(i in 1:100){
    intercept <- rnorm(1,-2,1)
    beta_min <- rnorm(1, 0, 0.1)
    beta_max <- rnorm(1, 0, 0.1)
    #beta_int <- rnorm(1, 0, 0.005)
    min_new <- age_min[age > age_min]
    lines(x = min_new, y = plogis(intercept + min_new*beta_min + 
                                  #min_new*age*beta_int + 
                                  age*beta_max), 
          col = rgb(0,0,1,0.2))
  }
}
rm(list = ls()) ; gc()
par(mfrow = c(1,1))

# load edge weight model
edge_binary <- cmdstan_model("models/edge_binary_with_dyadicregression.stan")
edge_binary

#### create data lists ####
counts_df <- read_csv('../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_shortwindows_impossiblepairsremoved.csv')

periods <- sort(unique(c(counts_df$period_start, counts_df$period_end)))

males <- readxl::read_excel('../data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()
males$id <- paste0('M',males$casename)
ids <- sort(unique(c(counts_df$id_1, counts_df$id_2)))
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

### identify older and younger of dyad
cdf_1$age_min <- NA ; cdf_1$age_max <- NA
for(i in 1:nrow(cdf_1)){
  x <- c(cdf_1$age_start_1[i],cdf_1$age_start_2[i])
  cdf_1$age_min[i] <- min(x)
  cdf_1$age_max[i] <- max(x)
}

### create data list
n_chains <- 4
n_samples <- 1000
n_dyads <- nrow(cdf_1)
counts_ls <- list(
  n_dyads    = n_dyads,                  # total number of times one or other of the dyad was observed
  dyad_ids   = cdf_1$dyad_id,            # identifier for each dyad
  together   = cdf_1$event_count,        # count number of sightings seen together
  count_dyad = cdf_1$period_count_dyad,  # count total number of times seen
  node_1     = cdf_1$node_1,
  node_2     = cdf_1$node_2,
  age_min    = cdf_1$age_min,
  age_max    = cdf_1$age_max
)
n_eles <- length(unique(c(cdf_1$id_1, cdf_1$id_2)))

### Fit model
fit_edges_anp1 <- edge_binary$sample(
  data = counts_ls, 
  chains = n_chains, 
  parallel_chains = n_chains)

### check model
fit_edges_anp1

# Extract posterior samples
posterior_samples <- fit_edges_anp1$draws()

# extract edge weights
edge_weights_matrix <- posterior_samples[,,2:(nrow(cdf_1)+1)]
posterior_samples <- posterior_samples[,,(nrow(cdf_1)+2):length(posterior_samples[1,1,])]

# extract dyadic regression slopes
sigma_weight <- posterior_samples[,,1]
intercept <- posterior_samples[,,2]
b_min <- posterior_samples[,,3]
b_max <- posterior_samples[,,4]
posterior_samples <- posterior_samples[,,5:length(posterior_samples[1,1,])]
#b_int <- posterior_samples[,,1]
#posterior_samples <- posterior_samples[,,2:length(posterior_samples[1,1,])]

# extract multimembership samples
mm_matrix <- posterior_samples[,,1:n_dyads]
posterior_samples <- posterior_samples[,,(n_dyads+1):length(posterior_samples[1,1,])]
sigma_mm <- posterior_samples[,,1]
rm(posterior_samples) ; gc()

# save edge samples
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
saveRDS(edges, '../data_processed/anp1_edgedistributions_conditionalprior_dyadreg.RDS')
#edges <- readRDS('../data_processed/anp1_edgedistributions_conditionalprior_dyadreg.RDS')

# save parameter values
extract_slopes <- function(draws, n_samples = 1000, n_chains = 4){
  draws <- as.data.frame(draws) %>% 
    pivot_longer(everything(), names_to = 'chain', values_to = 'slope_draw') %>% 
    separate(chain, sep = c(1,2), remove = T, into = c('chain','.', 'parameter')) %>% 
    select(-.)
  draws$chain <- as.numeric(draws$chain)
  draws$position <- rep(1:n_samples, each = n_chains)
  draws$draw_id <- draws$position + (draws$chain-1)*n_samples
  return(draws)
}
b_min <- extract_slopes(b_min)
b_max <- extract_slopes(b_max)
#b_int <- extract_slopes(b_int)
intercept <- extract_slopes(intercept)
parameters <- rbind(b_min, b_max, #b_int,
                    intercept)

### save data 
saveRDS(parameters, '../data_processed/anp1_dyadicregression_slopeparameters.RDS')

#### check outputs ####
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
plot(NULL, xlim = c(0,1), ylim = c(0,30), las = 1, xlab = 'edge weight', ylab = 'density')
for(i in 1:length(plot_dyads)){
  x <- plot_edges[plot_edges$dyad == plot_dyads[i],]
  lines(density(x$edge_draw), col = ifelse(x$seen_together == 1, rgb(0,0,1,0.1), rgb(1,0,0,0.1)))
}

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

### save image
save.image('anp1_edgeweights_stanmodel.RData')
#load('anp1_edgeweights_stanmodel.RData')

### create mean data frame to plot average values over full distribution
mean_edges <- data.frame(dyad = cdf_1$dyad_id,
                         node_1 = cdf_1$node_1, node_2 = cdf_1$node_2,
                         together = cdf_1$event_count, count_dyad = cdf_1$period_count_dyad,
                         age_1 = cdf_1$age_start_1, age_2 = cdf_1$age_start_2,
                         age_min = cdf_1$age_min, age_max = cdf_1$age_max, age_diff = cdf_1$age_diff,
                         edge_mean = NA)
for(i in 1:nrow(mean_edges)) {
  x <- edges[edges$dyad == mean_edges$dyad[i],]
  mean_edges$edge_mean[i] <- mean(x$edge_draw)
}

### add categories for box plots
mean_edges$age_cat_min <- ifelse(mean_edges$age_min < 15, '10-15',
                                 ifelse(mean_edges$age_min < 20, '15-20',
                                        ifelse(mean_edges$age_min < 25, '20-25',
                                               ifelse(mean_edges$age_min < 40, '25-40', '40+'))))
mean_edges$age_cat_max <- ifelse(mean_edges$age_max < 15, '10-15',
                                 ifelse(mean_edges$age_max < 20, '15-20',
                                        ifelse(mean_edges$age_max < 25, '20-25',
                                               ifelse(mean_edges$age_max < 40, '25-40', '40+'))))

### plot raw data
theme_set(theme_classic())
ggplot()+
  geom_point(data = mean_edges, aes(x = age_min, y = edge_mean),
             shape = 19, colour = rgb(0,0,1,0.1)
  )+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')
ggplot()+
  geom_point(data = mean_edges, aes(x = age_min, y = edge_mean,
                                    colour = age_max), # this makes it look more interesting than it is -- there is a definite pattern at first glance, but actually appears to be nothing more than that the older the youngest male is, the older the oldest must be to be older than the youngest, whereas when the youngest is very young the oldest doesn't have to be very old to still be the oldest
             shape = 19)+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')
ggplot()+
  geom_point(data = mean_edges, aes(x = age_min, y = edge_mean, colour = age_diff), # much less exciting now
             shape = 19)+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')

edges <- left_join(edges, mean_edges, by = 'dyad')
ggplot()+
  geom_point(data = edges, aes(x = age_min, y = edge_draw),
             shape = 19, colour = rgb(0,0,1,0.01))+
  geom_point(data = mean_edges, aes(x = age_min, y = edge_mean),
             shape = 19, colour = 'red')+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')
ggplot()+
  geom_violin(data = edges, aes(x = age_cat_min, y = edge_draw))+
  geom_jitter(data = mean_edges, aes(x = age_cat_min, y = edge_mean, size = count_dyad),
              shape = 1, width = 0.3, colour = rgb(0,0,1,0.1))+
  scale_x_discrete('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')

ggplot()+
  geom_point(data = edges, aes(x = count_dyad, y = edge_draw), colour = rgb(0,0,1,0.01))+
  geom_point(data = mean_edges, aes(x = count_dyad, y = edge_mean), colour = 'red')+
  scale_x_continuous('number of sightings of dyad')+
  scale_y_continuous('edge weight')

ggplot()+
  geom_point(data = mean_edges, aes(x = age_min, y = age_max, colour = edge_mean))+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('age of older dyad member')

# compute posterior means
age_min <- seq(5,46,1)
age_max <- seq(5,46,1)
#mu <- array(NA, dim = c(length(age_max), length(age_min), n_samples), dimnames = list(age_max, age_min, 1:n_samples))
#for(min in 1:nrow(mu)){
#  for(max in 1:ncol(mu)){
#    for(sample in 1:n_samples){
#      if( min <= max ){
#      mu[min,max,sample] <- b_min$slope_draw[sample]*age_min[min] + b_max$slope_draw[sample]*age_max[max] + b_int$slope_draw[sample]*age_min[min]*age_max[max]
#      }
#    }
#  }
#}
mu <- array(NA, dim = c(n_samples, length(age_min), length(age_max)), dimnames = list(1:n_samples, age_min, age_max))
for(sample in 1:n_samples){
  for(min in 1:length(age_min)){
    for(max in 1:length(age_max)){
      if( min <= max ){
        mu[sample,min,max] <- intercept$slope_draw[sample] + b_min$slope_draw[sample]*age_min[min] + b_max$slope_draw[sample]*age_max[max] + b_int$slope_draw[sample]*age_min[min]*age_max[max]
      }
    }
  }
}
mu <- plogis(mu)

mu_max_age <- function(edge_prediction_array, age_threshold){
  minimum_ages <- as.numeric(colnames(edge_prediction_array[1,,]))
  maximum_ages <- as.numeric(names(edge_prediction_array[1,1,]))
  predictions <- edge_prediction_array[,which(minimum_ages <= age_threshold),which(maximum_ages == age_threshold)]
  if(is.null(dim(predictions)) == TRUE) { mean_edge <- mean(predictions) } else { mean_edge <- apply(predictions, 2, mean) }
  return(mean_edge)
}
mu_mean_max20 <- mu_max_age(mu, 20)
mu_mean_max25 <- mu_max_age(mu, 25)
mu_mean_max30 <- mu_max_age(mu, 30)
mu_mean_max35 <- mu_max_age(mu, 35)
mu_mean_max41 <- mu_max_age(mu, 41)
mu_mean_max46 <- mu_max_age(mu, 46)

pi_max_age <- function(edge_prediction_array, age_threshold){
  minimum_ages <- as.numeric(colnames(edge_prediction_array[1,,]))
  maximum_ages <- as.numeric(names(edge_prediction_array[1,1,]))
  predictions <- edge_prediction_array[,which(minimum_ages <= age_threshold),which(maximum_ages == age_threshold)]
  if(is.null(dim(predictions)) == TRUE) { 
    mean_edge <- rethinking::HPDI(predictions) 
    } else { 
    mean_edge <- apply(predictions, 2, rethinking::HPDI) 
    }
  return(mean_edge)
}
pi_mean_max20 <- pi_max_age(mu, 20)
pi_mean_max25 <- pi_max_age(mu, 25)
pi_mean_max30 <- pi_max_age(mu, 30)
pi_mean_max35 <- pi_max_age(mu, 35)
pi_mean_max41 <- pi_max_age(mu, 41)
pi_mean_max46 <- pi_max_age(mu, 46)

# simulate from posterior -- NOT YET WORKED THIS BIT OUT
#sim_edges <- matrix(ncol = length(age_min), nrow = n_samples)
#for(i in 1:nrow(sim_edges)){
#  for(j in 1:ncol(sim_edges)){
#    mean_sim <- intercept$slope_draw[i] + b_min$slope_draw[i]*age_min[j] + b_max$slope_draw[i]*age_max[j] + b_int$slope_draw[i]*age_min[j]*age_max[j]
#    sigma_sim <- ????????
#    sim_edges[i,j] <- rnorm(1, mu = mean_sim, sigma = sigma_sim)
#  }
#}
#hist(sim_edges) # far too many getting a centrality score of 1
#sim_pi <- apply(sim_edges, 2, rethinking::HPDI, prob = 0.95)

# plot
par(mfrow = c(3,2))
ages <- c(20,25,30,35,41,46)
mu_list <- list(mu_mean_max20, mu_mean_max25, mu_mean_max30, mu_mean_max35, mu_mean_max41, mu_mean_max46)
pi_list <- list(pi_mean_max20, pi_mean_max25, pi_mean_max30, pi_mean_max35, pi_mean_max41, pi_mean_max46)
for(i in 1:6){
  plot(edge_mean ~ age_min, data = mean_edges[mean_edges$age_max == ages[i],], col = rgb(0,0,1,0.5), # raw sightings
     pch = 19, las = 1, xlim = c(5,45), ylim = c(0,1),
     xlab = 'age of younger dyad member', ylab = 'edge weight',
     main = paste0('older = ',ages[i],' years old'))
  rethinking::shade(pi_list[[i]], age_min[which(age_min <= ages[i])], lwd = 2, col = rgb(0.5,0,1,0.2))      # add mean line
  lines(age_min[which(age_min <= ages[i])], mu_list[[i]], lwd = 2, col = 'purple')                          # add mean shading
  #rethinking::shade(sim_pi, age_min)                                                                   # add predictions
}
par(mfrow = c(1,1))

### heat map: x = min age, y = max age, colour = edge weight
ggplot()+
  geom_tile(data = mean_edges, mapping = aes(x = age_min, y = age_max, fill = edge_mean))+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('age of older dyad member')

mean_predict <- matrix(NA, nrow = length(age_max), ncol = length(age_min))
for(i in 1:nrow(mean_predict)){
  for(j in 1:ncol(mean_predict)){
    if(i <= j){
      x <- mu_max_age(mu, age_max[j])
      mean_predict[,j] <- c(x ,rep(NA, length(mean_predict[,j]) - length(x)))
    }
  }
}
colnames(mean_predict) <- paste0('max.',age_max)
rownames(mean_predict) <- paste0('min.',age_min)
mean_predict_long <- as.data.frame(mean_predict) %>% 
  pivot_longer(everything(), names_to = 'age_max', values_to = 'mean_edge_weight') %>% 
  separate(age_max, into = c('max.min','age_max'), sep = 4) %>% 
  select(-max.min) %>% 
  mutate(age_min = rep(rownames(mean_predict), each = ncol(mean_predict))) %>% 
  separate(age_min, into = c('max.min','age_min'), sep = 4) %>% 
  select(-max.min) %>% 
  mutate(age_diff = as.numeric(age_max) - as.numeric(age_min))
  #filter(!is.na(mean_edge_weight)) %>% 
mean_predict_long$age_max <- as.numeric(mean_predict_long$age_max)
mean_predict_long$age_min <- as.numeric(mean_predict_long$age_min)
ggplot()+
  geom_tile(data = mean_predict_long, mapping = aes(x = age_min, y = age_max, fill = mean_edge_weight))
ggplot()+
  geom_line(data = mean_predict_long, mapping = aes(x = age_diff, y = mean_edge_weight, colour = age_min, group = age_min))+
  scale_x_continuous(limits = c(0,40))

### end pdf
dev.off()

### clear workspace
rm(cdf_1, n_dyads, counts_ls, n_eles, fit_edges_anp1, edges, plot_dyads, plot_edges, nodes, edge_weights_matrix, edge_samples)

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
