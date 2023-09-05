#### set up ####
# load packages
#library(cmdstanr, lib.loc = '../packages/')
#library(tidyverse, lib.loc = '../packages/')
#library(dplyr, lib.loc = '../packages/')
#library(igraph, lib.loc = '../packages/')
#library(janitor, lib.loc = '../packages/')
#library(lubridate, lib.loc = '../packages/')
#library(readxl, lib.loc = '../packages/')

library(cmdstanr)
library(tidyverse)
library(dplyr)
library(igraph)
library(janitor)
library(lubridate)
library(readxl)

# set stan path
#set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0')

# load edge weight model
edge_binary_model <- cmdstan_model("models/edge_binary_basic.stan")
#edge_binary <- cmdstan_model("edge_binary_basic.stan")
edge_binary_model

# set seed
set.seed(12345)

#### create data lists ####
cdf_1 <- read_csv('../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_shortwindows_impossiblepairsremoved.csv') %>% 
  filter(period == 1)

### create nodes data frame for period 1
nodes <- data.frame(id = sort(unique(c(cdf_1$id_1,cdf_1$id_2))),          # all unique individuals
                    node = NA, age = NA, sightings = NA)                  # data needed on each
for(i in 1:nrow(nodes)){
  # extract data about individual from cdf_1 data frame
  if(nodes$id[i] %in% cdf_1$id_1) {
    x <- cdf_1[cdf_1$id_1 == nodes$id[i], c('id_1','node_1','count_period_1','age_start_1')] %>% distinct()
  } else { x <- cdf_1[cdf_1$id_2 == nodes$id[i], c('id_2','node_2','count_period_2','age_start_2')] %>% distinct() }
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
n_eles <- length(unique(c(cdf_1$id_1, cdf_1$id_2)))
counts_ls <- list(
  n_dyads    = n_dyads,                  # total number of times one or other of the dyad was observed
  dyad_ids   = cdf_1$dyad_id,            # identifier for each dyad
  together   = cdf_1$event_count,        # count number of sightings seen together
  count_dyad = cdf_1$count_period_dyad)  # count total number of times seen

### Fit model
fit_edges_anp1 <- edge_binary_model$sample(
  data = counts_ls, 
  chains = n_chains, 
  parallel_chains = n_chains)

### save output
save.image('anp1_edges_fit.RData')