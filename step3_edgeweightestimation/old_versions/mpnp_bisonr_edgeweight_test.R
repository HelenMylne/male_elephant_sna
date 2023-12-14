#### Bayesian analysis of EfA data ####
# Script to process association data from Makgadikgadi Pans National Park, Botswana
# Data collected: 2012-2021 by Elephants for Africa
# Data supplied by: Kate Evans

#### Set up ####
# load packages
library(cmdstanr, lib.loc = '../packages/')    # library(cmdstanr)
library(tidyverse, lib.loc = '../packages/')   # library(tidyverse)
library(dplyr, lib.loc = '../packages/')       # library(dplyr)
#library(rstan, lib.loc = '../packages/')      # library(rstan)
library(bisonR, lib.loc = '../packages/')      # library(bisonR)
library(igraph, lib.loc = '../packages/')     # library(igraph)
library(asnipe, lib.loc = '../packages/')      # library(asnipe)
library(sna, lib.loc = '../packages/')         # library(sna)
#library(raster, lib.loc = '../packages/')      # library(raster)

# information
sessionInfo()
R.Version()

# set seed
set.seed(12345)

## set up ####
### add time marker
print(paste0('start window 5 at ', Sys.time()))

### create output pdf
pdf(file = paste0('../outputs/mpnpshort5_testrun.pdf'))

### read in data frame for edge weight model
filename <- paste0('../data_processed/mpnp_period5_pairwiseevents.csv')
counts_df <- read_csv(filename)

### read in ages for filtering
filename <- paste0('../data_processed/mpnp5_ageestimates_mcmcoutput.rds')
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

### create priors
priors <- get_default_priors('binary')
priors$edge <- 'normal(-2.5, 1.5)'    # slightly right skewed so that more density is towards 0 but still very weak and all values are possible
priors$fixed <- 'normal(0, 1.5)'  # slightly less wide that before, but zero centred so allows age to have a positive or negative impact on centrality and edge weight
prior_check(priors, 'binary')

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

## clean up ####
### save pdf
dev.off()

# save workspace image
save.image(paste0('mpnp_edgecalculations/mpnpshort5_testrun.RData'))

### clear environment
rm(list = ls()[!(ls() %in% c('time_window','priors',
                             'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])

### add time marker
print(paste0('time window ', time_window, ' completed at ', Sys.time()))
