## script to extract calculations of population size per time window for ANP data
library(tidyverse, lib.loc = '../packages/')

## import data
anp_all <- read_csv('../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_shortwindows_impossiblepairsremoved.csv')

## create empty data frame
population_size <- data.frame(window = unique(anp_all$period),
                              population = NA)

## fill data frame
for(window in 1:length(unique(anp_all$period))){
  load(paste0('anp_edgecalculations/anpshort',window,'_bisonr_edgescalculated.RData')) # load data frame after model run
  population_size$population[window] <- length(unique(counts_df_model$node_1_id)) + 1  # population size = number of elephants in node_1 position plus one elephant always classed as node_2
}

## check output
population_size

## save output
write_csv(population_size, '../data_processed/anp_population_per_timewindow.csv')
