library(tidyverse, lib.loc = '../packages/')

anp_all <- read_csv('../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_shortwindows_impossiblepairsremoved.csv')

population_size <- data.frame(window = unique(anp_all$period),
                              population = NA)

for(window in 1:length(unique(anp_all$period))){
  load(paste0('anp_edgecalculations/anpshort',window,'_bisonr_edgescalculated.RData'))
  population_size$population[window] <- length(unique(counts_df_model$node_1_id)) + 1
}

population_size