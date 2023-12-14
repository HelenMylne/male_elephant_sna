library(tidyverse)

sample_size <- data.frame(window = c('short1','short2','short3','short4','short5','long'),
                          count_nodes = NA,
                          count_dyads = NA)
for(time_window in 1:nrow(sample_size)){
  if(sample_size$window[time_window] == 'long'){
    df <- read_csv('../data_processed/step1_dataprocessing/mpnp_longtimewindow_pairwiseevents.csv') 
  } else {
    df <- read_csv(paste0('../data_processed/step1_dataprocessing/mpnp_period',time_window,'_pairwiseevents.csv'))
  }
  sample_size$count_nodes[time_window] <- length(unique(c(df$id_1, df$id_2)))
  sample_size$count_dyads[time_window] <- nrow(df)
}

sample_size



