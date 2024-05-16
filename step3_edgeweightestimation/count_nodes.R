library(tidyverse)

nodes <- read_csv('../data_processed/step4_nodalregression/anp_allnodes.csv')

counts <- data.frame(window = sort(unique(nodes$window)),
                     count = NA)
for(i in 1:nrow(counts)){
  x <- nodes %>% 
    filter(window == counts$window[i])
  counts$count[i] <- length(unique(x$id))
}

rm(list = ls()) ; gc()

counts <- data.frame(window = 1:5,
                     count = NA)

for(time_window in 1:5){
  ages <- readRDS(paste0('../data_processed/step2_ageestimation/mpnp',
                         time_window,'_ageestimates_mcmcoutput.rds'))
  counts$count[time_window] <- ncol(ages)
}

rm(list = ls()) ; gc()

counts <- data.frame(window = 1:7,
                     count = NA)

for(time_window in 1:7){
  eigen <- readRDS(paste0('../data_processed/step4_nodalregression/anplong',
                         time_window,'eigenvectorestimates.rds'))
  counts$count[time_window] <- length(unique(eigen$node))
}

