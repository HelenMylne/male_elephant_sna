## Information ##
# script to check for if there is evidence in the ANP data for males associating repeatedly over many years
library(tidyverse)

load('step5_dyadicregression/anpshort_dyadicregression_dataimported.RData')

head(edge_summary)

node_data <- node_data %>% 
  rename(node_1 = node) %>% 
  left_join(nodes_random[,c('node_1','node_1_random')], by = 'node_1') %>% 
  rename(node = node_1, node_random = node_1_random)

rm(nodes_random, edge_weights, edges, edge_samples, n_chains, n_samples, time_window) ; gc()

check <- data.frame(dyad = unique(edge_summary$dyad),
                    dyad_id = unique(edge_summary$dyad_id)) %>% 
  separate(dyad, into = c('node_1','node_2'), sep = '_', remove = F) %>% 
  mutate(w1 = NA, w2 = NA, w3 = NA, w4 = NA, w5 = NA, w6 = NA, w7 = NA, w8 = NA, w9 = NA,
         w10 = NA, w11 = NA, w12 = NA, w13 = NA, w14 = NA, w15 = NA, w16 = NA, w17 = NA, w18 = NA,
         w19 = NA, w20 = NA, w21 = NA, w22 = NA, w23 = NA, w24 = NA, w25 = NA, w26 = NA, w27 = NA,
         w28 = NA, w29 = NA, w30 = NA, w31 = NA, w32 = NA, w33 = NA, w34 = NA, w35 = NA, w36 = NA)
check <- as.matrix(check)
for(j in 1:n_windows){
  window <- edge_summary %>% 
    dplyr::select(dyad, dyad_id, node_1, node_2, period, period_count_1, period_count_2, period_count_dyad, event_count) %>% 
    filter(period == j)
  for(i in 1:nrow(check)){
    if(check[i,1] %in% window$dyad){
      check[i,j+4] <- ifelse(window$event_count[window$dyad == check[i,1]] > 0, 1, 0)
    }
  }
  print(paste0('window ',j,' completed at ',Sys.time()))
}

save.image('check_repeat_interactions_anp.RData')

check2 <- as.data.frame(check)

table(check2)







