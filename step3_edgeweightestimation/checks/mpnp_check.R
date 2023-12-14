library(tidyverse)
counts_df <- read_csv('../data_processed/mpnp_longtimewindow_pairwiseevents.csv') %>%
  distinct()
length(unique(counts_df$dyad))
check <- as.data.frame(table(counts_df$dyad))
check[check$Freq > 1,]
unique(check$Freq)
check <- check[check$Freq > 1,]
check2 <- counts_df[counts_df$dyad %in% check$Var1,]
rm(eles, mean_ages, missing_age, nodes, x, i, edge_weights_matrix) ; gc()

