#### Information ####
# script takes data input from edge weight estimation for MOTNP population (input = 1000 draws for each of 4 chains per dyad posterior distribution)
# runs through dyadic regression as specified by Jordan Hart in BISoN examples (https://github.com/JHart96/bison_examples/blob/main/examples/dyadic_regression_stan.md)

#### Set up ####
library(tidyverse)
library(rstan)
library(car)
library(bisonR)

#library(tidyverse, lib.loc = 'packages/')
#library(rstan, lib.loc = 'packages/')
#library(car, lib.loc = 'packages/')
#library(bisonR, lib.loc = 'packages/')

load('motnp_bisonr_edgescalculated_strongprior.RData')
rm(random_networks, counts_df_model, gbi_males, m_mat, motnp_edges_null_strongpriors) ; gc()

#### create data frame of male ages ####
ids <- counts_df[,c('id_1','node_1_males')] %>% distinct()
colnames(ids) <- c('id_2','node_2_males')
ids <- rbind(ids, counts_df[nrow(counts_df),c('id_2','node_2_males')])
colnames(ids) <- c('id','node_males')
motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age') %>% 
  filter(id %in% ids$id) %>% 
  left_join(ids, by = 'id')
motnp_ages$draw <- rep(1:8000, each = length(ids$id))

# create data frame of ages and age differences
counts_df_dyadic <- counts_df[,c('dyad_males','node_1_males','node_2_males')]
colnames(motnp_ages)[3] <- 'node_1_males'
counts_df_dyadic <- left_join(counts_df_dyadic, motnp_ages, by = 'node_1_males', multiple = 'all')
colnames(counts_df_dyadic)[4:5] <- c('id_1','age_1')
colnames(motnp_ages)[3] <- 'node_2_males'
counts_df_dyadic <- left_join(counts_df_dyadic, motnp_ages, by = c('node_2_males','draw'))
colnames(counts_df_dyadic)[7:8] <- c('id_2','age_2')
counts_df_dyadic <- counts_df_dyadic[,c(1:3,6,4:5,7:8)]
colnames(counts_df_dyadic)[2:3] <- c('node_1_id','node_2_id')
head(counts_df_dyadic)

#counts_df_dyadic$age_diff <- counts_df_dyadic$age_1 - counts_df_dyadic$age_2
#counts_df_dyadic$age_mean <- mean(counts_df_dyadic$age_1, counts_df_dyadic$age_2)

#counts_df_dyadic <- counts_df_dyadic[,c('node_1_males','node_2_males','age_diff','age_mean')]

cdf_dyadic <- counts_df_dyadic[,c('dyad_males','node_1_id','node_2_id')] %>% distinct()
cdf_dyadic$age_1 <- NA ; cdf_dyadic$age_2 <- NA
for(i in 1:nrow(cdf_dyadic)){
  dyad <- counts_df_dyadic[counts_df_dyadic$dyad_males == cdf_dyadic$dyad_males[i],]
  cdf_dyadic$age_1 <- mean(dyad$age_1)
  cdf_dyadic$age_2 <- mean(dyad$age_2)
}

#### fit model to mean age ####
mean_age_dyadic <- bison_brm (
  bison(edge_weight(node_1_id, node_2_id)) ~ age_1 + age_2 + (1 | mm(node_1_id, node_2_id)),
  motnp_edge_weights_strongpriors,
  cdf_dyadic,
  num_draws = 5, # Small sample size for demonstration purposes
  refresh = 0
)
summary(mean_age_dyadic)


