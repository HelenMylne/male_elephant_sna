##### set up ####
# library(tidyverse) ; library(cmdstanr) ; library(brms) ; library(Rcpp) ; library(ggdist) ; library(bayesplot) ; library(igraph) ; library(LaplacesDemon) ; library(bisonR) ; library(janitor) ; library(mice)

library(tidyverse, lib.loc = '../packages/')       # library(tidyverse)
library(cmdstanr, lib.loc = '../packages/')        # library(cmdstanr)
library(brms, lib.loc = '../packages/')            # library(brms)
library(Rcpp, lib.loc = '../packages/')            # library(Rcpp)
library(ggdist, lib.loc = '../packages/')          # library(ggdist)
#library(posterior, lib.loc = '../packages/')       # library(posterior)
library(bayesplot, lib.loc = '../packages/')       # library(bayesplot)
#library(rstan, lib.loc = '../packages/')           # library(rstan)
library(igraph, lib.loc = '../packages/')          # library(igraph)
library(LaplacesDemon, lib.loc = '../packages/')   # library(LaplacesDemon)
library(bisonR, lib.loc = '../packages/')          # library(bisonR)
library(janitor, lib.loc = '../packages/')         # library(janitor)

# load edge weight model and data frames
load('motnp_bisonr_edgescalculated_strongprior.RData')
#rm(counts_df_model, motnp_edges_null_strongpriors) ; gc()

##### read in data ####
df_nodal <- distinct(counts_df[,c('node_1_males','id_1')])
colnames(df_nodal) <- c('node_2_males','id_2')
df_nodal <- rbind(df_nodal, counts_df[nrow(counts_df),c('node_2_males','id_2')])
colnames(df_nodal) <- c('node','id')

motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>% 
  dplyr::select(sort(unique(c(counts_df$id_1, counts_df$id_2)))) %>% 
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')
motnp_ages <- left_join(motnp_ages, df_nodal, by = 'id')
motnp_ages$draw <- rep(1:8000, length(unique(motnp_ages$id)))

#### convert to mids (multiply imputed data set) object ####
ages_test <- motnp_ages[motnp_ages$draw <= 8,]

ages_mice <- mice(motnp_ages)
ages_mice
ages_mice$imp$age
ages_comp <- mice::complete(ages_mice, action = 'long', include = T)

num_eles <- length(unique(df_nodal$id))
num_draws <- 8
num_info <- 3
ages_comp <- array(data = NA, dim = c(num_eles,num_info,num_draws),
                   dimnames = list(df_nodal$id,
                                   colnames(ages_test[c(3,2)]),
                                   1:num_draws))
for(elephant in 1:num_eles){
  ages_elephant <- ages_test[ages_test$node == elephant,c('node','age')]
  for(column in 1:num_info){
    for(draw in 1:num_draws){
      ages_comp[elephant,column,draw] <- as.numeric(ages_elephant[draw,column])
    }
  }
}

num_info <- 4
ages_comp <- array(data = NA, dim = c(num_eles,num_info,num_draws),
                   dimnames = list(df_nodal$id,
                                   colnames(ages_test[c(1,3,2,4)]),
                                   1:num_draws))
for(elephant in 1:num_eles){
  ages_elephant <- ages_test[ages_test$node == elephant,]
  for(draw in 1:num_draws){
    ages_comp[elephant,1,draw] <- as.character(ages_elephant$id[draw])
    ages_comp[elephant,2,draw] <- as.numeric(ages_elephant$node[draw])
    ages_comp[elephant,3,draw] <- as.numeric(ages_elephant$age[draw])
    ages_comp[elephant,4,draw] <- as.numeric(ages_elephant$draw[draw])
  }
}

ages_mids <- as.mids(long = ages_comp, where = which(is.na(motnp_ages$age) == F),
                     .imp = 'draw', .id = 'id')



ages_list <- vector("list", num_draws)
for(i in 1:num_draws){
  ages_list[[i]] <- ages_test[ages_test$draw == i,1:3]
}

ages_mids <- miceadds::datalist2mids(ages_list)


load('motnp_bisonr_edgescalculated_strongprior.RData')
mean_motnp_eigen <- bison_brm(
  formula = bison(node_eigen(node)) ~ age,
  edgemodel_list = motnp_edge_weights_strongpriors,
  data_list = ages_list,
  chains = 4,
  iter = 10000,
  thin = 2
)

# bison_brm -- actually uses bison_mice within itself, not something to necessarily be used before running model. need brms_multiple to use the imputed data
bison_brm <- function(formula, edgemodel_list, data_list, num_draws=100, z_score=FALSE, cores=4, chains=2, ...) {
  # Parse formula
  parsed_formula <- parse_bison_brms_formula(formula)
  brms_formula <- parsed_formula$brms_formula
  
  # Generate mice object
  mice_obj <- bison_mice(
    edgemodel_list,
    data_list,
    parsed_formula$param_names,
    parsed_formula$target_name,
    parsed_formula$metric_name,
    num_draws,
    z_score
  )
  
  # Run brms imputation.
  brms::brm_multiple(brms_formula, mice_obj, backend="cmdstanr", cores=cores, chains=chains, ...)
}












