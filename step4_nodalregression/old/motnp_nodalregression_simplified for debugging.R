#### set up ####
options(future.globals.maxSize = 10000*(1024^2))   # running model with full age distributions = 5.01GB of globals, which exceeds default maximum allowed size. Set to allow up to 10 GB to allow model to run

library(tidyverse, lib.loc = '../packages/')       # library(tidyverse)
library(cmdstanr, lib.loc = '../packages/')        # library(cmdstanr)
library(brms, lib.loc = '../packages/')            # library(brms)
library(Rcpp, lib.loc = '../packages/')            # library(Rcpp)
library(ggdist, lib.loc = '../packages/')          # library(ggdist)
library(posterior, lib.loc = '../packages/')       # library(posterior)
library(bayesplot, lib.loc = '../packages/')       # library(bayesplot)
#library(rstan, lib.loc = '../packages/')           # library(rstan)
library(igraph, lib.loc = '../packages/')          # library(igraph)
library(LaplacesDemon, lib.loc = '../packages/')   # library(LaplacesDemon)
library(bisonR, lib.loc = '../packages/')          # library(bisonR)

# load edge weight model and data frames
load('motnp_bisonr_edgescalculated_strongprior.RData') # currently all saved versions of this use prior N(0,1.5) for fixed, but running regression at the moment with normal(0,1) instead

#### nodal data frame ####
df_nodal <- distinct(counts_df[,c('node_1_males','id_1')])
colnames(df_nodal) <- c('node_2_males','id_2')
df_nodal <- rbind(df_nodal, counts_df[nrow(counts_df),c('node_2_males','id_2')])
colnames(df_nodal) <- c('node','id')

motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>% 
  select(sort(unique(c(counts_df$id_1, counts_df$id_2)))) %>% 
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')
motnp_ages <- left_join(motnp_ages, df_nodal, by = 'id')
motnp_ages$draw <- rep(1:8000, length(unique(motnp_ages$id)))

## eigenvector only ####
pdf('../outputs/motnp_nodalregression_eigen_agedistribution.pdf')

## reduce age data to manageable number of draws
test_ages <- motnp_ages[motnp_ages$draw %in% sample(motnp_ages$draw, 2000, replace = F),]

## run model
motnp_eigen <- bison_brm(
  bison(node_eigen(node)) ~ age,
  motnp_edge_weights_strongpriors,
  #motnp_ages, # actual model was planning to use all 8000 draws from posterior
  test_ages,   # 2000 draws instead of 8000 as taking too long
  chains = 4,
  iter = 1000, # previously running for 10000 but was taking too long
  cores = 4
)
summary(motnp_eigen)                  # rhat values very large
hist(motnp_eigen$rhats[,2], las = 1,  # rhat values look fine here
     main = 'Rhat values for 100 imputed model runs', xlab = 'Rhat')

eigen_values <- motnp_eigen$data
plot(eigen_values$bison_node_eigen ~ eigen_values$age, # basic plot of model data
     las = 1, pch = 19, col = rgb(0,0,1,0.2),
     xlab = 'mean age estimate', ylab = 'eigenvector centrality',
     main = 'effect of age on eigenvector centrality')

# compare to null model
motnp_eigen_null <- bison_brm(
  bison(node_eigen(node)) ~ 1,
  motnp_edge_weights,
  #motnp_ages,
  test_ages,
  chains = 4,
  cores = 4
)
model_comparison(list(non_random_model = motnp_eigen, random_model = motnp_eigen_null))

save.image('motnp_nodalregression_eigenvector.RData')

rm(motnp_eigen, eigen_values, mod_summary, motnp_eigen_null) ; gc()
dev.off()
