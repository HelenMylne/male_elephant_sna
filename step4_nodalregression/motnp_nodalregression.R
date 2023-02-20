#### information #####
# 1) to estimate the ages of individuals in other populations, first came up with a survival curve for the Amboseli elephants which have good birth and death data.
# 2) selected the Gompertz bathtub distribution
# 3) use distribution to convert MOTNP ages in categories to probability distributions of true age (first part of this script)
# 4) extract node centrality distribution from draws of edge weight
# 5) use MOTNP age distributions to predict network centrality (third part of this script)

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

# define PDF output
pdf('../outputs/motnp_nodalregression_plots.pdf')

# load edge weight model and data frames
load('motnp_bisonr_edgescalculated_strongprior.RData')

#### nodal regression -- mean age estimate only, not full distribution ####
df_nodal <- distinct(counts_df[,c('node_1_males','id_1')])
colnames(df_nodal) <- c('node_2_males','id_2')
df_nodal <- rbind(df_nodal, counts_df[nrow(counts_df),c('node_2_males','id_2')])
colnames(df_nodal) <- c('node','id')

motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>% 
  select(sort(unique(c(counts_df$id_1, counts_df$id_2)))) %>% 
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')
motnp_ages <- left_join(motnp_ages, df_nodal, by = 'id')

mean_motnp_ages <- df_nodal
mean_motnp_ages$age <- NA
for(i in 1:nrow(mean_motnp_ages)){
  x <- motnp_ages[motnp_ages$id == mean_motnp_ages$id[i],]
  mean_motnp_ages$age[i] <- mean(x$age)
  rm(x)
}

rm(list= ls()[!(ls() %in% c('mean_motnp_ages','motnp_edge_weights_strongpriors','motnp_ages')) ]) ; gc()

mean_motnp_nodal <- bison_brm(
  bison(node_eigen(node)) ~ age,
  motnp_edge_weights_strongpriors,
  mean_motnp_ages,
  chains = 4,
  iter = 100000,
  thin = 2
)

summary(mean_motnp_nodal) # FIT 100 IMPUTED MODELS, 4 CHAINS PER MODEL, EACH 1000 DRAWS LONG (+1000 DRAWS WARM UP). WARNING AT END OF MODEL RUN THAT CHAINS <3 DRAWS LONG AS ACTUAL CHAIN IS ONLY 1 WARMUP AND 1 SAMPLE. ONLY IMPUTED CHAINS FOLLOW THE SPECIFIED ITERATIONS AND THINNING

mean_eigen_values <- mean_motnp_nodal$data
plot(mean_eigen_values$bison_node_eigen ~ mean_eigen_values$age, 
     las = 1, pch = 19, col = rgb(0,0,1,0.2),
     xlab = 'mean age estimate', ylab = 'eigenvector centrality',
     main = 'no effect of age on eigenvector centrality')

mean_mod_summary <- mean_motnp_nodal$fit

hist(mean_motnp_nodal$rhats[,2], las = 1, main = 'Rhat values for 100 imputed model runs', xlab = 'Rhat')

post <- as_draws_df(mean_motnp_nodal)
hist(post$b_age) # what scale is this on? should it go through plogis() (aka inverse logit) or not -- what scale does eigenvector read as?
plot(data = post[post$.chain == 1,], b_age ~ .draw, type = 'l')

# save workspace image for reloading at a later date that doesn't require running model again
save.image('motnp_bisonr_nodalregression_meanage.RData')

# remove nodal mean age
rm(mean_motnp_nodal, mean_eigen_values, mean_mod_summary) ; gc()

#### nodal regression -- age distribution ####
motnp_nodal <- bison_brm(
  bison(node_eigen(node)) ~ age,
  motnp_edge_weights_strongpriors,
  motnp_ages,
  chains = 4
)
summary(motnp_nodal)

eigen_values <- motnp_nodal$data
plot(eigen_values$bison_node_eigen ~ eigen_values$age, las = 1, pch = 19, col = rgb(0,0,1,0.2),
     xlab = 'mean age estimate', ylab = 'eigenvector centrality',
     main = 'no effect of age on eigenvector centrality')

mod_summary <- motnp_nodal$fit

hist(motnp_nodal$rhats[,2], las = 1, main = 'Rhat values for 100 imputed model runs', xlab = 'Rhat')

save.image('motnp_bisonr_nodalregressionrun.RData')

#### compare to null model ####
motnp_nodal_null <- bison_brm(
  bison(node_eigen(node)) ~ 1,
  motnp_edge_weights,
  motnp_ages,
  chains = 4
)

model_comparison(list(non_random_model = motnp_nodal, random_model = motnp_nodal_null))

#### clean up ####
dev.off()