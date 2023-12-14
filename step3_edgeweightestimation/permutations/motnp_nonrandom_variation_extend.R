#### asnipe non-random variation test ####
# library(tidyverse) ; library(plyr) ; library(aninet) ; library(rstan) ; library(asnipe)
library(cmdstanr, lib.loc = '../packages/')
library(cli, lib.loc = '../packages/')
library(rstan, lib.loc = '../packages/')
library(plyr, lib.loc = '../packages/')
library(tidyverse, lib.loc = '../packages/')
library(asnipe, lib.loc = '../packages/')
#library(aninet, lib.loc = '../packages/')
library(data.table, lib.loc = '../packages/')
library(spatsoc, lib.loc = '../packages/')
library(raster, lib.loc = '../packages/')
library(bisonR, lib.loc = '../packages/')
library(terra, lib.loc = '../packages/')

load('motnp_permutations.RData')

current_seed <- seed
for(seed in asnipe_seeds[current_seed:length(asnipe_seeds)]){
  permutation <- chain_length*seed
  set.seed(seed)
  ### run network permutations
  random_networks <- asnipe::network_permutation(association_data = gbi_males,   # permute network -- raw data = gbi_matrix
                                                 association_matrix = mat_new,   # SRI matrix
                                                 permutations = chain_length,    # 100 times
                                                 days = sightings$date,          # permute within day only
                                                 within_day = TRUE)              # permute within day only
  
  ### calculate coefficient of variation for final randomised network
  cv_random_networks <- c(cv_random_networks, rep(0,chain_length))               # join onto previous outputs
  for (i in c(1:chain_length)) {                                                 # set up for loop
    net_rand <- random_networks[i,,]                                             # select random network
    cv_random_networks[(permutation+i)-chain_length] <- raster::cv(net_rand)     # calculate CV and save into vector
  }
  
  # extract final network for continuation of chain
  mat_new <- random_networks[chain_length,,]
  
  # fill in outputs data frame with information on 1000th permutation
  permute_num <- which(outputs$permutation == permutation)
  outputs$cv[permute_num]    <- raster::cv(mat_new)
  outputs$mean[permute_num]  <- mean(mat_new)
  outputs$stdev[permute_num] <- sd(mat_new)
  outputs$rhat[permute_num]  <- rhat(mat_new)
  outputs$p_value[permute_num] <- length(which(outputs$cv[permute_num] >= cv_random_networks))/outputs$permutation[permute_num]
  
 ### save workspace every 10000th network
  if(seed %% 250 == 0) {
    saveRDS(cv_random_networks, '../data_processed/motnp_permutations_cv.RDS')
    saveRDS(outputs, '../data_processed/motnp_permutations_allmeasures.RDS')
    save.image('motnp_permutations.RData')
    print(paste0(permutation,' permutations done'))
    }
}

hist(outputs$p_value, breaks = 100)
abline(v = 0.05, col = 'red')

hist(outputs$cv, breaks = 100)
abline(v = outputs$cv[1], col = 'red')

outputs2 <- outputs[!is.na(outputs$cv),]
outputs2$p_value[nrow(outputs2)]
rm(outputs2) ; gc()
