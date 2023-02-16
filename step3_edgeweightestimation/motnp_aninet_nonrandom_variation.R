#### aninet non-random variation test ####
library(tidyverse, lib.loc = '../packages/')

#### generate gbi_matrix ####
eles <- read_csv('../data_processed/motnp_eles_long.csv')
eles$location <- paste(eles$gps_s, eles$gps_e, sep = '_')       # make single variable for unique locations
eles <- eles[,c(1,16,2,3,17,4,5,14,7,8,10,13)]                  # rearrange variables
nodes <- read_csv('../data_processed/motnp_elenodes.csv')       # read in node data
eles_asnipe <- eles[,c(3,4,2,5)]                                # date, time, elephant, location
eles_asnipe$Date <- as.integer(eles_asnipe$date)                # make date numeric
eles_asnipe$Date <- 1+eles_asnipe$Date - min(eles_asnipe$Date)  # start from 1, not 1st January 1970
eles_asnipe$Time <- ifelse(eles_asnipe$time > 1, NA, eles_asnipe$time) # time = proportion of day so anything >1 has to be wrong
eles_asnipe$Time <- eles_asnipe$Time*(24*60*60)                        # convert time values to seconds through day
eles_asnipe <- eles_asnipe[,c(5,6,3,4)]                                # create data frame to produce gbi matrix from
colnames(eles_asnipe) <- c('Date','Time','ID','Location')              # rename variables for get_gbi
eles_asnipe$ID <- as.character(eles_asnipe$ID)                         # correct data type
eles_asnipe$d_pad <- str_pad(eles_asnipe$Date, 3, pad = '0')           # 0-pad dates
eles_asnipe$encounter <- paste(eles_asnipe$d_pad, eles_asnipe$Time, eles_asnipe$Location, sep = '_') # unique value for each sighting
eles_asnipe$group <- as.integer(as.factor(eles_asnipe$encounter)) # unique factor for every sighting
eles_asnipe <- eles_asnipe[,c(3,7)]            # create data table for gbi matrix
eles_asnipe <- data.table::setDT(eles_asnipe)  # create data table for gbi matrix
gbi_matrix <- spatsoc::get_gbi(DT = eles_asnipe, group = 'group', id = 'ID')  # create gbi matrix
motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds')
ids <- colnames(motnp_ages)
motnp_ages <- pivot_longer(motnp_ages, cols = everything(), names_to = 'id', values_to = 'age')
gbi_males <- gbi_matrix[,colnames(gbi_matrix) %in% ids]
gbi_males <- gbi_males[rowSums(gbi_males) > 0,] # remove male-only sightings
rm(eles, eles_asnipe, motnp_ages, nodes, ids) ; gc()

( num_permutations <- aninet::gbi_MCMC_iters(gbi = gbi_males, target_samples = 1000, quiet = FALSE) )
num_permutations
#num_chains <- 4

#permute_gbi <- gbi_MCMC(
#  data = gbi_males,
#  ind_constraint = NULL,
#  group_constraint = NULL,
#  samples = num_permutations/num_chains,
#  thin = 100,
#  burnin = 1000,
#  chains = num_chains,
#  FUN = NULL      # leaving as NULL here will mean that permute_gbi$FUN(gbi_males) below gives CV, Mean, SD and SRI
#)
#str(permute_gbi)
#permute_gbi$FUN(gbi_males)
#permute_gbi$permuted_data

# make parallel
library(doParallel)
(num_chains <- detectCores())

cl <- makeCluster(num_chains)
registerDoParallel(cl)

seeds <- 1:num_chains
library(aninet, lib.loc = '../packages/')
permute_gbi <- foreach(seeds = seeds) %dopar% { # , .packages = 'aninet'
  set.seed(seeds)
  aninet::gbi_MCMC(data = gbi_males,
                   ind_constraint = NULL, group_constraint = NULL,
                   samples = ceiling(num_permutations/num_chains),
                   thin = 100, burnin = 1000,
                   chains = num_chains,
                   FUN = NULL      # NULL = CV, Mean, SD and SRI in permute_gbi$FUN(gbi_males)
  )
}
str(permute_gbi)
permute_gbi$FUN(gbi_males)
permute_gbi$permuted_data

save.image('motnp_nonrandomtest_gbiMCMCpermutations.RData')

#### ignore all this -- I thought the permutations needed to be on SRI data but nope it just wants the GBI matrix ####
#counts_df$sri <- counts_df$event_count / counts_df$count_dyad
#sri_males <- matrix(data = NA, nrow = length(ids), ncol = length(ids), dimnames = list(x = ids, y = ids))
#for(i in 1:length(ids)){
#  for(j in 1:length(ids)){
#    if( i == j ) { sri_males[i,j] <- 0 }
#    else {
#      if( i < j ){
#        sri_males[i,j] <- counts_df$sri[which(counts_df$id_1 == rownames(sri_males)[i] &
#                                                counts_df$id_2 == colnames(sri_males)[j])]
#      }
#      else {
#        sri_males[i,j] <- counts_df$sri[which(counts_df$id_1 == rownames(sri_males)[j] &
#                                                counts_df$id_2 == colnames(sri_males)[i])]
#      }
#    }
#  }
#}

#( num_permutations <- aninet::gbi_MCMC_iters(gbi = sri_males, target_samples = 1000, quiet = FALSE) )

#permute_network <- gbi_MCMC(
#  data = sri_males,
#  ind_constraint = NULL,
#  group_constraint = NULL,
#  samples = num_permutations,
#  thin = 100,
#  burnin = 1000,
#  chains = 4,
#  FUN = NULL      # leaving as NULL here will mean that permute_network$FUN(gbi_males) below gives CV, Mean, SD and SRI
#)
#str(permute_network)
#permute_network$FUN(gbi_males)
#permute_network$permuted_data

