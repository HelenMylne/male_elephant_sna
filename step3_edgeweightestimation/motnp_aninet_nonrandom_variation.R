#### aninet non-random variation test ####
library(tidyverse, lib.loc = '../packages/')
library(aninet, lib.loc = '../packages/')
library(rstan, lib.loc = '../packages/')

#### generate gbi_matrix ####
load('motnp_bisonr_edgescalculated_strongprior.RData')
rm(list = ls()[ls() != 'counts_df'])
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
gbi_males <- gbi_matrix[,colnames(gbi_matrix) %in% sort(unique(c(counts_df$id_1, counts_df$id_2)))]
gbi_males <- gbi_males[rowSums(gbi_males) > 0,] # remove male-only sightings
rm(eles, eles_asnipe, nodes) ; gc()

( num_permute_hypo <- aninet::gbi_MCMC_iters(gbi = gbi_males, target_samples = 1000, quiet = FALSE) )    # hypothetical minimum number of permutations necessary to reach stability
num_permute_real <- 1000 # test how long this takes, but can't actually run parallel because needs to be as close as possible to a single long chain -- measuring burn in, not sampling from chain

#permute_gbi <- gbi_MCMC(
#  data = gbi_males,
#  ind_constraint = NULL,
#  group_constraint = NULL,
#  samples = 10,#num_permutations/num_chains,
#  thin = 2, #100,
#  burnin = 1, #1000,
#  chains = num_chains,
#  FUN = NULL      # leaving as NULL here will mean that permute_gbi$FUN(gbi_males) below gives CV, Mean, SD and SRI
#)
seeds <- 1:num_permute_real
outputs <- data.frame(permutation = 1:num_permute_real,
                      cv = NA, mean = NA, stdev = NA, sri = NA, rhat = NA)
N <- ncol(gbi_males)
obs <- nrow(gbi_males)
num_chains <- 2
permute_gbi <- array(NA, c(N,obs,num_permute_real,num_chains),
                     dimnames = list(colnames(gbi_males),
                                     rownames(gbi_males),
                                     1:num_permute_real,
                                     1:num_chains))

rhat_rfun <- function(sims) {
  chains <- ncol(sims)
  n_samples <- nrow(sims)
  chain_mean <- numeric(chains)
  chain_var <- numeric(chains)
  for (i in seq_len(chains)) {
    chain_mean[i] <- mean(sims[, i])
    chain_var[i] <- var(sims[, i])
  }
  var_between <- n_samples * var(chain_mean)
  var_within <- mean(chain_var)
  sqrt((var_between / var_within + n_samples - 1) / n_samples)
}

z_scale <- function(x) {
  S <- length(x)
  r <- rank(x, ties.method = 'average')
  z <- qnorm((r - 1 / 2) / S)
  z[is.na(x)] <- NA
  if (!is.null(dim(x))) {
    # output should have the input dimension
    z <- array(z, dim = dim(x), dimnames = dimnames(x))
  }
  z
}

output_function <- function(gbi) {
  x <- get_numerator(gbi, data_format = "GBI", return = "vector")
  d <- get_denominator(gbi, data_format = "GBI", return = "vector")
  sri <- x/d
  bulk_rhat <- rhat_rfun(z_scale(gbi))
  gbi_folded <- abs(gbi - median(gbi))
  tail_rhat <- rhat_rfun(z_scale(gbi_folded))
  res <- c(mean(sri, na.rm = T),
           stats::sd(sri, na.rm = T), 
           stats::sd(sri, na.rm = T)/mean(sri, na.rm = T), 
           mean(sri > 0, na.rm = T),
           max(bulk_rhat, tail_rhat))
  names(res) <- c("Mean", "SD", "CV", "Non-zero", 'Rhat')
  return(res)
}

for(seed in seeds){
  set.seed(seed)
  print(paste0('start permutation number ', seed, ' at ', Sys.time()))
  permute <- aninet::gbi_MCMC(data = gbi_males,
                              ind_constraint = NULL, group_constraint = NULL,
                              samples = 1000, # short chains that can repeatedly extended and see how long they take -- previously: ceiling(num_permutations/num_chains),
                              thin = 100, burnin = 1000,
                              chains = num_chains,
                              FUN = output_function # NULL = CV, Mean, SD and SRI in permute_gbi$FUN(gbi_males)
  )
  print(paste0('initial 1000 samples drawn for permutation ', seed, ' at ', Sys.time()))
  permute <- aninet::extend.gbi_MCMC(permute, samples = 1000)
  print(paste0('extended 1000 samples drawn for permutation ', seed, ' at ', Sys.time()))
  outputs[which(outputs$permutation == seed),2:6] <- permute$FUN(gbi_males)
  permute_gbi[,,seed,1] <- permute$permuted_data[[1]]
  permute_gbi[,,seed,2] <- permute$permuted_data[[2]]
}

write_csv(outputs, '../data_processed/motnp_simulation_cv_rhat.csv')
saveRDS(permute_gbi, '../data_processed/motnp_simulation_permutations.RDS')

# make parallel
#library(doParallel)
#num_chains <- 1000 # test how long this takes, but can't actually run parallel because needs to be as close as possible to a single long chain -- measuring burn in, not sampling from chain

#cl <- makeCluster(detectCores())
#registerDoParallel(cl)

#seeds <- 1:num_chains
#library(aninet, lib.loc = '../packages/')
#permute_gbi <- foreach(seed = seeds) %dopar% { # , .packages = 'aninet'
#  set.seed(seed)
#  print(paste0('start permutation number ', seed, ' at ', Sys.time()))
#  permute <- aninet::gbi_MCMC(data = gbi_males,
#                              ind_constraint = NULL, group_constraint = NULL,
#                              samples = 1000, # short chains that can repeatedly extended and see how long they take -- previously: ceiling(num_permutations/num_chains),
#                              thin = 100, burnin = 100, # previously burn in of 1000 but then that wouldn't register anything
#                              chains = num_chains,
#                              FUN = NULL # NULL = CV, Mean, SD and SRI in permute_gbi$FUN(gbi_males)
#  )
#  print(paste0('initial 1000 samples drawn for chain ', seed, ' at ', Sys.time()))
#  permute <- aninet::extend.gbi_MCMC(permute, samples = 1000)
#  print(paste0('extended 1000 samples drawn for chain ', seed, ' at ', Sys.time()))
#}
#str(permute_gbi)
#permute_gbi$FUN(gbi_males)
#permute_gbi$permuted_data

#save.image('motnp_nonrandomtest_gbiMCMCpermutations.RData')

