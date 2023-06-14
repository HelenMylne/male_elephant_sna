#### asnipe non-random variation test ####
# library(tidyverse) ; library(plyr) ; library(aninet) ; library(rstan) ; library(asnipe)
library(cli, lib.loc = '../packages/')
library(cmdstanr, lib.loc = '../packages/')
library(tidyverse, lib.loc = '../packages/')
library(rstan, lib.loc = '../packages/')
library(plyr, lib.loc = '../packages/')
library(asnipe, lib.loc = '../packages/')
#library(aninet, lib.loc = '../packages/')     # can't install runjags on Viking so can't install aninet
library(data.table, lib.loc = '../packages/')
library(spatsoc, lib.loc = '../packages/')
library(raster, lib.loc = '../packages/')
library(bisonR, lib.loc = '../packages/')
library(terra, lib.loc = '../packages/')

load('motnp_bisonr_edgescalculated_strongprior.RData')
rm(list = ls()[! ls() %in% c('motnp_edge_weights_strongpriors', 'counts_df')]) ; gc()

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
gbi_males <- gbi_matrix[,colnames(gbi_matrix) %in% sort(unique(c(counts_df$id_1, counts_df$id_2)))]
gbi_males <- gbi_males[rowSums(gbi_males) > 0,] # remove female-only sightings

#### calculate number of permtuations to run ####
#( num_permute <- aninet::gbi_MCMC_iters(gbi = gbi_males, target_samples = 1000, quiet = FALSE) )
num_permute <- 455216741
num_permute <- round_any(num_permute, 1000, f = ceiling)

#### functions to print out Rhat at the end of the permutation ####
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

get_numerator <- function(data, data_format = "SP", return = "vector"){                    # lifted from aninet to make it work
  if(data_format == "SP"){
    X <- apply(data,c(2,3),sum)
  }
  if(data_format == "GBI"){
    X <- t(data) %*% data
    diag(X) <- 0
  }
  if(return == "vector" | is.null(return)){
    return(X[lower.tri(X)])
  }
  if(return == "matrix"){
    diag(X) = 0
    colnames(X) = colnames(data)
    row.names(X) = colnames(data)
    return(X)
  }
}

get_denominator <- function(data, index = "SRI", data_format = "SP", return = "vector"){   # lifted from aninet to make it work
  if(!index %in% c("SRI", "HWI", "BII")) stop("Invalid Association Index")
  if(data_format == "SP"){
    occur <- ifelse(apply(data,c(1,2),sum) > 0, 1, 0)
    if(index == "BII") denominator <- t(occur) %*% occur
    if(index == "SRI") denominator <- nrow(occur) - t(1-occur) %*% (1-occur)
    if(index == "HWI"){
      ya <- t(occur) %*% (1-occur)
      yb <- t(ya)
      yab <- t(occur) %*% (occur)
      denominator <- 0.5*(ya+yb) + yab
    }
  }
  if(data_format =="GBI"){
    if(index == "BII") denominator <- t(data) %*% data
    if(index == "SRI") denominator <- nrow(data) - t(1-data) %*% (1- data)
    if(index == "HWI"){
      ya <- t(data) %*% (1-data)
      yb <- t(ya)
      yab <- t(data) %*% data
      denominator <- 0.5*(ya+yb) + yab
    }
  }
  if(return == "vector" | is.null(return)){
    return(denominator[lower.tri(denominator)])
  }
  if(return == "matrix"){
    D = denominator
    D[upper.tri(D)] = t(D)[upper.tri(D)]
    colnames(D) = colnames(data)
    row.names(D) = colnames(data)
    diag(D) = 0
    return(D)
  }
}

rhat <- function(gbi) {                                                                 # lifted from aninet to make it work
  x <- get_numerator(gbi, data_format = "GBI", return = "vector")
  d <- get_denominator(gbi, data_format = "GBI", return = "vector")
  sri <- x/d
  bulk_rhat <- rhat_rfun(z_scale(gbi))
  gbi_folded <- abs(gbi - median(gbi))
  tail_rhat <- rhat_rfun(z_scale(gbi_folded))
  rhat <- max(bulk_rhat, tail_rhat)
  return(rhat)
}

#### run initial permutation ####
#set.seed(1)
#permute <- aninet::gbi_MCMC(data = gbi_males,
#                            ind_constraint = NULL, group_constraint = NULL,
#                            samples = 1000, # short chains to be repeatedly extended
#                            thin = 100, burnin = 1000,
#                            chains = num_chains,
#                            FUN = output_function) # CV, Mean, SD, SRI and Rhat

#outputs[which(outputs$permutation == 1),2:6] <- permute$FUN(gbi_males)
##permute_gbi[,,1,1] <- permute$permuted_data[[1]]
##permute_gbi[,,1,2] <- permute$permuted_data[[2]]

#### extend as far as necessary ####
#for(seed in seeds){
#  print(paste0('start seed number ', seed, ' at ', Sys.time()))
#  set.seed(seed)
#  permute <- aninet::extend.gbi_MCMC(permute, samples = 1000)
#  outputs[which(outputs$permutation == seed),2:6] <- permute$FUN(gbi_males)
#  #permute_gbi[,,seed,1] <- permute$permuted_data[[1]]
#  #permute_gbi[,,seed,2] <- permute$permuted_data[[2]]
#  # write out data every 100 extensions
#  if(seed %% 100 == 0){
#    write_csv(outputs, '../data_processed/motnp_permutation_cv_rhat.csv')
#    saveRDS(permute$permuted_data[[1]], '../data_processed/motnp_permutations_chain1.RDS')
#    saveRDS(permute$permuted_data[[2]], '../data_processed/motnp_permutations_chain2.RDS')
#    save.image('motnp_nonrandomtest_gbiMCMCpermutations.RData')
#  }
#}

# plots and p-values
#pdf('../outputs/')

#### set up asnipe run ####
# generate matrix
ids <- unique(c(counts_df$id_1, counts_df$id_2))            # IDs of elephants
num_eles <- length(ids)                                     # number of elephants
m_mat <- diag(nrow = num_eles)                              # matrix of males -- NxN to fill with SRI
colnames(m_mat) <- ids ; rownames(m_mat) <- ids             # dimnames = elephant IDs
print('matrix generated')

# obtain SRI
counts_df$sri <- counts_df$event_count / counts_df$count_dyad
print('SRI calculated')

# populate matrix with SRI values
for( i in 1:num_eles ) {                                    # rows
  for( j in 1:num_eles ) {                                  # columns
    if( i >= j ) { m_mat[i,j] <- m_mat[j,i] }               # make symmetrical about diagonal
    else {
      id1 <- colnames(m_mat)[i]                             # identify elephant 1 of dyad
      id2 <- rownames(m_mat)[j]                             # identify elephant 2 of dyad
      m_mat[i,j] <- counts_df$sri[which(counts_df$id_1 == id1 & counts_df$id_2 == id2)]      # add edge values -- match IDs and dyad
    }
  }
}
print('SRI matrix calculated')

### set up inputs
chain_length <- 1000
asnipe_seeds <- 2:(num_permute/chain_length)
obs <- nrow(gbi_males)

### set up output data frame
outputs <- data.frame(permutation = c(0,seq(from = chain_length, to = num_permute, by = chain_length)),
                      cv = NA, mean = NA, stdev = NA, rhat = NA, p = NA)

### fill in data for first row
outputs$cv[outputs$permutation == 0]    <- raster::cv(m_mat)
outputs$mean[outputs$permutation == 0]  <- mean(m_mat)
outputs$stdev[outputs$permutation == 0] <- sd(m_mat)
outputs$rhat[outputs$permutation == 0]  <- rhat(m_mat)
outputs$p[outputs$permutation == 0]     <- length(which(outputs$cv[outputs$permutation == 0] >= outputs$cv[outputs$permutation == 0]))/1 # for all others will be divided by outputs$permutation

### create vector of days for each sighting
sightings <- eles[,c('elephant','encounter','date')] %>%  # set up dataframe of actual encounters
  filter(elephant %in% ids) %>%                           # remove encounters that only include females, youngsters, or dead males
  dplyr::select(-elephant) %>%                                   # remove ID column
  distinct()                                              # cut down to one row per encounter

### run network permutations
set.seed(1)
random_networks <- asnipe::network_permutation(association_data = gbi_males, # permute network -- raw data = gbi_matrix
                                               association_matrix = m_mat,   # SRI matrix
                                               permutations = chain_length,  # 1000 times
                                               days = sightings$date,        # permute within day only
                                               within_day = TRUE)            # permute within day only

# calculate coefficient of variation for final randomised network
cv_random_networks <- rep(0,chain_length) ; for (i in c(1:chain_length)) {   # set up for loop
  net_rand <- random_networks[i,,]                                           # select random network
  cv_random_networks[i] <- raster::cv(net_rand)                              # calculate CV and save into vector
}

mat_new <- random_networks[chain_length,,]
outputs$cv[outputs$permutation == chain_length]    <- raster::cv(mat_new)
outputs$mean[outputs$permutation == chain_length]  <- mean(mat_new)
outputs$stdev[outputs$permutation == chain_length] <- sd(mat_new)
outputs$rhat[outputs$permutation == chain_length]  <- rhat(mat_new)
outputs$p[outputs$permutation == chain_length]     <- length(which(outputs$cv[outputs$permutation == chain_length] >= cv_random_networks))/outputs$permutation[outputs$permutation == chain_length]

# save every 1000th network
saved_networks <- array(NA, dim = c(nrow(outputs)-1, num_eles, num_eles),
                        dimnames = list(NULL, ids, ids))
saved_networks[1,,] <- mat_new

for(seed in asnipe_seeds){
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
  
  mat_new <- random_networks[chain_length,,]
  outputs$cv[outputs$permutation == permutation]    <- raster::cv(mat_new)
  outputs$mean[outputs$permutation == permutation]  <- mean(mat_new)
  outputs$stdev[outputs$permutation == permutation] <- sd(mat_new)
  outputs$rhat[outputs$permutation == permutation]  <- rhat(mat_new)
  outputs$p[outputs$permutation == permutation]     <- length(which(outputs$cv[outputs$permutation == permutation] >= cv_random_networks))/outputs$permutation[outputs$permutation == permutation]
  
  ### save every 1000th network
  saved_networks <- array(NA, dim = c(nrow(outputs)-1, num_eles, num_eles), dimnames = list(NULL, ids, ids))
  saved_networks[seed,,] <- mat_new
  
  ### save workspace every 10000th network
  if(seed %% 10 == 0) {
    saveRDS(cv_random_networks, '../data_processed/motnp_cv_permutations.csv')
    save.image('motnp_permutations.RData')
    }
}
