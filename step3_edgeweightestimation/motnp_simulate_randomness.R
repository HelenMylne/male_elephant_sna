#### set up ####
library(tidyverse)
library(asnipe)
library(bisonR)

load('motnp_bisonr_edgescalculated_strongprior.RData')
rm(list=setdiff(ls(), c("motnp_edge_weights_strongpriors","motnp_edges_null_strongpriors","counts_df")))

#### model comparison ####
model_comparison(models = list(random = motnp_edges_null_strongpriors, non_random = motnp_edge_weights_strongpriors) )
# Method: stacking
#           weight
# random     0.954 
# non_random 0.046 
# Warning message: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.

#### simulate population to test permutation sensitivity ####
# model-comparison shows non-random (but only 4.6% of dyads show "significantly" higher than average edge weight). Additional test: run simulation of network for 4.6% network structure -- weight[i,j] = Beta(alpha+friends_prob*gamma), alpha = mean of non-zero connections = Beta(alpha, error), friends_prob = Bernoulli(probability of being a strong connection == 4.6%), gamma = increase created by friendship --> loop through i and j (where j starts at i so only doing top half of matrix) { if w[i,j] ≤ alpha { if Bernoulli(0.8) = 0 then w[i,j] = 0 }} --> then transpose to bottom half of matrix --> run permutations and see if it says network is non-random

# simulate population
set.seed(12345)
N <- 186
ew <- data.frame(node_1 = rep(1:N, each = N),
                 node_2 = rep(1:N, N),
                 weight = NA) %>% 
  filter(node_1 != node_2)
ew$dyad <- ifelse(ew$node_1 > ew$node_2, paste0(ew$node_2,'_',ew$node_1), paste0(ew$node_1,'_',ew$node_2))
ew <- ew[,c(4,3)] %>% distinct()
ew <- ew %>% separate(dyad, into = c('node_1','node_2'), remove = F)
ew$node_1 <- as.numeric(ew$node_1)
ew$node_2 <- as.numeric(ew$node_2)

# alpha = mean of non-zero connections = Beta(alpha, error)
hist(LaplacesDemon::invlogit(motnp_edge_weights_strongpriors$edge_samples))
counts_df$sri <- counts_df$event_count / counts_df$count_dyad
hist(counts_df$sri)
a <- mean(counts_df$sri[counts_df$sri > 0] )
e <-  sd( counts_df$sri[counts_df$sri > 0] )
alpha <- rbeta(nrow(ew), a, e)
hist(alpha)

# friends_prob = Bernoulli(probability of being a strong connection == 4.6%)
ew$friends <- rbinom(n = nrow(ew), size = 1, prob = 0.046)
sum(ew$friends)/nrow(ew)
#ew$not_friends <- 1-ew$friends
ew$not_friends <- rbinom(n = nrow(ew), size = 1, prob = (1-0.046))

# gamma = edge weight increase created by friendship
gamma <- 0.1

# weight[i,j] = Beta(alpha+friends_prob*gamma)
hist(rbeta(10000, 0.07, 0.0), xlim = c(0,1))
hist(rbeta(10000, 0.07, 0.1), xlim = c(0,1))
hist(rbeta(10000, 0.07, 0.2), xlim = c(0,1))
hist(rbeta(10000, 0.07, 0.3), xlim = c(0,1))
hist(rbeta(10000, 0.1, 0.1), xlim = c(0,1))
hist(rbeta(10000, 0.2, 0.1), xlim = c(0,1))
hist(rbeta(10000, 0.3, 0.1), xlim = c(0,1))
hist(rbeta(10000, 0.5, 0.5), xlim = c(0,1))
hist(rbeta(10000, 1, 1), xlim = c(0,1))
hist(rbeta(10000, 3, 1), xlim = c(0,1))
hist(rbeta(10000, 3, 3), xlim = c(0,1))
hist(rbeta(10000, 1, 3), xlim = c(0,1))
hist(rbeta(10000, 0.1, 3), xlim = c(0,1))
hist(rbeta(10000, 1, 30), xlim = c(0,1))
#ew$weight <- rbeta(n = nrow(ew), shape1 = alpha, shape2 = gamma*ew$friends)
#ew$weight <- rbeta(n = nrow(ew), shape1 = a+gamma*ew$friends, shape2 = e)
#ew$weight <- rbeta(n = nrow(ew), shape1 = alpha+gamma*ew$friends, shape2 = e)
# Beta(times x did happen, times x didn't happen) --> Beta(a+Bern(0.046),a+Bern(0.956))
shapes <- simstudy::betaGetShapes(a+ew$friends*gamma, a+ew$not_friends*gamma)
ew$weight <- rbeta(n = nrow(ew), shape1 = shapes$shape1, shape2 = shapes$shape2)
hist(ew$weight, ylim = c(0,1000))

# loop through i and j (where j starts at i so only doing top half of matrix) { if w[i,j] ≤ alpha { if Bernoulli(0.8) = 0 then w[i,j] = 0 }} --> then transpose to bottom half of matrix


# run permutations and see if it says network is non-random



#### side note: what is CV output when SRI ≠ 1 or 0? #####
counts_df$sri
sri <- counts_df$sri[counts_df$sri > 0 & counts_df$sri < 1]
raster::cv(sri)

#### aninet non-random variation test ####
library(aninet)

### generate gbi_matrix
eles <- read_delim(file = '../data_processed/motnp_eles_long.csv', delim = ',')
eles$location <- paste(eles$gps_s, eles$gps_e, sep = '_') # make single variable for unique locations
eles <- eles[,c(1,16,2,3,17,4,5,14,7,8,10,13)]            # rearrange variables
nodes <- read_delim(file = '../data_processed/motnp_elenodes.csv', delim = ',') # read in node data
eles_asnipe <- eles[,c(3,4,2,5)]                                       # date, time, elephant, location
eles_asnipe$Date <- as.integer(eles_asnipe$date)                       # make date numeric
eles_asnipe$Date <- 1+eles_asnipe$Date - min(eles_asnipe$Date)         # start from 1, not 1st January 1970
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
motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')
ids <- unique(c(counts_df$id_1, counts_df$id_2))
gbi_males <- gbi_matrix[,colnames(gbi_matrix) %in% ids]
gbi_males <- gbi_males[rowSums(gbi_males) > 0,] # remove male-only sightings

( num_permutations <- aninet::gbi_MCMC_iters(gbi = gbi_males, target_samples = 1000, quiet = FALSE) )

permute_gbi <- gbi_MCMC(
  data = gbi_males,
  ind_constraint = NULL,
  group_constraint = NULL,
  samples = num_permutations,
  thin = 100,
  burnin = 1000,
  chains = 4,
  FUN = NULL      # leaving as NULL here will mean that permute_gbi$FUN(gbi_males) below gives CV, Mean, SD and SRI
)
str(permute_gbi)
permute_gbi$FUN(gbi_males)
permute_gbi$permuted_data


#### ignore all this -- I thought the permutations needed to be on SRI data but nope it just wants the GBI matrix ####
counts_df$sri <- counts_df$event_count / counts_df$count_dyad
sri_males <- matrix(data = NA, nrow = length(ids), ncol = length(ids), dimnames = list(x = ids, y = ids))
for(i in 1:length(ids)){
  for(j in 1:length(ids)){
    if( i == j ) { sri_males[i,j] <- 0 }
    else {
      if( i < j ){
        sri_males[i,j] <- counts_df$sri[which(counts_df$id_1 == rownames(sri_males)[i] &
                                                counts_df$id_2 == colnames(sri_males)[j])]
      }
        else {
          sri_males[i,j] <- counts_df$sri[which(counts_df$id_1 == rownames(sri_males)[j] &
                                                  counts_df$id_2 == colnames(sri_males)[i])]
      }
    }
  }
}

( num_permutations <- aninet::gbi_MCMC_iters(gbi = sri_males, target_samples = 1000, quiet = FALSE) )

permute_network <- gbi_MCMC(
  data = sri_males,
  ind_constraint = NULL,
  group_constraint = NULL,
  samples = num_permutations,
  thin = 100,
  burnin = 1000,
  chains = 4,
  FUN = NULL      # leaving as NULL here will mean that permute_network$FUN(gbi_males) below gives CV, Mean, SD and SRI
)
str(permute_network)
permute_network$FUN(gbi_males)
permute_network$permuted_data








