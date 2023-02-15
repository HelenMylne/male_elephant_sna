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
(a <- mean(counts_df$sri[counts_df$sri > 0] ))
(e <-  sd( counts_df$sri[counts_df$sri > 0] ))
alpha <- rbeta(nrow(ew), a, e)
hist(alpha)

# friends_prob = Bernoulli(probability of being a strong connection == 4.6%)
ew$friends <- rbinom(n = nrow(ew), size = 1, prob = 0.046)
sum(ew$friends)/nrow(ew)
#ew$not_friends <- 1-ew$friends
#ew$not_friends <- rbinom(n = nrow(ew), size = 1, prob = (1-0.046))

# gamma = edge weight increase created by friendship
gamma <- 0.1

### weight[i,j] = Beta(alpha+friends_prob*gamma) ####
## trying to find the shape parameters that give me what I want
#hist(rbeta(10000, 0.07, 0.0), xlim = c(0,1))
#hist(rbeta(10000, 0.07, 0.1), xlim = c(0,1))
#hist(rbeta(10000, 0.07, 0.2), xlim = c(0,1))
#hist(rbeta(10000, 0.07, 0.3), xlim = c(0,1))
hist(rbeta(10000, 0.07, 1), xlim = c(0,1))
#hist(rbeta(10000, 0.1, 0.1), xlim = c(0,1))
#hist(rbeta(10000, 0.2, 0.1), xlim = c(0,1))
#hist(rbeta(10000, 0.3, 0.1), xlim = c(0,1))
#hist(rbeta(10000, 0.5, 0.5), xlim = c(0,1))
#hist(rbeta(10000, 1, 1), xlim = c(0,1))
#hist(rbeta(10000, 3, 1), xlim = c(0,1))
#hist(rbeta(10000, 3, 3), xlim = c(0,1))
#hist(rbeta(10000, 1, 3), xlim = c(0,1))
#hist(rbeta(10000, 0.1, 3), xlim = c(0,1))
#hist(rbeta(10000, 1, 30), xlim = c(0,1))

## theory: shape1 = mean = alpha, shape2 = source of variation = gamma*ew$friends --> (0.07,0.01) or (0.07,0.00) --> (0.07,0.00)=1 so not this one
#ew$weight <- rbeta(n = nrow(ew), shape1 = alpha, shape2 = gamma*ew$friends)

## theory: error in alpha = stdev in weight --> shape1 = mean = a+f*g, shape2 = e
ew$weight <- rbeta(n = nrow(ew), shape1 = a+gamma*ew$friends, shape2 = e)
#shapes <- simstudy::betaGetShapes(a+ew$friends*gamma, e)
#ew$weight <- rbeta(n = nrow(ew), shape1 = shapes$shape1, shape2 = shapes$shape2)
#hist(ew$weight, ylim = c(0,1000)) # too stong on the 1, not enough low but >0

## theory: shape1 = mean, shape2 = precision = mean+sources of variation
ew$weight <- rbeta(n = nrow(ew), shape1 = a, shape2 = a+ew$friends*gamma)
hist(ew$weight, ylim = c(0,1000)) # too strong on extremes
shapes <- simstudy::betaGetShapes(a, a+ew$friends*gamma)
ew$weight <- rbeta(n = nrow(ew), shape1 = shapes$shape1, shape2 = shapes$shape2)
hist(ew$weight, ylim = c(0,1000)) # WAY too strong on extremes
ew$weight <- rbeta(n = nrow(ew), shape1 = alpha, shape2 = alpha+ew$friends*gamma)
hist(ew$weight) # too strong on extremes, not enough variation in middle

# theory: not_friends = Beta(a, not_a), friends = Beta(a+g, not a+g)
shapes_f <- simstudy::betaGetShapes(a, 1-a)
shapes_n <- simstudy::betaGetShapes(a+0.1, 1-(a+0.1))
ew$weight <- ifelse(ew$friends == 0,
                    rbeta(n = 1, shape1 = shapes_f$shape1, shape2 = shapes_f$shape2),
                    rbeta(n = 1, shape1 = shapes_n$shape1, shape2 = shapes_n$shape2))
hist(ew$weight) # entirely at extremes, mostly 1

## theory: Beta(times Z did happen, times Z didn't happen) --> shape1 = mean+g*friends = a+Bern(0.046), shape2 = mean+g*not_friends = a+Bern(0.956))
ew$weight <- rbeta(n = nrow(ew),
                   shape1 = a+rbinom(1,1,0.046),
                   shape2 = a+rbinom(1,1,0.956))
hist(ew$weight)
# don't alter shapes here because these are actual shape parameters rather than mean/precision values
#shapes <- simstudy::betaGetShapes(a+rbinom(1,1,0.046), a+rbinom(1,1,0.956))
#ew$weight <- rbeta(n = nrow(ew), shape1 = shapes$shape1, shape2 = shapes$shape2)
#hist(ew$weight)

# loop through i and j (where j starts at i so only doing top half of matrix) { if w[i,j] ≤ alpha { if Bernoulli(0.8) = 0 then w[i,j] = 0 }} --> then transpose to bottom half of matrix
ew_mat <- matrix(data = NA, nrow = N, ncol = N, dimnames = list(x = 1:N, y = 1:N))
for(i in 1:N){
  for(j in i:N){
    if( i == j ){ ew_mat[i,j] <- 1 }
    else{
      weight <- ew$weight[which(ew$node_1 == colnames(ew_mat)[i] & ew$node_2 == rownames(ew_mat)[j])]
      if( weight <= a ){
        w <- weight*rbinom(n = 1, size = 1, prob = 0.8)
      }
      else{ w <- weight }
      ew_mat[i,j] <- w
    }
  }
}
for(i in 1:N) {
  for(j in 1:i) {
    ew_mat[i,j] = ew_mat[j,i]
  }
}

## build gbi matrix -- determine what is a reasonable sampling period (at what point is it roughly independent? half a day or a day?)  --> build a matrix which is NxNxSP. (check asnipe default sampling period because haven't specified -- if sampling period = group observation then only ever permuting within group)
motnp <- read_csv('../data_processed/motnp_eles_long.csv') %>% 
  filter(elephant %in% unique(c(counts_df$id_1,counts_df$id_2)))
SP <- length(unique(motnp$encounter))
assn <- array(NA, c(N, N, SP), dimnames = list(1:N, 1:N, 1:SP))
for(k in 1:SP){
  for(i in 1:N){
    for(j in i:N){
      if(i == j){ assn[i,j,k] <- NA }
      else { assn[i,j,k] <- rbinom(1,1,ew_mat[i,j]) }
    }
  }
}
for(k in 1:SP){
  for(i in 1:N) {
    for(j in 1:i) {
      assn[i,j,k] <- assn[j,i,k]
    }
  }
}
sum(assn[,,1], na.rm = T)
assn[,,1]

gbi_mat <- matrix(0, nrow = SP, ncol = N)
for( i in 1:N ){
  for( j in 1:N ){
    for( k in 1:SP ){
      if(assn[i,j,k] == 1){
        gbi_mat[k,j] <- 1
        gbi_mat[k,i] <- 1
      }
    }
  }
}

# run permutations and see if it says network is non-random
N_networks <- 10000
random_networks <- asnipe::network_permutation(association_data = gbi_sim,
                                               association_matrix = ew_mat,
                                               permutations = N_networks)
cv_random_networks <- rep(0,N_networks)   # generate empty vector to fill with cv values
for (i in c(1:N_networks)) {
  net_rand <- random_networks[i,,]
  cv_random_networks[i] <- raster::cv(net_rand)
  if(i %% 1000 == 0) {print(i)}
}

#### side note: what is CV output when SRI ≠ 1 or 0? #####
counts_df$sri
sri <- counts_df$sri[counts_df$sri > 0 & counts_df$sri < 1]
raster::cv(sri)

