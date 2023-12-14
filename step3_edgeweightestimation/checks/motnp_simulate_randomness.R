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
ew$id_1 <- paste0('M',ew$node_1)
ew$id_2 <- paste0('M',ew$node_2)

# alpha = mean of non-zero connections = Beta(alpha, error)
hist(LaplacesDemon::invlogit(motnp_edge_weights_strongpriors$edge_samples))
counts_df$sri <- counts_df$event_count / counts_df$count_dyad
hist(counts_df$sri)
(a <- mean(counts_df$sri[counts_df$sri > 0] ))
(e <-  sd( counts_df$sri[counts_df$sri > 0] ))

# friends_prob = Bernoulli(probability of being a strong connection == 4.6%)
ew$friends <- rbinom(n = nrow(ew), size = 1, prob = 0.046)
sum(ew$friends)/nrow(ew)
#ew$not_friends <- 1-ew$friends
#ew$not_friends <- rbinom(n = nrow(ew), size = 1, prob = (1-0.046))

# gamma = edge weight increase created by friendship
gamma <- 0.15

#### weight[i,j] = Beta(alpha+friends_prob*gamma) ####
## theory: Beta(times Z did happen, times Z didn't happen) --> shape1 = mean+g*friends = a+Bern(0.046), shape2 = mean+g*not_friends = a+Bern(0.956))
ew$weight <- rbeta(n = nrow(ew),
                   shape1 = a+rbinom(1,1,0.046),
                   shape2 = a+rbinom(1,1,0.956))
hist(ew$weight) ## this looks right but if you look at the values there are a few which are way too high and not necessarily those which are friends

#### weight[i,j] = random, 2% add 0.15 -- HAVE SUCCESSFULLY PRODUCED A BETA DISTRIBUTION WITH MEAN AND STDEV SIMILAR TO THAT OF ORIGINAL DATA, BUT NOT EXCLUDING SRI=0 IN CALCULATIONS OF RANDOM VALUE. IF I MAKE IT TOO SIMILAR TO ORIGINAL DATA, DOES IT THEN ACTUALLY SHOW US ANYTHING?? ####
hist(rbeta(nrow(counts_df), a, e), breaks = 10) # too strong on the 1s
shapes <- simstudy::betaGetShapes(a,e)
hist(rbeta(nrow(counts_df), shapes$shape1, shapes$shape2), breaks = 10) # still too strong on the 1s and not enough in the middle
hist(counts_df$sri[counts_df$sri > 0], breaks = 10)

# code to calculate parameters for a beta distribution of mean = a and stdev = e taken from https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
estBetaParams <- function(mu, stdev) {
  alpha <- ((1 - mu) / (stdev^2) - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = c(alpha = alpha, beta = beta))
}

new_parameters <- estBetaParams(mu = a, stdev = e)
hist(rbeta(nrow(counts_df), new_parameters[1], new_parameters[2]), breaks = 10)
#new_shapes <- simstudy::betaGetShapes(new_parameters[1], new_parameters[2])
#hist(rbeta(nrow(counts_df), new_shapes$shape1, new_shapes$shape2), breaks = 10)

ew$weight <- rbeta(n = nrow(ew),
                   shape1 = new_parameters[1],
                   shape2 = new_parameters[2])
hist(counts_df$sri[counts_df$sri > 0], breaks = 20)
hist(ew$weight, xlim = c(0,1)) # tail too short and fat
mean(ew$weight)

a_zerosinc <- mean(counts_df$sri)
e_zerosinc <- sd(counts_df$sri)
new_parameters <- estBetaParams(mu = a_zerosinc, stdev = e_zerosinc)
hist(rbeta(nrow(counts_df), new_parameters[1], new_parameters[2]), breaks = 10)
#new_shapes <- simstudy::betaGetShapes(new_parameters[1], new_parameters[2])
#hist(rbeta(nrow(counts_df), new_shapes$shape1, new_shapes$shape2), breaks = 10)

ew$weight <- rbeta(n = nrow(ew),
                   shape1 = new_parameters[1],
                   shape2 = new_parameters[2])
hist(counts_df$sri, breaks = 20)
hist(ew$weight, xlim = c(0,1)) # tail too short and fat
mean(ew$weight)

# add 0.15 to random 4.6%
head(ew)
ew$friends <- rbinom(n = nrow(ew), size = 1, prob = 0.046)
ew$new_weight <- ew$weight+ew$friends*0.15
hist(ew$new_weight, xlim = c(0,1))
mean(ew$new_weight)

#### create edge weight matrix ####
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

#### build gbi matrix ####
# NOTE: need to determine what is a reasonable sampling period (at what point is it roughly independent? half a day or a day?)  --> build a matrix which is NxNxSP. (check asnipe default sampling period because haven't specified -- if sampling period = group observation then only ever permuting within group)
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

gbi_mat <- matrix(NA, nrow = SP, ncol = N, dimnames = list())
for( k in 1:SP ){
  for( i in 1:N ){
    for( j in i:N ){
      if(i != j){
        if(is.na(gbi_mat[k,j]) == TRUE ){
          if(assn[i,j,k] == 1){
            gbi_mat[k,i] <- 1
            gbi_mat[k,j] <- 1
          } 
          else {
            gbi_mat[k,i] <- 0
            gbi_mat[k,j] <- 0
          }
        }
        else{
          if(assn[i,j,k] == 1){
            gbi_mat[k,i] <- 1
          }
          else{
            gbi_mat[k,i] <- 0
          }
        }
      }
    }
  }
}
which(gbi_mat != 1)

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

#### new method to build GBI matrix ####
# build data frame where every individual is observed at least once
motnp <- read_csv('../data_processed/motnp_eles_long.csv') %>% 
  filter(elephant %in% unique(c(counts_df$id_1,counts_df$id_2)))
SP <- length(unique(motnp$encounter))
ids <- colnames(ew_mat) <- rownames(ew_mat) <- unique(c(ew$id_1, ew$id_2))
s <- data.frame(sighting = 1:SP,
                num_eles = NA,
                id_1 = c(ids, ids,
                         sample(ids, replace = F, size = SP-(2*length(ids)))),
                id_2 = NA,  id_3 = NA,  id_4 = NA,  id_5 = NA,
                id_6 = NA,  id_7 = NA,  id_8 = NA,  id_9 = NA,  id_10 = NA)

# simulate group size using rpois(), average = average of MOTNP
motnp$num_males <- NA
for(i in 1:nrow(motnp)){
  e <- motnp[motnp$encounter == motnp$encounter[i],]
  motnp$num_males[i] <- nrow(e)
  rm(e)
}
encounters <- motnp[,c('encounter','num_males')] %>% distinct()
average_num_eles <- round(mean(encounters$num_males))
set.seed(1) ; s$num_eles <- rpois(SP, average_num_eles)
s$num_eles <- ifelse(s$num_eles == 0, rpois(1, average_num_eles), s$num_eles)
table(s$num_eles)
round(mean(s$num_eles))

# simulate group IDs based on combination of SRIs for all those already identified
for(i in 1:SP){
  # id2
  possible_friends <- ids[which(ids != s$id_1[i])]
  probs <- ew_mat[rownames(ew_mat) == s$id_1[i],
                  colnames(ew_mat) %in% possible_friends]
  s$id_2[i] <- sample(possible_friends, prob = probs, size = 1)
  # id3 -- multiply probability of hanging out with ID1 by probability of hanging out with ID2 such that most likely to join are those who like both of them
  possible_friends <- possible_friends[which(possible_friends != s$id_2[i])]
  probs <- ew_mat[rownames(ew_mat) == s$id_1[i] | rownames(ew_mat) == s$id_2[i],
                  colnames(ew_mat) %in% possible_friends]
  probs_combined <- probs[1,]*probs[2,]
  s$id_3[i] <- sample(possible_friends, prob = probs_combined, size = 1)
  # id4
  possible_friends <- possible_friends[which(possible_friends != s$id_2[i])]
  probs <- ew_mat[rownames(ew_mat) == s$id_1[i] | rownames(ew_mat) == s$id_2[i],
                  colnames(ew_mat) %in% possible_friends]
  probs_combined <- probs[1,]*probs[2,]
  s$id_3[i] <- sample(possible_friends, prob = probs_combined, size = 1)
}

identify_group_members <- function(all_individuals, current_group, edgeweight_matrix){
  possible_friends <- all_individuals[!(all_individuals %in% current_group)]
  probs <- edgeweight_matrix[rownames(edgeweight_matrix) %in% current_group,
                             colnames(edgeweight_matrix) %in% possible_friends]
  if(length(current_group) > 1) {probs <- apply(probs, 2, prod)} # multiply probability of hanging out with ID1 by probability of hanging out with ID2 such that most likely to join are those who like both of them
  return(sample(possible_friends, prob = probs, size = 1))
}

for(i in 1:SP){
  s$id_2[i] <- identify_group_members(all_individuals = ids,
                                      current_group = s$id_1[i],
                                      edgeweight_matrix = ew_mat)
  s$id_3[i] <- identify_group_members(all_individuals = ids,
                                      current_group = c(s$id_1[i],s$id_2[i]),
                                      edgeweight_matrix = ew_mat)
  s$id_4[i] <- identify_group_members(all_individuals = ids,
                                      current_group = c(s$id_1[i],s$id_2[i],s$id_3[i]),
                                      edgeweight_matrix = ew_mat)
  s$id_5[i] <- identify_group_members(all_individuals = ids,
                                      current_group = c(s$id_1[i],s$id_2[i],s$id_3[i],
                                                        s$id_4[i]),
                                      edgeweight_matrix = ew_mat)
  s$id_6[i] <- identify_group_members(all_individuals = ids,
                                      current_group = c(s$id_1[i],s$id_2[i],s$id_3[i],
                                                        s$id_4[i],s$id_5[i]),
                                      edgeweight_matrix = ew_mat)
  s$id_7[i] <- identify_group_members(all_individuals = ids,
                                      current_group = c(s$id_1[i],s$id_2[i],s$id_3[i],
                                                        s$id_4[i],s$id_5[i],s$id_6[i]),
                                      edgeweight_matrix = ew_mat)
  s$id_8[i] <- identify_group_members(all_individuals = ids,
                                      current_group = c(s$id_1[i],s$id_2[i],s$id_3[i],
                                                        s$id_4[i],s$id_5[i],s$id_6[i],
                                                        s$id_7[i]),
                                      edgeweight_matrix = ew_mat)
  s$id_9[i] <- identify_group_members(all_individuals = ids,
                                      current_group = c(s$id_1[i],s$id_2[i],s$id_3[i],
                                                        s$id_4[i],s$id_5[i],s$id_6[i],
                                                        s$id_7[i],s$id_8[i]),
                                      edgeweight_matrix = ew_mat)
  s$id_10[i] <- identify_group_members(all_individuals = ids,
                                       current_group = c(s$id_1[i],s$id_2[i],s$id_3[i],
                                                         s$id_4[i],s$id_5[i],s$id_6[i],
                                                         s$id_7[i],s$id_8[i],s$id_9[i]),
                                       edgeweight_matrix = ew_mat)
}

# remove elephants which are not in original group based on group size
for(i in 1:SP){
  for(j in 2:10){
    s[i,j+2] <- ifelse(s$num_eles[i] < j, NA, s[i,j+2])
  }
}

# convert data frame of individual sightings to long format
sightings <- pivot_longer(s, cols = c(3:12), names_to = 'group_member', values_to = 'id')
sightings <- sightings[!is.na(sightings$id),]

# convert to gbi_matrix
gbi_mat <- matrix(NA, nrow = SP, ncol = N)
colnames(gbi_mat) <- ids
rownames(gbi_mat) <- 1:SP
for(i in 1:SP) {
  for(j in 1:N) {
    sighting <- sightings[sightings$sighting == rownames(gbi_mat)[i],]
    eles <- unique(sighting$id)
    gbi_mat[i,j] <- ifelse((colnames(gbi_mat)[j] %in% eles) == TRUE, 1, 0)
  }
}

#### run permutations ####
# calculate CV for SRI matrix
raster::cv(ew_mat)                 # 267.21

# set up permutations
N_networks <- 100000

# run network permutations
random_networks <- asnipe::network_permutation(association_data = gbi_mat,
                                               association_matrix = ew_mat,
                                               permutations = N_networks)
print('random networks generated')

cv_random_networks <- rep(0,N_networks)
for (i in c(1:N_networks)) {
  net_rand <- random_networks[i,,]
  cv_random_networks[i] <- raster::cv(net_rand)
  if(i %% 1000 == 0) {print(i)}
}
cv_random_networks
print(paste0('network permutations for entire network completed at ', Sys.time()))

# compare permuted networks to actual network
hist(cv_random_networks, xlim = c(min(cv_random_networks)-50,max(cv_random_networks)+50),
     main = 'permutations for network of all male elephants')
abline(v = raster::cv(ew_mat), col = 'red')
text(round(raster::cv(ew_mat),3), col = 'red',
     x = raster::cv(ew_mat)+30, y = 20000)

hist(cv_random_networks, xlim = c(min(cv_random_networks)-50,max(cv_random_networks)+50),
     main = 'permutations for network of all male elephants',
     ylim = c(0,1000))
abline(v = raster::cv(ew_mat), col = 'red')
text(round(raster::cv(ew_mat),3), col = 'red',
     x = raster::cv(ew_mat)+30, y = 600)

save.image('motnp_simulate_randomness.RData')

#### ignore: attempt to permute using igraph ####
# Create the weighted, unipartite, undirected network and calculate the observed Louvain modularity
#nodes <- read.delim("nodes.txt")
#links <- read.delim("links.txt")
net <- igraph::as.igraph(ew_mat)
net

modularity1 = cluster_louvain(anurosnet)
modularity1$modularity #observed value 

obs=modularity1$modularity
obs
real<-data.frame(obs)
real
# create empty vector
Nperm = 9 #I am starting with a low n, but intend to use at least 1000 permutations
randomized.modularity=matrix(nrow=length(obs),ncol=Nperm+1)
row.names(randomized.modularity)=names(obs)
randomized.modularity[,1]=obs 
randomized.modularity
#Permute the original network preserving its characteristics, calculate the Louvain modularity for all randomized networks, and compile the results in the vector
i<-1
while(i<=Nperm){ 
  
  randomnet <- rewire(anurosnet, with=each_edge(0.5)) #rewire vertices with constant probability
  E(randomnet)$weight <- sample(E(anurosnet)$weight) #shuffle initial weights and assign them randomly to edges
  
  mod<-(cluster_louvain(randomnet))
  
  mod$modularity
  
  linha = mod$modularity
  
  randomized.modularity[,i+1]=linha
  print(i)
  i=i+1
}
randomized.modularity #Here the result is not as expected
#Plot the observed value against the distribution of randomized values
niveis<-row.names(randomized.modularity)
for(k in niveis)
{
  if(any(is.na(randomized.modularity[k,]) == TRUE))
  {
    print(c(k, "metrica tem NA"))
  } else {
    nome.arq<- paste("modularity",k,".png", sep="")
    png(filename= nome.arq, res= 300, height= 15, width=21, unit="cm")
    plot(density(randomized.modularity[k,]), main="Observed vs. randomized",)
    abline(v=obs[k], col="red", lwd=2, xlab="")
    dev.off()
    print(k)
    nome.arq<- paste("Patefield_Null_mean_sd_",k,".txt", sep="")
    write.table(cbind(mean(randomized.modularity[k,]),sd(randomized.modularity[k,])), file=paste(nome.arq,sep=""), 
                sep=" ",row.names=TRUE,col.names=FALSE)
  }
}
#Estimate the P-value (significance)
significance=matrix(nrow=nrow(randomized.modularity),ncol=3)
row.names(significance)=row.names(randomized.modularity)
colnames(significance)=c("p (rand <= obs)", "p (rand >= obs)", "p (rand=obs)")

signif.sup=function(x) sum(x>=x[1])/length(x)
signif.inf=function(x) sum(x<=x[1])/length(x)
signif.two=function(x) ifelse(min(x)*2>1,1,min(x)*2)

significance[,1]=apply(randomized.modularity,1,signif.inf)
significance[,2]=apply(randomized.modularity,1,signif.sup)
significance[,3]=apply(significance[,-3],1,signif.two)

significance

#### side note: what is CV output when SRI ≠ 1 or 0? #####
counts_df$sri
sri <- counts_df$sri[counts_df$sri > 0 & counts_df$sri < 1]
raster::cv(sri)


#### test sampling period from creation of GBI matrix ####
eles <- read_delim(file = '../data_processed/motnp_eles_long.csv', delim = ',')
eles$location <- paste(eles$gps_s, eles$gps_e, sep = '_') # make single variable for unique locations
eles <- eles[,c(1,16,2,3,17,4,5,14,7,8,10,13)]            # rearrange variables
eles_asnipe <- eles[,c(3,4,2,5)] ; rm(eles) ; gc()                     # date, time, elephant, location
eles_asnipe$Date <- as.integer(eles_asnipe$date)                       # make date numeric
eles_asnipe$Date <- 1+eles_asnipe$Date - min(eles_asnipe$Date)         # start from 1, not 1st January 1970
eles_asnipe$Time <- ifelse(eles_asnipe$time > 1, NA, eles_asnipe$time) # time = proportion of day so anything >1 has to be wrong
eles_asnipe$Time <- eles_asnipe$Time*(24*60*60)                        # convert time values to seconds through day
eles_asnipe <- eles_asnipe[,c(5,6,3,4)]                                # create data frame to produce gbi matrix from
colnames(eles_asnipe) <- c('Date','Time','ID','Location')              # rename variables for get_gbi
eles_asnipe$ID <- as.character(eles_asnipe$ID)                         # correct data type
eles_asnipe$d_pad <- str_pad(eles_asnipe$Date, 3, pad = '0')           # 0-pad dates
eles_asnipe$encounter <- paste(eles_asnipe$d_pad, eles_asnipe$Time, eles_asnipe$Location, sep = '_') # unique value for each sighting
nodes <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds')
plot(colMeans(nodes), pch = 19, las = 1, col = rgb(0,0,1,0.2))
nodes <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>%
  select(-M51,-M113)                                 # remove dead elephants
nodes <- nodes[,which(colMeans(nodes) > 15)] %>%     # remove elephants below age of independence
  colnames()
eles_asnipe <- eles_asnipe[eles_asnipe$ID %in% nodes,]
eles_asnipe$group <- as.integer(as.factor(eles_asnipe$encounter)) # unique factor for every sighting
encounters <- eles_asnipe[,c('encounter','Time')] %>% distinct()
eles_asnipe <- eles_asnipe[,c(3,7)]            # create data table for gbi matrix
eles_asnipe <- data.table::setDT(eles_asnipe)  # create data table for gbi matrix
gbi_matrix <- spatsoc::get_gbi(DT = eles_asnipe, group = 'group', id = 'ID')  # create gbi matrix

encounters$sp <- NA
for(i in 1:nrow(encounters)-1){
  encounters$sp[i] <- encounters$Time[i+1]-encounters$Time[i]
}
encounters$sp[nrow(encounters)] <- mean(encounters$sp, na.rm = T)
asnipe::get_sampling_periods(association_data = gbi_matrix,
                     association_times = encounters$Time,
                     sampling_period = encounters$sp,
                     identities = colnames(gbi_matrix),
                     data_format = 'gbi')
