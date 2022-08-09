#### information #####
# To estimate the ages of individuals in other populations, first come up with a survival curve for the Amboseli elephants which have good birth and death data.
# survival curve options = exponential/logistic/Gompertz/Weibull, where shape = simple/bathtub/Makeham
# follow through Mia's code to determine which curve will best fit Amboseli data. she ran it on a high performance cluster, so it’s set up to run in parallel, but the basta model can also be run in parallel on its own (just like Stan). we will need to implement it as a custom function in Stan for the final analysis, given that we need to define it in the prior anyway. 

#### load packages ####
library(tidyverse)
library(cmdstanr)
library(ggdist)
library(posterior)
library(bayesplot)
library(rstan)
library(igraph)
library(LaplacesDemon)

#library(ggthemes)
#library(survival)
#library(tidybayes)
#library(data.table)

#### load model ####
# read in Stan model to estimate ages based on Gompertz bathtub distribution from ANP
latent_age_ordinal_model <- cmdstan_model("models/motnp_age_ordinal_regression.stan")

#### simulate process of assigning age categories to elephants ####
gompertz_bt <- function(a0, a1, c, b0, b1, age){
  gompertz <- exp(b0 + b1*age)
  bathtub <- exp(a0 - a1*age) + c + gompertz
  return(bathtub)
}

# draw new output curves
max_age <- 60
mortality <- gompertz_bt(a0 = -5.13,
                         a1 = 3.00,
                         c  = 0.026,
                         b0 = -5.08,
                         b1 = 0.09, 1:max_age)
plot(mortality)
probs <- 1 - (mortality/max(mortality))
plot(probs)

# create a fictional population with ages selected from that distribution
N <- 200
K <- 8
elephants_ls <- list(
  N = N,
  K = K,
  age = sample(1:max_age, N, replace = TRUE, prob = probs)
)
hist(elephants_ls$age)

# simulate observing ages with error and then binning them into age groups
error <- 3 # Error (in either direction) of age estimates
elephants_ls$age_guess <- elephants_ls$age + sample(-error:error, N, replace = TRUE)
elephants_ls$age_category_index <- sapply(elephants_ls$age_guess,
                                          function(x) which.max(x < c(5, 10, 15, 20, 25, 40, 60, Inf)))
hist(elephants_ls$age_category_index)

# look at actual age vs biologist assigned age: chance of being mis-classified is lower if actual age is in middle of age category, and biased towards lower end of class.
data.frame(elephants_ls) %>%
  ggplot(aes(x=age, y=age_guess, col=factor(age_category_index))) +
  geom_point(size=4,alpha=0.6) +
  geom_vline(xintercept=c(5, 10, 15, 20, 25, 40, 60), col=factor(1:7), linetype="dashed", alpha=0.6) +
  theme_minimal() + 
  xlab("Assigned age") + ylab("Age")

#### fit model to simulated dataset ####
#Fit model with cmdstanr
age_estimation_fit <- latent_age_ordinal_model$sample(
  data = elephants_ls, 
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000
)

# Examine the estimates. We can plot the estimated ages against the biologist assigned ages and ...
age_est_mat <- age_estimation_fit$summary()[202:401, ]
summary(age_est_mat)
hist(age_est_mat$mean)
hist(age_est_mat$rhat, breaks = 20)

plot_data <- data.frame(age = elephants_ls$age,       # Biologists original age est
                        model_age = age_est_mat$mean) # Mean modelled age

plot_data %>%
  ggplot(aes(x=factor(age), y=model_age)) +
  geom_point(size=4,col = 'blue', alpha=0.6) +
  geom_vline(xintercept=c(5, 10, 15, 20, 25, 40, 60), linetype="dashed", alpha=0.6) +
  geom_hline(yintercept=c(5, 10, 15, 20, 25, 40, 60), linetype="dashed", alpha=0.6) +
  geom_abline(slope = 1, intercept = 0)+
  scale_y_continuous(limits = c(0,60))+
  theme_minimal() + 
  xlab("Assigned age") + ylab("Modelled age")

#### load MOTNP data ####
motnp_males <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_elenodes.csv') %>% 
  #filter(dem_class == 'AM' | dem_class == 'PM') %>% 
  filter(sex == 'M')
unique(motnp_males$age_category)
motnp_males$age_cat_id <- ifelse(motnp_males$age_category == "0-3", 1, 
                                 ifelse(motnp_males$age_category == "3-4", 1,
                                        ifelse(motnp_males$age_category == "1-2", 1,
                                               ifelse(motnp_males$age_category == "7-8", 1,
                                                      ifelse(motnp_males$age_category == "4-5", 1, 
                                                             ifelse(motnp_males$age_category ==  "6-7", 1,
                                                                    ifelse(motnp_males$age_category == "8-9", 1, 
                                                                           ifelse(motnp_males$age_category == "5-6", 1, NA))))))))
motnp_males$age_cat_id <- ifelse(motnp_males$age_category == '9-10', 2,
                                 ifelse(motnp_males$age_category == '10-15', 3,
                                        ifelse(motnp_males$age_category == '15-19', 4,
                                               ifelse(motnp_males$age_category == '20-25', 5,
                                                      ifelse(motnp_males$age_category == '25-40', 6,
                                                             ifelse(motnp_males$age_category == '40+', 7,
                                                                    motnp_males$age_cat_id))))))

#### create data list ####
N_motnp <- nrow(motnp_males)
K <- 8
motnp_ls <- list(
  N = N_motnp,
  K = K,
  age_category_index = motnp_males$age_cat_id)
hist(motnp_ls$age_category_index)
#hist(elephants_ls$age_category_index)

#### fit model to MOTNP data ####
# Fit model with cmdstanr
age_motnp_fit <- latent_age_ordinal_model$sample(
  data = motnp_ls, 
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000
)

# Examine the estimates
age_est_mat <- age_motnp_fit$summary()[(N_motnp+2):(N_motnp*2+1), ]
summary(age_est_mat)
hist(age_est_mat$mean)
hist(age_est_mat$rhat, breaks = 20)

plot_data <- data.frame(age = ifelse(motnp_ls$age == 1, 3, 
                                     ifelse(motnp_ls$age == 2, 8,
                                            ifelse(motnp_ls$age == 3, 12,
                                                   ifelse(motnp_ls$age == 4, 18,
                                                          ifelse(motnp_ls$age == 5, 22, 
                                                                 ifelse(motnp_ls$age == 6, 32, 45)))))),
                        model_age = age_est_mat$mean) # Mean modelled age

plot_data %>%
  ggplot(aes(x=factor(age), y=model_age)) +
  geom_point(size=4,col = 'blue', alpha=0.6) +
  geom_vline(xintercept = c(5, 10, 15, 20, 25, 40, 60), linetype = "dashed", alpha = 0.6) +
  geom_hline(yintercept = c(5, 10, 15, 20, 25, 40, 60), linetype = "dashed", alpha = 0.6) +
  #geom_abline(slope = 1, intercept = 0)+
  scale_y_continuous(limits = c(0,60))+
  theme_minimal() + 
  xlab("Assigned age") + ylab("Modelled age")

# posterior predictive plot using draws from distribution to show uncertainty around mean age
true_ages <- age_motnp_fit$draws("true_age", format="df")
#mcmc_dens(true_ages)
true_ages <- true_ages[,1:N_motnp]

df <- as.data.frame(do.call(rbind, true_ages)) %>%
  mutate(age_cat = motnp_ls$age) %>% relocate(age_cat) %>%
  mutate(ID = motnp_males$id) %>% relocate(ID)

df <- df %>% pivot_longer(cols = 3:102) %>% select(-name)

df$true_age <- ifelse(df$age_cat == 1, 3, 
                      ifelse(df$age_cat == 2, 8,
                             ifelse(df$age_cat == 3, 12,
                                    ifelse(df$age_cat == 4, 18,
                                           ifelse(df$age_cat == 5, 22, 
                                                  ifelse(df$age_cat == 6, 32, 45))))))

df %>% ggplot(aes(x=true_age, y=value, group=factor(ID))) +
  geom_point(size=2,col = 'blue', alpha=0.1) +
  #stat_halfeye() +
  geom_vline(xintercept=c(5,10,15,20,25,40,60), linetype="dashed", alpha=0.6) +
  geom_hline(yintercept=c(5,10,15,20,25,40,60), linetype="dashed", alpha=0.6) +
  theme_bw() + 
  xlab("Assigned age") + ylab("Modelled age")

df$age_cat_centre <- ifelse(df$age_cat == 1, (0+5)/2,
                            ifelse(df$age_cat == 2, (5+10)/2,
                                   ifelse(df$age_cat == 3, (10+15)/2,
                                          ifelse(df$age_cat == 4, (15+20/2),
                                                 ifelse(df$age_cat == 5, (20+25)/2,
                                                        ifelse(df$age_cat == 6, (25+40)/2, 45))))))
df %>% ggplot() +
  geom_violin(aes(x = true_age, y = value, group = factor(age_cat)), fill = rgb(0,0,1,0.8))+
  #geom_point(aes(x = true_age, y = value, group = factor(ID)), size = 2,col = 'blue', alpha = 0.1) +
  #stat_halfeye() +
  geom_vline(xintercept = 0, alpha = 0.6) +
  geom_vline(xintercept = c(5,10,15,20,25,40,60), linetype = "dashed", alpha = 0.6) +
  geom_hline(yintercept = 0, alpha = 0.6) +
  geom_hline(yintercept = c(5,10,15,20,25,40,60), linetype = "dashed", alpha = 0.6) +
  theme_bw() + 
  xlab("Assigned age") + ylab("Modelled age")

### save output
colnames(true_ages) <- motnp_males$id
saveRDS(true_ages, file = '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_ageestimates_mcmcoutput.rds')

#### recreate network plots with new ages ####
draws_motnp2.2 <- readRDS('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_edgeweightestimates_mcmcoutput.rds') %>% 
  data.matrix()
summaries <- data.frame(dyad = colnames(draws_motnp2.2[,2:106954]),
                        min = rep(NA, ncol(draws_motnp2.2)-1),
                        max = rep(NA, ncol(draws_motnp2.2)-1),
                        mean = rep(NA, ncol(draws_motnp2.2)-1),
                        median = rep(NA, ncol(draws_motnp2.2)-1),
                        sd = rep(NA, ncol(draws_motnp2.2)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_motnp2.2[,i+1])
  summaries$max[i]    <- max(draws_motnp2.2[,i+1])
  summaries$mean[i]   <- mean(draws_motnp2.2[,i+1])
  summaries$median[i] <- median(draws_motnp2.2[,i+1])
  summaries$sd[i]     <- sd(draws_motnp2.2[,i+1])
}

par(mai = c(0.2,0.2,0.2,0.2))

# create array of draws per dyad (distributions)
counts_df <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_bayesian_binomialpairwiseevents.csv')
adj_tensor <- array(0, c(NROW(unique(counts_df$id_1))+1,
                         NROW(unique(counts_df$id_2))+1,
                         NROW(draws_motnp2.2)),
                    dimnames = list(c(unique(counts_df$id_1),'U9'),
                                    c('F1',unique(counts_df$id_2)),
                                    NULL))
N <- nrow(counts_df)
half_draws1 <- draws_motnp2.2[,2:54000]
half_draws2 <- draws_motnp2.2[,54001:N+1]
rm(draws_motnp2.2)

for (i in 1:54000) {            # can cope with jumps of 50000 dyads at a time
  dyad_row <- counts_df[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- half_draws1[, i]
}
rm(half_draws1)
for (i in 54001:N) {
  dyad_row <- counts_df[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- half_draws2[, i-54000]
}
rm(half_draws2)
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0
rm(adj_tensor, adj_quantiles, adj_lower, adj_upper, dyad_row)

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ages <- as.matrix(true_ages)
motnp_males$age_mean <- NA
for(i in 1:nrow(motnp_males)){
  motnp_males$age_mean[i] <- mean(as.numeric(ages[i]))
}

# create variables for different degrees of node connectedness
motnp_males$degree_0.1 <- NA
motnp_males$degree_0.2 <- NA
motnp_males$degree_0.3 <- NA
motnp_males$degree_0.4 <- NA
motnp_males$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(motnp_males)){
  rows <- summaries[summaries$id_1 == motnp_males$id[i] | summaries$id_2 == motnp_males$id[i],]
  motnp_males$degree_0.1[i] <- length(which(rows$median > 0.1))
  motnp_males$degree_0.2[i] <- length(which(rows$median > 0.2))
  motnp_males$degree_0.3[i] <- length(which(rows$median > 0.3))
  motnp_males$degree_0.4[i] <- length(which(rows$median > 0.4))
  motnp_males$degree_0.5[i] <- length(which(rows$median > 0.5))
}

which(motnp_males$degree_0.1 < motnp_males$degree_0.2)
which(motnp_males$degree_0.2 < motnp_males$degree_0.3)
which(motnp_males$degree_0.3 < motnp_males$degree_0.4)
which(motnp_males$degree_0.4 < motnp_males$degree_0.5)

### create data frame of all nodes to subset out males
all_eles <- counts_df[,c(2,4,12,14,16,18,22)] %>% distinct()
u9 <- counts_df[counts_df$id_2 == 'U9',c(3,5,13,15,17,19,23)] %>% distinct()
colnames(all_eles) <- c("id","node","name","age_class","age_category","sex","dem_class")
colnames(u9) <- c("id","node","name","age_class","age_category","sex","dem_class")
all_eles <- rbind(all_eles, u9) ; rm(u9)

males <- all_eles[all_eles$id %in% motnp_males$id,]
females <- anti_join(all_eles, males)

### create small df for plotting
plot_males <- motnp_males[,c('id','degree_0.3','age_mean','count')]
plot_males <- left_join(males, plot_males, by = 'id')

### All males
g_mid_m <- delete.vertices(graph = g_mid,
                           v = all_eles$id[all_eles$id %in% females$id])
g_rng_m <- delete.vertices(graph = g_rng,
                           v = all_eles$id[all_eles$id %in% females$id])
coords_m <- layout_nicely(g_mid_m)
plot(g_mid_m,
     edge.color = rgb(0,0,0,0.25),
     edge.width = E(g_rng_m)$weight,
     vertex.size = 7,
     vertex.label = NA,
     layout = coords_m)
plot(g_mid_m,
     edge.width = E(g_mid_m)$weight,
     edge.color = 'black',
     vertex.size = 7,
     vertex.label = plot_males$id,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = plot_males$age_mean,
     layout = coords_m, add = TRUE)

### Only males degree > 0.3
males0.3 <- plot_males[plot_males$degree_0.3 > 0,]
g_mid_m0.3 <- delete.vertices(graph = g_mid_m,
                              v = plot_males$id[plot_males$degree_0.3 == 0])
g_rng_m0.3 <- delete.vertices(graph = g_rng_m,
                              v = plot_males$id[plot_males$degree_0.3 == 0])

coords_m0.3 <- layout_nicely(g_mid_m0.3)
plot(g_mid_m0.3,
     edge.width = ifelse(E(g_mid_m0.3)$weight > 0.3,
                         E(g_rng_m0.3)$weight,0),
     edge.color = rgb(0,0,0,0.25),
     vertex.size = 7,
     vertex.label = NA,
     layout = coords_m0.3)
plot(g_mid_m0.3,
     edge.width = ifelse(E(g_mid_m0.3)$weight > 0.3,
                         E(g_mid_m0.3)$weight, 0),
     edge.color = 'black',
     vertex.size = 7,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = males0.3$age_mean,
     layout = coords_m0.3, add = TRUE)

plot(g_mid_m0.3,
     edge.width = ifelse(E(g_mid_m0.3)$weight > 0.3,
                         E(g_rng_m0.3)$weight,0),
     edge.color = rgb(0,0,0,0.25),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_m0.3)
plot(g_mid_m0.3,
     edge.width = ifelse(E(g_mid_m0.3)$weight > 0.3,
                         E(g_mid_m0.3)$weight, 0),
     edge.color = 'black',
     vertex.size = males0.3$count,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = males0.3$age_mean,
     layout = coords_m0.3, add = TRUE)

