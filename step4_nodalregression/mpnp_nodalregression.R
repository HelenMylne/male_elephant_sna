#### information ####
# Makgadikgadi Pans National Park -- regression to test effect of individual age on eigenvector and betweenness centrality

#### set up ####
#library(tidyverse)
#library(cmdstanr)
#library(posterior)
#library(MASS)
#library(rstan)
#library(igraph)
#library(LaplacesDemon)

library(tidyverse, lib.loc = 'packages/')
library(cmdstanr, lib.loc = 'packages/')
#library(ggdist, lib.loc = 'packages/')
library(posterior, lib.loc = 'packages/')
#library(bayesplot, lib.loc = 'packages/')
library(janitor, lib.loc = 'packages/')
library(readxl, lib.loc = 'packages/')
library(MASS, lib.loc = 'packages/')
library(rstan, lib.loc = 'packages/')
library(igraph, lib.loc = 'packages/')
library(LaplacesDemon, lib.loc = 'packages/')

# define PDF output
#pdf('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/outputs/mpnp_nodalregression_plots.pdf')

#### define the model ####
eigen_nodal_cont <- stan_model("models/nodal_regression_eigenvector_agedistribution.stan")

########### Time window 5 ###########
#### load mpnp5 nodes, edges and interactions data ####
### import data for aggregated model (binomial)
df_agg_mpnp5 <- read_delim('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_period5_pairwiseevents.csv', delim = ',')
df_agg_mpnp5 <- df_agg_mpnp5[,c(1:13,20:21)] %>% distinct()
df_agg_mpnp5 <- df_agg_mpnp5 %>%
  separate(id_1, into = c('BTF1','num1'), sep = 1, remove = F) %>% 
  separate(id_2, into = c('BTF2','num2'), sep = 1, remove = F) %>% 
  filter(BTF1 == 'B' & BTF2 == 'B')
df_agg_mpnp5 <- df_agg_mpnp5[,c(1:3,6,9:19)]

df_agg_mpnp5$age_mid_round_1 <- floor(df_agg_mpnp5$age_median_1)
df_agg_mpnp5$age_mid_round_2 <- floor(df_agg_mpnp5$age_median_2)
df_agg_mpnp5$node_1_nogaps <- as.integer(as.factor(df_agg_mpnp5$node_1))
df_agg_mpnp5$node_2_nogaps <- as.integer(as.factor(df_agg_mpnp5$node_2))+1
df_agg_mpnp5$dyad_id_nogaps <- as.integer(as.factor(df_agg_mpnp5$dyad))

### load the edge weights
mpnp5 <- readRDS('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp5_edgeweightestimates_mcmcoutput.rds')
mpnp5 <- mpnp5[,c(2:ncol(mpnp5))]
logit_edge_samples_mpnp5 <- logit(mpnp5)
print('logit complete')

### load age data
mpnp5_ages <- readRDS('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp5_ageestimates_mcmcoutput.rds')

#### calculate posterior centralities ####
mpnp5 <- as.matrix(logit_edge_samples_mpnp5)
edge_samples_mpnp5 <- plogis(mpnp5)
print('plogis complete')

# Build adjacency tensor
num_nodes <- length(unique(c(df_agg_mpnp5$id_1, df_agg_mpnp5$id_2)))
num_samples <- 4000
adj_tensor <- array(0, c(num_samples, num_nodes, num_nodes),
                    dimnames = list(NULL,
                                    unique(c(df_agg_mpnp5$id_1, df_agg_mpnp5$id_2)),
                                    unique(c(df_agg_mpnp5$id_1, df_agg_mpnp5$id_2))))
for (dyad_id in 1:nrow(df_agg_mpnp5)) {
  dyad_row <- df_agg_mpnp5[df_agg_mpnp5$dyad_id_nogaps == dyad_id, ]
  adj_tensor[, dyad_row$id_1, dyad_row$id_2] <- edge_samples_mpnp5[, dyad_id]
}
rm(dyad_row)

# Calculate centrality and store posterior samples in a matrix
eigen_samples_mpnp5 <- matrix(0, num_samples, num_nodes)
eigen_samples_mpnp5_std <- matrix(0, num_samples, num_nodes)
#btwn_samples_mpnp5 <- matrix(0, num_samples, num_nodes)
#btwn_samples_mpnp5_std <- matrix(0, num_samples, num_nodes)
for (i in 1:num_samples) {
  g <- graph_from_adjacency_matrix(adj_tensor[i, , ], mode="undirected", weighted=TRUE)
  eigen_samples_mpnp5[i, ] <- eigen_centrality(g)$vector
  eigen_samples_mpnp5_std[i, ] <- (eigen_samples_mpnp5[i, ] - mean(eigen_samples_mpnp5[i, ]))/sd(eigen_samples_mpnp5[i, ])
#  btwn_samples_mpnp5[i, ] <- betweenness(g, normalize = T)
#  btwn_samples_mpnp5_std[i, ] <- (btwn_samples_mpnp5[i, ] - mean(btwn_samples_mpnp5[i, ]))/sd(btwn_samples_mpnp5[i, ])
  rm(g)
}

# convert betweenness to a Z-score (REMEMBER BETWEENNESS IS INVERSED IN IGRAPH -- WILL NEED TO CONVERT IT SOMEHOW (1-BETWEENNESS??))
#mu_btwn <- mean(btwn_samples_mpnp5)
#sd_btwn <- sd(btwn_samples_mpnp5)
#btwn_samples_mpnp5_std <- (btwn_samples_mpnp5 - mu_btwn)/sd_btwn

head(eigen_samples_mpnp5)     # Unstandardised eigenvector centrality
head(eigen_samples_mpnp5_std) # Standardised eigenvector centrality
#head(btwn_samples_mpnp5)      # Unstandardised eigenvector centrality
#head(btwn_samples_mpnp5_std)  # Standardised eigenvector centrality
print('centralities calculated')

# visualise eigenvector centralities
df_wide_mpnp5 <- data.frame(eigen_samples_mpnp5_std)
colnames(df_wide_mpnp5) <- 1:num_nodes
df_long_mpnp5 <- pivot_longer(df_wide_mpnp5, cols = 1:num_nodes, names_to="node_id", values_to="eigenvector_centrality")
ggplot(df_long_mpnp5, aes(x = eigenvector_centrality)) +
  geom_density(fill="#387780", alpha=0.7, size=0.4) +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

nodes1 <- df_agg_mpnp5[,c(3,16)] %>% distinct()
nodes2 <- df_agg_mpnp5[,c(4,17)] %>% distinct()
colnames(nodes1) <- c('id', 'age_mid_round')
colnames(nodes2) <- c('id', 'age_mid_round')
node_ages_mpnp5 <- rbind(nodes1, nodes2) %>% distinct()
rm(nodes1, nodes2)

node_ages_mpnp5$node_id <- as.integer(as.factor(node_ages_mpnp5$id))
df_long_mpnp5$node_id <- as.integer(df_long_mpnp5$node_id)
df_long_mpnp5 <- left_join(df_long_mpnp5, node_ages_mpnp5, by = 'node_id')

ggplot(df_long_mpnp5, aes(x = eigenvector_centrality)) +
  geom_density(aes(fill = as.factor(age_mid_round)), alpha=0.5, size=0.4) +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_void() + 
  #theme_light() + 
  facet_grid(rows=vars(as.factor(node_id)), scales="free") +
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) # no effect

# visualise betweenness centrality
#df_wide_mpnp5 <- data.frame(btwn_samples_mpnp5_std)
#colnames(df_wide_mpnp5) <- 1:num_nodes
#df_long_mpnp5 <- pivot_longer(df_wide_mpnp5, cols = 1:num_nodes, names_to="node_id", values_to="btwn_centrality")
#ggplot(df_long_mpnp5, aes(x = btwn_centrality)) +
#  geom_density(fill="#387780", alpha=0.7, size=0.4) +
#  labs(x="betweenness centrality (standardised)") + 
#  theme_light() + 
#  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
#        axis.title.x = element_text(size=12),
#        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

#node_ages_mpnp5$node_id <- as.integer(as.factor(node_ages_mpnp5$id))
#df_long_mpnp5$node_id <- as.integer(df_long_mpnp5$node_id)
#df_long_mpnp5 <- left_join(df_long_mpnp5, node_ages_mpnp5, by = 'node_id')

#ggplot(df_long_mpnp5, aes(x = btwn_centrality)) +
#  geom_density(aes(fill = as.factor(age_mid_round)), alpha=0.7, size=0.4) +
#  labs(x="betweenness centrality (standardised)") + 
#  theme_light() + 
#  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
#        axis.title.x = element_text(size=12),
#        strip.text.y = element_text(size=12),
#        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) # well that would definitely suggest that there is no effect...

print('graphs plotted')

#### compute normal approximation ####
# The posterior centralities can now be characterised by the multivariate normal distribution as an approximation to their true posterior distributions. To highlight the importance of sampling from a multivariate normal rather than a univariate normal, we can plot centrality samples from two nodes against each other to see how they co-vary.
plot(eigen_samples_mpnp5_std[, 1], eigen_samples_mpnp5_std[, 2])
#plot(btwn_samples_mpnp5_std[, 1], btwn_samples_mpnp5_std[, 2])

#plot(eigen_samples_mpnp5 ~ btwn_samples_mpnp5)
#plot(eigen_samples_mpnp5_std ~ btwn_samples_mpnp5_std) # don't worry that high eigenvector doesn't necessarily equate to high betweenness -- https://datasciencegenie.com/what-is-centrality-in-graphs/ -- OR MAYBE DO WORRY - BETWEENNESS IS CALCULATED WEIRDLY IN IGRAPH, NEED TO ACCOUNT FOR THAT 

# This joint distribution contains important information that we would ideally like to include in the model. The multivariate normal will allow us to do this. Parameterising the multivariate normal approximation is a relatively simple process and can be done by calculating the sample mean and sample covariance matrix of the posterior centrality samples:
eigen_mu_mpnp5 <- apply(eigen_samples_mpnp5_std, 2, mean)
plot(density(eigen_mu_mpnp5)) # double peak -- not normally distributed
eigen_cov_mpnp5 <- cov(eigen_samples_mpnp5_std)
#btwn_mu_mpnp5 <- apply(btwn_samples_mpnp5_std, 2, mean) 
#btwn_cov_mpnp5 <- cov(btwn_samples_mpnp5_std)

# These quantities will be given to the Stan model as data to model joint posteriors of centrality in the regression. We can run a few quick plots to see how well the approximation is working and the covariance for one of the nodes:
eigen_samples_mpnp5_sim <- MASS::mvrnorm(1e5, eigen_mu_mpnp5, eigen_cov_mpnp5)
plot(density(eigen_samples_mpnp5_std[, 1]), lwd=2, main="Estimated standardised eigenvector centrality vs normal approximation", xlab="Logit edge weight")
lines(density(eigen_samples_mpnp5_sim[, 1]), col=rgb(0, 0, 1, 0.5), lwd=2)

#### eigenvector prior predictive checks (univariate Gaussian) ####
range(eigen_mu_mpnp5)
age <- 1:60

# eigenvector linear effect only
plot(NULL, xlim = c(0,60), ylim = c(-4,4), las = 1, xlab = 'age', ylab = 'mean eigenvector centrality (std)')
abline(h = range(eigen_mu_mpnp5), lty = 2) # eigen_samples_mpnp5_std or eigen_mu_mpnp5??
for(i in 1:100){
  intercept <- rnorm(1, 0, 1)
  beta_age <- rnorm(1, 0, 0.02)
  lines(x = age, y = intercept + beta_age*age,
        col = rgb(0,0,1,0.4))
}

plot(NULL, xlim = c(0,60), ylim = c(-4,4), las = 1, xlab = 'age', ylab = 'mean eigenvector centrality (std)')
abline(h = range(eigen_samples_mpnp5_std), lty = 2) # eigen_samples_mpnp5_std or eigen_mu_mpnp5??
for(i in 1:100){
  intercept <- rnorm(1, 0, 1)
  beta_age <- rnorm(1, 0, 0.02)
  lines(x = age, y = intercept + beta_age*age,
        col = rgb(0,0,1,0.4))
}

# eigenvector allowing for extremes to be most/least central
plot(NULL, xlim = c(0,60), ylim = c(-4,4), las = 1, xlab = 'age', ylab = 'mean eigenvector centrality (std)')
abline(h = range(eigen_mu_mpnp5), lty = 2)
for(i in 1:100){
  intercept <- rnorm(1, 0, 1)
  beta_age <- rnorm(1, 0, 0.02)
  beta_age2 <- rnorm(1, 0, 0.003)
  lines(x = age, y = intercept + beta_age*age + beta_age2*(age^2),
        col = rgb(0,0,1,0.4))
}

plot(NULL, xlim = c(0,60), ylim = c(-4,4), las = 1, xlab = 'age', ylab = 'mean eigenvector centrality (std)')
abline(h = range(eigen_samples_mpnp5_std), lty = 2)
for(i in 1:100){
  intercept <- rnorm(1, 0, 1)
  beta_age <- rnorm(1, 0, 0.02)
  beta_age2 <- rnorm(1, 0, 0.003)
  lines(x = age, y = intercept + beta_age*age + beta_age2*(age^2),
        col = rgb(0,0,1,0.4))
}

#### eigenvector ~ age ####
# create data list
modeldata_ages <- as.matrix(mpnp5_ages) %>% t() # transposed to fit matrix input specified in model
modeldata_ages_sq <- (modeldata_ages)^2
eigen_data_mpnp5 <- list(
  num_nodes = num_nodes,               # Number of dyads
  centrality_mu  = eigen_mu_mpnp5,     # Sample means of logit edge weights
  centrality_cov = eigen_cov_mpnp5,    # Sample covariance of logit edge weights
  node_age = modeldata_ages,           # Age of individual -- are these in the wrong order? should be M1, M10, M100, M101?
  node_age2 = modeldata_ages_sq        # Age squared for quadratic term
)
str(eigen_data_mpnp5)

# plot data
plot(eigen_data_mpnp5$centrality_mu ~ eigen_data_mpnp5$node_age[,1], las = 1, pch = 19, col = rgb(0,0,1,0.3))

# fit model to data
fit_eigen_cont <- sampling(eigen_nodal_cont, data = eigen_data_mpnp5, 
                           #iter = 10000,         # added in because high Rhat recommended running for more iterations but did not improve anything and just made the time take longer
                           cores = 4, chains = 4)

#### diagnostics ####
traceplot(fit_eigen_cont) # these don't look that great -- several sections where it gets stuck in a very small area

params <- rstan::extract(fit_eigen_cont)
plot(density(eigen_samples_mpnp5_std[1, ]), main = "Posterior predictive density (standardised centrality)",
     col = rgb(0, 0, 0, 0.25), ylim = c(0, 1), las = 1)
for (i in 1:100) {
  j <- sample(1:num_samples, 1)                                         # select a network to plot
  lines(density(eigen_samples_mpnp5_std[j, ]), col=rgb(0, 0, 0, 0.25)) # plot centrality density for network j (black)
  mu <- params$beta_age[j]*apply(eigen_data_mpnp5$node_age,1,mean)     # extract age slope parameter for network j
  sigma <- eigen_cov_mpnp5 + diag(rep(params$sigma[j], num_nodes))     # extract variance parameter for network j
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))  # plot predicted centralities for network j (blue)
} # All quite messy but the right shape!

### plot age against predictions
node_age <- 1:60
beta_age_mean  <- summary(fit_eigen_cont)$summary[1,1]
beta_age_2.5   <- summary(fit_eigen_cont)$summary[1,4]
beta_age_97.5  <- summary(fit_eigen_cont)$summary[1,8]
beta_age2_mean <- summary(fit_eigen_cont)$summary[2,1]
beta_age2_2.5  <- summary(fit_eigen_cont)$summary[2,4]
beta_age2_97.5 <- summary(fit_eigen_cont)$summary[2,8]
intercept_mean <- summary(fit_eigen_cont)$summary[3,1]
intercept_2.5  <- summary(fit_eigen_cont)$summary[3,4]
intercept_97.5 <- summary(fit_eigen_cont)$summary[3,8]
mu <- intercept_mean  + beta_age_mean*node_age + beta_age2_mean*(node_age^2)
min <- intercept_2.5  + beta_age_2.5*node_age  + beta_age2_2.5*(node_age^2)
max <- intercept_97.5 + beta_age_97.5*node_age + beta_age2_97.5*(node_age^2)
summary(mu)
summary(min)
summary(max)
plot(NULL, xlim = c(0,60), ylim = c(-5, 5), lwd = 2, las = 1,
     xlab = 'age (years)', ylab = 'eigenvector centrality (std)')
for(i in 1:100){
  j <- sample(1:4000, 1)
  int <- params$intercept[j]
  bA1 <- params$beta_age[j]
  bA2 <- params$beta_age2[j]
  mean <- int + bA1*age + bA2*age
  lines(c(1:60), mean, col = rgb(0.5,0,1,0.2))
}
lines(mu ~ node_age, lwd = 2) 
lines(min ~ node_age, col = 'red')
lines(max ~ node_age, col = 'red')
abline(h = c(-2,2), lty = 2)
for(i in 1:num_nodes){
  s <- sample(x = 1:num_samples, size = 100, replace = F)
  for(j in 1:length(s)){
    points(eigen_mu_mpnp5[i] ~ eigen_data_mpnp5$node_age[i,s[j]], pch = 19, col = rgb(0,0,0,0.05))
  }
}

#### interpreting the model ####
# Calculate the 95% credible intervals of the model parameters
round(summary(fit_eigen_cont)$summary[1:3, c(1:4, 8)], 2) # basically no effect whatsoever

# clean environment
rm(adj_tensor,age,beta_age,beta_age2,btwn_cov_mpnp5,btwn_mu_mpnp5,btwn_samples_mpnp5,btwn_samples_mpnp5_std,df_agg_mpnp5,df_long_mpnp5,df_wide_mpnp5,dyad_id,edge_samples_mpnp5,eigen_cov_mpnp5,eigen_data_mpnp5,eigen_mu_mpnp5,eigen_samples_mpnp5,eigen_samples_mpnp5_sim,eigen_samples_mpnp5_std,i,intercept,logit_edge_samples_mpnp5,modeldata_ages,modeldata_ages_sq,mpnp5,mpnp5_ages,mu_btwn,node_ages_mpnp5,num_nodes,num_samples,sd_btwn,bA1,bA2,beta_age_2.5,beta_age_97.5,beta_age_mean,beta_age2_2.5,beta_age2_97.5,beta_age2_mean,fit_eigen_cont,int,intercept_2.5,intercept_97.5,intercept_mean,j,max,mean,min,mu,node_age,params,s,sigma)

########### Time window 4 ###########
#### load mpnp4 nodes, edges and interactions data ####
### import data for aggregated model (binomial)
df_agg_mpnp4 <- read_delim('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_period4_pairwiseevents.csv', delim = ',')
df_agg_mpnp4 <- df_agg_mpnp4[,c(1:13,20:21)] %>% distinct()
df_agg_mpnp4 <- df_agg_mpnp4 %>%
  separate(id_1, into = c('BTF1','num1'), sep = 1, remove = F) %>% 
  separate(id_2, into = c('BTF2','num2'), sep = 1, remove = F) %>% 
  filter(BTF1 == 'B' & BTF2 == 'B')
df_agg_mpnp4 <- df_agg_mpnp4[,c(1:3,6,9:19)]

df_agg_mpnp4$age_mid_round_1 <- floor(df_agg_mpnp4$age_median_1)
df_agg_mpnp4$age_mid_round_2 <- floor(df_agg_mpnp4$age_median_2)
df_agg_mpnp4$node_1_nogaps <- as.integer(as.factor(df_agg_mpnp4$node_1))
df_agg_mpnp4$node_2_nogaps <- as.integer(as.factor(df_agg_mpnp4$node_2))+1
df_agg_mpnp4$dyad_id_nogaps <- as.integer(as.factor(df_agg_mpnp4$dyad))

### load the edge weights
mpnp4 <- readRDS('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp4_edgeweightestimates_mcmcoutput.rds')
mpnp4 <- mpnp4[,c(2:ncol(mpnp4))]
mpnp4 <- mpnp4[, which(colnames(mpnp4) %in% df_agg_mpnp4$dyad)]
logit_edge_samples_mpnp4 <- logit(mpnp4)
print('logit complete')

### load age data
mpnp4_ages <- readRDS('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp4_ageestimates_mcmcoutput.rds')

#### calculate posterior centralities ####
mpnp4 <- as.matrix(logit_edge_samples_mpnp4)
edge_samples_mpnp4 <- plogis(mpnp4)
print('plogis complete')

# Build adjacency tensor
num_nodes <- length(unique(c(df_agg_mpnp4$id_1, df_agg_mpnp4$id_2)))
num_samples <- 4000
adj_tensor <- array(0, c(num_samples, num_nodes, num_nodes),
                    dimnames = list(NULL,
                                    unique(c(df_agg_mpnp4$id_1, df_agg_mpnp4$id_2)),
                                    unique(c(df_agg_mpnp4$id_1, df_agg_mpnp4$id_2))))
for (dyad_id in 1:nrow(df_agg_mpnp4)) {
  dyad_row <- df_agg_mpnp4[df_agg_mpnp4$dyad_id_nogaps == dyad_id, ]
  adj_tensor[, dyad_row$id_1, dyad_row$id_2] <- edge_samples_mpnp4[, dyad_id]
}
rm(dyad_row)

# Calculate centrality and store posterior samples in a matrix
eigen_samples_mpnp4 <- matrix(0, num_samples, num_nodes)
eigen_samples_mpnp4_std <- matrix(0, num_samples, num_nodes)
#btwn_samples_mpnp4 <- matrix(0, num_samples, num_nodes)
#btwn_samples_mpnp4_std <- matrix(0, num_samples, num_nodes)
for (i in 1:num_samples) {
  g <- graph_from_adjacency_matrix(adj_tensor[i, , ], mode="undirected", weighted=TRUE)
  eigen_samples_mpnp4[i, ] <- eigen_centrality(g)$vector
  eigen_samples_mpnp4_std[i, ] <- (eigen_samples_mpnp4[i, ] - mean(eigen_samples_mpnp4[i, ]))/sd(eigen_samples_mpnp4[i, ])
#  btwn_samples_mpnp4[i, ] <- betweenness(g, normalize = T)
#  btwn_samples_mpnp4_std[i, ] <- (btwn_samples_mpnp4[i, ] - mean(btwn_samples_mpnp4[i, ]))/sd(btwn_samples_mpnp4[i, ])
  rm(g)
}

# convert betweenness to a Z-score
#mu_btwn <- mean(btwn_samples_mpnp4)
#sd_btwn <- sd(btwn_samples_mpnp4)
#btwn_samples_mpnp4_std <- (btwn_samples_mpnp4 - mu_btwn)/sd_btwn

head(eigen_samples_mpnp4)     # Unstandardised eigenvector centrality
head(eigen_samples_mpnp4_std) # Standardised eigenvector centrality
#head(btwn_samples_mpnp4)      # Unstandardised eigenvector centrality
#head(btwn_samples_mpnp4_std)  # Standardised eigenvector centrality
print('centralities calculated')

# visualise eigenvector centralities
df_wide_mpnp4 <- data.frame(eigen_samples_mpnp4_std)
colnames(df_wide_mpnp4) <- 1:num_nodes
df_long_mpnp4 <- pivot_longer(df_wide_mpnp4, cols = 1:num_nodes, names_to="node_id", values_to="eigenvector_centrality")
ggplot(df_long_mpnp4, aes(x = eigenvector_centrality)) +
  geom_density(fill="#387780", alpha=0.7, size=0.4) +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) # NOT NICELY NORMAL!

nodes1 <- df_agg_mpnp4[,c(3,16)] %>% distinct()
nodes2 <- df_agg_mpnp4[,c(4,17)] %>% distinct()
colnames(nodes1) <- c('id', 'age_mid_round')
colnames(nodes2) <- c('id', 'age_mid_round')
node_ages_mpnp4 <- rbind(nodes1, nodes2) %>% distinct()
rm(nodes1, nodes2)

node_ages_mpnp4$node_id <- as.integer(as.factor(node_ages_mpnp4$id))
df_long_mpnp4$node_id <- as.integer(df_long_mpnp4$node_id)
df_long_mpnp4 <- left_join(df_long_mpnp4, node_ages_mpnp4, by = 'node_id')

ggplot(df_long_mpnp4, aes(x = eigenvector_centrality)) +
  geom_density(aes(fill = as.factor(age_mid_round)), alpha=0.5, size=0.4) +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) # NOT NICELY NORMAL!

# visualise betweenness centrality
#df_wide_mpnp4 <- data.frame(btwn_samples_mpnp4_std)
#colnames(df_wide_mpnp4) <- 1:num_nodes
#df_long_mpnp4 <- pivot_longer(df_wide_mpnp4, cols = 1:num_nodes, names_to="node_id", values_to="btwn_centrality")
#ggplot(df_long_mpnp4, aes(x = btwn_centrality)) +
#  geom_density(fill="#387780", alpha=0.7, size=0.4) +
#  labs(x="betweenness centrality (standardised)") + 
#  theme_light() + 
#  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
#        axis.title.x = element_text(size=12),
#        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#
#node_ages_mpnp4$node_id <- as.integer(as.factor(node_ages_mpnp4$id))
#df_long_mpnp4$node_id <- as.integer(df_long_mpnp4$node_id)
#df_long_mpnp4 <- left_join(df_long_mpnp4, node_ages_mpnp4, by = 'node_id')
#
#ggplot(df_long_mpnp4, aes(x = btwn_centrality)) +
#  geom_density(aes(fill = as.factor(age_mid_round)), alpha=0.7, size=0.4) +
#  labs(x="betweenness centrality (standardised)") + 
#  theme_light() + 
#  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
#        axis.title.x = element_text(size=12),
#        strip.text.y = element_text(size=12),
#        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) # ZERO BETWEENNESS
#
print('graphs plotted')

#### compute normal approximation ####
# The posterior centralities can now be characterised by the multivariate normal distribution as an approximation to their true posterior distributions. To highlight the importance of sampling from a multivariate normal rather than a univariate normal, we can plot centrality samples from two nodes against each other to see how they co-vary.
plot(eigen_samples_mpnp4_std[, 1], eigen_samples_mpnp4_std[, 2])
#plot(btwn_samples_mpnp4_std[, 1], btwn_samples_mpnp4_std[, 2])

#plot(eigen_samples_mpnp4 ~ btwn_samples_mpnp4)
#plot(eigen_samples_mpnp4_std ~ btwn_samples_mpnp4_std)

# This joint distribution contains important information that we would ideally like to include in the model. The multivariate normal will allow us to do this. Parameterising the multivariate normal approximation is a relatively simple process and can be done by calculating the sample mean and sample covariance matrix of the posterior centrality samples:
eigen_mu_mpnp4 <- apply(eigen_samples_mpnp4_std, 2, mean)
plot(density(eigen_mu_mpnp4)) # double peak -- not normally distributed
eigen_cov_mpnp4 <- cov(eigen_samples_mpnp4_std)
#btwn_mu_mpnp4 <- apply(btwn_samples_mpnp4_std, 2, mean) 
#btwn_cov_mpnp4 <- cov(btwn_samples_mpnp4_std)

# These quantities will be given to the Stan model as data to model joint posteriors of centrality in the regression. We can run a few quick plots to see how well the approximation is working and the covariance for one of the nodes:
eigen_samples_mpnp4_sim <- MASS::mvrnorm(1e5, eigen_mu_mpnp4, eigen_cov_mpnp4)
plot(density(eigen_samples_mpnp4_std[, 1]), lwd=2, main="Estimated standardised centrality vs normal approximation", xlab="Logit edge weight")
lines(density(eigen_samples_mpnp4_sim[, 1]), col=rgb(0, 0, 1, 0.5), lwd=2)

#### eigenvector prior predictive checks (univariate Gaussian) ####
range(eigen_mu_mpnp4)
age <- 1:60

# eigenvector linear effect only
plot(NULL, xlim = c(0,60), ylim = c(-4,4), las = 1, xlab = 'age', ylab = 'mean eigenvector centrality (std)')
abline(h = range(eigen_samples_mpnp4_std), lty = 2)
for(i in 1:100){
  intercept <- rnorm(1, -1, 1)
  beta_age <- rnorm(1, 0, 0.02)
  lines(x = age, y = intercept + beta_age*age,
        col = rgb(0,0,1,0.4))
}

# eigenvector allowing for extremes to be most/least central
plot(NULL, xlim = c(0,60), ylim = c(-4,4), las = 1, xlab = 'age', ylab = 'mean eigenvector centrality (std)')
abline(h = range(eigen_samples_mpnp4_std), lty = 2)
for(i in 1:100){
  intercept <- rnorm(1, -1, 1)
  beta_age <- rnorm(1, 0, 0.02)
  beta_age2 <- rnorm(1, 0, 0.003)
  lines(x = age, y = intercept + beta_age*age + beta_age2*(age^2),
        col = rgb(0,0,1,0.4))
}

#### fit the model ####
modeldata_ages <- as.matrix(mpnp4_ages) %>% t() # CHECK WHY DID I TRANSPOSE THIS???
modeldata_ages_sq <- (modeldata_ages)^2

#### eigenvector ~ age ####
eigen_data_mpnp4 <- list(
  num_nodes = num_nodes,               # Number of dyads
  centrality_mu  = eigen_mu_mpnp4,     # Sample means of logit edge weights
  centrality_cov = eigen_cov_mpnp4,    # Sample covariance of logit edge weights
  node_age = modeldata_ages,           # Age of individual -- are these in the wrong order? should be M1, M10, M100, M101?
  node_age2 = modeldata_ages_sq        # Age squared for quadratic term
)
str(eigen_data_mpnp4)

plot(eigen_data_mpnp4$centrality_mu ~ eigen_data_mpnp4$node_age[,1], las = 1, pch = 19, col = rgb(0,0,1,0.3))

fit_eigen_cont <- sampling(eigen_nodal_cont, data = eigen_data_mpnp4, cores = 4, chains = 4)

#### diagnostics ####
traceplot(fit_eigen_cont) # these look good

params <- rstan::extract(fit_eigen_cont)
plot(density(eigen_samples_mpnp4_std[1, ]), main = "Posterior predictive density (standardised centrality)",
     col = rgb(0, 0, 0, 0.25), ylim = c(0, 0.8), las = 1)
for (i in 1:100) {
  j <- sample(1:num_samples, 1)                                         # select a network to plot
  lines(density(eigen_samples_mpnp4_std[j, ]), col=rgb(0, 0, 0, 0.25)) # plot centrality density for network j (black)
  mu <- params$beta_age[j]*apply(eigen_data_mpnp4$node_age,1,mean)     # extract age slope parameter for network j
  sigma <- eigen_cov_mpnp4 + diag(rep(params$sigma[j], num_nodes))     # extract variance parameter for network j
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))  # plot predicted centralities for network j (blue)
} # not too bad -- this makes sense considering the main data were left skewed but the model is normal

### plot age against predictions
node_age <- 1:60
beta_age_mean  <- summary(fit_eigen_cont)$summary[1,1]
beta_age_2.5   <- summary(fit_eigen_cont)$summary[1,4]
beta_age_97.5  <- summary(fit_eigen_cont)$summary[1,8]
beta_age2_mean <- summary(fit_eigen_cont)$summary[2,1]
beta_age2_2.5  <- summary(fit_eigen_cont)$summary[2,4]
beta_age2_97.5 <- summary(fit_eigen_cont)$summary[2,8]
intercept_mean <- summary(fit_eigen_cont)$summary[3,1]
intercept_2.5  <- summary(fit_eigen_cont)$summary[3,4]
intercept_97.5 <- summary(fit_eigen_cont)$summary[3,8]
mu <- intercept_mean  + beta_age_mean*node_age + beta_age2_mean*(node_age^2)
min <- intercept_2.5  + beta_age_2.5*node_age  + beta_age2_2.5*(node_age^2)
max <- intercept_97.5 + beta_age_97.5*node_age + beta_age2_97.5*(node_age^2)
summary(mu)
summary(min)
summary(max)
plot(NULL, xlim = c(0,60), ylim = c(-5, 5), lwd = 2, las = 1,
     xlab = 'age (years)', ylab = 'eigenvector centrality (std)')
for(i in 1:100){
  j <- sample(1:4000, 1)
  int <- params$intercept[j]
  bA1 <- params$beta_age[j]
  bA2 <- params$beta_age2[j]
  mean <- int + bA1*age + bA2*age
  lines(c(1:60), mean, col = rgb(0.5,0,1,0.2))
}
lines(mu ~ node_age, lwd = 2) 
lines(min ~ node_age, col = 'red')
lines(max ~ node_age, col = 'red')
abline(h = c(-2,2), lty = 2)
for(i in 1:num_nodes){
  s <- sample(x = 1:num_samples, size = 100, replace = F)
  for(j in 1:length(s)){
    points(eigen_mu_mpnp4[i] ~ eigen_data_mpnp4$node_age[i,s[j]], pch = 19, col = rgb(0,0,0,0.05))
  }
}

#### interpreting the model ####
# Calculate the 95% credible intervals of the model parameters
round(summary(fit_eigen_cont)$summary[1:3, c(1:4, 8)], 2) # basically no effect whatsoever

# clean environment
rm(adj_tensor,age,beta_age,beta_age2,btwn_cov_mpnp4,btwn_mu_mpnp4,btwn_samples_mpnp4,btwn_samples_mpnp4_std,df_agg_mpnp4,df_long_mpnp4,df_wide_mpnp4,dyad_id,edge_samples_mpnp4,eigen_cov_mpnp4,eigen_data_mpnp4,eigen_mu_mpnp4,eigen_samples_mpnp4,eigen_samples_mpnp4_sim,eigen_samples_mpnp4_std,i,intercept,logit_edge_samples_mpnp4,modeldata_ages,modeldata_ages_sq,mpnp4,mpnp4_ages,mu_btwn,node_ages_mpnp4,num_nodes,num_samples,sd_btwn,bA1,bA2,beta_age_2.5,beta_age_97.5,beta_age_mean,beta_age2_2.5,beta_age2_97.5,beta_age2_mean,fit_eigen_cont,int,intercept_2.5,intercept_97.5,intercept_mean,j,max,mean,min,mu,node_age,params,s,sigma)

########### Time window 3 ###########
#### load mpnp3 nodes, edges and interactions data ####
### import data for aggregated model (binomial)
df_agg_mpnp3 <- read_delim('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_period3_pairwiseevents.csv', delim = ',')
df_agg_mpnp3 <- df_agg_mpnp3[,c(1:16,23:26)] %>% distinct()
df_agg_mpnp3 <- df_agg_mpnp3 %>%
  separate(id_1, into = c('BTF1','num1'), sep = 1, remove = F) %>% 
  separate(id_2, into = c('BTF2','num2'), sep = 1, remove = F) %>% 
  filter(BTF1 == 'B' & BTF2 == 'B')
df_agg_mpnp3 <- df_agg_mpnp3[,c(1:3,6,9:19)]

df_agg_mpnp3$age_mid_round_1 <- floor(df_agg_mpnp3$age_median_1)
df_agg_mpnp3$age_mid_round_2 <- floor(df_agg_mpnp3$age_median_2)
df_agg_mpnp3$node_1_nogaps <- as.integer(as.factor(df_agg_mpnp3$node_1))
df_agg_mpnp3$node_2_nogaps <- as.integer(as.factor(df_agg_mpnp3$node_2))+1
df_agg_mpnp3$dyad_id_nogaps <- as.integer(as.factor(df_agg_mpnp3$dyad))

### load the edge weights
mpnp3 <- readRDS('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp3_edgeweightestimates_mcmcoutput.rds')
mpnp3 <- mpnp3[,c(2:ncol(mpnp3))]
mpnp3 <- mpnp3[, which(colnames(mpnp3) %in% df_agg_mpnp3$dyad)]
logit_edge_samples_mpnp3 <- logit(mpnp3)
print('logit complete')

### load age data
mpnp3_ages <- readRDS('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp3_ageestimates_mcmcoutput.rds')

#### calculate posterior centralities ####
mpnp3 <- as.matrix(logit_edge_samples_mpnp3)
edge_samples_mpnp3 <- plogis(mpnp3)
print('plogis complete')

# Build adjacency tensor
num_nodes <- length(unique(c(df_agg_mpnp3$id_1, df_agg_mpnp3$id_2)))
num_samples <- 4000
adj_tensor <- array(0, c(num_samples, num_nodes, num_nodes),
                    dimnames = list(NULL,
                                    unique(c(df_agg_mpnp3$id_1, df_agg_mpnp3$id_2)),
                                    unique(c(df_agg_mpnp3$id_1, df_agg_mpnp3$id_2))))
for (dyad_id in 1:nrow(df_agg_mpnp3)) {
  dyad_row <- df_agg_mpnp3[df_agg_mpnp3$dyad_id_nogaps == dyad_id, ]
  adj_tensor[, dyad_row$id_1, dyad_row$id_2] <- edge_samples_mpnp3[, dyad_id]
}
rm(dyad_row)

# Calculate centrality and store posterior samples in a matrix
eigen_samples_mpnp3 <- matrix(0, num_samples, num_nodes)
eigen_samples_mpnp3_std <- matrix(0, num_samples, num_nodes)
btwn_samples_mpnp3 <- matrix(0, num_samples, num_nodes)
btwn_samples_mpnp3_std <- matrix(0, num_samples, num_nodes)
for (i in 1:num_samples) {
  g <- graph_from_adjacency_matrix(adj_tensor[i, , ], mode="undirected", weighted=TRUE)
  eigen_samples_mpnp3[i, ] <- eigen_centrality(g)$vector
  eigen_samples_mpnp3_std[i, ] <- (eigen_samples_mpnp3[i, ] - mean(eigen_samples_mpnp3[i, ]))/sd(eigen_samples_mpnp3[i, ])
  btwn_samples_mpnp3[i, ] <- betweenness(g, normalize = T)
  btwn_samples_mpnp3_std[i, ] <- (btwn_samples_mpnp3[i, ] - mean(btwn_samples_mpnp3[i, ]))/sd(btwn_samples_mpnp3[i, ])
  rm(g)
}

# convert betweenness to a Z-score
mu_btwn <- mean(btwn_samples_mpnp3)
sd_btwn <- sd(btwn_samples_mpnp3)
btwn_samples_mpnp3_std <- (btwn_samples_mpnp3 - mu_btwn)/sd_btwn

head(eigen_samples_mpnp3)     # Unstandardised eigenvector centrality
head(eigen_samples_mpnp3_std) # Standardised eigenvector centrality
head(btwn_samples_mpnp3)      # Unstandardised eigenvector centrality
head(btwn_samples_mpnp3_std)  # Standardised eigenvector centrality
print('centralities calculated')

# visualise eigenvector centralities
df_wide_mpnp3 <- data.frame(eigen_samples_mpnp3_std)
colnames(df_wide_mpnp3) <- 1:num_nodes
df_long_mpnp3 <- pivot_longer(df_wide_mpnp3, cols = 1:num_nodes, names_to="node_id", values_to="eigenvector_centrality")
ggplot(df_long_mpnp3, aes(x = eigenvector_centrality)) +
  geom_density(fill="#387780", alpha=0.7, size=0.4) +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) # NOT NICELY NORMAL!

nodes1 <- df_agg_mpnp3 %>% select(id_1, age_mid_round_1) %>% distinct()
nodes2 <- df_agg_mpnp3 %>% select(id_2, age_mid_round_2) %>% distinct()
colnames(nodes1) <- c('id', 'age_mid_round')
colnames(nodes2) <- c('id', 'age_mid_round')
node_ages_mpnp3 <- rbind(nodes1, nodes2) %>% distinct()
rm(nodes1, nodes2)

node_ages_mpnp3$node_id <- as.integer(as.factor(node_ages_mpnp3$id))
df_long_mpnp3$node_id <- as.integer(df_long_mpnp3$node_id)
df_long_mpnp3 <- left_join(df_long_mpnp3, node_ages_mpnp3, by = 'node_id')

ggplot(df_long_mpnp3, aes(x = eigenvector_centrality)) +
  geom_density(aes(fill = as.factor(age_mid_round)), alpha=0.5, size=0.4) +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) # NOT NICELY NORMAL!

# visualise betweenness centrality
df_wide_mpnp3 <- data.frame(btwn_samples_mpnp3_std)
colnames(df_wide_mpnp3) <- 1:num_nodes
df_long_mpnp3 <- pivot_longer(df_wide_mpnp3, cols = 1:num_nodes, names_to="node_id", values_to="btwn_centrality")
ggplot(df_long_mpnp3, aes(x = btwn_centrality)) +
  geom_density(fill="#387780", alpha=0.7, size=0.4) +
  labs(x="betweenness centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

node_ages_mpnp3$node_id <- as.integer(as.factor(node_ages_mpnp3$id))
df_long_mpnp3$node_id <- as.integer(df_long_mpnp3$node_id)
df_long_mpnp3 <- left_join(df_long_mpnp3, node_ages_mpnp3, by = 'node_id')

ggplot(df_long_mpnp3, aes(x = btwn_centrality)) +
  geom_density(aes(fill = as.factor(age_mid_round)), alpha=0.7, size=0.4) +
  labs(x="betweenness centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) # ZERO BETWEENNESS

print('graphs plotted')

#### compute normal approximation ####
# The posterior centralities can now be characterised by the multivariate normal distribution as an approximation to their true posterior distributions. To highlight the importance of sampling from a multivariate normal rather than a univariate normal, we can plot centrality samples from two nodes against each other to see how they co-vary.
plot(eigen_samples_mpnp3_std[, 1], eigen_samples_mpnp3_std[, 2])
plot(btwn_samples_mpnp3_std[, 1], btwn_samples_mpnp3_std[, 2])

plot(eigen_samples_mpnp3 ~ btwn_samples_mpnp3)
plot(eigen_samples_mpnp3_std ~ btwn_samples_mpnp3_std)

# This joint distribution contains important information that we would ideally like to include in the model. The multivariate normal will allow us to do this. Parameterising the multivariate normal approximation is a relatively simple process and can be done by calculating the sample mean and sample covariance matrix of the posterior centrality samples:
eigen_mu_mpnp3 <- apply(eigen_samples_mpnp3_std, 2, mean)
plot(density(eigen_mu_mpnp3)) # double peak -- not normally distributed
eigen_cov_mpnp3 <- cov(eigen_samples_mpnp3_std)
btwn_mu_mpnp3 <- apply(btwn_samples_mpnp3_std, 2, mean) 
btwn_cov_mpnp3 <- cov(btwn_samples_mpnp3_std)

# These quantities will be given to the Stan model as data to model joint posteriors of centrality in the regression. We can run a few quick plots to see how well the approximation is working and the covariance for one of the nodes:
eigen_samples_mpnp3_sim <- MASS::mvrnorm(1e5, eigen_mu_mpnp3, eigen_cov_mpnp3)
plot(density(eigen_samples_mpnp3_std[, 1]), lwd=2, main="Estimated standardised centrality vs normal approximation", xlab="Logit edge weight")
lines(density(eigen_samples_mpnp3_sim[, 1]), col=rgb(0, 0, 1, 0.5), lwd=2)

#### eigenvector prior predictive checks (univariate Gaussian) ####
range(eigen_mu_mpnp3)
age <- 1:60

# eigenvector linear effect only
plot(NULL, xlim = c(0,60), ylim = c(-4,1), las = 1, xlab = 'age', ylab = 'mean eigenvector centrality (std)')
abline(h = range(eigen_mu_mpnp3), lty = 2)
for(i in 1:100){
  intercept <- rnorm(1, -1, 1)
  beta_age <- rnorm(1, 0, 0.02)
  lines(x = age, y = intercept + beta_age*age,
        col = rgb(0,0,1,0.4))
}

# eigenvector allowing for extremes to be most/least central
plot(NULL, xlim = c(0,60), ylim = c(-4,1), las = 1, xlab = 'age', ylab = 'mean eigenvector centrality (std)')
abline(h = range(eigen_mu_mpnp3), lty = 2)
for(i in 1:100){
  intercept <- rnorm(1, -1, 1)
  beta_age <- rnorm(1, 0, 0.02)
  beta_age2 <- rnorm(1, 0, 0.003)
  lines(x = age, y = intercept + beta_age*age + beta_age2*(age^2),
        col = rgb(0,0,1,0.4))
}

#### fit the model ####
modeldata_ages <- as.matrix(mpnp3_ages) %>% t() # CHECK WHY DID I TRANSPOSE THIS???
modeldata_ages_sq <- (modeldata_ages)^2

#### eigenvector ~ age ####
eigen_data_mpnp3 <- list(
  num_nodes = num_nodes,               # Number of dyads
  centrality_mu  = eigen_mu_mpnp3,     # Sample means of logit edge weights
  centrality_cov = eigen_cov_mpnp3,    # Sample covariance of logit edge weights
  node_age = modeldata_ages,           # Age of individual -- are these in the wrong order? should be M1, M10, M100, M101?
  node_age2 = modeldata_ages_sq        # Age squared for quadratic term
)
str(eigen_data_mpnp3)

plot(eigen_data_mpnp3$centrality_mu ~ eigen_data_mpnp3$node_age[,1], las = 1, pch = 19, col = rgb(0,0,1,0.3))

fit_eigen_cont <- sampling(eigen_nodal_cont, data = eigen_data_mpnp3, cores = 4, chains = 4)

#### diagnostics ####
traceplot(fit_eigen_cont) # these look good

params <- rstan::extract(fit_eigen_cont)
plot(density(eigen_samples_mpnp3_std[1, ]), main = "Posterior predictive density (standardised centrality)",
     col = rgb(0, 0, 0, 0.25), ylim = c(0, 0.8), las = 1)
for (i in 1:100) {
  j <- sample(1:num_samples, 1)                                         # select a network to plot
  lines(density(eigen_samples_mpnp3_std[j, ]), col=rgb(0, 0, 0, 0.25)) # plot centrality density for network j (black)
  mu <- params$beta_age[j]*apply(eigen_data_mpnp3$node_age,1,mean)     # extract age slope parameter for network j
  sigma <- eigen_cov_mpnp3 + diag(rep(params$sigma[j], num_nodes))     # extract variance parameter for network j
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))  # plot predicted centralities for network j (blue)
} # not too bad -- this makes sense considering the main data were left skewed but the model is normal

### plot age against predictions
node_age <- 1:60
beta_age_mean  <- summary(fit_eigen_cont)$summary[1,1]
beta_age_2.5   <- summary(fit_eigen_cont)$summary[1,4]
beta_age_97.5  <- summary(fit_eigen_cont)$summary[1,8]
beta_age2_mean <- summary(fit_eigen_cont)$summary[2,1]
beta_age2_2.5  <- summary(fit_eigen_cont)$summary[2,4]
beta_age2_97.5 <- summary(fit_eigen_cont)$summary[2,8]
intercept_mean <- summary(fit_eigen_cont)$summary[3,1]
intercept_2.5  <- summary(fit_eigen_cont)$summary[3,4]
intercept_97.5 <- summary(fit_eigen_cont)$summary[3,8]
mu <- intercept_mean  + beta_age_mean*node_age + beta_age2_mean*(node_age^2)
min <- intercept_2.5  + beta_age_2.5*node_age  + beta_age2_2.5*(node_age^2)
max <- intercept_97.5 + beta_age_97.5*node_age + beta_age2_97.5*(node_age^2)
summary(mu)
summary(min)
summary(max)
plot(NULL, xlim = c(0,60), ylim = c(-5, 5), lwd = 2, las = 1,
     xlab = 'age (years)', ylab = 'eigenvector centrality (std)')
for(i in 1:100){
  j <- sample(1:4000, 1)
  int <- params$intercept[j]
  bA1 <- params$beta_age[j]
  bA2 <- params$beta_age2[j]
  mean <- int + bA1*age + bA2*age
  lines(c(1:60), mean, col = rgb(0.5,0,1,0.2))
}
lines(mu ~ node_age, lwd = 2) 
lines(min ~ node_age, col = 'red')
lines(max ~ node_age, col = 'red')
abline(h = c(-2,2), lty = 2)
for(i in 1:num_nodes){
  s <- sample(x = 1:num_samples, size = 100, replace = F)
  for(j in 1:length(s)){
    points(eigen_mu_mpnp3[i] ~ eigen_data_mpnp3$node_age[i,s[j]], pch = 19, col = rgb(0,0,0,0.05))
  }
}

#### interpreting the model ####
# Calculate the 95% credible intervals of the model parameters
round(summary(fit_eigen_cont)$summary[1:3, c(1:4, 8)], 2) # basically no effect whatsoever

# clean environment
rm(adj_tensor,age,beta_age,beta_age2,btwn_cov_mpnp3,btwn_mu_mpnp3,btwn_samples_mpnp3,btwn_samples_mpnp3_std,df_agg_mpnp3,df_long_mpnp3,df_wide_mpnp3,dyad_id,edge_samples_mpnp3,eigen_cov_mpnp3,eigen_data_mpnp3,eigen_mu_mpnp3,eigen_samples_mpnp3,eigen_samples_mpnp3_sim,eigen_samples_mpnp3_std,i,intercept,logit_edge_samples_mpnp3,modeldata_ages,modeldata_ages_sq,mpnp3,mpnp3_ages,mu_btwn,node_ages_mpnp3,num_nodes,num_samples,sd_btwn,bA1,bA2,beta_age_2.5,beta_age_97.5,beta_age_mean,beta_age2_2.5,beta_age2_97.5,beta_age2_mean,fit_eigen_cont,int,intercept_2.5,intercept_97.5,intercept_mean,j,max,mean,min,mu,node_age,params,s,sigma)

########### Time window 2 ###########
#### load mpnp2 nodes, edges and interactions data ####
### import data for aggregated model (binomial)
df_agg_mpnp2 <- read_delim('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_period2_pairwiseevents.csv', delim = ',')
df_agg_mpnp2 <- df_agg_mpnp2[,c(1:13,18:21)] %>% distinct()
df_agg_mpnp2 <- df_agg_mpnp2 %>%
  separate(id_1, into = c('BTF1','num1'), sep = 1, remove = F) %>% 
  separate(id_2, into = c('BTF2','num2'), sep = 1, remove = F) %>% 
  filter(BTF1 == 'B' & BTF2 == 'B')
df_agg_mpnp2 <- df_agg_mpnp2[,c(1:3,6,9:19)]

df_agg_mpnp2$age_mid_round_1 <- floor(df_agg_mpnp2$age_median_1)
df_agg_mpnp2$age_mid_round_2 <- floor(df_agg_mpnp2$age_median_2)
df_agg_mpnp2$node_1_nogaps <- as.integer(as.factor(df_agg_mpnp2$node_1))
df_agg_mpnp2$node_2_nogaps <- as.integer(as.factor(df_agg_mpnp2$node_2))+1
df_agg_mpnp2$dyad_id_nogaps <- as.integer(as.factor(df_agg_mpnp2$dyad))

### load the edge weights
mpnp2 <- readRDS('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp2_edgeweightestimates_mcmcoutput.rds')
mpnp2 <- mpnp2[,c(2:ncol(mpnp2))]
mpnp2 <- mpnp2[, which(colnames(mpnp2) %in% df_agg_mpnp2$dyad)]
logit_edge_samples_mpnp2 <- logit(mpnp2)
print('logit complete')

### load age data
mpnp2_ages <- readRDS('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp2_ageestimates_mcmcoutput.rds')

#### calculate posterior centralities ####
mpnp2 <- as.matrix(logit_edge_samples_mpnp2)
edge_samples_mpnp2 <- plogis(mpnp2)
print('plogis complete')

# Build adjacency tensor
num_nodes <- length(unique(c(df_agg_mpnp2$id_1, df_agg_mpnp2$id_2)))
num_samples <- 4000
adj_tensor <- array(0, c(num_samples, num_nodes, num_nodes),
                    dimnames = list(NULL,
                                    unique(c(df_agg_mpnp2$id_1, df_agg_mpnp2$id_2)),
                                    unique(c(df_agg_mpnp2$id_1, df_agg_mpnp2$id_2))))
for (dyad_id in 1:nrow(df_agg_mpnp2)) {
  dyad_row <- df_agg_mpnp2[df_agg_mpnp2$dyad_id_nogaps == dyad_id, ]
  adj_tensor[, dyad_row$id_1, dyad_row$id_2] <- edge_samples_mpnp2[, dyad_id]
}
rm(dyad_row)

# Calculate centrality and store posterior samples in a matrix
eigen_samples_mpnp2 <- matrix(0, num_samples, num_nodes)
eigen_samples_mpnp2_std <- matrix(0, num_samples, num_nodes)
btwn_samples_mpnp2 <- matrix(0, num_samples, num_nodes)
btwn_samples_mpnp2_std <- matrix(0, num_samples, num_nodes)
for (i in 1:num_samples) {
  g <- graph_from_adjacency_matrix(adj_tensor[i, , ], mode="undirected", weighted=TRUE)
  eigen_samples_mpnp2[i, ] <- eigen_centrality(g)$vector
  eigen_samples_mpnp2_std[i, ] <- (eigen_samples_mpnp2[i, ] - mean(eigen_samples_mpnp2[i, ]))/sd(eigen_samples_mpnp2[i, ])
  btwn_samples_mpnp2[i, ] <- betweenness(g, normalize = T)
  btwn_samples_mpnp2_std[i, ] <- (btwn_samples_mpnp2[i, ] - mean(btwn_samples_mpnp2[i, ]))/sd(btwn_samples_mpnp2[i, ])
  rm(g)
}

# convert betweenness to a Z-score
mu_btwn <- mean(btwn_samples_mpnp2)
sd_btwn <- sd(btwn_samples_mpnp2)
btwn_samples_mpnp2_std <- (btwn_samples_mpnp2 - mu_btwn)/sd_btwn

head(eigen_samples_mpnp2)     # Unstandardised eigenvector centrality
head(eigen_samples_mpnp2_std) # Standardised eigenvector centrality
head(btwn_samples_mpnp2)      # Unstandardised eigenvector centrality
head(btwn_samples_mpnp2_std)  # Standardised eigenvector centrality
print('centralities calculated')

# visualise eigenvector centralities
df_wide_mpnp2 <- data.frame(eigen_samples_mpnp2_std)
colnames(df_wide_mpnp2) <- 1:num_nodes
df_long_mpnp2 <- pivot_longer(df_wide_mpnp2, cols = 1:num_nodes, names_to="node_id", values_to="eigenvector_centrality")
ggplot(df_long_mpnp2, aes(x = eigenvector_centrality)) +
  geom_density(fill="#387780", alpha=0.7, size=0.4) +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) # NOT NICELY NORMAL!

nodes1 <- df_agg_mpnp2 %>% select(id_1, age_mid_round_1) %>% distinct()
nodes2 <- df_agg_mpnp2 %>% select(id_2, age_mid_round_2) %>% distinct()
colnames(nodes1) <- c('id', 'age_mid_round')
colnames(nodes2) <- c('id', 'age_mid_round')
node_ages_mpnp2 <- rbind(nodes1, nodes2) %>% distinct()
rm(nodes1, nodes2)

node_ages_mpnp2$node_id <- as.integer(as.factor(node_ages_mpnp2$id))
df_long_mpnp2$node_id <- as.integer(df_long_mpnp2$node_id)
df_long_mpnp2 <- left_join(df_long_mpnp2, node_ages_mpnp2, by = 'node_id')

ggplot(df_long_mpnp2, aes(x = eigenvector_centrality)) +
  geom_density(aes(fill = as.factor(age_mid_round)), alpha=0.5, size=0.4) +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) # NOT NICELY NORMAL!

# visualise betweenness centrality
df_wide_mpnp2 <- data.frame(btwn_samples_mpnp2_std)
colnames(df_wide_mpnp2) <- 1:num_nodes
df_long_mpnp2 <- pivot_longer(df_wide_mpnp2, cols = 1:num_nodes, names_to="node_id", values_to="btwn_centrality")
ggplot(df_long_mpnp2, aes(x = btwn_centrality)) +
  geom_density(fill="#387780", alpha=0.7, size=0.4) +
  labs(x="betweenness centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

node_ages_mpnp2$node_id <- as.integer(as.factor(node_ages_mpnp2$id))
df_long_mpnp2$node_id <- as.integer(df_long_mpnp2$node_id)
df_long_mpnp2 <- left_join(df_long_mpnp2, node_ages_mpnp2, by = 'node_id')

ggplot(df_long_mpnp2, aes(x = btwn_centrality)) +
  geom_density(aes(fill = as.factor(age_mid_round)), alpha=0.7, size=0.4) +
  labs(x="betweenness centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) # ZERO BETWEENNESS

print('graphs plotted')

#### compute normal approximation ####
# The posterior centralities can now be characterised by the multivariate normal distribution as an approximation to their true posterior distributions. To highlight the importance of sampling from a multivariate normal rather than a univariate normal, we can plot centrality samples from two nodes against each other to see how they co-vary.
plot(eigen_samples_mpnp2_std[, 1], eigen_samples_mpnp2_std[, 2])
plot(btwn_samples_mpnp2_std[, 1], btwn_samples_mpnp2_std[, 2])

plot(eigen_samples_mpnp2 ~ btwn_samples_mpnp2)
plot(eigen_samples_mpnp2_std ~ btwn_samples_mpnp2_std)

# This joint distribution contains important information that we would ideally like to include in the model. The multivariate normal will allow us to do this. Parameterising the multivariate normal approximation is a relatively simple process and can be done by calculating the sample mean and sample covariance matrix of the posterior centrality samples:
eigen_mu_mpnp2 <- apply(eigen_samples_mpnp2_std, 2, mean)
plot(density(eigen_mu_mpnp2)) # double peak -- not normally distributed
eigen_cov_mpnp2 <- cov(eigen_samples_mpnp2_std)
btwn_mu_mpnp2 <- apply(btwn_samples_mpnp2_std, 2, mean) 
btwn_cov_mpnp2 <- cov(btwn_samples_mpnp2_std)

# These quantities will be given to the Stan model as data to model joint posteriors of centrality in the regression. We can run a few quick plots to see how well the approximation is working and the covariance for one of the nodes:
eigen_samples_mpnp2_sim <- MASS::mvrnorm(1e5, eigen_mu_mpnp2, eigen_cov_mpnp2)
plot(density(eigen_samples_mpnp2_std[, 1]), lwd=2, main="Estimated standardised centrality vs normal approximation", xlab="Logit edge weight")
lines(density(eigen_samples_mpnp2_sim[, 1]), col=rgb(0, 0, 1, 0.5), lwd=2)

#### eigenvector prior predictive checks (univariate Gaussian) ####
range(eigen_mu_mpnp2)
age <- 1:60

# eigenvector linear effect only
plot(NULL, xlim = c(0,60), ylim = c(-4,1), las = 1, xlab = 'age', ylab = 'mean eigenvector centrality (std)')
abline(h = range(eigen_mu_mpnp2), lty = 2)
for(i in 1:100){
  intercept <- rnorm(1, -1, 1)
  beta_age <- rnorm(1, 0, 0.02)
  lines(x = age, y = intercept + beta_age*age,
        col = rgb(0,0,1,0.4))
}

# eigenvector allowing for extremes to be most/least central
plot(NULL, xlim = c(0,60), ylim = c(-4,1), las = 1, xlab = 'age', ylab = 'mean eigenvector centrality (std)')
abline(h = range(eigen_mu_mpnp2), lty = 2)
for(i in 1:100){
  intercept <- rnorm(1, -1, 1)
  beta_age <- rnorm(1, 0, 0.02)
  beta_age2 <- rnorm(1, 0, 0.003)
  lines(x = age, y = intercept + beta_age*age + beta_age2*(age^2),
        col = rgb(0,0,1,0.4))
}

#### fit the model ####
modeldata_ages <- as.matrix(mpnp2_ages) %>% t() # CHECK WHY DID I TRANSPOSE THIS???
modeldata_ages_sq <- (modeldata_ages)^2

#### eigenvector ~ age ####
eigen_data_mpnp2 <- list(
  num_nodes = num_nodes,               # Number of dyads
  centrality_mu  = eigen_mu_mpnp2,     # Sample means of logit edge weights
  centrality_cov = eigen_cov_mpnp2,    # Sample covariance of logit edge weights
  node_age = modeldata_ages,           # Age of individual -- are these in the wrong order? should be M1, M10, M100, M101?
  node_age2 = modeldata_ages_sq        # Age squared for quadratic term
)
str(eigen_data_mpnp2)

plot(eigen_data_mpnp2$centrality_mu ~ eigen_data_mpnp2$node_age[,1], las = 1, pch = 19, col = rgb(0,0,1,0.3))

fit_eigen_cont <- sampling(eigen_nodal_cont, data = eigen_data_mpnp2, cores = 4, chains = 4)

#### diagnostics ####
traceplot(fit_eigen_cont) # these look good

params <- rstan::extract(fit_eigen_cont)
plot(density(eigen_samples_mpnp2_std[1, ]), main = "Posterior predictive density (standardised centrality)",
     col = rgb(0, 0, 0, 0.25), ylim = c(0, 0.8), las = 1)
for (i in 1:100) {
  j <- sample(1:num_samples, 1)                                         # select a network to plot
  lines(density(eigen_samples_mpnp2_std[j, ]), col=rgb(0, 0, 0, 0.25)) # plot centrality density for network j (black)
  mu <- params$beta_age[j]*apply(eigen_data_mpnp2$node_age,1,mean)     # extract age slope parameter for network j
  sigma <- eigen_cov_mpnp2 + diag(rep(params$sigma[j], num_nodes))     # extract variance parameter for network j
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))  # plot predicted centralities for network j (blue)
} # not too bad -- this makes sense considering the main data were left skewed but the model is normal

### plot age against predictions
node_age <- 1:60
beta_age_mean  <- summary(fit_eigen_cont)$summary[1,1]
beta_age_2.5   <- summary(fit_eigen_cont)$summary[1,4]
beta_age_97.5  <- summary(fit_eigen_cont)$summary[1,8]
beta_age2_mean <- summary(fit_eigen_cont)$summary[2,1]
beta_age2_2.5  <- summary(fit_eigen_cont)$summary[2,4]
beta_age2_97.5 <- summary(fit_eigen_cont)$summary[2,8]
intercept_mean <- summary(fit_eigen_cont)$summary[3,1]
intercept_2.5  <- summary(fit_eigen_cont)$summary[3,4]
intercept_97.5 <- summary(fit_eigen_cont)$summary[3,8]
mu <- intercept_mean  + beta_age_mean*node_age + beta_age2_mean*(node_age^2)
min <- intercept_2.5  + beta_age_2.5*node_age  + beta_age2_2.5*(node_age^2)
max <- intercept_97.5 + beta_age_97.5*node_age + beta_age2_97.5*(node_age^2)
summary(mu)
summary(min)
summary(max)
plot(NULL, xlim = c(0,60), ylim = c(-5, 5), lwd = 2, las = 1,
     xlab = 'age (years)', ylab = 'eigenvector centrality (std)')
for(i in 1:100){
  j <- sample(1:4000, 1)
  int <- params$intercept[j]
  bA1 <- params$beta_age[j]
  bA2 <- params$beta_age2[j]
  mean <- int + bA1*age + bA2*age
  lines(c(1:60), mean, col = rgb(0.5,0,1,0.2))
}
lines(mu ~ node_age, lwd = 2) 
lines(min ~ node_age, col = 'red')
lines(max ~ node_age, col = 'red')
abline(h = c(-2,2), lty = 2)
for(i in 1:num_nodes){
  s <- sample(x = 1:num_samples, size = 100, replace = F)
  for(j in 1:length(s)){
    points(eigen_mu_mpnp2[i] ~ eigen_data_mpnp2$node_age[i,s[j]], pch = 19, col = rgb(0,0,0,0.05))
  }
}

#### interpreting the model ####
# Calculate the 95% credible intervals of the model parameters
round(summary(fit_eigen_cont)$summary[1:3, c(1:4, 8)], 2) # basically no effect whatsoever

########### Time window 1 ###########
#### load mpnp1 nodes, edges and interactions data ####
### import data for aggregated model (binomial)
df_agg_mpnp1 <- read_delim('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_period1_pairwiseevents.csv', delim = ',')
df_agg_mpnp1 <- df_agg_mpnp1[,c(1:13,18:21)] %>% distinct()
df_agg_mpnp1 <- df_agg_mpnp1 %>%
  separate(id_1, into = c('BTF1','num1'), sep = 1, remove = F) %>% 
  separate(id_2, into = c('BTF2','num2'), sep = 1, remove = F) %>% 
  filter(BTF1 == 'B' & BTF2 == 'B')
df_agg_mpnp1 <- df_agg_mpnp1[,c(1:3,6,9:19)]

df_agg_mpnp1$age_mid_round_1 <- floor(df_agg_mpnp1$age_median_1)
df_agg_mpnp1$age_mid_round_2 <- floor(df_agg_mpnp1$age_median_2)
df_agg_mpnp1$node_1_nogaps <- as.integer(as.factor(df_agg_mpnp1$node_1))
df_agg_mpnp1$node_2_nogaps <- as.integer(as.factor(df_agg_mpnp1$node_2))+1
df_agg_mpnp1$dyad_id_nogaps <- as.integer(as.factor(df_agg_mpnp1$dyad))

### load the edge weights
mpnp1 <- readRDS('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp1_edgeweightestimates_mcmcoutput.rds')
mpnp1 <- mpnp1[,c(2:ncol(mpnp1))]
mpnp1 <- mpnp1[, which(colnames(mpnp1) %in% df_agg_mpnp1$dyad)]
logit_edge_samples_mpnp1 <- logit(mpnp1)
print('logit complete')

### load age data
mpnp1_ages <- readRDS('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp1_ageestimates_mcmcoutput.rds')

#### calculate posterior centralities ####
mpnp1 <- as.matrix(logit_edge_samples_mpnp1)
edge_samples_mpnp1 <- plogis(mpnp1)
print('plogis complete')

# Build adjacency tensor
num_nodes <- length(unique(c(df_agg_mpnp1$id_1, df_agg_mpnp1$id_2)))
num_samples <- 4000
adj_tensor <- array(0, c(num_samples, num_nodes, num_nodes),
                    dimnames = list(NULL,
                                    unique(c(df_agg_mpnp1$id_1, df_agg_mpnp1$id_2)),
                                    unique(c(df_agg_mpnp1$id_1, df_agg_mpnp1$id_2))))
for (dyad_id in 1:nrow(df_agg_mpnp1)) {
  dyad_row <- df_agg_mpnp1[df_agg_mpnp1$dyad_id_nogaps == dyad_id, ]
  adj_tensor[, dyad_row$id_1, dyad_row$id_2] <- edge_samples_mpnp1[, dyad_id]
}
rm(dyad_row)

# Calculate centrality and store posterior samples in a matrix
eigen_samples_mpnp1 <- matrix(0, num_samples, num_nodes)
eigen_samples_mpnp1_std <- matrix(0, num_samples, num_nodes)
btwn_samples_mpnp1 <- matrix(0, num_samples, num_nodes)
btwn_samples_mpnp1_std <- matrix(0, num_samples, num_nodes)
for (i in 1:num_samples) {
  g <- graph_from_adjacency_matrix(adj_tensor[i, , ], mode="undirected", weighted=TRUE)
  eigen_samples_mpnp1[i, ] <- eigen_centrality(g)$vector
  eigen_samples_mpnp1_std[i, ] <- (eigen_samples_mpnp1[i, ] - mean(eigen_samples_mpnp1[i, ]))/sd(eigen_samples_mpnp1[i, ])
  btwn_samples_mpnp1[i, ] <- betweenness(g, normalize = T)
  btwn_samples_mpnp1_std[i, ] <- (btwn_samples_mpnp1[i, ] - mean(btwn_samples_mpnp1[i, ]))/sd(btwn_samples_mpnp1[i, ])
  rm(g)
}

# convert betweenness to a Z-score
mu_btwn <- mean(btwn_samples_mpnp1)
sd_btwn <- sd(btwn_samples_mpnp1)
btwn_samples_mpnp1_std <- (btwn_samples_mpnp1 - mu_btwn)/sd_btwn

head(eigen_samples_mpnp1)     # Unstandardised eigenvector centrality
head(eigen_samples_mpnp1_std) # Standardised eigenvector centrality
head(btwn_samples_mpnp1)      # Unstandardised eigenvector centrality
head(btwn_samples_mpnp1_std)  # Standardised eigenvector centrality
print('centralities calculated')

# visualise eigenvector centralities
df_wide_mpnp1 <- data.frame(eigen_samples_mpnp1_std)
colnames(df_wide_mpnp1) <- 1:num_nodes
df_long_mpnp1 <- pivot_longer(df_wide_mpnp1, cols = 1:num_nodes, names_to="node_id", values_to="eigenvector_centrality")
ggplot(df_long_mpnp1, aes(x = eigenvector_centrality)) +
  geom_density(fill="#387780", alpha=0.7, size=0.4) +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) # NOT NICELY NORMAL!

nodes1 <- df_agg_mpnp1 %>% select(id_1, age_mid_round_1) %>% distinct()
nodes2 <- df_agg_mpnp1 %>% select(id_2, age_mid_round_2) %>% distinct()
colnames(nodes1) <- c('id', 'age_mid_round')
colnames(nodes2) <- c('id', 'age_mid_round')
node_ages_mpnp1 <- rbind(nodes1, nodes2) %>% distinct()
rm(nodes1, nodes2)

node_ages_mpnp1$node_id <- as.integer(as.factor(node_ages_mpnp1$id))
df_long_mpnp1$node_id <- as.integer(df_long_mpnp1$node_id)
df_long_mpnp1 <- left_join(df_long_mpnp1, node_ages_mpnp1, by = 'node_id')

ggplot(df_long_mpnp1, aes(x = eigenvector_centrality)) +
  geom_density(aes(fill = as.factor(age_mid_round)), alpha=0.5, size=0.4) +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) # NOT NICELY NORMAL!

# visualise betweenness centrality
df_wide_mpnp1 <- data.frame(btwn_samples_mpnp1_std)
colnames(df_wide_mpnp1) <- 1:num_nodes
df_long_mpnp1 <- pivot_longer(df_wide_mpnp1, cols = 1:num_nodes, names_to="node_id", values_to="btwn_centrality")
ggplot(df_long_mpnp1, aes(x = btwn_centrality)) +
  geom_density(fill="#387780", alpha=0.7, size=0.4) +
  labs(x="betweenness centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

node_ages_mpnp1$node_id <- as.integer(as.factor(node_ages_mpnp1$id))
df_long_mpnp1$node_id <- as.integer(df_long_mpnp1$node_id)
df_long_mpnp1 <- left_join(df_long_mpnp1, node_ages_mpnp1, by = 'node_id')

ggplot(df_long_mpnp1, aes(x = btwn_centrality)) +
  geom_density(aes(fill = as.factor(age_mid_round)), alpha=0.7, size=0.4) +
  labs(x="betweenness centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) # ZERO BETWEENNESS

print('graphs plotted')

#### compute normal approximation ####
# The posterior centralities can now be characterised by the multivariate normal distribution as an approximation to their true posterior distributions. To highlight the importance of sampling from a multivariate normal rather than a univariate normal, we can plot centrality samples from two nodes against each other to see how they co-vary.
plot(eigen_samples_mpnp1_std[, 1], eigen_samples_mpnp1_std[, 2])
plot(btwn_samples_mpnp1_std[, 1], btwn_samples_mpnp1_std[, 2])

plot(eigen_samples_mpnp1 ~ btwn_samples_mpnp1)
plot(eigen_samples_mpnp1_std ~ btwn_samples_mpnp1_std)

# This joint distribution contains important information that we would ideally like to include in the model. The multivariate normal will allow us to do this. Parameterising the multivariate normal approximation is a relatively simple process and can be done by calculating the sample mean and sample covariance matrix of the posterior centrality samples:
eigen_mu_mpnp1 <- apply(eigen_samples_mpnp1_std, 2, mean)
plot(density(eigen_mu_mpnp1)) # double peak -- not normally distributed
eigen_cov_mpnp1 <- cov(eigen_samples_mpnp1_std)
btwn_mu_mpnp1 <- apply(btwn_samples_mpnp1_std, 2, mean) 
btwn_cov_mpnp1 <- cov(btwn_samples_mpnp1_std)

# These quantities will be given to the Stan model as data to model joint posteriors of centrality in the regression. We can run a few quick plots to see how well the approximation is working and the covariance for one of the nodes:
eigen_samples_mpnp1_sim <- MASS::mvrnorm(1e5, eigen_mu_mpnp1, eigen_cov_mpnp1)
plot(density(eigen_samples_mpnp1_std[, 1]), lwd=2, main="Estimated standardised centrality vs normal approximation", xlab="Logit edge weight")
lines(density(eigen_samples_mpnp1_sim[, 1]), col=rgb(0, 0, 1, 0.5), lwd=2)

#### eigenvector prior predictive checks (univariate Gaussian) ####
range(eigen_mu_mpnp1)
age <- 1:60

# eigenvector linear effect only
plot(NULL, xlim = c(0,60), ylim = c(-4,1), las = 1, xlab = 'age', ylab = 'mean eigenvector centrality (std)')
abline(h = range(eigen_mu_mpnp1), lty = 2)
for(i in 1:100){
  intercept <- rnorm(1, -1, 1)
  beta_age <- rnorm(1, 0, 0.02)
  lines(x = age, y = intercept + beta_age*age,
        col = rgb(0,0,1,0.4))
}

# eigenvector allowing for extremes to be most/least central
plot(NULL, xlim = c(0,60), ylim = c(-4,1), las = 1, xlab = 'age', ylab = 'mean eigenvector centrality (std)')
abline(h = range(eigen_mu_mpnp1), lty = 2)
for(i in 1:100){
  intercept <- rnorm(1, -1, 1)
  beta_age <- rnorm(1, 0, 0.02)
  beta_age2 <- rnorm(1, 0, 0.003)
  lines(x = age, y = intercept + beta_age*age + beta_age2*(age^2),
        col = rgb(0,0,1,0.4))
}

#### fit the model ####
modeldata_ages <- as.matrix(mpnp1_ages) %>% t() # CHECK WHY DID I TRANSPOSE THIS???
modeldata_ages_sq <- (modeldata_ages)^2

#### eigenvector ~ age ####
eigen_data_mpnp1 <- list(
  num_nodes = num_nodes,               # Number of dyads
  centrality_mu  = eigen_mu_mpnp1,     # Sample means of logit edge weights
  centrality_cov = eigen_cov_mpnp1,    # Sample covariance of logit edge weights
  node_age = modeldata_ages,           # Age of individual -- are these in the wrong order? should be M1, M10, M100, M101?
  node_age2 = modeldata_ages_sq        # Age squared for quadratic term
)
str(eigen_data_mpnp1)

plot(eigen_data_mpnp1$centrality_mu ~ eigen_data_mpnp1$node_age[,1], las = 1, pch = 19, col = rgb(0,0,1,0.3))

fit_eigen_cont <- sampling(eigen_nodal_cont, data = eigen_data_mpnp1, cores = 4, chains = 4)

#### diagnostics ####
traceplot(fit_eigen_cont) # these look good

params <- rstan::extract(fit_eigen_cont)
plot(density(eigen_samples_mpnp1_std[1, ]), main = "Posterior predictive density (standardised centrality)",
     col = rgb(0, 0, 0, 0.25), ylim = c(0, 0.8), las = 1)
for (i in 1:100) {
  j <- sample(1:num_samples, 1)                                         # select a network to plot
  lines(density(eigen_samples_mpnp1_std[j, ]), col=rgb(0, 0, 0, 0.25)) # plot centrality density for network j (black)
  mu <- params$beta_age[j]*apply(eigen_data_mpnp1$node_age,1,mean)     # extract age slope parameter for network j
  sigma <- eigen_cov_mpnp1 + diag(rep(params$sigma[j], num_nodes))     # extract variance parameter for network j
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))  # plot predicted centralities for network j (blue)
} # not too bad -- this makes sense considering the main data were left skewed but the model is normal

### plot age against predictions
node_age <- 1:60
beta_age_mean  <- summary(fit_eigen_cont)$summary[1,1]
beta_age_2.5   <- summary(fit_eigen_cont)$summary[1,4]
beta_age_97.5  <- summary(fit_eigen_cont)$summary[1,8]
beta_age2_mean <- summary(fit_eigen_cont)$summary[2,1]
beta_age2_2.5  <- summary(fit_eigen_cont)$summary[2,4]
beta_age2_97.5 <- summary(fit_eigen_cont)$summary[2,8]
intercept_mean <- summary(fit_eigen_cont)$summary[3,1]
intercept_2.5  <- summary(fit_eigen_cont)$summary[3,4]
intercept_97.5 <- summary(fit_eigen_cont)$summary[3,8]
mu <- intercept_mean  + beta_age_mean*node_age + beta_age2_mean*(node_age^2)
min <- intercept_2.5  + beta_age_2.5*node_age  + beta_age2_2.5*(node_age^2)
max <- intercept_97.5 + beta_age_97.5*node_age + beta_age2_97.5*(node_age^2)
summary(mu)
summary(min)
summary(max)
plot(NULL, xlim = c(0,60), ylim = c(-5, 5), lwd = 2, las = 1,
     xlab = 'age (years)', ylab = 'eigenvector centrality (std)')
for(i in 1:100){
  j <- sample(1:4000, 1)
  int <- params$intercept[j]
  bA1 <- params$beta_age[j]
  bA2 <- params$beta_age2[j]
  mean <- int + bA1*age + bA2*age
  lines(c(1:60), mean, col = rgb(0.5,0,1,0.2))
}
lines(mu ~ node_age, lwd = 2) 
lines(min ~ node_age, col = 'red')
lines(max ~ node_age, col = 'red')
abline(h = c(-2,2), lty = 2)
for(i in 1:num_nodes){
  s <- sample(x = 1:num_samples, size = 100, replace = F)
  for(j in 1:length(s)){
    points(eigen_mu_mpnp1[i] ~ eigen_data_mpnp1$node_age[i,s[j]], pch = 19, col = rgb(0,0,0,0.05))
  }
}

#### interpreting the model ####
# Calculate the 95% credible intervals of the model parameters
round(summary(fit_eigen_cont)$summary[1:3, c(1:4, 8)], 2) # basically no effect whatsoever
