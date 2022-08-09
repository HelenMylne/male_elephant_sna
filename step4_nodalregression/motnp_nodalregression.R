#### information #####
# 1) to estimate the ages of individuals in other populations, first came up with a survival curve for the Amboseli elephants which have good birth and death data.
# 2) selected the Gompertz bathtub distribution
# 3) use distribution to convert MOTNP ages in categories to probability distributions of true age (first part of this script)
# 4) extract node centrality distribution from draws of edge weight
# 5) use MOTNP age distributions to predict network centrality (third part of this script)

#### set up ####
library(tidyverse)
library(cmdstanr)
library(ggdist)
library(posterior)
library(bayesplot)
library(rstan)
library(igraph)
library(LaplacesDemon)

#library(tidyverse, lib.loc = 'packages/')
#library(cmdstanr, lib.loc = 'packages/')
#library(ggdist, lib.loc = 'packages/')
#library(posterior, lib.loc = 'packages/')
#library(bayesplot, lib.loc = 'packages/')
#library(rstan, lib.loc = 'packages/')
#library(igraph, lib.loc = 'packages/')
#library(LaplacesDemon, lib.loc = 'packages/')

# define PDF output
#pdf('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/outputs/motnp_nodalregression_plots.pdf')

############ Eigenvector ############
#### load MOTNP nodes, edges and interactions data ####
motnp_males <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_elenodes.csv') %>% 
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

### import data for aggregated model (binomial)
df_agg_motnp <- read_delim('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_bayesian_binomialpairwiseevents.csv', delim = ',') %>% 
  filter(dem_class_1 == 'AM' | dem_class_1 == 'PM') %>% 
  filter(dem_class_2 == 'AM' | dem_class_2 == 'PM')
df_agg_motnp$sex_1 <- 'M'
df_agg_motnp$age_cat_id_1 <- ifelse(df_agg_motnp$age_category_1 == '9-10', 2,
                                    ifelse(df_agg_motnp$age_category_1 == '10-15', 3,
                                           ifelse(df_agg_motnp$age_category_1 == '15-19', 4,
                                                  ifelse(df_agg_motnp$age_category_1 == '20-25', 5,
                                                         ifelse(df_agg_motnp$age_category_1 == '25-40', 6,
                                                                ifelse(df_agg_motnp$age_category_1 == '40+', 7, 1))))))
df_agg_motnp$age_class_1 <- ifelse(df_agg_motnp$age_cat_id_1 == 2, 'Juvenile',
                                   ifelse(df_agg_motnp$age_cat_id_1 > 4, 'Adult','Pubescent'))
df_agg_motnp$age_cat_id_2 <- ifelse(df_agg_motnp$age_category_2 == '9-10', 2,
                                    ifelse(df_agg_motnp$age_category_2 == '10-15', 3,
                                           ifelse(df_agg_motnp$age_category_2 == '15-19', 4,
                                                  ifelse(df_agg_motnp$age_category_2 == '20-25', 5,
                                                         ifelse(df_agg_motnp$age_category_2 == '25-40', 6,
                                                                ifelse(df_agg_motnp$age_category_2 == '40+', 7, 1))))))
df_agg_motnp$age_class_2 <- ifelse(df_agg_motnp$age_cat_id_2 == 2, 'Juvenile',
                                   ifelse(df_agg_motnp$age_cat_id_2 > 4, 'Adult','Pubescent'))
table(df_agg_motnp$age_cat_id_1,df_agg_motnp$age_category_1)
table(df_agg_motnp$age_cat_id_2,df_agg_motnp$age_category_2)

df_agg_motnp$dem_class_1 <- ifelse(df_agg_motnp$age_class_1 == 'Adult', 'AM',
                                   ifelse(df_agg_motnp$age_class_1 == 'Pubescent', 'PM', 'JM'))
df_agg_motnp$dem_class_2 <- ifelse(df_agg_motnp$age_class_2 == 'Adult', 'AM',
                                   ifelse(df_agg_motnp$age_class_2 == 'Pubescent', 'PM', 'JM'))
df_agg_motnp$dem_type <- ifelse(df_agg_motnp$age_cat_id_1 >= df_agg_motnp$age_cat_id_2,
                                paste(df_agg_motnp$dem_class_1, df_agg_motnp$dem_class_2, sep = '_'),
                                paste(df_agg_motnp$dem_class_2, df_agg_motnp$dem_class_1, sep = '_'))
df_agg_motnp$count_dyad <- (df_agg_motnp$count_1 + df_agg_motnp$count_2) - df_agg_motnp$event_count  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.
df_agg_motnp$node_1_nogaps <- as.integer(as.factor(df_agg_motnp$node_1))
df_agg_motnp$node_2_nogaps <- as.integer(as.factor(df_agg_motnp$node_2))+1
df_agg_motnp$dyad <- paste(df_agg_motnp$id_1, df_agg_motnp$id_2, sep = '_')
df_agg_motnp$dyad_id_nogaps <- as.integer(as.factor(df_agg_motnp$dyad))

### load the edge weights
motnp <- readRDS('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_edgeweightestimates_mcmcoutput.rds') %>% 
  select(-`1.lp__`)
motnp <- motnp[, which(colnames(motnp) %in% df_agg_motnp$dyad)]
logit_edge_samples_motnp <- logit(motnp)

#### calculate posterior centralities ####
motnp <- as.matrix(logit_edge_samples_motnp)
edge_samples_motnp <- plogis(motnp)

# Build adjacency tensor
num_nodes <- length(unique(c(df_agg_motnp$id_1, df_agg_motnp$id_2)))
num_samples <- 4000
adj_tensor <- array(0, c(num_samples, num_nodes, num_nodes),
                    dimnames = list(NULL,
                                    unique(c(df_agg_motnp$id_1, df_agg_motnp$id_2)),
                                    unique(c(df_agg_motnp$id_1, df_agg_motnp$id_2))))
for (dyad_id in 1:nrow(df_agg_motnp)) {
  dyad_row <- df_agg_motnp[df_agg_motnp$dyad_id_nogaps == dyad_id, ]
  adj_tensor[, dyad_row$id_1, dyad_row$id_2] <- edge_samples_motnp[, dyad_id]
}
rm(dyad_row)

# Calculate centrality and store posterior samples in a matrix
eigen_samples_motnp <- matrix(0, num_samples, num_nodes)
eigen_samples_motnp_std <- matrix(0, num_samples, num_nodes)
for (i in 1:num_samples) {
  g <- graph_from_adjacency_matrix(adj_tensor[i, , ], mode="undirected", weighted=TRUE)
  eigen_samples_motnp[i, ] <- eigen_centrality(g)$vector
  eigen_samples_motnp_std[i, ] <- (eigen_samples_motnp[i, ] - mean(eigen_samples_motnp[i, ]))/sd(eigen_samples_motnp[i, ])
  rm(g)
}

# check
head(eigen_samples_motnp)     # Unstandardised eigenvector centrality
head(eigen_samples_motnp_std) # Standardised eigenvector centrality

# visualise eigenvector centralities
df_wide_motnp <- data.frame(eigen_samples_motnp_std)
colnames(df_wide_motnp) <- 1:num_nodes
df_long_motnp <- pivot_longer(df_wide_motnp, cols = 1:num_nodes, names_to="node_id", values_to="eigenvector_centrality")
ggplot(df_long_motnp, aes(x = eigenvector_centrality)) +
  geom_density(fill="#387780", alpha=0.7, size=0.4) +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

nodes1 <- df_agg_motnp %>% select(id_1, age_class_1, age_cat_id_1) %>% distinct()
nodes2 <- df_agg_motnp %>% select(id_2, age_class_2, age_cat_id_2) %>% distinct()
colnames(nodes1) <- c('id', 'age_class', 'age_cat_id')
colnames(nodes2) <- c('id', 'age_class', 'age_cat_id')
node_ages_motnp <- rbind(nodes1, nodes2) %>% distinct()
rm(nodes1, nodes2)

node_ages_motnp$node_id <- as.integer(as.factor(node_ages_motnp$id))
df_long_motnp$node_id <- as.integer(df_long_motnp$node_id)
df_long_motnp <- left_join(df_long_motnp, node_ages_motnp, by = 'node_id')

ggplot(df_long_motnp, aes(x = eigenvector_centrality)) +
  geom_density(aes(fill = as.factor(age_cat_id)), alpha=0.7, size=0.4) +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) # 2,3,6 = fairly well spread. 7 mostly low centrality. 4,5 = mostly high centrality.

df_long_motnp$age_cat <- ifelse(df_long_motnp$age_cat_id == 2, '9-10',
                                ifelse(df_long_motnp$age_cat_id == 3, '10-15',
                                       ifelse(df_long_motnp$age_cat_id == 4, '15-20',
                                              ifelse(df_long_motnp$age_cat_id == 5, '20-25',
                                                     ifelse(df_long_motnp$age_cat_id == 6, '25-40','>40')))))
ggplot(df_long_motnp, aes(x = eigenvector_centrality, group = id)) +
  geom_density(colour = rgb(0,0,1,0.2)) +
  facet_wrap(. ~ factor(age_cat, levels = c('9-10','10-15','15-20','20-25','25-40','>40'))) +
  labs(x = "eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text = element_text(angle = 0, size = 10, debug = FALSE),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = rgb(0,0,1,0.6)),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), 'cm'),
        legend.position = 'none')

#### compute normal approximation ####
# The posterior centralities can now be characterised by the multivariate normal distribution as an approximation to their true posterior distributions. To highlight the importance of sampling from a multivariate normal rather than a univariate normal, we can plot centrality samples from two nodes against each other to see how they co-vary.
plot(eigen_samples_motnp_std[, 1], eigen_samples_motnp_std[, 2], las = 1, pch = 19, col = rgb(0,0,1,0.2))

# This joint distribution contains important information that we would ideally like to include in the model. The multivariate normal will allow us to do this. Parameterising the multivariate normal approximation is a relatively simple process and can be done by calculating the sample mean and sample covariance matrix of the posterior centrality samples:
eigen_mu_motnp <- apply(eigen_samples_motnp_std, 2, mean) 
eigen_cov_motnp <- cov(eigen_samples_motnp_std)

# These quantities will be given to the Stan model as data to model joint posteriors of centrality in the regression. We can run a few quick plots to see how well the approximation is working and the covariance for one of the nodes:
eigen_samples_motnp_sim <- MASS::mvrnorm(1e5, eigen_mu_motnp, eigen_cov_motnp)
plot(density(eigen_samples_motnp_std[, 1]), lwd=2, main="Estimated standardised eigenvector centrality vs normal approximation", xlab="Logit edge weight")
lines(density(eigen_samples_motnp_sim[, 1]), col=rgb(0, 0, 1, 0.5), lwd=2)

#### eigenvector prior predictive checks (univariate Gaussian) ####
range(eigen_mu_motnp)
age <- 1:60

# eigenvector linear effect only
plot(NULL, xlim = c(0,60), ylim = c(-3,3), las = 1, xlab = 'age', ylab = 'mean eigenvector centrality (std)')
abline(h = range(eigen_mu_motnp), lty = 2)
for(i in 1:100){
  intercept <- rnorm(1, 0, 1)
  beta_age <- rnorm(1, 0, 0.02)
  lines(x = age, y = intercept + beta_age*age,
        col = rgb(0,0,1,0.4))
}

# eigenvector allowing for extremes to be most/least central
plot(NULL, xlim = c(0,60), ylim = c(-3,3), las = 1, xlab = 'age', ylab = 'mean eigenvector centrality (std)')
abline(h = range(eigen_mu_motnp), lty = 2)
for(i in 1:100){
  intercept <- rnorm(1, 0, 1)
  beta_age <- rnorm(1, 0, 0.02)
  beta_age2 <- rnorm(1, 0, 0.003)
  lines(x = age, y = intercept + beta_age*age + beta_age2*(age^2),
        col = rgb(0,0,1,0.4))
}

#### load MOTNP ages ####
true_ages <- readRDS('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_ageestimates_mcmcoutput.rds')

colnames(true_ages) <- motnp_males$id # check this shouldn't be as.integer(as.factor(motnp_males$id)), but 90% sure this is correct -- goes into the model in order so should come out in same order?
motnp_ap <- motnp_males[motnp_males$id %in% node_ages_motnp$id, ]
true_ages_ap <- true_ages[, colnames(true_ages) %in% node_ages_motnp$id]

modeldata_ages <- as.matrix(true_ages_ap) %>% t() # transpose to fit model
modeldata_ages_sq <- (modeldata_ages)^2

mean_age <- rep(NA, 218)
for(i in 1:length(mean_age)){
  elephant <- as.data.frame(true_ages_ap[,i])
  mean_age[i] <- mean(elephant[,1])
}
plot(eigen_mu_motnp ~ mean_age, las = 1, pch = 19, col = rgb(0,0,1,0.2),
     xlab = 'mean estimate for male age', ylab = 'mean estimate for male eigenvector centrality (standardised)')

#### eigenvector ~ age ####
model_data_motnp <- list(
  num_nodes = num_nodes,               # Number of dyads
  centrality_mu  = eigen_mu_motnp,     # Sample means of logit edge weights
  centrality_cov = eigen_cov_motnp,    # Sample covariance of logit edge weights
  node_age = modeldata_ages,           # Age of individual -- are these in the wrong order? should be M1, M10, M100, M101?
  node_age2 = modeldata_ages_sq        # Age squared for quadratic term
)
str(model_data_motnp)

model_eigen_cent <- stan_model("models/nodal_regression_eigenvector_agedistribution.stan")
fit_nodal_eigen <- sampling(model_eigen_cent, data = model_data_motnp, cores = 4, chains = 4)

### check traceplot
traceplot(fit_nodal_eigen, pars = c('beta_age','beta_age2','sigma','intercept'))

### calculate leave-one-out cross validation
loo(fit_nodal_eigen)

### posterior predictive checks
params <- rstan::extract(fit_nodal_eigen)
plot(density(eigen_samples_motnp_std[1, ]), main = "Posterior predictive density (standardised centrality)",
     col = rgb(0, 0, 0, 0.25), ylim = c(0, 1), las = 1)
for (i in 1:100) {
  j <- sample(1:num_samples, 1)                                         # select a network to plot
  lines(density(eigen_samples_motnp_std[j, ]), col=rgb(0, 0, 0, 0.25)) # plot centrality density for network j (black)
  mu <- params$beta_age[j]*apply(model_data_motnp$node_age,1,mean)     # extract age slope parameter for network j
  sigma <- eigen_cov_motnp + diag(rep(params$sigma[j], num_nodes))     # extract variance parameter for network j
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))  # plot predicted centralities for network j (blue)
} # POOR FIT -- predicted are all falling as a bell curve, observed centrality density as multimodal

### plot age against predictions
node_age <- 1:60
beta_age_mean  <- summary(fit_nodal_eigen)$summary[1,1]
beta_age_2.5   <- summary(fit_nodal_eigen)$summary[1,4]
beta_age_97.5  <- summary(fit_nodal_eigen)$summary[1,8]
beta_age2_mean <- summary(fit_nodal_eigen)$summary[2,1]
beta_age2_2.5  <- summary(fit_nodal_eigen)$summary[2,4]
beta_age2_97.5 <- summary(fit_nodal_eigen)$summary[2,8]
intercept_mean <- summary(fit_nodal_eigen)$summary[3,1]
intercept_2.5  <- summary(fit_nodal_eigen)$summary[3,4]
intercept_97.5 <- summary(fit_nodal_eigen)$summary[3,8]
mu <- intercept_mean  + beta_age_mean*node_age + beta_age2_mean*(node_age^2)
min <- intercept_2.5  + beta_age_2.5*node_age  + beta_age2_2.5*(node_age^2)
max <- intercept_97.5 + beta_age_97.5*node_age + beta_age2_97.5*(node_age^2)
summary(mu)
summary(min)
summary(max)
plot(NULL, xlim = c(0,60), ylim = c(-4, 4), lwd = 2, las = 1,
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
    points(eigen_mu_motnp[i] ~ model_data_motnp$node_age[i,s[j]], pch = 19, col = rgb(0,0,0,0.05))
  }
}

# Calculate the 95% credible intervals of the model parameters
round(summary(fit_nodal_eigen)$summary[1:4, c(1:4, 8)], 5) # basically no effect whatsoever (remembering this is in standard deviations on the outcome scale but years of age on the exponent)

# clear environment
rm(list = ls())

############ Betweenness ############
#### load MOTNP nodes, edges and interactions data ####
motnp_males <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_elenodes.csv') %>% 
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

### import data for aggregated model (binomial)
df_agg_motnp <- read_delim('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_bayesian_binomialpairwiseevents.csv', delim = ',') %>% 
  filter(dem_class_1 == 'AM' | dem_class_1 == 'PM') %>% 
  filter(dem_class_2 == 'AM' | dem_class_2 == 'PM')
df_agg_motnp$sex_1 <- 'M'
df_agg_motnp$age_cat_id_1 <- ifelse(df_agg_motnp$age_category_1 == '9-10', 2,
                                    ifelse(df_agg_motnp$age_category_1 == '10-15', 3,
                                           ifelse(df_agg_motnp$age_category_1 == '15-19', 4,
                                                  ifelse(df_agg_motnp$age_category_1 == '20-25', 5,
                                                         ifelse(df_agg_motnp$age_category_1 == '25-40', 6,
                                                                ifelse(df_agg_motnp$age_category_1 == '40+', 7, 1))))))
df_agg_motnp$age_class_1 <- ifelse(df_agg_motnp$age_cat_id_1 == 2, 'Juvenile',
                                   ifelse(df_agg_motnp$age_cat_id_1 > 4, 'Adult','Pubescent'))
df_agg_motnp$age_cat_id_2 <- ifelse(df_agg_motnp$age_category_2 == '9-10', 2,
                                    ifelse(df_agg_motnp$age_category_2 == '10-15', 3,
                                           ifelse(df_agg_motnp$age_category_2 == '15-19', 4,
                                                  ifelse(df_agg_motnp$age_category_2 == '20-25', 5,
                                                         ifelse(df_agg_motnp$age_category_2 == '25-40', 6,
                                                                ifelse(df_agg_motnp$age_category_2 == '40+', 7, 1))))))
df_agg_motnp$age_class_2 <- ifelse(df_agg_motnp$age_cat_id_2 == 2, 'Juvenile',
                                   ifelse(df_agg_motnp$age_cat_id_2 > 4, 'Adult','Pubescent'))
table(df_agg_motnp$age_cat_id_1,df_agg_motnp$age_category_1)
table(df_agg_motnp$age_cat_id_2,df_agg_motnp$age_category_2)

df_agg_motnp$dem_class_1 <- ifelse(df_agg_motnp$age_class_1 == 'Adult', 'AM',
                                   ifelse(df_agg_motnp$age_class_1 == 'Pubescent', 'PM', 'JM'))
df_agg_motnp$dem_class_2 <- ifelse(df_agg_motnp$age_class_2 == 'Adult', 'AM',
                                   ifelse(df_agg_motnp$age_class_2 == 'Pubescent', 'PM', 'JM'))
df_agg_motnp$dem_type <- ifelse(df_agg_motnp$age_cat_id_1 >= df_agg_motnp$age_cat_id_2,
                                paste(df_agg_motnp$dem_class_1, df_agg_motnp$dem_class_2, sep = '_'),
                                paste(df_agg_motnp$dem_class_2, df_agg_motnp$dem_class_1, sep = '_'))
df_agg_motnp$count_dyad <- (df_agg_motnp$count_1 + df_agg_motnp$count_2) - df_agg_motnp$event_count  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.
df_agg_motnp$node_1_nogaps <- as.integer(as.factor(df_agg_motnp$node_1))
df_agg_motnp$node_2_nogaps <- as.integer(as.factor(df_agg_motnp$node_2))+1
df_agg_motnp$dyad <- paste(df_agg_motnp$id_1, df_agg_motnp$id_2, sep = '_')
df_agg_motnp$dyad_id_nogaps <- as.integer(as.factor(df_agg_motnp$dyad))

### load the edge weights
motnp <- readRDS('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_edgeweightestimates_mcmcoutput.rds') %>% 
  select(-`1.lp__`)
motnp <- motnp[, which(colnames(motnp) %in% df_agg_motnp$dyad)]
logit_edge_samples_motnp <- logit(motnp)

#### trim down data for debugging and running quickly ####
sample_males <- motnp_males[motnp_males$dem_class == 'AM' | motnp_males$dem_class == 'PM',]
set.seed(1) ; sample_males <- sample(sample_males$id, 30, FALSE)
test_males <- motnp_males[motnp_males$id %in% sample_males,]
df_agg_test <- df_agg_motnp[df_agg_motnp$id_1 %in% sample_males &
                              df_agg_motnp$id_2 %in% sample_males,]
df_agg_test$dyad_id_nogaps <- as.integer(as.factor(df_agg_test$dyad))
test <- motnp[,colnames(motnp) %in% df_agg_test$dyad]
logit_edge_samples_test <- logit_edge_samples_motnp[,colnames(logit_edge_samples_motnp) %in% df_agg_test$dyad]

#### calculate posterior centralities ####
test <- as.matrix(logit_edge_samples_test)
edge_samples_test <- plogis(test)

# Build adjacency tensor
num_nodes <- length(unique(c(df_agg_test$id_1, df_agg_test$id_2)))
num_samples <- 4000
adj_tensor <- array(0, c(num_samples, num_nodes, num_nodes),
                    dimnames = list(NULL,
                                    unique(c(df_agg_test$id_1, df_agg_test$id_2)),
                                    unique(c(df_agg_test$id_1, df_agg_test$id_2))))
for (dyad_id in 1:nrow(df_agg_test)) {
  dyad_row <- df_agg_test[df_agg_test$dyad_id_nogaps == dyad_id, ]
  adj_tensor[, dyad_row$id_1, dyad_row$id_2] <- edge_samples_test[, dyad_id]
}
rm(dyad_row)

# Calculate centrality and store posterior samples in a matrix
btwn_samples_test <- matrix(0, num_samples, num_nodes)
btwn_samples_test_std <- matrix(0, num_samples, num_nodes)
for (i in 1:num_samples) {
  g <- graph_from_adjacency_matrix(adj_tensor[i, , ], mode="undirected", weighted=TRUE)
  btwn_samples_test[i, ] <- 1 - (betweenness(g, normalized = TRUE))
  btwn_samples_test_std[i, ] <- (btwn_samples_test[i, ] - mean(btwn_samples_test[i, ]))/sd(btwn_samples_test[i, ])
  rm(g)
}

# convert betweenness to a Z-score
mu_btwn <- mean(btwn_samples_test)
sd_btwn <- sd(btwn_samples_test)

# check
hist(btwn_samples_test)      # Unstandardised betweenness centrality -- 0-1 bounded
hist(btwn_samples_test_std)  # betweenness centrality Z scores -- unbounded

# visualise betweenness centrality
df_wide_test <- data.frame(btwn_samples_test)
colnames(df_wide_test) <- 1:num_nodes
df_long_test <- pivot_longer(df_wide_test, cols = 1:num_nodes, names_to="node_id", values_to="btwn_centrality")
ggplot(df_long_test, aes(x = btwn_centrality)) +
  geom_density(fill="#387780", alpha=0.7, size=0.4) +
  labs(x="betweenness centrality") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

nodes1 <- df_agg_test %>% select(id_1, age_class_1, age_cat_id_1) %>% distinct()
nodes2 <- df_agg_test %>% select(id_2, age_class_2, age_cat_id_2) %>% distinct()
colnames(nodes1) <- c('id', 'age_class', 'age_cat_id')
colnames(nodes2) <- c('id', 'age_class', 'age_cat_id')
node_ages_test <- rbind(nodes1, nodes2) %>% distinct()
rm(nodes1, nodes2)

node_ages_test$node_id <- as.integer(as.factor(node_ages_test$id))
df_long_test$node_id <- as.integer(df_long_test$node_id)
df_long_test <- left_join(df_long_test, node_ages_test, by = 'node_id')

ggplot(df_long_test, aes(x = btwn_centrality, fill = as.factor(age_cat_id))) +
  geom_density(alpha=0.7, size=0.4) +
  labs(x="betweenness centrality") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

df_long_test$age_cat <- ifelse(df_long_test$age_cat_id == 2, '9-10',
                               ifelse(df_long_test$age_cat_id == 3, '10-15',
                                      ifelse(df_long_test$age_cat_id == 4, '15-20',
                                             ifelse(df_long_test$age_cat_id == 5, '20-25',
                                                    ifelse(df_long_test$age_cat_id == 6, '25-40','>40')))))
ggplot(df_long_test, aes(x = btwn_centrality, group = id)) +
  geom_density(colour = rgb(0,0,1,0.2)) +
  facet_wrap(. ~ factor(age_cat, levels = c('9-10','10-15','15-20','20-25','25-40','>40'))) +
  labs(x = "betweenness centrality") + 
  theme_light() + 
  theme(axis.text = element_text(angle = 0, size = 10, debug = FALSE),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = rgb(0,0,1,0.6)),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), 'cm'),
        legend.position = 'none')

#### compute normal approximation ####
# The posterior centralities can now be characterised by the multivariate normal distribution as an approximation to their true posterior distributions. To highlight the importance of sampling from a multivariate normal rather than a univariate normal, we can plot centrality samples from two nodes against each other to see how they co-vary.
plot(btwn_samples_test[, 1], btwn_samples_test[, 2], las = 1, pch = 19, col = rgb(0,0,1,0.2))
plot(btwn_samples_test[, 1], btwn_samples_test[, 3], las = 1, pch = 19, col = rgb(0,0,1,0.2))
plot(btwn_samples_test[, 3], btwn_samples_test[, 2], las = 1, pch = 19, col = rgb(0,0,1,0.2))
plot(btwn_samples_test[, 4], btwn_samples_test[, 3], las = 1, pch = 19, col = rgb(0,0,1,0.2))

plot(NULL, xlim = c(0,1), ylim = c(0,300), las = 1, xlab = 'betweenness', ylab = 'density')
for(i in 1:num_nodes){
  lines(density(btwn_samples_test[,i]), col = rgb(0,0,1,0.1))
}
plot(NULL, xlim = c(-6,2), ylim = c(0,3), las = 1, xlab = 'betweenness z-score', ylab = 'density')
for(i in 1:num_nodes){
  lines(density(btwn_samples_test_std[,i]), col = rgb(0,0,1,0.1))
}

# This joint distribution contains important information that we would ideally like to include in the model. The multivariate normal will allow us to do this. Parameterising the multivariate normal approximation is a relatively simple process and can be done by calculating the sample mean and sample covariance matrix of the posterior centrality samples:
btwn_mu_test <- apply(btwn_samples_test, 2, mean) 
btwn_cov_test <- cov(btwn_samples_test)

# These quantities will be given to the Stan model as data to model joint posteriors of centrality in the regression. We can run a few quick plots to see how well the approximation is working and the covariance for one of the nodes:
# betweenness (remembering these are based on z-scores of normalised measures)
btwn_samples_test_sim <- MASS::mvrnorm(1e5, btwn_mu_test, btwn_cov_test)
plot(density(btwn_samples_test_std[, 1]), lwd = 2, xlim = c(-2,6), ylim = c(0,5),
     main="Estimated standardised betweenness centrality vs normal approximation", xlab = "Logit edge weight")
lines(density(btwn_samples_test_sim[, 1]), col=rgb(0, 0, 1, 0.5), lwd=2)
for( i in 2:ncol(btwn_samples_test_sim)){
  lines(density(btwn_samples_test_std[, i]), col = rgb(0,0,0,0.2), lwd = 0.5)
  lines(density(btwn_samples_test_sim[, i]), col = rgb(0,0,1,0.2), lwd = 0.5)
}     # this DOES looks too bad! the black lines are the measured values, the blue lines are the simulated versions

# try again for betweenness using Dirichlet distribution instead of multi-normal -- HOW DO I CONVERT btwn_mu_test/btwn_cov_test TO SOMETHING THAT CAN FORM A VECTOR OF LENGTH 3?
plot(density(btwn_samples_test[, 1]), las = 1, main="Estimated standardised betweenness centrality vs normal approximation", xlab="Logit edge weight", xlim = c(0,1))
btwn_samples_test_sim <- rdirichlet(length(btwn_mu_test)*1000,
                                                   alpha = c(btwn_mu_test[1],btwn_cov_test[1,1],
                                                             btwn_cov_test[1,1]))
lines(density(btwn_samples_test_sim), col=rgb(0, 0, 1, 0.5), lwd=2) # this looks terrible, but might just be because I don't know how rdirichlet works

btwn_samples_test_sim <- array(NA, dim = c(num_nodes, num_nodes, 1000))
for( i in 1:num_nodes){
  for( j in 1:num_nodes){
  btwn_samples_test_sim[i,j,] <- rdirichlet(1000,alpha = c(btwn_mu_test[i],
                                                           btwn_cov_test[i,j],btwn_cov_test[i,j])) # doesn't work because some are -ve
  }
}

plot(density(btwn_samples_test[, 1]), las = 1, main="Estimated standardised betweenness centrality vs normal approximation", xlab="Logit edge weight", xlim = c(0,1))
btwn_samples_test_sim <- 1 - rdirichlet(num_nodes, alpha = btwn_mu_test) # only works when using unstandardised values
lines(density(btwn_samples_test_sim), col=rgb(0, 0, 1, 0.5), lwd=2) # this looks much better!

#### betweenness prior predictive checks (univariate Gaussian) ####
range(btwn_mu_test)
age <- 1:60

# betweenness linear effect only -- not normalised
range(btwn_mu_test)
plot(NULL, xlim = c(0,60), ylim = c(-1,6), las = 1, xlab = 'age', ylab = 'mean betweenness centrality (std)')
abline(h = range(btwn_mu_test), lty = 2)
for(i in 1:100){
  intercept <- rnorm(1, 2, 1.4)
  beta_age <- rnorm(1, 0, 0.02)
  lines(x = age, y = intercept + beta_age*age,
        col = rgb(0,0,1,0.4))
}

# betweenness allowing for extremes to be most/least central -- not normalised
plot(NULL, xlim = c(0,60),  ylim = c(-1,6), las = 1, xlab = 'age', ylab = 'mean betweenness centrality (std)')
abline(h = range(btwn_mu_test), lty = 2)
for(i in 1:100){
  intercept <- rnorm(1, 2, 1.4)
  beta_age <- rnorm(1, 0, 0.02)
  beta_age2 <- rnorm(1, 0, 0.001)
  lines(x = age, y = intercept + beta_age*age + beta_age2*(age^2),
        col = rgb(0,0,1,0.4))
}

# betweenness linear effect only -- normalised
plot(NULL, xlim = c(0,60), ylim = c(-0.2,0.2), las = 1,
     xlab = 'age', ylab = 'mean betweenness centrality (std)')
abline(h = range(btwn_mu_test), lty = 2)
for(i in 1:100){
  intercept <- rnorm(1, 0.05, 0.02)
  beta_age <- rnorm(1, 0, 0.002)
  lines(x = age, y = intercept + beta_age*age,
        col = rgb(0,0,1,0.4))
}

# betweenness allowing for extremes to be most/least central -- normalised
plot(NULL, xlim = c(0,60), ylim = c(-0.3,0.3), las = 1,
     xlab = 'age', ylab = 'mean betweenness centrality (std)')
abline(h = range(btwn_mu_test), lty = 2)
for(i in 1:100){
  intercept <- rnorm(1, 0.05, 0.02)
  beta_age <- rnorm(1, 0, 0.002)
  beta_age2 <- rnorm(1, 0, 0.0001)
  lines(x = age, y = intercept + beta_age*age + beta_age2*(age^2),
        col = rgb(0,0,1,0.4))
}

#### load test ages ####
true_ages <- readRDS('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_ageestimates_mcmcoutput.rds')

true_ages_ap <- true_ages[, colnames(true_ages) %in% node_ages_test$id]

modeldata_ages <- as.matrix(true_ages_ap) %>% t()
modeldata_ages_sq <- (modeldata_ages)^2

#### betweenness ~ age -- SIMULATE FIRST TO TEST: USE A DIRICHLET DISTRIBUTION (MULTIVARIATE BETA) WITH NORMALIZED BETWEENNESS SCORES ####
# work out what a dirichlet looks like -- order of alpha values doesn't matter, values do -- NUMBER OF VALUES DOES NOT HAVE TO BE 3!!
rethinking::dens(rdirichlet(10000, alpha = c(1,1,1)))       # right skewed, fairly triangular
rethinking::dens(rdirichlet(10000, alpha = c(5,1,1)))       # gains weight to the right
rethinking::dens(rdirichlet(10000, alpha = c(0.2,1,1)))     # narrow left peak
rethinking::dens(rdirichlet(10000, alpha = c(1,5,1)))       # same as 5,1,1
rethinking::dens(rdirichlet(10000, alpha = c(1,0.2,1)))     # same as 0.2,1,1
rethinking::dens(rdirichlet(10000, alpha = c(1,1,5)))       # same as 5,1,1
rethinking::dens(rdirichlet(10000, alpha = c(1,1,0.2)))     # same as 0.2,1,1
rethinking::dens(rdirichlet(10000, alpha = c(5,0.2,1)))     # narrow left peak with weight on right
rethinking::dens(rdirichlet(10000, alpha = c(0.2,5,1)))     # same as 5,0.2,1
rethinking::dens(rdirichlet(10000, alpha = c(1,0.2,5)))     # same as 5,0.2,1
rethinking::dens(rdirichlet(10000, alpha = c(2,2,2)))       # right skewed normal
rethinking::dens(rdirichlet(10000, alpha = c(10,10,10)))    # almost normal, narrower
rethinking::dens(rdirichlet(10000, alpha = c(0.2,0.2,0.2))) # double peak, almost at 0 and 1
rethinking::dens(rdirichlet(10000, alpha = c(5,0,0)))       # same no matter the number
rethinking::dens(rdirichlet(10000, alpha = c(0.01,1,1)))    # most density very close to zero
rethinking::dens(rdirichlet(10000, alpha = c(1,0.01,1)))    # same as 0.01,1,1
rethinking::dens(rdirichlet(10000, alpha = c(1,1,0.01)))    # same as 0.01,1,1
rethinking::dens(rdirichlet(10000, alpha = c(0.01,1,1)))    # same as 0.01,1,1
rethinking::dens(rdirichlet(10000, alpha = c(-2,-1,1)))     # can't include negatives

plot(density(btwn_mu_test))

# simulate from population similar to test
sim <- 1 - rdirichlet(num_nodes, alpha = c(0.01,1,1))
plot(density(sim))

#### run on test data ####
model_data_test <- list(
  num_nodes = 30,               # Number of dyads
  centrality_mu  = btwn_mu_test,      # Sample means of logit edge weights
  centrality_cov = btwn_cov_test,     # Sample covariance of logit edge weights
  node_age = modeldata_ages,           # Age of individual -- are these in the wrong order? should be M1, M10, M100, M101?
  node_age2 = modeldata_ages_sq        # Age squared for quadratic term
)
str(model_data_test)

model_btwn_cent <- stan_model("models/nodal_regression_betweenness_agedistribution.stan")
fit_nodal_btwn <- sampling(model_btwn_cent, data = model_data_test, cores = 1, chains = 1)

### check traceplot
traceplot(fit_nodal_btwn)

### calculate leave-one-out cross validation
fit_nodal_btwn

hist(summary(fit_nodal_btwn)$summary[,10], breaks = 20, xlim = c(0.99,1.02), xlab = 'rhat', main = 'rhat')

### posterior predictive checks
params <- rstan::extract(fit_nodal_btwn)
plot(density(btwn_samples_test_std[1, ]), main = "Posterior predictive density (standardised centrality)",
     col = rgb(0, 0, 0, 0.25), ylim = c(0, 2), las = 1)
for (i in 1:100) {
  j <- sample(1:num_samples, 1)                                         # select a network to plot
  lines(density(btwn_samples_test_std[j, ]), col=rgb(0, 0, 0, 0.25))   # plot centrality density for network j (black)
  mu <- params$beta_age[j]*apply(model_data_test$node_age,1,mean)      # extract age slope parameter for network j
  sigma <- btwn_cov_test + diag(rep(params$sigma[j], num_nodes))       # extract variance parameter for network j
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))   # plot predicted centralities for network j (blue)
} # POOR FIT -- predicted are all falling as a bell curve, observed centrality density as multimodal
abline(v = c(0, -0.5), lty = 2, lwd = 2, col = 'red')

### plot age against predictions
node_age <- 1:60
beta_age_mean  <- summary(fit_nodal_btwn)$summary[1,1]
beta_age_2.5   <- summary(fit_nodal_btwn)$summary[1,4]
beta_age_97.5  <- summary(fit_nodal_btwn)$summary[1,8]
#beta_age2_mean <- summary(fit_nodal_btwn)$summary[2,1]
#beta_age2_2.5  <- summary(fit_nodal_btwn)$summary[2,4]
#beta_age2_97.5 <- summary(fit_nodal_btwn)$summary[2,8]
intercept_mean <- summary(fit_nodal_btwn)$summary[3,1]
intercept_2.5  <- summary(fit_nodal_btwn)$summary[3,4]
intercept_97.5 <- summary(fit_nodal_btwn)$summary[3,8]
#intercept_mean <- summary(fit_nodal_btwn)$summary[2,1]
#intercept_2.5  <- summary(fit_nodal_btwn)$summary[2,4]
#intercept_97.5 <- summary(fit_nodal_btwn)$summary[2,8]
mu <- intercept_mean  + beta_age_mean*node_age# + beta_age2_mean*(node_age^2)
min <- intercept_2.5  + beta_age_2.5*node_age#  + beta_age2_2.5*(node_age^2)
max <- intercept_97.5 + beta_age_97.5*node_age# + beta_age2_97.5*(node_age^2)
summary(mu)
summary(min)
summary(max)
plot(NULL, xlim = c(0,60), ylim = c(-1,5), lwd = 2, las = 1,
     xlab = 'age (years)', ylab = 'betweenness centrality (std)')
for(i in 1:100){
  j <- sample(1:num_samples, 1)
  int <- params$intercept[j]
  bA1 <- params$beta_age[j]
  #bA2 <- params$beta_age2[j]
  mean <- int + bA1*age# + bA2*age
  lines(c(1:60), mean, col = rgb(0.5,0,1,0.2))
}
lines(mu ~ node_age, lwd = 2) 
lines(min ~ node_age, col = 'red')
lines(max ~ node_age, col = 'red')
abline(h = c(-2,2), lty = 2)
for(i in 1:num_nodes){
  s <- sample(x = 1:num_samples, size = 100, replace = F)
  for(j in 1:length(s)){
    points(btwn_samples_test_std[i,j] ~ model_data_test$node_age[i,s[j]], pch = 19, cex = 0.5, col = rgb(0,0,1,0.05))
 }
}

# Calculate the 95% credible intervals of the model parameters
round(summary(fit_nodal_btwn)$summary[1:4, c(1:4, 8)], 5) # basically no effect whatsoever (remembering this is in standard deviations on the outcome scale but years of age on the exponent)

#### run on MOTNP data ####
#model_data_motnp <- list(
#  num_nodes = num_nodes,               # Number of dyads
#  centrality_mu  = btwn_mu_motnp,      # Sample means of logit edge weights
#  centrality_cov = btwn_cov_motnp,     # Sample covariance of logit edge weights
#  node_age = modeldata_ages#,           # Age of individual -- are these in the wrong order? should be M1, M10, M100, M101?
#node_age2 = modeldata_ages_sq        # Age squared for quadratic term
#)
#str(model_data_motnp)

#model_btwn_cent <- stan_model("models/nodal_regression_betweenness_agedistribution.stan")
#fit_nodal_btwn <- sampling(model_btwn_cent, data = model_data_motnp, cores = 1, chains = 1)

### check traceplot
#traceplot(fit_nodal_btwn)

### calculate leave-one-out cross validation
#fit_nodal_btwn

#hist(summary(fit_nodal_btwn)$summary[,10], breaks = 20, xlim = c(0.99,1.02), xlab = 'rhat', main = 'rhat')

### posterior predictive checks
#params <- rstan::extract(fit_nodal_btwn)
#plot(density(btwn_samples_motnp_std[1, ]), main = "Posterior predictive density (standardised centrality)",
#     col = rgb(0, 0, 0, 0.25), ylim = c(0, 2), las = 1)
#for (i in 1:100) {
#  j <- sample(1:num_samples, 1)                                         # select a network to plot
#  lines(density(btwn_samples_motnp_std[j, ]), col=rgb(0, 0, 0, 0.25))   # plot centrality density for network j (black)
#  mu <- params$beta_age[j]*apply(model_data_motnp$node_age,1,mean)      # extract age slope parameter for network j
#  sigma <- btwn_cov_motnp + diag(rep(params$sigma[j], num_nodes))       # extract variance parameter for network j
#  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))   # plot predicted centralities for network j (blue)
#} # POOR FIT -- predicted are all falling as a bell curve, observed centrality density as multimodal
#abline(v = c(0, -0.5), lty = 2, lwd = 2, col = 'red')

### plot age against predictions
#node_age <- 1:60
#beta_age_mean  <- summary(fit_nodal_btwn)$summary[1,1]
#beta_age_2.5   <- summary(fit_nodal_btwn)$summary[1,4]
#beta_age_97.5  <- summary(fit_nodal_btwn)$summary[1,8]
##beta_age2_mean <- summary(fit_nodal_btwn)$summary[2,1]
##beta_age2_2.5  <- summary(fit_nodal_btwn)$summary[2,4]
##beta_age2_97.5 <- summary(fit_nodal_btwn)$summary[2,8]
#intercept_mean <- summary(fit_nodal_btwn)$summary[3,1]
#intercept_2.5  <- summary(fit_nodal_btwn)$summary[3,4]
#intercept_97.5 <- summary(fit_nodal_btwn)$summary[3,8]
##intercept_mean <- summary(fit_nodal_btwn)$summary[2,1]
##intercept_2.5  <- summary(fit_nodal_btwn)$summary[2,4]
##intercept_97.5 <- summary(fit_nodal_btwn)$summary[2,8]
#mu <- intercept_mean  + beta_age_mean*node_age# + beta_age2_mean*(node_age^2)
#min <- intercept_2.5  + beta_age_2.5*node_age#  + beta_age2_2.5*(node_age^2)
#max <- intercept_97.5 + beta_age_97.5*node_age# + beta_age2_97.5*(node_age^2)
#summary(mu)
#summary(min)
#summary(max)
#plot(NULL, xlim = c(0,60), ylim = c(-1,5), lwd = 2, las = 1,
#     xlab = 'age (years)', ylab = 'betweenness centrality (std)')
#for(i in 1:100){
#  j <- sample(1:num_samples, 1)
#  int <- params$intercept[j]
#  bA1 <- params$beta_age[j]
#  #bA2 <- params$beta_age2[j]
#  mean <- int + bA1*age# + bA2*age
#  lines(c(1:60), mean, col = rgb(0.5,0,1,0.2))
#}
#lines(mu ~ node_age, lwd = 2) 
#lines(min ~ node_age, col = 'red')
#lines(max ~ node_age, col = 'red')
#abline(h = c(-2,2), lty = 2)
#for(i in 1:num_nodes){
#  s <- sample(x = 1:num_samples, size = 100, replace = F)
#  for(j in 1:length(s)){
#    points(btwn_samples_motnp_std[i,j] ~ model_data_motnp$node_age[i,s[j]], pch = 19, cex = 0.5, col = rgb(0,0,1,0.05))
# }
#}

# Calculate the 95% credible intervals of the model parameters
#round(summary(fit_nodal_btwn)$summary[1:4, c(1:4, 8)], 5) # basically no effect whatsoever (remembering this is in standard deviations on the outcome scale but years of age on the exponent)
