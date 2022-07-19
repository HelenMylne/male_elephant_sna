#### information #####
# 1) to estimate the ages of individuals in other populations, first came up with a survival curve for the Amboseli elephants which have good birth and death data.
# 2) selected the Gompertz bathtub distribution
# 3) use distribution to convert MOTNP ages in categories to probability distributions of true age (first part of this script)
# 4) extract node centrality distribution from draws of edge weight
# 5) use MOTNP age distributions to predict network centrality (third part of this script)

############ 1) Obtain age estimates ############
#### load packages ####
library(tidyverse)
library(cmdstanr)
library(ggdist)
library(posterior)
library(bayesplot)
library(rstan)
library(igraph)
library(LaplacesDemon)

#### load model ####
# read in Stan model to estimate ages based on Gompertz bathtub distribution from ANP
latent_age_ordinal_model <- cmdstan_model("models/age_estimation/motnp_elephant_latent_age_ordinal_regression_hkm_22.07.07.stan")

#### load MOTNP data ####
motnp_males <- read_csv('data_processed/motnp_elenodes_22.01.13.csv') %>% 
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

#### create data list
N_motnp <- nrow(motnp_males)
K <- 8
motnp_ls <- list(
  N = N_motnp,
  K = K,
  age_category_index = motnp_males$age_cat_id)
hist(motnp_ls$age_category_index)

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

############ 2) Extract node centralities ############
#### read in MOTNP data ########
### import data for aggregated model (binomial)
df_agg_motnp <- read_delim('data_processed/motnp_bayesian_allpairwiseevents_splitbygrouptype_22.01.13.csv', delim = ',') %>% 
  filter(dem_class_1 == 'AM' | dem_class_1 == 'PM') %>% 
  filter(dem_class_2 == 'AM' | dem_class_2 == 'PM')
df_agg_motnp$sex_1 <- 'M'
df_agg_motnp$age_cat_id_1 <- ifelse(df_agg_motnp$age_category_1 == '9-10', 2,
                                    ifelse(df_agg_motnp$age_category_1 == '10-15', 3,
                                           ifelse(df_agg_motnp$age_category_1 == '15-19', 4,
                                                  ifelse(df_agg_motnp$age_category_1 == '20-25', 5,
                                                         ifelse(df_agg_motnp$age_category_1 == '25-40', 6,
                                                                ifelse(df_agg_motnp$age_category_1 == '40+', 7,
                                                                       df_agg_motnp$age_category_1))))))
df_agg_motnp$age_class_1 <- ifelse(df_agg_motnp$age_cat_id_1 == 2, 'Juvenile',
                                   ifelse(df_agg_motnp$age_cat_id_1 > 4, 'Adult','Pubescent'))
df_agg_motnp$age_cat_id_2 <- ifelse(df_agg_motnp$age_category_2 == '9-10', 2,
                                    ifelse(df_agg_motnp$age_category_2 == '10-15', 3,
                                           ifelse(df_agg_motnp$age_category_2 == '15-19', 4,
                                                  ifelse(df_agg_motnp$age_category_2 == '20-25', 5,
                                                         ifelse(df_agg_motnp$age_category_2 == '25-40', 6,
                                                                ifelse(df_agg_motnp$age_category_2 == '40+', 7,
                                                                       df_agg_motnp$age_category_2))))))
df_agg_motnp$age_class_2 <- ifelse(df_agg_motnp$age_cat_id_2 == 2, 'Juvenile',
                                   ifelse(df_agg_motnp$age_cat_id_2 > 4, 'Adult','Pubescent'))
df_agg_motnp$dem_class_1 <- ifelse(df_agg_motnp$age_class_1 == 'Adult', 'AM',
                                   ifelse(df_agg_motnp$age_class_1 == 'Pubescent', 'PM', 'JM'))
df_agg_motnp$dem_class_2 <- ifelse(df_agg_motnp$age_class_2 == 'Adult', 'AM',
                                   ifelse(df_agg_motnp$age_class_2 == 'Pubescent', 'PM', 'JM'))
df_agg_motnp$dem_type <- ifelse(df_agg_motnp$age_cat_id_1 >= df_agg_motnp$age_cat_id_2,
                                paste(df_agg_motnp$dem_class_1, df_agg_motnp$dem_class_2, sep = '_'),
                                paste(df_agg_motnp$dem_class_2, df_agg_motnp$dem_class_1, sep = '_'))
df_agg_motnp$count_dyad <- (df_agg_motnp$count_1 + df_agg_motnp$count_2) - df_agg_motnp$all_events  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.
df_agg_motnp$node_1_nogaps <- as.integer(as.factor(df_agg_motnp$node_1))
df_agg_motnp$node_2_nogaps <- as.integer(as.factor(df_agg_motnp$node_2))+1
df_agg_motnp$dyad_id_nogaps <- as.integer(as.factor(df_agg_motnp$dyad))

### load the edge weights
motnp <- readRDS('data_processed/motnp_bayesian_edgedistributions_a2.b2_22.02.07.rds') %>% 
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
btwn_samples_motnp <- matrix(0, num_samples, num_nodes)
btwn_samples_motnp_std <- matrix(0, num_samples, num_nodes)
for (i in 1:num_samples) {
  g <- graph_from_adjacency_matrix(adj_tensor[i, , ], mode="undirected", weighted=TRUE)
  eigen_samples_motnp[i, ] <- eigen_centrality(g)$vector
  eigen_samples_motnp_std[i, ] <- (eigen_samples_motnp[i, ] - mean(eigen_samples_motnp[i, ]))/sd(eigen_samples_motnp[i, ])
  btwn_samples_motnp[i, ] <- betweenness(g)
  btwn_samples_motnp_std[i, ] <- (btwn_samples_motnp[i, ] - mean(btwn_samples_motnp[i, ]))/sd(btwn_samples_motnp[i, ])
  rm(g)
}
head(eigen_samples_motnp)     # Unstandardised eigenvector centrality
head(eigen_samples_motnp_std) # Standardised eigenvector centrality
head(btwn_samples_motnp)      # Unstandardised betweenness centrality
head(btwn_samples_motnp_std)  # Standardised betweenness centrality

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
  geom_density(aes(fill = age_cat_id), alpha=0.7, size=0.4) +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) # 2,3,6 = fairly well spread. 7 mostly low centrality. 4,5 = mostly high centrality.

# visualise betweenness centrality
df_wide_motnp <- data.frame(btwn_samples_motnp_std)
colnames(df_wide_motnp) <- 1:num_nodes
df_long_motnp <- pivot_longer(df_wide_motnp, cols = 1:num_nodes, names_to="node_id", values_to="btwn_centrality")
ggplot(df_long_motnp, aes(x = btwn_centrality)) +
  geom_density(fill="#387780", alpha=0.7, size=0.4) +
  labs(x="betweenness centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

node_ages_motnp$node_id <- as.integer(as.factor(node_ages_motnp$id))
df_long_motnp$node_id <- as.integer(df_long_motnp$node_id)
df_long_motnp <- left_join(df_long_motnp, node_ages_motnp, by = 'node_id')

ggplot(df_long_motnp, aes(x = btwn_centrality)) +
  geom_density(aes(fill = age_cat_id), alpha=0.7, size=0.4) +
  labs(x="betweenness centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) # well that would definitely suggest that there is no effect...

############ 3) Run nodal regression model ############
#### compute normal approximation ####
# The posterior centralities can now be characterised by the multivariate normal distribution as an approximation to their true posterior distributions. To highlight the importance of sampling from a multivariate normal rather than a univariate normal, we can plot centrality samples from two nodes against each other to see how they co-vary.
plot(eigen_samples_motnp_std[, 1], eigen_samples_motnp_std[, 2])
plot(btwn_samples_motnp_std[, 1], btwn_samples_motnp_std[, 2])

# This joint distribution contains important information that we would ideally like to include in the model. The multivariate normal will allow us to do this. Parameterising the multivariate normal approximation is a relatively simple process and can be done by calculating the sample mean and sample covariance matrix of the posterior centrality samples:
eigen_mu_motnp <- apply(eigen_samples_motnp_std, 2, mean) 
eigen_cov_motnp <- cov(eigen_samples_motnp_std)
btwn_mu_motnp <- apply(btwn_samples_motnp_std, 2, mean) 
btwn_cov_motnp <- cov(btwn_samples_motnp_std)

# These quantities will be given to the Stan model as data to model joint posteriors of centrality in the regression. We can run a few quick plots to see how well the approximation is working and the covariance for one of the nodes:
eigen_samples_motnp_sim <- MASS::mvrnorm(1e5, eigen_mu_motnp, eigen_cov_motnp)
plot(density(eigen_samples_motnp_std[, 1]), lwd=2, main="Estimated standardised eigenvector centrality vs normal approximation", xlab="Logit edge weight")
lines(density(eigen_samples_motnp_sim[, 1]), col=rgb(0, 0, 1, 0.5), lwd=2)

btwn_samples_motnp_sim <- MASS::mvrnorm(1e5, btwn_mu_motnp, btwn_cov_motnp)
plot(density(btwn_samples_motnp_std[, 1]), lwd=2, main="Estimated standardised betweenness centrality vs normal approximation", xlab="Logit edge weight")
lines(density(btwn_samples_motnp_sim[, 1]), col=rgb(0, 0, 1, 0.5), lwd=2)

#### prior predictive checks (univariate Gaussian) ####
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

# betweenness linear effect only
range(btwn_mu_motnp)
plot(NULL, xlim = c(0,60), ylim = c(-1,6), las = 1, xlab = 'age', ylab = 'mean betweenness centrality (std)')
abline(h = range(btwn_mu_motnp), lty = 2)
for(i in 1:100){
  intercept <- rnorm(1, 2, 1.4)
  beta_age <- rnorm(1, 0, 0.02)
  lines(x = age, y = intercept + beta_age*age,
        col = rgb(0,0,1,0.4))
}

# betweenness allowing for extremes to be most/least central
plot(NULL, xlim = c(0,60),  ylim = c(-1,6), las = 1, xlab = 'age', ylab = 'mean betweenness centrality (std)')
abline(h = range(btwn_mu_motnp), lty = 2)
for(i in 1:100){
  intercept <- rnorm(1, 2, 1.4)
  beta_age <- rnorm(1, 0, 0.02)
  beta_age2 <- rnorm(1, 0, 0.001)
  lines(x = age, y = intercept + beta_age*age + beta_age2*(age^2),
        col = rgb(0,0,1,0.4))
}

#### load MOTNP data ####
colnames(true_ages) <- motnp_males$id # check this shouldn't be as.integer(as.factor(motnp_males$id)), but 90% sure this is correct -- goes into the model in order so should come out in same order?
motnp_ap <- motnp_males[motnp_males$id %in% node_ages_motnp$id, ]
true_ages_ap <- true_ages[, colnames(true_ages) %in% node_ages_motnp$id]

modeldata_ages <- as.matrix(true_ages_ap) %>% t()
modeldata_ages_sq <- (modeldata_ages)^2

#### eigenvector ~ age ####
model_data_motnp <- list(
  num_nodes = num_nodes,               # Number of dyads
  centrality_mu  = eigen_mu_motnp,     # Sample means of logit edge weights
  centrality_cov = eigen_cov_motnp,    # Sample covariance of logit edge weights
  node_age = modeldata_ages,           # Age of individual -- are these in the wrong order? should be M1, M10, M100, M101?
  node_age2 = modeldata_ages_sq        # Age squared for quadratic term
)
str(model_data_motnp)

model_eigen_cent <- stan_model("models/nodal_regression/nodal_regression_hkm_distributionage_22.07.11.stan")
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
    points(eigen_mu_motnp[i] ~ model_data_motnp$node_age[i,s[j]], pch = 19, col = rgb(0,0,0,0.05))
  }
}

# Calculate the 95% credible intervals of the model parameters
round(summary(fit_nodal_eigen)$summary[1:4, c(1:4, 8)], 5) # basically no effect whatsoever (remembering this is in standard deviations on the outcome scale but years of age on the exponent)

#### betweenness ~ age ####
model_data_motnp <- list(
  num_nodes = num_nodes,               # Number of dyads
  centrality_mu  = btwn_mu_motnp,      # Sample means of logit edge weights
  centrality_cov = btwn_cov_motnp,     # Sample covariance of logit edge weights
  node_age = modeldata_ages,           # Age of individual -- are these in the wrong order? should be M1, M10, M100, M101?
  node_age2 = modeldata_ages_sq        # Age squared for quadratic term
)
str(model_data_motnp)

model_btwn_cent <- stan_model("models/nodal_regression/nodal_regression_hkm_distributionage_betweenness_22.07.13.stan")
fit_nodal_btwn <- sampling(model_nodal_cent, data = model_data_motnp, cores = 1, chains = 1)

### check traceplot
traceplot(fit_nodal_btwn)

### calculate leave-one-out cross validation
fit_nodal_btwn

hist(summary(fit_nodal_btwn)$summary[,10], breaks = 20, xlim = c(0.99,1.02))

### posterior predictive checks
params <- rstan::extract(fit_nodal_btwn)
plot(density(btwn_samples_motnp_std[1, ]), main = "Posterior predictive density (standardised centrality)",
     col = rgb(0, 0, 0, 0.25), ylim = c(0, 1), las = 1)
for (i in 1:100) {
  j <- sample(1:num_samples, 1)                                         # select a network to plot
  lines(density(btwn_samples_motnp_std[j, ]), col=rgb(0, 0, 0, 0.25)) # plot centrality density for network j (black)
  mu <- params$beta_age[j]*apply(model_data_motnp$node_age,1,mean)     # extract age slope parameter for network j
  sigma <- btwn_cov_motnp + diag(rep(params$sigma[j], num_nodes))     # extract variance parameter for network j
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))  # plot predicted centralities for network j (blue)
} # POOR FIT -- predicted are all falling as a bell curve, observed centrality density as multimodal

### plot age against predictions
node_age <- 1:60
beta_age_mean  <- summary(fit_nodal_btwn)$summary[1,1]
beta_age_2.5   <- summary(fit_nodal_btwn)$summary[1,4]
beta_age_97.5  <- summary(fit_nodal_btwn)$summary[1,8]
beta_age2_mean <- summary(fit_nodal_btwn)$summary[2,1]
beta_age2_2.5  <- summary(fit_nodal_btwn)$summary[2,4]
beta_age2_97.5 <- summary(fit_nodal_btwn)$summary[2,8]
intercept_mean <- summary(fit_nodal_btwn)$summary[3,1]
intercept_2.5  <- summary(fit_nodal_btwn)$summary[3,4]
intercept_97.5 <- summary(fit_nodal_btwn)$summary[3,8]
mu <- intercept_mean  + beta_age_mean*node_age + beta_age2_mean*(node_age^2)
min <- intercept_2.5  + beta_age_2.5*node_age  + beta_age2_2.5*(node_age^2)
max <- intercept_97.5 + beta_age_97.5*node_age + beta_age2_97.5*(node_age^2)
summary(mu)
summary(min)
summary(max)
plot(NULL, xlim = c(0,60), ylim = c(-5, 5), lwd = 2, las = 1,
     xlab = 'age (years)', ylab = 'betweenness centrality (std)')
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
    points(btwn_mu_motnp[i] ~ model_data_motnp$node_age[i,s[j]], pch = 19, col = rgb(0,0,0,0.05))
  }
}

# Calculate the 95% credible intervals of the model parameters
round(summary(fit_nodal_btwn)$summary[1:4, c(1:4, 8)], 5) # basically no effect whatsoever (remembering this is in standard deviations on the outcome scale but years of age on the exponent)
