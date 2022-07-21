############ Second half of script: Use in nodal regression ############
#### load packages ####
library(tidyverse, lib.loc = 'packages/')
library(cmdstanr, lib.loc = 'packages/')
#library(ggdist, lib.loc = 'packages/')
library(posterior, lib.loc = 'packages/')
#library(bayesplot, lib.loc = 'packages/')
library(janitor, lib.loc = 'packages/')
library(readxl, lib.loc = 'packages/')
library(MASS, lib.loc = 'packages/')
library(rstan)
library(igraph, lib.loc = 'packages/')
library(tidyverse, lib.loc = 'packages/')
library(LaplacesDemon, lib.loc = 'packages/')

#### read in MPNP data ########
### import data for aggregated model (binomial)
df_agg_mpnp1 <- read_delim('data_processed/mpnp_period1_pairwiseevents_22.05.30.csv', delim = ',')
df_agg_mpnp1 <- df_agg_mpnp1[,c(1:13,20:21)]
df_agg_mpnp1$sex_1 <- 'M'
df_agg_mpnp1$age_mid_round_1 <- floor(df_agg_mpnp1$age_median_1)
df_agg_mpnp1$age_mid_round_2 <- floor(df_agg_mpnp1$age_median_2)
df_agg_mpnp1$node_1_nogaps <- as.integer(as.factor(df_agg_mpnp1$node_1))
df_agg_mpnp1$node_2_nogaps <- as.integer(as.factor(df_agg_mpnp1$node_2))+1
df_agg_mpnp1$dyad_id_nogaps <- as.integer(as.factor(df_agg_mpnp1$dyad))

### load the edge weights
rm(df, mpnp1_long, age_est_mat)
mpnp1 <- readRDS('data_processed/mpnp1_bayesian_edgedistributions_a2.b2_period1_22.05.30.rds')
mpnp1 <- mpnp1[,c(2:ncol(mpnp1))]
mpnp1 <- mpnp1[, which(colnames(mpnp1) %in% df_agg_mpnp1$dyad)]
logit_edge_samples_mpnp1 <- logit(mpnp1)

### load age data
true_ages <- readRDS('22.07.15_mpnp_agedistributions.rds')

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
for (i in 1:num_samples) {
  g <- graph_from_adjacency_matrix(adj_tensor[i, , ], mode="undirected", weighted=TRUE)
  eigen_samples_mpnp1[i, ] <- eigen_centrality(g)$vector
  eigen_samples_mpnp1_std[i, ] <- (eigen_samples_mpnp1[i, ] - mean(eigen_samples_mpnp1[i, ]))/sd(eigen_samples_mpnp1[i, ])
}
head(eigen_samples_mpnp1)     # Unstandardised eigenvector centrality
head(eigen_samples_mpnp1_std) # Standardised eigenvector centrality
rm(g)
print('centralities calculated')

# Computing these matrices is essential for model fitting, but it’s hard to understand differences in centralities intuitively from these raw values. Instead we can visualise the centralities of each node using the following code, using tidyverse as it simplifies the process:

df_wide_mpnp1 <- data.frame(eigen_samples_mpnp1_std)
colnames(df_wide_mpnp1) <- 1:num_nodes
df_long_mpnp1 <- pivot_longer(df_wide_mpnp1, cols = 1:num_nodes, names_to="node_id", values_to="eigenvector_centrality")
ggplot(df_long_mpnp1, aes(x = eigenvector_centrality)) +
  geom_density(fill="#387780", alpha=0.7, size=0.4) +
  #facet_grid(rows=vars(as.factor(node_id)), scales="free") +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        #strip.text.y = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

nodes1 <- df_agg_mpnp1[,c(3,17)] %>% distinct()
nodes2 <- df_agg_mpnp1[,c(4,18)] %>% distinct()
colnames(nodes1) <- c('id', 'age_mid_round')
colnames(nodes2) <- c('id', 'age_mid_round')
node_ages_mpnp1 <- rbind(nodes1, nodes2) %>% distinct()
rm(nodes1, nodes2)

node_ages_mpnp1$node_id <- as.integer(as.factor(node_ages_mpnp1$id))
df_long_mpnp1$node_id <- as.integer(df_long_mpnp1$node_id)
df_long_mpnp1 <- left_join(df_long_mpnp1, node_ages_mpnp1, by = 'node_id')

ggplot(df_long_mpnp1, aes(x = eigenvector_centrality)) +
  geom_density(aes(fill = age_mid_round), alpha=0.7, size=0.4) +
  facet_grid(rows=vars(as.factor(age_mid_round)), scales="free") +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
print('graphs plotted')

#### compute normal approximation ####
# The posterior centralities can now be characterised by the multivariate normal distribution as an approximation to their true posterior distributions. To highlight the importance of sampling from a multivariate normal rather than a univariate normal, we can plot centrality samples from two nodes against each other to see how they co-vary.
plot(eigen_samples_mpnp1_std[, 1], eigen_samples_mpnp1_std[, 2])

# This joint distribution contains important information that we would ideally like to include in the model. The multivariate normal will allow us to do this. Parameterising the multivariate normal approximation is a relatively simple process and can be done by calculating the sample mean and sample covariance matrix of the posterior centrality samples:
eigen_mu_mpnp1 <- apply(eigen_samples_mpnp1_std, 2, mean) 
eigen_cov_mpnp1 <- cov(eigen_samples_mpnp1_std)

# These quantities will be given to the Stan model as data to model joint posteriors of centrality in the regression. We can run a few quick plots to see how well the approximation is working and the covariance for one of the nodes:
eigen_samples_mpnp1_sim <- MASS::mvrnorm(1e5, eigen_mu_mpnp1, eigen_cov_mpnp1)
plot(density(eigen_samples_mpnp1_std[, 1]), lwd=2, main="Estimated standardised centrality vs normal approximation", xlab="Logit edge weight")
lines(density(eigen_samples_mpnp1_sim[, 1]), col=rgb(0, 0, 1, 0.5), lwd=2)

# If we’re happy with the approximation, we can now get ready to fit the model.

#### define the model ####
model_nodal_cont <- stan_model("models/nodal_regression/nodal_regression_hkm_distributionage_22.07.11.stan")

#### prior predictive check -- WORK OUT HOW TO DO THIS ####













#### fit the model ####
### read in age data
mpnp_long <- readxl::read_excel('data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211019.xlsx') %>%
  janitor::clean_names()
mpnp_long <- mpnp_long %>% 
  filter(elephant_id != '-') %>% 
  separate(elephant_id, into = c('BTF','num'), sep = 1, remove = F) %>% 
  filter(BTF == 'B')
mpnp_long <- mpnp_long[,c(1:3,6:33)]
periods <- seq(from = min(mpnp_long$date, na.rm = T),
               to = max(mpnp_long$date, na.rm = T),
               length.out = 7)
periods[7] <- periods[7]+1 # set to be one higher than the final higher otherwise it takes the last date and creates a whole new period
mpnp1_long <- mpnp_long[mpnp_long$date < periods[2],]
mpnp1_long$age_range_id <- as.numeric(mpnp1_long$age_range_id)
mpnp1_long$age_mid <- NA
for(i in 1:nrow(mpnp1_long)){
  ele <- mpnp1_long[mpnp1_long$elephant_id == mpnp1_long$elephant_id[i],]
  ele <- ele[!is.na(ele$age_range_id),]
  ele <- ele[ele$age_range_id != 10,]
  mpnp1_long$age_mid[i] <- ifelse(nrow(ele) == 0, 10, median(ele$age_range_id, na.rm = T))
  rm(ele)
}
mpnp1_long$age_mid_round <- floor(mpnp1_long$age_mid)
mpnp1_long <- mpnp1_long[mpnp1_long$age_mid_round < 10 , c(3,5,6,7,32:33)]

### use raw age data to give ID numbers for true age data
colnames(true_ages) <- mpnp1_long$elephant_id # check this shouldn't be as.integer(as.factor(mpnp1_long$id)), but 90% sure this is correct -- goes into the model in order so should come out in same order?
mpnp1_ap <- mpnp1_long[mpnp1_long$elephant_id %in% node_ages_mpnp1$id, ]
true_ages_ap <- true_ages[, colnames(true_ages) %in% node_ages_mpnp1$elephant_id]

### create data list
model_data_mpnp1 <- list(
  num_nodes = num_nodes,                          # Number of dyads
  centrality_mu  = eigen_mu_mpnp1,                # Sample means of logit edge weights
  centrality_cov = eigen_cov_mpnp1,               # Sample covariance of logit edge weights
  node_age = mpnp1_long$age_mid_round,            # Age of individual -- these qre in the wrong order? should be M1, M10, M100, M101...?
  node_age2 = (node_ages_mpnp1$age_mean)^2        # Age squared for quadratic term
)
str(model_data_mpnp1)

plot(model_data_mpnp1$centrality_mu ~ model_data_mpnp1$node_age, las = 1, pch = 19, col = rgb(0,0,1,0.3))

fit_nodal_cont <- sampling(model_nodal_cont, data = model_data_mpnp1, cores = 4, chains = 4)

#### diagnostics ####
traceplot(fit_nodal_cont)

### posterior predictive checks -- I think this is probably a sign that I've written the model wrong because I can't generate an mvnorm from this output -> params$beta_age is 1-dimensional so can't select j'th network and also ages
params <- rstan::extract(fit_nodal_cont)
plot(density(eigen_samples_mpnp1_std[1, ]), main = "Posterior predictive density (standardised centrality)",
     col = rgb(0, 0, 0, 0.25), ylim = c(0, 1), las = 1)
for (i in 1:100) {
  j <- sample(1:num_samples, 1)                                        # select a network to plot
  lines(density(eigen_samples_mpnp1_std[j, ]), col=rgb(0, 0, 0, 0.25)) # plot centrality density for network j (black)
  mu <- params$beta_age[j]#, model_data_mpnp1$node_age]                  # extract age slope parameter for network j
  sigma <- eigen_cov_mpnp1 + diag(rep(params$sigma[j], num_nodes))     # extract variance parameter for network j
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))  # plot predicted centralities for network j (blue)
}

#### interpreting the model ####
# Calculate the 95% credible intervals of the model parameters
round(summary(fit_nodal_cont)$summary[1:3, c(1:4, 8)], 2) # basically no effect whatsoever
