#### information ####
# script to extract betweenness and eigenvector centrality measures from each set of output draws
# follows workflow of Jordan Hart BISoN examples nodal regression

#### set up ####
library(rstan)
library(igraph)
library(tidyverse)
library(LaplacesDemon)

######## read in data for motnp -- first run with categories ########
# Nodal regression is a regression-based analysis often used in social network analysis to determine factors that may be associated with node-level network properties, for example if age or sex predicts network centrality. We will use edge weight posteriors from a previously-run edge weight model to conduct a nodal regression of eigenvector centrality against node type.
# BISoN adopts a fully Bayesian philosophy, so not only does it use edge weight models to estimate uncertainty over edge weights, but it also propagates this uncertainty through downstream analyses, such as nodal regression in this case.

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
# To parameterise the multivariate normal approximation, we use the sample mean and covariance matrix, calculated from the posterior centrality samples. But first we need to actually calculate the posterior centralities. We will be using eigenvector centrality in this example, which is bounded between 0 and 1. We plan on using a Gaussian family regression model as they are generally easier to interpret and fit. The response of a Gaussian regression is expected to be unbounded, so we will standardise eigenvector centrality. Centrality posteriors can be calculated and standardised using the following code:
  
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
}
head(eigen_samples_motnp)     # Unstandardised eigenvector centrality
head(eigen_samples_motnp_std) # Standardised eigenvector centrality
rm(g)

# Computing these matrices is essential for model fitting, but it’s hard to understand differences in centralities intuitively from these raw values. Instead we can visualise the centralities of each node using the following code, using tidyverse as it simplifies the process:
  
df_wide_motnp <- data.frame(eigen_samples_motnp_std)
colnames(df_wide_motnp) <- 1:num_nodes
df_long_motnp <- pivot_longer(df_wide_motnp, cols = 1:num_nodes, names_to="node_id", values_to="eigenvector_centrality")
ggplot(df_long_motnp, aes(x = eigenvector_centrality)) +
  geom_density(fill="#387780", alpha=0.7, size=0.4) +
  #facet_grid(rows=vars(as.factor(node_id)), scales="free") +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        #strip.text.y = element_text(size=12),
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

eles2 <- df_long_motnp[df_long_motnp$age_cat_id == 2,]
eles3 <- df_long_motnp[df_long_motnp$age_cat_id == 3,]
eles3 <- eles3[eles3$node_id %in% sample(eles3$node_id, size = 5, replace = F),]
eles4 <- df_long_motnp[df_long_motnp$age_cat_id == 4,]
eles4 <- eles4[eles4$node_id %in% sample(eles4$node_id, size = 5, replace = F),]
eles5 <- df_long_motnp[df_long_motnp$age_cat_id == 5,]
eles5 <- eles5[eles5$node_id %in% sample(eles5$node_id, size = 5, replace = F),]
eles6 <- df_long_motnp[df_long_motnp$age_cat_id == 6,]
eles6 <- eles6[eles6$node_id %in% sample(eles6$node_id, size = 5, replace = F),]
eles7 <- df_long_motnp[df_long_motnp$age_cat_id == 7,]
eles7 <- eles7[eles7$node_id %in% sample(eles7$node_id, size = 5, replace = F),]
plot_eles <- rbind(eles2, eles3, eles4, eles5, eles6, eles7)
ggplot(plot_eles, aes(x = eigenvector_centrality)) +
  geom_density(aes(fill = as.factor(age_cat_id)), alpha=0.7, size=0.4) +
  facet_grid(rows=vars(as.factor(node_id)), scales="free") +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
rm(eles2, eles3, eles4, eles5, eles6, eles7)
rm(plot_eles)

#### compute normal approximation ####
# The posterior centralities can now be characterised by the multivariate normal distribution as an approximation to their true posterior distributions. To highlight the importance of sampling from a multivariate normal rather than a univariate normal, we can plot centrality samples from two nodes against each other to see how they co-vary.
plot(eigen_samples_motnp_std[, 1], eigen_samples_motnp_std[, 2])

# This joint distribution contains important information that we would ideally like to include in the model. The multivariate normal will allow us to do this. Parameterising the multivariate normal approximation is a relatively simple process and can be done by calculating the sample mean and sample covariance matrix of the posterior centrality samples:
eigen_mu_motnp <- apply(eigen_samples_motnp_std, 2, mean) 
eigen_cov_motnp <- cov(eigen_samples_motnp_std)

# These quantities will be given to the Stan model as data to model joint posteriors of centrality in the regression. We can run a few quick plots to see how well the approximation is working and the covariance for one of the nodes:
eigen_samples_motnp_sim <- MASS::mvrnorm(1e5, eigen_mu_motnp, eigen_cov_motnp)
plot(density(eigen_samples_motnp_std[, 1]), lwd=2, main="Estimated standardised centrality vs normal approximation", xlab="Logit edge weight")
lines(density(eigen_samples_motnp_sim[, 1]), col=rgb(0, 0, 1, 0.5), lwd=2)

# If we’re happy with the approximation, we can now get ready to fit the model.

#### define the model ####
# The nodal regression model we’ll be using will predict node centrality using a Gaussian family model where node type is the main and only effect. Since node centralities will almost always co-vary, we need to model the joint posterior distributions over network centrality as the response variable in the regression model. This can be achieved by modelling the posterior mean centralities y_mu as outcomes of a multivariate normal with a covariance matrix y_sigma, calculated from the centrality posteriors. In Stan this looks like:
  
#  centrality_mu ~ multi_normal(predictor, centrality_cov + diag_matrix(rep_vector(square(sigma), N)));

# where predictor is the predictor term (like a + b * x in simple linear regression). Modelling centralities with a multivariate normal allows the joint uncertainty over centralities to be taken into account by the model. Weakly informative priors are used in this example, but in any real analysis they should be determined by domain knowledge and predictive checks.

model_nodal_cat <- stan_model("models/nodal_regression_hkm_categoricalage_22.07.05.stan")

#### prior predictive check -- WORK OUT HOW TO DO THIS ####
beta_nodetype <- rnorm(0, 1)
sigma <- rnorm(0, 1)
mv <- rmvnorm()










#### fit the model ####
model_data_motnp <- list(
  num_nodes = num_nodes,                                   # Number of dyads
  num_node_types = length(unique(node_ages_motnp$age_cat_id)),   # number of node types
  centrality_mu  = eigen_mu_motnp,                               # Sample means of logit edge weights
  centrality_cov = eigen_cov_motnp,                              # Sample covariance of logit edge weights
  node_age = as.integer(node_ages_motnp$age_cat_id)-1            # Integer node types corresponding to age categories (-1 because starts at 2 so now 1 = <10, 2 = 10-15, 3 = 15-20, 4 = 20-25, 5 = 25-40, and 6 = 40+)
)
str(model_data_motnp)

fit_nodal <- sampling(model_nodal_cat, data = model_data_motnp, cores = 4, chains = 4)

#### diagnostics ####
traceplot(fit_nodal)
# The chains are stationary and well-mixed so we should be able to proceed with model diagnostics.

### posterior predictive checks
# We will run a brief diagnostic check by comparing the density of expected centralities (draws of which are shown in black) against the density of predicted centralities from the regression model (draws of which are shown in blue).
params <- rstan::extract(fit_nodal)
plot(density(eigen_samples_motnp_std[1, ]), main="Posterior predictive density (standardised centrality)",
     col=rgb(0, 0, 0, 0.25), ylim=c(0, 1), las = 1)
for (i in 1:100) {
  j <- sample(1:num_samples, 1)                                        # select a network to plot
  lines(density(eigen_samples_motnp_std[j, ]), col=rgb(0, 0, 0, 0.25)) # plot centrality density for network j (black)
  mu <- params$beta_nodetype[j, model_data_motnp$node_age]                   # extract nodetype slope parameter for network j
  sigma <- eigen_cov_motnp + diag(rep(params$sigma[j], num_nodes))           # extract variance parameter for network j
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))  # plot predicted centralities for network j (blue)
}
# If almost all of the expected (black) densities fall within the range expected from the predicted (blue) densities, then the distribution of centralities appears to have been captured reasonably by the regression model. If the regression model indicates a lot of variance around predicted values, this can suggest potential scope for improving the regression model.

#### interpreting the model ####
# Calculate the 95% credible intervals of the model parameters
round(summary(fit_nodal)$summary[1:7, c(1, 4, 8)], 2)

# At this point it becomes clear that the regression we’ve conducted is not exactly the same as what might be expected from standard frequentist regressions, where categories are interpreted relative to a reference category. Instead, a parameter is estimated for each category, and we can use contrasts to calculate the magnitude of differences between categories of interest. Contrasts are easily calculated as the statistic of interest from the posteriors of the model. NOTE: THESE VALUES ARE IN STANDARD DEVIATIONS
beta_diff_1_2 <- params$beta_nodetype[, 1] - params$beta_nodetype[, 2]
plot(density(beta_diff_1_2), main="Posterior difference between nodes of age <10 and 10-15")
abline(v=0, lty=2)
(beta_diff_summary_1_2 <- round(quantile(beta_diff_1_2, probs=c(0.5, 0.025, 0.975)), 2))

beta_diff_1_3 <- params$beta_nodetype[, 1] - params$beta_nodetype[, 3]
plot(density(beta_diff_1_3), main="Posterior difference between nodes of age <10 and 15-20")
abline(v=0, lty=2)
(beta_diff_summary_1_3 <- round(quantile(beta_diff_1_3, probs=c(0.5, 0.025, 0.975)), 2))

beta_diff_1_4 <- params$beta_nodetype[, 1] - params$beta_nodetype[, 4]
plot(density(beta_diff_1_4), main="Posterior difference between nodes of age <10 and 20-25")
abline(v=0, lty=2)
(beta_diff_summary_1_4 <- round(quantile(beta_diff_1_4, probs=c(0.5, 0.025, 0.975)), 2))

beta_diff_1_5 <- params$beta_nodetype[, 1] - params$beta_nodetype[, 5]
plot(density(beta_diff_1_5), main="Posterior difference between nodes of age <10 and 25-40")
abline(v=0, lty=2)
(beta_diff_summary_1_5 <- round(quantile(beta_diff_1_5, probs=c(0.5, 0.025, 0.975)), 2))

beta_diff_1_6 <- params$beta_nodetype[, 1] - params$beta_nodetype[, 6]
plot(density(beta_diff_1_6), main="Posterior difference between nodes of age <10 and 40+")
abline(v=0, lty=2)
(beta_diff_summary_1_6 <- round(quantile(beta_diff_1_6, probs=c(0.5, 0.025, 0.975)), 2))

beta_diff_2_3 <- params$beta_nodetype[, 2] - params$beta_nodetype[, 3]
plot(density(beta_diff_2_3), main="Posterior difference between nodes of age 10-15 and 15-20")
abline(v=0, lty=2)
(beta_diff_summary_2_3 <- round(quantile(beta_diff_2_3, probs=c(0.5, 0.025, 0.975)), 2))

beta_diff_2_4 <- params$beta_nodetype[, 2] - params$beta_nodetype[, 4]
plot(density(beta_diff_2_4), main="Posterior difference between nodes of age 10-15 and 20-25")
abline(v=0, lty=2)
(beta_diff_summary_2_4 <- round(quantile(beta_diff_2_4, probs=c(0.5, 0.025, 0.975)), 2))

beta_diff_2_5 <- params$beta_nodetype[, 2] - params$beta_nodetype[, 5]
plot(density(beta_diff_2_5), main="Posterior difference between nodes of age 10-15 and 25-40")
abline(v=0, lty=2)
(beta_diff_summary_2_5 <- round(quantile(beta_diff_2_5, probs=c(0.5, 0.025, 0.975)), 2))

beta_diff_2_6 <- params$beta_nodetype[, 2] - params$beta_nodetype[, 6]
plot(density(beta_diff_2_6), main="Posterior difference between nodes of age 10-15 and 40+")
abline(v=0, lty=2)
(beta_diff_summary_2_6 <- round(quantile(beta_diff_2_6, probs=c(0.5, 0.025, 0.975)), 2))

beta_diff_3_4 <- params$beta_nodetype[, 3] - params$beta_nodetype[, 4]
plot(density(beta_diff_3_4), main="Posterior difference between nodes of age 15-20 and 20-25")
abline(v=0, lty=2)
(beta_diff_summary_3_4 <- round(quantile(beta_diff_3_4, probs=c(0.5, 0.025, 0.975)), 2))

beta_diff_3_5 <- params$beta_nodetype[, 3] - params$beta_nodetype[, 5]
plot(density(beta_diff_3_5), main="Posterior difference between nodes of age 15-20 and 25-40", xlim = c(0,2))
abline(v=0, lty=2)
(beta_diff_summary_3_5 <- round(quantile(beta_diff_3_5, probs=c(0.5, 0.025, 0.975)), 2))

beta_diff_3_6 <- params$beta_nodetype[, 3] - params$beta_nodetype[, 6]
plot(density(beta_diff_3_6), main="Posterior difference between nodes of age 15-20 and 40+")
abline(v=0, lty=2)
(beta_diff_summary_3_6 <- round(quantile(beta_diff_3_6, probs=c(0.5, 0.025, 0.975)), 2))

beta_diff_4_5 <- params$beta_nodetype[, 4] - params$beta_nodetype[, 5]
plot(density(beta_diff_4_5), main="Posterior difference between nodes of age 20-25 and 25-40")
abline(v=0, lty=2)
(beta_diff_summary_4_5 <- round(quantile(beta_diff_4_5, probs=c(0.5, 0.025, 0.975)), 2))

beta_diff_4_6 <- params$beta_nodetype[, 4] - params$beta_nodetype[, 6]
plot(density(beta_diff_4_6), main="Posterior difference between nodes of age 20-25 and 40+")
abline(v=0, lty=2)
(beta_diff_summary_4_6 <- round(quantile(beta_diff_4_6, probs=c(0.5, 0.025, 0.975)), 2))

beta_diff_5_6 <- params$beta_nodetype[, 5] - params$beta_nodetype[, 6]
plot(density(beta_diff_5_6), main="Posterior difference between nodes of age 25-40 and 40+")
abline(v=0, lty=2)
(beta_diff_summary_5_6 <- round(quantile(beta_diff_5_6, probs=c(0.5, 0.025, 0.975)), 2))

######## read in data for motnp -- now run with continuous ########
### clean environment ###
rm(df_long_motnp, df_wide_motnp, )

### import data for aggregated model (binomial)
df_agg_motnp

### load the edge weights
logit_edge_samples_motnp

#### calculate posterior centralities ####
edge_samples_motnp

# Build adjacency tensor
num_nodes
num_samples
adj_tensor

# Calculate centrality and store posterior samples in a matrix
eigen_samples_motnp
eigen_samples_motnp_std



ggplot(plot_eles, aes(x = eigenvector_centrality)) +
  geom_density(aes(fill = AGE), alpha=0.7, size=0.4) +
  facet_grid(rows=vars(as.factor(node_id)), scales="free") +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

#### compute normal approximation ####
eigen_mu_motnp
eigen_cov_motnp
eigen_samples_motnp_sim

#### define the model ####
model_nodal_cont <- stan_model("models/nodal_regression_hkm_continuousage_22.07.06.stan")

#### prior predictive check ####













#### fit the model -- FOR NOW, HAVE SELECTED A SINGLE VALUE FROM THE CENTRE OF EACH AGE CATEGORY. THESE VALUES ARE ABSOLUTELY NOT WHAT SHOULD GO INTO THE REAL MODEL ####
node_ages_motnp$age_continuous <- ifelse(node_ages_motnp$age_cat_id == 2, 7,
                                         ifelse(node_ages_motnp$age_cat_id == 3, 12,
                                                ifelse(node_ages_motnp$age_cat_id == 4, 18,
                                                       ifelse(node_ages_motnp$age_cat_id == 5, 22,
                                                              ifelse(node_ages_motnp$age_cat_id == 6, 32, 48)))))
node_ages_motnp$age_continuous <- NA
for(i in 1:nrow(node_ages_motnp)) {
  node_ages_motnp$age_continuous[i] <- ifelse(node_ages_motnp$age_cat_id[i] == 2, runif(1,1,10),
                                         ifelse(node_ages_motnp$age_cat_id[i] == 3, runif(1,10,15),
                                                ifelse(node_ages_motnp$age_cat_id[i] == 4, runif(1,15,20),
                                                       ifelse(node_ages_motnp$age_cat_id[i] == 5, runif(1,20,25),
                                                              ifelse(node_ages_motnp$age_cat_id[i] == 6, runif(1,25,40),
                                                                     runif(1,40,60))))))
}
model_data_motnp <- list(
  num_nodes = num_nodes,                               # Number of dyads
  centrality_mu  = eigen_mu_motnp,                     # Sample means of logit edge weights
  centrality_cov = eigen_cov_motnp,                    # Sample covariance of logit edge weights
  node_age = node_ages_motnp$age_continuous,            # Age of individual (NOTE -- NOT CALCULATED HOW THEY SHOULD BE!! PURELY FOR TESTING AND DEBUGGING MODEL)
  node_age2 = (node_ages_motnp$age_continuous)^2
)
str(model_data_motnp)

plot(model_data_motnp$centrality_mu ~ model_data_motnp$node_age, las = 1, pch = 19, col = rgb(0,0,1,0.3))

fit_nodal_cont <- sampling(model_nodal_cont, data = model_data_motnp, cores = 4, chains = 4)

#### diagnostics ####
traceplot(fit_nodal_cont)

### posterior predictive checks
params <- rstan::extract(fit_nodal_cont)
plot(density(eigen_samples_motnp_std[1, ]), main = "Posterior predictive density (standardised centrality)",
     col = rgb(0, 0, 0, 0.25), ylim = c(0, 1), las = 1)
for (i in 1:100) {
  j <- sample(1:num_samples, 1)                                        # select a network to plot
  lines(density(eigen_samples_motnp_std[j, ]), col=rgb(0, 0, 0, 0.25)) # plot centrality density for network j (black)
  mu <- params$beta_age[j, model_data_motnp$node_age]                  # extract age slope parameter for network j
  sigma <- eigen_cov_motnp + diag(rep(params$sigma[j], num_nodes))     # extract variance parameter for network j
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))  # plot predicted centralities for network j (blue)
}

#### interpreting the model ####
# Calculate the 95% credible intervals of the model parameters
round(summary(fit_nodal_cont)$summary[1:3, c(1:4, 8)], 2)
