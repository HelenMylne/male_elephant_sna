#### set up ####
# load packages
library(dplyr)
library(cmdstanr)

# set seed
set.seed(12345)

# set cmdstan path
set_cmdstan_path('../../../packages/.cmdstan/cmdstan-2.31.0/')

#### create model inputs ####
### import data for aggregated model (binomial) -- counts of positive associations and total sightings
counts_df <- readr::read_csv('../../../data_processed/motnp_binomialpairwiseevents_malesonly.csv')

counts_df_model <- counts_df[, c('node_1_males','node_2_males','event_count','count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

counts_df <- counts_df[1:500,] # reduce the data for the purposes of testing only

### create data list for Stan model
edge_list <- list(
  N = nrow(counts_df),
  #dyad_ids = counts_df$dyad_males,
  together = counts_df$event_count,     # count of sightings seen together
  count_dyad = counts_df$count_dyad     # count of sightings seen
)

#### fit edge model ####
### load edge model
edge_model <- cmdstan_model("edge_binary.stan")

### fit edge model
fit <- edge_model$sample(
  data = edge_list,
  iter_warmup = 1000,
  iter_sampling = 2000,
  chains = 1,
  parallel_chains = 8
)

### extract edge draws
edge_draws <- fit$draws("edge_weight", format = "df")

#### compute normal approximation ####
# To parameterise the multivariate normal approximation, we use the sample mean and covariance matrix, calculated from the posterior edge weight samples using the following code:

### remove .chain .iteration .draw to just leave draws
edge_draws_cols <- as.data.frame(edge_draws[, 1:(ncol(edge_draws)-3)])
N <- ncol(edge_draws_cols)

### get the weights on the logit scale (they are not currently because we used a beta prior and identity link here rather than logistic link)
logit_weights <- apply(edge_draws_cols, 2, qlogis)

### fit a multivariate normal dist to the edges -- These quantities will be given to the Stan model as data to model joint posteriors of edge weight in the regression.
logit_edge_draws_mu <- apply(logit_weights, 2, mean)
logit_edge_draws_cov <- cov(logit_weights)

#### plot to see how well the approximation is working ####
### Randomly selecting samples to examine
num_samples <- 20
selected_samples <- sample(1:N, num_samples, replace = FALSE)

### Setting grid layout
rows <- floor(sqrt(num_samples))
cols <- ceiling(num_samples / rows)
par(mfrow=c(rows, cols), mar=c(2,2,2,1))

### plot
for (i in selected_samples) {
  mu <- logit_edge_draws_mu[i]
  sd <- sqrt(logit_edge_draws_cov[i,i])
  
  fitted_values_logit <- rnorm(1e5, mean=mu, sd=sd)
  fitted_values_original <- plogis(fitted_values_logit)
  
  hist(unlist(edge_draws_cols[,i]), probability=TRUE, main=paste("Dyad", i), xlab="Value", breaks=50)
  lines(density(fitted_values_original), col="blue", lw=1.5)
}

### reset plot window
par(mfrow=c(1,1))

#### dyadic regression ####
# Edge weight values (mu and cov) go into the regression model on the logit scale.
### Define the model
# The dyadic regression model we'll be using will predict the edge weight using a Gaussian family model where dyad type is the main effect, and multi-membership terms are included as random effects to account for non-independence between edges due to nodes. Since edge weights can co-vary, we need to model the joint posterior distributions over edge weights as the response variable in the regression model. This can be achieved by modelling the mean edge weights `y_mu` as a multivariate normal with a covariance matrix `y_sigma` calculated from the edge weight posterior. In Stan this looks like:
#     logit_edge_mu ~ multi_normal(predictor, logit_edge_cov + diag_matrix(rep_vector(square(sigma), N)));
# where `predictor` is the predictor term (like `a + b * x` in simple linear regression). Modelling edge weights with a multivariate normal allows the joint uncertainty over edge weights to be taken into account by the model. Weakly informative priors are used in this example, but in any real analysis they should be determined by domain knowledge and predictive checks.

### load dyadic regression model
model_dyadic <- cmdstan_model("dyadic_regression.stan")

### create data list for dyadic regression
all_node_IDs <- unique(c(counts_df$id_1, counts_df$id_2))
edge_model_list <- list(
  num_dyads = length(unique(counts_df$dyad_id)),   # Number of dyads
  K = length(all_node_IDs),                # Number of nodes
  logit_edge_mu = logit_edge_draws_mu,     # Sample means of the logit edge weights
  logit_edge_cov = logit_edge_draws_cov,   # Sample covariance of logit edge weights
  age_diff = counts_df$age_diff,           # designed as a factor variable, this is a continuous
  id_1 = as.integer(factor(counts_df$id_1, levels = all_node_IDs)), # Node IDs for multimembership effects
  id_2 = as.integer(factor(counts_df$id_2, levels = all_node_IDs)),
  jitter = 1e-6                            # jitter to add to the diag of the cov matrix for numerical stability
)

### Fit the dyadic regression model
fit_dyadic <- model_dyadic$sample(
  data = edge_model_list,
  iter_warmup = 1000,
  iter_sampling = 2000,
  chains = 4,
  parallel_chains = 8
)

### check model fit to logit(edge_weight) -- when looking at marginal effects, etc., you might want to transform back to the original scale for better interpretability.
fit_dyadic
