#### Information ####
# script takes data input from edge weight estimation for MOTNP population (input = 1000 draws for each of 4 chains per dyad posterior distribution)
# runs through dyadic regression as specified by Jordan Hart in BISoN examples (https://github.com/JHart96/bison_examples/blob/main/examples/dyadic_regression_stan.md)

#### Set up ####
library(tidyverse)
library(rstan)
library(car)
#library(rethinking)

#library(tidyverse, lib.loc = 'packages/')
#library(rstan, lib.loc = 'packages/')
#library(car, lib.loc = 'packages/')
#library(rethinking, lib.loc = 'packages/')

#### Read in data ####
#df_agg <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_bayesian_allpairwiseevents_splitbygrouptype_22.01.13.csv')
df_agg <- read_csv('data_processed/motnp_bayesian_binomialpairwiseevents.csv')
str(df_agg)
df_agg$dem_type <- paste(sort(df_agg$dem_class_1, dem_class_2), sep = '_')
#df_agg$dem_type <- ifelse(df_agg$dem_type == 'PM_AM', 'AM_PM', df_agg$dem_type)
df_agg <- df_agg[df_agg$dem_type == 'AM_AM' | df_agg$dem_type == 'AM_PM' | df_agg$dem_type == 'PM_PM',]

#ages <- readRDS('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_ageestimates_mcmcoutput.rds')
ages <- readRDS('data_processed/motnp_ageestimates_mcmcoutput.rds')
males <- unique(c(df_agg$id_1,df_agg$id_2))
ages <- as.data.frame(ages[,males])
rm(males)

#edge_samples <- readRDS('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_edgeweightestimates_mcmcoutput.rds') %>% 
#  select(-`1.lp__`)
edge_samples <- readRDS('data_processed/motnp_edgeweightestimates_mcmcoutput.rds') %>% 
  select(-`1.lp__`)
edge_males <- edge_samples[, colnames(edge_samples) %in% df_agg$dyad]
str(edge_males)
logit_edge_males <- car::logit(edge_males)
rm(edge_samples)#, edge_males)

#### define ages and age differences ####
# old version using categories:
#df_agg$age1 <- ifelse(df_agg$age_category_1 == '40+', 7,
#                      ifelse(df_agg$age_category_1 == '25-40', 6,
#                             ifelse(df_agg$age_category_1 == '20-25', 5,
#                                    ifelse(df_agg$age_category_1 == '15-19', 4,
#                                           ifelse(df_agg$age_category_1 == '10-15', 3, 
#                                                  ifelse(df_agg$age_category_1 == '9-10', 2, 'error'))))))
#df_agg$age2 <- ifelse(df_agg$age_category_2 == '40+', 7,
#                      ifelse(df_agg$age_category_2 == '25-40', 6,
#                             ifelse(df_agg$age_category_2 == '20-25', 5,
#                                    ifelse(df_agg$age_category_2 == '15-19', 4,
#                                           ifelse(df_agg$age_category_2 == '10-15', 3,
#                                                  ifelse(df_agg$age_category_2 == '9-10', 2, 'error'))))))
#df_agg$dyad_type <- paste0(df_agg$age1, df_agg$age2)
#table(df_agg$dyad_type)

#df_agg$age1 <- as.numeric(df_agg$age1)
#df_agg$age2 <- as.numeric(df_agg$age2)
#df_agg$age_diff <- abs(df_agg$age1 - df_agg$age2)

#df_agg$age_diff_cat <- ifelse(df_agg$age_diff < 2,
#                              ifelse(mean(c(df_agg$age1, df_agg$age2)) <= 4, 'YY',
#                                     ifelse(mean(c(df_agg$age1, df_agg$age2)) >= 6, 'OO', 'MM')),
#                              ifelse(df_agg$age_diff >= 4, 'distant', 'mid'))

# new version using posterior draws
df_agg$age_mean_1 <- NA
df_agg$age_mean_2 <- NA
for(i in 1:nrow(df_agg)) {
  age1 <- ages[,df_agg$id_1[i]]
  age2 <- ages[,df_agg$id_2[i]]
  df_agg$age_mean_1[i] <- mean(age1)
  df_agg$age_mean_2[i] <- mean(age2)
  rm(age1,age2)
}
df_agg$mean_diff <- abs(df_agg$age_mean_1 - df_agg$age_mean_2)

#### plot -- do not run unless necessary (takes a LOT of RAM and will usually need to kill all other processes) ####
plot_df <- pivot_longer(edge_males, cols = everything(), names_to = 'dyad', values_to = 'draw_value')
num_dyads <- length(unique(plot_df$dyad))
plot_df$draw_id <- rep(1:4000, each = num_dyads)
df_agg_short <- df_agg[,c(27,1:3,16:17,22:23,28:30)] # this still seems to have the mo_events/mx_events/bh_events columns so these values may change
plot_df2 <- left_join(plot_df, df_agg_short, by = 'dyad')

plot_df2$age_category_1 <- factor(plot_df2$age_category_1, levels = c("9-10","10-15","15-19","20-25","25-40","40+"))
plot_df2$age_category_2 <- factor(plot_df2$age_category_2, levels = c("9-10","10-15","15-19","20-25","25-40","40+"))

ggplot(plot_df2[plot_df2$draw_id < 1000,], aes(x = draw_value, group = dyad))+
  geom_density(colour = rgb(0,0,1,0.1))+
  facet_grid(age_category_1 ~ age_category_2)+
  labs(x = 'edge weight')+
  theme_light()+
  theme(strip.background = element_rect(fill = rgb(0,0,1,0.6)))

#### Defining the model ####
# The dyadic regression model we’ll be using will predict the edge weight using a Gaussian family model where dyad type is the main effect, and multi-membership terms are included as random effects to account for non-independence between edges due to nodes. Since edge weights can co-vary, we need to model the joint posterior distributions over edge weights as the response variable in the regression model. This can be achieved by modelling the mean edge weights y_mu as a multivariate normal with a covariance matrix y_sigma calculated from the edge weight posterior. In Stan this looks like:
  
#  logit_edge_mu ~ multi_normal(predictor, logit_edge_cov + diag_matrix(rep_vector(square(sigma), N)));

# where predictor is the predictor term (like a + b * x in simple linear regression). Modelling edge weights with a multivariate normal allows the joint uncertainty over edge weights to be taken into account by the model. Weakly informative priors are used in this example, but in any real analysis they should be determined by domain knowledge and predictive checks.

### define the priors:
# Hypothesis 0: Age has no impact on association.
# Hypothesis 1: Elephants will prefer to associate with age mates because they are more similar and are available as sparring partners. Effect of age difference will therefore be negative so a greater age gap = smaller probability of association.
# Hypothesis 2: Young elephants will prefer to associate with older males because they can act as sources of information. Old elephants will prefer to associate with young males because they can control their musth and don't act as sources of competition. Effect of age difference will therefore be positive so a greater age gap = greater probability of association.
# Hypothesis 3: Some combination of Hypotheses 1 and 2 is true. Effect of age difference will therefore have quite high variance and uncertainty, likely spanning 0 due to preferences in both directions.
# Hypothesis 4: Young elephants will prefer to associate with other young elephants for sparring and because they may be more familiar with individuals closer in age through their maternal herds (esp. in MOTNP). Older males will show no preference as they are already established in the environment and may not require sparring partners. There will therefore be an interaction effect between average age of the dyad and the age difference between them.

df_agg$age_diff_std <- standardize(df_agg$mean_diff)
df_agg$E <- apply(logit_edge_males, 2, mean)
a <- rnorm(100, 0.5, 0.2)
bA <- rnorm(100, 0, 0.1)

plot(NULL, las = 1,
     xlim = c(min(df_agg$age_diff_std),max(df_agg$age_diff_std)), xlab = 'age category difference',
     ylim = c(-0.5,1.5), ylab = 'edge weight')
abline(h = c(0,1), lty = 2)
for( i in 1:100 ) curve(expr = a[i] + bA[i]*x, add = T, col = rgb(0,0,1,0.5),
                        from = min(df_agg$age_diff_std), to = max(df_agg$age_diff_std))

# load model:
#model_dyadic_cat <- stan_model('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/models/dyadic_regression_agecategorical.stan')
#model_dyadic_cat <- stan_model(source('https://raw.githubusercontent.com/HelenMylne/test/main/models/dyadic_regression_agecategorical.stan'))
model_dyadic_cat <- stan_model('models/dyadic_regression_agecategorical.stan')
#model_dyadic_cont <- stan_model('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/models/dyadic_regression.stan')
model_dyadic_cont <- stan_model('models/dyadic_regression.stan')

#### trim down data for debugging ####
test_eles <- sample(x = df_agg$id_1, size = 30, replace = F)
test_agg <- df_agg[df_agg$id_1 %in% test_eles | df_agg$id_2 %in% test_eles,]
test_ages <- ages[,colnames(ages) %in% test_eles]
test_edge <- edge_males[,colnames(edge_males) %in% test_agg$dyad]
test_logit <- logit_edge_males[,colnames(logit_edge_males) %in% test_agg$dyad]

#### Compute normal approximation ####
# To parameterise the multivariate normal approximation, we use the sample mean and covariance matrix, calculated from the posterior edge weight samples using the following code:
test_edge_mu <- apply(test_logit, 2, mean)
test_edge_cov <- cov(test_logit)

# These quantities will be given to the Stan model as data to model joint posteriors of edge weight in the regression. We can run a few quick plots to see how well the approximation is working and the covariance between one pair of edge weights.

logit_edge_newsamples <- MASS::mvrnorm(1e5, test_edge_mu, test_edge_cov)
par(mfrow=c(1, 2))
plot(density(logit_edge_newsamples[, 1]), lwd = 2, main = "Estimated logit edge weight vs normal approximation", xlab = "logit edge weight")
lines(density(logit_edge_newsamples[, 1]), col = rgb(0, 0, 1, 0.5), lwd = 2)

plot(logit_edge_newsamples[, 1], logit_edge_newsamples[, 2], col=rgb(0, 0, 1, 0.05), main="Covariance between edges 1 & 2", xlab="Edge 1 samples", ylab="Edge 2 samples")

#### Fit model -- categorical version ####
model_data <- list(
  num_dyads = nrow(test_agg),                                     # Number of dyads
  num_nodes = length(unique(c(test_agg$node_1, test_agg$node_2))),  # Number of nodes
  logit_edge_mu = test_edge_mu,                                # Sample means of edge weights
  logit_edge_cov = test_edge_cov,                              # Sample covariance of edge weights
  dyad_type = as.integer(as.factor(test_agg$dem_type))#,           # Integer dyad types corresponding to ages of pairs
  #age_1 = as.integer(as.factor(test_agg$age_class_1)),            # Mean of age estimates for individual 1
  #age_2 = as.integer(as.factor(test_agg$age_class_2)),            # Mean of age estimates for individual 2
  #node_1_id = test_agg$node_1,                                    # Node IDs for multimembership effects
  #node_2_id = test_agg$node_2                                     # Node IDs for multimembership effects
)
fit_dyadic <- sampling(model_dyadic_cat, data = model_data, cores = 1, chains = 1) # THIS WILL RUN WHEN DON'T INCLUDE THE MULITMEMBERSHIP TERMS -- INCLUDING THEM CONFUSES IT BECAUSE THE VECTOR SIZES DON'T MATCH

#### Diagnostics ####
traceplot(fit_dyadic)

#### Posterior predictive check ####
# We will run a brief diagnostic check by comparing the density of expected edge weights (draws of which are shown in black) against the density of predicted edge weights from the regression model (draws of which are shown in blue).

params <- rstan::extract(fit_dyadic)
plot(density(edge_samples[1, ]), main="Posterior predictive density of responses (edge weights)",
     ylim = c(0, 0.8), col = rgb(0, 0, 0, 0.25))
for (i in 1:100) {
  j <- sample(1:4000, 1)
  lines(density(edge_samples[j, ]), col = rgb(0, 0, 0, 0.25))
  mu <- params$beta_dyadtype[j, model_data$dyad_type] + params$mm[j, model_data$node_1_id] + params$mm[j, model_data$node_2_id]
  sigma <- model_data$edge_cov + diag(rep(params$sigma[j], model_data$N))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}

# Almost all of the expected (black) densities should fall within the range expected from the predicted (blue) densities. There are many other types of diagnostic check that could be carried out, but we won’t go into detail here. See the github repository page for more information.

