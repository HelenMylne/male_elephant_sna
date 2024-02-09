## set up #### 
library(LaplacesDemon)
library(tidyverse)
library(cmdstanr)
#set_cmdstan_path('../packages/.cmdstan2/cmdstan-2.33.1/')

set.seed(1)

## simulate population ####
## define population parameters
min_age <- 11                                                # youngest individual
max_age <- 60                                                #  oldest  individual
n_nodes <- 2*((max_age+1) - min_age)                         # total nodes = average 2 per age

## simulate population for first time window
sim <- data.frame(node = 1:n_nodes,                          # create data frame of individuals
                  age = sample(x = min_age:max_age, size = n_nodes,
                               replace = T, prob = 1/(min_age:max_age)),
                  age_cat = NA,
                  mu = NA, sd = NA) %>% 
  mutate(age_cat = ifelse(age <= 15, 1,
                          ifelse(age <= 20, 2,
                                 ifelse(age <= 25, 3,
                                        ifelse(age <= 40, 4, 5)))))

## simulate centralities ####
## simulate age effect
sim_slope <- 0.1 #c(-1.5,-0.5,0.5,1.5,2.5)                    # set age effect -- smaller = bigger impact on invlogit scale as values large
sim_intcp <- -4
sim$mu <- sim$age * sim_slope + sim_intcp #sim$mu <- sim_slope[sim$age_cat] + sim_intcp       # simulate mean centrality on normal scale
plot(sim$mu ~ sim$age) #plot(sim$mu ~ sim$age_cat)           # plot
#sim$mu_std <- ( sim$mu - mean(sim$mu) ) / sd(sim$mu)

## simulate full distribution of samples per node
sim$sd <- 1#abs(sim_slope/3)           # make small to start with to be sure model should be able to detect difference
sim_dat <- matrix(data = NA, nrow = 4000, ncol = n_nodes, dimnames = list(NULL, sim$node_window))    # create matrix
for(j in 1:n_nodes){
  sim_dat[,j] <- rnorm(n = nrow(sim_dat), mean = sim$mu[j], sd = sim$sd[j])  # simulate distribution
}
plot(sim_dat[1,] ~ sim$age)          # plot simulated values against age

# ## standardise
# sim_dat_std <- sim_dat               # create matrix to fill
# for(i in 1:nrow(sim_dat_std)){
#   sim_dat_std[i,] <- (sim_dat[i,] - mean(sim_dat[i,]) ) / sd(sim_dat[i,]) # standardise values
# }
# plot(sim_dat_std[1,] ~ sim$age)      # plot simulated values against age

## visualise
data.frame(sim_dat) %>% 
  pivot_longer(cols = everything(),
               names_to = "node", values_to = "centrality") %>% 
  separate(node, into = c('X','node'), remove = T, sep = 1) %>% 
  dplyr::select(-X) %>% 
  mutate(node = as.integer(node)) %>% 
  left_join(sim[,c('node','age')], by = 'node') %>% 
  filter(node %in% seq(1, n_nodes, by = 2)) %>% 
  ggplot(aes(x = centrality, fill = age)) +
  geom_density(linewidth = 0.4) +
  facet_grid(rows = vars(as.factor(node)), scales = "free") +
  labs(x = "Eigenvector centrality (standardised)") + 
  theme_void() + 
  theme(strip.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

data.frame(sim_dat) %>% 
  pivot_longer(cols = everything(),
               names_to = "node", values_to = "centrality") %>% 
  separate(node, into = c('X','node'), remove = T, sep = 1) %>% 
  dplyr::select(-X) %>% 
  mutate(node = as.integer(node)) %>% 
  left_join(sim[,c('node','age_cat')], by = 'node') %>% 
  filter(node %in% seq(1, n_nodes, by = 2)) %>% 
  ggplot(aes(x = centrality, fill = as.factor(age_cat))) +
  geom_density(linewidth = 0.4) +
  facet_grid(rows = vars(as.factor(node)), scales = "free") +
  labs(x = "Eigenvector centrality (standardised)") + 
  scale_fill_viridis_d() +
  theme_void() + 
  theme(strip.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

## normal approximation ####
## normal approximation
plot(sim_dat[,which(sim$age == min(sim$age))[1]],
     sim_dat[,which(sim$age == max(sim$age))[1]])   # plot covariance (oldest and youngest to be sure it works for all pairs)
sim_cent_mu <- apply(sim_dat, 2, mean)     # calculate means per node
sim_cent_cov <- cov(sim_dat)               # calculate covariance matrix

## check normal approximation
sim_cent_samples <- MASS::mvrnorm(1e5, sim_cent_mu, sim_cent_cov)     # simulate from multivariate normal
plot(density(sim_dat[, 1]), lwd = 2, las = 1,                     # plot true density curve
     main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
lines(density(sim_cent_samples[, 1]), col = rgb(0,0,1,0.5), lwd = 2)  # overlay normal approximation

## prior predictive check ####
n <- 100
beta_age <- rnorm(n, 0, 1)
intercept  <- rnorm(n, 0, 1.5)
age_dirichlet <- rdirichlet(n, c(1,1,1,1,1))
plot(NULL, las = 1, xlab = 'age (standardised)', ylab = 'eigenvector (standardised)',
     ylim = c(min(sim_cent_mu)-1, max(sim_cent_mu)+1), xlim = c(min(sim$age_cat), max(sim$age_cat)))
abline(h = min(sim_cent_mu), lty = 2) ; abline(h = max(sim_cent_mu), lty = 2)
x <- min(sim$age_cat):max(sim$age_cat)
for(i in 1:n){
  y <- rep(NA, length(x))
  for(j in 1:length(x)){
    y[j] <- intercept[i] + beta_age[i]*sum(age_dirichlet[i,][1:x[j]])
  }
  lines(x = x, y = y, col = rgb(0,0,1,0.4))
}
rm(n, beta_age, intercept, age_dirichlet, x, y) ; gc()

## run model -- age as an ordered categorical variable with 1 window ####
## create data
n_age_cat <- length(unique(sim$age_cat))
eigen_list <- list(num_nodes = n_nodes,
                   num_age_cat = n_age_cat,
                   length_dirichlet = n_age_cat + 1,
                   centrality_mu = sim_cent_mu,
                   centrality_cov = sim_cent_cov,
                   node_age = sim$age_cat,
                   prior_age = rep(1, n_age_cat))

## check inputs
plot(sim_cent_mu ~ sim$age_cat)

## load model
nodal_regression <- cmdstan_model('models/eigen_regression_motnp.stan')

## run model
n_chains <- 4
n_samples <- 1000
fit_sim <- nodal_regression$sample(data = eigen_list,
                                   chains = n_chains, parallel_chains = n_chains,
                                   iter_warmup = n_samples, iter_sampling = n_samples)

## check outputs ####
## view summary
(summary <- fit_sim$summary())
par(mfrow = c(3,1))
hist(summary$rhat, breaks = 50)
hist(summary$ess_bulk, breaks = 50)
hist(summary$ess_tail, breaks = 50)
par(mfrow = c(1,1))

## extract posterior
params <- fit_sim$draws(format = 'draws_df')
delta <- params %>% 
  select(`delta[1]`,`delta[2]`,`delta[3]`,`delta[4]`,`delta[5]`,
         `.chain`,`.iteration`,`.draw`)
delta_j <- params %>% 
  select(`delta_j[1]`,`delta_j[2]`,`delta_j[3]`,`delta_j[4]`,`delta_j[5]`,`delta_j[6]`,
         `.chain`,`.iteration`,`.draw`)

## traceplot linear effect size
parameters_to_check <- c('intercept','beta_age','sigma','predictor[1]','predictor[20]','predictor[50]')
params %>% 
  select(all_of(parameters_to_check), .chain, .iteration, .draw) %>%
  pivot_longer(cols = all_of(parameters_to_check), 
               names_to = 'parameter', values_to = 'value') %>%
  rename(chain_position = .iteration,
         chain = .chain,
         draw = .draw) %>% 
  #filter(chain == 4) %>% 
  ggplot(aes(x = chain_position, y = value, colour = as.factor(chain)))+ # beta_age has fully mixed, but more stretch on one side than the other
  geom_line()+
  facet_wrap(. ~ parameter, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'none')
delta %>% 
  pivot_longer(cols = c(`delta[1]`,`delta[2]`,`delta[3]`,`delta[4]`,`delta[5]`), 
               names_to = 'parameter', values_to = 'value') %>% 
  rename(chain_position = .iteration,
         chain = .chain,
         draw = .draw) %>% 
  #filter(chain == 4) %>% 
  ggplot(aes(x = chain_position, y = value, colour = as.factor(chain)))+
  geom_line()+
  facet_wrap(. ~ parameter, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'none')
delta_j %>% 
  pivot_longer(cols = c(`delta_j[1]`,`delta_j[2]`,`delta_j[3]`,`delta_j[4]`,`delta_j[5]`,`delta_j[6]`), 
               names_to = 'parameter', values_to = 'value') %>% 
  rename(chain_position = .iteration,
         chain = .chain,
         draw = .draw) %>% 
  #filter(chain == 4) %>% 
  ggplot(aes(x = chain_position, y = value, colour = as.factor(chain)))+
  geom_line()+
  facet_wrap(. ~ parameter, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'none')

## posterior predictive check ####
plot(density(sim_dat[1, ]), las = 1, ylim = c(0,0.4),
     main = "Posterior predictive check (standardised centrality):\nblack = data, blue = predicted",
     col=rgb(0, 0, 0, 0.25))
for (i in 1:100) {
  j <- sample(1:length(params$beta_age), 1)
  lines(density(sim_dat[j, ]), col=rgb(0, 0, 0, 0.25))
  mu <- rep(NA, length(eigen_list$node_age))
  for(k in 1:length(eigen_list$node_age)){
    mu[k] <- params$intercept[j] + params$beta_age[j]*sum(delta_j[j,(1:eigen_list$node_age[k])])
  }
  sigma <- sim_cent_cov + diag(rep(params$sigma[j], n_nodes))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}

## predict from model ####
# create functions for predictions
predict_mean_centrality <- function(params, delta_j, pred_data){
  mu_matrix <- matrix(data = NA,
                      nrow = length(params$beta_age),
                      ncol = nrow(pred_data),
                      dimnames = list(1:length(params$beta_age),
                                      pred_data$node))
  for(i in 1:nrow(mu_matrix)){
    for(j in 1:ncol(mu_matrix)){
      mu_matrix[i,j] <- params$intercept[i] + params$beta_age[i] * sum(delta_j[i,(1:pred_data$age_cat[j])])
    }
  }
  return(mu_matrix)
}
predict_full_centrality <- function(params, mu_matrix, cov_matrix){
  full_matrix <- mu_matrix
  for(i in 1:nrow(full_matrix)){
    full_matrix[i,] <- MASS::mvrnorm(n = 1, mu = mu_matrix[i,],
                                      Sigma = cov_matrix + diag(rep(params$sigma[i], ncol(full_matrix))))
  }
  return(full_matrix)
}

## create predictions data frame
pred_data <- sim[order(sim$age_cat),]

## predicted means
predict_mu <- predict_mean_centrality(params = params, delta_j = delta_j, pred_data = pred_data)
pred_data$mu_mean_predict <- apply(predict_mu, 2, mean)
pred_data$mu_lwr_predict <- apply(predict_mu, 2, quantile, probs = 0.025)
pred_data$mu_upr_predict <- apply(predict_mu, 2, quantile, probs = 0.975)

## full predictions
predict_full <- predict_full_centrality(params = params, mu_matrix = predict_mu, cov_matrix = sim_cent_cov)
pred_data$full_lwr_predict <- apply(predict_full, 2, quantile, probs = 0.025)
pred_data$full_upr_predict <- apply(predict_full, 2, quantile, probs = 0.975)

## compare to raw data
ggplot(pred_data)+
  geom_ribbon(aes(x = age_cat, ymin = full_lwr_predict, ymax = full_upr_predict),
              alpha = 0.3, fill = 'purple')+
  geom_ribbon(aes(x = age_cat, ymin = mu_lwr_predict, ymax = mu_upr_predict),
              alpha = 0.3, fill = 'blue')+
  geom_point(aes(x = age_cat, y = mu))+
  geom_line(aes(x = age_cat, y = mu_mean_predict))+
  theme_classic()
ggplot(pred_data,
       aes(y = mu_mean_predict, x = mu))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)

# ## unstandardise predictions
# predict_mu_ustd <- predict_mu * sd(sim$mu) + mean(sim$mu)
# predict_full_ustd <- predict_full * sd(sim$mu) + mean(sim$mu)
# pred_data$mu_mean_predict_ustd <- apply(predict_mu_ustd, 2, mean)
# pred_data$mu_lwr_predict_ustd <- apply(predict_mu_ustd, 2, quantile, probs = 0.025)
# pred_data$mu_upr_predict_ustd <- apply(predict_mu_ustd, 2, quantile, probs = 0.975)
# pred_data$full_lwr_predict_ustd <- apply(predict_full_ustd, 2, quantile, probs = 0.025)
# pred_data$full_upr_predict_ustd <- apply(predict_full_ustd, 2, quantile, probs = 0.975)

## check
plot(pred_data$mu_mean_predict ~ pred_data$mu)
abline(a = 0, b = 1)

# ## compare to unstandardised raw data
# plot(NULL, las = 1,  type = 'p', ylim = c(-15,2), xlim = c(1,5),
#      xlab = 'age category', ylab = 'eigenvector centrality (unstd)')    # set up plot
# polygon(y = c(pred_data$full_lwr_predict_ustd, rev(pred_data$full_upr_predict_ustd)),
#         x = c(pred_data$age_cat, rev(pred_data$age_cat)),               # plot predicted full distribution
#         col = rgb(1,1,0,0.4), border = NA)
# polygon(y = c(pred_data$mu_lwr_predict_ustd, rev(pred_data$mu_upr_predict_ustd)),
#         x = c(pred_data$age_cat, rev(pred_data$age_cat)),               # plot predicted means distribution
#         col = rgb(1,0,0,0.4), border = NA)
# points(sim$mu ~ sim$age_cat, col = rgb(0,0,1,0.2))                      # add raw points
# points(pred_data$mu_mean_predict_ustd ~ pred_data$age_cat, pch = 19)    # add mean predicted points

## convert to invlogit scale
predict_mu_invlogit <- invlogit(predict_mu)
predict_full_invlogit <- invlogit(predict_full)
pred_data$mu_mean_predict_invlogit <- apply(predict_mu_invlogit, 2, mean)
pred_data$mu_lwr_predict_invlogit <- apply(predict_mu_invlogit, 2, quantile, probs = 0.025)
pred_data$mu_upr_predict_invlogit <- apply(predict_mu_invlogit, 2, quantile, probs = 0.975)
pred_data$full_lwr_predict_invlogit <- apply(predict_full_invlogit, 2, quantile, probs = 0.025)
pred_data$full_upr_predict_invlogit <- apply(predict_full_invlogit, 2, quantile, probs = 0.975)

## compare to invlogit raw data
pred_data$mu_raw_invlogit <- invlogit(pred_data$mu)
ggplot(pred_data)+
  geom_ribbon(aes(x = age_cat, ymin = full_lwr_predict_invlogit, ymax = full_upr_predict_invlogit),
              alpha = 0.3, fill = 'purple')+
  geom_ribbon(aes(x = age_cat, ymin = mu_lwr_predict_invlogit, ymax = mu_upr_predict_invlogit),
              alpha = 0.3, fill = 'blue')+
  geom_point(aes(x = age_cat, y = LaplacesDemon::invlogit(mu)))+
  geom_line(aes(x = age_cat, y = mu_mean_predict_invlogit))+
  theme_classic()

## clean up
rm(eigen_list, nodal_regression, predict_full_invlogit, predict_mu_invlogit,#predict_full,  predict_mu,
   sigma, sim_cent_samples, summary, i, j, k, parameters_to_check) ; gc()

## extract original values from output -- NOTE: CURRENT METHOD WILL NOT WORK FOR REAL DATA AS I DON'T HAVE THE AGE VALUES IN YEARS TO GO ABOUT ADDING 1 YEAR TO EVERY ELEPHANT. WILL NEED TO COME UP WITH AN ALTERNATIVE OR JUST SHOW CONTRASTS BETWEEN AGE CATEGORIES FOR REAL DATA, BUT HERE CHECK THE EFFECT ON REAL AGE SCALE TO BE SURE THAT THE MODEL IS WORKING. ####
## set up objects to store predictions
pred_means <- list()
pred_fulls <- list()

## predict for ages from 10 years younger than current up to 10 years older: multiple opportunities to pass through age category thresholds 
for( i in -10:10 ){
  pred_data_new <- sim[order(sim$age_cat),] %>% 
    mutate(age = age + i) %>% 
    mutate(age_cat = ifelse(age <= 15, 1,
                            ifelse(age <= 20, 2,
                                   ifelse(age <= 25, 3,
                                          ifelse(age <= 40, 4, 5)))))
  pred_mu_new <- predict_mean_centrality(params = params, delta_j = delta_j, pred_data = pred_data_new)
  pred_full_new <- predict_full_centrality(params = params, mu_matrix = pred_mu_new, cov_matrix = sim_cent_cov)
  pred_means[[(i+11)]] <- pred_mu_new
  pred_fulls[[(i+11)]] <- pred_full_new
  print(i)
}
warnings() # should just be "Dropping 'draws_df' class as required metadata was removed." over and over

## create data frame to store outputs
contrasts_all <- data.frame(older = rep(-10:10, each = 21),
                            younger = rep(-10:10, 21),
                            mean_contrast = NA)

## for each combination of predictions, calculate contrast and divide difference by number of years changed
for(i in 1:21){
  for(j in 1:21){
    if(i > j){
      contrast <- (pred_fulls[[i]] - pred_fulls[[j]]) / (i-j)
      contrasts_all$mean_contrast[contrasts_all$older == (i - 11) & 
                                    contrasts_all$younger == (j - 11)] <- mean(contrast)
    }
  }
}

## remove impossible age combinations
contrasts_all <- contrasts_all %>%
  filter(is.na(mean_contrast) == FALSE)

## compare to input
sim_slope
mean(contrasts_all$mean_contrast)
quantile(contrasts_all$mean_contrast, prob = c(0.025, 0.975))

#### ignore below unless above method is actually wrong and need to come back to this
# ## get mean predictions for altered age categories
# sim2 <- sim %>% 
#   mutate(age_cat = age_cat + 1) %>%                   # shift to new age category
#   mutate(age_cat = ifelse(age_cat == 6, 1, age_cat))  # won't work with age category 6 so make these category 1 to fill in the gaps (can't actually use these)
# sim2 <- sim2[order(sim2$age_cat),]
# predict_mu2 <- predict_mean_centrality(params = params, delta_j = delta_j, pred_data = sim2)
# sim2$mu_mean_plus1 <- apply(predict_mu2, 2, mean)
# 
# ## full distribution of predictions using age_cat + 1
# predict_full2 <- predict_full_centrality(params = params, mu_matrix = predict_mu2, cov_matrix = sim_cent_cov)
# 
# ## contrast predictions on standardised scale -- check that this returns the marginal effect presented in the summary
# contrast <- predict_full2 - predict_full              # contrast between predicted values
# head(contrast[,1:5])                                  # check matrix looks right
# sim_slope                                             # original input parameter
# mean(contrast)                                        # that's not right
# quantile(contrast, prob = c(0.025, 0.975))            # very wide
# 
# ## contrast per category
# contrast_1_2 <- contrast[,sim2$age_cat == 2]
# contrast_2_3 <- contrast[,sim2$age_cat == 3]
# contrast_3_4 <- contrast[,sim2$age_cat == 4]
# contrast_4_5 <- contrast[,sim2$age_cat == 5]
# 
# mean(contrast_1_2)/5                                      # difference per year between 10-15 and 16-20
# quantile(contrast_1_2, prob = c(0.025, 0.975))            # very wide
# 
# mean(contrast_2_3)/5                                      # difference per year between 16-20 and 21-15
# quantile(contrast_2_3, prob = c(0.025, 0.975))            # very wide
# 
# mean(contrast_3_4)/10                                     # difference per year between 21-25 and 26-40
# quantile(contrast_3_4, prob = c(0.025, 0.975))            # very wide
# 
# mean(contrast_4_5)/10                                     # difference between 26-40 and 40+
# quantile(contrast_4_5, prob = c(0.025, 0.975))            # very wide
