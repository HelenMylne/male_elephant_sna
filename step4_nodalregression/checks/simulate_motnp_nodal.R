## set up #### 
library(LaplacesDemon)
library(tidyverse)
library(cmdstanr)

## simulate population ####
## define population parameters
min_age <- 11                                                # youngest individual
max_age <- 60                                                #  oldest  individual
n_nodes <- ((max_age+1) - min_age)                           # total nodes = 2 per age

## simulate population for first time window
sim <- data.frame(node = 1:n_nodes,                          # create data frame of individuals
                  age = min_age:max_age,
                  mu = NA, sd = NA)

## redefine parameters
n_nodes <- length(unique(sim$node))

## simulate centralities ####
## simulate age effect
sim_slope <- -0.3                     # set age effect -- smaller = bigger impact on invlogit scale as values large
sim_intcp <- 2
sim$mu <- sim$age * sim_slope + sim_intcp       # simulate mean centrality on normal scale
plot(sim$mu ~ sim$age)               # plot
sim$mu_std <- ( sim$mu - mean(sim$mu) ) / sd(sim$mu)

## simulate full distribution of samples per node
sim$sd <- 1#abs(sim_slope/3)           # make small to start with to be sure model should be able to detect difference
sim_dat <- matrix(data = NA, nrow = 4000, ncol = n_data, dimnames = list(NULL, sim$node_window))    # create matrix
for(j in 1:n_data){
  sim_dat[,j] <- rnorm(n = nrow(sim_dat), mean = sim$mu[j], sd = sim$sd[j])  # simulate distribution
}
plot(sim_dat[1,] ~ sim$age)          # plot simulated values against age

## standardise
sim_dat_std <- sim_dat               # create matrix to fill
for(i in 1:nrow(sim_dat_std)){
  sim_dat_std[i,] <- (sim_dat[i,] - mean(sim_dat[i,]) ) / sd(sim_dat[i,]) # standardise values
}
plot(sim_dat_std[1,] ~ sim$age)      # plot simulated values against age

## visualise
df_wide <- data.frame(sim_dat_std)
df_plot <- df_wide %>% 
  pivot_longer(cols = everything(),
               names_to = "node", values_to = "centrality") %>% 
  separate(node, into = c('X','node'), remove = T, sep = 1) %>% 
  dplyr::select(-X) %>% 
  mutate(node = as.integer(node)) %>% 
  left_join(sim[,c('node','age')], by = 'node') %>% 
  filter(node %in% seq(1, n_data, by = 2))
ggplot(data = df_plot, aes(x = centrality, fill = age)) +
  geom_density(linewidth = 0.4) +
  facet_grid(rows = vars(as.factor(node)), scales = "free") +
  labs(x = "Eigenvector centrality (standardised)") + 
  theme_void() + 
  theme(strip.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

## normal approximation ####
## normal approximation
plot(sim_dat_std[,1], sim_dat_std[,n_nodes])   # plot covariance (oldest and youngest to be sure it works for all pairs)
sim_cent_mu <- apply(sim_dat_std, 2, mean)     # calculate means per node
sim_cent_cov <- cov(sim_dat_std)               # calculate covariance matrix

## check normal approximation
sim_cent_samples <- MASS::mvrnorm(1e5, sim_cent_mu, sim_cent_cov)     # simulate from multivariate normal
plot(density(sim_dat_std[, 1]), lwd = 2, las = 1,                     # plot true density curve
     main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
lines(density(sim_cent_samples[, 1]), col = rgb(0,0,1,0.5), lwd = 2)  # overlay normal approximation

## run model -- age as an ordered categorical variable with 1 window ####
## categorise age
sim$age_cat <- ifelse(sim$age <= 15, 1,
                      ifelse(sim$age <= 20, 2,
                             ifelse(sim$age <= 25, 3,
                                    ifelse(sim$age <= 40, 4, 5))))

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

## prior predictive check
n <- 100
beta_age <- rnorm(n, 0, 1)
intercept  <- rnorm(n, 0, 0.8)
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
rm(n, beta_age, intercept, age_dirichlet, sigma, x, y, df_plot, df_wide) ; gc()

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
fit_sim$summary()

## extract posterior
params <- fit_sim$draws(format = 'draws_df')
delta <- params %>% 
  select(`delta[1]`,`delta[2]`,`delta[3]`,`delta[4]`,`delta[5]`)
delta_j <- params %>% 
  select(`delta_j[1]`,`delta_j[2]`,`delta_j[3]`,`delta_j[4]`,`delta_j[5]`,`delta_j[6]`)

## traceplot linear effect size
parameters_to_check <- c('intercept','beta_age','sigma','predictor[1]','predictor[20]','predictor[50]')
params %>% 
  select(all_of(parameters_to_check)) %>%
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'draw') %>% 
  mutate(chain_position = rep(rep(1:n_samples, each = length(parameters_to_check)), n_chains),
         chain = rep(1:n_chains, each = length(parameters_to_check)*n_samples)) %>% 
  ggplot(aes(x = chain_position, y = draw, colour = as.factor(chain)))+
  geom_line()+
  facet_wrap(. ~ parameter, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'none')
delta %>% 
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'draw') %>% 
  mutate(chain_position = rep(rep(1:n_samples, each = length(unique(parameter))), n_chains),
         chain = rep(1:n_chains, each = length(unique(parameter))*n_samples)) %>% 
  ggplot(aes(x = chain_position, y = draw, colour = as.factor(chain)))+
  geom_line()+
  facet_wrap(. ~ parameter, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'none')
delta_j %>% 
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'draw') %>% 
  mutate(chain_position = rep(rep(1:n_samples, each = length(unique(parameter))), n_chains),
         chain = rep(1:n_chains, each = length(unique(parameter))*n_samples)) %>% 
  ggplot(aes(x = chain_position, y = draw, colour = as.factor(chain)))+
  geom_line()+
  facet_wrap(. ~ parameter, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'none')

## posterior predictive check
plot(density(sim_dat_std[1, ]), las = 1, ylim = c(0,0.4),
     main = "Posterior predictive check (standardised centrality):\nblack = data, blue = predicted",
     col=rgb(0, 0, 0, 0.25))
for (i in 1:100) {
  j <- sample(1:length(params$beta_age), 1)
  lines(density(sim_dat_std[j, ]), col=rgb(0, 0, 0, 0.25))
  mu <- rep(NA, length(eigen_list$node_age))
  for(k in 1:length(eigen_list$node_age)){
    mu[k] <- params$intercept[j] + params$beta_age[j]*sum(delta_j[j,(1:eigen_list$node_age[k])])
  }
  sigma <- sim_cent_cov + diag(rep(params$sigma[j], n_data))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}

## predict from model
new_data <- data.frame(node = 1:length(unique(sim$age_cat)),
                       age_cat = rep(sort(unique(sim$age_cat))))
predict <- matrix(data = NA, nrow = length(params$beta_age), ncol = nrow(new_data),
                  dimnames = list(1:length(params$beta_age),
                                  1:length(unique(sim$age_cat))))
for(i in 1:nrow(predict)){
  for(j in 1:ncol(predict)){
    predict[i,j] <- params$intercept[i] + params$beta_age[i] * sum(delta_j[i,(1:j)])
  }
}
new_data$mean_predict <- apply(predict, 2, mean)
new_data$lwr_predict <- apply(predict, 2, quantile, probs = 0.025)
new_data$upr_predict <- apply(predict, 2, quantile, probs = 0.975)

## compare to standardised raw data
plot(new_data$mean_predict ~ new_data$age_cat, las = 1, type = 'l', ylim = c(-2,2))        # plot against age
points(sim$mu_std ~ sim$age_cat)                                                           # add raw points
polygon(y = c(new_data$lwr_predict, rev(new_data$upr_predict)), x = c(new_data$age_cat, rev(new_data$age_cat)),
        col = rgb(1,1,0,0.5), border = NA)

## unstandardise predictions
sim$mu
(new_data$predict_ustd <- new_data$mean_predict * sd(sim$mu) + mean(sim$mu))
sim <- sim %>% 
  left_join(new_data, by = 'age_cat')
plot(sim$predict_ustd ~ sim$mu)
abline(a = 0, b = 1)

## compare to unstandardised raw data
sim$lwr_ustd <- sim$lwr_predict * sd(sim$mu) + mean(sim$mu)
sim$upr_ustd <- sim$upr_predict * sd(sim$mu) + mean(sim$mu)
plot(sim$predict_ustd ~ sim$age_cat, las = 1, type = 'l')        # plot against age
points(sim$mu ~ sim$age_cat)                                     # add raw points
polygon(y = c(sim$lwr_ustd, rev(sim$upr_ustd)), x = c(sim$age_cat, rev(sim$age_cat)),
        col = rgb(1,1,0,0.5), border = NA)

## convert to invlogit scale
sim$mean_predict_invlogit <- invlogit(sim$predict_ustd)
sim$lwr_invlogit <- invlogit(sim$lwr_ustd)
sim$upr_invlogit <- invlogit(sim$upr_ustd)
sim$mu_invlogit <- invlogit(sim$mu)

## compare to invlogit raw data
plot(sim$mu_invlogit ~ sim$age_cat, ylim = c(0,1))
lines(sim$mean_predict_invlogit ~ sim$age_cat)
polygon(y = c(sim$lwr_invlogit, rev(sim$upr_invlogit)), x = c(sim$age_cat, rev(sim$age_cat)),
        col = rgb(1,1,0,0.5), border = NA)

## extract slope estimates -- none of these match the original input, but the predictions work??!
# original (true) slope value
sim_slope
# predicted (modelled) slope value = mean of differences bewteen categories, divided by the number of years each category represents
mean( c( (new_data$predict_ustd[2] - new_data$predict_ustd[1]) / (16-10), 
         (new_data$predict_ustd[3] - new_data$predict_ustd[2]) / (21-16), 
         (new_data$predict_ustd[4] - new_data$predict_ustd[3]) / (26-21),
         (new_data$predict_ustd[5] - new_data$predict_ustd[4]) / (40-26) )
)

## extract model fit
summary <- fit_sim$summary()
par(mfrow = c(3,1))
hist(summary$rhat, breaks = 50)
hist(summary$ess_bulk, breaks = 50)
hist(summary$ess_tail, breaks = 50)
par(mfrow = c(1,1))

