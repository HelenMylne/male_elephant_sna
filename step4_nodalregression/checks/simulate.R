library(LaplacesDemon)
library(tidyverse)
library(rstan)

## simulate population
min_age <- 11                                                # youngest individual
max_age <- 60                                                #  oldest  individual
n_nodes <- ((max_age+1) - min_age)*2                         # total nodes = 2 per age
sim <- data.frame(node = 1:n_nodes,                          # create data frame of individuals
                  age = rep(min_age:max_age, each = 2 ),
                  mu = NA, sd = NA)

## simulate age effect
sim_slope <- 0.8                     # set age effect -- smaller = bigger impact on invlogit scale as values large
sim_intcp <- -3
sim$mu <- sim$age * sim_slope + sim_intcp       # simulate mean centrality on normal scale
plot(sim$mu ~ sim$age)               # plot

## simulate full distribution of samples per node
sim$sd <- abs(sim_slope/3)           # make small to start with to be sure model should be able to detect difference
sim_dat <- matrix(data = NA, nrow = 4000, ncol = n_nodes)    # create matrix
for(j in 1:n_nodes){
  sim_dat[,j] <- rnorm(n = nrow(sim_dat), mean = sim$mu[j], sd = sim$sd[j])  # simulate distribution
}
plot(sim_dat[1,] ~ sim$age)          # plot simulated values against age

## convert to 0-1 bounded scale
sim_invlogit <- invlogit(sim_dat)    # convert to probability scale
plot(sim_invlogit[1,] ~ sim$age)     # plot simulated values against age

## standardise
sim_dat_std <- sim_dat               # create matrix to fill
for(i in 1:nrow(sim_dat_std)){
  sim_dat_std[i,] <- (sim_invlogit[i,] - mean(sim_invlogit[i,]) ) / sd(sim_invlogit[i,]) # standardise values
}
plot(sim_dat_std[1,] ~ sim$age)      # plot simulated values against age

## visualise
df_wide <- data.frame(sim_dat_std)
colnames(df_wide) <- 1:n_nodes
df_wide %>% 
  pivot_longer(cols = everything(),
               names_to = "node", values_to = "centrality") %>% 
  mutate(node = as.integer(node)) %>% 
  left_join(sim[,c('node','age')], by = 'node') %>% 
  filter(node %in% seq(2, 100, by = 2)) %>% 
  ggplot(aes(x = centrality, fill = age)) +
  geom_density(linewidth = 0.4) +
  facet_grid(rows = vars(as.factor(node)), scales = "free") +
  labs(x = "Eigenvector centrality (standardised)") + 
  theme_void() + 
  theme(strip.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

## normal approximation
plot(sim_dat_std[,1], sim_dat_std[,n_nodes])   # plot covariance (oldest and youngest to be sure it works for all pairs)
sim_cent_mu <- apply(sim_dat_std, 2, mean)     # calculate means per node
sim_cent_cov <- cov(sim_dat_std)               # calculate covariance matrix

## check normal approximation
sim_cent_samples <- MASS::mvrnorm(1e5, sim_cent_mu, sim_cent_cov)     # simulate from multivariate normal
plot(density(sim_dat_std[, 1]), lwd = 2, las = 1,                     # plot true density curve
     main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
lines(density(sim_cent_samples[, 1]), col = rgb(0,0,1,0.5), lwd = 2)  # overlay normal approximation

## create data
eigen_list <- list(num_nodes = n_nodes,
                   nodes = sim$node,
                   centrality_mu = sim_cent_mu,
                   centrality_cov = sim_cent_cov,
                   node_age = sim$age)

# load model
nodal_regression <- stan_model('models/eigen_regression_intercept.stan')

## run model
n_chains <- 4
fit_sim <- sampling(nodal_regression, data = eigen_list,
                    cores = n_chains, chains = n_chains)

## traceplot linear effect size
traceplot(fit_sim, pars = c('intercept','beta_age','sigma','predictor[1]','predictor[2]','predictor[3]'))

## posterior predictive check
params <- rstan::extract(fit_sim)
plot(density(sim_dat_std[1, ]), las = 1, ylim = c(0,1),
     main = "Posterior predictive check (standardised centrality):\nblack = data, blue = predicted",
     col=rgb(0, 0, 0, 0.25))
for (i in 1:100) {
  j <- sample(1:length(params$beta_age), 1)
  lines(density(sim_dat_std[j, ]), col=rgb(0, 0, 0, 0.25))
  mu <- params$beta_age[j]*eigen_list$node_age + params$intercept[j]
  sigma <- sim_cent_cov + diag(rep(params$sigma[j], n_nodes))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}

## extract slope parameter
plot(density(params$beta_age), las = 1,                # plot slope draws
     main = "Posterior slope draws", xlim = c(-0.2, 0.2))
lines(density(rnorm(1e5,0,0.1)), col = 'blue')         # prior slope as defined in model code = normal(0,0.1)
abline(v = sim_slope, lty = 2)                         # true slope
( age_slope <- round(quantile(params$beta_age, probs = c(0.5, 0.025, 0.975)), 2) )

## compare to raw data
sim$mean_predict <- invlogit(sim$age * mean(params$beta_age) + mean(params$intercept))           # mean prediction = age * mean slope
sim$lwr_predict <- invlogit(sim$age * quantile(params$beta_age, probs = 0.025) + quantile(params$intercept, probs = 0.025))
sim$upr_predict <- invlogit(sim$age * quantile(params$beta_age, probs = 0.975) + quantile(params$intercept, probs = 0.975))
plot(sim$mean_predict ~ sim$age, type = 'l', ylim = c(0,1), las = 1)    # plot against age
points(invlogit(sim$mu) ~ sim$age)                                 # add raw points
polygon(y = c(sim$lwr_predict, rev(sim$upr_predict)), x = c(sim$age, rev(sim$age)), col = rgb(1,1,0,0.5), border = NA)


