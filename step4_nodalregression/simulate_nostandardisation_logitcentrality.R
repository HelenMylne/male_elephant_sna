library(LaplacesDemon)
library(tidyverse)
library(cmdstanr)

## simulate population
min_age <- 11                                                # youngest individual
max_age <- 60                                                #  oldest  individual
n_nodes <- ((max_age+1) - min_age)*2                         # total nodes = 2 per age
sim <- data.frame(node = 1:n_nodes,                          # create data frame of individuals
                  age = rep(min_age:max_age, each = 2 ),
                  mu = NA, sd = NA)

## simulate age effect
sim_slope <- 0.2                     # set age effect -- smaller = bigger impact on invlogit scale as values large
sim_intcp <- 2
sim$mu <- sim$age * sim_slope + sim_intcp       # simulate mean centrality on normal scale
plot(sim$mu ~ sim$age)               # plot

## simulate full distribution of samples per node
sim$sd <- 1#abs(sim_slope/3)           # make small to start with to be sure model should be able to detect difference
sim_dat <- matrix(data = NA, nrow = 4000, ncol = n_nodes)    # create matrix
for(j in 1:n_nodes){
  sim_dat[,j] <- rnorm(n = nrow(sim_dat), mean = sim$mu[j], sd = sim$sd[j])  # simulate distribution
}
plot(sim_dat[1,] ~ sim$age)          # plot simulated values against age

## visualise
df_wide <- data.frame(sim_dat)
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
plot(sim_dat[,1], sim_dat[,n_nodes])   # plot covariance (oldest and youngest to be sure it works for all pairs)
sim_cent_mu <- apply(sim_dat, 2, mean)     # calculate means per node
sim_cent_cov <- cov(sim_dat)               # calculate covariance matrix

## check normal approximation
sim_cent_samples <- MASS::mvrnorm(1e5, sim_cent_mu, sim_cent_cov)     # simulate from multivariate normal
plot(density(sim_dat[, 1]), lwd = 2, las = 1,                     # plot true density curve
     main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
lines(density(sim_cent_samples[, 1]), col = rgb(0,0,1,0.5), lwd = 2)  # overlay normal approximation

## create data
eigen_list <- list(num_nodes = n_nodes,
                   nodes = sim$node,
                   centrality_mu = sim_cent_mu,
                   centrality_cov = sim_cent_cov,
                   node_age = sim$age)

## check inputs
plot(sim_cent_mu ~ sim$age)

## prior predictive check
n <- 100
beta_age <- rnorm(n, 0, 0.5)
intercept  <- rnorm(n, 0, 5)
plot(NULL, las = 1, xlab = 'age (years)', ylab = 'eigenvector (logit transformed)',
     ylim = c(min(sim_cent_mu)-1, max(sim_cent_mu)+1), xlim = c(min(sim$age), max(sim$age)))
abline(h = min(sim_cent_mu), lty = 2) ; abline(h = max(sim_cent_mu), lty = 2)
for(i in 1:n){
  lines(x = seq(min(sim$age), max(sim$age), length.out = 2),
        y = intercept[i] + beta_age[i]*seq(min(sim$age), max(sim$age), length.out = 2),
        col = rgb(0,0,1,0.4))
} # looks really bad, but I think that's because of the intercept -- the intercept should be centred on zero because I don't know if an elephant of age 0 should have an average centrality score greater or less than 0.5 (0 once logit transformed), but then this will start the lines way off centre depending on the particular values I set for sim_intcp. Generally looking at it going yes the intercept value covers most of the required range of parameter space, and the slopes coming out of it would be sufficient to cover all space depending on where intercept started them

## load model
#nodal_regression <- stan_model('models/eigen_regression_intercept.stan')
nodal_regression <- cmdstan_model('models/eigen_regression_intercept_invlogitscale.stan')

## run model
n_chains <- 4
n_samples <- 1000
#fit_sim <- sampling(nodal_regression, data = eigen_list, chains = n_chains, cores = n_chains)
fit_sim <- nodal_regression$sample(data = eigen_list,
                                   chains = n_chains, parallel_chains = n_chains,
                                   iter_warmup = n_samples, iter_sampling = n_samples)

## view summary
fit_sim$summary()

## extract posterior
#params <- rstan::extract(fit_sim)
params <- fit_sim$draws(format = 'draws_df')

## traceplot linear effect size
#traceplot(fit_sim, pars = c('intercept','beta_age','sigma','predictor[1]','predictor[50]','predictor[100]'))
params %>% 
  select(intercept,beta_age,sigma,`predictor[1]`,`predictor[50]`,`predictor[100]`) %>% 
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'draw') %>% 
  mutate(chain_position = rep(rep(1:n_samples, each = 6), n_chains),
         chain = rep(1:n_chains, each = 6*n_samples)) %>% 
  #filter(chain == 4) %>% 
  ggplot(aes(x = chain_position, y = draw, colour = as.factor(chain)))+
  geom_line()+
  facet_wrap(. ~ parameter, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'none') # chains 1-3 are good, chain 4 poorly mixed after about 600 draws

## posterior predictive check
plot(density(sim_dat[1, ]), las = 1, ylim = c(0,0.4),
     main = "Posterior predictive check (standardised centrality):\nblack = data, blue = predicted",
     col=rgb(0, 0, 0, 0.25))
for (i in 1:100) {
  j <- sample(1:length(params$beta_age), 1)
  lines(density(sim_dat[j, ]), col=rgb(0, 0, 0, 0.25))
  mu <- params$beta_age[j]*eigen_list$node_age + params$intercept[j]
  sigma <- sim_cent_cov + diag(rep(params$sigma[j], n_nodes))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}

## extract slope parameter
plot(density(params$beta_age), las = 1,                # plot slope draws
     main = "Posterior slope draws")
abline(v = sim_slope, lty = 2)                         # true slope
( age_slope <- round(quantile(params$beta_age, probs = c(0.5, 0.025, 0.975)), 2) )

## extract intercept parameter
plot(density(params$intercept), las = 1,                # plot slope draws
     main = "Posterior slope draws")
abline(v = sim_intcp, lty = 2)                         # true slope
( intercept <- round(quantile(params$intercept, probs = c(0.5, 0.025, 0.975)), 2) )

## compare to raw data
sim$mean_predict <- sim$age * mean(params$beta_age) + mean(params$intercept)           # mean prediction = age * mean slope
sim$lwr_predict <- sim$age * quantile(params$beta_age, probs = 0.025) + quantile(params$intercept, probs = 0.025)
sim$upr_predict <- sim$age * quantile(params$beta_age, probs = 0.975) + quantile(params$intercept, probs = 0.975)
plot(sim$mean_predict ~ sim$age, type = 'l', las = 1)        # plot against age
points(sim$mu ~ sim$age)                                 # add raw points -- seems to be very good at predicting when the value true value is not 0
polygon(y = c(sim$lwr_predict, rev(sim$upr_predict)), x = c(sim$age, rev(sim$age)),
        col = rgb(1,1,0,0.5), border = NA)

## compare to original values
paste0('true slope value = ', sim_slope, '; model slope value = ', round(mean(params$beta_age),2) )
paste0('true intercept value = ', sim_intcp, '; model intercept value = ', round(mean(params$intercept),2) )

## extract model fit
summary <- fit_sim$summary()
par(mfrow = c(3,1))
hist(summary$rhat, breaks = 50)
hist(summary$ess_bulk, breaks = 50)
hist(summary$ess_tail, breaks = 50)
par(mfrow = c(1,1))

## convert to invlogit scale
sim$mean_predict_invlogit <- invlogit(sim$mean_predict)
sim$lwr_invlogit <- invlogit(sim$lwr_predict)
sim$upr_invlogit <- invlogit(sim$upr_predict)
sim$mu_invlogit <- invlogit(sim$mu)
plot(sim$mu_invlogit ~ sim$age, ylim = c(0.95,1.01))
lines(sim$mean_predict_invlogit ~ sim$age)
polygon(y = c(sim$lwr_invlogit, rev(sim$upr_invlogit)), x = c(sim$age, rev(sim$age)),
        col = rgb(1,1,0,0.5), border = NA)
