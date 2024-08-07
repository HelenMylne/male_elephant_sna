#### set up ####
# script to simulate dyadic regression and check actually working
#library(tidyverse) ; library(LaplacesDemon) ; library(car) ; library(cmdstanr) ; library(bisonR) ; library(brms)
library(cmdstanr, lib.loc = '../packages/')  # library(cmdstanr)
library(tidyverse, lib.loc = '../packages/') # library(tidyverse)
library(car, lib.loc = '../packages/')       # library(car)
#library(bisonR, lib.loc = '../packages/')    # library(bisonR)
#library(brms, lib.loc = '../packages/')      # library(brms)
library(LaplacesDemon, lib.loc = '../packages/')

set_cmdstan_path('R:/rsrch/df525/phd/hkm513/packages/.cmdstan2/cmdstan-2.33.1/')
theme_set(theme_classic())

#### simulate data ####
## simulate population
min_age <- 10                               # youngest individual
max_age <- 60                               #  oldest  individual
n_nodes <- ((max_age - min_age)/5)+1        # total nodes = 2 per age
sim_nodes1 <- data.frame(node = 1:n_nodes,   # create data frame of individuals
                  age = seq(from = min_age, 
                            to = max_age, by = 5)) #%>% 
  #mutate(age_std = ( age - mean(age) ) / sd(age))  # should age be standardised across both min and max variables or just within variable? when I ran this without them done together, max could be lower than min, but when together the predictions are thrown off and the wrong coefficients returned

## convert to dyadic dataframe
sim_dyads1 <- data.frame(node = rep(sim_nodes1$node, each = n_nodes),
                        node_2 = rep(sim_nodes1$node, n_nodes)) %>% 
  filter(node != node_2) %>% filter(node < node_2) %>% 
  left_join(sim_nodes1, by = 'node') %>% 
  rename(node_1 = node,
         #age1_std = age_std,
         age_1 = age) %>% 
  rename(node = node_2) %>% left_join(sim_nodes1, by = 'node') %>% 
  rename(node_2 = node,
         #age2_std = age_std,
         age_2 = age) %>% 
  mutate(age_min = NA, age_max = NA, mu = NA, sd = NA) %>% 
  mutate(window = 1)
for(i in 1:nrow(sim_dyads1)){
  sim_dyads1$age_min[i] <- min(sim_dyads1$age_1[i], sim_dyads1$age_2[i])
  sim_dyads1$age_max[i] <- max(sim_dyads1$age_1[i], sim_dyads1$age_2[i])
  #sim_dyads1$age_min_std[i] <- min(sim_dyads1$age1_std[i], sim_dyads1$age2_std[i])
  #sim_dyads1$age_max_std[i] <- max(sim_dyads1$age1_std[i], sim_dyads1$age2_std[i])
}

## simulate population 2
sim_nodes2 <- data.frame(node = c(sample(1:n_nodes, 5, replace = F),
                                  (n_nodes+1):(n_nodes*2)),
                         age = NA)
for(i in 1:nrow(sim_nodes2)){
  if(sim_nodes2$node[i] %in% sim_nodes1$node) {
    sim_nodes2$age[i] <- sim_nodes1$age[sim_nodes1$node == sim_nodes2$node[i]] + 2
  }
}
sim_nodes2$age[which(is.na(sim_nodes2$age) == TRUE)] <- seq(min_age, max_age, length.out = length(which(is.na(sim_nodes2$age) == TRUE)))

## convert to dyadic dataframe
sim_dyads2 <- data.frame(node = rep(sim_nodes2$node, each = n_nodes),
                         node_2 = rep(sim_nodes2$node, n_nodes)) %>% 
  filter(node != node_2) %>% filter(node < node_2) %>% 
  left_join(sim_nodes2, by = 'node') %>% 
  rename(node_1 = node,
         age_1 = age) %>% 
  rename(node = node_2) %>% left_join(sim_nodes2, by = 'node') %>% 
  rename(node_2 = node,
         age_2 = age) %>% 
  mutate(age_min = NA, age_max = NA, mu = NA, sd = NA) %>% 
  mutate(window = 2)
for(i in 1:nrow(sim_dyads2)){
  sim_dyads2$age_min[i] <- min(sim_dyads2$age_1[i], sim_dyads2$age_2[i])
  sim_dyads2$age_max[i] <- max(sim_dyads2$age_1[i], sim_dyads2$age_2[i])
}

## combine
sim_dyads <- rbind(sim_dyads1, sim_dyads2)
sim_dyads <- sim_dyads %>% 
  mutate(dyad_all = paste0(node_1, '_', node_2),
         dyad_window = paste0(node_1, '_', node_2, '_', window)) %>% 
  mutate(dyad_id_all = as.integer(as.factor(dyad_all)),
         dyad_id_window = as.integer(as.factor(dyad_window)))
sim_nodes <- rbind(sim_nodes1, sim_nodes2)
n_dyads <- length(unique(sim_dyads$dyad_id_all))
n_data <- nrow(sim_dyads)
n_nodes <- length(unique(sim_nodes$node))

## standardise ages -- CHECK THIS SHOULD ACTUALLY BE DONE SEPARATELY AND NOT ALTOGETHER
ages <- data.frame(age_min = unique(c(sim_dyads$age_1, sim_dyads$age_2)),
                   age_max = unique(c(sim_dyads$age_1, sim_dyads$age_2))) %>% 
  mutate(age_min_std = (age_min - mean(age_min) ) / sd(age_min),
         age_max_std = (age_max - mean(age_max) ) / sd(age_max))
sim_dyads <- sim_dyads %>% 
  left_join(ages[,c('age_min','age_min_std')], by = 'age_min') %>% 
  left_join(ages[,c('age_max','age_max_std')], by = 'age_max')

## simulate age effect
sim_min <- 1
sim_max <- 0.5
sim_int <- 10
sim_dyads$mu <- sim_dyads$age_min * sim_min + sim_dyads$age_max * sim_max   # simulate mean centrality on normal scale
plot(sim_dyads$mu ~ sim_dyads$age_min)               # plot
plot(sim_dyads$mu ~ sim_dyads$age_max)               # plot

## standardise edges
sim_dyads$mu_std <- ( sim_dyads$mu - mean(sim_dyads$mu) ) / sd(sim_dyads$mu)

## simulate full distribution of samples per node
sim_dyads$sd <- abs(sim_max/3)           # make small to start with to be sure model should be able to detect difference
sim_dat <- matrix(data = NA, nrow = 1000, ncol = n_data, dimnames = list(1:1000, sim_dyads$dyad_id_window))    # create matrix
for(j in 1:n_data){
  sim_dat[,j] <- rnorm(n = nrow(sim_dat), mean = sim_dyads$mu_std[j], sd = sim_dyads$sd[j])  # simulate distribution
}
plot(sim_dat[1,] ~ sim_dyads$age_min)          # plot simulated values against minimum age
plot(sim_dat[1,] ~ sim_dyads$age_max)          # plot simulated values against maximum age

#### plot raw data ####
ggplot()+
  geom_point(data = sim_dyads, aes(x = age_min_std, y = mu_std, colour = age_max_std), shape = 19)+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')

sim_dat_long <- sim_dat %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'dyad_id_window', values_to = 'edge_draw') %>% 
  mutate(dyad_id_window = as.numeric(dyad_id_window)) %>% 
  left_join(sim_dyads, by = 'dyad_id_window')
ggplot()+
  geom_point(data = sim_dat_long, aes(x = age_min_std, y = edge_draw, colour = as.factor(window)),
             shape = 19, alpha = 0.05)+
  geom_point(data = sim_dyads[sim_dyads$window == 1,], aes(x = age_min_std, y = mu_std),
             shape = 19, colour = 'black', size = 1)+
  geom_point(data = sim_dyads[sim_dyads$window == 2,], aes(x = age_min_std, y = mu_std),
             shape = 19, colour = 'white', size = 0.5)+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')

## heat map: x = min age, y = max age, colour = edge weight
ggplot()+
  geom_tile(data = sim_dyads, mapping = aes(x = age_min_std, y = age_max_std, fill = mu_std))+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('age of older dyad member')+
  labs(fill = 'logit(mean edge weight)')

#### fit multivariate Gaussian distribution to output of edge weight model ####
### fit a multivariate normal dist to the edges -- quantities will be given to Stan model as data to model joint posteriors of edge weight in the regression
logit_edge_draws_mu <- apply(sim_dat, 2, mean)
logit_edge_draws_cov <- cov(sim_dat)

#### plot to see how well the approximation is working -- pretty unnecessary when using simulated values! ####
### Randomly selecting samples to examine
num_check <- 20
selected_samples <- sample(1:(n_dyads-1), num_check, replace = FALSE)

### Setting grid layout
rows <- floor(sqrt(num_check))
cols <- ceiling(num_check / rows)
par(mfrow=c(rows, cols), mar=c(2,2,2,1))

### plot
for (i in selected_samples) {
  mu <- logit_edge_draws_mu[i]
  sd <- sqrt(logit_edge_draws_cov[i,i])
  
  fitted_values <- rnorm(1e5, mean = mu, sd = sd)
  
  hist(unlist(sim_dat[,i]), probability = TRUE, las = 1,
       main = paste("Dyad", i), xlab = "Value", breaks = 50)
  lines(density(fitted_values), col="blue", lw=1.5)
}

for (i in selected_samples) {
  mu <- logit_edge_draws_mu[i]
  sd <- sqrt(logit_edge_draws_cov[i,i])
  
  fitted_values <- rnorm(1e5, mean = mu, sd = sd)
  
  plot(unlist(sim_dat[,i]), unlist(sim_dat[,i+1]), col = rgb(0,0,1,0.05), las = 1,
       main = paste("cov ", i ,"&",i+1))
}

### reset plot window and clean up
par(mfrow=c(1,1), mai = c(1,1,0.5,0.5))
rm(ages, sim_dyads1, sim_dyads2, sim_nodes1, sim_nodes2, cols, fitted_values, i, j, max_age, min_age, mu, num_check, rows, sd, selected_samples) ; gc()

######### Standardised data without intercept -- single window only ########
#### prior predictive check ####
n <- 100
beta_age_min <- rnorm(n, 0, 1)
beta_age_max <- rnorm(n, 0, 1)
plot(NULL, las = 1, xlab = 'minimum age (standardised)', ylab = 'logit edge weight (standardised)',
     ylim = c(min(logit_edge_draws_mu)-5, max(logit_edge_draws_mu)+5), xlim = c(min(sim_dyads$age_min_std), max(sim_dyads$age_max_std)))
abline(h = min(logit_edge_draws_mu), lty = 2) ; abline(h = max(logit_edge_draws_mu), lty = 2)
for(i in 1:n){
  lines(x = seq(min(sim_dyads$age_min_std), max(sim_dyads$age_min_std), length.out = 2),
        y = beta_age_min[i]*seq(min(sim_dyads$age_min_std), max(sim_dyads$age_min_std), length.out = 2) + beta_age_max[i]*seq(min(sim_dyads$age_max_std), max(sim_dyads$age_max_std), length.out = 2),
        col = rgb(0,0,1,0.4))
} # looks mad when I set true slope very small, but decent when larger -- given that I have no idea what the real version should be, I've kept them relatively weak

#### fit dyadic regression ####
## create data list
dyad_data <- list(
  num_dyads = n_dyads,                      # number of dyads
  num_nodes = n_nodes,                      # number of nodes
  logit_edge_mu = logit_edge_draws_mu,      # sample means of the logit edge weights
  logit_edge_cov = logit_edge_draws_cov,    # sample covariance of logit edge weights
  age_min = sim_dyads$age_min_std,              # age of younger dyad member
  age_max = sim_dyads$age_max_std,              # age of  older  dyad member
  node_1 = sim_dyads$node_1,                # node IDs for multimembership effects
  node_2 = sim_dyads$node_2                 # node IDs for multimembership effects
)

## load dyadic regression model
dyadic_regression <- cmdstan_model('models/dyadic_regression.stan')
n_chains <- 4

## fit dyadic regression
fit_dyadreg_sim <- dyadic_regression$sample(
  data = dyad_data,
  chains = n_chains,
  parallel_chains = n_chains)

#### check outputs ####
# obtain summary
fit_dyadreg_sim$summary() %>% 
  filter(variable %in% c('beta_age_max','beta_age_min'))

## extract draws
draws <- fit_dyadreg_sim$draws(format = 'df')

## extract dyadic regression slopes
b_max <- draws$beta_age_max
b_min <- draws$beta_age_min
sigma <- draws$sigma
parameters <- data.frame(beta_age_max = b_max,
                         beta_age_min = b_min,
                         sigma = sigma) %>% 
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>% 
  pivot_longer(cols = c('beta_age_max','beta_age_min','sigma'), names_to = 'parameter', values_to = 'slope_draw')

## traceplots
ggplot(data = parameters)+
  geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')

#### plot predictions ####
## posterior predictive check
plot(density(as.numeric(sim_dat[1, ])), main = "Posterior predictive density of edge weights:\nblack = measured edge, red = predicted edge",
     ylim = c(0, 1), col = rgb(0, 0, 0, 0.25), las = 1)
for (i in 1:100) {
  j <- sample(1:1000, 1)
  mu_plot <- b_min[j]*dyad_data$age_min + b_max[j]*dyad_data$age_max
  sigma_plot <- dyad_data$logit_edge_cov + diag(rep(sigma[j], n_dyads))
  mv_norm <- MASS::mvrnorm(1, mu_plot, sigma_plot)
  
  lines(density(as.numeric(sim_dat[j, ])), col = rgb(0, 0, 0, 0.25)) # black lines for edge samples
  lines(density(mv_norm), col = rgb(1, 0, 0, 0.25))                  # red lines for predictions

}

## obtain median and 95% CI
plot(density(b_max)) ; abline(v = 0, lty = 2)
plot(density(b_min)) ; abline(v = 0, lty = 2)
plot(density(sigma)) ; abline(v = 0, lty = 2)
(b_max_summary <- round(quantile(b_max, probs=c(0.5, 0.025, 0.975)), 2)) ; sim_max
(b_min_summary <- round(quantile(b_min, probs=c(0.5, 0.025, 0.975)), 2)) ; sim_min

## convert to unstandardised scale
b_max_unstd <- b_max * (sd(sim_dyads$mu) / sd(sim_dyads$age_max))
b_min_unstd <- b_min * (sd(sim_dyads$mu) / sd(sim_dyads$age_min))
plot(density(b_max_unstd)) ; abline(v = sim_max, lty = 2)
plot(density(b_min_unstd)) ; abline(v = sim_min, lty = 2)
round(quantile(b_max_unstd, probs=c(0.5, 0.025, 0.975)), 2) ; sim_max
round(quantile(b_min_unstd, probs=c(0.5, 0.025, 0.975)), 2) ; sim_min

## calculate predictions
pred <- data.frame(age_min = rep(seq(from = min(sim_dyads$age_min_std), to = max(sim_dyads$age_min_std), length.out = 10),
                                     each = 10),
                   age_max = rep(seq(from = min(sim_dyads$age_max_std), to = max(sim_dyads$age_max_std), length.out = 10),
                                 10),
                   b_min_mid = b_min_summary[1],
                   b_min_lwr = b_min_summary[2],
                   b_min_upr = b_min_summary[3],
                   b_max_mid = b_max_summary[1],
                   b_max_lwr = b_max_summary[2],
                   b_max_upr = b_max_summary[3]) %>%
  mutate(mu_mid = b_min_mid*age_min + b_max_mid*age_max,
         mu_lwr = b_min_lwr*age_min + b_max_lwr*age_max,
         mu_upr = b_min_upr*age_min + b_max_upr*age_max) %>%
  mutate(invlogit_mu_mid = LaplacesDemon::invlogit(mu_mid),
         invlogit_mu_lwr = LaplacesDemon::invlogit(mu_lwr),
         invlogit_mu_upr = LaplacesDemon::invlogit(mu_upr))

## convert to unstandardised scale
pred$mu_mid_unstd <- pred$mu_mid*sd(sim_dyads$mu) + mean(sim_dyads$mu)
pred$mu_lwr_unstd <- pred$mu_lwr*sd(sim_dyads$mu) + mean(sim_dyads$mu)
pred$mu_upr_unstd <- pred$mu_upr*sd(sim_dyads$mu) + mean(sim_dyads$mu)

## plot
ggplot()+
  geom_ribbon(data = pred,#[pred$age_max %in% seq(from = min(sim_dyads$age_min_std), to = max(sim_dyads$age_min_std), length.out = 10)[c(1,3,5,8,10)],],
              aes(x = age_min, group = age_max, fill = age_max,
                  #ymin = invlogit_mu_lwr, ymax = invlogit_mu_upr),
                  ymin = mu_lwr, ymax = mu_upr),
              alpha = 0.3)+
  geom_line(data = pred,#[pred$age_max %in% seq(from = min(sim_dyads$age_min_std), to = max(sim_dyads$age_min_std), length.out = 5),],
            aes(x = age_min, y = mu_mid, #y = invlogit_mu_mid,
                colour = age_max, group = age_max),
            linewidth = 1)+
  geom_point(data = sim_dyads, aes(x = age_min_std, y = mu_std, colour = age_max_std),
             alpha = 0.5)+
  scale_colour_viridis_c()+ scale_fill_viridis_c()+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+
  labs(colour = 'maximum age', fill = 'maximum age')

pred$age_min_unstd <- pred$age_min * sd(sim_dyads$age_min) + mean(sim_dyads$age_min)
pred$age_max_unstd <- pred$age_max * sd(sim_dyads$age_max) + mean(sim_dyads$age_max)
ggplot()+
  geom_ribbon(data = pred,#[pred$age_max %in% seq(from = min(sim_dyads$age_min_std), to = max(sim_dyads$age_min_std), length.out = 10)[c(1,3,5,8,10)],],
              aes(x = age_min_unstd, group = age_max_unstd, fill = age_max_unstd,
                  #ymin = invlogit_mu_lwr, ymax = invlogit_mu_upr),
                  ymin = mu_lwr_unstd, ymax = mu_upr_unstd),
              alpha = 0.3)+
  geom_line(data = pred,#[pred$age_max %in% seq(from = min(sim_dyads$age_min_std), to = max(sim_dyads$age_min_std), length.out = 5),],
            aes(x = age_min_unstd, y = mu_mid_unstd, #y = invlogit_mu_mid,
                colour = age_max_unstd, group = age_max_unstd),
            linewidth = 1)+
  geom_point(data = sim_dyads, aes(x = age_min, y = mu, colour = age_max),
             alpha = 0.5)+
  scale_colour_viridis_c()+ scale_fill_viridis_c()+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+
  labs(colour = 'maximum age', fill = 'maximum age')

ggplot()+
  geom_ribbon(data = pred,#[pred$age_max %in% seq(from = min(sim_dyads$age_min_std), to = max(sim_dyads$age_min_std), length.out = 10)[c(1,3,5,8,10)],],
              aes(x = age_min_unstd, group = age_max_unstd, fill = age_max_unstd,
                  #ymin = invlogit_mu_lwr, ymax = invlogit_mu_upr),
                  ymin = invlogit_mu_lwr, ymax = invlogit_mu_upr),
              alpha = 0.3)+
  geom_line(data = pred,#[pred$age_max %in% seq(from = min(sim_dyads$age_min_std), to = max(sim_dyads$age_min_std), length.out = 5),],
            aes(x = age_min_unstd, y = invlogit_mu_mid, #y = invlogit_mu_mid,
                colour = age_max_unstd, group = age_max_unstd),
            linewidth = 1)+
  geom_point(data = sim_dyads, aes(x = age_min, y = LaplacesDemon::invlogit(mu_std), colour = age_max),
             alpha = 0.5)+
  scale_colour_viridis_c()+ scale_fill_viridis_c()+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+
  labs(colour = 'maximum age', fill = 'maximum age')

#### clean up ####
rm(list = ls()[! ls() %in% c('sim_dat','sim_dyads','sim_nodes','max_age','min_age','n_dyads','n_nodes','sim_max','sim_min','sim_int')]) ; gc()

######### Standardised data with intercept -- combine two time windows ########
#### prior predictive check ####
n <- 100
beta_age_min <- rnorm(n, 0, 1)
beta_age_max <- rnorm(n, 0, 1)
intercept <- rnorm(n, 0, 1)
plot(NULL, las = 1, xlab = 'minimum age (standardised)', ylab = 'logit edge weight (standardised)',
     ylim = c(min(logit_edge_draws_mu)-5, max(logit_edge_draws_mu)+5), xlim = c(min(sim_dyads$age_min_std), max(sim_dyads$age_max_std)))
abline(h = min(logit_edge_draws_mu), lty = 2) ; abline(h = max(logit_edge_draws_mu), lty = 2)
for(i in 1:n){
  lines(x = seq(min(sim_dyads$age_min_std), max(sim_dyads$age_min_std), length.out = 2),
        y = intercept[i] + beta_age_min[i]*seq(min(sim_dyads$age_min_std), max(sim_dyads$age_min_std), length.out = 2) + beta_age_max[i]*seq(min(sim_dyads$age_max_std), max(sim_dyads$age_max_std), length.out = 2),
        col = rgb(0,0,1,0.4))
} # looks mad when I set true slope very small, but decent when larger -- given that I have no idea what the real version should be, I've kept them relatively weak
rm(n, beta_age_max, beta_age_min, intercept, i) ; gc()

#### fit dyadic regression ####
## create data list
n_windows <- length(unique(sim_dyads$window))
dyad_data <- list(
  num_data = n_data,                        # number of edges in total data
  num_dyads = n_dyads,                      # number of dyads
  num_nodes = n_nodes,                      # number of nodes
  num_windows = n_windows,                  # Number of time windows
  logit_edge_mu = logit_edge_draws_mu,      # sample means of the logit edge weights
  logit_edge_cov = logit_edge_draws_cov,    # sample covariance of logit edge weights
  age_min = sim_dyads$age_min_std,          # age of younger dyad member
  age_max = sim_dyads$age_max_std,          # age of  older  dyad member
  node_1 = sim_dyads$node_1,                # node IDs for multimembership effects
  node_2 = sim_dyads$node_2,                # node IDs for multimembership effects
  window = sim_dyads$window,                # ID of time window (random effect)
  dyad_id = sim_dyads$dyad_id_all           # ID of dyad (random effect, can have the same value where same dyad seen multiple times)
)

## load dyadic regression model
dyadic_regression <- cmdstan_model('models/dyadic_regression_combinewindows.stan')
n_chains <- 4

## fit dyadic regression
fit_dyadreg_sim <- dyadic_regression$sample(
  data = dyad_data,
  chains = n_chains,
  parallel_chains = n_chains)

#### check outputs ####
# obtain summary
fit_dyadreg_sim$summary() %>% 
  filter(variable %in% c('beta_age_max','beta_age_min','intercept'))

## extract draws
draws <- fit_dyadreg_sim$draws(format = 'df')

## extract dyadic regression slopes
b_max <- draws$beta_age_max
b_min <- draws$beta_age_min
intercept <- draws$intercept
sigma <- draws$sigma
rand_window <- draws[,c('rand_window[1]','rand_window[2]')]
colnames(rand_window) <- 1:2
parameters <- data.frame(beta_age_max = b_max,
                         beta_age_min = b_min,
                         intercept = intercept,
                         sigma = sigma) %>% 
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>% 
  pivot_longer(cols = c('beta_age_max','beta_age_min','sigma','intercept'),
               names_to = 'parameter', values_to = 'slope_draw')

## traceplots -- quite wandery, but well mixed with one another -- could just be because there isn't much data?
ggplot(data = parameters)+
  geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')

#### plot predictions ####
## posterior predictive check
plot(density(as.numeric(sim_dat[1, ])), main = "Posterior predictive density of edge weights:\nblack = measured edge, red = predicted edge",
     ylim = c(0, 1), col = rgb(0, 0, 0, 0.25), las = 1)
for (i in 1:100) {
  j <- sample(1:1000, 1)
  mu_plot <- b_min[j]*dyad_data$age_min + b_max[j]*dyad_data$age_max + intercept[j]
  sigma_plot <- dyad_data$logit_edge_cov + diag(rep(sigma[j], n_data))
  mv_norm <- MASS::mvrnorm(1, mu_plot, sigma_plot)
  
  lines(density(as.numeric(sim_dat[j, ])), col = rgb(0, 0, 0, 0.25)) # black lines for edge samples
  lines(density(mv_norm), col = rgb(1, 0, 0, 0.25))                  # red lines for predictions
  
}

## obtain median and 95% CI
plot(density(b_max)) ; abline(v = 0, lty = 2)
plot(density(b_min)) ; abline(v = 0, lty = 2)
plot(density(intercept)) ; abline(v = 0, lty = 2)
(b_max_summary <- round(quantile(b_max, probs=c(0.5, 0.025, 0.975)), 2)) ; sim_max
(b_min_summary <- round(quantile(b_min, probs=c(0.5, 0.025, 0.975)), 2)) ; sim_min
(intcp_summary <- round(quantile(intercept, probs=c(0.5, 0.025, 0.975)), 2)) ; sim_int

## convert to unstandardised scale
b_max_unstd <- b_max * (sd(sim_dyads$mu) / sd(sim_dyads$age_max))
b_min_unstd <- b_min * (sd(sim_dyads$mu) / sd(sim_dyads$age_min))
intcp_unstd <- mean(sim_dyads$mu) - intercept * mean(sim_dyads$age_min) # this is missing age_max_std
plot(density(b_max_unstd)) ; abline(v = sim_max, lty = 2)
plot(density(b_min_unstd)) ; abline(v = sim_min, lty = 2)
plot(density(intcp_unstd)) ; abline(v = sim_int, lty = 2)
round(quantile(b_max_unstd, probs=c(0.5, 0.025, 0.975)), 2) ; sim_max
round(quantile(b_min_unstd, probs=c(0.5, 0.025, 0.975)), 2) ; sim_min
round(quantile(intcp_unstd, probs=c(0.5, 0.025, 0.975)), 2) ; sim_int

## calculate predictions
pred <- data.frame(age_min = rep(seq(from = min(sim_dyads$age_min_std), to = max(sim_dyads$age_min_std), length.out = 10),
                                     each = 10*2*5),
                   age_max = rep(rep(seq(from = min(sim_dyads$age_max_std), to = max(sim_dyads$age_max_std), length.out = 10),
                                     2*5),
                                 each = 10),
                   window = rep(rep(c(1,2), 10*10*2), each = 5),
                   dyad_id = rep(sample(dyad_data$dyad_id, 5, replace = F),10*10*2)) %>%
  mutate(mu_mid = NA, mu_lwr = NA, mu_upr = NA)
for(i in 1:nrow(pred)){
  window_summary <- round(quantile(unlist(rand_window[,pred$window[i]]), probs=c(0.5, 0.025, 0.975)), 2)
  rand_dyad <- draws[,paste0('rand_dyad[',pred$dyad_id[i],']')]
  dyadid_summary <- round(quantile(unlist(rand_dyad), probs=c(0.5, 0.025, 0.975)), 2)
  pred$mu_mid[i] <-  b_min_summary[1]*pred$age_min[i] + b_max_summary[1]*pred$age_max[i] + intcp_summary[1] + window_summary[1] + dyadid_summary[1]
  pred$mu_lwr[i] <-  b_min_summary[2]*pred$age_min[i] + b_max_summary[2]*pred$age_max[i] + intcp_summary[2] + window_summary[2] + dyadid_summary[2]
  pred$mu_upr[i] <-  b_min_summary[3]*pred$age_min[i] + b_max_summary[3]*pred$age_max[i] + intcp_summary[3] + window_summary[3] + dyadid_summary[3]
}
pred <- pred %>% 
  mutate(mu_mid = unlist(mu_mid),
         mu_lwr = unlist(mu_lwr),
         mu_upr = unlist(mu_upr)) %>% 
  mutate(invlogit_mu_mid = invlogit(mu_mid),
         invlogit_mu_lwr = invlogit(mu_lwr),
         invlogit_mu_upr = invlogit(mu_upr))

## convert to unstandardised scale
pred$mu_mid_unstd <- pred$mu_mid*sd(sim_dyads$mu) + mean(sim_dyads$mu)
pred$mu_lwr_unstd <- pred$mu_lwr*sd(sim_dyads$mu) + mean(sim_dyads$mu)
pred$mu_upr_unstd <- pred$mu_upr*sd(sim_dyads$mu) + mean(sim_dyads$mu)

## plot
ggplot()+
  geom_ribbon(data = pred[pred$age_max %in% seq(from = min(sim_dyads$age_max_std), to = max(sim_dyads$age_max_std), length.out = 10)[c(1,5,10)],],
              aes(x = age_min, group = as.factor(round(age_max,2)), fill = as.factor(round(age_max, 2)),
                  #ymin = invlogit_mu_lwr, ymax = invlogit_mu_upr),
                  ymin = mu_lwr, ymax = mu_upr),
              alpha = 0.3)+
  geom_line(data = pred[pred$age_max %in% seq(from = min(sim_dyads$age_max_std), to = max(sim_dyads$age_max_std), length.out = 10)[c(1,5,10)],],
            aes(x = age_min, y = mu_mid, #y = invlogit_mu_mid,
                colour = as.factor(round(age_max,2)), group = as.factor(round(age_max,2))),
            linewidth = 1)+
  scale_colour_viridis_d()+ scale_fill_viridis_d()+
  # geom_point(data = sim_dyads, aes(x = age_min_std, y = mu_std, colour = age_max_std),
  #            alpha = 0.5)+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+
  labs(colour = 'maximum age', fill = 'maximum age')

pred$age_min_unstd <- pred$age_min * sd(sim_dyads$age_min) + mean(sim_dyads$age_min)
pred$age_max_unstd <- pred$age_max * sd(sim_dyads$age_max) + mean(sim_dyads$age_max)
ggplot(data = pred[pred$age_max %in% seq(from = min(sim_dyads$age_max_std), to = max(sim_dyads$age_max_std), length.out = 10)[c(1,5,10)],],)+
  geom_ribbon(aes(x = age_min_unstd, group = age_max_unstd, fill = age_max_unstd,
                  #ymin = invlogit_mu_lwr, ymax = invlogit_mu_upr),
                  ymin = mu_lwr_unstd, ymax = mu_upr_unstd),
              alpha = 0.3)+
  geom_line(aes(x = age_min_unstd, y = mu_mid_unstd, #y = invlogit_mu_mid,
                colour = age_max_unstd, group = age_max_unstd),
            linewidth = 1)+
  # geom_point(data = sim_dyads, aes(x = age_min, y = mu, colour = age_max),
  #            alpha = 0.5)+
  scale_colour_viridis_c()+ scale_fill_viridis_c()+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+
  labs(colour = 'maximum age', fill = 'maximum age')

ggplot()+
  geom_ribbon(data = pred,#[pred$age_max %in% seq(from = min(sim_dyads$age_min_std), to = max(sim_dyads$age_min_std), length.out = 10)[c(1,3,5,8,10)],],
              aes(x = age_min_unstd, group = age_max_unstd, fill = age_max_unstd,
                  #ymin = invlogit_mu_lwr, ymax = invlogit_mu_upr),
                  ymin = invlogit_mu_lwr, ymax = invlogit_mu_upr),
              alpha = 0.3)+
  geom_line(data = pred,#[pred$age_max %in% seq(from = min(sim_dyads$age_min_std), to = max(sim_dyads$age_min_std), length.out = 5),],
            aes(x = age_min_unstd, y = invlogit_mu_mid, #y = invlogit_mu_mid,
                colour = age_max_unstd, group = age_max_unstd),
            linewidth = 1)+
  geom_point(data = sim_dyads, aes(x = age_min, y = LaplacesDemon::invlogit(mu_std), colour = age_max),
             alpha = 0.5)+
  scale_colour_viridis_c()+ scale_fill_viridis_c()+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+
  labs(colour = 'maximum age', fill = 'maximum age')

unique(sim_dyads$age_min_std) ; unique(pred$age_min)
unique(sim_dyads$age_max_std) ; unique(pred$age_max)
sim_dyads$age_pair <- paste0(sim_dyads$age_min_std, '_', sim_dyads$age_max_std)
compare <- pred %>% 
  rename(dyad_id_all = dyad_id) %>% 
  mutate(age_pair = paste0(age_min, '_', age_max)) %>% 
  left_join(sim_dyads[,c('mu','mu_std','age_pair')], by = c('age_pair'), multiple = 'all') %>% 
  mutate(invlogit_mu_raw = LaplacesDemon::invlogit(mu)) %>%
  rename(mu_raw = mu, mu_raw_std = mu_std) %>% 
  select(age_min, age_max, mu_raw, mu_mid_unstd, mu_raw_std, mu_mid, invlogit_mu_raw, invlogit_mu_mid, mu_lwr, mu_lwr_unstd, invlogit_mu_lwr, mu_upr, mu_upr_unstd, invlogit_mu_upr)
ggplot(compare)+
  # geom_ribbon(data = compare,
  #             aes(x = mu_raw, group = mu_mid_unstd,
  #                 ymin = mu_lwr_unstd, ymax = mu_upr_unstd))+
  # geom_line(aes(x = mu_raw, y = mu_mid_unstd),
  #           linewidth = 1)+
  geom_point(aes(x = mu_raw, y = mu_mid_unstd),
             alpha = 0.1)+
  geom_line(data = data.frame(x = c(15,85), y = c(15,85)),
            aes(x = x, y = y),
            linewidth = 0.5)+
  scale_x_continuous('raw mean')+
  scale_y_continuous('predicted mean')+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

#### clean up ####
rm(list = ls()[! ls() %in% c('sim_dat','sim_dyads','sim_nodes','max_age','min_age','n_dyads','n_nodes','sim_max','sim_min','sim_int')]) ; gc()

######### Standardised data with intercept -- single window, ordered categorical predictor with 1 window: run all above as is but change definition of sim_dyads1 to sim_dyads and sim_nodes1 to sim_nodes then don't rbind to sim2 ########
sim_nodes$age_cat <- ifelse(sim_nodes$age <= 15, 1,
                            ifelse(sim_nodes$age <= 20, 2,
                                   ifelse(sim_nodes$age <= 25, 3,
                                          ifelse(sim_nodes$age <= 40, 4, 5))))
sim_dyads$age_max_cat <- ifelse(sim_dyads$age_max <= 15, 1,
                                ifelse(sim_dyads$age_max <= 20, 2,
                                       ifelse(sim_dyads$age_max <= 25, 3,
                                              ifelse(sim_dyads$age_max <= 40, 4, 5))))
sim_dyads$age_min_cat <- ifelse(sim_dyads$age_min <= 15, 1,
                                ifelse(sim_dyads$age_min <= 20, 2,
                                       ifelse(sim_dyads$age_min <= 25, 3,
                                              ifelse(sim_dyads$age_min <= 40, 4, 5))))

#### prior predictive check ####
n <- 100
beta_age_min <- rnorm(n, 0, 1.5)
beta_age_max <- rnorm(n, 0, 1.5)
intercept <- rnorm(n, 0, 1)
min_dirichlet <- rdirichlet(n, c(1,1,1,1,1))
max_dirichlet <- rdirichlet(n, c(1,1,1,1,1))
plot(NULL, las = 1, xlab = 'minimum age (standardised)', ylab = 'logit edge weight (standardised)',
     ylim = c(min(logit_edge_draws_mu)-5, max(logit_edge_draws_mu)+5), xlim = c(min(sim_dyads$age_min_cat), max(sim_dyads$age_max_cat)))
abline(h = min(logit_edge_draws_mu), lty = 2) ; abline(h = max(logit_edge_draws_mu), lty = 2)
x <- min(sim_nodes$age_cat):max(sim_nodes$age_cat)
for(i in 1:n){
  y <- rep(NA, length(x))
  for(j in 1:length(x)){
    y[j] <- intercept[i] + beta_age_min[i]*sum(min_dirichlet[i,][1:x[j]]) + beta_age_max[i]*sum(max_dirichlet[i,][1:x[j]])
  }
  lines(x = x, y = y, col = rgb(0,0,1,0.4))
}
rm(n, beta_age_max, beta_age_min, intercept, min_dirichlet, max_dirichlet, i, y, j, x) ; gc()

#### fit dyadic regression ####
## create data list
#n_windows <- length(unique(sim_dyads$window))
n_age_cat <- length(unique(sim_dyads$age_max_cat))
dyad_data <- list(
  #num_data = n_data,                        # number of edges in total data
  num_dyads = n_dyads,                      # number of dyads
  num_nodes = n_nodes,                      # number of nodes
  num_age_cat = n_age_cat,                  # number of unique age categories
  length_dirichlet = n_age_cat + 1,         # number of unique age categories + 1
  #num_windows = n_windows,                  # Number of time windows
  logit_edge_mu = logit_edge_draws_mu,      # sample means of the logit edge weights
  logit_edge_cov = logit_edge_draws_cov,    # sample covariance of logit edge weights
  age_min_cat = sim_dyads$age_min_cat,      # age of younger dyad member
  age_max_cat = sim_dyads$age_max_cat,      # age of  older  dyad member
  node_1 = sim_dyads$node_1,                # node IDs for multimembership effects
  node_2 = sim_dyads$node_2,                # node IDs for multimembership effects
  #window = sim_dyads$window,                # ID of time window (random effect)
  #dyad_id = sim_dyads$dyad_id_all           # ID of dyad (random effect, can have the same value where same dyad seen multiple times)
  prior_min = rep(1, n_age_cat),
  prior_max = rep(1, n_age_cat)
)

## load dyadic regression model
dyadic_regression <- cmdstan_model('models/dyadic_regression_motnp.stan')
n_chains <- 4
n_samples <- 1000

## fit dyadic regression
fit_dyadreg_sim <- dyadic_regression$sample(
  data = dyad_data,
  iter_warmup = n_samples/2,
  iter_sampling = n_samples,
  chains = n_chains,
  parallel_chains = n_chains)

#### check outputs ####
# obtain summary
fit_dyadreg_sim$summary() %>% 
  filter(variable %in% c('beta_age_max','beta_age_min','intercept',
                         'delta_min[1]','delta_min[2]','delta_min[3]','delta_min[4]','delta_min[5]',
                         'delta_max[1]','delta_max[2]','delta_max[3]','delta_max[4]','delta_max[5]'))

## extract draws
draws <- fit_dyadreg_sim$draws(format = 'df')

## extract dyadic regression slopes
b_max <- draws$beta_age_max
b_min <- draws$beta_age_min
intercept <- draws$intercept
sigma <- draws$sigma
#rand_window <- draws[,c('rand_window[1]','rand_window[2]')]
#colnames(rand_window) <- 1:2
delta_min <- draws[,c('delta_min[1]','delta_min[2]','delta_min[3]','delta_min[4]','delta_min[5]')] ; colnames(delta_min) <- 1:n_age_cat
delta_max <- draws[,c('delta_max[1]','delta_max[2]','delta_max[3]','delta_max[4]','delta_max[5]')] ; colnames(delta_min) <- 1:n_age_cat
delta_j_min <- draws[,c('delta_j_min[1]','delta_j_min[2]','delta_j_min[3]','delta_j_min[4]','delta_j_min[5]','delta_j_min[6]')] ; colnames(delta_j_min) <- 1:(n_age_cat+1)
delta_j_max <- draws[,c('delta_j_max[1]','delta_j_max[2]','delta_j_max[3]','delta_j_max[4]','delta_j_max[5]','delta_j_max[6]')] ; colnames(delta_j_max) <- 1:(n_age_cat+1)
parameters <- data.frame(beta_age_max = b_max,
                         beta_age_min = b_min,
                         intercept = intercept,
                         sigma = sigma) %>% 
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>% 
  pivot_longer(cols = c('beta_age_max','beta_age_min','sigma','intercept'),
               names_to = 'parameter', values_to = 'slope_draw')

## traceplots -- quite wandery, but well mixed with one another -- could just be because there isn't much data?
ggplot(data = parameters)+
  geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')
delta_min %>%  as.data.frame() %>% 
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>% 
  pivot_longer(cols = all_of(1:n_age_cat),
               names_to = 'parameter', values_to = 'slope_draw') %>% 
  ggplot()+
  geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')
delta_max %>% as.data.frame() %>% 
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>% 
  pivot_longer(cols = all_of(1:n_age_cat),
               names_to = 'parameter', values_to = 'slope_draw') %>% 
  ggplot()+
  geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')
delta_j_min %>% as.data.frame() %>% 
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>% 
  pivot_longer(cols = all_of(2:(n_age_cat+1)),
               names_to = 'parameter', values_to = 'slope_draw') %>% 
  ggplot()+
  geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')
delta_j_max %>% as.data.frame(delta_j_max) %>% 
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>% 
  pivot_longer(cols = all_of(2:(n_age_cat+1)),
               names_to = 'parameter', values_to = 'slope_draw') %>% 
  ggplot()+
  geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')

#### plot predictions ####
## posterior predictive check
plot(density(as.numeric(sim_dat[1, ])),
     main = "Posterior predictive density of edge weights:\nblack = measured edge, red = predicted edge",
     ylim = c(0, 1), col = rgb(0, 0, 0, 0.25), las = 1)
for (i in 1:100) {
  j <- sample(1:1000, 1)
  
  mu_plot <- rep(NA, n_dyads)
  for(k in 1:n_dyads){
    mu_plot[k] <- intercept[j] + b_min[j]*sum(delta_j_min[j,(1:dyad_data$age_min_cat[k])]) + b_max[j]*sum(delta_j_max[j,(1:dyad_data$age_max_cat[k])]) #+ as.numeric(rand_window[j,eigen_list$window[k]]) + as.numeric(rand_node[j,eigen_list$nodes[k]])
  }
  
  sigma_plot <- dyad_data$logit_edge_cov + diag(rep(sigma[j], n_data))
  mv_norm <- MASS::mvrnorm(1, mu_plot, sigma_plot)
  
  lines(density(as.numeric(sim_dat[j, ])), col = rgb(0, 0, 0, 0.25)) # black lines for edge samples
  lines(density(mv_norm), col = rgb(1, 0, 0, 0.25))                  # red lines for predictions
  
}

# ## obtain median and 95% CI
# plot(density(b_max)) ; abline(v = 0, lty = 2)
# plot(density(b_min)) ; abline(v = 0, lty = 2)
# plot(density(intercept)) ; abline(v = 0, lty = 2)
# (b_max_summary <- round(quantile(b_max, probs=c(0.5, 0.025, 0.975)), 2)) ; sim_max
# (b_min_summary <- round(quantile(b_min, probs=c(0.5, 0.025, 0.975)), 2)) ; sim_min
# (intcp_summary <- round(quantile(intercept, probs=c(0.5, 0.025, 0.975)), 2)) ; sim_int
 
# ## convert to unstandardised scale
# b_max_unstd <- b_max * (sd(sim_dyads$mu) / sd(sim_dyads$age_max))
# b_min_unstd <- b_min * (sd(sim_dyads$mu) / sd(sim_dyads$age_min))
# intcp_unstd <- mean(sim_dyads$mu) - intercept * mean(sim_dyads$age_min) # this is missing age_max_std
# plot(density(b_max_unstd)) ; abline(v = sim_max, lty = 2)
# plot(density(b_min_unstd)) ; abline(v = sim_min, lty = 2)
# plot(density(intcp_unstd)) ; abline(v = sim_int, lty = 2)
# round(quantile(b_max_unstd, probs=c(0.5, 0.025, 0.975)), 2) ; sim_max
# round(quantile(b_min_unstd, probs=c(0.5, 0.025, 0.975)), 2) ; sim_min
# round(quantile(intcp_unstd, probs=c(0.5, 0.025, 0.975)), 2) ; sim_int
 
## calculate predictions
pred <- data.frame(age_min = rep(rep(unique(sim_dyads$age_min_cat), each = length(unique(sim_dyads$age_max_cat))),
                                 n_chains*n_samples),
                   age_max = rep(rep(unique(sim_dyads$age_max_cat), length(unique(sim_dyads$age_min_cat))),
                                 n_chains*n_samples),
                   #window = rep(rep(c(1,2), 10*10*2), each = 5),
                   #dyad_id = rep(sample(dyad_data$dyad_id, 5, replace = F),10*10*2)
                   intcp = rep(intercept, each = length(unique(sim_dyads$age_max_cat))*length(unique(sim_dyads$age_min_cat))),
                   beta_min = rep(b_min, each = length(unique(sim_dyads$age_max_cat))*length(unique(sim_dyads$age_min_cat))),
                   beta_max = rep(b_max, each = length(unique(sim_dyads$age_max_cat))*length(unique(sim_dyads$age_min_cat))),
                   delta_j_min = NA, delta_j_max = NA, pred = NA) %>%
  filter(age_min <= age_max)
for(i in 1:ncol(delta_j_max)){
  pred$delta_j_min[pred$age_min == i] <- rowSums(delta_j_min[,1:i])
  pred$delta_j_max[pred$age_max == i] <- rowSums(delta_j_max[,1:i])
}
pred$pred <- pred$intcp + pred$beta_min * pred$delta_j_min + pred$beta_max * pred$delta_j_max

## convert to unstandardised scale
pred$pred_unstd <- pred$pred*sd(sim_dyads$mu) + mean(sim_dyads$mu)

## summarise
pred_summary <- pred %>% 
  group_by(age_min, age_max) %>% 
  mutate(pred_lwr = quantile(pred, probs = 0.025),
         pred_mean = mean(pred),
         pred_upr = quantile(pred, probs = 0.975),
         pred_unstd_lwr = quantile(pred_unstd, probs = 0.025),
         pred_unstd_mean = mean(pred_unstd),
         pred_unstd_upr = quantile(pred_unstd, probs = 0.975)) %>% 
  ungroup() %>% 
  select(age_min, age_max,
         pred_lwr, pred_mean, pred_upr,
         pred_unstd_lwr, pred_unstd_mean, pred_unstd_upr) %>% 
  distinct()
  
## plot
ggplot()+
  geom_ribbon(data = pred_summary,
              aes(x = age_min,
                  #ymin = pred_lwr, ymax = pred_upr,
                  ymin = invlogit(pred_lwr), ymax = invlogit(pred_upr),
                  group = as.factor(age_max), fill = as.factor(age_max)),
              alpha = 0.3)+
  geom_line(data = pred_summary,
            aes(x = age_min,
                #y = pred_mean,
                y = invlogit(pred_mean),
                colour = as.factor(age_max), group = as.factor(age_max)),
            linewidth = 1)+
  scale_colour_viridis_d(direction = -1)+ scale_fill_viridis_d(direction = -1)+
  geom_point(data = sim_dyads, aes(x = age_min_cat, y = invlogit(mu_std), colour = as.factor(age_max_cat)))+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+
  labs(colour = 'maximum age', fill = 'maximum age')

raw_summary <- sim_dyads %>% 
  select(age_min_cat, age_max_cat, mu, mu_std) %>% 
  group_by(age_min_cat, age_max_cat) %>% 
  mutate(mu_mean = mean(mu),
         mu_std_mean = mean(mu_std)) %>% 
  ungroup() %>% 
  select(-mu, -mu_std) %>% 
  distinct() %>% 
  mutate(age_pair = paste0(age_min_cat, '_', age_max_cat))
compare <- pred_summary %>% 
  mutate(age_pair = paste0(age_min, '_', age_max)) %>% 
  left_join(raw_summary[,c('mu_mean','mu_std_mean','age_pair')], by = 'age_pair') %>% 
  rename(raw_unstd = mu_mean, raw_std = mu_std_mean) %>% 
  select(age_min, age_max, raw_unstd, pred_unstd_mean, raw_std, pred_mean, pred_lwr, pred_unstd_lwr, pred_upr, pred_unstd_upr)
ggplot(compare)+
  geom_ribbon(aes(x = raw_unstd,
                  ymin = pred_unstd_lwr, ymax = pred_unstd_upr),
              alpha = 0.3)+
  geom_line(aes(x = raw_unstd, y = pred_unstd_mean),
            linewidth = 1)+
  # geom_point(aes(x = mu_raw, y = mu_mid_unstd),
  #            alpha = 0.1)+
  geom_line(data = data.frame(x = c(15,85), y = c(15,85)),
            aes(x = x, y = y),
            linewidth = 0.5)+
  scale_x_continuous('raw mean')+
  scale_y_continuous('predicted mean')+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

#### clean up ####
rm(list = ls()[! ls() %in% c('sim_dat','sim_dyads','sim_nodes','max_age','min_age','n_dyads','n_nodes','sim_max','sim_min','sim_int')]) ; gc()

