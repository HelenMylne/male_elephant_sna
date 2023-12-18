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
n_nodes <- ((max_age - min_age)/2)+1        # total nodes = 2 per age
sim_nodes <- data.frame(node = 1:n_nodes,   # create data frame of individuals
                         age = seq(from = min_age, 
                                   to = max_age, by = 2)) %>% 
  mutate(age_cat = ifelse(age <= 15, 1,
                          ifelse(age <= 20, 2,
                                 ifelse(age <= 25, 3,
                                        ifelse(age <= 40, 4, 5)))))

## convert to dyadic dataframe
sim_dyads <- data.frame(node = rep(sim_nodes$node, each = n_nodes),
                         node_2 = rep(sim_nodes$node, n_nodes)) %>% 
  filter(node != node_2) %>% filter(node < node_2) %>% 
  left_join(sim_nodes, by = 'node') %>% 
  rename(node_1 = node,
         age_1 = age,
         age_cat_1 = age_cat) %>% 
  rename(node = node_2) %>%
  left_join(sim_nodes, by = 'node') %>% 
  rename(node_2 = node,
         age_2 = age,
         age_cat_2 = age_cat) %>% 
  mutate(age_min = NA, age_max = NA, mu = NA, sd = NA) %>% 
  mutate(dyad = factor(paste0(node_1, '_', node_2)),
         dyad_id = as.integer(dyad))
for(i in 1:nrow(sim_dyads)){
  sim_dyads$age_min[i] <- min(sim_dyads$age_cat_1[i], sim_dyads$age_cat_2[i])
  sim_dyads$age_max[i] <- max(sim_dyads$age_cat_1[i], sim_dyads$age_cat_2[i])
}
n_dyads <- nrow(sim_dyads)

## simulate age effect
sim_min <- 1
sim_max <- 0.5
sim_int <- -3
sim_dyads$mu <- sim_int + sim_dyads$age_min * sim_min + sim_dyads$age_max * sim_max   # simulate mean centrality on normal scale
plot(sim_dyads$mu ~ sim_dyads$age_min)               # plot
plot(sim_dyads$mu ~ sim_dyads$age_max)               # plot

## standardise edges
sim_dyads$mu_std <- ( sim_dyads$mu - mean(sim_dyads$mu) ) / sd(sim_dyads$mu)

## simulate full distribution of samples per node
sim_dyads$sd <- abs(sim_max/3)           # make small to start with to be sure model should be able to detect difference
sim_dat <- matrix(data = NA, nrow = 1000, ncol = n_dyads, dimnames = list(1:1000, sim_dyads$dyad_id))    # create matrix
for(j in 1:n_dyads){
  sim_dat[,j] <- rnorm(n = nrow(sim_dat), mean = sim_dyads$mu_std[j], sd = sim_dyads$sd[j])  # simulate distribution
}
plot(sim_dat[1,] ~ sim_dyads$age_min)          # plot simulated values against minimum age
plot(sim_dat[1,] ~ sim_dyads$age_max)          # plot simulated values against maximum age

#### plot raw data ####
ggplot()+
  geom_point(data = sim_dyads, aes(x = age_min, y = mu_std, colour = as.factor(age_max)), shape = 19)+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')

sim_dat_long <- sim_dat %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'dyad_id', values_to = 'edge_draw') %>% 
  mutate(dyad_id = as.numeric(dyad_id)) %>% 
  left_join(sim_dyads, by = 'dyad_id')
ggplot()+
  geom_point(data = sim_dat_long, aes(x = age_min, y = edge_draw, colour = as.factor(age_max)),
             shape = 19, alpha = 0.05)+
  geom_point(data = sim_dyads, aes(x = age_min, y = mu_std),
             shape = 19, colour = 'white', size = 1)+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')

## heat map: x = min age, y = max age, colour = edge weight
ggplot()+
  geom_tile(data = sim_dyads, mapping = aes(x = age_min, y = age_max, fill = mu_std))+
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
rm(cols, fitted_values, i, j, max_age, min_age, mu, num_check, rows, sd, selected_samples) ; gc()

#### prior predictive check ####
n <- 100
n_age_cat <- length(unique(sim_dyads$age_max))
beta_age_min <- rnorm(n, 0, 1.5)
beta_age_max <- rnorm(n, 0, 1.5)
intercept <- rnorm(n, 0, 1)
min_dirichlet <- rdirichlet(n, rep(1, n_age_cat))
max_dirichlet <- rdirichlet(n, rep(1, n_age_cat))
plot(NULL, las = 1, xlab = 'minimum age (standardised)', ylab = 'logit edge weight (standardised)',
     ylim = c(min(logit_edge_draws_mu)-5, max(logit_edge_draws_mu)+5), xlim = c(1,n_age_cat))
abline(h = min(logit_edge_draws_mu), lty = 2) ; abline(h = max(logit_edge_draws_mu), lty = 2)
x <- 1:n_age_cat
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
dyad_data <- list(
  num_dyads = n_dyads,                      # number of dyads
  num_nodes = n_nodes,                      # number of nodes
  num_age_cat = n_age_cat,                  # number of unique age categories
  length_dirichlet = n_age_cat + 1,         # number of unique age categories + 1
  logit_edge_mu = logit_edge_draws_mu,      # sample means of the logit edge weights
  logit_edge_cov = logit_edge_draws_cov,    # sample covariance of logit edge weights
  age_min_cat = sim_dyads$age_min,          # age of younger dyad member
  age_max_cat = sim_dyads$age_max,          # age of  older  dyad member
  node_1 = sim_dyads$node_1,                # node IDs for multimembership effects
  node_2 = sim_dyads$node_2,                # node IDs for multimembership effects
  prior_min = rep(1, n_age_cat),            # prior for minimum age slope
  prior_max = rep(1, n_age_cat)             # prior for maximum age slope
)

## load dyadic regression model
dyadic_regression <- cmdstan_model('models/dyadic_regression_motnp.stan')
n_chains <- 4
n_samples <- 1000

## fit dyadic regression
fit_dyadreg_sim <- dyadic_regression$sample(
  data = dyad_data,
  iter_warmup = n_samples,
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
    mu_plot[k] <- intercept[j] + b_min[j]*sum(delta_j_min[j,(1:dyad_data$age_min_cat[k])]) + b_max[j]*sum(delta_j_max[j,(1:dyad_data$age_max_cat[k])])
  }
  
  sigma_plot <- dyad_data$logit_edge_cov + diag(rep(sigma[j], n_dyads))
  mv_norm <- MASS::mvrnorm(1, mu_plot, sigma_plot)
  
  lines(density(as.numeric(sim_dat[j, ])), col = rgb(0, 0, 0, 0.25)) # black lines for edge samples
  lines(density(mv_norm), col = rgb(1, 0, 0, 0.25))                  # red lines for predictions
  
}

## create predictive data frame
pred <- data.frame(age_min = rep(rep(unique(sim_dyads$age_min), each = length(unique(sim_dyads$age_max))),
                                 n_chains*n_samples),
                   age_max = rep(rep(unique(sim_dyads$age_max), length(unique(sim_dyads$age_min))),
                                 n_chains*n_samples),
                   intcp = rep(intercept, each = length(unique(sim_dyads$age_max))*length(unique(sim_dyads$age_min))),
                   beta_min = rep(b_min, each = length(unique(sim_dyads$age_max))*length(unique(sim_dyads$age_min))),
                   beta_max = rep(b_max, each = length(unique(sim_dyads$age_max))*length(unique(sim_dyads$age_min))),
                   delta_j_min = NA, delta_j_max = NA) #%>%
  #filter(age_min <= age_max)
for(i in 1:(n_samples*n_chains)){
  for(j in 1:ncol(delta_j_max)){
    pred$delta_j_min[pred$age_min == j] <- rowSums(delta_j_min[i,1:j])
    pred$delta_j_max[pred$age_max == j] <- rowSums(delta_j_max[i,1:j])
  }
}

## calculate predictions
pred$pred <- pred$intcp + pred$beta_min * pred$delta_j_min + pred$beta_max * pred$delta_j_max

## convert to unstandardised scale
pred$pred_unstd <- pred$pred*sd(sim_dyads$mu) + mean(sim_dyads$mu)

## summarise
pred_summary_all <- pred %>% 
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
pred_summary_possible <- pred %>% 
  filter(age_min <= age_max) %>% 
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
  geom_ribbon(data = pred_summary_possible,
              aes(x = age_min,
                  #ymin = pred_lwr, ymax = pred_upr,
                  ymin = invlogit(pred_lwr), ymax = invlogit(pred_upr),
                  group = as.factor(age_max), fill = as.factor(age_max)),
              alpha = 0.3)+
  geom_line(data = pred_summary_possible,
            aes(x = age_min,
                #y = pred_mean,
                y = invlogit(pred_mean),
                colour = as.factor(age_max), group = as.factor(age_max)),
            linewidth = 1)+
  scale_colour_viridis_d(direction = -1)+ scale_fill_viridis_d(direction = -1)+
  geom_point(data = sim_dyads, aes(x = age_min, y = invlogit(mu_std), colour = as.factor(age_max)))+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+
  labs(colour = 'maximum age', fill = 'maximum age')

raw_summary <- sim_dyads %>% 
  select(age_min, age_max, mu, mu_std) %>% 
  group_by(age_min, age_max) %>% 
  mutate(mu_mean = mean(mu),
         mu_std_mean = mean(mu_std)) %>% 
  ungroup() %>% 
  select(-mu, -mu_std) %>% 
  distinct() %>% 
  mutate(age_pair = paste0(age_min, '_', age_max))
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
  geom_line(data = data.frame(x = c(min(sim_dyads$mu),max(sim_dyads$mu)),
                              y = c(min(sim_dyads$mu),max(sim_dyads$mu))),
             aes(x = x, y = y),
             linewidth = 0.5, colour = 'red')+
  scale_x_continuous('raw mean')+
  scale_y_continuous('predicted mean')+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

## extract slope estimates -- this is definitely wrong but I can't work out how to do it properly
contrasts <- data.frame(intercept = intercept, b_min = b_min, b_max = b_max) %>% 
  mutate(diff_1_2_min = NA, diff_2_3_min = NA, diff_3_4_min = NA, diff_4_5_min = NA,
         diff_1_2_max = NA, diff_2_3_max = NA, diff_3_4_max = NA, diff_4_5_max = NA)
for(i in 1:nrow(contrasts)){
  x <- pred %>%
    filter(intcp == contrasts$intercept[i] & 
             beta_min == contrasts$b_min[i] &
             beta_max == contrasts$b_max[i])
  a <- x %>% filter(age_min == 1) ; b <- x %>% filter(age_min == 2) ; c <- x %>% filter(age_min == 3) ; d <- x %>% filter(age_min == 4) ; e <- x %>% filter(age_min == 5)
  contrasts$diff_1_2_min[i] <- mean(a$pred_unstd - b$pred_unstd)
  contrasts$diff_2_3_min[i] <- mean(a$pred_unstd - b$pred_unstd)
  contrasts$diff_3_4_min[i] <- mean(a$pred_unstd - b$pred_unstd)
  contrasts$diff_4_5_min[i] <- mean(a$pred_unstd - b$pred_unstd)
  a <- x %>% filter(age_max == 1) ; b <- x %>% filter(age_max == 2) ; c <- x %>% filter(age_max == 3) ; d <- x %>% filter(age_max == 4) ; e <- x %>% filter(age_max == 5)
  contrasts$diff_1_2_max[i] <- mean(a$pred_unstd - b$pred_unstd)
  contrasts$diff_2_3_max[i] <- mean(a$pred_unstd - b$pred_unstd)
  contrasts$diff_3_4_max[i] <- mean(a$pred_unstd - b$pred_unstd)
  contrasts$diff_4_5_max[i] <- mean(a$pred_unstd - b$pred_unstd)
}
rm(a,b,c,d,e,x) ; gc()
contrasts$diff_1_2_min_year <- contrasts$diff_1_2_min/(16-10)
contrasts$diff_1_2_max_year <- contrasts$diff_1_2_max/(16-10)
contrasts$diff_2_3_min_year <- contrasts$diff_2_3_min/(21-16)
contrasts$diff_2_3_max_year <- contrasts$diff_2_3_max/(21-16)
contrasts$diff_3_4_min_year <- contrasts$diff_3_4_min/(26-21)
contrasts$diff_3_4_max_year <- contrasts$diff_3_4_max/(26-21)
contrasts$diff_4_5_min_year <- contrasts$diff_4_5_min/(41-26)
contrasts$diff_4_5_max_year <- contrasts$diff_4_5_max/(41-26)

sim_min  # original (true) slope value
mean(c(contrasts$diff_1_2_min_year,contrasts$diff_2_3_min_year,contrasts$diff_3_4_min_year, contrasts$diff_4_5_min_year))

sim_max  # original (true) slope value
mean(c(contrasts$diff_1_2_max_year,contrasts$diff_2_3_max_year,contrasts$diff_3_4_max_year, contrasts$diff_4_5_max_year))

## extract model fit
summary <- fit_dyadreg_sim$summary()
par(mfrow = c(3,1))
hist(summary$rhat, breaks = 50)
hist(summary$ess_bulk, breaks = 50)
hist(summary$ess_tail, breaks = 50)
par(mfrow = c(1,1))

## summarise predictions
sum_pred <- pred %>% 
  group_by(age_min, age_max) %>% 
  mutate(mean_std = mean(pred),
         mean_unstd = mean(pred_unstd)) %>% 
  select(age_min, age_max, mean_std, mean_unstd) %>% 
  distinct()
ggplot(sum_pred)+
  geom_point(aes(x = age_min, y = mean_std, colour = as.factor(age_max)))

