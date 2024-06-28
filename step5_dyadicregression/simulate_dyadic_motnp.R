#### set up ####
# script to simulate dyadic regression and check actually working
#library(tidyverse) ; library(LaplacesDemon) ; library(car) ; library(StanHeaders) ; library(rstan)
# library(cmdstanr) ; library(bisonR) ; library(brms)
#library(cmdstanr, lib.loc = '../packages/')
library(StanHeaders, lib.loc = '../packages/')
library(rstan, lib.loc = '../packages/')
library(tidyverse, lib.loc = '../packages/')
library(car, lib.loc = '../packages/')
library(LaplacesDemon, lib.loc = '../packages/')
library(float, lib.loc = '../packages/')

#set_cmdstan_path('R:/rsrch/df525/phd/hkm513/packages/.cmdstan2/cmdstan-2.33.1/')
theme_set(theme_bw())

pdf('step5_dyadicregression/simulate_dyadic_motnp_prep.pdf')

#### load raw data as basis for simulation ####
load('motnp_edgeweights_conditionalprior.RData')
sim <- counts_df %>%
  dplyr::select(dyad_males, node_1, node_2, node_1_males, node_2_males, id_1, id_2, count_dyad, event_count, apart, count_1, count_2, age_category_1, age_category_2, age_cat_id_1, age_cat_id_2) %>%
  rename(dyad_id = dyad_males, #node_1 = node_1_males, node_2 = node_2_males,
         total_sightings = count_dyad, together = event_count,
         age_cat_1 = age_category_1, age_cat_2 = age_category_2,
         age_num_1 = age_cat_id_1, age_num_2 = age_cat_id_2) %>%
  mutate(age_num_1 = as.numeric(age_num_1) - 2,
         age_num_2 = as.numeric(age_num_2) - 2)

## simulate actual ages from categorical
nodes <- nodes %>%
  mutate(age_cat = ifelse(age < 15, 1,
                          ifelse(age < 20, 2,
                                 ifelse(age < 25, 3,
                                        ifelse(age < 40, 4, 5))))) %>%
  mutate(min_age = ifelse(age_cat == 1, 10,
                            ifelse(age_cat == 2, 16,
                                   ifelse(age_cat == 3, 21,
                                          ifelse(age_cat == 4, 26, 40)))),
         max_age = ifelse(age_cat == 1, 15,
                            ifelse(age_cat == 2, 20,
                                   ifelse(age_cat == 3, 25,
                                          ifelse(age_cat == 4, 40, 60)))))
nodes$age_sim <- NA
for(i in 1:nrow(nodes)){
  sd <- (nodes$max_age[i] - nodes$min_age[i])/2
  age <- round(rnorm(1000, mean = nodes$age[i], sd = sd),0)
  age <- age[which(age < nodes$max_age[i] & age > nodes$min_age[i])[1]]
  nodes$age_sim[i] <- age
}
rm(list = ls()[! ls() %in% c('nodes', 'sim','counts_df')]) ; gc()

# ## filter down to only 50 elephants as proof of concept before trying to work out how much memory this requires
# percent_1 <- ( length(which(nodes$age_cat == 1))/nrow(nodes) ) * 100 # what proportion of elephants are in age category 1
# percent_2 <- ( length(which(nodes$age_cat == 2))/nrow(nodes) ) * 100 # proportion in cat 2
# percent_3 <- ( length(which(nodes$age_cat == 3))/nrow(nodes) ) * 100 # proportion in cat 3
# percent_4 <- ( length(which(nodes$age_cat == 4))/nrow(nodes) ) * 100 # proportion in cat 4
# percent_5 <- ( length(which(nodes$age_cat == 5))/nrow(nodes) ) * 100 # proportion in cat 5
# nodes_to_include <- c(sample(x = nodes$id[nodes$age_cat == 1],            # sample individuals from population proportional to natural ages
#                              replace = F, size = round(percent_1/2,0)),
#                       sample(x = nodes$id[nodes$age_cat == 2],
#                              replace = F, size = round(percent_2/2,0)),
#                       sample(x = nodes$id[nodes$age_cat == 3],
#                              replace = F, size = round(percent_3/2,0)),
#                       sample(x = nodes$id[nodes$age_cat == 4],
#                              replace = F, size = round(percent_4/2,0)),
#                       sample(x = nodes$id[nodes$age_cat == 5],
#                              replace = F, size = round(percent_5/2,0)))
# nodes <- nodes %>%
#   filter(id %in% nodes_to_include)
# sim <- sim %>%
#   filter(id_1 %in% nodes_to_include) %>%
#   filter(id_2 %in% nodes_to_include)

## randomise node IDs
nodes$node_rand_1 <- sample(size = nrow(nodes), x = 1:nrow(nodes), replace = F)
nodes$node_rand_2 <- nodes$node_rand_1

## randomise dyad IDs
sim$dyad_rand <- sample(size = nrow(sim), x = 1:nrow(sim), replace = F)

## combine with dyad data
nodes <- nodes %>%
  mutate(node_1 = node, node_2 = node,
         age_sim_1 = age_sim, age_sim_2 = age_sim)
sim <- sim %>%
  left_join(nodes[,c('node_1', 'age_sim_1', 'node_rand_1')], by = 'node_1') %>%
  left_join(nodes[,c('node_2', 'age_sim_2', 'node_rand_2')], by = 'node_2')

## define population parameters
n_dyads <- nrow(sim)
n_nodes <- nrow(nodes)

## define oldest and youngest
sim$age_max <- NA ; sim$age_min <- NA
for(i in 1:nrow(sim)){
  sim$age_max[i] <- max(c(sim$age_sim_1[i], sim$age_sim_2[i]))
  sim$age_min[i] <- min(c(sim$age_sim_1[i], sim$age_sim_2[i]))
}

## print progress marker
print('data loaded')

#### simulate age effect ####
sim_min <- -0.2
sim_max <- 0.8
sim_int <- 5
sim$mu <- sim_int + sim$age_min * sim_min + sim$age_max * sim_max   # simulate mean centrality on normal scale
ggplot(sim)+
  geom_point(aes(x = age_max, y = mu, colour = age_min))

## standardise edges
sim$mu_std <- ( sim$mu - mean(sim$mu) ) / sd(sim$mu)
ggplot(sim)+
  geom_point(aes(x = age_max, y = mu_std, colour = age_min))

## simulate full distribution of samples per node
sim$sd <- abs(sim_max/3)           # make small to start with to be sure model should be able to detect difference
sim_dat <- matrix(data = NA, nrow = 1000, ncol = n_dyads, dimnames = list(1:1000, sim$dyad_id))    # create matrix
for(j in 1:n_dyads){
  sim_dat[,j] <- rnorm(n = nrow(sim_dat), mean = sim$mu_std[j], sd = sim$sd[j])  # simulate distribution
}
plot(sim_dat[1,] ~ sim$age_min)          # plot simulated values against minimum age
plot(sim_dat[1,] ~ sim$age_max)          # plot simulated values against maximum age

# ## convert to float to reduce memory requirements
# x <- object.size(sim_dat) ; format(x, units = 'Mb')         # 9.5 Mb
# sim_dat_fl <- float::fl(sim_dat)
# x <- object.size(sim_dat_fl) ; format(x, units = 'Mb')      # 4.8 Mb
# x <- object.size(sim_dat_fl@Data) ; format(x, units = 'Mb') # 4.8 Mb

## print progress marker
print('edges simulated')

#### plot raw data ####
sim$age_cat_max <- NA ; sim$age_cat_min <- NA
for(i in 1:nrow(sim)){
  sim$age_cat_max[i] <- max(sim$age_num_1[i], sim$age_num_2[i])
  sim$age_cat_min[i] <- min(sim$age_num_1[i], sim$age_num_2[i])
}

ggplot()+
  geom_point(data = sim, aes(x = age_min,
                             y = mu_std,
                             colour = as.factor(age_cat_max)),
             shape = 19)+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')

ggplot()+
  geom_boxplot(data = sim, aes(x = as.factor(age_cat_min),
                               y = mu_std,
                               fill = as.factor(age_cat_max)),
             shape = 19)+
  scale_x_discrete('age category of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  theme(legend.position = 'bottom')+
  labs(fill = 'age category of older male')

# sim_dat_fl_long <- sim_dat_fl %>%
sim_dat_long <- sim_dat %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = 'dyad_id', values_to = 'edge_draw') %>%
  mutate(dyad_id = as.numeric(dyad_id)) %>%
  left_join(sim, by = 'dyad_id')
ggplot()+
  geom_violin(data = sim_dat_long, #sim_dat_fl_long,
              aes(x = as.factor(age_cat_min),
                  y = edge_draw,
                  fill = as.factor(age_cat_max)),
              #shape = 19,
              alpha = 0.5)+
  geom_point(data = sim,
             aes(x = as.factor(age_cat_min), y = mu_std,
                 group = as.factor(age_cat_max)),
             shape = 19,
             #colour = 'white',
             size = 1)+
  scale_x_discrete('age category of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  theme(legend.position = 'bottom')+
  labs(fill = 'age category of older elephant')
rm(sim_dat_fl_long, sim_dat_long) ; gc()

## print progress marker
print('edges plotted')

#### fit multivariate Gaussian distribution to output of edge weight model ####
## fit a multivariate normal dist to the edges -- quantities will be given to Stan model as data to model joint posteriors of edge weight in the regression
# logit_edge_draws_mu_fl <- apply(sim_dat_fl@Data, 2, mean)
# logit_edge_draws_cov_fl <- as.matrix(float::dbl(cov(sim_dat_fl@Data)))
logit_edge_draws_mu <- apply(sim_dat, 2, mean)
logit_edge_draws_cov <- cov(sim_dat)

# logit_edge_draws_mu_fl <- float::fl(logit_edge_draws_mu)
# object.size(logit_edge_draws_mu) - object.size(logit_edge_draws_mu_fl)
# logit_edge_draws_cov_fl <- float::fl(logit_edge_draws_cov)
# x <- object.size(logit_edge_draws_cov) - object.size(logit_edge_draws_cov_fl)
# format(x, units = 'Mb')      # 4.8 Mb

## randomly select samples to examine
num_check <- 20
selected_samples <- sample(1:(n_dyads-1), num_check, replace = FALSE)

## set grid layout
rows <- floor(sqrt(num_check))
cols <- ceiling(num_check / rows)
par(mfrow=c(rows, cols), mar=c(2,2,2,1))

## plot
for (i in selected_samples) {
  mu <- logit_edge_draws_mu[i]
  sd <- sqrt(logit_edge_draws_cov[i,i])

  fitted_values <- rnorm(1e5, mean = mu, sd = sd)

  hist(sim_dat[,i], probability = TRUE, las = 1,
       main = paste("Dyad", i), xlab = "Value", breaks = 50)
  lines(density(fitted_values), col="blue", lw=1.5)
}

# for (i in selected_samples) {
#   mu <- logit_edge_draws_mu_fl[i]
#   sd <- sqrt(logit_edge_draws_cov_fl[i,i])
#
#   fitted_values <- rnorm(1e5, mean = mu, sd = sd)
#
#   hist(sim_dat_fl[,i], probability = TRUE, las = 1,
#        main = paste("Dyad", i), xlab = "Value", breaks = 50)
#   lines(density(fitted_values), col="blue", lw=1.5)
# }

for (i in selected_samples) {
  mu <- logit_edge_draws_mu[i]
  sd <- sqrt(logit_edge_draws_cov[i,i])

  fitted_values <- rnorm(1e5, mean = mu, sd = sd)

  plot(unlist(sim_dat[,i]), unlist(sim_dat[,i+1]), col = rgb(0,0,1,0.05), las = 1,
       main = paste("cov ", i ,"&",i+1))
}

## reset plot window and clean up
par(mfrow=c(1,1), mai = c(1,1,0.5,0.5))
rm(cols, fitted_values, i, j, mu, num_check, rows, sd, selected_samples) ; gc()

## print progress marker
print('normal approximation complete')

#### prior predictive check ####
n <- 100
n_age_cat <- length(unique(sim$age_cat_max))
beta_age_min <- rnorm(n, 0, 1.5)
beta_age_max <- rnorm(n, 0, 1.5)
intercept <- rnorm(n, 0, 2)
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

## print progress marker
print('prior predictive complete')

#### fit dyadic regression ####
## create data list
dyad_data <- list(
  num_dyads = n_dyads,                      # number of dyads
  num_nodes = n_nodes,                      # number of nodes
  num_age_cat = n_age_cat,                  # number of unique age categories
  length_dirichlet = n_age_cat + 1,         # number of unique age categories + 1
  logit_edge_mu = logit_edge_draws_mu,      # sample means of the logit edge weights
  logit_edge_cov = logit_edge_draws_cov,    # sample covariance of logit edge weights
  age_min_cat = sim$age_cat_min,            # age of younger dyad member
  age_max_cat = sim$age_cat_max,            # age of  older  dyad member
  node_1 = sim$node_rand_1,                 # node IDs for multimembership effects
  node_2 = sim$node_rand_2,                 # node IDs for multimembership effects
  prior_min = rep(1, n_age_cat),            # prior for minimum age slope
  prior_max = rep(1, n_age_cat)             # prior for maximum age slope
)

## print progress marker
print('data list created')

## load dyadic regression model
dyadic_regression <- stan_model('models/dyadic_regression_motnp.stan')
#dyadic_regression <- cmdstan_model('models/dyadic_regression_motnp.stan')

## print progress marker
print('model loaded in')

## define parameters
n_chains <- 4
n_samples <- 1000

## fit dyadic regression
print('start model')
fit_dyadreg_sim <- sampling(dyadic_regression,
                            data = dyad_data,
                            iter = n_samples*2, warmup = n_samples,
                            chains = n_chains, cores = n_chains)
# fit_dyadreg_sim <- dyadic_regression$sample(
#   data = dyad_data,
#   iter_warmup = n_samples,
#   iter_sampling = n_samples,
#   chains = n_chains,
#   parallel_chains = n_chains)

## save output
save.image('step5_dyadicregression/motnp_dyadic_simulation.RData')

## print progress marker
print('model run')

#### check outputs ####
dev.off()
pdf('step5_dyadicregression/simulate_dyadic_motnp_checkoutputs.pdf')
#load('step5_dyadicregression/motnp_dyadic_simulation.RData')

## extract draws
posterior_samples <- fit_dyadreg_sim %>%
  as.data.frame()
# draws <- fit_dyadreg_sim$draws(format = 'df')

## extract model fit
# summary <- fit_dyadreg_sim$summary()
s <- summary(fit_dyadreg_sim)
summary <- s$summary %>%
  as.data.frame() %>%
  filter(is.nan(se_mean) == FALSE)
summary
# par(mfrow = c(3,1))
# hist(summary$rhat, breaks = 50)
# hist(summary$ess_bulk, breaks = 50)
# hist(summary$ess_tail, breaks = 50)
par(mfrow = c(2,1))
hist(summary$Rhat,  breaks = 50)
hist(summary$n_eff, breaks = 50)
par(mfrow = c(1,1))

save.image('step5_dyadicregression/motnp_dyadic_simulation.RData')

## extract dyadic regression slopes
b_max <- posterior_samples$beta_age_max
b_min <- posterior_samples$beta_age_min
intercept <- posterior_samples$intercept
sigma <- posterior_samples$sigma
delta_min <- posterior_samples[,c('delta_min[1]','delta_min[2]','delta_min[3]','delta_min[4]','delta_min[5]')] ; colnames(delta_min) <- 1:n_age_cat
delta_max <- posterior_samples[,c('delta_max[1]','delta_max[2]','delta_max[3]','delta_max[4]','delta_max[5]')] ; colnames(delta_min) <- 1:n_age_cat
delta_j_min <- posterior_samples[,c('delta_j_min[1]','delta_j_min[2]','delta_j_min[3]','delta_j_min[4]','delta_j_min[5]','delta_j_min[6]')] ; colnames(delta_j_min) <- 1:(n_age_cat+1)
delta_j_max <- posterior_samples[,c('delta_j_max[1]','delta_j_max[2]','delta_j_max[3]','delta_j_max[4]','delta_j_max[5]','delta_j_max[6]')] ; colnames(delta_j_max) <- 1:(n_age_cat+1)
mm_nodes <- posterior_samples[,grep(pattern = 'mm_nodes', x = colnames(posterior_samples))]
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
  geom_line(aes(x = position,
                y = slope_draw,
                colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')
delta_min %>%  as.data.frame() %>%
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>%
  pivot_longer(cols = all_of(1:n_age_cat),
               names_to = 'parameter',
               values_to = 'slope_draw') %>%
  ggplot()+
  geom_line(aes(x = position,
                y = slope_draw,
                colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')
delta_max %>% as.data.frame() %>%
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>%
  pivot_longer(cols = all_of(1:n_age_cat),
               names_to = 'parameter',
               values_to = 'slope_draw') %>%
  ggplot()+
  geom_line(aes(x = position,
                y = slope_draw,
                colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')
delta_j_min %>% as.data.frame() %>%
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>%
  pivot_longer(cols = all_of(2:(n_age_cat+1)),
               names_to = 'parameter',
               values_to = 'slope_draw') %>%
  ggplot()+
  geom_line(aes(x = position,
                y = slope_draw,
                colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')
delta_j_max %>% as.data.frame() %>%
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>%
  pivot_longer(cols = all_of(2:(n_age_cat+1)),
               names_to = 'parameter',
               values_to = 'slope_draw') %>%
  ggplot()+
  geom_line(aes(x = position,
                y = slope_draw,
                colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')

#### posterior predictive check ####
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

## save output
save.image('step5_dyadicregression/motnp_dyadic_simulation.RData')
dev.off()

## print progress marker
print('outputs checked')

#### predict from raw data ####
pdf('step5_dyadicregression/simulate_dyadic_motnp_predictions.pdf')
load('step5_dyadicregression/motnp_dyadic_simulation.RData')

# ## MADE UP VALUES SO THAT I CAN CHECK THE CODE WORKS WITHOUT LOADING THE FULL WORKSPASCE ONTO THE DESKTOP BECAUSE THEN IT THROWS A WOBBLY. DO NOT KEEP THIS SECTION!
# nodes <- nodes[1:10,]
# sim <- sim %>%
#   filter(id_1 %in% nodes$id & id_2 %in% nodes$id)
# b_max <- rnorm(n = n_chains*n_samples, mean = 4, sd = 0.2)
# b_min <- rnorm(n = n_chains*n_samples, mean = -0.7, sd = 0.1)
# intercept <- rnorm(n = n_chains*n_samples, mean = -1.8, sd = 0.1)
# sigma <- rnorm(n = n_chains*n_samples, mean = 0.022, sd = 0.005)
# delta_min <- data.frame(
#   `1` = abs(rnorm(n = n_chains*n_samples, mean = 0.25, sd = 0.1)),
#   `2` = abs(rnorm(n = n_chains*n_samples, mean = 0.3, sd = 0.1)),
#   `3` = abs(rnorm(n = n_chains*n_samples, mean = 0.15, sd = 0.05)),
#   `4` = abs(rnorm(n = n_chains*n_samples, mean = 0.2, sd = 0.1)),
#   `5` = abs(rnorm(n = n_chains*n_samples, mean = 0.2, sd = 0.1))) %>%
#   as.matrix()
# delta_max <- data.frame(
#   `1` = abs(rnorm(n = n_chains*n_samples, mean = 0.12, sd = 0.02)),
#   `2` = abs(rnorm(n = n_chains*n_samples, mean = 0.16, sd = 0.02)),
#   `3` = abs(rnorm(n = n_chains*n_samples, mean = 0.32, sd = 0.03)),
#   `4` = abs(rnorm(n = n_chains*n_samples, mean = 0.35, sd = 0.02)),
#   `5` = abs(rnorm(n = n_chains*n_samples, mean = 0.05, sd = 0.01))) %>%
#   as.matrix()
# delta_j_min <- delta_min
# delta_j_max <- delta_max
# mm_nodes_mean <- rnorm(n = nrow(nodes), mean = 0.1, sd = 0.5)
# mm_nodes <- matrix(data = NA, nrow = n_samples*n_chains, ncol = nrow(nodes))
# for(i in 1:ncol(mm_nodes)){
#   mm_nodes[,i] <- rnorm(n = nrow(mm_nodes), mean = mm_nodes_mean[i], sd = 0.1)
# }
# n_samples <- 50

## create predictive data frame
pred <- sim
# nodes$node_test_1 <- 1:nrow(nodes) ; nodes$node_test_2 <- 1:nrow(nodes) ; nodes$node_1 <- nodes$node ; nodes$node_2 <- nodes$node
# pred <- pred %>%
#   left_join(nodes[,c('node_1','node_test_1')], by = 'node_1') %>%
#   left_join(nodes[,c('node_2','node_test_2')], by = 'node_2') %>%
#   mutate(age_cat_min = NA, age_cat_max = NA)
# for(i in 1:nrow(pred)){
#   pred$age_cat_min[i] <- min(pred$age_num_1[i], pred$age_num_2[i])
#   pred$age_cat_max[i] <- max(pred$age_num_1[i], pred$age_num_2[i])
# }

str(mm_nodes)
dim(mm_nodes)

## calculate mean predictions
pred$dyad_rank <- 1:nrow(sim)
pred_mu <- matrix(NA, nrow = n_samples*n_chains, ncol = nrow(sim),
                  dimnames = list(1:(n_samples*n_chains),
                                  pred$dyad_rank))
for(dyad in 1:ncol(pred_mu)){
  node_effects <- mm_nodes[,pred$node_1[dyad]] + mm_nodes[,pred$node_2[dyad]]
  # node_effects <- mm_nodes[, pred$node_test_1[dyad]] + mm_nodes[, pred$node_test_2[dyad]]
  for(draw in 1:nrow(pred_mu)){
    pred_mu[draw, dyad] <- intercept[draw] + b_min[draw]*sum(delta_j_min[draw,(1:pred$age_cat_min[dyad])]) + b_max[draw]*sum(delta_j_max[draw,(1:pred$age_cat_max[dyad])]) + node_effects[draw]
  }
}
pred$pred_mu <- apply(pred_mu, 2, mean)
pred$pred_mu_lwr <- apply(pred_mu, 2, quantile, prob = 0.025)
pred$pred_mu_upr <- apply(pred_mu, 2, quantile, prob = 0.975)

## calculate full prediction distribution
pred_full <- pred_mu
for(dyad in 1:ncol(pred_mu)){
  for(draw in 1:nrow(pred_mu)){
    pred_full[draw, dyad] <- MASS::mvrnorm(n = 1, mu = pred_mu[draw, dyad], Sigma = sigma[draw])
  }
}
pred$pred_full_lwr <- apply(pred_full, 2, quantile, prob = 0.025)
pred$pred_full_upr <- apply(pred_full, 2, quantile, prob = 0.975)

## convert to long format
pred_full_df <- pred_full %>% 
  as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = 'dyad_id',
               values_to = 'prediction') %>% 
  mutate(dyad_id = as.integer(dyad_id)) %>% 
  left_join(pred, by = 'dyad_id')
pred_full_df <- pred_full_df %>% 
  rename(age_num_min = age_cat_min,
         age_num_max = age_cat_max) %>% 
  mutate(age_cat_min = ifelse(age_num_min == 1, '10-15 yrs',
                              ifelse(age_num_min == 2, '16-20 yrs',
                                     ifelse(age_num_min == 3, '21-25 yrs',
                                            ifelse(age_num_min == 4, '26-40 yrs',
                                                   '>40 yrs')))),
         age_cat_max = ifelse(age_num_max == 1, '10-15 yrs',
                              ifelse(age_num_max == 2, '16-20 yrs',
                                     ifelse(age_num_max == 3, '21-25 yrs',
                                            ifelse(age_num_max == 4, '26-40 yrs',
                                                   '>40 yrs'))))) %>% 
  mutate(age_cat_max = factor(age_cat_max,
                              levels = c('10-15 yrs','16-20 yrs',
                                         '21-25 yrs','26-40 yrs',
                                         '>40 yrs')),
         age_cat_max = factor(age_cat_max,
                              levels = c('10-15 yrs','16-20 yrs',
                                         '21-25 yrs','26-40 yrs',
                                         '>40 yrs')))

## sort age categories in pred data frame too
pred <- pred %>% 
  rename(age_num_min = age_cat_min,
         age_num_max = age_cat_max) %>% 
  mutate(age_cat_min = ifelse(age_num_min == 1, '10-15 yrs',
                              ifelse(age_num_min == 2, '16-20 yrs',
                                     ifelse(age_num_min == 3, '21-25 yrs',
                                            ifelse(age_num_min == 4, '26-40 yrs',
                                                   '>40 yrs')))),
         age_cat_max = ifelse(age_num_max == 1, '10-15 yrs',
                              ifelse(age_num_max == 2, '16-20 yrs',
                                     ifelse(age_num_max == 3, '21-25 yrs',
                                            ifelse(age_num_max == 4, '26-40 yrs',
                                                   '>40 yrs'))))) %>% 
  mutate(age_cat_max = factor(age_cat_max,
                              levels = c('10-15 yrs','16-20 yrs',
                                         '21-25 yrs','26-40 yrs',
                                         '>40 yrs')),
         age_cat_max = factor(age_cat_max,
                              levels = c('10-15 yrs','16-20 yrs',
                                         '21-25 yrs','26-40 yrs',
                                         '>40 yrs')))

## plot
ggplot()+
  geom_violin(data = pred_full_df,
              mapping = aes(x = age_cat_min,
                            y = prediction,
                            fill = age_cat_max),
              alpha = 0.4)+
  geom_boxplot(data = pred,
               mapping = aes(x = age_cat_min,
                             y = pred_mu,
                             fill = age_cat_max),
               width = 0.2,
               position = position_dodge(0.9))+
  geom_point(data = sim,
             mapping = aes(x = age_cat_min,
                           y = mu_std,
                           fill = age_cat_max,
                           size = total_sightings,
                           group = age_cat_max),
             pch = 21,
             position = position_dodge(0.9))+
  scale_fill_viridis_d()+
  labs(x = 'age category of younger elephant',
       y = 'logit edge weight',
       fill = 'age category of older elephant')

ggplot()+
  geom_violin(data = pred_full_df,
              mapping = aes(x = age_cat_min,
                            y = invlogit(prediction),
                            fill = age_cat_max),
              alpha = 0.4)+
  geom_boxplot(data = pred,
               mapping = aes(x = age_cat_min,
                             y = invlogit(pred_mu),
                             fill = age_cat_max),
               width = 0.2,
               position = position_dodge(0.9))+
  geom_point(data = sim,
             mapping = aes(x = age_cat_min,
                           y = invlogit(mu_std),
                           fill = age_cat_max,
                           size = total_sightings,
                           group = age_cat_max),
             pch = 21,
             position = position_dodge(0.9))+
  scale_fill_viridis_d()+
  labs(x = 'age category of younger elephant',
       y = 'logit edge weight',
       fill = 'age category of older elephant')







# ## convert to unstandardised scale
# pred$pred_unstd <- pred$pred*sd(sim$mu) + mean(sim$mu)
# 
# ## summarise
# pred_summary_all <- pred %>%
#   group_by(age_min, age_max) %>%
#   mutate(pred_lwr = quantile(pred, probs = 0.025),
#          pred_mean = mean(pred),
#          pred_upr = quantile(pred, probs = 0.975),
#          pred_unstd_lwr = quantile(pred_unstd, probs = 0.025),
#          pred_unstd_mean = mean(pred_unstd),
#          pred_unstd_upr = quantile(pred_unstd, probs = 0.975)) %>%
#   ungroup() %>%
#   select(age_min, age_max,
#          pred_lwr, pred_mean, pred_upr,
#          pred_unstd_lwr, pred_unstd_mean, pred_unstd_upr) %>%
#   distinct()
# pred_summary_possible <- pred %>%
#   filter(age_min <= age_max) %>%
#   group_by(age_min, age_max) %>%
#   mutate(pred_lwr = quantile(pred, probs = 0.025),
#          pred_mean = mean(pred),
#          pred_upr = quantile(pred, probs = 0.975),
#          pred_unstd_lwr = quantile(pred_unstd, probs = 0.025),
#          pred_unstd_mean = mean(pred_unstd),
#          pred_unstd_upr = quantile(pred_unstd, probs = 0.975)) %>%
#   ungroup() %>%
#   select(age_min, age_max,
#          pred_lwr, pred_mean, pred_upr,
#          pred_unstd_lwr, pred_unstd_mean, pred_unstd_upr) %>%
#   distinct()

dev.off()
save.image('step5_dyadicregression/motnp_dyadic_simulation.RData')

## print progress marker
print('predictions complete')

##################### GO FROM HERE ONCE IT'S RERUN ####################################
# ## extract slope estimates -- this is definitely wrong but I can't work out how to do it properly
# contrasts <- data.frame(intercept = intercept, b_min = b_min, b_max = b_max) %>% 
#   mutate(diff_1_2_min = NA, diff_2_3_min = NA, diff_3_4_min = NA, diff_4_5_min = NA,
#          diff_1_2_max = NA, diff_2_3_max = NA, diff_3_4_max = NA, diff_4_5_max = NA)
# for(i in 1:nrow(contrasts)){
#   x <- pred %>%
#     filter(intcp == contrasts$intercept[i] & 
#              beta_min == contrasts$b_min[i] &
#              beta_max == contrasts$b_max[i])
#   a <- x %>% filter(age_min == 1) ; b <- x %>% filter(age_min == 2) ; c <- x %>% filter(age_min == 3) ; d <- x %>% filter(age_min == 4) ; e <- x %>% filter(age_min == 5)
#   contrasts$diff_1_2_min[i] <- mean(a$pred_unstd - b$pred_unstd)
#   contrasts$diff_2_3_min[i] <- mean(a$pred_unstd - b$pred_unstd)
#   contrasts$diff_3_4_min[i] <- mean(a$pred_unstd - b$pred_unstd)
#   contrasts$diff_4_5_min[i] <- mean(a$pred_unstd - b$pred_unstd)
#   a <- x %>% filter(age_max == 1) ; b <- x %>% filter(age_max == 2) ; c <- x %>% filter(age_max == 3) ; d <- x %>% filter(age_max == 4) ; e <- x %>% filter(age_max == 5)
#   contrasts$diff_1_2_max[i] <- mean(a$pred_unstd - b$pred_unstd)
#   contrasts$diff_2_3_max[i] <- mean(a$pred_unstd - b$pred_unstd)
#   contrasts$diff_3_4_max[i] <- mean(a$pred_unstd - b$pred_unstd)
#   contrasts$diff_4_5_max[i] <- mean(a$pred_unstd - b$pred_unstd)
# }
# rm(a,b,c,d,e,x) ; gc()
# contrasts$diff_1_2_min_year <- contrasts$diff_1_2_min/(16-10)
# contrasts$diff_1_2_max_year <- contrasts$diff_1_2_max/(16-10)
# contrasts$diff_2_3_min_year <- contrasts$diff_2_3_min/(21-16)
# contrasts$diff_2_3_max_year <- contrasts$diff_2_3_max/(21-16)
# contrasts$diff_3_4_min_year <- contrasts$diff_3_4_min/(26-21)
# contrasts$diff_3_4_max_year <- contrasts$diff_3_4_max/(26-21)
# contrasts$diff_4_5_min_year <- contrasts$diff_4_5_min/(41-26)
# contrasts$diff_4_5_max_year <- contrasts$diff_4_5_max/(41-26)
# 
# sim_min  # original (true) slope value
# mean(c(contrasts$diff_1_2_min_year,contrasts$diff_2_3_min_year,contrasts$diff_3_4_min_year, contrasts$diff_4_5_min_year))
# 
# sim_max  # original (true) slope value
# mean(c(contrasts$diff_1_2_max_year,contrasts$diff_2_3_max_year,contrasts$diff_3_4_max_year, contrasts$diff_4_5_max_year))
# 
# ## summarise predictions
# sum_pred <- pred %>% 
#   group_by(age_min, age_max) %>% 
#   mutate(mean_std = mean(pred),
#          mean_unstd = mean(pred_unstd)) %>% 
#   select(age_min, age_max, mean_std, mean_unstd) %>% 
#   distinct()
# ggplot(sum_pred)+
#   geom_point(aes(x = age_min, y = mean_std, colour = as.factor(age_max)))
# 
# 
# 
# 
# 
# 
# 
# n = 4
# (mat <- matrix(data = invlogit(rnorm(n = 10*n, mean = -3, sd = 5)),
#               ncol = n, byrow = F))
# cov(mat)
# 
# for(i in 1:ncol(mat)){
#   for(j in 2:(ncol(mat)-1)){
#     print(cov(mat[,c(i,j)]))
#   }
# }
# 
# cov <- matrix(NA, ncol = n, nrow = n)
# for(i in 1:nrow(cov)){
#   for(j in 1:ncol(cov)){
#     section <- cov(mat[,c(i,j)])
#     if(i < j){
#       cov[i,j] <- section[]
#     }
#   }
# }
# 
# 
# 
# 
# 
# 
