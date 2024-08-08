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
library(patchwork, lib.loc = '../packages/')

#set_cmdstan_path('R:/rsrch/df525/phd/hkm513/packages/.cmdstan2/cmdstan-2.33.1/')
theme_set(theme_bw())

pdf('step5_dyadicregression/simulate_dyadic_sd_motnp_prep.pdf')

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
rm(sim_dat_long) ; gc()

## print progress marker
print('edges plotted')

#### fit Gaussian distribution to output of edge weight model ####
## fit a multivariate normal dist to the edges -- quantities will be given to Stan model as data to model joint posteriors of edge weight in the regression
logit_edge_draws_mu <- apply(sim_dat, 2, mean)
logit_edge_draws_sd <- apply(sim_dat, 2, sd)

## randomly select samples to examine
num_check <- 20
selected_samples <- sample(1:(n_dyads-1), num_check, replace = FALSE)

## set grid layout
rows <- floor(sqrt(num_check))
cols <- ceiling(num_check / rows)
par(mfrow=c(rows, cols), mar=c(2,2,2,1))

## plot normal approximation
for (i in selected_samples) {
  fitted_values <- rnorm(1e5,
                         mean = logit_edge_draws_mu[i],
                         sd = logit_edge_draws_sd[i])

  hist(sim_dat[,i], probability = TRUE, las = 1,
       main = paste("Dyad", i), xlab = "Value", breaks = 50)
  lines(density(fitted_values), col="blue", lw=1.5)
}

## plot covariance (though we aren't using that now)
for (i in selected_samples) {
  fitted_values <- rnorm(1e5,
                         mean = logit_edge_draws_mu[i],
                         sd = logit_edge_draws_sd[i])

  plot(unlist(sim_dat[,i]), unlist(sim_dat[,i+1]), col = rgb(0,0,1,0.05), las = 1,
       main = paste("cov ", i ,"&",i+1))
}

## reset plot window and clean up
par(mfrow=c(1,1), mai = c(1,1,0.5,0.5))
rm(cols, fitted_values, i, j, num_check, rows, selected_samples) ; gc()

## print progress marker
print('normal approximation complete')
dev.off()

##################### Minimum age #####################
pdf('step5_dyadicregression/simulate_dyadic_sd_motnp_checkoutputs_min.pdf')

#### prior predictive check ####
n <- 100
n_age_cat <- length(unique(sim$age_cat_max))
beta_age_min <- rnorm(n, 0, 2)
intercept <- rnorm(n, 0, 2)
min_dirichlet <- rdirichlet(n, rep(1, n_age_cat))
mm_nodes <- rnorm(n, 0, rexp(1, 1))
plot(NULL, las = 1, xlab = 'minimum age (standardised)', ylab = 'logit edge weight (standardised)',
     ylim = c(min(logit_edge_draws_mu)-5, max(logit_edge_draws_mu)+5),
     xlim = c(1,n_age_cat))
abline(h = min(logit_edge_draws_mu), lty = 2) ; abline(h = max(logit_edge_draws_mu), lty = 2)
x <- 1:n_age_cat
for(i in 1:n){
  y <- rep(NA, length(x))
  for(j in 1:length(x)){
    y[j] <- intercept[i] + beta_age_min[i]*sum(min_dirichlet[i,][1:x[j]]) + mm_nodes[i]
  }
  lines(x = x, y = y, col = rgb(0,0,1,0.4))
}
rm(n, beta_age_min, intercept, min_dirichlet, i, y, j, x, mm_nodes) ; gc()

## print progress marker
print('prior predictive complete')

#### fit dyadic regression ####
## create data list
dyad_data <- list(
  num_dyads = n_dyads,                      # number of dyads
  num_nodes = n_nodes,                      # number of nodes
  num_age_cat = n_age_cat,                  # number of unique age categories
  length_dirichlet = n_age_cat + 1,         # number of unique age categories + 1
  logit_edge_mu = logit_edge_draws_mu,      # means of the logit edge weights
  logit_edge_sd = logit_edge_draws_sd,      # standard deviation of logit edge weights
  age_min_cat = sim$age_cat_min,            # age of younger dyad member
  node_1 = sim$node_rand_1,                 # node IDs for multimembership effects
  node_2 = sim$node_rand_2,                 # node IDs for multimembership effects
  prior_min = rep(1, n_age_cat)             # prior for minimum age slope
)

## print progress marker
print('data list created')

## save output
save.image('step5_dyadicregression/motnp_dyadic_simulation_agemin.RData')

## load dyadic regression model
dyadic_regression <- stan_model('models/dyadic_regression_motnp_min.stan')
#dyadic_regression <- cmdstan_model('models/dyadic_regression_motnp_min.stan')

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
save.image('step5_dyadicregression/motnp_dyadic_simulation_agemin.RData')

## print progress marker
print('model run')

#### check outputs ####
#load('step5_dyadicregression/motnp_dyadic_simulation_agemin.RData')

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

rownames(summary)[summary$n_eff < 1000]
summary[summary$n_eff < 1000,c(1,3,9,10)]

save.image('step5_dyadicregression/motnp_dyadic_simulation_agemin.RData')

## extract dyadic regression slopes
b_min <- posterior_samples$beta_age_min
intercept <- posterior_samples$intercept
sigma <- posterior_samples$sigma
delta_min <- posterior_samples[,c('delta_min[1]','delta_min[2]','delta_min[3]','delta_min[4]','delta_min[5]')] ; colnames(delta_min) <- 1:n_age_cat
delta_j_min <- posterior_samples[,c('delta_j_min[1]','delta_j_min[2]','delta_j_min[3]','delta_j_min[4]','delta_j_min[5]','delta_j_min[6]')] ; colnames(delta_j_min) <- 1:(n_age_cat+1)
mm_nodes <- posterior_samples[,grep(pattern = 'mm_nodes', x = colnames(posterior_samples))]
parameters <- data.frame(beta_age_min = b_min,
                         intercept = intercept,
                         sigma = sigma) %>%
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>%
  pivot_longer(cols = c('beta_age_min','sigma','intercept'),
               names_to = 'parameter', values_to = 'slope_draw')

## traceplots
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

#### posterior predictive check ####
plot(density(as.numeric(sim_dat[1, ])),
     main = "Posterior predictive density of edge weights:\nblack = measured edge, red = predicted edge",
     ylim = c(0, 1), col = rgb(0, 0, 0, 0.25), las = 1)
for (i in 1:100) {
  j <- sample(1:1000, 1)

  mu_plot <- rep(NA, n_dyads)
  for(k in 1:n_dyads){
    mu_plot[k] <- intercept[j] + b_min[j]*sum(delta_j_min[j,(1:dyad_data$age_min_cat[k])])
  }

  sigma_plot <- dyad_data$logit_edge_sd + rep(sigma[j], n_dyads)

  norm <- rep(NA, n_dyads)
  for(k in 1:n_dyads){
    norm[k] <- rnorm(1, mu_plot[k], sigma_plot[k])
  }

  lines(density(as.numeric(sim_dat[j, ])), col = rgb(0, 0, 0, 0.25)) # black lines for edge samples
  lines(density(norm), col = rgb(1, 0, 0, 0.25))                  # red lines for predictions

}
rm(norm) ; gc()

## save output
save.image('step5_dyadicregression/motnp_dyadic_simulation_agemin.RData')
dev.off()

## print progress marker
print('outputs checked')

#### predict from raw data ####
pdf('step5_dyadicregression/simulate_dyadic_sd_motnp_predictions_agemin.pdf')
#load('step5_dyadicregression/motnp_dyadic_simulation_agemin.RData')

## create predictive data frame
pred <- sim
nodes$node_test_1 <- 1:nrow(nodes) ; nodes$node_test_2 <- 1:nrow(nodes) ; nodes$node_1 <- nodes$node ; nodes$node_2 <- nodes$node
pred <- pred %>%
  left_join(nodes[,c('node_1','node_test_1')], by = 'node_1') %>%
  left_join(nodes[,c('node_2','node_test_2')], by = 'node_2')

str(mm_nodes)
dim(mm_nodes)

## calculate mean predictions
pred$dyad_rank <- 1:nrow(pred)
pred_mu <- matrix(NA, nrow = n_samples*n_chains, ncol = nrow(sim),
                  dimnames = list(1:(n_samples*n_chains),
                                  pred$dyad_rank))
for(dyad in 1:ncol(pred_mu)){
  # node_effects <- mm_nodes[,paste0('mm_nodes[',pred$node_1_males[dyad],']')] + mm_nodes[,paste0('mm_nodes[',pred$node_2_males[dyad],']')]
  node_effects <- mm_nodes[,paste0('mm_nodes[',pred$node_test_1[dyad],']')] + mm_nodes[,paste0('mm_nodes[',pred$node_test_2[dyad],']')]
  for(draw in 1:nrow(pred_mu)){
    pred_mu[draw, dyad] <- intercept[draw] + b_min[draw]*sum(delta_j_min[draw,(1:pred$age_cat_min[dyad])]) + node_effects[draw]
  }
  if(dyad %% 1000 == 0){ print(dyad) }
}
pred$pred_mu <- apply(pred_mu, 2, mean)
pred$pred_mu_lwr <- apply(pred_mu, 2, quantile, prob = 0.025)
pred$pred_mu_upr <- apply(pred_mu, 2, quantile, prob = 0.975)

save.image('step5_dyadicregression/motnp_dyadic_simulation_agemin.RData')

## plot predicted vs original
ggplot(pred)+
  geom_point(aes(x = mu_std, y = pred_mu,
                 colour = as.factor(age_cat_min)))+
  geom_abline(slope = 1, intercept = 0)+
  scale_colour_viridis_d()

## calculate full prediction distribution
pred_full <- pred_mu
for(dyad in 1:ncol(pred_mu)){
  for(draw in 1:nrow(pred_mu)){
    pred_full[draw, dyad] <- rnorm(n = 1, mean = pred_mu[draw, dyad], sd = posterior_samples$sigma[draw])
  }
  if(dyad %% 1000 == 0){ print(dyad) }
}
save.image('step5_dyadicregression/motnp_dyadic_simulation_agemin.RData')
pred$pred_full_lwr <- apply(pred_full, 2, quantile, prob = 0.025)
pred$pred_full_upr <- apply(pred_full, 2, quantile, prob = 0.975)

## convert to long format
colnames(pred_full) <- pred$dyad_id
pred_full_df <- pred_full %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = 'dyad_id',
               values_to = 'prediction') %>%
  mutate(dyad_id = as.integer(dyad_id)) %>%
  left_join(pred, by = 'dyad_id')
pred_full_df <- pred_full_df %>%
  rename(age_num_min = age_cat_min) %>%
  mutate(age_cat_min = ifelse(age_num_min == 1, '10-15 yrs',
                              ifelse(age_num_min == 2, '15-19 yrs',
                                     ifelse(age_num_min == 3, '21-25 yrs',
                                            ifelse(age_num_min == 4, '26-40 yrs',
                                                   '>40 yrs'))))) %>%
  mutate(age_cat_min = factor(age_cat_min,
                              levels = c('10-15 yrs','15-19 yrs',
                                         '21-25 yrs','26-40 yrs',
                                         '>40 yrs')))

## sort age categories in pred data frame too
pred <- pred %>%
  rename(age_num_min = age_cat_min) %>%
  mutate(age_cat_min = ifelse(age_num_min == 1, '10-15 yrs',
                              ifelse(age_num_min == 2, '15-19 yrs',
                                     ifelse(age_num_min == 3, '21-25 yrs',
                                            ifelse(age_num_min == 4, '26-40 yrs',
                                                   '>40 yrs'))))) %>%
  mutate(age_cat_min = factor(age_cat_min,
                              levels = c('10-15 yrs','15-19 yrs',
                                         '21-25 yrs','26-40 yrs',
                                         '>40 yrs')))

## ...and sim data frame...
sim <- sim %>%
  rename(age_num_min = age_cat_min) %>%
  mutate(age_cat_min = ifelse(age_num_min == 1, '10-15 yrs',
                              ifelse(age_num_min == 2, '15-19 yrs',
                                     ifelse(age_num_min == 3, '21-25 yrs',
                                            ifelse(age_num_min == 4, '26-40 yrs',
                                                   '>40 yrs'))))) %>%
  mutate(age_cat_min = factor(age_cat_min,
                              levels = c('10-15 yrs','15-19 yrs',
                                         '21-25 yrs','26-40 yrs',
                                         '>40 yrs')))

## plot
ggplot()+
  geom_violin(data = pred_full_df,
              mapping = aes(x = age_cat_min,
                            y = prediction),
              alpha = 0.4)+
  geom_boxplot(data = pred,
               mapping = aes(x = age_cat_min,
                             y = pred_mu),
               width = 0.2,
               position = position_dodge(0.9))+
  geom_point(data = sim,
             mapping = aes(x = age_cat_min,
                           y = mu_std,
                           size = total_sightings),
             pch = 21,
             position = position_dodge(0.9))+
  scale_fill_viridis_d()+
  labs(x = 'age category of younger elephant',
       y = 'logit edge weight')

ggplot()+
  geom_violin(data = pred_full_df,
              mapping = aes(x = age_cat_min,
                            y = invlogit(prediction)),
              alpha = 0.4)+
  geom_boxplot(data = pred,
               mapping = aes(x = age_cat_min,
                             y = invlogit(pred_mu)),
               width = 0.2,
               position = position_dodge(0.9))+
  geom_point(data = sim,
             mapping = aes(x = age_cat_min,
                           y = invlogit(mu_std),
                           size = total_sightings),
             pch = 21,
             position = position_dodge(0.9))+
  scale_fill_viridis_d()+
  labs(x = 'age category of younger elephant',
       y = 'logit edge weight')

## convert to unstandardised scale
pred_full_df$pred_unstd <- pred_full_df$prediction*sd(sim$mu) + mean(sim$mu)
pred$pred_mu_unstd <- pred$pred_mu*sd(sim$mu) + mean(sim$mu)
pred$pred_mulwr_unstd <- pred$pred_mu_lwr*sd(sim$mu) + mean(sim$mu)
pred$pred_muupr_unstd <- pred$pred_mu_upr*sd(sim$mu) + mean(sim$mu)

## summarise
pred_summary_all <- pred_full_df %>%
  group_by(age_cat_min) %>%
  mutate(pred_lwr = quantile(prediction, probs = 0.025),
         pred_mean = mean(prediction),
         pred_upr = quantile(prediction, probs = 0.975),
         pred_unstd_lwr = quantile(pred_unstd, probs = 0.025),
         pred_unstd_mean = mean(pred_unstd),
         pred_unstd_upr = quantile(pred_unstd, probs = 0.975)) %>%
  ungroup() %>%
  select(age_cat_min,
         pred_lwr, pred_mean, pred_upr,
         pred_unstd_lwr, pred_unstd_mean, pred_unstd_upr) %>%
  distinct()

ggplot()+
  geom_violin(data = pred_full_df,
              mapping = aes(x = age_cat_min,
                            y = prediction),
              alpha = 0.4)+
  geom_boxplot(data = sim,
               mapping = aes(x = age_cat_min,
                             y = mu_std),
               width = 0.2,
               position = position_dodge(0.9))+
  # geom_point(data = pred_summary_all,
  #            mapping = aes(x = age_cat_min,
  #                          y = pred_mean),
  #            pch = 21,
  #            size = 2,
  #            fill = 'white',
  #            position = position_dodge(0.9))+
  scale_fill_viridis_d()+
  labs(x = 'age category of younger elephant',
       y = 'logit edge weight')

dev.off()
save.image('step5_dyadicregression/motnp_dyadic_simulation_agemin.RData')

## print progress marker
print('predictions complete')

#### extract original slope values ####
load('step5_dyadicregression/motnp_dyadic_simulation_agemin.RData')
pdf('step5_dyadicregression/simulate_dyadic_sd_motnp_contrasts_agemin.pdf')

## calculate predictions from altered data frames
pred_new <- pred %>%
  rename(age_cat_1_org = age_cat_1,
         age_cat_2_org = age_cat_2,
         age_num_1_org = age_num_1,
         age_num_2_org = age_num_2,
         age_cat_min_org = age_cat_min,
         age_num_min_org = age_num_min)

## create function for predicting mean
get_mean_predictions <- function(prediction_df){
  prediction_df$dyad_rank <- 1:nrow(prediction_df)
  mu <- matrix(NA, nrow = n_samples*n_chains, ncol = nrow(sim),
                    dimnames = list(1:(n_samples*n_chains),
                                    prediction_df$dyad_rank))
  for(dyad in 1:ncol(mu)){
    # node_effects <- mm_nodes[,paste0('mm_nodes[',prediction_df$node_1_males[dyad],']')] + mm_nodes[,paste0('mm_nodes[',prediction_df$node_2_males[dyad],']')]
    node_effects <- mm_nodes[,paste0('mm_nodes[',prediction_df$node_test_1[dyad],']')] + mm_nodes[,paste0('mm_nodes[',prediction_df$node_test_2[dyad],']')]
    for(draw in 1:nrow(mu)){
      mu[draw, dyad] <- intercept[draw] + b_min[draw]*sum(delta_j_min[draw,(1:prediction_df$age_num_min[dyad])])
    }
  }
  return(mu)
  }

## create function for predicting full distribution
get_full_predictions <- function(mu_matrix, sigma){
  predictions <- mu_matrix
  for(dyad in 1:ncol(predictions)){
    for(draw in 1:nrow(predictions)){
      predictions[draw, dyad] <- rnorm(n = 1,
                                       mean = mu_matrix[draw, dyad],
                                       sd = sigma[draw])
    }
  }
  return(predictions)
}

## add progress marker
print('functions created')

## calculate predictions for each age combination
predictions_mean <- list()
predictions_full <- list()
for(i in c(1:5)){
  ## create data frame to predict from
  pred_new_age <- pred_new %>%
    mutate(age_cat_min = ifelse(i == 1, '10-15',
                                ifelse(i == 2, '15-19',
                                       ifelse(i == 3, '20-25',
                                              ifelse(i == 4, '25-40',
                                                     '40+')))),
           age_num_min = i)
  
  ## extract mean predictions
  predictions_mean[[i]] <- get_mean_predictions(prediction_df = pred_new_age)

  ## extract full predictions
  predictions_full[[i]] <- get_full_predictions(mu_matrix = predictions_mean[[i]],
                                                sigma = posterior_samples$sigma )

  ## save outputs in case R crashes
  saveRDS(object = as.data.frame(predictions_mean[[i]]),
          file = paste0('../data_processed/step5_dyadicregression/simulate_motnp_meanpredictions_minage',i,'.RDS'))
  saveRDS(object = as.data.frame(predictions_full[[i]]),
          file = paste0('../data_processed/step5_dyadicregression/simulate_motnp_fullpredictions_minage',i,'.RDS'))
}

# ## if image does not write properly, will need to rebuild lists from written predictions
# predictions_mean <- list()
# predictions_full <- list()
# for(i in c(1:5)){
#   ## save outputs in case R crashes
#   mean <- readRDS(file = paste0('../data_processed/step5_dyadicregression/simulate_motnp_meanpredictions_agecombo',i,'.RDS'))
#   full <- readRDS(file = paste0('../data_processed/step5_dyadicregression/simulate_motnp_fullpredictions_agecombo',i,'.RDS'))
# 
#   ## extract mean predictions
#   predictions_mean[[i]] <- mean
# 
#   ## extract full predictions
#   predictions_full[[i]] <- full
#   }

## save workspace
save.image('step5_dyadicregression/motnp_dyadic_simulation_agemin.RData')

## add progress marker
print('predictions complete')

## calculate contrasts -- age category
contrasts_1vs2 <- predictions_full[[2]] - predictions_full[[1]]
print(paste0('contrast age 1 becomes age 2: mean = ', mean(contrasts_1vs2),
             ', stdv = ', sd(contrasts_1vs2),
             ', 95% CI = [', quantile(contrasts_1vs2, prob = 0.025),':',
             quantile(contrasts_1vs2, prob = 0.975),']'))
contrasts_2vs3 <- predictions_full[[3]] - predictions_full[[2]]
print(paste0('contrast age 2 becomes age 3: mean = ', mean(contrasts_2vs3),
             ', stdv = ', sd(contrasts_2vs3),
             ', 95% CI = [', quantile(contrasts_2vs3, prob = 0.025),':',
             quantile(contrasts_2vs3, prob = 0.975),']'))
contrasts_3vs4 <- predictions_full[[4]] - predictions_full[[3]]
print(paste0('contrast age 3 becomes age 4: mean = ', mean(contrasts_3vs4),
             ', stdv = ', sd(contrasts_3vs4),
             ', 95% CI = [', quantile(contrasts_3vs4, prob = 0.025),':',
             quantile(contrasts_3vs4, prob = 0.975),']'))
contrasts_4vs5 <- predictions_full[[5]] - predictions_full[[4]]
print(paste0('contrast age 4 becomes age 5: mean = ', mean(contrasts_4vs5),
             ', stdv = ', sd(contrasts_4vs5),
             ', 95% CI = [', quantile(contrasts_4vs5, prob = 0.025),':',
             quantile(contrasts_4vs5, prob = 0.975),']'))

## add progress marker
print('contrasts calculated per age category')

## create data frame to summarise contrasts
contrasts <- expand.grid(org_min = c(1:5),
                         alt_min = c(1:5)) %>% 
  filter(org_min < alt_min) %>% 
  mutate(min_change = alt_min - org_min) %>% 
  mutate(org_min_lwr = ifelse(org_min == 1, 10,
                              ifelse(org_min == 2, 15,
                                     ifelse(org_min == 3, 20,
                                            ifelse(org_min == 4, 25,
                                                   40)))),
         alt_min_lwr = ifelse(alt_min == 1, 10,
                              ifelse(alt_min == 2, 15,
                                     ifelse(alt_min == 3, 20,
                                            ifelse(alt_min == 4, 25,
                                                   40))))) %>% 
  mutate(mean_contrast = NA,
         stdv_contrast = NA,
         lwr_contrast = NA,
         upr_contrast = NA,
         num_years_change = NA,
         mean_per_year = NA)

## add progress marker
print('contrasts data frame created')

for(i in 1:nrow(contrasts)){
  contrast_matrix <- predictions_full[[as.numeric(as.character(contrasts$alt_min[i]))]] - predictions_full[[as.numeric(as.character(contrasts$org_min[i]))]]
  contrast_matrix <- as.matrix(contrast_matrix)
  contrasts$mean_contrast[i] <- mean(contrast_matrix)
  contrasts$stdv_contrast[i] <- sd(contrast_matrix)
  contrasts$lwr_contrast[i]  <- quantile(contrast_matrix, prob = 0.025)
  contrasts$upr_contrast[i]  <- quantile(contrast_matrix, prob = 0.975)
  contrasts$num_years_change[i] <- ifelse(contrasts$min_change[i] > 0, 
                                          ifelse(contrasts$max_change[i] > 0, NA,
                                                 contrasts$alt_min_lwr[i] - contrasts$org_min_lwr[i]),
                                          contrasts$alt_max_lwr[i] - contrasts$org_max_lwr[i])
  contrasts$mean_per_year[i] <- mean(contrast_matrix) / contrasts$num_years_change[i]
}

## add progress marker
print('contrasts data frame populated')

## calculate contrats -- year-by-year
print(paste0('original (true) slope value: ', sim_min))
print(paste0('contrast per year: mean = ', mean(contrasts$mean_per_year),
             ', stdv = ', sd(contrasts$mean_per_year),
             ', 95% CI = [', quantile(contrasts$mean_per_year, prob = 0.025),':',
             quantile(contrasts$mean_per_year, prob = 0.975),']'))

## add progress marker
print('contrasts calculated per year')

## plot
contrasts <- contrasts %>% 
  mutate(org_min_cat = ifelse(org_min == 1, '10-15',
                              ifelse(org_min == 2, '15-19',
                                     ifelse(org_min == 3, '20-25',
                                            ifelse(org_min == 4, '25-40',
                                                   '40+')))),
         alt_min_cat = ifelse(alt_min == 1, '10-15',
                              ifelse(alt_min == 2, '15-19',
                                     ifelse(alt_min == 3, '20-25',
                                            ifelse(alt_min == 4, '25-40',
                                                   '40+'))))) %>% 
  mutate(org_min_cat = factor(org_min_cat,
                              levels = c('10-15','15-19','20-25',
                                         '25-40','40+')),
         alt_min_cat = factor(alt_min_cat,
                              levels = c('10-15','15-19','20-25',
                                         '25-40','40+')))
(min_plot <- ggplot(data = contrasts)+
    geom_tile(aes(x = org_min_cat, y = alt_min_cat,
                  fill = mean_contrast))+
    scale_colour_viridis_c()+
    scale_x_discrete(drop = T)+
    scale_y_discrete(drop = T)+
    labs(x = 'original minimum age category',
         y = 'altered minimum age category',
         fill = 'contrast in\nlogit(probability\nof associating)'))

save.image('step5_dyadicregression/motnp_dyadic_simulation_agemin.RData')
dev.off()

##################### Maximum age #####################
pdf('step5_dyadicregression/simulate_dyadic_sd_motnp_checkoutputs_max.pdf')

#### prior predictive check ####
n <- 100
n_age_cat <- length(unique(sim$age_cat_max))
beta_age_max <- rnorm(n, 0, 2)
intercept <- rnorm(n, 0, 2)
max_dirichlet <- rdirichlet(n, rep(1, n_age_cat))
mm_nodes <- rnorm(n, 0, rexp(1, 1))
plot(NULL, las = 1, xlab = 'maximum age (standardised)', ylab = 'logit edge weight (standardised)',
     ylim = c(min(logit_edge_draws_mu)-5, max(logit_edge_draws_mu)+5),
     xlim = c(1,n_age_cat))
abline(h = max(logit_edge_draws_mu), lty = 2) ; abline(h = max(logit_edge_draws_mu), lty = 2)
x <- 1:n_age_cat
for(i in 1:n){
  y <- rep(NA, length(x))
  for(j in 1:length(x)){
    y[j] <- intercept[i] + beta_age_max[i]*sum(max_dirichlet[i,][1:x[j]]) + mm_nodes[i]
  }
  lines(x = x, y = y, col = rgb(0,0,1,0.4))
}
rm(n, beta_age_max, intercept, max_dirichlet, i, y, j, x, mm_nodes) ; gc()

## print progress marker
print('prior predictive complete')

#### fit dyadic regression ####
## create data list
dyad_data <- list(
  num_dyads = n_dyads,                      # number of dyads
  num_nodes = n_nodes,                      # number of nodes
  num_age_cat = n_age_cat,                  # number of unique age categories
  length_dirichlet = n_age_cat + 1,         # number of unique age categories + 1
  logit_edge_mu = logit_edge_draws_mu,      # means of the logit edge weights
  logit_edge_sd = logit_edge_draws_sd,      # standard deviation of logit edge weights
  age_max_cat = sim$age_cat_max,            # age of younger dyad member
  node_1 = sim$node_rand_1,                 # node IDs for multimembership effects
  node_2 = sim$node_rand_2,                 # node IDs for multimembership effects
  prior_max = rep(1, n_age_cat)             # prior for maximum age slope
)

## print progress marker
print('data list created')

## save output
save.image('step5_dyadicregression/motnp_dyadic_simulation_agemax.RData')

## load dyadic regression model
dyadic_regression <- stan_model('models/dyadic_regression_motnp_max.stan')
#dyadic_regression <- cmdstan_model('models/dyadic_regression_motnp_sd.stan')

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
save.image('step5_dyadicregression/motnp_dyadic_simulation_agemax.RData')

## print progress marker
print('model run')

#### check outputs ####
#load('step5_dyadicregression/motnp_dyadic_simulation_agemax.RData')

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

rownames(summary)[summary$n_eff < 1000]
summary[summary$n_eff < 1000,c(1,3,9,10)]

save.image('step5_dyadicregression/motnp_dyadic_simulation_agemax.RData')

## extract dyadic regression slopes
b_max <- posterior_samples$beta_age_max
intercept <- posterior_samples$intercept
sigma <- posterior_samples$sigma
delta_max <- posterior_samples[,c('delta_max[1]','delta_max[2]','delta_max[3]','delta_max[4]','delta_max[5]')] ; colnames(delta_max) <- 1:n_age_cat
delta_j_max <- posterior_samples[,c('delta_j_max[1]','delta_j_max[2]','delta_j_max[3]','delta_j_max[4]','delta_j_max[5]','delta_j_max[6]')] ; colnames(delta_j_max) <- 1:(n_age_cat+1)
mm_nodes <- posterior_samples[,grep(pattern = 'mm_nodes', x = colnames(posterior_samples))]
parameters <- data.frame(beta_age_max = b_max,
                         intercept = intercept,
                         sigma = sigma) %>%
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>%
  pivot_longer(cols = c('beta_age_max','sigma','intercept'),
               names_to = 'parameter', values_to = 'slope_draw')

## traceplots
ggplot(data = parameters)+
  geom_line(aes(x = position,
                y = slope_draw,
                colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')
delta_max %>%  as.data.frame() %>%
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
    mu_plot[k] <- intercept[j] + b_max[j]*sum(delta_j_max[j,(1:dyad_data$age_max_cat[k])])
  }
  
  sigma_plot <- dyad_data$logit_edge_sd + rep(sigma[j], n_dyads)
  
  norm <- rep(NA, n_dyads)
  for(k in 1:n_dyads){
    norm[k] <- rnorm(1, mu_plot[k], sigma_plot[k])
  }
  
  lines(density(as.numeric(sim_dat[j, ])), col = rgb(0, 0, 0, 0.25)) # black lines for edge samples
  lines(density(norm), col = rgb(1, 0, 0, 0.25))                  # red lines for predictions
  
}
rm(norm) ; gc()

## save output
save.image('step5_dyadicregression/motnp_dyadic_simulation_agemax.RData')
dev.off()

## print progress marker
print('outputs checked')

#### predict from raw data ####
pdf('step5_dyadicregression/simulate_dyadic_sd_motnp_predictions_agemax.pdf')
#load('step5_dyadicregression/motnp_dyadic_simulation_agemax.RData')

## create predictive data frame
pred <- sim
nodes$node_test_1 <- 1:nrow(nodes) ; nodes$node_test_2 <- 1:nrow(nodes) ; nodes$node_1 <- nodes$node ; nodes$node_2 <- nodes$node
pred <- pred %>%
  left_join(nodes[,c('node_1','node_test_1')], by = 'node_1') %>%
  left_join(nodes[,c('node_2','node_test_2')], by = 'node_2')

str(mm_nodes)
dim(mm_nodes)

## calculate mean predictions
pred$dyad_rank <- 1:nrow(pred)
pred_mu <- matrix(NA, nrow = n_samples*n_chains, ncol = nrow(sim),
                  dimnames = list(1:(n_samples*n_chains),
                                  pred$dyad_rank))
for(dyad in 1:ncol(pred_mu)){
  # node_effects <- mm_nodes[,paste0('mm_nodes[',pred$node_1_males[dyad],']')] + mm_nodes[,paste0('mm_nodes[',pred$node_2_males[dyad],']')]
  node_effects <- mm_nodes[,paste0('mm_nodes[',pred$node_test_1[dyad],']')] + mm_nodes[,paste0('mm_nodes[',pred$node_test_2[dyad],']')]
  for(draw in 1:nrow(pred_mu)){
    pred_mu[draw, dyad] <- intercept[draw] + b_max[draw]*sum(delta_j_max[draw,(1:pred$age_cat_max[dyad])]) + node_effects[draw]
  }
  if(dyad %% 1000 == 0){ print(dyad) }
}
pred$pred_mu <- apply(pred_mu, 2, mean)
pred$pred_mu_lwr <- apply(pred_mu, 2, quantile, prob = 0.025)
pred$pred_mu_upr <- apply(pred_mu, 2, quantile, prob = 0.975)

save.image('step5_dyadicregression/motnp_dyadic_simulation_agemax.RData')

## plot predicted vs original
ggplot(pred)+
  geom_point(aes(x = mu_std, y = pred_mu,
                 colour = as.factor(age_cat_max)))+
  geom_abline(slope = 1, intercept = 0)+
  scale_colour_viridis_d()

## calculate full prediction distribution
pred_full <- pred_mu
for(dyad in 1:ncol(pred_mu)){
  for(draw in 1:nrow(pred_mu)){
    pred_full[draw, dyad] <- rnorm(n = 1, mean = pred_mu[draw, dyad], sd = posterior_samples$sigma[draw])
  }
  if(dyad %% 1000 == 0){ print(dyad) }
}
save.image('step5_dyadicregression/motnp_dyadic_simulation_agemax.RData')
pred$pred_full_lwr <- apply(pred_full, 2, quantile, prob = 0.025)
pred$pred_full_upr <- apply(pred_full, 2, quantile, prob = 0.975)

## convert to long format
colnames(pred_full) <- pred$dyad_id
pred_full_df <- pred_full %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = 'dyad_id',
               values_to = 'prediction') %>%
  mutate(dyad_id = as.integer(dyad_id)) %>%
  left_join(pred, by = 'dyad_id')
pred_full_df <- pred_full_df %>%
  rename(age_num_max = age_cat_max) %>%
  mutate(age_cat_max = ifelse(age_num_max == 1, '10-15 yrs',
                              ifelse(age_num_max == 2, '15-19 yrs',
                                     ifelse(age_num_max == 3, '21-25 yrs',
                                            ifelse(age_num_max == 4, '26-40 yrs',
                                                   '>40 yrs'))))) %>%
  mutate(age_cat_max = factor(age_cat_max,
                              levels = c('10-15 yrs','15-19 yrs',
                                         '21-25 yrs','26-40 yrs',
                                         '>40 yrs')))

## sort age categories in pred data frame too
pred <- pred %>%
  rename(age_num_max = age_cat_max) %>%
  mutate(age_cat_max = ifelse(age_num_max == 1, '10-15 yrs',
                              ifelse(age_num_max == 2, '15-19 yrs',
                                     ifelse(age_num_max == 3, '21-25 yrs',
                                            ifelse(age_num_max == 4, '26-40 yrs',
                                                   '>40 yrs'))))) %>%
  mutate(age_cat_max = factor(age_cat_max,
                              levels = c('10-15 yrs','15-19 yrs',
                                         '21-25 yrs','26-40 yrs',
                                         '>40 yrs')))

## ...and sim data frame...
sim <- sim %>%
  rename(age_num_max = age_cat_max) %>%
  mutate(age_cat_max = ifelse(age_num_max == 1, '10-15 yrs',
                              ifelse(age_num_max == 2, '15-19 yrs',
                                     ifelse(age_num_max == 3, '21-25 yrs',
                                            ifelse(age_num_max == 4, '26-40 yrs',
                                                   '>40 yrs'))))) %>%
  mutate(age_cat_max = factor(age_cat_max,
                              levels = c('10-15 yrs','15-19 yrs',
                                         '21-25 yrs','26-40 yrs',
                                         '>40 yrs')))

## plot
ggplot()+
  geom_violin(data = pred_full_df,
              mapping = aes(x = age_cat_max,
                            y = prediction),
              alpha = 0.4)+
  geom_boxplot(data = pred,
               mapping = aes(x = age_cat_max,
                             y = pred_mu),
               width = 0.2,
               position = position_dodge(0.9))+
  geom_point(data = sim,
             mapping = aes(x = age_cat_max,
                           y = mu_std,
                           size = total_sightings),
             pch = 21,
             position = position_dodge(0.9))+
  scale_fill_viridis_d()+
  labs(x = 'age category of younger elephant',
       y = 'logit edge weight')

ggplot()+
  geom_violin(data = pred_full_df,
              mapping = aes(x = age_cat_max,
                            y = invlogit(prediction)),
              alpha = 0.4)+
  geom_boxplot(data = pred,
               mapping = aes(x = age_cat_max,
                             y = invlogit(pred_mu)),
               width = 0.2,
               position = position_dodge(0.9))+
  geom_point(data = sim,
             mapping = aes(x = age_cat_max,
                           y = invlogit(mu_std),
                           size = total_sightings),
             pch = 21,
             position = position_dodge(0.9))+
  scale_fill_viridis_d()+
  labs(x = 'age category of younger elephant',
       y = 'logit edge weight')

## convert to unstandardised scale
pred_full_df$pred_unstd <- pred_full_df$prediction*sd(sim$mu) + mean(sim$mu)
pred$pred_mu_unstd <- pred$pred_mu*sd(sim$mu) + mean(sim$mu)
pred$pred_mulwr_unstd <- pred$pred_mu_lwr*sd(sim$mu) + mean(sim$mu)
pred$pred_muupr_unstd <- pred$pred_mu_upr*sd(sim$mu) + mean(sim$mu)

## summarise
pred_summary_all <- pred_full_df %>%
  group_by(age_cat_max) %>%
  mutate(pred_lwr = quantile(prediction, probs = 0.025),
         pred_mean = mean(prediction),
         pred_upr = quantile(prediction, probs = 0.975),
         pred_unstd_lwr = quantile(pred_unstd, probs = 0.025),
         pred_unstd_mean = mean(pred_unstd),
         pred_unstd_upr = quantile(pred_unstd, probs = 0.975)) %>%
  ungroup() %>%
  select(age_cat_max,
         pred_lwr, pred_mean, pred_upr,
         pred_unstd_lwr, pred_unstd_mean, pred_unstd_upr) %>%
  distinct()

ggplot()+
  geom_violin(data = pred_full_df,
              mapping = aes(x = age_cat_max,
                            y = prediction),
              alpha = 0.4)+
  geom_boxplot(data = sim,
               mapping = aes(x = age_cat_max,
                             y = mu_std),
               width = 0.2,
               position = position_dodge(0.9))+
  # geom_point(data = pred_summary_all,
  #            mapping = aes(x = age_cat_max,
  #                          y = pred_mean),
  #            pch = 21,
  #            size = 2,
  #            fill = 'white',
  #            position = position_dodge(0.9))+
  scale_fill_viridis_d()+
  labs(x = 'age category of younger elephant',
       y = 'logit edge weight')

dev.off()
save.image('step5_dyadicregression/motnp_dyadic_simulation_agemax.RData')

## print progress marker
print('predictions complete')

#### extract original slope values ####
load('step5_dyadicregression/motnp_dyadic_simulation_agemax.RData')
pdf('step5_dyadicregression/simulate_dyadic_sd_motnp_contrasts_agemax.pdf')

## calculate predictions from altered data frames
pred_new <- pred %>%
  rename(age_cat_1_org = age_cat_1,
         age_cat_2_org = age_cat_2,
         age_num_1_org = age_num_1,
         age_num_2_org = age_num_2,
         age_cat_max_org = age_cat_max,
         age_num_max_org = age_num_max)

## create function for predicting mean
get_mean_predictions <- function(prediction_df){
  prediction_df$dyad_rank <- 1:nrow(prediction_df)
  mu <- matrix(NA, nrow = n_samples*n_chains, ncol = nrow(sim),
               dimnames = list(1:(n_samples*n_chains),
                               prediction_df$dyad_rank))
  for(dyad in 1:ncol(mu)){
    # node_effects <- mm_nodes[,paste0('mm_nodes[',prediction_df$node_1_males[dyad],']')] + mm_nodes[,paste0('mm_nodes[',prediction_df$node_2_males[dyad],']')]
    node_effects <- mm_nodes[,paste0('mm_nodes[',prediction_df$node_test_1[dyad],']')] + mm_nodes[,paste0('mm_nodes[',prediction_df$node_test_2[dyad],']')]
    for(draw in 1:nrow(mu)){
      mu[draw, dyad] <- intercept[draw] + b_max[draw]*sum(delta_j_max[draw,(1:prediction_df$age_num_max[dyad])])
    }
  }
  return(mu)
}

## create function for predicting full distribution
get_full_predictions <- function(mu_matrix, sigma){
  predictions <- mu_matrix
  for(dyad in 1:ncol(predictions)){
    for(draw in 1:nrow(predictions)){
      predictions[draw, dyad] <- rnorm(n = 1,
                                       mean = mu_matrix[draw, dyad],
                                       sd = sigma[draw])
    }
  }
  return(predictions)
}

## add progress marker
print('functions created')

## calculate predictions for each age combination
predictions_mean <- list()
predictions_full <- list()
for(i in c(1:5)){
  ## create data frame to predict from
  pred_new_age <- pred_new %>%
    mutate(age_cat_max = ifelse(i == 1, '10-15',
                                ifelse(i == 2, '15-19',
                                       ifelse(i == 3, '20-25',
                                              ifelse(i == 4, '25-40',
                                                     '40+')))),
           age_num_max = i)
  
  ## extract mean predictions
  predictions_mean[[i]] <- get_mean_predictions(prediction_df = pred_new_age)
  
  ## extract full predictions
  predictions_full[[i]] <- get_full_predictions(mu_matrix = predictions_mean[[i]],
                                                sigma = posterior_samples$sigma )
  
  ## save outputs in case R crashes
  saveRDS(object = as.data.frame(predictions_mean[[i]]),
          file = paste0('../data_processed/step5_dyadicregression/simulate_motnp_meanpredictions_maxage',i,'.RDS'))
  saveRDS(object = as.data.frame(predictions_full[[i]]),
          file = paste0('../data_processed/step5_dyadicregression/simulate_motnp_fullpredictions_maxage',i,'.RDS'))
}

# ## if image does not write properly, will need to rebuild lists from written predictions
# predictions_mean <- list()
# predictions_full <- list()
# for(i in c(1:5)){
#   ## save outputs in case R crashes
#   mean <- readRDS(file = paste0('../data_processed/step5_dyadicregression/simulate_motnp_meanpredictions_agecombo',i,'.RDS'))
#   full <- readRDS(file = paste0('../data_processed/step5_dyadicregression/simulate_motnp_fullpredictions_agecombo',i,'.RDS'))
# 
#   ## extract mean predictions
#   predictions_mean[[i]] <- mean
# 
#   ## extract full predictions
#   predictions_full[[i]] <- full
#   }

## save workspace
save.image('step5_dyadicregression/motnp_dyadic_simulation_agemax.RData')

## add progress marker
print('predictions complete')

## calculate contrasts -- age category
contrasts_1vs2 <- predictions_full[[2]] - predictions_full[[1]]
print(paste0('contrast age 1 becomes age 2: mean = ', mean(contrasts_1vs2),
             ', stdv = ', sd(contrasts_1vs2),
             ', 95% CI = [', quantile(contrasts_1vs2, prob = 0.025),':',
             quantile(contrasts_1vs2, prob = 0.975),']'))
contrasts_2vs3 <- predictions_full[[3]] - predictions_full[[2]]
print(paste0('contrast age 2 becomes age 3: mean = ', mean(contrasts_2vs3),
             ', stdv = ', sd(contrasts_2vs3),
             ', 95% CI = [', quantile(contrasts_2vs3, prob = 0.025),':',
             quantile(contrasts_2vs3, prob = 0.975),']'))
contrasts_3vs4 <- predictions_full[[4]] - predictions_full[[3]]
print(paste0('contrast age 3 becomes age 4: mean = ', mean(contrasts_3vs4),
             ', stdv = ', sd(contrasts_3vs4),
             ', 95% CI = [', quantile(contrasts_3vs4, prob = 0.025),':',
             quantile(contrasts_3vs4, prob = 0.975),']'))
contrasts_4vs5 <- predictions_full[[5]] - predictions_full[[4]]
print(paste0('contrast age 4 becomes age 5: mean = ', mean(contrasts_4vs5),
             ', stdv = ', sd(contrasts_4vs5),
             ', 95% CI = [', quantile(contrasts_4vs5, prob = 0.025),':',
             quantile(contrasts_4vs5, prob = 0.975),']'))

## add progress marker
print('contrasts calculated per age category')

## create data frame to summarise contrasts
contrasts <- expand.grid(org_max = c(1:5),
                         alt_max = c(1:5)) %>% 
  filter(org_max < alt_max) %>% 
  mutate(max_change = alt_max - org_max) %>% 
  mutate(org_max_lwr = ifelse(org_max == 1, 10,
                              ifelse(org_max == 2, 15,
                                     ifelse(org_max == 3, 20,
                                            ifelse(org_max == 4, 25,
                                                   40)))),
         alt_max_lwr = ifelse(alt_max == 1, 10,
                              ifelse(alt_max == 2, 15,
                                     ifelse(alt_max == 3, 20,
                                            ifelse(alt_max == 4, 25,
                                                   40))))) %>% 
  mutate(mean_contrast = NA,
         stdv_contrast = NA,
         lwr_contrast = NA,
         upr_contrast = NA,
         num_years_change = NA,
         mean_per_year = NA)

## add progress marker
print('contrasts data frame created')

for(i in 1:nrow(contrasts)){
  contrast_matrix <- predictions_full[[as.numeric(as.character(contrasts$alt_max[i]))]] - predictions_full[[as.numeric(as.character(contrasts$org_max[i]))]]
  contrast_matrix <- as.matrix(contrast_matrix)
  contrasts$mean_contrast[i] <- mean(contrast_matrix)
  contrasts$stdv_contrast[i] <- sd(contrast_matrix)
  contrasts$lwr_contrast[i]  <- quantile(contrast_matrix, prob = 0.025)
  contrasts$upr_contrast[i]  <- quantile(contrast_matrix, prob = 0.975)
  contrasts$num_years_change[i] <- ifelse(contrasts$max_change[i] > 0, 
                                          ifelse(contrasts$max_change[i] > 0, NA,
                                                 contrasts$alt_max_lwr[i] - contrasts$org_max_lwr[i]),
                                          contrasts$alt_max_lwr[i] - contrasts$org_max_lwr[i])
  contrasts$mean_per_year[i] <- mean(contrast_matrix) / contrasts$num_years_change[i]
}

## add progress marker
print('contrasts data frame populated')

## calculate contrats -- year-by-year
print(paste0('original (true) slope value: ', sim_max))
print(paste0('contrast per year: mean = ', mean(contrasts$mean_per_year),
             ', stdv = ', sd(contrasts$mean_per_year),
             ', 95% CI = [', quantile(contrasts$mean_per_year, prob = 0.025),':',
             quantile(contrasts$mean_per_year, prob = 0.975),']'))

## add progress marker
print('contrasts calculated per year')

## plot
contrasts <- contrasts %>% 
  mutate(org_max_cat = ifelse(org_max == 1, '10-15',
                              ifelse(org_max == 2, '15-19',
                                     ifelse(org_max == 3, '20-25',
                                            ifelse(org_max == 4, '25-40',
                                                   '40+')))),
         alt_max_cat = ifelse(alt_max == 1, '10-15',
                              ifelse(alt_max == 2, '15-19',
                                     ifelse(alt_max == 3, '20-25',
                                            ifelse(alt_max == 4, '25-40',
                                                   '40+'))))) %>% 
  mutate(org_max_cat = factor(org_max_cat,
                              levels = c('10-15','15-19','20-25',
                                         '25-40','40+')),
         alt_max_cat = factor(alt_max_cat,
                              levels = c('10-15','15-19','20-25',
                                         '25-40','40+')))
(max_plot <- ggplot(data = contrasts)+
    geom_tile(aes(x = org_max_cat, y = alt_max_cat,
                  fill = mean_contrast))+
    scale_colour_viridis_c()+
    scale_x_discrete(drop = T)+
    scale_y_discrete(drop = T)+
    labs(x = 'original maximum age category',
         y = 'altered maximum age category',
         fill = 'contrast in\nlogit(probability\nof associating)'))
(min_plot + max_plot) +
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect")

save.image('step5_dyadicregression/motnp_dyadic_simulation_agemax.RData')
dev.off()

