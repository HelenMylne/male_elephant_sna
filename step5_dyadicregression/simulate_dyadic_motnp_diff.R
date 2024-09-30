#### set up ####
# script to simulate dyadic regression and check actually working
#library(tidyverse) ; library(LaplacesDemon) ; library(car) ; library(StanHeaders) ; library(rstan)
library(StanHeaders, lib.loc = '../packages/')
library(rstan, lib.loc = '../packages/')
library(tidyverse, lib.loc = '../packages/')
library(car, lib.loc = '../packages/')
library(LaplacesDemon, lib.loc = '../packages/')
library(patchwork, lib.loc = '../packages/')

theme_set(theme_bw())

pdf('step5_dyadicregression/simulate_dyadic_sd_motnp_prep.pdf')

#### load raw data as basis for simulation ####
load('motnp_edgeweights_conditionalprior.RData')
sim <- counts_df %>%
  dplyr::select(dyad_males, node_1, node_2, node_1_males, node_2_males, id_1, id_2, count_dyad, event_count, apart, count_1, count_2, age_category_1, age_category_2, age_cat_id_1, age_cat_id_2) %>%
  rename(dyad_id = dyad_males,
         total_sightings = count_dyad,
         age_cat_1 = age_category_1, age_cat_2 = age_category_2,
         age_num_1 = age_cat_id_1, age_num_2 = age_cat_id_2) %>%
  mutate(age_num_1 = as.numeric(age_num_1) - 2,
         age_num_2 = as.numeric(age_num_2) - 2) %>%
  mutate(age_diff = abs(age_num_1 - age_num_2),
         together = ifelse(event_count == 0, 1, 2))

## simulate actual ages from categorical
nodes <- nodes %>%
  mutate(age_cat = ifelse(age < 15, 1,
                          ifelse(age < 20, 2,
                                 ifelse(age < 25, 3,
                                        ifelse(age < 40, 4, 5))))) %>%
  mutate(diff_age = ifelse(age_cat == 1, 10,
                            ifelse(age_cat == 2, 16,
                                   ifelse(age_cat == 3, 21,
                                          ifelse(age_cat == 4, 26, 40)))),
         max_age = ifelse(age_cat == 1, 15,
                            ifelse(age_cat == 2, 20,
                                   ifelse(age_cat == 3, 25,
                                          ifelse(age_cat == 4, 40, 60)))))
nodes$age_sim <- NA
for(i in 1:nrow(nodes)){
  sd <- (nodes$max_age[i] - nodes$diff_age[i])/2
  age <- round(rnorm(2500, mean = nodes$age[i], sd = sd),0)
  age <- age[which(age < nodes$max_age[i] & age > nodes$diff_age[i])[1]]
  nodes$age_sim[i] <- age
}

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

rm(list = ls()[! ls() %in% c('nodes', 'sim','counts_df')]) ; gc()

save.image('step5_dyadicregression/motnp_dyadic_simulation_dataprep.RData')
# load('motnp_dyadic_simulation_dataprep.RData') # load('step5_dyadicregression/motnp_dyadic_simulation_dataprep.RData')

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

## trim down to desired columns and rearrange
sim <- sim %>% 
  dplyr::select(dyad_id, dyad_rand,
                id_1, id_2, node_1, node_2, node_rand_1, node_rand_2,
                together,
                age_cat_1, age_cat_2, age_num_1, age_num_2, age_diff,
                age_sim_1,age_sim_2)

## print progress marker
print('data loaded')

#### simulate age effect ####
## effect of age difference on Pr(together)
sim_diff_tg <- -0.2
sim_nodes <- rnorm(n = n_nodes, mean = logit(0.01), sd = 0.4)

## simulate probability of being together
sim$pr_together <- NA ; sim$together <- NA
for(i in 1:nrow(sim)){
  sim$pr_together[i] <- sim$age_diff[i] * sim_diff_tg + 
    sim_nodes[sim$node_rand_1[i]] + sim_nodes[sim$node_rand_2[i]]
  sim$together[i] <- rbern(1, 0.1 + invlogit(sim$pr_together[i]))
}
table(sim$together)

## effect of age difference on edge weight
sim_int <- -5
sim_diff_ew <- -0.3
sim_togt_ew <- 0.6

## simulate mean centrality on normal scale
sim$mu <- NA
for(i in 1:nrow(sim)){
  sim$mu[i] <- sim_int + sim$age_diff[i] * sim_diff_ew + sim$together[i] * sim_togt_ew
}

sim$age_diff_cat <- as.factor(sim$age_diff + 1)

ggplot(sim)+
  geom_violin(aes(x = age_diff_cat, y = mu), colour = rgb(0,0,1,0.2))

# ## standardise edges
# sim$mu_std <- ( sim$mu - mean(sim$mu) ) / sd(sim$mu)
# ggplot(sim)+
#   geom_violin(aes(x = age_diff_cat, y = mu_std), colour = rgb(0,0,1,0.2))

## simulate full distribution of samples per node
sim$sd <- NA
for(i in 1:nrow(sim)){
  sim$sd[i] <- abs(sim_nodes[sim$node_rand_1[i]] + sim_nodes[sim$node_rand_2[i]])/8
}

sim_dat <- matrix(data = NA, nrow = 1000, ncol = n_dyads,
                  dimnames = list(1:1000, sim$dyad_id))    # create matrix
for(j in 1:n_dyads){
  sim_dat[,j] <- rnorm(n = nrow(sim_dat), mean = sim$mu[j], sd = sim$sd[j])  # simulate distribution
}
sim %>% 
  mutate(sim_dat1 = sim_dat[1,]) %>% 
  ggplot()+
  geom_boxplot(aes(x = age_diff_cat, y = sim_dat1, fill = as.factor(together)))

## print progress marker
print('edges simulated')

#### plot raw data ####
ggplot()+
  geom_violin(data = sim, aes(x = age_diff_cat,
                              y = mu, fill = age_diff_cat))+
  scale_x_discrete('age difference between dyad')+
  scale_y_continuous('mean estimated edge weight')+
  theme(legend.position = 'none')

sim_dat_long <- sim_dat %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = 'dyad_id', values_to = 'edge_draw') %>%
  mutate(dyad_id = as.numeric(dyad_id)) %>%
  left_join(sim, by = 'dyad_id')
ggplot()+
  geom_violin(data = sim_dat_long,
              aes(x = age_diff_cat,
                  y = edge_draw),
              fill = rgb(1,1,1,0),
              colour = rgb(0,0,0,1))+
  geom_violin(data = sim,
              aes(x = age_diff_cat, y = mu),
              fill = rgb(0.5,0,1,0.4),
              colour = rgb(1,1,1,0),
              linewidth = 1)+
  scale_x_discrete('age difference between dyad')+
  scale_y_continuous('mean edge weight')
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
pdf('step5_dyadicregression/simulate_dyadic_sd_motnp_checkoutputs_diff.pdf')

#### prior predictive check ####
n <- 100
n_age_diffs <- length(unique(sim$age_diff_cat))

## global parameters
# beta_age_diff <- rnorm(n,0,3)
beta_diff_together <- rnorm(n,0,2)
beta_together_weight <- cbind(rep(0,n),rnorm(n,1,3))
beta_diff_weight <- rnorm(n,0,3)

diff_dirichlet <- rdirichlet(n, rep(1, n_age_diffs))
diff_dirichlet <- cbind(rep(0, n), diff_dirichlet)

# global sigma
tau_sigma_raw <- rnorm(n,0,1)
tau_sigma <- exp(tau_sigma_raw)

## multimembership mean
mu_mm <- rnorm(n,logit(0.1),1)
rand_mm <- rnorm(n,0,1)
tau_mm <- rnorm(n,0,0.5)
node_mean <- mu_mm + rand_mm * tau_mm

## multimembership variance
raw_sigma <- rnorm(n,0,0.5)
node_sigma <- exp(raw_sigma)

## multimembership
mm_nodes <- matrix(NA, nrow = n, ncol = n)
for(i in 1:ncol(mm_nodes)){
  mm_nodes[,i] <- rnorm(n, mean = node_mean[i], sd = node_sigma[i])
}

## plot
par(mfrow = c(2,2))
x_age <- rep(1:n_age_diffs, 2)
x_tgt <- rep(1:2, each = n_age_diffs)
for(bounds in c(100, 50, 20, 5)){
  
  plot(NULL, las = 1, xlab = 'difference age categeory', 
       ylab = 'logit edge weight (standardised)',
       ylim = c(min(logit_edge_draws_mu)-bounds,
                max(logit_edge_draws_mu)+bounds),#c(-8,10),
       xlim = c(1,n_age_diffs))
  for(i in 1:n){
    y <- rep(NA, length(x_age))#length(x))
    nodes_sample <- sample(1:n, size = 2, replace = F)
    for(j in 1:length(x_age)){#length(x)){
      # y[j] <- beta_age_diff[i]*sum(diff_dirichlet[i,][1:x[j]]) + mm_nodes[i,nodes_sample[1]] + mm_nodes[i,nodes_sample[2]]
      y[j] <- beta_diff_weight[i] * sum(diff_dirichlet[i,][1:x_age[j]]) +
        beta_together_weight[i,x_tgt[j]] + 
        mm_nodes[i,nodes_sample[1]] + mm_nodes[i,nodes_sample[2]]
    }
    lines(x = x_age[1:n_age_diffs],
          y = y[1:n_age_diffs],
          col = rgb(0,0,1,0.4))
    lines(x = x_age[(n_age_diffs+1):length(x_age)],
          y = y[1:n_age_diffs],
          col = rgb(0,0,1,0.4))
  }
    abline(h = min(logit_edge_draws_mu),
           lty = 2) ;
    abline(h = max(logit_edge_draws_mu),
           lty = 2)
}
par(mfrow = c(1,1))

# # regression equation for Pr(ever together)
# vector[num_dyads] predictor_together;
# for (i in 1:num_dyads) {
#   predictor_together[i] = beta_diff_together * sum(dj_diff_together[1:age_diff[i]]) + mm_nodes[node_1[i]] + mm_nodes[node_2[i]];
# }
# 
# # regression equation for edge weight
# vector[num_dyads] predictor_edge;
# for (i in 1:num_dyads) {
#   predictor_edge[i] = beta_diff_weight * sum(dj_diff_weight[1:age_diff[i]]) + together * beta_together_weight + mm_nodes[node_1[i]] + mm_nodes[node_2[i]];
# }
# 
# together ~ bernoulli(predictor_together)
# logit_edge_mu ~ normal(predictor_edge, logit_edge_sd + rep_vector(tau_sigma, num_dyads));

# rm(n, beta_age_diff, intercept, diff_dirichlet, i, y, j, x, mm_nodes) ; gc()
rm(n, beta_diff_together, beta_together_weight, beta_diff_weight, diff_dirichlet, tau_sigma_raw, tau_sigma, mu_mm, rand_mm, tau_mm, node_mean, raw_sigma, node_sigma, mm_nodes, x_age, x_tgt, y, nodes_sample) ; gc()

## print progress marker
print('prior predictive complete')

#### fit dyadic regression ####
## create data list
dyad_data <- list(
  # global data size
  num_dyads = n_dyads,                      # number of dyads
  num_nodes = n_nodes,                      # number of nodes
  num_age_diffs = n_age_diffs,              # number of unique age categories
  length_dirichlet = n_age_diffs + 1,       # number of unique age categories + 1
  
  # Gaussian approximation of edge weights
  logit_edge_mu = logit_edge_draws_mu,      # means of the logit edge weights
  logit_edge_sd = logit_edge_draws_sd,      # standard deviation of logit edge weights
  
  # explanatory variables
  age_diff = as.integer(sim$age_diff_cat),  # age of younger dyad member
  together = as.integer(sim$together),
  
  # node IDs
  node_1 = sim$node_rand_1,                 # node IDs for multimembership effects
  node_2 = sim$node_rand_2,                 # node IDs for multimembership effects
  
  # simplex prior
  prior_diff = rep(1, n_age_diffs)          # prior for difference age slope
)

## print progress marker
print('data list created')

## save output
save.image('step5_dyadicregression/motnp_dyadic_simulation_agediff_together.RData')

## load dyadic regression model
dyadic_regression <- stan_model('models/dyadic_regression_motnp_agediff_together.stan')
dyadic_regression

## print progress marker
print('model loaded in')

## define parameters
n_chains <- 4
n_samples <- 2500

## fit dyadic regression
print(paste0('start model at ', Sys.time()))
fit_dyadreg_sim <- sampling(dyadic_regression,
                            data = dyad_data,
                            iter = n_samples*3, warmup = n_samples*2,
                            chains = n_chains, cores = n_chains,
                            control = list(adapt_delta = 0.9,
                                           max_treedepth = 15))

## save output
save.image('step5_dyadicregression/motnp_dyadic_simulation_agediff_together.RData')

## print progress marker
print('model run')

#### check outputs ####
# load('step5_dyadicregression/motnp_dyadic_simulation_agediff_together.RData')

## extract draws
posterior_samples <- fit_dyadreg_sim %>%
  as.data.frame()

## extract model fit
s <- summary(fit_dyadreg_sim)
summary <- s$summary %>%
  as.data.frame() %>%
  filter(is.nan(se_mean) == FALSE)
summary
par(mfrow = c(2,1))
hist(summary$Rhat,  breaks = 50)
hist(summary$n_eff, breaks = 50)
par(mfrow = c(1,1))

rownames(summary)[summary$n_eff < 600]
summary[summary$n_eff < 600,c(1,3,9,10)]

save.image('step5_dyadicregression/motnp_dyadic_simulation_agediff_together.RData')

## extract dyadic regression slopes
# b_diff <- posterior_samples$beta_age_diff
# sigma <- posterior_samples$tau_sigma
# delta_j_diff <- posterior_samples[,c('delta_j_diff[1]','delta_j_diff[2]','delta_j_diff[3]','delta_j_diff[4]','delta_j_diff[5]','delta_j_diff[6]')] ; colnames(delta_j_diff) <- 1:(n_age_diffs+1)
# mm_nodes <- posterior_samples[,grep(pattern = 'mm_nodes', x = colnames(posterior_samples))]
# parameters <- data.frame(beta_age_diff = b_diff,
#                          # intercept = intercept,
#                          sigma = sigma) %>%
#   mutate(chain = rep(1:4, each = n_samples),
#          position = rep(1:n_samples, 4)) %>%
#   pivot_longer(cols = c('beta_age_diff','sigma'),#,'intercept'),
#                names_to = 'parameter', values_to = 'slope_draw')
b_diff_tg <- posterior_samples$beta_diff_together
b_togt_ew <- posterior_samples$beta_together_weight
b_diff_ew <- posterior_samples$beta_diff_weight
tau_sigma <- posterior_samples$tau_sigma
dj_diff_tg <- posterior_samples[,c('dj_diff_together[1]','dj_diff_together[2]','dj_diff_together[3]','dj_diff_together[4]','dj_diff_together[5]','dj_diff_together[6]')] ; colnames(dj_diff_tg) <- 1:(n_age_diffs+1)
dj_diff_ew <- posterior_samples[,c('dj_diff_weight[1]','dj_diff_weight[2]','dj_diff_weight[3]','dj_diff_weight[4]','dj_diff_weight[5]','dj_diff_weight[6]')] ; colnames(dj_diff_ew) <- 1:(n_age_diffs+1)
mm_nodes <- posterior_samples[,grep(pattern = 'mm_nodes', x = colnames(posterior_samples))]
parameters <- data.frame(b_diff_tg = b_diff_tg,
                         b_togt_ew = b_togt_ew,
                         b_diff_ew = b_diff_ew,
                         tau_sigma = tau_sigma) %>%
  mutate(chain = rep(1:4, each = n_samples),
         position = rep(1:n_samples, 4)) %>%
  pivot_longer(cols = c('b_diff_tg','b_togt_ew','b_diff_ew','tau_sigma'),
               names_to = 'parameter', values_to = 'slope_draw')

## traceplots
ggplot(data = parameters)+
  geom_line(aes(x = position,
                y = slope_draw,
                colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')
dj_diff_tg %>% as.data.frame() %>%
  mutate(chain = rep(1:4, each = n_samples),
         position = rep(1:n_samples, 4)) %>%
  pivot_longer(cols = all_of(2:(n_age_diffs+1)),
               names_to = 'parameter',
               values_to = 'slope_draw') %>%
  ggplot()+
  geom_line(aes(x = position,
                y = slope_draw,
                colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')
dj_diff_ew %>% as.data.frame() %>%
  mutate(chain = rep(1:4, each = n_samples),
         position = rep(1:n_samples, 4)) %>%
  pivot_longer(cols = all_of(2:(n_age_diffs+1)),
               names_to = 'parameter',
               values_to = 'slope_draw') %>%
  ggplot()+
  geom_line(aes(x = position,
                y = slope_draw,
                colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')

sample_nodes <- sample(colnames(mm_nodes), 25, replace = F)
mm_nodes %>%
  as.data.frame() %>%
  select(all_of(sample_nodes)) %>%
  mutate(chain = rep(1:4, each = n_samples),
         position = rep(1:n_samples, 4)) %>%
  pivot_longer(cols = all_of(1:25),
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
     main = "Posterior predictive density of edge weights:\nblack = measured, colours = predicted (col by chain)",
     ylim = c(0, 1), col = rgb(0, 0, 0, 0.25), las = 1)
for (i in 1:100) {
  j <- sample(1:nrow(sim_dat), 1)

  mu_togt <- matrix(NA, ncol = n_dyads, nrow = n_chains)
  mu_edge <- matrix(NA, ncol = n_dyads, nrow = n_chains)

  for(l in c(1:n_chains)){
    n <- (l*n_samples) - n_samples
    for(k in 1:n_dyads){
      # mu_edge[l,k] <- b_diff[j+n]*sum(delta_j_diff[j,(1:dyad_data$age_diff[k])]) +
      #   mm_nodes[j+n,dyad_data$node_1[k]] + mm_nodes[j+n,dyad_data$node_2[k]]

      mu_togt[l,k] <- b_diff_tg[j+n] * sum(dj_diff_tg[j+n,(1:dyad_data$age_diff[k])]) +
        mm_nodes[j+n,dyad_data$node_1[k]] + mm_nodes[j+n,dyad_data$node_2[k]]

      mu_edge[l,k] <- b_diff_ew[j+n] * sum(dj_diff_ew[j+n,(1:dyad_data$age_diff[k])]) +
        b_togt_ew[j+n] * mu_togt[l,k]
    }
    sigma_edge <- dyad_data$logit_edge_sd + diag(rep(tau_sigma[j+n], n_dyads))
    norm <- rnorm(n_samples, mu_edge[l,], sigma_edge)
    lines(density(as.numeric(sim_dat[j, ])), col = rgb(0, 0, 0, 0.25))
    lines(density(norm),
          col = ifelse(l == 1, rgb(1, 0, 0, 0.25),
                       ifelse(l == 2, rgb(0, 1, 0, 0.25),
                              ifelse(l == 3, rgb(0, 0, 1, 0.25),
                                     rgb(1, 1, 0, 0.25)))) )
  }
}
rm(mu_edge, mu_togt, sigma_edge, i,j,k,l,n) ; gc()

## save output
save.image('step5_dyadicregression/motnp_dyadic_simulation_agediff_together.RData')
dev.off()

## print progress marker
print('outputs checked')

#### predict from raw data ####
pdf('step5_dyadicregression/simulate_dyadic_sd_motnp_predictions_agediff_together.pdf')
#load('step5_dyadicregression/motnp_dyadic_simulation_agediff_together.RData')

## create predictive data frame
pred <- sim %>% 
  mutate(dyad_rank = 1:n_dyads) %>% 
  relocate(dyad_rank, .after = 'dyad_rand') %>% 
  relocate(age_diff_cat, .after = 'age_diff')
# nodes$node_test_1 <- 1:nrow(nodes) ; nodes$node_test_2 <- 1:nrow(nodes) ; nodes$node_1 <- nodes$node ; nodes$node_2 <- nodes$node
# pred <- pred %>%
#   left_join(nodes[,c('node_1','node_test_1')], by = 'node_1') %>%
#   left_join(nodes[,c('node_2','node_test_2')], by = 'node_2')

str(mm_nodes)
dim(mm_nodes)

## calculate mean predictions
pred_mu <- matrix(NA, nrow = 100,#n_samples*n_chains,
                  ncol = nrow(sim),
                  dimnames = list(1:100,#(n_samples*n_chains),
                                  pred$dyad_rank))
for(dyad in 1:ncol(pred_mu)){
  node_effects <- mm_nodes[,paste0('mm_nodes[',pred$node_rand_1[dyad],']')] + mm_nodes[,paste0('mm_nodes[',pred$node_rand_2[dyad],']')]
  for(draw in 1:nrow(pred_mu)){
    # pred_mu[draw, dyad] <- b_diff[draw]*sum(delta_j_diff[draw,(1:pred$age_diff_cat[dyad])]) + node_effects[draw]
    
    pred_mu[draw,dyad] <- ( b_diff_tg[draw] * sum(dj_diff_tg[draw,(1:dyad_data$age_diff[dyad])]) + 
                            mm_nodes[draw,dyad_data$node_1[dyad]] + 
                              mm_nodes[draw,dyad_data$node_2[dyad]] ) * b_togt_ew[draw] +
      b_diff_ew[draw] * sum(dj_diff_ew[draw,(1:dyad_data$age_diff[dyad])])
    
    }
  if(dyad %% 100 == 0){ print(dyad) }
}

pred$pred_mu <- apply(pred_mu, 2, mean)
pred$pred_mu_lwr <- apply(pred_mu, 2, quantile, prob = 0.025)
pred$pred_mu_upr <- apply(pred_mu, 2, quantile, prob = 0.975)

save.image('step5_dyadicregression/motnp_dyadic_simulation_agediff_together.RData')

## plot predicted vs original
ggplot(pred)+
  geom_point(aes(x = mu, y = pred_mu,
                 colour = as.factor(age_diff_cat)))+
  geom_abline(slope = 1, intercept = 0)+
  scale_colour_viridis_d()+
  labs(colour = 'age of younger\ndyad member')
ggplot(pred)+
  geom_point(aes(x = mu, y = pred_mu,
                 colour = as.factor(together)),
             alpha = 0.1)+
  geom_abline(slope = 1, intercept = 0)+
  scale_colour_viridis_d(begin = 0.5, end = 0)+
  labs(colour = 'ever together\n1=no,2=yes')

## calculate full prediction distribution
pred_full <- pred_mu
for(dyad in 1:ncol(pred_mu)){
  for(draw in 1:nrow(pred_mu)){
    pred_full[draw, dyad] <- rnorm(n = 1, mean = pred_mu[draw, dyad], sd = logit_edge_draws_sd[dyad] + posterior_samples$tau_sigma[draw])
  }
  if(dyad %% 100 == 0){ print(dyad) }
}
save.image('step5_dyadicregression/motnp_dyadic_simulation_agediff_together.RData')
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

## plot
ggplot()+
  geom_violin(data = pred_full_df,
              mapping = aes(x = age_diff_cat,
                            y = prediction),
              alpha = 0.4)+
  geom_boxplot(data = pred,
               mapping = aes(x = age_diff_cat,
                             y = pred_mu),
               width = 0.2,
               position = position_dodge(0.9))+
  geom_point(data = sim,
             mapping = aes(x = age_diff_cat,
                           # size = total_sightings,
                           y = mu),
             pch = 21,
             position = position_dodge(0.9))+
  scale_fill_viridis_d()+
  labs(x = 'difference in age categories',
       y = 'logit edge weight')

ggplot()+ # purple = predicted mean, clear = true means
  geom_violin(data = pred,
               mapping = aes(x = age_diff_cat,
                             y = pred_mu),
              fill = rgb(0.5,0,1,0.4),
              colour = rgb(1,1,1,0))+
  geom_violin(data = sim,
             mapping = aes(x = age_diff_cat,
                           y = mu),
             fill = rgb(1,1,1,0),
             colour = rgb(0,0,0,1))+
  scale_fill_viridis_d()+
  labs(x = 'difference in age categories',
       y = 'logit edge weight')

ggplot()+
  geom_violin(data = pred_full_df,
              mapping = aes(x = age_diff_cat,
                            y = invlogit(prediction)),
              fill = rgb(0.5,0,1,0.4), colour = 'transparent')+
  geom_boxplot(data = pred,
               mapping = aes(x = age_diff_cat,
                             y = invlogit(pred_mu)),
               width = 0.2,
               position = position_dodge(0.9))+
  # geom_point(data = sim,
  #            mapping = aes(x = age_diff_cat,
  #                          # size = total_sightings,
  #                          y = invlogit(mu_std)),
  #            pch = 21,
  #            position = position_dodge(0.9))+
  geom_violin(data = sim,
              mapping = aes(x = age_diff_cat,
                            y = invlogit(mu)),
              fill = rgb(1,1,1,0))+
  scale_fill_viridis_d()+
  labs(x = 'difference in age categories',
       y = 'edge weight')

dev.off()
save.image('step5_dyadicregression/motnp_dyadic_simulation_agediff_together.RData')

## print progress marker
print('predictions complete')

# #### extract original slope values ####
# #load('step5_dyadicregression/motnp_dyadic_simulation_agediff_together.RData')
# pdf('step5_dyadicregression/simulate_dyadic_sd_motnp_contrasts_agediff_together.pdf')
# 
# ## calculate predictions from altered data frames
# pred_new <- pred %>%
#   mutate(age_diff_org = as.integer(age_diff_cat),
#          age_diff_new = ifelse(age_diff_cat == 5, 1, as.integer(age_diff_cat) + 1))
# 
# ## create function for predicting mean
# get_mean_predictions <- function(prediction_df, new){
#   prediction_df$dyad_rank <- 1:nrow(prediction_df)
#   mu <- matrix(NA,
#                nrow = 100,#n_samples*n_chains,
#                ncol = nrow(sim),
#                     dimnames = list(1:100,#(n_samples*n_chains),
#                                     prediction_df$dyad_rank))
#   if(new == TRUE){
#     prediction_df <- prediction_df %>% 
#       rename(age_num_diff = age_diff_new)
#   }
#   if(new == FALSE){
#     prediction_df <- prediction_df %>%
#       rename(age_num_diff = age_diff_org)
#   }
#   
#   for(dyad in 1:ncol(mu)){
#     node_effects <- mm_nodes[,paste0('mm_nodes[',prediction_df$node_rand_1[dyad],']')] + 
#       mm_nodes[,paste0('mm_nodes[',prediction_df$node_rand_2[dyad],']')]
#     for(draw in 1:nrow(mu)){
#       mu[draw, dyad] <- #intercept[draw] + 
#         b_diff[draw]*sum(delta_j_diff[draw,(1:prediction_df$age_num_diff[dyad])]) + node_effects[draw]
#     }
#   }
#   return(mu)
#   }
# 
# ## create function for predicting full distribution
# get_full_predictions <- function(mu_matrix, sigma){
#   predictions <- mu_matrix
#   for(dyad in 1:ncol(predictions)){
#     for(draw in 1:nrow(predictions)){
#       predictions[draw, dyad] <- rnorm(n = 1,
#                                        mean = mu_matrix[draw, dyad],
#                                        sd = sigma[draw])
#     }
#   }
#   return(predictions)
# }
# 
# ## add progress marker
# print('functions created')
# 
# # ## calculate predictions for each age combination
# # predictions_mean <- list()
# # predictions_full <- list()
# # for(i in c(1:5)){
# #   ## create data frame to predict from
# #   pred_new_age <- pred_new %>%
# #     mutate(age_diff_cat = ifelse(i == 1, '10-15',
# #                                 ifelse(i == 2, '16-20',
# #                                        ifelse(i == 3, '20-25',
# #                                               ifelse(i == 4, '25-40',
# #                                                      '40+')))),
# #            age_num_diff = i)
# # 
# #   ## extract mean predictions
# #   predictions_mean[[i]] <- get_mean_predictions(prediction_df = pred_new_age)
# # 
# #   ## extract full predictions
# #   predictions_full[[i]] <- get_full_predictions(mu_matrix = predictions_mean[[i]],
# #                                                 sigma = posterior_samples$sigma )
# # 
# #   ## save outputs in case R crashes
# #   saveRDS(object = as.data.frame(predictions_mean[[i]]),
# #           file = paste0('../data_processed/step5_dyadicregression/simulate_motnp_meanpredictions_diffage',i,'.RDS'))
# #   saveRDS(object = as.data.frame(predictions_full[[i]]),
# #           file = paste0('../data_processed/step5_dyadicregression/simulate_motnp_fullpredictions_diffage',i,'.RDS'))
# # }
# 
# # ## if image does not write properly, will need to rebuild lists from written predictions
# # predictions_mean <- list()
# # predictions_full <- list()
# # for(i in c(1:5)){
# #   ## save outputs in case R crashes
# #   mean <- readRDS(file = paste0('../data_processed/step5_dyadicregression/simulate_motnp_meanpredictions_agecombo',i,'.RDS'))
# #   full <- readRDS(file = paste0('../data_processed/step5_dyadicregression/simulate_motnp_fullpredictions_agecombo',i,'.RDS'))
# #
# #   ## extract mean predictions
# #   predictions_mean[[i]] <- mean
# #
# #   ## extract full predictions
# #   predictions_full[[i]] <- full
# #   }
# 
# ## calculate mean predictions
# pred_mu_org <- get_mean_predictions(prediction_df = pred_new, new = FALSE)
# pred_mu_new <- get_mean_predictions(prediction_df = pred_new, new = TRUE)
# 
# ## calculate full predictions
# pred_full_org <- get_full_predictions(mu_matrix = pred_mu_org, sigma = posterior_samples$tau_sigma)
# pred_full_new <- get_full_predictions(mu_matrix = pred_mu_new, sigma = posterior_samples$tau_sigma)
# 
# ## save workspace
# save.image('step5_dyadicregression/motnp_dyadic_simulation_agediff_together.RData')
# # load('step5_dyadicregression/motnp_dyadic_simulation_agediff_together.RData')
# # rm(list = ls()[! ls() %in% c('predictions_full', 'predictions_mean','sim_diff')]) ; gc()
# 
# ## add progress marker
# print('predictions complete')
# 
# # ## calculate contrasts -- age category
# # contrasts_1vs2 <- predictions_full[[2]] - predictions_full[[1]]
# # print(paste0('contrast age 1 becomes age 2: mean = ', mean(contrasts_1vs2),
# #              ', stdv = ', sd(contrasts_1vs2),
# #              ', 95% CI = [', quantile(contrasts_1vs2, prob = 0.025),':',
# #              quantile(contrasts_1vs2, prob = 0.975),']'))
# # contrasts_2vs3 <- predictions_full[[3]] - predictions_full[[2]]
# # print(paste0('contrast age 2 becomes age 3: mean = ', mean(contrasts_2vs3),
# #              ', stdv = ', sd(contrasts_2vs3),
# #              ', 95% CI = [', quantile(contrasts_2vs3, prob = 0.025),':',
# #              quantile(contrasts_2vs3, prob = 0.975),']'))
# # contrasts_3vs4 <- predictions_full[[4]] - predictions_full[[3]]
# # print(paste0('contrast age 3 becomes age 4: mean = ', mean(contrasts_3vs4),
# #              ', stdv = ', sd(contrasts_3vs4),
# #              ', 95% CI = [', quantile(contrasts_3vs4, prob = 0.025),':',
# #              quantile(contrasts_3vs4, prob = 0.975),']'))
# # contrasts_4vs5 <- predictions_full[[5]] - predictions_full[[4]]
# # print(paste0('contrast age 4 becomes age 5: mean = ', mean(contrasts_4vs5),
# #              ', stdv = ', sd(contrasts_4vs5),
# #              ', 95% CI = [', quantile(contrasts_4vs5, prob = 0.025),':',
# #              quantile(contrasts_4vs5, prob = 0.975),']'))
# # 
# # ## add progress marker
# # print('contrasts calculated per age category')
# # 
# # ## create data frame to summarise contrasts
# # contrasts <- expand.grid(org_diff = c(1:5),
# #                          alt_diff = c(1:5)) %>%
# #   filter(org_diff < alt_diff) %>%
# #   mutate(diff_change = alt_diff - org_diff) %>%
# #   mutate(org_diff_lwr = ifelse(org_diff == 1, 10,
# #                               ifelse(org_diff == 2, 15,
# #                                      ifelse(org_diff == 3, 20,
# #                                             ifelse(org_diff == 4, 25,
# #                                                    40)))),
# #          alt_diff_lwr = ifelse(alt_diff == 1, 10,
# #                               ifelse(alt_diff == 2, 15,
# #                                      ifelse(alt_diff == 3, 20,
# #                                             ifelse(alt_diff == 4, 25,
# #                                                    40))))) %>%
# #   mutate(mean_contrast = NA,
# #          stdv_contrast = NA,
# #          lwr_contrast = NA,
# #          upr_contrast = NA,
# #          num_years_change = NA,
# #          mean_per_year = NA)
# # 
# # ## add progress marker
# # print('contrasts data frame created')
# # 
# # for(i in 1:nrow(contrasts)){
# #   contrast_matrix <- predictions_full[[as.numeric(as.character(contrasts$alt_diff[i]))]] - predictions_full[[as.numeric(as.character(contrasts$org_diff[i]))]]
# #   contrast_matrix <- as.matrix(contrast_matrix)
# #   contrasts$mean_contrast[i] <- mean(contrast_matrix)
# #   contrasts$stdv_contrast[i] <- sd(contrast_matrix)
# #   contrasts$lwr_contrast[i]  <- quantile(contrast_matrix, prob = 0.025)
# #   contrasts$upr_contrast[i]  <- quantile(contrast_matrix, prob = 0.975)
# #   contrasts$num_years_change[i] <- contrasts$alt_diff_lwr[i] - contrasts$org_diff_lwr[i]
# #   contrasts$mean_per_year[i] <- mean(contrast_matrix) / contrasts$num_years_change[i]
# # }
# 
# small_age_gap <- which(pred_new$age_diff_org != 5)
# big_age_gap <- which(pred_new$age_diff_org == 5)
# contrasts <- pred_full_new[,small_age_gap] - pred_full_org[,small_age_gap]
# contrasts5 <- pred_full_new[,big_age_gap] - pred_full_org[,big_age_gap] * (-1)
# 
# ## add progress marker
# print('contrasts data frame populated')
# 
# ## calculate contrasts
# print(paste0('contrast per category: mean = ', round(mean(contrasts),4),
#              ', stdv = ', round(sd(contrasts),4),
#              ', 95% CI = [', round(quantile(contrasts, prob = 0.025, na.rm = T),4),':',
#              round(quantile(contrasts, prob = 0.975, na.rm = T),4),']'))
# print(paste0('contrast big vs small: mean = ', round(mean(contrasts5),4),
#              ', stdv = ', round(sd(contrasts5),4),
#              ', 95% CI = [', round(quantile(contrasts5, prob = 0.025, na.rm = T),4),':',
#              round(quantile(contrasts5, prob = 0.975, na.rm = T),4),']'))
# print(paste0('original (true) slope value: ', sim_diff))
# # can't do contrast per year because contrast in age difference is meaningless -- shifting from a difference of 1 category to 2 categories could mean anything in terms of actual age difference between them
# 
# ## add progress marker
# print('contrasts calculated')
# 
# ## plot
# # contrasts5 <- contrasts5 %>%
# #   as.data.frame() %>%
# #   pivot_longer(cols = everything(), names_to = 'dyad_rank', values_to = 'contrast') %>%
# #   mutate(dyad_rank = as.numeric(dyad_rank)) %>%
# #   left_join(pred_new, by = 'dyad_rank')
# contrasts_long <- contrasts %>% 
#   as.data.frame() %>% 
#   pivot_longer(cols = everything(), names_to = 'dyad_rank', values_to = 'contrast') %>% 
#   mutate(dyad_rank = as.numeric(dyad_rank)) %>% 
#   left_join(pred_new, by = 'dyad_rank')
# # contrasts_long <- rbind(contrasts_long, contrasts5)
# 
# contrast_summary <- contrasts_long %>% 
#   group_by(dyad_rand) %>% 
#   summarise(mean_contrast = mean(contrast)) %>% 
#   left_join(pred_new[,c('dyad_rand','age_diff_org','age_diff_new')]) %>% 
#   mutate(change = paste0(age_diff_org-1, ' -> ', age_diff_new-1))
# # (diff_plot <- ggplot(data = contrast_summary)+
# #     geom_density(aes(x = mean_contrast, fill = change))+
# #     scale_fill_viridis_d()+
# #     facet_wrap(. ~ change, scales = 'free_x')+
# #     labs(x = 'contrast in logit(probability of associating)',
# #          fill = 'change in age\ndifference'))
# ggplot(data = contrast_summary)+
#     geom_density_ridges(aes(x = mean_contrast, y = change, fill = change))+
#     scale_fill_viridis_d()+
#     labs(x = 'contrast in logit(probability of associating)',
#          y = 'change in age difference',
#          fill = 'change in age\ndifference')
# 
# save.image('step5_dyadicregression/motnp_dyadic_simulation_agediff_together.RData')
# dev.off()
# 
