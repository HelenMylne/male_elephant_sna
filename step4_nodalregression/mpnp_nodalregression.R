#### information ####
# Makgadikgadi Pans National Park -- regression to test effect of individual age on eigenvector centrality

#### set up ####
#options(future.globals.maxSize = 10000*(1024^2))   # running model with full age distributions = 5.01GB of globals, which exceeds default maximum allowed size. Set to allow up to 10 GB to allow model to run

# library(rstan) ; library(igraph) ; library(tidyverse) ; library(LaplacesDemon) ; library(MASS)
library(StanHeaders, lib.loc = '../packages/')
library(rstan, lib.loc = '../packages/')
library(igraph, lib.loc = '../packages/')          # library(igraph)
library(tidyverse, lib.loc = '../packages/')       # library(tidyverse)
library(LaplacesDemon, lib.loc = '../packages/')   # library(LaplacesDemon)
library(MASS, lib.loc = '../packages/')            # library(MASS)
#library(cmdstanr, lib.loc = '../packages/')        # library(cmdstanr)
#library(brms, lib.loc = '../packages/')            # library(brms)
#library(Rcpp, lib.loc = '../packages/')            # library(Rcpp)
#library(ggdist, lib.loc = '../packages/')          # library(ggdist)
#library(posterior, lib.loc = '../packages/')       # library(posterior)
#library(bayesplot, lib.loc = '../packages/')       # library(bayesplot)
#library(bisonR, lib.loc = '../packages/')          # library(bisonR)
#library(janitor, lib.loc = '../packages/')         # library(janitor)

## set cmdstan path
#set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0/')

## set seed for reproducibility
set.seed(12345)

##### run model with rstan -- long window ####
## load data and remove additional data
load('mpnp_edgecalculations/mpnplong_edgeweights_conditionalprior.RData')
rm(eles, edgelist, summary, edge_binary, fit_edges_mpnp1, mean_ages, missing_age, make_edgelist, plot_network_threshold_mpnp) ; gc()

## set up pdf
pdf('../outputs/mpnpshort1_nodalregression.pdf')

### prior predictive check ####
age <- 10:60
mu_mean <- 0
mu_stdv <- 0.1
mu <- rnorm(100, mu_mean, mu_stdv)
mu <- sort(mu)
mu_mtrx <- matrix(NA, nrow = length(mu), ncol = length(age))
for(i in 1:nrow(mu_mtrx)){
  mu_mtrx[i,] <- mu[i]*age
}
sigma_range <- rexp(25, 2)
sigma_range <- sort(sigma_range)
par(mfrow = c(5,5), mai = c(0.2,0.2,0.2,0.2))
for(j in 1:25){
  plot(NULL, xlim = c(10,60), ylim = c(0,1), las = 1,
       xlab = '', ylab = '')
  sigma <- diag(rep(sigma_range[j], length(age)))
  for(i in 1:length(mu)){
    predictor <- mu_mtrx[i,]
    y <- MASS::mvrnorm(1, predictor, sigma)
    lines(x = age, y = invlogit(y), col = rgb(0,0,1,1))
  }
}
par(mfrow = c(1,1))

### filter down to only elephants with known age categories ####
## ages
hist(nodes$age, breaks = 50)
nodes$age_cat <- ifelse(nodes$age < 10, 1,
                        ifelse(nodes$age < 16, 2,
                               ifelse(nodes$age < 21, 3,
                                      ifelse(nodes$age < 26, 4,
                                             ifelse(nodes$age < 36, 5, 6)))))
## nodes data frame
nodes_to_include <- nodes %>% 
  filter(! is.na(age_cat))
n_age_cat <- length(unique(nodes_to_include$age_cat))
ele_ids <- nodes_to_include$id
n_eles <- length(ele_ids)

## dyads data frame
dyads_to_include <- counts_df %>% 
  filter(id_1 %in% ele_ids) %>% 
  filter(id_2 %in% ele_ids)
n_dyads <- nrow(dyads_to_include)

### extract centralities ####
## build adjacency tensor
ncol(edge_samples) ; nrow(dyads_to_include)
node_join <- nodes_to_include %>% 
  mutate(node_rank = as.integer(as.factor(node))) %>% 
  dplyr::select(node_rank, node, id) %>% 
  rename(node_1 = node, id_1 = id, node_rank_1 = node_rank) %>% 
  mutate(node_2 = node_1, id_2 = id_1, node_rank_2 = node_rank_1)
dyads_to_include <- dyads_to_include %>% 
  left_join(node_join[,c('node_rank_1','node_1','id_1')], by = c('node_1', 'id_1')) %>% 
  left_join(node_join[,c('node_rank_2','node_2','id_2')], by = c('node_2', 'id_2'))
adj_tensor <- array(0, c(n_chains*n_samples, n_eles, n_eles))
for (dyad_id in 1:n_dyads) {
  dyad_row <- dyads_to_include[dyad_id, ]
  adj_tensor[, dyad_row$node_rank_1, dyad_row$node_rank_2] <- edge_samples[, dyad_id]
}

## calculate centrality and store posterior samples in a matrix
centrality_samples <- matrix(0, n_chains*n_samples, n_eles)
centrality_samples_std <- matrix(0, n_chains*n_samples, n_eles)
for (i in 1:(n_chains*n_samples)) {
  g <- graph_from_adjacency_matrix(adj_tensor[i, , ], mode="undirected", weighted=TRUE)
  centrality_samples[i, ] <- eigen_centrality(g)$vector
  centrality_samples_std[i, ] <- (centrality_samples[i, ] - mean(centrality_samples[i, ]))/sd(centrality_samples[i, ])
}
head(centrality_samples) # Unstandardised eigenvector centrality
head(centrality_samples_std)

## visualise centralities
nodes_to_include$node_rank <- as.integer(as.factor(nodes_to_include$node))
df_wide <- data.frame(centrality_samples_std)
colnames(df_wide) <- 1:n_eles
df_long <- pivot_longer(df_wide, cols = everything(),
                        names_to = "node_rank", values_to = "centrality") %>% 
  mutate(node_rank = as.integer(node_rank)) %>% 
  left_join(nodes_to_include[,c('node_rank','age')], by = 'node_rank')
df_long %>% 
  filter(node_rank <= 50) %>% 
  mutate(nodes_reordered = fct_reorder(.f = as.factor(node_rank), .x = age, .desc = T)) %>% 
  ggplot(aes(x = centrality, fill = age)) +
  geom_density(linewidth = 0.4) +
  facet_grid(rows = vars(as.factor(nodes_reordered)), scales = "free") +
  labs(x = "Eigenvector centrality (standardised)") + 
  theme_void() + 
  theme(strip.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

### compute normal approximation ####
## check covariance
plot(centrality_samples_std[, 1], centrality_samples_std[, 2],
     xlab = 'standardised eigenvectors ID1',
     ylab = 'standardised eigenvectors ID2',
     las = 1, pch = 19, col = rgb(0,0,1,0.2))

## check variance
par(mfrow = c(5,5), mai = c(0.2,0.2,0.2,0.2))
to_plot <- sample(1:ncol(centrality_samples), 25, replace = F)
for(i in 1:length(to_plot)){
  hist(centrality_samples[,to_plot[i]], main = '')
}

## compute normal approximation
centrality_mu <- apply(centrality_samples_std, 2, mean)
centrality_cov <- cov(centrality_samples_std)

centrality_samples_sim <- MASS::mvrnorm(1e5, centrality_mu, centrality_cov)

for(i in 1:length(to_plot)){
  plot(density(centrality_samples_std[, i]), lwd = 2, main = "", xlab = "")
  lines(density(centrality_samples_sim[, i]), col = rgb(0,0,1,0.5), lwd = 2)
}
par(mfrow = c(1,1), mai = c(1,1,1,1))
rm(g, edge_samples, adj_tensor, i) ; gc()

save.image('mpnp_nodalregression/mpnp_long_nodalregression.RData')

### run model ####
# ## extract age distributions (using a normal for now)
# mpnp_ages <- readRDS('../data_processed/step2_ageestimation/mpnp_longwindow_ageestimates_mcmcoutput.rds')
# mpnp_ages <- motnp_ages[,which(colnames(mpnp_ages) %in% ele_ids)]
# age_mu <- apply(motnp_ages, 2, mean)
# age_sd <- apply(motnp_ages, 2, sd)

## create data list
eigen_list <- list(num_nodes = n_eles,
                   num_age_cat = n_age_cat,
                   length_dirichlet = n_age_cat + 1,
                   #nodes = nodes_to_include$node_rank,
                   centrality_mu = centrality_mu,
                   centrality_cov = centrality_cov,
                   #age_mu = age_mu,
                   #age_sd = age_sd,
                   node_age = as.integer(nodes_to_include$age_cat),
                   prior_age = rep(1, n_age_cat))

# load model
nodal_regression <- stan_model('models/eigen_regression_motnp.stan') # MOTNP version for long window because single window with ordered categorical predictor

## run model
fit_mpnp_eigen <- sampling(nodal_regression,
                           data = eigen_list,
                           cores = n_chains,
                           chains = n_chains)

## save output
save.image('mpnp_nodalregression/mpnp_long_nodalregression.RData')

### posterior check ####
# load('mpnp_nodalregression/mpnpshort1_nodalregression_conditionaledge_rstan.RData')
## traceplot linear effect size
traceplot(fit_mpnp_eigen, pars = c('beta_age','sigma','predictor[1]','predictor[2]','predictor[3]','predictor[4]','predictor[5]','predictor[6]','predictor[7]','predictor[8]','predictor[9]','predictor[10]'))

traceplot(fit_mpnp_eigen, pars = c('delta[1]','delta[1]','delta[1]','delta[1]','delta[1]','delta[1]'))

traceplot(fit_mpnp_eigen, pars = c('delta_j[1]','delta_j[1]','delta_j[1]','delta_j[1]','delta_j[1]','delta_j[1]','delta_j[1]'))

## posterior predictive check
params <- rstan::extract(fit_mpnp_eigen)
plot(density(centrality_samples_std[1, ]), main="Posterior predictive density of responses (standardised centrality)", col=rgb(0, 0, 0, 0.25), ylim=c(0, 1))
for (i in 1:100) {
  j <- sample(1:(n_chains*n_samples), 1)
  lines(density(centrality_samples_std[j, ]), col=rgb(0, 0, 0, 0.25))
  mu <- rep(NA, n_eles)
  for(k in 1:n_eles){
    mu[k] <- params$beta_age[j]*sum(params$delta_j[j, 1:nodes_to_include$age_cat[k]]) + params$intercept[j]
  }
  sigma <- centrality_cov + diag(rep(params$sigma[j], n_eles))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}

## interpret model
par(mfrow = c(2,2))
beta_diff_12 <- params$delta[, 1] - params$delta[, 2]
plot(density(beta_diff_12), main="Posterior difference:\n<10 and 10-15") ; abline(v=0, lty=2)
beta_diff_23 <- params$delta[, 2] - params$delta[, 3]
plot(density(beta_diff_12), main="Posterior difference:\n10-15 and 15-20") ; abline(v=0, lty=2)
beta_diff_34 <- params$delta[, 3] - params$delta[, 4]
plot(density(beta_diff_23), main="Posterior difference:\n15-20 and 21-25") ; abline(v=0, lty=2)
beta_diff_45 <- params$delta[, 4] - params$delta[, 5]
plot(density(beta_diff_34), main="Posterior difference:\n21-25 and 26-35") ; abline(v=0, lty=2)
beta_diff_56 <- params$delta[, 5] - params$delta[, 6]
plot(density(beta_diff_45), main="Posterior difference:\n26-35 and >35") ; abline(v=0, lty=2)

par(mfrow = c(1,1))
plot(density(beta_diff_12), main="Posterior differences")
lines(density(beta_diff_23), col = 'red')
lines(density(beta_diff_34), col = 'blue')
lines(density(beta_diff_45), col = 'purple')
abline(v=0, lty=2)

plot(density(params$beta_age),
     main = "Posterior age effect estimate")
abline(v = 0, lty = 2)

### predict from model ####
## summarise
(summary <- as.data.frame(round(summary(fit_mpnp_eigen)$summary[1:14, c(1, 4, 8)], 3)))
summary$parameter <- rownames(summary)

## simulate full predictions for model
ages <- as.numeric(unique(nodes_to_include$age_cat))
sim_mean <- matrix(NA, nrow = n_chains*n_samples, ncol = n_age_cat,
                   dimnames = list(1:(n_chains*n_samples), ages))
sim_full <- matrix(NA, nrow = n_chains*n_samples, ncol = n_age_cat,
                   dimnames = list(1:(n_chains*n_samples), ages))
for(i in 1:(n_chains*n_samples)){
  for(j in 1:(n_age_cat)){
    sim_mean[i,j] <- params$beta_age[i]*params$delta[ages[j]] + params$intercept[i]
    sim_full[i,j] <- MASS::mvrnorm(n = 1,
                                   mu = sim_mean[i,j],
                                   Sigma = params$sigma[i])
  }
}

## summarise mean predictions
sum_mean_pred <- data.frame(age = sort(ages),
                            lwr = apply(sim_mean, 2, quantile, probs = 0.025),
                            mean = apply(sim_mean, 2, mean),
                            mid = apply(sim_mean, 2, quantile, probs = 0.5),
                            upr = apply(sim_mean, 2, quantile, probs = 0.975))
plot(sum_mean_pred$mean ~ sum_mean_pred$age, type = 'b', col = 'blue', ylim = c(-5,5))
lines(sum_mean_pred$lwr ~ sum_mean_pred$age, lty = 2)
lines(sum_mean_pred$upr ~ sum_mean_pred$age, lty = 2)

## summarise total predictions
sum_full_pred <- data.frame(age = sort(ages),
                            lwr = apply(sim_full, 2, quantile, probs = 0.025),
                            mean = apply(sim_full, 2, mean),
                            mid = apply(sim_full, 2, quantile, probs = 0.5),
                            upr = apply(sim_full, 2, quantile, probs = 0.975))
lines(sum_full_pred$lwr ~ sum_full_pred$age, lty = 3)
lines(sum_full_pred$upr ~ sum_full_pred$age, lty = 3)

## plot simulations
sim_df <- as.data.frame(sim_full) %>% 
  pivot_longer(cols = everything(), names_to = 'age', values_to = 'eigen_sim')
points(sim_df$eigen_sim ~ sim_df$age, col = rgb(0,0,0,0.01), pch = 19, cex = 0.5)

## plot raw with model output
df_long_ <- df_long %>%
  group_by(node_rank) %>% 
  mutate(mean_eigen = mean(centrality)) %>% 
  ungroup() %>% 
  left_join(nodes_to_include, by = 'node_rank') %>% 
  dplyr::select(-age.x, -age.y)
ggplot()+
  geom_ribbon(data = sum_full_pred,
              aes(x = age, ymin = lwr, ymax = upr),          # shade simulations
              colour = 'transparent', fill = rgb(0,0,0,0.1))+
  geom_ribbon(data = sum_mean_pred,
              aes(x = age, ymin = lwr, ymax = upr),          # shade mean distribution
              colour = 'transparent', 
              fill = rgb(33/255, 145/255, 140/255, 0.5))+
  geom_point(data = df_long, 
             aes(x = as.numeric(age_cat), 
                 y = centrality),                  # all eigenvector draws
             colour = rgb(253/255, 231/255, 37/255, 0.01))+
  geom_point(data = distinct(df_long[,c('node','age_cat','mean_eigen','sightings')]),
             aes(x = as.numeric(age_cat), 
                 y = mean_eigen, 
                 size = sightings),  # mean eigenvector
             colour = rgb(68/255, 1/255, 84/255))+
  geom_line(data = sum_mean_pred,
            aes(x = age, y = mid),                                # mean line
            colour = rgb(33/255, 145/255, 140/255),
            linewidth = 1)+
  # geom_errorbar(data = nodes, aes(xmin = age_lwr, xmax = age_upr,                # age distribution
  #                                   y = mean_eigen, group = node_rank),
  #           colour = rgb(68/255, 1/255, 84/255), linewidth = 0.5, width = 0)+
  scale_x_continuous('age (years)')+
  scale_y_continuous('eigenvector centrality (standardised)')+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22))
ggsave(filename = '../outputs/mpnp_nodalregression_line.png', device = 'png',
       plot = last_plot(), width = 16.6, height = 11.6)

ggplot()+
  geom_violin(data = sim_df, aes(x = as.factor(age), y = eigen_sim),          # shade simulations
              #colour = 'transparent', 
              fill = rgb(0,0,0,0.1))+
  geom_errorbar(data = sum_mean_pred, aes(x = age, ymin = lwr, ymax = upr),               # shade mean distribution
                colour = rgb(33/255, 145/255, 140/255),
                linewidth = 1.5, width = 0.5)+
  geom_line(data = sum_mean_pred, aes(x = age, y = mid),                                # mean line
            colour = rgb(33/255, 145/255, 140/255), linewidth = 1.5)+
  geom_point(data = df_long, aes(x = as.numeric(age_cat_fct), y = centrality),                  # all eigenvector draws
             colour = rgb(253/255, 231/255, 37/255, 0.01))+
  geom_point(data = distinct(df_long[,c('node','age_cat_fct','mean_eigen','sightings')]),
             aes(x = as.numeric(age_cat_fct), y = mean_eigen, size = sightings),  # mean eigenvector
             colour = rgb(68/255, 1/255, 84/255))+
  # geom_errorbar(data = nodes, aes(xmin = age_lwr, xmax = age_upr,                # age distribution
  #                                   y = mean_eigen, group = node_rank),
  #           colour = rgb(68/255, 1/255, 84/255), linewidth = 0.5, width = 0)+
  scale_x_discrete('age category')+
  scale_y_continuous('eigenvector centrality (standardised)')+
  theme_classic()+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22))
ggsave(filename = '../outputs/mpnp_nodalregression_violin.png', device = 'png',
       plot = last_plot(), width = 16.6, height = 11.6)

## save output
save.image('mpnp_nodalregression/mpnp_long_nodalregression.RData')
dev.off()
rm(list = ls()[!ls() %in% c('n_samples','n_age_cat','n_chains')])

# ##### run model with rstan -- time window 1-5 (short) ####
# set.seed(12345)
# 
# for(time_window in 1:5){
#   # load data and remove additional data
#   load(paste0('mpnp_nodalregression/mpnpshort',time_window,'_edgeweights_conditionaledge.RData'))
#   rm(counts_df, adj_mat, edges, ele_obs, obs, summary, x, make_edgelist, plot_network_threshold_mpnp) ; gc()
#   pdf(paste0('../outputs/mpnplong',time_window,'_nodalregression.pdf'))
#   
#   ## add progress marker
#   print(paste0('start window ',time_window,' at ',Sys.time()))
#   
#   ### extract centralities ####
#   ## build adjacency tensor
#   ncol(edge_samples)
#   cdf$node_1_id <- as.integer(as.factor(cdf$node_1))
#   cdf$node_2_id <- as.integer(as.factor(cdf$node_2))+1
#   adj_tensor <- array(0, c(n_chains*n_samples, n_eles, n_eles))
#   for (dyad_id in 1:n_dyads) {
#     dyad_row <- cdf[dyad_id, ]
#     adj_tensor[, dyad_row$node_1_id, dyad_row$node_2_id] <- edge_samples[, dyad_id]
#   }
#   
#   ## calculate centrality and store posterior samples in a matrix
#   centrality_samples <- matrix(0, n_chains*n_samples, n_eles)
#   centrality_samples_std <- matrix(0, n_chains*n_samples, n_eles)
#   for (i in 1:(n_chains*n_samples)) {
#     g <- graph_from_adjacency_matrix(adj_tensor[i, , ], mode="undirected", weighted=TRUE)
#     centrality_samples[i, ] <- eigen_centrality(g)$vector
#     centrality_samples_std[i, ] <- (centrality_samples[i, ] - mean(centrality_samples[i, ]))/sd(centrality_samples[i, ])
#   }
#   head(centrality_samples) # Unstandardised eigenvector centrality
#   head(centrality_samples_std)
#   
#   ## visualise centralities
#   nodes_to_include$node_rank <- as.integer(as.factor(nodes_to_include$node))
#   df_wide <- data.frame(centrality_samples_std)
#   colnames(df_wide) <- 1:n_eles
#   df_long <- pivot_longer(df_wide, cols = everything(),
#                           names_to = "node_rank", values_to = "centrality") %>% 
#     mutate(node_rank = as.integer(node_rank)) %>% 
#     left_join(nodes[,c('node_rank','age')], by = 'node_rank')
#   df_long %>% 
#     mutate(nodes_reordered = fct_reorder(.f = as.factor(node_rank), .x = age, .desc = T)) %>% 
#     ggplot(aes(x = centrality, fill = age)) +
#     geom_density(size = 0.4) +
#     facet_grid(rows = vars(as.factor(nodes_reordered)), scales = "free") +
#     labs(x = "Eigenvector centrality (standardised)") + 
#     theme_void() + 
#     theme(strip.text.y = element_text(size = 12),
#           axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
#           axis.title.x = element_text(size = 12),
#           plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
#   
#   ## add progress marker
#   print('eigen values extracted')
#   
#   ### compute normal approximation ####
#   ## check covariance
#   plot(centrality_samples_std[, 1], centrality_samples_std[, 2],
#        xlab = 'standardised eigenvectors ID1',
#        ylab = 'standardised eigenvectors ID2',
#        las = 1, pch = 19, col = rgb(0,0,1,0.2))
#   
#   ## compute normal approximation
#   centrality_mu <- apply(centrality_samples_std, 2, mean)
#   centrality_cov <- cov(centrality_samples_std)
#   
#   centrality_samples_sim <- MASS::mvrnorm(1e5, centrality_mu, centrality_cov)
#   
#   plot(density(centrality_samples_std[, 1]), lwd = 2, main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
#   lines(density(centrality_samples_sim[, 1]), col = rgb(0,0,1,0.5), lwd = 2)
#   
#   ## add progress marker
#   print('normal approximation computed')
#   
#   ### run model ####
#   ## create data list
#   eigen_list <- list(num_nodes = n_eles,
#                      nodes = nodes_to_include$node_rank,
#                      centrality_mu = centrality_mu,
#                      centrality_cov = centrality_cov,
#                      node_age = nodes_to_include$age)
#   
#   # load model
#   nodal_regression <- stan_model('models/eigen_regression.stan')
#   
#   ## run model
#   fit_mpnp_eigen <- sampling(nodal_regression,
#                             data = eigen_list,
#                             cores = n_chains,
#                             chains = n_chains)
#   
#   ## save output
#   rm(g, edge_samples, adj_tensor, i) ; gc()
#   save.image(paste0('mpnp_nodalregression/mpnplong',time_window,'_nodalregression_conditionaledge_rstan.RData'))
#   
#   ## add progress marker
#   print('model run')
#   
#   ### posterior check ####
#   # load(paste0('mpnp_nodalregression/mpnplong',time_window,'_nodalregression_conditionaledge_rstan.RData'))
#   ## traceplot linear effect size
#   traceplot(fit_mpnp_eigen, pars = c('beta_age','sigma','predictor[1]','predictor[2]','predictor[3]','predictor[4]','predictor[5]','predictor[6]','predictor[7]','predictor[8]','predictor[9]','predictor[10]'))
#   
#   ## posterior predictive check
#   params <- rstan::extract(fit_mpnp_eigen)
#   plot(density(centrality_samples_std[1, ]), main="Posterior predictive density of responses (standardised centrality)", col=rgb(0, 0, 0, 0.25), ylim=c(0, 1))
#   for (i in 1:100) {
#     j <- sample(1:(n_chains*n_samples), 1)
#     lines(density(centrality_samples_std[j, ]), col=rgb(0, 0, 0, 0.25))
#     mu <- params$beta_age[j]*eigen_list$node_age
#     sigma <- centrality_cov + diag(rep(params$sigma[j], n_eles))
#     lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
#   }
#   
#   ## interpret model
#   plot(density(params$beta_age), # don't need to do the contrast because continuous -- is there an equivalent using the do operator? by standardising have we already dealt with the logit transformation?
#        main = "Posterior difference between node types")
#   abline(v = 0, lty = 2)
#   
#   ## add progress marker
#   print('posterior predictive check complete')
#   
#   ### predict from model ####
#   (summary <- as.data.frame(round(summary(fit_mpnp_eigen)$summary[1:2, c(1, 4, 8)], 3)))
#   summary$parameter <- rownames(summary)
#   
#   ## calculate mean predictions for model
#   mod_mu <- data.frame(age = min(nodes_to_include$age):max(nodes_to_include$age)) %>% 
#     mutate(lwr = age*summary$`2.5%`[1],
#            mid = age*summary$mean[1],
#            upr = age*summary$`97.5%`[1])
#   
#   ## plot mean predictions
#   plot(mid ~ age, data = mod_mu, type = 'l', las = 1, col = 'blue', lwd = 2,
#        ylim = c(min(df_long$centrality), max(df_long$centrality)),
#        xlab = 'age (years)', ylab = 'eigenvector centrality (standardised)')
#   lines(lwr ~ age, data = mod_mu, lty = 2)
#   lines(upr ~ age, data = mod_mu, lty = 2)
#   
#   ## simulate full predictions for model
#   sim <- matrix(NA, nrow = n_chains*n_samples, ncol = nrow(mod_mu),
#                 dimnames = list(1:(n_chains*n_samples), mod_mu$age))
#   for(i in 1:nrow(sim)){
#     for(j in 1:ncol(sim)){
#       #sim[i,j] <- params$beta_age[i]*x$age[j]
#       sim[i,j] <- MASS::mvrnorm(n = 1, mu = params$beta_age[i]*mod_mu$age[j],
#                                 Sigma = params$sigma[i])
#     }
#   }
#   
#   ## plot simulations
#   sim_df <- as.data.frame(sim) %>% 
#     pivot_longer(cols = everything(), names_to = 'age', values_to = 'eigen_sim')
#   points(sim_df$eigen_sim ~ sim_df$age, col = rgb(0,0,0,0.01), pch = 19, cex = 0.5)
#   
#   ## summarise simulations
#   sim_summary <- data.frame(age = as.numeric(unique(sim_df$age)),
#                             lwr = NA, mid = NA, upr = NA)
#   for(i in 1:nrow(sim_summary)){
#     x <- sim_df %>% filter(age == sim_summary$age[i])
#     sim_summary$lwr[i] <- quantile(x$eigen_sim, 0.025)
#     sim_summary$mid[i] <- quantile(x$eigen_sim, 0.5)
#     sim_summary$upr[i] <- quantile(x$eigen_sim, 0.975)
#   }
#   
#   ## plot raw with model output
#   df_long <- df_long %>%
#     group_by(node_rank) %>% 
#     mutate(mean_eigen = mean(centrality)) %>% 
#     ungroup() %>% 
#     left_join(nodes[,c('node_rank','sightings')], by = 'node_rank')
#   ggplot()+
#     geom_ribbon(data = sim_summary, aes(x = age, ymin = lwr, ymax = upr),       # shade simulations
#                 colour = 'transparent', fill = rgb(0,0,0,0.1))+
#     geom_ribbon(data = mod_mu, aes(x = age, ymin = lwr, ymax = upr),            # shade mean distribution
#                 colour = 'transparent', fill = rgb(33/255, 145/255, 140/255, 0.5))+
#     geom_point(data = df_long, aes(x = age, y = centrality),                    # all eigenvector draws
#                colour = rgb(253/255, 231/255, 37/255, 0.01))+
#     geom_point(data = df_long, aes(x = age, y = mean_eigen, size = sightings),  # mean eigenvector
#                colour = rgb(68/255, 1/255, 84/255))+
#     geom_line(data = mod_mu, aes(x = age, y = mid),                             # mean line
#               colour = rgb(33/255, 145/255, 140/255), linewidth = 1)+
#     scale_x_continuous('age (years)')+
#     scale_y_continuous('eigenvector centrality (standardised)')+
#     theme_classic()+
#     theme(axis.text = element_text(size = 18),
#           axis.title = element_text(size = 22),
#           legend.text = element_text(size = 18),
#           legend.title = element_text(size = 22))
#   ggsave(filename = paste0('../outputs/mpnplong',time_window,'_nodalregression.png'), device = 'png',
#          plot = last_plot())
#   
#   ## save output ####
#   save.image(paste0('mpnp_nodalregression/mpnplong',time_window,'_nodalregression_conditionaledge_rstan.RData'))
#   dev.off()
#   rm(list = ls()[! ls() %in% c('time_window','nodal_regression')])
#   
#   ## add progress marker
#   print(paste0('time window ',time_window,' complete'))
#   
# }
# 
