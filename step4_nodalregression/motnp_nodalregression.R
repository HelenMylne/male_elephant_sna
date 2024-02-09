##### information #####
# 1) to estimate the ages of individuals in other populations, first came up with a survival curve for the Amboseli elephants which have good birth and death data.
# 2) selected the Gompertz bathtub distribution
# 3) use distribution to convert MOTNP ages in categories to probability distributions of true age (first part of this script)
# 4) extract node centrality distribution from draws of edge weight
# 5) use MOTNP age distributions to predict network centrality (third part of this script)

##### set up ####
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

## set up pdf
pdf('../outputs/motnp_nodalregression.pdf')

## set up theme for plots
theme_set(theme_classic())

##### read in data ####
# load edge weight model and data frames
load('motnp_edgeweights_conditionalprior.RData')
rm(edgelist, x, make_edgelist, plot_network_threshold, i, motnp_ages) ; gc()

# df_nodal <- distinct(counts_df[,c('node_1_males','id_1')])
# colnames(df_nodal) <- c('node_2_males','id_2')
# df_nodal <- rbind(df_nodal, counts_df[nrow(counts_df),c('node_2_males','id_2')])
# colnames(df_nodal) <- c('node','id')
# 
# motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>% 
#   dplyr::select(sort(unique(c(counts_df$id_1, counts_df$id_2)))) %>% 
#   pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')
# motnp_ages <- left_join(motnp_ages, df_nodal, by = 'id')
# motnp_ages$draw <- rep(1:8000, length(unique(motnp_ages$id)))
# 
# mean_motnp_ages <- df_nodal
# mean_motnp_ages$age <- NA
# for(i in 1:nrow(mean_motnp_ages)){
#   x <- motnp_ages[motnp_ages$id == mean_motnp_ages$id[i],]
#   mean_motnp_ages$age[i] <- mean(x$age)
#   rm(x)
# }

## convert age data back to categories
nodes <- nodes %>% 
  rename(id_1 = id) %>% 
  left_join(distinct(counts_df[,c('id_1','age_category_1','age_cat_id_1')]),
            by = 'id_1') %>% 
  rename(id = id_1, age_cat_chr = age_category_1, age_cat_num = age_cat_id_1)
nodes$age_cat_chr[nrow(nodes)] <- counts_df$age_category_2[nrow(counts_df)]
nodes$age_cat_num[nrow(nodes)] <- counts_df$age_cat_id_2[nrow(counts_df)]
nodes$age_cat_num <- as.numeric(nodes$age_cat_num)
nodes$age_cat_fct <- as.factor(nodes$age_cat_num - 2)

## clean up data frame
ele_ids <- nodes$id
n_eles <- length(ele_ids)
edges$chain <- ifelse(edges$chain == 'chain1', 1,
                      ifelse(edges$chain == 'chain2', 2,
                             ifelse(edges$chain == 'chain3', 3, 4)))
edges$draw_id <- edges$position + (edges$chain-1) * 1000
edges <- edges %>% 
  rename(dyad_id = dyad) %>% 
  left_join(counts_df[,c('dyad_id','node_1','node_2','id_1','id_2',
                         'count_1','count_2','event_count','count_dyad')],
            by = 'dyad_id')

## build adjacency tensor
adj_tensor <- array(NA, c(n_eles, n_eles, n_samples*n_chains),
                    dimnames = list(ele_ids, ele_ids, NULL))
for(i in 1:n_dyads) {
  dyad_row <- counts_df[i,]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- edge_samples[,i]
}
adj_tensor[,,1]

#### extract centralities ####
## calculate centrality and store posterior samples in a matrix
centrality_samples <- matrix(0, n_chains*n_samples, n_eles, dimnames = list(NULL,nodes$id))
centrality_samples_std <- matrix(0, n_chains*n_samples, n_eles, dimnames = list(NULL,nodes$id))
for (draw in 1:(n_chains*n_samples)) {
  g <- graph_from_adjacency_matrix(adj_tensor[,,draw], mode = "undirected",
                                   diag = FALSE, # TRUE in Jordan's example and in ANP -- shouldn't make a difference as NA anyway
                                   weighted = TRUE)
  centrality_samples[draw, ] <- eigen_centrality(g)$vector
  centrality_samples_std[draw, ] <- (centrality_samples[draw, ] - mean(centrality_samples[draw, ]))/sd(centrality_samples[draw, ])
}
head(centrality_samples)      # unstandardised eigenvector centrality
head(centrality_samples_std)  #  standardised  eigenvector centrality

## visualise centralities
nodes$node_rank <- as.integer(as.factor(nodes$node))
df_wide <- data.frame(centrality_samples_std)
colnames(df_wide) <- 1:n_eles
df_long <- pivot_longer(df_wide, cols = everything(),
                        names_to = "node_rank", values_to = "centrality") %>% 
  mutate(node_rank = as.integer(node_rank)) %>% 
  left_join(nodes[,c('node_rank','age')], by = 'node_rank')
df_long %>% 
  filter(node_rank <= 30) %>% 
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

#### compute normal approximation ####
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

save.image('motnp_nodalregression.RData')

##### prior predictive check ####
## prior predictive check
n <- 100
n_age_cat <- length(unique(nodes$age_cat_fct))
beta_age <- rnorm(n, 0, 1)
intercept  <- rnorm(n, 0, 1)
age_dirichlet <- rdirichlet(n, c(1,1,1,1,1))
plot(NULL, las = 1, xlab = 'age (standardised)', ylab = 'eigenvector (standardised)',
     ylim = c(min(centrality_mu)-1, max(centrality_mu)+1), xlim = c(1,n_age_cat))
abline(h = min(centrality_mu), lty = 2) ; abline(h = max(centrality_mu), lty = 2)
x <- 1:n_age_cat
for(i in 1:n){
  y <- rep(NA, length(x))
  for(j in 1:length(x)){
    y[j] <- intercept[i] + beta_age[i]*sum(age_dirichlet[i,][1:x[j]])
  }
  lines(x = x, y = y, col = rgb(0,0,1,0.4))
}
rm(n, beta_age, intercept, age_dirichlet, sigma, x, y, df_plot, df_wide) ; gc()

#### run model -- using ordered categorical exposure ####
#load('motnp_nodalregression.RData')

# ## extract age distributions (using a normal for now)
# motnp_ages <- readRDS('../data_processed/step2_ageestimation/motnp_ageestimates_mcmcoutput.rds')
# motnp_ages <- motnp_ages[,which(colnames(motnp_ages) %in% ele_ids)]
# age_mu <- apply(motnp_ages, 2, mean)
# age_sd <- apply(motnp_ages, 2, sd)
# 
## create data list
eigen_list <- list(num_nodes = n_eles,
                   num_age_cat = n_age_cat,
                   length_dirichlet = n_age_cat+1,
                   #nodes = nodes$node_rank,
                   centrality_mu = centrality_mu,
                   centrality_cov = centrality_cov,
                   #age_mu = age_mu,
                   #age_sd = age_sd)
                   node_age = as.integer(nodes$age_cat_fct),
                   prior_age = c(1,1,1,1,1))

## check inputs
boxplot(centrality_mu ~ nodes$age_cat_fct, notch = T)

## load model
nodal_regression <- stan_model('models/eigen_regression_motnp.stan')

## run model
fit_motnp_eigen <- sampling(nodal_regression,
                            data = eigen_list,
                            cores = n_chains,
                            chains = n_chains)

## save output (get it saved, then clean it up once it hasn't crashed, then save the cleaner version!)
save.image('motnp_nodalregression.RData')
rm(g, edge_samples, adj_tensor, i, to_plot, draw,summary,centrality_samples_sim, dyad_row) ; gc()
save.image('motnp_nodalregression.RData')

#### posterior check ####
# load('motnp_nodalregression.RData')
## traceplot linear effect size
traceplot(fit_motnp_eigen, pars = c('beta_age','sigma','predictor[1]','predictor[2]','predictor[3]','predictor[4]','predictor[5]','predictor[6]','predictor[7]','predictor[8]','predictor[9]','predictor[10]'))

## posterior predictive check
params <- rstan::extract(fit_motnp_eigen)
plot(density(centrality_samples_std[1, ]), main="Posterior predictive density of responses (standardised centrality)", col=rgb(0, 0, 0, 0.25), ylim=c(0, 0.6))
for (i in 1:100) {
  j <- sample(1:(n_chains*n_samples), 1)
  lines(density(centrality_samples_std[j, ]), col=rgb(0, 0, 0, 0.25))
  mu <- rep(NA, n_eles)
  for(k in 1:n_eles){
    mu[k] <- params$beta_age[j]*sum(params$delta_j[j,1:nodes$age_cat_fct[k]]) + params$intercept[j]
  }
  sigma <- centrality_cov + diag(rep(params$sigma[j], n_eles))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}

## interpret model
par(mfrow = c(2,2))
beta_diff_12 <- params$delta[, 1] - params$delta[, 2]
plot(density(beta_diff_12), main="Posterior difference:\n10-15 and 15-20") ; abline(v=0, lty=2)
beta_diff_23 <- params$delta[, 2] - params$delta[, 3]
plot(density(beta_diff_23), main="Posterior difference:\n15-20 and 21-25") ; abline(v=0, lty=2)
beta_diff_34 <- params$delta[, 3] - params$delta[, 4]
plot(density(beta_diff_34), main="Posterior difference:\n21-25 and 26-40") ; abline(v=0, lty=2)
beta_diff_45 <- params$delta[, 4] - params$delta[, 5]
plot(density(beta_diff_45), main="Posterior difference:\n26-40 and 40+") ; abline(v=0, lty=2)

par(mfrow = c(1,1))
plot(density(beta_diff_12), main="Posterior differences")
lines(density(beta_diff_23), col = 'red')
lines(density(beta_diff_34), col = 'blue')
lines(density(beta_diff_45), col = 'purple')
abline(v=0, lty=2)

plot(density(params$beta_age),
     main = "Posterior age effect estimate")
abline(v = 0, lty = 2)

# #### predict from model -- dummy ####
# ## summarise
# (summary <- as.data.frame(round(summary(fit_motnp_eigen)$summary[1:14, c(1, 4, 8)], 3)))
# summary$parameter <- rownames(summary)
# 
# ## simulate full predictions for model
# ages <- as.numeric(unique(nodes$age_cat_fct))
# sim_mean <- matrix(NA, nrow = n_chains*n_samples, ncol = n_age_cat,
#                    dimnames = list(1:(n_chains*n_samples), ages))
# sim_full <- matrix(NA, nrow = n_chains*n_samples, ncol = n_age_cat,
#                    dimnames = list(1:(n_chains*n_samples), ages))
# for(i in 1:(n_chains*n_samples)){
#   for(j in 1:(n_age_cat)){
#     sim_mean[i,j] <- params$beta_age[i]*params$delta[ages[j]] + params$intercept[i]
#     sim_full[i,j] <- MASS::mvrnorm(n = 1,
#                               mu = sim_mean[i,j],
#                               Sigma = params$sigma[i])
#   }
# }
# 
# ## summarise mean predictions
# sum_mean_pred <- data.frame(age = sort(as.numeric(unique(nodes$age_cat_fct))),
#                        lwr = apply(sim_mean, 2, quantile, probs = 0.025),
#                        mean = apply(sim_mean, 2, mean),
#                        mid = apply(sim_mean, 2, quantile, probs = 0.5),
#                        upr = apply(sim_mean, 2, quantile, probs = 0.975))
# plot(sum_mean_pred$mean ~ sum_mean_pred$age, type = 'b', col = 'blue', ylim = c(-5,5))
# lines(sum_mean_pred$lwr ~ sum_mean_pred$age, lty = 2)
# lines(sum_mean_pred$upr ~ sum_mean_pred$age, lty = 2)
# 
# ## summarise total predictions
# sum_full_pred <- data.frame(age = sort(as.numeric(unique(nodes$age_cat_fct))),
#                             lwr = apply(sim_full, 2, quantile, probs = 0.025),
#                             mean = apply(sim_full, 2, mean),
#                             mid = apply(sim_full, 2, quantile, probs = 0.5),
#                             upr = apply(sim_full, 2, quantile, probs = 0.975))
# lines(sum_full_pred$lwr ~ sum_full_pred$age, lty = 3)
# lines(sum_full_pred$upr ~ sum_full_pred$age, lty = 3)
# 
# ## plot simulations
# sim_df <- as.data.frame(sim_full) %>% 
#   pivot_longer(cols = everything(), names_to = 'age', values_to = 'eigen_sim')
# points(sim_df$eigen_sim ~ sim_df$age, col = rgb(0,0,0,0.01), pch = 19, cex = 0.5)
# 
# ## plot raw with model output
# df_long_mean <- df_long %>%
#   group_by(node_rank) %>% 
#   mutate(mean_eigen = mean(centrality)) %>% 
#   ungroup() %>% 
#   left_join(nodes, by = 'node_rank') %>% 
#   dplyr::select(age_cat_fct, mean_eigen, sightings, node) %>% 
#   distinct()
# df_long_plot <- df_long %>%
#   left_join(nodes, by = 'node_rank') %>% 
#   dplyr::select(-age.x, -age.y)
# ggplot()+
#   geom_ribbon(data = sum_full_pred,
#               aes(x = age, ymin = lwr, ymax = upr),          # shade simulations
#               colour = 'transparent', fill = rgb(0,0,0,0.1))+
#   geom_ribbon(data = sum_mean_pred,
#               aes(x = age, ymin = lwr, ymax = upr),          # shade mean distribution
#               colour = 'transparent', 
#               fill = rgb(33/255, 145/255, 140/255, 0.5))+
#   geom_point(data = df_long_plot, 
#              aes(x = as.numeric(age_cat_fct), 
#                  y = centrality),                  # all eigenvector draws
#              colour = rgb(253/255, 231/255, 37/255, 0.01))+
#   geom_point(data = df_long_mean,
#              aes(x = as.numeric(age_cat_fct), 
#                  y = mean_eigen, 
#                  size = sightings),  # mean eigenvector
#              colour = rgb(68/255, 1/255, 84/255))+
#   geom_line(data = sum_mean_pred,
#             aes(x = age, y = mid),                                # mean line
#             colour = rgb(33/255, 145/255, 140/255),
#             linewidth = 1)+
#   # geom_errorbar(data = nodes, aes(xmin = age_lwr, xmax = age_upr,                # age distribution
#   #                                   y = mean_eigen, group = node_rank),
#   #           colour = rgb(68/255, 1/255, 84/255), linewidth = 0.5, width = 0)+
#   scale_x_continuous('age (years)')+
#   scale_y_continuous('eigenvector centrality (standardised)')+
#   theme(axis.text = element_text(size = 18),
#         axis.title = element_text(size = 22),
#         legend.text = element_text(size = 18),
#         legend.title = element_text(size = 22))
# ggsave(filename = '../outputs/motnp_nodalregression_line.png', device = 'png',
#        plot = last_plot(), width = 16.6, height = 11.6)
# 
# ggplot()+
#   geom_violin(data = sim_df, aes(x = as.factor(age), y = eigen_sim),          # shade simulations
#               #colour = 'transparent', 
#               fill = rgb(0,0,0,0.1))+
#   geom_errorbar(data = sum_mean_pred, aes(x = age, ymin = lwr, ymax = upr),               # shade mean distribution
#               colour = rgb(33/255, 145/255, 140/255),
#               linewidth = 1.5, width = 0.5)+
#   geom_line(data = sum_mean_pred,
#             aes(x = age,
#                 y = mid),                                # mean line
#             colour = rgb(33/255, 145/255, 140/255), linewidth = 1.5)+
#   geom_point(data = df_long_plot,
#              aes(x = as.numeric(age_cat_fct),
#                  y = centrality),                  # all eigenvector draws
#              colour = rgb(253/255, 231/255, 37/255, 0.01))+
#   geom_point(data = df_long_mean,
#              aes(x = as.numeric(age_cat_fct),
#                  y = mean_eigen,
#                  size = sightings),  # mean eigenvector
#              colour = rgb(68/255, 1/255, 84/255))+
#   # geom_errorbar(data = nodes, aes(xmin = age_lwr, xmax = age_upr,                # age distribution
#   #                                   y = mean_eigen, group = node_rank),
#   #           colour = rgb(68/255, 1/255, 84/255), linewidth = 0.5, width = 0)+
#   scale_x_discrete('age category')+
#   scale_y_continuous('eigenvector centrality (standardised)')+
#   theme_classic()+
#   theme(axis.text = element_text(size = 18),
#         axis.title = element_text(size = 22),
#         legend.text = element_text(size = 18),
#         legend.title = element_text(size = 22))
# ggsave(filename = '../outputs/motnp_nodalregression_violin.png', device = 'png',
#        plot = last_plot(), width = 16.6, height = 11.6)
# 
# ## save output
# save.image('motnp_nodalregression.RData')
# dev.off()
# 
#### predict from model -- raw data ####
# load('motnp_nodalregression.RData')

## summarise
(summary <- as.data.frame(round(summary(fit_motnp_eigen)$summary[1:14, c(1, 4, 8)], 3)))
summary$parameter <- rownames(summary)

## simulate full predictions for model
sim_mean <- matrix(NA, nrow = n_chains*n_samples, ncol = n_eles,
                   dimnames = list(1:(n_chains*n_samples), nodes$id))
sim_full <- matrix(NA, nrow = n_chains*n_samples, ncol = n_eles,
                   dimnames = list(1:(n_chains*n_samples), nodes$id))
for(i in 1:(n_chains*n_samples)){
  for(j in 1:n_eles){
    sim_mean[i,j] <- params$beta_age[i]*params$delta[nodes$age_cat_fct[j]] + params$intercept[i]
    sim_full[i,j] <- MASS::mvrnorm(n = 1,
                                   mu = sim_mean[i,j],
                                   Sigma = params$sigma[i])
  }
}
boxplot(sim_full[,1])

## summarise mean predictions -- all the same for all nodes in the same category because the only predictor
nodes_mean <- sim_mean %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'node_mean') %>% 
  mutate(draw = rep(1:n_eles, each = nrow(sim_mean))) %>% 
  group_by(id) %>% 
  left_join(nodes[,c('id','node','age_cat_chr','age_cat_fct')]) %>% 
  mutate(pred_mu_lwr = quantile(node_mean, probs = 0.025),
         pred_mu_mean = mean(node_mean),
         pred_mu_mid = quantile(node_mean, probs = 0.5),
         pred_mu_upr = quantile(node_mean, probs = 0.975))
ggplot(nodes_mean)+
  geom_boxplot(aes(x = age_cat_fct, y = pred_mu_mean))

## summarise total predictions
nodes_full <- sim_full %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'node_pred') %>% 
  mutate(draw = rep(1:n_eles, each = nrow(sim_full))) %>% 
  group_by(id) %>% 
  left_join(nodes[,c('id','node','age_cat_chr','age_cat_fct')]) %>% 
  mutate(pred_full_lwr = quantile(node_pred, probs = 0.025),
         pred_full_mean = mean(node_pred),
         pred_full_mid = quantile(node_pred, probs = 0.5),
         pred_full_upr = quantile(node_pred, probs = 0.975))
ggplot()+
  geom_violin(data = nodes_full,
              aes(x = age_cat_fct, y = node_pred,
                  fill = factor(age_cat_chr, levels = c('10-15', '15-19',
                                                        '20-25','25-40','40+'))),
              alpha = 0.5)+
  geom_boxplot(data = nodes_full,
               aes(x = age_cat_fct, y = pred_full_mean, fill = age_cat_chr))+
  #geom_boxplot(data = nodes_mean, aes(x = age_cat_fct, y = pred_mu_mean), colour = 'white')+
  scale_fill_viridis_d()+
  labs(fill = 'age category', x = 'age category', y = 'predicted eigenvector centrality')

## plot raw with model output -- GOT TO HERE WITH ATTEMPTING TO FIX PREDICTIONS
df_long_mean <- df_long %>%
  group_by(node_rank) %>% 
  mutate(mean_eigen = mean(centrality)) %>% 
  ungroup() %>% 
  left_join(nodes, by = 'node_rank') %>% 
  dplyr::select(age_cat_fct, mean_eigen, sightings, node) %>% 
  distinct()
df_long_plot <- df_long %>%
  left_join(nodes, by = 'node_rank') %>% 
  dplyr::select(-age.x, -age.y)
ggplot()+
  geom_ribbon(data = nodes_full,
              aes(x = as.numeric(age_cat_fct),
                  ymin = pred_full_lwr, ymax = pred_full_upr),  # shade simulations
              colour = 'transparent', fill = rgb(0,0,0,0.1))+
  geom_ribbon(data = nodes_mean,
              aes(x = as.numeric(age_cat_fct),
                  ymin = pred_mu_lwr, ymax = pred_mu_upr),      # shade mean distribution
              colour = 'transparent', 
              fill = rgb(33/255, 145/255, 140/255, 0.5))+
  geom_point(data = df_long_plot, 
             aes(x = as.numeric(age_cat_fct), 
                 y = centrality),                               # all eigenvector draws
             colour = rgb(253/255, 231/255, 37/255, 0.01))+
  geom_point(data = df_long_mean,
             aes(x = as.numeric(age_cat_fct), 
                 y = mean_eigen, 
                 size = sightings),                             # mean eigenvector
             colour = rgb(68/255, 1/255, 84/255))+
  geom_line(data = nodes_mean,
            aes(x = as.numeric(age_cat_fct),
                y = pred_mu_mean),                              # mean line
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
ggsave(filename = '../outputs/motnp_nodalregression_line.png', device = 'png',
       plot = last_plot(), width = 16.6, height = 11.6)

ggplot()+
  geom_violin(data = sim_df, aes(x = as.factor(age), y = eigen_sim),          # shade simulations
              #colour = 'transparent', 
              fill = rgb(0,0,0,0.1))+
  geom_errorbar(data = sum_mean_pred, aes(x = age, ymin = lwr, ymax = upr),               # shade mean distribution
                colour = rgb(33/255, 145/255, 140/255),
                linewidth = 1.5, width = 0.5)+
  geom_line(data = sum_mean_pred,
            aes(x = age,
                y = mid),                                # mean line
            colour = rgb(33/255, 145/255, 140/255), linewidth = 1.5)+
  geom_point(data = df_long_plot,
             aes(x = as.numeric(age_cat_fct),
                 y = centrality),                  # all eigenvector draws
             colour = rgb(253/255, 231/255, 37/255, 0.01))+
  geom_point(data = df_long_mean,
             aes(x = as.numeric(age_cat_fct),
                 y = mean_eigen,
                 size = sightings),  # mean eigenvector
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
ggsave(filename = '../outputs/motnp_nodalregression_violin.png', device = 'png',
       plot = last_plot(), width = 16.6, height = 11.6)

## save output
save.image('motnp_nodalregression.RData')
dev.off()
