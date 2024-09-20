#### information ####
# script takes data input from edge weight estimation for MOTNP population (input = 1000 draws for each of 4 chains per dyad posterior distribution)
# runs through dyadic regression as specified by Jordan Hart in BISoN examples (https://github.com/JHart96/bison_examples/blob/main/examples/dyadic_regression_stan.md)

#### set up ####
# library(StanHeaders) ; library(rstan) ; library(tidyverse) ; library(car) ; library(LaplacesDemon)
library(StanHeaders, lib.loc = '../packages/')         # library(rstan)
library(rstan, lib.loc = '../packages/')         # library(rstan)
library(tidyverse, lib.loc = '../packages/')     # library(tidyverse)
library(car, lib.loc = '../packages/')           # library(car)
library(LaplacesDemon, lib.loc = '../packages/') # library(LaplacesDemon)
theme_set(theme_classic())

pdf('../outputs/motnp_dyadicregression_dataprep.pdf')

load('motnp_edgeweights_conditionalprior.RData')
rm(edge_binary, fit_edges_motnp, edgelist, edges, motnp_ages, x, summary, i, make_edgelist, plot_network_threshold) ; gc()

set.seed(15)

#### create analysis data frame ####
## calculate integer values
n_dyads <- nrow(counts_df)
n_nodes <- length(unique(c(counts_df$id_1, counts_df$id_2)))

## randomise node id
nodes$node_2 <- nodes$node_1 <- nodes$node
nodes$node_random_2 <- nodes$node_random_1 <- sample(x = 1:n_nodes, size = n_nodes, replace = F)

## clean up counts_df
counts_df <- counts_df %>%
  mutate(age_num_1 = as.integer(age_cat_id_1) - 2,
         age_num_2 = as.integer(age_cat_id_2) - 2) %>%
  mutate(age_diff = abs(age_num_1 - age_num_2) + 1) %>%
  left_join(nodes[c('node_1','node_random_1')],
            by = c('node_1')) %>%
  left_join(nodes[c('node_2','node_random_2')],
            by = c('node_2')) %>%
  dplyr::select('dyad_id','dyad_males',
                'node_1','node_2','node_random_1','node_random_2',
                'id_1','id_2','id_pad_1','id_pad_2',
                'age_category_1','age_category_2',
                'age_num_1','age_num_2','age_diff',
                'count_1','count_2',
                'count_dyad','event_count','apart')

## print progress marker
print('data prepped')

#### summarise edges ####
counts_df$mean_edge <- apply(edge_samples, 2, mean)

## plot raw data
ggplot()+
  geom_jitter(data = counts_df,
              aes(x = age_diff, y = mean_edge),
              shape = 19, alpha = 0.2)+
  scale_x_continuous('difference in age categories')+
  scale_y_continuous('mean estimated edge weight')

edges <- edge_samples %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = 'dyad_id', values_to = 'edge_draw') %>%
  mutate(dyad_id = as.numeric(dyad_id)) %>%
  left_join(counts_df, by = 'dyad_id')
ggplot()+
  geom_point(data = edges,
             aes(x = age_diff, y = edge_draw),
             shape = 19, alpha = 0.05)+
  geom_point(data = counts_df, aes(x = age_diff, y = mean_edge),
             shape = 19, colour = 'white', size = 1)+
  scale_x_continuous('difference in dyad age categories')+
  scale_y_continuous('mean estimated edge weight')

## print progress marker
print('raw data plotted')

#### fit multivariate Gaussian distribution to output of age model ####
### quantities will be given to Stan model as data to model joint posteriors of edge weight in the regression
logit_edge_draws <- logit(edge_samples)
logit_edge_draws_mu <- apply(logit_edge_draws, 2, mean)
logit_edge_draws_sd <- apply(logit_edge_draws, 2, sd)

#### plot to check how well Gaussian approximation is working -- not particularly brilliant to be honest... ####
### Randomly selecting samples to examine
num_check <- 50
selected_samples <- sample(1:(n_dyads-1), num_check, replace = FALSE)

### Setting grid layout
rows <- floor(sqrt(num_check))
cols <- ceiling(num_check / rows)
par(mfrow=c(rows, cols), mar=c(2,2,2,1))

### plot
for (i in selected_samples) {
  mu <- logit_edge_draws_mu[i]
  sd <- logit_edge_draws_sd[i]

  fitted_values <- rnorm(1e5, mean = mu, sd = sd)

  hist(unlist(logit_edge_draws[,i]), probability = TRUE, las = 1,
       main = paste("Dyad", i), xlab = "Value", breaks = 50)
  lines(density(fitted_values), col="blue", lwd=1.5)
}

for (i in selected_samples) {
  mu <- logit_edge_draws_mu[i]
  sd <- logit_edge_draws_sd[i]

  fitted_values <- rnorm(1e5, mean = mu, sd = sd)

  plot(unlist(logit_edge_draws[,i]), unlist(logit_edge_draws[,i+1]),
       col = rgb(0,0,1,0.05), las = 1, main = paste("cov ", i ,"&",i+1))
}

### reset plot window and clean up
par(mfrow=c(1,1), mai = c(1,1,0.5,0.5))
rm(cols, fitted_values, i, num_check, rows, sd, selected_samples) ; gc()

## save workspace for later
save.image('step5_dyadicregression/motnp_dyadicregression.RData')

## print progress marker
print('Gaussian approximation complete')

#### plot raw data ####
counts_df$raw_mu <- logit_edge_draws_mu

## plot
ggplot(counts_df)+
  geom_violin(aes(x = as.factor(age_diff),
                  y = raw_mu),
              fill = rgb(0.2,0,1,0.4))

#### prior predictive check ####
n <- 100
n_age_diffs <- length(unique(counts_df$age_diff))

## global parameters
# intercept <- rnorm(n,0,2)
beta_age_diff <- rnorm(n,0,3)
diff_dirichlet <- rdirichlet(n, rep(1, n_age_diffs))
diff_dirichlet <- cbind(rep(0, n), diff_dirichlet)
tau_sigma_raw <- rnorm(n,0,1)
tau_sigma <- exp(tau_sigma_raw)

## multimembership mean
mu_mm <- rnorm(n,logit(0.1),1) # taken from Chiyo 2011 as upper end of most association strengths
rand_mm <- rnorm(n,0,1)
tau_mm <- rnorm(n,0,1)
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
x <- 1:n_age_diffs
for(bounds in c(30, 15, 5, 2)){

  plot(NULL, las = 1, xlab = 'difference age categeory',
       ylab = 'logit edge weight (standardised)',
       ylim = c(min(logit_edge_draws_mu)-bounds,
                max(logit_edge_draws_mu)+bounds),#c(-8,10),
       xlim = c(1,n_age_diffs))
  for(i in 1:n){
    y <- rep(NA, length(x))
    nodes_sample <- sample(1:n, size = 2, replace = F)
    for(j in 1:length(x)){
      y[j] <- #intercept[i] +
        beta_age_diff[i]*sum(diff_dirichlet[i,][1:x[j]]) +
        mm_nodes[i,nodes_sample[1]] + mm_nodes[i,nodes_sample[2]]
    }
    lines(x = x, y = y, col = rgb(0,0,1,0.4))
    abline(h = min(logit_edge_draws_mu), #-2,
           lty = 2) ;
    abline(h = max(logit_edge_draws_mu), #3.5,
           lty = 2)
  }
}

rm(n, beta_age_diff, diff_dirichlet, tau_sigma_raw, tau_sigma, mu_mm, rand_mm, tau_mm, node_mean, raw_sigma, node_sigma, mm_nodes, bounds, i, j, x, y, nodes_sample, mu, sd) ; gc()

## print progress marker
print('prior predictive complete')

## save pdf and start new one
dev.off()
pdf('../outputs/step5_dyadicregression/motnp_dyadicregression_checkoutputs.pdf')

#### fit dyadic regression ####
## create data list
dyad_data <- list(
  num_dyads = n_dyads,                        # number of dyads
  num_nodes = n_nodes,                        # number of nodes
  num_age_diffs = n_age_diffs,                # number of unique age categories
  length_dirichlet = n_age_diffs + 1,         # number of unique age categories + 1
  logit_edge_mu = logit_edge_draws_mu,        # means of the logit edge weights
  logit_edge_sd = logit_edge_draws_sd,        # standard deviation of logit edge weights
  age_diff = as.integer(counts_df$age_diff),  # age of younger dyad member
  node_1 = counts_df$node_random_1,           # node IDs for multimembership effects
  node_2 = counts_df$node_random_2,           # node IDs for multimembership effects
  prior_diff = rep(1, n_age_diffs)            # prior for difference age slope
)

## load dyadic regression model
#dyadic_regression <- cmdstan_model('models/dyadic_regression_motnp.stan')
dyadic_regression <- stan_model('models/dyadic_regression_motnp_agediff.stan')
n_chains <- 4
n_samples <- 1000

## fit dyadic regression
fit_dyadreg_motnp <- sampling(dyadic_regression,
                              data = dyad_data,
                              iter = n_samples*3, warmup = n_samples*2,
                              chains = n_chains, cores = n_chains,
                              control = list(adapt_delta = 0.9,
                                             max_treedepth = 15),
                              seed = 15)

## save output
save.image('step5_dyadicregression/motnp_dyadicregression.RData')

#### check outputs ####
# load('step5_dyadicregression/motnp_dyadicregression.RData')

## extract draws
posterior_samples <- fit_dyadreg_motnp %>%
  as.data.frame()
# draws <- fit_dyadreg_motnp$draws(format = 'df')

## extract model fit
s <- summary(fit_dyadreg_motnp)
summary <- s$summary %>%
  as.data.frame() %>%
  filter(is.nan(se_mean) == FALSE)
# summary <- fit_dyadreg_motnp$summary()
summary
par(mfrow = c(2,1))
hist(summary$Rhat,  breaks = 50)
hist(summary$n_eff, breaks = 50)
# par(mfrow = c(3,1))
# hist(summary$rhat, breaks = 50)
# hist(summary$ess_bulk, breaks = 50)
# hist(summary$ess_tail, breaks = 50)
par(mfrow = c(1,1))

rownames(summary)[summary$n_eff < 600]
summary[summary$n_eff < 600,c(1,3,9,10)]

save.image('step5_dyadicregression/motnp_dyadicregression.RData')

## extract dyadic regression slopes
b_diff <- posterior_samples$beta_age_diff
# intercept <- posterior_samples$intercept
sigma <- posterior_samples$tau_sigma
delta_j_diff <- posterior_samples[,c('delta_j_diff[1]','delta_j_diff[2]','delta_j_diff[3]','delta_j_diff[4]','delta_j_diff[5]','delta_j_diff[6]')] ; colnames(delta_j_diff) <- 1:(n_age_diffs+1)
mm_nodes <- posterior_samples[,grep(pattern = 'mm_nodes', x = colnames(posterior_samples))]
parameters <- data.frame(beta_age_diff = b_diff,
                         # intercept = intercept,
                         sigma = sigma) %>%
  mutate(chain = rep(1:4, each = n_samples),
         position = rep(1:n_samples, 4)) %>%
  pivot_longer(cols = c('beta_age_diff','sigma'),#,'intercept'),
               names_to = 'parameter', values_to = 'slope_draw')

## traceplots
ggplot(data = parameters)+
  geom_line(aes(x = position,
                y = slope_draw,
                colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')
delta_j_diff %>% as.data.frame() %>%
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
plot(density(as.numeric(logit_edge_draws[1, ])),
     main = "Posterior predictive density of edge weights:\nblack = measured, colours = predicted (col by chain)",
     ylim = c(0, 0.6), col = rgb(0, 0, 0, 0.25), las = 1)
for (i in 1:100) {
  j <- sample(1:n_samples, 1)#nrow(logit_edge_draws), 1)

  mu_plot <- matrix(NA, ncol = n_dyads, nrow = n_chains)

  for(l in c(1:n_chains)){
    n <- (l*n_samples) - n_samples
    for(k in 1:n_dyads){
      mu_plot[l,k] <- #intercept[j+n] +
        b_diff[j+n]*sum(delta_j_diff[j,(1:dyad_data$age_diff[k])]) +
        mm_nodes[j+n,dyad_data$node_1[k]] + mm_nodes[j+n,dyad_data$node_2[k]]
    }
    sigma_plot <- dyad_data$logit_edge_sd + rep(sigma[j+n], n_dyads)
    norm <- rnorm(n_samples, mu_plot[l,], sigma_plot)
    lines(density(as.numeric(logit_edge_draws[j, ])), col = rgb(0, 0, 0, 0.25))
    lines(density(norm),
          col = ifelse(l == 1, rgb(1, 0, 0, 0.25),
                       ifelse(l == 2, rgb(0, 1, 0, 0.25),
                              ifelse(l == 3, rgb(0, 0, 1, 0.25),
                                     rgb(1, 1, 0, 0.25)))) )
  }
}

rm(i,j,mu_plot,l,n,k,sigma_plot,norm,sample_nodes) ; gc()

## save output
save.image('step5_dyadicregression/motnp_dyadicregression.RData')

## print progress marker
print('outputs checked')

#### predict from raw data ####
dev.off()
pdf('../outputs/step5_dyadicregression/motnp_dyadicregression_predictions.pdf')
#load('step5_dyadicregression/motnp_dyadicregression.RData')

## create predictive data frame
pred <- counts_df %>%
  mutate(dyad_rank = 1:n_dyads,
         age_diff_num = as.numeric(age_diff)) %>%
  relocate(dyad_rank, .after = 'dyad_males') %>% 
  relocate(age_diff_num, .after = 'age_diff')
pred$raw_mu <- logit_edge_draws_mu

## create function for predicting mean
get_mean_predictions <- function(prediction_df, new){
  prediction_df$dyad_rank <- 1:nrow(prediction_df)
  mu <- matrix(NA,
               nrow = n_samples*n_chains,
               ncol = nrow(counts_df),
               dimnames = list(1:(n_samples*n_chains),
                               prediction_df$dyad_rank))
  
  if(is.null(new) == FALSE){
    if(new == TRUE){
      prediction_df <- prediction_df %>% 
        rename(age_num_diff = age_diff_new)
    } else {
      prediction_df <- prediction_df %>%
        rename(age_num_diff = age_diff_org)
    }
  }
  
  for(dyad in 1:ncol(mu)){
    node_effects <- mm_nodes[,paste0('mm_nodes[',prediction_df$node_random_1[dyad],']')] + 
      mm_nodes[,paste0('mm_nodes[',prediction_df$node_random_2[dyad],']')]
    for(draw in 1:nrow(mu)){
      mu[draw, dyad] <- #intercept[draw] + 
        b_diff[draw]*sum(delta_j_diff[draw,(1:prediction_df$age_num_diff[dyad])]) + node_effects[draw]
    }
  }
  return(mu)
}

## calculate mean predictions
pred_mu <- get_mean_predictions(prediction_df = pred, new = NULL)
pred$pred_mu <- apply(pred_mu, 2, mean)
pred$pred_mu_lwr <- apply(pred_mu, 2, quantile, prob = 0.025)
pred$pred_mu_upr <- apply(pred_mu, 2, quantile, prob = 0.975)

save.image('step5_dyadicregression/motnp_dyadicregression.RData')

## plot predicted vs original
ggplot(pred)+
  geom_point(aes(x = raw_mu, y = pred_mu,
                 colour = as.factor(age_diff)))+
  geom_abline(slope = 1, intercept = 0)+
  scale_colour_viridis_d()+
  labs(colour = 'age difference between dyad')

## create function for predicting full distribution
get_full_predictions <- function(mu_matrix, sigma, sd){
  predictions <- mu_matrix
  for(dyad in 1:ncol(predictions)){
    for(draw in 1:nrow(predictions)){
      predictions[draw, dyad] <- rnorm(n = 1,
                                       mean = mu_matrix[draw, dyad],
                                       sd = sigma[draw] + sd[dyad])
    }
  }
  return(predictions)
}

## calculate full prediction distribution
pred_full <- get_full_predictions(mu_matrix = pred_mu,
                                  sigma = posterior_samples$tau_sigma,
                                  sd = logit_edge_draws_sd)
pred$pred_full_lwr <- apply(pred_full, 2, quantile, prob = 0.025)
pred$pred_full_upr <- apply(pred_full, 2, quantile, prob = 0.975)
save.image('step5_dyadicregression/motnp_dyadicregression.RData')

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
categorise_age_diff <- function(df){
  if('age_diff_cat' %in% colnames(df)){
    df <- df %>% dplyr::select(-age_diff_cat)
  }
  df <- df %>%
    mutate(age_diff_cat = ifelse(age_diff == 1, 'same category',
                                 ifelse(age_diff == 2, '1 apart',
                                        ifelse(age_diff == 3, '2 apart',
                                               ifelse(age_diff == 4, '3 apart',
                                                      '4 apart'))))) %>%
    mutate(age_diff_cat = factor(age_diff_cat,
                                 levels = c('same category',
                                            '1 apart','2 apart','3 apart',
                                            '4 apart'))) %>%
    relocate(age_diff_cat, .after = age_diff)
  return(df)
}
counts_df <- categorise_age_diff(counts_df)
pred <- categorise_age_diff(pred)
pred_full_df <- categorise_age_diff(pred_full_df)
ggplot()+
  geom_violin(data = pred_full_df,
              mapping = aes(x = age_diff_cat,
                            y = prediction,
                            fill = age_diff_cat),
              #fill = rgb(0.5,0,1,0.4),
              alpha = 0.4,
              colour = 'transparent')+
  geom_boxplot(data = pred,
               mapping = aes(x = age_diff_cat,
                             y = pred_mu),
               width = 0.2,
               position = position_dodge(0.9))+
  # geom_point(data = counts_df,
  #            mapping = aes(x = age_diff,
  #                          # size = total_sightings,
  #                          y = raw_mu),
  #            pch = 21,
  #            position = position_dodge(0.9))+
  geom_violin(data = counts_df,
              mapping = aes(x = age_diff_cat,
                            y = raw_mu),
              fill = rgb(1,1,1,0))+
  scale_fill_viridis_d()+
  labs(x = 'difference in age categories',
       y = 'logit edge weight',
       fill = 'difference in\nage categories')
ggsave(plot = last_plot(), device = 'svg', width = 2400, height = 1800, units = 'px',
       filename = 'motnp_dyadic_logit.svg',
       path = '../outputs/step5_dyadicregression/')
ggsave(plot = last_plot(), device = 'png', width = 2400, height = 1800, units = 'px',
       filename = 'motnp_dyadic_logit.png',
       path = '../outputs/step5_dyadicregression/')

## convert to invlogit scale
pred_mu_invlogit <- invlogit(pred_mu)
pred_full_invlogit <- invlogit(pred_full)
pred$pred_mu_invlogit <- apply(pred_mu_invlogit, 2, mean)
pred$pred_mu_lwr_invlogit <- apply(pred_mu_invlogit, 2, quantile, prob = 0.025)
pred$pred_mu_upr_invlogit <- apply(pred_mu_invlogit, 2, quantile, prob = 0.975)
pred$pred_full_lwr_invlogit <- apply(pred_full_invlogit, 2, quantile, prob = 0.025)
pred$pred_full_upr_invlogit <- apply(pred_full_invlogit, 2, quantile, prob = 0.975)
counts_df$raw_mu_invlogit <- apply(edge_samples, 2, mean)

## convert to long format
pred_full_df_invlogit <- pred_full_invlogit %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(),
               names_to = 'dyad_id',
               values_to = 'prediction') %>%
  mutate(dyad_id = as.integer(dyad_id)) %>%
  left_join(pred, by = 'dyad_id')

## plot on invlogit scale
ggplot()+
  geom_violin(data = pred_full_df_invlogit,
              mapping = aes(x = age_diff_cat,
                            y = prediction),
              #fill = rgb(0.5,0,1,0.4),
              alpha = 0.4,
              colour = 'transparent')+
  geom_boxplot(data = pred,
               mapping = aes(x = age_diff_cat,
                             y = pred_mu_invlogit),
               width = 0.2,
               position = position_dodge(0.9))+
  geom_violin(data = counts_df,
              mapping = aes(x = age_diff_cat,
                            y = raw_mu_invlogit),
              fill = rgb(1,1,1,0))+
  #scale_y_continuous(limits = c(0,0.5))+
  scale_fill_viridis_d()+
  labs(x = 'difference in age categories',
       y = 'edge weight',
       fill = 'difference in\nage categories')
ggsave(plot = last_plot(), device = 'svg', width = 2400, height = 1800, units = 'px',
       filename = 'motnp_dyadic_invlogit.svg',
       path = '../outputs/step5_dyadicregression/')
ggsave(plot = last_plot(), device = 'png', width = 2400, height = 1800, units = 'px',
       filename = 'motnp_dyadic_invlogit.png',
       path = '../outputs/step5_dyadicregression/')

save.image('step5_dyadicregression/motnp_dyadicregression.RData')

## print progress marker
print('predictions complete')

#### extract original slope values ####
load('step5_dyadicregression/motnp_dyadicregression.RData')
dev.off()
pdf('../outputs/step5_dyadicregression/motnp_dyadicregression_contrasts.pdf')

## calculate predictions from altered data frames
pred_new <- pred %>%
  mutate(age_diff_org = as.integer(age_diff),
         age_diff_new = ifelse(age_diff == 5, 1, as.integer(age_diff) + 1))

# ## create function for predicting mean
# get_mean_predictions <- function(prediction_df, new){
#   prediction_df$dyad_rank <- 1:nrow(prediction_df)
#   mu <- matrix(NA,
#                nrow = 100,#n_samples*n_chains,
#                ncol = nrow(counts_df),
#                dimnames = list(1:100,#(n_samples*n_chains),
#                                prediction_df$dyad_rank))
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
#     node_effects <- mm_nodes[,paste0('mm_nodes[',prediction_df$node_random_1[dyad],']')] + 
#       mm_nodes[,paste0('mm_nodes[',prediction_df$node_random_2[dyad],']')]
#     for(draw in 1:nrow(mu)){
#       mu[draw, dyad] <- #intercept[draw] + 
#         b_diff[draw]*sum(delta_j_diff[draw,(1:prediction_df$age_num_diff[dyad])]) + node_effects[draw]
#     }
#   }
#   return(mu)
# }
# 
# ## create function for predicting full distribution
# get_full_predictions <- function(mu_matrix, sigma, sd){
#   predictions <- mu_matrix
#   for(dyad in 1:ncol(predictions)){
#     for(draw in 1:nrow(predictions)){
#       predictions[draw, dyad] <- rnorm(n = 1,
#                                        mean = mu_matrix[draw, dyad],
#                                        sd = sigma[draw] + sd[dyad])
#     }
#   }
#   return(predictions)
# }
# 
# ## add progress marker
# print('functions created')

## calculate mean predictions
pred_mu_org <- get_mean_predictions(prediction_df = pred_new, new = FALSE)
pred_mu_new <- get_mean_predictions(prediction_df = pred_new, new = TRUE)

## calculate full predictions
pred_full_org <- get_full_predictions(mu_matrix = pred_mu_org,
                                      sigma = posterior_samples$tau_sigma,
                                      sd = logit_edge_draws_sd)
pred_full_new <- get_full_predictions(mu_matrix = pred_mu_new,
                                      sigma = posterior_samples$tau_sigma,
                                      sd = logit_edge_draws_sd)

## save workspace
save.image('step5_dyadicregression/motnp_dyadicregression.RData')
# load('step5_dyadicregression/motnp_dyadicregression.RData')
# rm(list = ls()[! ls() %in% c('predictions_full', 'predictions_mean','sim_diff')]) ; gc()

## add progress marker
print('predictions complete')

# ## calculate contrasts -- age category
small_age_gap <- which(pred_new$age_diff_org != 5)
big_age_gap <- which(pred_new$age_diff_org == 5)
contrasts <- pred_full_new[,small_age_gap] - pred_full_org[,small_age_gap]
contrasts5 <- pred_full_new[,big_age_gap] - pred_full_org[,big_age_gap] * (-1)

## add progress marker
print('contrasts data frame populated')

## calculate contrasts
print(paste0('contrast per category: mean = ', round(mean(contrasts),4),
             ', stdv = ', round(sd(contrasts),4),
             ', 95% CI = [', round(quantile(contrasts, prob = 0.025, na.rm = T),4),':',
             round(quantile(contrasts, prob = 0.975, na.rm = T),4),']'))
print(paste0('contrast big vs small: mean = ', round(mean(contrasts5),4),
             ', stdv = ', round(sd(contrasts5),4),
             ', 95% CI = [', round(quantile(contrasts5, prob = 0.025, na.rm = T),4),':',
             round(quantile(contrasts5, prob = 0.975, na.rm = T),4),']'))
# can't do contrast per year because contrast in age difference is meaningless -- shifting from a difference of 1 category to 2 categories could mean anything in terms of actual age difference between them

## add progress marker
print('contrasts calculated')

## plot
contrasts5 <- contrasts5 %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = 'dyad_rank', values_to = 'contrast') %>%
  mutate(dyad_rank = as.numeric(dyad_rank)) %>%
  left_join(pred_new, by = 'dyad_rank')
contrasts_long <- contrasts %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'dyad_rank', values_to = 'contrast') %>% 
  mutate(dyad_rank = as.numeric(dyad_rank)) %>% 
  left_join(pred_new, by = 'dyad_rank')
# contrasts_long <- rbind(contrasts_long, contrasts5)

contrast_summary <- contrasts_long %>% 
  group_by(dyad_rank) %>% 
  summarise(mean_contrast = mean(contrast)) %>% 
  left_join(pred_new[,c('dyad_rank','age_diff_org','age_diff_new')]) %>% 
  mutate(change = paste0(age_diff_org-1, ' -> ', age_diff_new-1))
ggplot(data = contrast_summary)+
  geom_density_ridges(aes(x = mean_contrast, y = change, fill = change))+
  scale_fill_viridis_d()+
  labs(x = 'contrast in logit(probability of associating)',
       y = 'change in age difference',
       fill = 'change in age\ndifference')

## finish up
save.image('step5_dyadicregression/motnp_dyadicregression.RData')
dev.off()
