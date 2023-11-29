##### information ####
# 1) to estimate the ages of individuals in other populations, first came up with a survival curve for the Amboseli elephants which have good birth and death data.
# 2) selected the Gompertz bathtub distribution
# 3) use distribution to convert MOTNP ages in categories to probability distributions of true age (first part of this script)
# 4) extract node centrality distribution from draws of edge weight
# 5) use MOTNP age distributions to predict network centrality (third part of this script)

##### set up ####
#options(future.globals.maxSize = 10000*(1024^2))   # running model with full age distributions = 5.01GB of globals, which exceeds default maximum allowed size. Set to allow up to 10 GB to allow model to run

# library(rstan) ; library(sna) ; library(tidyverse) ; library(LaplacesDemon) ; library(MASS)
library(StanHeaders, lib.loc = '../packages/')     # library(StanHeaders)
library(rstan, lib.loc = '../packages/')           # library(rstan)
library(sna, lib.loc = '../packages/')             # library(sna)
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

##### prior predictive check ####
## define PDF output
#pdf('../outputs/anp_nodalregression_priorcheck.pdf')

# ## simulate
# age <- 10:60
# mu_mean <- 0
# mu_stdv <- 0.1
# mu <- rnorm(100, mu_mean, mu_stdv)
# mu <- sort(mu)
# mu_mtrx <- matrix(NA, nrow = length(mu), ncol = length(age))
# for(i in 1:nrow(mu_mtrx)){
#   mu_mtrx[i,] <- mu[i]*age
# }
# sigma_range <- rexp(25, 2)
# sigma_range <- sort(sigma_range)
# par(mfrow = c(5,5), mai = c(0.2,0.2,0.2,0.2))
# for(j in 1:25){
#   plot(NULL, xlim = c(10,60), ylim = c(0,1), las = 1,
#        xlab = '', ylab = '')
#   sigma <- diag(rep(sigma_range[j], length(age)))
#   for(i in 1:length(mu)){
#   predictor <- mu_mtrx[i,]
#   y <- MASS::mvrnorm(1, predictor, sigma)
#   lines(x = age, y = invlogit(y), col = rgb(0,0,1,1))
#   }
# }
# par(mfrow = c(1,1))
# dev.off()
# rm(list = ls())

##### run model with rstan -- time window 1 ####
## load data and remove additional data
load('anp_nodalregression/anpshort1_nodalregression_conditionaledge.RData')
rm(counts_df, adj_mat, edges, ele_obs, obs, summary, make_edgelist, plot_network_threshold_anp) ; gc()

## set up pdf
pdf('../outputs/anpshort1_nodalregression.pdf')

### extract centralities ####
## build adjacency tensor
ncol(edge_samples)
cdf_1$node_1_id <- as.integer(as.factor(cdf_1$node_1))
cdf_1$node_2_id <- as.integer(as.factor(cdf_1$node_2))+1
adj_tensor <- array(0, c(n_chains*n_samples, n_eles, n_eles))
for (dyad_id in 1:n_dyads) {
  dyad_row <- cdf_1[dyad_id, ]
  adj_tensor[, dyad_row$node_1_id, dyad_row$node_2_id] <- edge_samples[, dyad_id]
  adj_tensor[, dyad_row$node_2_id, dyad_row$node_1_id] <- edge_samples[, dyad_id]
}
adj_tensor[1,,]

## calculate centrality and store posterior samples in a matrix
centrality_samples_invlogit <- matrix(0, n_chains*n_samples, n_eles)
for(i in 1:(n_chains*n_samples)){
  centrality_samples_invlogit[i, ] <- sna::evcent(adj_tensor[i,,], gmode = 'graph')
}
head(centrality_samples_invlogit) # Unstandardised eigenvector centrality, 0-1 bounded

## convert to logit scale
centrality_samples <- logit(centrality_samples_invlogit)
head(centrality_samples)          # unstandardised eigenvector centrality, natural scale
centrality_mu_ustd <- apply(centrality_samples, 2, mean)

## add mean estimate per ID to nodes data frame
df <- data.frame(node_1_id = 1:n_eles,
                 mean_eigen_ustd = centrality_mu_ustd) %>% 
  left_join(distinct(cdf_1[,c('node_1_id','id_1')]), by = 'node_1_id') %>% 
  rename(node_rank = node_1_id,
         id = id_1)
df$id <- ifelse(is.na(df$id) == FALSE, df$id,
                nodes$id[! nodes$id %in% df$id])
nodes <- nodes %>% 
  mutate(node_rank = as.integer(as.factor(node))) %>% 
  left_join(df, by = c('node_rank','id'))
rm(df); gc()

## standardise
centrality_samples_std <- centrality_samples
for(i in 1:nrow(centrality_samples_std)){
  centrality_samples_std[i,] <- (centrality_samples[i,] - mean(centrality_samples[i,])) / sd(centrality_samples[i,])
}
head(centrality_samples_std)      # standardised eigenvector centrality, natural scale

## visualise centralities
df_wide <- data.frame(centrality_samples_std)
colnames(df_wide) <- 1:n_eles
df_long <- pivot_longer(df_wide, cols = everything(),
                        names_to = "node_rank", values_to = "centrality") %>% 
  mutate(node_rank = as.integer(node_rank)) %>% 
  left_join(nodes[,c('node_rank','age')], by = 'node_rank')
df_long %>% 
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

## add mean estimate per ID to nodes data frame
df <- df_long %>% 
  group_by(node_rank) %>% 
  mutate(mean_eigen_std = mean(centrality)) %>% 
  dplyr::select(node_rank, mean_eigen_std) %>% 
  distinct()
nodes <- nodes %>% 
  left_join(df, by = 'node_rank')
rm(df); gc()

### compute normal approximation ####
## check covariance
plot(centrality_samples_std[, 1], centrality_samples_std[, 2],
     xlab = 'standardised eigenvectors ID1',
     ylab = 'standardised eigenvectors ID2',
     las = 1, pch = 19, col = rgb(0,0,1,0.2))

## compute normal approximation
centrality_mu <- apply(centrality_samples_std, 2, mean)
centrality_cov <- cov(centrality_samples_std)

centrality_samples_sim <- MASS::mvrnorm(1e5, centrality_mu, centrality_cov)

plot(density(centrality_samples_std[, 1]), lwd = 2, main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
lines(density(centrality_samples_sim[, 1]), col = rgb(0,0,1,0.5), lwd = 2)

### run model ####
## standardise age variable
nodes$age_std <- (nodes$age - mean(nodes$age))/sd(nodes$age)

## create data list
eigen_list <- list(num_nodes = n_eles,
                   nodes = nodes$node_rank,
                   centrality_mu = centrality_mu,
                   centrality_cov = centrality_cov,
                   node_age = nodes$age_std)

# load model
nodal_regression <- stan_model('models/eigen_regression_intercept_standardiseall.stan') # nodal_regression <- cmdstan_model('models/eigen_regression_intercept_standardiseall.stan')

## run model
fit_anp1_eigen <- sampling(nodal_regression, data = eigen_list,
                           cores = n_chains, chains = n_chains,
                           iter = n_samples*2, warmup = n_samples) # fit_anp1_eigen <- nodal_regression$sample(data = eigen_list, chains = n_chains, parallel_chains = n_chains, iter_warmup = n_samples, iter_sampling = n_samples)

## save output
rm(edge_samples, adj_tensor, i) ; gc()
save.image('anp_nodalregression/anpshort1_nodalregression_conditionaledge_rstan.RData')

### posterior check ####
# load('anp_nodalregression/anpshort1_nodalregression_conditionaledge.RData')
## summarise model
#fit_anp1_eigen$summary()
summary(fit_anp1_eigen)

## extract posterior
params <- rstan::extract(fit_anp1_eigen)

## traceplot linear effect size
traceplot(fit_anp1_eigen, pars = c('beta_age','sigma','predictor[1]','predictor[10]','predictor[20]'))
# params %>% 
#   select(intercept,beta_age,sigma,`predictor[1]`,`predictor[50]`,`predictor[100]`) %>% 
#   pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'draw') %>% 
#   mutate(chain_position = rep(rep(1:n_samples, each = 6), n_chains),
#          chain = rep(1:n_chains, each = 6*n_samples)) %>% 
#   #filter(chain == 4) %>% 
#   ggplot(aes(x = chain_position, y = draw, colour = as.factor(chain)))+
#   geom_line()+
#   facet_wrap(. ~ parameter, scales = 'free_y')+
#   theme_bw()+
#   theme(legend.position = 'none') # chains 1-3 are good, chain 4 poorly mixed after about 600 draws

## check model fit
# summary <- summary(fit_anp1_eigen)
# par(mfrow = c(3,1))
# hist(summary$rhat, breaks = 50)
# hist(summary$ess_bulk, breaks = 50)
# hist(summary$ess_tail, breaks = 50)
# par(mfrow = c(1,1))

## posterior predictive check
plot(density(centrality_samples_std[1, ]), main="Posterior predictive check (standardised centrality):\nblack = data, blue = predicted", col=rgb(0, 0, 0, 0.25))
for (i in 1:100) {
  j <- sample(1:length(params$beta_age), 1)
  lines(density(centrality_samples_std[j, ]), col=rgb(0, 0, 0, 0.25))
  mu <- params$beta_age[j]*eigen_list$node_age + params$intercept[j]
  sigma <- centrality_cov + diag(rep(params$sigma[j], n_eles))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}

## extract original values from output
( fit_slope <- params$beta_age * (sd(centrality_mu_ustd)/sd(nodes$age)) )
( fit_intcp <- mean(nodes$age) - fit_slope * mean(nodes$age) )
plot((nodes$age*mean(fit_slope) + mean(fit_intcp)) ~ nodes$age,
     main = 'mean predictions vs age')

## plot density curves for each
plot(density(fit_slope)) ; abline(v = 0, lty = 2)
plot(density(fit_intcp)) ; abline(v = 0, lty = 2)

## save image
save.image('anp_nodalregression/anpshort1_nodalregression_conditionaledge_rstan.RData')

### predict from model ####
#load('anp_nodalregression/anpshort1_nodalregression_conditionaledge_rstan.RData')
#(summary <- as.data.frame(round(summary(fit_anp1_eigen)$summary[1:2, c(1, 4, 8)], 3)))
#summary$parameter <- rownames(summary)

## calculate mean predictions for model -- THERE'S SOMETHING WRONG WITH MY LOWER AND UPPER PREDICTIVE BOUNDS HERE, AS THEY FLIP WHICH IS LARGER IN A COUPLE OF PLACES, BUT i HAVE NO IDEA WHAT SO COME BACK TO THIS
nodes$mean_predict <- nodes$age_std * mean(params$beta_age) + mean(params$intercept)           # mean prediction = age * mean slope
nodes$lwr_predict <- nodes$age_std * quantile(params$beta_age, probs = 0.025) + quantile(params$intercept, probs = 0.025)
nodes$upr_predict <- nodes$age_std * quantile(params$beta_age, probs = 0.975) + quantile(params$intercept, probs = 0.975)

## compare to raw data
plot(mean_predict ~ age_std, data = nodes, type = 'l', las = 1,
     ylim = c(min(c(nodes$mean_eigen_std,nodes$mean_predict)),
              max(c(nodes$mean_eigen_std,nodes$mean_predict))))
points(nodes$mean_eigen_std ~ nodes$age_std)
polygon(y = c(nodes$lwr_predict, rev(nodes$upr_predict)), 
        x = c(nodes$age_std, rev(nodes$age_std)),
        col = rgb(1,1,0,0.5), border = NA)

## reverse standardisation
nodes$mean_eigen_ustd
(nodes$predict_ustd <- nodes$mean_predict * sd(nodes$mean_eigen_ustd) + mean(nodes$mean_eigen_ustd))
plot(nodes$predict_ustd ~ nodes$mean_eigen_ustd,
     ylim = c(min(nodes$mean_eigen_ustd), max(nodes$mean_eigen_ustd)))
abline(a = 0, b = 1)

## compare unstandardised to raw data -- unstandardising lwr and upr does not work like this!
nodes$lwr_ustd <- nodes$age * quantile(fit_slope, probs = 0.025) + quantile(fit_intcp, probs = 0.025)
nodes$upr_ustd <- nodes$age * quantile(fit_slope, probs = 0.975) + quantile(fit_intcp, probs = 0.975)
plot(nodes$predict_ustd ~ nodes$age, type = 'l', las = 1,
     ylim = c(min(nodes$mean_eigen_ustd), max(nodes$mean_eigen_ustd)))    # plot against age
points(nodes$mean_eigen_ustd ~ nodes$age)                                 # add raw points
polygon(y = c(nodes$lwr_ustd, rev(nodes$upr_ustd)), x = c(nodes$age, rev(nodes$age)),
        col = rgb(1,1,0,0.5), border = NA) # THIS DOES NOT WORK CURRENTLY!

## convert to invlogit scale
nodes$mean_predict_invlogit <- invlogit(nodes$predict_ustd)
nodes$lwr_invlogit <- invlogit(nodes$lwr_ustd)
nodes$upr_invlogit <- invlogit(nodes$upr_ustd)
nodes$mu_invlogit <- invlogit(nodes$mean_eigen_ustd)
plot(nodes$mu_invlogit ~ nodes$age, ylim = c(0,0.3))
lines(nodes$mean_predict_invlogit ~ nodes$age)
polygon(y = c(nodes$lwr_invlogit, rev(nodes$upr_invlogit)), x = c(nodes$age, rev(nodes$age)),
        col = rgb(1,1,0,0.5), border = NA)

## simulate full predictions for model
mod_mu <- data.frame(id = 1:100,
                     age = seq(from = min(nodes$age),
                               to = max(nodes$age),
                               length.out = 100))
sim <- matrix(NA, nrow = n_chains*n_samples, ncol = nrow(mod_mu),
              dimnames = list(1:(n_chains*n_samples), mod_mu$age))
for(i in 1:nrow(sim)){
  for(j in 1:ncol(sim)){
    sim[i,j] <- MASS::mvrnorm(n = 1, mu = params$beta_age[i]*mod_mu$age[j],
                              Sigma = params$sigma[i])
  }
}

## plot simulations
sim_df <- as.data.frame(sim) %>% 
  pivot_longer(cols = everything(), names_to = 'age', values_to = 'eigen_sim')
plot(sim_df$eigen_sim ~ sim_df$age, col = rgb(0,0,0,0.01), pch = 19, cex = 0.5)

## summarise simulations
sim_summary <- data.frame(age = as.numeric(unique(sim_df$age)),
                          lwr = NA, mid = NA, upr = NA)
for(i in 1:nrow(sim_summary)){
  x <- sim_df %>% filter(age == sim_summary$age[i])
  sim_summary$lwr[i] <- quantile(x$eigen_sim, 0.025)
  sim_summary$mid[i] <- quantile(x$eigen_sim, 0.5)
  sim_summary$upr[i] <- quantile(x$eigen_sim, 0.975)
}

## plot raw with model output
df_long <- df_long %>%
  group_by(node_rank) %>% 
  mutate(mean_eigen = mean(centrality)) %>% 
  ungroup() %>% 
  left_join(nodes[,c('node_rank','sightings')], by = 'node_rank')
ggplot()+
  # geom_ribbon(data = sim_summary, aes(x = age, ymin = lwr, ymax = upr),       # shade simulations
  #             colour = 'transparent', fill = rgb(0,0,0,0.1))+
  # geom_ribbon(data = mod_mu, aes(x = age, ymin = lwr, ymax = upr),            # shade mean distribution
  #             colour = 'transparent', fill = rgb(33/255, 145/255, 140/255, 0.5))+
  geom_point(data = df_long, aes(x = age, y = centrality),                    # all eigenvector draws
             colour = rgb(253/255, 231/255, 37/255, 0.01))+
  geom_point(data = df_long, aes(x = age, y = mean_eigen, size = sightings),  # mean eigenvector
             colour = rgb(68/255, 1/255, 84/255))+
  # geom_line(data = mod_mu, aes(x = age, y = mid),                             # mean line
  #           colour = rgb(33/255, 145/255, 140/255), linewidth = 1)+
  scale_x_continuous('age (years)')+
  scale_y_continuous('eigenvector centrality (standardised)')+
  theme_classic()+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22))
ggsave(filename = '../outputs/anpshort1_nodalregression_raw.png', device = 'png',
       plot = last_plot())

## save output
save.image('anp_nodalregression/anpshort1_nodalregression_conditionaledge_predictions.RData')
dev.off()
rm(list = ls()[!ls() %in% 'nodal_regression'])

##### run model with rstan -- time window 2-36 ####
set.seed(12345)

for(time_window in 2:36){
  ### load data ####
  load(paste0('anp_nodalregression/anpshort',time_window,'_nodalregression_conditionaledge.RData'))
  rm(counts_df, adj_mat, edges, ele_obs, obs, summary, make_edgelist, plot_network_threshold_anp) ; gc()
  
  ## set up pdf
  pdf(paste0('../outputs/anpshort',time_window,'_nodalregression.pdf'))
  
  ## add progress marker
  print(paste0('start window ',time_window,' at ',Sys.time()))
  
  ### extract centralities ####
  ## build adjacency tensor
  ncol(edge_samples)
  cdf_1$node_1_id <- as.integer(as.factor(cdf_1$node_1))
  cdf_1$node_2_id <- as.integer(as.factor(cdf_1$node_2))+1
  adj_tensor <- array(0, c(n_chains*n_samples, n_eles, n_eles))
  for (dyad_id in 1:n_dyads) {
    dyad_row <- cdf_1[dyad_id, ]
    adj_tensor[, dyad_row$node_1_id, dyad_row$node_2_id] <- edge_samples[, dyad_id]
    adj_tensor[, dyad_row$node_2_id, dyad_row$node_1_id] <- edge_samples[, dyad_id]
  }
  head(adj_tensor[1,,])
  
  ## calculate centrality and store posterior samples in a matrix
  centrality_samples_invlogit <- matrix(0, n_chains*n_samples, n_eles)
  for(i in 1:(n_chains*n_samples)){
    centrality_samples_invlogit[i, ] <- sna::evcent(adj_tensor[i,,], gmode = 'graph')
  }
  head(centrality_samples_invlogit) # Unstandardised eigenvector centrality, 0-1 bounded
  
  ## convert to logit scale
  centrality_samples <- logit(centrality_samples_invlogit)
  head(centrality_samples)          # unstandardised eigenvector centrality, natural scale
  centrality_mu_ustd <- apply(centrality_samples, 2, mean)
  
  ## add mean estimate per ID to nodes data frame
  df <- data.frame(node_1_id = 1:n_eles,
                   mean_eigen_ustd = centrality_mu_ustd) %>% 
    left_join(distinct(cdf_1[,c('node_1_id','id_1')]), by = 'node_1_id') %>% 
    rename(node_rank = node_1_id,
           id = id_1)
  df$id <- ifelse(is.na(df$id) == FALSE, df$id,
                  nodes$id[! nodes$id %in% df$id])
  nodes <- nodes %>% 
    mutate(node_rank = as.integer(as.factor(node))) %>% 
    left_join(df, by = c('node_rank','id'))
  rm(df); gc()
  
  ## standardise
  centrality_samples_std <- centrality_samples
  for(i in 1:nrow(centrality_samples_std)){
    centrality_samples_std[i,] <- (centrality_samples[i,] - mean(centrality_samples[i,])) / sd(centrality_samples[i,])
  }
  head(centrality_samples_std)      # standardised eigenvector centrality, natural scale
  
  ## visualise centralities
  df_wide <- data.frame(centrality_samples_std)
  colnames(df_wide) <- 1:n_eles
  df_long <- pivot_longer(df_wide, cols = everything(),
                          names_to = "node_rank", values_to = "centrality") %>% 
    mutate(node_rank = as.integer(node_rank)) %>% 
    left_join(nodes[,c('node_rank','age')], by = 'node_rank')
  df_long %>% 
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
  
  ## add mean estimate per ID to nodes data frame
  df <- df_long %>% 
    group_by(node_rank) %>% 
    mutate(mean_eigen_std = mean(centrality)) %>% 
    dplyr::select(node_rank, mean_eigen_std) %>% 
    distinct()
  nodes <- nodes %>% 
    left_join(df, by = 'node_rank')
  rm(df); gc()
  
    ## add progress marker
  print('eigen values extracted')
  
  ### compute normal approximation ####
  ## check covariance
  plot(centrality_samples_std[, 1], centrality_samples_std[, 2],
       xlab = 'standardised eigenvectors ID1',
       ylab = 'standardised eigenvectors ID2',
       las = 1, pch = 19, col = rgb(0,0,1,0.2))
  
  ## compute normal approximation
  centrality_mu <- apply(centrality_samples_std, 2, mean)
  centrality_cov <- cov(centrality_samples_std)
  
  centrality_samples_sim <- MASS::mvrnorm(1e5, centrality_mu, centrality_cov)
  
  plot(density(centrality_samples_std[, 1]), lwd = 2, main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
  lines(density(centrality_samples_sim[, 1]), col = rgb(0,0,1,0.5), lwd = 2)
  
  ## add progress marker
  print('normal approximation computed')

  ### run model ####
  ## standardise age variable
  nodes$age_std <- (nodes$age - mean(nodes$age))/sd(nodes$age)
  
  ## create data list
  eigen_list <- list(num_nodes = n_eles,
                     nodes = nodes$node_rank,
                     centrality_mu = centrality_mu,
                     centrality_cov = centrality_cov,
                     node_age = nodes$age_std)
  
  ## run model
  fit_anp_eigen <- sampling(nodal_regression, data = eigen_list,
                             cores = n_chains, chains = n_chains,
                             iter = n_samples*2, warmup = n_samples)
  # fit_anp_eigen <- nodal_regression$sample(data = eigen_list, chains = n_chains, parallel_chains = n_chains, iter_warmup = n_samples, iter_sampling = n_samples)
  
  ## save output
  rm(edge_samples, adj_tensor, i) ; gc()
  save.image(paste0('anp_nodalregression/anpshort',time_window,'_nodalregression_conditionaledge_rstan.RData'))
  
  ## add progress marker
  print('model run')

  ### posterior check -- new ####
  ## summarise model
  #fit_anp_eigen$summary()
  summary(fit_anp_eigen)
  
  ## extract posterior
  params <- rstan::extract(fit_anp_eigen)
  
  ## traceplot linear effect size
  traceplot(fit_anp_eigen, pars = c('beta_age','sigma','predictor[1]','predictor[10]','predictor[20]'))
  # params %>% 
  #   select(intercept,beta_age,sigma,`predictor[1]`,`predictor[50]`,`predictor[100]`) %>% 
  #   pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'draw') %>% 
  #   mutate(chain_position = rep(rep(1:n_samples, each = 6), n_chains),
  #          chain = rep(1:n_chains, each = 6*n_samples)) %>% 
  #   #filter(chain == 4) %>% 
  #   ggplot(aes(x = chain_position, y = draw, colour = as.factor(chain)))+
  #   geom_line()+
  #   facet_wrap(. ~ parameter, scales = 'free_y')+
  #   theme_bw()+
  #   theme(legend.position = 'none') # chains 1-3 are good, chain 4 poorly mixed after about 600 draws
  
  ## check model fit
  # summary <- summary(fit_anp_eigen)
  # par(mfrow = c(3,1))
  # hist(summary$rhat, breaks = 50)
  # hist(summary$ess_bulk, breaks = 50)
  # hist(summary$ess_tail, breaks = 50)
  # par(mfrow = c(1,1))
  
  ## posterior predictive check
  plot(density(centrality_samples_std[1, ]), main="Posterior predictive check (standardised centrality):\nblack = data, blue = predicted", col=rgb(0, 0, 0, 0.25))
  for (i in 1:100) {
    j <- sample(1:length(params$beta_age), 1)
    lines(density(centrality_samples_std[j, ]), col=rgb(0, 0, 0, 0.25))
    mu <- params$beta_age[j]*eigen_list$node_age + params$intercept[j]
    sigma <- centrality_cov + diag(rep(params$sigma[j], n_eles))
    lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
  }
  
  ## add progress marker
  print('posterior predictive check complete')
  
    ## extract original values from output
  ( fit_slope <- params$beta_age * (sd(centrality_mu_ustd)/sd(nodes$age)) )
  ( fit_intcp <- mean(nodes$age) - fit_slope * mean(nodes$age) )
  plot((nodes$age*mean(fit_slope) + mean(fit_intcp)) ~ nodes$age,
       main = 'mean predictions vs age')
  
  ## plot density curves for each
  plot(density(fit_slope)) ; abline(v = 0, lty = 2)
  plot(density(fit_intcp)) ; abline(v = 0, lty = 2)
  
  ## save image
  save.image('anp_nodalregression/anpshort1_nodalregression_conditionaledge_rstan.RData')
  
  ## save output ####
  save.image(paste0('anp_nodalregression/anpshort',time_window,'_nodalregression_conditionaledge_rstan.RData'))
  dev.off()
  rm(list = ls()[! ls() %in% c('time_window','nodal_regression')])

  ## add progress marker
  print(paste0('time window ',time_window,' complete'))

}

##### predict from models -- time windows 2-36 (short) ####

##### run model with rstan -- time window 1-7 (long) ####
set.seed(12345)

for(time_window in 1:7){
  ### load data ####
  load(paste0('anp_nodalregression/anplong',time_window,'_nodalregression_conditionaledge.RData'))
  rm(counts_df, adj_mat, edges, ele_obs, obs, summary, make_edgelist, plot_network_threshold_anp) ; gc()
  
  ## set up pdf
  pdf(paste0('../outputs/anplong',time_window,'_nodalregression.pdf'))
  
  ## add progress marker
  print(paste0('start window ',time_window,' at ',Sys.time()))
  
  ### extract centralities ####
  ## build adjacency tensor
  ncol(edge_samples)
  cdf_1$node_1_id <- as.integer(as.factor(cdf_1$node_1))
  cdf_1$node_2_id <- as.integer(as.factor(cdf_1$node_2))+1
  adj_tensor <- array(0, c(n_chains*n_samples, n_eles, n_eles))
  for (dyad_id in 1:n_dyads) {
    dyad_row <- cdf_1[dyad_id, ]
    adj_tensor[, dyad_row$node_1_id, dyad_row$node_2_id] <- edge_samples[, dyad_id]
    adj_tensor[, dyad_row$node_2_id, dyad_row$node_1_id] <- edge_samples[, dyad_id]
  }
  head(adj_tensor[1,,])
  
  ## calculate centrality and store posterior samples in a matrix
  centrality_samples_invlogit <- matrix(0, n_chains*n_samples, n_eles)
  for(i in 1:(n_chains*n_samples)){
    centrality_samples_invlogit[i, ] <- sna::evcent(adj_tensor[i,,], gmode = 'graph')
  }
  head(centrality_samples_invlogit) # Unstandardised eigenvector centrality, 0-1 bounded
  
  ## convert to logit scale
  centrality_samples <- logit(centrality_samples_invlogit)
  head(centrality_samples)          # unstandardised eigenvector centrality, natural scale
  centrality_mu_ustd <- apply(centrality_samples, 2, mean)
  
  ## add mean estimate per ID to nodes data frame
  df <- data.frame(node_1_id = 1:n_eles,
                   mean_eigen_ustd = centrality_mu_ustd) %>% 
    left_join(distinct(cdf_1[,c('node_1_id','id_1')]), by = 'node_1_id') %>% 
    rename(node_rank = node_1_id,
           id = id_1)
  df$id <- ifelse(is.na(df$id) == FALSE, df$id,
                  nodes$id[! nodes$id %in% df$id])
  nodes <- nodes %>% 
    mutate(node_rank = as.integer(as.factor(node))) %>% 
    left_join(df, by = c('node_rank','id'))
  rm(df); gc()
  
  ## standardise
  centrality_samples_std <- centrality_samples
  for(i in 1:nrow(centrality_samples_std)){
    centrality_samples_std[i,] <- (centrality_samples[i,] - mean(centrality_samples[i,])) / sd(centrality_samples[i,])
  }
  head(centrality_samples_std)      # standardised eigenvector centrality, natural scale
  
  ## visualise centralities
  df_wide <- data.frame(centrality_samples_std)
  colnames(df_wide) <- 1:n_eles
  df_long <- pivot_longer(df_wide, cols = everything(),
                          names_to = "node_rank", values_to = "centrality") %>% 
    mutate(node_rank = as.integer(node_rank)) %>% 
    left_join(nodes[,c('node_rank','age')], by = 'node_rank')
  df_long %>% 
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
  
  ## add mean estimate per ID to nodes data frame
  df <- df_long %>% 
    group_by(node_rank) %>% 
    mutate(mean_eigen_std = mean(centrality)) %>% 
    dplyr::select(node_rank, mean_eigen_std) %>% 
    distinct()
  nodes <- nodes %>% 
    left_join(df, by = 'node_rank')
  rm(df); gc()
  
  ## add progress marker
  print('eigen values extracted')
  
  ### compute normal approximation ####
  ## check covariance
  plot(centrality_samples_std[, 1], centrality_samples_std[, 2],
       xlab = 'standardised eigenvectors ID1',
       ylab = 'standardised eigenvectors ID2',
       las = 1, pch = 19, col = rgb(0,0,1,0.2))
  
  ## compute normal approximation
  centrality_mu <- apply(centrality_samples_std, 2, mean)
  centrality_cov <- cov(centrality_samples_std)
  
  centrality_samples_sim <- MASS::mvrnorm(1e5, centrality_mu, centrality_cov)
  
  plot(density(centrality_samples_std[, 1]), lwd = 2, main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
  lines(density(centrality_samples_sim[, 1]), col = rgb(0,0,1,0.5), lwd = 2)
  
  ## add progress marker
  print('normal approximation computed')
  
  ### run model ####
  ## standardise age variable
  nodes$age_std <- (nodes$age - mean(nodes$age))/sd(nodes$age)
  
  ## create data list
  eigen_list <- list(num_nodes = n_eles,
                     nodes = nodes$node_rank,
                     centrality_mu = centrality_mu,
                     centrality_cov = centrality_cov,
                     node_age = nodes$age_std)
  
  ## run model
  fit_anp_eigen <- sampling(nodal_regression, data = eigen_list,
                            cores = n_chains, chains = n_chains,
                            iter = n_samples*2, warmup = n_samples)
  # fit_anp_eigen <- nodal_regression$sample(data = eigen_list, chains = n_chains, parallel_chains = n_chains, iter_warmup = n_samples, iter_sampling = n_samples)
  
  ## save output
  rm(edge_samples, adj_tensor, i) ; gc()
  save.image(paste0('anp_nodalregression/anplong',time_window,'_nodalregression_conditionaledge_rstan.RData'))
  
  ## add progress marker
  print('model run')
  
  ### posterior check -- new ####
  ## summarise model
  #fit_anp_eigen$summary()
  summary(fit_anp_eigen)
  
  ## extract posterior
  params <- rstan::extract(fit_anp_eigen)
  
  ## traceplot linear effect size
  traceplot(fit_anp_eigen, pars = c('beta_age','sigma','predictor[1]','predictor[10]','predictor[20]'))
  # params %>% 
  #   select(intercept,beta_age,sigma,`predictor[1]`,`predictor[50]`,`predictor[100]`) %>% 
  #   pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'draw') %>% 
  #   mutate(chain_position = rep(rep(1:n_samples, each = 6), n_chains),
  #          chain = rep(1:n_chains, each = 6*n_samples)) %>% 
  #   #filter(chain == 4) %>% 
  #   ggplot(aes(x = chain_position, y = draw, colour = as.factor(chain)))+
  #   geom_line()+
  #   facet_wrap(. ~ parameter, scales = 'free_y')+
  #   theme_bw()+
  #   theme(legend.position = 'none') # chains 1-3 are good, chain 4 poorly mixed after about 600 draws
  
  ## check model fit
  # summary <- summary(fit_anp_eigen)
  # par(mfrow = c(3,1))
  # hist(summary$rhat, breaks = 50)
  # hist(summary$ess_bulk, breaks = 50)
  # hist(summary$ess_tail, breaks = 50)
  # par(mfrow = c(1,1))
  
  ## posterior predictive check
  plot(density(centrality_samples_std[1, ]), main="Posterior predictive check (standardised centrality):\nblack = data, blue = predicted", col=rgb(0, 0, 0, 0.25))
  for (i in 1:100) {
    j <- sample(1:length(params$beta_age), 1)
    lines(density(centrality_samples_std[j, ]), col=rgb(0, 0, 0, 0.25))
    mu <- params$beta_age[j]*eigen_list$node_age + params$intercept[j]
    sigma <- centrality_cov + diag(rep(params$sigma[j], n_eles))
    lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
  }
  
  ## add progress marker
  print('posterior predictive check complete')
  
  ## extract original values from output
  ( fit_slope <- params$beta_age * (sd(centrality_mu_ustd)/sd(nodes$age)) )
  ( fit_intcp <- mean(nodes$age) - fit_slope * mean(nodes$age) )
  plot((nodes$age*mean(fit_slope) + mean(fit_intcp)) ~ nodes$age,
       main = 'mean predictions vs age')
  
  ## plot density curves for each
  plot(density(fit_slope)) ; abline(v = 0, lty = 2)
  plot(density(fit_intcp)) ; abline(v = 0, lty = 2)
  
  ## save image
  save.image('anp_nodalregression/anplong1_nodalregression_conditionaledge_rstan.RData')
  
  ## save output ####
  save.image(paste0('anp_nodalregression/anplong',time_window,'_nodalregression_conditionaledge_rstan.RData'))
  dev.off()
  rm(list = ls()[! ls() %in% c('time_window','nodal_regression')])
  
  ## add progress marker
  print(paste0('time window ',time_window,' complete'))
  
}

##### predict from models -- time windows 1-7 (long) ####

