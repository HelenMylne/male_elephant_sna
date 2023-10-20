##### information ####
# 1) to estimate the ages of individuals in other populations, first came up with a survival curve for the Amboseli elephants which have good birth and death data.
# 2) selected the Gompertz bathtub distribution
# 3) use distribution to convert MOTNP ages in categories to probability distributions of true age (first part of this script)
# 4) extract node centrality distribution from draws of edge weight
# 5) use MOTNP age distributions to predict network centrality (third part of this script)

##### set up ####
options(future.globals.maxSize = 10000*(1024^2))   # running model with full age distributions = 5.01GB of globals, which exceeds default maximum allowed size. Set to allow up to 10 GB to allow model to run

# library(tidyverse) ; library(cmdstanr) ; library(brms) ; library(Rcpp) ; library(ggdist) ; library(posterior) ; library(bayesplot) ; library(igraph) ; library(LaplacesDemon) ; library(bisonR) ; library(janitor)
library(tidyverse, lib.loc = '../packages/')       # library(tidyverse)
library(cmdstanr, lib.loc = '../packages/')        # library(cmdstanr)
library(brms, lib.loc = '../packages/')            # library(brms)
library(Rcpp, lib.loc = '../packages/')            # library(Rcpp)
library(ggdist, lib.loc = '../packages/')          # library(ggdist)
library(posterior, lib.loc = '../packages/')       # library(posterior)
library(bayesplot, lib.loc = '../packages/')       # library(bayesplot)
library(igraph, lib.loc = '../packages/')          # library(igraph)
library(LaplacesDemon, lib.loc = '../packages/')   # library(LaplacesDemon)
library(bisonR, lib.loc = '../packages/')          # library(bisonR)
library(janitor, lib.loc = '../packages/')         # library(janitor)

## set cmdstan path
set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0/')

## load work space from calculating edge weights
load('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')

##### prior predictive check ####
## define PDF output
pdf('../outputs/anpshort1_nodalregression_conditionalprior_bisonR.pdf')

## simulate
age <- 1:60
beta_mu <- 0
beta_sigma <- 0.1
mean_age <- mean(nodes$age)
plot(NULL, xlim = c(10,60), ylim = c(0,1), las = 1,
     xlab = 'age', ylab = 'eigenvector centrality')
for(i in 1:100){
  beta <- rnorm(1, beta_mu, beta_sigma)
  y <- age*beta
  lines(x = age, y = LaplacesDemon::invlogit(y), col = rgb(0,0,1,0.5)) # vast majority come out somewhere sensible, and those that don't would if they started at a different value for age 10 so fine in combo with posterior intercept
}

##### run model with rstan -- time window 1 ####
library(rstan) ; library(igraph) ; library(tidyverse) ; library(LaplacesDemon)
set.seed(12345)

# load data and remove additional data
load('anp_nodalregression/anpshort1_nodalregression_conditionaledge.RData')
rm(counts_df, adj_mat, edges, ele_obs, obs, summary, x, make_edgelist, plot_network_threshold_anp) ; gc()

### extract centralities ####
## build adjacency tensor
ncol(edge_samples)
cdf_1$node_1_id <- as.integer(as.factor(cdf_1$node_1))
cdf_1$node_2_id <- as.integer(as.factor(cdf_1$node_2))+1
adj_tensor <- array(0, c(n_chains*n_samples, n_eles, n_eles))
for (dyad_id in 1:n_dyads) {
  dyad_row <- cdf_1[dyad_id, ]
  adj_tensor[, dyad_row$node_1_id, dyad_row$node_2_id] <- edge_samples[, dyad_id]
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
nodes$node_rank <- as.integer(as.factor(nodes$node))
df_wide <- data.frame(centrality_samples_std)
colnames(df_wide) <- 1:n_eles
df_long <- pivot_longer(df_wide, cols = everything(),
                        names_to = "node_rank", values_to = "centrality") %>% 
  mutate(node_rank = as.integer(node_rank)) %>% 
  left_join(nodes[,c('node_rank','age')], by = 'node_rank')
df_long %>% 
  mutate(nodes_reordered = fct_reorder(.f = as.factor(node_rank), .x = age, .desc = T)) %>% 
  ggplot(aes(x = centrality, fill = age)) +
  geom_density(size = 0.4) +
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

## compute normal approximation
centrality_mu <- apply(centrality_samples_std, 2, mean)
centrality_cov <- cov(centrality_samples_std)

centrality_samples_sim <- MASS::mvrnorm(1e5, centrality_mu, centrality_cov)

plot(density(centrality_samples_std[, 1]), lwd = 2, main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
lines(density(centrality_samples_sim[, 1]), col = rgb(0,0,1,0.5), lwd = 2)

### run model ####
## create data list
eigen_list <- list(num_nodes = n_eles,
                   nodes = nodes$node_rank,
                   centrality_mu = centrality_mu,
                   centrality_cov = centrality_cov,
                   node_age = nodes$age)

# load model
nodal_regression <- stan_model('models/eigen_regression.stan')

## run model
fit_anp1_eigen <- sampling(nodal_regression,
                           data = eigen_list,
                           cores = n_chains,
                           chains = n_chains)

## save output
rm(g, edge_samples, adj_tensor, i) ; gc()
save.image('anp_nodalregression/anpshort1_nodalregression_conditionaledge_rstan.RData')

### posterior check ####
# load('anp_nodalregression/anpshort1_nodalregression_conditionaledge_rstan.RData')
## traceplot linear effect size
traceplot(fit_anp1_eigen, pars = c('beta_age','sigma','predictor[1]','predictor[2]','predictor[3]','predictor[4]','predictor[5]','predictor[6]','predictor[7]','predictor[8]','predictor[9]','predictor[10]'))

## posterior predictive check
params <- rstan::extract(fit_anp1_eigen)
plot(density(centrality_samples_std[1, ]), main="Posterior predictive density of responses (standardised centrality)", col=rgb(0, 0, 0, 0.25), ylim=c(0, 1))
for (i in 1:100) {
  j <- sample(1:(n_chains*n_samples), 1)
  lines(density(centrality_samples_std[j, ]), col=rgb(0, 0, 0, 0.25))
  mu <- params$beta_age[j]*eigen_list$node_age
  sigma <- centrality_cov + diag(rep(params$sigma[j], n_eles))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}

## interpret model
plot(density(params$beta_age), # don't need to do the contrast because continuous -- is there an equivalent using the do operator? by standardising have we already dealt with the logit transformation?
     main = "Posterior difference between node types")
abline(v = 0, lty = 2)

### predict from model ####
(summary <- as.data.frame(round(summary(fit_anp1_eigen)$summary[1:2, c(1, 4, 8)], 3)))
summary$parameter <- rownames(summary)

## calculate mean predictions for model
mod_mu <- data.frame(age = min(nodes$age):max(nodes$age)) %>% 
  mutate(lwr = age*summary$`2.5%`[1],
         mid = age*summary$mean[1],
         upr = age*summary$`97.5%`[1])

## plot mean predictions
plot(mid ~ age, data = mod_mu, type = 'l', las = 1, col = 'blue', lwd = 2,
     ylim = c(min(df_long$centrality), max(df_long$centrality)),
     xlab = 'age (years)', ylab = 'eigenvector centrality (standardised)')
lines(lwr ~ age, data = mod_mu, lty = 2)
lines(upr ~ age, data = mod_mu, lty = 2)

## simulate full predictions for model
sim <- matrix(NA, nrow = n_chains*n_samples, ncol = nrow(mod_mu),
              dimnames = list(1:(n_chains*n_samples), mod_mu$age))
for(i in 1:nrow(sim)){
  for(j in 1:ncol(sim)){
    #sim[i,j] <- params$beta_age[i]*x$age[j]
    sim[i,j] <- MASS::mvrnorm(n = 1, mu = params$beta_age[i]*mod_mu$age[j],
                              Sigma = params$sigma[i])
  }
}

## plot simulations
sim_df <- as.data.frame(sim) %>% 
  pivot_longer(cols = everything(), names_to = 'age', values_to = 'eigen_sim')
points(sim_df$eigen_sim ~ sim_df$age, col = rgb(0,0,0,0.01), pch = 19, cex = 0.5)

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
  geom_ribbon(data = sim_summary, aes(x = age, ymin = lwr, ymax = upr),       # shade simulations
              colour = 'transparent', fill = rgb(0,0,0,0.1))+
  geom_ribbon(data = mod_mu, aes(x = age, ymin = lwr, ymax = upr),            # shade mean distribution
              colour = 'transparent', fill = rgb(33/255, 145/255, 140/255, 0.5))+
  geom_point(data = df_long, aes(x = age, y = centrality),                    # all eigenvector draws
             colour = rgb(253/255, 231/255, 37/255, 0.01))+
  geom_point(data = df_long, aes(x = age, y = mean_eigen, size = sightings),  # mean eigenvector
             colour = rgb(68/255, 1/255, 84/255))+
  geom_line(data = mod_mu, aes(x = age, y = mid),                             # mean line
            colour = rgb(33/255, 145/255, 140/255), linewidth = 1)+
  scale_x_continuous('age (years)')+
  scale_y_continuous('eigenvector centrality (standardised)')+
  theme_classic()+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22))
ggsave(filename = '../outputs/anpshort1_nodalregression.png', device = 'png',
       plot = last_plot())

## save output
save.image('anp_nodalregression/anpshort1_nodalregression_conditionaledge_rstan.RData')
dev.off()
rm(list = ls()[!ls() %in% 'nodal_regression'])

##### run model with rstan -- time window 2-36 ####
set.seed(12345)

for(time_window in 2:36){
  # load data and remove additional data
  load(paste0('anp_nodalregression/anpshort',time_window,'_nodalregression_conditionaledge.RData'))
  rm(counts_df, adj_mat, edges, ele_obs, obs, summary, x, make_edgelist, plot_network_threshold_anp) ; gc()
  
  ### extract centralities ####
  ## build adjacency tensor
  ncol(edge_samples)
  cdf$node_1_id <- as.integer(as.factor(cdf$node_1))
  cdf$node_2_id <- as.integer(as.factor(cdf$node_2))+1
  adj_tensor <- array(0, c(n_chains*n_samples, n_eles, n_eles))
  for (dyad_id in 1:n_dyads) {
    dyad_row <- cdf[dyad_id, ]
    adj_tensor[, dyad_row$node_1_id, dyad_row$node_2_id] <- edge_samples[, dyad_id]
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
  nodes$node_rank <- as.integer(as.factor(nodes$node))
  df_wide <- data.frame(centrality_samples_std)
  colnames(df_wide) <- 1:n_eles
  df_long <- pivot_longer(df_wide, cols = everything(),
                          names_to = "node_rank", values_to = "centrality") %>% 
    mutate(node_rank = as.integer(node_rank)) %>% 
    left_join(nodes[,c('node_rank','age')], by = 'node_rank')
  df_long %>% 
    mutate(nodes_reordered = fct_reorder(.f = as.factor(node_rank), .x = age, .desc = T)) %>% 
    ggplot(aes(x = centrality, fill = age)) +
    geom_density(size = 0.4) +
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
  
  ## compute normal approximation
  centrality_mu <- apply(centrality_samples_std, 2, mean)
  centrality_cov <- cov(centrality_samples_std)
  
  centrality_samples_sim <- MASS::mvrnorm(1e5, centrality_mu, centrality_cov)
  
  plot(density(centrality_samples_std[, 1]), lwd = 2, main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
  lines(density(centrality_samples_sim[, 1]), col = rgb(0,0,1,0.5), lwd = 2)
  
  ### run model ####
  ## create data list
  eigen_list <- list(num_nodes = n_eles,
                     nodes = nodes$node_rank,
                     centrality_mu = centrality_mu,
                     centrality_cov = centrality_cov,
                     node_age = nodes$age)
  
  # load model
  nodal_regression <- stan_model('models/eigen_regression.stan')
  
  ## run model
  fit_anp_eigen <- sampling(nodal_regression,
                             data = eigen_list,
                             cores = n_chains,
                             chains = n_chains)
  
  ## save output
  rm(g, edge_samples, adj_tensor, i) ; gc()
  save.image(paste0('anp_nodalregression/anpshort',time_window,'_nodalregression_conditionaledge_rstan.RData'))
  
  ### posterior check ####
  # load(paste0('anp_nodalregression/anpshort',time_window,'_nodalregression_conditionaledge_rstan.RData'))
  ## traceplot linear effect size
  traceplot(fit_anp_eigen, pars = c('beta_age','sigma','predictor[1]','predictor[2]','predictor[3]','predictor[4]','predictor[5]','predictor[6]','predictor[7]','predictor[8]','predictor[9]','predictor[10]'))
  
  ## posterior predictive check
  params <- rstan::extract(fit_anp_eigen)
  plot(density(centrality_samples_std[1, ]), main="Posterior predictive density of responses (standardised centrality)", col=rgb(0, 0, 0, 0.25), ylim=c(0, 1))
  for (i in 1:100) {
    j <- sample(1:(n_chains*n_samples), 1)
    lines(density(centrality_samples_std[j, ]), col=rgb(0, 0, 0, 0.25))
    mu <- params$beta_age[j]*eigen_list$node_age
    sigma <- centrality_cov + diag(rep(params$sigma[j], n_eles))
    lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
  }
  
  ## interpret model
  plot(density(params$beta_age), # don't need to do the contrast because continuous -- is there an equivalent using the do operator? by standardising have we already dealt with the logit transformation?
       main = "Posterior difference between node types")
  abline(v = 0, lty = 2)
  
  ### predict from model ####
  (summary <- as.data.frame(round(summary(fit_anp_eigen)$summary[1:2, c(1, 4, 8)], 3)))
  summary$parameter <- rownames(summary)
  
  ## calculate mean predictions for model
  mod_mu <- data.frame(age = min(nodes$age):max(nodes$age)) %>% 
    mutate(lwr = age*summary$`2.5%`[1],
           mid = age*summary$mean[1],
           upr = age*summary$`97.5%`[1])
  
  ## plot mean predictions
  plot(mid ~ age, data = mod_mu, type = 'l', las = 1, col = 'blue', lwd = 2,
       ylim = c(min(df_long$centrality), max(df_long$centrality)),
       xlab = 'age (years)', ylab = 'eigenvector centrality (standardised)')
  lines(lwr ~ age, data = mod_mu, lty = 2)
  lines(upr ~ age, data = mod_mu, lty = 2)
  
  ## simulate full predictions for model
  sim <- matrix(NA, nrow = n_chains*n_samples, ncol = nrow(mod_mu),
                dimnames = list(1:(n_chains*n_samples), mod_mu$age))
  for(i in 1:nrow(sim)){
    for(j in 1:ncol(sim)){
      #sim[i,j] <- params$beta_age[i]*x$age[j]
      sim[i,j] <- MASS::mvrnorm(n = 1, mu = params$beta_age[i]*mod_mu$age[j],
                                Sigma = params$sigma[i])
    }
  }
  
  ## plot simulations
  sim_df <- as.data.frame(sim) %>% 
    pivot_longer(cols = everything(), names_to = 'age', values_to = 'eigen_sim')
  points(sim_df$eigen_sim ~ sim_df$age, col = rgb(0,0,0,0.01), pch = 19, cex = 0.5)
  
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
    geom_ribbon(data = sim_summary, aes(x = age, ymin = lwr, ymax = upr),       # shade simulations
                colour = 'transparent', fill = rgb(0,0,0,0.1))+
    geom_ribbon(data = mod_mu, aes(x = age, ymin = lwr, ymax = upr),            # shade mean distribution
                colour = 'transparent', fill = rgb(33/255, 145/255, 140/255, 0.5))+
    geom_point(data = df_long, aes(x = age, y = centrality),                    # all eigenvector draws
               colour = rgb(253/255, 231/255, 37/255, 0.01))+
    geom_point(data = df_long, aes(x = age, y = mean_eigen, size = sightings),  # mean eigenvector
               colour = rgb(68/255, 1/255, 84/255))+
    geom_line(data = mod_mu, aes(x = age, y = mid),                             # mean line
              colour = rgb(33/255, 145/255, 140/255), linewidth = 1)+
    scale_x_continuous('age (years)')+
    scale_y_continuous('eigenvector centrality (standardised)')+
    theme_classic()+
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 22),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 22))
  ggsave(filename = '../outputs/anpshort',time_window,'_nodalregression.png', device = 'png',
         plot = last_plot())
  
  ## save output ####
  save.image(paste0('anp_nodalregression/anpshort',time_window,'_nodalregression_conditionaledge_rstan.RData'))
  dev.off()
  rm(list = ls()[! ls() %in% c('time_window','nodal_regression')])
  
}
  
##### run model with rstan -- time window 1-7 (long) ####
set.seed(12345)

for(time_window in 1:7){
  # load data and remove additional data
  load(paste0('anp_nodalregression/anplong',time_window,'_nodalregression_conditionaledge.RData'))
  rm(counts_df, adj_mat, edges, ele_obs, obs, summary, x, make_edgelist, plot_network_threshold_anp) ; gc()
  
  ### extract centralities ####
  ## build adjacency tensor
  ncol(edge_samples)
  cdf$node_1_id <- as.integer(as.factor(cdf$node_1))
  cdf$node_2_id <- as.integer(as.factor(cdf$node_2))+1
  adj_tensor <- array(0, c(n_chains*n_samples, n_eles, n_eles))
  for (dyad_id in 1:n_dyads) {
    dyad_row <- cdf[dyad_id, ]
    adj_tensor[, dyad_row$node_1_id, dyad_row$node_2_id] <- edge_samples[, dyad_id]
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
  nodes$node_rank <- as.integer(as.factor(nodes$node))
  df_wide <- data.frame(centrality_samples_std)
  colnames(df_wide) <- 1:n_eles
  df_long <- pivot_longer(df_wide, cols = everything(),
                          names_to = "node_rank", values_to = "centrality") %>% 
    mutate(node_rank = as.integer(node_rank)) %>% 
    left_join(nodes[,c('node_rank','age')], by = 'node_rank')
  df_long %>% 
    mutate(nodes_reordered = fct_reorder(.f = as.factor(node_rank), .x = age, .desc = T)) %>% 
    ggplot(aes(x = centrality, fill = age)) +
    geom_density(size = 0.4) +
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
  
  ## compute normal approximation
  centrality_mu <- apply(centrality_samples_std, 2, mean)
  centrality_cov <- cov(centrality_samples_std)
  
  centrality_samples_sim <- MASS::mvrnorm(1e5, centrality_mu, centrality_cov)
  
  plot(density(centrality_samples_std[, 1]), lwd = 2, main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
  lines(density(centrality_samples_sim[, 1]), col = rgb(0,0,1,0.5), lwd = 2)
  
  ### run model ####
  ## create data list
  eigen_list <- list(num_nodes = n_eles,
                     nodes = nodes$node_rank,
                     centrality_mu = centrality_mu,
                     centrality_cov = centrality_cov,
                     node_age = nodes$age)
  
  # load model
  nodal_regression <- stan_model('models/eigen_regression.stan')
  
  ## run model
  fit_anp_eigen <- sampling(nodal_regression,
                            data = eigen_list,
                            cores = n_chains,
                            chains = n_chains)
  
  ## save output
  rm(g, edge_samples, adj_tensor, i) ; gc()
  save.image(paste0('anp_nodalregression/anplong',time_window,'_nodalregression_conditionaledge_rstan.RData'))
  
  ### posterior check ####
  # load(paste0('anp_nodalregression/anplong',time_window,'_nodalregression_conditionaledge_rstan.RData'))
  ## traceplot linear effect size
  traceplot(fit_anp_eigen, pars = c('beta_age','sigma','predictor[1]','predictor[2]','predictor[3]','predictor[4]','predictor[5]','predictor[6]','predictor[7]','predictor[8]','predictor[9]','predictor[10]'))
  
  ## posterior predictive check
  params <- rstan::extract(fit_anp_eigen)
  plot(density(centrality_samples_std[1, ]), main="Posterior predictive density of responses (standardised centrality)", col=rgb(0, 0, 0, 0.25), ylim=c(0, 1))
  for (i in 1:100) {
    j <- sample(1:(n_chains*n_samples), 1)
    lines(density(centrality_samples_std[j, ]), col=rgb(0, 0, 0, 0.25))
    mu <- params$beta_age[j]*eigen_list$node_age
    sigma <- centrality_cov + diag(rep(params$sigma[j], n_eles))
    lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
  }
  
  ## interpret model
  plot(density(params$beta_age), # don't need to do the contrast because continuous -- is there an equivalent using the do operator? by standardising have we already dealt with the logit transformation?
       main = "Posterior difference between node types")
  abline(v = 0, lty = 2)
  
  ### predict from model ####
  (summary <- as.data.frame(round(summary(fit_anp_eigen)$summary[1:2, c(1, 4, 8)], 3)))
  summary$parameter <- rownames(summary)
  
  ## calculate mean predictions for model
  mod_mu <- data.frame(age = min(nodes$age):max(nodes$age)) %>% 
    mutate(lwr = age*summary$`2.5%`[1],
           mid = age*summary$mean[1],
           upr = age*summary$`97.5%`[1])
  
  ## plot mean predictions
  plot(mid ~ age, data = mod_mu, type = 'l', las = 1, col = 'blue', lwd = 2,
       ylim = c(min(df_long$centrality), max(df_long$centrality)),
       xlab = 'age (years)', ylab = 'eigenvector centrality (standardised)')
  lines(lwr ~ age, data = mod_mu, lty = 2)
  lines(upr ~ age, data = mod_mu, lty = 2)
  
  ## simulate full predictions for model
  sim <- matrix(NA, nrow = n_chains*n_samples, ncol = nrow(mod_mu),
                dimnames = list(1:(n_chains*n_samples), mod_mu$age))
  for(i in 1:nrow(sim)){
    for(j in 1:ncol(sim)){
      #sim[i,j] <- params$beta_age[i]*x$age[j]
      sim[i,j] <- MASS::mvrnorm(n = 1, mu = params$beta_age[i]*mod_mu$age[j],
                                Sigma = params$sigma[i])
    }
  }
  
  ## plot simulations
  sim_df <- as.data.frame(sim) %>% 
    pivot_longer(cols = everything(), names_to = 'age', values_to = 'eigen_sim')
  points(sim_df$eigen_sim ~ sim_df$age, col = rgb(0,0,0,0.01), pch = 19, cex = 0.5)
  
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
    geom_ribbon(data = sim_summary, aes(x = age, ymin = lwr, ymax = upr),       # shade simulations
                colour = 'transparent', fill = rgb(0,0,0,0.1))+
    geom_ribbon(data = mod_mu, aes(x = age, ymin = lwr, ymax = upr),            # shade mean distribution
                colour = 'transparent', fill = rgb(33/255, 145/255, 140/255, 0.5))+
    geom_point(data = df_long, aes(x = age, y = centrality),                    # all eigenvector draws
               colour = rgb(253/255, 231/255, 37/255, 0.01))+
    geom_point(data = df_long, aes(x = age, y = mean_eigen, size = sightings),  # mean eigenvector
               colour = rgb(68/255, 1/255, 84/255))+
    geom_line(data = mod_mu, aes(x = age, y = mid),                             # mean line
              colour = rgb(33/255, 145/255, 140/255), linewidth = 1)+
    scale_x_continuous('age (years)')+
    scale_y_continuous('eigenvector centrality (standardised)')+
    theme_classic()+
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 22),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 22))
  ggsave(filename = '../outputs/anplong',time_window,'_nodalregression.png', device = 'png',
         plot = last_plot())
  
  ## save output ####
  save.image(paste0('anp_nodalregression/anplong',time_window,'_nodalregression_conditionaledge_rstan.RData'))
  dev.off()
  rm(list = ls()[! ls() %in% c('time_window','nodal_regression')])
  
}
