#### information ####
# script for prior predicitve checks for nodal regressions -- supplementary information

#### set up ####
library(LaplacesDemon)
library(tidyverse)
theme_set(theme_bw())
library(patchwork)

set.seed(12345)

n <- 100
window <- rnorm(n, 0, 1)
node <- rnorm(n, 0, 1)

#### motnp ####
load('step4_nodalregression/motnp_nodalregression.RData')
rm(list = ls()[!ls() %in% c('nodes','centrality_mu','n','beta_age','window','node')]) ; gc()

intercept  <- rnorm(n, logit(0.05), 1)
n_age_cat <- length(unique(nodes$age_cat_fct))
age_dirichlet <- rdirichlet(n, c(1,1,1,1,1))
beta_age <- rnorm(n, 0, 2)

## plot
plot(NULL, las = 1, xlab = 'age category', ylab = 'logit eigenvector',
     ylim = c(min(centrality_mu)-5, max(centrality_mu)+5), xlim = c(1,n_age_cat))
abline(h = min(centrality_mu), lty = 2) ; abline(h = max(centrality_mu), lty = 2)
x <- 1:n_age_cat
for(i in 1:n){
  y <- rep(NA, length(x))
  for(j in 1:length(x)){
    y[j] <- intercept[i] + beta_age[i]*sum(age_dirichlet[i,][1:x[j]])
  }
  lines(x = x, y = y, col = rgb(0,0,1,0.4))
}

## long format simulation
sim_mot <- data.frame(row_num = rep(1:n, each = n_age_cat),
                      x = rep(1:n_age_cat, n),
                      intercept = rep(intercept, each = n_age_cat),
                      beta_age = rep(beta_age, each = n_age_cat),
                      y = NA)
for(i in 1:nrow(sim_mot)){
  sim_mot$y[i] <- sim_mot$intercept[i] + 
    sim_mot$beta_age[i]*sum(age_dirichlet[sim_mot$row_num[i],1:sim_mot$x[i]])
}
(motnp <- ggplot(sim_mot)+
  geom_line(aes(x = x, y = y, group = row_num),
            colour = '#FDE725FF')+
  geom_hline(yintercept = min(centrality_mu),
             linetype = 2, linewidth = 1)+
  geom_hline(yintercept = max(centrality_mu),
             linetype = 2, linewidth = 1)+
  labs(x = 'age category',
       y = 'logit eigenvector centrality'))
sim_mot <- sim_mot %>% 
  mutate(population = 'Mosi-Oa-Tunya',
         min_raw = min(centrality_mu),
         max_raw = max(centrality_mu))

rm(list = ls()[!ls() %in% c('n','window','node','motnp','sim_mot')]) ; gc()

#### anp short ####
load('anp_nodalregression/anp_short_nodal_modelprep.RData')
beta_age <- rnorm(n, 0, 0.8)
n_data <- nrow(nodes_all)
n_nodes <- length(unique(nodes_all$node_random))
n_windows <- length(unique(nodes_all$window))

model_data_list <- list(
  # global data size
  num_data = n_data,
  num_nodes = n_nodes,
  num_windows = n_windows,
  # per time window data size
  num_nodes_window1 = length(unique(nodes_all$node_random[nodes_all$window == 1])),
  num_nodes_window2 = length(unique(nodes_all$node_random[nodes_all$window == 2])),
  num_nodes_window3 = length(unique(nodes_all$node_random[nodes_all$window == 3])),
  num_nodes_window4 = length(unique(nodes_all$node_random[nodes_all$window == 4])),
  num_nodes_window5 = length(unique(nodes_all$node_random[nodes_all$window == 5])),
  num_nodes_window6 = length(unique(nodes_all$node_random[nodes_all$window == 6])),
  num_nodes_window7 = length(unique(nodes_all$node_random[nodes_all$window == 7])),
  num_nodes_window8 = length(unique(nodes_all$node_random[nodes_all$window == 8])),
  num_nodes_window9 = length(unique(nodes_all$node_random[nodes_all$window == 9])),
  num_nodes_window10 = length(unique(nodes_all$node_random[nodes_all$window == 10])),
  num_nodes_window11 = length(unique(nodes_all$node_random[nodes_all$window == 11])),
  num_nodes_window12 = length(unique(nodes_all$node_random[nodes_all$window == 12])),
  num_nodes_window13 = length(unique(nodes_all$node_random[nodes_all$window == 13])),
  num_nodes_window14 = length(unique(nodes_all$node_random[nodes_all$window == 14])),
  num_nodes_window15 = length(unique(nodes_all$node_random[nodes_all$window == 15])),
  num_nodes_window16 = length(unique(nodes_all$node_random[nodes_all$window == 16])),
  num_nodes_window17 = length(unique(nodes_all$node_random[nodes_all$window == 17])),
  num_nodes_window18 = length(unique(nodes_all$node_random[nodes_all$window == 18])),
  num_nodes_window19 = length(unique(nodes_all$node_random[nodes_all$window == 19])),
  num_nodes_window20 = length(unique(nodes_all$node_random[nodes_all$window == 20])),
  num_nodes_window21 = length(unique(nodes_all$node_random[nodes_all$window == 21])),
  num_nodes_window22 = length(unique(nodes_all$node_random[nodes_all$window == 22])),
  num_nodes_window23 = length(unique(nodes_all$node_random[nodes_all$window == 23])),
  num_nodes_window24 = length(unique(nodes_all$node_random[nodes_all$window == 24])),
  num_nodes_window25 = length(unique(nodes_all$node_random[nodes_all$window == 25])),
  num_nodes_window26 = length(unique(nodes_all$node_random[nodes_all$window == 26])),
  num_nodes_window27 = length(unique(nodes_all$node_random[nodes_all$window == 27])),
  num_nodes_window28 = length(unique(nodes_all$node_random[nodes_all$window == 28])),
  num_nodes_window29 = length(unique(nodes_all$node_random[nodes_all$window == 29])),
  num_nodes_window30 = length(unique(nodes_all$node_random[nodes_all$window == 30])),
  num_nodes_window31 = length(unique(nodes_all$node_random[nodes_all$window == 31])),
  num_nodes_window32 = length(unique(nodes_all$node_random[nodes_all$window == 32])),
  num_nodes_window33 = length(unique(nodes_all$node_random[nodes_all$window == 33])),
  num_nodes_window34 = length(unique(nodes_all$node_random[nodes_all$window == 34])),
  num_nodes_window35 = length(unique(nodes_all$node_random[nodes_all$window == 35])),
  num_nodes_window36 = length(unique(nodes_all$node_random[nodes_all$window == 36])),
  # number of nodes in all preceding time windows for node age indexing
  num_nodes_prev_windows = c(0,
                             length(which(nodes_all$window < 2)),
                             length(which(nodes_all$window < 3)),
                             length(which(nodes_all$window < 4)),
                             length(which(nodes_all$window < 5)),
                             length(which(nodes_all$window < 6)),
                             length(which(nodes_all$window < 7)),
                             length(which(nodes_all$window < 8)),
                             length(which(nodes_all$window < 9)),
                             length(which(nodes_all$window < 10)),
                             length(which(nodes_all$window < 11)),
                             length(which(nodes_all$window < 12)),
                             length(which(nodes_all$window < 13)),
                             length(which(nodes_all$window < 14)),
                             length(which(nodes_all$window < 15)),
                             length(which(nodes_all$window < 16)),
                             length(which(nodes_all$window < 17)),
                             length(which(nodes_all$window < 18)),
                             length(which(nodes_all$window < 19)),
                             length(which(nodes_all$window < 20)),
                             length(which(nodes_all$window < 21)),
                             length(which(nodes_all$window < 22)),
                             length(which(nodes_all$window < 23)),
                             length(which(nodes_all$window < 24)),
                             length(which(nodes_all$window < 25)),
                             length(which(nodes_all$window < 26)),
                             length(which(nodes_all$window < 27)),
                             length(which(nodes_all$window < 28)),
                             length(which(nodes_all$window < 29)),
                             length(which(nodes_all$window < 30)),
                             length(which(nodes_all$window < 31)),
                             length(which(nodes_all$window < 32)),
                             length(which(nodes_all$window < 33)),
                             length(which(nodes_all$window < 34)),
                             length(which(nodes_all$window < 35)),
                             length(which(nodes_all$window < 36))),
  # centrality means per time window
  centrality_mu_1 = nodes_all$mean_eigen[nodes_all$window == 1],
  centrality_mu_2 = nodes_all$mean_eigen[nodes_all$window == 2],
  centrality_mu_3 = nodes_all$mean_eigen[nodes_all$window == 3],
  centrality_mu_4 = nodes_all$mean_eigen[nodes_all$window == 4],
  centrality_mu_5 = nodes_all$mean_eigen[nodes_all$window == 5],
  centrality_mu_6 = nodes_all$mean_eigen[nodes_all$window == 6],
  centrality_mu_7 = nodes_all$mean_eigen[nodes_all$window == 7],
  centrality_mu_8 = nodes_all$mean_eigen[nodes_all$window == 8],
  centrality_mu_9 = nodes_all$mean_eigen[nodes_all$window == 9],
  centrality_mu_10 = nodes_all$mean_eigen[nodes_all$window == 10],
  centrality_mu_11 = nodes_all$mean_eigen[nodes_all$window == 11],
  centrality_mu_12 = nodes_all$mean_eigen[nodes_all$window == 12],
  centrality_mu_13 = nodes_all$mean_eigen[nodes_all$window == 13],
  centrality_mu_14 = nodes_all$mean_eigen[nodes_all$window == 14],
  centrality_mu_15 = nodes_all$mean_eigen[nodes_all$window == 15],
  centrality_mu_16 = nodes_all$mean_eigen[nodes_all$window == 16],
  centrality_mu_17 = nodes_all$mean_eigen[nodes_all$window == 17],
  centrality_mu_18 = nodes_all$mean_eigen[nodes_all$window == 18],
  centrality_mu_19 = nodes_all$mean_eigen[nodes_all$window == 19],
  centrality_mu_20 = nodes_all$mean_eigen[nodes_all$window == 20],
  centrality_mu_21 = nodes_all$mean_eigen[nodes_all$window == 21],
  centrality_mu_22 = nodes_all$mean_eigen[nodes_all$window == 22],
  centrality_mu_23 = nodes_all$mean_eigen[nodes_all$window == 23],
  centrality_mu_24 = nodes_all$mean_eigen[nodes_all$window == 24],
  centrality_mu_25 = nodes_all$mean_eigen[nodes_all$window == 25],
  centrality_mu_26 = nodes_all$mean_eigen[nodes_all$window == 26],
  centrality_mu_27 = nodes_all$mean_eigen[nodes_all$window == 27],
  centrality_mu_28 = nodes_all$mean_eigen[nodes_all$window == 28],
  centrality_mu_29 = nodes_all$mean_eigen[nodes_all$window == 29],
  centrality_mu_30 = nodes_all$mean_eigen[nodes_all$window == 30],
  centrality_mu_31 = nodes_all$mean_eigen[nodes_all$window == 31],
  centrality_mu_32 = nodes_all$mean_eigen[nodes_all$window == 32],
  centrality_mu_33 = nodes_all$mean_eigen[nodes_all$window == 33],
  centrality_mu_34 = nodes_all$mean_eigen[nodes_all$window == 34],
  centrality_mu_35 = nodes_all$mean_eigen[nodes_all$window == 35],
  centrality_mu_36 = nodes_all$mean_eigen[nodes_all$window == 36],
  # covariance matrix per time window
  centrality_cov_1 = covs_all[[1]],
  centrality_cov_2 = covs_all[[2]],
  centrality_cov_3 = covs_all[[3]],
  centrality_cov_4 = covs_all[[4]],
  centrality_cov_5 = covs_all[[5]],
  centrality_cov_6 = covs_all[[6]],
  centrality_cov_7 = covs_all[[7]],
  centrality_cov_8 = covs_all[[8]],
  centrality_cov_9 = covs_all[[9]],
  centrality_cov_10 = covs_all[[10]],
  centrality_cov_11 = covs_all[[11]],
  centrality_cov_12 = covs_all[[12]],
  centrality_cov_13 = covs_all[[13]],
  centrality_cov_14 = covs_all[[14]],
  centrality_cov_15 = covs_all[[15]],
  centrality_cov_16 = covs_all[[16]],
  centrality_cov_17 = covs_all[[17]],
  centrality_cov_18 = covs_all[[18]],
  centrality_cov_19 = covs_all[[19]],
  centrality_cov_20 = covs_all[[20]],
  centrality_cov_21 = covs_all[[21]],
  centrality_cov_22 = covs_all[[22]],
  centrality_cov_23 = covs_all[[23]],
  centrality_cov_24 = covs_all[[24]],
  centrality_cov_25 = covs_all[[25]],
  centrality_cov_26 = covs_all[[26]],
  centrality_cov_27 = covs_all[[27]],
  centrality_cov_28 = covs_all[[28]],
  centrality_cov_29 = covs_all[[29]],
  centrality_cov_30 = covs_all[[30]],
  centrality_cov_31 = covs_all[[31]],
  centrality_cov_32 = covs_all[[32]],
  centrality_cov_33 = covs_all[[33]],
  centrality_cov_34 = covs_all[[34]],
  centrality_cov_35 = covs_all[[35]],
  centrality_cov_36 = covs_all[[36]],
  # node IDs for all time windows
  node_id = nodes_all$node_random,
  # exposure variable
  node_age = nodes_all$age_std)

intercept  <- rnorm(n, LaplacesDemon::logit(0.05), 2) # taking the intercept from the logit of results from Chiyo 2011 (doesn't state mean/median centrality so estimated from graph based on where correlation line would cross x = 0)
window <- rnorm(n, 0, 1)
node <- rnorm(n, 0, 1)
min_raw <- min(unlist(model_data_list[grep('centrality_mu', names(model_data_list))]))
max_raw <- max(unlist(model_data_list[grep('centrality_mu', names(model_data_list))]))
plot(NULL, las = 1, xlab = 'age (standardised)', ylab = 'eigenvector',
     #ylim = c(min_raw-2, max_raw+2),
     ylim = c(-15, 5),
     xlim = c(min(model_data_list$node_age),
              max(model_data_list$node_age)))
abline(h = min_raw, lty = 2) ; abline(h = max_raw, lty = 2)
for(i in 1:n){
  lines(x = seq(min(nodes_all$age_std), max(nodes_all$age_std), length.out = 2),
        y = intercept[i] + beta_age[i]*c(min(nodes_all$age_std),
                                         max(nodes_all$age_std)) + window[i] + node[i],
        col = rgb(0,0,1,0.4))
}

## long format simulation
sim_ash <- data.frame(row_num = rep(1:n, each = 2),
                      x = rep(c(min(nodes_all$age_std), max(nodes_all$age_std)), n),
                      intercept = rep(intercept, each = 2),
                      beta_age = rep(beta_age, each = 2)) %>% 
  mutate(y = intercept + beta_age * x)
(short <- ggplot(sim_ash)+
    geom_line(aes(x = x, y = y, group = row_num),
              colour = '#1F968BFF')+
    geom_hline(yintercept = min_raw,
               linetype = 2, linewidth = 1)+
    geom_hline(yintercept = max_raw,
               linetype = 2, linewidth = 1)+
    labs(x = 'age (standardised)',
         y = 'logit eigenvector centrality'))
sim_ash <- sim_ash %>% 
  mutate(population = 'Amboseli (short)',
         min_raw = min_raw,
         max_raw = max_raw)
rm(list = ls()[!ls() %in% c('n','beta_age','window','node','motnp','short','sim_mot','sim_ash')]) ; gc()

#### anp long ####
load('anp_nodalregression/anp_long_nodal_modelprep.RData')
n_data <- nrow(nodes_all)
n_nodes <- length(unique(nodes_all$node_random))
n_windows <- length(unique(nodes_all$window))

model_data_list <- list(
  # global data size
  num_data = n_data,
  num_nodes = n_nodes,
  num_windows = n_windows,
  # per time window data size
  num_nodes_window1 = length(unique(nodes_all$node_random[nodes_all$window == 1])),
  num_nodes_window2 = length(unique(nodes_all$node_random[nodes_all$window == 2])),
  num_nodes_window3 = length(unique(nodes_all$node_random[nodes_all$window == 3])),
  num_nodes_window4 = length(unique(nodes_all$node_random[nodes_all$window == 4])),
  num_nodes_window5 = length(unique(nodes_all$node_random[nodes_all$window == 5])),
  num_nodes_window6 = length(unique(nodes_all$node_random[nodes_all$window == 6])),
  num_nodes_window7 = length(unique(nodes_all$node_random[nodes_all$window == 7])),
  # number of nodes in all preceding time windows for node age indexing
  num_nodes_prev_windows = c(0,
                             length(which(nodes_all$window < 2)),
                             length(which(nodes_all$window < 3)),
                             length(which(nodes_all$window < 4)),
                             length(which(nodes_all$window < 5)),
                             length(which(nodes_all$window < 6)),
                             length(which(nodes_all$window < 7))),
  # centrality means per time window
  centrality_mu_1 = nodes_all$mean_eigen[nodes_all$window == 1],
  centrality_mu_2 = nodes_all$mean_eigen[nodes_all$window == 2],
  centrality_mu_3 = nodes_all$mean_eigen[nodes_all$window == 3],
  centrality_mu_4 = nodes_all$mean_eigen[nodes_all$window == 4],
  centrality_mu_5 = nodes_all$mean_eigen[nodes_all$window == 5],
  centrality_mu_6 = nodes_all$mean_eigen[nodes_all$window == 6],
  centrality_mu_7 = nodes_all$mean_eigen[nodes_all$window == 7],
  # covariance matrix per time window
  centrality_cov_1 = covs_all[[1]],
  centrality_cov_2 = covs_all[[2]],
  centrality_cov_3 = covs_all[[3]],
  centrality_cov_4 = covs_all[[4]],
  centrality_cov_5 = covs_all[[5]],
  centrality_cov_6 = covs_all[[6]],
  centrality_cov_7 = covs_all[[7]],
  # node IDs for all time windows
  node_id = nodes_all$node_random,
  # exposure variable
  node_age = nodes_all$age_std)

intercept  <- rnorm(n, LaplacesDemon::logit(0.05), 2)
min_raw <- min(unlist(model_data_list[grep('centrality_mu', names(model_data_list))]))
max_raw <- max(unlist(model_data_list[grep('centrality_mu', names(model_data_list))]))
plot(NULL, las = 1, xlab = 'age (standardised)', ylab = 'eigenvector',
     #ylim = c(min_raw-2, max_raw+2),
     ylim = c(-15, 5),
     xlim = c(min(model_data_list$node_age),
              max(model_data_list$node_age)))
abline(h = min_raw, lty = 2) ; abline(h = max_raw, lty = 2)
for(i in 1:n){
  lines(x = seq(min(nodes_all$age_std), max(nodes_all$age_std), length.out = 2),
        y = intercept[i] + beta_age[i]*c(min(nodes_all$age_std),
                                         max(nodes_all$age_std)) + window[i] + node[i],
        col = rgb(0,0,1,0.4))
}

## long format simulation
sim_alg <- data.frame(row_num = rep(1:n, each = 2),
                      x = rep(c(min(nodes_all$age_std), max(nodes_all$age_std)), n),
                      intercept = rep(intercept, each = 2),
                      beta_age = rep(beta_age, each = 2)) %>% 
  mutate(y = intercept + beta_age * x)
(long <- ggplot(sim_alg)+
    geom_line(aes(x = x, y = y, group = row_num),
              colour = '#440154FF')+
    geom_hline(yintercept = min_raw,
               linetype = 2, linewidth = 1)+
    geom_hline(yintercept = max_raw,
               linetype = 2, linewidth = 1)+
    labs(x = 'age (standardised)',
         y = 'logit eigenvector centrality'))
sim_alg <- sim_alg %>% 
  mutate(population = 'Amboseli (long)',
         min_raw = min_raw,
         max_raw = max_raw)
rm(list = ls()[!ls() %in% c('n','beta_age','window','node','motnp','short','long','sim_mot','sim_ash','sim_alg')]) ; gc()

#### combine ####
(motnp + short + long) +
  plot_annotation(tag_levels = 'a')
ggsave(plot = last_plot(),
       filename = 'all_priors.png',
       path = '../outputs/step4_nodalregression/',
       device = 'png', height = 800, width = 2400, unit = 'px')

motnp2 <- motnp+
  labs(title = 'Mosi-Oa-Tunya')
short2 <- short+
  labs(title = 'Amboseli (short)')
long2 <- long+
  labs(title = 'Amboseli (long)')
(motnp2 + short2 + long2) +
  plot_annotation(tag_levels = 'a')
ggsave(plot = last_plot(),
       filename = 'all_priors.png',
       path = '../outputs/step4_nodalregression/',
       device = 'png', height = 800, width = 2400, unit = 'px')

ggsave(plot = motnp2,
       filename = 'motnp_prior.png',
       path = '../outputs/step4_nodalregression/',
       device = 'png', height = 800, width = 750, unit = 'px')
ggsave(plot = short2,
       filename = 'anpshort_prior.png',
       path = '../outputs/step4_nodalregression/',
       device = 'png', height = 800, width = 750, unit = 'px')
ggsave(plot = long2,
       filename = 'anplong_prior.png',
       path = '../outputs/step4_nodalregression/',
       device = 'png', height = 800, width = 750, unit = 'px')


sim_dat <- rbind(sim_mot, sim_ash, sim_alg) %>% 
  mutate(population = factor(population,
                             levels = c('Mosi-Oa-Tunya',
                                        'Amboseli (short)',
                                        'Amboseli (long)')))

raw_data <- sim_dat %>% 
  dplyr::select(population, min_raw, max_raw) %>% 
  distinct() %>% 
  mutate(population = factor(population,
                             levels = c('Mosi-Oa-Tunya',
                                        'Amboseli (short)',
                                        'Amboseli (long)')))

ggplot()+
    geom_line(data = sim_dat,
              mapping = aes(x = x, y = y,
                            group = row_num,
                            colour = population),
              linewidth = 0.2)+
    geom_hline(data = raw_data,
               mapping = aes(yintercept = min_raw),
               linetype = 2, linewidth = 1)+
    geom_hline(data = raw_data,
               mapping = aes(yintercept = max_raw),
               linetype = 2, linewidth = 1)+
    labs(x = 'age (standardised)',
         y = 'logit eigenvector centrality')+
    facet_wrap(. ~ population)+
    scale_colour_viridis_d(direction = -1)+
    theme(legend.position = 'bottom')
ggsave(plot = last_plot(),
       filename = 'all_priors2.png',
       path = '../outputs/step4_nodalregression/',
       device = 'png', height = 1200, width = 2400, unit = 'px')

