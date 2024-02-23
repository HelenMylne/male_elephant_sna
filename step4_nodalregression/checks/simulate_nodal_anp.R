#### set up ####
# library(StanHeaders) ; library(rstan) ; library(LaplacesDemon) ; library(tidyverse)
library(StanHeaders, lib.loc = '../packages/')
library(rstan, lib.loc = '../packages/')
#library(cmdstanr, lib.loc = '../packages/')
#set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0/')
library(LaplacesDemon, lib.loc = '../packages/')
library(tidyverse, lib.loc = '../packages/')

# pdf('step4_nodalregression/checks/simulation_anp_modelprep.pdf')
theme_set(theme_bw())
set.seed(2)

#### derive population information directly from ANP data ####
sim <- read_csv('../data_processed/step4_nodalregression/anp_allnodes.csv') %>%
  select(-mean_eigen) %>%
  relocate(age_std, .after = age) %>%
  relocate(node_random, .after = node)

## define parameters
n_data <- nrow(sim)
n_nodes <- length(unique(sim$node))
n_windows <- length(unique(sim$window))

## view ages
hist(sim$age)

#### simulate centralities ####
# simulate age effect
sim_slope <- -0.02

## simulate intercept
sim_intcp <- -5

## simulate random effects -- these are just random draws from a distribution centred on zero, but they do not have to have a mean of 0 between them. Should they? (i.e. if random effect of windows 1 and 2 are both positive, does window 3 need to be negative to counteract it?)
sim_window_unique <- rnorm(n_windows, mean = 0, sd = 1)
sim_node_unique <- rnorm(n_nodes, mean = 0, sd = 0.2)

## calculate mean centrality per elephant
sim$mu <- sim$age * sim_slope + sim_intcp + sim_window_unique[sim$window] + sim_node_unique[sim$node_random]    # simulate mean centrality on outcome scale

## plot mean centrality against age
plot(mu ~ age, data = sim[sim$window == 1,], col = 'red', pch = 19, ylim = c(-5,1), las = 1)
points(mu ~ age, data = sim[sim$window == 2,], col = 'blue', pch = 19)
points(mu ~ age, data = sim[sim$window == 3,], col = 'green', pch = 19)

## check inputs
colour <- rep(c('red','orange','yellow','green','blue','purple'), each = 6)
shapes <- rep(c(3,4,15,16,17,18),6)
plot(mu ~ age_std, data = sim[sim$window == 1,],
     pch = shapes[1], col = colour[1], las = 1,
     ylim = c(-10,0), xlim = c(-2.5,5))
for(i in 2:n_windows){
  points(mu ~ age_std, data = sim[sim$window == i,],
         pch = shapes[i], col = colour[i])
}

## simulate full distribution of samples per node
sim$sd <- 1#(rpois(nrow(sim),lambda = 1)+1)/10         # can be made to vary amongst nodes, currently all elephants have equal variance in centrality (not realistic given that in the real data some are seen more often than others)
sim_dat <- matrix(data = NA, nrow = 4000, ncol = n_data, dimnames = list(NULL, sim$node_window))    # create matrix for full centrality distribution
for(j in 1:n_data){
  sim_dat[,j] <- rnorm(n = nrow(sim_dat), mean = sim$mu[j], sd = sim$sd[j])  # simulate distribution for each elephant
}

## visualise
data.frame(sim_dat) %>%
  select(sample(sim$node_random, 40, replace = F)) %>%
  pivot_longer(cols = everything(),
               names_to = "node_random", values_to = "centrality") %>%
  separate(node_random, into = c('X','node_window'), remove = T, sep = 1) %>%
  dplyr::select(-X) %>%
  separate(node_window, into = c('node_random','window'), sep = '_', remove = F) %>%
  mutate(node = as.integer(node_random)) %>%
  left_join(sim[,c('node_window','age')], by = 'node_window') %>%
  ggplot(aes(x = centrality, fill = age, group = node_window)) +
  geom_density(linewidth = 0.4, alpha = 0.6) +
  facet_grid(rows = vars(as.factor(node)), scales = "free") +
  labs(x = "Eigenvector centrality (standardised)") +
  theme_void() +
  theme(strip.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

#### prep inputs ####
eigen_means <- list()
eigen_covs <- list()

# compute normal approximation by window -- calculate means per node
for(time_window in 1:n_windows){
  eigen_means[[time_window]] <- apply(sim_dat[,which(sim$window == time_window)], 2, mean)
  eigen_covs[[time_window]] <- cov(sim_dat[,which(sim$window == time_window)])
}

## check normal approximation -- simulate from combined mean and covariance, plot curve against relative node ID
par(mfrow = c(6,6), mai = c(0.1,0.1,0.1,0.1))
for(time_window in 1:n_windows){
  ## simulate from multivariate normal
  sim_cent_samples <- MASS::mvrnorm(1e5, eigen_means[[time_window]], eigen_covs[[time_window]])

  ## identify random node of interest
  node_id_sample <- sample(which(sim$window == time_window),1)
  nodes_timewindow <- sim %>%
    filter(window == time_window)
  node_id_check <- which(nodes_timewindow$node == sim$node[node_id_sample])

  ## plot true density curve
  plot(density(sim_dat[,node_id_sample]), lwd = 2, las = 1,
       main = paste0('time window ', time_window),
       xlab = '', ylab = '')

  ## plot normal approximation
  lines(density(sim_cent_samples[, node_id_check]),
        col = rgb(0,0,1,0.5), lwd = 2)
}
par(mfrow = c(1,1), mai = c(0.5,0.5,0.5,0.5))
rm(sim_cent_samples, node_id_sample, node_id_check, nodes_timewindow) ; gc()

## create data
nodes_per_window <- as.data.frame(table(sim$window)) %>%
  rename(window = Var1,
         node_count = Freq)

eigen_list <- list(
  # global data size
  num_data = n_data,
  num_nodes = n_nodes,
  num_windows = n_windows,
  # per time window data size
  num_nodes_window1 = nodes_per_window$node_count[1],
  num_nodes_window2 = nodes_per_window$node_count[2],
  num_nodes_window3 = nodes_per_window$node_count[3],
  num_nodes_window4 = nodes_per_window$node_count[4],
  num_nodes_window5 = nodes_per_window$node_count[5],
  num_nodes_window6 = nodes_per_window$node_count[6],
  num_nodes_window7 = nodes_per_window$node_count[7],
  num_nodes_window8 = nodes_per_window$node_count[8],
  num_nodes_window9 = nodes_per_window$node_count[9],
  num_nodes_window10 = nodes_per_window$node_count[10],
  num_nodes_window11 = nodes_per_window$node_count[11],
  num_nodes_window12 = nodes_per_window$node_count[12],
  num_nodes_window13 = nodes_per_window$node_count[13],
  num_nodes_window14 = nodes_per_window$node_count[14],
  num_nodes_window15 = nodes_per_window$node_count[15],
  num_nodes_window16 = nodes_per_window$node_count[16],
  num_nodes_window17 = nodes_per_window$node_count[17],
  num_nodes_window18 = nodes_per_window$node_count[18],
  num_nodes_window19 = nodes_per_window$node_count[19],
  num_nodes_window20 = nodes_per_window$node_count[20],
  num_nodes_window21 = nodes_per_window$node_count[21],
  num_nodes_window22 = nodes_per_window$node_count[22],
  num_nodes_window23 = nodes_per_window$node_count[23],
  num_nodes_window24 = nodes_per_window$node_count[24],
  num_nodes_window25 = nodes_per_window$node_count[25],
  num_nodes_window26 = nodes_per_window$node_count[26],
  num_nodes_window27 = nodes_per_window$node_count[27],
  num_nodes_window28 = nodes_per_window$node_count[28],
  num_nodes_window29 = nodes_per_window$node_count[29],
  num_nodes_window30 = nodes_per_window$node_count[30],
  num_nodes_window31 = nodes_per_window$node_count[31],
  num_nodes_window32 = nodes_per_window$node_count[32],
  num_nodes_window33 = nodes_per_window$node_count[33],
  num_nodes_window34 = nodes_per_window$node_count[34],
  num_nodes_window35 = nodes_per_window$node_count[35],
  num_nodes_window36 = nodes_per_window$node_count[36],
  # number of nodes in all preceding time windows for node age indexing
  num_nodes_prev_windows = c(0,
                             length(which(sim$window < 2)),
                             length(which(sim$window < 3)),
                             length(which(sim$window < 4)),
                             length(which(sim$window < 5)),
                             length(which(sim$window < 6)),
                             length(which(sim$window < 7)),
                             length(which(sim$window < 8)),
                             length(which(sim$window < 9)),
                             length(which(sim$window < 10)),
                             length(which(sim$window < 11)),
                             length(which(sim$window < 12)),
                             length(which(sim$window < 13)),
                             length(which(sim$window < 14)),
                             length(which(sim$window < 15)),
                             length(which(sim$window < 16)),
                             length(which(sim$window < 17)),
                             length(which(sim$window < 18)),
                             length(which(sim$window < 19)),
                             length(which(sim$window < 20)),
                             length(which(sim$window < 21)),
                             length(which(sim$window < 22)),
                             length(which(sim$window < 23)),
                             length(which(sim$window < 24)),
                             length(which(sim$window < 25)),
                             length(which(sim$window < 26)),
                             length(which(sim$window < 27)),
                             length(which(sim$window < 28)),
                             length(which(sim$window < 29)),
                             length(which(sim$window < 30)),
                             length(which(sim$window < 31)),
                             length(which(sim$window < 32)),
                             length(which(sim$window < 33)),
                             length(which(sim$window < 34)),
                             length(which(sim$window < 35)),
                             length(which(sim$window < 36))),
  # centrality means per time window
  centrality_mu_1 = eigen_means[[1]],
  centrality_mu_2 = eigen_means[[2]],
  centrality_mu_3 = eigen_means[[3]],
  centrality_mu_4 = eigen_means[[4]],
  centrality_mu_5 = eigen_means[[5]],
  centrality_mu_6 = eigen_means[[6]],
  centrality_mu_7 = eigen_means[[7]],
  centrality_mu_8 = eigen_means[[8]],
  centrality_mu_9 = eigen_means[[9]],
  centrality_mu_10 = eigen_means[[10]],
  centrality_mu_11 = eigen_means[[11]],
  centrality_mu_12 = eigen_means[[12]],
  centrality_mu_13 = eigen_means[[13]],
  centrality_mu_14 = eigen_means[[14]],
  centrality_mu_15 = eigen_means[[15]],
  centrality_mu_16 = eigen_means[[16]],
  centrality_mu_17 = eigen_means[[17]],
  centrality_mu_18 = eigen_means[[18]],
  centrality_mu_19 = eigen_means[[19]],
  centrality_mu_20 = eigen_means[[20]],
  centrality_mu_21 = eigen_means[[21]],
  centrality_mu_22 = eigen_means[[22]],
  centrality_mu_23 = eigen_means[[23]],
  centrality_mu_24 = eigen_means[[24]],
  centrality_mu_25 = eigen_means[[25]],
  centrality_mu_26 = eigen_means[[26]],
  centrality_mu_27 = eigen_means[[27]],
  centrality_mu_28 = eigen_means[[28]],
  centrality_mu_29 = eigen_means[[29]],
  centrality_mu_30 = eigen_means[[30]],
  centrality_mu_31 = eigen_means[[31]],
  centrality_mu_32 = eigen_means[[32]],
  centrality_mu_33 = eigen_means[[33]],
  centrality_mu_34 = eigen_means[[34]],
  centrality_mu_35 = eigen_means[[35]],
  centrality_mu_36 = eigen_means[[36]],
  # covariance matrix per time window
  centrality_cov_1 = eigen_covs[[1]],
  centrality_cov_2 = eigen_covs[[2]],
  centrality_cov_3 = eigen_covs[[3]],
  centrality_cov_4 = eigen_covs[[4]],
  centrality_cov_5 = eigen_covs[[5]],
  centrality_cov_6 = eigen_covs[[6]],
  centrality_cov_7 = eigen_covs[[7]],
  centrality_cov_8 = eigen_covs[[8]],
  centrality_cov_9 = eigen_covs[[9]],
  centrality_cov_10 = eigen_covs[[10]],
  centrality_cov_11 = eigen_covs[[11]],
  centrality_cov_12 = eigen_covs[[12]],
  centrality_cov_13 = eigen_covs[[13]],
  centrality_cov_14 = eigen_covs[[14]],
  centrality_cov_15 = eigen_covs[[15]],
  centrality_cov_16 = eigen_covs[[16]],
  centrality_cov_17 = eigen_covs[[17]],
  centrality_cov_18 = eigen_covs[[18]],
  centrality_cov_19 = eigen_covs[[19]],
  centrality_cov_20 = eigen_covs[[20]],
  centrality_cov_21 = eigen_covs[[21]],
  centrality_cov_22 = eigen_covs[[22]],
  centrality_cov_23 = eigen_covs[[23]],
  centrality_cov_24 = eigen_covs[[24]],
  centrality_cov_25 = eigen_covs[[25]],
  centrality_cov_26 = eigen_covs[[26]],
  centrality_cov_27 = eigen_covs[[27]],
  centrality_cov_28 = eigen_covs[[28]],
  centrality_cov_29 = eigen_covs[[29]],
  centrality_cov_30 = eigen_covs[[30]],
  centrality_cov_31 = eigen_covs[[31]],
  centrality_cov_32 = eigen_covs[[32]],
  centrality_cov_33 = eigen_covs[[33]],
  centrality_cov_34 = eigen_covs[[34]],
  centrality_cov_35 = eigen_covs[[35]],
  centrality_cov_36 = eigen_covs[[36]],
  # node IDs for all time windows
  nodes_window1 = sim$node_random[sim$window == 1],
  nodes_window2 = sim$node_random[sim$window == 2],
  nodes_window3 = sim$node_random[sim$window == 3],
  nodes_window4 = sim$node_random[sim$window == 4],
  nodes_window5 = sim$node_random[sim$window == 5],
  nodes_window6 = sim$node_random[sim$window == 6],
  nodes_window7 = sim$node_random[sim$window == 7],
  nodes_window8 = sim$node_random[sim$window == 8],
  nodes_window9 = sim$node_random[sim$window == 9],
  nodes_window10 = sim$node_random[sim$window == 10],
  nodes_window11 = sim$node_random[sim$window == 11],
  nodes_window12 = sim$node_random[sim$window == 12],
  nodes_window13 = sim$node_random[sim$window == 13],
  nodes_window14 = sim$node_random[sim$window == 14],
  nodes_window15 = sim$node_random[sim$window == 15],
  nodes_window16 = sim$node_random[sim$window == 16],
  nodes_window17 = sim$node_random[sim$window == 17],
  nodes_window18 = sim$node_random[sim$window == 18],
  nodes_window19 = sim$node_random[sim$window == 19],
  nodes_window20 = sim$node_random[sim$window == 20],
  nodes_window21 = sim$node_random[sim$window == 21],
  nodes_window22 = sim$node_random[sim$window == 22],
  nodes_window23 = sim$node_random[sim$window == 23],
  nodes_window24 = sim$node_random[sim$window == 24],
  nodes_window25 = sim$node_random[sim$window == 25],
  nodes_window26 = sim$node_random[sim$window == 26],
  nodes_window27 = sim$node_random[sim$window == 27],
  nodes_window28 = sim$node_random[sim$window == 28],
  nodes_window29 = sim$node_random[sim$window == 29],
  nodes_window30 = sim$node_random[sim$window == 30],
  nodes_window31 = sim$node_random[sim$window == 31],
  nodes_window32 = sim$node_random[sim$window == 32],
  nodes_window33 = sim$node_random[sim$window == 33],
  nodes_window34 = sim$node_random[sim$window == 34],
  nodes_window35 = sim$node_random[sim$window == 35],
  nodes_window36 = sim$node_random[sim$window == 36],
  # exposure variable
  node_age = sim$age_std)

#### prior predictive check ####
n <- 100
intercept  <- rnorm(n, -6, 2) #rnorm(n, LaplacesDemon::logit(0.05), 3) # taking the intercept from the logit of results from Chiyo 2011 (doesn't state mean/median centrality so estimated from graph based on where correlation line would cross x = 0)
beta_age <- rnorm(n, 0, 0.8)   # beta_age <- rnorm(n, 0, 0.8)
sigma_window <- rexp(n, 3)
sigma_node <- rexp(n, 3)
rand_window <- rnorm(n, 0, sigma_window)
rand_node <- rnorm(n, 0,sigma_node)
min_raw <- min(unlist(eigen_list[grep('centrality_mu', names(eigen_list))]))
max_raw <- max(unlist(eigen_list[grep('centrality_mu', names(eigen_list))]))
plot(NULL, las = 1, xlab = 'age (standardised)', ylab = 'eigenvector',
     ylim = c(-15, 5),
     xlim = c(min(eigen_list$node_age),
              max(eigen_list$node_age)))
abline(h = min_raw, lty = 2) ; abline(h = max_raw, lty = 2)
for(i in 1:n){
  lines(x = seq(min(sim$age_std), max(sim$age_std), length.out = 2),
        y = intercept[i] + beta_age[i]*c(min(sim$age_std), max(sim$age_std)) + rand_window[i] + rand_node[i],
        col = rgb(0,0,1,0.4))
}
rm(n, max_raw, min_raw, beta_age, intercept) ; gc()

dev.off()

#### run model -- age as a continuous variable with 2 windows ####
## load model
nodal_regression <- stan_model('models/eigen_regression_anp.stan')
#nodal_regression <- cmdstan_model('models/eigen_regression_anp.stan')

## run model
n_chains <- 4
n_samples <- 1000
fit_sim <- sampling(nodal_regression, data = eigen_list,
                    chains = n_chains, cores = n_chains)
# fit_sim <- nodal_regression$sample(data = eigen_list,
#                                    chains = n_chains, parallel_chains = n_chains,
#                                    iter_warmup = n_samples, iter_sampling = n_samples)
# temp_rds_file <- tempfile(fileext = '.RDS')
# fit_sim$save.object(file = temp_rds_file)

save.image('anp_nodalregression/simulation.RData')

#### check outputs ####
#load('anp_nodalregression/simulation.RData')
#fit_sim <- readRDS(temp_rds_file)
pdf('step4_nodalregression/checks/simulation_anp_modelchecks.pdf')

## extract model fit -- very good!
# ( summary <- fit_sim$summary() )
# par(mfrow = c(3,1))
# hist(summary$rhat,#[grep(pattern = 'chol_',
#      #x = summary$variable,
#      #invert = T)],
#      breaks = 50)
# hist(summary$ess_bulk,#[grep(pattern = 'chol_',
#      #x = summary$variable,
#      #invert = T)],
#      breaks = 50)
# hist(summary$ess_tail,#[grep(pattern = 'chol_',
#      #x = summary$variable,
#      #invert = T)],
#      breaks = 50)
summary <- summary(fit_sim)
fit_summary <- as.data.frame(summary$summary)
par(mfrow = c(2,1))
hist(fit_summary$Rhat[grep(pattern = 'predictor_window', x = rownames(fit_summary),invert = T)],
     breaks = 50, main = 'Rhat distribution')
hist(fit_summary$n_eff[grep(pattern = 'predictor_window', x = rownames(fit_summary), invert = T)],
     breaks = 50, main = 'n_eff distribution')
par(mfrow = c(1,1))

## extract posterior
# params_std <- fit_sim$draws(format = 'draws_df')
params_std <- rstan::extract(fit_sim) %>%
  as.data.frame()

## separate random effects from global parameters
rand_window <- params_std %>%
  dplyr::select(grep('rand_window', colnames(params_std), value=TRUE))
check_windows <- apply(rand_window, 2, mean) %>%
  as.data.frame() %>%
  mutate(node = colnames(rand_window)) %>%
  relocate(node) %>%
  mutate(sim_input = sim_window_unique)
colnames(check_windows)[2] <- 'model_output'
plot(check_windows$model_output ~ check_windows$sim_input,
     main = 'model calculated effect of window vs\nactual simulated effect of window')

rand_node <- params_std %>%
  dplyr::select(grep('rand_node', colnames(params_std), value=TRUE))
check_nodes <- apply(rand_node, 2, mean) %>%
  as.data.frame() %>%
  mutate(node = colnames(rand_node)) %>%
  relocate(node) %>%
  mutate(sim_input = sim_node_unique)
colnames(check_nodes)[2] <- 'model_output'
plot(check_nodes$model_output ~ check_nodes$sim_input,
     main = 'model calculated effect of node ID vs\nactual simulated effect of node ID')   # wouldn't expect this to be an actual 1:1 relationship because of the scales, but there should be a positive relationship

rm(check_windows, check_nodes) ; gc()

## traceplot all parameters
# plot_params <- c('intercept','beta_age','sigma',
#                  'rand_node[1]','rand_node[50]','rand_node[100]','rand_node[150]',
#                  'rand_window[1]','rand_window[2]','rand_window[3]','rand_window[4]',
#                  'predictor_window1[1]','predictor_window1[16]',
#                  'predictor_window2[1]','predictor_window2[2]')
# params_std %>%
#   select(all_of(plot_params),`.draw`,`.chain`,`.iteration`) %>%
#   pivot_longer(cols = all_of(plot_params), names_to = 'parameter', values_to = 'draw') %>%
#   rename(chain = .chain,
#          chain_position = .iteration,
#          draw_id = .draw) %>%
#   #filter(chain == 1) %>% # inspect individual chains -- check for wandery sections that might be hidden by other chains when all plotted together
#   ggplot(aes(x = chain_position, y = draw, colour = as.factor(chain)))+
#   geom_line()+
#   facet_wrap(. ~ parameter, scales = 'free_y')+
#   theme_bw()+
#   theme(legend.position = 'none')
# rm(plot_params) ; gc()

traceplot(fit_sim, pars = 'intercept')
traceplot(fit_sim, pars = 'beta_age')
traceplot(fit_sim, pars = 'sigma')

#plot_rand_nodes <- sample(colnames(rand_node), size = 12, replace = F)
plot_rand_nodes <- c('rand_node[1]','rand_node[60]','rand_node[119]','rand_node[179]','rand_node[238]','rand_node[298]','rand_node[357]','rand_node[417]','rand_node[476]','rand_node[536]','rand_node[595]','rand_node[655]')
traceplot(fit_sim, pars = plot_rand_nodes)

plot_rand_windows1 <- c('rand_window[1]','rand_window[2]','rand_window[3]','rand_window[4]','rand_window[5]','rand_window[6]','rand_window[7]','rand_window[8]','rand_window[9]')
traceplot(fit_sim, pars = plot_rand_windows1)
plot_rand_windows2 <- c('rand_window[10]','rand_window[11]','rand_window[12]','rand_window[13]','rand_window[14]','rand_window[15]','rand_window[16]','rand_window[17]','rand_window[18]')
traceplot(fit_sim, pars = plot_rand_windows2)
plot_rand_windows3 <- c('rand_window[19]','rand_window[20]','rand_window[21]','rand_window[22]','rand_window[23]','rand_window[24]','rand_window[25]','rand_window[26]','rand_window[27]')
traceplot(fit_sim, pars = plot_rand_windows3)
plot_rand_windows4 <- c('rand_window[28]','rand_window[29]','rand_window[30]','rand_window[31]','rand_window[32]','rand_window[33]','rand_window[34]','rand_window[35]','rand_window[36]')
traceplot(fit_sim, pars = plot_rand_windows4)
rm(plot_params, plot_rand_nodes, plot_rand_windows1, plot_rand_windows2, plot_rand_windows3, plot_rand_windows4) ; gc()

#### posterior predictive check ####
par(mfrow = c(3,3))

## create check function
ppcheck <- function(eigen_mat, eigen_df, cent_cov, window, params, rand_node, rand_window){
  plot(density(eigen_mat[1, which(eigen_df$window == window)]), las = 1, ylim = c(0,1),
       main = "Posterior predictive check:\nblack = data, blue = predicted",
       #main = "Posterior predictive check (standardised centrality):\nblack = data, blue = predicted",
       col=rgb(0, 0, 0, 0.25))
  eigen_data <- list(num_nodes_window = length(which(eigen_df$window == window)),
                     node_age = eigen_df$age_std[eigen_df$window == window],
                     window = window,
                     nodes = eigen_df$node_random[eigen_df$window == window],
                     nodes_window = eigen_df$node_random[eigen_df$window == window])
  for (i in 1:100) {
    j <- sample(1:length(params$beta_age), 1)
    lines(density(eigen_mat[j, which(eigen_df$window == window)]), col=rgb(0, 0, 0, 0.25))
    mu <- params$beta_age[j]*eigen_data$node_age + params$intercept[j]
    for(k in 1:length(mu)) {
      mu[k] <- mu[k] + as.numeric(rand_window[j,eigen_data$window]) + as.numeric(rand_node[j,eigen_data$nodes[k]])
    }
    sigma <- cent_cov + diag(rep(params$sigma[j], eigen_data$num_nodes_window))
    lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
  }
}

## check per time window
for(time_window in 1:n_windows){
  ppcheck(eigen_mat = sim_dat, eigen_df = sim,
          cent_cov = eigen_covs[[time_window]], window = time_window,
          params = params_std, rand_node = rand_node, rand_window = rand_window)
}
par(mfrow = c(1,1))
dev.off()

#### predict from model -- standardised scale ####
pdf('step4_nodalregression/checks/simulation_anp_predictions.pdf')

## create mean prediction function
get_mean_predictions <- function(predict_df, parameters, include_node = TRUE, include_window = TRUE){
  ## create empty matrix to fill with predictions
  mean_matrix <- matrix(NA, nrow = nrow(parameters), ncol = nrow(predict_df),
                        dimnames = list(NULL, predict_df$node_window))

  ## populate matrix = mean centrality values per node, predicting for real data
  for(i in 1:nrow(mean_matrix)){
    mean_matrix[i,] <- parameters$beta_age[i] * predict_df$age_std + parameters$intercept[i]
  }

  if(include_window == TRUE){
    window_effect <- parameters %>% dplyr::select(grep('rand_window', colnames(parameters), value=TRUE))
    for(i in 1:nrow(mean_matrix)){
      mean_matrix[i,] <- mean_matrix[i,] + as.numeric(window_effect[i,predict_df$window])
    }
  }

  if(include_node == TRUE){
    node_effect <- parameters %>% dplyr::select(grep('rand_node', colnames(parameters), value=TRUE))
    for(i in 1:nrow(mean_matrix)){
      mean_matrix[i,] <- mean_matrix[i,] + as.numeric(node_effect[i, predict_df$node_random])
    }
  }

  return(mean_matrix)

}

## get mean predictions
mu_std <- get_mean_predictions(predict_df = sim, parameters = params_std, include_window = TRUE, include_node = TRUE)

## add mean and CI of predicted means to input data frame for comparison
sim$mu_mean_std <- apply(mu_std, 2, mean)
sim$mu_lwr_std <- apply(mu_std, 2, quantile, prob = 0.025)
sim$mu_upr_std <- apply(mu_std, 2, quantile, prob = 0.975)

## plot mean of model vs mean of raw data
ggplot()+
  geom_point(data = sim, aes(x = mu, y = mu_mean_std, colour = as.factor(window)))+
  scale_colour_viridis_d()+
  labs(colour = 'window', x = 'simulated mean (standardised)', y = 'predicted mean (standardised)')+
  geom_abline(slope = 1, intercept = 0) # add line showing where points would lie if model fit was perfect

## put together sigma arrays, separated by time window
num_nodes_windows <- unlist(eigen_list[grep('num_nodes_window', names(eigen_list))])
sigma_list <- list()
for(time_window in 1:n_windows) {
  sigma_window <- array(NA, dim = c(num_nodes_windows[time_window],
                                    num_nodes_windows[time_window],
                                    nrow(params_std)),
                        dimnames = list(sim$node_random[sim$window == time_window],
                                        sim$node_random[sim$window == time_window],
                                        NULL))
  cent_cov <- eigen_covs[[time_window]]
  for(i in 1:nrow(params_std)){
    sigma_window[,,i] <- cent_cov + diag(rep(params_std$sigma[i], num_nodes_windows[time_window]))
  }
  sigma_list[[time_window]] <- sigma_window
}

## create empty matrix to take full set of predicted values per elephant
predictions_std <- matrix(NA, nrow = nrow(params_std), ncol = nrow(sim), dimnames = list(NULL, sim$node_window))

## populate matrix using mean values in matrix mu_std, and sigma values based on time window
for(time_window in 1:n_windows){
  sigma_window <- sigma_list[[time_window]]
  columns <- which(sim$window == time_window)
  for(i in 1:nrow(predictions_std)){
    predictions_std[i,columns] <- MASS::mvrnorm(1,
                                                mu_std[i,columns],
                                                sigma_window[,,i])
  }
}

## add CI of predicted data points to input data frame for comparison
sim$predict_lwr_std <- apply(predictions_std, 2, quantile, prob = 0.025)
sim$predict_upr_std <- apply(predictions_std, 2, quantile, prob = 0.025)

## plot predictions
ggplot(sim)+
  geom_ribbon(aes(x = age, ymin = predict_lwr_std, ymax = predict_upr_std, fill = as.factor(window)),
              alpha = 0.2)+                      # background layer showing the 95% CI of all predictions
  geom_ribbon(aes(x = age, ymin = mu_lwr_std, ymax = mu_upr_std, fill = as.factor(window)),
              alpha = 0.4)+                      # mid layer showing the 95% CI of predicted means
  geom_line(aes(x = age, y = mu_mean_std, colour = as.factor(window)))+  # line showing mean of predicted means
  geom_point(aes(x = age, y = mu))+              # original data points (standardised centrality, actual age)
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  facet_wrap(. ~ as.factor(window))+             # separate plots per window
  theme(legend.position = 'bottom')+
  labs(colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'age (years)')

save.image('anp_nodalregression/simulation.RData')

#### extract original values from output -- simulated slope value originally used produces an effect on the unstandardised scale. The model works on the standardised scale. Convert predictions to unstandardised scale and then run contrasts to calculate the slope coefficient. ####
pdf('step4_nodalregression/checks/simulation_anp_extractoriginal.pdf')
#load('anp_nodalregression/simulation.RData')

## predict means from model again, now using age_std + 1 stdev -- for now working with age_std+1SD not age+1yr because that's how I originally simulated the data. if change the simulation so that mean centrality = age*slope rather age_std*slope, then change this to extract the correct slope value
## get mean predictions
sim2 <- sim %>% 
  mutate(age_std = age_std+1)
mu2_std <- get_mean_predictions(predict_df = sim2, parameters = params_std, include_window = TRUE, include_node = TRUE)
sim$mu_mean_plus1 <- apply(mu2_std, 2, mean)

## full distribution of predictions using age_std + 1 stdev
predictions2_std <- matrix(NA, nrow = nrow(params_std), ncol = nrow(sim2), dimnames = list(NULL, sim$node_window))
for(time_window in 1:n_windows){
  sigma_window <- sigma_list[[time_window]]
  columns <- which(sim$window == time_window)
  for(i in 1:nrow(predictions2_std)){
    predictions2_std[i,columns] <- MASS::mvrnorm(1, mu2_std[i,columns], sigma_window[,,i])
  }
}

## contrast predictions on standardised scale -- check that this returns the marginal effect presented in the summary
contrast_std <- predictions2_std - predictions_std    # contrast between predicted values for raw data and all same data but add 1 to age
head(contrast_std[,1:5])                              # check matrix looks right
mean(params_std$beta_age)                             # output parameter direct from model
mean(contrast_std)                                    # should be very similar to mean of output parameter
quantile(contrast_std, prob = c(0.025, 0.975))        # very wide

## contrast predictions on unstandardised scale
contrast <- contrast_std / sd(sim$age)                # contrast between predicted values for raw data and all same data but add 1 to age
head(contrast[,1:5])                                  # check matrix looks right
sim_slope                                             # output parameter direct from model
mean(contrast)                                        # should be very similar to mean of output parameter
quantile(contrast, prob = c(0.025, 0.975))            # very wide

#### final "clean" plots using hypothetical data ####
## create dummy dataset
# newdat <- expand.grid(window = 1:n_windows,
#                       age = sort(unique(sim$age))) %>% 
#   left_join(distinct(sim[,c('age','age_std')]), by = 'age')
newdat <- sim %>% 
  select(node_random, age, age_std, window)

## get mean predictions
fake_mean <- get_mean_predictions(predict_df = newdat, parameters = params_std, include_window = TRUE, include_node = FALSE)
newdat$predict_mean <- apply(fake_mean, 2, mean)
newdat$predict_mean_lwr <- apply(fake_mean, 2, quantile, prob = 0.025)
newdat$predict_mean_upr <- apply(fake_mean, 2, quantile, prob = 0.975)

## create prediction matrix
fake_pred <- matrix(NA, nrow = nrow(params_std), ncol = nrow(newdat))
for(time_window in 1:n_windows){
  sigma_window <- sigma_list[[time_window]]
  columns <- which(newdat$window == time_window)
  for(i in 1:nrow(fake_pred)){
    fake_pred[i,columns] <- MASS::mvrnorm(1, fake_mean[i,columns], sigma_window[,,i])
  }
}

## add CI of predicted data points to input data frame for comparison
newdat$predict_pred_lwr <- apply(fake_pred, 2, quantile, prob = 0.025)
newdat$predict_pred_upr <- apply(fake_pred, 2, quantile, prob = 0.975)

## plot predictions
newdat_summary <- newdat %>% 
  group_by(age, window) %>% 
  mutate(predict_pred_lwr = mean(predict_pred_lwr),
         predict_pred_upr = mean(predict_pred_upr)) %>% 
  select(age, predict_pred_lwr, predict_pred_upr,window) %>% 
  distinct()

## plot mean values
ggplot()+
  geom_ribbon(data = newdat_summary,
              aes(x = age, ymin = predict_pred_lwr, ymax = predict_pred_upr, fill = as.factor(window)),
              alpha = 0.2)+                      # background layer showing the 95% CI of all predictions
  geom_ribbon(data = newdat,
              aes(x = age, ymin = predict_mean_lwr, ymax = predict_mean_upr, fill = as.factor(window)),
              alpha = 0.4)+                      # mid layer showing the 95% CI of predicted means
  geom_line(data = newdat,
            aes(x = age, y = predict_mean, colour = as.factor(window)))+  # line showing mean of predicted means
  geom_point(data = sim,
             aes(x = age, y = mu))+              # original data points (standardised centrality, actual age)
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  facet_wrap(. ~ as.factor(window))+             # separate plots per window
  theme(legend.position = 'bottom')+
  labs(colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'age (years)')

## convert full predictions to data frame
sim_dat_df <- sim_dat %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'node_window', values_to = 'eigen')

sim$node_window <- paste0(sim$node_random, '_', sim$window)
sim_dat_df <- sim_dat_df %>% 
  left_join(sim[,c('node_window','age','window')], by = 'node_window')

## plot full distribution
ggplot()+
  geom_ribbon(data = newdat_summary,
              aes(x = age, ymin = predict_pred_lwr, ymax = predict_pred_upr, fill = as.factor(window)),
              alpha = 0.2)+                      # background layer showing the 95% CI of all predictions
  geom_ribbon(data = newdat,
              aes(x = age, ymin = predict_mean_lwr, ymax = predict_mean_upr, fill = as.factor(window)),
              alpha = 0.4)+                      # mid layer showing the 95% CI of predicted means
  geom_point(data = sim_dat_df,
             aes(x = age, y = eigen),
             alpha = 0.01, size = 0.5)+         # original data points (standardised centrality, actual age)
  geom_line(data = newdat,
            aes(x = age, y = predict_mean, colour = as.factor(window)))+  # line showing mean of predicted means
  geom_point(data = sim,
             aes(x = age, y = mu),
             size = 0.5, colour = 'white')+      # original mean data points (standardised centrality, actual age)
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  facet_wrap(. ~ as.factor(window))+             # separate plots per window
  theme(legend.position = 'bottom')+
  labs(colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'age (years)')

dev.off()
save.image('anp_nodalregression/simulation.RData')
