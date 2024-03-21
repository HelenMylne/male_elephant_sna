## set up #### 
library(LaplacesDemon)
library(tidyverse)
library(cmdstanr)
set_cmdstan_path('../packages/.cmdstan2/cmdstan-2.33.1/')

set.seed(1)
pdf('step4_nodalregression/checks/simulation_motnp_nodal.pdf')

## simulate population ####
# ## define population parameters
# min_age <- 11                                                # youngest individual
# max_age <- 60                                                #  oldest  individual
# n_nodes <- 2*((max_age+1) - min_age)                         # total nodes = average 2 per age
# 
# ## simulate population for first time window
# sim <- data.frame(node = 1:n_nodes,                          # create data frame of individuals
#                   age = sample(x = min_age:max_age, size = n_nodes,
#                                replace = T, prob = 1/(min_age:max_age))) %>% 
#   mutate(age_cat = ifelse(age <= 15, 1,
#                           ifelse(age <= 20, 2,
#                                  ifelse(age <= 25, 3,
#                                         ifelse(age <= 40, 4, 5)))),
#          mu = NA,
#          sd = NA)

load('motnp_nodalregression.RData')
sim <- nodes %>% 
  mutate(age_cat_chr = ifelse(age_cat_chr == '40+', '40-60', age_cat_chr)) %>% 
  separate(age_cat_chr, into = c('min_age','max_age'),
           remove = F, sep = '-') %>% 
  mutate(min_age = as.numeric(min_age),
         max_age = as.numeric(max_age))
rm(list = ls()[! ls() %in% 'sim'])

sim$age_sim <- NA
for(i in 1:nrow(sim)){
  sd <- (sim$max_age[i] - sim$min_age[i])/2
  age <- round(rnorm(1000, mean = sim$age[i], sd = sd),0)
  age <- age[which(age < sim$max_age[i] & age > sim$min_age[i])[1]]
  sim$age_sim[i] <- age
}
rm(age, sd, i) ; gc()

sim <- sim %>% 
  select(-age, -min_age, -max_age) %>% 
  rename(age = age_sim) %>% 
  relocate(age, .after = node)

n_nodes <- nrow(sim)

## simulate centralities ####
## simulate age effect
sim_slope <- 0.1 #c(-1.5,-0.5,0.5,1.5,2.5)
sim_intcp <- -4
sim$mu <- sim$age * sim_slope + sim_intcp #sim$mu <- sim_slope[sim$age_cat] + sim_intcp       # simulate mean centrality on normal scale
plot(sim$mu ~ sim$age) #plot(sim$mu ~ sim$age_cat)

## simulate full distribution of samples per node
sim$sd <- 1            # allow to vary amongst nodes once confident that it's working with the simple version
sim_dat <- matrix(data = NA, nrow = 4000, ncol = n_nodes,
                  dimnames = list(NULL, sim$node_rank))    # create matrix
for(j in 1:n_nodes){
  sim_dat[,j] <- rnorm(n = nrow(sim_dat), mean = sim$mu[j], sd = sim$sd[j])  # simulate distribution
}
plot(sim_dat[1,] ~ sim$age)          # plot simulated values against age

## visualise -- centrality ~ age (years)
data.frame(sim_dat) %>% 
  pivot_longer(cols = everything(),
               names_to = "node", values_to = "centrality") %>% 
  separate(node, into = c('X','node_rank'), remove = T, sep = 1) %>% 
  dplyr::select(-X) %>% 
  mutate(node_rank = as.integer(node_rank)) %>% 
  left_join(sim[,c('node_rank','age')], by = 'node_rank') %>% 
  filter(node_rank %in% seq(1, n_nodes, by = 3)) %>% 
  ggplot(aes(x = centrality, fill = age)) +
  geom_density(linewidth = 0.4) +
  facet_grid(rows = vars(as.factor(node_rank)), scales = "free") +
  labs(x = "Eigenvector centrality (standardised)") + 
  theme_void() + 
  theme(strip.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

## visualise -- centrality ~ age (category)
data.frame(sim_dat) %>% 
  pivot_longer(cols = everything(),
               names_to = "node", values_to = "centrality") %>% 
  separate(node, into = c('X','node_rank'), remove = T, sep = 1) %>% 
  dplyr::select(-X) %>% 
  mutate(node_rank = as.integer(node_rank)) %>% 
  left_join(sim[,c('node_rank','age_cat_num')], by = 'node_rank') %>% 
  filter(node_rank %in% seq(1, n_nodes, by = 3)) %>% 
  ggplot(aes(x = centrality, fill = as.factor(age_cat_num))) +
  geom_density(linewidth = 0.4) +
  facet_grid(rows = vars(as.factor(node_rank)), scales = "free") +
  labs(x = "Eigenvector centrality (standardised)") + 
  scale_fill_viridis_d() +
  theme_void() + 
  theme(strip.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

## normal approximation ####
## plot covariance
plot(sim_dat[,which(sim$age == min(sim$age))[1]],
     sim_dat[,which(sim$age == max(sim$age))[1]])   # oldest and youngest

## calculate mean and covariance for all
sim_cent_mu <- apply(sim_dat, 2, mean)     # calculate means per node
sim_cent_cov <- cov(sim_dat)               # calculate covariance matrix

## check normal approximation
sim_cent_samples <- MASS::mvrnorm(1e5, sim_cent_mu, sim_cent_cov)     # simulate from multivariate normal
plot(density(sim_dat[, 1]), lwd = 2, las = 1,                     # plot true density curve
     main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
lines(density(sim_cent_samples[, 1]), col = rgb(0,0,1,0.5), lwd = 2)  # overlay normal approximation

## prior predictive check ####
n <- 100
beta_age <- rnorm(n, 0, 1)
intercept  <- rnorm(n, LaplacesDemon::logit(0.05), 2)
age_dirichlet <- rdirichlet(n, c(1,1,1,1,1))
x <- (min(sim$age_cat_num):max(sim$age_cat_num)) - (min(sim$age_cat_num) - 1)
plot(NULL, las = 1, xlab = 'age category', ylab = 'eigenvector (standardised)',
     ylim = c(min(sim_cent_mu)-5, max(sim_cent_mu)+5),
     xlim = c(min(x), max(x)))
abline(h = min(sim_cent_mu), lty = 2) ; abline(h = max(sim_cent_mu), lty = 2)
for(i in 1:n){
  y <- rep(NA, length(x))
  for(j in 1:length(x)){
    y[j] <- intercept[i] + beta_age[i]*sum(age_dirichlet[i,][1:x[j]])
  }
  lines(x = x, y = y, col = rgb(0,0,1,0.4))
}
rm(n, beta_age, intercept, age_dirichlet, x, y) ; gc()

## run model -- age as an ordered categorical variable with 1 window ####
## create data
n_age_cat <- length(unique(sim$age_cat_num))
eigen_list <- list(num_nodes = n_nodes,
                   num_age_cat = n_age_cat,
                   length_dirichlet = n_age_cat + 1,
                   centrality_mu = sim_cent_mu,
                   centrality_cov = sim_cent_cov,
                   node_age = sim$age_cat_fct,
                   prior_age = rep(1, n_age_cat))

## check inputs
plot(eigen_list$centrality_mu ~ eigen_list$node_age)

## load model
nodal_regression <- cmdstan_model('models/eigen_regression_motnp.stan')

## run model
n_chains <- 4
n_samples <- 1000
fit_sim <- nodal_regression$sample(data = eigen_list,
                                   chains = n_chains, parallel_chains = n_chains,
                                   iter_warmup = n_samples, iter_sampling = n_samples)

## check outputs ####
## view summary
(summary <- fit_sim$summary())
par(mfrow = c(3,1))
hist(summary$rhat, breaks = 50)
hist(summary$ess_bulk, breaks = 50)
hist(summary$ess_tail, breaks = 50)
par(mfrow = c(1,1))

## extract posterior
params <- fit_sim$draws(format = 'draws_df')
delta <- params %>% 
  select(`delta[1]`,`delta[2]`,`delta[3]`,`delta[4]`,`delta[5]`,
         `.chain`,`.iteration`,`.draw`)
delta_j <- params %>% 
  select(`delta_j[1]`,`delta_j[2]`,`delta_j[3]`,`delta_j[4]`,`delta_j[5]`,`delta_j[6]`,
         `.chain`,`.iteration`,`.draw`)

## traceplot linear effect size
parameters_to_check <- c('intercept','beta_age','sigma','predictor[1]','predictor[20]','predictor[50]')
params %>% 
  select(all_of(parameters_to_check), .chain, .iteration, .draw) %>%
  pivot_longer(cols = all_of(parameters_to_check), 
               names_to = 'parameter', values_to = 'value') %>%
  rename(chain_position = .iteration,
         chain = .chain,
         draw = .draw) %>% 
  #filter(chain == 4) %>% 
  ggplot(aes(x = chain_position, y = value, colour = as.factor(chain)))+ # beta_age has fully mixed, but more stretch on one side than the other
  geom_line()+
  facet_wrap(. ~ parameter, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'none')
delta %>% 
  pivot_longer(cols = c(`delta[1]`,`delta[2]`,`delta[3]`,`delta[4]`,`delta[5]`), 
               names_to = 'parameter', values_to = 'value') %>% 
  rename(chain_position = .iteration,
         chain = .chain,
         draw = .draw) %>% 
  #filter(chain == 4) %>% 
  ggplot(aes(x = chain_position, y = value, colour = as.factor(chain)))+
  geom_line()+
  facet_wrap(. ~ parameter, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'none')
delta_j %>% 
  pivot_longer(cols = c(`delta_j[1]`,`delta_j[2]`,`delta_j[3]`,`delta_j[4]`,`delta_j[5]`,`delta_j[6]`), 
               names_to = 'parameter', values_to = 'value') %>% 
  rename(chain_position = .iteration,
         chain = .chain,
         draw = .draw) %>% 
  #filter(chain == 4) %>% 
  ggplot(aes(x = chain_position, y = value, colour = as.factor(chain)))+
  geom_line()+
  facet_wrap(. ~ parameter, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'none')

## posterior predictive check ####
plot(density(sim_dat[1, ]), las = 1, ylim = c(0,0.4),
     main = "Posterior predictive check (standardised centrality):\nblack = data, blue = predicted",
     col=rgb(0, 0, 0, 0.25))
for (i in 1:100) {
  j <- sample(1:length(params$beta_age), 1)
  lines(density(sim_dat[j, ]), col=rgb(0, 0, 0, 0.25))
  mu <- rep(NA, length(eigen_list$node_age))
  for(k in 1:length(eigen_list$node_age)){
    mu[k] <- params$intercept[j] + params$beta_age[j]*sum(delta_j[j,(1:eigen_list$node_age[k])])
  }
  sigma <- sim_cent_cov + diag(rep(params$sigma[j], n_nodes))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}
warnings()

## predict from model ####
# create functions for predictions
predict_mean_centrality <- function(params, delta_j, pred_data){
  mu_matrix <- matrix(data = NA,
                      nrow = length(params$beta_age),
                      ncol = nrow(pred_data),
                      dimnames = list(1:length(params$beta_age),
                                      pred_data$node))
  for(i in 1:nrow(mu_matrix)){
    for(j in 1:ncol(mu_matrix)){
      mu_matrix[i,j] <- params$intercept[i] + params$beta_age[i] * sum(delta_j[i,(1:pred_data$age_cat_fct[j])])
    }
  }
  return(mu_matrix)
}
predict_full_centrality <- function(params, mu_matrix, cov_matrix){
  full_matrix <- mu_matrix
  for(i in 1:nrow(full_matrix)){
    full_matrix[i,] <- MASS::mvrnorm(n = 1, mu = mu_matrix[i,],
                                      Sigma = cov_matrix + diag(rep(params$sigma[i], ncol(full_matrix))))
  }
  return(full_matrix)
}

## create predictions data frame
pred_data <- sim[order(sim$age_cat_num),] # ordered because otherwise a later plotting function was being super annoying about forming logical lines and which bits it connected

## predicted means
predict_mu <- predict_mean_centrality(params = params, delta_j = delta_j, pred_data = pred_data)
pred_data$mu_mean_predict <- apply(predict_mu, 2, mean)
pred_data$mu_lwr_predict <- apply(predict_mu, 2, quantile, probs = 0.025)
pred_data$mu_upr_predict <- apply(predict_mu, 2, quantile, probs = 0.975)

## full predictions
predict_full <- predict_full_centrality(params = params, mu_matrix = predict_mu, cov_matrix = sim_cent_cov)
pred_data$full_lwr_predict <- apply(predict_full, 2, quantile, probs = 0.025)
pred_data$full_upr_predict <- apply(predict_full, 2, quantile, probs = 0.975)

## compare to raw data
ggplot(pred_data)+
  geom_ribbon(aes(x = age_cat_num, ymin = full_lwr_predict, ymax = full_upr_predict),
              alpha = 0.3, fill = 'purple')+
  geom_ribbon(aes(x = age_cat_num, ymin = mu_lwr_predict, ymax = mu_upr_predict),
              alpha = 0.3, fill = 'blue')+
  geom_point(aes(x = age_cat_num, y = mu))+
  geom_line(aes(x = age_cat_num, y = mu_mean_predict))+
  theme_classic()
ggplot(pred_data,
       aes(y = mu_mean_predict, x = mu))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)

## convert to invlogit scale
predict_mu_invlogit <- invlogit(predict_mu)
predict_full_invlogit <- invlogit(predict_full)
pred_data$mu_mean_predict_invlogit <- apply(predict_mu_invlogit, 2, mean)
pred_data$mu_lwr_predict_invlogit <- apply(predict_mu_invlogit, 2, quantile, probs = 0.025)
pred_data$mu_upr_predict_invlogit <- apply(predict_mu_invlogit, 2, quantile, probs = 0.975)
pred_data$full_lwr_predict_invlogit <- apply(predict_full_invlogit, 2, quantile, probs = 0.025)
pred_data$full_upr_predict_invlogit <- apply(predict_full_invlogit, 2, quantile, probs = 0.975)

## compare to invlogit raw data
pred_data$mu_raw_invlogit <- invlogit(pred_data$mu)
ggplot(pred_data)+
  geom_ribbon(aes(x = age_cat_num, ymin = full_lwr_predict_invlogit, ymax = full_upr_predict_invlogit),
              alpha = 0.3, fill = 'purple')+
  geom_ribbon(aes(x = age_cat_num, ymin = mu_lwr_predict_invlogit, ymax = mu_upr_predict_invlogit),
              alpha = 0.3, fill = 'blue')+
  geom_point(aes(x = age_cat_num, y = LaplacesDemon::invlogit(mu)))+
  geom_line(aes(x = age_cat_num, y = mu_mean_predict_invlogit))+
  theme_classic()

## clean up
rm(eigen_list, nodal_regression, predict_full_invlogit, predict_mu_invlogit,#predict_full,  predict_mu,
   sigma, sim_cent_samples, summary, i, j, k, parameters_to_check) ; gc()

## extract original values from output -- NOTE: CURRENT METHOD WILL NOT WORK FOR REAL DATA AS I DON'T HAVE THE AGE VALUES IN YEARS TO GO ABOUT ADDING 1 YEAR TO EVERY ELEPHANT. WILL NEED TO COME UP WITH AN ALTERNATIVE OR JUST SHOW CONTRASTS BETWEEN AGE CATEGORIES FOR REAL DATA, BUT HERE CHECK THE EFFECT ON REAL AGE SCALE TO BE SURE THAT THE MODEL IS WORKING. ####
## set up objects to store predictions
pred_means <- list()
pred_fulls <- list()

## predict for ages from 10 years younger than current up to 10 years older: multiple opportunities to pass through age category thresholds 
for( i in -10:10 ){
  pred_data_new <- sim[order(sim$age_cat_num),] %>% 
    mutate(age = age + i) %>% 
    mutate(age_cat_num = ifelse(age <= 15, 1,
                                ifelse(age <= 20, 2,
                                       ifelse(age <= 25, 3,
                                              ifelse(age <= 40, 4, 5))))) %>% 
    mutate(age_cat_fct = as.integer(as.factor(age_cat_num)))
  pred_mu_new <- predict_mean_centrality(params = params,
                                         delta_j = delta_j,
                                         pred_data = pred_data_new)
  pred_full_new <- predict_full_centrality(params = params,
                                           mu_matrix = pred_mu_new,
                                           cov_matrix = sim_cent_cov)
  pred_means[[(i+11)]] <- pred_mu_new
  pred_fulls[[(i+11)]] <- pred_full_new
  print(i)
}
warnings() # should just be "Dropping 'draws_df' class as required metadata was removed." over and over

## create data frame to store outputs
contrasts_all <- data.frame(older = rep(-10:10, each = 21),
                            younger = rep(-10:10, 21),
                            mean_contrast = NA)

## for each combination of predictions, calculate contrast and divide difference by number of years changed
contrast_matrices <- list()
for(i in 1:21){
  for(j in 1:21){
    if(i > j){
      contrast <- (pred_fulls[[i]] - pred_fulls[[j]]) / (i-j)
      row_to_fill <- which(contrasts_all$older == (i - 11) & 
                             contrasts_all$younger == (j - 11))
      contrast_matrices[[row_to_fill]] <- contrast
      contrasts_all$mean_contrast[row_to_fill] <- mean(contrast)
    }
  }
}

## remove impossible age combinations
contrasts_all <- contrasts_all %>%
  filter(is.na(mean_contrast) == FALSE)

## take only the times where the elephant changed age category
contrast <- pred_fulls[[2]]-pred_fulls[[1]]
mean(contrast)

## compare to input
sim_slope
mean(contrasts_all$mean_contrast)
quantile(contrasts_all$mean_contrast, prob = c(0.025, 0.975))

## save
save.image('step4_nodalregression/checks/simulation_motnp_nodal.RData')
dev.off()
