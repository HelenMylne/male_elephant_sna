## set up #### 
library(LaplacesDemon)
library(tidyverse)
library(cmdstanr)

#set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0/')

## simulate population ####
## define population parameters
min_age <- 11                                                # youngest individual
max_age <- 60                                                #  oldest  individual
n_nodes <- ((max_age+1) - min_age)                           # total nodes = 2 per age

## simulate population for first time window
sim <- data.frame(node = 1:n_nodes,                          # create data frame of individuals
                  age = c(seq(min_age, max_age, by = 2),     # some present right up to max
                          sample(min_age:max_age, n_nodes/2, prob = 1/(min_age:max_age), # more at lower ages
                                 replace = T)),
                  mu = NA, sd = NA)
sim$age_std <- ( sim$age - mean(sim$age) ) / sd(sim$age)

## redefine parameters
n_data <- nrow(sim)
n_nodes <- length(unique(sim$node))

## simulate centralities ####
## simulate age effect
sim_slope <- -1                     # set age effect -- smaller = bigger impact on invlogit scale as values large
sim_intcp <- 2
sim$true_mean <- sim$age_std * sim_slope + sim_intcp   # simulate mean centrality on normal scale
plot(true_mean ~ age, data = sim, col = 'red', pch = 19)            # plot
sim$true_mean_std <- ( sim$true_mean - mean(sim$true_mean) ) / sd(sim$true_mean)
plot(true_mean_std ~ age_std, data = sim, col = 'red', pch = 19)   # plot

## simulate full distribution of samples per node
sim$sd <- 0.01
sim_dat <- matrix(data = NA, nrow = 4000, ncol = n_data, dimnames = list(NULL, sim$node_window))    # create matrix
for(j in 1:n_data){
  sim_dat[,j] <- rnorm(n = nrow(sim_dat), mean = sim$true_mean[j], sd = sim$sd[j])  # simulate distribution
}
plot(sim_dat[1,] ~ sim$age, col = 'red', pch = 19)       # plot simulated values against age for window 1
lwr <- apply(sim_dat, 2, quantile, prob = 0.025)
upr <- apply(sim_dat, 2, quantile, prob = 0.975)
for(i in 1:ncol(sim_dat)){
  lines(x = c(sim$age[i], sim$age[i]), y = c(lwr[i],upr[i]))
}

## standardise
sim_dat_std <- sim_dat               # create matrix to fill
for(i in 1:nrow(sim_dat_std)){
  sim_dat_std[i,] <- (sim_dat[i,] - mean(sim_dat[i,]) ) / sd(sim_dat[i,]) # standardise values
}
plot(sim_dat_std[1,] ~ sim$age, col = 'red', pch = 19)       # plot simulated values against age for window 1

## visualise
df_wide <- data.frame(sim_dat_std)
df_plot <- df_wide %>% 
  pivot_longer(cols = everything(),
               names_to = "node", values_to = "centrality") %>% 
  separate(node, into = c('X','node'), remove = T, sep = 1) %>% 
  dplyr::select(-X) %>% 
  #separate(node_window, into = c('node','window'), sep = '_', remove = F) %>% 
  mutate(node = as.integer(node)) %>% 
  left_join(sim[,c('node','age')], by = 'node')
ggplot(data = df_plot, aes(x = centrality, fill = age, group = node)) +
  geom_density(linewidth = 0.4) +
  facet_grid(rows = vars(as.factor(node)), scales = "free") +
  labs(x = "Eigenvector centrality (standardised)") + 
  theme_void() + 
  theme(strip.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

## normal approximation ####
# compute normal approximation by window -- calculate means per node
sim_cent_mu_1 <- apply(sim_dat_std, 2, mean)

# compute normal approximation by window -- calculate covariance matrix
sim_cent_cov_1 <- cov(sim_dat_std)

## check normal approximation -- simulate from combined mean and covariance, plot curve against relative node ID 
sim_cent_samples <- MASS::mvrnorm(1e5, sim_cent_mu_1, sim_cent_cov_1)     # simulate from multivariate normal
plot(density(sim_dat_std[,1]), lwd = 2, las = 1,                         # plot true density curve
     main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
lines(density(sim_cent_samples[, 1]), col = rgb(0,0,1,0.5), lwd = 2)      # overlay normal approximation

## prep inputs ####
## create data
eigen_list <- list(num_data = n_data,
                   num_nodes = n_nodes,
                   #num_windows = length(unique(sim$window)),
                   num_nodes_window1 = length(which(sim$window == 1)),
                   centrality_mu_1 = sim_cent_mu_1,
                   centrality_cov_1 = sim_cent_cov_1,
                   node_age = sim$age_std,
                   #window = sim$window,
                   #nodes = sim$node,
                   nodes_window1 = sim$node[sim$window == 1])

## check inputs
plot(eigen_list$centrality_mu_1 ~ eigen_list$node_age, pch = 19, col = 'red')

## prior predictive check ####
n <- 100
beta_age <- rnorm(n, -1, 0.1)
intercept  <- rnorm(n, 2, 0.1)
plot(NULL, las = 1, xlab = 'age (standardised)', ylab = 'eigenvector (standardised)',
     ylim = c(min(sim_cent_mu_1)-2, max(sim_cent_mu_1)+2),
     xlim = c(min(sim$age_std), max(sim$age_std)))
abline(h = min(sim_cent_mu_1), lty = 2)
abline(h = max(sim_cent_mu_1), lty = 2)
for(i in 1:n){
  lines(x = seq(min(sim$age_std), max(sim$age_std), length.out = 2),
        y = ( intercept[i] + beta_age[i] * c(min(sim$age_std),max(sim$age_std)) )# * sd(sim_dat) + mean(sim_dat), # this is predicting on unstandardised scale, not standardised
        col = rgb(0,0,1,0.4))
} # looks mad when I set true slope very small, but decent when larger -- given that I have no idea what the real version should be, I've kept them relatively weak

## run model -- age as a continuous variable with multiple windows ####
## load model
#nodal_regression <- stan_model('models/eigen_regression_intercept.stan')
nodal_regression <- cmdstan_model('models/eigen_regression_singlewindow.stan')

## run model
n_chains <- 4
n_samples <- 1000
#fit_sim <- sampling(nodal_regression, data = eigen_list, chains = n_chains, cores = n_chains)
fit_sim <- nodal_regression$sample(data = eigen_list,
                                   chains = n_chains, parallel_chains = n_chains,
                                   iter_warmup = n_samples, iter_sampling = n_samples)

## check outputs ####
## extract model fit -- very good!
fit_sim$summary()
summary <- fit_sim$summary()
par(mfrow = c(3,1))
hist(summary$rhat, breaks = 50)
hist(summary$ess_bulk, breaks = 50)
hist(summary$ess_tail, breaks = 50)
par(mfrow = c(1,1))

## extract posterior
#params <- rstan::extract(fit_sim)
params <- fit_sim$draws(format = 'draws_df')
colnames(params)
rand_window <- params %>% 
  dplyr::select(`rand_window[1]`,`rand_window[2]`,`rand_window[3]`)
rand_node <- params[,8:(n_nodes+7)]

## traceplot linear effect size
#traceplot(fit_sim, pars = c('intercept','beta_age','sigma','predictor[1]','predictor[50]','predictor[100]'))
params %>% 
  select(intercept,beta_age,sigma,
         `rand_node[1]`,`rand_node[50]`,`rand_node[100]`,
         `rand_window[1]`,`rand_window[2]`,
         `predictor_window1[1]`,`predictor_window1[16]`,`predictor_window1[38]`,`predictor_window1[50]`,
         `predictor_window2[1]`,`predictor_window2[25]`,`predictor_window2[50]`,`predictor_window2[75]`) %>% 
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'draw') %>% 
  mutate(chain_position = rep(rep(1:n_samples, each = 16,
  ), n_chains),
  chain = rep(1:n_chains, each = 16*n_samples)) %>% 
  #filter(chain == 4) %>% # inspect individual chains -- some are really bad and haven't explored at all or are very wandery
  ggplot(aes(x = chain_position, y = draw, colour = as.factor(chain)))+
  geom_line()+
  facet_wrap(. ~ parameter, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'none')

## posterior predictive check
par(mfrow = c(3,1))
plot(density(sim_dat_std[1, which(sim$window == 1)]), las = 1, ylim = c(0,0.6),
     main = "Posterior predictive check (standardised centrality):\nblack = data, blue = predicted",
     col=rgb(0, 0, 0, 0.25))
eigen_data1 <- list(num_nodes_window1 = length(which(sim$window == 1)),
                    centrality_mu_1 = sim_cent_mu_1,
                    centrality_cov_1 = sim_cent_cov_1,
                    node_age = sim$age_std[sim$window == 1],
                    window = 1,
                    nodes = sim$node[sim$window == 1],
                    nodes_window1 = sim$node[sim$window == 1])
for (i in 1:100) {
  j <- sample(1:length(params$beta_age), 1)
  lines(density(sim_dat_std[j, which(sim$window == 1)]), col=rgb(0, 0, 0, 0.25))
  mu <- params$beta_age[j]*eigen_data1$node_age + params$intercept[j]
  for(k in 1:length(mu)) {
    mu[k] <- mu[k] + as.numeric(rand_window[j,eigen_data1$window]) + as.numeric(rand_node[j,eigen_data1$nodes[k]])
  }
  sigma <- sim_cent_cov_1 + diag(rep(params$sigma[j], eigen_list$num_nodes_window1))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}


plot(density(sim_dat_std[1, which(sim$window == 2)]), las = 1, ylim = c(0,0.6),
     main = "Posterior predictive check (standardised centrality):\nblack = data, blue = predicted",
     col=rgb(0, 0, 0, 0.25))
eigen_data2 <- list(num_nodes_window2 = length(which(sim$window == 2)),
                    centrality_mu_2 = sim_cent_mu_2,
                    centrality_cov_2 = sim_cent_cov_2,
                    node_age = sim$age_std[sim$window == 2],
                    window = 2,
                    nodes = sim$node[sim$window == 2],
                    nodes_window2 = sim$node[sim$window == 2])
for (i in 1:100) {
  j <- sample(1:length(params$beta_age), 1)
  lines(density(sim_dat_std[j, which(sim$window == 2)]), col=rgb(0, 0, 0, 0.25))
  mu <- params$beta_age[j]*eigen_data2$node_age + params$intercept[j]
  for(k in 1:length(mu)) {
    mu[k] <- mu[k] + as.numeric(rand_window[j,eigen_data2$window]) + as.numeric(rand_node[j,eigen_data2$nodes[k]])
  }
  sigma <- sim_cent_cov_2 + diag(rep(params$sigma[j], eigen_list$num_nodes_window2))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}

plot(density(sim_dat_std[1, which(sim$window == 3)]), las = 1, ylim = c(0,0.6),
     main = "Posterior predictive check (standardised centrality):\nblack = data, blue = predicted",
     col=rgb(0, 0, 0, 0.25))
eigen_data3 <- list(num_nodes_window3 = length(which(sim$window == 3)),
                    centrality_mu_3 = sim_cent_mu_3,
                    centrality_cov_3 = sim_cent_cov_3,
                    node_age = sim$age_std[sim$window == 3],
                    window = 3,
                    nodes = sim$node[sim$window == 3],
                    nodes_window3 = sim$node[sim$window == 3])
for (i in 1:100) {
  j <- sample(1:length(params$beta_age), 1)
  lines(density(sim_dat_std[j, which(sim$window == 3)]), col=rgb(0, 0, 0, 0.25))
  mu <- params$beta_age[j]*eigen_data3$node_age + params$intercept[j]
  for(k in 1:length(mu)) {
    mu[k] <- mu[k] + as.numeric(rand_window[j,eigen_data3$window]) + as.numeric(rand_node[j,eigen_data3$nodes[k]])
  }
  sigma <- sim_cent_cov_3 + diag(rep(params$sigma[j], eigen_list$num_nodes_window3))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}
par(mfrow = c(1,1))

## predict from model -- predict from raw data, not from hypothetical/counterfactual data ####
mu <- matrix(NA, nrow = nrow(params), ncol = nrow(sim), dimnames = list(NULL, sim$node_window))
for(i in 1:nrow(mu)){
  mu[i,] <- params$beta_age[i] * sim$age_std + params$intercept[i] + as.numeric(rand_window[i,sim$window]) + as.numeric(rand_node[i, sim$node])
}
sim$mu_mean <- apply(mu, 2, mean)

sigma1 <- array(NA, dim = c(eigen_list$num_nodes_window1, eigen_list$num_nodes_window1, nrow(params)),
                dimnames = list(eigen_list$nodes_window1, eigen_list$nodes_window1, NULL))
sigma2 <- array(NA, dim = c(eigen_list$num_nodes_window2, eigen_list$num_nodes_window2, nrow(params)),
                dimnames = list(eigen_list$nodes_window2, eigen_list$nodes_window2, NULL))
sigma3 <- array(NA, dim = c(eigen_list$num_nodes_window3, eigen_list$num_nodes_window3, nrow(params)),
                dimnames = list(eigen_list$nodes_window3, eigen_list$nodes_window3, NULL))
for(i in 1:nrow(params)){
  sigma1[,,i] <- sim_cent_cov_1 + diag(rep(params$sigma[i], eigen_list$num_nodes_window1))
  sigma2[,,i] <- sim_cent_cov_2 + diag(rep(params$sigma[i], eigen_list$num_nodes_window2))
  sigma3[,,i] <- sim_cent_cov_3 + diag(rep(params$sigma[i], eigen_list$num_nodes_window3))
}

predictions <- matrix(NA, nrow = nrow(params), ncol = nrow(sim), dimnames = list(NULL, sim$node_window))
for(i in 1:nrow(predictions)){
  predictions[i,sim$window == 1] <- MASS::mvrnorm(1, mu[i,sim$window == 1], sigma1[,,i])
  predictions[i,sim$window == 2] <- MASS::mvrnorm(1, mu[i,sim$window == 2], sigma2[,,i])
  predictions[i,sim$window == 3] <- MASS::mvrnorm(1, mu[i,sim$window == 3], sigma3[,,i])
}

pred_data <- sim
pred_data$mu_model <- NA
pred_data$prediction <- NA
for(i in 1:nrow(pred_data)){
  pred_data$mu_model[i] <- list(mu[,i])
  pred_data$prediction[i] <- list(predictions[,i])
}
pred_data <- unnest(cols = prediction, data = pred_data)

sim$mu_lwr <- NA ; sim$mu_upr <- NA
sim$predict_lwr <- NA ; sim$predict_upr <- NA
for( i in 1:nrow(sim) ){
  sim$mu_lwr[i] <- rethinking::HPDI(mu[,i], prob = 0.95)[1]
  sim$mu_upr[i] <- rethinking::HPDI(mu[,i], prob = 0.95)[2]
  sim$predict_lwr[i] <- rethinking::HPDI(predictions[,i], prob = 0.95)[1]
  sim$predict_upr[i] <- rethinking::HPDI(predictions[,i], prob = 0.95)[2]
}

## plot mean of model vs mean of raw data
plot(sim$mu_mean ~ sim$mu_std, pch = 19, col = rgb(0,0,1,0.2)) ; abline(a = 0, b = 1, lwd = 2)

## plot all predictions with points = mean of raw data
ggplot(sim)+
  geom_ribbon(aes(x = age, ymin = predict_lwr, ymax = predict_upr, fill = as.factor(window)),
              alpha = 0.2)+
  geom_ribbon(aes(x = age, ymin = mu_lwr, ymax = mu_upr, fill = as.factor(window)),
              alpha = 0.4)+
  geom_line(aes(x = age, y = mu_mean, colour = as.factor(window)))+
  geom_point(aes(x = age, y = mu_std, colour = as.factor(window)))+
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  facet_wrap(. ~ as.factor(window))+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'age (years)')

## plot on unstandardised scale ####
## revert predictions to unstandardised scale
mu_ustd <- mu * sd(sim$mu) + mean(sim$mu)
predictions_ustd <- predictions * sd(sim$mu) + mean(sim$mu)

sim$mu_mean_ustd <- apply(mu_ustd, 2, mean)
sim$mu_lwr_ustd <- NA ; sim$mu_upr_ustd <- NA
sim$predict_lwr_ustd <- NA ; sim$predict_upr_ustd <- NA
for( i in 1:nrow(sim) ){
  sim$mu_lwr_ustd[i] <- rethinking::HPDI(mu_ustd[,i], prob = 0.95)[1]
  sim$mu_upr_ustd[i] <- rethinking::HPDI(mu_ustd[,i], prob = 0.95)[2]
  sim$predict_lwr_ustd[i] <- rethinking::HPDI(predictions_ustd[,i], prob = 0.95)[1]
  sim$predict_upr_ustd[i] <- rethinking::HPDI(predictions_ustd[,i], prob = 0.95)[2]
}

ggplot(sim)+
  geom_ribbon(aes(x = age, ymin = predict_lwr_ustd, ymax = predict_upr_ustd, fill = as.factor(window)),
              alpha = 0.2)+
  geom_ribbon(aes(x = age, ymin = mu_lwr_ustd, ymax = mu_upr_ustd, fill = as.factor(window)),
              alpha = 0.4)+
  geom_line(aes(x = age, y = mu_mean_ustd, colour = as.factor(window)))+
  geom_point(aes(x = age, y = mu, colour = as.factor(window)))+
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  facet_wrap(. ~ as.factor(window))+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'standardised age')

## extract original values from output -- I CAN'T WORK OUT WHY CALCULATING THE MARGNAL EFFECTS FROM THE PREDICTIONS WOULD BE ANY DIFFERENT TO THOSE PRESENTED IN THE SUMMARY, GIVEN THAT THERE IS ONLY ONE PREDICTOR VARIABLE. ALSO, IF I'M MAKING PREDICTIONS BASED ON THE RAW DATA RATHER THAN A HYPOTHETICAL/COUNTERFACTUAL DATASET, HOW DO YOU CALCULATE THE EQUIVALENT OF THE DO OPERATOR? WOULD I BASICALLY JUST ADD 1 YEAR TO THE AGE OF ALL INDIVIDUALS AND REDO THE PREDICTIONS, THEN LOOK AT THE DIFFERENCE? ####

## add additional time windows ####
## simulate population for first time window
sim1 <- sim

## simulate population for second time window -- half the nodes resampled from window 1, half new
sim2 <- data.frame(node = c(sample(1:n_nodes, n_nodes/2, replace = F),
                            (n_nodes+1):(n_nodes*2)),                          # create data frame of individuals
                   age = NA,
                   mu = NA, sd = NA)
for(i in 1:nrow(sim2)){
  if(sim2$node[i] %in% sim1$node){
    sim2$age[i] <- sim1$age[sim1$node == sim2$node[i]] + 2 # if present in 1st window, age = age in t1 + 2 y
  } else {
    sim2$age[i] <- sample(min_age:max_age, 1, prob = 1/(min_age:max_age)) # if not present, randomly assign a new age, more likely to be young than old
  }
}

## simulate population for third time window -- as above, half the nodes are resampled from oervious two time windows, so some may be represented in all three windows, some twice, and some only once
n_nodes <- length(unique(sim2$node))
sim3 <- data.frame(node = c(sample(1:n_nodes, n_nodes/2, replace = F),
                            (n_nodes+1):(n_nodes*2)),                          # create data frame of individuals
                   age = NA,
                   mu = NA, sd = NA)
for(i in 1:nrow(sim3)){
  if(sim3$node[i] %in% sim2$node){
    sim3$age[i] <- sim2$age[sim2$node == sim3$node[i]] + 2
  } else {
    if(sim3$node[i] %in% sim1$node){
      sim3$age[i] <- sim1$age[sim1$node == sim3$node[i]] + 2
    } else {
      sim3$age[i] <- sample(min_age:max_age, 1, prob = 1/(min_age:max_age)) #(min_age:max_age)[i-(n_nodes/2)]
    }
  }
}

## combine populations into single data frame
sim1$window <- 1 ; sim2$window <- 2 ; sim3$window <- 3 # create window ID variable
sim <- rbind(sim1, sim2, sim3)
sim$node_window <- paste0(sim$node, '_', sim$window)
sim$age_std <- ( sim$age - mean(sim$age) ) / sd(sim$age)
rm(sim1, sim2, sim3) ; gc()

## redefine parameters
n_data <- nrow(sim)
n_nodes <- length(unique(sim$node))
n_windows <- length(unique(sim$window))

## view ages
hist(sim$age)

## simulate centralities ####
## simulate age effect
sim_slope <- -1                     # set age effect -- smaller = bigger impact on invlogit scale as values large
sim_intcp <- 2
sim_window_unique <- rnorm(n_windows, mean = 0, sd = 1)
sim_node_unique <- rnorm(n_nodes, mean = 0, sd = 0.2)
sim$mu <- sim$age_std * sim_slope + sim_intcp + sim_window_unique[sim$window] + sim_node_unique[sim$node]    # simulate mean centrality on normal scale
plot(mu ~ age, data = sim[sim$window == 1,], col = 'red', pch = 19, ylim = c(0,5))            # plot
points(mu ~ age, data = sim[sim$window == 2,], col = 'blue', pch = 19)                        # plot
points(mu ~ age, data = sim[sim$window == 3,], col = 'green', pch = 19)                       # plot
sim$mu_std <- ( sim$mu - mean(sim$mu) ) / sd(sim$mu)
plot(mu_std ~ age_std, data = sim[sim$window == 1,], col = 'red', pch = 19, ylim = c(-5,5))   # plot
points(mu_std ~ age_std, data = sim[sim$window == 2,], col = 'blue', pch = 19)                # plot
points(mu_std ~ age_std, data = sim[sim$window == 3,], col = 'green', pch = 19)               # plot

## simulate full distribution of samples per node
sim$sd <- 1
sim_dat <- matrix(data = NA, nrow = 4000, ncol = n_data, dimnames = list(NULL, sim$node_window))    # create matrix
for(j in 1:n_data){
  sim_dat[,j] <- rnorm(n = nrow(sim_dat), mean = sim$mu[j], sd = sim$sd[j])  # simulate distribution
}
plot(sim_dat[1,1:(length(which(sim$window == 1)))] ~ sim$age[1:(length(which(sim$window == 1)))],
     col = 'red', pch = 19, ylim = c(-1,6))       # plot simulated values against age for window 1
points(sim_dat[1,(length(which(sim$window == 1))+1):(length(which(sim$window == 1 | sim$window == 2)))] ~ sim$age[(length(which(sim$window == 1))+1):(length(which(sim$window == 1 | sim$window == 2)))],
       col = 'blue', pch = 19)    # plot simulated values against age for window 2
points(sim_dat[1,(length(which(sim$window == 1 | sim$window == 2))+1):nrow(sim)] ~ sim$age[(length(which(sim$window == 1 | sim$window == 2))+1):nrow(sim)],
       col = 'green', pch = 19)   # plot simulated values against age for window 3

## standardise
sim_dat_std <- sim_dat               # create matrix to fill
for(i in 1:nrow(sim_dat_std)){
  sim_dat_std[i,] <- (sim_dat[i,] - mean(sim_dat[i,]) ) / sd(sim_dat[i,]) # standardise values
}
plot(sim_dat_std[1,which(sim$window == 1)] ~ sim$age[which(sim$window == 1)],
     col = 'red', pch = 19, ylim = c(-5,5))       # plot simulated values against age for window 1
points(sim_dat_std[1,which(sim$window == 2)] ~ sim$age[which(sim$window == 2)],
       col = 'blue', pch = 19)    # plot simulated values against age for window 2
points(sim_dat_std[1,which(sim$window == 3)] ~ sim$age[which(sim$window == 3)],
       col = 'green', pch = 19)   # plot simulated values against age for window 3

## visualise
df_wide <- data.frame(sim_dat_std)
df_plot <- df_wide %>% 
  pivot_longer(cols = everything(),
               names_to = "node", values_to = "centrality") %>% 
  separate(node, into = c('X','node_window'), remove = T, sep = 1) %>% 
  dplyr::select(-X) %>% 
  separate(node_window, into = c('node','window'), sep = '_', remove = F) %>% 
  mutate(node = as.integer(node)) %>% 
  left_join(sim[,c('node_window','age')], by = 'node_window') %>% 
  filter(node %in% seq(1, n_data, by = 5))
ggplot(data = df_plot, aes(x = centrality, fill = age, group = node_window)) +
  geom_density(linewidth = 0.4, alpha = 0.6) +
  facet_grid(rows = vars(as.factor(node)), scales = "free") +
  labs(x = "Eigenvector centrality (standardised)") + 
  theme_void() + 
  theme(strip.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

## normal approximation ####
# compute normal approximation by window -- calculate means per node
sim_cent_mu_1 <- apply(sim_dat_std[,which(sim$window == 1)], 2, mean)
sim_cent_mu_2 <- apply(sim_dat_std[,which(sim$window == 2)], 2, mean)
sim_cent_mu_3 <- apply(sim_dat_std[,which(sim$window == 3)], 2, mean)

# compute normal approximation by window -- calculate covariance matrix
sim_cent_cov_1 <- cov(sim_dat_std[,which(sim$window == 1)])
sim_cent_cov_2 <- cov(sim_dat_std[,which(sim$window == 2)])
sim_cent_cov_3 <- cov(sim_dat_std[,which(sim$window == 3)])

## check normal approximation -- simulate from combined mean and covariance, plot curve against relative node ID 
sim_cent_samples <- MASS::mvrnorm(1e5, sim_cent_mu_1, sim_cent_cov_1)     # simulate from multivariate normal
node_id_sample <- sample(which(sim$window == 1),1)
plot(density(sim_dat_std[,node_id_sample]), lwd = 2, las = 1,                         # plot true density curve
     main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
lines(density(sim_cent_samples[, node_id_sample]), col = rgb(0,0,1,0.5), lwd = 2)      # overlay normal approximation

sim_cent_samples <- MASS::mvrnorm(1e5, sim_cent_mu_2, sim_cent_cov_2)     # simulate from multivariate normal
node_id_sample <- sample(which(sim$window == 2),1)
plot(density(sim_dat_std[,node_id_sample]),          # plot true density curve
     lwd = 2, las = 1,
     main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
lines(density(sim_cent_samples[, node_id_sample - length(which(sim$window == 1))]), col = rgb(0,0,1,0.5), lwd = 2)      # overlay normal approximation

sim_cent_samples <- MASS::mvrnorm(1e5, sim_cent_mu_3, sim_cent_cov_3)     # simulate from multivariate normal
node_id_sample <- sample(which(sim$window == 3),1)
plot(density(sim_dat_std[,node_id_sample]),          # plot true density curve
     lwd = 2, las = 1,
     main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
lines(density(sim_cent_samples[, node_id_sample - length(which(sim$window < 3))]), col = rgb(0,0,1,0.5), lwd = 2)      # overlay normal approximation

## prep inputs ####
## create data
eigen_list <- list(num_data = n_data,
                   num_nodes = n_nodes,
                   num_windows = length(unique(sim$window)),
                   num_nodes_window1 = length(which(sim$window == 1)),
                   num_nodes_window2 = length(which(sim$window == 2)),
                   num_nodes_window3 = length(which(sim$window == 3)),
                   centrality_mu_1 = sim_cent_mu_1,
                   centrality_mu_2 = sim_cent_mu_2,
                   centrality_mu_3 = sim_cent_mu_3,
                   centrality_cov_1 = sim_cent_cov_1,
                   centrality_cov_2 = sim_cent_cov_2,
                   centrality_cov_3 = sim_cent_cov_3,
                   node_age = sim$age_std,
                   window = sim$window,
                   nodes = sim$node,
                   nodes_window1 = sim$node[sim$window == 1],
                   nodes_window2 = sim$node[sim$window == 2],
                   nodes_window3 = sim$node[sim$window == 3])

## check inputs
plot(eigen_list$centrality_mu_1 ~ eigen_list$node_age[eigen_list$window == 1], pch = 19, col = 'red')
points(eigen_list$centrality_mu_2 ~ eigen_list$node_age[eigen_list$window == 2], pch = 19, col = 'blue')
points(eigen_list$centrality_mu_3 ~ eigen_list$node_age[eigen_list$window == 3], pch = 19, col = 'green')

## prior predictive check ####
n <- 100
beta_age <- rnorm(n, 0, 0.8)
intercept  <- rnorm(n, 0, 0.8)
plot(NULL, las = 1, xlab = 'age (standardised)', ylab = 'eigenvector (standardised)',
     ylim = c(min(c(sim_cent_mu_1, sim_cent_mu_2))-2, max(c(sim_cent_mu_1, sim_cent_mu_2))+2),
     xlim = c(min(sim$age_std), max(sim$age_std)))
abline(h = min(c(sim_cent_mu_1, sim_cent_mu_2, sim_cent_mu_3)), lty = 2)
abline(h = max(c(sim_cent_mu_1, sim_cent_mu_2, sim_cent_mu_3)), lty = 2)
for(i in 1:n){
  lines(x = seq(min(sim$age_std), max(sim$age_std), length.out = 2),
        y = intercept[i] + beta_age[i]*c(min(sim$age_std), max(sim$age_std)),
        col = rgb(0,0,1,0.4))
} # looks mad when I set true slope very small, but decent when larger -- given that I have no idea what the real version should be, I've kept them relatively weak

## run model -- age as a continuous variable with multiple windows ####
## load model
#nodal_regression <- stan_model('models/eigen_regression_intercept.stan')
nodal_regression <- cmdstan_model('models/eigen_regression_combinewindows.stan')

## run model
n_chains <- 4
n_samples <- 1000
#fit_sim <- sampling(nodal_regression, data = eigen_list, chains = n_chains, cores = n_chains)
fit_sim <- nodal_regression$sample(data = eigen_list,
                                   chains = n_chains, parallel_chains = n_chains,
                                   iter_warmup = n_samples, iter_sampling = n_samples)

## check outputs ####
## extract model fit -- very good!
fit_sim$summary()
summary <- fit_sim$summary()
par(mfrow = c(3,1))
hist(summary$rhat, breaks = 50)
hist(summary$ess_bulk, breaks = 50)
hist(summary$ess_tail, breaks = 50)
par(mfrow = c(1,1))

## extract posterior
#params <- rstan::extract(fit_sim)
params <- fit_sim$draws(format = 'draws_df')
colnames(params)
rand_window <- params %>% 
  dplyr::select(`rand_window[1]`,`rand_window[2]`,`rand_window[3]`)
rand_node <- params[,8:(n_nodes+7)]

## traceplot linear effect size
#traceplot(fit_sim, pars = c('intercept','beta_age','sigma','predictor[1]','predictor[50]','predictor[100]'))
params %>% 
  select(intercept,beta_age,sigma,
         `rand_node[1]`,`rand_node[50]`,`rand_node[100]`,
         `rand_window[1]`,`rand_window[2]`,
         `predictor_window1[1]`,`predictor_window1[16]`,`predictor_window1[38]`,`predictor_window1[50]`,
         `predictor_window2[1]`,`predictor_window2[25]`,`predictor_window2[50]`,`predictor_window2[75]`) %>% 
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'draw') %>% 
  mutate(chain_position = rep(rep(1:n_samples, each = 16,
  ), n_chains),
  chain = rep(1:n_chains, each = 16*n_samples)) %>% 
  #filter(chain == 4) %>% # inspect individual chains -- some are really bad and haven't explored at all or are very wandery
  ggplot(aes(x = chain_position, y = draw, colour = as.factor(chain)))+
  geom_line()+
  facet_wrap(. ~ parameter, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'none')

## posterior predictive check
par(mfrow = c(3,1))
plot(density(sim_dat_std[1, which(sim$window == 1)]), las = 1, ylim = c(0,0.6),
     main = "Posterior predictive check (standardised centrality):\nblack = data, blue = predicted",
     col=rgb(0, 0, 0, 0.25))
eigen_data1 <- list(num_nodes_window1 = length(which(sim$window == 1)),
                    centrality_mu_1 = sim_cent_mu_1,
                    centrality_cov_1 = sim_cent_cov_1,
                    node_age = sim$age_std[sim$window == 1],
                    window = 1,
                    nodes = sim$node[sim$window == 1],
                    nodes_window1 = sim$node[sim$window == 1])
for (i in 1:100) {
  j <- sample(1:length(params$beta_age), 1)
  lines(density(sim_dat_std[j, ]), col=rgb(0, 0, 0, 0.25))
  mu <- params$beta_age[j]*eigen_data1$node_age + params$intercept[j]
  for(k in 1:length(mu)) {
    mu[k] <- mu[k] + as.numeric(rand_window[j,eigen_data1$window]) + as.numeric(rand_node[j,eigen_data1$nodes[k]])
  }
  sigma <- sim_cent_cov_1 + diag(rep(params$sigma[j], eigen_list$num_nodes_window1))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}


plot(density(sim_dat_std[1, which(sim$window == 2)]), las = 1, ylim = c(0,0.6),
     main = "Posterior predictive check (standardised centrality):\nblack = data, blue = predicted",
     col=rgb(0, 0, 0, 0.25))
eigen_data2 <- list(num_nodes_window2 = length(which(sim$window == 2)),
                    centrality_mu_2 = sim_cent_mu_2,
                    centrality_cov_2 = sim_cent_cov_2,
                    node_age = sim$age_std[sim$window == 2],
                    window = 2,
                    nodes = sim$node[sim$window == 2],
                    nodes_window2 = sim$node[sim$window == 2])
for (i in 1:100) {
  j <- sample(1:length(params$beta_age), 1)
  lines(density(sim_dat_std[j, ]), col=rgb(0, 0, 0, 0.25))
  mu <- params$beta_age[j]*eigen_data2$node_age + params$intercept[j]
  for(k in 1:length(mu)) {
    mu[k] <- mu[k] + as.numeric(rand_window[j,eigen_data2$window]) + as.numeric(rand_node[j,eigen_data2$nodes[k]])
  }
  sigma <- sim_cent_cov_2 + diag(rep(params$sigma[j], eigen_list$num_nodes_window2))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}

plot(density(sim_dat_std[1, which(sim$window == 3)]), las = 1, ylim = c(0,0.6),
     main = "Posterior predictive check (standardised centrality):\nblack = data, blue = predicted",
     col=rgb(0, 0, 0, 0.25))
eigen_data3 <- list(num_nodes_window3 = length(which(sim$window == 3)),
                    centrality_mu_3 = sim_cent_mu_3,
                    centrality_cov_3 = sim_cent_cov_3,
                    node_age = sim$age_std[sim$window == 3],
                    window = 3,
                    nodes = sim$node[sim$window == 3],
                    nodes_window3 = sim$node[sim$window == 3])
for (i in 1:100) {
  j <- sample(1:length(params$beta_age), 1)
  lines(density(sim_dat_std[j, ]), col=rgb(0, 0, 0, 0.25))
  mu <- params$beta_age[j]*eigen_data3$node_age + params$intercept[j]
  for(k in 1:length(mu)) {
    mu[k] <- mu[k] + as.numeric(rand_window[j,eigen_data3$window]) + as.numeric(rand_node[j,eigen_data3$nodes[k]])
  }
  sigma <- sim_cent_cov_3 + diag(rep(params$sigma[j], eigen_list$num_nodes_window3))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}
par(mfrow = c(1,1))

## predict from model -- predict from raw data, not from hypothetical/counterfactual data ####
mu <- matrix(NA, nrow = nrow(params), ncol = nrow(sim), dimnames = list(NULL, sim$node_window))
for(i in 1:nrow(mu)){
  mu[i,] <- params$beta_age[i] * sim$age_std + params$intercept[i] + as.numeric(rand_window[i,sim$window]) + as.numeric(rand_node[i, sim$node])
}
sim$mu_mean <- apply(mu, 2, mean)

sigma1 <- array(NA, dim = c(eigen_list$num_nodes_window1, eigen_list$num_nodes_window1, nrow(params)),
                dimnames = list(eigen_list$nodes_window1, eigen_list$nodes_window1, NULL))
sigma2 <- array(NA, dim = c(eigen_list$num_nodes_window2, eigen_list$num_nodes_window2, nrow(params)),
                dimnames = list(eigen_list$nodes_window2, eigen_list$nodes_window2, NULL))
sigma3 <- array(NA, dim = c(eigen_list$num_nodes_window3, eigen_list$num_nodes_window3, nrow(params)),
                dimnames = list(eigen_list$nodes_window3, eigen_list$nodes_window3, NULL))
for(i in 1:nrow(params)){
  sigma1[,,i] <- sim_cent_cov_1 + diag(rep(params$sigma[i], eigen_list$num_nodes_window1))
  sigma2[,,i] <- sim_cent_cov_2 + diag(rep(params$sigma[i], eigen_list$num_nodes_window2))
  sigma3[,,i] <- sim_cent_cov_3 + diag(rep(params$sigma[i], eigen_list$num_nodes_window3))
}

predictions <- matrix(NA, nrow = nrow(params), ncol = nrow(sim), dimnames = list(NULL, sim$node_window))
for(i in 1:nrow(predictions)){
  predictions[i,sim$window == 1] <- MASS::mvrnorm(1, mu[i,sim$window == 1], sigma1[,,i])
  predictions[i,sim$window == 2] <- MASS::mvrnorm(1, mu[i,sim$window == 2], sigma2[,,i])
  predictions[i,sim$window == 3] <- MASS::mvrnorm(1, mu[i,sim$window == 3], sigma3[,,i])
}

pred_data <- sim
pred_data$mu_model <- NA
pred_data$prediction <- NA
for(i in 1:nrow(pred_data)){
  pred_data$mu_model[i] <- list(mu[,i])
  pred_data$prediction[i] <- list(predictions[,i])
}
pred_data <- unnest(cols = prediction, data = pred_data)

sim$mu_lwr <- NA ; sim$mu_upr <- NA
sim$predict_lwr <- NA ; sim$predict_upr <- NA
for( i in 1:nrow(sim) ){
  sim$mu_lwr[i] <- rethinking::HPDI(mu[,i], prob = 0.95)[1]
  sim$mu_upr[i] <- rethinking::HPDI(mu[,i], prob = 0.95)[2]
  sim$predict_lwr[i] <- rethinking::HPDI(predictions[,i], prob = 0.95)[1]
  sim$predict_upr[i] <- rethinking::HPDI(predictions[,i], prob = 0.95)[2]
}

## plot mean of model vs mean of raw data
plot(sim$mu_mean ~ sim$mu_std, pch = 19, col = rgb(0,0,1,0.2)) ; abline(a = 0, b = 1, lwd = 2)

## plot all predictions with points = mean of raw data
ggplot(sim)+
  geom_ribbon(aes(x = age, ymin = predict_lwr, ymax = predict_upr, fill = as.factor(window)),
              alpha = 0.2)+
  geom_ribbon(aes(x = age, ymin = mu_lwr, ymax = mu_upr, fill = as.factor(window)),
              alpha = 0.4)+
  geom_line(aes(x = age, y = mu_mean, colour = as.factor(window)))+
  geom_point(aes(x = age, y = mu_std, colour = as.factor(window)))+
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  facet_wrap(. ~ as.factor(window))+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'age (years)')

## plot on unstandardised scale ####
## revert predictions to unstandardised scale
mu_ustd <- mu * sd(sim$mu) + mean(sim$mu)
predictions_ustd <- predictions * sd(sim$mu) + mean(sim$mu)

sim$mu_mean_ustd <- apply(mu_ustd, 2, mean)
sim$mu_lwr_ustd <- NA ; sim$mu_upr_ustd <- NA
sim$predict_lwr_ustd <- NA ; sim$predict_upr_ustd <- NA
for( i in 1:nrow(sim) ){
  sim$mu_lwr_ustd[i] <- rethinking::HPDI(mu_ustd[,i], prob = 0.95)[1]
  sim$mu_upr_ustd[i] <- rethinking::HPDI(mu_ustd[,i], prob = 0.95)[2]
  sim$predict_lwr_ustd[i] <- rethinking::HPDI(predictions_ustd[,i], prob = 0.95)[1]
  sim$predict_upr_ustd[i] <- rethinking::HPDI(predictions_ustd[,i], prob = 0.95)[2]
}

ggplot(sim)+
  geom_ribbon(aes(x = age, ymin = predict_lwr_ustd, ymax = predict_upr_ustd, fill = as.factor(window)),
              alpha = 0.2)+
  geom_ribbon(aes(x = age, ymin = mu_lwr_ustd, ymax = mu_upr_ustd, fill = as.factor(window)),
              alpha = 0.4)+
  geom_line(aes(x = age, y = mu_mean_ustd, colour = as.factor(window)))+
  geom_point(aes(x = age, y = mu, colour = as.factor(window)))+
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  facet_wrap(. ~ as.factor(window))+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'standardised age')

## extract original values from output -- I CAN'T WORK OUT WHY CALCULATING THE MARGNAL EFFECTS FROM THE PREDICTIONS WOULD BE ANY DIFFERENT TO THOSE PRESENTED IN THE SUMMARY, GIVEN THAT THERE IS ONLY ONE PREDICTOR VARIABLE. ALSO, IF I'M MAKING PREDICTIONS BASED ON THE RAW DATA RATHER THAN A HYPOTHETICAL/COUNTERFACTUAL DATASET, HOW DO YOU CALCULATE THE EQUIVALENT OF THE DO OPERATOR? WOULD I BASICALLY JUST ADD 1 YEAR TO THE AGE OF ALL INDIVIDUALS AND REDO THE PREDICTIONS, THEN LOOK AT THE DIFFERENCE? ####
