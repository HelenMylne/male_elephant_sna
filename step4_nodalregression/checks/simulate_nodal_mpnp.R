#### set up #### 
library(LaplacesDemon)
library(tidyverse)
library(cmdstanr) #set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0/')

theme_set(theme_bw())
set.seed(2)

#### simulate population ####
## define population parameters
min_age <- 11                                                # youngest individual
max_age <- 60                                                #  oldest  individual
n_nodes <- ((max_age+1) - min_age)                         # total nodes = 2 per possible age

all_nodes <- data.frame(node = 1:n_nodes,
                        age = sample(min_age:max_age, n_nodes,
                                     prob = 1/(min_age:max_age), replace = T)) # more at lower ages

nodes1 <- sample(1:n_nodes, n_nodes/2, replace = F)
nodes2 <- sample(1:n_nodes, n_nodes/2, replace = F)
nodes3 <- sample(1:n_nodes, n_nodes/2, replace = F)

length(unique(c(nodes1,nodes2,nodes3)))

## simulate first half of each population -- all ages drawn from the same baseline data
sim1.1 <- data.frame(node = nodes1) %>% 
  left_join(all_nodes, by = 'node') %>% 
  mutate(age = age - 2) %>% 
  mutate(window = 1)
sim2.1 <- data.frame(node = nodes2) %>% 
  left_join(all_nodes, by = 'node') %>% 
  mutate(window = 2)
sim3.1 <- data.frame(node = nodes3) %>% 
  left_join(all_nodes, by = 'node') %>% 
  mutate(age = age + 2) %>% 
  mutate(window = 3)

## add second part of each population -- all ages drawn in the same way, but shifted so that average age of all time windows should be the same
sim1.2 <- data.frame(node = (n_nodes+1):(2*n_nodes),
                     age = sample(min_age:max_age, n_nodes,
                                  prob = 1/(min_age:max_age), replace = T),
                     window = 1) %>% 
  mutate(age = age + 2)
sim2.2 <- data.frame(node = ((2*n_nodes)+1):(3*n_nodes),
                     age = sample(min_age:max_age, n_nodes,
                                  prob = 1/(min_age:max_age), replace = T),
                     window = 2)
sim3.2 <- data.frame(node = ((3*n_nodes)+1):(4*n_nodes),
                     age = sample(min_age:max_age, n_nodes,
                                  prob = 1/(min_age:max_age), replace = T),
                     window = 3) %>% 
  mutate(age = age - 2)

## combine populations into single data frame
sim <- rbind(sim1.1, sim1.2,
             sim2.1, sim2.2,
             sim3.1, sim3.2)

## check ages have similar average across time windows
mean(sim$age[sim$window == 1])
mean(sim$age[sim$window == 2])
mean(sim$age[sim$window == 3])

## randomise node ID to ensure there is no correlation between node and age
random_nodes <- data.frame(node = unique(sim$node)) %>% 
  mutate(node_random = sample(1:length(unique(sim$node)), replace = F))
sim <- sim %>% left_join(random_nodes, by = 'node')
sim$node_window <- paste0(sim$node_random, '_', sim$window)
rm(sim1.1, sim1.2, sim2.1, sim2.2, sim3.1, sim3.2, nodes1, nodes2, nodes3, random_nodes) ; gc()

## redefine parameters
n_data <- nrow(sim)
n_nodes <- length(unique(sim$node))
n_windows <- length(unique(sim$window))

## view ages
hist(sim$age)

## calculate age categories
sim$age_cat <- ifelse(sim$age < 15, 1,
                      ifelse(sim$age < 20, 2,
                             ifelse(sim$age < 25, 3,
                                    ifelse(sim$age < 35, 4, 5))))

#### simulate centralities ####
## simulate age effect
sim_slope <- 0.2

## simulate intercept
sim_intcp <- -1

## simulate random effects
sim_window_unique <- rnorm(n_windows, mean = 0, sd = 1)
sim_node_unique <- rnorm(n_nodes, mean = 0, sd = 0.2)

## calculate mean centrality per elephant
sim$mu <- sim$age * sim_slope + sim_intcp + sim_window_unique[sim$window] + sim_node_unique[sim$node_random]    # simulate mean centrality on outcome scale

## plot mean centrality against age
plot(mu ~ age, data = sim[sim$window == 1,], col = 'red', pch = 19, ylim = c(-3,12), las = 1)           # plot
points(mu ~ age, data = sim[sim$window == 2,], col = 'blue', pch = 19)                                 # plot
points(mu ~ age, data = sim[sim$window == 3,], col = 'green', pch = 19)                                # plot

## simulate full distribution of samples per node
sim$sd <- 1#(rpois(nrow(sim),lambda = 1)+1)/10         # can be made to vary amongst nodes, currently all elephants have equal variance in centrality (not realistic given that in the real data some are seen more often than others)
sim_dat <- matrix(data = NA, nrow = 4000, ncol = n_data, dimnames = list(NULL, sim$node_window))    # create matrix for full centrality distribution
for(j in 1:n_data){
  sim_dat[,j] <- rnorm(n = nrow(sim_dat), mean = sim$mu[j], sd = sim$sd[j])  # simulate distribution for each elephant
}

## plot full unstandardised distribution
plot(sim_dat[1,which(sim$window == 1)] ~ sim$age[which(sim$window == 1)],
     col = 'red', pch = 19, ylim = c(-20,20))       # plot simulated values against age for window 1 (first row of simulated centralities only)
points(sim_dat[1,which(sim$window == 2)] ~ sim$age[which(sim$window == 2)],
       col = 'blue', pch = 19)                    # plot simulated values against age for window 2 (first row of simulated centralities only)
points(sim_dat[1,which(sim$window == 3)] ~ sim$age[which(sim$window == 3)],
       col = 'green', pch = 19)                   # plot simulated values against age for window 3 (first row of simulated centralities only)

## visualise
data.frame(sim_dat) %>% 
#data.frame(sim_dat_std) %>% 
  pivot_longer(cols = everything(),
               names_to = "node_random", values_to = "centrality") %>% 
  separate(node_random, into = c('X','node_window'), remove = T, sep = 1) %>% 
  dplyr::select(-X) %>% 
  separate(node_window, into = c('node_random','window'), sep = '_', remove = F) %>% 
  mutate(node = as.integer(node_random)) %>% 
  left_join(sim[,c('node_window','age')], by = 'node_window') %>% 
  filter(node_random %in% sample(node_random, 40)) %>% 
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
# compute normal approximation by window -- calculate means per node
sim_cent_mu_1 <- apply(sim_dat[,which(sim$window == 1)], 2, mean) # sim_cent_mu_1 <- apply(sim_dat_std[,which(sim$window == 1)], 2, mean)
sim_cent_mu_2 <- apply(sim_dat[,which(sim$window == 2)], 2, mean) # sim_cent_mu_2 <- apply(sim_dat_std[,which(sim$window == 2)], 2, mean)
sim_cent_mu_3 <- apply(sim_dat[,which(sim$window == 3)], 2, mean) # sim_cent_mu_3 <- apply(sim_dat_std[,which(sim$window == 3)], 2, mean)

# compute normal approximation by window -- calculate covariance matrix
sim_cent_cov_1 <- cov(sim_dat[,which(sim$window == 1)]) # sim_cent_cov_1 <- cov(sim_dat_std[,which(sim$window == 1)])
sim_cent_cov_2 <- cov(sim_dat[,which(sim$window == 2)]) # sim_cent_cov_2 <- cov(sim_dat_std[,which(sim$window == 2)])
sim_cent_cov_3 <- cov(sim_dat[,which(sim$window == 3)]) # sim_cent_cov_3 <- cov(sim_dat_std[,which(sim$window == 3)])

## check normal approximation -- simulate from combined mean and covariance, plot curve against relative node ID 
par(mfrow = c(3,1))
sim_cent_samples <- MASS::mvrnorm(1e5, sim_cent_mu_1, sim_cent_cov_1)     # simulate from multivariate normal
node_id_sample <- sample(which(sim$window == 1),1)
plot(density(sim_dat[,node_id_sample]), lwd = 2, las = 1,                         # plot true density curve
#plot(density(sim_dat_std[,node_id_sample]), lwd = 2, las = 1,                         # plot true density curve
     main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
lines(density(sim_cent_samples[, node_id_sample]), col = rgb(0,0,1,0.5), lwd = 2)      # overlay normal approximation

sim_cent_samples <- MASS::mvrnorm(1e5, sim_cent_mu_2, sim_cent_cov_2)     # simulate from multivariate normal
node_id_sample <- sample(which(sim$window == 2),1)
plot(density(sim_dat[,node_id_sample]),          # plot true density curve
#plot(density(sim_dat_std[,node_id_sample]),          # plot true density curve
     lwd = 2, las = 1,
     main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
lines(density(sim_cent_samples[, node_id_sample - length(which(sim$window == 1))]), col = rgb(0,0,1,0.5), lwd = 2)      # overlay normal approximation

sim_cent_samples <- MASS::mvrnorm(1e5, sim_cent_mu_3, sim_cent_cov_3)     # simulate from multivariate normal
node_id_sample <- sample(which(sim$window == 3),1)
plot(density(sim_dat[,node_id_sample]),          # plot true density curve
#plot(density(sim_dat_std[,node_id_sample]),          # plot true density curve
     lwd = 2, las = 1,
     main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
lines(density(sim_cent_samples[, node_id_sample - length(which(sim$window < 3))]), col = rgb(0,0,1,0.5), lwd = 2)      # overlay normal approximation

par(mfrow = c(1,1))
rm(sim_cent_samples, node_id_sample) ; gc()

## create data
eigen_list <- list(num_data = n_data,
                   num_nodes = n_nodes,
                   num_windows = length(unique(sim$window)),
                   num_age_cat = length(unique(sim$age_cat)),
                   length_dirichlet = length(unique(sim$age_cat)) + 1,
                   prior_age = rep(1,length(unique(sim$age_cat))),
                   num_nodes_window1 = length(which(sim$window == 1)),
                   num_nodes_window2 = length(which(sim$window == 2)),
                   num_nodes_window3 = length(which(sim$window == 3)),
                   num_nodes_prev_windows = c(0,
                                              length(which(sim$window < 2)),
                                              length(which(sim$window < 3))),
                   centrality_mu_1 = sim_cent_mu_1,
                   centrality_mu_2 = sim_cent_mu_2,
                   centrality_mu_3 = sim_cent_mu_3,
                   centrality_cov_1 = sim_cent_cov_1,
                   centrality_cov_2 = sim_cent_cov_2,
                   centrality_cov_3 = sim_cent_cov_3,
                   node_age = sim$age_cat,
                   nodes_window1 = sim$node_random[sim$window == 1],
                   nodes_window2 = sim$node_random[sim$window == 2],
                   nodes_window3 = sim$node_random[sim$window == 3])

## check inputs
plot(eigen_list$centrality_mu_1 ~ sim$age_cat[sim$window == 1], pch = 19, col = 'red', ylim = c(-1,12), las = 1)
points(eigen_list$centrality_mu_2 ~ sim$age_cat[sim$window == 2], pch = 19, col = 'blue')
points(eigen_list$centrality_mu_3 ~ sim$age_cat[sim$window == 3], pch = 19, col = 'green')

#### prior predictive check ####
n <- 100
beta_age <- rnorm(n, 0, 3)
intercept <- rnorm(n, 6, 2.5) # intercept  <- rnorm(n, logit(0.05), 2)
age_dirichlet <- rdirichlet(n, c(1,1,1,1,1))
plot(NULL, las = 1, xlab = 'age category', ylab = 'logit(eigenvector)',
     ylim = c(min(c(sim_cent_mu_1, sim_cent_mu_2))-2,
              max(c(sim_cent_mu_1, sim_cent_mu_2))+2),
     xlim = c(min(sim$age_cat), max(sim$age_cat)))
abline(h = min(c(sim_cent_mu_1, sim_cent_mu_2, sim_cent_mu_3)), lty = 2)
abline(h = max(c(sim_cent_mu_1, sim_cent_mu_2, sim_cent_mu_3)), lty = 2)
x <- min(sim$age_cat):max(sim$age_cat)
for(i in 1:n){
  y <- rep(NA, length(x))
  for(j in 1:length(x)){
    y[j] <- intercept[i] + beta_age[i]*sum(age_dirichlet[i,][1:x[j]])
  }
  lines(x = x, y = y, col = rgb(0,0,1,0.4))
}
rm(n, beta_age, intercept, age_dirichlet, x, y) ; gc()

#### run model -- age as a continuous variable with 2 windows ####
## load model
nodal_regression <- cmdstan_model('models/eigen_regression_mpnp.stan')

## run model
n_chains <- 4
n_samples <- 1000
#fit_sim <- sampling(nodal_regression, data = eigen_list, chains = n_chains, cores = n_chains)
fit_sim <- nodal_regression$sample(data = eigen_list,
                                   chains = n_chains, parallel_chains = n_chains,
                                   iter_warmup = n_samples, iter_sampling = n_samples)

#### check outputs ####
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
params_std <- fit_sim$draws(format = 'draws_df')

## extract delta and delta_j parameters
delta <- params_std %>% 
  select(`delta[1]`,`delta[2]`,`delta[3]`,`delta[4]`,`delta[5]`,
         `.chain`,`.iteration`,`.draw`)
delta_j <- params_std %>% 
  select(`delta_j[1]`,`delta_j[2]`,`delta_j[3]`,`delta_j[4]`,`delta_j[5]`,`delta_j[6]`,
         `.chain`,`.iteration`,`.draw`)

## separate random effects from global parameters
rand_window <- params_std %>% 
  dplyr::select(grep('rand_window', colnames(params_std), value=TRUE))
check_windows <- apply(rand_window, 2, mean) %>% 
  as.data.frame() %>% 
  mutate(node = colnames(rand_window)) %>% 
  relocate(node) %>% 
  mutate(sim_input = sim_window_unique)
colnames(check_windows)[2] <- 'model_output'
plot(check_windows$model_output ~ check_windows$sim_input)

rand_node <- params_std %>% 
  dplyr::select(grep('rand_node', colnames(params_std), value=TRUE))
check_nodes <- apply(rand_node, 2, mean) %>% 
  as.data.frame() %>% 
  mutate(node = colnames(rand_node)) %>% 
  relocate(node) %>% 
  mutate(sim_input = sim_node_unique)
colnames(check_nodes)[2] <- 'model_output'
plot(check_nodes$model_output ~ check_nodes$sim_input)   # wouldn't expect this to be an actual 1:1 relationship because of the scales, but there should be a positive relationship

rm(check_windows, check_nodes) ; gc()

## traceplot all parameters
#traceplot(fit_sim, pars = c('intercept','beta_age','sigma','predictor[1]','predictor[50]','predictor[100]'))
plot_params <- c('intercept','beta_age','sigma',
                 'rand_node[1]','rand_node[50]','rand_node[100]',
                 'rand_window[1]','rand_window[2]',
                 'predictor_window1[1]','predictor_window1[16]','predictor_window1[38]','predictor_window1[50]',
                 'predictor_window2[1]','predictor_window2[25]','predictor_window2[50]','predictor_window2[75]')
params_std %>% 
  select(all_of(plot_params),`.draw`,`.chain`,`.iteration`) %>% 
  pivot_longer(cols = all_of(plot_params), names_to = 'parameter', values_to = 'draw') %>% 
  rename(chain = .chain,
         chain_position = .iteration,
         draw_id = .draw) %>% 
  #filter(chain == 1) %>% # inspect individual chains -- check for wandery sections that might be hidden by other chains when all plotted together
  ggplot(aes(x = chain_position, y = draw, colour = as.factor(chain)))+
  geom_line()+
  facet_wrap(. ~ parameter, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'none')
rm(plot_params) ; gc()
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


#### posterior predictive check ####
par(mfrow = c(3,1))

## check on standardised scale
ppcheck <- function(eigen_mat, eigen_df, cent_cov, window, params, rand_node, rand_window){
  plot(density(eigen_mat[1, which(eigen_df$window == window)]), las = 1, ylim = c(0,0.2),
       main = "Posterior predictive check:\nblack = data, blue = predicted",
       col = rgb(0, 0, 0, 0.25))
  n_nodes <- length(which(eigen_df$window == window))
  eigen_data <- data.frame(node_age = eigen_df$age_cat[which(eigen_df$window == window)],
                           window = rep(window, n_nodes),
                           nodes = eigen_df$node_random[eigen_df$window == window],
                           nodes_window = eigen_df$node_random[eigen_df$window == window])
  for (i in 1:100) {
    j <- sample(1:length(params$beta_age), 1)
    lines(density(eigen_mat[j, which(eigen_df$window == window)]), col=rgb(0, 0, 0, 0.25))
    mu <- rep(NA, length(eigen_data$node_age))
    for(k in 1:length(mu)) {
      mu[k] <- params$intercept[j] + params$beta_age[j]*sum(delta_j[j,(1:eigen_data$node_age[k])]) + as.numeric(rand_window[j,eigen_data$window[k]]) + as.numeric(rand_node[j,eigen_data$nodes[k]])
    }
    sigma <- cent_cov + diag(rep(params$sigma[j], n_nodes))
    lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
  }
}
ppcheck(eigen_mat = sim_dat, eigen_df = sim,
        cent_cov = sim_cent_cov_1, window = 1, 
        params = params_std, rand_node = rand_node, rand_window = rand_window)
ppcheck(eigen_mat = sim_dat, eigen_df = sim,
        cent_cov = sim_cent_cov_1, window = 2, 
        params = params_std, rand_node = rand_node, rand_window = rand_window)
ppcheck(eigen_mat = sim_dat, eigen_df = sim,
        cent_cov = sim_cent_cov_1, window = 3, 
        params = params_std, rand_node = rand_node, rand_window = rand_window)
par(mfrow = c(1,1))

#### predict from model -- standardised scale ####
## create mean prediction function
get_mean_predictions <- function(predict_df, delta_j, parameters, include_node = TRUE, include_window = TRUE){
  ## create empty matrix to fill with predictions
  mean_matrix <- matrix(NA, nrow = nrow(parameters), ncol = nrow(predict_df),
                           dimnames = list(NULL, predict_df$node_window))
  
  ## populate matrix = mean centrality values per node, predicting for real data
  for(i in 1:nrow(mean_matrix)){
    for(j in 1:ncol(mean_matrix)){
      mean_matrix[i,j] <- params$intercept[i] + params$beta_age[i] * sum(delta_j[i,(1:predict_df$age_cat[j])])
    }
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
mu_std <- get_mean_predictions(predict_df = sim, parameters = params_std, delta_j = delta_j,
                               include_window = TRUE, include_node = TRUE)

## add mean and CI of predicted means to input data frame for comparison
sim$mu_mean_std <- apply(mu_std, 2, mean)
sim$mu_lwr_std <- NA ; sim$mu_upr_std <- NA
for( i in 1:nrow(sim) ){
  sim$mu_lwr_std[i] <- rethinking::HPDI(mu_std[,i], prob = 0.95)[1]
  sim$mu_upr_std[i] <- rethinking::HPDI(mu_std[,i], prob = 0.95)[2]
}

## plot mean of model vs mean of raw data
ggplot()+
  geom_point(data = sim, aes(x = mu, y = mu_mean_std, colour = as.factor(window)))+
  scale_colour_viridis_d()+
  labs(colour = 'window', x = 'simulated mean (standardised)', y = 'predicted mean (standardised)')+
  geom_abline(slope = 1, intercept = 0) # add line showing where points would lie if model fit was perfect

## put together sigma arrays, separated by time window
sigma_all <- list()
for(time_window in 1:n_windows){
  cent_cov <- eigen_list[[grep('centrality_cov', names(eigen_list))[time_window] ]]
  n_nodes_window <- eigen_list[[grep('num_nodes_window', names(eigen_list))[time_window] ]]
  nodes_window <- eigen_list[[grep('nodes_window', names(eigen_list))[time_window + n_windows] ]]
  sigma_array <- array(NA, dim = c(n_nodes_window,
                                   n_nodes_window,
                                   nrow(params_std)),
                       dimnames = list(nodes_window,
                                       nodes_window,
                                       NULL))
  
  for(i in 1:nrow(params_std)){
    sigma_array[,,i] <- cent_cov + diag(rep(params_std$sigma[i], n_nodes_window))
  }
  sigma_all[[time_window]] <- sigma_array
}

## create empty matrix to take full set of predicted values per elephant
predictions_std <- matrix(NA, nrow = nrow(params_std), ncol = nrow(sim), dimnames = list(NULL, sim$node_window))

## populate matrix using mean values in matrix mu_std, and sigma values based on time window
for(time_window in 1:n_windows){
  sigma_array <- sigma_all[[time_window]]
  for(i in 1:nrow(predictions_std)){
    predictions_std[i,sim$window == time_window] <- MASS::mvrnorm(1, mu_std[i,sim$window == time_window], sigma_array[,,i])
  }
}

## add CI of predicted data points to input data frame for comparison
sim$predict_lwr_std <- apply(predictions_std, 2, quantile, prob = 0.025)
sim$predict_upr_std <- apply(predictions_std, 2, quantile, prob = 0.975)

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

#### extract original values from output -- simulated slope value originally used produces an effect on the unstandardised scale. The model works on the standardised scale. Convert predictions to unstandardised scale and then run contrasts to calculate the slope coefficient. ####
## get mean predictions
sim2 <- sim %>% 
  mutate(age_cat = ifelse(age_cat == 5, 1, age_cat + 1))
mu2_std <- get_mean_predictions(predict_df = sim2, parameters = params_std, delta_j = delta_j,
                                include_window = TRUE, include_node = TRUE)
sim$mu_mean_plus1 <- apply(mu2_std, 2, mean)

## full distribution of predictions using age_std + 1 stdev
predictions2_std <- matrix(NA, nrow = nrow(params_std), ncol = nrow(sim2), dimnames = list(NULL, sim$node_window))
for(time_window in 1:n_windows){
  sigma_array <- sigma_all[[time_window]]
  for(i in 1:nrow(predictions2_std)){
    predictions2_std[i,sim$window == time_window] <- MASS::mvrnorm(1, mu2_std[i,sim$window == time_window], sigma_array[,,i])
  }
}

## contrast predictions on standardised scale -- check that this returns the marginal effect presented in the summary
contrast_std <- predictions2_std - predictions_std    # contrast between predicted values for raw data and all same data but add 1 to age
head(contrast_std[,1:5])                              # check matrix looks right
mean(params_std$beta_age)                                 # output parameter direct from model
mean(contrast_std)                                    # should be very similar to mean of output parameter
quantile(contrast_std, prob = c(0.025, 0.975))        # very wide

## contrast predictions on unstandardised scale
contrast <- contrast_std / sd(sim$age)                # contrast between predicted values for raw data and all same data but add 1 to age
head(contrast[,1:5])                                  # check matrix looks right
sim_slope                                             # output parameter direct from model
mean(contrast)                                        # should be very similar to mean of output parameter
quantile(contrast, prob = c(0.025, 0.975))            # very wide

## extract original values from output -- NOTE: CURRENT METHOD WILL NOT WORK FOR REAL DATA AS I DON'T HAVE THE AGE VALUES IN YEARS TO GO ABOUT ADDING 1 YEAR TO EVERY ELEPHANT. WILL NEED TO COME UP WITH AN ALTERNATIVE OR JUST SHOW CONTRASTS BETWEEN AGE CATEGORIES FOR REAL DATA, BUT HERE CHECK THE EFFECT ON REAL AGE SCALE TO BE SURE THAT THE MODEL IS WORKING. ####
## set up objects to store predictions
pred_means <- list()
pred_fulls <- list()

## predict for ages from 10 years younger than current up to 10 years older: multiple opportunities to pass through age category thresholds 
for( age_difference in -10:10 ){
  pred_data_new <- sim[order(sim$age_cat),] %>% 
    mutate(age = age + age_difference) %>% 
    mutate(age_cat = ifelse(age <= 15, 1,
                            ifelse(age <= 20, 2,
                                   ifelse(age <= 25, 3,
                                          ifelse(age <= 40, 4, 5)))))
  pred_mu_new <- get_mean_predictions(parameters = params, delta_j = delta_j, predict_df = pred_data_new,
                                      include_node = TRUE, include_window = TRUE)
  
  pred_full_new <- matrix(NA, nrow = nrow(params_std), ncol = nrow(pred_data_new), dimnames = list(NULL, sim$node_window))
  for(time_window in 1:n_windows){
    sigma_array <- sigma_all[[time_window]]
    for(i in 1:nrow(pred_full_new)){
      pred_full_new[i,pred_data_new$window == time_window] <- MASS::mvrnorm(1,
                                                                            pred_mu_new[i,pred_data_new$window == time_window],
                                                                            sigma_array[,,i])
    }
  }
  
  pred_means[[(age_difference+11)]] <- pred_mu_new
  pred_fulls[[(age_difference+11)]] <- pred_full_new
  print(age_difference)
}
warnings() # should just be "Dropping 'draws_df' class as required metadata was removed." over and over

## create data frame to store outputs
contrasts_all <- data.frame(older = rep(-10:10, each = 21),
                            younger = rep(-10:10, 21),
                            mean_contrast = NA)

## for each combination of predictions, calculate contrast and divide difference by number of years changed
for(i in 1:21){
  for(j in 1:21){
    if(i > j){
      contrast <- (pred_fulls[[i]] - pred_fulls[[j]]) / (i-j)
      contrasts_all$mean_contrast[contrasts_all$older == (i - 11) & 
                                    contrasts_all$younger == (j - 11)] <- mean(contrast)
    }
  }
}

## remove impossible age combinations
contrasts_all <- contrasts_all %>%
  filter(is.na(mean_contrast) == FALSE)

## compare to input
sim_slope
mean(contrasts_all$mean_contrast)
quantile(contrasts_all$mean_contrast, prob = c(0.025, 0.975))

#### final "clean" plots using hypothetical data ####
## create dummy dataset
newdat <- sim %>% 
  select(node_random, age, age_cat, window)

## get mean predictions
fake_mean <- get_mean_predictions(predict_df = newdat, delta_j = delta_j, parameters = params_std,
                                  include_window = TRUE, include_node = FALSE)
newdat$predict_mean <- apply(fake_mean, 2, mean)
newdat$predict_mean_lwr <- apply(fake_mean, 2, quantile, prob = 0.025)
newdat$predict_mean_upr <- apply(fake_mean, 2, quantile, prob = 0.975)

## create prediction matrix
fake_pred <- matrix(NA, nrow = nrow(params_std), ncol = nrow(newdat))
for(time_window in 1:n_windows){
  sigma_array <- sigma_all[[time_window]]
  for(i in 1:nrow(fake_pred)){
    fake_pred[i,newdat$window == time_window] <- MASS::mvrnorm(1,
                                                               fake_mean[i,newdat$window == time_window],
                                                               sigma_array[,,i])
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
  pivot_longer(cols = everything(), names_to = 'node_window', values_to = 'eigen') %>% 
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
