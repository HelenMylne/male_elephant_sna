#### information ####
# script takes data input from edge weight estimation for ANP population (input = 1000 draws for each of 4 chains per dyad posterior distribution)
# runs through dyadic regression as specified by Jordan Hart in BISoN examples (https://github.com/JHart96/bison_examples/blob/main/examples/dyadic_regression_stan.md)

#### set up ####
#library(tidyverse) ; library(LaplacesDemon) ; library(car) ; library(cmdstanr) ; library(bisonR) ; library(brms)
library(cmdstanr, lib.loc = '../packages/')  # library(cmdstanr)
library(tidyverse, lib.loc = '../packages/') # library(tidyverse)
library(car, lib.loc = '../packages/')       # library(car)
#library(bisonR, lib.loc = '../packages/')    # library(bisonR)
#library(brms, lib.loc = '../packages/')      # library(brms)
library(LaplacesDemon, lib.loc = '../packages/')

set_cmdstan_path('R:/rsrch/df525/phd/hkm513/packages/.cmdstan2/cmdstan-2.33.1/')

load('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')
# rm(edgelist, edge_binary, nodes) ; gc()

pdf('../outputs/anp1_dyadicregression_plots.pdf')

#### prior predictive check ####
age_min <- 5:50
age_max <- seq(10, 50, length.out = 12)
par(mfrow = c(4,3))
for(age in age_max){
 plot(NULL, xlim = c(5,50), ylim = c(0,1), las = 1,
      xlab = 'age of younger', ylab = 'edge weight', main = paste0('oldest = ',round(age)))
 for(i in 1:100){
   beta_min <- rnorm(1, 0, 0.1)
   beta_max <- rnorm(1, 0, 0.1)
   #beta_int <- rnorm(1, 0, 0.005)
   min_new <- age_min[age > age_min]
   lines(x = min_new, y = plogis(min_new*beta_min +
                                   #min_new*age*beta_int +
                                   age*beta_max),
         col = rgb(0,0,1,0.2))
 }
}
rm(age, age_min, age_max, beta_max, beta_min, i, min_new) ; gc()
par(mfrow = c(1,1))
# 
# age_min <- 5:50
# age_diff <- c(0,2,5,10,15,20,25,30)
# par(mfrow = c(4,2), mai = c(0.3,0.3,0.3,0.1))
# for(age in age_diff){
#   plot(NULL, xlim = range(age_min), ylim = c(0,1), las = 1,
#        xlab = 'age of younger', ylab = 'edge weight',
#        main = paste0('diff = ',age))
#   for(i in 1:100){
#     beta_min <- rnorm(1, 0, 0.1)
#     beta_diff <- rnorm(1, 0, 0.1)
#     lines(x = age_min, y = plogis(age_min*beta_min + 
#                                     age*beta_diff), 
#           col = rgb(0,0,1,0.2))
#     lines(x = age_min, y = age_min*beta_min + age*beta_diff, 
#           col = rgb(1,0,0,0.2))
#   }
# }
# rm(age, age_min, age_max, beta_max, beta_min, i, min_new) ; gc()
# par(mfrow = c(1,1))

# add time marker
print(paste0('prior predictive check completed at ', Sys.time()))

#### fit multivariate Gaussian distribution to output of edge weight model ####
# To parameterise the multivariate normal approximation, we use the sample mean and covariance matrix, calculated from the posterior edge weight samples using the following code:
### get the weights on the logit scale (they are not currently because we used a beta prior and identity link here rather than logistic link)
logit_weights <- apply(edge_samples, 2, qlogis)

### fit a multivariate normal dist to the edges -- quantities will be given to Stan model as data to model joint posteriors of edge weight in the regression
logit_edge_draws_mu <- apply(logit_weights, 2, mean)
logit_edge_draws_cov <- cov(logit_weights)

#### plot to see how well the approximation is working ####
### Randomly selecting samples to examine
num_check <- 20
selected_samples <- sample(1:n_dyads, num_check, replace = FALSE)

### Setting grid layout
rows <- floor(sqrt(num_check))
cols <- ceiling(num_check / rows)
par(mfrow=c(rows, cols), mar=c(2,2,2,1))

### plot
for (i in selected_samples) {
  mu <- logit_edge_draws_mu[i]
  sd <- sqrt(logit_edge_draws_cov[i,i])
  
  fitted_values_logit <- rnorm(1e5, mean = mu, sd = sd)
  fitted_values_original <- plogis(fitted_values_logit)
  
  hist(unlist(edge_samples[,i]), probability = TRUE, las = 1,
       main = paste("Dyad", i), xlab = "Value", breaks = 50)
  lines(density(fitted_values_original), col="blue", lw=1.5)
}

for (i in selected_samples) {
  mu <- logit_edge_draws_mu[i]
  sd <- sqrt(logit_edge_draws_cov[i,i])
  
  fitted_values_logit <- rnorm(1e5, mean = mu, sd = sd)
  fitted_values_original <- plogis(fitted_values_logit)
  
  plot(unlist(edge_samples[,i]), unlist(edge_samples[,i+1]), col = rgb(0,0,1,0.05), las = 1,
       main = paste("cov ", i ,"&",i+1))
}

### reset plot window
par(mfrow=c(1,1))

# save image so far
save.image('anpshort1_dyadicregression.RData')

# add time marker
print(paste0('multivariate Gaussian approximation fitted at ', Sys.time()))

#### fit dyadic regression ####
### identify older and younger of dyad
cdf_1$age_min <- NA ; cdf_1$age_max <- NA
for(i in 1:nrow(cdf_1)){
  x <- c(cdf_1$age_start_1[i],cdf_1$age_start_2[i])
  cdf_1$age_min[i] <- min(x)
  cdf_1$age_max[i] <- max(x)
}

## convert node IDs to values 1:52 not numbers based on casename
all_node_IDs <- unique(c(cdf_1$id_1, cdf_1$id_2))
n_nodes <- length(all_node_IDs)
cdf_1 <- cdf_1 %>% 
  rename(node_1_original = node_1,
         node_2_original = node_2) %>% 
  mutate(node_1_period = as.integer(factor(id_1, levels = all_node_IDs)),
         node_2_period = as.integer(factor(id_2, levels = all_node_IDs)),)

## create data list
dyad_data <- list(
  num_dyads = n_dyads,                      # number of dyads
  num_nodes = n_nodes,                      # number of nodes
  logit_edge_mu = logit_edge_draws_mu,      # sample means of the logit edge weights
  logit_edge_cov = logit_edge_draws_cov,    # sample covariance of logit edge weights
  age_min = cdf_1$age_min,                  # age of younger dyad member
  age_max = cdf_1$age_max,                  # age of  older  dyad member
  #age_diff = cdf_1$age_diff,                # age difference between dyad members
  node_1 = cdf_1$node_1_period,             # node IDs for multimembership effects
  node_2 = cdf_1$node_2_period              # node IDs for multimembership effects
  #jitter = 1e-6                             # jitter to add to the diag of the cov matrix for numerical stability -- COME BACK TO THIS, MAY NEED TO ALTER THE MODEL TO ADD SIGMA BACK IN AND REMOVE JITTER
)

## load dyadic regression model
dyadic_regression <- cmdstan_model('models/dyadic_regression.stan')

# add time marker
print(paste0('start model run at ', Sys.time()))

## fit dyadic regression
fit_dyadreg_anp1 <- dyadic_regression$sample(
  data = dyad_data,
  chains = n_chains,
  parallel_chains = n_chains)

# save image so far
save.image('anpshort1_dyadicregression_minmax.RData')

# add time marker
print(paste0('finish model run at ', Sys.time()))

#### check outputs ####
#load('anpshort1_dyadicregression_minmax.RData')
rm(counts_df, dyad_data, edge_binary, edgelist, edges, i, x, make_edgelist, plot_network_threshold_anp) ; gc()

# obtain summary
fit_dyadreg_anp1$summary()

## extract draws
draws <- fit_dyadreg_anp1$draws(format = 'df')

## extract dyadic regression slopes
b_diff <- draws[,,'beta_age_diff']
b_min <- draws[,,'beta_age_min']
sigma <- draws[,,'sigma']

#b_max <- draws[,,'b_max']
#b_int <- draws[,,'b_int']

## extract overall age effects
dimnames(draws)$variable
draws <- draws[,,4:length(draws[1,1,])]   # WILL CHANGE IF YOU ADD ANY MORE REGRESSION SLOPES

## extract multi-membership samples -- CHECK THIS IS IN THE RIGHT ORDER
mm_matrix <- draws[,,1:n_nodes]
sigma_mm <- draws[,,n_nodes+2]
predictor <- draws[,,(n_nodes+3):length(draws[1,1,])]
length(draws[1,1,]) == length(mm_matrix[1,1,]) + length(sigma_mm[1,1,]) + length(sigma[1,1,]) + length(predictor[1,1,])
rm(draws) ; gc()

# add time marker
print(paste0('parameters extracted at ', Sys.time()))

## save parameter values
extract_slopes <- function(draws, n_samples = 1000, n_chains = 4){
  draws <- as.data.frame(draws) %>% 
    pivot_longer(everything(), names_to = 'chain', values_to = 'slope_draw') %>% 
    separate(chain, sep = c(1,2), remove = T, into = c('chain','.', 'parameter')) %>% 
    select(-.)
  draws$chain <- as.numeric(draws$chain)
  draws$position <- rep(1:n_samples, each = n_chains)
  draws$draw_id <- draws$position + (draws$chain-1)*n_samples
  return(draws)
}
b_diff <- extract_slopes(b_diff)
b_min <- extract_slopes(b_min)
#b_max <- extract_slopes(b_max)
sigma <- extract_slopes(sigma)
#b_int <- extract_slopes(b_int)
parameters <- rbind(b_min, b_diff, #b_max, #b_int,
                    sigma)

## save data 
saveRDS(parameters, '../data_processed/step5_dyadicregression/anp1_dyadicregression_slopeparameters.RDS')
#parameters <- readRDS('../data_processed/step5_dyadicregression/anp1_dyadicregression_slopeparameters.RDS')

# add time marker
print(paste0('parameters saved to file at ', Sys.time()))

pdf('../outputs/anp1_dyadicregression_plots2.pdf')

## traceplots
ggplot(data = parameters)+
  geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
  theme_classic()+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap( . ~ parameter , scales = 'free_y')

# add time marker
print(paste0('traceplots run at ', Sys.time()))

#### plot edges against age values ####
#edges <- readRDS('../data_processed/step3_edgeweightestimation/anpshort1_edgedistributions_conditionalprior.RDS')
#parameters <- readRDS('../data_processed/step5_dyadicregression/anp1_dyadicregression_slopeparameters.RDS')

### create mean data frame to plot average values over full distribution
mean_edges <- data.frame(dyad = cdf_1$dyad_id,
                         node_1 = cdf_1$node_1, node_2 = cdf_1$node_2,
                         together = cdf_1$event_count, count_dyad = cdf_1$period_count_dyad,
                         age_1 = cdf_1$age_start_1, age_2 = cdf_1$age_start_2,
                         age_min = NA, age_max = NA,
                         age_diff = cdf_1$age_diff)
for(i in 1:nrow(mean_edges)){
  mean_edges$age_min[i] <- min(c(cdf_1$age_start_1[i], cdf_1$age_start_2[i]))
  mean_edges$age_max[i] <- max(c(cdf_1$age_start_1[i], cdf_1$age_start_2[i]))
}
length(which( (mean_edges$dyad == sort(unique(mean_edges$dyad)))  == FALSE ))
edge_means <- as.data.frame(apply(edge_samples, 2, mean))
colnames(edge_means) <- 'edge_mean'
mean_edges$edge_mean <- edge_means$edge_mean

### add categories for box plots
mean_edges$age_cat_min <- ifelse(mean_edges$age_min < 15, '10-15',
                                 ifelse(mean_edges$age_min < 20, '15-20',
                                        ifelse(mean_edges$age_min < 25, '20-25',
                                               ifelse(mean_edges$age_min < 40, '25-40', '40+'))))
mean_edges$age_cat_max <- ifelse(mean_edges$age_max < 15, '10-15',
                                 ifelse(mean_edges$age_max < 20, '15-20',
                                        ifelse(mean_edges$age_max < 25, '20-25',
                                               ifelse(mean_edges$age_max < 40, '25-40', '40+'))))

#### plot raw data ####
theme_set(theme_classic())
ggplot()+
  geom_point(data = mean_edges, aes(x = age_min, y = edge_mean),
             shape = 19, colour = rgb(0,0,1,0.1)
  )+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')
ggplot()+
  geom_point(data = mean_edges, aes(x = age_min, y = edge_mean,
                                    colour = age_max), # this makes it look more interesting than it is -- there is a definite pattern at first glance, but actually appears to be nothing more than that the older the youngest male is, the older the oldest must be to be older than the youngest, whereas when the youngest is very young the oldest doesn't have to be very old to still be the oldest
             shape = 19)+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')
ggplot()+
  geom_point(data = mean_edges, aes(x = age_min, y = edge_mean, colour = age_diff), # much less exciting now
             shape = 19)+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')

edges <- left_join(edges, mean_edges, by = 'dyad')
ggplot()+
  geom_point(data = edges, aes(x = age_min, y = edge_draw),
             shape = 19, colour = rgb(0,0,1,0.01))+
  geom_point(data = mean_edges, aes(x = age_min, y = edge_mean),
             shape = 19, colour = 'red')+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')
ggplot()+
  geom_violin(data = edges, aes(x = age_cat_min, y = edge_draw))+
  geom_jitter(data = mean_edges, aes(x = age_cat_min, y = edge_mean, size = count_dyad),
              shape = 1, width = 0.3, colour = rgb(0,0,1,0.1))+
  scale_x_discrete('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')

ggplot()+
  geom_point(data = edges, aes(x = count_dyad, y = edge_draw), colour = rgb(0,0,1,0.01))+
  geom_point(data = mean_edges, aes(x = count_dyad, y = edge_mean), colour = 'red')+
  scale_x_continuous('number of sightings of dyad')+
  scale_y_continuous('edge weight')

ggplot()+
  geom_point(data = mean_edges, aes(x = age_min, y = age_max, colour = edge_mean))+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('age of older dyad member')

## heat map: x = min age, y = max age, colour = edge weight
ggplot()+
  geom_tile(data = mean_edges, mapping = aes(x = age_min, y = age_max, fill = edge_mean))+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('age of older dyad member')

# add time marker
print(paste0('age vs edge predictions completed at ', Sys.time()))

#### plot predictions ####
## posterior predictive check
b_min <- parameters %>% filter(parameter == 'beta_age_min')
b_diff <- parameters %>% filter(parameter == 'beta_age_diff')
sigma <- parameters %>% filter(parameter == 'sigma')
plot(density(as.numeric(edge_samples[1, ])), main = "Posterior predictive density of responses (edge weights)",
     ylim = c(0, 20), col = rgb(0, 0, 0, 0.25))
for (i in 1:100) {
  j <- sample(1:1000, 1)
  mu_plot <- b_min$slope_draw[j]*dyad_data$age_min + b_diff$slope_draw[j]*dyad_data$age_diff
  sigma_plot <- dyad_data$logit_edge_cov + diag(rep(sigma$slope_draw[j], n_dyads))
  mv_norm <- MASS::mvrnorm(1, mu_plot, sigma_plot)
  
  lines(density(as.numeric(edge_samples[j, ])), col = rgb(0, 0, 0, 0.25)) # black lines for edge samples
  lines(density(plogis(mv_norm)), col = rgb(1, 0, 0, 0.25))   # red lines for predictions
  # 
  # sigma_plot <- diag(rep(sigma$slope_draw[j], n_dyads))
  # mv_norm <- MASS::mvrnorm(1, mu_plot, sigma_plot)
  # lines(density(plogis(mv_norm)), col = rgb(0, 0, 1, 0.25))   # blue lines for predictions without covariance
  # 
} # definitely distinct but following the same trends at least

## obtain median and 95% CI
plot(density(b_diff$slope_draw), xlim = c(min(b_diff$slope_draw),0)) ; abline(v = 0, lty = 2)
plot(density(b_min$slope_draw), xlim = c(min(b_min$slope_draw),0)) ; abline(v = 0, lty = 2)

(b_diff_summary <- round(quantile(b_diff$slope_draw, probs=c(0.5, 0.025, 0.975)), 2))
(b_min_summary <- round(quantile(b_min$slope_draw, probs=c(0.5, 0.025, 0.975)), 2))
(sigma_summary <- round(quantile(sigma$slope_draw, probs=c(0.5, 0.025, 0.975)), 2))

## calculate predictions
# pred <- data.frame(age_min = rep(rep(5:45, each = 7),1000),
#                    age_diff = rep(seq(0,30, by = 5),41000),
#                    b_min = b_min$slope_draw[b_min$chain == 1],
#                    b_diff = b_diff$slope_draw[b_diff$chain == 1],
#                    sigma = sigma$slope_draw[sigma$chain == 1]) %>%
#   mutate(mu = b_min*age_min + b_diff*age_diff) %>%
#   mutate(sim = rnorm(1,mu,sigma))#MASS::mvrnorm(1,pred$mu[1], (logit_edge_draws_cov + diag(rep(pred$sigma[1], n_dyads))) ) )

pred <- data.frame(age_min = rep(rep(seq(from = min(c(cdf_1$age_start_1, cdf_1$age_start_2)),
                                         to = max(c(cdf_1$age_start_1, cdf_1$age_start_2)),
                                         by = 1), each = 20),3),
                   age_difference = rep(seq(from = min(cdf_1$age_diff),
                                      to = max(cdf_1$age_diff),
                                      by = 2),(39*3)),
                   b_min_mid = b_min_summary[1],
                   b_min_lwr = b_min_summary[2],
                   b_min_upr = b_min_summary[3],
                   b_dif_mid = b_diff_summary[1],
                   b_dif_lwr = b_diff_summary[2],
                   b_dif_upr = b_diff_summary[3]) %>% #,sigma = sigma$slope_draw[sigma$chain == 1]) %>%
  mutate(mu_mid = b_min_mid*age_min + b_dif_mid*age_difference,
         mu_lwr = b_min_lwr*age_min + b_dif_lwr*age_difference,
         mu_upr = b_min_upr*age_min + b_dif_upr*age_difference) %>%
  #mutate(sim = rnorm(1,mu,sigma))#MASS::mvrnorm(1,pred$mu[1], (logit_edge_draws_cov + diag(rep(pred$sigma[1], n_dyads))) ) )
  mutate(invlogit_mu_mid = invlogit(mu_mid),
         invlogit_mu_lwr = invlogit(mu_lwr),
         invlogit_mu_upr = invlogit(mu_upr))

ggplot()+
  geom_ribbon(data = pred[pred$age_difference %in% c(0,8,16,24,30,38),],
              aes(x = age_min, group = age_difference, fill = age_difference,
                  ymin = invlogit_mu_lwr, ymax = invlogit_mu_upr),
              alpha = 0.3)+
  geom_line(data = pred[pred$age_difference %in% c(0,8,16,24,30,38),],
            aes(x = age_min, y = invlogit_mu_mid,
                colour = age_difference, group = age_difference),
            linewidth = 1)+
  geom_point(data = mean_edges, aes(x = age_min, y = edge_mean, colour = age_diff), # much less exciting now
            #shape = 19,
            alpha = 0.5)+
  scale_colour_viridis_c()+ scale_fill_viridis_c()+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+
  labs(colour = 'age difference', fill = 'age difference')
ggsave('../outputs/anpshort1_dyadicregression_conditionalprior.png')


save.image('anpshort1_dyadicregression_plotting_minmax.RData')













#################
load('TESTING_DYADIC_REGRESSION_17THAUGUST2023_minmax.RData') # delete this once you've loaded it back in tomorrow!
rm(counts_df, dyad_data, edge_binary, edgelist, fit_edges_anp1, parameters, sigma_plot, summary, x, all_node_IDs, i, j, mm_matrix, mu_plot, mv_norm, n_windows, periods, predictor, sigma_mm, extract_slopes, make_edgelist, plot_network_threshold_anp) ; gc










# compute posterior means
ages <- min(nodes$age):max(nodes$age)
age_differences <- 0:30
mu <- array(NA, dim = c(n_samples, length(ages), length(age_differences)),
            dimnames = list(1:n_samples, ages, age_differences))
for(draw in 1:n_samples){
  for(minimum in 1:length(ages)){
    for(difference in 1:length(age_differences)){
      #if( min <= diff ){
        mu[draw,minimum,difference] <- b_min$slope_draw[draw]*ages[minimum] + b_diff$slope_draw[draw]*age_differences[difference]# + b_int$slope_draw[sample]*age_min[min]*age_max[max]
      #}
    }
  }
}
logit_mu <- plogis(mu)     # IS THIS ONLY NECESSARY IF PRODUCING VALUES ON LOGIT NORMAL SCALE RATHER THAN BETA?? DON'T THINK SO BECAUSE WORKING FROM THE APPROXIMATION NOT THE ACTUAL EDGE DISTRIBUTIONS

mu[10,,] # 10th draw for all pairs of min and diff
mu[,16,] # 16 years minimum age for all age differences (max age = 16-46)
mu[,,25] # 25 years age difference (max age = 30-45)

#mu_max_age <- function(edge_prediction_array, age_threshold){
#  minimum_ages <- as.numeric(colnames(edge_prediction_array[1,,]))
#  maximum_ages <- as.numeric(names(edge_prediction_array[1,1,]))
#  predictions <- edge_prediction_array[,which(minimum_ages <= age_threshold),which(maximum_ages == age_threshold)]
#  if(is.null(dim(predictions)) == TRUE) { mean_edge <- mean(predictions) } else { mean_edge <- apply(predictions, 2, mean) }
#  return(mean_edge)
#}

#pi_max_age <- function(edge_prediction_array, age_threshold){
#  minimum_ages <- as.numeric(colnames(edge_prediction_array[1,,]))
#  maximum_ages <- as.numeric(names(edge_prediction_array[1,1,]))
#  predictions <- edge_prediction_array[,which(minimum_ages <= age_threshold),which(maximum_ages == age_threshold)]
#  if(is.null(dim(predictions)) == TRUE) { 
#    mean_edge <- rethinking::HPDI(predictions) 
#  } else { 
#    mean_edge <- apply(predictions, 2, rethinking::HPDI) 
#  }
#  return(mean_edge)
#}

#mu_mean_max20 <- mu_max_age(mu, 20)
#mu_mean_max25 <- mu_max_age(mu, 25)
#mu_mean_max30 <- mu_max_age(mu, 30)
#mu_mean_max35 <- mu_max_age(mu, 35)
#mu_mean_max41 <- mu_max_age(mu, 41)
#mu_mean_max46 <- mu_max_age(mu, 46)

#pi_mean_max20 <- pi_max_age(mu, 20)
#pi_mean_max25 <- pi_max_age(mu, 25)
#pi_mean_max30 <- pi_max_age(mu, 30)
#pi_mean_max35 <- pi_max_age(mu, 35)
#pi_mean_max41 <- pi_max_age(mu, 41)
#pi_mean_max46 <- pi_max_age(mu, 46)

quant_95 <- function(x) { quantile(x, probs = c(0.025,0.975), na.rm = T) }
mu_diff_age <- function(edge_prediction_array, min_threshold, diff_threshold){
  minimum_ages <- as.numeric(rownames(edge_prediction_array[1,,]))
  difference_ages <- as.numeric(names(edge_prediction_array[1,1,]))
  predictions <- edge_prediction_array[,which(minimum_ages <= min_threshold),
                                       which(difference_ages == diff_threshold)]
  if(is.null(dim(predictions)) == TRUE) { 
    mean_edge <- mean(predictions) 
  } else { 
      mean_edge <- apply(predictions, 2, mean) 
      }
  return(mean_edge)
}
pi_diff_age <- function(edge_prediction_array, min_threshold, diff_threshold){
  minimum_ages <- as.numeric(rownames(edge_prediction_array[1,,]))
  difference_ages <- as.numeric(names(edge_prediction_array[1,1,]))
  predictions <- edge_prediction_array[,which(minimum_ages <= min_threshold),
                                       which(difference_ages == diff_threshold)]
  if(is.null(dim(predictions)) == TRUE) { 
    pi_edge <- quant_95(as.numeric(predictions))
  } else { 
    pi_edge <- apply(predictions, 2, quant_95) 
  }
  return(pi_edge)
}

mu_means_list <- data.frame(min_age = rep(ages, each = length(age_differences)),
                            diff_age = rep(age_differences, length(ages)),
                            mean = NA, hpdi = NA)
for(i in 1:nrow(mu_means_list)){
  mu_means_list$mean[i] <- list(mu_diff_age(edge_prediction_array = mu,
                                            min_threshold = mu_means_list$min_age[i],
                                            diff_threshold = mu_means_list$diff_age[i]))
  mu_means_list$hpdi[i] <- list(pi_diff_age(edge_prediction_array = mu,
                                            min_threshold = mu_means_list$min_age[i],
                                            diff_threshold = mu_means_list$diff_age[i]))
}
pi <- as.data.frame(mu_means_list$hpdi[i])
pi_lower <- pi[1,] %>% 
  pivot_longer(everything(), names_to = 'x_min', values_to = 'lower_bound') %>% 
  separate(x_min, into = c('x','min_age'), remove = T, sep = 1) %>% 
  mutate(min_age = as.numeric(min_age)) %>% 
  select(-x)
pi_upper <- pi[2,] %>% 
  pivot_longer(everything(), names_to = 'x_min', values_to = 'upper_bound') %>% 
  separate(x_min, into = c('x','min_age'), remove = T, sep = 1) %>% 
  mutate(min_age = as.numeric(min_age)) %>% 
  select(-x)
pi <- left_join(pi_lower, pi_upper, by = 'min_age') %>% 
  mutate(diff_age = mu_means_list$diff_age[i]) %>% 
  relocate(diff_age, .after = min_age)

for(i in 1:nrow(mu_means_list)){
  pi_loop <- as.data.frame(mu_means_list$hpdi[i])
  if(ncol(pi_loop) > 1){
    pi_loop_lower <- pi_loop[1,] %>% 
      pivot_longer(everything(), names_to = 'x_min', values_to = 'lower_bound') %>% 
      separate(x_min, into = c('x','min_age'), remove = T, sep = 1) %>% 
      mutate(min_age = as.numeric(min_age)) %>% 
      select(-x)
    pi_loop_upper <- pi_loop[2,] %>% 
      pivot_longer(everything(), names_to = 'x_min', values_to = 'upper_bound') %>% 
      separate(x_min, into = c('x','min_age'), remove = T, sep = 1) %>% 
      mutate(min_age = as.numeric(min_age)) %>% 
      select(-x)
    pi_loop <- left_join(pi_loop_lower, pi_loop_upper, by = 'min_age') %>% 
      mutate(diff_age = mu_means_list$diff_age[i])
  } else {
    pi_loop <- data.frame(min_age = mu_means_list$min_age[i],
                            diff_age = mu_means_list$diff_age[i],
                            lower_bound = pi_loop[1,1],
                            upper_bound = pi_loop[2,1])
  }
  pi <- rbind(pi, pi_loop)
}

mu_means <- pi %>% 
  left_join(mu_means_list[,c('min_age','diff_age','mean')],
            by = c('min_age','diff_age')) %>% 
  unnest(cols = mean) %>% 
  distinct()

rm(pi, pi_loop, pi_loop_lower, pi_loop_upper, pi_lower, pi_upper) ; gc()

for(difference in 1:length(age_differences)){
  
  ## simulate from posterior
  sim_edges <- matrix(nrow = n_samples, ncol = length(ages))
  for(draw in 1:nrow(sim_edges)){
    for(minimum in 1:ncol(sim_edges)){
      mean_sim <- b_min$slope_draw[draw]*ages[minimum] + b_diff$slope_draw[draw]*age_differences[difference]
      sigma_sim <- sigma$slope_draw[draw]
      sim_edges[draw,minimum] <- MASS::mvrnorm(n = 1, mu = mean_sim, Sigma = sigma_sim)
      #sim_edges[i,j] <- LaplacesDemon::rmvnc(n = 1, mu = mean_sim, U = sigma_sim)
    }
  }
  
  hist(sim_edges)  # almost entirely equal
  sim_mu <- apply(sim_edges, 2, mean)
  sim_pi <- apply(sim_edges, 2, quant_95)
  
  # plot
  plot(data = mean_edges[mean_edges$age_diff == age_differences[difference],],
       edge_mean ~ age_min, col = rgb(0,0,1,0.5), # raw sightings
       pch = 19, las = 1, xlim = c(5,45), ylim = c(0,1),
       xlab = 'minimum age', ylab = 'edge weight',
       main = paste0('difference = ',age_differences[difference],' years'))
  #rethinking::shade(pi_list[[i]],
  #                  age_min[which(age_min <= ages[i])],
  #                  lwd = 2, col = rgb(0.5,0,1,0.2))      # add mean shading
  #polygon(y = plogis(sim_pi),
  #        x = c(ages,ages),
  #        lwd = 1, col = rgb(0.5,0,1,0.2))      # add mean shading
  polygon(y = c(plogis(sim_pi[1,]),plogis(sim_pi[2,])),
          x = c(rev(ages),ages),
          lwd = 1, col = rgb(0.5,0,1,0.2))       # add mean shading
  lines(x = ages,
        y = plogis(sim_mu), lwd = 2, col = 'purple')   # add mean line
  #rethinking::shade(sim_pi, age_min)             # add predictions
}




## simulate from posterior -- IS THIS RIGHT??
# sim_edges <- matrix(nrow = n_samples, ncol = length(age_min))
# for(i in 1:nrow(sim_edges)){
#   for(j in 1:ncol(sim_edges)){
#     mean_sim <- b_min$slope_draw[i]*age_min[j] + b_diff$slope_draw[i]*age_diff[j]
#     sigma_sim <- sigma$slope_draw[i]
#     sim_edges[i,j] <- MASS::mvrnorm(n = 1, mu = mean_sim, Sigma = sigma_sim)
#     #sim_edges[i,j] <- LaplacesDemon::rmvnc(n = 1, mu = mean_sim, U = sigma_sim)
#   }
# }
# hist(sim_edges)  # almost entirely equal
# sim_pi <- apply(sim_edges, 2, quant_95)

sim_edges <- array(NA, dim = c(n_samples, length(age_min), length(age_diff)))
for(i in 1:n_samples){
  for(j in 1:length(age_min)){
    for(k in 1:length(age_diff)){
      mean_sim <- b_min$slope_draw[i]*age_min[j] + b_diff$slope_draw[i]*age_diff[k]
      sigma_sim <- sigma$slope_draw[i]
      sim_edges[i,j,k] <- MASS::mvrnorm(n = 1, mu = mean_sim, Sigma = sigma_sim)
    }
  }
}
par(mfrow = c(4,2))
sim_pi <- array(NA, dim = c(2, length(age_min), length(age_diff)), dimnames = list(c('2.5%','97.5%'), NULL,NULL))
for(k in 1:length(age_diff)){
  hist(sim_edges[,,k])  # almost entirely equal
  sim_pi[,,k] <- apply(sim_edges[,,k], 2, quant_95)
}

sim_pi_df <- sim_pi[1,,1] %>% 
  as.data.frame() %>% 
  mutate(age_min = age_min,
         age_diff = age_diff[1])
colnames(sim_pi_df)[1] <- 'p2.5'
sim_pi_df$p97.5 <- sim_pi[2,,1]
for(i in 2:length(age_diff)){
  sim_df <- sim_pi[1,,i] %>% 
    as.data.frame() %>% 
    mutate(age_min = age_min,
           age_diff = age_diff[i])
  colnames(sim_df)[1] <- 'p2.5'
  sim_df$p97.5 <- sim_pi[2,,i]
  sim_pi_df <- rbind(sim_pi_df, sim_df)
  rm(sim_df)
}
sim_pi_df <- sim_pi_df %>% 
  relocate(p2.5, .before = p97.5) %>% 
  mutate(invlogit_2.5 = LaplacesDemon::invlogit(p2.5),
         invlogit_97.5 = LaplacesDemon::invlogit(p97.5))

# plot
#par(mfrow = c(length(unique(mu_means$min_age)),length(unique(mu_means$diff_age))))
par(mai = c(1,1,1,0.6))
for(i in 1:length(age_diff)){
  plot(data = mean_edges[mean_edges$age_min == age_min[i],],
       edge_mean ~ age_diff, col = rgb(0,0,1,0.5), # raw sightings
       pch = 19, las = 1, xlim = c(0,50), ylim = c(0,1),
       xlab = 'age difference', ylab = 'edge weight',
       main = paste0('younger = ',ages[i],' years old'))
  polygon(y = c(rev(sim_pi_df$invlogit_2.5[age_min == age_min[i]]),
                sim_pi_df$invlogit_97.5[age_min == age_min[i]]),
          x = c(rev(unique(sim_pi_df$age_diff)),unique(sim_pi_df$age_diff)),
          lwd = 2, col = rgb(0.5,0,1,0.2))      # add mean shading
  #polygon(x = age_diff, y = 
  #                  age_min[which(age_min <= ages[i])],
  #                  lwd = 2, col = rgb(0.5,0,1,0.2))      # add mean shading
  lines(data = mu_means[mu_means$min_age <= age_min[i],],
        LaplacesDemon::invlogit(mean) ~ diff_age, lwd = 2, col = 'purple')            # add mean line
  #rethinking::shade(sim_pi, age_min)                      # add predictions
}
par(mfrow = c(1,1))

## heat map: x = min age, y = max age, colour = edge weight
mean_predict <- matrix(NA, nrow = length(age_max), ncol = length(age_min))
for(i in 1:nrow(mean_predict)){
  for(j in 1:ncol(mean_predict)){
    if(i <= j){
      x <- mu_max_age(mu, age_max[j])
      mean_predict[,j] <- c(x ,rep(NA, length(mean_predict[,j]) - length(x)))
    }
  }
}
colnames(mean_predict) <- paste0('max.',age_max)
rownames(mean_predict) <- paste0('min.',age_min)
mean_predict_long <- as.data.frame(mean_predict) %>% 
  pivot_longer(everything(), names_to = 'age_max', values_to = 'mean_edge_weight') %>% 
  separate(age_max, into = c('max.min','age_max'), sep = 4) %>% 
  select(-max.min) %>% 
  mutate(age_min = rep(rownames(mean_predict), each = ncol(mean_predict))) %>% 
  separate(age_min, into = c('max.min','age_min'), sep = 4) %>% 
  select(-max.min) %>% 
  mutate(age_diff = as.numeric(age_max) - as.numeric(age_min))
#filter(!is.na(mean_edge_weight)) %>% 
mean_predict_long$age_max <- as.numeric(mean_predict_long$age_max)
mean_predict_long$age_min <- as.numeric(mean_predict_long$age_min)
ggplot()+
  geom_tile(data = mean_predict_long, mapping = aes(x = age_min, y = age_max, fill = mean_edge_weight))
ggplot()+
  geom_line(data = mean_predict_long, mapping = aes(x = age_diff, y = mean_edge_weight, colour = age_min, group = age_min))+
  scale_x_continuous(limits = c(0,40))

### end pdf
dev.off()

# save image
save.image('anpshort1_dyadicregression_minmax.RData')

### clear workspace
rm(cdf_1, n_dyads, counts_ls, n_eles, fit_edges_anp1, edges, plot_dyads, plot_edges, nodes, edge_weights_matrix, edge_samples)

# add time marker
print(paste0('completed at ', Sys.time()))
