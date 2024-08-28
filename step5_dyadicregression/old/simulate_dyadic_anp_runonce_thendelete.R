#### set up ####
# script to simulate dyadic regression and check actually working
# library(StanHeaders) ; library(rstan) ; library(tidyverse) ; library(car) ; library(LaplacesDemon) ; library(patchwork)
library(StanHeaders, lib.loc = '../packages/')    # library(StanHeaders)
library(rstan, lib.loc = '../packages/')          # library(rstan)
library(tidyverse, lib.loc = '../packages/')      # library(tidyverse)
library(car, lib.loc = '../packages/')            # library(car)
library(LaplacesDemon, lib.loc = '../packages/')  # library(LaplacesDemon)
library(patchwork, lib.loc = '../packages/')  # library(LaplacesDemon)

theme_set(theme_bw())

pdf('step5_dyadicregression/simulate_dyadic_anp_allplots.pdf')

set.seed(12345)

print('ready to start')

#### load workspace ####
load('step5_dyadicregression/anp_dyadic_simulation.RData')

#### simulate edges ####
## set mean edge weight
ggplot(edge_summary)+
  geom_point(aes(x = age_max, y = mu, colour = age_min))

## standardise edges
ggplot(edge_summary)+
  geom_point(aes(x = age_max, y = mu_std, colour = age_min))

## simulate full distribution of samples per node
plot(sim_dat[1,] ~ edge_summary$age_min)          # plot simulated values against minimum age
plot(sim_dat[1,] ~ edge_summary$age_max)          # plot simulated values against maximum age

## print progress marker
print('edges simulated')

#### plot raw data ####
ggplot()+
  geom_point(data = edge_summary, aes(x = age_min_std,
                             y = mu_std,
                             colour = age_cat_max),
             shape = 19)+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  scale_colour_viridis_d()+
  labs(colour = 'age category of\nolder male')

ggplot()+
  geom_boxplot(data = edge_summary, aes(x = age_cat_min,
                               y = mu_std,
                               fill = age_cat_max),
               shape = 19)+
  scale_x_discrete('age category of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  theme(legend.position = 'bottom')+
  scale_fill_viridis_d()+
  labs(fill = 'age category of older male')

sim_dat_long <- sim_dat %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = 'dyad_window', values_to = 'edge_draw') %>%
  left_join(edge_summary, by = 'dyad_window')
ggplot()+
  geom_violin(data = sim_dat_long,
              aes(x = as.factor(age_cat_min),
                  y = edge_draw,
                  fill = as.factor(age_cat_max)),
              #shape = 19,
              alpha = 0.5)+
  geom_point(data = edge_summary,
             aes(x = as.factor(age_cat_min), y = mu_std,
                 group = as.factor(age_cat_max)),
             shape = 19,
             #colour = 'white',
             size = 1)+
  scale_x_discrete('age category of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  scale_fill_viridis_d()+
  theme(legend.position = 'bottom')+
  labs(fill = 'age category of older elephant')

## print progress marker
print('edges plotted')

#### fit multivariate Gaussian distribution to output of edge weight model ####
## set up plotting
num_check <- 64
selected_samples <- sample(1:n_dyads, num_check, replace = FALSE)

### Setting grid layout
rows <- floor(sqrt(num_check))
cols <- ceiling(num_check / rows)
par(mfrow=c(rows, cols), mar=c(2,2,2,1))

### plot
for (i in selected_samples) {
  mu <- logit_edge_draws_mu[i]
  sd <- logit_edge_draws_sd[i] #sqrt(logit_edge_draws_cov[i,i])

  fitted_values <- rnorm(1e5, mean = mu, sd = sd)

  hist(unlist(sim_dat[,i]), probability = TRUE, las = 1,
       main = paste("Dyad", i), xlab = "Value", breaks = 50)
  lines(density(fitted_values), col="red", lw=1.5)
}

for (i in selected_samples) {
  plot(unlist(sim_dat[,i]), unlist(sim_dat[,i+1]),
       col = rgb(0,0,1,0.05), las = 1,
       main = paste("cov ", i ,"&",i+1))
}

## reset plot window and clean up
par(mfrow=c(1,1), mai = c(1,1,0.5,0.5))
rm(cols, fitted_values, i, j, mu, num_check, rows, sd, selected_samples) ; gc()

## print progress marker
print('normal approximation complete')

#### prior predictive check ####
n <- 100
beta_age_min <- rnorm(n, 0, 0.8)
beta_age_max <- rnorm(n, 0, 0.8)
intercept <- rnorm(n, 0, 2)
plot(NULL, las = 1, xlab = 'minimum age', ylab = 'logit edge weight (standardised)',
     ylim = c(min(logit_edge_draws_mu)-5, max(logit_edge_draws_mu)+5),
     xlim = c(min(edge_summary$age_min_std), max(edge_summary$age_max_std)))
abline(h = min(logit_edge_draws_mu), lty = 2)
abline(h = max(logit_edge_draws_mu), lty = 2)
x <- c(min(edge_summary$age_min_std), max(edge_summary$age_max_std))
for(i in 1:n){
  y <- rep(NA, length(x))
  for(j in 1:length(x)){
    y[j] <- intercept[i] + beta_age_min[i]*x[j] + beta_age_max[i]*x[j]
  }
  lines(x = x, y = y, col = rgb(0,0,1,0.4))
}
rm(n, beta_age_max, beta_age_min, intercept, i, y, j, x) ; gc()

## print progress marker
print('prior predictive complete')

#### check outputs ####
## extract model fit
summary <- rstan::summary(fit_dyadreg_sim)
(summary <- as.data.frame(summary$summary))
par(mfrow = c(3,1))
hist(summary$Rhat, breaks = 50)
hist(summary$n_eff, breaks = 50)
par(mfrow = c(1,1))

## extract draws
draws <- rstan::extract(fit_dyadreg_sim)

## extract dyadic regression parameters
b_max <- draws$beta_age_max
b_min <- draws$beta_age_min
intercept <- draws$intercept
sigma <- draws$sigma
parameters <- data.frame(beta_age_max = b_max,
                         beta_age_min = b_min,
                         intercept = intercept,
                         sigma = sigma) %>%
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>%
  pivot_longer(cols = c('beta_age_max','beta_age_min','sigma','intercept'),
               names_to = 'parameter', values_to = 'slope_draw')

## traceplot function
traceplot <- function(parameter_data){
  plot <- ggplot(data = parameter_data)+
    geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
    theme(legend.position = 'none')+
    scale_colour_viridis_d()+
    facet_wrap(. ~ parameter , scales = 'free_y')
  return(plot)
}

## global parameter traceplots
traceplot(parameter_data = parameters)

print('global parameters plotted')

## random effect extraction function
long_random <- function(random_effect, n){
  random_effect_long <- random_effect %>%
    mutate(chain = rep(1:4, each = 1000),
           position = rep(1:1000, 4)) %>%
    pivot_longer(cols = all_of(1:n),
                 names_to = 'parameter',
                 values_to = 'slope_draw')
  return(random_effect_long)
}

## extract multimembership draws
mm_nodes <- as.data.frame(draws$mm_nodes)
colnames(mm_nodes) <- unique(c(dyad_data$node_1, dyad_data$node_2)) # I don't know that this is right, but I can't see any other order it would be in?
mm_nodes_long <- long_random(mm_nodes, n_nodes)

## multimembership traceplots
mm_nodes_long %>%
  filter(parameter %in% sample(unique(mm_nodes_long$parameter), 25, replace = F)) %>%
  traceplot()

print('multimembership parameters plotted')

## extract window random effect draws
rand_window <- as.data.frame(draws$rand_window)
colnames(rand_window) <- unique(dyad_data$window)
rand_window_long <- long_random(rand_window, n_windows)

## window random effect traceplots
traceplot(rand_window_long)

print('window parameters plotted')

## extract dyad random effect draws
rand_dyad <- as.data.frame(draws$rand_dyad)
colnames(rand_dyad) <- unique(dyad_data$dyad)
rand_dyad_long <- long_random(rand_dyad, n_dyads)

## dyad random effect traceplots
rand_dyad_long %>%
  filter(parameter %in% sample(unique(rand_dyad_long$parameter), 25, replace = F)) %>%
  traceplot()

## print progress marker
print('outputs checked')

#### posterior predictive check ####
plot(density(as.numeric(sim_dat[1, ])),
     main = "Posterior predictive density of edge weights:\nblack = measured edge, red = predicted edge",
     ylim = c(0, 1), col = rgb(0, 0, 0, 0.25), las = 1)
for (i in 1:100) {
  j <- sample(1:nrow(sim_dat), 1)

  mu_plot <- rep(NA, n_dyads)
  for(k in 1:n_dyads){
    mu_plot[k] <- intercept[j] + b_min[j]*dyad_data$age_min[k] + b_max[j]*dyad_data$age_max[k] + mm_nodes[j,dyad_data$node_1[k]] + mm_nodes[j,dyad_data$node_2[k]]
  }

  sigma_plot <- dyad_data$logit_edge_sd + diag(rep(sigma[j], n_data))
  norm <- rnorm(1000, mu_plot, sigma_plot)

  lines(density(as.numeric(sim_dat[j, ])), col = rgb(0, 0, 0, 0.25)) # black lines for edge samples
  lines(density(norm), col = rgb(1, 0, 0, 0.25))                  # red lines for predictions

}

print('posterior predictive check complete')

#### plot predictions ####
ggplot()+
  geom_ribbon(data = edge_summary,
              aes(x = age_min,
                  ymin = pred_full_lwr, ymax = pred_full_upr,
                  group = as.factor(age_max), fill = as.factor(age_max)),
              alpha = 0.3)+
  geom_line(data = edge_summary,
            aes(x = age_min,
                y = pred_mu,
                colour = as.factor(age_max), group = as.factor(age_max)),
            linewidth = 1)+
  scale_colour_viridis_d(direction = -1)+ scale_fill_viridis_d(direction = -1)+
  geom_point(data = edge_summary,
             aes(x = age_min,
                 y = mu_std,
                 colour = as.factor(age_max))
             )+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  theme(axis.text = element_text(size = 14),
        legend.position = 'bottom', #c(0.8,0.2),
        # legend.title = element_text(size = 14),
        # legend.text = element_text(size = 12)
        axis.title = element_text(size = 18))+
  labs(colour = 'maximum age', fill = 'maximum age')+
  guides(colour = guide_legend(ncol = 6),
         fill = guide_legend(ncol = 6))

ggplot()+
  geom_point(data = edge_summary,
             mapping = aes(x = mu_std, y = pred_mu))+
  geom_abline(slope = 1, intercept = 0)

ggplot()+
  geom_point(data = edge_summary,
             mapping = aes(x = mu, y = pred_mu_ustd))+
  geom_abline(slope = 1, intercept = 0)

ggplot()+
  geom_point(data = edge_summary,
             mapping = aes(x = invlogit(mu), y = pred_mu_invlogit))+
  geom_abline(slope = 1, intercept = 0)


## print progress marker
print('predictions complete')
dev.off()

#### extract original parameters ####
## calculate contrasts: minimum age increases by 1, standardised scale
contrast_min <- pred_new_min_full - pred_org_min_full
print(paste0('contrast for age_min + 1 (std): ',
             round(mean(contrast_min), 5), ' ± ', round(sd(contrast_min), 5),
             ' [', round(quantile(contrast_min, prob = 0.025), 5),
             ' - ',round(quantile(contrast_min, prob = 0.975), 5), ']'))

## calculate contrasts: minimum age increases by 1, outcome scale
contrast_min_ustd <- pred_new_min_full_ustd - pred_org_min_full_ustd
print(paste0('contrast for age_min + 1 (ustd): ',
             round(mean(contrast_min_ustd), 5), ' ± ', round(sd(contrast_min_ustd), 5),
             ' [', round(quantile(contrast_min_ustd, prob = 0.025), 5),
             ' - ',round(quantile(contrast_min_ustd, prob = 0.975), 5), ']'))
print(paste0('true minimum effect: ', sim_min))

## calculate contrasts: maximum age increases by 1, standardised scale
contrast_max <- pred_new_max_full - pred_full
print(paste0('contrast for age_max + 1 (std): ',
             round(mean(contrast_max), 5), ' ± ', round(sd(contrast_max), 5),
             ' [', round(quantile(contrast_max, prob = 0.025), 5),
             ' - ',round(quantile(contrast_max, prob = 0.975), 5), ']'))

## calculate contrasts: maximum age increases by 1, outcome scale
contrast_max_ustd <- pred_new_max_full_ustd - pred_full_ustd
print(paste0('contrast for age_max + 1 (ustd): ',
             round(mean(contrast_max_ustd), 5), ' ± ', round(sd(contrast_max_ustd), 5),
             ' [', round(quantile(contrast_max_ustd, prob = 0.025), 5),
             ' - ',round(quantile(contrast_max_ustd, prob = 0.975), 5), ']'))
print(paste0('true minimum effect: ', sim_max))

## calculate contrasts: both ages increase by 1, standardised scale
contrast_both <- pred_new_both_full - pred_full
print(paste0('contrast for both ages + 1 (std): ',
             round(mean(contrast_both), 5), ' ± ', round(sd(contrast_both), 5),
             ' [', round(quantile(contrast_both, prob = 0.025), 5),
             ' - ',round(quantile(contrast_both, prob = 0.975), 5), ']'))

## calculate contrasts: both ages increase by 1, outcome scale
contrast_both_ustd <- pred_new_both_full_ustd - pred_full_ustd
print(paste0('contrast for both ages + 1 (ustd): ',
             round(mean(contrast_both_ustd), 5), ' ± ', round(sd(contrast_both_ustd), 5),
             ' [', round(quantile(contrast_both_ustd, prob = 0.025), 5),
             ' - ',round(quantile(contrast_both_ustd, prob = 0.975), 5), ']'))
