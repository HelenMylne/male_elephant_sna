#### information #####
# To estimate the ages of individuals in other populations, first come up with a survival curve for the Amboseli elephants which have good birth and death data.
# survival curve options = exponential/logistic/Gompertz/Weibull, where shape = simple/bathtub/Makeham
# follow through Mia's code to determine which curve will best fit Amboseli data. she ran it on a high performance cluster, so it’s set up to run in parallel, but the basta model can also be run in parallel on its own (just like Stan). we will need to implement it as a custom function in Stan for the final analysis, given that we need to define it in the prior anyway. 

#### load packages ####
library(tidyverse)
library(cmdstanr)
library(ggdist)
library(posterior)
library(bayesplot)

#library(ggthemes)
#library(survival)
#library(tidybayes)
#library(data.table)

#### simulate process of assigning age categories to elephants ####
gompertz_bt <- function(a0, a1, c, b0, b1, age){
  gompertz <- exp(b0 + b1*age)
  bathtub <- exp(a0 - a1*age) + c + gompertz
  return(bathtub)
}

# draw new output curves
max_age <- 60
mortality <- gompertz_bt(a0 = -5.13,
                         a1 = 3.00,
                         c  = 0.026,
                         b0 = -5.08,
                         b1 = 0.09, 1:max_age)
plot(mortality)
probs <- 1 - (mortality/max(mortality))
plot(probs)

# create a fictional population with ages selected from that distribution
N <- 200
K <- 8
elephants_ls <- list(
  N = N,
  K = K,
  age = sample(1:max_age, N, replace = TRUE, prob = probs)
)
hist(elephants_ls$age)

# simulate observing ages with error and then binning them into age groups
error <- 3 # Error (in either direction) of age estimates
elephants_ls$age_guess <- elephants_ls$age + sample(-error:error, N, replace = TRUE)
elephants_ls$age_category_index <- sapply(elephants_ls$age_guess,
                                          function(x) which.max(x < c(5, 10, 15, 20, 25, 40, 60, Inf)))
hist(elephants_ls$age_category_index)

# look at actual age vs biologist assigned age: chance of being mis-classified is lower if actual age is in middle of age category, and biased towards lower end of class.
data.frame(elephants_ls) %>%
  ggplot(aes(x=age, y=age_guess, col=factor(age_category_index))) +
  geom_point(size=4,alpha=0.6) +
  geom_vline(xintercept=c(5, 10, 15, 20, 25, 40, 60), col=factor(1:7), linetype="dashed", alpha=0.6) +
  theme_minimal() + 
  xlab("Assigned age") + ylab("Age")

#### fit model to simulated dataset ####
# read in Stan model to estimate ages based on Gompertz bathtub distribution from ANP
latent_age_ordinal_model <- cmdstan_model("models/age_estimation/elephant_latent_age_ordinal_regression_hkm_22.07.07.stan")

#Fit model with cmdstanr
age_estimation_fit <- latent_age_ordinal_model$sample(
  data = elephants_ls, 
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000
)

# Examine the estimates. We can plot the estimated ages against the biologist assigned ages and ...
age_est_mat <- age_estimation_fit$summary()[202:401, ]
summary(age_est_mat)
hist(age_est_mat$mean)
hist(age_est_mat$rhat, breaks = 20)

plot_data <- data.frame(age = elephants_ls$age,       # Biologists original age est
                        model_age = age_est_mat$mean) # Mean modelled age

plot_data %>%
  ggplot(aes(x=factor(age), y=model_age)) +
  geom_point(size=4,col = 'blue', alpha=0.6) +
  geom_vline(xintercept=c(5, 10, 15, 20, 25, 40, 60), linetype="dashed", alpha=0.6) +
  geom_hline(yintercept=c(5, 10, 15, 20, 25, 40, 60), linetype="dashed", alpha=0.6) +
  geom_abline(slope = 1, intercept = 0)+
  scale_y_continuous(limits = c(0,60))+
  theme_minimal() + 
  xlab("Assigned age") + ylab("Modelled age")

#### load MOTNP data ####
motnp_males <- read_csv('data_processed/motnp_elenodes_22.01.13.csv') %>% 
  #filter(dem_class == 'AM' | dem_class == 'PM') %>% 
  filter(sex == 'M')
unique(motnp_males$age_category)
motnp_males$age_cat_id <- ifelse(motnp_males$age_category == "0-3", 1, 
                                 ifelse(motnp_males$age_category == "3-4", 1,
                                        ifelse(motnp_males$age_category == "1-2", 1,
                                               ifelse(motnp_males$age_category == "7-8", 1,
                                                      ifelse(motnp_males$age_category == "4-5", 1, 
                                                             ifelse(motnp_males$age_category ==  "6-7", 1,
                                                                    ifelse(motnp_males$age_category == "8-9", 1, 
                                                                           ifelse(motnp_males$age_category == "5-6", 1, NA))))))))
motnp_males$age_cat_id <- ifelse(motnp_males$age_category == '9-10', 2,
                                 ifelse(motnp_males$age_category == '10-15', 3,
                                        ifelse(motnp_males$age_category == '15-19', 4,
                                               ifelse(motnp_males$age_category == '20-25', 5,
                                                      ifelse(motnp_males$age_category == '25-40', 6,
                                                             ifelse(motnp_males$age_category == '40+', 7,
                                                                    motnp_males$age_cat_id))))))

#### create data list ####
N_motnp <- nrow(motnp_males)
K <- 8
motnp_ls <- list(
  N = N_motnp,
  K = K,
  age_category_index = motnp_males$age_cat_id)
hist(motnp_ls$age_category_index)
#hist(elephants_ls$age_category_index)

#### fit model to MOTNP data ####
# Fit model with cmdstanr
age_motnp_fit <- latent_age_ordinal_model$sample(
  data = motnp_ls, 
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000
)

# Examine the estimates. We can plot the estimated ages against the biologist assigned ages and ...
age_est_mat <- age_motnp_fit$summary()[(N_motnp+2):(N_motnp*2+1), ]
summary(age_est_mat)
hist(age_est_mat$mean)
hist(age_est_mat$rhat, breaks = 20)

plot_data <- data.frame(age = ifelse(motnp_ls$age == 1, 3, 
                                     ifelse(motnp_ls$age == 2, 8,
                                            ifelse(motnp_ls$age == 3, 12,
                                                   ifelse(motnp_ls$age == 4, 18,
                                                          ifelse(motnp_ls$age == 5, 22, 
                                                                 ifelse(motnp_ls$age == 6, 32, 45)))))),
                        model_age = age_est_mat$mean) # Mean modelled age

plot_data %>%
  ggplot(aes(x=factor(age), y=model_age)) +
  geom_point(size=4,col = 'blue', alpha=0.6) +
  geom_vline(xintercept = c(5, 10, 15, 20, 25, 40, 60), linetype = "dashed", alpha = 0.6) +
  geom_hline(yintercept = c(5, 10, 15, 20, 25, 40, 60), linetype = "dashed", alpha = 0.6) +
  #geom_abline(slope = 1, intercept = 0)+
  scale_y_continuous(limits = c(0,60))+
  theme_minimal() + 
  xlab("Assigned age") + ylab("Modelled age")

# posterior predictive plot using draws from distribution to show uncertainty around mean age
true_ages <- age_motnp_fit$draws("true_age", format="df")
#mcmc_dens(true_ages)
true_ages <- true_ages[,1:N_motnp]

df <- as.data.frame(do.call(rbind, true_ages)) %>%
  mutate(age_cat = motnp_ls$age) %>% relocate(age_cat) %>%
  mutate(ID = motnp_males$id) %>% relocate(ID)

df <- df %>% pivot_longer(cols = 3:102) %>% select(-name)

df$true_age <- ifelse(df$age_cat == 1, 3, 
                      ifelse(df$age_cat == 2, 8,
                             ifelse(df$age_cat == 3, 12,
                                    ifelse(df$age_cat == 4, 18,
                                           ifelse(df$age_cat == 5, 22, 
                                                  ifelse(df$age_cat == 6, 32, 45))))))

df %>% ggplot(aes(x=true_age, y=value, group=factor(ID))) +
  geom_point(size=2,col = 'blue', alpha=0.1) +
  #stat_halfeye() +
  geom_vline(xintercept=c(5,10,15,20,25,40,60), linetype="dashed", alpha=0.6) +
  geom_hline(yintercept=c(5,10,15,20,25,40,60), linetype="dashed", alpha=0.6) +
  theme_bw() + 
  xlab("Assigned age") + ylab("Modelled age")
