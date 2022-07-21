#### information #####
# To estimate the ages of individuals in other populations, first come up with a survival curve for the Amboseli elephants which have good birth and death data. -- this identified the best curve as being the Gompertz-Bathtub.
# Next fit this to the MPNP dataset

#### load packages ####
#library(tidyverse)
#library(cmdstanr)
#library(ggdist)
#library(posterior)
#library(bayesplot)

library(tidyverse, lib.loc = 'packages/')
library(cmdstanr, lib.loc = 'packages/')
#library(ggdist, lib.loc = 'packages/')
library(posterior, lib.loc = 'packages/')
#library(bayesplot, lib.loc = 'packages/')
library(janitor, lib.loc = 'packages/')
library(readxl, lib.loc = 'packages/')
library(MASS, lib.loc = 'packages/')

#### set up graph outputs ####
pdf('mpnp_age_ouputs.pdf', width = 20, height = 10)

#### load model -- adapted from MOTNP to use MPNP thresholds ####
# read in Stan model to estimate ages based on Gompertz bathtub distribution from ANP
latent_age_ordinal_model <- cmdstan_model("models/mpnp_age_ordinal_regression.stan")
# Age_Range_ID
# 1	>1
# 2	1-4
# 3	5-9
# 4	10-15
# 5	16-20
# 6	21-25
# 7	26-35
# 8	36+
#10	UK
print('model loaded')

#### simulate data with full age ranges in population -- working fine (divergent transitions = 0%) ####
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
N <- 500
K <- 9
elephants_ls <- list(
  id = 1:N,
  N = N,
  K = K,
  age = sample(1:max_age, N, replace = TRUE, prob = probs)
)
hist(elephants_ls$age)

# simulate observing ages with error and then binning them into age groups
error <- 3 # Error (in either direction) of age estimates
elephants_ls$age_guess <- elephants_ls$age + sample(-error:error, N, replace = TRUE)
elephants_ls$age_category_index <- sapply(elephants_ls$age_guess,
                                          function(x) which.max(x < c(1, 4, 9, 15, 20, 25, 35, 70)))
hist(elephants_ls$age_category_index)

# look at actual age vs biologist assigned age: chance of being mis-classified is lower if actual age is in middle of age category, and biased towards lower end of class.
data.frame(elephants_ls) %>%
  ggplot(aes(x=age, y=age_guess, col=factor(age_category_index))) +
  geom_point(size=4,alpha=0.6) +
  geom_vline(xintercept=c(1, 4, 9, 15, 20, 25, 35, 60), col=factor(1:8), linetype="dashed", alpha=0.6) +
  theme_minimal() + 
  xlab("Assigned age") + ylab("Age")

#### fit model to simulated dataset
#Fit model with cmdstanr
age_estimation_fit <- latent_age_ordinal_model$sample(
  data = elephants_ls, 
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000
)

# Examine the estimates
age_est_mat <- age_estimation_fit$summary()[(N+2):(2*N+1), ]
summary(age_est_mat)
hist(age_est_mat$mean)
hist(age_est_mat$rhat, breaks = 20)

plot_data <- data.frame(age = elephants_ls$age,       # Biologists original age est
                        model_age = age_est_mat$mean) # Mean modelled age

plot_data %>%
  ggplot(aes(x=factor(age), y=model_age)) +
  geom_point(size=4,col = 'blue', alpha=0.1) +
  geom_vline(xintercept=c(5, 10, 15, 20, 25, 40, 60), linetype="dashed", alpha=0.6) +
  geom_hline(yintercept=c(5, 10, 15, 20, 25, 40, 60), linetype="dashed", alpha=0.6) +
  geom_abline(slope = 1, intercept = 0)+
  scale_y_continuous(limits = c(0,60))+
  theme_minimal() + 
  xlab("Assigned age") + ylab("Modelled age")

# posterior predictive plot using draws from distribution to show uncertainty around mean age
sim_ages <- age_estimation_fit$draws("true_age", format="df")
sim_ages <- sim_ages[,1:N]

df <- as.data.frame(do.call(rbind, sim_ages)) %>%
  mutate(age_cat = elephants_ls$age_category_index) %>% relocate(age_cat) %>%
  mutate(ID = elephants_ls$id) %>% relocate(ID)

df <- df %>% pivot_longer(cols = 3:102) %>% select(-name)

df$true_age <- ifelse(df$age_cat == 1, 1, 
                      ifelse(df$age_cat == 2, 3,
                             ifelse(df$age_cat == 3, 7,
                                    ifelse(df$age_cat == 4, 12,
                                           ifelse(df$age_cat == 5, 18, 
                                                  ifelse(df$age_cat == 6, 22, 
                                                         ifelse(df$age_cat == 7, 30, 40)))))))

df %>% ggplot(aes(x=true_age, y=value, group=factor(ID))) +
  geom_point(size=2,col = 'blue', alpha=0.1) +
  geom_vline(xintercept=c(1,4,9,15,20,25,35,60), linetype="dashed", alpha=0.6) +
  geom_hline(yintercept=c(1,4,9,15,20,25,35,60), linetype="dashed", alpha=0.6) +
  theme_bw() + 
  scale_x_continuous(breaks = c(1,4,9,15,20,35,60)) +
  xlab("Assigned age") + ylab("Modelled age")

print('initial simulation with complete age range complete')

#### repeat but now simulating absence of all individuals younger than category 3, and most from category 3 -- working almost fine (divergent transitions = 4%) ####
# create a fictional population with ages selected from that distribution
N2 <- 1000
age <- sample(1:max_age, N2, replace = TRUE, prob = probs)
age <- age[age > 8]
hist(age)
elephants_ls <- list(
  id = 1:N,
  N = N,
  K = K,
  age = sample(age, N, replace = FALSE)
)

hist(elephants_ls$age)

# simulate observing ages with error and then binning them into age groups
elephants_ls$age_guess <- elephants_ls$age + sample(-error:error, N, replace = TRUE)
elephants_ls$age_category_index <- sapply(elephants_ls$age_guess,
                                          function(x) which.max(x < c(1, 4, 9, 15, 20, 25, 35, 70)))
hist(elephants_ls$age_category_index)

# look at actual age vs biologist assigned age: chance of being mis-classified is lower if actual age is in middle of age category, and biased towards lower end of class.
data.frame(elephants_ls) %>%
  ggplot(aes(x=age, y=age_guess, col=factor(age_category_index))) +
  geom_point(size=4,alpha=0.6) +
  geom_vline(xintercept=c(1, 4, 9, 15, 20, 25, 35, 60), col=factor(1:8), linetype="dashed", alpha=0.6) +
  theme_minimal() + 
  xlab("Assigned age") + ylab("Age")

#### fit model to simulated dataset
#Fit model with cmdstanr
age_estimation_fit <- latent_age_ordinal_model$sample(
  data = elephants_ls, 
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000
)

# Examine the estimates
age_est_mat <- age_estimation_fit$summary()[(N+2):(2*N+1), ]
summary(age_est_mat)
hist(age_est_mat$mean)
hist(age_est_mat$rhat, breaks = 20)

plot_data <- data.frame(age = elephants_ls$age,       # Biologists original age est
                        model_age = age_est_mat$mean) # Mean modelled age

plot_data %>%
  ggplot(aes(x=factor(age), y=model_age)) +
  geom_point(size=4,col = 'blue', alpha=0.1) +
  #geom_vline(xintercept=c(5, 10, 15, 20, 25, 40, 60), linetype="dashed", alpha=0.3) +
  geom_hline(yintercept=c(5, 10, 15, 20, 25, 40, 60), alpha=0.1) +
  geom_abline(slope = 1, intercept = 8)+
  scale_y_continuous(limits = c(0,60))+
  theme_minimal() + 
  xlab("Assigned age") + ylab("Modelled age")

# posterior predictive plot using draws from distribution to show uncertainty around mean age
sim_ages <- age_estimation_fit$draws("true_age", format="df")
sim_ages <- sim_ages[,1:N]

df <- as.data.frame(do.call(rbind, sim_ages)) %>%
  mutate(age_cat = elephants_ls$age_category_index) %>% relocate(age_cat) %>%
  mutate(ID = elephants_ls$id) %>% relocate(ID)

df <- df %>% pivot_longer(cols = 3:102) %>% select(-name)

df$true_age <- ifelse(df$age_cat == 1, 1, 
                      ifelse(df$age_cat == 2, 3,
                             ifelse(df$age_cat == 3, 7,
                                    ifelse(df$age_cat == 4, 12,
                                           ifelse(df$age_cat == 5, 18, 
                                                  ifelse(df$age_cat == 6, 22, 
                                                         ifelse(df$age_cat == 7, 30, 40)))))))

df %>% ggplot(aes(x=true_age, y=value, group=factor(ID))) +
  geom_point(size=2,col = 'blue', alpha=0.1) +
  geom_vline(xintercept=c(1,4,9,15,20,25,35,60), linetype="dashed", alpha=0.6) +
  geom_hline(yintercept=c(1,4,9,15,20,25,35,60), linetype="dashed", alpha=0.6) +
  theme_bw() + 
  scale_x_continuous(breaks = c(1,4,9,15,20,35,60)) +
  xlab("Assigned age") + ylab("Modelled age")

print('simulation without young males complete')

#### load MPNP data ####
mpnp_long <- readxl::read_excel('../data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211019.xlsx') %>%
  janitor::clean_names()

# identify period 1 only
periods <- seq(from = min(mpnp_long$date, na.rm = T),
               to = max(mpnp_long$date, na.rm = T),
               length.out = 7)
periods[7] <- periods[7]+1 # set to be one higher than the final higher otherwise it takes the last date and creates a whole new period

# remove non-B-numbered elephants
mpnp_long <- mpnp_long %>% 
  filter(elephant_id != '-') %>% 
  separate(elephant_id, into = c('BTF','num'), sep = 1, remove = F) %>% 
  filter(BTF == 'B') %>% 
  select(-BTF, -num)

# select only period 1
mpnp1_long <- mpnp_long[mpnp_long$date < periods[2],]
max(mpnp_long$date, na.rm = T) - min(mpnp_long$date, na.rm = T)

# get range of ages
mpnp1_long$age_range_id <- as.numeric(mpnp1_long$age_range_id)
mpnp1_long$age_mid <- NA
mpnp1_long$age_range <- NA
for(i in 1:nrow(mpnp1_long)){
  ele <- mpnp1_long[mpnp1_long$elephant_id == mpnp1_long$elephant_id[i],]
  ele <- ele[!is.na(ele$age_range_id),]
  ele <- ele[ele$age_range_id != 10,]
  mpnp1_long$age_mid[i] <- ifelse(nrow(ele) == 0, 10, median(ele$age_range_id, na.rm = T))
  mpnp1_long$age_range[i] <- ifelse(nrow(ele) == 0, 0, max(ele$age_range_id, na.rm = T) - min(ele$age_range_id, na.rm = T))
  rm(ele)
}

# take median category and then round down to nearest integer
mpnp1_long$age_mid_round <- floor(mpnp1_long$age_mid)
mpnp1_long <- mpnp1_long[mpnp1_long$age_mid_round < 10 , c(3,5,6,7,32:34)]

rm(mpnp_long, periods)

#### create MPNP data list ####
mpnp1_males <- mpnp1_long %>%
  select(elephant_id, age_mid_round) %>% 
  distinct()
N_mpnp1 <- nrow(mpnp1_males)
mpnp1_ls <- list(
  N = N_mpnp1,
  K = K,
  age_category_index = mpnp1_males$age_mid_round)
hist(mpnp1_ls$age_category_index)

print('model data prepped')

#### fit model to MPNP data -- divergent transitions = 3% ####
# Fit model with cmdstanr
age_mpnp1_fit <- latent_age_ordinal_model$sample(
  data = mpnp1_ls, 
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000
)
print ('model fitted')

# Examine the estimates. We can plot the estimated ages against the biologist assigned ages
age_est_mat <- age_mpnp1_fit$summary()[(N_mpnp1+2):(N_mpnp1*2+1), ]
summary(age_est_mat)
hist(age_est_mat$mean)
hist(age_est_mat$rhat, breaks = 20)

plot_data <- data.frame(age = ifelse(mpnp1_ls$age == 1, 1, 
                                     ifelse(mpnp1_ls$age == 2, 3,
                                            ifelse(mpnp1_ls$age == 3, 7,
                                                   ifelse(mpnp1_ls$age == 4, 12,
                                                          ifelse(mpnp1_ls$age == 5, 18,
                                                                 ifelse(mpnp1_ls$age == 6, 22, 
                                                                        ifelse(mpnp1_ls$age == 7, 30, 45))))))),
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
true_ages <- age_mpnp1_fit$draws("true_age", format="df")
true_ages <- true_ages[,1:N_mpnp1]

df <- as.data.frame(do.call(rbind, true_ages)) %>%
  mutate(age_cat = mpnp1_ls$age) %>% relocate(age_cat) %>%
  mutate(ID = mpnp1_males$elephant_id) %>% relocate(ID)

df <- df %>% pivot_longer(cols = 3:102) %>% select(-name)

df$true_age <- ifelse(df$age_cat == 1, 1, 
                      ifelse(df$age_cat == 2, 3,
                             ifelse(df$age_cat == 3, 7,
                                    ifelse(df$age_cat == 4, 12,
                                           ifelse(df$age_cat == 5, 18, 
                                                  ifelse(df$age_cat == 6, 22, 
                                                         ifelse(df$age_cat == 7, 30, 40)))))))

df %>% ggplot(aes(x=true_age, y=value, group=factor(ID))) +
  geom_point(size=2,col = 'blue', alpha=0.1) +
  #stat_halfeye() +
  geom_vline(xintercept=c(1,4,9,15,20,25,35,60), linetype="dashed", alpha=0.6) +
  geom_hline(yintercept=c(1,4,9,15,20,25,35,60), linetype="dashed", alpha=0.6) +
  theme_bw() + 
  scale_x_continuous(breaks = c(1,4,9,15,20,35,60)) +
  xlab("Assigned age") + ylab("Modelled age")

#### save output ####
saveRDS(true_ages, '../data_processed/mpnp1_bayesian_agedistribution.rds')
