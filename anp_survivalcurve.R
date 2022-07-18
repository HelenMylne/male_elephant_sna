### Information ####
# Data collected by Amboseli Trust for Elephants (ATE) 1972-2021
# Data supplied by Vicki Fishlock, 24th February 2022
# model to assess the probability distribution of male age using the age category assigned and the survival curve for elephants in the Amboseli National Park population

#I've made good progress with the age estimation and it looks like it will work. I've decided to use a Weibull distribution for the survival because it's already in Stan, plus it allows higher infant mortality.

#Fit a Weibull model to the elephant data where you have ages. The way you need to do this is to have the data as age reached by the elephant, and then an event column that is 1 if it died at that age (known max age for that elephant) or 0 if it is still living (it lives at *least* to that age).
#Here is an example of how to do this in brms: https://statwonk.com/bayesian-right-censored-weibull-model.html
#Also: https://rdrr.io/cran/brms/man/kidney.html

### Set up ####
library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(survival)
library(brms)
library(tidybayes)
library(ggthemes)

### Import nodes data ####
all_nodes <- readxl::read_excel('data_raw/Raw_ATE_AllElephants_Lee220118.xlsx')
boxplot(as.numeric(all_nodes$Age_death) ~ all_nodes$Death_accuracy) # adds NA values for anything classed as "living"
str(all_nodes)

### calculate oldest age each observed ####
# create dataframe that only considers node ages
ages <- all_nodes[,c(1:3,15,16,18:20)] %>% 
  janitor::clean_names()
ages$id_no <- ifelse(ages$sex == 'unknown', paste('U', ages$id, sep = ''),
                     ifelse(ages$sex == 'Female', paste('F', ages$id, sep = ''),
                            ifelse(ages$sex == 'Male', paste('M', ages$id, sep = ''),
                                   'check')))

# generate numeric variable for age in 2021 of living elephants or age at death
ages$age_death_num <- as.numeric(ifelse(ages$age_death == 'living', 2021 - ages$birth_year, ages$age_death))
ages$age_death_round <- round(ages$age_death_num, 0)
table(ages$age_death_round, ages$sex)
hist(ages$age_death_round, breaks = 75)                               # overall age distribution
hist(ages$age_death_round[which(ages$sex == 'Male')], breaks = 70)    # male age distribution
hist(ages$age_death_round[which(ages$sex == 'Female')], breaks = 75)  # female age distribution
ages <- ages[,c(9,1,3,8,10)]

# create new 
ages$censor <- ifelse(ages$age_death == 'living', 'TRUE', 'FALSE')    # censor = TRUE when elephant is still alive

# clean up
age_cens <- ages[!is.na(ages$age_death_num),c('age_death_num','censor','sex')] # single data frame for age at death/2021
rm(ages, all_nodes)

### Weibull example - Statwonk (https://statwonk.com/bayesian-right-censored-weibull-model.html) ####
# example from Statwonk
rweibull_cens <- function(n, shape, scale) {
  a_random_death_time <- rweibull(n, shape = shape, scale = scale) 
  a_random_censor_time <- rweibull(n, shape = shape, scale = scale)
  observed_time <- pmin(a_random_censor_time, a_random_death_time)
  censor <- observed_time == a_random_death_time
  tibble(time = observed_time, censor = censor)
}

rweibull_cens(1e4, shape = 1, scale = 10) -> d
d %>%
  mutate(censor = if_else(censor == 0, 1, 0)) %>%
  brm(time | cens(censor) ~ 1, data = ., family = "weibull", 
      refresh = 2e3, cores = 4) -> bfit
print(bfit, digits = 3)
as_sample_tibble(bfit) %>%
  mutate_at(vars(b_Intercept), exp) %>%
  ggplot(aes(x = b_Intercept, y = shape)) +
  geom_density_2d(colour = "red", size = 1) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 1, colour = "blue") +
  geom_vline(xintercept = 10, colour = "blue") +
  xlab("Scale") +
  ylab("Shape") +
  ggtitle("Right-censored weibull(shape, scale)") +
  theme_fivethirtyeight()

### adapt to ATE data ####
rethinking::dens(age_cens$age_death_num[which(age_cens$sex == 'Male')], col = 'red')
rethinking::dens(age_cens$age_death_num[which(age_cens$sex == 'Female')], add = T, col = 'blue')
rethinking::dens(age_cens$age_death_num[which(age_cens$sex == 'unknown')], add = T, col = 'purple', lty = 3)
abline(v = 0)

colnames(age_cens)[1] <- 'age'
age_cens$age_non0 <- age_cens$age+0.01

age_cens %>%
  mutate(censor = if_else(censor == 0, 1, 0)) %>%
  brm(age_non0 | cens(censor) ~ 1, data = ., family = "weibull", 
      refresh = 2e3, cores = 4) -> bfit_all
print(bfit_all)
as_sample_tibble(bfit_all) %>%
  mutate_at(vars(b_Intercept), exp) %>%
  ggplot(aes(x = b_Intercept, y = shape)) +
  geom_density_2d(colour = "red", size = 1) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 0.89, colour = "blue") +
  geom_vline(xintercept = 16.2, colour = "blue") +
  xlab("Scale") +
  ylab("Shape") +
  ggtitle("Right-censored weibull all elephants (shape, scale)") +
  theme_fivethirtyeight()

age_cens[age_cens$sex == 'Male',] %>%
  mutate(censor = if_else(censor == 0, 1, 0)) %>%
  brm(age_non0 | cens(censor) ~ 1, data = ., family = "weibull", 
      refresh = 2e3, cores = 4) -> bfit_m
print(bfit_m)
as_sample_tibble(bfit_m) %>%
  mutate_at(vars(b_Intercept), exp) %>%
  ggplot(aes(x = b_Intercept, y = shape)) +
  geom_density_2d(colour = "red", size = 1) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 0.87, colour = "blue") +
  geom_vline(xintercept = 14.7, colour = "blue") +
  xlab("Scale") +
  ylab("Shape") +
  ggtitle("Right-censored weibull males (shape, scale)") +
  theme_fivethirtyeight()

age_cens[age_cens$sex == 'Female',] %>%
  mutate(censor = if_else(censor == 0, 1, 0)) %>%
  brm(age_non0 | cens(censor) ~ 1, data = ., family = "weibull", 
      refresh = 2e3, cores = 4) -> bfit_f
print(bfit_f)
as_sample_tibble(bfit_f) %>%
  mutate_at(vars(b_Intercept), exp) %>%
  ggplot(aes(x = b_Intercept, y = shape)) +
  geom_density_2d(colour = "red", size = 1) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 0.96, colour = "blue") +
  geom_vline(xintercept = 18, colour = "blue") +
  xlab("Scale") +
  ylab("Shape") +
  ggtitle("Right-censored weibull females (shape, scale)") +
  theme_fivethirtyeight()

# extract estimates for priors
bfit_all$fit  # shape = 0.89 ± 0.01. lwr = 0.87, upr = 0.92. Scale from graph = 16.2 (15.7 - 16.7)
bfit_m$fit    # shape = 0.87 ± 0.02. lwr = 0.84, upr = 0.91. Scale from graph = 14.7 (13.9 - 15.4)
bfit_f$fit    # shape = 0.96 ± 0.02. lwr = 0.92, upr = 1.    Scale from graph = 18.0 (17.1 - 18.9)

### Weibull example - DWF, 18th March 2022 ####
# simulate an age distribution based on a Weibull distribution where we set the shape and scale parameters.
probs <- dweibull(1:100, shape = 1.2, scale = 30)
probs <- (probs/sum(probs))
plot(probs) 

# create a fictional population with ages selected from that distribution
N <- 100 # Number of individuals
K <- 6   # Number of age bin categories
elephants_ls_example <- list(
  N = N,
  K = K,
  age = sample(1:100, N, replace=TRUE, prob=probs)) # Simulated ages drawn from the Weibull distribution plus an error in age estimation

# simulate observing ages with error and then binning them into age groups
E <- 3 # Error (in either direction) of age estimates
elephants_ls_example$age_guess <- elephants_ls_example$age + sample(-E:E, N, replace=TRUE)
elephants_ls_example$age_category_index <- sapply(elephants_ls_example$age_guess, function(x) which.max(x < c(15, 30, 45, 60, Inf)))
hist(elephants_ls_example$age_category)

# look at the actual age verses the biologist assigned age and look at the thresholds.
data.frame(elephants_ls_example) %>%
  ggplot(aes(x=age, y=age_guess, col=factor(age_category_index))) +
  geom_point(size=4,alpha=0.6) +
  geom_vline(xintercept=c(15,30,45,60,75), col=factor(2:6), linetype="dashed", alpha=0.6) +
  theme_minimal() + 
  xlab("Assigned age") + ylab("Age") # chance of being mis-classified is lower if an actual age is in the middle of the age category. Also, there is mostly a bias towards ages within the age class actually being towards the lower end of the class.
# fit Stan model to estimate ages from the categories based on a Weibull distribution (for real data these will be centred on parameters from fitting to age estimates based on another elephant population)
latent_age_ordinal_model <- cmdstan_model("models/elephant_latent_age_ordinal_regression_dwf_22.03.18.stan")
age_estimation_fit_example <- latent_age_ordinal_model$sample(
  data = elephants_ls_example, 
  chains = 4, 
  parallel_chains = 1,
  iter_sampling = 2000
)

# Examine the estimates.
age_est_mat_example <- age_estimation_fit_example$summary()[102:201,]  # true ages
plot_data <- data.frame(age = elephants_ls_example$age,       # Biologists original age est
                        model_age = age_est_mat_example$mean) # Mean modelled age
plot_data %>%
  ggplot(aes(x = factor(age), y = model_age)) +
  geom_point(size = 4, alpha = 0.6) +
  geom_vline(xintercept = c(15,30,45,60,75), linetype = "dashed", alpha = 0.6) +
  theme_minimal() + 
  xlab("Assigned age") + ylab("Modelled age") # plot the estimated ages against the biologist assigned ages

# Now an equivalent posterior predictive plot from draws from the distribution to show uncertainty around the mean age (then need to do it with full uncertainty).
true_ages <- age_estimation_fit$draws("true_age", format="df")
true_ages <- true_ages[1:100,1:N]
df <- as.data.frame(do.call(rbind, true_ages)) %>%
  mutate(true_age = elephants_ls$age) %>% relocate(true_age) %>%
  mutate(ID = 1:nrow(true_ages)) %>% relocate(ID)
df <- df %>% pivot_longer(cols = 3:102) %>% select(-name)
df %>% ggplot(aes(x=true_age, y=value, group=factor(ID))) +
  geom_point(size=2,alpha=0.1) +
  geom_vline(xintercept=c(15,30,45,60,75), linetype="dashed", alpha=0.6) +
  geom_hline(yintercept=c(15,30,45,60,75), linetype="dashed", alpha=0.6) +
  theme_bw() + 
  xlab("Assigned age") + ylab("Modelled age")

### estimate Weibull distribution -- DWF 28th March 2022 ####
colnames(age_cens)[1] <- 'age'
age_cens$age_non0 <- age_cens$age+0.01 # wouldn't run when I allowed age = 0, so increased all values by 0.01 to make it run

age_cens$censor <- as.integer(as.logical(age_cens$censor))  # 1647 elephants LIVING = 1647 elephants need to be censor=1

# Examine default brms priors
# prior_summary(bfit_m)

# run model for male elephants with slightly better priors
bfit_m <- age_cens %>% filter(sex == 'Male') %>%
  brm(age_non0 | cens(censor) ~ 1,
      prior = c(
        prior(student_t(3, 5, 5), class = Intercept),
        prior(exponential(0.5), class = shape)
      ),
      data = ., family = "weibull", 
      chains = 4, cores = 4)

plot(bfit_m)

# extract scale and shape from posterior draws
# to convert to scale we need to both undo the link function by taking the exponent
# and then refer to the brms documentation to understand how the mean relates to the scale
bfit_m_draws <- as_draws_df(bfit_m) %>%
  mutate(scale = exp(b_Intercept) / (gamma(1 + 1 / shape)))

est_shape <- mean(bfit_m_draws$shape)
est_scale <- mean(bfit_m_draws$scale)

# Let's do a posterior predictive plot based on the mean
probs <- dweibull(1:100, shape = est_shape, scale = est_scale)
probs <- probs / sum(probs)
plot(probs) 

# Lets simulate ages for 3k elephants
hist(rweibull(3000,est_shape,est_scale))

# Shape = 1.29, scale = 29.5 (so we will use shape = 1.3 and scale = 30)

### estimate Weibull distribution -- HKM 29th March 2022 ####
colnames(age_cens)[1] <- 'age'
age_cens$age_non0 <- ifelse(age_cens$age == 0, 0.001, age_cens$age) # wouldn't run when I allowed age = 0, so increased all values by 0.01 to make it run

age_cens$censor <- as.integer(!as.logical(age_cens$censor)) # 1 = FALSE = DEAD

# Examine default brms priors
# prior_summary(bfit_m)

# run model for male elephants with slightly better priors
bfit_m <- age_cens %>% filter(sex == 'Male') %>%
  brm(age_non0 | cens(censor) ~ 1,
      prior = c(
        prior(student_t(3,5,5), class = Intercept),
        prior(exponential(0.5), class = shape)
      ),
      data = ., family = "weibull", 
      chains = 4, cores = 4)
plot(bfit_m)

# extract scale and shape from posterior draws
# to convert to scale we need to both undo the link function by taking the exponent
# and then refer to the brms documentation to understand how the mean relates to the scale
bfit_m_draws <- as_draws_df(bfit_m) %>%
  mutate(scale = exp(b_Intercept) / (gamma(1 + 1 / shape)))

est_shape <- mean(bfit_m_draws$shape)
est_scale <- mean(bfit_m_draws$scale)

# Let's do a posterior predictive plot based on the mean
probs <- dweibull(1:100, shape = est_shape, scale = est_scale)
probs <- probs / sum(probs)
plot(probs) 

# Lets simulate ages for 3k elephants
hist(rweibull(3000,est_shape,est_scale))

# Shape = 1.29, scale = 29.5 (so we will use shape = 1.3 and scale = 30)

# DAN CODE REVIEW DONE UP TO HERE
















#https://discourse.mc-stan.org/t/fitting-time-to-event-data-with-weibull-hazard-using-brm-function/4638/8
#https://rileyking.netlify.app/post/bayesian-modeling-of-censored-and-uncensored-fatigue-data-in-r/#use-brm()-to-generate-a-posterior-distribution-for-shape-and-scale


### test on same population ####
# read in male elephant data
males <- readxl::read_excel('data_raw/Raw_ATE_Males_Lee220121.xlsx', sheet = 1) %>% janitor::clean_names()
males$age_2020 <- ifelse(males$dyr < 20, 2020 - males$byr, NA)
males$age_2020_cat <- ifelse(males$age_2020 < 15, 3,
                             ifelse(males$age_2020 > 14 & males$age_2020 < 20, 4,
                                    ifelse(males$age_2020 > 19 & males$age_2020 < 25, 5,
                                           ifelse(males$age_2020 > 24 & males$age_2020 < 40, 6,7))))
model_data <- males[!is.na(males$age_2020_cat),]

# create data list
N <- length(unique(model_data$casename))   # Number of individuals
K <- 7                                     # Number of age bin categories
elephants_ls <- list(
  N = N,
  K = K,
  age_category_index = model_data$age_2020_cat)           # Age category if observed in 2020

# fit Stan model to estimate ages from the categories based on a Weibull distribution (for real data these will be centred on parameters from fitting to age estimates based on another elephant population)
anp_age_ord_mod <- cmdstan_model("models/elephant_latent_age_ordinal_regression_dwf_22.03.29.stan")

age_estimation_fit <- anp_age_ord_mod$sample(
  data = elephants_ls, 
  chains = 4, 
  parallel_chains = 1,
  iter_sampling = 2000
)

# Examine the estimates.
age_est_mat <- age_estimation_fit$summary()[458:913,]
plot_data <- data.frame(age_cat = elephants_ls$age,       # Biologists original age est
                        age_true = model_data$age_2020,   # real age
                        model_age = age_est_mat$mean)     # Mean modelled age
plot_data %>%
  ggplot(aes(x = factor(age_cat), y = model_age)) +
  geom_point(size = 4, alpha = 0.6) +
  geom_hline(yintercept = c(5,10,15,20,25,40), linetype = "dashed", alpha = 0.6) +
  theme_minimal() + 
  xlab("Assigned age") + ylab("Modelled age") # plot the estimated ages against the biologist assigned ages

plot_data %>%
  ggplot(aes(x = age_true, y = model_age)) +
  geom_point(size = 4, alpha = 0.2) +
  geom_vline(xintercept = c(5,10,15,20,25,40), linetype = "dashed", alpha = 0.6) +
  theme_minimal() + 
  xlab("Actual age") + ylab("Modelled age") # plot the estimated ages against the biologist assigned ages

# Now an equivalent posterior predictive plot from draws from the distribution to show uncertainty around the mean age (then need to do it with full uncertainty).
true_ages <- age_estimation_fit$draws("true_age", format="df")
true_ages <- true_ages[1:100,1:N]
means <- apply(true_ages, 2, mean)
stdevs <- apply(true_ages, 2, sd)

summary(means)
summary(stdevs)

# I can't get th

df <- as.data.frame(do.call(rbind, true_ages)) %>%
  mutate(true_age = elephants_ls$age) %>%
  relocate(true_age) %>%
  mutate(id = 1:ncol(true_ages)) %>% # changed nrow to ncol, else it had to be the same number of draws as elephants
  relocate(id)
df <- df %>% pivot_longer(cols = 3:102) %>% dplyr::select(-name)
df %>% ggplot(aes(x=jitter(true_age), y=value, group=factor(id))) +
  geom_point(size=2,alpha=0.1) +
  #geom_vline(xintercept=c(5,10,15,20,25,40), linetype="dashed", alpha=0.6) +
  geom_hline(yintercept=c(5,10,15,20,25,40), linetype="dashed", alpha=0.6) +
  theme_bw() + 
  xlab("Assigned age") + ylab("Modelled age")
summary(df$value[which(df$true_age == 7)])
quantile(df$value[which(df$true_age == 7)], seq(0,1,0.01))


