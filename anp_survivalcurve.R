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
library(dplyr)
### Import nodes data ####
#setwd("~/Downloads/ElephantAges")
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

### estimate Weibull distribution ####
colnames(age_cens)[1] <- 'age'
age_cens$age_non0 <- age_cens$age+0.01 # wouldn't run when I allowed age = 0, so increased all values by 0.01 to make it run

age_cens$censor <- as.integer(!as.logical(age_cens$censor))

# Examine default brms priors (male elephants)
bfit_m_default <- age_cens %>% filter(sex == 'Male') %>%
  brm(age_non0 | cens(censor) ~ 1,
      data = ., family = "weibull", 
      chains = 4, cores = 4)
prior_summary(bfit_m_default)
plot(bfit_m_default)
# Defaults
#intercept - student_t(3, 2.3, 2.5)
#shape - gamma(0.01, 0.01)

# Prior predictive check
num_samples <- 500

# Default priors
prior_intercept <- rstudent_t(num_samples, 3, 2.3, 2.5)
prior_shape <- rgamma(num_samples, 0.01, 0.01)

# Priors used
prior_intercept <- rstudent_t(num_samples, 3, 5, 5)
prior_shape <- rexp(num_samples, 0.5)
prior_scale <- exp(prior_intercept) / (gamma(1 + 1 / prior_shape))

hist(prior_shape)
hist(prior_scale[prior_scale<100])

probs <- dweibull(1:100, shape = prior_shape[1], scale = prior_scale[1])
probs <- probs / sum(probs)
plot(probs, type="l", lwd=0.5, ylim = c(0,1))
for( i in 2:num_samples) {
  probs <- dweibull(1:100, shape = prior_shape[i], scale = prior_scale[i])
  probs <- probs / sum(probs)
  lines(probs, lwd=0.5)
}

# Run model for male elephants with slightly better priors
bfit_m <- age_cens %>% filter(sex == 'Male') %>%
  brm(age_non0 | cens(censor) ~ 1,
      prior = c(
        prior(student_t(3,5,5), class = Intercept),
        prior(exponential(0.5), class = shape)
      ),
      data = ., family = "weibull", 
      chains = 4, cores = 4)
plot(bfit_m)  # actually looks very similar to the default one

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
anp_age_ord_mod <- cmdstan_model("models/age_estimation/elephant_latent_age_ordinal_regression_dwf_22.03.30.stan")

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

### Helen -- repeat for females ####
colnames(age_cens)[1] <- 'age'
age_cens$age_non0 <- age_cens$age+0.01 # wouldn't run when I allowed age = 0, so increased all values by 0.01 to make it run

age_cens$censor <- as.integer(!as.logical(age_cens$censor))
hist(age_cens$age_non0[age_cens$sex == 'Female'])
hist(age_cens$age_non0[age_cens$sex == 'Female' & age_cens$censor == 0])
hist(rweibull(3000,est_shape,est_scale), breaks = 40)

# Examine default brms priors
bfit_f_default <- age_cens %>% filter(sex == 'Female') %>%
  brm(age_non0 | cens(censor) ~ 1,
      data = ., family = "weibull", 
      chains = 4, cores = 4)
prior_summary(bfit_f_default)
plot(bfit_f_default)
# Defaults
#intercept - student_t(3, 2.6, 2.5)
#shape - gamma(0.01, 0.01)

# Prior predictive check
num_samples <- 500

# Default priors
prior_intercept <- rstudent_t(num_samples, 3, 2.6, 2.5)
prior_shape <- rgamma(num_samples, 0.01, 0.01)
prior_scale <- exp(prior_intercept) / (gamma(1 + 1 / prior_shape))

hist(prior_shape)
hist(prior_scale[prior_scale<100])

probs <- dweibull(1:100, shape = prior_shape[1], scale = prior_scale[1])
probs <- probs / sum(probs)
plot(probs, type="l", lwd=0.5, ylim = c(0,1))
for( i in 2:num_samples) {
  probs <- dweibull(1:100, shape = prior_shape[i], scale = prior_scale[i])
  probs <- probs / sum(probs)
  lines(probs, lwd=0.5)
}

# Run model for male elephants with slightly better priors
bfit_f <- age_cens %>% filter(sex == 'Female') %>%
  brm(age_non0 | cens(censor) ~ 1,
      prior = c(
        prior(student_t(3,5,5), class = Intercept),
        prior(exponential(0.5), class = shape)
      ),
      data = ., family = "weibull", 
      chains = 4, cores = 4)
plot(bfit_f)  # again very similar to the default one

# extract scale and shape from posterior draws
# to convert to scale we need to both undo the link function by taking the exponent
# and then refer to the brms documentation to understand how the mean relates to the scale
bfit_f_draws <- as_draws_df(bfit_f) %>%
  mutate(scale = exp(b_Intercept) / (gamma(1 + 1 / shape)))

est_shape <- mean(bfit_f_draws$shape)
est_scale <- mean(bfit_f_draws$scale)

# Let's do a posterior predictive plot based on the mean
probs <- dweibull(1:100, shape = est_shape, scale = est_scale)
probs <- probs / sum(probs)
plot(probs) 

# Lets simulate ages for 3k elephants
hist(rweibull(3000,est_shape,est_scale), breaks = 15)

# Shape = 1.25, scale = 32.9 (so we will use shape = 1.3 and scale = 33)

### test on same population ####
# read in male elephant data
females <- readxl::read_excel('data_raw/Raw_ATE_AllElephants_Lee220118.xlsx', sheet = 1) %>% janitor::clean_names()
females <- females[females$sex == 'Female',]
females$age_2020 <- ifelse(females$age_death == 'living', 2020 - females$birth_year, NA)
females$age_2020_cat <- ifelse(females$age_2020 < 5,1,
                             ifelse(females$age_2020 > 4 & females$age_2020 < 10, 2,
                                    ifelse(females$age_2020 > 9 & females$age_2020 < 15, 3,
                                           ifelse(females$age_2020 > 14 & females$age_2020 < 20, 4,
                                                  ifelse(females$age_2020 > 19 & females$age_2020 < 35, 5, 
                                                         ifelse(females$age_2020 > 34 & females$age_2020 < 50, 6,7))))))
model_data <- females[!is.na(females$age_2020_cat),]

# create data list
N <- length(unique(model_data$id))         # Number of individuals
K <- 7                                     # Number of age bin categories
elephants_ls <- list(
  N = N,
  K = K,
  age_category_index = model_data$age_2020_cat)           # Age category if observed in 2020

# fit Stan model to estimate ages from the categories based on a Weibull distribution (for real data these will be centred on parameters from fitting to age estimates based on another elephant population)
anp_age_ord_mod <- cmdstan_model("models/age_estimation/elephant_latent_age_ordinal_regression_dwf_22.03.30.stan")

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


### compare categorical distributions amongst populations ####
anp_nodes <- janitor::clean_names(all_nodes)
anp_nodes$age_2020 <- round(as.numeric(ifelse(anp_nodes$age_death == 'living', 2020 - anp_nodes$birth_year,
                                              anp_nodes$age_death)),0)
anp_nodes$age_2020_cat <- ifelse(anp_nodes$age_2020 < 5, 1,
                               ifelse(anp_nodes$age_2020 < 10, 2,
                                      ifelse(anp_nodes$age_2020 < 15, 3,
                                             ifelse(anp_nodes$age_2020 < 20, 4,
                                                    ifelse(anp_nodes$age_2020 < 35, 5,
                                                           ifelse(anp_nodes$age_2020 < 50, 6, 7))))))
females_anp <- anp_nodes[anp_nodes$sex == 'Female',]
males_anp <- anp_nodes[anp_nodes$sex == 'Male',]

motnp_nodes <- read_csv('data_processed/motnp_elenodes_22.01.13.csv')
motnp_nodes$k <- ifelse(motnp_nodes$age_category == '0-3' | motnp_nodes$age_category == '1-2' | 
                          motnp_nodes$age_category == '3-4' |motnp_nodes$age_category == '4-5' | 
                          motnp_nodes$age_category == '0-5', 1,
                        ifelse(motnp_nodes$age_category == '5-6' | motnp_nodes$age_category == '6-7' | 
                                 motnp_nodes$age_category == '7-8' |motnp_nodes$age_category == '8-9' | 
                                 motnp_nodes$age_category == '9-10', 2,
                               ifelse(motnp_nodes$age_category == '10-15', 3,
                                      ifelse(motnp_nodes$age_category == '15-19', 4,
                                             ifelse(motnp_nodes$age_category == '20-25' | 
                                                      motnp_nodes$age_category == '20-35', 5,
                                                    ifelse(motnp_nodes$age_category == '25-40' | 
                                                             motnp_nodes$age_category == '35-50', 6, 7))))))
females_motnp <- motnp_nodes[motnp_nodes$sex == 'F',]
males_motnp <- motnp_nodes[motnp_nodes$sex == 'M',]

par(mfrow = c(3,3))
hist(anp_nodes$age_2020_cat, main = 'ANP, all elephants',
     col = 'grey', xlab = 'age category')
hist(females_anp$age_2020_cat, main = 'ANP, all females',
     col = 'magenta', xlab = 'age category')
hist(males_anp$age_2020_cat, main = 'ANP, all males',
     col = 'blue', xlab = 'age category')
hist(anp_nodes$age_2020_cat[anp_nodes$age_death == 'living'], main = 'ANP, living elephants',
     col = 'grey', xlab = 'age category')
hist(females_anp$age_2020_cat[females_anp$age_death == 'living'], main = 'ANP, living females',
     col = 'magenta', xlab = 'age category')
hist(males_anp$age_2020_cat[males_anp$age_death == 'living'], main = 'ANP, living males',
     col = 'blue', xlab = 'age category')
hist(motnp_nodes$k, main = 'MOTNP, living elephants',
     col = 'grey', xlab = 'age category')
hist(females_motnp$k, main = 'MOTNP, living females',
     col = 'magenta', xlab = 'age category')
hist(males_motnp$k, main = 'MOTNP, living males',
     col = 'blue', xlab = 'age category')





