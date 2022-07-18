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
library(rstan)
library(rethinking)
library(cmdstanr)
library(lubridate)
library(posterior)
library(bayesplot)
library(survival)
library(brms)
library(tidybayes)
library(ggthemes)

### Import sightings data ####
info <- readxl::read_excel(path = 'data_raw/Raw_ATE_Sightings_Fishlock220224.xlsx', sheet = 2)
vars <- info[1:22,1:2]
activities <- info[1:10, 11:12]; colnames(activities) <- c('code','definition')
habitats <- info[14:30, 11:12] ; colnames(habitats) <- c('code','definition')
codes <- data.frame(type = c(rep('activity', nrow(activities)), rep('habitat',nrow(habitats))),
                    code = c(activities$code, habitats$code),
                    definition = c(activities$definition, habitats$definition))
rm(info, activities, habitats)

ate <- readxl::read_excel(path = 'data_raw/Raw_ATE_Sightings_Fishlock220224.xlsx', sheet = 1)
ate <- janitor::clean_names(ate)
ate$id <- paste('M',ate$casename, sep = '')
ate$node_id <- as.integer(as.factor(ate$casename))
ate$obs_date <- lubridate::as_date(ate$obs_date)
ate <- separate(ate, obs_time, into = c('wrong_date','correct_time'), remove = F, sep = ' ')
ate$correct_time_hms <- hms::as_hms(ate$correct_time)
ate$corrected_time <- lubridate::hour(ate$correct_time_hms)*60*60 + lubridate::minute(ate$correct_time_hms) + lubridate::second(ate$correct_time_hms)

lu <- function(x) { length(unique(x)) }
ate_nums <- tapply(X = ate$obs_num, INDEX = ate$obs_date, FUN = lu )

ate$obs_num <- ifelse(ate$obs_num == '0','00',
                      ifelse(ate$obs_num == '0a','0A',
                             ifelse(ate$obs_num == '0b','0B',
                                    ifelse(ate$obs_num == '1','01', ate$obs_num))))
ate$obs_num_std <- NA
for(i in 1:length(ate)){
  date_row <- ate[ate$obs_date == ate$obs_date[i],]
  date_row$obs_num_std <- as.integer(as.factor(sort(date_row$obs_num)))
  ate$obs_num_std[i] <- date_row$obs_num_std[which(date_row$obs_id == ate$obs_id[i])[1]]
}
## rearrange
ate <- ate[,c(1:3,25,26,4,7,27:28,8,29,9,12:24)]
colnames(ate)

## clean environment
rm(ate_nums, date_row, i)

### Import nodes data ####
males <- readxl::read_excel('data_raw/Raw_ATE_Males_Lee220121.xlsx')
males <- janitor::clean_names(males)
males$id <- paste('M', males$casename, sep = '')

all_nodes <- readxl::read_excel('data_raw/Raw_ATE_Males_Lee220118.xlsx')
boxplot(as.numeric(all_nodes$Age_death) ~ all_nodes$Death_accuracy) # adds NA values for anything classed as "living"
str(all_nodes)

### calculate oldest age each observed ####
colnames(all_nodes)
ages <- all_nodes[,c(1:3,15,16,18:20)]
ages <- janitor::clean_names(ages)
unique(ages$death_accuracy)
ages$dacc <- ifelse(ages$death_accuracy == 'within 1mo', (1/12)*3,
                    ifelse(ages$death_accuracy == 'within 3mo', (3/12)*2,
                           ifelse(ages$death_accuracy == 'within 6mo', (6/12)*1.5,
                                  ifelse(ages$death_accuracy == 'within 1yr', 1.5, NA))))
ages$bacc <- ifelse(ages$birth_accuracy == '+/- 1mo', (1/12)*3,
                    ifelse(ages$birth_accuracy == '+/- 6mo', (6/12)*1.5,
                           ifelse(ages$birth_accuracy == '+/- 1yr', 1.5,
                                  ifelse(ages$birth_accuracy == '+/- 2.5yrs', 2.5, 5))))
sort(unique(ages$age_death))
sort(unique(ages$dacc))
sort(unique(ages$bacc))
#length(which(is.na(ages$age_death) == TRUE))
#ages[is.na(ages$age_death),]
ages$age_death_num <- as.numeric(ifelse(ages$age_death == 'living', 2021 - ages$birth_year, ages$age_death))

ages$youngest_death <- ifelse(ages$death_accuracy == 'living',
                              ages$age_death_num + ages$bacc,
                              ages$age_death_num - ages$dacc + ages$bacc)
ages$oldest_death <- ifelse(ages$death_accuracy == 'living',
                            ages$age_death_num - ages$bacc,
                            ages$age_death_num + ages$dacc - ages$bacc)
str(ages)

ages$id_no <- ifelse(ages$sex == 'unknown', paste('U', ages$id, sep = ''),
                     ifelse(ages$sex == 'Female', paste('F', ages$id, sep = ''),
                            ifelse(ages$sex == 'Male', paste('M', ages$id, sep = ''),
                                   'check')))
ages$age_death_round <- round(ages$age_death_num, 0)
table(ages$age_death_round, ages$sex)
hist(ages$age_death_round, breaks = 75)
hist(ages$age_death_round[which(ages$sex == 'Male')], breaks = 70)
hist(ages$age_death_round[which(ages$sex == 'Female')], breaks = 75)
ages <- ages[,c(14,1,3,8,11,15,12:13)]
unique(ages2$age_death)

ages$youngest_death

age_df <- data.frame(id = rep(ages$id, each = 75),
                     age = rep(1:75, length(unique(ages$id))),
                     died = NA,
                     censor = NA)
for(i in 1:nrow(age_df)){
  age_row <- ages[ages$id == age_df$id[i],]
  age_df$died[i] <- ifelse(age_row$age_death_round[1] > age_df$age[i], 0, 1)
  age_df$censor[i] <- ifelse(age_row$age_death == 'living', 'TRUE','FALSE')
}

ages$censor <- ifelse(ages$age_death == 'living', 'TRUE', 'FALSE')
table(ages$censor)

age_cens <- ages[!is.na(ages$age_death_num),c('age_death_num','censor')]

rm(age_row, codes, vars)
rm(age_df, ages, all_nodes, ate, males)

### Weibull example - DWF, 18th March 2022 ####
# simulate an age distribution based on a Weibull distribution where we set the shape and scale parameters.
probs <- dweibull(1:100, shape = 1.2, scale = 30)
probs <- (probs/sum(probs))
plot(probs) 

# create a fictional population with ages selected from that distribution
N <- 100 # Number of individuals
K <- 6   # Number of age bin categories
elephants_ls <- list(
  N = N,
  K = K,
  age = sample(1:100, N, replace=TRUE, prob=probs)) # Simulated ages drawn from the Weibull distribution plus an error in age estimation

# simulate observing ages with error and then binning them into age groups
E <- 3 # Error (in either direction) of age estimates
elephants_ls$age_guess <- elephants_ls$age + sample(-E:E, N, replace=TRUE)
elephants_ls$age_category_index <- sapply(elephants_ls$age_guess, function(x) which.max(x < c(15, 30, 45, 60, Inf)))
hist(elephants_ls$age_category)

# look at the actual age verses the biologist assigned age and look at the thresholds.
data.frame(elephants_ls) %>%
  ggplot(aes(x=age, y=age_guess, col=factor(age_category_index))) +
  geom_point(size=4,alpha=0.6) +
  geom_vline(xintercept=c(15,30,45,60,75), col=factor(2:6), linetype="dashed", alpha=0.6) +
  theme_minimal() + 
  xlab("Assigned age") + ylab("Age") # chance of being mis-classified is lower if an actual age is in the middle of the age category. Also, there is mostly a bias towards ages within the age class actually being towards the lower end of the class.
# fit Stan model to estimate ages from the categories based on a Weibull distribution (for real data these will be centred on parameters from fitting to age estimates based on another elephant population)
latent_age_ordinal_model <- cmdstan_model("models/elephant_latent_age_ordinal_regression_dwf_22.03.18.stan")
age_estimation_fit <- latent_age_ordinal_model$sample(
  data = elephants_ls, 
  chains = 4, 
  parallel_chains = 1,
  iter_sampling = 2000
)

# Examine the estimates.
age_est_mat <- age_estimation_fit$summary()[102:201, ]
plot_data <- data.frame(age = elephants_ls$age,       # Biologists original age est
                        model_age = age_est_mat$mean) # Mean modelled age
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

### Weibull with real values ####
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

# play with shapes of Weibull just to see how they look
play <- data.frame(draws = rweibull(n = 1000, shape = 1, scale = 10),
                   position = NA)
play$position <- as.integer(as.factor(play$draws))
plot(play$draws ~ play$position)

play <- data.frame(draws = rweibull(n = 1000, shape = 3, scale = 10),
                   position = NA)
play$position <- as.integer(as.factor(play$draws))
plot(play$draws ~ play$position)

play <- data.frame(draws = rweibull(n = 1000, shape = 5, scale = 10),
                   position = NA)
play$position <- as.integer(as.factor(play$draws))
plot(play$draws ~ play$position)

play <- data.frame(draws = rweibull(n = 1000, shape = 10, scale = 10),
                   position = NA)
play$position <- as.integer(as.factor(play$draws))
plot(play$draws ~ play$position)

play <- data.frame(draws = rweibull(n = 1000, shape = 1, scale = 1),
                   position = NA)
play$position <- as.integer(as.factor(play$draws))
plot(play$draws ~ play$position)

play <- data.frame(draws = rweibull(n = 1000, shape = 3, scale = 1),
                   position = NA)
play$position <- as.integer(as.factor(play$draws))
plot(play$draws ~ play$position)

play <- data.frame(draws = rweibull(n = 1000, shape = 5, scale = 1),
                   position = NA)
play$position <- as.integer(as.factor(play$draws))
plot(play$draws ~ play$position)

play <- data.frame(draws = rweibull(n = 1000, shape = 10, scale = 1),
                   position = NA)
play$position <- as.integer(as.factor(play$draws))
plot(play$draws ~ play$position) # shape = rate of initial incline, scale = only seems to affect y-axis

# adapt to ATE data
head(ages)
dens(ages$age_death_round[which(ages$sex == 'Male')], col = 'red')
dens(ages$age_death_round[which(ages$sex == 'Female')], add = T, col = 'blue')
dens(ages$age_death_round[which(ages$sex == 'unknown')], add = T, col = 'purple', lty = 3)
abline(v = 0)

play <- data.frame(draws = rweibull(n = 1000, shape = 10, scale = 1),
                   draws_rev = NA,
                   position = NA)
play$draws_rev <- 1 - play$draws
play$position <- as.integer(as.factor(play$draws))
plot(play$draws_rev ~ play$position) # shape = rate of initial incline, scale = only seems to affect y-axis
abline(h = 0)

play <- data.frame(draws = rweibull(n = 1000, shape = 4, scale = 0.8),
                   draws_rev = NA,
                   position = NA)
play$draws_rev <- 1 - play$draws
play$position <- as.integer(as.factor(play$draws))
plot(play$draws_rev ~ play$position) # shape = rate of initial incline, scale = only seems to affect y-axis
abline(h = 0)

# Need to fit a Weibull distribution to age data, and that will then go into upcoming model
#m <- ulam(alist(
#  age_death_round ~ dweibull(shape, scale),
#  shape ~ dexp(1),
#  scale ~ dexp(1)
#), data = ages)

#m <- cmdstan_model("models/weibull.stan")
#ls <- list(
#  N = nrow(ages),
#  age_death_round = ifelse(ages$age_death_round == 0, 0.01, ages$age_death_round)
#)
#summary(ls$age_death_round)
#weibull <- m$sample(
#  data = ls,
#  chains = 1,
#  iter_sampling = 2000
#)

colnames(age_cens)[1] <- 'age'
age_cens$age_non0 <- age_cens$age+0.01
age_cens %>%
  mutate(censor = if_else(censor == 0, 1, 0)) %>%
  brm(age_non0 | cens(censor) ~ 1, data = ., family = "weibull", 
      refresh = 2e3, cores = 4) -> bfit

print(bfit)

as_sample_tibble(bfit) %>%
  mutate_at(vars(b_Intercept), exp) %>%
  ggplot(aes(x = b_Intercept, y = shape)) +
  geom_density_2d(colour = "red", size = 1) +
  geom_point(alpha = 0.1) +
  #geom_hline(yintercept = 1, colour = "blue") +
  #geom_vline(xintercept = 10, colour = "blue") +
  xlab("Scale") +
  ylab("Shape") +
  ggtitle("Right-censored weibull(shape, scale)") +
  theme_fivethirtyeight()












# test on same population
males$age_2020 <- ifelse(males$dyr < 20, 2020 - males$byr, NA)
males$age_2020_cat <- ifelse(males$age_2020 < 15, 3,
                             ifelse(males$age_2020 > 14 & males$age_2020 < 20, 4,
                                    ifelse(males$age_2020 > 19 & males$age_2020 < 25, 5,
                                           ifelse(males$age_2020 > 24 & males$age_2020 < 40, 6,7))))
model_data <- males[!is.na(males$age_2020_cat),]
N <- length(unique(model_data$casename))   # Number of individuals
K <- 7                                     # Number of age bin categories
elephants_ls <- list(
  N = N,
  K = K,
  age = model_data$age_2020_cat)           # Age category if observed in 2020

# fit Stan model to estimate ages from the categories based on a Weibull distribution (for real data these will be centred on parameters from fitting to age estimates based on another elephant population)
latent_age_ordinal_model <- cmdstan_model("models/elephant_latent_age_ordinal_regression_dwf_22.03.18.stan")
age_estimation_fit <- latent_age_ordinal_model$sample(
  data = elephants_ls, 
  chains = 4, 
  parallel_chains = 1,
  iter_sampling = 2000
)

# Examine the estimates.
age_est_mat <- age_estimation_fit$summary()[102:201, ]
plot_data <- data.frame(age = elephants_ls$age,       # Biologists original age est
                        model_age = age_est_mat$mean) # Mean modelled age
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





