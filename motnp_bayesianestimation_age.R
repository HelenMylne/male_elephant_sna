### Information ####
# Inputs = draws from MOTNP edge weight estimation model, node characteristics from MOTNP database, Weibull distribution parameters estimated from ANP ages of death.
# Fit the ANP Weibull distribution to MOTNP individuals, then those age distributions can be used in regression models

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
motnp <- read_csv('data_processed/motnp_elenodes_22.01.13.csv')
motnp$id_no[is.na(motnp$age_category)]
motnp$age_category <- ifelse(motnp$id == 'U8', '0-5', motnp$age_category)
motnp$k <- ifelse(motnp$age_category == '0-3' | motnp$age_category == '1-2' | 
                    motnp$age_category == '3-4' |motnp$age_category == '4-5' | 
                    motnp$age_category == '0-5', 1,
                  ifelse(motnp$age_category == '5-6' | motnp$age_category == '6-7' | 
                           motnp$age_category == '7-8' |motnp$age_category == '8-9' | 
                           motnp$age_category == '9-10', 2,
                         ifelse(motnp$age_category == '10-15', 3,
                                ifelse(motnp$age_category == '15-19', 4,
                                       ifelse(motnp$age_category == '20-25' | motnp$age_category == '20-35', 5,
                                              ifelse(motnp$age_category == '25-40' | 
                                                       motnp$age_category == '35-50', 6, 7))))))
unique(motnp$k)

### create data lists ####
males <- motnp[motnp$sex == 'M',]
females <- motnp[motnp$sex == 'F',]

# Males
N_m <- length(unique(males$id_no))     # Number of individuals
K <- 7                                 # Number of age bin categories
males_ls <- list(
  N = N_m,
  K = K,
  age_category_index = males$k)        # age category

# Females
N_f <- length(unique(females$id_no))   # Number of individuals
females_ls <- list(
  N = N_f,
  K = K,
  age_category_index = females$k)      # age category

### 1) Fit male model of ages ####
# Shape = 1.29, scale = 29.5 (so we will use shape = 1.3 and scale = 30)
motnp_male_age_ord_mod <- cmdstan_model("models/age_estimation/male_age_ordreg_22.03.31.stan")
male_age_fit <- motnp_male_age_ord_mod$sample(
  data = males_ls, 
  chains = 4, 
  parallel_chains = 1,
  iter_sampling = 2000
)

# Examine the estimates.
age_est_m <- male_age_fit$summary()[247:491,]
plot_male <- data.frame(age_cat = males_ls$age,         # Biologists original age est
                        model_age = age_est_m$mean)     # Mean modelled age
plot_male %>%
  ggplot(aes(x = factor(age_cat), y = model_age)) +
  geom_point(size = 4, alpha = 0.6) +
  geom_hline(yintercept = c(5,10,15,20,25,40), linetype = "dashed", alpha = 0.6) +
  theme_minimal() + 
  xlab("Age category") + ylab("Modelled age") # plot the estimated ages against the biologist assigned ages

# Now an equivalent posterior predictive plot from draws from the distribution to show uncertainty around the mean age (then need to do it with full uncertainty).
true_ages_m <- male_age_fit$draws("true_age", format="df")
true_ages_m <- true_ages_m[1:100,1:N_m]
means_m <- apply(true_ages_m, 2, mean)
stdevs_m <- apply(true_ages_m, 2, sd)

summary(means_m)
summary(stdevs_m)

df_m <- as.data.frame(do.call(rbind, true_ages_m)) %>%
  mutate(true_age = males_ls$age) %>%
  relocate(true_age) %>%
  mutate(id = 1:ncol(true_ages_m)) %>% # changed nrow to ncol, else it had to be the same number of draws as elephants
  relocate(id)
df_m <- df_m %>% pivot_longer(cols = 3:102) %>% dplyr::select(-name)
df_m %>% ggplot(aes(x=jitter(true_age), y=value, group=factor(id))) +
  geom_point(size=2,alpha=0.1) +
  geom_boxplot(aes(group = factor(true_age)), col = 'skyblue', fill = 'transparent') +
  geom_hline(yintercept=c(5,10,15,20,25,40), linetype="dashed", alpha=0.6) +
  theme_bw() + 
  xlab("Assigned age") + ylab("Modelled age")
summary(df_m$value[which(df_m$true_age == 7)])
quantile(df_m$value[which(df_m$true_age == 7)], seq(0,1,0.01)) # top 3% over 70 years

### 2) Fit female model of ages ####
# Shape = 1.25, scale = 32.9 (so we will use shape = 1.3 and scale = 33)
motnp_female_age_ord_mod <- cmdstan_model("models/age_estimation/female_age_ordreg_22.03.31.stan")
female_age_fit <- motnp_female_age_ord_mod$sample(
  data = females_ls, 
  chains = 4, 
  parallel_chains = 1,
  iter_sampling = 2000
)

# Examine the estimates.
age_est_f <- female_age_fit$summary()[153:303,]
plot_female <- data.frame(age_cat = females_ls$age,         # Biologists original age est
                        model_age = age_est_f$mean)     # Mean modelled age
plot_female %>%
  ggplot(aes(x = factor(age_cat), y = model_age)) +
  geom_point(size = 4, alpha = 0.6) +
  geom_hline(yintercept = c(5,10,15,20,35,50), linetype = "dashed", alpha = 0.6) +
  theme_minimal() + 
  xlab("Age category") + ylab("Modelled age") # plot the estimated ages against the biologist assigned ages

# Now an equivalent posterior predictive plot from draws from the distribution to show uncertainty around the mean age (then need to do it with full uncertainty).
true_ages_f <- female_age_fit$draws("true_age", format="df")
true_ages_f <- true_ages_f[1:100,1:N_f]
means_f <- apply(true_ages_f, 2, mean)
stdevs_f <- apply(true_ages_f, 2, sd)

summary(means_f)
summary(stdevs_f)

df_f <- as.data.frame(do.call(rbind, true_ages_f)) %>%
  mutate(true_age = females_ls$age) %>%
  relocate(true_age) %>%
  mutate(id = 1:ncol(true_ages_f)) %>% # changed nrow to ncol, else it had to be the same number of draws as elephants
  relocate(id)
df_f <- df_f %>% pivot_longer(cols = 3:102) %>% dplyr::select(-name)
df_f %>% ggplot(aes(x=jitter(true_age), y=value, group=factor(id))) +
  geom_point(size=2,alpha=0.1) +
  geom_boxplot(aes(group = factor(true_age)), col = 'magenta', fill = 'transparent') +
  geom_hline(yintercept=c(5,10,15,20,35,50), linetype="dashed", alpha=0.6) +
  theme_bw() + 
  xlab("Assigned age") + ylab("Modelled age")
summary(df_f$value[which(df_f$true_age == 7)])
quantile(df_f$value[which(df_f$true_age == 7)], seq(0,1,0.01)) # top 6% over 75 years


