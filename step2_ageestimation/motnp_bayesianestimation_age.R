#### information #####
# To estimate the ages of individuals in other populations, first come up with a survival curve for the Amboseli elephants which have good birth and death data.
# survival curve options = exponential/logistic/Gompertz/Weibull, where shape = simple/bathtub/Makeham
# follow through Mia's code to determine which curve will best fit Amboseli data. she ran it on a high performance cluster, so it’s set up to run in parallel, but the basta model can also be run in parallel on its own (just like Stan). we will need to implement it as a custom function in Stan for the final analysis, given that we need to define it in the prior anyway. 

#### load packages ####
#library(tidyverse) ; library(cmdstanr) ; library(ggdist) ; library(posterior) ; library(bayesplot) ; library(rstan) ; library(igraph) ; library(LaplacesDemon)
library(tidyverse, lib.loc = '../packages')
library(cmdstanr, lib.loc = '../packages')
library(ggdist, lib.loc = '../packages')
library(posterior, lib.loc = '../packages')
library(bayesplot, lib.loc = '../packages')
library(rstan, lib.loc = '../packages')
library(igraph, lib.loc = '../packages')
library(LaplacesDemon, lib.loc = '../packages')

set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0/')

#### load model ####
# read in Stan model to estimate ages based on Gompertz bathtub distribution from ANP
latent_age_ordinal_model <- cmdstan_model("models/motnp_age_ordinal_regression.stan")

#### simulate process of assigning age categories to elephants ####
gompertz_bt <- function(a0, a1, c, b0, b1, age){  # create custom function for mortality distribution
  gompertz <- exp(b0 + b1*age)
  bathtub <- exp(a0 - a1*age) + c + gompertz
  return(bathtub)
}

# draw new output curves
max_age <- 60                                     # maximum realistic age of males in the population
mortality <- gompertz_bt(a0 = -5.13,              # set priors
                         a1 = 3.00,
                         c  = 0.026,
                         b0 = -5.08,
                         b1 = 0.09, 1:max_age)
plot(mortality)                                   # plot mortality curve
probs <- 1 - (mortality/max(mortality))           # survival = 1-mortality
plot(probs)                                       # plot probability of survival

# create a fictional population with ages selected from that distribution
N <- 200                # number of eles in fictional population
K <- 8                  # number of thresholds for categories
elephants_ls <- list(   # make list of elephants with death ages sampled from survival curve
  N = N,
  K = K,
  age = sample(1:max_age, N, replace = TRUE, prob = probs)
)
hist(elephants_ls$age)  # histogram

# simulate observing ages with error and then binning them into age groups
error <- 3 # Error (in either direction) of age estimates
elephants_ls$age_guess <- elephants_ls$age + sample(-error:error, N, replace = TRUE)  # simulate observations of true age
elephants_ls$age_category_index <- sapply(elephants_ls$age_guess,                     # simulate categories from observations
                                          function(x) which.max(x < c(5, 10, 15, 20, 25, 40, 60, Inf)))
hist(elephants_ls$age_category_index)                                                 # bar chart age categories

# look at actual age vs biologist assigned age: chance of being mis-classified is lower if actual age is in middle of age category, and biased towards lower end of class.
data.frame(elephants_ls) %>%       
  ggplot(aes(x = age, y = age_guess, col = factor(age_category_index))) +             # plot assigned vs true age
  geom_point(size = 4, alpha = 0.6) +                                                 # plot points
  geom_vline(xintercept = c(5,10,15,20,25,40,60), col = factor(1:7), linetype = "dashed", alpha = 0.6) +        # add vertical lines at category boundaries
  theme_minimal() + 
  xlab("Assigned age") + ylab("Age")

#### fit model to simulated dataset ####
# Fit model with cmdstanr
age_estimation_fit <- latent_age_ordinal_model$sample(
  data = elephants_ls, 
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000
)

# Examine the estimates
age_est_mat <- age_estimation_fit$summary()[202:401, ] # true age estimates for each elephant
summary(age_est_mat)
hist(age_est_mat$mean)                                 # plot histogram of mean age values
hist(age_est_mat$rhat, breaks = 20)                    # check rhat values

plot_data <- data.frame(age = elephants_ls$age,        # Biologists original age est
                        model_age = age_est_mat$mean)  # Mean modelled age
plot_data %>%
  ggplot(aes(x = factor(age), y = model_age)) +        # true age vs modelled age
  geom_point(size = 4,col = 'blue', alpha = 0.6) +     # add points
  geom_vline(xintercept = c(5,10,15,20,25,40,60), linetype = "dashed", alpha = 0.6) +    # add vertical lines at category boundaries
  geom_hline(yintercept = c(5,10,15,20,25,40,60), linetype = "dashed", alpha = 0.6) +    # add horizontal lines at category boundaries
  geom_abline(slope = 1, intercept = 0)+               # x=y line to show where values should fall if perfectly estimated
  scale_y_continuous(limits = c(0,60))+
  theme_minimal() + 
  xlab("Assigned age") + ylab("Modelled age")

#### load MOTNP data ####
motnp_males <- read_csv('../data_processed/step1_dataprocessing/motnp_elenodes.csv') %>% filter(sex == 'M')  # select male nodes only 
unique(motnp_males$age_category)
motnp_males$age_cat_id <- ifelse(motnp_males$age_category == "0-3", 1,     # standardise calf categories
                                 ifelse(motnp_males$age_category == "3-4", 1,
                                        ifelse(motnp_males$age_category == "1-2", 1,
                                               ifelse(motnp_males$age_category == "7-8", 1,
                                                      ifelse(motnp_males$age_category == "4-5", 1, 
                                                             ifelse(motnp_males$age_category ==  "6-7", 1,
                                                                    ifelse(motnp_males$age_category == "8-9", 1, 
                                                                           ifelse(motnp_males$age_category == "5-6", 1, NA))))))))
motnp_males$age_cat_id <- ifelse(motnp_males$age_category == '9-10', 2,    # convert to numeric categories
                                 ifelse(motnp_males$age_category == '10-15', 3,
                                        ifelse(motnp_males$age_category == '15-19', 4,
                                               ifelse(motnp_males$age_category == '20-25', 5,
                                                      ifelse(motnp_males$age_category == '25-40', 6,
                                                             ifelse(motnp_males$age_category == '40+', 7,
                                                                    motnp_males$age_cat_id))))))

#### create data list ####
N_motnp <- nrow(motnp_males)  # number of males
K <- 8                        # number of age thresholds
motnp_ls <- list(             # generate list of males to be included
  N = N_motnp,
  K = K,
  age_category_index = motnp_males$age_cat_id)
hist(motnp_ls$age_category_index)

#### fit model to MOTNP data ####
# Fit model with cmdstanr
age_motnp_fit <- latent_age_ordinal_model$sample(
  data = motnp_ls, 
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000
)

# Examine the estimates
age_est_mat <- age_motnp_fit$summary()[(N_motnp+2):(N_motnp*2+1), ] # true ages of elephants
summary(age_est_mat)
hist(age_est_mat$mean)                                              # plot histogram of mean age values
hist(age_est_mat$rhat, breaks = 20)                                 # check rhat values

plot_data <- data.frame(age = ifelse(motnp_ls$age == 1, 3,          # set dummy "true" age in centre of each age category (ONLY for plotting purposes)
                                     ifelse(motnp_ls$age == 2, 8,
                                            ifelse(motnp_ls$age == 3, 12,
                                                   ifelse(motnp_ls$age == 4, 18,
                                                          ifelse(motnp_ls$age == 5, 22, 
                                                                 ifelse(motnp_ls$age == 6, 32, 45)))))),
                        model_age = age_est_mat$mean)               # Mean modelled age

plot_data %>%
  ggplot(aes(x = factor(age), y = model_age)) +        # true age vs modelled age
  geom_point(size = 4,col = 'blue', alpha = 0.6) +     # add points
  #geom_vline(xintercept = c(5,10,15,20,25,40,60), linetype = "dashed", alpha = 0.6) +     # add vertical lines at category boundaries
  geom_hline(yintercept = c(5,10,15,20,25,40,60), linetype = "dashed", alpha = 0.6) +      # add horizontal lines at category boundaries
  #geom_abline(slope = 1, intercept = 0)+              # x=y line to show where values should fall if perfectly estimated
  scale_y_continuous(limits = c(0,60))+
  theme_minimal() + 
  xlab("Assigned age") + ylab("Modelled age")

# posterior predictive plot using draws from distribution to show uncertainty around mean age
true_ages <- age_motnp_fit$draws("true_age", format="df")    # data frame of true ages produced by model
#mcmc_dens(true_ages)
true_ages <- true_ages[,1:N_motnp]                           # select only columns relevant to individuals

df <- as.data.frame(do.call(rbind, true_ages)) %>%           # create small data frame of just age and ID
  mutate(age_cat = motnp_ls$age) %>% relocate(age_cat) %>%
  mutate(ID = motnp_males$id) %>% relocate(ID)

df <- df %>% pivot_longer(cols = 3:ncol(df)) %>% select(-name)     # convert to long format

df$age_cat_centre <- ifelse(df$age_cat == 1, (0+5)/2,              # for plot ONLY, set "age" as central value in each category
                            ifelse(df$age_cat == 2, (5+10)/2,
                                   ifelse(df$age_cat == 3, (10+15)/2,
                                          ifelse(df$age_cat == 4, (15+20)/2,
                                                 ifelse(df$age_cat == 5, (20+25)/2,
                                                        ifelse(df$age_cat == 6, (25+40)/2, 45))))))

df %>% ggplot(aes(x=age_cat_centre, y=value, group=factor(ID))) +  # plot intervals for each category against values set above
  geom_point(size=2,col = 'blue', alpha=0.1) +
  geom_vline(xintercept = c(5,10,15,20,25,40,60), linetype = "dashed", alpha = 0.6) +     # add vertical lines at category boundaries
  geom_hline(yintercept = c(5,10,15,20,25,40,60), linetype = "dashed", alpha = 0.6) +     # add horizontal lines at category boundaries
  theme_bw() + 
  xlab("Assigned age") + ylab("Modelled age")

df %>% ggplot() +                                                  # plot intervals for each category against values set above
  geom_violin(aes(x = age_cat_centre, y = value, group = factor(age_cat)), fill = rgb(0,0,1,0.8))+
  #geom_point(aes(x = true_age, y = value, group = factor(ID)), size = 2, col = 'red', alpha = 0.1) +
  geom_vline(xintercept = 0, alpha = 0.6) +                                               # add vertical lines at 0
  geom_vline(xintercept = c(5,10,15,20,25,40,60), linetype = "dashed", alpha = 0.6) +     # add vertical lines at category boundaries
  geom_hline(yintercept = 0, alpha = 0.6) +                                               # add horizontal lines at 0
  geom_hline(yintercept = c(5,10,15,20,25,40,60), linetype = "dashed", alpha = 0.6) +     # add horizontal lines at category boundaries
  theme_bw() + 
  xlab("Assigned age") + ylab("Modelled age")+
  theme(axis.text = element_text(size = 14))

df %>% ggplot() +                                                 # plot intervals for each category against values set above
  geom_violin(aes(x = age_cat_centre, y = value, group = factor(age_cat), fill = factor(age_cat, levels = c(7:3,NA,NA))))+
  geom_vline(xintercept = 0, alpha = 0.6) +                                               # add vertical lines at 0
  geom_vline(xintercept = c(5,10,15,20,25,40,60), linetype = "dashed", alpha = 0.6) +     # add vertical lines at category boundaries
  geom_hline(yintercept = 0, alpha = 0.6) +                                               # add horizontal lines at 0
  geom_hline(yintercept = c(5,10,15,20,25,40,60), linetype = "dashed", alpha = 0.6) +     # add horizontal lines at category boundaries
  theme_bw() + 
  xlab("Assigned age") + ylab("Modelled age")+
  theme(axis.text = element_text(size = 14), legend.position = 'none', axis.title = element_text(size = 18))+
  scale_fill_viridis_d()

### save output
colnames(true_ages) <- motnp_males$id
saveRDS(true_ages, file = '../data_processed/motnp_ageestimates_mcmcoutput.rds')          # save output for next steps
save.image('motnp_ageestimation.RData')

#### extract probability distributions ####
age_probs <- age_motnp_fit$draws(format = 'df')
colnames(age_probs)
age_probs <- age_probs[,c('a0','a0_std','a1','a1_std',
                          'b0','b0_std','b1','b1_std',
                          'c','c_std','sigma_age')]
age_probs

## plot parameters
par(mfrow = c(3,4))
for(i in 1:ncol(age_probs)){
  hist(as.matrix(age_probs[,i]),
       main = colnames(age_probs)[i],
       xlab = 'parameter value')
}
par(mfrow = c(1,1))

# draw new output curves
max_age <- 60                             # maximum realistic age of males in the population
gompertz_bt <- function(a0, a1, c, b0, b1, age){  # create custom function for mortality distribution
  gompertz <- exp(b0 + b1*age)
  bathtub <- exp(a0 - a1*age) + c + gompertz
  return(bathtub)
}
mortality <- gompertz_bt(a0 = mean(age_probs$a0),
                         a1 = mean(age_probs$a1),
                         c  = mean(age_probs$c),
                         b0 = mean(age_probs$b0),
                         b1 = mean(age_probs$b1),
                         1:max_age)
plot(mortality)                           # plot mortality curve
probs <- 1 - (mortality/max(mortality))   # survival = 1-mortality
plot(probs)                               # plot probability of survival

#### plot distribution per elephant ####
rm(list = ls()) ; gc() ; load('motnp_edgeweights_conditionalprior.RData') ; rm(list = ls()[ ! ls() %in% 'counts_df']) ; gc() ; load('motnp_ageestimation.RData')
eles <- unique(c(counts_df$id_1, counts_df$id_2))

df <- df %>% 
  filter(ID %in% eles)
par(mfrow = c(5,4))
for(i in 1:length(eles)){
  x <- df %>% 
    filter(ID == eles[i])
  hist(x$value)
}
par(mfrow = c(1,1))
