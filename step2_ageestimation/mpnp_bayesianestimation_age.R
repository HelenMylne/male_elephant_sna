#### information #####
# Script to estimate ages of MPNP elephants from categorical ages.
# Categories based on first time window only: take the largest dataset and then add time difference for next time windows. -- HOW DOES THIS WORK IF THERE ARE ELEPHANTS ONLY PRESENT IN LATER TIME WINDOWS?!!
# OPTION 1: RUN ALL TIME WINDOWS AS A MIXTURE MODEL WITH POOLING ACROSS TIME WINDOWS
# OPTION 2: SELECT A RANDOM POSTERIOR DISTRIBUTION FROM THE POOL OF OTHER INDIVIDUALS WITH THE SAME CATEGORICAL ESTIMATE
# OPTION 3: AMALGAMATE POSTERIOR DISTRIBUTION DRAWS OF ALL INDIVIDUALS SHARING A CATEGORY AND THEN TKAE 8000 RANDOM DRAWS FROM THAT AMALGAMATED POSTERIOR

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

#### create MPNP1 data list ####
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

#### fit model to MPNP1 data -- divergent transitions = 3% ####
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

#### save MPNP1 output ####
colnames(true_ages) <- mpnp1_males$elephant_id
saveRDS(true_ages, '../data_processed/mpnp1_ageestimates_mcmcoutput.rds')

#### produce estimate distributions for MPNP2-5 ####
time_diff <- periods[2]-periods[1]         # duration of time window in days
time_diff <- as.numeric(time_diff)/365.25  # duration of time window in years
mpnp1 <- as.data.frame(readRDS('../data_processed/mpnp1_bayesian_agedistribution.rds'))

## time window 2 ####
mpnp2_long <- mpnp_long[mpnp_long$date >= periods[2] & mpnp_long$date < periods[3],]

# identify category
mpnp2_long$age_range_id <- as.numeric(mpnp2_long$age_range_id)
mpnp2_long$age_mid <- NA
mpnp2_long$age_range <- NA
for(i in 1:nrow(mpnp2_long)){
  ele <- mpnp2_long[mpnp2_long$elephant_id == mpnp2_long$elephant_id[i],]
  ele <- ele[!is.na(ele$age_range_id),]
  ele <- ele[ele$age_range_id != 10,]
  mpnp2_long$age_mid[i] <- ifelse(nrow(ele) == 0, 10, median(ele$age_range_id, na.rm = T))
  mpnp2_long$age_range[i] <- ifelse(nrow(ele) == 0, 0, max(ele$age_range_id, na.rm = T) - min(ele$age_range_id, na.rm = T))
  rm(ele)
}

# take rounded-down-median
mpnp2_long$age_mid_round <- floor(mpnp2_long$age_mid)
mpnp2_long <- mpnp2_long[mpnp2_long$age_mid_round < 10 , c(3,5,6,7,32:34)]

# identify individuals present in both
length(unique(mpnp1_long$elephant_id)) ; sort(unique(mpnp1_long$elephant_id))
length(unique(mpnp2_long$elephant_id)) ; sort(unique(mpnp2_long$elephant_id))
estimated <- mpnp2_long[mpnp2_long$elephant_id %in% mpnp1_long$elephant_id,]
sort(unique(estimated$elephant_id))
missing <- anti_join(mpnp2_long, estimated)
sort(unique(missing$elephant_id))

# create a variable that randomly selects another elephant from time window 1 which shares an age category
missing_elephants <- missing[,c('elephant_id','age_mid_round')] %>% distinct()
sort(unique(missing_elephants$age_mid_round)) # 3:8
missing_elephants$shared_age <- NA
for(age_cat in 3:8){
  missing_ids <- missing_elephants[missing_elephants$age_mid_round == age_cat,]  # select all elephants in age category i that have not had their age probability distribution estimated in first time window 
  missing_elephants <- anti_join(missing_elephants, missing_ids)           # remove these elephants from current list of elephants lacking an age probability distribution
  age_mates <- mpnp1_males[mpnp1_males$age_mid_round == age_cat,]                # select males who WERE in age category i during first time window
  
  if(nrow(missing_ids) <= nrow(age_mates)){                                # where possible, draw a different random elephant for every missing individual from the pooled of males who shared that age distribution
    missing_ids$shared_age <- sample(x = age_mates$elephant_id, size = nrow(missing_ids), replace = F)
  } else {                                                                 # where not possible to always draw a different one, select entirely at random
    missing_ids$shared_age <- sample(x = age_mates$elephant_id, size = nrow(missing_ids), replace = T)
  }
  
  missing_elephants <- rbind(missing_elephants, missing_ids)               # return elephants from category i to list of those lacking an age 
  rm(missing_ids, age_mates)                                               # remove that age category
}
table(missing_elephants$shared_age)

# recombine data for second time window to include all elephants, present in MPNP1 or not
missing <- left_join(missing, missing_elephants, by = c('elephant_id','age_mid_round'))
head(missing)
head(estimated)
estimated$shared_age <- estimated$elephant_id
mpnp2_long <- rbind(estimated, missing)

# create data distribution for all 
colnames(mpnp1)                                                                           # check that mpnp1 data are named by the individual
mpnp2_males <- mpnp2_long[, c('elephant_id','age_mid_round','shared_age')] %>% distinct() # identify all individuals (and age-mate where appropriate)
mpnp2 <- mpnp1[,mpnp2_males$shared_age]                                                   # select mpnp1 columns matching mpnp2 age-mates
nrow(mpnp2_males) ; ncol(mpnp2)                                                           # check correctly selected

for(i in 1:nrow(mpnp2)){
  for(j in 1:ncol(mpnp2)){
    mpnp2[i,j] <- mpnp2[i,j] + time_diff
  }
}

# save output
saveRDS(mpnp2, '../data_processed/mpnp2_ageestimates_mcmcoutput.rds')

# clean environment
rm(estimated, missing, missing_elephants, mpnp2_long, mpnp2_males, mpnp2, age_cat, i, j)

## time window 3 ####
mpnp3_long <- mpnp_long[mpnp_long$date >= periods[3] & mpnp_long$date < periods[4],]

# identify category
mpnp3_long$age_range_id <- as.numeric(mpnp3_long$age_range_id)
mpnp3_long$age_mid <- NA
mpnp3_long$age_range <- NA
for(i in 1:nrow(mpnp3_long)){
  ele <- mpnp3_long[mpnp3_long$elephant_id == mpnp3_long$elephant_id[i],]
  ele <- ele[!is.na(ele$age_range_id),]
  ele <- ele[ele$age_range_id != 10,]
  mpnp3_long$age_mid[i] <- ifelse(nrow(ele) == 0, 10, median(ele$age_range_id, na.rm = T))
  mpnp3_long$age_range[i] <- ifelse(nrow(ele) == 0, 0, max(ele$age_range_id, na.rm = T) - min(ele$age_range_id, na.rm = T))
  rm(ele)
}

# take rounded-down-median
mpnp3_long$age_mid_round <- floor(mpnp3_long$age_mid)
mpnp3_long <- mpnp3_long[mpnp3_long$age_mid_round < 10 , c(3,5,6,7,32:34)]

# identify individuals present in both
length(unique(mpnp1_long$elephant_id)) ; length(unique(mpnp3_long$elephant_id))
estimated <- mpnp3_long[mpnp3_long$elephant_id %in% mpnp1_long$elephant_id,]
sort(unique(estimated$elephant_id))
missing <- anti_join(mpnp3_long, estimated)
sort(unique(missing$elephant_id))

# create a variable that randomly selects another elephant from time window 1 which shares an age category
missing_elephants <- missing[,c('elephant_id','age_mid_round')] %>% distinct()
sort(unique(missing_elephants$age_mid_round)) # 4:8
missing_elephants$shared_age <- NA
for(age_cat in 4:8){
  missing_ids <- missing_elephants[missing_elephants$age_mid_round == age_cat,]  # select all elephants in age category i that have not had their age probability distribution estimated in first time window 
  missing_elephants <- anti_join(missing_elephants, missing_ids)           # remove these elephants from current list of elephants lacking an age probability distribution
  age_mates <- mpnp1_males[mpnp1_males$age_mid_round == age_cat,]                # select males who WERE in age category i during first time window
  
  if(nrow(missing_ids) <= nrow(age_mates)){                                # where possible, draw a different random elephant for every missing individual from the pooled of males who shared that age distribution
    missing_ids$shared_age <- sample(x = age_mates$elephant_id, size = nrow(missing_ids), replace = F)
  } else {                                                                 # where not possible to always draw a different one, select entirely at random
    missing_ids$shared_age <- sample(x = age_mates$elephant_id, size = nrow(missing_ids), replace = T)
  }
  
  missing_elephants <- rbind(missing_elephants, missing_ids)               # return elephants from category i to list of those lacking an age 
  rm(missing_ids, age_mates)                                               # remove that age category
}
table(missing_elephants$shared_age)

# recombine data for second time window to include all elephants, present in MPNP1 or not
missing <- left_join(missing, missing_elephants, by = c('elephant_id','age_mid_round'))
head(missing)
head(estimated)
estimated$shared_age <- estimated$elephant_id
mpnp3_long <- rbind(estimated, missing)

# create data distribution for all 
mpnp3_males <- mpnp3_long[, c('elephant_id','age_mid_round','shared_age')] %>% distinct() # identify all individuals (and age-mate where appropriate)
mpnp3 <- mpnp1[,mpnp3_males$shared_age]                                                   # select mpnp1 columns matching mpnp3 age-mates
nrow(mpnp3_males) ; ncol(mpnp3)                                                           # check correctly selected

for(i in 1:nrow(mpnp3)){
  for(j in 1:ncol(mpnp3)){
    mpnp3[i,j] <- mpnp3[i,j] + time_diff
  }
}

# save output
saveRDS(mpnp3, '../data_processed/mpnp3_ageestimates_mcmcoutput.rds')

# clean environment
rm(estimated, missing, missing_elephants, mpnp3_long, mpnp3_males, mpnp3, age_cat, i, j)

## time window 4 ####
mpnp4_long <- mpnp_long[mpnp_long$date >= periods[4] & mpnp_long$date < periods[5],]

# identify category
mpnp4_long$age_range_id <- as.numeric(mpnp4_long$age_range_id)
mpnp4_long$age_mid <- NA
mpnp4_long$age_range <- NA
for(i in 1:nrow(mpnp4_long)){
  ele <- mpnp4_long[mpnp4_long$elephant_id == mpnp4_long$elephant_id[i],]
  ele <- ele[!is.na(ele$age_range_id),]
  ele <- ele[ele$age_range_id != 10,]
  mpnp4_long$age_mid[i] <- ifelse(nrow(ele) == 0, 10, median(ele$age_range_id, na.rm = T))
  mpnp4_long$age_range[i] <- ifelse(nrow(ele) == 0, 0, max(ele$age_range_id, na.rm = T) - min(ele$age_range_id, na.rm = T))
  rm(ele)
}

# take rounded-down-median
mpnp4_long$age_mid_round <- floor(mpnp4_long$age_mid)
mpnp4_long <- mpnp4_long[mpnp4_long$age_mid_round < 10 , c(3,5,6,7,32:34)]

# identify individuals present in both
length(unique(mpnp1_long$elephant_id)) ; length(unique(mpnp4_long$elephant_id))
estimated <- mpnp4_long[mpnp4_long$elephant_id %in% mpnp1_long$elephant_id,]
sort(unique(estimated$elephant_id))
missing <- anti_join(mpnp4_long, estimated)
sort(unique(missing$elephant_id))

# create a variable that randomly selects another elephant from time window 1 which shares an age category
missing_elephants <- missing[,c('elephant_id','age_mid_round')] %>% distinct()
sort(unique(missing_elephants$age_mid_round)) # 4:8
missing_elephants$shared_age <- NA
for(age_cat in 4:8){
  missing_ids <- missing_elephants[missing_elephants$age_mid_round == age_cat,]  # select all elephants in age category i that have not had their age probability distribution estimated in first time window 
  missing_elephants <- anti_join(missing_elephants, missing_ids)           # remove these elephants from current list of elephants lacking an age probability distribution
  age_mates <- mpnp1_males[mpnp1_males$age_mid_round == age_cat,]                # select males who WERE in age category i during first time window
  
  if(nrow(missing_ids) <= nrow(age_mates)){                                # where possible, draw a different random elephant for every missing individual from the pooled of males who shared that age distribution
    missing_ids$shared_age <- sample(x = age_mates$elephant_id, size = nrow(missing_ids), replace = F)
  } else {                                                                 # where not possible to always draw a different one, select entirely at random
    missing_ids$shared_age <- sample(x = age_mates$elephant_id, size = nrow(missing_ids), replace = T)
  }
  
  missing_elephants <- rbind(missing_elephants, missing_ids)               # return elephants from category i to list of those lacking an age 
  rm(missing_ids, age_mates)                                               # remove that age category
}
table(missing_elephants$shared_age)

# recombine data for second time window to include all elephants, present in MPNP1 or not
missing <- left_join(missing, missing_elephants, by = c('elephant_id','age_mid_round'))
head(missing)
head(estimated)
estimated$shared_age <- estimated$elephant_id
mpnp4_long <- rbind(estimated, missing)

# create data distribution for all 
mpnp4_males <- mpnp4_long[, c('elephant_id','age_mid_round','shared_age')] %>% distinct() # identify all individuals (and age-mate where appropriate)
mpnp4 <- mpnp1[,mpnp4_males$shared_age]                                                   # select mpnp1 columns matching mpnp4 age-mates
nrow(mpnp4_males) ; ncol(mpnp4)                                                           # check correctly selected

for(i in 1:nrow(mpnp4)){
  for(j in 1:ncol(mpnp4)){
    mpnp4[i,j] <- mpnp4[i,j] + time_diff
  }
}

# save output
saveRDS(mpnp4, '../data_processed/mpnp4_ageestimates_mcmcoutput.rds')

# clean environment
rm(estimated, missing, missing_elephants, mpnp4_long, mpnp4_males, mpnp4, age_cat, i, j)

## time window 5 ####
mpnp5_long <- mpnp_long[mpnp_long$date >= periods[5] & mpnp_long$date < periods[6],]

# identify category
mpnp5_long$age_range_id <- as.numeric(mpnp5_long$age_range_id)
mpnp5_long$age_mid <- NA
mpnp5_long$age_range <- NA
for(i in 1:nrow(mpnp5_long)){
  ele <- mpnp5_long[mpnp5_long$elephant_id == mpnp5_long$elephant_id[i],]
  ele <- ele[!is.na(ele$age_range_id),]
  ele <- ele[ele$age_range_id != 10,]
  mpnp5_long$age_mid[i] <- ifelse(nrow(ele) == 0, 10, median(ele$age_range_id, na.rm = T))
  mpnp5_long$age_range[i] <- ifelse(nrow(ele) == 0, 0, max(ele$age_range_id, na.rm = T) - min(ele$age_range_id, na.rm = T))
  rm(ele)
}

# take rounded-down-median
mpnp5_long$age_mid_round <- floor(mpnp5_long$age_mid)
mpnp5_long <- mpnp5_long[mpnp5_long$age_mid_round < 10 , c(3,5,6,7,32:34)]

# identify individuals present in both
length(unique(mpnp1_long$elephant_id)) ; length(unique(mpnp5_long$elephant_id))
estimated <- mpnp5_long[mpnp5_long$elephant_id %in% mpnp1_long$elephant_id,]
sort(unique(estimated$elephant_id))
missing <- anti_join(mpnp5_long, estimated)
sort(unique(missing$elephant_id))

# create a variable that randomly selects another elephant from time window 1 which shares an age category
missing_elephants <- missing[,c('elephant_id','age_mid_round')] %>% distinct()
sort(unique(missing_elephants$age_mid_round)) # 4:7
missing_elephants$shared_age <- NA
for(age_cat in 4:7){
  missing_ids <- missing_elephants[missing_elephants$age_mid_round == age_cat,]  # select all elephants in age category i that have not had their age probability distribution estimated in first time window 
  missing_elephants <- anti_join(missing_elephants, missing_ids)           # remove these elephants from current list of elephants lacking an age probability distribution
  age_mates <- mpnp1_males[mpnp1_males$age_mid_round == age_cat,]                # select males who WERE in age category i during first time window
  
  if(nrow(missing_ids) <= nrow(age_mates)){                                # where possible, draw a different random elephant for every missing individual from the pooled of males who shared that age distribution
    missing_ids$shared_age <- sample(x = age_mates$elephant_id, size = nrow(missing_ids), replace = F)
  } else {                                                                 # where not possible to always draw a different one, select entirely at random
    missing_ids$shared_age <- sample(x = age_mates$elephant_id, size = nrow(missing_ids), replace = T)
  }
  
  missing_elephants <- rbind(missing_elephants, missing_ids)               # return elephants from category i to list of those lacking an age 
  rm(missing_ids, age_mates)                                               # remove that age category
}
table(missing_elephants$shared_age)

# recombine data for second time window to include all elephants, present in MPNP1 or not
missing <- left_join(missing, missing_elephants, by = c('elephant_id','age_mid_round'))
head(missing)
head(estimated)
estimated$shared_age <- estimated$elephant_id
mpnp5_long <- rbind(estimated, missing)

# create data distribution for all 
mpnp5_males <- mpnp5_long[, c('elephant_id','age_mid_round','shared_age')] %>% distinct() # identify all individuals (and age-mate where appropriate)
mpnp5 <- mpnp1[,mpnp5_males$shared_age]                                                   # select mpnp1 columns matching mpnp5 age-mates
nrow(mpnp5_males) ; ncol(mpnp5)                                                           # check correctly selected

for(i in 1:nrow(mpnp5)){
  for(j in 1:ncol(mpnp5)){
    mpnp5[i,j] <- mpnp5[i,j] + time_diff
  }
}

# save output
saveRDS(mpnp5, '../data_processed/mpnp5_ageestimates_mcmcoutput.rds')

# clean environment
rm(estimated, missing, missing_elephants, mpnp5_long, mpnp5_males, mpnp5, age_cat, i, j)
