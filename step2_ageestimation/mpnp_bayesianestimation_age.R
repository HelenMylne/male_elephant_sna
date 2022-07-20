#### information #####
# To estimate the ages of individuals in other populations, first come up with a survival curve for the Amboseli elephants which have good birth and death data. -- this identified the best curve as being the Gompertz-Bathtub.
# Next fit this to the MPNP dataset

#### load packages ####
library(tidyverse)
library(cmdstanr)
library(ggdist)
library(posterior)
library(bayesplot)

#### load model -- adapted from MOTNP to use MPNP thresholds ####
# read in Stan model to estimate ages based on Gompertz bathtub distribution from ANP
latent_age_ordinal_model <- cmdstan_model("models/age_estimation/mpnp_elephant_latent_age_ordinal_regression_hkm_22.07.07.stan")
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

#### load MPNP data ####
mpnp_long <- readxl::read_excel('data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211019.xlsx') %>%
  janitor::clean_names()

# identify period 1 only
periods <- seq(from = min(mpnp_long$date, na.rm = T),
               to = max(mpnp_long$date, na.rm = T),
               length.out = 7)
periods[7] <- periods[7]+1 # set to be one higher than the final higher otherwise it takes the last date and creates a whole new period
#as.numeric(max(mpnp_groups$date, na.rm = T)) - as.numeric(min(mpnp_groups$date, na.rm = T))
#mpnp1_groups <- mpnp_groups[mpnp_groups$date < periods[2],]

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
mpnp1_long$age_min <- NA
mpnp1_long$age_max <- NA
mpnp1_long$age_mid <- NA
mpnp1_long$age_mean <- NA
mpnp1_long$age_range <- NA
for(i in 1:nrow(mpnp1_long)){
  ele <- mpnp1_long[mpnp1_long$elephant_id == mpnp1_long$elephant_id[i],]
  ele <- ele[!is.na(ele$age_range_id),]
  ele <- ele[ele$age_range_id != 10,]
  mpnp1_long$age_min[i] <- ifelse(nrow(ele) == 0, 10, min(ele$age_range_id))
  mpnp1_long$age_max[i] <- ifelse(nrow(ele) == 0, 10, max(ele$age_range_id, na.rm = T))
  mpnp1_long$age_mid[i] <- ifelse(nrow(ele) == 0, 10, median(ele$age_range_id, na.rm = T))
  mpnp1_long$age_mean[i] <- ifelse(nrow(ele) == 0, 10, round(mean(ele$age_range_id, na.rm = T),0))
  mpnp1_long$age_range[i] <- ifelse(nrow(ele) == 0, 0, max(ele$age_range_id, na.rm = T) - min(ele$age_range_id, na.rm = T))
  rm(ele)
}

# take median category and then round down to nearest integer
mpnp1_long$age_mid_round <- floor(mpnp1_long$age_mid)
mpnp1_long <- mpnp1_long[mpnp1_long$age_mid_round < 10 , c(3,5,6,7,32:34)]

rm(mpnp_long, periods)

#### create data list ####
N_mpnp1 <- nrow(mpnp1_males)
K <- 9
mpnp1_ls <- list(
  N = N_mpnp1,
  K = K,
  age_category_index = mpnp1_males$age_cat_id)
hist(mpnp1_ls$age_category_index)
#hist(elephants_ls$age_category_index)

#### fit model to MPNP data ####
# Fit model with cmdstanr
age_mpnp1_fit <- latent_age_ordinal_model$sample(
  data = mpnp1_ls, 
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000
)

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
#mcmc_dens(true_ages)
true_ages <- true_ages[,1:N_mpnp1]

df <- as.data.frame(do.call(rbind, true_ages)) %>%
  mutate(age_cat = mpnp1_ls$age) %>% relocate(age_cat) %>%
  mutate(ID = mpnp1_males$id) %>% relocate(ID)

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

#### save output ####
saveRDS(true_ages, 'data_processed/mpnp1_age')
