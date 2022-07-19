################# Information ##################
# To estimate the ages of individuals in other populations, first come up with a survival curve for the Amboseli elephants which have good birth and death data.
# survival curve options = exponential/logistic/Gompertz/Weibull, where shape = simple/bathtub/Makeham
# follow through Mia's code to determine which curve will best fit Amboseli data. she ran it on a high performance cluster, so it’s set up to run in parallel, but the basta model can also be run in parallel on its own (just like Stan). we will need to implement it as a custom function in Stan for the final analysis, given that we need to define it in the prior anyway. 

################# Load packages #################
library(BaSTA)
library(tidyverse)
library(cmdstanr)
library(rstan)

################# Run Gompertz Bathtub basta model to obtain prior estimates #################
sightings <- read_csv('data_processed/anp_sightings_updated_22.06.22.csv')
sightings$year <- lubridate::year(sightings$obs_date)
sightings <- sightings[,c('casename','year')]

males <- readxl::read_excel('data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% 
  janitor::clean_names() %>% 
  select(casename, byr, dyr)
str(males)
males <- males[males$casename %in% sightings$casename,]

# Create columns for sightings in each year
males$obs_1972 <- NA ; males$obs_1973 <- NA ; males$obs_1974 <- NA ; males$obs_1975 <- NA ; males$obs_1976 <- NA ; males$obs_1977 <- NA ; males$obs_1978 <- NA ; males$obs_1979 <- NA ; males$obs_1980 <- NA ; males$obs_1981 <- NA ; males$obs_1982 <- NA ; males$obs_1983 <- NA ; males$obs_1984 <- NA ; males$obs_1985 <- NA ; males$obs_1986 <- NA ; males$obs_1987 <- NA ; males$obs_1988 <- NA ; males$obs_1989 <- NA ; males$obs_1990 <- NA ; males$obs_1991 <- NA ; males$obs_1992 <- NA ; males$obs_1993 <- NA ; males$obs_1994 <- NA ; males$obs_1995 <- NA ; males$obs_1996 <- NA ; males$obs_1997 <- NA ; males$obs_1998 <- NA ; males$obs_1999 <- NA ; males$obs_2000 <- NA ; males$obs_2001 <- NA ; males$obs_2002 <- NA ; males$obs_2003 <- NA ; males$obs_2004 <- NA ; males$obs_2005 <- NA ; males$obs_2006 <- NA ; males$obs_2007 <- NA ; males$obs_2008 <- NA ; males$obs_2009 <- NA ; males$obs_2010 <- NA ; males$obs_2011 <- NA ; males$obs_2012 <- NA ; males$obs_2013 <- NA ; males$obs_2014 <- NA ; males$obs_2015 <- NA ; males$obs_2016 <- NA ; males$obs_2017 <- NA ; males$obs_2018 <- NA ; males$obs_2019 <- NA ; males$obs_2020 <- NA ; males$obs_2021 <- NA
(years <- 2022-1972)
for(i in 1:nrow(males)){ # 1 = sighting in that year, 0 = no sighting
  for(j in 1:years){
    casename <- sightings[sightings$casename == males$casename[i],]
    males[i,j+3] <- ifelse(length(which(casename$year == j+1971)) > 0, 1, 0)
    rm(casename)
  }
}

# prep data
basta_data <- males
basta_data$count_years <- rowSums(basta_data[4:53])
summary(basta_data$count_years) # no individuals never observed, as it should be. No NAs.

# set birth and death years
basta_data$birth_year <- ifelse(basta_data$byr < 1972, 0, basta_data$byr)
basta_data$death_year <- ifelse(basta_data$dyr < 100,  0, basta_data$dyr)
basta_data <- basta_data[,c(1,55,56,4:53,2:3,54)]

# remove elephants with incorrect sightings
check_birth <- basta_data[c(638,648,665,682,686),]
check_death <- basta_data[c(3,5,7,13,15,16,17,44,54,55,57,61,65,67,71,82,87,89,92,95,104,106,107,110,113,117,126,127,131,132,134,135,137,145,150,154,158,161,178,181,186,189,190,207,209,218,228,242,246,252,254,256,263,288,291,294,323,335,338,365,378,410,411,429,653),]
basta_data <- anti_join(basta_data, check_birth, by = 'casename')
basta_data <- anti_join(basta_data, check_death, by = 'casename')
rm(check_birth, check_death)

# convert birth years to year first sighted
basta_data$year_first_sighted <- NA
for(i in 1:nrow(basta_data)){
  for(j in 4:53){
    if(is.na(basta_data$year_first_sighted[i]) == TRUE){
      basta_data$year_first_sighted[i] <- ifelse(basta_data[i,j] == 1,
                                                 (j+1968)-1, # j+1968 = year of first independent sighting. -1 because won't allow sightings in "birth" year
                                                 basta_data$year_first_sighted[i]) # no sighting, leave as NA
    }
    else { basta_data$year_first_sighted[i] <- basta_data$year_first_sighted[i] }
  }
}
basta_data <- basta_data[,c(1,57,2,54,3:53,55:56)]
basta_data$byr_apparent <- ifelse(basta_data$byr < 1962, 0,
                                  ifelse(basta_data$year_first_sighted == 1971, 0,
                                         basta_data$year_first_sighted))
basta_data <- basta_data[,c(1,58,5:55,2:4,56,57)]
basta_data <- basta_data[,c(1:53)]

# Set the arguments for the iterations
sim <- 4
iter <- 50000
warmup <- iter/2
thin <- 100
start_year <- 1972
end_year <- 2021

# run model
a0 <- -4
a1 <- 3
c <- 1
b0 <- -6
b1 <- 0.05
pM.Gb <-  c(a0,a1,c,b0,b1)        # prior means
pSD.Gb <- c(1,0.1,0.4,2,0.05)     # prior standard deviations
Go.Bt <- basta(basta_data, studyStart = start_year,
               studyEnd = end_year, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM.Gb, thetaPriorSd = pSD.Gb,
               model = "GO", shape = c("bathtub"))
summary(Go.Bt)                             # posterior values to add to Stan model
plot(Go.Bt, las = 1)                       # Plot traces -- looks a bit better than simple, but still a bit wandery
plot(Go.Bt, plot.trace = FALSE, las = 1)   # Plot survival and mortality curves -- better than simple
par(mfrow = c(1,1))

# clean environment
rm(list = ls())
################# Run on MOTNP data #################
# read in Stan model
anp_age_ord_mod <- cmdstan_model("models/age_estimation/male_age_ordreg_22.06.27.stan")

# read in male elephant data
motnp <- read_csv('data_processed/motnp_elenodes_22.01.13.csv') %>% 
  filter(sex == 'M')

sort(unique(motnp$age_category))
motnp$age_cat <- ifelse(motnp$age_category == '40+', '7',
                        ifelse(motnp$age_category == '25-40','6',
                               ifelse(motnp$age_category == '20-25', '5',
                                      ifelse(motnp$age_category == '15-19', '4',
                                             ifelse(motnp$age_category == '10-15', '3',
                                                    ifelse(motnp$age_category == '8-9', '2',
                                                           ifelse(motnp$age_category == '7-8', '2',
                                                                  ifelse(motnp$age_category == '6-7', '2',
                                                                         ifelse(motnp$age_category == '5-6', '2',
                                                                                '1')))))))))
motnp$age_cat

independents <- motnp[motnp$age_cat > 2,]

# create data list
N <- nrow(motnp)                      # Number of individuals -- should we dropping anything <10 years old?
K <- 7                                # Number of age bin categories -- should this be 5 and we just drop anything in 1/2?
elephants_ls <- list(
  N = N,
  K = K,
  age_category_index = motnp$age_cat) # Age category

# fit Stan model to estimate ages from the categories
age_estimation_fit <- anp_age_ord_mod$sample(
  data = elephants_ls, 
  chains = 4, 
  parallel_chains = 4,
  iter_sampling = 2000
)

# examine estimates
output <- age_estimation_fit$summary()
observed_age_std <- output[2:246,]
observed_age <- output[261:505,]
true_age <- output[247,]
sigma_age <- output[248,]
gompertz <- output[249:253,]
independence_age <- output[254,]

plot(observed_age$mean,observed_age_std$mean)

plot_data <- data.frame(age_cat = elephants_ls$age,     # Biologists original age est
                        model_age = observed_age$mean)  # Mean modelled age
plot_data %>%
  ggplot(aes(x = factor(age_cat), y = model_age)) +
  geom_point(size = 4, alpha = 0.6) +
  geom_hline(yintercept = c(5,10,15,20,25,40), linetype = "dashed", alpha = 0.6) +
  theme_minimal() + 
  xlab("Assigned age") + ylab("Modelled age") # plot the estimated ages against the biologist assigned ages

# Now an equivalent posterior predictive plot from draws from the distribution to show uncertainty around the mean age (then need to do it with full uncertainty).
true_ages <- age_estimation_fit$draws("true_age", format = "df")
true_ages <- true_ages[1:1000,1:N]
means <- apply(true_ages, 2, mean)
stdevs <- apply(true_ages, 2, sd)

summary(means)
summary(stdevs)

df <- as.data.frame(do.call(rbind, true_ages)) %>%
  mutate(true_age = elephants_ls$age) %>%
  relocate(true_age) %>%
  mutate(id = 1:ncol(true_ages)) %>% # changed nrow to ncol, else it had to be the same number of draws as elephants
  relocate(id)
df <- df %>% pivot_longer(cols = 3:1002) %>% dplyr::select(-name) %>% distinct()
df %>% ggplot(aes(x = jitter(as.numeric(true_age)), y = value, group = factor(id))) +
  geom_point(size=2,alpha=0.1) +
  #geom_vline(xintercept=c(5,10,15,20,25,40), linetype="dashed", alpha=0.6) +
  geom_hline(yintercept=c(5,10,15,20,25,40), linetype="dashed", alpha=0.6) +
  theme_bw() + 
  xlab("Assigned age") + ylab("Modelled age")
summary(df$value[which(df$true_age == 7)])
quantile(df$value[which(df$true_age == 7)], seq(0,1,0.01))

