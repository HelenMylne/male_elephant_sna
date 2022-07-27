#### information #####
# To estimate the ages of individuals in other populations, first come up with a survival curve for the Amboseli elephants which have good birth and death data.
# survival curve options = exponential/logistic/Gompertz/Weibull, where shape = simple/bathtub/Makeham
# follow through Mia's code to determine which curve will best fit Amboseli data. she ran it on a high performance cluster, so it’s set up to run in parallel, but the basta model can also be run in parallel on its own (just like Stan). we will need to implement it as a custom function in Stan for the final analysis, given that we need to define it in the prior anyway. 

#### load packages ####
library(BaSTA)
library(tidyverse)
library(cmdstanr)
library(ggdist)
library(posterior)
library(bayesplot)

#library(ggthemes)
#library(survival)
#library(tidybayes)
#library(data.table)

#### read in Amboseli data ####
sightings <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_sightings_updated.csv')
sightings$year <- lubridate::year(sightings$obs_date)
sightings <- sightings[,c('casename','year')]

males <- readxl::read_excel('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% 
  janitor::clean_names() %>% 
  select(casename, byr, dyr)
str(males)
males <- males[males$casename %in% sightings$casename,]

#### create basta data frame ####
# Create columns for sightings in each year
males$obs_1972 <- NA ; males$obs_1973 <- NA ; males$obs_1974 <- NA ; males$obs_1975 <- NA ; males$obs_1976 <- NA ; males$obs_1977 <- NA ; males$obs_1978 <- NA ; males$obs_1979 <- NA ; males$obs_1980 <- NA ; males$obs_1981 <- NA ; males$obs_1982 <- NA ; males$obs_1983 <- NA ; males$obs_1984 <- NA ; males$obs_1985 <- NA ; males$obs_1986 <- NA ; males$obs_1987 <- NA ; males$obs_1988 <- NA ; males$obs_1989 <- NA ; males$obs_1990 <- NA ; males$obs_1991 <- NA ; males$obs_1992 <- NA ; males$obs_1993 <- NA ; males$obs_1994 <- NA ; males$obs_1995 <- NA ; males$obs_1996 <- NA ; males$obs_1997 <- NA ; males$obs_1998 <- NA ; males$obs_1999 <- NA ; males$obs_2000 <- NA ; males$obs_2001 <- NA ; males$obs_2002 <- NA ; males$obs_2003 <- NA ; males$obs_2004 <- NA ; males$obs_2005 <- NA ; males$obs_2006 <- NA ; males$obs_2007 <- NA ; males$obs_2008 <- NA ; males$obs_2009 <- NA ; males$obs_2010 <- NA ; males$obs_2011 <- NA ; males$obs_2012 <- NA ; males$obs_2013 <- NA ; males$obs_2014 <- NA ; males$obs_2015 <- NA ; males$obs_2016 <- NA ; males$obs_2017 <- NA ; males$obs_2018 <- NA ; males$obs_2019 <- NA ; males$obs_2020 <- NA ; males$obs_2021 <- NA

# fill columns with 1 = sighting, 0 = no sighting
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

#### remove elephants with incorrect sightings ####
check_birth <- basta_data[c(638,648,665,682,686),] # elephants recorded before birth
check_death <- basta_data[c(3,5,7,13,15,16,17,44,54,55,57,61,65,67,71,82,87,89,92,95,104,106,107,110,113,117,126,127,131,132,134,135,137,145,150,154,158,161,178,181,186,189,190,207,209,218,228,242,246,252,254,256,263,288,291,294,323,335,338,365,378,410,411,429,653),]    # elephants recorded after death
basta_data <- anti_join(basta_data, check_birth, by = 'casename')
basta_data <- anti_join(basta_data, check_death, by = 'casename')
rm(check_birth, check_death)

#### add minimum age of independence to all birth years ####
basta_data$year_first_sighted <- NA
for(i in 1:nrow(basta_data)){ # loop through all rows and all years, identify first column with a "1" per row -- apparent "birth" year will be first sighting -1, as basta won't allow an individual to be sighted in the year it is born.
  for(j in 4:53){
    if(is.na(basta_data$year_first_sighted[i]) == TRUE){
      basta_data$year_first_sighted[i] <- ifelse(basta_data[i,j] == 1,
                                                 (j+1968)-1, # j+1968 = year of first independent sighting. -1 because won't allow sightings in "birth" year
                                                 basta_data$year_first_sighted[i]) # no sighting, leave as NA
    }
    else { basta_data$year_first_sighted[i] <- basta_data$year_first_sighted[i] }
  }
}

# look at ages each individual was first seen
basta_data$age_first_sighted <- basta_data$year_first_sighted - basta_data$byr
hist(basta_data$age_first_sighted, breaks = 50)

# remove all individuals first observed at <8 years old
table(basta_data$age_first_sighted)
basta_data <- basta_data[basta_data$age_first_sighted > 7,]

# set "birth year" to 8 years after birth -- 8 years because minimum age reported for independence in literature
basta_data$byr_apparent <- basta_data$byr + 8

# set birth and death years
basta_data$birth_year <- ifelse(basta_data$byr_apparent < 1972, 0, basta_data$byr_apparent)
basta_data$death_year <- ifelse(basta_data$dyr < 100,  0, basta_data$dyr)

# create data frame for basta
basta_anp <- basta_data[,c(1,58:59,4:53)]

#### prior predictive check ####
gompertz_bt <- function(a0, a1, c, b0, b1, age){
  gompertz <- exp(b0 + b1*age)
  bathtub <- exp(a0 - a1*age) + c + gompertz
  return(bathtub)
}

a0 <- -4
a1 <- 3
c <- 1
b0 <- -6
b1 <- 0.05
mortality <- gompertz_bt(a0 = a0, a1 = a1, c = c, b0 = b0 , b1 = b1, 1:50) ; plot(mortality)
to.plot <- data.frame(a0 = c(rep(a0, 50), rep(NA,4950)),
                      a1 = c(rep(a1, 50), rep(NA,4950)),
                      c  = c(rep(c,  50), rep(NA,4950)),
                      b0 = c(rep(b0, 50), rep(NA,4950)),
                      b1 = c(rep(b1, 50), rep(NA,4950)),
                      age = rep(1:50,100),
                      y = rep(NA,5000),
                      iter = rep(1:100, each = 50))
for(i in 1:99){ # populate data frame with curves indicating the probability of survival for 100 draws from the prior distribution
  to.plot$a0[(i*50+1):(i*50+50)] <- rep(rnorm(1, a0, 1),    50)
  to.plot$a1[(i*50+1):(i*50+50)] <- rep(rnorm(1, a1, 0.1),  50)
  to.plot$c[ (i*50+1):(i*50+50)] <- rep(rnorm(1, c,  0.4), 50)
  to.plot$b0[(i*50+1):(i*50+50)] <- rep(rnorm(1, b0, 2),    50)
  to.plot$b1[(i*50+1):(i*50+50)] <- rep(rnorm(1, b1, 0.01), 50)
}
to.plot$y  <- gompertz_bt(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1,
                          age = to.plot$age)
to.plot$y_std <- to.plot$y / max(to.plot$y)
plot(NULL, xlab = 'age', ylab = 'mortality probability', xlim = c(0,50), ylim = c(0,1), las = 1)
for(i in 1:100){
  data <- to.plot[to.plot$iter == i,]
  lines(data$y_std ~ data$age, col = rgb(0,0,1,0.2))
}

pM.Gb <-  c(a0,a1,c,b0,b1)        # prior means
pSD.Gb <- c(1,0.1,0.4,2,0.05)     # prior standard deviations

#### set iteration arguments ####
start_year <- 1972
end_year <- 2021
iter <- 50000      # number of iterations
sim <- 4           # number of chains
warmup <- iter/2   # number of iterations to use during warmup/burn-in period
thin <- 100        # thinning number

#### run model ####
Go.Bt <- basta(basta_anp, studyStart = start_year,
               studyEnd = end_year, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM.Gb, thetaPriorSd = pSD.Gb,
               model = "GO", shape = c("bathtub"))
summary(Go.Bt)
#     Jump.sds Prior.means Prior.sds
#a0 3.31117128       -4.00      1.00
#a1 0.49077852        3.00      0.10
#c  0.01031378        1.00      0.40
#b0 0.32489720       -6.00      2.00
#b1 0.01175845        0.05      0.05
#KLDC was not calculated due to lack of convergence, or because covariates were not included in the model.
#         Estimate   StdErr Lower95%CI Upper95%CI SerAutocor UpdateRate PotScaleReduc
#a0      -5.12213 0.718990  -6.601708   -3.83542  0.0005584     0.2544        0.9991
#a1       3.00076 0.105410   2.796268    3.21278 -0.0103880     0.2445        1.0021
#c        0.02583 0.006686   0.009942    0.03669  0.6750486     0.2547        1.0059
#b0      -5.06816 0.566933  -6.237723   -3.98137  0.7651080     0.2483        1.0094
#b1       0.09482 0.017827   0.061215    0.13199  0.7687475     0.2525        1.0110
#pi.1972  0.78884 0.004540   0.779725    0.79752 -0.0142274     1.0000        1.0007
# Appropriate convergence reached for all parameters.
plot(Go.Bt)                      # chains, esp. b0 and b1, a little wander-y
plot(Go.Bt, plot.trace = FALSE)  #very similar to above

#### simulate from output ####
a0.post <- rnorm(100, Go.Bt$coefficients[1,1], Go.Bt$coefficients[1,2])
a1.post <- rnorm(100, Go.Bt$coefficients[2,1], Go.Bt$coefficients[2,2])
c.post  <- rnorm(100, Go.Bt$coefficients[3,1], Go.Bt$coefficients[3,2])
b0.post <- rnorm(100, Go.Bt$coefficients[4,1], Go.Bt$coefficients[4,2])
b1.post <- rnorm(100, Go.Bt$coefficients[5,1], Go.Bt$coefficients[5,2])

gompertz_bt <- function(a0, a1, c, b0, b1, age){
  gompertz <- exp(b0 + b1*age)
  bathtub <- exp(a0 - a1*age) + c + gompertz
  return(bathtub)
}

mortality <- gompertz_bt(a0 = Go.Bt$coefficients[1,1],
                         a1 = Go.Bt$coefficients[2,1],
                         c  = Go.Bt$coefficients[3,1],
                         b0 = Go.Bt$coefficients[4,1],
                         b1 = Go.Bt$coefficients[5,1], 1:50) ; plot(mortality)

to.plot <- data.frame(a0 = rep(a0.post, each = 50),
                      a1 = rep(a1.post, each = 50),
                      c  = rep(c.post,  each = 50),
                      b0 = rep(b0.post, each = 50),
                      b1 = rep(b1.post, each = 50),
                      age = rep(1:50,100),
                      y = rep(NA,5000),
                      iter = rep(1:100, each = 50))
to.plot$y  <- gompertz_bt(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1,
                          age = to.plot$age)
plot(NULL, xlab = 'age', ylab = 'mortality probability', xlim = c(0,50), ylim = c(0,1), las = 1)
for(i in 1:100){
  data <- to.plot[to.plot$iter == i,]
  lines(data$y ~ data$age, col = rgb(0,0,1,0.2))
}

# draw new output curve
mortality <- gompertz_bt(a0 = Go.Bt$coefficients[1,1],
                         a1 = Go.Bt$coefficients[2,1],
                         c  = Go.Bt$coefficients[3,1],
                         b0 = Go.Bt$coefficients[4,1],
                         b1 = Go.Bt$coefficients[5,1], 1:70)
probs <- 1 - (mortality/max(mortality))
plot(probs)

# create a fictional population with ages selected from that distribution
N <- 700
elephants_ls <- list(
  N = N,
  age = sample(1:70, N, replace = TRUE, prob = probs)
)
hist(elephants_ls$age)

basta_data$age_death <- ifelse(basta_data$death_year == 0, NA, basta_data$death_year - basta_data$byr)
hist(basta_data$age_death)

#### simulate process of assigning age categories to elephants -- THIS WORKS PERFECTLY, just missing double 0-5 threshold ####
# draw new output curves
max_age <- 60
mortality <- gompertz_bt(a0 = Go.Bt$coefficients[1,1],
                         a1 = Go.Bt$coefficients[2,1],
                         c  = Go.Bt$coefficients[3,1],
                         b0 = Go.Bt$coefficients[4,1],
                         b1 = Go.Bt$coefficients[5,1], 1:max_age)
plot(mortality)
probs <- 1 - (mortality/max(mortality))
plot(probs)

# create a fictional population with ages selected from that distribution
N <- 200
K <- 7
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
                                          function(x) which.max(x < c(10, 15, 20, 25, 40, 60, Inf)))
hist(elephants_ls$age_category_index)

# look at actual age vs biologist assigned age: chance of being mis-classified is lower if actual age is in middle of age category, and biased towards lower end of class.
data.frame(elephants_ls) %>%
  ggplot(aes(x=age, y=age_guess, col=factor(age_category_index))) +
    geom_point(size=4,alpha=0.6) +
    geom_vline(xintercept=c(15, 20, 25, 40, 60), col=factor(3:7), linetype="dashed", alpha=0.6) +
    theme_minimal() + 
    xlab("Assigned age") + ylab("Age")

# read in Stan model to estimate ages based on Gompertz bathtub distribution from ANP
latent_age_ordinal_model <- cmdstan_model("models/motnp_age_ordinal_regression.stan")

#Fit model with cmdstanr
age_estimation_fit <- latent_age_ordinal_model$sample(
  data = elephants_ls, 
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000
)

#Examine the estimates. We can plot the estimated ages against the biologist assigned ages and ...
age_est_mat <- age_estimation_fit$summary()[202:401, ]
summary(age_est_mat)
hist(age_est_mat$mean)
hist(age_est_mat$rhat, breaks = 20)

plot_data <- data.frame(age = elephants_ls$age,       # Biologists original age est
                        model_age = age_est_mat$mean) # Mean modelled age

plot_data %>%
  ggplot(aes(x=factor(age), y=model_age)) +
    geom_point(size=4,col = 'blue', alpha=0.6) +
    geom_vline(xintercept=c(10, 15, 20, 25, 40, 60), linetype="dashed", alpha=0.6) +
    geom_hline(yintercept=c(10, 15, 20, 25, 40, 60), linetype="dashed", alpha=0.6) +
    geom_abline(slope = 1, intercept = 0)+
    scale_y_continuous(limits = c(0,60))+
    theme_minimal() + 
    xlab("Assigned age") + ylab("Modelled age")

#### simulate process of assigning age categories to elephants -- THIS WORKS PERFECTLY ####
# draw new output curves
max_age <- 60
mortality <- gompertz_bt(a0 = Go.Bt$coefficients[1,1],
                         a1 = Go.Bt$coefficients[2,1],
                         c  = Go.Bt$coefficients[3,1],
                         b0 = Go.Bt$coefficients[4,1],
                         b1 = Go.Bt$coefficients[5,1], 1:max_age)
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

# read in Stan model to estimate ages based on Gompertz bathtub distribution from ANP
latent_age_ordinal_model <- cmdstan_model("models/motnp_age_ordinal_regression.stan")

#Fit model with cmdstanr
age_estimation_fit <- latent_age_ordinal_model$sample(
  data = elephants_ls, 
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000
)

#Examine the estimates. We can plot the estimated ages against the biologist assigned ages and ...
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

#### run on MOTNP elephants -- THIS WORKS PERFECTLY ####
motnp_males <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_elenodes.csv') %>% 
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
N_motnp <- nrow(motnp_males)
K <- 8
motnp_ls <- list(
  N = N_motnp,
  K = K,
  age_category_index = motnp_males$age_cat_id)
hist(motnp_ls$age_category_index)
#hist(elephants_ls$age_category_index)

# read in Stan model to estimate ages based on Gompertz bathtub distribution from ANP
latent_age_ordinal_model <- cmdstan_model("models/motnp_age_ordinal_regression.stan")

#Fit model with cmdstanr
age_motnp_fit <- latent_age_ordinal_model$sample(
  data = motnp_ls, 
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000
)

#Examine the estimates. We can plot the estimated ages against the biologist assigned ages and ...
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
  geom_vline(xintercept=c(5,10,15, 20, 25, 40, 60), linetype="dashed", alpha=0.6) +
  geom_hline(yintercept=c(5,10,15, 20, 25, 40, 60), linetype="dashed", alpha=0.6) +
  #geom_abline(slope = 1, intercept = 0)+
  scale_y_continuous(limits = c(0,60))+
  theme_minimal() + 
  xlab("Assigned age") + ylab("Modelled age")

# Now an equivalent posterior predictive plot from draws from the distribution to show uncertainty around the mean age (then need to do it with full uncertainty).
true_ages <- age_motnp_fit$draws("true_age", format="df")
#mcmc_dens(true_ages)
true_ages <- true_ages[1:100,1:N_motnp]

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
