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
sightings <- read_csv('../data_processed/step1_dataprocessing/anp_sightings_updated.csv') # updated sightings = new data sheet sent over by Vicki Fishlock with corrections to old observations
sightings$year <- lubridate::year(sightings$obs_date)      # fix date variable
sightings <- sightings[,c('casename','year')]              # select only ID and year of sighting variables

males <- readxl::read_excel('../data_raw/Raw_ATE_Males_Lee220121.xlsx') %>%                                                # information about individual males
  janitor::clean_names() %>%                               # clean up 
  select(casename, byr, dyr)                               # select birth and death years for each individual
str(males)
males <- males[males$casename %in% sightings$casename,]    # filter for birth and death years only of males present in sightings data

#### create basta data frame ####
# Create columns for sightings in each year
males$obs_1972 <- NA ; males$obs_1973 <- NA ; males$obs_1974 <- NA ; males$obs_1975 <- NA ; males$obs_1976 <- NA ; males$obs_1977 <- NA ; males$obs_1978 <- NA ; males$obs_1979 <- NA ; males$obs_1980 <- NA ; males$obs_1981 <- NA ; males$obs_1982 <- NA ; males$obs_1983 <- NA ; males$obs_1984 <- NA ; males$obs_1985 <- NA ; males$obs_1986 <- NA ; males$obs_1987 <- NA ; males$obs_1988 <- NA ; males$obs_1989 <- NA ; males$obs_1990 <- NA ; males$obs_1991 <- NA ; males$obs_1992 <- NA ; males$obs_1993 <- NA ; males$obs_1994 <- NA ; males$obs_1995 <- NA ; males$obs_1996 <- NA ; males$obs_1997 <- NA ; males$obs_1998 <- NA ; males$obs_1999 <- NA ; males$obs_2000 <- NA ; males$obs_2001 <- NA ; males$obs_2002 <- NA ; males$obs_2003 <- NA ; males$obs_2004 <- NA ; males$obs_2005 <- NA ; males$obs_2006 <- NA ; males$obs_2007 <- NA ; males$obs_2008 <- NA ; males$obs_2009 <- NA ; males$obs_2010 <- NA ; males$obs_2011 <- NA ; males$obs_2012 <- NA ; males$obs_2013 <- NA ; males$obs_2014 <- NA ; males$obs_2015 <- NA ; males$obs_2016 <- NA ; males$obs_2017 <- NA ; males$obs_2018 <- NA ; males$obs_2019 <- NA ; males$obs_2020 <- NA ; males$obs_2021 <- NA

# fill columns with 1 = sighting, 0 = no sighting
(years <- 2022-1972)         # number of years of study = 50
for(i in 1:nrow(males)){     # fill year columns with 1 = sighting in that year, 0 = no sighting
  for(j in 1:years){
    casename <- sightings[sightings$casename == males$casename[i],]
    males[i,j+3] <- ifelse(length(which(casename$year == j+1971)) > 0, 1, 0)
    rm(casename)
  }
}

# prep data
basta_data <- males                                   # new data frame
basta_data$count_years <- rowSums(basta_data[4:53])   # count total number of years in which an individual was observed
summary(basta_data$count_years)                       # no individuals never observed, as it should be. No NAs.

#### remove elephants with incorrect sightings ####
check_birth <- basta_data[c(638,648,665,682,686),]                          # elephants recorded before birth
check_death <- basta_data[c(3,5,7,13,15,16,17,44,54,55,57,61,65,67,71,82,87,89,92,95,104,106,107,110,113,117,126,127,131,132,
                            134,135,137,145,150,154,158,161,178,181,186,189,190,207,209,218,228,242,246,252,254,256,263,288,
                            291,294,323,335,338,365,378,410,411,429,653),]  # elephants recorded after death
basta_data <- anti_join(basta_data, check_birth, by = 'casename')           # remove individuals recorded as present before they were born
basta_data <- anti_join(basta_data, check_death, by = 'casename')           # remove individuals recorded as present after they died
rm(check_birth, check_death)                                                # clean environment

#### add minimum age of independence to all birth years ####
basta_data$year_first_sighted <- NA     # new variable for the first year that an elephant was recorded as sighted separate from family
for(i in 1:nrow(basta_data)){           # loop through all rows and all years, identify first column with a "1" per row -- apparent "birth" year will be first sighting -1, as basta won't allow an individual to be sighted in the year it is born.
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
basta_data$age_first_sighted <- basta_data$year_first_sighted - basta_data$byr  # age first recorded as independent of family
hist(basta_data$age_first_sighted, breaks = 50)                                 # histogram

# remove all individuals first observed at <8 years old
table(basta_data$age_first_sighted)                            # check ages
basta_data <- basta_data[basta_data$age_first_sighted > 7,]    # assume any dispersing below 8 is likely an error (literature)

# set "birth year" to 8 years after birth -- 8 years because minimum age reported for independence in literature
basta_data$byr_apparent <- basta_data$byr + 8                  # alter the apparent birth year to byr+8 -- survival curve starting from first age at which male elephants start to become likely to disperse

# set birth and death years
basta_data$birth_year <- ifelse(basta_data$byr_apparent < 1972, 0, basta_data$byr_apparent) # anything over 8 years old at start of study set as unknown birth year
basta_data$death_year <- ifelse(basta_data$dyr < 100,  0, basta_data$dyr)                   # anything with unknown death year set to still alive

# create data frame for basta
basta_anp <- basta_data[,c(1,58:59,4:53)]            # select only necessary variables

#### prior predictive check ####
gompertz_bt <- function(a0, a1, c, b0, b1, age){          # create function for visualising gompertz bathtub distribution
  gompertz <- exp(b0 + b1*age)
  bathtub <- exp(a0 - a1*age) + c + gompertz
  return(bathtub)
}

a0 <- -5 ; a1 <- 1 ; c <- 2 ; b0 <- -4 ; b1 <- 0       # set prior values
max_age <- 70
mortality <- gompertz_bt(a0 = a0, a1 = a1, c = c, b0 = b0 , b1 = b1, 1:max_age) ; plot(mortality/max(mortality))  # plot distribution
to.plot <- data.frame(a0 = c(rep(a0, max_age), rep(NA,99*max_age)),  # create data frame with 100x1:50 year old elephants and priors 
                      a1 = c(rep(a1, max_age), rep(NA,99*max_age)),
                      c  = c(rep(c,  max_age), rep(NA,99*max_age)),
                      b0 = c(rep(b0, max_age), rep(NA,99*max_age)),
                      b1 = c(rep(b1, max_age), rep(NA,99*max_age)),
                      age = rep(1:max_age,100),
                      y = rep(NA,100*max_age),
                      iter = rep(1:100, each = max_age))
for(i in 1:99){ # populate rest of data frame with other prior values around those set out above
  to.plot$a0[(i*max_age+1):(i*max_age+max_age)] <- rep(rnorm(1, a0, 4),   max_age)
  to.plot$a1[(i*max_age+1):(i*max_age+max_age)] <- rep(rnorm(1, a1, 0.5), max_age)
  to.plot$c[ (i*max_age+1):(i*max_age+max_age)] <- rep(rnorm(1, c,  1), max_age)
  to.plot$b0[(i*max_age+1):(i*max_age+max_age)] <- rep(rnorm(1, b0, 2),   max_age)
  to.plot$b1[(i*max_age+1):(i*max_age+max_age)] <- rep(rnorm(1, b1, 0.1),max_age)
}

to.plot$y  <- gompertz_bt(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1,
                          age = to.plot$age)  # calculate curves indicating the probability of survival for 100 draws from the prior distribution
to.plot$y_std <- to.plot$y / max(to.plot$y)   # standardised output
plot(NULL, xlab = 'age', ylab = 'mortality probability', xlim = c(0,max_age), ylim = c(0,1), las = 1) # create plot window
for(i in 1:100){                              # plot distribution lines
  data <- to.plot[to.plot$iter == i,]
  lines(data$y_std ~ data$age, col = rgb(0,0,1,0.2))
}

pM.Gb <-  c(a0,a1,c,b0,b1)        # prior means
pSD.Gb <- c(1,0.1,0.4,2,0.05)     # prior standard deviations

#### set iteration arguments ####
start_year <- 1972    # year of first observations
end_year <- 2021      # year of last observations
iter <- 50000         # number of iterations
sim <- 4              # number of chains
warmup <- iter/2      # number of iterations to use during warmup period
thin <- 100           # thinning number - reduce data output

#### run model ####
Go.Bt <- basta(basta_anp, studyStart = start_year, # run basta function to estimate best GoBt distribution for ANP males
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
plot(Go.Bt, plot.trace = FALSE)  # very similar to above

#### simulate from output ####
a0.post <- rnorm(100, Go.Bt$coefficients[1,1], Go.Bt$coefficients[1,2])  # create a0 distribution from posterior
a1.post <- rnorm(100, Go.Bt$coefficients[2,1], Go.Bt$coefficients[2,2])  # create a1 distribution from posterior
c.post  <- rnorm(100, Go.Bt$coefficients[3,1], Go.Bt$coefficients[3,2])  # create c  distribution from posterior
b0.post <- rnorm(100, Go.Bt$coefficients[4,1], Go.Bt$coefficients[4,2])  # create b0 distribution from posterior
b1.post <- rnorm(100, Go.Bt$coefficients[5,1], Go.Bt$coefficients[5,2])  # create b1 distribution from posterior

gompertz_bt <- function(a0, a1, c, b0, b1, age){                         # repeat function definition
  gompertz <- exp(b0 + b1*age)
  bathtub <- exp(a0 - a1*age) + c + gompertz
  return(bathtub)
}

mortality <- gompertz_bt(a0 = Go.Bt$coefficients[1,1],                   # plot distribution
                         a1 = Go.Bt$coefficients[2,1],
                         c  = Go.Bt$coefficients[3,1],
                         b0 = Go.Bt$coefficients[4,1],
                         b1 = Go.Bt$coefficients[5,1], 1:50) ; plot(mortality)

to.plot <- data.frame(a0 = rep(a0.post, each = 50),                      # create data frame for plotting posterior
                      a1 = rep(a1.post, each = 50),
                      c  = rep(c.post,  each = 50),
                      b0 = rep(b0.post, each = 50),
                      b1 = rep(b1.post, each = 50),
                      age = rep(1:50,100),
                      y = rep(NA,5000),
                      iter = rep(1:100, each = 50))
to.plot$y  <- gompertz_bt(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1,
                          age = to.plot$age)                             # calculate output for each value of posterior
plot(NULL, xlab = 'age', ylab = 'mortality probability', xlim = c(0,50), ylim = c(0,1), las = 1)  # create empty plot window
for(i in 1:100){                                                         # plot output lines
  data <- to.plot[to.plot$iter == i,]
  lines(data$y ~ data$age, col = rgb(0,0,1,0.2))
}

# draw new output curve
mortality <- gompertz_bt(a0 = Go.Bt$coefficients[1,1],
                         a1 = Go.Bt$coefficients[2,1],
                         c  = Go.Bt$coefficients[3,1],
                         b0 = Go.Bt$coefficients[4,1],
                         b1 = Go.Bt$coefficients[5,1], 1:70)             # mean output mortality curve
probs <- 1 - (mortality/max(mortality))                                  # convert mortality curve to survival curve
plot(probs)                                                              # plot survival curve

# create a fictional population with ages selected from that distribution
N <- 700                      # 700 elephants in population
elephants_ls <- list(         # make a list of all elephants and their true ages
  N = N,
  age = sample(1:70, N, replace = TRUE, prob = probs)
)
hist(elephants_ls$age)        # histogram of elephant age

basta_data$age_death <- ifelse(basta_data$death_year == 0, NA, basta_data$death_year - basta_data$byr) # age at death for elephants with known death year
hist(basta_data$age_death)    # histogram of years died

#### simulate process of assigning age categories to elephants -- THIS WORKS PERFECTLY, just missing double 0-5 threshold ####
# draw new output curves
max_age <- 60  # no male elephants over 60 years old
mortality <- gompertz_bt(a0 = Go.Bt$coefficients[1,1],
                         a1 = Go.Bt$coefficients[2,1],
                         c  = Go.Bt$coefficients[3,1],
                         b0 = Go.Bt$coefficients[4,1],
                         b1 = Go.Bt$coefficients[5,1], 1:max_age) # run mortality function on fictional population
plot(mortality)                                                   # plot mortality output
probs <- 1 - (mortality/max(mortality))                           # convert to survival
plot(probs)                                                       # plot survival output

# create a fictional population with ages selected from that distribution
N <- 200                   # number of elephants
K <- 7                     # number of age categories for elephants to be sorted into
elephants_ls <- list(      # list of elephants
  N = N,
  K = K,
  age = sample(1:max_age, N, replace = TRUE, prob = probs)
)
hist(elephants_ls$age)     # histogram of ages

# simulate observing ages with error and then binning them into age groups
error <- 3                 # Error (in either direction) of age estimates
elephants_ls$age_guess <- elephants_ls$age + sample(-error:error, N, replace = TRUE)  # add error for age observations
elephants_ls$age_category_index <- sapply(elephants_ls$age_guess,                     # categorise based on estimated ages
                                          function(x) which.max(x < c(10, 15, 20, 25, 40, 60, Inf)))
hist(elephants_ls$age_category_index)                                                 # histogram categories

# look at actual age vs biologist assigned age: chance of being mis-classified is lower if actual age is in middle of age category, and biased towards lower end of class.
data.frame(elephants_ls) %>%
  ggplot(aes(x = age, y = age_guess, col = factor(age_category_index))) + # compare estimated age against true age
    geom_point(size = 4, alpha = 0.6) +
    geom_vline(xintercept = c(15, 20, 25, 40, 60), col = factor(3:7), linetype = "dashed", alpha = 0.6) +
    theme_minimal() + 
    xlab("Assigned age") + ylab("Age")

# read in Stan model to estimate ages based on Gompertz bathtub distribution from ANP
latent_age_ordinal_model <- cmdstan_model("models/motnp_age_ordinal_regression.stan")  # Gompertz bathtub model

# Fit model with cmdstanr
age_estimation_fit <- latent_age_ordinal_model$sample(  # fit model
  data = elephants_ls,                                  # data list
  chains = 4,                                           # number of MCMC chains
  parallel_chains = 4,                                  # number of chains to run at the same time
  iter_sampling = 2000)                                 # length of each chain

# Examine the estimates 
age_est_mat <- age_estimation_fit$summary()[202:401, ]  # rows containing estimates for each elephant
summary(age_est_mat)                                    # summarise
hist(age_est_mat$mean)                                  # plot
hist(age_est_mat$rhat, breaks = 20)                     # check convergence

plot_data <- data.frame(age = elephants_ls$age,         # Biologists original age est
                        model_age = age_est_mat$mean)   # Mean modelled age

plot_data %>%
  ggplot(aes(x = factor(age), y = model_age)) +         # plot true age against modelled age estimate
    geom_point(size = 4,col = 'blue', alpha = 0.6) +
    geom_vline(xintercept = c(10, 15, 20, 25, 40, 60), linetype = "dashed", alpha = 0.6) +
    geom_hline(yintercept = c(10, 15, 20, 25, 40, 60), linetype = "dashed", alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0)+
    scale_y_continuous(limits = c(0,60))+
    theme_minimal() + 
    xlab("Assigned age") + ylab("Modelled age")

#### simulate process of assigning age categories to elephants -- THIS WORKS PERFECTLY ####
# draw new output curves
max_age <- 60  # no male elephants over 60 years old
mortality <- gompertz_bt(a0 = Go.Bt$coefficients[1,1],
                         a1 = Go.Bt$coefficients[2,1],
                         c  = Go.Bt$coefficients[3,1],
                         b0 = Go.Bt$coefficients[4,1],
                         b1 = Go.Bt$coefficients[5,1], 1:max_age) # run mortality function on fictional population
plot(mortality)                                                   # plot mortality output
probs <- 1 - (mortality/max(mortality))                           # convert to survival
plot(probs)                                                       # plot survival output

# create a fictional population with ages selected from that distribution
N <- 200                   # number of elephants
K <- 8                     # number of age categories for elephants to be sorted into PLUS ONE -- UPPER LIMIT INCLUDED
elephants_ls <- list(      # list of elephants
  N = N,
  K = K,
  age = sample(1:max_age, N, replace = TRUE, prob = probs)
)
hist(elephants_ls$age)     # histogram of ages

# simulate observing ages with error and then binning them into age groups
error <- 3                 # Error (in either direction) of age estimates
elephants_ls$age_guess <- elephants_ls$age + sample(-error:error, N, replace = TRUE)  # add error for age observations
elephants_ls$age_category_index <- sapply(elephants_ls$age_guess,                     # categorise based on estimated ages
                                          function(x) which.max(x < c(5, 10, 15, 20, 25, 40, 60, Inf)))
hist(elephants_ls$age_category_index)

# look at actual age vs biologist assigned age: chance of being mis-classified is lower if actual age is in middle of age category, and biased towards lower end of class.
data.frame(elephants_ls) %>%
  ggplot(aes(x = age, y = age_guess, col = factor(age_category_index))) + # compare estimated age against true age
  geom_point(size = 4, alpha = 0.6) +
  geom_vline(xintercept = c(5, 10, 15, 20, 25, 40, 60), col = factor(1:7), linetype = "dashed", alpha = 0.6) + # MORE INTERCEPTS
  theme_minimal() + 
  xlab("Assigned age") + ylab("Age")

# read in Stan model to estimate ages based on Gompertz bathtub distribution from ANP
latent_age_ordinal_model <- cmdstan_model("models/motnp_age_ordinal_regression.stan")  # Gompertz bathtub model

#Fit model with cmdstanr
age_estimation_fit <- latent_age_ordinal_model$sample(  # fit model
  data = elephants_ls,                                  # data list
  chains = 4,                                           # number of MCMC chains
  parallel_chains = 4,                                  # number of chains to run at the same time
  iter_sampling = 2000)                                 # length of each chain

#Examine the estimates. We can plot the estimated ages against the biologist assigned ages and ...
age_est_mat <- age_estimation_fit$summary()[202:401, ]  # rows containing estimates for each elephant
summary(age_est_mat)                                    # summarise
hist(age_est_mat$mean)                                  # plot
hist(age_est_mat$rhat, breaks = 20)                     # check convergence

plot_data <- data.frame(age = elephants_ls$age,         # Biologists original age est
                        model_age = age_est_mat$mean)   # Mean modelled age

plot_data %>%
  ggplot(aes(x = factor(age), y = model_age)) +         # plot true age against modelled age estimate
  geom_point(size = 4,col = 'blue', alpha = 0.6) +
  geom_vline(xintercept = c(5, 10, 15, 20, 25, 40, 60), linetype = "dashed", alpha = 0.6) +
  geom_hline(yintercept = c(5, 10, 15, 20, 25, 40, 60), linetype = "dashed", alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0)+
  scale_y_continuous(limits = c(0,60))+
  theme_minimal() + 
  xlab("Assigned age") + ylab("Modelled age")

#### run on MOTNP elephants -- THIS WORKS PERFECTLY ####
motnp_males <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_elenodes.csv') %>% 
  #filter(dem_class == 'AM' | dem_class == 'PM') %>% 
  filter(sex == 'M')               # read in MOTNP individual data and select only males
unique(motnp_males$age_category)   # check age categories
motnp_males$age_cat_id <- ifelse(motnp_males$age_category == "0-3", 1,  # standardise age categories -- all <6 = category 1
                                 ifelse(motnp_males$age_category == "3-4", 1,
                                        ifelse(motnp_males$age_category == "1-2", 1,
                                               ifelse(motnp_males$age_category == "7-8", 1,
                                                      ifelse(motnp_males$age_category == "4-5", 1, 
                                                             ifelse(motnp_males$age_category ==  "6-7", 1,
                                                                    ifelse(motnp_males$age_category == "8-9", 1, 
                                                                           ifelse(motnp_males$age_category == "5-6", 1, NA))))))))
motnp_males$age_cat_id <- ifelse(motnp_males$age_category == '9-10', 2, # standardise age categories
                                 ifelse(motnp_males$age_category == '10-15', 3,
                                        ifelse(motnp_males$age_category == '15-19', 4,
                                               ifelse(motnp_males$age_category == '20-25', 5,
                                                      ifelse(motnp_males$age_category == '25-40', 6,
                                                             ifelse(motnp_males$age_category == '40+', 7,
                                                                    motnp_males$age_cat_id))))))

N_motnp <- nrow(motnp_males)        # number of males in data set
K <- 8                              # number of age categories + 1
motnp_ls <- list(                   # list of elephants
  N = N_motnp,
  K = K,
  age_category_index = motnp_males$age_cat_id)
hist(motnp_ls$age_category_index)   # plot ages
#hist(elephants_ls$age_category_index)

# read in Stan model to estimate ages based on Gompertz bathtub distribution from ANP
latent_age_ordinal_model <- cmdstan_model("models/motnp_age_ordinal_regression.stan")  # read in model

#Fit model with cmdstanr
age_motnp_fit <- latent_age_ordinal_model$sample(   # fit model
  data = motnp_ls,                                  # data list
  chains = 4,                                       # number of MCMC chains
  parallel_chains = 4,                              # run all chains at same time
  iter_sampling = 2000)                             # length of chains

# Examine the estimates
age_est_mat <- age_motnp_fit$summary()[(N_motnp+2):(N_motnp*2+1), ] # select individuals
summary(age_est_mat)                                                # summarise
hist(age_est_mat$mean)                                              # plot mean estimates
hist(age_est_mat$rhat, breaks = 20)                                 # plot convergence estimates

plot_data <- data.frame(age = ifelse(motnp_ls$age == 1, 3,          # set "true age" to value in the middle of each category for comparison -- NOT TO BE USED FOR AY ACTUAL ANALYSIS, JUST CHECKING MODEL FIT
                                     ifelse(motnp_ls$age == 2, 8,
                                            ifelse(motnp_ls$age == 3, 12,
                                                   ifelse(motnp_ls$age == 4, 18,
                                                          ifelse(motnp_ls$age == 5, 22, 
                                                                 ifelse(motnp_ls$age == 6, 32, 45)))))),
                        model_age = age_est_mat$mean)               # mean modelled age

plot_data %>%
  ggplot(aes(x = factor(age), y = model_age)) +     # plot "true age" assigned above against modelled age
  geom_point(size = 4, col = 'blue', alpha = 0.6) +
  geom_vline(xintercept = c(5, 10, 15, 20, 25, 40, 60), linetype = "dashed", alpha = 0.6) +
  geom_hline(yintercept = c(5, 10, 15, 20, 25, 40, 60), linetype = "dashed", alpha = 0.6) +
  #geom_abline(slope = 1, intercept = 0)+
  scale_y_continuous(limits = c(0,60))+
  theme_minimal() + 
  xlab("Assigned age") + ylab("Modelled age")

# Now an equivalent posterior predictive plot from draws from the distribution to show uncertainty around the mean age (then need to do it with full uncertainty).
true_ages <- age_motnp_fit$draws("true_age", format = "df")  # extract actual draw values from model output
#mcmc_dens(true_ages)
true_ages <- true_ages[1:100,1:N_motnp]                      # first 100 draws for all elephants

df <- as.data.frame(do.call(rbind, true_ages)) %>%           # convert data frame format
  mutate(age_cat = motnp_ls$age) %>% relocate(age_cat) %>%
  mutate(ID = motnp_males$id) %>% relocate(ID)
df <- df %>% pivot_longer(cols = 3:102) %>% select(-name)    # make long format

df$true_age <- ifelse(df$age_cat == 1, 3,                    # set "true age" to value in the middle of each category for comparison -- NOT TO BE USED FOR AY ACTUAL ANALYSIS, JUST CHECKING MODEL FIT
                      ifelse(df$age_cat == 2, 8,
                             ifelse(df$age_cat == 3, 12,
                                    ifelse(df$age_cat == 4, 18,
                                           ifelse(df$age_cat == 5, 22, 
                                                  ifelse(df$age_cat == 6, 32, 45))))))

df %>% ggplot(aes(x = true_age, y = value, group = factor(ID))) +   # plot "true age" assigned above against modelled age
  geom_point(size = 2,col = 'blue', alpha = 0.1) +
  #stat_halfeye() +
  geom_vline(xintercept = c(5,10,15,20,25,40,60), linetype = "dashed", alpha = 0.6) +
  geom_hline(yintercept = c(5,10,15,20,25,40,60), linetype = "dashed", alpha = 0.6) +
  theme_bw() + 
  xlab("Assigned age") + ylab("Modelled age")

#### run on ANP elephants to test if it can accurately predict the ages of ANP ####
motnp_males <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_elenodes.csv') %>% 
  #filter(dem_class == 'AM' | dem_class == 'PM') %>% 
  filter(sex == 'M')               # read in MOTNP individual data and select only males
unique(motnp_males$age_category)   # check age categories
motnp_males$age_cat_id <- ifelse(motnp_males$age_category == "0-3", 1,  # standardise age categories -- all <6 = category 1
                                 ifelse(motnp_males$age_category == "3-4", 1,
                                        ifelse(motnp_males$age_category == "1-2", 1,
                                               ifelse(motnp_males$age_category == "7-8", 1,
                                                      ifelse(motnp_males$age_category == "4-5", 1, 
                                                             ifelse(motnp_males$age_category ==  "6-7", 1,
                                                                    ifelse(motnp_males$age_category == "8-9", 1, 
                                                                           ifelse(motnp_males$age_category == "5-6", 1, NA))))))))
motnp_males$age_cat_id <- ifelse(motnp_males$age_category == '9-10', 2, # standardise age categories
                                 ifelse(motnp_males$age_category == '10-15', 3,
                                        ifelse(motnp_males$age_category == '15-19', 4,
                                               ifelse(motnp_males$age_category == '20-25', 5,
                                                      ifelse(motnp_males$age_category == '25-40', 6,
                                                             ifelse(motnp_males$age_category == '40+', 7,
                                                                    motnp_males$age_cat_id))))))

N_motnp <- nrow(motnp_males)        # number of males in data set
K <- 8                              # number of age categories + 1
motnp_ls <- list(                   # list of elephants
  N = N_motnp,
  K = K,
  age_category_index = motnp_males$age_cat_id)
hist(motnp_ls$age_category_index)   # plot ages
#hist(elephants_ls$age_category_index)

# read in Stan model to estimate ages based on Gompertz bathtub distribution from ANP
latent_age_ordinal_model <- cmdstan_model("models/motnp_age_ordinal_regression.stan")  # read in model

#Fit model with cmdstanr
age_motnp_fit <- latent_age_ordinal_model$sample(   # fit model
  data = motnp_ls,                                  # data list
  chains = 4,                                       # number of MCMC chains
  parallel_chains = 4,                              # run all chains at same time
  iter_sampling = 2000)                             # length of chains

# Examine the estimates
age_est_mat <- age_motnp_fit$summary()[(N_motnp+2):(N_motnp*2+1), ] # select individuals
summary(age_est_mat)                                                # summarise
hist(age_est_mat$mean)                                              # plot mean estimates
hist(age_est_mat$rhat, breaks = 20)                                 # plot convergence estimates

plot_data <- data.frame(age = ifelse(motnp_ls$age == 1, 3,          # set "true age" to value in the middle of each category for comparison -- NOT TO BE USED FOR AY ACTUAL ANALYSIS, JUST CHECKING MODEL FIT
                                     ifelse(motnp_ls$age == 2, 8,
                                            ifelse(motnp_ls$age == 3, 12,
                                                   ifelse(motnp_ls$age == 4, 18,
                                                          ifelse(motnp_ls$age == 5, 22, 
                                                                 ifelse(motnp_ls$age == 6, 32, 45)))))),
                        model_age = age_est_mat$mean)               # mean modelled age

plot_data %>%
  ggplot(aes(x = factor(age), y = model_age)) +     # plot "true age" assigned above against modelled age
  geom_point(size = 4, col = 'blue', alpha = 0.6) +
  geom_vline(xintercept = c(5, 10, 15, 20, 25, 40, 60), linetype = "dashed", alpha = 0.6) +
  geom_hline(yintercept = c(5, 10, 15, 20, 25, 40, 60), linetype = "dashed", alpha = 0.6) +
  #geom_abline(slope = 1, intercept = 0)+
  scale_y_continuous(limits = c(0,60))+
  theme_minimal() + 
  xlab("Assigned age") + ylab("Modelled age")

# Now an equivalent posterior predictive plot from draws from the distribution to show uncertainty around the mean age (then need to do it with full uncertainty).
true_ages <- age_motnp_fit$draws("true_age", format = "df")  # extract actual draw values from model output
#mcmc_dens(true_ages)
true_ages <- true_ages[1:100,1:N_motnp]                      # first 100 draws for all elephants

df <- as.data.frame(do.call(rbind, true_ages)) %>%           # convert data frame format
  mutate(age_cat = motnp_ls$age) %>% relocate(age_cat) %>%
  mutate(ID = motnp_males$id) %>% relocate(ID)
df <- df %>% pivot_longer(cols = 3:102) %>% select(-name)    # make long format

df$true_age <- ifelse(df$age_cat == 1, 3,                    # set "true age" to value in the middle of each category for comparison -- NOT TO BE USED FOR AY ACTUAL ANALYSIS, JUST CHECKING MODEL FIT
                      ifelse(df$age_cat == 2, 8,
                             ifelse(df$age_cat == 3, 12,
                                    ifelse(df$age_cat == 4, 18,
                                           ifelse(df$age_cat == 5, 22, 
                                                  ifelse(df$age_cat == 6, 32, 45))))))

df %>% ggplot(aes(x = true_age, y = value, group = factor(ID))) +   # plot "true age" assigned above against modelled age
  geom_point(size = 2,col = 'blue', alpha = 0.1) +
  #stat_halfeye() +
  geom_vline(xintercept = c(5,10,15,20,25,40,60), linetype = "dashed", alpha = 0.6) +
  geom_hline(yintercept = c(5,10,15,20,25,40,60), linetype = "dashed", alpha = 0.6) +
  theme_bw() + 
  xlab("Assigned age") + ylab("Modelled age")
