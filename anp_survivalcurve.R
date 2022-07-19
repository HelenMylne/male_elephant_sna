#### Information #####
# To estimate the ages of individuals in other populations, first come up with a survival curve for the Amboseli elephants which have good birth and death data.
# survival curve options = exponential/logistic/Gompertz/Weibull, where shape = simple/bathtub/Makeham
# follow through Mia's cde to determine which curve will best fit Amboseli data. she ran it on a high performance cluster, so it’s set up to run in parallel, but the basta model can also be run in parallel on its own (just like Stan). we will need to implement it as a custom function in Stan for the final analysis, given that we need to define it in the prior anyway. 

#### Load packages ####
lapply(c("BaSTA"), require, character.only = TRUE) # required packages for the running the bayesian hierarchichal function
library(tidyverse)

#### Read in Amboseli data ####
males <- readxl::read_excel('data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% 
  select(CASENAME, BYR, DYR) %>% 
  janitor::clean_names()
str(males)

sightings <- readxl::read_excel('data_raw/Raw_ATE_Sightings_Fishlock220224.xlsx') %>% janitor::clean_names()
sightings$year <- lubridate::year(sightings$obs_date)
sightings <- sightings[,c('casename','year')]
colnames(sightings)[1] <- 'id'
tapply(X = sightings$id, INDEX = sightings$year, FUN = table)

#### Create data frame for basta ####
# Create columns for sightings in each year
males$obs_1972 <- NA ; males$obs_1973 <- NA ; males$obs_1974 <- NA ; males$obs_1975 <- NA
males$obs_1976 <- NA ; males$obs_1977 <- NA ; males$obs_1978 <- NA ; males$obs_1979 <- NA
males$obs_1980 <- NA ; males$obs_1981 <- NA ; males$obs_1982 <- NA ; males$obs_1983 <- NA
males$obs_1984 <- NA ; males$obs_1985 <- NA ; males$obs_1986 <- NA ; males$obs_1987 <- NA
males$obs_1988 <- NA ; males$obs_1989 <- NA ; males$obs_1990 <- NA ; males$obs_1991 <- NA
males$obs_1992 <- NA ; males$obs_1993 <- NA ; males$obs_1994 <- NA ; males$obs_1995 <- NA
males$obs_1996 <- NA ; males$obs_1997 <- NA ; males$obs_1998 <- NA ; males$obs_1999 <- NA
males$obs_2000 <- NA ; males$obs_2001 <- NA ; males$obs_2002 <- NA ; males$obs_2003 <- NA
males$obs_2004 <- NA ; males$obs_2005 <- NA ; males$obs_2006 <- NA ; males$obs_2007 <- NA
males$obs_2008 <- NA ; males$obs_2009 <- NA ; males$obs_2010 <- NA ; males$obs_2011 <- NA
males$obs_2012 <- NA ; males$obs_2013 <- NA ; males$obs_2014 <- NA ; males$obs_2015 <- NA
males$obs_2016 <- NA ; males$obs_2017 <- NA ; males$obs_2018 <- NA ; males$obs_2019 <- NA
males$obs_2020 <- NA ; males$obs_2021 <- NA

(N <- 2022-1972)
for(i in 1:nrow(males)){
  for(j in 1:N){
    casename <- sightings[sightings$id == males$casename[i],]
    males[i,j+3] <- ifelse(length(which(casename$year == j+1971)) > 0, 1, 0)
  }
}

#### Prepare input data ####
# Getting start and end times from the input datasets - will depend on your input if this works, but can also manually assign the year the observations began and ended
Ss <- 1972  # start of study
Se <- 2021  # end of study

# prep data
basta_data <- males
basta_data$casename <- paste0('M', basta_data$casename)
basta_data <- basta_data[-c(1:4,6,8,11,14:18,22,23,26:30,35,38,44:46,52,55,56,58,59,62,63,66,68,71,72,76,83,87,90,93,94,96,105,107,108,111,114,115,118,127,128,132,133,135,136,138,145,146,151,155,159,162,174,179,182,186,187,190,191,193,197,201,208,210,213,220,224,230,235,238,244,248,250,254,256,259,261,263,266,268,288,291,292,295,298,306,307,309,319,327,339,342,346,370,384,416,417,435,449,472,717,693,708,735,769,774),]  # elephants that have sightings before birth/after death
basta_data <- basta_data[basta_data$casename %in% sort(unique(paste0('M',sightings$id))),]

basta_data$dyr <- ifelse(basta_data$dyr < 100,  0, basta_data$dyr)
basta_data$byr <- ifelse(basta_data$byr < 1972, 0, basta_data$byr)

################# Single chain, 500 iterations only #################
# Set the arguments for the iterations
iter <- 500        # number of iterations
sim <- 1           # number of chains
warmup <- iter/2   # number of iterations to use during warmup/burn-in period
thin <- 10         # thinning number

#### Gompertz simple ####
gompertz <- function(b0, b1, age){
  y <- exp(b0 + b1*age)
  return(y)
}
mortality <- gompertz(b0 = 0.3 , b1 = 0.6, 1:70) ; plot(mortality)
mortality <- gompertz(b0 = 0.1 , b1 = 0.05, 1:70) ; plot(mortality)

pM <-  c(0.1,0.05)  # prior means for parameters in a Gompertz model with a simple shape
pSD <- c(0.1,0.05)  # prior standard deviation for parameters in a Gompertz model with a simple shape

Go.Si <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               model = "GO", shape = c("simple"))
summary(Go.Si)
plot(Go.Si)                       # Plot traces
plot(Go.Si, plot.trace = FALSE)   # Plot survival and mortality curves

par(mfrow = c(1,1))

#### Gompertz bathtub ####
gompertz_bt <- function(a0, a1, c, b0, b1, age){
  gompertz <- gompertz(b0, b1, age)
  bathtub <- exp(a0 - a1*age) + c + gompertz
  return(bathtub)
}

mortality <- gompertz_bt(a0 = -3, a1 = 0.2, c = 0.001, b0 = -4 , b1 = 0.05, 1:100) ; plot(mortality) # original prior taken from Killer Whale paper
mortality <- gompertz_bt(a0 = 0, a1 = 0.2, c = 0.001, b0 = -4 , b1 = 0.05, 1:100)  ; plot(mortality) # adapted prior allowing greater infant mortality in first 10 years

pM <-  c(0,0.2,0.001,-4,0.05)  # prior means for parameters in a Gompertz model with a bathtub shape
pSD <- c(1,0.05,0.001,1,0.01)  # prior standard deviation for parameters in a Gompertz model with a bathtub shape
Go.Bt <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               model = "GO", shape = c("bathtub"))
summary(Go.Bt)
plot(Go.Bt)                       # Plot traces
plot(Go.Bt, plot.trace = FALSE)   # Plot survival and mortality curves

par(mfrow = c(1,1))

#### Weibull simple ####
weibull <- function(b0, b1, age){
  y <- (b0) * (b1^b0) * (age^(b0-1))
  return(y)
}
mortality <- weibull(b0 = 4, b1 = 2, age = 1:70)      ; plot(mortality)

pM <-  c(4,2)    # prior means for parameters in a Weibull model with a simple shape
pSD <- c(1.5,1)  # prior standard deviation for parameters in a Weibull model with a simple shape -- smaller values makes no visible difference to final curves
Wb.Si <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = 1, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               model = "WE", shape = c("simple"))
summary(Wb.Si)
plot(Wb.Si)                       # Plot traces
plot(Wb.Si, plot.trace = FALSE)   # Plot survival and mortality curves

par(mfrow = c(1,1))

#### Weibull bathtub ####
weibull_bt <- function(a0, a1, c, b0, b1, age) {
  f1 <- b0 * b1^(b0) * age^(b0-1)
  f2 <- exp(a0 - a1*age)
  y <- f1 + f2 + c
  return(y)
}

mortality <- weibull_bt(a0 = 1, a1 = 0.5, c = 1, b0 = 1.5, b1 = 0.1, 1:70)  ; plot(mortality)
mortality <- weibull_bt(a0 = 1, a1 = 0.3, c = 1, b0 = 1.5, b1 = 0.1, 1:70)  ; plot(mortality)
mortality <- weibull_bt(a0 = 2, a1 = 0.2, c = 1, b0 = 1, b1 = 1.2, 1:70)  ; plot(mortality)
mortality <- weibull_bt(a0 = 0.5, a1 = 0.2, c = 1, b0 = 1.2, b1 = 0.5, 1:70)  ; plot(mortality, las = 1)

to.plot <- data.frame(a0 = c(rep(0.5, 70), rep(NA,6930)),
                      a1 = c(rep(0.2, 70), rep(NA,6930)),
                      c  = c(rep(1,   70), rep(NA,6930)),
                      b0 = c(rep(1.2, 70), rep(NA,6930)),
                      b1 = c(rep(0.5, 70), rep(NA,6930)),
                      age = rep(1:70, 100),
                      y = rep(NA,7000), iter = rep(1:100,each = 70))
to.plot$y <- weibull_bt(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1, age = to.plot$age)
for(i in 1:99){
  to.plot$a0[(i*70+1):(i*70+70)] <- rep(rnorm(1, 0.5, 0.3),70)
  to.plot$a1[(i*70+1):(i*70+70)] <- rep(rnorm(1, 0.2, 0.1),70)
  to.plot$c[(i*70+1):(i*70+70)]  <- rep(rnorm(1, 1, 1),70)
  to.plot$b0[(i*70+1):(i*70+70)] <- rep(rnorm(1, 1.2, 0.5),70)
  to.plot$b1[(i*70+1):(i*70+70)] <- rep(rnorm(1, 0.5, 0.3),70)
}
to.plot$y  <- weibull_bt(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1, 1:70)
plot(NULL, xlab = 'age', ylab = 'mortality', xlim = c(0,70), ylim = c(0,20), las = 1)
for(i in 1:100){
  data <- to.plot[to.plot$iter == i, 6:7]
  lines(data$y ~ data$age, col = rgb(0,0,1,0.2))
}

pM <-  c(0.5,0.2,1,1.2,0.5) # prior means for parameters in a Weibull model with a bathtub shape
pSD <- c(0.3,0.1,1,0.5,0.3)  # prior standard deviation for parameters in a Weibull model with a bathtub shape

Wb.Bt <- basta(basta_data, studyStart = Ss,
                studyEnd = Se, nsim = sim, niter = iter,
                burnin = warmup, thinning = thin,
                thetaPriorMean = pM,thetaPriorSd = pSD,
                model = "WE", shape = c("bathtub"))
summary(Wb.Bt)
plot(Wb.Bt)                       # Plot traces
plot(Wb.Bt, plot.trace = FALSE)   # Plot survival and mortality curves

################# Repeat with 4 chains and more iterations #################
# Set the arguments for the iterations
iter <- 50000      # number of iterations
sim <- 4           # number of chains
warmup <- iter/2   # number of iterations to use during warmup/burn-in period
thin <- 100        # thinning number

#### Gompertz simple -- DIC = uncalculated, lack of convergence ####
pM <-  c(0.1,0.05)  # prior means for parameters in a Gompertz model with a simple shape
pSD <- c(0.1,0.05)  # prior standard deviation for parameters in a Gompertz model with a simple shape
Go.Si <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               model = "GO", shape = c("simple"))
summary(Go.Si)
plot(Go.Si)                       # Plot traces
plot(Go.Si, plot.trace = FALSE)   # Plot survival and mortality curves

#### Gompertz bathtub -- DIC = 10256.84 ####
pM <-  c(0,0.2,0.001,-4,0.05)  # prior means for parameters in a Gompertz model with a bathtub shape
pSD <- c(1,0.05,0.001,1,0.01)  # prior standard deviation for parameters in a Gompertz model with a bathtub shape

Go.Bt <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               model = "GO", shape = c("bathtub"))
summary(Go.Bt)
plot(Go.Bt)                       # Plot traces
plot(Go.Bt, plot.trace = FALSE)   # Plot survival and mortality curves

#### Weibull simple -- DIC = 10167.53 ####
pM <-  c(4,2)    # prior means for parameters in a Weibull model with a simple shape
pSD <- c(1.5,1)  # prior standard deviation for parameters in a Weibull model with a simple shape

Wb.Si <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               model = "WE", shape = c("simple"))
summary(Wb.Si)
plot(Wb.Si)                       # Plot traces
plot(Wb.Si, plot.trace = FALSE)   # Plot survival and mortality curves

#### Weibull bathtub -- DIC = uncalculated, lack of convergence ####
pM <-  c(0.5,0.2,1,1.2,0.5) # prior means for parameters in a Weibull model with a bathtub shape
pSD <- c(0.3,0.1,1,0.5,0.3)  # prior standard deviation for parameters in a Weibull model with a bathtub shape

Wb.Bt <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               model = "WE", shape = c("bathtub"))
summary(Wb.Bt)
plot(Wb.Bt)                       # Plot traces
plot(Wb.Bt, plot.trace = FALSE)   # Plot survival and mortality curves

################# Improve priors using AEDD ####
#### Load packages ####
library(cmdstanr)
library(posterior)
library(bayesplot)
library(survival)
library(brms)
library(tidybayes)
library(ggthemes)
library(dplyr)
library(data.table)

#### Read in AEDD -- old ####
#aedd <- read_csv('data_raw/AfricanElephantDemographicDatabase_aedd_22.04.04.csv') %>% 
  filter(sex == 'male' | sex == 'dispersed_male' | sex == 'family_male')
table(aedd$metric)
aedd <- aedd[!is.na(aedd$age_reported), c(1,3,5,6,7,8,13,14,15,22)]
boxplot(as.numeric(aedd$age_reported) ~ aedd$sex)
par(mai = c(1,5,0.2,0.2)); boxplot(as.numeric(aedd$age_reported) ~ aedd$sitename,
                                   horizontal = T, las = 1, ylab = '', xlab = 'age') ; par(mai = c(1,1,1,1))
str(aedd)

# create age data frame
table(aedd$age_reported)          # every study uses a different format!
length(unique(aedd$age_reported)) # 877 different values

aedd <- separate(aedd, age_reported, sep = '>', into = c('min_age','value'), remove = F) # if age is reported as >x, then x is minimum age it might have lived to -- censor at x
aedd$min_age <- ifelse(aedd$min_age == '', aedd$value, aedd$min_age)
aedd <- aedd %>% select(-value)

aedd <- separate(aedd, min_age, sep = '<', into = c('max_age','value'), remove = F)      # if age is reported as <x, then x is maximum age and can't assume it actually reached that age -- subtract 1
aedd$max_age <- ifelse(aedd$max_age == '', as.numeric(aedd$value) - 1, aedd$max_age)
aedd <- aedd %>% select(-value)

aedd <- separate(aedd, max_age, sep = '-', into = c('lwr_age','upr_age'), remove = F)    # split age ranges up to give a maximum and minimum possible.
aedd$upr_age <- ifelse(is.na(aedd$upr_age) == TRUE, aedd$lwr_age, aedd$upr_age)
length(unique(aedd$upr_age))

str(aedd$upr_age)
aedd$age_all <- ifelse(aedd$upr_age == "calf",                 median(x = c(0,5)),
                ifelse(aedd$upr_age == "infant",               median(x = c(0,5)),
                ifelse(aedd$upr_age == "juvenile; infant",     median(x = c(0,10)),
                ifelse(aedd$upr_age == "offspring",            median(x = c(0,10)),
                ifelse(aedd$upr_age == "young",                median(x = c(0,25)),
                ifelse(aedd$upr_age == "juvenile",             median(x = c(5,10)),
                ifelse(aedd$upr_age == "halfgrown",            median(x = c(5,15)),
                ifelse(aedd$upr_age == "juvenile; subadult",   median(x = c(5,15)),
                ifelse(aedd$upr_age == "sexually_immature",    median(x = c(5,20)),
                ifelse(aedd$upr_age == "subadult",             median(x = c(10,20)),
                ifelse(aedd$upr_age == "intermediate",         median(x = c(10,25)),
                ifelse(aedd$upr_age == "subadult; adult",      median(x = c(15,25)),
                ifelse(aedd$upr_age == "sexually_mature",      median(x = c(15,30)),
                ifelse(aedd$upr_age == "fullgrown",            median(x = c(20,30)),
                ifelse(aedd$upr_age == "adult",                median(x = c(20,50)),
                ifelse(aedd$upr_age == "old, nonreproductive", median(x = c(50,70)),
                ifelse(aedd$upr_age == "unknown",              NA,    aedd$upr_age)))))))))))))))))
table(aedd$age_all)
aedd$lwr_age <- as.numeric(aedd$lwr_age)
aedd$upr_age <- as.numeric(aedd$upr_age)
aedd$age <- NA
for(i in 1:nrow(aedd)){
  aedd$age[i] <- median(x = c(aedd$lwr_age[i], aedd$upr_age[i]))
}

# add censorship (all are alive, censor = 1)
aedd$censor <- 1

# clean up data
age_cens <- aedd[!is.na(aedd$age),c('age','censor','sex')]

#### estimate Weibull distribution -- old, 22.04.04_Weibull.R ####
colnames(age_cens)[1] <- 'age'
length(which(age_cens$age == 0))
age_cens$age_non0 <- ifelse(age_cens$age == 0, 0.01, age_cens$age) # wouldn't run when age = 0

# Examine default brms priors (male elephants)
bfit_m_default <- age_cens_m %>%
  brm(age_non0 | cens(censor) ~ 1,
      data = ., family = "weibull", 
      chains = 4, cores = 4)
prior_summary(bfit_m_default)
# Defaults
#intercept - student_t(3, 3.1, 2.5)
#shape - gamma(0.01, 0.01)
plot(bfit_m_default)

# Default prior predictive check
prior_intercept <- rstudent_t(num_samples, 3, 3.1, 2.5)
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

# Prior predictive check
num_samples <- 500
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
bfit_m <- age_cens_m %>%
  brm(age_non0 | cens(censor) ~ 1,
      prior = c(
        prior(NEW_PRIOR(NEW_PRIOR_PARAMETERS), class = Intercept),
        prior(NEW_PRIOR(NEW_PRIOR_PARAMETERS), class = shape)
      ),
      data = ., family = "weibull", 
      chains = 4, cores = 4)
plot(bfit_m)

# extract scale and shape from posterior draws
bfit_m_draws <- as_draws_df(bfit_m) %>%
  mutate(scale = exp(b_Intercept) / (gamma(1 + 1 / shape)))

est_shape <- mean(bfit_m_draws$shape)  # 
est_scale <- mean(bfit_m_draws$scale)  # 

# posterior predictive plot based on the mean
probs <- dweibull(1:100, shape = est_shape, scale = est_scale)
probs <- probs / sum(probs)
plot(probs) 

# simulate ages for 3k elephants
hist(rweibull(3000,est_shape,est_scale)) # still horrendous


### start over -- need to use the rest of the data too: e.g. for population size, if estimate = 8 that's 8 individuals of that age, not 1, but for age at first calving there is no estimate for age_reported but they must have reached at least the calving age, and for survivorship they provide an estimate for the proportion reaching those ages ####
aedd <- read_csv('data_raw/AfricanElephantDemographicDatabase_aedd_22.04.04.csv')
aedd <- aedd[aedd$species == 'L. africana' & !is.na(aedd$species),]
sort(unique(aedd$metric))

## studies of interest
ages <- aedd[aedd$metric == 'first calf age min' |
               aedd$metric == 'first calf age mean' |
               aedd$metric == 'first calf age median' |
               aedd$metric == 'first calf age log_mean' |
               aedd$metric == 'population size' |
               aedd$metric == 'life span max',
             c(1,6,7,8,13,14,22)]
colnames(ages)
#rm(aedd)

## minimum age first calf: estimate = censored age, 31 elephants made it to at least 9 -- all individuals in study birthed, minimum age was 9, therefore all reached at least 9 years old
fca_min <- ages[ages$metric == 'first calf age min',]
fca_min
fca_min_cens <- data.frame(age = rep(NA, nrow(fca_min)),
                           sex = rep(NA, nrow(fca_min)),
                           sid = rep(NA, nrow(fca_min)),
                           loc = rep(NA, nrow(fca_min)))
for(i in 1:nrow(fca_min)){
  fca_min_cens$age[i] <- list(rep(fca_min$estimate[i], fca_min$sample_size[i]))
  fca_min_cens$sex[i] <- list(rep(fca_min$sex[i], fca_min$sample_size[i]))
  fca_min_cens$sid[i] <- list(rep(fca_min$study_id[i], fca_min$sample_size[i]))
  fca_min_cens$loc[i] <- list(rep(fca_min$sitename[i], fca_min$sample_size[i]))
}
fca_min_cens_unlist <- data.frame(age = unlist(fca_min_cens$age),
                                  alive = TRUE,
                                  cens = 1,
                                  sex = unlist(fca_min_cens$sex),
                                  study_id = unlist(fca_min_cens$sid),
                                  location = unlist(fca_min_cens$loc))
rm(fca_min, fca_min_cens)

## median age first calf: estimate = censored age, sample size/2 = number made it to at least that age -- median age for first calf, so at least half must have reached this age
fca_mid <- ages[ages$metric == 'first calf age median',]
fca_mid
fca_mid_cens <- data.frame(age = rep(NA, nrow(fca_mid)),
                           sex = rep(NA, nrow(fca_mid)),
                           sid = rep(NA, nrow(fca_mid)),
                           loc = rep(NA, nrow(fca_mid)))
for(i in 1:nrow(fca_mid)){
  fca_mid_cens$age[i] <- list(rep(fca_mid$estimate[i], fca_mid$sample_size[i]/2))
  fca_mid_cens$sex[i] <- list(rep(fca_mid$sex[i], fca_mid$sample_size[i]/2))
  fca_mid_cens$sid[i] <- list(rep(fca_mid$study_id[i], fca_mid$sample_size[i]/2))
  fca_mid_cens$loc[i] <- list(rep(fca_mid$sitename[i], fca_mid$sample_size[i]/2))
}
fca_mid_cens_unlist <- data.frame(age = unlist(fca_mid_cens$age),
                                  alive = TRUE,
                                  cens = 1,
                                  sex = unlist(fca_mid_cens$sex),
                                  study_id = unlist(fca_mid_cens$sid),
                                  location = unlist(fca_mid_cens$loc))

rm(fca_mid, fca_mid_cens)

## mean age first calf: estimate = censored age, sample size 1 if not indicated, not sure if it is. remove study which is also present in fca_mid to avoid including same elephants twice. approximately half will have birthed above average age, therefore reasonable to assume that at least half survived to this age. if no sample size given, only include as a single record
fca_mean <- ages[ages$metric == 'first calf age mean' & ages$study_id != 3891,]
fca_mean
fca_mean$sex <- 'female'
fca_mean$sample_size[is.na(fca_mean$sample_size)] <- 1
fca_mean_cens <- data.frame(age = rep(NA, nrow(fca_mean)),
                            sex = rep(NA, nrow(fca_mean)),
                            sid = rep(NA, nrow(fca_mean)),
                            loc = rep(NA, nrow(fca_mean)))
for(i in 1:nrow(fca_mean)){
  fca_mean_cens$age[i] <- list(rep(fca_mean$estimate[i],
                                   ifelse(fca_mean$sample_size[i] > 1,
                                          fca_mean$sample_size[i]/2,
                                          fca_mean$sample_size[i])))
  fca_mean_cens$sex[i] <- list(rep(fca_mean$sex[i],
                                   ifelse(fca_mean$sample_size[i] > 1,
                                          fca_mean$sample_size[i]/2,
                                          fca_mean$sample_size[i])))
  fca_mean_cens$sid[i] <- list(rep(fca_mean$study_id[i],
                                   ifelse(fca_mean$sample_size[i] > 1,
                                          fca_mean$sample_size[i]/2,
                                          fca_mean$sample_size[i])))
  fca_mean_cens$loc[i] <- list(rep(fca_mean$sitename[i],
                                   ifelse(fca_mean$sample_size[i] > 1,
                                          fca_mean$sample_size[i]/2,
                                          fca_mean$sample_size[i])))
}
fca_mean_cens_unlist <- data.frame(age = unlist(fca_mean_cens$age),
                                   alive = TRUE,
                                   cens = 1,
                                   sex = unlist(fca_mean_cens$sex),
                                   study_id = unlist(fca_mean_cens$sid),
                                   location = unlist(fca_mean_cens$loc))
rm(fca_mean, fca_mean_cens)

## log means age first calf: as above, but anti-log estimate first
fca_lmean <- ages[ages$metric == 'first calf age log_mean',]
fca_lmean$exp_estimate <- round(exp(fca_lmean$estimate),1)
fca_lmean
fca_lmean_cens <- data.frame(age = rep(NA, nrow(fca_lmean)),
                             sex = rep(NA, nrow(fca_lmean)),
                             sid = rep(NA, nrow(fca_lmean)),
                             loc = rep(NA, nrow(fca_lmean)))
for(i in 1:nrow(fca_lmean)){
  fca_lmean_cens$age[i] <- list(rep(fca_lmean$exp_estimate[i], fca_lmean$sample_size[i]/2))
  fca_lmean_cens$sex[i] <- list(rep(fca_lmean$sex[i], fca_lmean$sample_size[i]/2))
  fca_lmean_cens$sid[i] <- list(rep(fca_lmean$study_id[i], fca_lmean$sample_size[i]/2))
  fca_lmean_cens$loc[i] <- list(rep(fca_lmean$sitename[i], fca_lmean$sample_size[i]/2))
}
fca_lmean_cens_unlist <- data.frame(age = unlist(fca_lmean_cens$age),
                                    alive = TRUE,
                                    cens = 1,
                                    sex = unlist(fca_lmean_cens$sex),
                                    study_id = unlist(fca_lmean_cens$sid),
                                    location = unlist(fca_lmean_cens$loc))
rm(fca_lmean, fca_lmean_cens)

## maximum life span: estimate = uncensored age, but only for a single individual (199 animals in study, but only one reached this age)
lsm <- ages[ages$metric == 'life span max',]
lsm
lsm_cens_unlist <- data.frame(age = lsm$estimate,
                              alive = FALSE,
                              cens = 0,
                              sex = lsm$sex,
                              study_id = lsm$study_id,
                              location = lsm$sitename)
rm(lsm)

## population size: estimate = number observed to reach minimum age (censor at min or median?)
ps <- ages[ages$metric == 'population size' & !is.na(ages$age_reported),]
unique(ps$age_reported)
ps$age <- with(ps,
               case_when(#age_reported == 'calf'    ~ min(x = c(0,5)),
                 #age_reported == 'juvenile' ~ min(x = c(5,10)),
                 #age_reported == 'halfgrown' ~ min(x = c(10,20)),
                 #age_reported == 'sexually_immature' ~ min(x = c(10,20)),
                 #age_reported == 'sexually_mature' ~ min(x = c(20,30)),
                 #age_reported == 'fullgrown' ~ min(x = c(20,50)),
                 #age_reported == 'adult'   ~ min(x = c(20,50)),
                 age_reported == '<1'      ~ min(x = c(0.1,1)),
                 age_reported == '0-1'     ~ min(x = c(0.1,1)),
                 age_reported == '0-4.9'   ~ min(x = c(0.1,4.9)),
                 age_reported == '0-5'     ~ min(x = c(0.1,5)),
                 age_reported == '<10'     ~ min(x = c(0.1,10)),
                 age_reported == '1-11'    ~ min(x = c(1,11)),
                 age_reported == '1-12'    ~ min(x = c(1,12)),
                 age_reported == '4-10'    ~ min(x = c(4,10)),
                 age_reported == '5-9.9'   ~ min(x = c(5,9.9)),
                 age_reported == '6-15'    ~ min(x = c(6,15)),
                 age_reported == '10-14.9' ~ min(x = c(10,14.9)),
                 age_reported == '10-20'   ~ min(x = c(10,20)),
                 age_reported == '15-19.9' ~ min(x = c(15,19.9)),
                 age_reported == '16-30'   ~ min(x = c(16,30)),
                 age_reported == '20-24.9' ~ min(x = c(20,24.9)),
                 age_reported == '20-42'   ~ min(x = c(20,42)),
                 age_reported == '25-34.9' ~ min(x = c(25,34.9)),
                 age_reported == '35-39.9' ~ min(x = c(35,39.9)),
                 age_reported == '25-29.9' ~ min(x = c(10,20)),
                 age_reported == '30-34.9' ~ min(x = c(30,34.9)),
                 age_reported == '20-34.9' ~ min(x = c(20,34.9)),
                 age_reported == '35-49.9' ~ min(x = c(35,49.9)),
                 age_reported == '>30'     ~ min(x = c(30,65)),
                 age_reported == '5-10'    ~ min(x = c(5,10)),
                 age_reported == '10-15'   ~ min(x = c(10,15)),
                 age_reported == '15-20'   ~ min(x = c(15,20)),
                 age_reported == '20-25'   ~ min(x = c(20,25)),
                 age_reported == '25-30'   ~ min(x = c(25,30)),
                 age_reported == '30-35'   ~ min(x = c(30,35)),
                 age_reported == '35-40'   ~ min(x = c(35,40)),
                 age_reported == '40-45'   ~ min(x = c(40,45)),
                 age_reported == '>40'     ~ min(x = c(40,65)),
                 age_reported == '45-50'   ~ min(x = c(45,50)),
                 age_reported == '50-55'   ~ min(x = c(50,55)),
                 age_reported == '>50'     ~ min(x = c(50,65)),
                 age_reported == '55-60'   ~ min(x = c(55,60)),
                 age_reported == '60-65'   ~ min(x = c(60,65))))
ps <- ps[!is.na(ps$age),]
ps
ps_cens <- data.frame(age = rep(NA, nrow(ps)),
                      sex = rep(NA, nrow(ps)),
                      sid = rep(NA, nrow(ps)),
                      loc = rep(NA, nrow(ps)))
for(i in 1:nrow(ps)){
  ps_cens$age[i] <- list(rep(ps$age[i], ifelse(ps$estimate[i] < 1, 0, ps$estimate[i])))
  ps_cens$sex[i] <- list(rep(ps$sex[i], ifelse(ps$estimate[i] < 1, 0, ps$estimate[i])))
  ps_cens$sid[i] <- list(rep(ps$study_id[i], ifelse(ps$estimate[i] < 1, 0, ps$estimate[i])))
  ps_cens$loc[i] <- list(rep(ps$sitename[i], ifelse(ps$estimate[i] < 1, 0, ps$estimate[i])))
}
ps_cens_unlist <- data.frame(age = unlist(ps_cens$age),
                             alive = TRUE,
                             cens = 1,
                             sex = unlist(ps_cens$sex),
                             study_id = unlist(ps_cens$sid),
                             location = unlist(ps_cens$loc))
rm(ps, ps_cens)

age_cens <- rbind(fca_min_cens_unlist, fca_mid_cens_unlist, fca_mean_cens_unlist, fca_lmean_cens_unlist, lsm_cens_unlist, ps_cens_unlist)
rm(fca_min_cens_unlist, fca_mid_cens_unlist, fca_mean_cens_unlist, fca_lmean_cens_unlist, lsm_cens_unlist, ps_cens_unlist, i)

#####
## annual mortality rate: estimate = proportion of elephants in that age category found dead that year -- can only include those with a recorded sample size, else any kind of simulated individuals would contribute a weird amount to the final total
amr <-  ages[!is.na(ages$age_reported) & ages$metric == 'annual mortality rate',]
amr
amr$sim_sample <- round(10-10*amr$estimate,0)





s <-    ages[!is.na(ages$age_reported) & ages$metric == 'survivorship',]           # estimate = proportion of the population that survives up to a given age
sd <-   ages[!is.na(ages$age_reported) & ages$metric == 'survivorship: drought in first two years',]      # estimate = proportion of the population that survives up to a given age for elephants that experienced a drought year in the first two years of their life
snd <-  ages[!is.na(ages$age_reported) & ages$metric == 'survivorship: no drought in first two years',]   # estimate = proportion of the population that survives up to a given age for elephants that did not experience a drought year in the first two years of their life -- same study as that looking at drought in first two years, but no crossover with survivorship studies that ignore drought. All three can be included?
pgb <-  ages[!is.na(ages$age_reported) & ages$metric == 'probability of giving birth',]                   # sample size indicates the number of individuals in the study that reached each age category -- censor at minimum?
cp <-   ages[!is.na(ages$age_reported) & ages$metric == 'carcass proportion',]   # estimate = number of carcasses of a given age(-class) and/or sex divided by the number of all carcasses
##### old ####
# including reported ages
ages_reported <- aedd[!is.na(aedd$age_reported), c(1,3,5,6,7,8,13,14,15,22)]
str(ages_reported)
table(ages_reported$metric)
table(ages_reported$age_reported)          # every study uses a different format depending on focus
length(unique(ages_reported$age_reported)) # 877 different values

ages_reported <- separate(ages_reported, age_reported, sep = '>', into = c('min_age','value'), remove = F) # if age is reported as >x, then x is minimum age it might have lived to -- censor at x
ages_reported$min_age <- ifelse(ages_reported$min_age == '', ages_reported$value, ages_reported$min_age)
ages_reported <- ages_reported %>% select(-value)

ages_reported <- separate(ages_reported, min_age, sep = '<', into = c('max_age','value'), remove = F)      # if age is reported as <x, then x is maximum age and can't assume it actually reached that age -- subtract 1
ages_reported$max_age <- ifelse(ages_reported$max_age == '', as.numeric(ages_reported$value) - 1, ages_reported$max_age)
ages_reported <- ages_reported %>% select(-value)

ages_reported <- separate(ages_reported, max_age, sep = '-', into = c('lwr_age','upr_age'), remove = F)    # split age ranges up to give a maximum and minimum possible.
ages_reported$upr_age <- ifelse(is.na(ages_reported$upr_age) == TRUE, ages_reported$lwr_age, ages_reported$upr_age)
length(unique(ages_reported$upr_age))

str(ages_reported$upr_age)
ages_reported$age_all <- ifelse(ages_reported$upr_age == "calf",                 median(x = c(0,5)),
                                ifelse(ages_reported$upr_age == "infant",               median(x = c(0,5)),
                                       ifelse(ages_reported$upr_age == "juvenile; infant",     median(x = c(0,10)),
                                              ifelse(ages_reported$upr_age == "offspring",            median(x = c(0,10)),
                                                     ifelse(ages_reported$upr_age == "young",                median(x = c(0,25)),
                                                            ifelse(ages_reported$upr_age == "juvenile",             median(x = c(5,10)),
                                                                   ifelse(ages_reported$upr_age == "halfgrown",            median(x = c(5,15)),
                                                                          ifelse(ages_reported$upr_age == "juvenile; subadult",   median(x = c(5,15)),
                                                                                 ifelse(ages_reported$upr_age == "sexually_immature",    median(x = c(5,20)),
                                                                                        ifelse(ages_reported$upr_age == "subadult",             median(x = c(10,20)),
                                                                                               ifelse(ages_reported$upr_age == "intermediate",         median(x = c(10,25)),
                                                                                                      ifelse(ages_reported$upr_age == "subadult; adult",      median(x = c(15,25)),
                                                                                                             ifelse(ages_reported$upr_age == "sexually_mature",      median(x = c(15,30)),
                                                                                                                    ifelse(ages_reported$upr_age == "fullgrown",            median(x = c(20,30)),
                                                                                                                           ifelse(ages_reported$upr_age == "adult",                median(x = c(20,50)),
                                                                                                                                  ifelse(ages_reported$upr_age == "old, nonreproductive", median(x = c(50,70)),
                                                                                                                                         ifelse(ages_reported$upr_age == "unknown",              NA,
                                                                                                                                                ages_reported$upr_age)))))))))))))))))
table(ages_reported$age_all)
ages_reported$lwr_age <- as.numeric(ages_reported$lwr_age)
ages_reported$upr_age <- as.numeric(ages_reported$upr_age)
ages_reported$age <- NA
for(i in 1:nrow(ages_reported)){
  ages_reported$age[i] <- median(x = c(ages_reported$lwr_age[i], ages_reported$upr_age[i]))
}

ages_reported$censor <- 1    # censor = 1 when elephant is still alive -- database contains no information about whether individual is still alive or not, so assume all are alive

# ages unreported
ages_unrep <- aedd[is.na(aedd$age_reported), c(1,3,5,6,7,8,13,14,15,22)]
str(ages_unrep)
table(ages_unrep$metric)

ages_unrep$age <- NA
for(i in 1:nrow(ages_unrep)){
  ages_unrep$age[i]
}


ages_unrep <- separate(ages_unrep, age_reported, sep = '>', into = c('min_age','value'), remove = F) # if age is reported as >x, then x is minimum age it might have lived to -- censor at x
ages_unrep$min_age <- ifelse(ages_unrep$min_age == '', ages_unrep$value, ages_unrep$min_age)
ages_unrep <- ages_unrep %>% select(-value)

ages_unrep <- separate(ages_unrep, min_age, sep = '<', into = c('max_age','value'), remove = F)      # if age is reported as <x, then x is maximum age and can't assume it actually reached that age -- subtract 1
ages_unrep$max_age <- ifelse(ages_unrep$max_age == '', as.numeric(ages_unrep$value) - 1, ages_unrep$max_age)
ages_unrep <- ages_unrep %>% select(-value)

ages_unrep <- separate(ages_unrep, max_age, sep = '-', into = c('lwr_age','upr_age'), remove = F)    # split age ranges up to give a maximum and minimum possible.
ages_unrep$upr_age <- ifelse(is.na(ages_unrep$upr_age) == TRUE, ages_unrep$lwr_age, ages_unrep$upr_age)
length(unique(ages_unrep$upr_age))

str(ages_unrep$upr_age)
ages_unrep$age_all <- ifelse(ages_unrep$upr_age == "calf",                 median(x = c(0,5)),
                             ifelse(ages_unrep$upr_age == "infant",               median(x = c(0,5)),
                                    ifelse(ages_unrep$upr_age == "juvenile; infant",     median(x = c(0,10)),
                                           ifelse(ages_unrep$upr_age == "offspring",            median(x = c(0,10)),
                                                  ifelse(ages_unrep$upr_age == "young",                median(x = c(0,25)),
                                                         ifelse(ages_unrep$upr_age == "juvenile",             median(x = c(5,10)),
                                                                ifelse(ages_unrep$upr_age == "halfgrown",            median(x = c(5,15)),
                                                                       ifelse(ages_unrep$upr_age == "juvenile; subadult",   median(x = c(5,15)),
                                                                              ifelse(ages_unrep$upr_age == "sexually_immature",    median(x = c(5,20)),
                                                                                     ifelse(ages_unrep$upr_age == "subadult",             median(x = c(10,20)),
                                                                                            ifelse(ages_unrep$upr_age == "intermediate",         median(x = c(10,25)),
                                                                                                   ifelse(ages_unrep$upr_age == "subadult; adult",      median(x = c(15,25)),
                                                                                                          ifelse(ages_unrep$upr_age == "sexually_mature",      median(x = c(15,30)),
                                                                                                                 ifelse(ages_unrep$upr_age == "fullgrown",            median(x = c(20,30)),
                                                                                                                        ifelse(ages_unrep$upr_age == "adult",                median(x = c(20,50)),
                                                                                                                               ifelse(ages_unrep$upr_age == "old, nonreproductive", median(x = c(50,70)),
                                                                                                                                      ifelse(ages_unrep$upr_age == "unknown",              NA,
                                                                                                                                             ages_unrep$upr_age)))))))))))))))))
table(ages_unrep$age_all)
ages_unrep$lwr_age <- as.numeric(ages_unrep$lwr_age)
ages_unrep$upr_age <- as.numeric(ages_unrep$upr_age)
ages_unrep$age <- NA
for(i in 1:nrow(ages_unrep)){
  ages_unrep$age[i] <- median(x = c(ages_unrep$lwr_age[i], ages_unrep$upr_age[i]))
}

ages_unrep$censor <- 1    # censor = 1 when elephant is still alive -- database contains no information about whether individual is still alive or not, so assume all are alive


##### continue ####
# clean up
colnames(age_cens)[3] <- 'censor'
age_cens$sex <- ifelse(age_cens$sex == 'female', 'F', 'M')
age_cens$location_id <- as.integer(as.factor(age_cens$location))
sort(unique(age_cens$location)) # many different studies based on Amboseli: a) don't add in extra Amboseli data and b) should it be just one or is this biasing it towards the Amboseli shape?

### estimate Weibull distribution -- new, but still shit ####
length(which(age_cens$age == 0))

# Examine default brms priors (all elephants together)
bfit_default <- age_cens %>%
  brm(age | cens(censor) ~ 1,
      data = ., family = "weibull", 
      chains = 4, cores = 4)
prior_summary(bfit_default)
# Defaults
#intercept - student_t(3, 1.6, 2.5)
#shape - gamma(0.01, 0.01)
plot(bfit_default) # shape max >600

# Prior predictive check
num_samples <- 500

# Default priors
prior_intercept <- rstudent_t(num_samples, 3, 1.6, 2.5)
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
bfit <- age_cens %>% 
  brm(age | cens(censor) ~ 1,
      prior = c(
        prior(student_t(3,5,5), class = Intercept),
        prior(exponential(0.5), class = shape)
      ),
      data = ., family = "weibull", 
      chains = 4, cores = 4)
plot(bfit)  # much shorter shape tail, but a lot wider distribution

# extract scale and shape from posterior draws
# to convert to scale we need to both undo the link function by taking the exponent
# and then refer to the brms documentation to understand how the mean relates to the scale
bfit_draws <- as_draws_df(bfit) %>%
  mutate(scale = exp(b_Intercept) / (gamma(1 + 1 / shape)))

est_shape <- mean(bfit_draws$shape)  # 6.72
est_scale <- mean(bfit_draws$scale)  # 4279.97?!!

# posterior predictive plot based on the mean
probs <- dweibull(1:100, shape = est_shape, scale = est_scale)
probs <- probs / sum(probs)
plot(probs) 

# simulate ages for 3k elephants
hist(rweibull(3000,est_shape,est_scale)) # horrendous -- left skewed not right, peaking ay 4500 years old!)

#### Age estimation for mega data set based on Amboseli distribution ####
## population size: estimate = number observed to reach minimum age (censor at min or median?)
ps <- ages[ages$metric == 'population size' & !is.na(ages$age_reported),]
unique(ps$age_reported)
ps$min <- with(ps, case_when(age_reported == '<1'      ~ min(x = c(0.1,1)),
                             age_reported == '0-1'     ~ min(x = c(0.1,1)),
                             age_reported == '0-4.9'   ~ min(x = c(0.1,4.9)),
                             age_reported == '0-5'     ~ min(x = c(0.1,5)),
                             age_reported == '<10'     ~ min(x = c(0.1,10)),
                             age_reported == '1-11'    ~ min(x = c(1,11)),
                             age_reported == '1-12'    ~ min(x = c(1,12)),
                             age_reported == '4-10'    ~ min(x = c(4,10)),
                             age_reported == '5-9.9'   ~ min(x = c(5,9.9)),
                             age_reported == '6-15'    ~ min(x = c(6,15)),
                             age_reported == '10-14.9' ~ min(x = c(10,14.9)),
                             age_reported == '10-20'   ~ min(x = c(10,20)),
                             age_reported == '15-19.9' ~ min(x = c(15,19.9)),
                             age_reported == '16-30'   ~ min(x = c(16,30)),
                             age_reported == '20-24.9' ~ min(x = c(20,24.9)),
                             age_reported == '20-42'   ~ min(x = c(20,42)),
                             age_reported == '25-34.9' ~ min(x = c(25,34.9)),
                             age_reported == '35-39.9' ~ min(x = c(35,39.9)),
                             age_reported == '25-29.9' ~ min(x = c(10,20)),
                             age_reported == '30-34.9' ~ min(x = c(30,34.9)),
                             age_reported == '20-34.9' ~ min(x = c(20,34.9)),
                             age_reported == '35-49.9' ~ min(x = c(35,49.9)),
                             age_reported == '>30'     ~ min(x = c(30,65)),
                             age_reported == '5-10'    ~ min(x = c(5,10)),
                             age_reported == '10-15'   ~ min(x = c(10,15)),
                             age_reported == '15-20'   ~ min(x = c(15,20)),
                             age_reported == '20-25'   ~ min(x = c(20,25)),
                             age_reported == '25-30'   ~ min(x = c(25,30)),
                             age_reported == '30-35'   ~ min(x = c(30,35)),
                             age_reported == '35-40'   ~ min(x = c(35,40)),
                             age_reported == '40-45'   ~ min(x = c(40,45)),
                             age_reported == '>40'     ~ min(x = c(40,65)),
                             age_reported == '45-50'   ~ min(x = c(45,50)),
                             age_reported == '50-55'   ~ min(x = c(50,55)),
                             age_reported == '>50'     ~ min(x = c(50,65)),
                             age_reported == '55-60'   ~ min(x = c(55,60)),
                             age_reported == '60-65'   ~ min(x = c(60,65))))
ps$max <- with(ps, case_when(age_reported == '<1'      ~ max(x = c(0.1,1)),
                             age_reported == '0-1'     ~ max(x = c(0.1,1)),
                             age_reported == '0-4.9'   ~ max(x = c(0.1,4.9)),
                             age_reported == '0-5'     ~ max(x = c(0.1,5)),
                             age_reported == '<10'     ~ max(x = c(0.1,10)),
                             age_reported == '1-11'    ~ max(x = c(1,11)),
                             age_reported == '1-12'    ~ max(x = c(1,12)),
                             age_reported == '4-10'    ~ max(x = c(4,10)),
                             age_reported == '5-9.9'   ~ max(x = c(5,9.9)),
                             age_reported == '6-15'    ~ max(x = c(6,15)),
                             age_reported == '10-14.9' ~ max(x = c(10,14.9)),
                             age_reported == '10-20'   ~ max(x = c(10,20)),
                             age_reported == '15-19.9' ~ max(x = c(15,19.9)),
                             age_reported == '16-30'   ~ max(x = c(16,30)),
                             age_reported == '20-24.9' ~ max(x = c(20,24.9)),
                             age_reported == '20-42'   ~ max(x = c(20,42)),
                             age_reported == '25-34.9' ~ max(x = c(25,34.9)),
                             age_reported == '35-39.9' ~ max(x = c(35,39.9)),
                             age_reported == '25-29.9' ~ max(x = c(10,20)),
                             age_reported == '30-34.9' ~ max(x = c(30,34.9)),
                             age_reported == '20-34.9' ~ max(x = c(20,34.9)),
                             age_reported == '35-49.9' ~ max(x = c(35,49.9)),
                             age_reported == '>30'     ~ max(x = c(30,65)),
                             age_reported == '5-10'    ~ max(x = c(5,10)),
                             age_reported == '10-15'   ~ max(x = c(10,15)),
                             age_reported == '15-20'   ~ max(x = c(15,20)),
                             age_reported == '20-25'   ~ max(x = c(20,25)),
                             age_reported == '25-30'   ~ max(x = c(25,30)),
                             age_reported == '30-35'   ~ max(x = c(30,35)),
                             age_reported == '35-40'   ~ max(x = c(35,40)),
                             age_reported == '40-45'   ~ max(x = c(40,45)),
                             age_reported == '>40'     ~ max(x = c(40,65)),
                             age_reported == '45-50'   ~ max(x = c(45,50)),
                             age_reported == '50-55'   ~ max(x = c(50,55)),
                             age_reported == '>50'     ~ max(x = c(50,65)),
                             age_reported == '55-60'   ~ max(x = c(55,60)),
                             age_reported == '60-65'   ~ max(x = c(60,65))))
ps <- ps[!is.na(ps$min),]
unique(ps$age_reported)
# calves:  "0-4.9"   "0-1"   "0-5"   "<1" 
# juveniles: "5-9.9"   "5-10"   "4-10"  
# general young: "1-12"   "<10"   "1-11"  
# young pubescents:   "10-14.9"  "10-15" 
# old pubescents: "15-19.9" "15-20"  
# general pubescent:  "10-20" 
# young adults:"20-24.9"  "20-25"  "25-29.9, females"  "25-30, females"   "30-35, females" 
# mid adults:"25-34.9"   "35-39.9"   "30-34.9"  "20-34.9, females" "35-49.9, females"  "25-29.9, males"   "25-30, males"   "30-35 males"  "35-40"  "40-45 females"  "45-50 females" 
# old adults: ">50"  "40-45 males"  "45-50 males"  "50-55"  "55-60"   "60-65"  ">40 males" 
# confused:  "20-42"  "20-34.9, males"  "35-49.9, males"  "6-15"   "16-30"    ">30"   ">40 females" 

ps$sex <- ifelse(ps$sex == 'male', 'M', 'F')
ps$age <- with(ps, case_when(age_reported == '<1'      ~ '0-5',
                             age_reported == '0-1'     ~ '0-5',
                             age_reported == '0-4.9'   ~ '0-5',
                             age_reported == '0-5'     ~ '0-5',
                             age_reported == '4-10'    ~ '5-10',
                             age_reported == '5-9.9'   ~ '5-10',
                             age_reported == '5-10'    ~ '5-10',
                             age_reported == '10-14.9' ~ '10-15',
                             age_reported == '10-15'   ~ '10-15',
                             age_reported == '15-19.9' ~ '15-20',
                             age_reported == '15-20'   ~ '15-20',
                             age_reported == '20-24.9' & sex == 'M' ~ '20-25',
                             age_reported == '20-25'   & sex == 'M' ~ '20-25',
                             age_reported == '25-29.9' & sex == 'M' ~ '25-40',
                             age_reported == '25-30'   & sex == 'M' ~ '25-40',
                             age_reported == '35-39.9' & sex == 'M' ~ '25-40',
                             age_reported == '35-40'   & sex == 'M' ~ '25-40',
                             age_reported == '35-49.9' & sex == 'M' ~ '40+',    # not ideal
                             age_reported == '30-34.9' & sex == 'M' ~ '25-40',
                             age_reported == '20-34.9' & sex == 'M' ~ '25-40',  # not ideal
                             age_reported == '25-34.9' & sex == 'M' ~ '25-40',
                             age_reported == '30-35'   & sex == 'M' ~ '25-40',
                             age_reported == '40-45'   & sex == 'M' ~ '40+',
                             age_reported == '>40'     & sex == 'M' ~ '40+',
                             age_reported == '45-50'   & sex == 'M' ~ '40+',
                             age_reported == '50-55'   & sex == 'M' ~ '40+',
                             age_reported == '55-60'   & sex == 'M' ~ '40+',
                             age_reported == '60-65'   & sex == 'M' ~ '40+',
                             age_reported == '>50'     & sex == 'M' ~ '40+',
                             age_reported == '20-24.9' & sex == 'F' ~ '20-35',
                             age_reported == '20-25'   & sex == 'F' ~ '20-35',
                             age_reported == '25-29.9' & sex == 'F' ~ '25-40',
                             age_reported == '25-30'   & sex == 'F' ~ '25-40',
                             age_reported == '35-39.9' & sex == 'F' ~ '35-50',
                             age_reported == '35-40'   & sex == 'F' ~ '35-50',
                             age_reported == '35-49.9' & sex == 'F' ~ '35-50',
                             age_reported == '30-34.9' & sex == 'F' ~ '20-35',
                             age_reported == '20-34.9' & sex == 'F' ~ '20-35',
                             age_reported == '25-34.9' & sex == 'F' ~ '20-35',
                             age_reported == '30-35'   & sex == 'F' ~ '20-35',
                             age_reported == '40-45'   & sex == 'F' ~ '35-50',
                             age_reported == '>40'     & sex == 'F' ~ '35-50',  # not ideal
                             age_reported == '45-50'   & sex == 'F' ~ '35-50',
                             age_reported == '50-55'   & sex == 'F' ~ '50+',
                             age_reported == '55-60'   & sex == 'F' ~ '50+',
                             age_reported == '60-65'   & sex == 'F' ~ '50+',
                             age_reported == '>50'     & sex == 'F' ~ '50+',
                             age_reported == '10-20'   ~ 'NA', # can't differentiate pubescent age groups at all
                             age_reported == '<10'     ~ 'NA', # can't differentiate calf-juvenile age groups at all
                             age_reported == '1-11'    ~ 'NA', # can't differentiate calf-juvenile age groups at all
                             age_reported == '1-12'    ~ 'NA', # can't differentiate calf-juvenile age groups at all
                             age_reported == '20-42'   ~ 'NA', # can't differentiate adult age groups at all
                             age_reported == '16-30'   ~ 'NA', # can't differentiate adult-pubescent age groups at all
                             age_reported == '6-15'    ~ 'NA', # can't differentiate juvenile-pubescent age groups at all
                             age_reported == '>30'     ~ 'NA', # can't differentiate adult age groups at all
                             TRUE ~ as.character(age_reported)))
unique(ps$age)
ps <- ps[ps$age != 'NA',]
ps
ps_cens <- data.frame(age = rep(NA, nrow(ps)),
                      sex = rep(NA, nrow(ps)),
                      sid = rep(NA, nrow(ps)),
                      loc = rep(NA, nrow(ps)))
for(i in 1:nrow(ps)){
  ps_cens$age[i] <- list(rep(ps$age[i], ifelse(ps$estimate[i] < 1, 0, ps$estimate[i])))
  ps_cens$sex[i] <- list(rep(ps$sex[i], ifelse(ps$estimate[i] < 1, 0, ps$estimate[i])))
  ps_cens$sid[i] <- list(rep(ps$study_id[i], ifelse(ps$estimate[i] < 1, 0, ps$estimate[i])))
  ps_cens$loc[i] <- list(rep(ps$sitename[i], ifelse(ps$estimate[i] < 1, 0, ps$estimate[i])))
}
ps_cens_unlist <- data.frame(age = unlist(ps_cens$age),
                             alive = TRUE,
                             cens = 1,
                             sex = unlist(ps_cens$sex),
                             study_id = unlist(ps_cens$sid),
                             location = unlist(ps_cens$loc))
rm(ps, ps_cens)

ps_cens_unlist$k <- ifelse(ps_cens_unlist$age == '0-5', 1,
                           ifelse(ps_cens_unlist$age == '5-10', 2,
                                  ifelse(ps_cens_unlist$age == '10-15', 3,
                                         ifelse(ps_cens_unlist$age == '15-20', 4,
                                                ifelse(ps_cens_unlist$age == '20-25' | 
                                                         ps_cens_unlist$age == '20-35', 5,
                                                       ifelse(ps_cens_unlist$age == '25-40' | 
                                                                ps_cens_unlist$age == '35-50', 6, 7))))))
unique(ps_cens_unlist$k)

### create data lists 
males <- ps_cens_unlist[ps_cens_unlist$sex == 'M',]   ; males   <-   males[!is.na(males$age),]
females <- ps_cens_unlist[ps_cens_unlist$sex == 'F',] ; females <- females[!is.na(females$age),]

# Males
N_m <- nrow(males)                     # Number of individuals
K <- 7                                 # Number of age bin categories
males_ls <- list(
  N = N_m,
  K = K,
  age_category_index = males$k)        # age category

# Females
N_f <- nrow(females)                   # Number of individuals
females_ls <- list(
  N = N_f,
  K = K,
  age_category_index = females$k)      # age category

### 1) Fit male model of ages 
# Shape = 1.29, scale = 29.5 (so we will use shape = 1.3 and scale = 30)
aedd_male_age_ord_mod <- cmdstan_model("models/age_estimation/male_age_ordreg_22.03.31.stan")
male_age_fit <- aedd_male_age_ord_mod$sample(
  data = males_ls, 
  chains = 4, 
  parallel_chains = 1,
  iter_sampling = 2000
)

# Examine the estimates.
age_est_m <- male_age_fit$summary()[906:1809,]
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

### 2) Fit female model of ages
# Shape = 1.25, scale = 32.9 (so we will use shape = 1.3 and scale = 33)
aedd_female_age_ord_mod <- cmdstan_model("models/age_estimation/female_age_ordreg_22.03.31.stan")
female_age_fit <- aedd_female_age_ord_mod$sample(
  data = females_ls, 
  chains = 4, 
  parallel_chains = 1,
  iter_sampling = 2000
)

# Examine the estimates.
age_est_f <- female_age_fit$summary()[N_f+2:2*N_f+3,]
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











# Examine default brms priors (all elephants together)
bfit_default <- age_cens %>%
  brm(age_non0 | cens(censor) ~ 1,
      data = ., family = "weibull", 
      chains = 4, cores = 4)
prior_summary(bfit_default)
# Defaults
#intercept - student_t(3, 2.9, 2.5)
#shape - gamma(0.01, 0.01)
plot(bfit_default) # shape max >600


# Run model for male elephants with slightly better priors
bfit <- age_cens %>% 
  brm(age_non0 | cens(censor) ~ 1,
      prior = c(
        prior(student_t(3,5,5), class = Intercept),
        prior(exponential(0.5), class = shape)
      ),
      data = ., family = "weibull", 
      chains = 4, cores = 4)
plot(bfit)  # looks very similar shape to the default one, but much shorter tail

# extract scale and shape from posterior draws
# to convert to scale we need to both undo the link function by taking the exponent
# and then refer to the brms documentation to understand how the mean relates to the scale
bfit_draws <- as_draws_df(bfit) %>%
  mutate(scale = exp(b_Intercept) / (gamma(1 + 1 / shape)))

est_shape <- mean(bfit_draws$shape)  # 3.47
est_scale <- mean(bfit_draws$scale)  # 2.53

# posterior predictive plot based on the mean
probs <- dweibull(1:100, shape = est_shape, scale = est_scale)
probs <- probs / sum(probs)
plot(probs) 

# simulate ages for 3k elephants
hist(rweibull(3000,est_shape,est_scale)) # horrendous -- looks like a normal and the age values are insane (peaks at 2500000000000000000000000000000000000000 years old!)

# Shape = 3.47, scale = 2.53 BUT histogram and distribution look terrible

### Run on MOTNP data ####
# read in male elephant data
motnp <- read_csv('data_processed/motnp_elenodes_22.01.13.csv') %>% 
  filter(sex == 'M')

motnp$age_cat <- ifelse(motnp$age_category == '40+', '7',
                        ifelse(motnp$age_category == '25-40', '6',
                               ifelse(motnp$age_category == '20-25', '5',
                                      ifelse(motnp$age_category == '16-20', '4',
                                             ifelse(motnp$age_category == '10-15', '3',
                                                    ifelse(motnp$age_category == '6-9', '2', '1'))))))

# create data list
N <- nrow(motnp)                      # Number of individuals
K <- 7                                # Number of age bin categories
elephants_ls <- list(
  N = N,
  K = K,
  age_category_index = motnp$age_cat) # Age category

# fit Stan model to estimate ages from the categories based on a Weibull distribution (for real data these will be centred on parameters from fitting to age estimates based on another elephant population)
anp_age_ord_mod <- cmdstan_model("models/age_estimation/male_age_ordreg_test_22.04.04.stan")

age_estimation_fit <- anp_age_ord_mod$sample(
  data = elephants_ls, 
  chains = 4, 
  parallel_chains = 1,
  iter_sampling = 2000
)

# Examine the estimates.
age_est_mat <- age_estimation_fit$summary()[1640:3277,]
plot_data <- data.frame(age_cat = elephants_ls$age,                       # Biologists original age est
                        age_true = males$age[!is.na(males$age_cat)],      # real age
                        model_age = age_est_mat$mean)                     # Mean modelled age
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
df %>% ggplot(aes(x=jitter(as.numeric(true_age)), y=value, group=factor(id))) +
  geom_point(size=2,alpha=0.1) +
  #geom_vline(xintercept=c(5,10,15,20,25,40), linetype="dashed", alpha=0.6) +
  geom_hline(yintercept=c(5,10,15,20,25,40), linetype="dashed", alpha=0.6) +
  theme_bw() + 
  xlab("Assigned age") + ylab("Modelled age")
summary(df$value[which(df$true_age == 7)])
quantile(df$value[which(df$true_age == 7)], seq(0,1,0.01))

### males only -- old ####