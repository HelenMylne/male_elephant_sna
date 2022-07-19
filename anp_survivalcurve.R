#### information #####
# To estimate the ages of individuals in other populations, first come up with a survival curve for the Amboseli elephants which have good birth and death data.
# survival curve options = exponential/logistic/Gompertz/Weibull, where shape = simple/bathtub/Makeham
# follow through Mia's code to determine which curve will best fit Amboseli data. she ran it on a high performance cluster, so it’s set up to run in parallel, but the basta model can also be run in parallel on its own (just like Stan). we will need to implement it as a custom function in Stan for the final analysis, given that we need to define it in the prior anyway. 

#### load packages ####
lapply(c("BaSTA"), require, character.only = TRUE) # required packages for the running the bayesian hierarchichal function
library(tidyverse)
library(cmdstanr)
library(ggdist)
library(posterior)
library(bayesplot)
library(survival)
library(tidybayes)
library(ggthemes)
library(data.table)

#### read in Amboseli data ####
sightings <- read_csv('data_processed/anp_sightings_updated_22.06.22.csv')
sightings$year <- lubridate::year(sightings$obs_date)
sightings <- sightings[,c('casename','year')]

males <- readxl::read_excel('data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% 
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

#### prepare input data ####
# prep data
basta_data <- males
basta_data$count_years <- rowSums(basta_data[4:53])
summary(basta_data$count_years) # no individuals never observed, as it should be. No NAs.

# set birth and death years
basta_data$birth_year <- ifelse(basta_data$byr < 1972, 0, basta_data$byr)
basta_data$death_year <- ifelse(basta_data$dyr < 100,  0, basta_data$dyr)

basta_data <- basta_data[,c(1,55,56,4:53,2:3,54)]

#### remove elephants with incorrect sightings ####
check_birth <- basta_data[c(638,648,665,682,686),] # elephants recorded before birth
check_death <- basta_data[c(3,5,7,13,15,16,17,44,54,55,57,61,65,67,71,82,87,89,92,95,104,106,107,110,113,117,126,127,131,132,134,135,137,145,150,154,158,161,178,181,186,189,190,207,209,218,228,242,246,252,254,256,263,288,291,294,323,335,338,365,378,410,411,429,653),]    # elephants recorded after death
basta_data <- anti_join(basta_data, check_birth, by = 'casename')
basta_data <- anti_join(basta_data, check_death, by = 'casename')
rm(check_birth, check_death)

#### convert birth years to year first sighted ####
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
hist(basta_data$age_first_sighted)

# remove all individuals first observed at <8 years old or > 20 years old, but retaining those already over 20 at the strat of the study to avoid truncating maximum age data.
basta_data <- basta_data[,c(1,57,2,54,3:53,55:56,58)]
basta_data$remove <- ifelse(basta_data$byr > 1962 & basta_data$age_first_sighted > 20 | 
                              basta_data$age_first_sighted < 8, 'remove', 'keep')
table(basta_data$age_first_sighted[which(basta_data$age_first_sighted < 8)])
table(basta_data$remove)
basta_data <- basta_data[basta_data$remove == 'keep',c(1:58)]

# set "birth year" to the year first sighted, with 0 if already of independent age in 1972.
basta_data$byr_apparent <- ifelse(basta_data$byr < 1962, 0,
                                  ifelse(basta_data$year_first_sighted == 1971, 0,
                                         basta_data$year_first_sighted))

basta_data <- basta_data[,c(1,59,5:55,2:4,56:58)]
basta_anp <- basta_data[,c(1:53)]

#### prior predictive check ####
gompertz_bt <- function(a0, a1, c, b0, b1, age){
  gompertz <- gompertz(b0, b1, age)
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

#### run model -- DIC = 10162 / 9601.075 (without cutting out very early/late dispersers / cut out) ####
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

# following Dan's previous method for simulating population
mortality <- gompertz_bt(a0 = Go.Bt$coefficients[1,1],
                         a1 = Go.Bt$coefficients[2,1],
                         c  = Go.Bt$coefficients[3,1],
                         b0 = Go.Bt$coefficients[4,1],
                         b1 = Go.Bt$coefficients[5,1], 1:50)
probs <- 1 - (mortality/max(mortality))
plot(probs)

# Now we create a fictional population with ages selected from that distribution. We also add an error for age estimation.
N <- 700
elephants_ls <- list(
  N = N,
  age = sample(1:50, N, replace = TRUE, prob = probs)
)
hist(elephants_ls$age)

# Next simulate observing ages with error and then binning them into age groups
E <- 3 # Error (in either direction) of age estimates
elephants_ls$age_guess <- elephants_ls$age + sample(-E:E, N, replace=TRUE)
elephants_ls$age_category_index <- sapply(elephants_ls$age_guess, function(x) which.max(x < c(5, 10, 15, 20, 25, 40, Inf)))
hist(elephants_ls$age_guess)
hist(elephants_ls$age_category_index)

# Let's look at the actual age verses the biologist assigned age and look at the thresholds. You can see that the chance of being mis-classified is lower if an actual age is in the middle of the age category. Also, there is mostly a bias towards ages within the age class actually being towards the lower end of the class.
data.frame(elephants_ls) %>%
  ggplot(aes(x=age, y=age_guess, col=factor(age_category_index))) +
  geom_point(size=4,alpha=0.6) +
  geom_vline(xintercept=c(5, 10, 15, 20, 25, 40, 60), col=factor(1:7), linetype="dashed", alpha=0.6) +
  theme_minimal() + 
  theme(legend.position = 'none')+
  xlab("Assigned age") + ylab("Age")

