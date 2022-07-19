#### Information #####
# To estimate the ages of individuals in other populations, first come up with a survival curve for the Amboseli elephants which have good birth and death data.
# survival curve options = exponential/logistic/Gompertz/Weibull, where shape = simple/bathtub/Makeham
# follow through Mia's cde to determine which curve will best fit Amboseli data. she ran it on a high performance cluster, so it’s set up to run in parallel, but the basta model can also be run in parallel on its own (just like Stan). we will need to implement it as a custom function in Stan for the final analysis, given that we need to define it in the prior anyway. 

#### Load packages ####
lapply(c("BaSTA"), require, character.only = TRUE) # required packages for the running the bayesian hierarchichal function
library(tidyverse)

#### Read in Amboseli data ####
males <- readxl::read_excel('data_raw/Raw_ATE_AllElephants_Lee220118.xlsx') %>% 
  janitor::clean_names() %>% 
  filter(sex == 'Male') %>% 
  mutate('casename' = id) %>% 
  select(casename, birth_year, death_year, age_death)
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
    males[i,j+4] <- ifelse(length(which(casename$year == j+1971)) > 0, 1, 0)
    rm(casename)
  }
}

#### Prepare input data ####
# Study start and end years
Ss <- 1972  # start of study
Se <- 2021  # end of study

# prep data
basta_data <- males
basta_data$count_years <- rowSums(basta_data[5:54])
sum(table(basta_data$count_years)[2:35])

basta_data$dyr <- ifelse(basta_data$dyr < 100,  0, basta_data$dyr)
basta_data$byr <- ifelse(basta_data$byr < 1972, 0, basta_data$byr)

write_csv(basta_data, 'data_processed/anp_agesightings_basta_22.05.30.csv')

old <- readxl::read_excel('data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% 
  janitor::clean_names() %>% 
  select(casename, byr, dyr)
old$casename <- paste0('M', old$casename)
old <- old[-c(1:4,6,8,11,14:18,22,23,26:30,35,38,44:46,52,55,56,58,59,62,63,66,68,71,72,76,83,87,90,93,94,96,105,107,108,111,114,115,118,127,128,132,133,135,136,138,145,146,151,155,159,162,174,179,182,186,187,190,191,193,197,201,208,210,213,220,224,230,235,238,244,248,250,254,256,259,261,263,266,268,288,291,292,295,298,306,307,309,319,327,339,342,346,370,384,416,417,435,449,472,717,693,708,735,769,774),]  # elephants that have sightings before birth/after death
old_inc <- old[old$casename %in% sort(unique(paste0('M',sightings$id))),]
old_exc <- anti_join(old, old_inc)
old_exc$age_death <- ifelse(old_exc$dyr < 10,
                            paste0('living ', 2021 - old_exc$byr),
                            old_exc$dyr - old_exc$byr)
table(old_exc$age_death)

old_readin <- read_csv('data_processed/anp_agesightings_basta_22.05.24.csv')
old_readin$count_years <- rowSums(old_readin[4:53])
table(old_readin$count_years)

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

mortality <- weibull(b0 = 1.3, b1 = 30, age = 1:70)      ; plot(mortality)

pM <-  c(1.3,30)    # prior means for parameters in a Weibull model with a simple shape
pSD <- c(0.2, 5)    # prior standard deviation for parameters in a Weibull model with a simple shape -- smaller values makes no visible difference to final curves
Wb.Si2 <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = 1, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               model = "WE", shape = c("simple"))
summary(Wb.Si2)
plot(Wb.Si2)                       # Plot traces
plot(Wb.Si2, plot.trace = FALSE)   # Plot survival and mortality curves

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

pM <-  c(1.3,30)    # prior means for parameters in a Weibull model with a simple shape
pSD <- c(0.2, 5)    # prior standard deviation for parameters in a Weibull model with a simple shape -- smaller values makes no visible difference to final curves
Wb.Si2 <- basta(basta_data, studyStart = Ss,
                studyEnd = Se, nsim = sim, niter = iter,
                burnin = warmup, thinning = thin,
                thetaPriorMean = pM,thetaPriorSd = pSD,
                model = "WE", shape = c("simple"))
summary(Wb.Si2)                    # DIC lower than for priors given above, and output estimates look no better 
plot(Wb.Si2)                       # Plot traces
plot(Wb.Si2, plot.trace = FALSE)   # Plot survival and mortality curves

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
rm(ages, amr, basta_data, casename, males, sightings, iter, j, N, pM, pSD, Se, sim, Ss, thin, warmup)
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

#### read in AEDD data -- need to use just population size, if estimate = 8 that's 8 individuals of that age, not 1 ####
aedd <- read_csv('data_raw/AfricanElephantDemographicDatabase_aedd_22.04.04.csv')
aedd <- aedd[aedd$species == 'L. africana' & !is.na(aedd$species),]
sort(unique(aedd$metric))

## studies of interest
aedd <- aedd[aedd$metric == 'population size', c(1,6,7,8,13,14,22)]
colnames(aedd)

## population size: estimate = number observed to reach minimum age (censor at min or median?)
ps <- aedd[!is.na(aedd$age_reported),]
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
ps$age_censor <- ifelse(ps$age == '0-5', 1,
                        ifelse(ps$age == '5-10', 5,
                               ifelse(ps$age == '10-15', 10,
                                      ifelse(ps$age == '15-20', 15,
                                             ifelse(ps$age == '20-25', 20,
                                                    ifelse(ps$age == '25-40', 25,
                                                           ifelse(ps$age == '40+', 40,
                                                                  ifelse(ps$age == '20-35', 20,
                                                                         ifelse(ps$age == '35-50', 35,
                                                                                ifelse(ps$age == '50+', 50, NA))))))))))

ps_cens <- data.frame(age = rep(NA, nrow(ps)),
                      sex = rep(NA, nrow(ps)),
                      sid = rep(NA, nrow(ps)),
                      loc = rep(NA, nrow(ps)))
for(i in 1:nrow(ps)){
  ps_cens$age[i] <- list(rep(ps$age_censor[i], ifelse(ps$estimate[i] < 1, 0, ps$estimate[i])))
  ps_cens$sex[i] <- list(rep(ps$sex[i], ifelse(ps$estimate[i] < 1, 0, ps$estimate[i])))
  ps_cens$sid[i] <- list(rep(ps$study_id[i], ifelse(ps$estimate[i] < 1, 0, ps$estimate[i])))
  ps_cens$loc[i] <- list(rep(ps$sitename[i], ifelse(ps$estimate[i] < 1, 0, ps$estimate[i])))
}
age_cens <- data.frame(age = unlist(ps_cens$age),
                       alive = TRUE,
                       cens = 1,
                       sex = unlist(ps_cens$sex),
                       study_id = unlist(ps_cens$sid),
                       location = unlist(ps_cens$loc))
rm(ps, ps_cens)

age_cens$k <- ifelse(age_cens$age == 1, 1,
                     ifelse(age_cens$age == 5, 2,
                            ifelse(age_cens$age == 10, 3,
                                   ifelse(age_cens$age == 15, 4,
                                          ifelse(age_cens$age == 20, 5,
                                                 ifelse(age_cens$age == 25 | 
                                                          age_cens$age == 35, 6, 7))))))
unique(age_cens$k)

# clean up
colnames(age_cens)[3] <- 'censor'
age_cens$sex <- ifelse(is.na(age_cens$sex), 'U', ifelse(age_cens$sex == 'female', 'F', 'M'))
age_cens$location_id <- as.integer(as.factor(age_cens$location))
sort(unique(age_cens$location)) # should we remove the Amboseli study from this given that that is already included in the prior, or retain it?

#### improve Weibull distribution parameters for prior -- https://statwonk.com/bayesian-right-censored-weibull-model.html ####
length(which(age_cens$age == 0))

# Priors used
summary(Wb.Si)
num_samples <- 500
prior_intercept <- rnorm(num_samples, Wb.Si$coefficients[3,1], Wb.Si$coefficients[3,2]) # parameter pi.1972
prior_shape <- rnorm(num_samples, Wb.Si$coefficients[2,1], Wb.Si$coefficients[2,2])     # parameter b1
prior_scale1 <- exp(prior_intercept) / (gamma(1 + 1 / prior_shape))                     # Formula Dan previously used. Output values are REALLY tiny
prior_scale2 <- rnorm(num_samples, Wb.Si$coefficients[1,1], Wb.Si$coefficients[1,2])    # parameter b0

hist(prior_intercept)
hist(prior_shape)
hist(prior_scale1)
hist(prior_scale2)

probs <- dweibull(1:100, shape = prior_shape[1], scale = prior_scale2[1])
probs <- probs / sum(probs)
plot(probs, type = "l", lwd = 0.5, ylim = c(0,1), las = 1)
for( i in 2:num_samples) {
  probs <- dweibull(1:100, shape = prior_shape[i], scale = prior_scale2[i])
  probs <- probs / sum(probs)
  lines(probs, lwd=0.5)
} # this just shows that everything will be dead by age 5? no variation at all -- coefficients are SE not SD --> tried multiplying by sqrt(iter/2) but gave same shape graphs

# Run model for male elephants with Amboseli priors
exp(Wb.Si$coefficients[3,1]) # pi.1972 mean -- exponentiated because one of Dan's emails mentioned that these would be on the log scale, but I actually think he was walking about something else! These don't make sense as reasonable parameters though... log scale or not!
exp(Wb.Si$coefficients[3,2]) # pi.1972 SE
exp(Wb.Si$coefficients[2,1]) # b1 mean
exp(Wb.Si$coefficients[2,2]) # b1 SD

hist(age_cens$age)

bfit <- age_cens %>% 
  brm(age | cens(censor) ~ 1,
      prior = c(
        prior(normal(0.4928, 0.00476), class = Intercept), # values from pi.1972
        prior(normal(0.0318, 0.000585), class = shape)     # values from b1
      ),
      data = ., family = "weibull", 
      chains = 4, cores = 4)
plot(bfit)

as_sample_tibble(bfit) %>%
  mutate_at(vars(b_Intercept), exp) %>%
  ggplot(aes(x = b_Intercept, y = shape)) +
  geom_density_2d(colour = "red", size = 1) +
  geom_point(alpha = 0.1) +
  #geom_hline(yintercept = 3.35, colour = "blue") +
  #geom_vline(xintercept = 170, colour = "blue") +
  geom_hline(yintercept = 0.0835, colour = "blue") +
  geom_vline(xintercept = 1.6875, colour = "blue") +
  xlab("Scale") +
  ylab("Shape") +
  ggtitle("Right-censored weibull(shape, scale)") +
  theme_fivethirtyeight()

# extract scale and shape from posterior draws
# to convert to scale we need to both undo the link function by taking the exponent
# and then refer to the brms documentation to understand how the mean relates to the scale
bfit_draws <- as_draws_df(bfit) %>%
  mutate(scale = exp(b_Intercept) / (gamma(1 + 1 / shape)))

est_shape <- mean(bfit_draws$shape)  # 0.40977
est_scale <- mean(bfit_draws$scale)  # 0.57311

# posterior predictive plot based on the mean
probs <- dweibull(1:100, shape = est_shape, scale = est_scale)
probs <- probs / sum(probs)
plot(probs, type = 'l') ; abline(h = 0, lty = 2) # all dead by 40, almost all by 20

# simulate ages for 3k elephants
hist(rweibull(3000,est_shape,est_scale), ylim = c(0,20),
     #xlim = c(0,50),
     breaks = 1000) # early death rate much too high, up to 150 years

#### I don't think it works to use basta again to refine the parameters, because these are also categories -- ignore this section ####
rethinking::dens(rweibull(1000, Wb.Si$coefficients[2,1], Wb.Si$coefficients[1,1]))
rethinking::dens(rweibull(1000, Wb.Si$coefficients[1,1], Wb.Si$coefficients[2,1]))

weibull <- function(b0, b1, age){
  y <- (b0) * (b1^b0) * (age^(b0-1))
  return(y)
}
mortality <- weibull(b0 = Wb.Si$coefficients[1,1],
                     b1 = Wb.Si$coefficients[2,1],
                     age = 1:70)
plot(mortality)

max(age_cens$age)
aedd_basta <- data.frame(id = 1:nrow(age_cens),
                         byr = 2000-age_cens$age,
                         dyr = rep(0,nrow(age_cens)))
aedd_basta$obs_1950 <- NA ; aedd_basta$obs_1951 <- NA ; aedd_basta$obs_1952 <- NA ; aedd_basta$obs_1953 <- NA ; aedd_basta$obs_1954 <- NA ; aedd_basta$obs_1955 <- NA ; aedd_basta$obs_1956 <- NA ; aedd_basta$obs_1957 <- NA ; aedd_basta$obs_1958 <- NA ; aedd_basta$obs_1959 <- NA ; aedd_basta$obs_1960 <- NA ; aedd_basta$obs_1961 <- NA ; aedd_basta$obs_1962 <- NA ; aedd_basta$obs_1963 <- NA ; aedd_basta$obs_1964 <- NA ; aedd_basta$obs_1965 <- NA ; aedd_basta$obs_1966 <- NA ; aedd_basta$obs_1967 <- NA ; aedd_basta$obs_1968 <- NA ; aedd_basta$obs_1969 <- NA ; aedd_basta$obs_1970 <- NA ; aedd_basta$obs_1971 <- NA ; aedd_basta$obs_1972 <- NA ; aedd_basta$obs_1973 <- NA ; aedd_basta$obs_1974 <- NA ; aedd_basta$obs_1975 <- NA ; aedd_basta$obs_1976 <- NA ; aedd_basta$obs_1977 <- NA ; aedd_basta$obs_1978 <- NA ; aedd_basta$obs_1979 <- NA ; aedd_basta$obs_1980 <- NA ; aedd_basta$obs_1981 <- NA ; aedd_basta$obs_1982 <- NA ; aedd_basta$obs_1983 <- NA ; aedd_basta$obs_1984 <- NA ; aedd_basta$obs_1985 <- NA ; aedd_basta$obs_1986 <- NA ; aedd_basta$obs_1987 <- NA ; aedd_basta$obs_1988 <- NA ; aedd_basta$obs_1989 <- NA ; aedd_basta$obs_1990 <- NA ; aedd_basta$obs_1991 <- NA ; aedd_basta$obs_1992 <- NA ; aedd_basta$obs_1993 <- NA ; aedd_basta$obs_1994 <- NA ; aedd_basta$obs_1995 <- NA ; aedd_basta$obs_1996 <- NA ; aedd_basta$obs_1997 <- NA ; aedd_basta$obs_1998 <- NA ; aedd_basta$obs_1999 <- NA ; aedd_basta$obs_2000 <- NA

for(i in 1:nrow(aedd_basta)){
  for(j in 4:ncol(aedd_basta)){
    aedd_basta[i,j] <- ifelse(aedd_basta$byr[i] < j + 1946, 1, 0)
  }
}

pM <-  c(Wb.Si$coefficients[1,1], Wb.Si$coefficients[2,1])
pSD <- c(Wb.Si$coefficients[1,2], Wb.Si$coefficients[2,2])

sim <- 1 ; iter <- 500 ; warmup <- iter/2 ; thin <- 10

Wb.Si.updated <- basta(aedd_basta, studyStart = 1950,
                       studyEnd = 2000, nsim = sim, niter = iter,
                       burnin = warmup, thinning = thin,
                       thetaPriorMean = pM,thetaPriorSd = pSD,
                       model = "WE", shape = c("simple"))
summary(Wb.Si.updated)
plot(Wb.Si.updated)                       # Plot traces
plot(Wb.Si.updated, plot.trace = FALSE)   # Plot survival and mortality curves

#### improve Gompertz bathtub distribution parameters for prior ####
# Priors used
summary(Go.Bt)
num_samples <- 500
prior_a0 <- rnorm(num_samples, Go.Bt$coefficients[1,1], Go.Bt$coefficients[1,2])        # parameter a1
prior_a1 <- rnorm(num_samples, Go.Bt$coefficients[2,1], Go.Bt$coefficients[2,2])        # parameter a0
prior_c  <- rnorm(num_samples, Go.Bt$coefficients[3,1], Go.Bt$coefficients[3,2])        # parameter c
prior_b0 <- rnorm(num_samples, Go.Bt$coefficients[4,1], Go.Bt$coefficients[4,2])        # parameter b1
prior_b1 <- rnorm(num_samples, Go.Bt$coefficients[5,1], Go.Bt$coefficients[5,2])        # parameter b0
prior_intercept <- rnorm(num_samples, Go.Bt$coefficients[6,1], Go.Bt$coefficients[6,2]) # parameter pi.1972

hist(prior_a0)
hist(prior_a1)
hist(prior_c)
hist(prior_b0)
hist(prior_b1)
hist(prior_intercept)

gompertz <- function(b0, b1, age){
  y <- exp(b0 + b1*age)
  return(y)
}
gompertz_bt <- function(a0, a1, c, b0, b1, age){
  gompertz <- gompertz(b0, b1, age)
  bathtub <- exp(a0 - a1*age) + c + gompertz
  return(bathtub)
}
mortality <- gompertz_bt(a0 = Go.Bt$coefficients[1,1],
                         a1 = Go.Bt$coefficients[2,1],
                         c  = Go.Bt$coefficients[3,1],
                         b0 = Go.Bt$coefficients[4,1],
                         b1 = Go.Bt$coefficients[5,1],
                         1:60)  ; plot(mortality) # no infant mortality, could this be because ATE only supplied sightings data for post-independence?

# Run model for male elephants with Amboseli priors
Go.Bt$coefficients[1,1]
Go.Bt$coefficients[1,2]
Go.Bt$coefficients[2,1]
Go.Bt$coefficients[2,2]
Go.Bt$coefficients[3,1]
Go.Bt$coefficients[3,2]
Go.Bt$coefficients[4,1]
Go.Bt$coefficients[4,2]
Go.Bt$coefficients[5,1]
Go.Bt$coefficients[5,2]

hist(age_cens$age)

bfit <- age_cens %>% 
  brm(age | cens(censor) ~ 1,
      prior = c(
        prior(normal(), class = Intercept), # values from pi.1972
        prior(normal(), class = shape)      # values from b0
      ),
      data = ., family = "SIMPLE GOMPERTZ AND BOTH BATHTUBS DON'T APPEAR TO BE OPTIONS", 
      chains = 4, cores = 4)
plot(bfit)

as_sample_tibble(bfit) %>%
  mutate_at(vars(b_Intercept), exp) %>%
  ggplot(aes(x = b_Intercept, y = shape)) +
  geom_density_2d(colour = "red", size = 1) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = 3.35, colour = "blue") +
  geom_vline(xintercept = 170, colour = "blue") +
  xlab("Scale") +
  ylab("Shape") +
  ggtitle("Right-censored weibull(shape, scale)") +
  theme_fivethirtyeight()

# extract scale and shape from posterior draws
# to convert to scale we need to both undo the link function by taking the exponent
# and then refer to the brms documentation to understand how the mean relates to the scale
bfit_draws <- as_draws_df(bfit) %>%
  mutate(scale = exp(b_Intercept) / (gamma(1 + 1 / shape)))

est_shape <- mean(bfit_draws$shape)  # 0.40977
est_scale <- mean(bfit_draws$scale)  # 0.57311

# posterior predictive plot based on the mean
probs <- dweibull(1:100, shape = est_shape, scale = est_scale)
probs <- probs / sum(probs)
plot(probs, type = 'l') ; abline(h = 0, lty = 2) # all dead by 40, almost all by 20

# simulate ages for 3k elephants
hist(rweibull(3000,est_shape,est_scale)) # early death rate much too high, up to 150 years

#### Age estimation for mega data set based on Amboseli distribution -- this is just to test quality of prior, not actually change anything in the values. Run this method for MOTNP after testing. ####
### create data lists 
males <- age_cens[age_cens$sex == 'M',]   ; males   <-   males[!is.na(males$age),]

# Males
N_m <- nrow(males)                     # Number of individuals
K <- 7                                 # Number of age bin categories
males_ls <- list(
  N = N_m,
  K = K,
  age_category_index = males$k)        # age category

### 1) Fit male model of ages 
# Shape = 0.41, scale = 0.57
aedd_male_age_ord_mod <- cmdstan_model("models/age_estimation/male_age_ordreg_22.05.24.stan")
male_age_fit <- aedd_male_age_ord_mod$sample(
  data = males_ls, 
  chains = 4, 
  parallel_chains = 1,
  iter_sampling = 2000
)

# Examine the estimates.
#age_est_m <- male_age_fit$summary()[906:1809,]
age_est_m <- male_age_fit$summary()[2051:4099,]
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
  scale_y_continuous(limits = c(0,100)) +
  theme_bw() + 
  xlab("Assigned age") + ylab("Modelled age")
summary(df_m$value[which(df_m$true_age == 7)])
quantile(df_m$value[which(df_m$true_age == 7)], seq(0,1,0.01))

### Run on MOTNP data ####
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

# create data list
N <- nrow(motnp)                      # Number of individuals
K <- 7                                # Number of age bin categories
elephants_ls <- list(
  N = N,
  K = K,
  age_category_index = motnp$age_cat) # Age category

# fit Stan model to estimate ages from the categories based on a Weibull distribution (for real data these will be centred on parameters from fitting to age estimates based on another elephant population)
anp_age_ord_mod <- cmdstan_model("models/age_estimation/male_age_ordreg_test_22.05.25.stan")

age_estimation_fit <- anp_age_ord_mod$sample(
  data = elephants_ls, 
  chains = 4, 
  parallel_chains = 1,
  iter_sampling = 2000
)

# Examine the estimates.
age_est_mat <- age_estimation_fit$summary()[247:491,]
plot_data <- data.frame(age_cat = elephants_ls$age,                       # Biologists original age est
                        model_age = age_est_mat$mean)                     # Mean modelled age
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
