#### Load packages ####
lapply(c("foreach", "doParallel"), require, character.only = TRUE) # required if you want to run the model using parallel computing
lapply(c("BaSTA"), require, character.only = TRUE) # required packages for the running the bayesian hierarchichal function
library(tidyverse)

## read in Amboseli data ####
males <- readxl::read_excel('data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% 
  select(CASENAME, BYR, DYR) %>% 
  janitor::clean_names()
str(males)

# generate numeric variable for age in 2021 of living elephants or age at death
males$age <- as.numeric(ifelse(males$dyr < 10, 2021 - males$byr, males$dyr - males$byr))
table(males$age)
hist(males$age, breaks = 75)    # overall age distribution

# create new variable
males$censor <- ifelse(males$dyr < 10, 'TRUE', 'FALSE')    # censor = TRUE when elephant is still alive
table(males$censor)
nrow(males) - length(unique(males$casename))               # 0 -- every elephant has a unique number

## create data frame for basta ####
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

sightings <- readxl::read_excel('data_raw/Raw_ATE_Sightings_Fishlock220224.xlsx') %>% janitor::clean_names()
sightings$year <- lubridate::year(sightings$obs_date)
sightings <- sightings[,c('casename','year')]
colnames(sightings)[1] <- 'id'
tapply(X = sightings$id, INDEX = sightings$year, FUN = table)

(N <- 2022-1972)
for(i in 1:nrow(males)){
  for(j in 1:N){
    casename <- sightings[sightings$id == males$casename[i],]
    males[i,j+5] <- ifelse(length(which(casename$year == j+1971)) > 0, 1, 0)
  }
}

## Preparing input data ####
# Getting start and end times from the input datasets - will depend on your input if this works, but can also manually assign the year the observations began and ended
Ss <- 1972  # start of study
Se <- 2021  # end of study

# Set the arguments for the iterations
iter <- 500   # number of iterations
sim <- 1         # number of chains -- when this is set to 1, it works fine. When set to 4, the entire thing produces an error of "Error in rep(0, nbd) : invalid 'times' argument". I know this has something to do with updating the jump standard deviations in the MCMC chain, but I'm not 100% sure I know what that means (is this like the step size/initial energy of the ball on the posterior distribution?) and therefore why Mia had it set to 4. 
warmup <- iter/2 # number of iterations to use during warmup/burn-in period
thin <- 10      # thinning number

# Set the prior values
pM <-  c(-3,0.2,0.001,-4,0.05) # prior means for parameters in a Gompertz model with a bathtub shape
pSD <- c(1,0.05,0.001,1,0.01)  # prior standard deviation for parameters in a Gompertz model with a bathtub shape

basta_data <- males
basta_data$dyr <- ifelse(basta_data$censor == 'TRUE', 0, basta_data$dyr)
basta_data$casename <- paste0('M', basta_data$casename)
str(basta_data)
basta_data <- basta_data[-c(1:4,6,8,11,14:18,22,23,26:30,35,38,44:46,52,55,56,58,59,62,63,66,68,71,72,76,83,87,90,93,94,96,105,107,108,111,114,115,118,127,128,132,133,135,136,138,145,146,151,155,159,162,174,179,182,186,187,190,191,193,197,201,208,210,213,220,224,230,235,238,244,248,250,254,256,259,261,263,266,268,288,291,292,295,298,306,307,309,319,327,339,342,346,370,384,416,417,435,449,472,717,693,708,735,769,774),  # elephants that have sightings before birth/after death
                         c(1:3,6:55)]

basta_data <- basta_data[basta_data$casename %in% sort(unique(paste0('M',sightings$id))),] # remove any elephants without sightings

Out <- basta(basta_data, studyStart = Ss,
             studyEnd = Se, nsim = sim, niter = iter, 
             burnin = warmup, thinning = thin, 
             thetaPriorMean = pM,thetaPriorSd = pSD,
             models = c("GO"), shape = c("bathtub"))
str(basta_data)
summary(Out)
plot(Out)
plot(Out, plot.trace = FALSE)   # Plot survival and mortality curves
