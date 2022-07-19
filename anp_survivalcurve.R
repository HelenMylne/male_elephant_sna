#### Information #####
# To estimate the ages of individuals in other populations, first come up with a survival curve for the Amboseli elephants which have good birth and death data.
# survival curve options = exponential/logistic/Gompertz/Weibull, where shape = simple/bathtub/Makeham
# follow through Mia's cde to determine which curve will best fit Amboseli data. she ran it on a high performance cluster, so it’s set up to run in parallel, but the basta model can also be run in parallel on its own (just like Stan). we will need to implement it as a custom function in Stan for the final analysis, given that we need to define it in the prior anyway. 

#### Load packages ####
lapply(c("BaSTA"), require, character.only = TRUE) # required packages for the running the bayesian hierarchichal function
library(tidyverse)

#### Mia Nielsen example code ####
# prep data and model ####
### Model survival based on census data using a Bayesian framework (BaSTA package)
input_path = input_path   # sighting matrix of the population
output_path = output_path # where you want to save the output

### load data
inputMat <- read.csv(input_path, header = TRUE)

### Preparing input data
# Getting start and end times from the input datasets - will depend on your input if this works, but can also manually assign the year the observations began and ended
Ss = colnames(inputMat[4]) # start of study
Se = colnames(inputMat[length(inputMat)-3]) # end of study (not the last row as the three last rows are male, female and unknown)

### Assigning either male or female to individuals of unknown sex
# Probability of unknowns being male, with assumption: the probability of unknowns being male is 1-adult male ratio
numUnknown <- sum(inputMat[ , length(inputMat)])
numMales <- sum(inputMat[ , length(inputMat)-1])
numFemales <- sum(inputMat[ , length(inputMat)-2])
howManyMoreMalesThanFemales <- numMales - numFemales
moreFemalesNeeded <- howManyMoreMalesThanFemales + (0.5 * (numUnknown - howManyMoreMalesThanFemales))
moreMalesNeeded <- numUnknown - moreFemalesNeeded
p_male <- moreMalesNeeded / numUnknown # probability of being male if assigned as unknown sex

### Function to create the input list to be able to run in parallel where sex was assigned for unknowns
MakeInput <- function(n,input,pm) {
  # Create a list of dataframes 
  InputList <- list()
  i <- 1
  
  while (i<=n) {
    
    print(i) 
    
    ## Assign either M or F if sex is U (remove column w U) ####
    male <- rbinom(n = nrow(input), size = 1, prob = pm)
    
    # is it a female? (then set female to 1 and male to 0)
    input[,length(input)-1][input[,length(input)] == 1 & male == 0] <- 0
    input[,length(input)-2][input[,length(input)] == 1 & male == 0] <- 1
    # is it a male?  (then set female to 0 and male to 1)
    input[,length(input)-1][input[,length(input)] == 1 & male == 1] <- 1
    input[,length(input)-2][input[,length(input)] == 1 & male == 1] <- 0
    
    ## Remove the "unknowns" columns in dataframe
    NewInput <- input[1:nrow(input),1:length(input)-1]
    InputList[[i]] <- NewInput
    
    i <- i+1
  }
  return(InputList)
}

### Values for the permutations and clusters
nperm = 10 # number of lists you want in each inputList
ncl = 10 # number of cpus used in cluster - check how many cores you have available/want to use

# run model ####
### Start the function in parallel -- unnecessary for me
cl <- makeCluster(ncl) # create the cluster element
registerDoParallel(cl)

InputList <- MakeInput(nperm,inputMat,p_male) # create input from function above

# Save input to be able to check it later
save(InputList, file=paste0(output_path,"input",".RData")) # Saved my input files to be able to go back and check

# Set the arguments for the iterations
iter = 500000 # number of iterations
sim = 4 # number of chains
warmup = iter/2 # number of iterations to use during warmup/burn-in period
thin = 101 # thinning number

# Set the prior values
pM = c(-3,0.2,0.001,-4,0.05) # prior means for parameters in a Gompertz model with a bathtub shape
pSD = c(1,0.05,0.001,1,0.01) # prior standard deviation for parameters in a Gompertz model with a bathtub shape

# begin the parallel loop
Perm.Out <- foreach(i=1:length(InputList), .packages='BaSTA',.verbose = T) %dopar% { # start the parallel processing - change 'BaSTA' with the package you will be using in the parallel loop
  
  # Run survival model
  Out <- basta(InputList[[i]],studyStart = Ss, # the model function will run through all lists in the input
               studyEnd = Se, nsim = sim, niter = iter, 
               burnin = warmup, thinning = thin, 
               thetaPriorMean=pM,thetaPriorSd=pSD,
               models = c("GO"), shape = c("bathtub"))
  
  return(Out)
}
# save the output
save(Perm.Out, file=paste0(output_path,"output.RData"))

stopCluster(cl)

#### Try with Amboseli data ####
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
iter <- 500      # number of iterations
sim <- 4         # number of chains
warmup <- iter/2 # number of iterations to use during warmup/burn-in period
thin <- 10       # thinning number

# Set the prior values
pM <-  c(-3,0.2,0.001,-4,0.05) # prior means for parameters in a Gompertz model with a bathtub shape
pSD <- c(1,0.05,0.001,1,0.01)  # prior standard deviation for parameters in a Gompertz model with a bathtub shape

basta_data <- males
basta_data$dyr <- ifelse(basta_data$censor == 'TRUE', 0, basta_data$dyr)
basta_data$casename <- paste0('M', basta_data$casename)
basta_data <- basta_data[-c(1:4,6,8,11,14:18,22,23,26:30,35,38,44:46,52,55,56,58,59,62,63,66,68,71,72,76,83,87,90,93,94,96,105,107,108,111,114,115,118,127,128,132,133,135,136,138,145,146,151,155,159,162,174,179,182,186,187,190,191,193,197,201,208,210,213,220,224,230,235,238,244,248,250,254,256,259,261,263,266,268,288,291,292,295,298,306,307,309,319,327,339,342,346,370,384,416,417,435,449,472,717,693,708,735,769,774),  # elephants that have sightings before birth/after death
                         c(1:3,6:55)]

#Go.Bt <- basta(basta_data, studyStart = Ss,
#               studyEnd = Se, nsim = sim, niter = iter,
#               burnin = warmup, thinning = thin,
#               thetaPriorMean = pM,thetaPriorSd = pSD,
#               models = c("GO"), shape = c("bathtub"))
#str(basta_data)
#summary(Go.Bt)
#plot(Go.Bt)
#plot(Go.Bt, plot.trace = FALSE)   # Plot survival and mortality curves

basta_data$age_2021 <- ifelse(basta_data$dyr > 100, basta_data$dyr - basta_data$byr, 2021 - basta_data$byr)
basta_data$censor_2021 <- ifelse(basta_data$dyr > 100, 'CENSORED', 'DEAD')
hist(basta_data$age_2021)
rethinking::dens(basta_data$age_2021, ylim = c(0, 0.07))
rethinking::dens(basta_data$age_2021[basta_data$censor_2021 == 'DEAD'], add = T, col = 'red')
rethinking::dens(basta_data$age_2021[basta_data$censor_2021 == 'CENSORED'], add = T, col = 'blue')

males$total_years_sighted <- rowSums(males[,6:55])
table(males$total_years_sighted)

basta_data$total_years_sighted <- rowSums(basta_data[,4:54])
table(basta_data$total_years_sighted)

basta_data <- basta_data[basta_data$casename %in% sort(unique(paste0('M',sightings$id))),]

sighted_males <- males[males$casename %in% sightings$id,]
unsighted_males <- anti_join(males, sighted_males)
head(unsighted_males)
summary(unsighted_males$age)


#### lets just start trying shit and see what happens! ####
# Set the arguments for the iterations
iter <- 500      # number of iterations
sim <- 4         # number of chains
warmup <- iter/2 # number of iterations to use during warmup/burn-in period
thin <- 11       # thinning number

# Gompertz Bathtub -- run using priors from killer whale paper ####
pM <-  c(-3,0.2,0.001,-4,0.05) # prior means for parameters in a Gompertz model with a bathtub shape
pSD <- c(1,0.05,0.001,1,0.01)  # prior standard deviation for parameters in a Gompertz model with a bathtub shape
Go.Bt <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               models = c("GO"), shape = c("bathtub"))
summary(Go.Bt)
plot(Go.Bt)                       # Plot traces
plot(Go.Bt, plot.trace = FALSE)   # Plot survival and mortality curves

# try the prior predictive plot for these parameters which we've agreed upon and see how they look
gompertz_bathtub <- function(b0, b1, a0, a1, c, age){
  gompertz <- exp(b0 + b1*age)
  makeham <- c
  bathtub <- exp(a0 - a1*age)
  y <- gompertz + makeham + bathtub
  print(y)
}

prior_test <- data.frame(b_1 = rnorm(6000, -3, 1),
                         b_0 = rnorm(6000, 0.2, 0.05),
                         a_1 = rnorm(6000, 0.001, 0.001),
                         a_0 = rnorm(6000, -4, 1),
                          c  = rnorm(6000, -0.05, 0.01),
                         age = rep(seq(10, 60, length.out = 6),1000),
                         iter= rep(1:1000, each = 6))
prior_test$y <- gompertz_bathtub(b0 = prior_test$b_0, b1 = prior_test$b_1,
                                 a0 = prior_test$a_0, a1 = prior_test$a_1,
                                 c = prior_test$c,   age = prior_test$age)
plot(prior_test$y,pch = 16, col = rgb(0,0,1,0.1))
plot(prior_test$y ~ prior_test$age, pch = 16, col = rgb(0,0,1,0.1))
plot(NULL, xlim = c(10,60), ylim = c(0,1), xlab = 'age', ylab = 'mortality')
for(i in 1:1000){
  plot_lines <- prior_test[prior_test$iter == i,]
  lines(plot_lines$y ~ plot_lines$age, col = rgb(0,0,1,0.1))
}
rethinking::dens(prior_test$y)

prior_test[23,]
b_0 <- prior_test[23,1]
b_1 <- prior_test[23,2]
a_0 <- prior_test[23,3]
a_1 <- prior_test[23,4]
c <- prior_test[23,5]
age <- prior_test[23,6]

y <- exp(b_0 + b_1*age) + c + exp(a_0 - a_1*age)

exp(250)
rethinking::dens(prior_test$y)


gompertz <- function(b0, b1, age){
  y <- exp(b0 + b1*age)
  print(y)
}
rethinking::dens(gompertz(1, 0.1, c(0:5)))

rethinking::dens(VGAM::rgompertz(1000000, scale = 1, shape = 0.1), ylim = c(0,1), xlim = c(0,5))


# Gompertz simple - Only requires 2 parameters -- run ####
rethinking::dens(VGAM::rgompertz(1000, scale = 1, shape = 1), ylim = c(0,5), xlim = c(0,2))
rethinking::dens(VGAM::rgompertz(1000, scale = 1, shape = 2), add = T, col = 'red')
rethinking::dens(VGAM::rgompertz(1000, scale = 1, shape = 3), add = T, col = 'orange')
rethinking::dens(VGAM::rgompertz(1000, scale = 1, shape = 4), add = T, col = 'yellow')
rethinking::dens(VGAM::rgompertz(1000, scale = 1, shape = 5), add = T, col = 'green')
rethinking::dens(VGAM::rgompertz(1000, scale = 2, shape = 1), add = T, col = 'blue')
rethinking::dens(VGAM::rgompertz(1000, scale = 3, shape = 1), add = T, col = 'purple')
rethinking::dens(VGAM::rgompertz(1000, scale = 4, shape = 1), add = T, col = 'magenta')
rethinking::dens(VGAM::rgompertz(1000, scale = 5, shape = 1), add = T, col = 'chocolate')
rethinking::dens(VGAM::rgompertz(1000, scale = 2, shape = 2), add = T, col = 'tan')
rethinking::dens(VGAM::rgompertz(1000, scale = 3, shape = 3), add = T, col = 'skyblue')
rethinking::dens(VGAM::rgompertz(1000, scale = 4, shape = 4), add = T, col = 'seagreen1')
rethinking::dens(VGAM::rgompertz(1000, scale = 5, shape = 5), add = T, col = 'darkgreen')
# Increasing shape makes it steeper and more right skewed. Increasing scale maintains peak position but reduces skewness.

rethinking::dens(VGAM::rgompertz(1000, scale = 1, shape = 1), xlim = c(0,5))
rethinking::dens(VGAM::rgompertz(1000, scale = 1, shape = 0.8), add = T, col = 'red')
rethinking::dens(VGAM::rgompertz(1000, scale = 1, shape = 0.5), add = T, col = 'orange')
rethinking::dens(VGAM::rgompertz(1000, scale = 1, shape = 0.2), add = T, col = 'yellow')
rethinking::dens(VGAM::rgompertz(1000, scale = 1, shape = 0.05), add = T, col = 'green')
rethinking::dens(VGAM::rgompertz(1000, scale = 0.8, shape = 1), add = T, col = 'blue')
rethinking::dens(VGAM::rgompertz(1000, scale = 0.5, shape = 1), add = T, col = 'purple')
rethinking::dens(VGAM::rgompertz(1000, scale = 0.2, shape = 1), add = T, col = 'magenta')
rethinking::dens(VGAM::rgompertz(1000, scale = 0.05, shape = 1), add = T, col = 'chocolate')
# Decreasing shape makes it broader and less skewed. Decreasing scale maintains peak position but increases tail
# want it to peak on the left but tail off slowly
rethinking::dens(VGAM::rgompertz(1000, scale = 0.3, shape = 0.6), add = T, col = 'seagreen1', lwd = 3)

pM <-  c(0.3,0.6)  # prior means for parameters in a Gompertz model with a simple shape
pSD <- c(0.1,0.2)  # prior standard deviation for parameters in a Gompertz model with a simple shape

rethinking::dens(VGAM::rgompertz(1000,
                                 scale = abs(rnorm(1,0.3,0.1)),
                                 shape = abs(rnorm(1,0.6,0.2))),
                 xlim = c(0,8), ylim = c(0,1))
for(i in 1:1000){
  rethinking::dens(VGAM::rgompertz(1000,
                                   scale = abs(rnorm(1,0.3,0.1)),
                                   shape = abs(rnorm(1,0.6,0.2))),
                   add = T, col = rethinking::col.alpha('blue',0.2))
}

Go.Si <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               models = c("GO"), shape = c("simple"))
summary(Go.Si)
plot(Go.Si)                       # Plot traces
plot(Go.Si, plot.trace = FALSE)   # Plot survival and mortality curves

# WEIBULL BATHTUB -- not yet worked out priors ####

prior_test <- data.frame(b_0 = rnorm(6000, -3, 1),
                         b_1 = rnorm(6000, 0.2, 0.05),
                         a_0 = rnorm(6000, 0.001, 0.001),
                         a_1 = rnorm(6000, -4, 1),
                          c  = rnorm(6000, -0.05, 0.01),
                         age = rep(seq(10, 60, length.out = 6),1000),
                         iter= rep(1:1000, each = 6))
weibull_bathtub <- function(b_0, b_1, a_0, a_1, c, age){
  f1 <- b_0 * b_1^(b_0) * age^(b_0-1)
  f2 <- exp(a_0 + a_1*age)
  y <- f1 + f2 + c
  print(y)
}
prior_test$y <- weibull_bathtub(prior_test$b_0, prior_test$b_1, prior_test$a_0, prior_test$a_1, prior_test$c, prior_test$age)
plot(prior_test$y,pch = 16, col = rgb(0,0,1,0.1))
plot(prior_test$y ~ prior_test$age, pch = 16, col = rgb(0,0,1,0.1))
plot(NULL, xlim = c(10,60), ylim = c(-4,1), xlab = 'age', ylab = 'mortality')
for(i in 1:1000){
  plot_lines <- prior_test[prior_test$iter == i,]
  lines(plot_lines$y ~ plot_lines$age, col = rgb(0,0,1,0.1))
}
rethinking::dens(prior_test$y)


prior_test <- data.frame(b_0 = rnorm(6000, 0, 5),
                         b_1 = rnorm(6000, 0, 1),
                         a_0 = rnorm(6000, 0, 1),
                         a_1 = rnorm(6000, 0, 5),
                          c  = rnorm(6000, 0, 1),
                         age = rep(seq(10, 60, length.out = 6),1000),
                         iter= rep(1:1000, each = 6))
weibull_bathtub <- function(b_0, b_1, a_0, a_1, c, age){
  f1 <- b_0 * b_1^(b_0) * age^(b_0-1)
  f2 <- exp(a_0 + a_1*age)
  y <- f1 + f2 + c
  print(y)
}
prior_test$y <- weibull_bathtub(prior_test$b_0, prior_test$b_1, prior_test$a_0, prior_test$a_1, prior_test$c, prior_test$age)
plot(prior_test$y,pch = 16, col = rgb(0,0,1,0.1))
plot(prior_test$y ~ prior_test$age, pch = 16, col = rgb(0,0,1,0.1))
plot(NULL, xlim = c(10,60), ylim = c(-50000000000,50000000000), xlab = 'age', ylab = 'mortality')
for(i in 1:1000){
  plot_lines <- prior_test[prior_test$iter == i,]
  lines(plot_lines$y ~ plot_lines$age, col = rgb(0,0,1,0.1))
}
rethinking::dens(prior_test$y)

prior_test <- data.frame(b_0 = rnorm(6000, -3, 1),
                         b_1 = rnorm(6000, 0.2, 0.05),
                         a_0 = rnorm(6000, 0.001, 0.001),
                         a_1 = rnorm(6000, -4, 1),
                         c  = rnorm(6000, -0.05, 0.01),
                         age = rep(seq(10, 60, length.out = 6),1000),
                         iter= rep(1:1000, each = 6))
weibull_bathtub <- function(b_0, b_1, a_0, a_1, c, age){
  f1 <- b_0 * b_1^(b_0) * age^(b_0-1)
  f2 <- exp(a_0 + a_1*age)
  y <- f1 + f2 + c
  print(y)
}
prior_test$y <- weibull_bathtub(prior_test$b_0, prior_test$b_1, prior_test$a_0, prior_test$a_1, prior_test$c, prior_test$age)
plot(prior_test$y,pch = 16, col = rgb(0,0,1,0.1))
plot(prior_test$y ~ prior_test$age, pch = 16, col = rgb(0,0,1,0.1))
plot(NULL, xlim = c(10,60), ylim = c(-4,1), xlab = 'age', ylab = 'mortality')
for(i in 1:1000){
  plot_lines <- prior_test[prior_test$iter == i,]
  lines(plot_lines$y ~ plot_lines$age, col = rgb(0,0,1,0.1))
}
rethinking::dens(prior_test$y)

prior_test <- data.frame(b_0 = rnorm(6000, 3, 1),
                         b_1 = rnorm(6000, 3, 2),
                         a_0 = rnorm(6000, 0, 2),
                         a_1 = rnorm(6000, 4, 1),
                          c  = rnorm(6000, 0, 0.05),
                         age = rep(seq(10, 60, length.out = 6),1000),
                         iter= rep(1:1000, each = 6))
weibull_bathtub <- function(b_0, b_1, a_0, a_1, c, age){
  f1 <- b_0 * b_1^(b_0) * age^(b_0-1)
  f2 <- exp(a_0 + a_1*age)
  y <- f1 + f2 + c
  print(y)
}
prior_test$y <- weibull_bathtub(prior_test$b_0, prior_test$b_1, prior_test$a_0, prior_test$a_1, prior_test$c, prior_test$age)
plot(prior_test$y,pch = 16, col = rgb(0,0,1,0.1))
plot(prior_test$y ~ prior_test$age, pch = 16, col = rgb(0,0,1,0.1))
plot(NULL, xlim = c(10,60), ylim = c(-4,1), xlab = 'age', ylab = 'mortality')
for(i in 1:1000){
  plot_lines <- prior_test[prior_test$iter == i,]
  lines(plot_lines$y ~ plot_lines$age, col = rgb(0,0,1,0.1))
}
rethinking::dens(prior_test$y)























pM <-  c(5,5,5,-1,5) # prior means for parameters in a Weibull model with a bathtub shape
pSD <- c(0.1,0.1,0.1,0.1,0.1)  # prior standard deviation for parameters in a Weibull model with a bathtub shape

Wb.Bt <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               models = c("WE"), shape = c("bathtub"))
summary(Wb.Bt)
plot(Wb.Bt)                       # Plot traces
plot(Wb.Bt, plot.trace = FALSE)   # Plot survival and mortality curves

# Weibull simple -- run ####
rethinking::dens(rweibull(1000, scale = 1, shape = 1), ylim = c(0,2), xlim = c(0,8))
rethinking::dens(rweibull(1000, scale = 1, shape = 2), add = T, col = 'red')
rethinking::dens(rweibull(1000, scale = 1, shape = 3), add = T, col = 'orange')
rethinking::dens(rweibull(1000, scale = 1, shape = 4), add = T, col = 'yellow')
rethinking::dens(rweibull(1000, scale = 1, shape = 5), add = T, col = 'green')
rethinking::dens(rweibull(1000, scale = 1, shape = 0.5), add = T, col = 'blue')
rethinking::dens(rweibull(1000, scale = 1, shape = 0.1), add = T, col = 'purple')
rethinking::dens(rweibull(1000, scale = 1, shape = 0.005), add = T, col = 'magenta')
# Increasing shape makes it narrower and less right skewed. sub 1 quickly moves to crazy narrow distributions of tiny values
rethinking::dens(rweibull(1000, scale = 1, shape = 1), ylim = c(0,2), xlim = c(0,20))
rethinking::dens(rweibull(1000, scale = 2, shape = 1), add = T, col = 'red')
rethinking::dens(rweibull(1000, scale = 3, shape = 1), add = T, col = 'orange')
rethinking::dens(rweibull(1000, scale = 4, shape = 1), add = T, col = 'yellow')
rethinking::dens(rweibull(1000, scale = 5, shape = 1), add = T, col = 'green')
rethinking::dens(rweibull(1000, scale = 0.5, shape = 1), add = T, col = 'blue')
rethinking::dens(rweibull(1000, scale = 0.1, shape = 1), add = T, col = 'purple')
rethinking::dens(rweibull(1000, scale = 0.05, shape = 1), add = T, col = 'magenta')
# Increasing scale makes it wider and more right skewed. sub 1 quickly moves to crazy narrow distributions of tiny values

rethinking::dens(rweibull(1000, scale = 1, shape = 1), ylim = c(0,2), xlim = c(0,15))
rethinking::dens(rweibull(1000, scale = 2, shape = 2), add = T, col = 'red')
rethinking::dens(rweibull(1000, scale = 3, shape = 3), add = T, col = 'orange')
rethinking::dens(rweibull(1000, scale = 4, shape = 4), add = T, col = 'yellow')
rethinking::dens(rweibull(1000, scale = 5, shape = 5), add = T, col = 'green')
rethinking::dens(rweibull(1000, scale = 0.5, shape = 0.5), add = T, col = 'blue')
rethinking::dens(rweibull(1000, scale = 0.1, shape = 0.1), add = T, col = 'purple')
rethinking::dens(rweibull(1000, scale = 0.05, shape = 0.05), add = T, col = 'magenta')
rethinking::dens(rweibull(1000, scale = 10, shape = 10), add = T, col = 'chocolate')
# Increasing together moves probability density to the right, and makes whole thing more left skewed

rethinking::dens(rweibull(1000, scale = 3, shape = 2), add = T, col = 'seagreen1', lwd = 3)


pM <-  c(3,2) # prior means for parameters in a Weibull model with a simple shape
pSD <- c(0.5,0.5)  # prior standard deviation for parameters in a Weibull model with a simple shape
Wb.Si <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               models = c("WE"), shape = c("simple"))
summary(Wb.Si)
plot(Wb.Si)                       # Plot traces
plot(Wb.Si, plot.trace = FALSE)   # Plot survival and mortality curves


# prior predictive check attempt?
data <- data.frame(scale = rnorm(6100, 3, 0.5),
                   shape = rnorm(6100, 2, 0.5),
                   age   = rep(0:60,100),
                   iteration = rep(1:100, each = 61))
weibull <- function(data){
  k <- data[,1]
  lambda <- data[,2]
  x <- data[,3]
  a <- k/lambda
  b <- (x/lambda)^(k-1)
  c <- -1*(x/lambda)^(k)
  y <- a*b*exp(c)
  print(y)
}
data$y <- weibull(data = data)
head(data)

plot(data$y ~ data$age)
rethinking::dens(data$y)

data[1:62,]
plot(NULL, xlim = c(0,60), ylim = c(0,1), xlab = 'age', ylab = 'proportion surviving')
for(i in 1:100){
  iteration <- data[data$iteration == i,]
  lines(iteration$y ~ iteration$age, col = rgb(0,0,1,0.2))
}

set.seed(1)
data <- data.frame(scale = rnorm(6100, 3, 2),
                   shape = rnorm(6100, 2, 2),
                   age   = rep(0:60,100),
                   iteration = rep(1:100, each = 61))
data$y <- weibull(data = data)
head(data, 20)
data$y_0 <- ifelse(is.nan(data$y) == TRUE, 0, data$y) # should only get NaN when lambda = 0, but this has apparently happened way too many times for that!
# scale = 2.389223, shape = -3.0614537, age = 9 --> NaN
scale = 2.389223; shape =  -3.0614537; age = 9
a <- shape/scale              # -1.281
b <- (age/scale)^(shape - 1)  #  0.005
c <- -1*(age/scale)^shape     # -0.017
y <- a*b*exp(c)               # -0.006 --> so it IS a number... so why did it come out as NaN?

plot(NULL, xlim = c(0,60), ylim = c(0,1), xlab = 'age', ylab = 'proportion surviving')
for(i in 1:100){
  iteration <- data[data$iteration == i,]
  lines(iteration$y_0 ~ iteration$age, col = rgb(0,0,1,0.2))
}

set.seed(1)
data <- data.frame(scale = rnorm(6100, 15, 2),
                   shape = rnorm(6100, 10, 2),
                   age   = rep(0:60,100),
                   iteration = rep(1:100, each = 61))
data$y <- weibull(data = data)
head(data, 20)
data$y_0 <- ifelse(is.nan(data$y) == TRUE, 0, data$y) # should only get NaN when lambda = 0, but this has apparently happened way too many times for that!
plot(y_0 ~ age, data = data[data$iteration == 1,], lwd = 2, type = 'l',
     xlim = c(0,60), ylim = c(0,1), xlab = 'age', ylab = 'proportion surviving', las = 1)
for(i in 2:100){
  iteration <- data[data$iteration == i,]
  lines(iteration$y_0 ~ iteration$age, col = rgb(0,0,1,0.2))
}

set.seed(1)
data <- data.frame(scale = rnorm(6100, 10, 2),
                   shape = rnorm(6100, 40, 2),
                   age   = rep(0:60,100),
                   iteration = rep(1:100, each = 61))
data$y <- weibull(data = data)
head(data, 20)
data$y_0 <- ifelse(is.nan(data$y) == TRUE, 0, data$y) # should only get NaN when lambda = 0, but this has apparently happened way too many times for that!
plot(y_0 ~ age, data = data[data$iteration == 1,], lwd = 2, type = 'l',
     xlim = c(0,60), ylim = c(0,1), xlab = 'age', ylab = 'proportion surviving', las = 1)
for(i in 2:100){
  iteration <- data[data$iteration == i,]
  lines(iteration$y_0 ~ iteration$age, col = rgb(0,0,1,0.2))
}


set.seed(1)
data <- data.frame(scale = rnorm(6100, 1.5, 0.5),
                   shape = rnorm(6100, 30, 8),
                   age   = rep(0:60,100),
                   iteration = rep(1:100, each = 61))
data$y <- weibull(data = data)
head(data, 20)
data$y_0 <- ifelse(is.nan(data$y) == TRUE, 0, data$y) # should only get NaN when lambda = 0, but this has apparently happened way too many times for that!
plot(y_0 ~ age, data = data[data$iteration == 1,], lwd = 2, type = 'l',
     xlim = c(0,60), ylim = c(0,0.2), xlab = 'age', ylab = 'proportion surviving', las = 1)
for(i in 2:100){
  iteration <- data[data$iteration == i,]
  lines(iteration$y_0 ~ iteration$age, col = rgb(0,0,1,0.2))
}
rethinking::dens(rweibull(1000, scale = 1.5, shape = 30)) # this doesn't look like how I would expect it to be, but maybe I've got the wrong idea in my head?

iter <- 500 ; warmup <- iter/2 
pM <-  c(1.5,30) # prior means for parameters in a Weibull model with a simple shape
pSD <- c(0.5,8)  # prior standard deviation for parameters in a Weibull model with a simple shape
Wb.Si <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               models = c("WE"), shape = c("simple"))
summary(Wb.Si)
plot(Wb.Si)                       # Plot traces
plot(Wb.Si, plot.trace = FALSE)   # Plot survival and mortality curves -- looks almost exactly the same as the one from the other prior I came up with

##### Dan debugging #####
males <- readxl::read_excel('data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% 
  select(CASENAME, BYR, DYR) %>% 
  janitor::clean_names()
males$age <- as.numeric(ifelse(males$dyr < 10, 2021 - males$byr, males$dyr - males$byr))
table(males$age)
hist(males$age, breaks = 75)    # overall age distribution
males$censor <- ifelse(males$dyr < 10, 'TRUE', 'FALSE')    # censor = TRUE when elephant is still alive
table(males$censor)
males$X1972 <- NA ; males$X1973 <- NA ; males$X1974 <- NA ; males$X1975 <- NA
males$X1976 <- NA ; males$X1977 <- NA ; males$X1978 <- NA ; males$X1979 <- NA
males$X1980 <- NA ; males$X1981 <- NA ; males$X1982 <- NA ; males$X1983 <- NA
males$X1984 <- NA ; males$X1985 <- NA ; males$X1986 <- NA ; males$X1987 <- NA
males$X1988 <- NA ; males$X1989 <- NA ; males$X1990 <- NA ; males$X1991 <- NA
males$X1992 <- NA ; males$X1993 <- NA ; males$X1994 <- NA ; males$X1995 <- NA
males$X1996 <- NA ; males$X1997 <- NA ; males$X1998 <- NA ; males$X1999 <- NA
males$X2000 <- NA ; males$X2001 <- NA ; males$X2002 <- NA ; males$X2003 <- NA
males$X2004 <- NA ; males$X2005 <- NA ; males$X2006 <- NA ; males$X2007 <- NA
males$X2008 <- NA ; males$X2009 <- NA ; males$X2010 <- NA ; males$X2011 <- NA
males$X2012 <- NA ; males$X2013 <- NA ; males$X2014 <- NA ; males$X2015 <- NA
males$X2016 <- NA ; males$X2017 <- NA ; males$X2018 <- NA ; males$X2019 <- NA
males$X2020 <- NA ; males$X2021 <- NA

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

basta_data <- males
basta_data$dyr <- ifelse(basta_data$censor == 'TRUE', 0, basta_data$dyr)
basta_data$casename <- paste0('M', basta_data$casename)
basta_data <- basta_data[-c(1:4,6,8,11,14:18,22,23,26:30,35,38,44:46,52,55,56,58,59,62,63,66,68,71,72,76,83,87,90,93,94,96,105,107,108,111,114,115,118,127,128,132,133,135,136,138,145,146,151,155,159,162,174,179,182,186,187,190,191,193,197,201,208,210,213,220,224,230,235,238,244,248,250,254,256,259,261,263,266,268,288,291,292,295,298,306,307,309,319,327,339,342,346,370,384,416,417,435,449,472,717,693,708,735,769,774),  # elephants that have sightings before birth/after death
                         c(1:3,6:55)]
basta_data <- basta_data[basta_data$casename %in% paste0('M',sightings$id),]

basta_data$age_2021 <- ifelse(basta_data$dyr > 100, basta_data$dyr - basta_data$byr, 2021 - basta_data$byr)
basta_data$censor_2021 <- ifelse(basta_data$dyr > 100, 'CENSORED', 'DEAD')
hist(basta_data$age_2021)
rethinking::dens(basta_data$age_2021, ylim = c(0, 0.07))
rethinking::dens(basta_data$age_2021[basta_data$censor_2021 == 'DEAD'], add = T, col = 'red')
rethinking::dens(basta_data$age_2021[basta_data$censor_2021 == 'CENSORED'], add = T, col = 'blue')

males$total_years_sighted <- rowSums(males[,6:55])
table(males$total_years_sighted)

basta_data$total_years_sighted <- rowSums(basta_data[,4:53])
table(basta_data$total_years_sighted)
basta_data[basta_data$total_years_sighted == 42,]

sighted_males <- males[males$casename %in% sightings$id,]
unsighted_males <- anti_join(males, sighted_males)
head(unsighted_males)
summary(unsighted_males$age)

source('basta_source_codes/basta.R')
source('basta_source_codes/basta.default.R')
source('basta_source_codes/CensusToCaptHist.R')
source('basta_source_codes/DataCheck.R')
source('basta_source_codes/MakeCovMat.R')
source('basta_source_codes/MakeLifeTable.R')
source('basta_source_codes/multibasta.R')
source('basta_source_codes/plot.basta.R')
source('basta_source_codes/print.basta.R')
source('basta_source_codes/summary.basta.R')

debugonce(basta)
library(snowfall)

Ss <- colnames(basta_data[,4])
Se <- colnames(basta_data[,53])

Go.Bt <- basta(basta_data,
               studyStart = 1972, studyEnd = 2021,
               nsim = 4, niter = 500, burnin = 250, thinning = 10,
               thetaPriorMean = c(-3,0.2,0.001,-4,0.05),
               thetaPriorSd = c(1,0.05,0.001,1,0.01),
               models = c("GO"), shape = c("bathtub"), parallel = T)













