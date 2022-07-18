#### Information #####
# To estimate the ages of individuals in other populations, first come up with a survival curve for the Amboseli elephants which have good birth and death data.
# survival curve options = exponential/logistic/Gompertz/Weibull, where shape = simple/bathtub/Makeham
# follow through Mia's cde to determine which curve will best fit Amboseli data. she ran it on a high performance cluster, so it’s set up to run in parallel, but the basta model can also be run in parallel on its own (just like Stan). we will need to implement it as a custom function in Stan for the final analysis, given that we need to define it in the prior anyway. 

#### Load packages ####
lapply(c("foreach", "doParallel"), require, character.only = TRUE) # required if you want to run the model using parallel computing
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
### Start the function in parallel
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
?basta # object = A data.frame to be used as an input data file for BaSTA. The first column is a vector of individual unique IDs, the second and third columns are birth and death years respectively. Columns 4, …, n_t-1 represent the observation window of n_t years. This is followed (optionally) by columns for categorical and continuous covariates (see details).
# Example:
data("sim1", package = "BaSTA")
new.dat  <- DataCheck(sim1, studyStart = 51, studyEnd = 70, autofix = rep(1,7)) # Check data consistency
Out <- basta(sim1, studyStart = 51, studyEnd = 70, niter = 200, 
             burnin = 11, thinning = 10, 
             thetaPriorMean = pM,thetaPriorSd = pSD, nsim = 4,  # thetaPriorMean = mean prior, thetaPriorSd = SD prior, nsim = number chains
             models = c("GO"), shape = c("bathtub"))
summary(Out, digits = 3) # Print results
plot(Out)                # Plot traces for survival parameters
plot(Out, trace.name = "gamma") # Plot traces for proportional hazards parameter
plot(Out, plot.trace = FALSE)   # Plot survival and mortality curves

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

### Preparing input data
# Getting start and end times from the input datasets - will depend on your input if this works, but can also manually assign the year the observations began and ended
Ss <- 1972  # start of study
Se <- 2021  # end of study

# Set the arguments for the iterations
iter <- 5000   # number of iterations
sim <- 1         # number of chains
warmup <- iter/2 # number of iterations to use during warmup/burn-in period
thin <- 101      # thinning number

# Set the prior values
pM <-  c(-3,0.2,0.001,-4,0.05) # prior means for parameters in a Gompertz model with a bathtub shape
pSD <- c(1,0.05,0.001,1,0.01)  # prior standard deviation for parameters in a Gompertz model with a bathtub shape

basta_data <- males
basta_data$dyr <- ifelse(basta_data$censor == 'TRUE', 0, basta_data$dyr)
basta_data <- basta_data[-c(1:4,6,8,11,14:18,22,23,26:30,35,38,44:46,52,55,56,58,59,62,63,66,68,71,72,76,83,87,90,93,94,96,105,107,108,111,114,115,118,127,128,132,133,135,136,138,145,146,151,155,159,162,174,179,182,186,187,190,191,193,197,201,208,210,213,220,224,230,235,238,244,248,250,254,256,259,261,263,266,268,288,291,292,295,298,306,307,309,319,327,339,342,346,370,384,416,417,435,449,472,717,693,708,735,769,774),
                         c(1:3,6:55)]
Out <- basta(basta_data, studyStart = Ss,
             studyEnd = Se, nsim = sim, niter = iter, 
             burnin = warmup, thinning = thin, 
             thetaPriorMean = pM,thetaPriorSd = pSD,
             models = c("GO"), shape = c("bathtub"))
summary(Out)
plot(Out)










### Assigning either male or female to individuals of unknown sex -- THIS DOESN'T WORK WHEN THE ONLY SIGHTINGS DATA I HAVE ARE FOR MALES!
numUnknown <- sum(males[ , ncol(males)])
numMales <- sum(males[ , length(males)-1])
numFemales <- sum(males[ , length(males)-2])
howManyMoreFemalesThanMales <- numFemales - numMales
moreMalesNeeded <- howManyMoreFemalesThanMales + (0.5 * (numUnknown - howManyMoreFemalesThanMales))
moreFemalesNeeded <- numUnknown - moreMalesNeeded
p_female <- moreFemalesNeeded / numUnknown # probability of being male if assigned as unknown sex

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
nperm  <- 10 # number of lists you want in each inputList
ncl <- 4     # number of cpus used in cluster - check how many cores you have available/want to use

# run model ####
### Start the function in parallel
cl <- makeCluster(ncl) # create the cluster element
doParallel::registerDoParallel(cl)

InputList <- MakeInput(nperm,males,p_female) # create input from function above
table(InputList[[10]]$male) ; table(InputList[[10]]$female)


# Set the arguments for the iterations
iter <- 500000   # number of iterations
sim <- 4         # number of chains
warmup <- iter/2 # number of iterations to use during warmup/burn-in period
thin <- 101      # thinning number

# Set the prior values
pM <-  c(-3,0.2,0.001,-4,0.05) # prior means for parameters in a Gompertz model with a bathtub shape
pSD <- c(1,0.05,0.001,1,0.01)  # prior standard deviation for parameters in a Gompertz model with a bathtub shape

# begin the parallel loop
Perm.Out <- foreach(i = 1:length(InputList), .packages='BaSTA', .verbose = T) %dopar% { # start the parallel processing - change 'BaSTA' with the package you will be using in the parallel loop
  
  # Run survival model
  Out <- basta(InputList[[i]],studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter, 
               burnin = warmup, thinning = thin, 
               thetaPriorMean = pM,thetaPriorSd = pSD,
               models = c("GO"), shape = c("bathtub"))
  
  return(Out)
}
# save the output -- ???

stopCluster(cl)
