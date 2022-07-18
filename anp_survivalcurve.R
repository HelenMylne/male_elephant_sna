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
## read in AEDD data ####
all_nodes <- read_csv('data_raw/AfricanElephantDemographicDatabase_aedd_22.04.04.csv')
all_nodes <- all_nodes[all_nodes$species == 'L. africana' & !is.na(all_nodes$species) &
                         all_nodes$metric == 'population size' & !is.na(all_nodes$age_reported),
                       c(1,5,7,8,13:17,22)]

## population size: estimate = number observed to reach minimum age (censor at min or median?)
sort(unique(all_nodes$age_reported))
all_nodes <- all_nodes[all_nodes$age_reported != 'adult' & all_nodes$age_reported != 'calf' & 
                         all_nodes$age_reported != 'fullgrown' & all_nodes$age_reported != 'halfgrown' & 
                         all_nodes$age_reported != 'juvenile' & all_nodes$age_reported != 'sexually_immature' & 
                         all_nodes$age_reported != 'sexually_mature',]
all_nodes$min <- with(all_nodes,
                      case_when(age_reported == '<1'   ~ min(x = c(0.1,1)),
                                age_reported == '<10'     ~ min(x = c(0.1,10)),
                                age_reported == '>30'     ~ min(x = c(30,65)),
                                age_reported == '>40'     ~ min(x = c(40,65)),
                                age_reported == '>50'     ~ min(x = c(50,65)),
                                age_reported == '0-1'     ~ min(x = c(0.1,1)),
                                age_reported == '0-4.9'   ~ min(x = c(0.1,4.9)),
                                age_reported == '0-5'     ~ min(x = c(0.1,5)),
                                age_reported == '1-11'    ~ min(x = c(1,11)),
                                age_reported == '1-12'    ~ min(x = c(1,12)),
                                age_reported == '10-14.9' ~ min(x = c(10,14.9)),
                                age_reported == '10-15'   ~ min(x = c(10,15)),
                                age_reported == '10-20'   ~ min(x = c(10,20)),
                                age_reported == '15-19.9' ~ min(x = c(15,19.9)),
                                age_reported == '15-20'   ~ min(x = c(15,20)),
                                age_reported == '16-30'   ~ min(x = c(16,30)),
                                age_reported == '20-24.9' ~ min(x = c(20,24.9)),
                                age_reported == '20-25'   ~ min(x = c(20,25)),
                                age_reported == '20-34.9' ~ min(x = c(20,34.9)),
                                age_reported == '20-42'   ~ min(x = c(20,42)),
                                age_reported == '25-29.9' ~ min(x = c(10,20)),
                                age_reported == '25-30'   ~ min(x = c(25,30)),
                                age_reported == '25-34.9' ~ min(x = c(25,34.9)),
                                age_reported == '30-34.9' ~ min(x = c(30,34.9)),
                                age_reported == '30-35'   ~ min(x = c(30,35)),
                                age_reported == '35-39.9' ~ min(x = c(35,39.9)),
                                age_reported == '35-40'   ~ min(x = c(35,40)),
                                age_reported == '35-49.9' ~ min(x = c(35,49.9)),
                                age_reported == '4-10'    ~ min(x = c(4,10)),
                                age_reported == '40-45'   ~ min(x = c(40,45)),
                                age_reported == '45-50'   ~ min(x = c(45,50)),
                                age_reported == '5-10'    ~ min(x = c(5,10)),
                                age_reported == '5-9.9'   ~ min(x = c(5,9.9)),
                                age_reported == '50-55'   ~ min(x = c(50,55)),
                                age_reported == '55-60'   ~ min(x = c(55,60)),
                                age_reported == '6-15'    ~ min(x = c(6,15)),
                                age_reported == '60-65'   ~ min(x = c(60,65))))
all_nodes$max <- with(all_nodes,
                      case_when(age_reported == '<1'      ~ max(x = c(0.1,1)),
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
unique(all_nodes$age_reported)
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

all_nodes$sex <- ifelse(all_nodes$sex == 'male', 'M', 'F')
all_nodes$age <- with(all_nodes,
                      case_when(age_reported == '<1'      ~ '0-5',
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
                                age_reported == '10-20'   ~ 'NA', # can't differentiate pubescent age groups
                                age_reported == '<10'     ~ 'NA', # can't differentiate calf-juvenile age groups
                                age_reported == '1-11'    ~ 'NA', # can't differentiate calf-juvenile age groups
                                age_reported == '1-12'    ~ 'NA', # can't differentiate calf-juvenile age groups
                                age_reported == '20-42'   ~ 'NA', # can't differentiate adult age groups
                                age_reported == '16-30'   ~ 'NA', # can't differentiate adult-pubescent age groups
                                age_reported == '6-15'    ~ 'NA', # can't differentiate juvenile-pubescent age groups
                                age_reported == '>30'     ~ 'NA', # can't differentiate adult age groups
                                TRUE ~ as.character(age_reported)))
unique(all_nodes$age)
all_nodes <- all_nodes[all_nodes$age != 'NA',]
all_nodes$year_of_investigation <- ifelse(all_nodes$year == '2001-2002', 2001,
                                          ifelse(all_nodes$year == '1997-1999', 1998,
                                                 ifelse(all_nodes$year == '1968-1969', 1968,
                                                        ifelse(all_nodes$year == '1998-2003', 2000,
                                                               as.numeric(all_nodes$year)))))
all_nodes$byr_earliest <- all_nodes$year_of_investigation - round(all_nodes$max,0)
all_nodes$byr_latest   <- all_nodes$year_of_investigation - round(all_nodes$min,0)
all_nodes

all_nodes_cens <- data.frame(age = rep(NA, nrow(all_nodes)),
                             max = rep(NA, nrow(all_nodes)),
                             min = rep(NA, nrow(all_nodes)),
                             bye = rep(NA, nrow(all_nodes)),
                             byl = rep(NA, nrow(all_nodes)),
                             sex = rep(NA, nrow(all_nodes)),
                             sid = rep(NA, nrow(all_nodes)),
                             loc = rep(NA, nrow(all_nodes)),
                             syr = rep(NA, nrow(all_nodes)))
for(i in 1:nrow(all_nodes)){
  all_nodes_cens$age[i] <- list(rep(all_nodes$age[i], ifelse(all_nodes$estimate[i] < 1, 0, all_nodes$estimate[i])))
  all_nodes_cens$max[i] <- list(rep(all_nodes$max[i], ifelse(all_nodes$estimate[i] < 1, 0, all_nodes$estimate[i])))
  all_nodes_cens$min[i] <- list(rep(all_nodes$min[i], ifelse(all_nodes$estimate[i] < 1, 0, all_nodes$estimate[i])))
  all_nodes_cens$bye[i] <- list(rep(all_nodes$byr_earliest[i],
                                    ifelse(all_nodes$estimate[i] < 1, 0, all_nodes$estimate[i])))
  all_nodes_cens$byl[i] <- list(rep(all_nodes$byr_latest[i],
                                    ifelse(all_nodes$estimate[i] < 1, 0, all_nodes$estimate[i])))
  all_nodes_cens$sex[i] <- list(rep(all_nodes$sex[i], ifelse(all_nodes$estimate[i] < 1, 0, all_nodes$estimate[i])))
  all_nodes_cens$sid[i] <- list(rep(all_nodes$study_id[i], ifelse(all_nodes$estimate[i] < 1, 0, all_nodes$estimate[i])))
  all_nodes_cens$loc[i] <- list(rep(all_nodes$sitename[i], ifelse(all_nodes$estimate[i] < 1, 0, all_nodes$estimate[i])))
  all_nodes_cens$syr[i] <- list(rep(all_nodes$`time-period_studied`[i], ifelse(all_nodes$estimate[i] < 1, 0, all_nodes$estimate[i])))
}
nodes_cens <- data.frame(age = unlist(all_nodes_cens$age),
                         alive = TRUE,
                         cens = 1,
                         max_age = unlist(all_nodes_cens$max),
                         min_age = unlist(all_nodes_cens$min),
                         byr_earliest = unlist(all_nodes_cens$bye),
                         byr_latest = unlist(all_nodes_cens$byl),
                         sex = unlist(all_nodes_cens$sex),
                         study_id = unlist(all_nodes_cens$sid),
                         location = unlist(all_nodes_cens$loc),
                         study_year = unlist(all_nodes_cens$syr))

nodes_cens$k <- ifelse(nodes_cens$age == '0-5', 1,
                       ifelse(nodes_cens$age == '5-10', 2,
                              ifelse(nodes_cens$age == '10-15', 3,
                                     ifelse(nodes_cens$age == '15-20', 4,
                                            ifelse(nodes_cens$age == '20-25' |
                                                     nodes_cens$age == '20-35', 5,
                                                   ifelse(nodes_cens$age == '25-40' | 
                                                            nodes_cens$age == '35-50', 6, 7))))))
unique(nodes_cens$k)
#rm(all_nodes_cens)

## create data frame for basta ####
?basta # object = A data.frame to be used as an input data file for BaSTA. The first column is a vector of individual unique IDs, the second and third columns are birth and death years respectively. Columns 4, …, n_t-1 represent the observation window of n_t years. This is followed (optionally) by columns for categorical and continuous covariates (see details).
# Example:
data("sim1", package = "BaSTA")
new.dat  <- DataCheck(sim1, studyStart = 51, studyEnd = 70, autofix = rep(1,7)) # Check data consistency
out <- basta(sim1, studyStart = 51, studyEnd = 70, niter = 200, burnin = 11, 
             thinning = 10, updateJumps = FALSE) # Run short version of BaSTA on the data
summary(out, digits = 3) # Print results
plot(out)                # Plot traces for survival parameters
plot(out, trace.name = "gamma") # Plot traces for proportional hazards parameter
plot(out, plot.trace = FALSE)   # Plot survival and mortality curves








nodes_cens$id <- 1:nrow(nodes_cens)
nodes_cens$byr_sample <- NA # not sure if this gives a better estimation of the true ages, or if this makes the entire thing circular because I'm assuming a uniform distribution within age category and then using that to predict true age within category?
set.seed(1)
for(i in 1:nrow(nodes_cens)) {
  possible_years <- seq(from = nodes_cens$byr_earliest[i], to = nodes_cens$byr_latest[i], by = 1)
  nodes_cens$byr_sample[i] <- sample(x = possible_years, size = 1)
}

basta_df <- data.frame(id = nodes_cens$id,
                       byr_sample = nodes_cens$byr_sample,
                       dyr = NA,
                       obs_window = nodes_cens$study_year,
                       sex = nodes_cens$sex)
unique(basta_df$obs_window) # "1995-2002" "1964" "1977-1989" "1997-1999" "1998-2003" "1931-1998" "1978-1997" "1996" "1968-1969"
basta_df <- separate(basta_df, obs_window, sep = '-', into = c('studyStart','studyEnd'), remove = FALSE)
#test <- basta_df[,4:6] %>% distinct()
#rm(test)
basta_df$studyEnd[is.na(basta_df$studyEnd)] <- basta_df$studyStart[is.na(basta_df$studyEnd)]





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
ncl = 10   # number of cpus used in cluster - check how many cores you have available/want to use

# run model ####
### Start the function in parallel
cl <- makeCluster(ncl) # create the cluster element
registerDoParallel(cl)

InputList <- MakeInput(nperm,inputMat,p_male) # create input from function above

# Save input to be able to check it later
save(InputList, file=paste0(output_path,"input",".RData")) # Saved my input files to be able to go back and check

# Set the arguments for the iterations
iter = 500000   # number of iterations
sim = 4         # number of chains
warmup = iter/2 # number of iterations to use during warmup/burn-in period
thin = 101      # thinning number

# Set the prior values
pM = c(-3,0.2,0.001,-4,0.05) # prior means for parameters in a Gompertz model with a bathtub shape
pSD = c(1,0.05,0.001,1,0.01) # prior standard deviation for parameters in a Gompertz model with a bathtub shape

basta_ls <- list(id = basta_df$id,
                 byr_sample = basta_df$byr_sample,
                 dyr = NA, # test with -- rep(2020, nrow(basta_df)),
                 studyStart = basta_df$studyStart,
                 studyEnd = basta_df$studyEnd)#,
                 sex = basta_df$sex)


# begin the parallel loop
Perm.Out <- foreach(i=1:length(basta_df), .packages='BaSTA',.verbose = T) %dopar% { # start the parallel processing - change 'BaSTA' with the package you will be using in the parallel loop
  
  # Run survival model
  Out <- basta(basta_df[[i]],studyStart = basta_df$studyStart, # the model function will run through all lists in the input
               studyEnd = basta_df$studyEnd, nsim = sim, niter = iter, 
               burnin = warmup, thinning = thin, 
               thetaPriorMean=pM,thetaPriorSd=pSD,
               models = c("GO"), shape = c("bathtub"))
  
  return(Out)
}
# save the output
save(Perm.Out, file=paste0(output_path,"output.RData"))

stopCluster(cl)
