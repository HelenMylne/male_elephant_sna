#### Information #####
# To estimate the ages of individuals in other populations, first come up with a survival curve for the Amboseli elephants which have good birth and death data.
# survival curve options = exponential/logistic/Gompertz/Weibull, where shape = simple/bathtub/Makeham
# follow through Mia's code to determine which curve will best fit Amboseli data. she ran it on a high performance cluster, so it’s set up to run in parallel, but the basta model can also be run in parallel on its own (just like Stan). we will need to implement it as a custom function in Stan for the final analysis, given that we need to define it in the prior anyway. 

#### Load packages ####
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

################# Simulate population #################
#### simulate population ####
# using Dan's example from Weibull model (March 2022)
# First let's simulate an age distribution based on a Weibull distribution where we set the shape and scale parameters.
probs <- dweibull(1:60, shape = 1.2, scale = 30)
probs <- (probs/sum(probs))
plot(probs) 

# Now we create a fictional population with ages selected from that distribution. We also add an error for age estimation.
N <- 700 # Number of individuals
K <- 7   # Number of age bin categories
start_year <- 1972
end_year <- 2021
elephants_ls <- list(
  N = N,
  K = K, 
  age = sample(1:60, N, replace = TRUE, prob = probs) # Simulated ages drawn from the Weibull distribution plus an error in age estimation
)

# Next simulate observing ages with error and then binning them into age groups
E <- 3 # Error (in either direction) of age estimates
elephants_ls$age_guess <- elephants_ls$age + sample(-E:E, N, replace=TRUE)
elephants_ls$age_category_index <- sapply(elephants_ls$age_guess, function(x) which.max(x < c(5, 10, 15, 20, 25, 40, Inf)))
hist(elephants_ls$age_category)

# Let's look at the actual age verses the biologist assigned age and look at the thresholds. You can see that the chance of being mis-classified is lower if an actual age is in the middle of the age category. Also, there is mostly a bias towards ages within the age class actually being towards the lower end of the class.
data.frame(elephants_ls) %>%
  ggplot(aes(x=age, y=age_guess, col=factor(age_category_index))) +
    geom_point(size=4,alpha=0.6) +
    geom_vline(xintercept=c(5, 10, 15, 20, 25, 40, 60), col=factor(1:7), linetype="dashed", alpha=0.6) +
    theme_minimal() + 
    theme(legend.position = 'none')+
    xlab("Assigned age") + ylab("Age")

#### create simulated population data for basta ####
N_all <- 10000              # number of individuals in entire population from 1930 to 2022
N_sample <- 700             # number of individuals to sample from 1972 to 2021
popn <- data.frame(id = 1:N_all,
                   byr = sample(1930:2021, N_all, replace = TRUE),                 # sample birth years for population
                   dyr = rep(NA,N_all),
                   true_age = sample(1:60, N_all, replace = TRUE, prob = probs))   # sample ages at death for population from survival curve probabilities
popn$dyr <- popn$byr + popn$true_age
popn$alive_during_study <- ifelse(popn$dyr > start_year, 'yes', 'no')
popn <- popn[popn$alive_during_study == 'yes', 1:4]
popn$ind <- popn$byr + round(rnorm(1, mean = 14, sd = 1.5), 0)
popn$ind <- ifelse(popn$ind > popn$dyr, NA, popn$ind)

# sample sightings per year per elephant
popn$obs_1972 <- NA ; popn$obs_1973 <- NA ; popn$obs_1974 <- NA ; popn$obs_1975 <- NA ; popn$obs_1976 <- NA ; popn$obs_1977 <- NA ; popn$obs_1978 <- NA ; popn$obs_1979 <- NA ; popn$obs_1980 <- NA ; popn$obs_1981 <- NA ; popn$obs_1982 <- NA ; popn$obs_1983 <- NA ; popn$obs_1984 <- NA ; popn$obs_1985 <- NA ; popn$obs_1986 <- NA ; popn$obs_1987 <- NA ; popn$obs_1988 <- NA ; popn$obs_1989 <- NA ; popn$obs_1990 <- NA ; popn$obs_1991 <- NA ; popn$obs_1992 <- NA ; popn$obs_1993 <- NA ; popn$obs_1994 <- NA ; popn$obs_1995 <- NA ; popn$obs_1996 <- NA ; popn$obs_1997 <- NA ; popn$obs_1998 <- NA ; popn$obs_1999 <- NA ; popn$obs_2000 <- NA ; popn$obs_2001 <- NA ; popn$obs_2002 <- NA ; popn$obs_2003 <- NA ; popn$obs_2004 <- NA ; popn$obs_2005 <- NA ; popn$obs_2006 <- NA ; popn$obs_2007 <- NA ; popn$obs_2008 <- NA ; popn$obs_2009 <- NA ; popn$obs_2010 <- NA ; popn$obs_2011 <- NA ; popn$obs_2012 <- NA ; popn$obs_2013 <- NA ; popn$obs_2014 <- NA ; popn$obs_2015 <- NA ; popn$obs_2016 <- NA ; popn$obs_2017 <- NA ; popn$obs_2018 <- NA ; popn$obs_2019 <- NA ; popn$obs_2020 <- NA ; popn$obs_2021 <- NA

(years <- 2022-1972)    # total number of possible years of observation
for(i in 1:nrow(popn)){ 
  for(j in 1:years){    # simulate a string of sightings (if reach independence then sample every year while non-independent, then randomly with poisson distribution afterwards, if don't reach independence then just sample every year birth until death)
    if(is.na(popn$ind[i]) == FALSE){
      sighting_years <- c((popn$byr[i]+1):popn$ind[i],
                          sample(x = (popn$ind[i]-1):popn$dyr[i], # i-1 rather than i to prevent random values being selected when ind = dyr. -1 rather than +1 because will have been seen in previous year anyway.
                                 size = rpois(1,lambda = 8), replace = T))
    } else {sighting_years <- (popn$byr[i]+1):popn$dyr[i] }
    popn[i,j+5] <- ifelse(length(which(sighting_years == j+1971)) > 0, 1, 0)
    rm(sighting_years)
  }
}
popn <- popn[rowSums(popn[6:55]) > 0,] # trim down to only individuals that were sighted

# separate out individuals that reach independence before the end of the study
indes <- popn[popn$dyr > popn$ind & popn$ind < 2022,]
indes <- indes[!is.na(indes$id),]

# sample 700 individuals from this population -- similar to size of ANP population
indes_sampled <- indes[sort(sample(1:nrow(indes), N_sample, replace = F)),]
indes_sampled <- indes_sampled[!is.na(indes_sampled$id),]

# calculate date individuals were first observed independently -- this will become "birth" date
indes_sampled$first_inde_sighting <- NA
for(i in 1:nrow(indes_sampled)){
  for(j in 6:55){
    indes_sampled$first_inde_sighting[i] <- ifelse(is.na(indes_sampled$first_inde_sighting[i]) == TRUE,
                                                   ifelse(indes_sampled[i,j] == 1 & 
                                                            j+1966 > indes_sampled$ind[i],
                                                          j+1966, NA),
                                                   indes_sampled$first_inde_sighting[i])
  }
}
summary(indes_sampled$first_inde_sighting - indes_sampled$ind)
plot(indes_sampled$first_inde_sighting - indes_sampled$ind, cex = 0.5, pch = 19)
indes_sampled <- indes_sampled[!is.na(indes_sampled$first_inde_sighting),c(1:3,5,56,6:55)] # remove all elephants not sampled post-independence

# set up data set so that year of first sighting = "birth" year
indes_basta <- indes_sampled[,c(1,5,3,6:55,2)] # rearrange data so that "birth year" column is actually year first sighted post-independence
colnames(indes_basta)[2] <- 'byr_apparent'
for(i in 1:nrow(indes_basta)){ # write over basta sightings so that sightings pre-independence are removed
  for(j in 1:years){
    if(indes_basta[i,j+3] == 1 & indes_basta$byr_apparent[i] >= j+1971 ){
      indes_basta[i,j+3] <- 0
    }
    else { indes_basta[i,j+3] <- indes_basta[i,j+3] }
  }
}

# set "birth year" to 0 for anything born more than 10 years prior to the start or observed as already independent in 1972, otherwise "birth" year = first sighting
indes_basta$byr_apparent <- ifelse(indes_basta$byr < 1962, 0,
                                   ifelse(indes_basta$byr_apparent < 1973, 0, indes_basta$byr_apparent))
#indes_basta$dyr <- ifelse(indes_basta$dyr > end_year,   0, indes_basta$dyr) # previously having death dates after the end of the study caused problems, but now it is objecting to having 0s in these rows and looking at the sim1 data stored in BaSTA, this step appears to be unnecessary
indes_basta <- indes_basta[,1:53]  # remove byr column

# create equivalent data frame for if we have all of the information about the elephants
all_basta <- popn[sort(sample(1:nrow(popn), nrow(indes_basta), replace = F)), c(1:3, 6:55)]
all_basta$byr <- ifelse(all_basta$byr < start_year, 0, all_basta$byr)
#all_basta$dyr <- ifelse(all_basta$dyr > end_year,   0, all_basta$dyr) # same as above -- removed as appears to be unnecessary in sim1 data and want this to be comparable to the data using only independent elephants

################# test basta on simulated population #################
sim <- 4 ; iter <- 500 ; warmup <- iter/2 ; thin <- 10 # set basta arguments
pM <-  c(2, 25)
pSD <- c(1.5, 2)
all_data <- basta(all_basta,
                  studyStart = start_year, studyEnd = end_year,
                  nsim = sim, niter = iter,
                  burnin = warmup, thinning = thin,
                  thetaPriorMean = pM,thetaPriorSd = pSD,
                  model = "WE", shape = "simple")
summary(all_data)                    # close-ish for the shape, scale it returns precision = 1/variance so 1/b1 is not too bad
plot(all_data)                       # Plot traces -- slow to converge, but good
plot(all_data, plot.trace = FALSE)   # Plot survival and mortality curves -- don't look good, but try again with more samples and see they improve

sim_inde <- basta(indes_basta,
                  studyStart = start_year, studyEnd = end_year,
                  nsim = sim, niter = iter,
                  burnin = warmup, thinning = thin,
                  thetaPriorMean = pM,thetaPriorSd = pSD,
                  model = "WE", shape = "simple")
summary(sim_inde)                    # actually better than with all the data! 1/b1 is a bit further away
plot(sim_inde)                       # Plot traces -- converges better
plot(sim_inde, plot.trace = FALSE)   # Plot survival and mortality curves -- look much better

#### Repeat with longer chain and multiple chains
sim <- 4 ; iter <- 10000 ; warmup <- iter/2 ; thin <- 100
all_data <- basta(all_basta,
                  studyStart = start_year, studyEnd = end_year,
                  nsim = sim, niter = iter,
                  burnin = warmup, thinning = thin,
                  thetaPriorMean = pM,thetaPriorSd = pSD,
                  model = "WE", shape = "simple")
summary(all_data)                    # same as above, slightly overestimates b0, slightly underestimates 1/b1. DIC = 8793
plot(all_data)                       # Plot traces -- much better
plot(all_data, plot.trace = FALSE)   # Plot survival and mortality curves -- look good

sim_inde <- basta(indes_basta,
                  studyStart = start_year, studyEnd = end_year,
                  nsim = sim, niter = iter,
                  burnin = warmup, thinning = thin,
                  thetaPriorMean = pM,thetaPriorSd = pSD,
                  model = "WE", shape = "simple")
summary(sim_inde)                    # same as above, better for shape, slightly worse for scale. DIC = 5428 (but much fewer data points)
plot(sim_inde)                       # Plot traces -- good
plot(sim_inde, plot.trace = FALSE)   # Plot survival and mortality curves -- good

rm(all_basta, elephants_ls, indes, indes_sampled, E, K, N, N_all, N_sample, probs)

################# Prepare real data for modelling #################
#### Read in Amboseli data ####
sightings <- read_csv('data_processed/anp_sightings_updated_22.06.22.csv')
sightings$year <- lubridate::year(sightings$obs_date)
sightings <- sightings[,c('casename','year')]

males <- readxl::read_excel('data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% 
  janitor::clean_names() %>% 
  select(casename, byr, dyr)
str(males)
males <- males[males$casename %in% sightings$casename,]

#### Create data frame for basta ####
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

#### Prepare input data ####
# prep data
basta_data <- males
basta_data$count_years <- rowSums(basta_data[4:53])
summary(basta_data$count_years) # no individuals never observed, as it should be. No NAs.

#### outdated code chunk but which might be useful later ####
#basta_data$count_sightings_old <- NA
#basta_data$count_sightings_new <- NA
#basta_data$count_years_w_sightings_old <- NA
#basta_data$count_years_w_sightings_new <- NA
#for(i in 1:nrow(basta_data)){
#  old <- initial[initial$casename == basta_data$casename[i],]
#  new <- updated[updated$casename == basta_data$casename[i],]
#  basta_data$count_sightings_old[i] <- nrow(old)
#  basta_data$count_sightings_new[i] <- nrow(new)
#  basta_data$count_years_w_sightings_old[i] <- length(unique(lubridate::year(old$obs_date)))
#  basta_data$count_years_w_sightings_new[i] <- length(unique(lubridate::year(new$obs_date)))
#}
#table(basta_data$count_sightings_old)
#table(basta_data$count_sightings_new)

#updated <- read_csv('data_processed/anp_sightings_updated_22.06.22.csv')
#initial <- readxl::read_xlsx('data_raw/Raw_ATE_Sightings_Fishlock220224.xlsx') %>% janitor::clean_names()
#check <- basta_data[1,] # full row never sighted
#upd266 <- updated[updated$casename == check266$casename[1],]   # no sightings in updated data
#ini266 <- initial[initial$casename == check266$casename[1],]   # no sightings in initial data

#### set birth and death years ####
basta_data$birth_year <- ifelse(basta_data$byr < 1972, 0, basta_data$byr)
basta_data$death_year <- ifelse(basta_data$dyr < 100,  0, basta_data$dyr)

basta_data <- basta_data[,c(1,55,56,4:53,2:3,54)]

#### remove elephants with incorrect sightings ####
### code chunk sent to Vicki to highlight males with problems
#old_sightings <- readxl::read_excel('data_raw/Raw_ATE_Sightings_Fishlock220224.xlsx') %>% janitor::clean_names()
#old_males <- readxl::read_excel('data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()
#old_males <- old_males[old_males$casename %in% old_sightings$casename,]
#old_males$first <- NA ; old_males$last <- NA
#for(i in 1:nrow(males)){
#  casename <- old_sightings[old_sightings$casename == old_males$casename[i],]
#  old_males$first[i] <- lubridate::year(casename$obs_date[1])
#  old_males$last[i] <- lubridate::year(casename$obs_date[nrow(casename)])
#}
#m <- old_males[,c(2,6,9,21,22)]
#m$dyr_living <- ifelse(m$dyr < 100, 3000, m$dyr) # give living elephants impossibly high death year so don't get counted as dying before their last sighting
#m$error <- ifelse(m$byr > m$first, 'before_birth',
#                  ifelse(m$dyr_living < m$last, 'after_death','good'))
#table(m$error)
#before <- m[m$error == 'before_birth',]
#after <- m[m$error == 'after_death',]

check_birth <- basta_data[c(638,648,665,682,686),]
#which(check_birth$casename != before$casename)# these are the same individuals I previously had problems with -- Vicki has not managed to correct these. REMOVE FROM DATA.
#rm(before)

check_death <- basta_data[c(3,5,7,13,15,16,17,44,54,55,57,61,65,67,71,82,87,89,92,95,104,106,107,110,113,117,126,127,131,132,134,135,137,145,150,154,158,161,178,181,186,189,190,207,209,218,228,242,246,252,254,256,263,288,291,294,323,335,338,365,378,410,411,429,653),]
#nrow(after) - nrow(check_death)          # 46 elephants corrected from previous
#check_death$casename %in% after$casename # all current issues were there previously - no new ones have arisen -- REMOVE JUST THESE ONES FROM DATA, NOT ALL OF THE PREVIOUS ONES
#rm(after, m, old_males, old_sightings)

basta_data <- anti_join(basta_data, check_birth, by = 'casename')
basta_data <- anti_join(basta_data, check_death, by = 'casename')
rm(check_birth, check_death)

#### convert birth years to year first sighted ####
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

basta_data$age_first_sighted <- basta_data$year_first_sighted - basta_data$byr
hist(basta_data$age_first_sighted)
basta_data <- basta_data[,c(1,57,2,54,3:53,55:56,58)]

#### remove elephants dispersing at abnormal ages ####
young <- basta_data[basta_data$byr > 1962,]
mean_age_ind <- mean(young$age_first_sighted)
stdv_age_ind <-   sd(young$age_first_sighted)

basta_data$remove <- ifelse(basta_data$byr > 1962 & basta_data$age_first_sighted > (mean_age_ind+stdv_age_ind) | 
                            basta_data$byr > 1962 & basta_data$age_first_sighted < (mean_age_ind-stdv_age_ind),
                            'remove', 'keep')
table(basta_data$remove)
basta_data <- basta_data[basta_data$remove == 'keep',c(1:58)]

basta_data$byr_apparent <- ifelse(basta_data$byr < 1962, 0,
                                  ifelse(basta_data$year_first_sighted == 1971, 0,
                                         basta_data$year_first_sighted))

basta_data <- basta_data[,c(1,59,5:55,2:4,56:58)]
basta_anp <- basta_data[,c(1:53)]

################# Single chain, 500 iterations only #################
# Set the arguments for the iterations
sim <- 1; iter <- 500 ; warmup <- iter/2 ; thin <- 10
start_year <- 1972 ; end_year <- 2021

#### Gompertz simple ####
gompertz <- function(b0, b1, age){
  y <- exp(b0 + b1*age)
  return(y)
}

b0 <- -4
b1 <- 0.05
mortality <- gompertz(b0 = b0 , b1 = b1, 1:50) ; plot(mortality, las = 1)
to.plot <- data.frame(b0 = c(rep(b0, 50), rep(NA,4950)),
                      b1 = c(rep(b1, 50), rep(NA,4950)),
                      age = rep(1:50, 100),
                      y = rep(NA,5000),
                      iter = rep(1:100,each = 50))
for(i in 1:99){
  to.plot$b0[(i*50+1):(i*50+50)] <- rep(rnorm(1, b0, 2),50)
  to.plot$b1[(i*50+1):(i*50+50)] <- rep(rnorm(1, b1, 0.2),50)
}
to.plot$y  <- gompertz(b0 = to.plot$b0, b1 = to.plot$b1, age = to.plot$age)
plot(NULL, xlab = 'age', ylab = 'mortality', xlim = c(0,50), ylim = c(0,10), las = 1)
for(i in 1:100){
  data <- to.plot[to.plot$iter == i,]
  lines(data$y ~ data$age, col = rgb(0,0,1,0.2))
}

pM.Gs <-  c(b0,b1)  # prior means
pSD.Gs <- c(2,0.2)  # prior standard deviation
Go.Si <- basta(basta_anp, studyStart = start_year,
               studyEnd = end_year, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM.Gs,thetaPriorSd = pSD.Gs,
               model = "GO", shape = "simple")
summary(Go.Si)
plot(Go.Si)                       # Plot traces -- not great. might just be because exploring only a small area, but even so these are pretty bad. check with the longer curves
plot(Go.Si, plot.trace = FALSE)   # Plot survival and mortality curves -- mortality might be a bit too low 

#### Gompertz bathtub ####
gompertz_bt <- function(a0, a1, c, b0, b1, age){
  gompertz <- gompertz(b0, b1, age)
  bathtub <- exp(a0 - a1*age) + c + gompertz
  return(bathtub)
}

a0 <- -3
a1 <- 0.2
c <- 0.001
b0 <- -4
b1 <- 0.05
mortality <- gompertz_bt(a0 = a0, a1 = a1, c = c, b0 = b0 , b1 = b1, 1:50) ; plot(mortality) # original prior taken from Killer Whale paper
to.plot <- data.frame(a0 = c(rep(a0, 50), rep(NA,4950)),
                      a1 = c(rep(a1, 50), rep(NA,4950)),
                      c  = c(rep(c,  50), rep(NA,4950)),
                      b0 = c(rep(b0, 50), rep(NA,4950)),
                      b1 = c(rep(b1, 50), rep(NA,4950)),
                      age = rep(1:50,100),
                      y = rep(NA,5000),
                      iter = rep(1:100, each = 50))
for(i in 1:99){
  to.plot$a0[(i*50+1):(i*50+50)] <- rep(rnorm(1, a0, 1),    50)
  to.plot$a1[(i*50+1):(i*50+50)] <- rep(rnorm(1, a1, 0.1),  50)
  to.plot$c[ (i*50+1):(i*50+50)] <- rep(rnorm(1, c,  0.01), 50)
  to.plot$b0[(i*50+1):(i*50+50)] <- rep(rnorm(1, b0, 2),    50)
  to.plot$b1[(i*50+1):(i*50+50)] <- rep(rnorm(1, b1, 0.05), 50)
}
to.plot$y  <- gompertz_bt(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1,
                          age = to.plot$age)
plot(NULL, xlab = 'age', ylab = 'mortality', xlim = c(0,50), ylim = c(0,5), las = 1)
for(i in 1:100){
  data <- to.plot[to.plot$iter == i,]
  lines(data$y ~ data$age, col = rgb(0,0,1,0.2))
}

pM.Gb_kw <-  c(a0,a1,c,b0,b1)        # prior means
pSD.Gb_kw <- c(1,0.1,0.01,2,0.05)    # prior standard deviations
Go.Bt_kw <- basta(basta_anp, studyStart = start_year,
                  studyEnd = end_year, nsim = sim, niter = iter,
                  burnin = warmup, thinning = thin,
                  thetaPriorMean = pM.Gb_kw, thetaPriorSd = pSD.Gb_kw,
                  model = "GO", shape = "bathtub")
summary(Go.Bt_kw)
plot(Go.Bt_kw)                       # Plot traces -- looks a bit better than simple, but still a bit wandery
plot(Go.Bt_kw, plot.trace = FALSE)   # Plot survival and mortality curves -- better than simple

par(mfrow = c(1,1))

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
for(i in 1:99){
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
Go.Bt <- basta(basta_anp, studyStart = start_year,
               studyEnd = end_year, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM.Gb,thetaPriorSd = pSD.Gb,
               model = "GO", shape = c("bathtub"))
summary(Go.Bt)
plot(Go.Bt)                       # Plot traces -- looks a bit better than simple, but still a bit wandery
plot(Go.Bt, plot.trace = FALSE)   # Plot survival and mortality curves -- better than simple
par(mfrow = c(1,1))

#### Weibull simple ####
weibull <- function(b0, b1, age){
  y <- (b0) * (b1^b0) * (age^(b0-1))
  return(y)
}

b0 <- 2
b1 <- 0.2
mortality <- weibull(b0 = b0 , b1 = b1, 1:50) ; plot(mortality, las = 1)
to.plot <- data.frame(b0 = c(rep(b0, 50), rep(NA,4950)),
                      b1 = c(rep(b1, 50), rep(NA,4950)),
                      age = rep(1:50, 100),
                      y = rep(NA,5000),
                      iter = rep(1:100,each = 50))
for(i in 1:99){
  to.plot$b0[(i*50+1):(i*50+50)] <- rep(rnorm(1, b0, 0.4),50)
  to.plot$b1[(i*50+1):(i*50+50)] <- rep(rnorm(1, b1, 0.2),50)
}
to.plot$y  <- weibull(b0 = to.plot$b0, b1 = to.plot$b1, age = to.plot$age)
plot(NULL, xlab = 'age', ylab = 'mortality', xlim = c(0,50), ylim = c(0,35), las = 1)
for(i in 1:100){
  data <- to.plot[to.plot$iter == i,]
  lines(data$y ~ data$age, col = rgb(0,0,1,0.2))
}

pM.Ws <-  c(b0, b1)   # prior means
pSD.Ws <- c(0.4,0.2)  # prior standard deviations
Wb.Si <- basta(basta_anp, studyStart = start_year,
               studyEnd = end_year, nsim = 1, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM.Ws,thetaPriorSd = pSD.Ws,
               model = "WE", shape = c("simple"))
summary(Wb.Si)
plot(Wb.Si)                       # Plot traces -- take a while to stabilise but then fine
plot(Wb.Si, plot.trace = FALSE)   # Plot survival and mortality curves -- look good

### priors Dan suggested
b0 <- 1.3
b1 <- 30
mortality <- weibull(b0 = b0 , b1 = b1, 1:50) ; plot(mortality, las = 1)
to.plot <- data.frame(b0 = c(rep(b0, 50), rep(NA,4950)),
                      b1 = c(rep(b1, 50), rep(NA,4950)),
                      age = rep(1:50, 100),
                      y = rep(NA,5000),
                      iter = rep(1:100,each = 50))
for(i in 1:99){
  to.plot$b0[(i*50+1):(i*50+50)] <- rep(rnorm(1, b0, 0.4),50)
  to.plot$b1[(i*50+1):(i*50+50)] <- rep(rnorm(1, b1, 0.2),50)
}
to.plot$y  <- weibull(b0 = to.plot$b0, b1 = to.plot$b1, age = to.plot$age)
plot(NULL, xlab = 'age', ylab = 'mortality', xlim = c(0,50), ylim = c(0,3500), las = 1)
for(i in 1:100){
  data <- to.plot[to.plot$iter == i,]
  lines(data$y ~ data$age, col = rgb(0,0,1,0.2))
}

pM.Ws <-  c(b0, b1)   # prior means
pSD.Ws <- c(0.2, 5)   # prior standard deviation for parameters in a Weibull model with a simple shape -- smaller values makes no visible difference to final curves
Wb.Si2 <- basta(basta_anp, studyStart = start_year,
               studyEnd = end_year, nsim = 1, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM.Ws,thetaPriorSd = pSD.Ws,
               model = "WE", shape = c("simple"))
summary(Wb.Si2)                    # very similar the one above -- same overall shape just on a different scale?
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

a0 <- 0.5
a1 <- 0.2
c <- 1
b0 <- 1.4
b1 <- 0.5
to.plot <- data.frame(a0 = c(rep(a0, 50), rep(NA,4950)),
                      a1 = c(rep(a1, 50), rep(NA,4950)),
                      c  = c(rep(c,  50), rep(NA,4950)),
                      b0 = c(rep(b0, 50), rep(NA,4950)),
                      b1 = c(rep(b1, 50), rep(NA,4950)),
                      age = rep(1:50, 100),
                      y = rep(NA,5000), iter = rep(1:100,each = 50))
mortality <- weibull_bt(a0 = a0, a1 = a1, c = c, b0 = b0 , b1 = b1, 1:50)  ; plot(mortality, las = 1)
to.plot$y <- weibull_bt(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1, age = to.plot$age)
for(i in 1:99){
  to.plot$a0[(i*50+1):(i*50+50)] <- rep(rnorm(1, a0, 0.3),50)
  to.plot$a1[(i*50+1):(i*50+50)] <- rep(rnorm(1, a1, 0.1),50)
  to.plot$c[ (i*50+1):(i*50+50)] <- rep(rnorm(1, c,  1),  50)
  to.plot$b0[(i*50+1):(i*50+50)] <- rep(rnorm(1, b0, 0.5),50)
  to.plot$b1[(i*50+1):(i*50+50)] <- rep(rnorm(1, b1, 0.3),50)
}
to.plot$y  <- weibull_bt(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1, 1:50)
plot(NULL, xlab = 'age', ylab = 'mortality', xlim = c(0,50), ylim = c(0,20), las = 1)
for(i in 1:100){
  data <- to.plot[to.plot$iter == i,]
  lines(data$y ~ data$age, col = rgb(0,0,1,0.2))
}

pM.Wb <-  c(a0, a1, c, b0, b1)  # prior means
pSD.Wb <- c(0.3,0.1,1,0.5,0.3)  # prior standard deviations
Wb.Bt <- basta(basta_anp, studyStart = start_year,
                studyEnd = end_year, nsim = sim, niter = iter,
                burnin = warmup, thinning = thin,
                thetaPriorMean = pM.Wb,thetaPriorSd = pSD.Wb,
                model = "WE", shape = c("bathtub"))
summary(Wb.Bt)
plot(Wb.Bt)                       # Plot traces
plot(Wb.Bt, plot.trace = FALSE)   # Plot survival and mortality curves

################# Repeat with 4 chains and more iterations #################
# Set the arguments for the iterations
start_year <- 1972
end_year <- 2021
iter <- 50000      # number of iterations
sim <- 4           # number of chains
warmup <- iter/2   # number of iterations to use during warmup/burn-in period
thin <- 100        # thinning number

#### Gompertz simple -- DIC = 10248.61 with all elephants / 9694.605 with only those that were first seen independently between 8 and 20 years old ####
Go.Si <- basta(basta_anp, studyStart = start_year,
               studyEnd = end_year, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM.Gs,thetaPriorSd = pSD.Gs,
               model = "GO", shape = "simple")
summary(Go.Si)
#     Jump.sds Prior.means Prior.sds
#b0 0.19077737       -4.00       2.0
#b1 0.00911661        0.05       0.2
# KLDC was not calculated due to lack of convergence, or because covariates were not included in the model.
#        Estimate   StdErr Lower95%CI Upper95%CI SerAutocor UpdateRate PotScaleReduc
#b0      -3.62296 0.092487    -3.8017   -3.44899    0.23762     0.2519        0.9982
#b1       0.05189 0.005245     0.0422    0.06223    0.38390     0.2503        0.9998
#pi.1972  0.78836 0.004515     0.7794    0.79721   -0.02033     1.0000        0.9981
# Appropriate convergence reached for all parameters.
plot(Go.Si)
plot(Go.Si, plot.trace = FALSE)

#### Gompertz bathtub -- DIC = 10195.55 / 9644.401 (Killer Whale prior) or 10162 / 9601.075 (adapted prior) ####
Go.Bt_kw <- basta(basta_anp, studyStart = start_year,
                  studyEnd = end_year, nsim = sim, niter = iter,
                  burnin = warmup, thinning = thin,
                  thetaPriorMean = pM.Gb_kw, thetaPriorSd = pSD.Gb_kw,
                  model = "GO", shape = c("bathtub"))
summary(Go.Bt_kw)
#     Jump.sds Prior.means Prior.sds
#a0 1.63915656      -3.000      1.00
#a1 0.48229931       0.200      0.10
#c  0.01152920       0.001      0.01
#b0 0.23094553      -4.000      2.00
#b1 0.01019753       0.050      0.05
#KLDC was not calculated due to lack of convergence, or because covariates were not included in the model.
#         Estimate   StdErr Lower95%CI Upper95%CI SerAutocor UpdateRate PotScaleReduc
#a0      -4.35505 0.527073 -5.5372849   -3.51693   0.036405     0.2486        0.9981
#a1       0.17545 0.091076  0.0211328    0.36399   0.017531     0.2410        1.0017
#c        0.01000 0.005971  0.0007024    0.02284   0.269089     0.2504        1.0029
#b0      -4.39121 0.364851 -5.1761041   -3.82022   0.664460     0.2496        1.0078
#b1       0.07744 0.012500  0.0561452    0.10384   0.682329     0.2583        1.0166
#pi.1972  0.78820 0.004443  0.7795289    0.79657  -0.008825     1.0000        1.0006
# Appropriate convergence reached for all parameters.
plot(Go.Bt_kw)                      # chains, esp. b0 and b1, a little wander-y
plot(Go.Bt_kw, plot.trace = FALSE)  # mortality rate seems much higher (probably makes more sense), overall fairly linear survival

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

#### Weibull simple -- DIC = 10419.9 / 9898.886 ####
Wb.Si <- basta(basta_anp, studyStart = start_year,
               studyEnd = end_year, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM.Ws,thetaPriorSd = pSD.Ws,
               model = "WE", shape = c("simple"))
summary(Wb.Si)
#      Jump.sds Prior.means Prior.sds
#b0 0.213910536         2.0       0.4
#b1 0.006648478         0.2       0.2
# KLDC was not calculated due to lack of convergence, or because covariates were not included in the model.
#        Estimate   StdErr Lower95%CI Upper95%CI SerAutocor UpdateRate PotScaleReduc
#b0       1.46039 0.060343    1.34456    1.57901    0.32192     0.2507         1.017
#b1       0.05148 0.001723    0.04817    0.05472    0.21411     0.2552         1.007
#pi.1972  0.78856 0.004381    0.78006    0.79667    0.07894     1.0000         1.003
# Appropriate convergence reached for all parameters.
plot(Wb.Si)                     # good
plot(Wb.Si, plot.trace = FALSE) # maybe a bit steep on the initial mortality given that we are starting from age 14ish

#### Weibull bathtub -- DIC = 10283.44 / 9721.131 ####
Wb.Bt <- basta(basta_anp, studyStart = start_year,
               studyEnd = end_year, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM.Wb,thetaPriorSd = pSD.Wb,
               model = "WE", shape = c("bathtub"))
summary(Wb.Bt)
#      Jump.sds Prior.means Prior.sds
#a0 0.430899474         0.5       0.3
#a1 0.184761855         0.2       0.1
#c  0.009965117         1.0       1.0
#b0 0.393429881         1.4       0.5
#b1 0.004332755         0.5       0.3
# KLDC was not calculated due to lack of convergence, or because covariates were not included in the model.
#        Estimate   StdErr Lower95%CI Upper95%CI SerAutocor UpdateRate PotScaleReduc
#a0      -2.46590 0.128153 -2.722e+00  -2.229259   0.026728     0.2640        1.0023
#a1       0.32732 0.067257  2.003e-01   0.463862   0.078014     0.2500        1.0017
#c        0.00267 0.002402  6.901e-05   0.009221   0.003755     0.2605        1.0012
#b0       2.24187 0.174338  1.913e+00   2.596075   0.358325     0.2570        1.0066
#b1       0.04189 0.001658  3.879e-02   0.044929   0.183821     0.2483        0.9984
#pi.1972  0.78891 0.004390  7.805e-01   0.797469   0.042629     1.0000        1.0061
# Appropriate convergence reached for all parameters.
plot(Wb.Bt)                     # good
plot(Wb.Bt, plot.trace = FALSE) # more uncertain, again possibly too steep on the initial mortality


