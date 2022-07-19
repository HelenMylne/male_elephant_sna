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

# Set the arguments for the iterations
iter <- 500      # number of iterations
sim <- 1           # number of chains
warmup <- iter/2   # number of iterations to use during warmup/burn-in period
thin <- 10        # thinning number

# prep data
basta_data <- males
basta_data$casename <- paste0('M', basta_data$casename)
basta_data <- basta_data[-c(1:4,6,8,11,14:18,22,23,26:30,35,38,44:46,52,55,56,58,59,62,63,66,68,71,72,76,83,87,90,93,94,96,105,107,108,111,114,115,118,127,128,132,133,135,136,138,145,146,151,155,159,162,174,179,182,186,187,190,191,193,197,201,208,210,213,220,224,230,235,238,244,248,250,254,256,259,261,263,266,268,288,291,292,295,298,306,307,309,319,327,339,342,346,370,384,416,417,435,449,472,717,693,708,735,769,774),]  # elephants that have sightings before birth/after death
basta_data <- basta_data[basta_data$casename %in% sort(unique(paste0('M',sightings$id))),]

basta_data$dyr <- ifelse(basta_data$dyr < 100,  0, basta_data$dyr)
basta_data$byr <- ifelse(basta_data$byr < 1972, 0, basta_data$byr)

#### Gompertz simple - Only requires 2 parameters -- run ####
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
               models = c("GO"), shape = c("simple"))
summary(Go.Si)
plot(Go.Si)                       # Plot traces
plot(Go.Si, plot.trace = FALSE)   # Plot survival and mortality curves

#### Gompertz Bathtub -- run using priors from killer whale paper ####
gompertz_bt <- function(a0, a1, c, b0, b1, age){
  gompertz <- gompertz(b0, b1, age)
  bathtub <- exp(a0 - a1*age) + c + gompertz
  return(bathtub)
}

mortality <- gompertz_bt(a0 = -3, a1 = 0.2, c = 0.001, b0 = -4 , b1 = 0.05, 1:100) ; plot(mortality)
mortality <- gompertz_bt(a0 = 0, a1 = 0.2, c = 0.001, b0 = -4 , b1 = 0.05, 1:100)  ; plot(mortality)

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

pM <-  c(0,0.2,0.001,-4,0.05)  # prior means for parameters in a Gompertz model with a bathtub shape
pSD <- c(1,0.05,0.001,1,0.01)  # prior standard deviation for parameters in a Gompertz model with a bathtub shape
Go.Bt <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               models = c("GO"), shape = c("bathtub"))
summary(Go.Bt)
plot(Go.Bt)                       # Plot traces
plot(Go.Bt, plot.trace = FALSE)   # Plot survival and mortality curves

#### Weibull simple -- run ####
weibull <- function(b0, b1, age){
  y <- (b0) * (b1^b0) * (age^(b0-1))
  return(y)
}
mortality <- weibull(b0 = 3, b1 = 2, age = 1:70)  ; plot(mortality)
mortality <- weibull(b0 = 4, b1 = 2, age = 1:70)  ; plot(mortality)
mortality <- weibull(b0 = 0.5, b1 = 0.2, age = 1:70)  ; plot(mortality)
mortality <- weibull(b0 = 1.5, b1 = 2, age = 1:70)  ; plot(mortality)

pM <-  c(4,2) # prior means for parameters in a Weibull model with a simple shape
pSD <- c(1,0.5)  # prior standard deviation for parameters in a Weibull model with a simple shape
Wb.Si <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = 1, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               models = c("WE"), shape = c("simple"))
summary(Wb.Si)
plot(Wb.Si)                       # Plot traces
plot(Wb.Si, plot.trace = FALSE)   # Plot survival and mortality curves

#### WEIBULL BATHTUB -- not yet worked out priors ####
weibull_bt <- function(a0, a1, c, b0, b1, age) {
  f1 <- b0 * b1^(b0) * age^(b0-1)
  f2 <- exp(a0 - a1*age)
  y <- f1 + f2 + c
  return(y)
}

mortality <- weibull_bt(a0 = 1, a1 = 0.5, c = 1, b0 = 1.5, b1 = 0.1, 1:70)  ; plot(mortality)
mortality <- weibull_bt(a0 = 1, a1 = 0.3, c = 1, b0 = 1.5, b1 = 0.1, 1:70)  ; plot(mortality)
mortality <- weibull_bt(a0 = 2, a1 = 0.2, c = 1, b0 = 1, b1 = 1.2, 1:70)  ; plot(mortality)

to.plot <- data.frame(a0 = c(rep(1,70),rep(NA,6930)), a1 = c(rep(0.3,70),rep(NA,6930)), c = c(rep(1,70),rep(NA,6930)),
                      b0 = c(rep(1.5,70),rep(NA,6930)), b1 = c(rep(0.1,70),rep(NA,6930)), age = rep(1:70,100),
                      y = rep(NA,7000), iter = rep(1:100,each = 70))
to.plot$y <- weibull_bt(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1, age = to.plot$age)
for(i in 1:99){
  to.plot$a0[(i*70+1):(i*70+70)] <- rep(rnorm(1, 2, 0.5),70)
  to.plot$a1[(i*70+1):(i*70+70)] <- rep(rnorm(1, 0.2, 0.1),70)
  to.plot$c[(i*70+1):(i*70+70)]  <- rep(rnorm(1, 1, 0.2),70)
  to.plot$b0[(i*70+1):(i*70+70)] <- rep(rnorm(1, 1.2, 0.4),70)
  to.plot$b1[(i*70+1):(i*70+70)] <- rep(rexp(1),70)
}

to.plot$y  <- weibull_bt(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1, 1:70)

plot(NULL, xlab = 'age', ylab = 'mortality', xlim = c(0,70), ylim = c(0,20), las = 1)
for(i in 1:100){
  data <- to.plot[to.plot$iter == i, 6:7]
  lines(data$y ~ data$age, col = rgb(0,0,1,0.2))
}

pM <-  c(  2, 0.2,   1, 1.2, 4)  # prior means for parameters in a Weibull model with a bathtub shape
pSD <- c(0.5, 0.1, 0.2, 0.4, 1)  # prior standard deviation for parameters in a Weibull model with a bathtub shape

#Wb.Bt1 <- basta(basta_data, studyStart = Ss,
#               studyEnd = Se, nsim = sim, niter = iter,
#               burnin = warmup, thinning = thin,
#               thetaPriorMean = pM,thetaPriorSd = pSD,
#               models = c("WE"), shape = c("bathtub"))
#summary(Wb.Bt1)
#plot(Wb.Bt1)                       # Plot traces
#plot(Wb.Bt1, plot.trace = FALSE)   # Plot survival and mortality curves

#mortality <- weibull_bt(a0 = 5, a1 = 2, c = 1, b0 = 5, b1 = 2, 1:70)          ; plot(mortality)
#mortality <- weibull_bt(a0 = 0.5, a1 = 0.2, c = 1, b0 = 0.5, b1 = 0.2, 1:70)  ; plot(mortality)
#mortality <- weibull_bt(a0 = 5, a1 = 2, c = 1, b0 = 0.5, b1 = 0.2, 1:70)      ; plot(mortality)
#mortality <- weibull_bt(a0 = 0.5, a1 = 0.2, c = 1, b0 = 5, b1 = 2, 1:70)      ; plot(mortality)

mortality <- weibull_bt(a0 = 2, a1 = 0.2, c = 1, b0 = 1.2, b1 = 4, 1:70)      ; plot(mortality, las = 1)
mortality <- weibull_bt(a0 = 0.5, a1 = 0.2, c = 1, b0 = 1.2, b1 = 0.5, 1:70)  ; plot(mortality, las = 1)

pM <-  c(0.5,0.2,1,1.2,0.5) # prior means for parameters in a Weibull model with a bathtub shape
pSD <- c(0.3,0.1,1,0.5,0.3)  # prior standard deviation for parameters in a Weibull model with a bathtub shape

Wb.Bt3 <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               models = c("WE"), shape = c("bathtub"))
summary(Wb.Bt3)
plot(Wb.Bt3)                       # Plot traces
plot(Wb.Bt3, plot.trace = FALSE)   # Plot survival and mortality curves


plot(Wb.Bt1, plot.trace = FALSE)   # Plot survival and mortality curves
plot(Wb.Bt2, plot.trace = FALSE)   # Plot survival and mortality curves
plot(Wb.Bt3, plot.trace = FALSE)   # Plot survival and mortality curves



#### Weibull bathtub -- attempt 2 ####
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

to.plot <- data.frame(a0 = c(rep(1,   70), rep(NA,6930)),
                      a1 = c(rep(0.3, 70), rep(NA,6930)),
                      c  = c(rep(1,   70), rep(NA,6930)),
                      b0 = c(rep(1.5, 70), rep(NA,6930)),
                      b1 = c(rep(0.1, 70), rep(NA,6930)),
                      age = rep(1:70, 100),
                      y = rep(NA,7000), iter = rep(1:100,each = 70))
to.plot$y <- weibull_bt(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1, age = to.plot$age)
for(i in 1:99){
  to.plot$a0[(i*70+1):(i*70+70)] <- rep(rnorm(1, 1, 0.5),70)
  to.plot$a1[(i*70+1):(i*70+70)] <- rep(rnorm(1, 0.3, 0.1),70)
  to.plot$c[(i*70+1):(i*70+70)]  <- rep(rnorm(1, 1, 0.2),70)
  to.plot$b0[(i*70+1):(i*70+70)] <- rep(rnorm(1, 1.5, 0.4),70)
  to.plot$b1[(i*70+1):(i*70+70)] <- rep(rnorm(1, 0.1, 0.02),70)
}
to.plot$y  <- weibull_bt(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1, 1:70)
plot(NULL, xlab = 'age', ylab = 'mortality', xlim = c(0,70), ylim = c(0,10), las = 1)
for(i in 1:100){
  data <- to.plot[to.plot$iter == i, 6:7]
  lines(data$y ~ data$age, col = rgb(0,0,1,0.2))
}

pM <-  c(1,0.3,1,1.5,0.1) # prior means for parameters in a Weibull model with a bathtub shape
pSD <- c(0.5,0.1,1,0.4,0.02)  # prior standard deviation for parameters in a Weibull model with a bathtub shape
Wb.Bt1 <- basta(basta_data, studyStart = Ss,
                studyEnd = Se, nsim = sim, niter = iter,
                burnin = warmup, thinning = thin,
                thetaPriorMean = pM,thetaPriorSd = pSD,
                models = c("WE"), shape = c("bathtub"))
summary(Wb.Bt1)
plot(Wb.Bt1)                       # Plot traces
plot(Wb.Bt1, plot.trace = FALSE)   # Plot survival and mortality curves -- woah! no!

to.plot <- data.frame(a0 = c(rep(2,   70), rep(NA,6930)),
                      a1 = c(rep(0.2, 70), rep(NA,6930)),
                      c  = c(rep(1,   70), rep(NA,6930)),
                      b0 = c(rep(1.2, 70), rep(NA,6930)),
                      b1 = c(rep(0.1, 70), rep(NA,6930)),
                      age = rep(1:70, 100),
                      y = rep(NA,7000), iter = rep(1:100,each = 70))
to.plot$y <- weibull_bt(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1, age = to.plot$age)
for(i in 1:99){
  to.plot$a0[(i*70+1):(i*70+70)] <- rep(rnorm(1, 2, 0.5),70)
  to.plot$a1[(i*70+1):(i*70+70)] <- rep(rnorm(1, 0.2, 0.1),70)
  to.plot$c[(i*70+1):(i*70+70)]  <- rep(rnorm(1, 1, 0.2),70)
  to.plot$b0[(i*70+1):(i*70+70)] <- rep(rnorm(1, 1.2, 0.4),70)
  to.plot$b1[(i*70+1):(i*70+70)] <- rep(rnorm(1, 0.1, 0.02),70)
}
to.plot$y  <- weibull_bt(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1, 1:70)
plot(NULL, xlab = 'age', ylab = 'mortality', xlim = c(0,70), ylim = c(0,20), las = 1)
for(i in 1:100){
  data <- to.plot[to.plot$iter == i, 6:7]
  lines(data$y ~ data$age, col = rgb(0,0,1,0.2))
}

pM <-  c(2,0.2,1,1.2,0.1) # prior means for parameters in a Weibull model with a bathtub shape
pSD <- c(0.5,0.1,0.2,0.4,0.02)  # prior standard deviation for parameters in a Weibull model with a bathtub shape
Wb.Bt2 <- basta(basta_data, studyStart = Ss,
                studyEnd = Se, nsim = sim, niter = iter,
                burnin = warmup, thinning = thin,
                thetaPriorMean = pM,thetaPriorSd = pSD,
                models = c("WE"), shape = c("bathtub"))
summary(Wb.Bt2)
plot(Wb.Bt2)                       # Plot traces
plot(Wb.Bt2, plot.trace = FALSE)   # Plot survival and mortality curves -- very narrow 

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

mortality <- weibull_bt(a0 = 0.5, a1 = 0.2, c = 1, b0 = 1.2, b1 = 0.5, 1:70)  ; plot(mortality, las = 1)

pM <-  c(0.5,0.2,1,1.2,0.5) # prior means for parameters in a Weibull model with a bathtub shape
pSD <- c(0.3,0.1,1,0.5,0.3)  # prior standard deviation for parameters in a Weibull model with a bathtub shape

Wb.Bt3 <- basta(basta_data, studyStart = Ss,
                studyEnd = Se, nsim = sim, niter = iter,
                burnin = warmup, thinning = thin,
                thetaPriorMean = pM,thetaPriorSd = pSD,
                models = c("WE"), shape = c("bathtub"))
summary(Wb.Bt3)
plot(Wb.Bt3)                       # Plot traces
plot(Wb.Bt3, plot.trace = FALSE)   # Plot survival and mortality curves


plot(Wb.Bt1, plot.trace = FALSE)   # Plot survival and mortality curves
plot(Wb.Bt2, plot.trace = FALSE)   # Plot survival and mortality curves
plot(Wb.Bt3, plot.trace = FALSE)   # Plot survival and mortality curves
















