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

#basta_data$age_2021 <- ifelse(basta_data$dyr > 100, basta_data$dyr - basta_data$byr, 2021 - basta_data$byr)
#basta_data$censor_2021 <- ifelse(basta_data$dyr > 100, 'CENSORED', 'DEAD')
#hist(basta_data$age_2021)
#rethinking::dens(basta_data$age_2021, ylim = c(0, 0.07))
#rethinking::dens(basta_data$age_2021[basta_data$censor_2021 == 'DEAD'], add = T, col = 'red')
#rethinking::dens(basta_data$age_2021[basta_data$censor_2021 == 'CENSORED'], add = T, col = 'blue')

basta_data$dyr <- ifelse(basta_data$dyr < 100,  0, basta_data$dyr)
basta_data$byr <- ifelse(basta_data$byr < 1972, 0, basta_data$byr)

#basta_data$total_years_sighted <- rowSums(basta_data[,4:53])
#table(basta_data$total_years_sighted)

#### Gompertz simple - Only requires 2 parameters -- run ####
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

#### Gompertz Bathtub -- run using priors from killer whale paper ####
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
silerdist <- function(a0, a1, c, b0, b1, x) {  # The Siler (Gompertz-bathtub) distribution with priors for Basta 
  return( exp(a0-a1*x) + c + exp(b0 + b1*x) )
}

mortality <- silerdist(a0 = -3, a1 = 0.2, c = 0.001, b0 = -4 , b1 = 0.05, 1:100) ; plot(mortality)
mortality <- silerdist(a0 = 0, a1 = 0.2, c = 0.001, b0 = -4 , b1 = 0.05, 1:100)  ; plot(mortality)

gompertz_bathtub <- function(b0, b1, a0, a1, c, age){
  gompertz <- exp(b0 + b1*age)
  makeham <- c
  bathtub <- exp(a0 - a1*age)
  y <- gompertz + makeham + bathtub
  print(y)
}
prior_test <- data.frame(b1 = rep(-3,51),    b0 = rep(0.2,51),
                         a1 = rep(0.001,51), a0 = rep(-4,51),
                         c = rep(-0.05,51), age = 10:60)
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


#### Weibull simple -- run ####
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

weibull <- function(b0, b1, age){
  y <- (b0) * (b1^b0) * (age^(b0-1))
  return(y)
}
mortality <- weibull(b0 = 3, b1 = 2, age = 1:70)  ; plot(mortality)
mortality <- weibull(b0 = 5, b1 = 2, age = 1:70)  ; plot(mortality)
mortality <- weibull(b0 = 0.5, b1 = 0.2, age = 1:70)  ; plot(mortality)


pM <-  c(4,2) # prior means for parameters in a Weibull model with a simple shape
pSD <- c(0.5,0.5)  # prior standard deviation for parameters in a Weibull model with a simple shape
Wb.Si <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = 1, niter = iter,
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

#### WEIBULL BATHTUB -- not yet worked out priors ####
# first attempt ####
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

# second attempt ####
wb.dist <- function(a0, a1, c, b0, b1, age) {
  f1 <- b0 * b1^(b0) * age^(b0-1)
  f2 <- exp(a0 - a1*age)
  y <- f1 + f2 + c
  return(y)
}

mortality <- wb.dist(a0 = 1, a1 = 1, c = 1, b0 = 1, b1 = 1, 1:70)  ; plot(mortality)
mortality <- wb.dist(a0 = 10, a1 = 1, c = 1, b0 = 1, b1 = 1, 1:70)  ; plot(mortality)
mortality <- wb.dist(a0 = 1, a1 = 10, c = 1, b0 = 1, b1 = 1, 1:70)  ; plot(mortality)
mortality <- wb.dist(a0 = 1, a1 = 1, c = 10, b0 = 1, b1 = 1, 1:70)  ; plot(mortality)
mortality <- wb.dist(a0 = 1, a1 = 1, c = 1, b0 = 1, b1 = 10, 1:70)  ; plot(mortality)
mortality <- wb.dist(a0 = 0.1, a1 = 1, c = 1, b0 = 1, b1 = 1, 1:70)  ; plot(mortality)
mortality <- wb.dist(a0 = 1, a1 = 0.1, c = 1, b0 = 1, b1 = 1, 1:70)  ; plot(mortality)
mortality <- wb.dist(a0 = 1, a1 = 0.5, c = 1, b0 = 1, b1 = 1, 1:70)  ; plot(mortality)
mortality <- wb.dist(a0 = 1, a1 = 1, c = 1, b0 = 0.001, b1 = 1, 1:70)  ; plot(mortality)

mortality <- wb.dist(a0 = 1, a1 = 0.5, c = 1, b0 = 1, b1 = 0.1, 1:70)  ; plot(mortality)
mortality <- wb.dist(a0 = 10, a1 = 0.5, c = 1, b0 = 1, b1 = 0.1, 1:70)  ; plot(mortality)
mortality <- wb.dist(a0 = 1, a1 = 0.5, c = 10, b0 = 1, b1 = 0.1, 1:70)  ; plot(mortality)

mortality <- wb.dist(a0 = 1, a1 = 0.5, c = 1, b0 = 2, b1 = 0.1, 1:70)  ; plot(mortality)

mortality <- wb.dist(a0 = 1, a1 = 0.5, c = 1, b0 = 1.5, b1 = 0.05, 1:70)  ; plot(mortality)
mortality <- wb.dist(a0 = 1, a1 = 0.5, c = 1, b0 = 1.5, b1 = 0.05, 1:70)  ; plot(mortality)

mortality <- wb.dist(a0 = 1, a1 = 0.5, c = 10, b0 = 1, b1 = 0.1, 1:70)  ; plot(mortality)

# leave a0 and c at 1 as these only seem to affect the scale of the y axis
# leave a1 at 0.5 as this seems to set the bottom of the bathtub around 10 years old
# set b0 > 1 to obtain tail end of bathtub
# set b1 << 1 to avoid crazy rapid increase in mortality from age 10

mortality <- wb.dist(a0 = 1, a1 = 0.5, c = 1, b0 = 1.5, b1 = 0.1, 1:70)  ; plot(mortality)
mortality <- wb.dist(a0 = 1, a1 = 0.3, c = 1, b0 = 1.5, b1 = 0.1, 1:70)  ; plot(mortality)





to.plot <- data.frame(a0 = c(rep(1,70),rep(NA,6930)), a1 = c(rep(0.3,70),rep(NA,6930)), c = c(rep(1,70),rep(NA,6930)),
                      b0 = c(rep(1.5,70),rep(NA,6930)), b1 = c(rep(0.1,70),rep(NA,6930)), age = rep(1:70,100),
                      y = rep(NA,7000), iter = rep(1:100,each = 70))
to.plot$y <- wb.dist(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1, age = to.plot$age)
for(i in 1:99){
  to.plot$a0[(i*70+1):(i*70+70)] <- rep(rnorm(1, 1, 0.2),70)
  to.plot$a1[(i*70+1):(i*70+70)] <- rep(rnorm(1, 0.3, 0.1),70)
  to.plot$c[(i*70+1):(i*70+70)]  <- rep(rnorm(1, 1, 0.1),70)
  to.plot$b0[(i*70+1):(i*70+70)] <- rep(rnorm(1, 1.5, 0.3),70)
  to.plot$b1[(i*70+1):(i*70+70)] <- rep(rexp(1),70)
}

to.plot$y  <- wb.dist(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1, 1:70)

plot(NULL, xlab = 'age', ylab = 'mortality', xlim = c(0,70), ylim = c(0,50), las = 1)
for(i in 1:100){
  data <- to.plot[to.plot$iter == i, 6:7]
  lines(data$y ~ data$age, col = rgb(0,0,1,0.2))
}

to.plot <- data.frame(a0 = c(rep(1,70),rep(NA,6930)), a1 = c(rep(0.3,70),rep(NA,6930)), c = c(rep(1,70),rep(NA,6930)),
                      b0 = c(rep(1.5,70),rep(NA,6930)), b1 = c(rep(0.1,70),rep(NA,6930)), age = rep(1:70,100),
                      y = rep(NA,7000), iter = rep(1:100,each = 70))
to.plot$y <- wb.dist(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1, age = to.plot$age)
for(i in 1:99){
  to.plot$a0[(i*70+1):(i*70+70)] <- rep(rnorm(1, 2, 0.1),70)
  to.plot$a1[(i*70+1):(i*70+70)] <- rep(rnorm(1, 0.2, 0.1),70)
  to.plot$c[(i*70+1):(i*70+70)]  <- rep(rnorm(1, 1, 0.05),70)
  to.plot$b0[(i*70+1):(i*70+70)] <- rep(rnorm(1, 1, 0.2),70)
  to.plot$b1[(i*70+1):(i*70+70)] <- rep(rexp(1),70)
}

to.plot$y  <- wb.dist(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1, 1:70)

plot(NULL, xlab = 'age', ylab = 'mortality', xlim = c(0,70), ylim = c(0,20), las = 1)
for(i in 1:100){
  data <- to.plot[to.plot$iter == i, 6:7]
  lines(data$y ~ data$age, col = rgb(0,0,1,0.2))
}

to.plot <- data.frame(a0 = c(rep(1,70),rep(NA,6930)), a1 = c(rep(0.3,70),rep(NA,6930)), c = c(rep(1,70),rep(NA,6930)),
                      b0 = c(rep(1.5,70),rep(NA,6930)), b1 = c(rep(0.1,70),rep(NA,6930)), age = rep(1:70,100),
                      y = rep(NA,7000), iter = rep(1:100,each = 70))
to.plot$y <- wb.dist(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1, age = to.plot$age)
for(i in 1:99){
  to.plot$a0[(i*70+1):(i*70+70)] <- rep(rnorm(1, 2, 0.1),70)
  to.plot$a1[(i*70+1):(i*70+70)] <- rep(rnorm(1, 0.2, 0.1),70)
  to.plot$c[(i*70+1):(i*70+70)]  <- rep(rnorm(1, 1, 0.05),70)
  to.plot$b0[(i*70+1):(i*70+70)] <- rep(rnorm(1, 1.2, 0.2),70)
  to.plot$b1[(i*70+1):(i*70+70)] <- rep(rexp(1),70)
}

to.plot$y  <- wb.dist(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1, 1:70)

plot(NULL, xlab = 'age', ylab = 'mortality', xlim = c(0,70), ylim = c(0,20), las = 1)
for(i in 1:100){
  data <- to.plot[to.plot$iter == i, 6:7]
  lines(data$y ~ data$age, col = rgb(0,0,1,0.2))
}

to.plot <- data.frame(a0 = c(rep(1,70),rep(NA,6930)), a1 = c(rep(0.3,70),rep(NA,6930)), c = c(rep(1,70),rep(NA,6930)),
                      b0 = c(rep(1.5,70),rep(NA,6930)), b1 = c(rep(0.1,70),rep(NA,6930)), age = rep(1:70,100),
                      y = rep(NA,7000), iter = rep(1:100,each = 70))
to.plot$y <- wb.dist(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1, age = to.plot$age)
for(i in 1:99){
  to.plot$a0[(i*70+1):(i*70+70)] <- rep(rnorm(1, 2, 0.5),70)
  to.plot$a1[(i*70+1):(i*70+70)] <- rep(rnorm(1, 0.2, 0.1),70)
  to.plot$c[(i*70+1):(i*70+70)]  <- rep(rnorm(1, 1, 0.2),70)
  to.plot$b0[(i*70+1):(i*70+70)] <- rep(rnorm(1, 1.2, 0.4),70)
  to.plot$b1[(i*70+1):(i*70+70)] <- rep(rexp(1),70)
}

to.plot$y  <- wb.dist(a0 = to.plot$a0, a1 = to.plot$a1, c = to.plot$c, b0 = to.plot$b0, b1 = to.plot$b1, 1:70)

plot(NULL, xlab = 'age', ylab = 'mortality', xlim = c(0,70), ylim = c(0,20), las = 1)
for(i in 1:100){
  data <- to.plot[to.plot$iter == i, 6:7]
  lines(data$y ~ data$age, col = rgb(0,0,1,0.2))
}

pM <-  c(5,5,5,-1,5) # prior means for parameters in a Weibull model with a bathtub shape
pSD <- c(0.1,0.1,0.1,0.1,0.1)  # prior standard deviation for parameters in a Weibull model with a bathtub shape

Wb.Bt <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = 1, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               models = c("WE"), shape = c("bathtub"))
summary(Wb.Bt)
plot(Wb.Bt)                       # Plot traces
plot(Wb.Bt, plot.trace = FALSE)   # Plot survival and mortality curves




weibull <- function(b0, b1, age){
  y <- (b0) * (b1^b0) * (age^(b0-1))
  return(y)
}

wb.dist <- function(a0, a1, c, b0, b1, age) {
  f1 <- b0 * b1^(b0) * age^(b0-1)
  f2 <- exp(a0 - a1*age)
  y <- f1 + f2 + c
  return(y)
}

mortality <- weibull(b0 = 5, b1 = 2, age = 1:70)  ; plot(mortality)
mortality <- weibull(b0 = 0.5, b1 = 0.2, age = 1:70)  ; plot(mortality)

mortality <- wb.dist(a0 = 5, a1 = 2, c = 1, b0 = 5, b1 = 2, 1:70)  ; plot(mortality)
mortality <- wb.dist(a0 = 0.5, a1 = 0.2, c = 1, b0 = 0.5, b1 = 0.2, 1:70)  ; plot(mortality)
mortality <- wb.dist(a0 = 5, a1 = 2, c = 1, b0 = 0.5, b1 = 0.2, 1:70)  ; plot(mortality)
mortality <- wb.dist(a0 = 0.5, a1 = 0.2, c = 1, b0 = 5, b1 = 2, 1:70)  ; plot(mortality)


pM <-  c(0.5,0.2,1,5,2) # prior means for parameters in a Weibull model with a bathtub shape
pSD <- c(0.3,0.1,2,1,1)  # prior standard deviation for parameters in a Weibull model with a bathtub shape

Wb.Bt <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = 1, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               models = c("WE"), shape = c("bathtub"))
summary(Wb.Bt)
plot(Wb.Bt)                       # Plot traces
plot(Wb.Bt, plot.trace = FALSE)   # Plot survival and mortality curves

