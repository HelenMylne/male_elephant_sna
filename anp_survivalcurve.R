#### Information #####
# To estimate the ages of individuals in other populations, first come up with a survival curve for the Amboseli elephants which have good birth and death data.
# survival curve options = exponential/logistic/Gompertz/Weibull, where shape = simple/bathtub/Makeham
# follow through Mia's cde to determine which curve will best fit Amboseli data. she ran it on a high performance cluster, so it’s set up to run in parallel, but the basta model can also be run in parallel on its own (just like Stan). we will need to implement it as a custom function in Stan for the final analysis, given that we need to define it in the prior anyway. 

#emo::ji('elephant')
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


#### basta functions to try and run just part of the model and find the bug ####
# getAnywhere(.BuildAliveMatrix)
.BuildAliveMatrix <- function (f, l, dataObj) 
{
  Fm <- dataObj$Tm - f
  Fm[Fm >= 0] <- 1
  Fm[Fm < 0] <- 0
  Lm <- dataObj$Tm - l
  Lm[Lm <= 0] <- -1
  Lm[Lm > 0] <- 0
  return(Fm * (-Lm))
}

# getAnywhere(.CalcDiagnost)
.CalcDiagnost <- function (bastaOut, algObj, covObj, defTheta, fullParObj, dataObj, 
                           parsIni, parsCovIni, agesIni, postIni, CalcSurv, priorAgeObj) 
{
  thinned <- seq(algObj$burnin, algObj$niter, algObj$thinning)
  nthin <- length(thinned)
  parMat <- bastaOut[[1]]$par[thinned, ]
  fullParMat <- bastaOut[[1]]$par[algObj$burnin:algObj$niter, 
  ]
  posterior <- bastaOut[[1]]$post[thinned]
  if (dataObj$updB | dataObj$updD) {
    posterior <- posterior + bastaOut[[1]]$postXu[thinned]
  }
  if (algObj$minAge > 0) 
    lams <- bastaOut[[1]]$lambda[thinned]
  if (algObj$nsim > 1) {
    for (ii in 2:algObj$nsim) {
      parMat <- rbind(parMat, bastaOut[[ii]]$par[thinned, 
      ])
      fullParMat <- rbind(fullParMat, bastaOut[[ii]]$par[algObj$burnin:algObj$niter, 
      ])
      postii <- bastaOut[[ii]]$post[thinned]
      if (dataObj$updB | dataObj$updD) {
        postii <- postii + bastaOut[[ii]]$postXu[thinned]
      }
      posterior <- c(posterior, postii)
      if (algObj$minAge > 0) 
        lams <- c(lams, bastaOut[[ii]]$lambda[thinned])
    }
  }
  coef <- cbind(apply(parMat, 2, mean, na.rm = TRUE), apply(parMat, 
                                                            2, sd, na.rm = TRUE), t(apply(parMat, 2, quantile, c(0.025, 
                                                                                                                 0.975), na.rm = TRUE)), apply(parMat, 2, function(x) cor(x[-1], 
                                                                                                                                                                          x[-length(x)], use = "complete.obs")), apply(fullParMat, 
                                                                                                                                                                                                                       2, function(x) length(which(diff(x) != 0))/(length(x) - 
                                                                                                                                                                                                                                                                     1)))
  colnames(coef) <- c("Estimate", "StdErr", "Lower95%CI", "Upper95%CI", 
                      "SerAutocor", "UpdateRate")
  if (algObj$nsim > 1) {
    idSims <- rep(1:algObj$nsim, each = nthin)
    Means <- apply(parMat, 2, function(x) tapply(x, idSims, 
                                                 mean))
    Vars <- apply(parMat, 2, function(x) tapply(x, idSims, 
                                                var))
    meanall <- apply(Means, 2, mean)
    B <- nthin/(algObj$nsim - 1) * apply(t((t(Means) - meanall)^2), 
                                         2, sum)
    W <- 1/algObj$nsim * apply(Vars, 2, sum)
    Varpl <- (nthin - 1)/nthin * W + 1/nthin * B
    Rhat <- sqrt(Varpl/W)
    Rhat[Varpl == 0] <- 1
    conv <- cbind(B, W, Varpl, Rhat)
    rownames(conv) <- colnames(parMat)
    coef <- cbind(coef, conv[, "Rhat"])
    colnames(coef) <- c(colnames(coef)[-ncol(coef)], "PotScaleReduc")
    idnconv <- which(conv[, "Rhat"] > 1.1)
    if (length(idnconv) == 0) {
      Dave <- -2 * mean(posterior)
      parsMode <- parsIni
      parsMode$theta[1:fullParObj$theta$len] <- coef[1:fullParObj$theta$len, 
                                                     1]
      if (class(parsIni)[1] == "theGam") {
        idGam <- which(substr(fullParObj$allNames, 1, 
                              2) == "ga")
        parsMode$gamma <- coef[idGam, 1]
      }
      if (algObj$minAge > 0) 
        parsMode$lambda <- mean(lams)
      idPi <- which(substr(fullParObj$allNames, 1, 2) == 
                      "pi")
      parsMode$pi <- coef[idPi, 1]
      parsCovMode <- .CalcParCovObj(covObj, parsMode, parsCovIni)
      agesMode <- agesIni
      if (dataObj$updB | dataObj$updD) {
        cat("\nCalculating DIC...")
        nbd <- nrow(bastaOut[[1]]$birth)
        DmodeVec <- rep(0, nbd)
        for (dm in 1:nbd) {
          DmodeMean <- 0
          for (si in 1:algObj$nsim) {
            bi <- dataObj$bi
            di <- dataObj$di
            if (dataObj$updB) {
              bi[dataObj$idNoB] <- bastaOut[[si]]$birth[dm, 
              ]
            }
            if (dataObj$updD) {
              di[dataObj$idNoD] <- bastaOut[[si]]$death[dm, 
              ]
            }
            agesMode$ages[, "birth"] <- bi
            agesMode$ages[, "death"] <- di
            agesMode$ages[, "age"] <- agesMode$ages[, 
                                                    "death"] - agesMode$ages[, "birth"]
            if (algObj$minAge > 0) {
              agesMode$ages <- .SplitByMinAge(agesMode$ages, 
                                              algObj)
            }
            else {
              idtr <- which(agesMode$ages[, "birth"] < 
                              algObj$start)
              agesMode$ages[idtr, "ageTr"] <- algObj$start - 
                agesMode$ages[idtr, "birth"]
            }
            firstAlive <- c(apply(cbind(algObj$start, 
                                        agesMode$ages[, "birth"] + 1), 1, max))
            lastAlive <- c(apply(cbind(algObj$end, agesMode$ages[, 
                                                                 "death"]), 1, min))
            alive <- .BuildAliveMatrix(firstAlive, lastAlive, 
                                       dataObj)
            agesMode$alive <- alive
            postMode <- .CalcPostX(agesMode, parsMode, 
                                   postIni, parsCovMode, 1:dataObj$n, CalcSurv, 
                                   priorAgeObj, fullParObj, covObj, dataObj)
            DmodeMean <- DmodeMean + sum(postMode$mat[, 
                                                      "fx"] - postMode$mat[, "Sx"] + postMode$mat[, 
                                                                                                  "px"] + postMode$mat[, "lx"])
            if (dataObj$updB | dataObj$updD) {
              DmodeMean <- DmodeMean + sum(postMode$mat[dataObj$idNoA, 
                                                        "postX"])
            }
          }
          DmodeVec[dm] <- DmodeMean/algObj$nsim
        }
        Dmode <- -2 * mean(DmodeVec)
        cat(" done\n")
      }
      else {
        agesMode <- agesIni
        postMode <- .CalcPostX(agesMode, parsMode, postIni, 
                               parsCovMode, 1:dataObj$n, CalcSurv, priorAgeObj, 
                               fullParObj, covObj, dataObj)
        Dmode <- -2 * (sum(postMode$mat[, "fx"] - postMode$mat[, 
                                                               "Sx"] + postMode$mat[, "px"] + postMode$mat[, 
                                                                                                           "lx"]))
      }
      DIC <- 2 * Dave - Dmode
      pD <- Dave - Dmode
      runOld <- FALSE
      if (runOld) {
        L <- length(posterior)
        Dm <- -2 * posterior
        densD <- density(posterior)
        Dave <- mean(Dm)
        pD <- 1/2 * 1/(L - 1) * sum((Dm - Dave)^2)
        DIC <- Dave + pD
        Dmode <- abs(2 * Dave - DIC)
      }
      k <- ncol(parMat)
      modSel <- c(Dave, Dmode, pD, k, DIC)
      names(modSel) <- c("D.ave", "D.mode", "pD", "k", 
                         "DIC")
      cat("Survival parameters converged appropriately.", 
          "\nDIC was calculated.\n")
      kulLeib <- .CalcKulbackLeibler(coef, covObj, defTheta, 
                                     fullParObj, algObj, dataObj)
    }
    else {
      warning("Convergence not reached for some survival parameters.", 
              "\nDIC could not be calculated.\n", call. = FALSE)
      modSel <- "Not calculated"
      kulLeib <- "Not calculated"
    }
  }
  else {
    conv <- "Not calculated"
    modSel <- "Not calculated"
    kulLeib <- "Not calculated"
  }
  diagObj <- list(coefficients = coef, DIC = modSel, convergence = conv, 
                  KullbackLeibler = kulLeib, params = parMat)
  return(diagObj)
}


# getAnywhere(.CreateAlgObj)
.CreateAlgObj <- function (model, shape, studyStart, studyEnd, minAge, covarsStruct, 
                           recaptTrans, niter, burnin, thinning, updateJumps, nsim) 
{
  return(list(model = model, shape = shape, start = studyStart, 
              end = studyEnd, minAge = minAge, covStruc = covarsStruct, 
              recap = recaptTrans, niter = niter, burnin = burnin, 
              thinning = thinning, updJump = updateJumps, nsim = nsim))
}

# getAnywhere(.FindErrors)
.FindErrors <- function (object, algObj) 
{
  data.check <- DataCheck(object, algObj$start, algObj$end, 
                          silent = TRUE)
  if (!data.check[[1]]) {
    stop("You have an error in Dataframe 'object',\nplease use function ", 
         "'DataCheck'\n", call. = FALSE)
  }
  if (algObj$burnin > algObj$niter) {
    stop("Object 'algObj$burnin' larger than 'algObj$niter'.", 
         call. = FALSE)
  }
  if (algObj$thinning > algObj$niter) {
    stop("Object 'algObj$thinning' larger than 'algObj$niter'.", 
         call. = FALSE)
  }
  if (!is.element(algObj$model, c("EX", "GO", "WE", "LO"))) {
    stop("Model misspecification: specify available models", 
         " (i.e. 'EX', 'GO', 'WE' or 'LO')\n", call. = FALSE)
  }
  if (!is.element(algObj$shape, c("simple", "Makeham", "bathtub"))) {
    stop("algObj$shape misspecification. Appropriate arguments are:", 
         " 'simple', 'Makeham' or 'bathtub'.\n", call. = FALSE)
  }
  if (!is.element(algObj$covStruc, c("fused", "prop.haz", "all.in.mort"))) {
    stop("Covariate structure misspecification. Appropriate arguments are:", 
         " 'fused', 'prop.haz' or 'all.in.mort'.\n", call. = FALSE)
  }
  if (algObj$model == "EX" & algObj$shape != "simple") {
    stop("Model misspecification: EX algObj$model can only be fitted with a", 
         " simple algObj$shape", call. = FALSE)
  }
  if (algObj$model == "EX" & algObj$covStruc != "fused") {
    stop("Model misspecification: EX algObj$model can only be fitted with a", 
         " fused covariate structure", call. = FALSE)
  }
  if (algObj$covStruc == "all.in.mort" & sum(algObj$model == 
                                             "GO", algObj$shape == "simple") < 2) {
    stop("Model misspecification: all.in.mort is only available with", 
         " Gompertz (GO) models and simple algObj$shape.", 
         call. = FALSE)
  }
}

# getAnywhere(.PrepDataObj)
.PrepDataObj <- function (object, algObj) 
{
  dataObj <- list()
  dataObj$study <- algObj$start:algObj$end
  dataObj$studyLen <- length(dataObj$study)
  dataObj$n <- nrow(object)
  dataObj$Y <- as.matrix(object[, 1:dataObj$studyLen + 3])
  colnames(dataObj$Y) <- dataObj$study
  bd <- as.matrix(object[, 2:3])
  dataObj$bi <- bd[, 1]
  dataObj$di <- bd[, 2]
  bi0 <- which(dataObj$bi == 0)
  if (length(bi0) > 0) {
    dataObj$idNoB <- bi0
    dataObj$updB <- TRUE
  }
  else {
    dataObj$updB <- FALSE
  }
  di0 <- which(dataObj$di == 0)
  if (length(di0) > 0) {
    dataObj$idNoD <- di0
    dataObj$updD <- TRUE
  }
  else {
    dataObj$updD <- FALSE
  }
  if (!dataObj$updB & !dataObj$updD) {
    class(dataObj) <- "noAgeUpd"
  }
  else {
    dataObj$idNoA <- sort(unique(c(dataObj$idNoB, dataObj$idNoD)))
    ytemp <- t(t(dataObj$Y) * dataObj$study)
    dataObj$lastObs <- c(apply(ytemp, 1, max))
    ytemp[ytemp == 0] <- 10000
    dataObj$firstObs <- c(apply(ytemp, 1, min))
    dataObj$firstObs[dataObj$firstObs == 10000] <- 0
    dataObj$oi <- dataObj$Y %*% rep(1, dataObj$studyLen)
    dataObj$Tm <- matrix(dataObj$study, dataObj$n, dataObj$studyLen, 
                         byrow = TRUE)
    fii <- dataObj$firstObs
    id1 <- which(dataObj$bi > 0 & dataObj$bi >= algObj$start)
    fii[id1] <- dataObj$bi[id1] + 1
    fii[dataObj$bi > 0 & dataObj$bi < algObj$start] <- algObj$start
    lii <- dataObj$lastObs
    id2 <- which(dataObj$di > 0 & dataObj$di <= algObj$end)
    lii[id2] <- dataObj$di[id2] - 1
    lii[dataObj$di > 0 & dataObj$di > algObj$end] <- algObj$end
    dataObj$obsMat <- .BuildAliveMatrix(fii, lii, dataObj)
    dataObj$obsMat[lii == 0 | fii == 0, ] <- 0
    class(dataObj) <- "ageUpd"
  }
  dataObj$Dx <- 1
  return(dataObj)
}

# getAnywhere(.SetDefaultTheta)
.SetDefaultTheta <- function (algObj) 
{
  if (algObj$model == "EX") {
    nTh <- 1
    startTh <- 0.2
    jumpTh <- 0.1
    priorMean <- 0.06
    priorSd <- 1
    nameTh <- "b0"
    lowTh <- 0
    jitter <- 0.5
  }
  else if (algObj$model == "GO") {
    nTh <- 2
    startTh <- c(-2, 0.01)
    jumpTh <- c(0.1, 0.1)
    priorMean <- c(-3, 0.01)
    priorSd <- c(1, 1)
    nameTh <- c("b0", "b1")
    lowTh <- c(-Inf, -Inf)
    jitter <- c(0.5, 0.2)
    if (algObj$shape == "bathtub") {
      lowTh <- c(-Inf, 0)
    }
  }
  else if (algObj$model == "WE") {
    nTh <- 2
    startTh <- c(1.5, 0.2)
    jumpTh <- c(0.01, 0.1)
    priorMean <- c(1.5, 0.05)
    priorSd <- c(1, 1)
    nameTh <- c("b0", "b1")
    lowTh <- c(0, 0)
    jitter <- c(0.5, 0.2)
  }
  else if (algObj$model == "LO") {
    nTh <- 3
    startTh <- c(-2, 0.01, 1e-04)
    jumpTh <- c(0.1, 0.1, 0.1)
    priorMean <- c(-3, 0.01, 1e-10)
    priorSd <- c(1, 1, 1)
    nameTh <- c("b0", "b1", "b2")
    lowTh <- c(-Inf, 0, 0)
    jitter <- c(0.5, 0.2, 0.5)
  }
  if (algObj$shape == "Makeham") {
    nTh <- nTh + 1
    startTh <- c(0, startTh)
    jumpTh <- c(0.1, jumpTh)
    priorMean <- c(0, priorMean)
    priorSd <- c(1, priorSd)
    nameTh <- c("c", nameTh)
    lowTh <- c(0, lowTh)
    jitter <- c(0.25, jitter)
  }
  else if (algObj$shape == "bathtub") {
    nTh <- nTh + 3
    startTh <- c(-0.1, 0.6, 0, startTh)
    jumpTh <- c(0.1, 0.1, 0.1, jumpTh)
    priorMean <- c(-2, 0.01, 0, priorMean)
    priorSd <- c(1, 1, 1, priorSd)
    nameTh <- c("a0", "a1", "c", nameTh)
    lowTh <- c(-Inf, 0, 0, lowTh)
    jitter <- c(0.5, 0.2, 0.2, jitter)
  }
  defaultTheta <- list(length = nTh, start = startTh, jump = jumpTh, 
                       priorMean = priorMean, priorSd = priorSd, name = nameTh, 
                       low = lowTh, jitter = jitter)
  attr(defaultTheta, "algObj$model") = algObj$model
  attr(defaultTheta, "algObj$shape") = algObj$shape
  return(defaultTheta)
}

# getAnywhere(.DefineMort)
.DefineMort <- function (algObj) 
{
  if (algObj$model == "EX") {
    CalcMort <- function(x, theta) c(theta) * rep(1, length(x))
  }
  else if (algObj$model == "GO") {
    if (algObj$shape == "simple") {
      CalcMort <- function(x, theta) {
        exp(theta[, "b0"] + theta[, "b1"] * x)
      }
    }
    else if (algObj$shape == "Makeham") {
      CalcMort <- function(x, theta) {
        theta[, "c"] + exp(theta[, "b0"] + theta[, "b1"] * 
                             x)
      }
    }
    else {
      CalcMort <- function(x, theta) {
        exp(theta[, "a0"] - theta[, "a1"] * x) + theta[, 
                                                       "c"] + exp(theta[, "b0"] + theta[, "b1"] * 
                                                                    x)
      }
    }
  }
  else if (algObj$model == "WE") {
    if (algObj$shape == "simple") {
      CalcMort <- function(x, theta) {
        theta[, "b0"] * theta[, "b1"]^theta[, "b0"] * 
          x^(theta[, "b0"] - 1)
      }
    }
    else if (algObj$shape == "Makeham") {
      CalcMort <- function(x, theta) {
        theta[, "c"] + theta[, "b0"] * theta[, "b1"]^theta[, 
                                                           "b0"] * x^(theta[, "b0"] - 1)
      }
    }
    else {
      CalcMort <- function(x, theta) {
        exp(theta[, "a0"] - theta[, "a1"] * x) + theta[, 
                                                       "c"] + theta[, "b0"] * theta[, "b1"]^theta[, 
                                                                                                  "b0"] * x^(theta[, "b0"] - 1)
      }
    }
  }
  else if (algObj$model == "LO") {
    if (algObj$shape == "simple") {
      CalcMort <- function(x, theta) {
        exp(theta[, "b0"] + theta[, "b1"] * x)/(1 + theta[, 
                                                          "b2"] * exp(theta[, "b0"])/theta[, "b1"] * 
                                                  (exp(theta[, "b1"] * x) - 1))
      }
    }
    else if (algObj$shape == "Makeham") {
      CalcMort <- function(x, theta) {
        theta[, "c"] + exp(theta[, "b0"] + theta[, "b1"] * 
                             x)/(1 + theta[, "b2"] * exp(theta[, "b0"])/theta[, 
                                                                              "b1"] * (exp(theta[, "b1"] * x) - 1))
      }
    }
    else {
      CalcMort <- function(x, theta) {
        exp(theta[, "a0"] - theta[, "a1"] * x) + theta[, 
                                                       "c"] + exp(theta[, "b0"] + theta[, "b1"] * 
                                                                    x)/(1 + theta[, "b2"] * exp(theta[, "b0"])/theta[, 
                                                                                                                     "b1"] * (exp(theta[, "b1"] * x) - 1))
      }
    }
  }
  return(CalcMort)
}

# getAnywhere(.DefineSurv)
.DefineSurv <- function (algObj) 
{
  if (algObj$model == "EX") {
    CalcSurv <- function(x, theta) exp(-c(theta) * x)
  }
  else if (algObj$model == "GO") {
    if (algObj$shape == "simple") {
      CalcSurv <- function(x, theta) {
        exp(exp(theta[, "b0"])/theta[, "b1"] * (1 - exp(theta[, 
                                                              "b1"] * x)))
      }
    }
    else if (algObj$shape == "Makeham") {
      CalcSurv <- function(x, theta) {
        exp(-theta[, "c"] * x + exp(theta[, "b0"])/theta[, 
                                                         "b1"] * (1 - exp(theta[, "b1"] * x)))
      }
    }
    else {
      CalcSurv <- function(x, theta) {
        exp(exp(theta[, "a0"])/theta[, "a1"] * (exp(-theta[, 
                                                           "a1"] * x) - 1) - theta[, "c"] * x + exp(theta[, 
                                                                                                          "b0"])/theta[, "b1"] * (1 - exp(theta[, "b1"] * 
                                                                                                                                            x)))
      }
    }
  }
  else if (algObj$model == "WE") {
    if (algObj$shape == "simple") {
      CalcSurv <- function(x, theta) {
        exp(-(theta[, "b1"] * x)^theta[, "b0"])
      }
    }
    else if (algObj$shape == "Makeham") {
      CalcSurv <- function(x, theta) {
        exp(-theta[, "c"] * x - (theta[, "b1"] * x)^theta[, 
                                                          "b0"])
      }
    }
    else {
      CalcSurv <- function(x, theta) {
        exp(exp(theta[, "a0"])/theta[, "a1"] * (exp(-theta[, 
                                                           "a1"] * x) - 1) - theta[, "c"] * x - (theta[, 
                                                                                                       "b1"] * x)^theta[, "b0"])
      }
    }
  }
  else if (algObj$model == "LO") {
    if (algObj$shape == "simple") {
      CalcSurv <- function(x, theta) {
        (1 + theta[, "b2"] * exp(theta[, "b0"])/theta[, 
                                                      "b1"] * (exp(theta[, "b1"] * x) - 1))^(-1/theta[, 
                                                                                                      "b2"])
      }
    }
    else if (algObj$shape == "Makeham") {
      CalcSurv <- function(x, theta) {
        exp(-theta[, "c"] * x) * (1 + theta[, "b2"] * 
                                    exp(theta[, "b0"])/theta[, "b1"] * (exp(theta[, 
                                                                                  "b1"] * x) - 1))^(-1/theta[, "b2"])
      }
    }
    else {
      CalcSurv <- function(x, theta) {
        exp(exp(theta[, "a0"])/theta[, "a1"] * (exp(-theta[, 
                                                           "a1"] * x) - 1) - theta[, "c"] * x) * (1 + 
                                                                                                    theta[, "b2"] * exp(theta[, "b0"])/theta[, 
                                                                                                                                             "b1"] * (exp(theta[, "b1"] * x) - 1))^(-1/theta[, 
                                                                                                                                                                                             "b2"])
      }
    }
  }
  return(CalcSurv)
}

# getAnywhere(.CreateCovObj)
.CreateCovObj <- function (object, dataObj, algObj) 
{
  covObj <- list()
  covClass <- c("noCov", "noCovType")
  if (ncol(object) > dataObj$studyLen + 3) {
    covMat <- as.matrix(object[, (dataObj$studyLen + 4):ncol(object)])
    mode(covMat) <- "numeric"
    colnames(covMat) <- colnames(object)[(dataObj$studyLen + 
                                            4):ncol(object)]
    covType <- .FindCovType(covMat)
    if (algObj$covStruc == "fused") {
      covClass[1] <- "fused"
      if (!is.null(covType$cat)) {
        covObj$inMort <- covMat[, covType$cat]
        covObj$imLen <- ncol(covObj$inMort)
      }
      else {
        covClass[1] <- "propHaz"
      }
      if (!is.null(covType$cont)) {
        covObj$propHaz <- matrix(covMat[, c(covType$int, 
                                            covType$cont)], ncol = length(c(covType$int, 
                                                                            covType$cont)), dimnames = list(NULL, c(names(covType$int), 
                                                                                                                    names(covType$cont))))
        covObj$phLen <- ncol(covObj$propHaz)
      }
      else {
        covClass[1] <- "inMort"
      }
    }
    else if (algObj$covStruc == "all.in.mort") {
      if (is.null(covType$int) & is.null(covType$cat)) {
        covObj$inMort <- cbind(1, covMat)
        colnames(covObj$inMort) <- c("Intercept", colnames(covMat))
      }
      else {
        covObj$inMort <- covMat
      }
      covObj$imLen <- ncol(covObj$inMort)
      covClass[1] <- "inMort"
    }
    else {
      if (!is.null(covType$int)) {
        covObj$propHaz <- matrix(covMat[, -covType$int], 
                                 dataObj$n, ncol(covMat) - 1, dimnames = list(NULL, 
                                                                              colnames(covMat)[-covType$int]))
      }
      else if (!is.null(covType$cat)) {
        covObj$propHaz <- matrix(covMat[, -covType$cat[1]], 
                                 dataObj$n, ncol(covMat) - 1, dimnames = list(NULL, 
                                                                              colnames(covMat)[-covType$cat[1]]))
      }
      else {
        covObj$propHaz <- covMat
      }
      covObj$phLen <- ncol(covObj$propHaz)
      covClass[1] <- "propHaz"
    }
    if (!is.null(covType$cat) & !is.null(covType$cont)) {
      covClass[2] <- "bothCov"
      covObj$cat <- covType$cat
      covObj$cont <- covType$cont
    }
    else if (!is.null(covType$cat)) {
      covClass[2] <- "cateCov"
      covObj$cat <- covType$cat
    }
    else if (!is.null(covType$cont)) {
      covClass[2] <- "contCov"
      covObj$cont <- covType$cont
    }
  }
  else {
    covObj$covs <- NULL
  }
  class(covObj) <- covClass
  return(covObj)
}

# getAnywhere(.CreateUserPar)
.CreateUserPar <- function (covObj, argList) 
{
  userPars <- list()
  genParName <- c("theta", "gamma")
  parTypes <- c("Start", "Jumps", "PriorMean", "PriorSd")
  parTypesList <- c("start", "jump", "priorMean", "priorSd")
  for (genPp in 1:2) {
    if (all(genPp == 2 & class(covObj)[1] %in% c("fused", 
                                                 "propHaz")) | genPp == 1) {
      userPars[[genParName[genPp]]] <- list()
      for (pp in 1:4) {
        usrPar <- sprintf("%s%s", genParName[genPp], 
                          parTypes[pp])
        if (usrPar %in% names(argList)) {
          userPars[[genParName[genPp]]][[parTypesList[pp]]] <- argList[[usrPar]]
        }
        else {
          userPars[[genParName[genPp]]][[parTypesList[pp]]] <- NULL
        }
      }
    }
  }
  return(userPars)
}

# getAnywhere(.BuildFullParObj)
.BuildFullParObj <- function (covObj, defTheta, algObj, userPars, dataObj) 
{
  fullParObj <- list()
  fullParObj$theta <- list()
  parNames <- c("start", "priorMean", "priorSd", "jump")
  for (i in 1:4) {
    if (class(covObj)[1] %in% c("inMort", "fused")) {
      if (is.null(userPars$theta[[parNames[i]]])) {
        thetaMat <- matrix(defTheta[[parNames[i]]], covObj$imLen, 
                           defTheta$length, byrow = TRUE, dimnames = list(colnames(covObj$inMort), 
                                                                          defTheta$name))
        if (i %in% c(1, 2)) {
          if (class(covObj)[1] == "inMort" & class(covObj)[2] %in% 
              c("contCov", "bothCov")) {
            thetaMat[names(covObj$cont), ] <- 0
          }
        }
        fullParObj$theta[[parNames[[i]]]] <- thetaMat
      }
      else {
        if (is.element(length(userPars$theta[[parNames[i]]]), 
                       c(defTheta$length, defTheta$length * covObj$imLen))) {
          if (length(userPars$theta[[parNames[i]]]) == 
              defTheta$length) {
            fullParObj$theta[[parNames[[i]]]] <- matrix(userPars$theta[[parNames[i]]], 
                                                        covObj$imLen, defTheta$length, byrow = TRUE, 
                                                        dimnames = list(colnames(covObj$inMort), 
                                                                        defTheta$name))
          }
          else {
            fullParObj$theta[[parNames[[i]]]] <- userPars$theta[[parNames[i]]]
            dimnames(fullParObj$theta[[parNames[[i]]]]) <- list(colnames(covObj$inMort), 
                                                                defTheta$name)
          }
        }
        else {
          stop(paste("\nDimensions of theta ", parNames[i], 
                     " matrix are incorrect.\n", "Provide a single vector of length ", 
                     defTheta$length, "\nor a matrix of dimensions ", 
                     covObj$imLen, " times ", defTheta$length, 
                     ".\n(i.e. number of covariates times number", 
                     " of\n parameters for model ", algObj$model, 
                     " with ", algObj$shape, " shape).", sep = ""), 
               call. = FALSE)
        }
      }
      allParNames <- paste(rep(defTheta$name, each = ncol(covObj$inMort)), 
                           rep(colnames(covObj$inMort), defTheta$len), sep = ".")
    }
    else {
      if (is.null(userPars$theta[[parNames[i]]])) {
        fullParObj$theta[[parNames[[i]]]] <- matrix(defTheta[[parNames[i]]], 
                                                    1, defTheta$length, dimnames = list(NULL, defTheta$name))
      }
      else {
        if (length(userPars$theta[[parNames[i]]]) == 
            defTheta$length) {
          fullParObj$theta[[parNames[[i]]]] <- matrix(userPars$theta[[parNames[i]]], 
                                                      1, defTheta$length, dimnames = list(NULL, 
                                                                                          defTheta$name))
        }
        else {
          stop(paste("\nLength of theta ", parNames[i], 
                     " is incorrect.\n", "Provide a single vector of length ", 
                     defTheta$length, ".\n(i.e. number of parameters for model ", 
                     algObj$model, " with ", algObj$shape, " shape).", 
                     sep = ""), call. = FALSE)
        }
      }
      allParNames <- defTheta$name
    }
  }
  fullParObj$theta$low <- t(t(fullParObj$theta$start) * 0 + 
                              defTheta$low)
  fullParObj$theta$len <- length(fullParObj$theta$start)
  if (class(covObj)[1] %in% c("propHaz", "fused")) {
    fullParObj$gamma <- list()
    for (i in 1:4) {
      if (is.null(userPars$gamma[[parNames[i]]])) {
        fullParObj$gamma[[parNames[i]]] <- rep(c(0.01, 
                                                 0, 1, 0.1)[i], covObj$phLen)
        names(fullParObj$gamma[[parNames[i]]]) <- colnames(covObj$propHaz)
      }
      else {
        if (length(userPars$gamma[[parNames[i]]]) == 
            covObj$phLen) {
          fullParObj$gamma[[parNames[i]]] <- userPars$gamma[[parNames[i]]]
          names(fullParObj$gamma[[parNames[i]]]) <- colnames(covObj$propHaz)
        }
        else {
          stop(paste("\nLength of gamma parameters is incorrect.\n", 
                     "Provide a single vector of length ", covObj$phLen, 
                     ".\n(i.e. number of proportional hazards covariates).", 
                     sep = ""), call. = FALSE)
        }
      }
    }
    fullParObj$gamma$len <- length(fullParObj$gamma$start)
    allParNames <- c(allParNames, paste("gamma", colnames(covObj$propHaz), 
                                        sep = "."))
  }
  Classes <- ifelse(class(covObj)[1] %in% c("inMort", "noCov"), 
                    "theta", "theGam")
  if (algObj$minAge > 0) {
    fullParObj$lambda <- list(start = 0.01, priorMean = 0.01, 
                              priorSd = 1, jump = 0.01)
    Classes <- c(Classes, "lambda")
  }
  else {
    Classes <- c(Classes, "noLambda")
  }
  if (class(dataObj) == "ageUpd") {
    fullParObj$pi <- list()
    study <- algObj$start:algObj$end
    if (length(algObj$recap) == length(study)) {
      if (all(algObj$recap %in% study)) {
        idpi <- rep(1, length(study))
        for (i in 2:length(idpi)) {
          idpi[i] <- ifelse(algObj$recap[i] == algObj$recap[i - 
                                                              1], idpi[i - 1], ifelse(algObj$recap[i] %in% 
                                                                                        algObj$recap[1:(i - 1)], idpi[which(algObj$recap[1:(i - 
                                                                                                                                              1)] == algObj$recap[i])[1]], max(idpi[1:(i - 
                                                                                                                                                                                         1)] + 1)))
        }
      }
      else if (all(algObj$recap %in% 1:length(study))) {
        idpi <- algObj$recap
      }
      namespi <- unique(algObj$recap)
    }
    else {
      idpi <- findInterval(study, algObj$recap)
      namespi <- algObj$recap
    }
    names(idpi) <- study
    npi <- length(unique(idpi))
    fullParObj$pi$start <- rep(0.5, npi)
    names(fullParObj$pi$start) <- namespi
    fullParObj$pi$idpi <- idpi
    fullParObj$pi$n <- npi
    fullParObj$pi$prior2 <- 1
    fullParObj$pi$Prior1 <- tapply(1 + t(t(dataObj$Y) %*% 
                                           rep(1, dataObj$n)), idpi, sum)
    fullParObj$pi$len <- length(fullParObj$pi$start)
    fullParObj$allNames <- c(allParNames, paste("pi", namespi, 
                                                sep = "."))
    Classes <- c(Classes, "pi", "noEta")
  }
  else {
    fullParObj$allNames <- allParNames
    Classes <- c(Classes, "noPi", "noEta")
  }
  fullParObj$class <- Classes
  return(fullParObj)
}

# getAnywhere(.PrepAgeObj)
.PrepAgeObj <- function (dataObj, algObj) 
{
  ageObj <- list()
  birth <- dataObj$bi
  if (dataObj$updB) {
    idBi0Fi1 <- which(dataObj$bi == 0 & dataObj$firstObs > 
                        0)
    birth[idBi0Fi1] <- dataObj$firstObs[idBi0Fi1] - sample(1:6, 
                                                           length(idBi0Fi1), replace = TRUE)
    idBi0Fi0 <- which(dataObj$bi == 0 & dataObj$firstObs == 
                        0 & dataObj$di > 0)
    birth[idBi0Fi0] <- dataObj$di[idBi0Fi0] - sample(0:6, 
                                                     length(idBi0Fi0), replace = TRUE)
  }
  death <- dataObj$di
  if (dataObj$updD) {
    idDi0Li1 <- which(dataObj$di == 0 & dataObj$lastObs > 
                        0)
    death[idDi0Li1] <- dataObj$lastObs[idDi0Li1] + sample(0:6, 
                                                          length(idDi0Li1), replace = TRUE)
    idDi0Li0 <- which(dataObj$di == 0 & dataObj$lastObs == 
                        0)
    death[idDi0Li0] <- dataObj$bi[idDi0Li0] + sample(0:6, 
                                                     length(idDi0Li0), replace = TRUE)
    idDiNeg <- which(death < algObj$start)
    death[idDiNeg] <- algObj$start + sample(1:6, length(idDiNeg), 
                                            replace = TRUE)
  }
  age <- death - birth
  ageObj$ages <- cbind(birth, death, age)
  firstObs <- c(apply(cbind(algObj$start, birth + 1), 1, max))
  lastObs <- c(apply(cbind(algObj$end, death), 1, min))
  alive <- .BuildAliveMatrix(firstObs, lastObs, dataObj)
  ageObj$alive <- alive
  class(ageObj) <- c(class(dataObj), "noMinAge")
  if (algObj$minAge > 0) {
    indAd <- rep(0, dataObj$n)
    indJu <- indAd
    indAd[age >= algObj$minAge] <- 1
    indJu[age < algObj$minAge] <- 1
    ageJu <- age
    ageJu[age > algObj$minAge] <- algObj$minAge
    ageAd <- age - algObj$minAge
    ageAd[age < algObj$minAge] <- 0
    ageJuTr <- age * 0
    idtr <- which(birth < algObj$start & algObj$start - birth < 
                    algObj$minAge)
    ageJuTr[idtr] <- algObj$start - birth[idtr]
    ageAdTr <- age * 0
    idtr <- which(birth + algObj$minAge < algObj$start)
    ageAdTr[idtr] <- algObj$start - (birth[idtr] + algObj$minAge)
    ageObj$ages <- cbind(ageObj$ages, ageAd, ageAdTr, indAd, 
                         ageJu, ageJuTr, indJu)
    class(ageObj)[2] <- "minAge"
  }
  else {
    idtr <- which(birth < algObj$start)
    ageTr <- age * 0
    ageTr[idtr] <- algObj$start - birth[idtr]
    ageObj$ages <- cbind(ageObj$ages, ageTr)
  }
  return(ageObj)
}

# getAnywhere(.DefineIniParObj)
.DefineIniParObj <- function (fullParObj) 
{
  iniParObj <- list()
  iniParObj$theta <- fullParObj$theta$start
  if (fullParObj$class[1] %in% c("theGam")) {
    iniParObj$gamma <- fullParObj$gamma$start
  }
  if (fullParObj$class[2] == "lambda") {
    iniParObj$lambda <- fullParObj$lambda$start
  }
  if (fullParObj$class[3] %in% c("pi", "piEta")) {
    iniParObj$pi <- fullParObj$pi$start
  }
  class(iniParObj) <- fullParObj$class
  return(iniParObj)
}

# getAnywhere(.BuildParCovObj)
# print(methods('.BuildParCovObj'))
.BuildParCovObj <- function (covObj, ...) 
  UseMethod(".BuildParCovObj")

# getAnywhere(.SetPriorAgeDist)
.SetPriorAgeDist <- function (fullParObj, CalcSurv, dataObj, covObj, parsIni, parsCovIni) 
{
  dxx <- 0.001
  xx <- seq(0, 100, dxx)
  parsPrior <- list()
  parsPrior$theta <- fullParObj$theta$priorMean
  if (class(parsIni)[1] == "theGam") {
    parsPrior$gamma <- fullParObj$gamma$priorMean
  }
  parsPriorCov <- list()
  if (class(covObj)[1] == "noCov") {
    Ex <- sum(CalcSurv(xx, parsPrior$theta) * dxx)
    lifeExp <- rep(Ex, dataObj$n)
  }
  else if (class(covObj)[1] == "inMort") {
    if (class(covObj)[2] %in% c("bothCov", "contCov")) {
      meanCont <- apply(matrix(covObj$inMort[, names(covObj$cont)], 
                               ncol = length(covObj$cont)), 2, mean)
      thetaCont <- meanCont %*% matrix(parsPrior$theta[names(covObj$cont), 
      ], nrow = length(covObj$cont))
      colnames(thetaCont) <- colnames(fullParObj$theta$priorMean)
    }
    else {
      thetaCont <- t(as.matrix(parsPrior$theta[1, ] * 0))
    }
    Ex <- sapply(1:(ncol(covObj$inMort) - length(covObj$cont)), 
                 function(pp) sum(CalcSurv(xx, parsPrior$theta[pp, 
                 ] + thetaCont) * dxx))
    if (is.null(covObj$cat)) {
      lifeExp <- rep(Ex, dataObj$n)
    }
    else {
      names(Ex) <- names(covObj$cat)
      lifeExp <- covObj$inMor[, names(covObj$cat)] %*% 
        Ex
    }
  }
  else if (class(covObj)[1] == "fused") {
    meanCont <- apply(covObj$propHaz, 2, mean)
    gam <- sum(parsPrior$gamma * meanCont)
    Ex <- sapply(1:length(covObj$cat), function(pp) sum(CalcSurv(xx, 
                                                                 t(parsPrior$theta[pp, ]))^exp(gam) * dxx))
    names(Ex) <- names(covObj$cat)
    lifeExp <- covObj$inMor %*% Ex
  }
  else {
    if (class(covObj)[2] == "bothCov") {
      meanCont <- apply(matrix(covObj$propHaz[, names(covObj$cont)], 
                               ncol = length(covObj$cont)), 2, mean)
      gam <- c(0, parsPrior$gamma[names(covObj$cat)[-1]]) + 
        sum(parsPrior$gamma * meanCont)
      Ex <- sapply(1:(length(gam)), function(pp) sum(CalcSurv(xx, 
                                                              parsPrior$theta)^exp(gam[pp]) * dxx))
      names(Ex) <- names(covObj$cat)
      if (covObj$phLen - length(covObj$cont) == 1) {
        intercept <- 1 - covObj$propHaz[, names(covObj$cat)[-1]]
      }
      else {
        intercept <- 1 - apply(covObj$propHaz[, names(covObj$cat)[-1]], 
                               1, sum)
      }
      lifeExp <- cbind(intercept, covObj$propHaz[, names(covObj$cat)[-1]]) %*% 
        Ex
    }
    else if (class(covObj)[2] == "cateCov") {
      gam <- c(0, parsPrior$gamma[names(covObj$cat)[-1]])
      Ex <- sapply(1:(length(gam)), function(pp) sum(CalcSurv(xx, 
                                                              parsPrior$theta)^exp(gam[pp]) * dxx))
      names(Ex) <- names(covObj$cat)
      if (covObj$phLen == 1) {
        intercept <- 1 - covObj$propHaz[, names(covObj$cat)[-1]]
      }
      else {
        intercept <- 1 - apply(covObj$propHaz[, names(covObj$cat)[-1]], 
                               1, sum)
      }
      lifeExp <- cbind(intercept, covObj$propHaz[, names(covObj$cat)[-1]]) %*% 
        Ex
    }
    else {
      meanCont <- apply(matrix(covObj$propHaz[, names(covObj$cont)], 
                               ncol = length(covObj$cont)), 2, mean)
      gam <- sum(parsPrior$gamma * meanCont)
      Ex <- sum(CalcSurv(xx, parsPrior$theta)^exp(gam) * 
                  dxx)
      lifeExp <- rep(Ex, dataObj$n)
    }
  }
  class(parsPrior) <- class(parsIni)
  priorCov <- .CalcParCovObj(covObj, parsPrior, parsCovIni)
  priorAgeObj <- list(lifeExp = lifeExp, theta = priorCov$theta)
  if (class(covObj)[1] %in% c("fused", "propHaz")) {
    priorAgeObj$gamma <- priorCov$gamma
  }
  return(priorAgeObj)
}

# getAnywhere(.BuildPostObj)
.BuildPostObj <- function (ageObj, parObj, parCovObj, dataObj, CalcSurv, priorAgeObj, 
                           fullParObj, covObj) 
{
  postObj <- list()
  postObj$mat <- matrix(0, dataObj$n, 6, dimnames = list(NULL, 
                                                         c("fx", "Sx", "vx", "lx", "px", "postX")))
  postObj <- .CalcPostX(ageObj, parObj, postObj, parCovObj, 
                        1:dataObj$n, CalcSurv, priorAgeObj, fullParObj, covObj, 
                        dataObj)
  return(postObj)
}

# getAnywhere(.PrepJumpObj)
.PrepJumpObj <- function (fullParObj, covObj) 
{
  jump <- c(fullParObj$theta$jump)
  if (fullParObj$class[1] == "theGam") {
    jump <- c(jump, fullParObj$gamma$jump)
  }
  nPar <- length(jump)
  return(list(jump = jump, jumpsMat = matrix(NA, 0, nPar), 
              update = TRUE, updateRate = rep(0, nPar)))
}

# getAnywhere(.RunIniUpdJump)
.RunIniUpdJump <- function (argList, algObj, defTheta, CalcMort, CalcSurv, dataObj, 
                            covObj, userPars, fullParObj, agesIni, parsIni, priorAgeObj, 
                            parsCovIni, postIni, jumps, .jumpObjIni) 
{
  cat("Starting simulation to find jump sd's... ")
  jumpObj <- .jumpObjIni
  parsNow <- parsIni
  agesNow <- agesIni
  parsCovNow <- parsCovIni
  postNow <- postIni
  newJumps <- list()
  newJumps$theta <- fullParObj$theta$jump
  if (class(parsIni)[1] == "theGam") {
    newJumps$gamma <- fullParObj$gamma$jump
  }
  idMp <- which(substr(fullParObj$allNames, 1, 2) != "pi")
  if (algObj$shape != "simple") {
    idC <- which(substr(fullParObj$allNames, 1, 1) == "c")
    idMp <- c(idMp[-idC], idMp[idC])
  }
  idGam <- which(substr(fullParObj$allNames, 1, 2) == "ga")
  nMp <- length(idMp)
  updObj <- list(len = 50)
  updObj$targ <- ifelse("updateRate" %in% names(argList), argList$updateRate, 
                        0.25)
  niter <- updObj$len * 125
  updObj$int <- seq(updObj$len, niter, updObj$len)
  updObj$updVec <- matrix(0, niter, nMp)
  op <- options()
  options(warn = -1)
  for (m in 1:niter) {
    for (mp in idMp) {
      parsNew <- .ProposeMortPars(parsNow, agesNow, newJumps, 
                                  fullParObj, mp)
      newPar <- ifelse(c(parsNew$theta, parsNew$gamma)[mp] != 
                         c(parsNow$theta, parsNow$gamma)[mp], TRUE, FALSE)
      parsCovNew <- .CalcParCovObj(covObj, parsNew, parsCovNow)
      postNew <- .CalcLike(agesNow, parsNew, postNow, parsCovNew, 
                           1:dataObj$n, fullParObj, CalcSurv, dataObj)
      if (!is.na(postNew$mortPars) & newPar) {
        acceptCrit <- exp(postNew$mortPars - postNow$mortPars)
        acceptProb <- runif(1)
        if (acceptCrit > acceptProb) {
          parsNow <- parsNew
          parsCovNow <- parsCovNew
          postNow <- postNew
          updObj$updVec[m, mp] <- 1
        }
      }
    }
    parsNow <- .SamplePiPars(parsNow, agesNow, fullParObj, 
                             dataObj)
    postNow <- .CalcPostPi(parsNow, postNow, agesNow, 1:dataObj$n, 
                           fullParObj, dataObj)
    if (class(parsNow)[2] == "lambda") {
      parsNew <- .ProposeLambda(parsNow, fullParObj)
      postNew <- .CalcPostLambda(parsNew, postNow, agesNow, 
                                 1:dataObj$n, fullParObj)
      acceptCrit <- exp(postNew$lambda - postNow$lambda)
      acceptProb <- runif(1)
      if (acceptCrit > acceptProb) {
        parsNow <- parsNew
        postNow <- postNew
      }
    }
    postNow <- .SumPosts(postNow, 1:dataObj$n)
    if (class(dataObj) == "ageUpd") {
      agesNew <- .SampleAges(agesNow, dataObj, algObj)
      idNew <- which(agesNew$ages[, "birth"] != agesNow$ages[, 
                                                             "birth"] | agesNew$ages[, "death"] != agesNow$ages[, 
                                                                                                                "death"])
      postNew <- .CalcPostX(agesNew, parsNow, postNow, 
                            parsCovNow, dataObj$idNoA, CalcSurv, priorAgeObj, 
                            fullParObj, covObj, dataObj)
      idAccept <- .AcceptAges(postNow, postNew, idNew)
      postNow$mat[idAccept, ] <- postNew$mat[idAccept, 
      ]
      agesNow$ages[idAccept, ] <- agesNew$ages[idAccept, 
      ]
      agesNow$alive[idAccept, ] <- agesNew$alive[idAccept, 
      ]
      postNow$mortPars <- .CalcPostMortPars(parsNow, postNow, 
                                            fullParObj)
      if (class(parsNow)[2] == "lambda") {
        postNow$lambda <- sum(postNow$mat[, "lx"]) + 
          .dtnorm(parsNow$lambda, mean = fullParObj$lambda$priorMean, 
                  sd = fullParObj$lambda$priorSd, lower = 0)
      }
    }
    if (m %in% updObj$int) {
      jumpObj <- .UpdateJumps(jumpObj, updObj, m)
      newJumps$theta[1:fullParObj$theta$len] <- jumpObj$jump[1:fullParObj$theta$len]
      if (class(parsNow)[1] == "theGam") {
        newJumps$gamma <- jumpObj$jump[-c(1:fullParObj$theta$len)]
      }
    }
  }
  options(op)
  aveJumps <- apply(matrix(jumpObj$jumpsMat[nrow(jumpObj$jumpsMat) - 
                                              c(100:0), ], ncol = ncol(jumpObj$jumpsMat)), 2, mean)
  newJumps$theta[1:fullParObj$theta$len] <- aveJumps[1:fullParObj$theta$len]
  if (class(parsIni)[1] == "theGam") {
    newJumps$gamma <- aveJumps[idGam]
  }
  cat(" done.\n\n")
  return(list(updJumps = newJumps, jObject = jumpObj, m = m))
}

# getAnywhere(.CalcQuants)
.CalcQuants <- function (bastaOut, bastaResults, defTheta, fullParObj, algObj, 
                         dataObj, CalcSurv, CalcMort, covObj, agesIni) 
{
  if (class(agesIni)[1] == "ageUpd") {
    nthin <- ceiling((algObj$niter - algObj$burnin + 1)/algObj$thinning)
    bMat <- array(matrix(dataObj$bi, nthin, dataObj$n, byrow = TRUE), 
                  dim = list(nthin, dataObj$n, algObj$nsim))
    if (dataObj$updB) {
      for (sim in 1:algObj$nsim) {
        bMat[, dataObj$idNoB, sim] <- bastaOut[[sim]]$birth
      }
    }
    dMat <- array(matrix(dataObj$di, nthin, dataObj$n, byrow = TRUE), 
                  dim = list(nthin, dataObj$n, algObj$nsim))
    if (dataObj$updD) {
      for (sim in 1:algObj$nsim) {
        dMat[, dataObj$idNoD, sim] <- bastaOut[[sim]]$death
      }
    }
    qdMat <- apply(dMat, 2, quantile, c(0.5, 0.025, 0.975))
    qbMat <- apply(bMat, 2, quantile, c(0.5, 0.025, 0.975))
    xMat <- dMat - bMat
    qxMat <- rbind(round(apply(xMat, 2, mean)), apply(xMat, 
                                                      2, quantile, c(0.025, 0.975)))
  }
  else {
    qbMat <- matrix(dataObj$bi, nrow = 1)
    qdMat <- matrix(dataObj$di, nrow = 1)
    qxMat <- qdMat - qbMat
  }
  ageVec <- seq(0.001, max(qxMat[1, ]) * 5, 0.1)
  mxq <- Sxq <- list()
  if (is.null(covObj$cat)) {
    covNames <- c("noCov")
  }
  else {
    covNames <- paste(".", names(covObj$cat), sep = "")
  }
  for (cov in covNames) {
    if (class(covObj)[1] %in% c("noCov", "propHaz")) {
      thPars <- matrix(bastaResults$params[, defTheta$name], 
                       ncol = defTheta$length)
    }
    else if (class(covObj)[1] == "fused") {
      idTh <- grep(cov, fullParObj$allNames, fixed = TRUE)
      thPars <- matrix(bastaResults$params[, idTh], ncol = length(idTh))
    }
    else {
      if (class(covObj)[2] %in% c("cateCov", "bothCov")) {
        idTh <- grep(cov, fullParObj$allNames, fixed = TRUE)
      }
      else {
        idTh <- grep("Intercept", fullParObj$allNames, 
                     fixed = TRUE)
      }
      thPars <- matrix(bastaResults$params[, idTh], ncol = length(idTh))
      if (!is.null(covObj$cont)) {
        for (pp in names(covObj$cont)) {
          idCon <- which(substr(fullParObj$allNames, 
                                nchar(fullParObj$allNames) - (nchar(pp) - 
                                                                1), nchar(fullParObj$allNames)) == pp)
          thPars <- thPars + matrix(bastaResults$params[, 
                                                        idCon], ncol = length(idCon)) * mean(covObj$inMort[, 
                                                                                                           covObj$cont[pp]])
        }
      }
    }
    colnames(thPars) <- defTheta$name
    if (class(covObj)[1] %in% c("noCov", "inMort")) {
      gamPar <- rep(0, nrow(bastaResults$params))
    }
    else if (class(covObj)[1] == "fused") {
      idGam <- grep("gamma", fullParObj$allNames)
      gamPar <- matrix(bastaResults$params[, idGam], ncol = length(idGam)) %*% 
        c(apply(covObj$propHaz, 2, mean, na.rm = TRUE))
    }
    else {
      if (is.null(covObj$cat)) {
        gamPar <- rep(0, nrow(bastaResults$params))
      }
      else if (cov == sprintf(".%s", names(covObj$cat)[1])) {
        gamPar <- rep(0, nrow(bastaResults$params))
      }
      else {
        gamPar <- bastaResults$params[, grep(cov, fullParObj$allNames, 
                                             fixed = TRUE)]
      }
      if (!is.null(covObj$cont)) {
        idGam <- grep("gamma", fullParObj$allNames)
        idCon <- idGam[which(substr(fullParObj$allNames[idGam], 
                                    7, nchar(fullParObj$allNames[idGam])) %in% 
                               names(covObj$cont))]
        gamPar <- gamPar + t(t(as.matrix(bastaResults$params[, 
                                                             idCon])) * apply(as.matrix(covObj$propHaz[, 
                                                                                                       names(covObj$cont)]), 2, mean, na.rm = TRUE))
      }
    }
    Sxq[[cov]] <- apply(sapply(1:nrow(thPars), function(pp) CalcSurv(ageVec, 
                                                                     t(thPars[pp, ]))^exp(gamPar[pp])), 1, quantile, c(0.5, 
                                                                                                                       0.025, 0.975))
    colnames(Sxq[[cov]]) <- ageVec + algObj$minAge
    id01 <- which(Sxq[[cov]][1, ] < 0.01)[1]
    if (is.na(id01) | length(id01) == 0) 
      id01 <- length(ageVec)
    Sxq[[cov]] <- Sxq[[cov]][, 1:id01]
    mxq[[cov]] <- apply(sapply(1:nrow(thPars), function(pp) CalcMort(ageVec, 
                                                                     t(thPars[pp, ])) * exp(gamPar[pp])), 1, quantile, 
                        c(0.5, 0.025, 0.975))
    colnames(mxq[[cov]]) <- ageVec + algObj$minAge
    mxq[[cov]] <- mxq[[cov]][, 1:id01]
  }
  bastaResults$birthQuant <- qbMat
  bastaResults$deathQuant <- qdMat
  bastaResults$agesQuant <- qxMat
  bastaResults$mortQuant <- mxq
  bastaResults$survQuant <- Sxq
  return(bastaResults)
}

# getAnywhere(.CalcLifeTable)
.CalcLifeTable <- function (bastaResults, lifeTable, object, covObj, algObj) 
{
  if (lifeTable) {
    LT <- list()
    if (is.null(covObj$cat)) {
      covNames <- c("noCov")
    }
    else {
      covNames <- names(covObj$cat)
    }
    for (cov in covNames) {
      if (cov == "noCov") {
        x <- bastaResults$agesQuant[1, bastaResults$birthQuant[1, 
        ] >= algObj$start]
      }
      else {
        x <- bastaResults$agesQuant[1, object[, cov] == 
                                      1 & bastaResults$birthQuant[1, ] >= algObj$start]
      }
      tempLT <- MakeLifeTable(x, ax = 0.5, n = 1)
      tempLT <- subset(tempLT, tempLT$StartAge >= algObj$minAge)
      rownames(tempLT) <- NULL
      LT[[cov]] <- tempLT
    }
  }
  else {
    LT <- NULL
  }
  return(LT)
}

# getAnywhere(sfInit)
sfInit <- function (parallel = NULL, cpus = NULL, type = NULL, socketHosts = NULL, 
                    restore = NULL, slaveOutfile = NULL, nostart = FALSE, useRscript = FALSE) 
{
  reconnect <- FALSE
  initEnv <- new.env()
  if (nostart) 
    return(TRUE)
  if (length(.sfOption) == 0) {
    debug("Setup sfOption...")
    setOption("parallel", FALSE)
    setOption("session", NULL)
    setOption("priority", 1)
    setOption("nodes", 1)
    setOption("stopped", FALSE)
    setOption("init", FALSE)
    data("config", package = "snowfall", envir = initEnv)
    configM <- as.matrix(t(config))
    config <- as.list(configM)
    names(config) <- dimnames(configM)[[2]]
    if (.sfPresetCPUs > 0) 
      setOption("MAXNODES", .sfPresetCPUs)
    else setOption("MAXNODES", as.numeric(config[["MAXNODES"]]))
    setOption("LOCKFILE", "")
    if (as.character(config[["TMPDIR"]]) != "-") 
      setOption("TMPDIR", path.expand(as.character(config[["TMPDIR"]])))
    else {
      if (.Platform$OS.type == "unix") 
        setOption("TMPDIR", file.path(Sys.getenv("R_SESSION_TMPDIR"), 
                                      "sfCluster"))
      else setOption("TMPDIR", "")
    }
    setOption("RESTOREFILES", NULL)
    setOption("RESTOREUPDATE", 5)
    setOption("RESTORE", FALSE)
    setOption("CURRENT", NULL)
    setOption("type", "SOCK")
    setOption("sockHosts", NULL)
    if (as.character(config[["RESTDIR"]]) != "-") 
      setOption("RESTDIR", path.expand(as.character(config[["RESTDIR"]])))
    else setOption("RESTDIR", file.path(Sys.getenv("HOME"), 
                                        ".sfCluster", "restore"))
    rm(config, envir = initEnv)
  }
  else {
    reconnect <- TRUE
    if (.sfOption$stopped && !.sfOption$init) 
      debug("Irregluar init state (error on previous init)...")
    if (!.sfOption$stopped && .sfOption$init) {
      message("Explicit sfStop() is missing: stop now.")
      sfStop()
    }
  }
  searchCommandline(parallel, cpus = cpus, type = type, socketHosts = socketHosts, 
                    restore = restore)
  if (getOption("verbose") && !reconnect) 
    print(.sfOption)
  if (.sfOption$parallel && !nostart) {
    if (startedWithSfCluster() && is.null(.sfOption$session)) 
      stop("No session-ID but parallel run with sfCluster (something went wrong here?)...")
    if (is.null(.sfOption$nodes) || is.na(as.numeric(.sfOption$nodes))) 
      setOption("nodes", 2)
    else setOption("nodes", as.numeric(.sfOption$nodes))
    libList <- list(PVM = "rpvm", MPI = "Rmpi", NWS = "nws", 
                    SOCK = "")
    if (libList[[.sfOption$type]] != "") {
      if (!require(libList[[.sfOption$type]], character.only = TRUE)) {
        message(paste("Failed to load required library:", 
                      libList[[.sfOption$type]], "for parallel mode", 
                      .sfOption$type, "\nFallback to sequential execution"))
        return(sfInit(parallel = FALSE))
      }
      else message(paste("Library", libList[[.sfOption$type]], 
                         "loaded."))
    }
    if (startedWithSfCluster()) {
      tmp <- file.path(.sfOption$TMPDIR, paste("rout_", 
                                               .sfOption$session, sep = ""))
      if (!reconnect) 
        dirCreateStop(.sfOption$TMPDIR)
    }
    else tmp <- ifelse(is.null(slaveOutfile), "/dev/null", 
                       slaveOutfile)
    setDefaultClusterOptions(type = .sfOption$type)
    setDefaultClusterOptions(homogenous = FALSE)
    if (.sfOption$type == "SOCK") {
      if (is.null(.sfOption$sockHosts) || (length(.sfOption$sockHosts) == 
                                           0)) 
        setOption("sockHosts", c(rep("localhost", .sfOption$nodes)))
      else setOption("nodes", length(.sfOption$sockHosts))
      setOption("cluster", try(makeCluster(.sfOption$sockHosts, 
                                           type = "SOCK", outfile = tmp, homogenous = TRUE)))
    }
    else if (.sfOption$type == "PVM") {
      setOption("cluster", try(makeCluster(.sfOption$nodes, 
                                           outfile = tmp)))
    }
    else if (.sfOption$type == "NWS") {
      if (is.null(.sfOption$sockHosts) || (length(.sfOption$sockHosts) == 
                                           0)) 
        setOption("sockHosts", c(rep("localhost", .sfOption$nodes)))
      else setOption("nodes", length(.sfOption$sockHosts))
      setOption("cluster", try(makeNWScluster(.sfOption$sockHosts[1:.sfOption$nodes], 
                                              type = "NWS", outfile = tmp)))
    }
    else {
      setOption("cluster", try(makeMPIcluster(.sfOption$nodes, 
                                              outfile = tmp, homogenous = TRUE, useRscript = useRscript)))
    }
    if (is.null(.sfOption$cluster) || inherits(.sfOption$cluster, 
                                               "try-error")) 
      stop(paste("Starting of snow cluster failed!", geterrmessage(), 
                 .sfOption$cluster))
    setOption("init", TRUE)
    setOption("stopped", FALSE)
    if (!reconnect) {
      if (!is.null(.sfOption$LOCKFILE) && file.exists(.sfOption$LOCKFILE)) {
        if (unlink(.sfOption$LOCKFILE) != 0) 
          warning("Unable to remove startup lockfile: ", 
                  .sfOption$LOCKFILE)
        else message("Startup Lockfile removed: ", .sfOption$LOCKFILE)
      }
      if (getOption("verbose")) {
        if (tmp == "/dev/null") 
          message("Slave output suppressed. Use 'slaveOutfile' to activate.")
        else message(paste("Temporary log for STDOUT/STDERR (on each node): ", 
                           tmp, "\n", "Cluster started with", .sfOption$nodes, 
                           "CPUs.", "\n"))
      }
      else debug(paste("Temporary log for STDOUT/STDERR (on each node): ", 
                       tmp, "\n", "Cluster started with", .sfOption$nodes, 
                       "CPUs.", "\n"))
      .startInfo <- strsplit(Sys.info(), "\n")
      .startMsg <- paste(sep = "", "JOB STARTED AT ", date(), 
                         " ON ", .startInfo$nodename, " (OS", .startInfo$sysname, 
                         ") ", .startInfo$release, "\n")
      sfExport(".sfOption", ".startMsg", local = TRUE, 
               namespace = "snowfall", debug = DEBUG)
      sfCat(.startMsg, "\n", master = FALSE)
      sfCat(paste("R Version: ", R.version$version.string, 
                  "\n\n"))
      sfRemove(".startMsg")
    }
    else sfExport(".sfOption", local = FALSE, namespace = "snowfall")
  }
  else {
    setOption("init", TRUE)
    setOption("stopped", FALSE)
    setOption("cluster", NULL)
  }
  if (sfParallel()) {
    message(paste("snowfall ", packageDescription("snowfall")$Version, 
                  " initialized (using snow ", packageDescription("snow")$Version, 
                  "): parallel execution on ", sfCpus(), " CPUs.\n", 
                  sep = ""))
  }
  else {
    message(paste("snowfall", packageDescription("snowfall")$Version, 
                  "initialized: sequential execution, one CPU.\n"))
  }
  return(invisible(TRUE))
}

# getAnywhere(sfExport)
sfExport <- function (..., list = NULL, local = TRUE, namespace = NULL, debug = FALSE, 
                      stopOnError = TRUE) 
{
  sfCheck()
  if (!sfParallel()) {
    warning("sfExport() writes to global environment in sequential mode.\n")
  }
  names <- fetchNames(...)
  if (!is.null(list)) {
    if (!length(list) || !all(sapply(list, function(x) is.symbol(x) || 
                                     is.character(x)))) {
      if (stopOnError) 
        stop("'list' must contain names or character strings")
      else {
        warning("Error in sfExport: 'list' must contain names or character strings")
        return(invisible(FALSE))
      }
    }
    names <- c(names, list)
  }
  for (name in names) {
    if (!is.null(namespace) && is.character(namespace)) {
      val <- tryCatch(getFromNamespace(name, namespace), 
                      error = function(x) {
                        NULL
                      })
      if (!is.null(val) && !inherits(val, "try-error")) {
        res <- sfClusterCall(assign, name, val, env = globalenv(), 
                             stopOnError = FALSE)
        if (is.null(res) || !all(checkTryErrorAny(res))) {
          if (stopOnError) 
            stop(paste("Error exporting '", name, "': ", 
                       geterrmessage(), sep = ""))
          else {
            warning(paste("Error exporting '", name, 
                          "': ", geterrmessage(), sep = ""))
            return(invisible(FALSE))
          }
        }
        next
      }
    }
    if (local) {
      found <- FALSE
      for (pframe in seq(1, sys.nframe())) {
        if (exists(name, inherits = FALSE, envir = sys.frame(-pframe))) {
          found <- TRUE
          if (debug) {
            definedIn <- gsub("\n", "|", as.character(sys.call(-pframe)))
            cat("Export '", name, "' defined in '", definedIn, 
                "'", "\n", sep = "")
            print(get(name, envir = sys.frame(-pframe)))
          }
          res <- sfClusterCall(assign, name, get(name, 
                                                 envir = sys.frame(-pframe)), env = globalenv(), 
                               stopOnError = FALSE)
          if ((is.null(res) && !is.null(get(name, envir = sys.frame(-pframe)))) || 
              !all(checkTryErrorAny(res))) {
            if (stopOnError) 
              stop(paste("Error exporting '", name, "': ", 
                         geterrmessage(), sep = ""))
            else {
              message(paste("Error exporting '", name, 
                            "': ", geterrmessage(), sep = ""))
              return(invisible(FALSE))
            }
          }
          break
        }
      }
      if (!found) {
        if (stopOnError) 
          stop(paste("Unknown/unfound variable ", name, 
                     " in export. (local=", local, ")", sep = ""))
        else {
          message(paste("Unknown/unfound variable ", 
                        name, " in export. (local=", local, ")", 
                        sep = ""))
          return(invisible(FALSE))
        }
      }
    }
    else {
      if (exists(name, inherits = FALSE, envir = globalenv())) {
        res <- sfClusterCall(assign, name, get(name, 
                                               inherits = FALSE, envir = globalenv()), env = globalenv(), 
                             stopOnError = FALSE)
        if (is.null(res) || !all(checkTryErrorAny(res))) {
          if (stopOnError) 
            stop(paste("Error exporting global '", name, 
                       "': ", geterrmessage(), sep = ""))
          else {
            warning(paste("Error exporting global '", 
                          name, "': ", geterrmessage(), sep = ""))
            return(invisible(TRUE))
          }
        }
      }
      else {
        if (stopOnError) 
          stop(paste("Unknown variable ", name, " in export."))
        else {
          warning(paste("Unknown variable ", name, " in export."))
          return(invisible(TRUE))
        }
      }
    }
  }
  invisible(TRUE)
}

# getAnywhere(sfLibrary)
sfLibrary <- function (package, pos = 2, lib.loc = NULL, character.only = FALSE, 
                       warn.conflicts = TRUE, keep.source = NULL, verbose = getOption("verbose"), 
                       version, stopOnError = TRUE) 
{
  sfCheck()
  sfPars <- list()
  if (!missing(package)) {
    if (character.only) {
      if (is.character(package)) 
        sfPars$package <- package
      else stop(paste("Package", package, "is no character string."))
    }
    else sfPars$package <- deparse(substitute(package))
  }
  sfPars$character.only <- TRUE
  sfPars$pos <- pos
  sfPars$lib.loc <- lib.loc
  sfPars$warn.conflicts <- warn.conflicts
  if (as.integer(R.version$major) < 3) {
    if (is.null(keep.source)) 
      keep.source = getOption("keep.source.pkgs")
    sfPars$keep.source <- keep.source
  }
  sfPars$verbose <- verbose
  sfPars$logical.return <- TRUE
  if (!missing(version)) 
    sfPars$version <- version
  if (sfParallel()) {
    if (!is.null(sfPars$lib.loc)) 
      sfPars$lib.loc <- absFilePath(sfPars$lib.loc)
    setVar(".sfPars", sfPars)
    sfExport(".sfPars", local = FALSE, namespace = "snowfall")
    result <- try(sfClusterEval(do.call("library", .sfPars)))
    if (inherits(result, "try-error") || (length(result) != 
                                          sfCpus()) || !all(checkTryErrorAny(result)) || !all(unlist(result))) {
      if (stopOnError) 
        stop(paste("Stop: error loading library on slave(s):", 
                   .sfPars$package))
      else {
        warning(paste("Error loading library on slave(s):", 
                      package))
        return(invisible(FALSE))
      }
    }
    else {
      sfCat(paste("Library", .sfPars$package, "loaded.\n"))
      message(paste("Library", .sfPars$package, "loaded in cluster.\n"))
    }
  }
  result <- try(do.call("library", sfPars))
  sfRemove(".sfPars")
  if (inherits(result, "try-error") || !result) {
    if (stopOnError) {
      warning(paste("Unable to load library:", package))
      return(invisible(FALSE))
    }
    else stop(paste("Unable to load library:", package))
  }
  else {
    if (verbose) 
      message(paste("Library", package, "loaded.\n"))
    return(invisible(TRUE))
  }
}

# getAnywhere(sfClusterApplyLB)
sfClusterApplyLB <- function (x, fun, ...) 
{
  sfCheck()
  checkFunction(fun)
  if (sfParallel()) 
    return(clusterApplyLB(sfGetCluster(), x, fun, ...))
  else return(lapply(x, fun, ...))
}

# getAnywhere(sfRemoveAll)
sfRemoveAll <- function (except = NULL, debug = FALSE, hidden = TRUE) 
{
  sfCheck()
  if (sfParallel()) {
    if (hidden) 
      sfTmpAll <- sfClusterEval(ls(pos = globalenv(), all.names = TRUE))
    else sfTmpAll <- sfClusterEval(ls(pos = globalenv(), 
                                      all.names = FALSE))
    if (length(sfTmpAll) == 0) {
      message("sfRemoveAll: problems fetching variables from nodes (or none existant)...\n")
      return(invisible(FALSE))
    }
    sfTmp <- sfTmpAll[[which.max(sapply(sfTmpAll, length))]]
    if (length(sfTmp) > 0) {
      if (!is.null(except)) {
        if (is.list(except)) 
          except <- unlist(except)
        if (!is.vector(except)) {
          warning("sfRemoveAll: except is not a vector.\n")
          return(invisible(FALSE))
        }
        for (i in match(except, sfTmp)) sfTmp[i] <- NA
        sfTmp <- sort(sfTmp, na.last = NA)
      }
      if (debug) {
        message("sfRemoveAll: Remove variables from nodes:")
        message(paste(sfTmp, collapse = ", "))
      }
      sfExport("sfTmp", local = TRUE)
      sfClusterEval(rm(list = sfTmp, pos = globalenv()))
      sfClusterEval(rm("sfTmp", pos = globalenv()))
    }
    else {
      message("sfRemoveAll: no variables on nodes.\n")
      return(invisible(FALSE))
    }
    return(invisible(TRUE))
  }
  else {
    message("sfRemoveAll() ignored in sequential mode.\n")
    return(invisible(TRUE))
  }
}

# getAnywhere(sfStop)
sfStop <- function (nostop = FALSE) 
{
  if (nostop) 
    return(TRUE)
  if (exists(".sfOption") && (length(.sfOption) > 0)) {
    if (!.sfOption$stopped && .sfOption$init && .sfOption$parallel) {
      message("\nStopping cluster\n")
      stopCluster(.sfOption$cluster)
    }
    setOption("stopped", TRUE)
    setOption("parallel", FALSE)
    deleteRestoreFiles()
  }
  invisible(NULL)
}



# does it work with the full size function if I run it on a sub-sample of code?
# may be trying to log a negative number?

basta_test <- function (object, studyStart, studyEnd, minAge = 0, model = "GO", 
                        shape = "simple", covarsStruct = "fused", niter = 11000, 
                        burnin = 1001, thinning = 20, recaptTrans = studyStart, nsim = 1, 
                        parallel = FALSE, ncpus = 2, lifeTable = TRUE, updateJumps = TRUE, 
                        ...) 
{
  argList <- list(...)
  bastaIntVars <- c("algObj", "defTheta", "CalcMort", "CalcSurv", 
                    "dataObj", "covObj", "userPars", "fullParObj", "agesIni", 
                    "parsIni", "priorAgeObj", "parsCovIni", "postIni", "jumps", 
                    ".Random.seed")
  algObj <- .CreateAlgObj(model, shape, studyStart, studyEnd, 
                          minAge, covarsStruct, recaptTrans, niter, burnin, thinning, 
                          updateJumps, nsim)
  .FindErrors(object, algObj)
  dataObj <- .PrepDataObj(object, algObj)
  defTheta <- .SetDefaultTheta(algObj)
  CalcMort <- .DefineMort(algObj)
  CalcSurv <- .DefineSurv(algObj)
  covObj <- .CreateCovObj(object, dataObj, algObj)
  algObj$covStruc <- class(covObj)[1]
  userPars <- .CreateUserPar(covObj, argList)
  fullParObj <- .BuildFullParObj(covObj, defTheta, algObj, 
                                 userPars, dataObj)
  agesIni <- .PrepAgeObj(dataObj, algObj)
  parsIni <- .DefineIniParObj(fullParObj)
  parsCovIni <- .BuildParCovObj(covObj, parsIni)
  priorAgeObj <- .SetPriorAgeDist(fullParObj, CalcSurv, dataObj, 
                                  covObj, parsIni, parsCovIni)
  postIni <- .BuildPostObj(agesIni, parsIni, parsCovIni, dataObj, 
                           CalcSurv, priorAgeObj, fullParObj, covObj)
  jumps <- list()
  jumps$theta <- fullParObj$theta$jump
  if (fullParObj$class[1] == "theGam") {
    jumps$gamma <- fullParObj$gamma$jump
  }
  Start <- Sys.time()
  if (updateJumps) {
    .jumpObjIni <- .PrepJumpObj(fullParObj, covObj)
    updatedJumps <- .RunIniUpdJump(argList, algObj, defTheta, 
                                   CalcMort, CalcSurv, dataObj, covObj, userPars, fullParObj, 
                                   agesIni, parsIni, priorAgeObj, parsCovIni, postIni, 
                                   jumps, .jumpObjIni)
    jumps <- updatedJumps$updJumps
  }
  if (nsim > 1) {
    cat("Multiple simulations started...\n\n")
    if (parallel) {
      opp <- options()
      options(warn = -1)
      sfInit(parallel = TRUE, cpus = ncpus)
      sfExport(list = c(bastaIntVars, ".Random.seed"))
      sfLibrary("BaSTA", character.only = TRUE, warn.conflicts = FALSE)
      bastaOut <- sfClusterApplyLB(1:nsim, .RunBastaMCMC, 
                                   algObj, defTheta, CalcMort, CalcSurv, dataObj, 
                                   covObj, userPars, fullParObj, agesIni, parsIni, 
                                   priorAgeObj, parsCovIni, postIni, jumps)
      sfRemoveAll(hidden = TRUE)
      sfStop()
      options(opp)
    }
    else {
      bastaOut <- lapply(1:nsim, .RunBastaMCMC, algObj, 
                         defTheta, CalcMort, CalcSurv, dataObj, covObj, 
                         userPars, fullParObj, agesIni, parsIni, priorAgeObj, 
                         parsCovIni, postIni, jumps)
    }
  }
  else {
    cat("Simulation started...\n\n")
    bastaOut <- lapply(1:nsim, .RunBastaMCMC, algObj, defTheta, 
                       CalcMort, CalcSurv, dataObj, covObj, userPars, fullParObj, 
                       agesIni, parsIni, priorAgeObj, parsCovIni, postIni, 
                       jumps)
  }
  End <- Sys.time()
}



pM <-  c(-3,0.2,0.001,-4,0.05) # prior means for parameters in a Gompertz model with a bathtub shape
pSD <- c(1,0.05,0.001,1,0.01)  # prior standard deviation for parameters in a Gompertz model with a bathtub shape
Go.Bt <- basta_test(basta_data, studyStart = Ss,
                    studyEnd = Se, nsim = sim, niter = iter,
                    burnin = warmup, thinning = thin,
                    thetaPriorMean = pM,thetaPriorSd = pSD,
                    models = c("GO"), shape = c("bathtub"))





#### lets just start trying shit and see what happens! ####
### Start the function in parallel -- don't need to do multiple permutations as all male and the permutations were changing males and females around in the unknowns. Does this still need to use a cluster?
ncl <- 4               # number of cpus used in cluster - check how many cores you have available/want to use
cl <- makeCluster(ncl) # create the cluster element
doParallel::registerDoParallel(cl)

# Set the arguments for the iterations
iter <- 500   # number of iterations
sim <- 4         # number of chains
warmup <- iter/2 # number of iterations to use during warmup/burn-in period
thin <- 11      # thinning number

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

stopCluster(cl)

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

# Gompertz Makeham == u(age)=A*expB*age+C -- run ####
rethinking::dens(fmsb::GompertzMakeham(1,2,3,c(1:60)))
rethinking::dens(fmsb::GompertzMakeham(1,1,1,c(1:60)))
rethinking::dens(fmsb::GompertzMakeham(1,2,1,c(1:60)))
rethinking::dens(fmsb::GompertzMakeham(1,3,1,c(1:60)))
rethinking::dens(fmsb::GompertzMakeham(1,1,2,c(1:60)))
rethinking::dens(fmsb::GompertzMakeham(1,1,3,c(1:60)))
rethinking::dens(fmsb::GompertzMakeham(2,1,1,c(1:60)))
rethinking::dens(fmsb::GompertzMakeham(3,1,1,c(1:60)))
rethinking::dens(fmsb::GompertzMakeham(0.1,0.2,0.3,c(1:60)))
# haven't the faintest idea how these are working but none of them look even remotely like a mortality curve....  

rethinking::dens(fmsb::GompertzMakeham(0.1,0.2,3,c(1:60))) # final value compresses everything down to a much smaller set of x values, but don't know what that means.
rethinking::dens(fmsb::GompertzMakeham(0.1,0.5,1,c(1:60))) # middle value affects initial mortality
rethinking::dens(fmsb::GompertzMakeham(0.1,0.1,1,c(1:60))) # middle value affects initial mortality

pM <-  c(0.1,0.1,1)      # prior means for parameters in a Gompertz model with a Makeham shape -- these look pretty shit in the model but I don't know what would be better from plotting the priors
pSD <- c(0.1, 0.1, 0.5)  # prior standard deviation for parameters in a Gompertz model with a Makeham shape
Go.Mk <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               models = c("GO"), shape = c("Makeham"))
summary(Go.Mk)
plot(Go.Mk)                       # Plot traces
plot(Go.Mk, plot.trace = FALSE)   # Plot survival and mortality curves

# WEIBULL BATHTUB -- not yet worked out priors ####
pM <-  c() # prior means for parameters in a Weibull model with a bathtub shape
pSD <- c()  # prior standard deviation for parameters in a Weibull model with a bathtub shape
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

# WEIBULL MAKEHAM -- not yet worked out priors ####
pM <-  c() # prior means for parameters in a Weibull model with a Makeham shape
pSD <- c()  # prior standard deviation for parameters in a Weibull model with a Makeham shape
Wb.Mk <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               models = c("WE"), shape = c("Makeham"))
summary(Wb.Mk)
plot(Wb.Mk)                       # Plot traces
plot(Wb.Mk, plot.trace = FALSE)   # Plot survival and mortality curves

# Exponential simple -- requires 2 parameters (why 2 and not just 1?), no bathtub or Makeham options -- run ####
pM <-  c(1,1) # prior means for parameters in a Exponential model with a simple shape
pSD <- c(0.5,0.5)  # prior standard deviation for parameters in a Exponential model with a simple shape
Ex.Si <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               models = c("EX"), shape = c("simple"))
summary(Ex.Si)
plot(Ex.Si)                       # Plot traces
plot(Ex.Si, plot.trace = FALSE)   # Plot survival and mortality curves -- this actually looks pretty good, and almost identical to other curves.

# LOGISTIC BATHTUB -- not yet worked out priors ####
pM <-  c() # prior means for parameters in a Logistic model with a bathtub shape
pSD <- c()  # prior standard deviation for parameters in a Logistic model with a bathtub shape
Lg.Bt <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               models = c("LO"), shape = c("bathtub"))
summary(Lg.Bt)
plot(Lg.Bt)                       # Plot traces
plot(Lg.Bt, plot.trace = FALSE)   # Plot survival and mortality curves

stopCluster(cl)

# Logistic simple ####
rethinking::dens(rlogis(1000,0,1))
rethinking::dens(rlogis(1000,0,5))
rethinking::dens(rlogis(1000,10,8))
rethinking::dens(rlogis(1000,15,10))
rethinking::dens(rlogis(1000,16,9))

pM <-  c(16,10) # prior means for parameters in a Logistic model with a simple shape
pSD <- c(2,2)  # prior standard deviation for parameters in a Logistic model with a simple shape
Lg.Si <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               models = c("LO"), shape = c("simple"))
summary(Lg.Si)
plot(Lg.Si)                       # Plot traces
plot(Lg.Si, plot.trace = FALSE)   # Plot survival and mortality curves

# LOGISTIC MAKEHAM -- not yet worked out priors ####
pM <-  c() # prior means for parameters in a Logistic model with a Makeham shape
pSD <- c()  # prior standard deviation for parameters in a Logistic model with a Makeham shape
Lg.Mk <- basta(basta_data, studyStart = Ss,
               studyEnd = Se, nsim = sim, niter = iter,
               burnin = warmup, thinning = thin,
               thetaPriorMean = pM,thetaPriorSd = pSD,
               models = c("LO"), shape = c("Makeham"))
summary(Lg.Mk)
plot(Lg.Mk)                       # Plot traces
plot(Lg.Mk, plot.trace = FALSE)   # Plot survival and mortality curves


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













