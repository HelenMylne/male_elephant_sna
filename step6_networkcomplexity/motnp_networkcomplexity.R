## Bayesian analysis of ALERT data
#### Set up ####
# load packages
library(tidyverse)
library(rstan)
library(rethinking)
library(igraph)
library(cmdstanr)

################ 2) Create data frame ################
### import data for aggregated model (binomial)
counts_df <- read_delim('data_processed/motnp_bayesian_allpairwiseevents_splitbygrouptype_22.01.13.csv', delim = ',')

# correct sex_1, which has loaded in as a logical vector not a character/factor
unique(counts_df$sex_1) # FALSE or NA
sex_1 <- data.frame(sex_1 = counts_df$id_1)
sex_1 <- sex_1 %>% separate(sex_1, c("sex", "number"), sep = 1, remove = FALSE) ; unique(sex_1$sex) # F, M, U
counts_df$sex_1 <- as.character(sex_1$sex) ; rm(sex_1)
counts_df$sex_1 <- as.character(counts_df$sex_1)
str(counts_df)  # sex_1 still comes up as logical?

# create variable for age difference
unique(counts_df$age_category_1) # "0-3","1-2","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-15","15-19","20-25","20-35","25-40","35-50","40+","50+",NA 
counts_df$age_category_1 <- ifelse(counts_df$age_category_1 == '1-2','0-3',counts_df$age_category_1)
counts_df$age_cat_id_1 <- ifelse(counts_df$age_category_1 == '0-3', 1,
                                 ifelse(counts_df$age_category_1 == '3-4', 1,
                                        ifelse(counts_df$age_category_1 == '4-5', 1,
                                               ifelse(counts_df$age_category_1 == '5-6', 2,
                                                      ifelse(counts_df$age_category_1 == '6-7', 2,
                                                             ifelse(counts_df$age_category_1 == '7-8', 2,
                                                                    ifelse(counts_df$age_category_1 == '8-9', 2,
                                                                           ifelse(counts_df$age_category_1 == '9-10', 2,
                                                                                  ifelse(counts_df$age_category_1 == '10-15', 3,
                                                                                         ifelse(counts_df$age_category_1 == '15-19', 4,
                                                                                                ifelse(counts_df$age_category_1 == '20-25', 5,
                                                                                                       ifelse(counts_df$age_category_1 == '20-35', 5,
                                                                                                              ifelse(counts_df$age_category_1 == '25-40', 6,
                                                                                                                     ifelse(counts_df$age_category_1 == '35-50', 6,
                                                                                                                            ifelse(counts_df$age_category_1 == '40+', 7,
                                                                                                                                   ifelse(counts_df$age_category_1 == '50+', 7, counts_df$age_category_1))))))))))))))))
counts_df[is.na(counts_df$age_cat_id_1),]   # U8 doesn't have an age but not a problem -- won't be part of the main analysis. Is a calf so age_cat_id = 1
counts_df$age_cat_id_1[which(is.na(counts_df$age_cat_id_1))] <- 1

unique(counts_df$age_category_1[counts_df$age_class_1 == 'Calf'])      # shouldn't include any ages over 4-5
unique(counts_df$age_category_1[counts_df$age_class_1 == 'Juvenile'])  # shouldn't include any ages under 5-6
unique(counts_df$age_category_1[counts_df$age_class_1 == 'Pubescent']) # shouldn't include any ages under 9-10 or over 15-19
unique(counts_df$age_category_1[counts_df$age_class_1 == 'Adult'])     # shouldn't include any ages under 20-25

counts_df$age_class_1 <- ifelse(counts_df$age_cat_id_1 == 1, 'Calf',
                                ifelse(counts_df$age_cat_id_1 == 2, 'Juvenile',
                                       ifelse(counts_df$age_cat_id_1 > 4, 'Adult','Pubescent')))

unique(counts_df$age_category_2) # "0-3","1-2","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-15","15-19","20-25","20-35","25-40","35-50","40+","50+",NA 
counts_df$age_category_2 <- ifelse(counts_df$age_category_2 == '1-2','0-3',counts_df$age_category_2)
counts_df$age_cat_id_2 <- ifelse(counts_df$age_category_2 == '0-3', 1,
                                 ifelse(counts_df$age_category_2 == '3-4', 1,
                                        ifelse(counts_df$age_category_2 == '4-5', 1,
                                               ifelse(counts_df$age_category_2 == '5-6', 2,
                                                      ifelse(counts_df$age_category_2 == '6-7', 2,
                                                             ifelse(counts_df$age_category_2 == '7-8', 2,
                                                                    ifelse(counts_df$age_category_2 == '8-9', 2,
                                                                           ifelse(counts_df$age_category_2 == '9-10', 2,
                                                                                  ifelse(counts_df$age_category_2 == '10-15', 3,
                                                                                         ifelse(counts_df$age_category_2 == '15-19', 4,
                                                                                                ifelse(counts_df$age_category_2 == '20-25', 5,
                                                                                                       ifelse(counts_df$age_category_2 == '20-35', 5,
                                                                                                              ifelse(counts_df$age_category_2 == '25-40', 6,
                                                                                                                     ifelse(counts_df$age_category_2 == '35-50', 6,
                                                                                                                            ifelse(counts_df$age_category_2 == '40+', 7,
                                                                                                                                   ifelse(counts_df$age_category_2 == '50+', 7, counts_df$age_category_2))))))))))))))))
counts_df[is.na(counts_df$age_cat_id_2),]   # U8 doesn't have an age but not a problem -- won't be part of the main analysis. Is a calf so age_cat_id = 1
counts_df$age_cat_id_2[which(is.na(counts_df$age_cat_id_2))] <- 1

counts_df$age_class_2 <- ifelse(counts_df$age_cat_id_2 == 1, 'Calf',
                                ifelse(counts_df$age_cat_id_2 == 2, 'Juvenile',
                                       ifelse(counts_df$age_cat_id_2 > 4, 'Adult','Pubescent')))

unique(counts_df$age_category_2[counts_df$age_class_2 == 'Calf'])      # shouldn't include any ages over 4-5
unique(counts_df$age_category_2[counts_df$age_class_2 == 'Juvenile'])  # shouldn't include any ages under 5-6
unique(counts_df$age_category_2[counts_df$age_class_2 == 'Pubescent']) # shouldn't include any ages under 9-10 or over 15-19
unique(counts_df$age_category_2[counts_df$age_class_2 == 'Adult'])     # shouldn't include any ages under 20-25

### correct dem_class with corrected age classes
counts_df$dem_class_1 <- ifelse(counts_df$age_class_1 == 'Adult', paste('A',counts_df$sex_1, sep = ''),
                                ifelse(counts_df$age_class_1 == 'Pubescent', paste('P',counts_df$sex_1, sep = ''),
                                       ifelse(counts_df$age_class_1 == 'Juvenile', paste('J',counts_df$sex_1, sep = ''),
                                              paste('C',counts_df$sex_1, sep = ''))))
counts_df$dem_class_2 <- ifelse(counts_df$age_class_2 == 'Adult', paste('A',counts_df$sex_2, sep = ''),
                                ifelse(counts_df$age_class_2 == 'Pubescent', paste('P',counts_df$sex_2, sep = ''),
                                       ifelse(counts_df$age_class_2 == 'Juvenile', paste('J',counts_df$sex_2, sep = ''),
                                              paste('C',counts_df$sex_2, sep = ''))))

### add column for age difference between dyad
counts_df$age_diff <- abs(as.numeric(counts_df$age_cat_id_1) - as.numeric(counts_df$age_cat_id_2))

### add column for total number of sightings per pair
counts_df$count_dyad <- (counts_df$count_1 + counts_df$count_2) - counts_df$all_events  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.

### add column for total number of sightings per pair where they were NOT together
counts_df$apart <- counts_df$count_dyad - counts_df$all_events

### reassign dyad numbers to remove gaps
counts_df$node_1_nogaps <- as.integer(as.factor(counts_df$node_1))
counts_df$node_2_nogaps <- as.integer(as.factor(counts_df$node_2))+1
counts_df$dyad_id_nogaps <- as.integer(as.factor(counts_df$dyad_id))

################ 3) Create simulated data set ################
### Population
# 120 individuals: 50 male (10 each of ages 3-7), 50 female (10 each of ages 3-7), 20 unknown (10 each of ages 1-2)
# males no preference for any other individual -- edge weight 0.2
# females in 8 cliques -- edge weight 0.8, otherwise 0.2
# unknowns attached to a single female -- edge weight 0.95

# assign females to family units
cliques_F <- sample(1:8,50,replace=T)
table(cliques_F)

# create data frame of all individuals
population <- data.frame(id = rep(NA,120),
                         num = c(1:50,1:50,1:20),
                         sex = c(rep('M',50), rep('F',50), rep('U',20)),
                         age = c(rep(c(3:7), each = 10), rep(c(3:7), each = 10), rep(c(1:2), each = 10)),
                         clique = c(rep(NA, 50), cliques_F, rep(NA, 20)))
population$id <- paste(population$sex, population$num, sep = '')

# randomly assign calves to females
mothers <- data.frame(id = population$id[101:120],
                      family = sample(population$id[61:100], replace = F, size = 20))
mothers
population <- left_join(x = population, y = mothers, by = 'id')

# need a clique for U elephants and a family for F elephants
females <- population[population$sex == 'F', c(1,5,6)]
unknown <- population[population$sex == 'U', c(1,5,6)]
females$mother <- females$id ; unknown$mother <- unknown$family   # create variable to join by
families <- left_join(x = unknown, y = females, by = 'mother')    # join to give unknowns a clique
unknown <- families[,c(1,4,6)]

families <- left_join(x = females, y = unknown, by = 'mother')    # join to give mothers their calf
females <- families[,c(1,5,2)]

colnames(unknown) <- c('id','mother','clique')
colnames(females) <- c('id','mother','clique')
families <- rbind(females, unknown)                               # put mother and calf data into a single data frame

# join to combine male and female/unknown data together
population <- left_join(x = population, y = families, by = 'id')
population <- population[,c(1,3,4,7,8)]
colnames(population)[c(4,5)] <- c('family','clique')

### Create dyadic data frame
dyads <- data.frame(id_1 = rep(population$id, each = 120), id_2 = rep(population$id, 120))               # create new data frame of all dyads
dyads$same <- ifelse(dyads$id_1 == dyads$id_2, 'yes', 'no') ; dyads <- dyads[dyads$same == 'no', c(1,2)] # remove self-dyads (e.g. M1-M1)
dyads$node_1 <- as.integer(as.factor(dyads$id_1)) ; dyads$node_2 <- as.integer(as.factor(dyads$id_2))    # create factor of nodes
dyads$dyad <- ifelse(dyads$node_1 < dyads$node_2,                               # create variable of dyads
                     paste(dyads$id_1, dyads$id_2, sep = '_'),
                     paste(dyads$id_2, dyads$id_1, sep = '_'))

dyads <- data.frame(dyad = unique(dyads$dyad))                                  # create new data frame of only unique dyads
dyads <- separate(dyads, col = dyad, into = c('id_1','id_2'), remove = FALSE)   # nodes columns
dyads$node_1 <- as.integer(as.factor(dyads$id_1))                               # nodes columns
dyads$node_2 <- as.integer(as.factor(dyads$id_2))                               # nodes columns
dyads$dyad_id <- as.integer(as.factor(dyads$dyad))                              # dyad index variable

population$id_1 <- population$id ; population$id_2 <- population$id  # variable to join information on each individual
dyads <- left_join(x = dyads, y = population, by = 'id_1')           # add information about node 1
dyads <- dyads[,c(1:6,8:11)]
colnames(dyads)[c(3,7:10)] <- c('id_2','sex_1','age_1','family_1','clique_1')

dyads <- left_join(x = dyads, y = population, by = 'id_2')           # add information about node 2
dyads <- dyads[,c(1:10,12:15)]
colnames(dyads)[c(2,11:14)] <- c('id_1','sex_2','age_2','family_2','clique_2')

dyads <- dyads[,c(1:7,11,8,12,9,13,10,14)]

### Assign type of dyad
dyads$pair_type <- ifelse(is.na(dyads$clique_1) | is.na(dyads$clique_2), 'no_group',
                          ifelse(dyads$clique_1 == dyads$clique_2,
                                 ifelse(is.na(dyads$family_1) | is.na(dyads$family_2), 'group',
                                        ifelse(dyads$family_1 == dyads$id_2 | dyads$family_2 == dyads$id_1,
                                               'family', 'group')),
                                 'no_group'))

### Assign edge weight to dyad pairs
for(i in 1:nrow(dyads)){
  if(dyads$pair_type[i] == 'family')  {dyads$edge[i] <- rethinking::rbeta2(1,0.95,50)   # mother-calf: edge weight around 0.95
  } else 
    if(dyads$pair_type[i] == 'group') {dyads$edge[i] <- rethinking::rbeta2(1,0.8,50)    # same family, not mother-calf: edge ~ 0.8
    } else 
      dyads$edge[i] <- rethinking::rbeta2(1,0.20,50)                                        # different family: edge weight around 0.2
}
summary(dyads$edge)

boxplot(dyads$edge ~ dyads$pair_type)                                                   # check simulations are around what I expected
points(y = dyads$edge, x = rep(2, length(which(dyads$pair_type == 'group'))),           # had a few problems with this one, but looks good now
       col = col.alpha('black',0.2), pch = 19)

### clean up environment
rm(families, females, mothers, unknown, cliques_F)

### Sample observations
N <- 20
dyads$event_count <- rbinom(nrow(dyads), N, prob = dyads$edge)

### plot simulated observations against assigned edge weight to check model input is correlates to true value
plot(event_count ~ edge, data = dyads,
     las = 1, ylab = 'sightings together', xlab = 'assigned edge weight',
     pch = 16, col = col.alpha(rangi2, 0.2))

### add additional columns (dyadic sex, age and dem_class)
head(dyads)
# assign each node to age category (adult, pubescent, juvenile or calf) from age class
dyads$age_cat_1 <- ifelse(dyads$age_1 == 1, 'C', ifelse(dyads$age_1 == 2, 'J', ifelse(dyads$age_1 < 5, 'P', 'A')))
dyads$age_cat_2 <- ifelse(dyads$age_2 == 1, 'C', ifelse(dyads$age_2 == 2, 'J', ifelse(dyads$age_2 < 5, 'P', 'A')))

# combine age categories into single variable for dyad
dyads$age_dyad <- ifelse(dyads$age_cat_1 == 'A',
                         ifelse(dyads$age_cat_2 == 'A', 'A_A',
                                ifelse(dyads$age_cat_2 == 'P', 'A_P',
                                       ifelse(dyads$age_cat_2 == 'J', 'A_J', 'A_C'))),
                         ifelse(dyads$age_cat_1 == 'P',
                                ifelse(dyads$age_cat_2 == 'A', 'A_P',
                                       ifelse(dyads$age_cat_2 == 'P', 'P_P',
                                              ifelse(dyads$age_cat_2 == 'J', 'P_J', 'P_C'))),
                                ifelse(dyads$age_cat_1 == 'J',
                                       ifelse(dyads$age_cat_2 == 'A', 'A_J',
                                              ifelse(dyads$age_cat_2 == 'P', 'P_J',
                                                     ifelse(dyads$age_cat_2 == 'J', 'J_J', 'J_C'))),
                                       ifelse(dyads$age_cat_2 == 'A', 'A_C',
                                              ifelse(dyads$age_cat_2 == 'P', 'P_C',
                                                     ifelse(dyads$age_cat_2 == 'J', 'J_C', 'C_C'))))))
unique(dyads$age_dyad) # "P_P" "A_P" "P_C" "P_J" "A_A" "A_C" "A_J" "C_C" "J_C" "J_J"
dyads$age_diff <- abs(dyads$age_1 - dyads$age_2)

# create composite measure of demography per node (e.g. AM = adult male, CU = calf of unknown sex)
dyads$dem_class_1 <- paste(dyads$age_cat_1, dyads$sex_1, sep = '')
dyads$dem_class_2 <- paste(dyads$age_cat_2, dyads$sex_2, sep = '')

# combine demographies into single variable for dyad
dyads$demf_1 <- as.integer(as.factor(dyads$dem_class_1)) ; dyads$demf_2 <- as.integer(as.factor(dyads$dem_class_2))
dyads$dem_dyad <- ifelse(dyads$demf_1 < dyads$demf_2,
                         paste(dyads$dem_class_1, dyads$dem_class_2, sep = '_'),
                         paste(dyads$dem_class_2, dyads$dem_class_1, sep = '_'))
dyads$dem_dyad <- ifelse(dyads$dem_dyad == 'CU_JU', 'JU_CU',           # standardise dem_dyad so always in same order
                         ifelse(dyads$dem_dyad == 'CU_PF', 'PF_CU',
                                ifelse(dyads$dem_dyad == 'CU_PM','PM_CU',
                                       ifelse(dyads$dem_dyad == 'JU_PF','PF_JU',
                                              ifelse(dyads$dem_dyad == 'JU_PM','PM_JU',dyads$dem_dyad)))))
dyads$dem_dyad_f <- as.integer(as.factor(dyads$dem_dyad))              # make index variable

dyads$dem_diff <- ifelse(dyads$dem_class_1 == dyads$dem_class_2, 0, 1) # binary: same or different
dyads$sex_dyad <- paste(dyads$sex_1, dyads$sex_2, sep = '_')           # composite sex variable for dyad
dyads$sex_diff <- ifelse(dyads$sex_1 == dyads$sex_2, 0, 1)             # binary: same or different
dyads <- dyads[,c(1:23,26:30)]                                         # remove demf variables -- not needed

# together vs apart per dyad
summary(dyads$event_count)             # together, all out of 20
dyads$apart <- N - dyads$event_count   # apart = 20-together

# clear environment
rm(population, i, N)

################ 4) Social complexity ################
### create functions for measuring social complexity -- taken from Weiss et al. 2019 ####
#Function to fit a J component mixture model
bt = function(n,Y,J,maxiter=1000,maxrep=50,tolfun=1e-6,minprior=0){
  fit.rho = function(Q,A,n,Y){
    K = length(Q)
    nY = length(Y)
    nm = matrix(n,ncol=K,nrow=nY)
    Ym = matrix(Y,ncol=K,nrow=nY)
    Am = matrix(A,nrow=nY,ncol=K,byrow=T)
    Qm = matrix(Q,nrow=nY,ncol=K,byrow=T)
    ll.rho = function(r){
      b1 = Qm*(1/r - 1)
      b2 = ((Qm-1)*(r-1))/r
      ll = dbetabinom.ab(Ym,nm,b1,b2)*Am
      -sum(log(rowSums(ll)))
    }
    optimise(ll.rho,interval=c(0,1))$minimum
  } #function to fit overdispersion parameter
  mean = list() #lists to hold parameters
  freq = list()
  rho = list()
  mus = list()
  lllq = NULL 
  nY = length(Y) #number of dyads
  for(mr in 1:maxrep){
    ll0 = 0
    nm = matrix(n,nrow=nY,ncol=J) #a matrix for the n's
    Ym = matrix(Y,nrow=nY,ncol=J) #a matrix of Y's
    Q = 0.8*runif(J,0,1)
    A = rep(1/J,J)
    K = J
    for(j in 1:maxiter){
      Q = Q[A>minprior] #drop components with no weight (typically not an issue with fuzzy clustering)
      A = A[A>minprior]
      K = length(A) #get new number of components
      nm = matrix(n,ncol=K,nrow=nY)
      Ym = matrix(Y,ncol=K,nrow=nY)
      Am = matrix(A,nrow=nY,ncol=K,byrow=T) #turn component paramters into matrices
      Qm = matrix(Q,nrow=nY,ncol=K,byrow=T) 
      ll = dbinom(Ym,nm,Qm)*Am #get likelihoods
      mu = ll/rowSums(ll) #responsibilities
      A = colMeans(mu) #get fractions
      Q = colSums(mu*Y)/colSums(mu*n) #new parameters
      if(any(is.na(Q))|any(is.na(A))) break 
      lll = sum(log(rowSums(ll))) #log-likelihood
      if(abs(lll-ll0)<tolfun & j>(maxiter/10)) break #check for convergence
      ll0 = lll
    }
    rho[[mr]] = fit.rho(Q,A,n,Y)
    mean[[mr]] = Q #save parameters
    freq[[mr]] = A
    lllq = c(lllq,lll)
    lllm = max(lllq,na.rm=T)
    if(mr>4 & sum(abs(lllm-lllq)<tolfun,na.rm = T)>4) break
  }
  mean = mean[[which.max(lllq)]]
  freq = freq[[which.max(lllq)]]
  rho = rho[[which.max(lllq)]]
  freq = freq[order(mean)]
  mean = mean[order(mean)]
  K = length(mean)
  Qm = matrix(mean,nrow=nY,ncol=K,byrow=T)
  Am = matrix(freq,nrow=nY,ncol=K,byrow=T)
  Ym = matrix(Y, nrow=nY,ncol=K)
  nm = matrix(n, nrow = nY, ncol= K)
  ll = dbinom(Ym,nm,matrix(mean,nrow=nY,ncol=K,byrow=T))*matrix(freq,nrow=nY,ncol=K,byrow=T)
  mu = ll/rowSums(ll)
  qq = mean*freq
  qq = qq/sum(qq)
  qq = qq[qq>0] #prevents NAs in the entropy estimate
  S = -sum(qq*log(qq))
  lllm = sum(log(rowSums(ll)))
  cl = mu[mu!=0]
  vv = -sum(cl*log(cl))
  AIC = 2*(2*J-1) - 2*lllm
  BIC = (2*J-1)+log(nY)*(2*J-1) - 2*lllm
  ICL = BIC + 2*vv
  return(list(K.in = J, K.out = length(mean),Mean=mean,Frequency=freq, S = S, rho = rho, logLik=lllm, BIC = BIC, AIC = AIC, ICL = ICL, nrep=mr))
}

#Function to perform full fitting process. 
#Returns a summary table and the best model as determined by selected criterion.
binom_assoc_mixt = function(Den,Num,maxJ=9,nrep=20,criterion="ICL"){
  summary = matrix(ncol = 7,nrow=maxJ) #matrix for summary output
  models = list() #list to hold all models
  colnames(summary) = c("K.in", "K.out", "S", "rho", "AIC","BIC","ICL")
  worse = 0
  for(i in 1:maxJ){
    print(paste("Fitting",i,"Component(s)"))
    res = bt(Den,Num,i,maxrep=nrep)
    summary[i,"K.in"] = i
    summary[i,"K.out"] = res$K.out
    summary[i,"S"] = res$S
    summary[i, "rho"] = res$rho
    summary[i,"AIC"] = res$AIC
    summary[i,"BIC"] = res$BIC
    summary[i,"ICL"] = res$ICL
    if(i > 1){
      if(criterion == "AIC"){ #If the current model is worse than the previous, record this.
        worse = ifelse(summary[i,"AIC"]>summary[(i-1),"AIC"], worse+1, 0)
      }
      if(criterion == "BIC"){
        worse = ifelse(summary[i,"BIC"]>summary[(i-1),"BIC"], worse+1, 0)
      }
      if(criterion == "ICL"){
        worse = ifelse(summary[i,"ICL"]>summary[(i-1),"ICL"], worse+1, 0)
      }
    }
    models[[i]] = res
    if(worse > 1) break #If two models in a row have gotten worse, stop fitting
  }
  if(criterion == "BIC"){ #Get the best model, as chosen by the specified criteria
    best = models[[which.min(summary[,"BIC"])]]
  }
  if(criterion == "AIC"){
    best = models[[which.min(summary[,"AIC"])]]
  }
  if(criterion == "ICL"){
    best = models[[which.min(summary[,"ICL"])]]
  }
  summary = summary[1:i,]
  return(list(summary = summary, best_model = best)) #return list containing summary table and best model
}


### Weiss et al. 2019 example ####
#Generate simulated data
K = 3 #number of types
N = 20 #number of individuals
N.dyad = (N*(N-1))/2 #number of dyads
mean.d = 80 #sampling effort
bad.par = T
#generate valid parameters
while(bad.par){
  mu = runif(K,0,1)
  a = runif(K,0,1)
  a = a/sum(a)
  if(min(dist(mu))>=0.1 & min(a)>=0.1/K) bad.par = F
}
rho = runif(K,0,0.015) #overdispersion

b1 = mu*(1/rho - 1) #shape parameters from means and overdispersion
b2 = ((mu-1)*(rho-1))/rho

k = sample(K, size = N.dyad, rep = T, prob = a) #assign classes
p = rbeta(n=N.dyad,shape1=b1[k],shape2=b2[k]) #assign association probabilities
d = rpois(N.dyad,mean.d) #assign denominators
x = rbinom(n=N.dyad,size=d,prob=p) #assign numerators

#fit model
model.fit = binom_assoc_mixt(d,x,criterion="ICL")

#check summary table
model.fit$summary

#check best model
best.model = model.fit$best

#number of components
best.model$K.out

#means
best.model$Mean

#frequencies
best.model$Frequency

#overdispersion
best.model$rho

#complexity
best.model$S

### fit to data ####
# work out how each value in simulation fits to real data
# K = number of types -- to be estimated, already in there
# N = number of individuals = 218 for independent males or 463 for all elephants
# N.dyad = number of dyads = (N*(N-1))/2
# mean.d = sampling effort -- already in there
# bad.par = T --- NO IDEA

#d = rpois(N.dyad,mean.d) #assign denominators -- sightings total = counts_df$count_dyad
#x = rbinom(n=N.dyad,size=d,prob=p) #assign numerators -- sightings together = counts_df$all_events

motnp_fit <- binom_assoc_mixt(counts_df$count_dyad,counts_df$all_events,criterion="ICL") # if I've understood it all correctly, no complexity at all... this doesn't seem right, since this contains all of the data and surely the mother-calves should come into a different category than adult males??

#check summary table
motnp_fit$summary # model with only a single component by far the best according to all indices

#check best model
best.model = motnp_fit$best

#number of components
best.model$K.out # 1 component only

#means
best.model$Mean

#frequencies
best.model$Frequency

#overdispersion
best.model$rho

#complexity
best.model$S
