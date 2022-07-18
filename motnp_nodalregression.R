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

################ 4) DYADIC REGRESSION ################
# A common type of analysis in studies of social networks is to consider factors that affectedge weights.  This type of analysis is often called dyadic regression, where explanatoryfactors such as sex or age difference are regressed against edge weight. The edge weightsin dyadic regression are non-independent, as variables such as age difference or sex dif-ference are inherently linked to individual-level attributes.  This means that effects dueto age or sex in individualiaffect all dyads that connecti.  This non-independence canbe controlled for by including node-level effects in the regression (Tranmer et al., 2014).Using the count data example discussed earlier, we propose to conduct dyadic regres-sion using a standard regression model of the form:
#         X_ij ∼ Poisson(lambda_ij)
#         log(lambda_ij) = weight_ij + Location_ij
#         weight_ij ∼ Normal(beta_0 + beta_1*x_ij + u_i + u_j, sigma_sq)
# where beta_0 is the intercept parameter, beta_1 is the slope parameter, ui and uj are parameters accounting for the effect of node membership on edge weight, and sigma is the standard deviation of the residuals.  The dashed line indicates that the model can be fit in two separate parts:the first part being the core edge weight model, and the second part being the dyadic regression model.  It is not strictly necessary to separate out the model like this, but do-ing so makes it possible to fit the edge weight model once and conduct multiple types of analysis on it afterwards without fitting the entire model again.Theβ0andβ1parameters can be interpreted in the same way as usual regression coefficients, and will be accompanied with credible intervals describing the probability distribution of coefficient values.  There are several ways to gain inference from these values.  The overlap of the credible interval with zero is often used to indicate the presence of a ‘significant effect’, though this approach has been criticised for several reasons(Gelman et al., 2012; Ogle et al., 2019).  An alternative option is to use Bayes factorsto quantify support for or against the hypothesis H1:β16=0 (Kass and Raftery, 1995;Jeffreys, 1998). In regression settings, the Bayes factor for a point null hypothesis can be calculated quickly using the Savage-Dickey method (Wagenmakers et al., 2010).The model shown here is intended as a minimal example of a simple linear regression with node-level effects to account for the non-independence of edges. The regression can be extended to support different family distributions, hierarchical effects, non-linearities,and more. We have included examples of diagnostics and some extensions in the example code.

d_reg <- cmdstan_model("models/dyadic_regression_hkm_22.03.03.stan")
d_reg

simdat_ls <- list(
  n_dyads  = nrow(dyads),        # Number of dyads
  together = dyads$event_count,  # Number of sightings of each dyad together
  apart    = dyads$apart         # Number of sightings of each dyad apart
)

### Fit model
dyad_reg_sim <- d_reg$sample(
  data = simdat_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4
)

### check model
dyad_reg_sim$summary()
################ 5) MEASURE node centrality values -- Hart method -- given up on that bit ################
### read in data
draws_motnp2.2 <- read_csv('Random/draws1_motnp2.2_22.02.16')

### trim it down in the hopes that it might work -- first chain only, male-male links only
draws <- draws_motnp2.2[1:1000, which(counts_df$dem_type == 'AM_AM' | counts_df$dem_type == 'AM_PM' | counts_df$dem_type == 'PM_PM')]
rm(draws_motnp2.2)
males <- counts_df[which(counts_df$dem_type == 'AM_AM' | counts_df$dem_type == 'AM_PM' | counts_df$dem_type == 'PM_PM'),]
males$dyad_id_nogaps <- as.integer(as.factor(males$dyad_id))
males$node_1_nogaps <- as.integer(as.factor(males$node_1))
males$node_2_nogaps <- as.integer(as.factor(males$node_2))+1

### create count values
N <- length(unique(males$id_1))+1
num_iterations <- nrow(draws)

### create array of draws per dyad (distributions)
adj_tensor <- array(0, c(N, N, num_iterations))
for (i in min(males$dyad_id_nogaps):max(males$dyad_id_nogaps)) {
  dyad_row <- males[males$dyad_id_nogaps == i,]
  adj_tensor[dyad_row$node_1_nogaps,dyad_row$node_1_nogaps,1:1000] <- draws[,i]
}
adj_tensor[,,1]



adj_tensor <- array(0, c(N, N, num_iterations))#,
#dimnames = list(c(unique(counts_df$id_1),'U9'),
#                c('F1',unique(counts_df$id_2)),
#                NULL))
#for (i in 1:nrow(counts_df)) {
#  dyad_row <- counts_df[i, ]
#  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws[, i]
#}

for (i in 1:nrow(counts_df)) {
  dyad_row <- counts_df[counts_df$dyad_id_nogaps == 1, ]
  adj_tensor[dyad_row$node_1_nogaps, dyad_row$node_2_nogaps, ] <- draws[, 1]
}

### create centrality matrix
centrality_matrix <- matrix(0, nrow = num_iterations, ncol = N)
for (i in 1:num_iterations) {
  g <- graph_from_adjacency_matrix(adj_tensor[, , i], mode="undirected", weighted=TRUE)
  centrality_matrix[i, ] <- strength(g)
}

colnames(centrality_matrix) <- c(sort(unique(counts_df$id_1)),unique(counts_df$id_2)[462])
head(centrality_matrix)

centrality_quantiles <- t(apply(centrality_matrix, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975))))
centrality_quantiles

################ 6) NODAL REGRESSION ################
# The final common type of network analysis we will cover here is nodal regression, where a regression is performed to analyse the relationship between a nodal network metric(such as centrality) and nodal traits (such as age and sex).  These analyses are usually used to assess how network position depends on various biological factors, but can also be used where network position is a predictor.  Since node metrics are derivative measures of the network, uncertainty in edge weights should ideally propagate through to node metrics, and on to coefficients in regression analyses, giving an accurate estimate of the total uncertainty in inferred parameters.  The core of the nodal regression is similar to dyadic regression, taking the form:
#X(n)i j∼Poisson(λ(n)i jD(n)i j)
#logÄλ(n)i jä=ωi j+L(n)i j
#Mi(ω)∼Normal(β0+β1xi,σ2)
# where β0 is the intercept parameter, β1 is the slope parameter, and σ is the standard deviation of the residuals. Mi(ω) denotes the node metric estimates when applied to the edge weights estimated in the top part of the model. Calculating node metrics within the model may present a practical challenge when using standard model fitting software, as many node metrics can not be described purely in terms of closed-form equations (Freeman, 1978).  In this case, splitting up the model along the dashed line becomes an important step, as it allows the edge weights to be estimated using a standard piece of software, and leaves the regression part of the model to be fit using a sampler that supports numeric estimates of network centralities. As part of the code we have made available, we have written a custom Metropolis-Hastings sampler that allows numeric estimates of node metrics to be used when fitting the model. Importantly, our sampler samples from the joint distribution of edge weights, rather than sampling edge weights independently. This ensures that the effect of any structure in the observation data (such as location effects) is maintained and propagated through to the node metric estimates, and subsequently to regression coefficient estimates.

