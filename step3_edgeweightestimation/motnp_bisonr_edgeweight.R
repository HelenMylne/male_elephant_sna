#### set up ####
# load packages
# library(tidyverse) ; library(dplyr) ; library(cmdstanr) ; library(bisonR) ; library(asnipe) ; library(sna) ; library(raster)
library(tidyverse, lib.loc = '../packages/') # library(tidyverse)
library(dplyr, lib.loc = '../packages/')     # library(dplyr)
#library(rstan, lib.loc = '../packages/')     # library(rstan)
library(cmdstanr, lib.loc = '../packages/')  # library(cmdstanr)
library(bisonR, lib.loc = '../packages/')    # library(bisonR)
library(asnipe, lib.loc = '../packages/')    # library(asnipe)
library(sna, lib.loc = '../packages/')       # library(sna)
library(raster, lib.loc = '../packages/')    # library(raster)

# information
sessionInfo()
R.Version()
rstan::stan_version()

# set seed
set.seed(12345)

# set cmdstanr path
#set_cmdstan_path('H:/rlibs/4.2.1/')

#### create data lists ####
### import data for aggregated model (binomial) -- counts of positive associations and total sightings
counts_df <- read_csv('../data_processed/motnp_bayesian_binomialpairwiseevents.csv')

# correct sex_1, which has loaded in as a logical vector not a character/factor
unique(counts_df$sex_1) # FALSE or NA
sex_1 <- data.frame(sex_1 = counts_df$id_1)
sex_1 <- sex_1 %>% separate(sex_1, c("sex", "number"), sep = 1, remove = FALSE) ; unique(sex_1$sex) # F, M, U
counts_df$sex_1 <- as.character(sex_1$sex) ; rm(sex_1) ; gc()
str(counts_df)  # sex_1 still comes up as logical, but now contains the right levels

# create variable for age difference
unique(counts_df$age_category_1) # "0-3","1-2","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-15","15-19","20-25","20-35","25-40","35-50","40+","50+",NA --> convert anything <5 to cat 1, and anything 5-10 to cat 2
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
counts_df$age_cat_id_1[which(is.na(counts_df$age_cat_id_1))] <- 1 # U8 doesn't have an age but not a problem -- won't be part of the main analysis. Is a calf so age_cat_id = 1

counts_df$age_class_1 <- ifelse(counts_df$age_cat_id_1 == 1, 'Calf',
                                ifelse(counts_df$age_cat_id_1 == 2, 'Juvenile',
                                       ifelse(counts_df$age_cat_id_1 > 4, 'Adult','Pubescent')))

unique(counts_df$age_category_2) # "0-3","1-2","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-15","15-19","20-25","20-35","25-40","35-50","40+","50+",NA 
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
counts_df$age_cat_id_2[which(is.na(counts_df$age_cat_id_2))] <- 1   # U8 doesn't have an age but not a problem -- won't be part of the main analysis. Is a calf so age_cat_id = 1

counts_df$age_class_2 <- ifelse(counts_df$age_cat_id_2 == 1, 'Calf',
                                ifelse(counts_df$age_cat_id_2 == 2, 'Juvenile',
                                       ifelse(counts_df$age_cat_id_2 > 4, 'Adult','Pubescent')))

# correct dem_class with corrected age classes
counts_df$dem_class_1 <- ifelse(counts_df$age_class_1 == 'Adult', paste('A',counts_df$sex_1, sep = ''),
                                ifelse(counts_df$age_class_1 == 'Pubescent', paste('P',counts_df$sex_1, sep = ''),
                                       ifelse(counts_df$age_class_1 == 'Juvenile', paste('J',counts_df$sex_1, sep = ''),
                                              paste('C',counts_df$sex_1, sep = ''))))
counts_df$dem_class_2 <- ifelse(counts_df$age_class_2 == 'Adult', paste('A',counts_df$sex_2, sep = ''),
                                ifelse(counts_df$age_class_2 == 'Pubescent', paste('P',counts_df$sex_2, sep = ''),
                                       ifelse(counts_df$age_class_2 == 'Juvenile', paste('J',counts_df$sex_2, sep = ''),
                                              paste('C',counts_df$sex_2, sep = ''))))

# combined dem_class of dyad
counts_df$age_class_id_1 <- ifelse(counts_df$age_class_1 == 'Adult',4,
                                   ifelse(counts_df$age_class_1 == 'Pubescent',3,
                                          ifelse(counts_df$age_class_1 == 'Juvenile',2,1)))
counts_df$age_class_id_2 <- ifelse(counts_df$age_class_2 == 'Adult',4,
                                   ifelse(counts_df$age_class_2 == 'Pubescent',3,
                                          ifelse(counts_df$age_class_2 == 'Juvenile',2,1)))
counts_df$dem_type <- ifelse(counts_df$age_class_id_1 > counts_df$age_class_id_2,
                             paste(counts_df$dem_class_1, counts_df$dem_class_2, sep = '_'), # when 1 is older: dc1_dc2
                             ifelse(counts_df$age_class_id_1 < counts_df$age_class_id_2,
                                    paste(counts_df$dem_class_2, counts_df$dem_class_1, sep = '_'), # when 2 is older: dc2_dc1
                                    ifelse(counts_df$sex_1 == 'F',                                  # when age1 = age2...
                                           paste(counts_df$dem_class_1, counts_df$dem_class_2, sep = '_'), # ...when 1 is F: dc1_dc2
                                           ifelse(counts_df$sex_2 == 'F',
                                                  paste(counts_df$dem_class_2, counts_df$dem_class_1, sep = '_'), # ...when 2 is F: dc2_dc1
                                                  ifelse(counts_df$sex_1 == 'M', # ...when neither are F...
                                                         paste(counts_df$dem_class_1, counts_df$dem_class_2, sep = '_'), # ...when 1 is M: dc1_dc2
                                                         paste(counts_df$dem_class_2, counts_df$dem_class_1, sep = '_')))))) # ...when 2 is M or both are U: dc2_dc1
sort(unique(counts_df$dem_type))

# add column for age difference between dyad
counts_df$age_diff <- abs(as.numeric(counts_df$age_cat_id_1) - as.numeric(counts_df$age_cat_id_2))

# add column for total number of sightings per pair
counts_df$count_dyad <- (counts_df$count_1 + counts_df$count_2) - counts_df$event_count  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.

# add column for total number of sightings per pair where they were NOT together
counts_df$apart <- counts_df$count_dyad - counts_df$event_count

# reassign dyad numbers to remove gaps
counts_df$node_1_nogaps <- as.integer(as.factor(counts_df$node_1))
counts_df$node_2_nogaps <- as.integer(as.factor(counts_df$node_2))+1

### remove all except counts data frame
#rm(list = ls()[which(ls() != 'counts_df')]) ; gc()

### load in ages 
motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

### remove dead individuals (also removes all unknown sex calves)
counts_df <- counts_df %>% 
  filter(name_1 != 'Richard' & name_2 != 'Richard') %>% 
  filter(name_1 != 'Gabriel' & name_2 != 'Gabriel') %>% 
  filter(name_1 != 'Tori'    & name_2 != 'Tori')

### filter counts data frame down to males only
counts_df <- counts_df %>% 
  filter(id_1 %in% unique(motnp_ages$id)) %>% 
  filter(id_2 %in% unique(motnp_ages$id))

### remove young males
counts_df <- counts_df %>% 
  filter(age_cat_id_1 >= 3) %>% 
  filter(age_cat_id_2 >= 3)

### create nodes id variable which is 1:max(id) as currently covers all elephants, not just males
counts_df$node_1_males <- as.integer(as.factor(counts_df$node_1_nogaps))
counts_df$node_2_males <- as.integer(as.factor(counts_df$node_2_nogaps))+1

### standardise dyad_id for males only
counts_df$dyad_males <- as.integer(as.factor(counts_df$dyad_id))

### write out data as ready to go into model
# write_csv(counts_df, '../data_processed/motnp_binomialpairwiseevents_malesonly.csv')
# counts_df <- read_csv('../data_processed/motnp_binomialpairwiseevents_malesonly.csv')

### add time marker
print(paste0('data read in at ', Sys.time()))

#### edge weights -- stronger priors, males only ####
# add pdf output file
pdf(file = '../outputs/motnp_bisonr_edgeweight_strongprior.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[, c('node_1_males','node_2_males','event_count','count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### set priors
priors <- get_default_priors('binary') # obtain structure for bison model priors
priors$edge <- 'normal(-2.5,1.5)' 
#priors$edge <- "prior_string('
#real conditional_prior_lpdf(real mu, int seen_together) {
#    if (seen_together == 1) {
#      return normal_lpdf(mu | -2.5, 1.5);
#    } else {
#      return normal_lpdf(mu | -5, 3);
#    }
#  }
#')"     # slightly right skewed so that more density is towards 0 but still very weak and all values are possible
priors$fixed <- 'normal(-5,3)'       # slightly less wide than before, but zero centred so allows age to have a positive or negative impact on centrality and edge weight
prior_check(priors, 'binary')

# extract prior_parameters -- code copied directly from bisonR GitHub
extract_distribution <- function(prior_string) {
  parameter_string <- stringr::str_replace_all(prior_string, " ", "")
  parameter_split <- stringr::str_split(parameter_string, "\\(|\\)|,")[[1]]
  distribution_name <- parameter_split[1]
  parameter_values <- parameter_split[2:(length(parameter_split) - 1)]
  parameter_values <- as.numeric(parameter_values)
  return(list(distribution_name=distribution_name, parameter_values=parameter_values))
}
extract_prior_parameters <- function(priors) {
  prior_parameters <- list()
  for (parameter_name in names(priors)) {
    prior_distribution <- extract_distribution(priors[parameter_name])
    distribution_name <- prior_distribution$distribution_name
    parameter_values <- prior_distribution$parameter_values
    if (distribution_name == "normal") {
      prior_parameters[[paste("prior", parameter_name, "mu", sep="_")]] <- parameter_values[1]
      prior_parameters[[paste("prior", parameter_name, "sigma", sep="_")]] <- parameter_values[2]
    }
    if (distribution_name == "half-normal") {
      prior_parameters[[paste("prior", parameter_name, "sigma", sep="_")]] <- parameter_values[1]
    }
    if (distribution_name == "beta") {
      prior_parameters[[paste("prior", parameter_name, "alpha", sep="_")]] <- parameter_values[1]
      prior_parameters[[paste("prior", parameter_name, "beta", sep="_")]] <- parameter_values[2]
    }
    if (distribution_name == "gamma") {
      prior_parameters[[paste("prior", parameter_name, "alpha", sep="_")]] <- parameter_values[1]
      prior_parameters[[paste("prior", parameter_name, "beta", sep="_")]] <- parameter_values[2]
    }
  }
  # print(prior_parameters)
  prior_parameters
}
prior_parameters <- extract_prior_parameters(priors)
prior_parameters$prior_edge_mu <- 'c(-2.5, -5)'
prior_parameters$prior_edge_sigma <- list(1.5, 3)

### create dummy variable indicating if ever seen together or not -- use to select which prior to draw edge weight from
counts_df_model$seen_together <- ifelse(counts_df_model$event > 0, 1, 0)
counts_df_model$never_together <- 1 - counts_df_model$seen_together

### new prior for zero-inflated dyads (never seen together)
 plot(density( LaplacesDemon::invlogit(rnorm(1000, -5, 3)) ), col = 'blue')      # new -- dyads never sighted together
lines(density( LaplacesDemon::invlogit(rnorm(1000, -2.5, 1.5)) ), col = 'red')   # original -- dyads ≥ 1 seen together 

### run edge weight model
motnp_edge_weights_strongpriors <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id)#*seen_together,   # count of sightings together given count of total sightings as a result of the individuals contained within the dyad
  data = counts_df_model,
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

motnp_edge_weights_strongpriors$stan_model$code()
#  "data {"                                                                                                                 
#  "  int<lower=0> num_rows; // Number of data points"                                                                      
#  "  int<lower=0> num_edges; // Number of edge weights to estimate"                                                        
#  "  int<lower=0> num_fixed; // Number of fixed effect parameters"                                                         
#  "  int<lower=0> num_random; // Number of random effect parameters"                                                       
#  "  int<lower=0> num_random_groups; // Number of random effect groups"                                                    
#  ""                                                                                                                       
#  "  array[num_rows] int event; // Outcome for each data point (presence/absence)"                                         
#  "  array[num_rows] int divisor; // Duration of each observation"                                                         
#  "  array[num_rows] int dyad_ids; // Dyad IDs of each observation for indexing edge weights"                              
#  "  matrix[num_rows, num_fixed] design_fixed; // Design matrix for fixed effects"                                         
#  "  matrix[num_rows, num_random] design_random; // Design matrix for random effects."                                     
#  "  array[num_random] int<lower=0, upper=num_random_groups> random_group_index; // Index for groupings for random effects"
#  ""                                                                                                                       
#  "  real prior_edge_mu; // Prior mean for fixed effects"                                                                  
#  "  real<lower=0> prior_edge_sigma; // Prior standard deviation for fixed effects"                                        
#  "  real prior_fixed_mu; // Prior mean for fixed effects"                                                                 
#  "  real<lower=0> prior_fixed_sigma; // Prior standard deviation for fixed effects"                                       
#  "  real prior_random_mean_mu; // Prior mean on centralisation of random effects"                                         
#  "  real<lower=0> prior_random_mean_sigma; // Prior standard deviation on centralisation of random effects"               
#  "  real<lower=0> prior_random_std_sigma; // Prior standard deviation on dispersion of random effects"                    
#  "  real<lower=0> prior_zero_prob_alpha; // Prior alpha on zero inflation"                                                
#  "  real<lower=0> prior_zero_prob_beta; // Prior beta on zero inflation"                                                  
#  ""                                                                                                                       
#  "  int<lower=0, upper=1> priors_only; // Whether to sample from only the priors"                                         
#  "  int<lower=0, upper=1> partial_pooling; // Whether to pool edge weight estimates"                                      
#  "  int<lower=0, upper=1> zero_inflated; // Whether to use zero-inflated edge model"                                      
#  "}"                                                                                                                      
#  ""                                                                                                                       
#  "parameters {"                                                                                                           
#  "  vector[num_edges] edge_weight; // Parameters for edge weights."                                                       
#  "  vector[num_fixed] beta_fixed; // Parameters for fixed effects."                                                       
#  "  vector[num_random] beta_random; // Parameters for random effects."                                                    
#  "  vector[num_random_groups] random_group_mu; // Hyperpriors for random effects (mean)."                                 
#  "  vector<lower=0>[num_random_groups] random_group_sigma; // Hyperpriors for random effects (std. dev.)."                
#  "  vector<lower=0>[partial_pooling] edge_sigma; // Random effect for edge weight pooling."                               
#  "  vector<lower=0>[zero_inflated] zero_prob; // Zero inflated parameter for probability of zeroes."                      
#  "}"                                                                                                                      
#  ""                                                                                                                       
#  "transformed parameters {"                                                                                               
#  "  vector[num_rows] predictor;"                                                                                          
#  "  predictor = rep_vector(0, num_rows);"                                                                                 
#  "  if (num_edges > 0) {"                                                                                                 
#  "    predictor += edge_weight[dyad_ids];"                                                                                
#  "  }"                                                                                                                    
#  "  if (num_fixed > 0) {"                                                                                                 
#  "    predictor += design_fixed * beta_fixed;"                                                                            
#  "  }"                                                                                                                    
#  "  if (num_random > 0) {"                                                                                                
#  "    predictor += design_random * beta_random;"                                                                          
#  "  }"                                                                                                                    
#  "}"                                                                                                                      
#  ""                                                                                                                       
#  "model {"                                                                                                                
#  "  if (!priors_only) {"                                                                                                  
#  "    // Main model"                                                                                                      
#  "    if (zero_inflated == 0) {"                                                                                          
#  "      event ~ binomial(divisor, inv_logit(predictor));"                                                                 
#  "    } else {"                                                                                                           
#  "      for (i in 1:num_rows) {"                                                                                          
#  "        if (event[i] == 0) {"                                                                                           
#  "          target += log_sum_exp("                                                                                       
#  "            bernoulli_lpmf(1 | zero_prob[1]),"                                                                          
#  "            bernoulli_lpmf(0 | zero_prob[1]) + binomial_lpmf(event[i] | divisor[i], inv_logit(predictor[i]))"           
#  "          );"                                                                                                           
#  "        } else {"                                                                                                       
#  "          target += bernoulli_lpmf(0 | zero_prob[1]) + binomial_lpmf(event[i] | divisor[i], inv_logit(predictor[i]));"  
#  "        }"                                                                                                              
#  "      }"                                                                                                                
#  "      zero_prob ~ beta(prior_zero_prob_alpha, prior_zero_prob_beta);"                                                   
#  "    }"                                                                                                                  
#  "  }"                                                                                                                    
#  "  // Priors"                                                                                                            
#  "  if (num_edges > 0) {"                                                                                                 
#  "    if (partial_pooling == 0) {"                                                                                        
#  "      edge_weight ~ normal(prior_edge_mu, prior_edge_sigma);"                                                           
#  "    } else {"                                                                                                           
#  "      edge_weight ~ normal(prior_edge_mu, edge_sigma[1]);"                                                              
#  "      edge_sigma ~ normal(0, prior_edge_sigma);"                                                                        
#  "    }"                                                                                                                  
#  "  }"                                                                                                                    
#  ""                                                                                                                       
#  "  if (num_fixed > 0) {"                                                                                                 
#  "    beta_fixed ~ normal(prior_fixed_mu, prior_fixed_sigma);"                                                            
#  "  }"                                                                                                                    
#  ""                                                                                                                       
#  "  if (num_random > 0) {"                                                                                                
#  "    beta_random ~ normal(random_group_mu[random_group_index], random_group_sigma[random_group_index]);"                 
#  "    // Hyperpriors"                                                                                                     
#  "    random_group_mu ~ normal(prior_random_mean_mu, prior_random_mean_sigma);"                                           
#  "    random_group_sigma ~ normal(0, prior_random_std_sigma);"                                                            
#  "  }"                                                                                                                    
#  "}"                                                                                                                      
#  ""                                                                                                                       
#  "generated quantities {"                                                                                                 
#  "  array[num_rows] int event_pred;"                                                                                      
#  "  event_pred = binomial_rng(divisor, inv_logit(predictor));"                                                            
#  "  vector[num_rows] log_lik;"                                                                                            
#  "  for (i in 1:num_rows) {"                                                                                              
# "    log_lik[i] = binomial_lpmf(event[i] | divisor[i], inv_logit(predictor[i]));"                                        
# "  }"                                                                                                                    
# "}"

prior_parameters <- extract_prior_parameters(priors)
prior_parameters$prior_edge_mu_0 <- -5
prior_parameters$prior_edge_mu_1 <- -2.5
prior_parameters$prior_edge_sigma_0 <- 3
prior_parameters$prior_edge_sigma_1 <- 1.5
prior_parameters <- prior_parameters[c(10:13,3:9)]

motnp_edge_weights_strongpriors <- bison_model_hkm_2priors(
  ( event | duration ) ~ dyad(node_1_id, node_2_id),#*seen_together,
  #( event | duration ) ~ dyad(node_1_id, node_2_id) + seen_together*normal(-2.5,1.5) + never_together*normal(-5,3),
  data = counts_df_model,
  model_type = "binary",
  prior_parameters = prior_parameters
)

### run diagnostic plots
plot_trace(motnp_edge_weights_strongpriors, par_ids = 2)                             # trace plot
plot_predictions(motnp_edge_weights_strongpriors, num_draws = 20, type = "density")  # compare predictions to raw data -- predictions are more variable than observations, both overestimating the number of pairs only seen together a few times, but also allowing for together scores higher than observed
plot_predictions(motnp_edge_weights_strongpriors, num_draws = 20, type = "point")    # compare predictions to raw data -- anything below the line is predicting below SRI (e.g. upper end SHOULD be massively below line because we don't want scores of 1 from pairs seen once)

### extract edge weight summaries
edgelist <- get_edgelist(motnp_edge_weights_strongpriors, ci = 0.9, transform = TRUE)  # extract edge list from model (distribution summary for all dyads)
plot(density(edgelist$median))           # median edge strength estimate for all pairs -- mostly very low
summary(motnp_edge_weights_strongpriors)

### compare edge weights to prior
edges <- as.data.frame(motnp_edge_weights_strongpriors$chain) %>%               # extract chain of values from model
  pivot_longer(cols = everything(), names_to = 'dyad', values_to = 'draw') %>%  # convert to long format
  mutate(dyad_id = rep(counts_df$dyad_id, 4000),                                # add column to identify dyad per draw
         draw = plogis(draw))                                                   # convert draws to proportion
head(edges)
plot(NULL, xlim = c(0,1), ylim = c(0,50), las = 1,                              # prepare plot window
     main = 'edge distributions',ylab = 'density', xlab = 'edge weight')
plot_samples <- sample(counts_df$dyad_id, 10000, replace = F)                   # sample 10000 dyads to plot
for(dyad in plot_samples) {                                                     # plot probability density for sampled dyads
  x <- edges[edges$dyad_id == dyad,]
  lines(density(x$draw),
        col = ifelse(mean(x$draw < 0.15), rgb(0,0,1,0.1),
                     ifelse(mean(x$draw < 0.3), rgb(1,0,0,0.1),rgb(1,0,1,0.1))))
}
lines(density(edges$draw), lwd = 3)                                             # add probability density line for all draws from all dyads

edges$chain_position <- rep(1:4000, each = length(unique(edges$dyad)))          # add column for position in chain
edges <- edges[edges$chain_position < 1001,]                                    # only take first 1000 values to speed up next plot
edges$mean <- NA
for(dyad in unique(edges$dyad)) {                                               # calculate mean edge per dyad
  if(is.na(edges$mean[dyad]) == TRUE){
  edges$mean[edges$dyad == dyad] <- mean(edges$draw[edges$dyad == dyad])
  }
}
plot_low <- edges[edges$mean == min(edges$mean),]                               # dyad with minimum mean edge weight
plot_q25 <- edges[edges$mean == quantile(edges$mean, 0.25),]                    # dyad with 25th percentile mean edge weight
plot_mid <- edges[edges$mean == quantile(edges$mean, 0.501),]                   # dyad with median mean edge weight (0.5 exactly couldn't identify a single dyad, but 0.501 found it)
plot_q75 <- edges[edges$mean == quantile(edges$mean, 0.75),]                    # dyad with 75th percentile mean edge weight
plot_q98 <- edges[edges$mean == quantile(edges$mean, 0.98),]                    # dyad with 98th percentile mean edge weight
plot_high <- edges[edges$mean == max(edges$mean),]                              # dyad with maximum mean edge weight
plot_data <- rbind(plot_low, plot_q25, plot_mid, plot_q75, plot_q98, plot_high) # combine all to single data frame

ggplot(plot_data,                                                               # plot distributions of selected dyads
       aes(x = draw, colour = as.factor(mean), fill = as.factor(mean)))+
  geom_density(linewidth = 1)+
  theme_classic()+
  scale_fill_viridis_d(alpha = 0.2)+                                            # colour distributions by dyad (translucent)
  scale_color_viridis_d()+                                                      # colour distributions by dyad (solid)
  theme(legend.position = 'none',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_x_continuous(name = 'edge weight', limits = c(0,1))

#priors$edge <- 'normal(-2.5, 1.5)'                                              # set prior value (already set -- doesn't actually change anything, just in as a reminder of what it was)
prior_plot <- data.frame(dyad = 'prior',                                        # create prior dataframe with draws from prior distribution, formatted to match columns in plot_data
                         draw = plogis(rnorm(1000, -2.5, 1.5)),
                         dyad_id = 1000000,
                         mean = NA,
                         chain_position = 1:1000)
prior_plot$mean <- mean(prior_plot$draw)                                        # calculate mean of edge prior distribution
ggplot(plot_data,                                                               # plot distributions of selected dyads
       aes(x = draw, colour = as.factor(mean), fill = as.factor(mean)),
       #fill = viridis, colour = colours
       )+
  geom_density(linewidth = 1, alpha = 0.2)+
  geom_density(data = prior_plot, linewidth = 1, colour = 'black', linetype = 2, alpha = 0)+  # add probability density line for prior distribution -- no fill, dashed line not solid, different colour scale to others
  theme_classic()+
  scale_fill_viridis_d(option = 'plasma',                                       # use alternative viridis colour pallette to avoid confusion regarding ages (all other plots, viridis scale = age)
                       alpha = 0.2)+
  scale_color_viridis_d(
    option = 'plasma'
    )+
  theme(legend.position = 'none',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_x_continuous(name = 'edge weight', limits = c(0,1))

plot_chains <- rbind(plot_high, plot_mid, plot_low)                             # make data set for trace plots
plot_chains$order <- factor(plot_chains$dyad_id,
                            levels = c('73914','86520','85862'))                # reorder to colour by chain so matches plot above, and also has narrowest chain on top and most uncertain at the back
ggplot(#data = plot_chains,
       #aes(y = draw, x = chain_position, colour = order)
  )+
  geom_line(data = plot_high, aes(y = draw, x = chain_position),                # chain for most uncertain
            #colour = '#fde725',                                                # yellow --  viridis scale B
            colour = '#f0f921'                                                  # yellow -- viridis scale 'plasma'
            )+
  geom_line(data = plot_mid, aes(y = draw, x = chain_position),                 # chain for median
            #colour = '#21918c',                                                # turquoise -- viridis scale B
            colour = '#e16462'                                                  # orange -- viridis scale 'plasma'
            )+
  geom_line(data = plot_low, aes(y = draw, x = chain_position),                 # chain for minimum edge weight
            #colour = '#440154',                                                # purple -- viridis scale B
            colour = '#0d0887'                                                  # red -- viridis scale 'plasma'
            )+
  #scale_color_viridis_d()+
  theme_classic()+
  theme(#legend.position = 'none',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_x_continuous(name = 'chain position')+
  scale_y_continuous(name = 'chain position',
                     limits = c(0,1))

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
motnp_edges_null_strongpriors <- bison_model(
  (event | duration) ~ 1,                       # response variable does not vary with dyad, all dyads drawn from the same distribution
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = motnp_edge_weights_strongpriors, random_model = motnp_edges_null_strongpriors)) # compare fit for model allowed to vary by dyad vs model that draws all dyad strengths from the same distribution -- vast majority of network is best explained by random model, but 4.3% of dyads better explained by non-random. Network is therefore non-random, but with only a small proportion of dyads associating more strongly than random distribution will allow

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging <- function(models) {
  loos <- lapply(models, function(model) loo::loo(model$log_lik, r_eff=NA))
  
  results_matrix <- loo::loo_model_weights(loos, method="pseudobma")
  if (!is.null(names(models))) {
    names(results_matrix) <- names(models)
  }
  results_matrix
}  # produce alternative method for comparing models
model_averaging(models = list(non_random_model = motnp_edge_weights_strongpriors, random_model = motnp_edges_null_strongpriors))            # 100% confidence that random model is better

# save workspace image
save.image('motnp_bisonr_edgescalculated_strongprior.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## plot network ####
### adapt bisonR plot_network function to give more flexibility over plotting options
plot_network_threshold <- function (obj, ci = 0.9, lwd = 2, threshold = 0.3,  # threshold = median edge weight must be over X
                                    vertex.label.color1 = 'transparent', vertex.label.font1 = NA, 
                                    vertex.color1 = 'transparent',
                                    vertex.size1 = 1, edge.color1 = rgb(0, 0, 0, 1),
                                    vertex.label.color2 = 'black', vertex.label.font2 = 'Helvetica', 
                                    vertex.color2 = 'seagreen1',
                                    vertex.size2 = 8, edge.color2 = rgb(0, 0, 0, 0.3))
{
  edgelist <- get_edgelist(obj, ci = ci, transform = TRUE)
  net <- igraph::graph_from_edgelist(as.matrix(edgelist[, 1:2]), 
                                     directed = obj$directed)
  md <- edgelist[, 3]
  ub <- edgelist[, 5]
  coords <- igraph::layout_nicely(net)
  igraph::plot.igraph(net, layout = coords,
                      edge.width = ifelse(md < threshold, 0, md * lwd),
                      vertex.label.color = vertex.label.color1, vertex.color = vertex.color1,
                      vertex.size = vertex.size1,
                      edge.color = rgb(0, 0, 0, 1), edge.arrow.size = 0)
  igraph::plot.igraph(net, layout = coords,
                      edge.width = ifelse(md < threshold, 0, ub * lwd),
                      vertex.label.color = vertex.label.color2, vertex.color = vertex.color2,
                      vertex.size = vertex.size2,
                      edge.color = edge.color1, edge.arrow.size = 0, 
                      add = TRUE)
  if (obj$directed) {
    igraph::plot.igraph(net, edge.width = md * 0, layout = coords, 
                        vertex.label.color = vertex.label.color2, vertex.color = vertex.color2, 
                        edge.color = edge.color2, add = TRUE)
  }
}

### load in ages 
motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

### create nodes data frame
nodes <- data.frame(id = sort(unique(c(counts_df$id_1,counts_df$id_2))),  # all unique individuals
                    node = NA, age = NA, sightings = NA)                  # data needed on each
for(i in 1:nrow(nodes)){
  # add age data
  x <- motnp_ages[motnp_ages$id == nodes$id[i],]                          # vector of ages per individual
  nodes$age[i] <- mean(x$age)                                             # for initial test purposes, age = mean only
  # add node id and count data
  if(nodes$id[i] != 'M99') { y <- counts_df[counts_df$id_1 == nodes$id[i], c('id_1','count_1','node_1_males')] }   # add data for M99 -- only male who is only present in counts_df in ID2, all others are ID1
  else { y <- counts_df[counts_df$id_2 == nodes$id[i], c('id_2','count_2','node_2_males')] }  # add data for all other males
  nodes$sightings[i] <- y[1,2]
  nodes$node[i] <- y[1,3]
}
nodes$sightings <- as.numeric(nodes$sightings)   # convert to numeric
nodes$node <- as.character(nodes$node)           # convert to character
str(nodes)                                       # check structure

### plot network
plot_network_threshold(motnp_edge_weights_strongpriors, lwd = 15, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')

### adapt to remove unconnected nodes and simplify input requirements
plot_network_threshold2 <- function (obj, ci = 0.95, lwd = 2, threshold = 0.3,
                                     label.colour = 'transparent', label.font = 'Helvetica', 
                                     node.size = 4, node.colour = 'seagreen1',
                                     link.colour1 = 'black', link.colour2 = rgb(0, 0, 0, 0.3))
{
  edgelist <- get_edgelist(obj, ci = ci, transform = TRUE)
  threshold_edges <- edgelist[edgelist$median >= threshold,]
  net <- igraph::graph_from_edgelist(as.matrix(threshold_edges[, 1:2]))
  
  if(is.data.frame(node.size) == TRUE ) {
    nodes_list <- data.frame(node = as.numeric(names(net[1])),
                             sightings = NA)
    for(i in 1:nrow(nodes_list)){
      nodes_list$sightings[i] <- nodes$sightings[which(nodes$node == nodes_list$node[i])]
    }
    node_sightings <- log(nodes_list$sightings)*5
  } else { node_sightings <- node.size }
  
  if(is.data.frame(node.colour) == TRUE ) {
    nodes_list <- data.frame(node = as.numeric(names(net[1])),
                             age = NA)
    for(i in 1:nrow(nodes_list)){
      nodes_list$age[i] <- nodes$age[which(nodes$node == nodes_list$node[i])]
    }
    node_age <- nodes_list$age
  } else { node_age <- node.colour }
  
  md <- threshold_edges[, 3]
  ub <- threshold_edges[, 5]
  coords <- igraph::layout_nicely(net)
  igraph::plot.igraph(net, layout = coords,
                      vertex.label.color = ifelse(is.null(label.colour) == TRUE,
                                                  ifelse(node_age < 20, 'black', 'white'),
                                                  label.colour),
                      label.family = label.font,
                      vertex.color = ifelse(node_age < 15, '#FDE725FF',
                                            ifelse(node_age < 20, '#55C667FF',
                                                   ifelse(node_age < 30, '#1F968BFF', 
                                                          ifelse(node_age < 40, '#39568CFF', '#440154FF')))), 
                      vertex.size = node_sightings,
                      frame.color = NA, frame.width = 0,
                      edge.color = NA, edge.arrow.size = 0, edge.width = 0)
  igraph::plot.igraph(net, layout = coords, add = TRUE,
                      vertex.label = NA, vertex.color = 'transparent', vertex.size = 0, 
                      frame.color = NA, frame.width = 0,
                      edge.color = link.colour1, edge.arrow.size = 0, edge.width = md * lwd)
  igraph::plot.igraph(net, layout = coords, add = TRUE,
                      vertex.label = NA, vertex.color = 'transparent', vertex.size = 0, 
                      frame.color = NA, frame.width = 0,
                      edge.color = link.colour2, edge.arrow.size = 0, edge.width = ub * lwd)
}

### plot network
plot_network_threshold2(obj = motnp_edge_weights_strongpriors, threshold = 0.15,
                        node.size = nodes, node.colour = nodes, lwd = 15)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

## compare edge weight distributions to simple SRI ####
# load in workspace image with edge weights already calculated
#load('motnp_bisonr_edgescalculated_strongprior.RData')

### plot against SRI
head(edgelist)                                                                      # check structure of edgelist
colnames(edgelist)[1:2] <- colnames(counts_df_model)[1:2]                           # rename for join function
edgelist$node_1_id <- as.integer(edgelist$node_1_id)                                # convert to integers
edgelist$node_2_id <- as.integer(edgelist$node_2_id)                                # convert to integers
summary <- left_join(edgelist, counts_df_model, by = c('node_1_id','node_2_id'))    # combine distribution data with raw counts
counts <- counts_df[,c('node_1_males','node_2_males','count_1','count_2')]          # add individual count data
colnames(counts)[1:2] <- c('node_1_id','node_2_id')                                 # rename for join function
summary <- left_join(summary, counts, by = c('node_1_id','node_2_id'))              # combine summary data with individual counts
summary$sri <- summary$event / (summary$duration)                                   # calculate basic SRI

plot(density(summary$sri), main = 'SRI vs model output: blue=all,\nred=both seen 8 times, green=both 12 times') # plot SRI distribution
lines(density(summary$median), col = 'blue')                                                                    # distribution of median estimates
lines(density(summary$median[which(summary$count_1 >= 8 & summary$count_2 >= 8)]), col = 'red')                 # distribution of median estimates for frequently sighted males
lines(density(summary$median[which(summary$count_1 >= 12 & summary$count_2 >= 12)]), col = 'green')             # distribution of median estimates for frequently sighted males

# try plotting a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
dyads <- counts_df[,c('dyad_id','node_1_males','node_2_males')]                        # create data frame of dyads
colnames(dyads) <- c('dyad_id','node_1_id','node_2_id')                                # rename for joining
dyads <- left_join(dyads, counts_df_model, by = c('node_1_id','node_2_id'))            # join with counts data (creates counts_df_model but with dyad_id)
length(which(is.na(dyads$duration) == TRUE))

draws <- as.data.frame(motnp_edge_weights_strongpriors$chain) %>%                      # generate dataframe of draws per dyad
  pivot_longer(cols = everything(), names_to = 'dyad_model', values_to = 'edge')       # convert to a long format
draws$dyad_id <- rep(counts_df$dyad_id, 4000)                                          # add dyad_id data
draws$weight <- gtools::inv.logit(draws$edge)                                          # convert to 0-1 bounded values
draws$draw <- rep(1:4000,  each = nrow(counts_df_model))                               # add chain position number
draws <- left_join(draws, dyads, by = 'dyad_id')                                       # combine to add node data
draws$sri <- draws$event / draws$duration                                              # calculate basic SRI to compare with distributions

set.seed(15)                                                                           # make subsetting reproducible
subset_draws <- draws[draws$dyad_id %in% sample(draws$dyad_id, 150, replace = F),]     # sample 150 dyads at random to plot (takes a very long time)
subset_draws$median <- NA ; for(i in 1:nrow(subset_draws)){                            # set up for loop
  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]                  # generate dataset cut down to just the dyad in question
  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), subset_draws$median[i]) }   # record median draw value
head(subset_draws)                                                                     # check structure of data frame
which(is.na(subset_draws$median) == TRUE)[1]                                           # check that all are filled in

subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)           # order dataframe based on total sightings per pair
ggplot(data = subset_draws, mapping = aes(x = weight))+                                # plot dataframe
  geom_density(colour = 'blue')+                                                       # plot density plots of probability distirbtuions output per dyad
  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+                    # split by dyad, allow to vary height depending on dyad
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+            # add line showing where the median estimate is
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')                          # add line showing where the SRI value is

write_csv(subset_draws, '../data_processed/motnp_sampledyads_random_binary_vs_sri_strongprior.csv')    # save output for future reference

subset_draws <- draws[draws$sri > 0.2,]                                                # repeat subsetting but selecting only the dyads with the highest SRI values (generally some of the least sighted individuals)
subset_draws$median <- NA ; for(i in 1:nrow(subset_draws)){                            # set up for loop
  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]                  # generate dataset cut down to just the dyad in question
  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), subset_draws$median[i]) }   # record median draw value
head(subset_draws)                                                                     # check structure of data frame
which(is.na(subset_draws$median) == TRUE)[1]                                           # check that all are filled in

subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)           # order dataframe based on total sightings per pair
ggplot(data = subset_draws, mapping = aes(x = weight))+                                # plot dataframe
  geom_density(colour = 'blue')+                                                       # plot density plots of probability distirbtuions output per dyad
  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15)+                                       # split by dyad, allow to vary height depending on dyad
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+            # add line showing where the median estimate is
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')                          # add line showing where the SRI value is

ggplot(data = subset_draws, mapping = aes(x = weight))+                                # plot dataframe
  geom_density(colour = 'blue')+                                                       # plot density plots of probability distirbtuions output per dyad
  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+                    # split by dyad, allow to vary height depending on dyad
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+            # add line showing where the median estimate is
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')                          # add line showing where the SRI value is

write_csv(subset_draws, '../data_processed/motnp_sampledyads_sri0.2_binary_vs_sri_strongprior.csv')    # save output for future reference

# clean environment
rm(draws, dyads, priors, subset_draws, x) ; gc()
dev.off()

## coefficient of variation of edge weights (aka social differentiation) ####
#load('motnp_bisonr_edgescalculated_strongprior.RData')

# add pdf output file
pdf(file = '../outputs/motnp_bisonr_edgeweight_strongprior_cv.pdf')

# extract cv for model
global_cv_strongprior <- extract_metric(motnp_edge_weights_strongpriors, "global_cv", num_draws = 10000)
head(global_cv_strongprior)
hist(global_cv_strongprior)

# calculate SRI for all dyads
counts_df$sri <- counts_df$event_count / counts_df$count_dyad
summary(counts_df$sri)

# calculate CV of SRI for all dyads
raster::cv(counts_df$sri)         # very high, but lower than gbi_matrix
#raster::cv(m$sri[m$sri > 0])     # still massive even when I remove the 0s -- zero inflation is real

### create SRI matrix
# generate matrix
N <- length(unique(c(counts_df$id_1, counts_df$id_2)))      # number of elephants
ids <- unique(c(counts_df$id_1, counts_df$id_2))            # IDs of elephants
m_mat <- diag(nrow = N)                                     # matrix of males -- NxN to fill with SRI
colnames(m_mat) <- ids ; rownames(m_mat) <- ids             # dimnames = elephant IDs

# populate matrix with SRI values
for( i in 1:N ) {                                           # rows
  for( j in 1:N ) {                                         # columns
    if( i >= j ) { m_mat[i,j] <- m_mat[j,i] }               # make symmetrical about diagonal
    else {
      id1 <- colnames(m_mat)[i]                             # identify elephant 1 of dyad
      id2 <- rownames(m_mat)[j]                             # identify elephant 2 of dyad
      m_mat[i,j] <- counts_df$sri[which(counts_df$id_1 == id1 & counts_df$id_2 == id2)]      # add SRI values -- match IDs and dyad
    }
  }
}

# calculate CV for matrix
raster::cv(m_mat)                # still way too high
which(is.nan(m_mat) == TRUE)     # check if any are blank
sd(m_mat)/mean(m_mat)            # matches raster version -- check it is doing what I think it is doing

# create gbi_matrix
eles <- read_csv(file = '../data_processed/motnp_eles_long.csv')       # read in long format of raw data
eles$location <- paste(eles$gps_s, eles$gps_e, sep = '_')              # make single variable for unique locations
eles <- eles[,c(1,16,2,3,17,4,5,14,7,8,10,13)]                         # rearrange variables
nodes_data <- read_csv(file = '../data_processed/motnp_elenodes.csv' ) # read in node data
colnames(nodes_data)
str(nodes_data)
eles_asnipe <- eles[,c(3,4,2,5)]                                       # date, time, elephant, location
eles_asnipe$Date <- as.integer(eles_asnipe$date)                       # make date numeric
eles_asnipe$Date <- 1+eles_asnipe$Date - min(eles_asnipe$Date)         # start from 1, not 1st January 1970
eles_asnipe$Time <- ifelse(eles_asnipe$time > 1, NA, eles_asnipe$time) # time = proportion of day so anything >1 has to be wrong
eles_asnipe$Time <- eles_asnipe$Time*(24*60*60)                        # convert time values to seconds through day
which(is.na(eles_asnipe$Time))                                         # 161 698 1122 1469 1770
eles_asnipe[c(161,698,1122,1469,1770),]                                # all 1 sighting of B7+M44
eles_asnipe <- eles_asnipe[,c(5,6,3,4)]                                # create data frame to produce gbi matrix from
colnames(eles_asnipe) <- c('Date','Time','ID','Location')              # rename variables for get_gbi
eles_asnipe$ID <- as.character(eles_asnipe$ID)                         # correct data type
str(eles_asnipe)
eles_asnipe$d_pad <- str_pad(eles_asnipe$Date, 3, pad = '0')  # 0-pad dates
eles_asnipe$encounter <- paste(eles_asnipe$d_pad, eles_asnipe$Time, eles_asnipe$Location, sep = '_') # unique value for each sighting
eles_asnipe$group <- as.integer(as.factor(eles_asnipe$encounter)) # unique factor for every sighting
max(eles_asnipe$group)                         # 574 -- agrees with number of different sightings for which elephants were identified
eles_dt <- eles_asnipe[,c(3,7)]                # create data table for gbi matrix
eles_dt <- data.table::setDT(eles_dt)          # create data table for gbi matrix
gbi_matrix <- spatsoc::get_gbi(DT = eles_dt, group = 'group', id = 'ID')  # create gbi matrix
gbi_males <- gbi_matrix[,colnames(gbi_matrix) %in% ids]
gbi_males <- gbi_males[rowSums(gbi_males) > 0,] # remove male-only sightings

# set up permutations
N_networks <- 10000
# rm(counts_df_model, edgelist, eles_asnipe, eles_dt, females_df, gbi_matrix, motnp_ages, motnp_edge_weights_strongpriors, motnp_edges_null_strongpriors, i, id1, id2, j, model_averaging) ; gc()

# create vector of days for each sighting
sightings <- eles[,c('elephant','encounter','date')] %>%  # set up dataframe of actual encounters
  filter(elephant %in% ids) %>%                           # remove encounters that only include females, youngsters, or dead males
  select(-elephant) %>%                                   # remove ID column
  distinct()                                              # cut down to one row per encounter

# run network permutations
random_networks <- asnipe::network_permutation(association_data = gbi_males, # permute network -- raw data = gbi_matrix
                                               association_matrix = m_mat,   # SRI matrix
                                               permutations = N_networks,    # 10000 times
                                               days = sightings$date,        # permute within day only
                                               within_day = TRUE)            # permute within day only
print('random networks generated')

# calculate coefficient of variation per randomised network
cv_random_networks <- rep(0,N_networks) ; for (i in c(1:N_networks)) {       # set up for loop
  net_rand <- random_networks[i,,]                                           # select random network
  cv_random_networks[i] <- raster::cv(net_rand)                              # calculate CV and save into vector
  if(i %% 1000 == 0) {print(i)}                                              # track progress
}
print(paste0('network permutations for entire network completed at ', Sys.time()))

# compare permuted networks to actual network
hist(cv_random_networks, xlim = c(min(cv_random_networks)-50,max(cv_random_networks)+50),   # plot histogram of random netowkr CVs
     #ylim = c(0,10),  # there are 3 networks out of 10000 that have cv > 290, but then all are in the same bar as cv(m_mat)
     main = 'permutations for network of all male elephants')
abline(v = raster::cv(m_mat), col = 'red')                                                  # add line for measured network
text(round(raster::cv(m_mat),3), col = 'red', x = median(cv_random_networks)+50, y = 2000)  # add text to specify value for measured network

plot_cv <- as.data.frame(cv_random_networks)                # create dataframe
cv_network <- raster::cv(m_mat)                             # save value of measured CV
cv_random_networks <- plot_cv$cv_random_networks            # vector of randomised CV values
ggplot(data = plot_cv)+                                     # plot CV values
  geom_vline(xintercept = cv_network, linewidth = 1.5,      # add line for true value
             colour = rgb(68/255, 1/255, 84/255))+          # dark purple (viridis scale)
  #geom_text(x = cv_network - 50, y = 700,
  geom_histogram(aes(cv_random_networks),                   # random cv values
                 fill = rgb(253/255, 231/255, 37/255),      # yellow (viridis scale)
                 colour = 'black', bins = 100)+             # set bin width
  theme_classic()+
  scale_x_continuous(name = 'coefficient of variation',
                     expand = c(0,0),
                     limits = c(min(cv_random_networks)-10,
                                cv_network+10))+
  scale_y_continuous(name = 'frequency',
                     expand = c(0,0))+
  #          colour = rgb(68/255, 1/255, 84/255),
  #          label = 'coefficient of/nvariation for/nmeasured network')+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))

### write out outputs for future reference
write_csv(plot_cv, '../data_processed/motnp_networkpermutations_cv_strongprior.csv')

### combine with global_cv_strongpriors
plot_cv_model <- data.frame(cv = c(plot_cv$cv_random_networks, (global_cv_strongprior*100)),             # combine permutations and model values
                      iteration = rep(1:length(global_cv_strongprior), 2),
                      type = rep(c('permutation','model_draw'), each = length(global_cv_strongprior)))
write_csv(plot_cv_model, '../data_processed/motnp_networkpermutations_cv_strongprior.csv')               # save for future reference
ggplot()+
  geom_vline(xintercept = cv_network, linewidth = 1.5,                                                   # line for measured network from SRI values
             colour = rgb(68/255, 1/255, 84/255))+                                                       # dark purple (viridis)
  #annotate('text', x = cv_network - 50, y = 700), colour = rgb(68/255, 1/255, 84/255),
  #         label = 'coefficient of/nvariation for/nmeasured network')+
  theme_classic()+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))+
  geom_histogram(data = plot_cv_model[plot_cv_model$type == 'model_draw',], aes(x = cv),                 # histogram of CV draws from model
                 fill = rgb(33/255, 145/255, 140/255),                                                   # yellow (viridis)
                 bins = 100, colour = 'black')+
  geom_histogram(data = plot_cv_model[plot_cv_model$type == 'permutation',], aes(x = cv),                # histogram of CV draws from permutations
                 fill = rgb(253/255, 231/255, 37/255),                                                   # turquoise  (viridis)
                 bins = 100, colour = 'black')+
  scale_x_continuous(name = 'coefficient of variation',
                     #limits = c(min(cv_random_networks)-10,
                     #cv_network+10),
                     expand = c(0,0))+
  scale_y_continuous(name = 'frequency',
                     expand = c(0,0))

# repeat plot for poster -- larger fonts, fewer draws
plot_cv_model2 <- plot_cv_model[c(1:10000, sample(10001:20000, 4000, replace = F)),]                     # select only 4000 values from dataframe
ggplot()+
  geom_vline(xintercept = cv_network, linewidth = 1.5,                                                   # line for measured network from SRI values
             colour = rgb(68/255, 1/255, 84/255))+                                                       # dark purple (viridis)
  #annotate('text', x = cv_network - 50, y = 700), colour = rgb(68/255, 1/255, 84/255),
  #         label = 'coefficient of/nvariation for/nmeasured network')+
  theme_classic()+
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 22))+
  geom_histogram(data = plot_cv_model2[plot_cv_model2$type == 'model_draw',], aes(x = cv),               # histogram of CV draws from model
                 fill = rgb(33/255, 145/255, 140/255),                                                   # yellow (viridis)
                 bins = 100, colour = 'black')+
  geom_histogram(data = plot_cv_model2[plot_cv_model2$type == 'permutation',], aes(x = cv),              # histogram of CV draws from permutations
                 fill = rgb(253/255, 231/255, 37/255),                                                   # turquoise  (viridis)
                 bins = 100, colour = 'black')+
  scale_x_continuous(name = 'coefficient of variation',
                     #limits = c(min(cv_random_networks)-10,
                     #cv_network+10),
                     breaks = round(seq(100, 400, by = 50),-1),
                     expand = c(0,0))+
  scale_y_continuous(name = 'frequency',
                     expand = c(0,0))

### run statistical test to confirm that difference is significant
# plot_cv <- read_csv('../data_processed/motnp_networkpermutations_cv_strongprior.csv')
hist(plot_cv$cv_random_networks)    # plot cv values from random networks -- 1 tiny bar very high value and outside of the main distribution
cv_network                          # plot line for measured value -- falls over the tiny bar

mean(plot_cv$cv_random_networks)    # mean is far lower than measured value
sd(plot_cv$cv_random_networks)      # mean + 2*SD is still far lower than measured value

length(which(plot_cv$cv_random_networks >= cv_network))/length(plot_cv$cv_random_networks) # 3 / 10000 = 3e-04 = p-value (don't need an additional test because that will further randomise the data which are already randomised, just need to know proportion that are greater than or equal to your value)

mean(global_cv) ; sd(global_cv)     # totally outside either the random distribution and also doesn't overlap at all with the SRI measured network
max(global_cv) ; min(plot_cv$cv_random_networks)     # no overlap

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2% ####
counts_df_model <- counts_df[,c('node_1_males','node_2_males','event_count','count_dyad','id_1','id_2','dyad_males')]   # select variables of interest
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')                      # rename variables
ms_chain <- as.data.frame(motnp_edge_weights_strongpriors$chain)                                                        # extract chain values
colnames(ms_chain) <- counts_df_model$dyad_id                                                                           # rename with dyad IDs
ms_chain <- pivot_longer(ms_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')                       # pivot longer
ms_chain$chain_position <- rep(1:length(unique(ms_chain$dyad_id)), each = 4000)                                         # add chain positions
ms_chain$draw <- LaplacesDemon::invlogit(ms_chain$draw)                                                                 # convert values to 0-1 bounded
ms_chain$mean <- NA                                                                                                     # set up for loop for mean values
hist(ms_chain$draw)                                                                                                     # plot histogram of draw values
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior',           # empty plot for density lines
     xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sample(unique(ms_chain$dyad_id), size = 1000, replace = F)){                                                   # set up for loop using 1000 sample dyads
  x <- ms_chain[ms_chain$dyad_id == i,]                                                                                 # select values for dyad weights
  ms_chain$mean <- ifelse(ms_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ms_chain$mean)                             # calculate mean values
  lines(density(x$draw), col = rgb(0,0,1,0.1))                                                                          # plot sample line
}
lines(density(ms_chain$draw), lwd = 2)                                                                                  # plot average line

quantile(ms_chain$draw, 0.98)                                                                                           # identify 98th percentile draw
quantile(edgelist$median, 0.98)                                                                                         # identify 98th percentile median
counts_df_model$sri <- counts_df_model$event / counts_df_model$duration                                                 # calculate SRI values
quantile(counts_df_model$sri, 0.98)                                                                                     # identify 98th percentile SRI

## clean up ####
### save pdf
dev.off()

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

rm(list= ls()[!(ls() %in% c('motnp_edge_weights_strongpriors','motnp_edges_null_strongpriors',
                            'counts_df','counts_df_model',
                            'gbi_males','m_mat',
                            #'random_nodes',
                            'random_networks'))])

### write over saved workspace image -- only keep the things you might need later that take a long time to run
save.image('motnp_bisonr_edgescalculated_strongprior.RData')

