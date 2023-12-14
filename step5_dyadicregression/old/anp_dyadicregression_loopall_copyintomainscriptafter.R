#### information ####
# script takes data input from edge weight estimation for ANP population (input = 1000 draws for each of 4 chains per dyad posterior distribution)
# runs through dyadic regression as specified by Jordan Hart in BISoN examples (https://github.com/JHart96/bison_examples/blob/main/examples/dyadic_regression_stan.md)

#### set up ####
#library(tidyverse) ; library(LaplacesDemon) ; library(cmdstanr) ; library(bisonR) ; library(brms)
library(cmdstanr, lib.loc = '../packages/')  # library(cmdstanr)
library(tidyverse, lib.loc = '../packages/') # library(tidyverse)
#library(car, lib.loc = '../packages/')       # library(car)
#library(bisonR, lib.loc = '../packages/')    # library(bisonR)
#library(brms, lib.loc = '../packages/')      # library(brms)
library(LaplacesDemon, lib.loc = '../packages/')

set_cmdstan_path('R:/rsrch/df525/phd/hkm513/packages/.cmdstan/cmdstan-2.31.0')

## load dyadic regression model
dyadic_regression <- cmdstan_model('models/dyadic_regression.stan')

#### loop model ####
for(time_window in 1:36){
  #### import data ####
  load(file = paste0('anp_edgecalculations/anpshort',time_window,'_edgeweights_conditionalprior.RData'))
  # rm(edgelist, edge_binary, nodes) ; gc()
  
  pdf(file = paste0('../outputs/anpshort',time_window,'_dyadicregression_plots.pdf'))
  
  ## progress marker
  print(paste0('time window ',time_window,' data read in at ',Sys.time()))
  
  #### prior predictive check ####
  age_min <- 5:50
  age_max <- seq(10, 50, length.out = 12)
  par(mfrow = c(4,3))
  for(age in age_max){
    plot(NULL, xlim = c(5,50), ylim = c(0,1), las = 1,
         xlab = 'age of younger', ylab = 'edge weight', main = paste0('oldest = ',round(age)))
    for(i in 1:100){
      beta_min <- rnorm(1, 0, 0.1)
      beta_max <- rnorm(1, 0, 0.1)
      #beta_int <- rnorm(1, 0, 0.005)
      min_new <- age_min[age > age_min]
      lines(x = min_new, y = plogis(min_new*beta_min +
                                      #min_new*age*beta_int +
                                      age*beta_max),
            col = rgb(0,0,1,0.2))
    }
  }
  rm(age, age_min, age_max, beta_max, beta_min, i, min_new) ; gc()
  par(mfrow = c(1,1))
  
  # add time marker
  print(paste0('prior predictive check completed at ', Sys.time()))
  
  #### fit multivariate Gaussian distribution to output of edge weight model ####
  ### To parameterise the multivariate normal approximation, we use the sample mean and covariance matrix, calculated from the posterior edge weight samples using the following code:
  ### get the weights on the logit scale (they are not currently because we used a beta prior and identity link here rather than logistic link)
  logit_weights <- apply(edge_samples, 2, logit)
  
  ### fit a multivariate normal dist to the edges -- quantities will be given to Stan model as data to model joint posteriors of edge weight in the regression
  logit_edge_draws_mu <- apply(logit_weights, 2, mean)
  logit_edge_draws_cov <- cov(logit_weights)
  
  #### plot to see how well the approximation is working ####
  ### Randomly selecting samples to examine
  num_check <- 20
  selected_samples <- sample(1:n_dyads, num_check, replace = FALSE)
  
  ### Setting grid layout
  rows <- floor(sqrt(num_check))
  cols <- ceiling(num_check / rows)
  par(mfrow=c(rows, cols), mar=c(2,2,2,1))
  
  ### plot
  for (i in selected_samples) {
    mu <- logit_edge_draws_mu[i]
    sd <- sqrt(logit_edge_draws_cov[i,i])
    
    fitted_values_logit <- rnorm(1e5, mean = mu, sd = sd)
    fitted_values_original <- plogis(fitted_values_logit)
    
    hist(unlist(edge_samples[,i]), probability = TRUE, las = 1,
         main = paste("Dyad", i), xlab = "Value", breaks = 50)
    lines(density(fitted_values_original), col="blue", lw=1.5)
  }
  
  for (i in selected_samples) {
    mu <- logit_edge_draws_mu[i]
    sd <- sqrt(logit_edge_draws_cov[i,i])
    
    fitted_values_logit <- rnorm(1e5, mean = mu, sd = sd)
    fitted_values_original <- plogis(fitted_values_logit)
    
    plot(unlist(edge_samples[,i]), unlist(edge_samples[,i+1]), col = rgb(0,0,1,0.05), las = 1,
         main = paste("cov ", i ,"&",i+1))
  }
  
  ### reset plot window
  par(mfrow=c(1,1))
  
  # save image so far
  save.image(paste0('anpshort',time_window,'_dyadicregression.RData'))
  
  # add time marker
  print(paste0('multivariate Gaussian approximation fitted at ', Sys.time()))
  
  #### fit dyadic regression ####
  ### identify older and younger of dyad
  cdf$age_min <- NA ; cdf$age_max <- NA
  for(i in 1:nrow(cdf)){
    x <- c(cdf$age_start_1[i],cdf$age_start_2[i])
    cdf$age_min[i] <- min(x)
    cdf$age_max[i] <- max(x)
  }
  
  ## convert node IDs to values 1:52 not numbers based on casename
  all_node_IDs <- unique(c(cdf$id_1, cdf$id_2))
  n_nodes <- length(all_node_IDs)
  cdf <- cdf %>% 
    rename(node_1_original = node_1,
           node_2_original = node_2) %>% 
    mutate(node_1_period = as.integer(factor(id_1, levels = all_node_IDs)),
           node_2_period = as.integer(factor(id_2, levels = all_node_IDs)),)
  
  ## create data list
  dyad_data <- list(
    num_dyads = n_dyads,                      # number of dyads
    num_nodes = n_nodes,                      # number of nodes
    logit_edge_mu = logit_edge_draws_mu,      # sample means of the logit edge weights
    logit_edge_cov = logit_edge_draws_cov,    # sample covariance of logit edge weights
    age_min = cdf$age_min,                  # age of younger dyad member
    #age_max = cdf$age_max,                  # age of  older  dyad member
    age_diff = cdf$age_diff,                # age difference between dyad members
    node_1 = cdf$node_1_period,             # node IDs for multimembership effects
    node_2 = cdf$node_2_period              # node IDs for multimembership effects
    #jitter = 1e-6                             # jitter to add to the diag of the cov matrix for numerical stability -- COME BACK TO THIS, MAY NEED TO ALTER THE MODEL TO ADD SIGMA BACK IN AND REMOVE JITTER
  )
  
  ## add time marker
  print(paste0('start model run at ', Sys.time()))
  
  ## fit dyadic regression
  fit_dyadreg_anp <- dyadic_regression$sample(
    data = dyad_data,
    chains = n_chains,
    parallel_chains = n_chains)
  
  # save image so far
  save.image(paste0('anpshort',time_window,'_dyadicregression.RData'))
  
  # add time marker
  print(paste0('finish model run at ', Sys.time()))
  
  #### check outputs ####
  #load(paste0('anpshort',time_window,'_dyadicregression.RData'))
  rm(counts_df, dyad_data, edge_binary, edgelist, edges, i, x, make_edgelist, plot_network_threshold_anp) ; gc()
  
  # obtain summary
  fit_dyadreg_anp$summary()
  
  ## extract draws
  draws <- fit_dyadreg_anp$draws(format = 'df')
  
  ## extract dyadic regression slopes
  b_diff <- draws[,,'beta_age_diff']
  b_min <- draws[,,'beta_age_min']
  sigma <- draws[,,'sigma']
  
  ## extract overall age effects
  dimnames(draws)$variable
  draws <- draws[,,4:length(draws[1,1,])]   # WILL CHANGE IF YOU ADD ANY MORE REGRESSION SLOPES
  
  ## extract multi-membership samples -- CHECK THIS IS IN THE RIGHT ORDER
  mm_matrix <- draws[,,1:n_nodes]
  sigma_mm <- draws[,,n_nodes+2]
  predictor <- draws[,,(n_nodes+3):length(draws[1,1,])]
  length(draws[1,1,]) == length(mm_matrix[1,1,]) + length(sigma_mm[1,1,]) + length(sigma[1,1,]) + length(predictor[1,1,])
  rm(draws) ; gc()
  
  ## save parameter values
  extract_slopes <- function(draws, n_samples = 1000, n_chains = 4){
    draws <- as.data.frame(draws) %>% 
      pivot_longer(everything(), names_to = 'chain', values_to = 'slope_draw') %>% 
      separate(chain, sep = c(1,2), remove = T, into = c('chain','.', 'parameter')) %>% 
      select(-.)
    draws$chain <- as.numeric(draws$chain)
    draws$position <- rep(1:n_samples, each = n_chains)
    draws$draw_id <- draws$position + (draws$chain-1)*n_samples
    return(draws)
  }
  b_diff <- extract_slopes(b_diff)
  b_min <- extract_slopes(b_min)
  sigma <- extract_slopes(sigma)
  parameters <- rbind(b_min, b_diff, #b_max, #b_int,
                      sigma)
  
  ## save data 
  saveRDS(parameters, paste0('../data_processed/step5_dyadicregression/anpshort',time_window,'_dyadicregression_slopeparameters.RDS'))
  
  # add time marker
  print(paste0('parameters saved to file at ', Sys.time()))
  
  pdf(paste0('../outputs/anpshort',time_window,'_dyadicregression_plots.pdf'))
  
  ## traceplots
  ggplot(data = parameters)+
    geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
    theme_classic()+
    theme(legend.position = 'none')+
    scale_colour_viridis_d()+
    facet_wrap( . ~ parameter , scales = 'free_y')
  
  # add time marker
  print(paste0('finish time window ',time_window,' at ',Sys.time()))
  
  #### finish ####
  save.image(paste0('anpshort',time_window,'_dyadicregression.RData'))
  dev.off()
  rm(list = ls()[! ls() %in% c('time_window','dyadic_regression')])
  
}
