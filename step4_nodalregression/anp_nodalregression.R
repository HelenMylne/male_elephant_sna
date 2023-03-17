#### Information ####
# 1) to estimate the ages of individuals in other populations, first came up with a survival curve for the Amboseli elephants which have good birth and death data.
# 2) selected the Gompertz bathtub distribution
# 3) use distribution to convert MOTNP ages in categories to probability distributions of true age (first part of this script)
# 4) extract node centrality distribution from draws of edge weight
# 5) use MOTNP age distributions to predict network centrality (third part of this script)

#### set up ####
options(future.globals.maxSize = 10000*(1024^2))   # running model with full age distributions = 5.01GB of globals, which exceeds default maximum allowed size. Set to allow up to 10 GB to allow model to run

library(tidyverse, lib.loc = '../packages/')       # library(tidyverse)
library(cmdstanr, lib.loc = '../packages/')        # library(cmdstanr)
library(brms, lib.loc = '../packages/')            # library(brms)
library(Rcpp, lib.loc = '../packages/')            # library(Rcpp)
library(ggdist, lib.loc = '../packages/')          # library(ggdist)
library(posterior, lib.loc = '../packages/')       # library(posterior)
library(bayesplot, lib.loc = '../packages/')       # library(bayesplot)
library(igraph, lib.loc = '../packages/')          # library(igraph)
library(LaplacesDemon, lib.loc = '../packages/')   # library(LaplacesDemon)
library(bisonR, lib.loc = '../packages/')          # library(bisonR)
library(janitor, lib.loc = '../packages/')         # library(janitor)

#### prior predictive check ####
priors <- get_default_priors('binary')
priors$fixed
prior_check(priors, 'binary')

age <- 1:60
beta_mu <- 0
beta_sigma <- 0.005

load('anp_edgecalculations/anpshort1_bisonr_edgescalculated.RData')
rm(counts_df_model, counts_df_allwindows, df, edgelist, ew_chain, ew_edgelist, global_cv, anp_edges_null, draw98, i, ids, median98, sri98, plot_network_threshold, plot_network_threshold2) ; gc()

mean_age <- mean(nodes$age)
plot(NULL, xlim = c(10,60), ylim = c(0,1), las = 1,
     xlab = 'age', ylab = 'eigenvector centrality')
for(i in 1:100){
  intercept <- rbeta(1,2,2)   # this isn't right but I don't think there is a prior for the intercept?? I've gone for a symmetrical one here that in itself explores most of the parameter space and allows it to see whether some of the lines are steep enough to go from top to bottom, but on the assumption that when combined, they will explore only the space relevant to their starting position
  beta <- rnorm(1, beta_mu, beta_sigma)
  lines(x = age, y = intercept + (age - mean_age)*beta, col = rgb(0,0,1,0.5)) # vast majority come out somewhere sensible, and those that don't would if they started at a different value for age 10 so that comes down to my ability to work out what the intercept prior should actually be -- think this is a good prior for the slope
}

#### nodal regression -- only one age value ####
# load edge weight model and data frames
for( time_window in 1:36) {
  
  ## load work space from calculating edge weights
  load(file = paste0('anp_edgecalculations/anpshort',time_window,'_bisonr_edgescalculated.RData'))
  rm(counts_df_model, counts_df_allwindows, df, edgelist, ew_chain, ew_edgelist, global_cv, anp_edges_null, draw98, i, ids, median98, sri98, plot_network_threshold, plot_network_threshold2) ; gc()
  
  ## create nodal data frame from counts_df (data frame of dyad counts of sightings together vs apart)
  df_nodal <- counts_df[,c('node_1','id_1')] %>% distinct()
  colnames(df_nodal) <- c('node_2','id_2')
  df_nodal <- rbind(df_nodal, counts_df[,c('node_2','id_2')]) %>% distinct()
  colnames(df_nodal) <- c('node','id')
  colnames(nodes)[1] <- 'node'
  df_nodal <- left_join(df_nodal, nodes, by = 'node')
  
  # define PDF output
  pdf(paste0('../outputs/anpshort',time_window,'_nodalregression_plots_meanage.pdf'))

  ## eigenvector only ####
  anp_edge_weights$model_data$prior_fixed_mu    <- beta_mu
  anp_edge_weights$model_data$prior_fixed_sigma <- beta_sigma
  
  ## run model ####
  anp_eigen <- bison_brm(
    bison(node_eigen(node)) ~ age,
    anp_edge_weights,
    df_nodal,
    chains = 4,
    iter = 10000,
    thin = 2
  )
  summary(anp_eigen) # FIT 100 IMPUTED MODELS, 4 CHAINS PER MODEL, EACH 1000 DRAWS LONG (+1000 DRAWS WARM UP). WARNING AT END OF MODEL RUN THAT CHAINS <3 DRAWS LONG AS ACTUAL CHAIN IS ONLY 1 WARMUP AND 1 SAMPLE. ONLY IMPUTED CHAINS FOLLOW THE SPECIFIED ITERATIONS AND THINNING
  
  ## save output
  save.image(paste0('anp_edgecalculations/anpshort',time_window,'_nodalregression.RData'))
  
  ## posterior check ####
  summary(anp_eigen)
  
  # plot
  eigen_values <- anp_eigen$data
  plot(eigen_values$bison_node_eigen ~ eigen_values$age, 
       las = 1, pch = 19, col = rgb(0,0,1,0.2),
       xlab = 'mean age estimate', ylab = 'eigenvector centrality',
       main = 'effect of age on eigenvector centrality')
  
  hist(anp_eigen$rhats[,2], las = 1, main = 'Rhat values for 100 imputed model runs', xlab = 'Rhat')
  
  post_eigen <- as.data.frame(as_draws_df(anp_eigen)) %>% clean_names()
  plot(data = post_eigen[post_eigen$chain == 1,], b_age ~ draw, type = 'l', xlim = c(0,10000))
  lines(data = post_eigen[post_eigen$chain == 2,], b_age ~ draw, col = 'red')
  lines(data = post_eigen[post_eigen$chain == 3,], b_age ~ draw, col = 'blue')
  lines(data = post_eigen[post_eigen$chain == 4,], b_age ~ draw, col = 'green')
  lines(data = post_eigen[post_eigen$chain == 5,], b_age ~ draw, col = 'purple')
  unique(post_eigen$chain)
  
  post_eigen$imputation <- rep(1:100, each = 10000)
  post_eigen$chain_imp  <- as.factor(rep(1:4, each = length(post_eigen$chain[post_eigen$chain == 1])))
  ggplot(post_eigen, aes(x = iteration, y = b_age, colour = chain_imp))+
    geom_line()+
    facet_wrap(imputation ~ .)   # looks good, well mixed
  
  hist(post_eigen$b_age)         # natural scale?
  #hist(plogis(post_eigen$b_age)) # logit scale?
  
  #plot(anp_eigen)
  
  plot(conditional_effects(anp_eigen), points = TRUE) # older elephants have lower network centrality than younger
  
  if(length(which(eigen_values$age != nodes$age)) == 0 ){
    data <- cbind(eigen_values, nodes[,c('node','sightings')])
    
    ggplot(data = data, aes(x = age, y = bison_node_eigen,
                            colour = age))+
      geom_point(aes(size = sightings))+
      geom_smooth()+
      scale_x_continuous(name = 'age')+
      scale_y_continuous(name = 'mean eigenvector centrality')+
      scale_colour_viridis_c()+
      theme_classic()+
      theme(legend.position = 'none', axis.title = element_text(size = 18),
            axis.text = element_text(size = 14))
  }
  
  # compare empirical distribution to posterior predictive distribution
  y <- eigen_values$bison_node_eigen                         # extract eigenvector centralities
  yrep <- posterior_predict(anp_eigen, draws = 500)          # make predictions of eigenvector centtrality
  dim(yrep)
  ppc_dens_overlay(y, yrep[1:1000, ])                        # plot 1000 predictions over empirical distribution (is it ok that this slightly extends beyond x = 1??)
  ppc_hist(y, yrep[1:55, ])                                  # compare 55 predictions to empirical distribution
  ppc_pit_ecdf(y, yrep[1:50,])                               # no idea....
  ppc_ecdf_overlay(y, yrep[1:50,]) + xaxis_text()            # also no idea...
  
  ## save output
  save.image(paste0('anp_edgecalculations/anpshort',time_window,'_nodalregression.RData'))
  dev.off()
  
}

