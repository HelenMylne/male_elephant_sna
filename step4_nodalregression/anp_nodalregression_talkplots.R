##### information ####
# 1) to estimate the ages of individuals in other populations, first came up with a survival curve for the Amboseli elephants which have good birth and death data.
# 2) selected the Gompertz bathtub distribution
# 3) use distribution to convert MOTNP ages in categories to probability distributions of true age (first part of this script)
# 4) extract node centrality distribution from draws of edge weight
# 5) use MOTNP age distributions to predict network centrality (third part of this script)

##### set up ####
# library(rstan) ; library(igraph) ; library(tidyverse) ; library(LaplacesDemon) ; library(MASS)
library(rstan, lib.loc = '../packages/')
library(igraph, lib.loc = '../packages/')          # library(igraph)
library(tidyverse, lib.loc = '../packages/')       # library(tidyverse)
library(LaplacesDemon, lib.loc = '../packages/')   # library(LaplacesDemon)
library(MASS, lib.loc = '../packages/')            # library(MASS)

## set seed for reproducibility
set.seed(12345)

##### run model with rstan -- time window 1-36 ####
set.seed(12345)

for(time_window in 1:36){
  # load data and remove additional data
  load(paste0('anp_nodalregression/anpshort',time_window,'_nodalregression_conditionaledge_rstan.RData'))
  rm(list = ls()[! ls() %in% c('time_window','df_long','sim_summary','mod_mu')]) ; gc()
  
  ggplot()+
    geom_ribbon(data = sim_summary, aes(x = age, ymin = lwr, ymax = upr),       # shade simulations
                colour = 'transparent', fill = rgb(0,0,0,0.1))+
    geom_ribbon(data = mod_mu, aes(x = age, ymin = lwr, ymax = upr),            # shade mean distribution
                colour = 'transparent', fill = rgb(33/255, 145/255, 140/255, 0.5))+
    # geom_point(data = df_long, aes(x = age, y = centrality),                    # all eigenvector draws
    #            colour = rgb(253/255, 231/255, 37/255, 0.01))+
    # geom_point(data = df_long, aes(x = age, y = mean_eigen, size = sightings),  # mean eigenvector
    #            colour = rgb(68/255, 1/255, 84/255))+
    geom_line(data = mod_mu, aes(x = age, y = mid),                             # mean line
              colour = rgb(33/255, 145/255, 140/255), linewidth = 1)+
    scale_x_continuous('age (years)')+
    scale_y_continuous('eigenvector centrality (standardised)',
                       limits = c(-3,3))+
    theme_classic()+
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 22),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 22))
  ggsave(filename = paste0('../outputs/anpshort',time_window,'_nodalregression_nopoints.png'), 
         device = 'png',plot = last_plot(), width = 7, height = 8.5)
  
  ## add progress marker
  print(paste0('time window ',time_window,' complete'))
  
}
