#### set up ####
# library(cmdstanr) ; library(tidyverse) ; library(car) ; library(LaplacesDemon)
library(cmdstanr, lib.loc = '../packages/')  # library(cmdstanr)
library(tidyverse, lib.loc = '../packages/') # library(tidyverse)
library(car, lib.loc = '../packages/')       # library(car)
library(LaplacesDemon, lib.loc = '../packages/')

set_cmdstan_path('../packages/.cmdstan2/cmdstan-2.33.1/')

theme_set(theme_bw())

# #### import data for window 5 ####
# load('anp_edgecalculations/anpshort5_edgeweights_conditionalprior.RData')
# 
# edges <- edges %>% 
#   mutate(dyad_id = as.numeric(dyad_id)) %>% 
#   left_join(cdf, by = 'dyad_id') %>% 
#   mutate(together = ifelse(event_count == 0, 'never together', 'together at least once'))
# 
# sample <- sample(unique(edges$dyad_id), size = 100, replace = F)
# 
# edges_sample <- edges %>% 
#   filter(dyad_id %in% sample)
# 
# #### plot posterior distribution for window 5 ####
# (posterior <- ggplot(edges_sample)+
#   stat_density(aes(x = edge_draw,
#                    colour = as.factor(together),
#                    group = dyad_id),
#                geom = 'line',
#                position = 'identity',
#                linewidth = 0.2)+
#   scale_colour_viridis_d(end = 0, begin = 0.5)+
#   labs(x = 'edge weight', y = 'density',
#        colour = 'sightings together')+
#   scale_x_continuous(limits = c(0,1)))
# ggsave(posterior, device = 'png',
#        width = 1200, height = 600, units= 'px',
#        filename = 'anpshort5_posterior_sample.png',
#        path = '../outputs/step3_edgeweightestimation/')
# 
#### summarise posterior distribution medians ####
load('step5_dyadicregression/anpshort_dyadicregression_dataimported.RData')

edge_medians <- edge_summary$median

quantile(edge_medians, prob = c(0.025,0.05,0.1, 0.2, 0.25, 0.3, 0.4, 
                                0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.975))
max(edge_medians)
min(edge_medians)

