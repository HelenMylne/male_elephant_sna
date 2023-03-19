#### set up ####
# load packages
library(tidyverse) #, lib.loc = '../packages/') # library(tidyverse)
library(dplyr, lib.loc = '../packages/')     # library(dplyr)
#library(rstan, lib.loc = '../packages/')     # library(rstan)
library(cmdstanr, lib.loc = '../packages/')  # library(cmdstanr)
library(bisonR, lib.loc = '../packages/')    # library(bisonR)
library(asnipe, lib.loc = '../packages/')    # library(asnipe)
library(sna, lib.loc = '../packages/')       # library(sna)
library(raster, lib.loc = '../packages/')    # library(raster)
library(janitor, lib.loc = '../packages/')   # library(janitor)
library(gtools, lib.loc = '../packages/')    # library(gtools)
library(igraph, lib.loc = '../packages/')    # library(igraph)
library(loo, lib.loc = '../packages/')       # library(loo)
library(readxl, lib.loc = '../packages/')    # library(readxl)
library(LaplacesDemon, lib.loc = '../packages/')    # library(LaplacesDemon)

#library(tidyverse) ; library(igraph) ; library(bisonR)

# information
sessionInfo()
R.Version()
#rstan::stan_version()

# set seed
set.seed(12345)

# set cmdstanr path
#set_cmdstan_path('H:/rlibs/4.2.1/')

#### short time windows ####
max_time_window <- 5
for( time_window in 1:max_time_window){
  ## set up ####
  ### add time marker
  print(paste0('start window ', time_window ,' at ', Sys.time()))
  
  # create output pdf
  pdf(file = paste0('../outputs/mpnp',time_window,'_bisonr_plotnetwork.pdf'))
  
  ### load data
  load(paste0('mpnp_edgecalculations/mpnpshort',time_window,'_bisonr_edgescalculated.RData'))
  rm(counts, draws, dyads, ew_chain, ew_edgelist, global_cv, mpnp_edges_null, priors, summary, draw98, filename, median98, sri98, model_averaging) ; gc()
  
  ## plot network ####
  # create nodes data frame
  nodes <- data.frame(node = sort(unique(c(counts_df_model$node_1_id,
                                           counts_df_model$node_2_id))),
                      id = NA,
                      age = NA,
                      sightings = NA)
    for(i in 1:nrow(nodes)){
    df <- counts_df[(counts_df$period == time_window & counts_df$node_1 == nodes$node[i]),
                    c('id_1','count_period_1')] %>% distinct()
    if(nrow(df) > 0){
      nodes$sightings[i] <- df$count_period_1[1]
      nodes$id[i] <- df$id_1[1]
    }
    else {
      df <- counts_df[(counts_df$period == time_window & counts_df$node_2 == nodes$node[i]),
                      c('id_2','count_period_2')] %>% distinct()
      nodes$sightings[i] <- df$count_period_2[1]
      nodes$id[i] <- df$id_2[1]
    }
    nodes$age[i] <- mean(mpnp_ages$age[which(mpnp_ages$id == nodes$id[i])])
    rm(df)
  }
  nodes$sightings <- as.numeric(nodes$sightings)
  str(nodes)
  
  # plot network -- 0.05
  plot_network_threshold(mpnp_edge_weights, lwd = 2, ci = 0.9, threshold = 0.05,
                         vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                         vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                         vertex.size2 = nodes$sightings, edge.color2 = 'black') 
  plot_network_threshold2(obj = mpnp_edge_weights, threshold = 0.05,
                          node.size = nodes, node.colour = nodes, lwd = 10)
  
  # plot network -- 0.1
  plot_network_threshold(mpnp_edge_weights, lwd = 2, ci = 0.9, threshold = 0.1,
                         vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                         vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                         vertex.size2 = nodes$sightings, edge.color2 = 'black') 
  plot_network_threshold2(obj = mpnp_edge_weights, threshold = 0.1,
                          node.size = nodes, node.colour = nodes, lwd = 10)
  
  # plot network -- 0.15
  if(length(which(edgelist$median > 0.15)) > 0){
    plot_network_threshold(mpnp_edge_weights, lwd = 2, ci = 0.9, threshold = 0.15,
                           vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                           vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                           vertex.size2 = nodes$sightings, edge.color2 = 'black') 
    plot_network_threshold2(obj = mpnp_edge_weights, threshold = 0.15,
                            node.size = nodes, node.colour = nodes, lwd = 10)
  }
  
  # plot network -- 0.2
  if(length(which(edgelist$median > 0.2)) > 0){
    plot_network_threshold(mpnp_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                           vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                           vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                           vertex.size2 = nodes$sightings, edge.color2 = 'black') 
    plot_network_threshold2(obj = mpnp_edge_weights, threshold = 0.2,
                            node.size = nodes, node.colour = nodes, lwd = 10)
  }
  
  ### add time marker
  print(paste0('network plots completed at ', Sys.time()))
  
## clean up ####
### save pdf
dev.off()

### clear environment
rm(list = ls()[!(ls() %in% 'time_window')])

### add time marker
print(paste0('time window ', time_window, ' completed at ', Sys.time()))

}

#### long window ####
# create nodes data frame
nodes <- data.frame(bull = sort(unique(ages$id)),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  x <- ages[ages$id == nodes$bull[i],]
  nodes$age[i] <- mean(x$age)
  if(nodes$bull[i] != 'M99'){
    y <- counts_df[counts_df$id_1 == nodes$bull[i], c('id_1','count_period_1')]
    nodes$sightings[i] <- y[1,2]
  }
}
nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

# plot network
plot_network_threshold(mpnp_edge_weights, lwd = 2, ci = 0.9, threshold = 0.15,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')
plot_network_threshold2(obj = mpnp_edge_weights, threshold = 0.15,
                        node.size = nodes, node.colour = nodes, lwd = 10)

### add time marker
print(paste0('network plots completed at ', Sys.time()))
