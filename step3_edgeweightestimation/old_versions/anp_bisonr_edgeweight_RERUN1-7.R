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

# information
sessionInfo()
R.Version()
#rstan::stan_version()

# set seed
set.seed(12345)

# set cmdstanr path
#set_cmdstan_path('H:/rlibs/4.2.1/')

pdf(file = '../outputs/anp_edgeweight_setup.pdf')

#### create priors and functions for all ####
# set priors
priors <- get_default_priors('binary')
priors$edge <- 'normal(-2.5, 1.5)'    # slightly right skewed so that more density is towards 0 but still very weak and all values are possible
priors$fixed <- 'normal(0, 1.5)'  # slightly less wide that before, but zero centred so allows age to have a positive or negative impact on centrality and edge weight
prior_check(priors, 'binary')

### adapt bisonR plot_network function to give more flexibility over plotting options
plot_network_threshold <- function (obj, ci = 0.9, lwd = 2, threshold = 0.3,
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
                      vertex.label.color = vertex.label.color1, vertex.color = vertex.color1, vertex.size = vertex.size1,
                      edge.color = rgb(0, 0, 0, 1), edge.arrow.size = 0)
  igraph::plot.igraph(net, layout = coords,
                      edge.width = ifelse(md < threshold, 0, ub * lwd),
                      vertex.label.color = vertex.label.color2, vertex.color = vertex.color2, vertex.size = vertex.size2,
                      edge.color = edge.color1, edge.arrow.size = 0, 
                      add = TRUE)
  if (obj$directed) {
    igraph::plot.igraph(net, edge.width = md * 0, layout = coords, 
                        vertex.label.color = vertex.label.color2, vertex.color = vertex.color2, 
                        edge.color = edge.color2, add = TRUE)
  }
}

### adapt to remove unconnected nodes
plot_network_threshold2 <- function (obj, ci = 0.9, lwd = 2, threshold = 0.3,
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

# alternative method to compare models
model_averaging <- function(models) {
  loos <- lapply(models, function(model) loo::loo(model$log_lik, r_eff=NA))
  
  results_matrix <- loo::loo_model_weights(loos, method="pseudobma")
  if (!is.null(names(models))) {
    names(results_matrix) <- names(models)
  }
  results_matrix
}

dev.off()

#### short time windows ####
### create data lists ####
# read in dyad data
counts_df_allwindows <- read_csv('../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_shortwindows_impossiblepairsremoved.csv')

# read in sightings data
ate <- read_csv('../data_processed/anp_sightings_rawcombined.csv')
head(ate)
table(ate$num_bulls)
ate$num_bulls_recount <- NA
for(i in 1:nrow(ate)){
  x <- ate[ate$obs_id == ate$obs_id[i],]
  ate$num_bulls_recount[i] <- nrow(x)
}
mean(ate$num_bulls_recount) ; sd(ate$num_bulls_recount)
mean(ate$num_bulls) ; sd(ate$num_bulls)
rm(x, i) ; gc()

# count number of males per sighting
counts <- as.data.frame(table(ate$id)) ; colnames(counts) <- c('id','count')
summary(counts$count) ; sd(counts$count)

# narrow down to sightings data of interest
sightings <- ate[,c('obs_id','obs_date','correct_time_hms','grid_code','obs_type_old')] %>% distinct()
length(unique(sightings$obs_id))         # 24386 -- duplicates due to wrong obs_type
sightings$distinct <- NA
for(i in 1:nrow(sightings)) {
  x <- sightings[sightings$obs_id == sightings$obs_id[i], ]
  sightings$distinct[i] <- nrow(x)
  if(is.na(sightings$obs_type_old[i]) == TRUE ) {
    obs_type <- unique(x$obs_type_old[!is.na(x$obs_type_old)])
    sightings$obs_type_old[i] <- ifelse(length(obs_type) == 1, obs_type, 'U')
  }
}
sightings <- distinct(sightings) # now the right length
table(sightings$obs_type_old)

# read in individual male data
males <- readxl::read_excel('../data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% 
  clean_names()
males$id <- paste0('M',males$casename)
ids <- sort(unique(ate$id))
males <- males %>% dplyr::filter(id %in% ids)

# counts per male per time period
length(unique(counts_df_allwindows$period))
males$count_all <- NA
males$p1 <- NA ; males$p2 <- NA ; males$p3 <- NA ; males$p4 <- NA ; males$p5 <- NA ; males$p6 <- NA ; males$p7 <- NA
males$p8 <- NA ; males$p9 <- NA ; males$p10 <- NA ; males$p11 <- NA ; males$p12 <- NA ; males$p13 <- NA
males$p14 <- NA ; males$p15 <- NA ; males$p16 <- NA ; males$p17 <- NA ; males$p18 <- NA ; males$p19 <- NA
males$p20 <- NA ; males$p21 <- NA ; males$p22 <- NA ; males$p23 <- NA ; males$p24 <- NA ; males$p25 <- NA
males$p26 <- NA ; males$p27 <- NA ; males$p28 <- NA ; males$p29 <- NA ; males$p30 <- NA ; males$p31 <- NA
males$p32 <- NA ; males$p33 <- NA ; males$p34 <- NA ; males$p35 <- NA ; males$p36 <- NA
windows <- c(sort(unique(counts_df_allwindows$period_start)),max(counts_df_allwindows$period_end))

for(i in 1:nrow(males)){
  counts <- ate[ate$id == males$id[i],]
  males$count_all[i] <- nrow(counts)
  for(j in 1:max(counts_df_allwindows$period)){
    males[i,j+22] <- length(which(counts$obs_date >= windows[j] & counts$obs_date < windows[j+1]))
  }
}

# time windows
periods <- data.frame(period = 1:(length(windows)-1),
                      period_start = windows[1:(length(windows)-1)],
                      period_end = windows[2:(length(windows))])

# check out total sightings per individual/dyad per period and see if they are reasonable
table(counts_df_allwindows$count_period_1) ; hist(counts_df_allwindows$count_period_1)
table(counts_df_allwindows$count_period_2) ; hist(counts_df_allwindows$count_period_2)
table(counts_df_allwindows$count_period_dyad) ; hist(counts_df_allwindows$count_period_dyad, breaks = 30)
table(counts_df_allwindows$event_count) ; hist(counts_df_allwindows$event_count, breaks = 20, ylim = c(0,15000))
# many are not seen much, but repeat sightings of pairs do seem to be genuinely rare

rm(counts, sightings, i, x, j, obs_type)

### add time marker
print(paste0('data read in at ', Sys.time()))

### run models ####
for( time_window in 1:7){
  ## set up ####
  ### add time marker
  print(paste0('start window ', time_window ,' at ', Sys.time()))
  
  # create output pdf
  pdf(file = paste0('../outputs/anp',time_window,'_bisonr_edgeweight.pdf'))
  
  ### create data frame for edge weight model
  counts_df <- counts_df_allwindows[counts_df_allwindows$period == time_window,]
  counts_df_model <- counts_df[, c('node_1','node_2','event_count','count_period_dyad')] %>% distinct()
  colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')
  
  ## run model ####
  ### run edge weight model
  anp_edge_weights <- bison_model(
    ( event | duration ) ~ dyad(node_1_id, node_2_id), 
    data = counts_df_model, 
    model_type = "binary",
    #partial_pooling = TRUE,
    #mc_cores = 4,
    priors = priors
  )
  
  ### run diagnostic plots
  plot_trace(anp_edge_weights, par_ids = 2)
  plot_predictions(anp_edge_weights, num_draws = 20, type = "density")
  plot_predictions(anp_edge_weights, num_draws = 20, type = "point")
  
  ### extract edge weight summaries
  edgelist <- get_edgelist(anp_edge_weights, ci = 0.9, transform = TRUE)
  plot(density(edgelist$median))
  summary(anp_edge_weights)
  
  ### add time marker
  print(paste0('model completed at ', Sys.time()))
  
  ## non-random edge weights ####
  ### run null model
  anp_edges_null <- bison_model(
    (event | duration) ~ 1, 
    data = counts_df_model, 
    model_type = "binary",
    priors = priors
  )
  
  ### compare null model with fitted model -- bisonR model stacking
  #model_comparison(list(non_random_model = anp_edge_weights, random_model = anp_edges_null)) # NETWORK IS RANDOM
  
  ### compare null model with fitted model -- Jordan pseudo model averaging
  #model_averaging(models = list(non_random_model = anp_edge_weights, random_model = anp_edges_null))
  
  ### add time marker
  print(paste0('random network comparison completed at ', Sys.time()))
  
  ## plot network ####
  # create nodes data frame
  nodes <- data.frame(bull = sort(unique(c(counts_df_model$node_1_id,counts_df_model$node_2_id))),
                      age = NA,
                      sightings = NA)
  for(i in 1:nrow(nodes)){
    df <- counts_df[(counts_df$period == time_window & counts_df$node_1 == nodes$bull[i]),
                    c('id_1','age_start_1','count_period_1')] %>% distinct()
    if(nrow(df) > 0){
      nodes$age[i] <- df$age_start_1[1]
      nodes$sightings[i] <- df$count_period_1[1]
    }
    else {
      df <- counts_df[(counts_df$period == time_window & counts_df$node_2 == nodes$bull[i]),
                      c('id_2','age_start_2','count_period_2')] %>% distinct()
      nodes$age[i] <- df$age_start_2[1]
      nodes$sightings[i] <- df$count_period_2[1]
    }
  }
  nodes$sightings <- as.numeric(nodes$sightings)
  str(nodes)
  
  # plot network
  #plot_network_threshold(anp_edge_weights, lwd = 2, ci = 0.9, threshold = 0.15,
  #                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
  #                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
  #                       vertex.size2 = nodes$sightings, edge.color2 = 'black') 
  #plot_network_threshold2(obj = anp_edge_weights, threshold = 0.15,
  #                        node.size = nodes, node.colour = nodes, lwd = 10)
  
  ### add time marker
  print(paste0('network plots completed at ', Sys.time()))
  
  ## coefficient of variation of edge weights (aka social differentiation) ####
  # extract cv for model
  global_cv <- extract_metric(anp_edge_weights, "global_cv")
  head(global_cv)
  hist(global_cv)
  
  ### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
  counts_df_model <- counts_df[,c('node_1','node_2','event_count','count_period_dyad','id_1','id_2','dyad_id')]
  colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
  ew_chain <- as.data.frame(anp_edge_weights$chain)
  colnames(ew_chain) <- counts_df_model$dyad_id
  ew_chain <- pivot_longer(ew_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
  ew_chain$chain_position <- rep(1:length(unique(ew_chain$dyad_id)), each = 4000)
  ew_chain$draw <- LaplacesDemon::invlogit(ew_chain$draw)
  ew_chain$mean <- NA
  hist(ew_chain$draw)
  #plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'edge distribution', xlab = 'edge weight', ylab = 'density', las = 1)
  #for(i in sort(unique(ew_chain$dyad_id))){
  #  x <- ew_chain[ew_chain$dyad_id == i,]
  #  ew_chain$mean <- ifelse(ew_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ew_chain$mean)
  #  lines(density(x$draw), col = rgb(0,0,1,0.1))
  #}
  #lines(density(ew_chain$draw), lwd = 2)
  
  (draw98 <- quantile(ew_chain$draw, 0.98))
  
  ew_edgelist <- bisonR::get_edgelist(anp_edge_weights)
  (median98 <- quantile(ew_edgelist$median, 0.98))
  
  counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
  (sri98 <- quantile(counts_df_model$sri, 0.98))
  
  ## clean up ####
  ### save pdf
  dev.off()
  
  # save workspace image
  save.image(file = paste0('anp_edgecalculations/anpshort',time_window,'_bisonr_edgescalculated.RData'))
  
  ### clear environment
  rm(list = ls()[!(ls() %in% c('counts_df_allwindows','time_window','periods','windows','males',
                               'ate','sightings','priors',
                               'model_averaging', 'plot_network_threshold','plot_network_threshold2'))])
  
  ### add time marker
  print(paste0('time window ', time_window, ' completed at ', Sys.time()))
  
}

