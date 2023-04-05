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
females_df <- counts_df
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
priors$edge <- 'normal(-2.5, 1.5)'     # slightly right skewed so that more density is towards 0 but still very weak and all values are possible
priors$fixed <- 'normal(0, 1.5)'       # slightly less wide than before, but zero centred so allows age to have a positive or negative impact on centrality and edge weight
prior_check(priors, 'binary')

### run edge weight model
motnp_edge_weights_strongpriors <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id),   # count of sightings together given count of total sightings as a result of the individuals contained within the dyad
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
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

priors$edge <- 'normal(-2.5, 1.5)'                                              # set prior value (already set -- doesn't actually change anything, just in as a reminder of what it was)
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
nodes <- data.frame(id = sort(unique(c(counts_df$id_1,counts_df$id_2))),
                    node = NA,
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  x <- motnp_ages[motnp_ages$id == nodes$id[i],]
  nodes$age[i] <- mean(x$age)
  if(nodes$id[i] != 'M99') { y <- counts_df[counts_df$id_1 == nodes$id[i], c('id_1','count_1','node_1_males')] }
  else { y <- counts_df[counts_df$id_2 == nodes$id[i], c('id_2','count_2','node_2_males')] }
  nodes$sightings[i] <- y[1,2]
  nodes$node[i] <- y[1,3]
}
nodes$sightings <- as.numeric(nodes$sightings)
nodes$node <- as.character(nodes$node)
str(nodes)

### plot network
plot_network_threshold(motnp_edge_weights_strongpriors, lwd = 15, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')

### adapt to remove unconnected nodes
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
head(edgelist)
colnames(edgelist)[1:2] <- colnames(counts_df_model)[1:2]
edgelist$node_1_id <- as.integer(edgelist$node_1_id) ; edgelist$node_2_id <- as.integer(edgelist$node_2_id)
summary <- left_join(edgelist, counts_df_model, by = c('node_1_id','node_2_id'))
counts <- counts_df[,c('node_1_males','node_2_males','count_1','count_2')]
colnames(counts)[1:2] <- c('node_1_id','node_2_id')
summary <- left_join(summary, counts, by = c('node_1_id','node_2_id'))
summary$sri <- summary$event / (summary$duration)

plot(density(summary$sri), main = 'SRI vs model output: blue=all,\nred=both seen 8 times, green=both 12 times')
lines(density(summary$median), col = 'blue')
lines(density(summary$median[which(summary$count_1 >= 8 & summary$count_2 >= 8)]), col = 'red')
lines(density(summary$median[which(summary$count_1 >= 12 & summary$count_2 >= 12)]), col = 'green')

# try plotting a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
dyads <- counts_df[,c('dyad_id','node_1_males','node_2_males')]
colnames(dyads) <- c('dyad_id','node_1_id','node_2_id')
dyads <- left_join(dyads, counts_df_model, by = c('node_1_id','node_2_id'))
length(which(is.na(dyads$duration) == TRUE))

draws <- as.data.frame(motnp_edge_weights_strongpriors$chain) %>% pivot_longer(cols = everything(), names_to = 'dyad_model', values_to = 'edge')
draws$dyad_id <- rep(counts_df$dyad_id, 4000) ## IS THIS RIGHT??
draws$weight <- gtools::inv.logit(draws$edge)
draws$draw <- rep(1:4000,  each = nrow(counts_df_model))
draws <- left_join(draws, dyads, by = 'dyad_id')
draws$sri <- draws$event / draws$duration

set.seed(15)
subset_draws <- draws[draws$dyad_id %in% sample(draws$dyad_id, 150, replace = F),]
subset_draws$median <- NA
for(i in 1:nrow(subset_draws)){
  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]
  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), subset_draws$median[i])
}
head(subset_draws)
which(is.na(subset_draws$median) == TRUE)[1]

subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)
ggplot(data = subset_draws, mapping = aes(x = weight))+
  geom_density(colour = 'blue')+
  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

write_csv(subset_draws, '../data_processed/motnp_sampledyads_random_binary_vs_sri_strongprior.csv')

subset_draws <- draws[draws$sri > 0.2,]
subset_draws$median <- NA
for(i in 1:nrow(subset_draws)){
  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]
  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), subset_draws$median[i])
}
head(subset_draws)
which(is.na(subset_draws$median) == TRUE)[1]

subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)
ggplot(data = subset_draws, mapping = aes(x = weight))+
  geom_density(colour = 'blue')+
  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15)+
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

ggplot(data = subset_draws, mapping = aes(x = weight))+
  geom_density(colour = 'blue')+
  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

write_csv(subset_draws, '../data_processed/motnp_sampledyads_sri0.2_binary_vs_sri_strongprior.csv')

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
N <- length(unique(c(counts_df$id_1, counts_df$id_2)))
ids <- unique(c(counts_df$id_1, counts_df$id_2))
m_mat <- diag(nrow = N)
colnames(m_mat) <- ids
rownames(m_mat) <- ids

# populate matrix with SRI values
for( i in 1:N ) {
  for( j in 1:N ) {
    if( i >= j ) {
      m_mat[i,j] <- m_mat[j,i]
    }
    else {
      id1 <- colnames(m_mat)[i]
      id2 <- rownames(m_mat)[j]
      m_mat[i,j] <- counts_df$sri[which(counts_df$id_1 == id1 & counts_df$id_2 == id2)]
    }
  }
}

# calculate CV for matrix
raster::cv(m_mat) # still way too high
which(is.nan(m_mat) == TRUE)
sd(m_mat)/mean(m_mat) # matches raster version -- check it is doing what I think it is doing

# create gbi_matrix
eles <- read_csv(file = '../data_processed/motnp_eles_long.csv')
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
sightings <- eles[,c('elephant','encounter','date')] %>% 
  filter(elephant %in% ids) %>% 
  select(-elephant) %>% 
  distinct()

# run network permutations
random_networks <- asnipe::network_permutation(association_data = gbi_males,
                                               association_matrix = m_mat,
                                               permutations = N_networks,
                                               days = sightings$date,
                                               within_day = TRUE)
print('random networks generated')

cv_random_networks <- rep(0,N_networks)
for (i in c(1:N_networks)) {
  net_rand <- random_networks[i,,]
  cv_random_networks[i] <- raster::cv(net_rand)
  if(i %% 1000 == 0) {print(i)}
}
print(paste0('network permutations for entire network completed at ', Sys.time()))

# compare permuted networks to actual network
hist(cv_random_networks, xlim = c(min(cv_random_networks)-50,max(cv_random_networks)+50),
     #ylim = c(0,10),  # there are 3 networks out of 10000 that have cv > 290, but then all are in the same bar as cv(m_mat)
     main = 'permutations for network of all male elephants')
abline(v = raster::cv(m_mat), col = 'red')
text(round(raster::cv(m_mat),3), col = 'red', x = median(cv_random_networks)+50, y = 2000)

plot_cv <- as.data.frame(cv_random_networks)
cv_network <- raster::cv(m_mat)
cv_random_networks <- plot_cv$cv_random_networks
ggplot(data = plot_cv)+
  geom_vline(xintercept = cv_network, linewidth = 1.5,
             colour = rgb(68/255, 1/255, 84/255))+
  #geom_text(x = cv_network - 50, y = 700,
  geom_histogram(aes(cv_random_networks),
                 fill = rgb(253/255, 231/255, 37/255),
                 colour = 'black',
                 bins = 100)+
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
plot_cv_model <- data.frame(cv = c(plot_cv$cv_random_networks, (global_cv_strongprior*100)),
                      iteration = rep(1:length(global_cv_strongprior), 2),
                      type = rep(c('permutation','model_draw'), each = length(global_cv_strongprior)))
write_csv(plot_cv_model, '../data_processed/motnp_networkpermutations_cv_strongprior.csv')

hist(plot_cv_model$cv)

ggplot()+
  geom_vline(xintercept = cv_network, linewidth = 1.5,
             colour = rgb(68/255, 1/255, 84/255))+
  #annotate('text', x = cv_network - 50, y = 700), colour = rgb(68/255, 1/255, 84/255),
  #         label = 'coefficient of/nvariation for/nmeasured network')+
  theme_classic()+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))+
  geom_histogram(data = plot_cv_model[plot_cv_model$type == 'model_draw',], aes(x = cv),
                 fill = rgb(33/255, 145/255, 140/255),
                 bins = 100,
                 colour = 'black')+
  geom_histogram(data = plot_cv_model[plot_cv_model$type == 'permutation',], aes(x = cv),
                 fill = rgb(253/255, 231/255, 37/255),
                 bins = 100,
                 colour = 'black')+
  scale_x_continuous(name = 'coefficient of variation',
                     #limits = c(min(cv_random_networks)-10,
                     #cv_network+10),
                     expand = c(0,0))+
  scale_y_continuous(name = 'frequency',
                     expand = c(0,0))

plot_cv_model2 <- plot_cv_model[c(1:10000, sample(10001:20000, 4000, replace = F)),]
ggplot()+
  geom_vline(xintercept = cv_network, linewidth = 1.5,
             colour = rgb(68/255, 1/255, 84/255))+
  #annotate('text', x = cv_network - 50, y = 700), colour = rgb(68/255, 1/255, 84/255),
  #         label = 'coefficient of/nvariation for/nmeasured network')+
  theme_classic()+
  theme(axis.title = element_text(size = 28),
        axis.text = element_text(size = 22))+
  geom_histogram(data = plot_cv_model2[plot_cv_model2$type == 'model_draw',], aes(x = cv),
                 fill = rgb(33/255, 145/255, 140/255),
                 bins = 100,
                 colour = 'black')+
  geom_histogram(data = plot_cv_model2[plot_cv_model2$type == 'permutation',], aes(x = cv),
                 fill = rgb(253/255, 231/255, 37/255),
                 bins = 100,
                 colour = 'black')+
  scale_x_continuous(name = 'coefficient of variation',
                     #limits = c(min(cv_random_networks)-10,
                     #cv_network+10),
                     breaks = round(seq(100, 400, by = 50),-1),
                     expand = c(0,0))+
  scale_y_continuous(name = 'frequency',
                     expand = c(0,0))

### run statistical test to confirm that difference is significant
# plot_cv <- read_csv('../data_processed/motnp_networkpermutations_cv_strongprior.csv')
hist(plot_cv$cv_random_networks)
cv_network

mean(plot_cv$cv_random_networks)
sd(plot_cv$cv_random_networks)

length(which(plot_cv$cv_random_networks >= cv_network))/length(plot_cv$cv_random_networks) # 3 / 10000 = 3e-04 = p-value (don't need an additional test because that will further randomise the data which are already randomised, just need to know proportion that are greater than or equal to your value)

mean(global_cv)
sd(global_cv)

max(global_cv) ; min(plot_cv$cv_random_networks)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2% ####
counts_df_model <- counts_df[,c('node_1_males','node_2_males','event_count','count_dyad','id_1','id_2','dyad_males')]
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ms_chain <- as.data.frame(motnp_edge_weights_strongpriors$chain)
colnames(ms_chain) <- counts_df_model$dyad_id
ms_chain <- pivot_longer(ms_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ms_chain$chain_position <- rep(1:length(unique(ms_chain$dyad_id)), each = 4000)
ms_chain$draw <- LaplacesDemon::invlogit(ms_chain$draw)
ms_chain$mean <- NA
hist(ms_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sample(unique(ms_chain$dyad_id), size = 1000, replace = F)){
  x <- ms_chain[ms_chain$dyad_id == i,]
  ms_chain$mean <- ifelse(ms_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ms_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ms_chain$draw), lwd = 2)

quantile(ms_chain$draw, 0.98)

ms_edgelist <- bisonR::get_edgelist(motnp_edge_weights_strongpriors)
quantile(ms_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

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

