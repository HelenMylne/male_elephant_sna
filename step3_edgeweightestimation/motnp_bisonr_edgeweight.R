#### set up ####
# load packages
library(tidyverse)
library(dplyr)
library(rstan)
library(cmdstanr)
library(bisonR)
library(asnipe)
library(sna)

# information
sessionInfo()
R.Version()
rstan::stan_version()

# set seed
set.seed(12345)

#### create data lists ####
### import data for aggregated model (binomial)
#counts_df <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_bayesian_binomialpairwiseevents.csv')
counts_df <- read_csv('../data_processed/motnp_bayesian_binomialpairwiseevents.csv')

# correct sex_1, which has loaded in as a logical vector not a character/factor
unique(counts_df$sex_1) # FALSE or NA
sex_1 <- data.frame(sex_1 = counts_df$id_1)
sex_1 <- sex_1 %>% separate(sex_1, c("sex", "number"), sep = 1, remove = FALSE) ; unique(sex_1$sex) # F, M, U
counts_df$sex_1 <- as.character(sex_1$sex) ; rm(sex_1)
str(counts_df)  # sex_1 still comes up as logical, but now contains the right levels

# create variable for age difference
unique(counts_df$age_category_1) # "0-3","1-2","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-15","15-19","20-25","20-35","25-40","35-50","40+","50+",NA 
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
#rm(list = ls()[which(ls() != 'counts_df')])

### load in ages 
motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

### filter counts data frame down to males only
counts_df <- counts_df %>% 
  filter(id_1 %in% unique(motnp_ages$id)) %>% 
  filter(id_2 %in% unique(motnp_ages$id))

### remove young males
counts_df <- counts_df %>% 
  filter(age_cat_id_1 > 3) %>% 
  filter(age_cat_id_2 > 3)

### remove dead males
counts_df <- counts_df %>% 
  filter(name_1 != 'Richard' & name_2 != 'Richard') %>% 
  filter(name_1 != 'Gabriel' & name_2 != 'Gabriel')

### create nodes id variable which is 1:max(id) as currently covers all elephants, not just males
counts_df$node_1_males <- as.integer(as.factor(counts_df$node_1_nogaps))
counts_df$node_2_males <- as.integer(as.factor(counts_df$node_2_nogaps))+1

### standardise dyad_id for males only
counts_df$dyad_males <- as.integer(as.factor(counts_df$dyad_id))

#### edge weights ####
### create data frame for edge weight model
counts_df_model <- counts_df[, c('node_1_males','node_2_males','event_count','count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### set priors
priors <- get_default_priors('binary_conjugate')
priors$edge <- 'beta(2,2)'
#priors$edge <- 'normal(-1, 2)'    # slightly right skewed so that more density is towards 0 but still very weak and all values are possible
#priors$fixed <- 'normal(0, 1.5)'  # slightly less wide that before, but zero centred so allows age to have a positive or negative impact on centrality and edge weight
prior_check(priors, 'binary_conjugate')

### run edge weight model
motnp_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary_conjugate",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(motnp_edge_weights, par_ids = 2)
plot_predictions(motnp_edge_weights, num_draws = 20, type = "density")
plot_predictions(motnp_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(motnp_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(motnp_edge_weights)

### save workspace image for reloading at a later date that doesn't require running model again
save.image('motnp_bisonr_edgescalculated.RData')

#### non-random edge weights ####
### load in workspace image with edge weights already calculated
load('motnp_bisonr_edgescalculated.RData')

### run null model
motnp_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary_conjugate",
  priors = priors
)

### compare null model with fitted model
model_comparison(list(non_random_model = motnp_edge_weights, random_model = motnp_edges_null)) # NETWORK IS RANDOM
# Method: stacking
#                  weight
# non_random_model 0.020
# random_model     0.980
# Warning message:  Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.

# save workspace image
save.image('motnp_bisonr_randomnetwork.RData')

#### plot network ####
### adapt bisonr plot_netowrk function to give more flexibility over plotting options
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

### create nodes data frame
nodes <- data.frame(bull = sort(unique(motnp_ages$id)),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  x <- motnp_ages[motnp_ages$id == nodes$bull[i],]
  nodes$age[i] <- mean(x$age)
  if(nodes$bull[i] != 'M99'){
    y <- counts_df[counts_df$id_1 == nodes$bull[i], c('id_1','count_1')]
    nodes$sightings[i] <- y[1,2]
  }
}
nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

### plot network
plot_network_threshold(motnp_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')

### adapt to remove unconnected nodes
plot_network_threshold2 <- function (obj, ci = 0.9, lwd = 2, threshold = 0.3,
                                     vertex.label.color1 = 'transparent', vertex.label.font1 = NA, 
                                     vertex.color1 = 'transparent',
                                     vertex.size1 = 1, edge.color1 = rgb(0, 0, 0, 1),
                                     vertex.label.color2 = 'black', vertex.label.font2 = 'Helvetica', 
                                     vertex.color2 = 'seagreen1',
                                     vertex.size2 = 8, edge.color2 = rgb(0, 0, 0, 0.3))
{
  edgelist <- get_edgelist(obj, ci = ci, transform = TRUE)
  threshold_edges <- edgelist[edgelist$median >= threshold,]
  net <- igraph::graph_from_edgelist(as.matrix(threshold_edges[, 1:2]))
  md <- threshold_edges[, 3]
  ub <- threshold_edges[, 5]
  coords <- igraph::layout_nicely(net)
  igraph::plot.igraph(net, layout = coords,
                      vertex.label.color = vertex.label.color1, label.family = vertex.label.font1, 
                      vertex.color = vertex.color1, vertex.size = vertex.size1,
                      edge.color = rgb(0, 0, 0, 1), edge.arrow.size = 0, edge.width = md * lwd)
  igraph::plot.igraph(net, layout = coords,
                      vertex.label.color = vertex.label.color2, label.family = vertex.label.font2,
                      vertex.color = vertex.color2, vertex.size = vertex.size2,
                      edge.color = edge.color1, edge.arrow.size = 0, edge.width = ub * lwd,
                      add = TRUE)
  if (obj$directed) {
    igraph::plot.igraph(net, edge.width = md * 0, layout = coords, 
                        vertex.label.color = vertex.label.color2, vertex.color = vertex.color2, 
                        edge.color = edge.color2, add = TRUE)
  }
}

### plot network
plot_network_threshold2(obj = motnp_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)


#### check model output ####
### load in workspace image with edge weights already calculated
load('motnp_bisonr_edgescalculated.RData')

### compare edge weight distributions to simple SRI ####
head(edgelist)
colnames(edgelist)[1:2] <- colnames(counts_df_model)[1:2]
edgelist$node_1_id <- as.integer(edgelist$node_1_id) ; edgelist$node_2_id <- as.integer(edgelist$node_2_id)
summary <- left_join(edgelist, counts_df_model, by = c('node_1_id','node_2_id'))
summary$sri <- summary$event / (summary$duration)

plot(density(summary$sri))
lines(density(summary$median), col = 'blue')
lines(density(summary$`5%`), col = 'red')
lines(density(summary$`95%`), col = 'green')

plot(density(summary$sri[which(summary$duration >= 8)]), xlim = c(0,1))
lines(density(summary$median[which(summary$duration >= 8)]), col = 'blue')
lines(density(summary$`5%`[which(summary$duration >= 8)]), col = 'red')
lines(density(summary$`95%`[which(summary$duration >= 8)]), col = 'green')

plot(density(summary$sri[which(summary$duration >= 12)]), xlim = c(0,1))
lines(density(summary$median[which(summary$duration >= 12)]), col = 'blue')
lines(density(summary$`5%`[which(summary$duration >= 12)]), col = 'red')
lines(density(summary$`95%`[which(summary$duration >= 12)]), col = 'green')

certain_sri <- summary[summary$sri == 1,]
plot(density(certain_sri$sri), xlim = c(0,1))
lines(density(certain_sri$median), col = 'blue') # just above 0.5
lines(density(certain_sri$`5%`), col = 'red')    # very low 
lines(density(certain_sri$`95%`), col = 'green') # very high -- overall just very high uncertainty

# try plotting a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
dyads <- counts_df[,c('dyad_id','node_1_males','node_2_males')]
colnames(dyads) <- c('dyad_id','node_1_id','node_2_id')
dyads <- left_join(dyads, counts_df_model, by = c('node_1_id','node_2_id'))
length(which(is.na(dyads$duration) == TRUE))

draws <- as.data.frame(motnp_edge_weights$chain) %>% pivot_longer(cols = everything(), names_to = 'dyad_model', values_to = 'edge')
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

ggplot(data = subset_draws, mapping = aes(x = weight))+
  geom_density(colour = 'blue')+
  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

write_csv(subset_draws, '../data_processed/motnp_sampledyads_random_model_vs_sri.csv')

subset_draws <- draws[draws$sri > 0.2,]
subset_draws$median <- NA
for(i in 1:nrow(subset_draws)){
  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]
  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), subset_draws$median[i])
}
head(subset_draws)
which(is.na(subset_draws$median) == TRUE)[1]

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

write_csv(subset_draws, '../data_processed/motnp_sampledyads_sri0.2_model_vs_sri.csv')

# clean environment
rm(draws, dyads, priors, subset_draws, x)

### coefficient of variation of edge weights (aka social differentation) ####
global_cv <- extract_metric(motnp_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

# create gbi_matrix
eles <- read_delim(file = '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_eles_long.csv', delim = ',')
eles$location <- paste(eles$gps_s, eles$gps_e, sep = '_') # make single variable for unique locations
eles <- eles[,c(1,16,2,3,17,4,5,14,7,8,10,13)]            # rearrange variables
nodes <- read_delim(file = '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_elenodes.csv', delim = ',') # read in node data
colnames(nodes)
str(nodes)
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

# obtain 100,000 random permutations of network and compare distribution to global_cv
N <- 100000                                                                     # number of permutations
random_networks <- asnipe::network_permutation(association_data = gbi_matrix,          # permute network
                                               #association_matrix = gbi_matrix,
                                               permutations = N)
cv_random_networks <- rep(0,N)                                                  # generate empty vector to fill with cv values
for (i in c(1:N)) {
  net_rand <- random_networks[i,,]
  cv_random_networks[i] <- extract_metric(net_rand,'global_cv')
}

cv_random_nodes <- rep(0,N)
for (i in c(1:N)) {
  net_rand <- sna::rmperm(network)
  cv_random_nodes[i] <- extract_metric(net_rand,'global_cv')
}

ggplot()+
  geom_histogram(mapping = aes(xintercept = global_cv))+
  geom_density(mapping = aes(x = cv_random_networks), colour = 'blue')+
  geom_density(mapping = aes(x = cv_random_nodes), colour = 'red')
  
### standard deviation edge weight
edge_weight <- extract_metric(motnp_edge_weights, "edge_weight")
head(edge_weight)
summary$std <- NA
for(i in 1:nrow(summary)){
  summary$std[i] <- sd(edge_weight[,i])
}
hist(summary$std)

hist(summary$std/summary$median)



















