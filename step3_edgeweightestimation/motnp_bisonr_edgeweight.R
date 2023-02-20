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
### import data for aggregated model (binomial)
counts_df <- read_csv('../data_processed/motnp_bayesian_binomialpairwiseevents.csv')

# correct sex_1, which has loaded in as a logical vector not a character/factor
unique(counts_df$sex_1) # FALSE or NA
sex_1 <- data.frame(sex_1 = counts_df$id_1)
sex_1 <- sex_1 %>% separate(sex_1, c("sex", "number"), sep = 1, remove = FALSE) ; unique(sex_1$sex) # F, M, U
counts_df$sex_1 <- as.character(sex_1$sex) ; rm(sex_1) ; gc()
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
  filter(age_cat_id_1 > 3) %>% 
  filter(age_cat_id_2 > 3)

### create nodes id variable which is 1:max(id) as currently covers all elephants, not just males
counts_df$node_1_males <- as.integer(as.factor(counts_df$node_1_nogaps))
counts_df$node_2_males <- as.integer(as.factor(counts_df$node_2_nogaps))+1

### standardise dyad_id for males only
counts_df$dyad_males <- as.integer(as.factor(counts_df$dyad_id))

### add time marker
print(paste0('data read in at ', Sys.time()))

#### edge weights -- original: males only, weak priors ####
# add pdf output file
pdf(file = '../outputs/motnp_bisonr_edgeweight_weakprior.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[, c('node_1_males','node_2_males','event_count','count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### set priors
priors <- get_default_priors('binary')
priors$edge <- 'normal(-1, 2)'    # slightly right skewed so that more density is towards 0 but still very weak and all values are possible
priors$fixed <- 'normal(0, 1)'  # slightly less wide that before, but zero centred so allows age to have a positive or negative impact on centrality and edge weight
prior_check(priors, 'binary')

### run edge weight model
motnp_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### add time marker
print(paste0('basic model completed at ', Sys.time()))

### run diagnostic plots
plot_trace(motnp_edge_weights, par_ids = 2)
plot_predictions(motnp_edge_weights, num_draws = 20, type = "density")
plot_predictions(motnp_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(motnp_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(motnp_edge_weights)

## non-random edge weights ####
### run null model
motnp_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = motnp_edge_weights, random_model = motnp_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging <- function(models) {
  loos <- lapply(models, function(model) loo::loo(model$log_lik, r_eff=NA))
  
  results_matrix <- loo::loo_model_weights(loos, method="pseudobma")
  if (!is.null(names(models))) {
    names(results_matrix) <- names(models)
  }
  results_matrix
}
model_averaging(models = list(non_random_model = motnp_edge_weights, random_model = motnp_edges_null))

### save workspace image
save.image('motnp_bisonr_edgescalculated_weakprior.RData')

### add time marker
print(paste0('null model fitted and compared at ', Sys.time()))

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

### add time marker
print(paste0('networks plotted at ', Sys.time()))

## compare edge weight distributions to simple SRI ####
### load in workspace image with edge weights already calculated
#load('motnp_bisonr_edgescalculated_weakprior.RData')

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

subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)
#labels <- subset_draws[,c('dyad_id','node_1_id','node_2_id')] %>% distinct()
ggplot(data = subset_draws, mapping = aes(x = weight))+
  geom_density(colour = 'blue')+
  facet_wrap(. ~ dyad_id,
             #scales = 'free_y',
             nrow = 10, ncol = 15)+
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')#+
#geom_text(data = labels, aes(label = dyad_id))

write_csv(subset_draws, '../data_processed/motnp_sampledyads_random_binary_vs_sri_weakprior.csv')

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

write_csv(subset_draws, '../data_processed/motnp_sampledyads_sri0.2_binary_vs_sri_weakprior.csv')

# clean environment
rm(draws, dyads, priors, subset_draws, x) ; gc()
dev.off()

## coefficient of variation of edge weights (aka social differentiation) ####
# add pdf output file
pdf(file = '../outputs/motnp_bisonr_edgeweight_weakprior_cv.pdf')

# extract cv for model
global_cv_weakprior <- extract_metric(motnp_edge_weights, "global_cv")
head(global_cv_weakprior)
hist(global_cv_weakprior)

# calculate SRI for all dyads
counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
summary(counts_df_model$sri)

# identify id numbers from node id values
m_nodes <- counts_df[,c('node_1_males','node_2_males','id_1','id_2')]
colnames(m_nodes)[1:2] <- c('node_1_id','node_2_id')
m <- counts_df_model %>% left_join(m_nodes, by = c('node_1_id','node_2_id'))
rm(m_nodes) ; gc()

# calculate CV of SRI for all dyads
raster::cv(counts_df_model$sri)  # very high, but lower than gbi_matrix
#raster::cv(m$sri[m$sri > 0])     # still massive even when I remove the 0s -- zero inflation is real

### create SRI matrix
# generate matrix
N <- length(unique(c(m$id_1, m$id_2)))
ids <- unique(c(m$id_1, m$id_2))
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
      m_mat[i,j] <- m$sri[which(m$id_1 == id1 & m$id_2 == id2)]
    }
  }
}

# calculate CV for matrix
raster::cv(m_mat) # still way too high
which(is.nan(m_mat) == TRUE)
sd(m_mat)/mean(m_mat) # matches raster version -- check it is doing what I think it is doing

# create gbi_matrix
# create matrix of only males
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

# run network permutations
random_networks <- asnipe::network_permutation(association_data = gbi_males,
                                               association_matrix = m_mat,
                                               permutations = N_networks)
print('random networks generated')

cv_random_networks <- rep(0,N_networks)
for (i in c(1:N_networks)) {
  net_rand <- random_networks[i,,]
  cv_random_networks[i] <- raster::cv(net_rand)
  if(i %% 1000 == 0) {print(i)}
}
#cv_random_networks
print(paste0('network permutations for entire network completed at ', Sys.time()))

#cv_random_nodes <- rep(0,N)
#for (i in c(1:N)) {
#  net_rand <- sna::rmperm(gbi_matrix) --> TRY THIS AGAIN BUT RUN N TIMES IN ONE LOOP SO IT IS SEQUENTIALLY SWAPPING NODES RATHER THAN STARTING OVER EVERY TIME, AND THEN RUN LOOP AGAIN TO MEASURE CV
#  cv_random_nodes[i] <- raster::cv(net_rand)
#  if(i %% 1000 == 0) {print(i)}
#}
#print(paste0('node permutations for entire network completed at ', Sys.time()))

# compare permuted networks to actual network
hist(cv_random_networks, xlim = c(min(cv_random_networks)-50,max(cv_random_networks)+50),
     main = 'permutations for network of all male elephants')
abline(v = raster::cv(m_mat), col = 'red')
text(round(raster::cv(m_mat),3), col = 'red', x = median(cv_random_networks)+50, y = 2000)

#hist(cv_random_nodes, xlim = c(min(cv_random_networks)-50,max(cv_random_networks)+50),
#     main = 'permutations for network of all male elephants')
#abline(v = raster::cv(m_mat), col = 'red')
#text(round(raster::cv(m_mat),3), col = 'red', x = median(cv_random_networks)+50, y = 2000)

### write out outputs for future reference
cv_random_networks <- as.data.frame(cv_random_networks)
write_csv(cv_random_networks, '../data_processed/motnp_networkpermutations_cv_weakprior.csv')
#cv_random_nodes <- as.data.frame(cv_random_nodes)
#write_csv(cv_random_nodes,    '../data_processed/motnp_nodepermutations_cv_weakprior.csv')

### check for only elephants seen at least 5 times
# trim down data to remove any seen less than 5 times
seen5 <- counts_df[counts_df$count_1 > 4 & counts_df$count_2 > 4,]
ids5 <- unique(c(seen5$id_1, seen5$id_2))
m_mat5 <- m_mat[rownames(m_mat) %in% ids5, colnames(m_mat) %in% ids5]
gbi_m5 <- gbi_matrix[,colnames(gbi_matrix) %in% ids5]

# calculate CV for network
raster::cv(m_mat5) # still way too high
summary(m_mat5)
which(is.nan(m_mat5) == TRUE)

# run network permutations
random_networks <- asnipe::network_permutation(association_data = gbi_m5,
                                               association_matrix = m_mat5,
                                               permutations = N_networks)
cv_random_networks <- rep(0,N_networks)
for (i in c(1:N_networks)) {
  net_rand <- random_networks[i,,]
  cv_random_networks[i] <- raster::cv(net_rand)
  if(i %% 1000 == 0) {print(i)}
}
#cv_random_networks
raster::cv(m_mat5)
print(paste0('network permutations for all ids seen at least 5 times completed at ', Sys.time()))

#cv_random_nodes <- rep(0,N)
#for (i in c(1:N)) {
#  net_rand <- sna::rmperm(gbi_matrix)
#  cv_random_nodes[i] <- raster::cv(net_rand)
#  if(i %% 1000 == 0) {print(i)}
#}
#print(paste0('node permutations for all ids seen at least 5 times completed at ', Sys.time()))

# compare permuted networks to actual network
hist(cv_random_networks, xlim = c(min(cv_random_networks)-50,max(cv_random_networks)+50),
     main = 'permutations for network of all male elephants seen at least 5 times')
abline(v = raster::cv(m_mat5), col = 'red')
text(round(raster::cv(m_mat5),3), col = 'red', x = median(cv_random_networks)+50, y = 2000)

#hist(cv_random_nodes, xlim = c(min(cv_random_networks)-50,max(cv_random_networks)+50),
#     main = 'permutations for network of all male elephants seen at least 5 times')
#abline(v = raster::cv(m_mat5), col = 'red')
#text(round(raster::cv(m_mat5),3), col = 'red', x = median(cv_random_networks)+50, y = 2000)

### convert diagonal to 0s to check it doesn't change anything
#raster::cv(m_mat)
#for(i in 1:nrow(m_mat)){
#  m_mat[i,i] <- 0
#}
#raster::cv(m_mat)

### convert diagonal to NAs to check it doesn't change anything
#for(i in 1:nrow(m_mat)){
#  m_mat[i,i] <- NA
#}
#raster::cv(m_mat)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
mw_chain <- as.data.frame(motnp_edge_weights$chain)
colnames(mw_chain) <- counts_df$dyad_id
mw_chain <- pivot_longer(mw_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
mw_chain$chain_position <- rep(1:length(unique(mw_chain$dyad_id)), each = 4000)
mw_chain$draw <- LaplacesDemon::invlogit(mw_chain$draw)
mw_chain$mean <- NA
hist(mw_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using weak prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(mw_chain$dyad_id))){
  x <- mw_chain[mw_chain$dyad_id == i,]
  mw_chain$mean <- ifelse(mw_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), mw_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(mw_chain$draw), lwd = 2)

quantile(mw_chain$draw, 0.98)

mw_edgelist <- bisonR::get_edgelist(motnp_edge_weights)
quantile(mw_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### save pdf
dev.off()

### add time marker
print(paste0('original weak-priored model completed at ', Sys.time()))

## clean up ####
rm(list= ls()[!(ls() %in% c('motnp_edge_weights','motnp_edges_null',
                            'counts_df','counts_df_model',
                            'gbi_males','m_mat',
                            #'random_nodes',
                            'random_networks'))])

### write over saved workspace image -- only keep the things you might need later that take a long time to run
save.image('motnp_bisonr_edgescalculated_weakprior.RData')

#### edge weights -- stronger priors, males only ####
# add pdf output file
pdf(file = '../outputs/motnp_bisonr_edgeweight_strongprior.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[, c('node_1_males','node_2_males','event_count','count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### set priors
priors <- get_default_priors('binary')
priors$edge <- 'normal(-2.5, 1.5)'    # slightly right skewed so that more density is towards 0 but still very weak and all values are possible
priors$fixed <- 'normal(0, 1.5)'  # slightly less wide that before, but zero centred so allows age to have a positive or negative impact on centrality and edge weight
prior_check(priors, 'binary')

### run edge weight model
motnp_edge_weights_strongpriors <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(motnp_edge_weights_strongpriors, par_ids = 2)
plot_predictions(motnp_edge_weights_strongpriors, num_draws = 20, type = "density")
plot_predictions(motnp_edge_weights_strongpriors, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(motnp_edge_weights_strongpriors, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(motnp_edge_weights_strongpriors)

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
motnp_edges_null_strongpriors <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = motnp_edge_weights_strongpriors, random_model = motnp_edges_null_strongpriors)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging <- function(models) {
  loos <- lapply(models, function(model) loo::loo(model$log_lik, r_eff=NA))
  
  results_matrix <- loo::loo_model_weights(loos, method="pseudobma")
  if (!is.null(names(models))) {
    names(results_matrix) <- names(models)
  }
  results_matrix
}
model_averaging(models = list(non_random_model = motnp_edge_weights_strongpriors, random_model = motnp_edges_null_strongpriors))

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

### load in ages 
motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

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
plot_network_threshold(motnp_edge_weights_strongpriors, lwd = 2, ci = 0.9, threshold = 0.2,
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
plot_network_threshold2(obj = motnp_edge_weights_strongpriors, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)


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
# add pdf output file
pdf(file = '../outputs/motnp_bisonr_edgeweight_strongprior_cv.pdf')

# extract cv for model
global_cv_strongprior <- extract_metric(motnp_edge_weights_strongpriors, "global_cv")
head(global_cv_strongprior)
hist(global_cv_strongprior)

### standard deviation edge weight
#edge_weight <- extract_metric(motnp_edge_weights_strongpriors, "edge_weight")
#head(edge_weight)
#summary$std <- NA
#for(i in 1:nrow(summary)){
#  summary$std[i] <- sd(edge_weight[,i])
#}
#hist(summary$std)

#hist(summary$std/summary$median)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
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
for(i in sort(unique(ms_chain$dyad_id))){
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

### save pdf
dev.off()

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
rm(list= ls()[!(ls() %in% c('motnp_edge_weights_strongpriors','motnp_edges_null_strongpriors',
                            'counts_df','counts_df_model',
                            'gbi_males','m_mat',
                            #'random_nodes',
                            'random_networks'))])

### write over saved workspace image -- only keep the things you might need later that take a long time to run
save.image('motnp_bisonr_edgescalculated_strongprior.RData')

#### edge weights -- stronger priors, well sighted males only ####
# add pdf output file
pdf(file = '../outputs/motnp_bisonr_edgeweight_strongprior_mostsighted.pdf')

### create data frame for edge weight model
counts_df <- counts_df[counts_df$count_1 >= 10 & counts_df$count_2 >= 10,]
counts_df$node_1_wellsighted <- as.integer(as.factor(counts_df$id_1))
counts_df$node_2_wellsighted <- as.integer(as.factor(counts_df$id_2))+1
counts_df_model <- counts_df[, c('node_1_wellsighted','node_2_wellsighted','event_count','count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### set priors
priors <- get_default_priors('binary')
priors$edge <- 'normal(-2.5, 1.5)'    # slightly right skewed so that more density is towards 0 but still very weak and all values are possible
priors$fixed <- 'normal(0, 1.5)'  # slightly less wide that before, but zero centred so allows age to have a positive or negative impact on centrality and edge weight
prior_check(priors, 'binary')

### run edge weight model
motnp_edge_weights_strongprior_mostsighted <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(motnp_edge_weights_strongprior_mostsighted, par_ids = 2)
plot_predictions(motnp_edge_weights_strongprior_mostsighted, num_draws = 20, type = "density")
plot_predictions(motnp_edge_weights_strongprior_mostsighted, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(motnp_edge_weights_strongprior_mostsighted, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
motnp_edges_null_strongprior_mostsighted <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = motnp_edge_weights_strongprior_mostsighted, random_model = motnp_edges_null_strongprior_mostsighted))
# Method: stacking
#                weight
#non_random_model 0.227 
#random_model     0.773 
#Warning message: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging <- function(models) {
  loos <- lapply(models, function(model) loo::loo(model$log_lik, r_eff=NA))
  
  results_matrix <- loo::loo_model_weights(loos, method="pseudobma")
  if (!is.null(names(models))) {
    names(results_matrix) <- names(models)
  }
  results_matrix
}
model_averaging(models = list(non_random_model = motnp_edge_weights_strongprior_mostsighted, random_model = motnp_edges_null_strongprior_mostsighted)) # all weight still on random

# save workspace image
save.image('motnp_bisonr_edgescalculated_strongprior_mostsighted.RData')

### add time marker
print(paste0('random network comparison completed at ', Sys.time()))

## compare edge weight distributions to simple SRI ####
# load in workspace image with edge weights already calculated
#load('motnp_bisonr_edgescalculated_strongprior_mostsighted.RData')

### plot against SRI
head(edgelist)
colnames(edgelist)[1:2] <- colnames(counts_df_model)[1:2]
edgelist$node_1_id <- as.integer(edgelist$node_1_id) ; edgelist$node_2_id <- as.integer(edgelist$node_2_id)
summary <- left_join(edgelist, counts_df_model, by = c('node_1_id','node_2_id'))
counts <- counts_df[,c('node_1_males','node_2_males','count_1','count_2')]
colnames(counts)[1:2] <- c('node_1_id','node_2_id')
summary <- left_join(summary, counts, by = c('node_1_id','node_2_id'))
summary$sri <- summary$event / (summary$duration)

par(mar = c(4,4,3,1))
plot(density(summary$sri), las = 1, ylim = c(0,20),
     main = 'SRI vs model output for individuals seen ≥10 times:\nblack=sri, blue=grand median')
lines(density(summary$median), col = 'blue')

# try plotting a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
dyads <- counts_df[,c('dyad_id','node_1_wellsighted','node_2_wellsighted')]
colnames(dyads) <- c('dyad_id','node_1_id','node_2_id')
dyads <- left_join(dyads, counts_df_model, by = c('node_1_id','node_2_id'))
length(which(is.na(dyads$duration) == TRUE))

draws <- as.data.frame(motnp_edge_weights_strongprior_mostsighted$chain) %>% 
  pivot_longer(cols = everything(), names_to = 'dyad_model', values_to = 'edge')
draws$dyad_id <- rep(counts_df$dyad_id, 4000) ## IS THIS RIGHT??
draws$weight <- gtools::inv.logit(draws$edge)
draws$draw <- rep(1:4000,  each = nrow(counts_df_model))
draws <- left_join(draws, dyads, by = 'dyad_id')
draws$sri <- draws$event / draws$duration

set.seed(15)
subset_draws <- draws[draws$dyad_id %in% sample(draws$dyad_id, 64, replace = F),]
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
  facet_wrap(. ~ dyad_id, #scales = 'free_y',
             nrow = 8, ncol = 8)+
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

write_csv(subset_draws, '../data_processed/motnp_sampledyads_random_binary_vs_sri_strongprior_mostsighted.csv')

subset_draws <- draws[draws$sri > 0.15,]
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
  facet_wrap(. ~ dyad_id, nrow = 7, ncol = 6)+
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

write_csv(subset_draws, '../data_processed/motnp_sampledyads_sri0.15_binary_vs_sri_strongprior_mostsighted.csv')

# clean environment
rm(draws, dyads, priors, subset_draws, x) ; gc()
dev.off()

## coefficient of variation of edge weights (aka social differentiation) ####
# add pdf output file
pdf(file = '../outputs/motnp_bisonr_edgeweight_strongprior_mostsighted_cv.pdf')

# extract cv for model
global_cv_strongprior_mostsighted <- extract_metric(motnp_edge_weights_strongprior_mostsighted, "global_cv")
head(global_cv_strongprior_mostsighted)
hist(global_cv_strongprior_mostsighted)

### standard deviation edge weight
#edge_weight <- extract_metric(motnp_edge_weights_strongprior_mostsighted, "edge_weight")
#head(edge_weight)
#summary$std <- NA
#for(i in 1:nrow(summary)){
#  summary$std[i] <- sd(edge_weight[,i])
#}
#hist(summary$std)

#hist(summary$std/summary$median)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[,c('node_1_males','node_2_males','event_count','count_dyad','id_1','id_2','dyad_males')]
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ms_chain <- as.data.frame(motnp_edge_weights_strongprior_mostsighted$chain)
colnames(ms_chain) <- counts_df_model$dyad_id
ms_chain <- pivot_longer(ms_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ms_chain$chain_position <- rep(1:length(unique(ms_chain$dyad_id)), each = 4000)
ms_chain$draw <- LaplacesDemon::invlogit(ms_chain$draw)
ms_chain$mean <- NA
hist(ms_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(ms_chain$dyad_id))){
  x <- ms_chain[ms_chain$dyad_id == i,]
  ms_chain$mean <- ifelse(ms_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ms_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ms_chain$draw), lwd = 2)

quantile(ms_chain$draw, 0.98)

ms_edgelist <- bisonR::get_edgelist(motnp_edge_weights_strongprior_mostsighted)
quantile(ms_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### save pdf
dev.off()

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
rm(list= ls()[!(ls() %in% c('motnp_edge_weights_strongprior_mostsighted','motnp_edges_null_strongprior_mostsighted',
                            'counts_df','counts_df_model',
                            'gbi_males','m_mat',
                            #'random_nodes',
                            'random_networks'))])

### write over saved workspace image -- only keep the things you might need later that take a long time to run
save.image('motnp_bisonr_edgescalculated_strongprior_mostsighted.RData')



















## add pdf output file
pdf(file = '../outputs/motnp_bisonr_edgeweight_femalesonly.pdf')

## read in females data again ####
### import data for aggregated model (binomial)
counts_df <- read_csv('../data_processed/motnp_bayesian_binomialpairwiseevents.csv')
sex_1 <- data.frame(sex_1 = counts_df$id_1) %>% separate(sex_1, c("sex", "number"), sep = 1, remove = FALSE)
counts_df$sex_1 <- as.character(sex_1$sex) ; rm(sex_1) ; gc()
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

### remove dead individuals (also removes all unknown sex calves)
counts_df <- counts_df %>% 
  filter(name_1 != 'Richard' & name_2 != 'Richard') %>% 
  filter(name_1 != 'Gabriel' & name_2 != 'Gabriel') %>% 
  filter(name_1 != 'Tori'    & name_2 != 'Tori')

### filter counts data frame down to males only
females_df <- counts_df

### create data frame for edge weight model
motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')
motnp_ages <- motnp_ages[ ! motnp_ages$id %in% unique(females_df$id_1[females_df$age_cat_id_1 < 3]), ]
motnp_ages <- motnp_ages[ ! motnp_ages$id %in% females_df$id_2[females_df$age_class_2 < 3], ]
females_only <- females_df %>% 
  filter( ! id_1 %in% unique(motnp_ages$id)) %>% 
  filter( ! id_2 %in% unique(motnp_ages$id))
unique(females_only$dem_class_1)
unique(females_only$dem_class_2)
rm(counts_df, counts_df_model, females_df, motnp_ages, gbi_males, m_mat, motnp_edge_weights_strongprior_mostsighted, motnp_edges_null_strongprior_mostsighted) ; gc()
females_only_model <- females_only[, c('node_1','node_2','event_count','count_dyad')] %>% distinct()
colnames(females_only_model) <- c('node_1_id','node_2_id','event','duration')

### set priors
priors <- get_default_priors('binary')
priors$edge <- 'normal(-1, 2)'    # slightly right skewed so that more density is towards 0 but still very weak and all values are possible
priors$fixed <- 'normal(0, 1.5)'  # slightly less wide that before, but zero centred so allows age to have a positive or negative impact on centrality and edge weight
prior_check(priors, 'binary')

### run edge weight model
motnp_edge_weights_females <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = females_only_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(motnp_edge_weights_females, par_ids = 2)
plot_predictions(motnp_edge_weights_females, num_draws = 20, type = "density")
plot_predictions(motnp_edge_weights_females, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(motnp_edge_weights_females, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(motnp_edge_weights_females)

### add time marker
print(paste0('female-inclusive model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
motnp_edges_null_females <- bison_model(
  (event | duration) ~ 1, 
  data = females_only_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = motnp_edge_weights_females, random_model = motnp_edges_null_females)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging <- function(models) {
  loos <- lapply(models, function(model) loo::loo(model$log_lik, r_eff=NA))
  
  results_matrix <- loo::loo_model_weights(loos, method="pseudobma")
  if (!is.null(names(models))) {
    names(results_matrix) <- names(models)
  }
  results_matrix
}
model_averaging(models = list(non_random_model = motnp_edge_weights_females, random_model = motnp_edges_null_females))

### save workspace image
save.image('motnp_bisonr_edgescalculated_femalesonly.RData')
dev.off()

## plot network ####
### set up pdf
pdf(file = '../outputs/motnp_bisonr_edgeweight_femalesonly_networkplots_randomnesschecks.pdf')

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

### create nodes data frame
nodes <- data.frame(id_1 = sort(unique(c(females_only$id_1, females_only$id_2))))
nodes <- left_join(nodes, females_only[,c('id_1','age_category_1','count_1')]) %>% distinct()
nodes[nodes$id_1 == 'M89',] <- distinct(females_only[females_only$id_2 == 'M89',c('id_2','age_category_2','count_2')])
colnames(nodes) <- c('id','age','sightings')
nodes$age <- as.factor(as.integer(as.factor(nodes$age)))

### plot network
plot_network_threshold(motnp_edge_weights_females, lwd = 2, ci = 0.9, threshold = 0.2,
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
plot_network_threshold2(obj = motnp_edge_weights_females, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

### add time marker
print(paste0('network plots completed at ', Sys.time()))

### plot herds
females_chain <- as.data.frame(motnp_edge_weights_females$chain)
colnames(females_chain) <- females_only$dyad_id
females_chain <- pivot_longer(females_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
females_chain$chain_position <- rep(1:4000, each = nrow(females_only))
str(females_chain)
females_chain$dyad_id <- as.numeric(females_chain$dyad_id)
females_chain <- left_join(females_chain, females_only[,c('dyad_id','id_1','id_2','node_1','node_2','event_count','apart','count_dyad')], by = 'dyad_id')
head(females_chain)

unique(females_only$id_1)
test <- females_chain[females_chain$id_1 == 'F13' | females_chain$id_1 == 'F60' | 
                        females_chain$id_1 == 'F98' | females_chain$id_1 == 'F52' | 
                        females_chain$id_1 == 'F19' | females_chain$id_1 == 'F21' | 
                        females_chain$id_1 == 'M23' | females_chain$id_1 == 'M26',]
test <- test[test$id_2 == 'F60' | test$id_2 == 'F98' | test$id_2 == 'F52' | test$id_2 == 'F19' | test$id_2 == 'F21' | test$id_2 == 'M23' | test$id_2 == 'M26',]
head(test)

test$edge <- LaplacesDemon::invlogit(test$draw)
hist(test$edge)
test$dyad <- paste0(test$id_1,'_',test$id_2)
test$herd_1 <- ifelse(test$id_1 == 'F52' | test$id_1 == 'F60' | test$id_1 == 'F98', 'breeding herd 7', 'breeding herd 2')
test$herd_2 <- ifelse(test$id_2 == 'F52' | test$id_2 == 'F60' | test$id_2 == 'F98', 'breeding herd 7', 'breeding herd 2')
test$herd <- ifelse(test$herd_1 == test$herd_2,
                    ifelse(test$herd_1 == 'breeding herd 7', 'breeding herd 7','breeding herd 2'), 'different herds')
ggplot(test, aes(x = edge, colour = dyad))+
  geom_density()
ggplot(test, aes(x = edge, colour = herd, group = dyad))+
  geom_density()

rm(test, females_chain) ; gc()

## compare edge weight distributions to simple SRI ####
### plot against SRI
head(edgelist)
colnames(edgelist)[1:2] <- colnames(females_only_model)[1:2]
edgelist$node_1_id <- as.integer(edgelist$node_1_id) ; edgelist$node_2_id <- as.integer(edgelist$node_2_id)
summary <- left_join(edgelist, females_only_model, by = c('node_1_id','node_2_id'))
counts <- females_only[,c('node_1','node_2','count_1','count_2')]
colnames(counts)[1:2] <- c('node_1_id','node_2_id')
summary <- left_join(summary, counts, by = c('node_1_id','node_2_id'))
summary$sri <- summary$event / (summary$duration)

plot(density(summary$sri), main = 'SRI vs model output: blue=all,\nred=both seen 8 times, green=both 12 times')
lines(density(summary$median), col = 'blue')
lines(density(summary$median[which(summary$count_1 >= 8 & summary$count_2 >= 8)]), col = 'red')
lines(density(summary$median[which(summary$count_1 >= 12 & summary$count_2 >= 12)]), col = 'green')

# try plotting a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
dyads <- females_only[,c('dyad_id','node_1','node_2')]
colnames(dyads) <- c('dyad_id','node_1_id','node_2_id')
dyads <- left_join(dyads, females_only_model, by = c('node_1_id','node_2_id'))
length(which(is.na(dyads$duration) == TRUE))

draws <- as.data.frame(motnp_edge_weights_females$chain) %>% pivot_longer(cols = everything(), names_to = 'dyad_model', values_to = 'edge')
draws$dyad_id <- rep(females_only$dyad_id, 4000) ## IS THIS RIGHT??
draws$weight <- gtools::inv.logit(draws$edge)
draws$draw <- rep(1:4000,  each = nrow(females_only_model))
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

write_csv(subset_draws, '../data_processed/motnp_sampledyads_random_binary_vs_sri_femalesonly.csv')

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
  facet_wrap(. ~ dyad_id, ncol = 15)+
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

ggplot(data = subset_draws, mapping = aes(x = weight))+
  geom_density(colour = 'blue')+
  facet_wrap(. ~ dyad_id, ncol = 15, scales = 'free_y')+
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

write_csv(subset_draws, '../data_processed/motnp_sampledyads_sri0.2_binary_vs_sri_femalesonly.csv')

# clean environment
dev.off()
rm(draws, dyads, priors, subset_draws, x) ; gc()

## coefficient of variation of edge weights (aka social differentiation) ####
# add pdf output file
pdf(file = '../outputs/motnp_bisonr_edgeweight_femalesonly_cv.pdf')

# extract cv for model
global_cv_femalesonly <- extract_metric(motnp_edge_weights_females, "global_cv")
head(global_cv_femalesonly)
hist(global_cv_femalesonly)

# calculate SRI for all dyads
females_only_model$sri <- females_only_model$event / females_only_model$duration
summary(females_only_model$sri)

# identify id numbers from node id values
f_nodes <- females_only[,c('node_1','node_2','id_1','id_2')]
colnames(f_nodes)[1:2] <- c('node_1_id','node_2_id')
f <- females_only_model %>% left_join(f_nodes, by = c('node_1_id','node_2_id'))
rm(f_nodes) ; gc()

# calculate CV of SRI for all dyads
raster::cv(females_only_model$sri)  # way too high, but lower than gbi_matrix
raster::cv(f$sri[f$sri > 0])        # still massive even when I remove the 0s

# generate matrix
N <- length(unique(c(f$id_1, f$id_2)))
ids <- unique(c(f$id_1, f$id_2))
f_mat <- diag(nrow = N)
colnames(f_mat) <- ids
rownames(f_mat) <- ids

# populate matrix with SRI values
for( i in 1:N ) {
  for( j in 1:N ) {
    if( i >= j ) {
      f_mat[i,j] <- f_mat[j,i]
    }
    else {
      id1 <- colnames(f_mat)[i]
      id2 <- rownames(f_mat)[j]
      f_mat[i,j] <- f$sri[which(f$id_1 == id1 & f$id_2 == id2)]
    }
  }
}

# calculate CV for matrix
raster::cv(f_mat) # still way too high
which(is.nan(f_mat) == TRUE)
sd(f_mat)/mean(f_mat) # matches raster version -- check it is doing what I think it is doing

### permute matrix to compare to random
# create matrix of only females
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
gbi_females <- gbi_matrix[,colnames(gbi_matrix) %in% ids]
gbi_females <- gbi_females[rowSums(gbi_females) > 0,] # remove male-only sightings

# run network permutations
N_networks <- 10000
random_networks <- asnipe::network_permutation(association_data = gbi_females,
                                               association_matrix = f_mat,
                                               permutations = N_networks)
cv_random_networks <- rep(0,N_networks)   # generate empty vector to fill with cv values
for (i in c(1:N_networks)) {
  net_rand <- random_networks[i,,]
  cv_random_networks[i] <- raster::cv(net_rand)
  if(i %% 1000 == 0) {print(i)}
}
#cv_random_networks

#cv_random_nodes <- rep(0,N)
#for (i in c(1:N)) {
#  net_rand <- sna::rmperm(gbi_matrix)
#  cv_random_nodes[i] <- raster::cv(net_rand)
#  if(i %% 1000 == 0) {print(i)}
#}

# compare permuted networks to actual network
hist(cv_random_networks, xlim = c(min(cv_random_networks)-50,max(cv_random_networks)+50),
     main = 'permutations for network of all female elephants')
abline(v = raster::cv(f_mat), col = 'red')
text(round(raster::cv(f_mat),3), col = 'red', x = median(cv_random_networks)+50, y = 2000)

#hist(cv_random_nodes, xlim = c(min(cv_random_networks)-50,max(cv_random_networks)+50),
#     main = 'permutations for network of all female elephants')
#abline(v = raster::cv(f_mat), col = 'red')
#text(round(raster::cv(f_mat),3), col = 'red', x = median(cv_random_networks)+50, y = 2000)

cv_random_networks <- as.data.frame(cv_random_networks)
write_csv(cv_random_networks, '../data_processed/motnp_networkpermutations_cv_femalesonly.csv')
#cv_random_nodes <- as.data.frame(cv_random_nodes)
#write_csv(cv_random_nodes, '../data_processed/motnp_nodepermutations_cv_femalesonly.csv')

### check for only elephants seen at least 5 times
# trim down data to remove any seen less than 5 times
seen5 <- females_only[females_only$count_1 > 4 & females_only$count_2 > 4,]
ids5 <- unique(c(seen5$id_1, seen5$id_2))
f_mat5 <- f_mat[rownames(f_mat) %in% ids5, colnames(f_mat) %in% ids5]
gbi_f5 <- gbi_matrix[,colnames(gbi_matrix) %in% ids5]

# calculate CV for network
raster::cv(f_mat5) # still way too high
summary(f_mat5)
which(is.nan(f_mat5) == TRUE)

# run network permutations
random_networks <- asnipe::network_permutation(association_data = gbi_f5,
                                               association_matrix = f_mat5,
                                               permutations = N_networks)
cv_random_networks <- rep(0,N_networks)
for (i in c(1:N_networks)) {
  net_rand <- random_networks[i,,]
  cv_random_networks[i] <- raster::cv(net_rand)
  if(i %% 1000 == 0) {print(i)}
}
#cv_random_networks
raster::cv(f_mat5)

#cv_random_nodes <- rep(0,N)
#for (i in c(1:N)) {
#  net_rand <- sna::rmperm(gbi_matrix)
#  cv_random_nodes[i] <- raster::cv(net_rand)
#  if(i %% 1000 == 0) {print(i)}
#}

# compare permuted networks to actual network
hist(cv_random_networks, xlim = c(min(cv_random_networks)-50,max(cv_random_networks)+50),
     main = 'permutations for network of all female elephants seen at least 5 times')
abline(v = raster::cv(f_mat5), col = 'red')
text(round(raster::cv(f_mat5),3), col = 'red', x = median(cv_random_networks)+50, y = 2000)

#hist(cv_random_nodes, xlim = c(min(cv_random_networks)-50,max(cv_random_networks)+50),
#     main = 'permutations for network of all female elephants')
#abline(v = raster::cv(f_mat), col = 'red')
#text(round(raster::cv(f_mat),3), col = 'red', x = median(cv_random_networks)+50, y = 2000)

### convert diagonal to 0s to check it doesn't change anything
#raster::cv(f_mat)
#for(i in 1:nrow(f_mat)){
#  f_mat[i,i] <- 0
#}
#raster::cv(f_mat)

### convert diagonal to NAs to check it doesn't change anything
#for(i in 1:nrow(f_mat)){
#  f_mat[i,i] <- NA
#}
#raster::cv(f_mat)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
females_only_model <- females_only_model %>%
  mutate(node_1 = node_1_id,
         node_2 = node_2_id) %>% 
  left_join(y = females_only[,c('node_1','node_2','id_1','id_2','dyad_id')],
            by = c('node_1','node_2'))
#females_only_model <- females_only_model[,c(1:4,7:9)]
f_chain <- as.data.frame(motnp_edge_weights_females$chain)
colnames(f_chain) <- females_only_model$dyad_id
f_chain <- pivot_longer(f_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
f_chain$chain_position <- rep(1:length(unique(f_chain$dyad_id)), each = 4000)
f_chain$draw <- LaplacesDemon::invlogit(f_chain$draw)
f_chain$mean <- NA
hist(f_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for female network', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sort(unique(f_chain$dyad_id))){
  x <- f_chain[f_chain$dyad_id == i,]
  f_chain$mean <- ifelse(f_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), f_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(f_chain$draw), lwd = 2)

quantile(f_chain$draw, 0.98)

f_edgelist <- bisonR::get_edgelist(motnp_edge_weights_females)
quantile(f_edgelist$median, 0.98)

females_only_model$sri <- females_only_model$event / females_only_model$duration
quantile(females_only_model$sri, 0.98)

### save pdf
dev.off()

### add time marker
print(paste0('female-inclusive model completed at ', Sys.time()))

## clean up ####
rm(list= ls()[!(ls() %in% c('motnp_edge_weights_females','motnp_edges_null_females',
                            'females_df','females_df_model',
                            'gbi_females','f_mat',
                            #'random_nodes',
                            'random_networks'))])

### write over saved workspace image -- only keep the things you might need later that take a long time to run
save.image('motnp_bisonr_edgescalculated_females.RData')
