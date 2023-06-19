#### Bayesian analysis of ALERT data ####
# Script to process association data from Mosi-Oa-Tunya National Park, Zambia.
# Data collected: 19th May 2016-21st December 2017
# Collected by: Mr David Youldon, Mr Dabwiso Sakala, Miss Helen Mylne and other volunteers/interns/facilitated research students working with ALERT during this time
# Data supplied by: ALERT and Mr David Youldon (11th August 2021) and via Volunteer Encounter (Bex Saunders, 19th October 2021)
# Data input: raw data provided by ALERT processed using script 22.01.13_ALERT_bayesian.R

#### set up ####
### install packages
# install.packages('tidyverse')
# install.packages('dplyr')
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages(c("StanHeaders","rstan"),type="source")
# install.packages("remotes") ; remotes::install_github("stan-dev/cmdstanr") ## OR USE # install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages('rethinking')
# install.packages('igraph')
# install.packages('dagitty')
# install.packages('janitor')
# install.packages('lubridate')
# install.packages('hms')
# install.packages('readxl')

### load packages
# library(tidyverse) ; library(dplyr) ; library(cmdstanr) ; library(dagitty) ; library(igraph) ; library(janitor) ; library(lubridate) ; library(hms) ; library(readxl)
# set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0/')
library(cmdstanr, lib.loc = '../packages/')    # library(cmdstanr)
library(tidyverse, lib.loc = '../packages/')   # library(tidyverse)
library(dplyr, lib.loc = '../packages/')       # library(dplyr)
#library(rstan, lib.loc = '../packages/')      # library(rstan)
#library(dagitty)
library(igraph, lib.loc = '../packages/')      # library(igraph)
library(janitor, lib.loc = '../packages/')     # library(janitor)
library(lubridate, lib.loc = '../packages/')   # library(lubridate)
library(readxl, lib.loc = '../packages/')      # library(readxl)

### set stan path
#set_cmdstan_path('/Users/helen/.cmdstan/cmdstan-2.32.1')
#set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0')

### set seed
set.seed(12345)

### create file of output graphs
pdf('../outputs/motnp_edgeweights_conditionalprior.pdf', width = 20, height = 20)

#### draw dyadic regression DAGS ####
### plot with full names
binom <- dagitty::dagitty("dag{
                         age_1 [exposure];
                         age_2 [exposure];
                         edgeweight [outcome];
                         relatedness_dyad [unobserved];
                         age_1 -> edgeweight <- age_2;
                         edgeweight <- relatedness_dyad;
                         }")
dagitty::coordinates(binom) <- list(x = c(age_1 = 0, age_2 = 2, edgeweight = 1, relatedness_dyad = 1),
                                    y = c(age_1 = 0, age_2 = 0, edgeweight = 1, relatedness_dyad = 2))
rethinking::drawdag(binom)

### plot with letters
binom <- dagitty::dagitty("dag{
                         A_1 [exposure];
                         A_2 [exposure];
                         EW [outcome];
                         R_d [unobserved];
                         A_1 -> EW <- A_2;
                         EW <- R_d;
                         }")
dagitty::coordinates(binom) <- list(x = c(A_1 = 0, A_2 = 2, EW = 1, R_d = 1),
                                    y = c(A_1 = 0, A_2 = 0, EW = 1, R_d = 2))
rethinking::drawdag(binom, radius = 6, cex = 1.6)

### clear environment and reset plot window
rm(binom)

#### complete data processing ####
### load in age data
motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

### import data for aggregated model (binomial) -- counts of positive associations and total sightings
counts_df <- read_csv('../data_processed/motnp_bayesian_binomialpairwiseevents.csv')

### correct sex_1, which has loaded in as a logical vector not a character/factor
unique(counts_df$sex_1) # FALSE or NA
sex_1 <- data.frame(sex_1 = counts_df$id_1)
sex_1 <- sex_1 %>% separate(sex_1, c("sex", "number"), sep = 1, remove = FALSE) ; unique(sex_1$sex) # F, M, U
counts_df$sex_1 <- as.character(sex_1$sex) ; rm(sex_1) ; gc()
str(counts_df)  # sex_1 still comes up as logical, but now contains the right levels

### create variable for age difference
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

### correct dem_class with corrected age classes
counts_df$dem_class_1 <- ifelse(counts_df$age_class_1 == 'Adult', paste('A',counts_df$sex_1, sep = ''),
                                ifelse(counts_df$age_class_1 == 'Pubescent', paste('P',counts_df$sex_1, sep = ''),
                                       ifelse(counts_df$age_class_1 == 'Juvenile', paste('J',counts_df$sex_1, sep = ''),
                                              paste('C',counts_df$sex_1, sep = ''))))
counts_df$dem_class_2 <- ifelse(counts_df$age_class_2 == 'Adult', paste('A',counts_df$sex_2, sep = ''),
                                ifelse(counts_df$age_class_2 == 'Pubescent', paste('P',counts_df$sex_2, sep = ''),
                                       ifelse(counts_df$age_class_2 == 'Juvenile', paste('J',counts_df$sex_2, sep = ''),
                                              paste('C',counts_df$sex_2, sep = ''))))

### combined dem_class of dyad
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

### add column for age difference between dyad
counts_df$age_diff <- abs(as.numeric(counts_df$age_cat_id_1) - as.numeric(counts_df$age_cat_id_2))

### add column for total number of sightings per pair
counts_df$count_dyad <- (counts_df$count_1 + counts_df$count_2) - counts_df$event_count  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.

### add column for total number of sightings per pair where they were NOT together
counts_df$apart <- counts_df$count_dyad - counts_df$event_count

### reassign dyad numbers to remove gaps
counts_df$node_1_nogaps <- as.integer(as.factor(counts_df$node_1))
counts_df$node_2_nogaps <- as.integer(as.factor(counts_df$node_2))+1

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
write_csv(counts_df, '../data_processed/motnp_binomialpairwiseevents_malesonly.csv')
# counts_df <- read_csv('../data_processed/motnp_binomialpairwiseevents_malesonly.csv')

### add time marker
print(paste0('data read in at ', Sys.time()))

##### create data list ####
### create nodes data frame
nodes <- data.frame(id = sort(unique(c(counts_df$id_1,counts_df$id_2))),  # all unique individuals
                   node = NA, age = NA, sightings = NA)                  # data needed on each
for(i in 1:nrow(nodes)){
  # extract data about individual from counts_df data frame
  if(nodes$id[i] %in% counts_df$id_1) {
    x <- counts_df[counts_df$id_1 == nodes$id[i], c('id_1','node_1','count_1')] %>% distinct()
  } else { x <- counts_df[counts_df$id_2 == nodes$id[i], c('id_2','node_2','count_2')] %>% distinct() }
  colnames(x) <- c('id','node','count')
  # extract age data from motnp_ages data frame
  y <- motnp_ages[motnp_ages$id == nodes$id[i],]
  # add individual data
  nodes$node[i] <- x$node
  nodes$sightings[i] <- x$count
  nodes$age[i] <- mean(y$age)                                                   # for initial test purposes, age = mean only
}
rm(x,y) ; gc()

### create data list
n_chains <- 4
n_samples <- 1000
n_dyads <- nrow(counts_df)
counts_ls <- list(
  n_dyads    = n_dyads,                      # total number of times one or other of the dyad was observed
  dyad_ids   = counts_df$dyad_id,            # identifier for each dyad
  together   = counts_df$event_count,        # count number of sightings seen together
  count_dyad = counts_df$count_dyad          # count total number of times seen
)

#### compile Stan model ####
edge_binary <- cmdstan_model("models/edge_binary_basic.stan")
edge_binary

#### run model ####
### fit model
fit_edges_motnp <- edge_binary$sample(
  data = counts_ls, 
  chains = n_chains, parallel_chains = n_chains,
  iter_warmup = n_samples, iter_sampling = n_samples)

### check model
fit_edges_motnp

### extract posterior samples
posterior_samples <- fit_edges_motnp$draws()

### break down into parameters
edge_weights_matrix <- posterior_samples[,,2:(nrow(counts_df)+1)]
rm(posterior_samples) ; gc()

### save edge samples
edges <- as.data.frame(edge_weights_matrix[,,1])
colnames(edges) <- c('chain1','chain2','chain3','chain4')
edges <- pivot_longer(edges, everything(), values_to = 'edge_draw', names_to = 'chain')
edges$dyad <- counts_ls$dyad_ids[1]
edges$position <- rep(1:n_samples, each = n_chains)
for(i in 2:n_dyads){
  x <- as.data.frame(edge_weights_matrix[,,i])
  colnames(x) <- c('chain1','chain2','chain3','chain4')
  x <- pivot_longer(x, everything(), values_to = 'edge_draw', names_to = 'chain')
  x$dyad <- counts_ls$dyad_ids[i]
  x$position <- rep(1:n_samples, each = n_chains)
  edges <- rbind(edges, x)
}
saveRDS(edges, '../data_processed/motnp_edgedistributions_conditionalprior.RDS')
#edges <- readRDS('../data_processed/motnp_edgedistributions_conditionalprior.RDS')

##### check outputs: edge weights ####
### Assign random set of columns to check
if(length(which(counts_df$event_count >= 1)) >= 200){ n_test <- 200 } else { n_test <- length(which(counts_df$event_count >= 1)) }
plot_dyads <- c(sample(counts_df$dyad_id[counts_df$event_count >= 1], size = n_test, replace = F),
                sample(counts_df$dyad_id[counts_df$event_count == 0], size = n_test, replace = F))
plot_edges <- edges[edges$dyad %in% plot_dyads,]
plot_edges$seen_together <- NA ; for(i in 1:length(plot_dyads)){
  plot_edges$seen_together[plot_edges$dyad == plot_dyads[i]] <- ifelse(counts_df$event_count[counts_df$dyad_id == plot_dyads[i]] > 0, 1, 0)
}

### build traceplots
ggplot(data = plot_edges[plot_edges$seen_together == 1,], aes(y = edge_draw, x = position, colour = chain))+
  geom_line()+
  facet_wrap(. ~ dyad)+
  theme_classic()+
  theme(legend.position = 'none', strip.background = element_blank(), strip.text = element_blank())
ggplot(data = plot_edges[plot_edges$seen_together == 0,], aes(y = edge_draw, x = position, colour = chain))+
  geom_line()+
  facet_wrap(. ~ dyad)+
  theme_classic()+
  theme(legend.position = 'none', strip.background = element_blank(), strip.text = element_blank())

#### density plots
plot(NULL, xlim = c(0,1), ylim = c(0,30), las = 1, xlab = 'edge weight', ylab = 'density')
for(i in 1:length(plot_dyads)){
  x <- plot_edges[plot_edges$dyad == plot_dyads[i],]
  lines(density(x$edge_draw), col = ifelse(x$seen_together == 1, rgb(0,0,1,0.1), rgb(1,0,0,0.1)))
}

#### check outputs: plot network ####
### create plotting function
plot_network_threshold <- function (edge_samples, dyad_data, lwd = 2, threshold = 0.3,
                                    label.colour = 'transparent', label.font = 'Helvetica', 
                                    node.size = 4, node.colour = 'seagreen1',
                                    link.colour1 = 'black', link.colour2 = rgb(0, 0, 0, 0.3))
{
  dyad_name <- do.call(paste, c(dyad_data[c("node_1", "node_2")], sep=" <-> "))
  edge_lower <- apply(edge_samples, 2, function(x) quantile(x, probs=0.025))
  edge_upper <- apply(edge_samples, 2, function(x) quantile(x, probs=0.975))
  edge_median <- apply(edge_samples, 2, function(x) quantile(x, probs=0.5))
  edge_list <- cbind(
    "median"=round(edge_median, 3), 
    "2.5%"=round(edge_lower, 3), 
    "97.5%"=round(edge_upper, 3)
  )
  rownames(edge_list) <- dyad_name
  edgelist <- as.data.frame(edge_list)
  edgelist$node_1 <- as.character(dyad_data$node_1)
  edgelist$node_2 <- as.character(dyad_data$node_2)
  edgelist <- edgelist[,c(4:5,1:3)]
  threshold_edges <- edgelist[edgelist$median >= threshold,]
  if(nrow(threshold_edges) == 0) { stop('No edges above threshold') }
  net <- igraph::graph_from_edgelist(as.matrix(threshold_edges[, 1:2]), directed = F)
  
  if(is.data.frame(node.size) == TRUE ) {
    nodes_list <- data.frame(node = rep(NA, length(unique(c(threshold_edges$node_1, threshold_edges$node_2)))),
                             sightings = NA)
    for(i in 1:nrow(nodes_list)){
      nodes_all <- rep(NA, 2*nrow(threshold_edges))  
      for(a in 1:2){
        for(b in 1:nrow(threshold_edges)){
          nodes_all[a + (b-1)*2] <- threshold_edges[b,a]
        }
      }
      nodes_list$node <- unique(nodes_all)
      nodes_list$sightings[i] <- nodes$sightings[which(nodes$node == nodes_list$node[i])]
    }
    node_sightings <- nodes_list$sightings*8 #log(nodes_list$sightings)*5
  } else { node_sightings <- node.size }
  
  if(is.data.frame(node.colour) == TRUE ) {
    nodes_list <- data.frame(node = rep(NA, length(unique(c(threshold_edges$node_1, threshold_edges$node_2)))),
                             age = NA)
    for(i in 1:nrow(nodes_list)){
      nodes_all <- rep(NA, 2*nrow(threshold_edges))  
      for(a in 1:2){
        for(b in 1:nrow(threshold_edges)){
          nodes_all[a + (b-1)*2] <- threshold_edges[b,a]
        }
      }
      nodes_list$node <- unique(nodes_all)
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
                      vertex.size = log(node_sightings)*5,
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

### create single matrix of edge samples
edge_samples <- matrix(data = NA, nrow = n_samples*n_chains, ncol = n_dyads)
for(j in 1:n_dyads){
  edge_samples[,j] <- edge_weights_matrix[,,j]
}
colnames(edge_samples) <- counts_df$dyad_id

### plot network
plot_network_threshold(edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.05,
                       node.size = nodes, node.colour = nodes, lwd = 15)
plot_network_threshold(edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.10,
                       node.size = nodes, node.colour = nodes, lwd = 15)
plot_network_threshold(edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.15,
                       node.size = nodes, node.colour = nodes, lwd = 15)
plot_network_threshold(edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.20,
                       node.size = nodes, node.colour = nodes, lwd = 15)
plot_network_threshold(edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.25,
                       node.size = nodes, node.colour = nodes, lwd = 15)
plot_network_threshold(edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.30,
                       node.size = nodes, node.colour = nodes, lwd = 15)

### clean workspace
rm(counts_ls, x, edge_weights_matrix, i, j) ; gc()

### save image
save.image('motnp_edgeweights_conditionalpriors.RData')
#load('motnp_edgeweights_conditionalpriors.RData')

#### check outputs: compare to simple SRI ####
### extract edgelist and calculate SRI
make_edgelist <- function (edge_samples, dyad_data)
{
  dyad_name <- do.call(paste, c(dyad_data[c("node_1", "node_2")], sep=" <-> "))
  edge_lower <- apply(edge_samples, 2, function(x) quantile(x, probs=0.025))
  edge_upper <- apply(edge_samples, 2, function(x) quantile(x, probs=0.975))
  edge_median <- apply(edge_samples, 2, function(x) quantile(x, probs=0.5))
  edge_list <- cbind(
    "median"=round(edge_median, 3), 
    "2.5%"=round(edge_lower, 3), 
    "97.5%"=round(edge_upper, 3)
  )
  rownames(edge_list) <- dyad_name
  edgelist <- as.data.frame(edge_list)
  edgelist$node_1 <- as.character(dyad_data$node_1)
  edgelist$node_2 <- as.character(dyad_data$node_2)
  edgelist <- edgelist[,c(4:5,1:3)]
  return(edgelist)
}
edgelist <- make_edgelist(edge_samples = edge_samples, dyad_data = counts_df)       # obtain edge list
head(edgelist)                                                                      # check structure of edgelist
colnames(edgelist)[1:2] <- c('node_1','node_2')                                     # rename for join function
edgelist$node_1 <- as.integer(edgelist$node_1)                                      # convert to integers
edgelist$node_2 <- as.integer(edgelist$node_2)                                      # convert to integers
summary <- left_join(edgelist, counts_df, by = c('node_1','node_2'))              # combine distribution data with raw counts
summary$sri <- summary$event_count / summary$count_dyad                             # calculate basic SRI

### plot median vs SRI
plot(density(summary$sri), main = 'red = SRI, blue = BISoN', col = 'red', lwd = 2, las = 1) # plot SRI distribution
lines(density(summary$median), col = 'blue', lwd = 2)                                       # distribution of median estimates
plot(density(summary$sri), main = 'SRI vs model output: blue=all,\nred=both seen 8 times, green=both 12 times') # plot SRI distribution
lines(density(summary$median), col = 'blue')                                                                    # distribution of median estimates
lines(density(summary$median[which(summary$count_1 >= 8 & summary$count_2 >= 8)]), col = 'red')                 # distribution of median estimates for frequently sighted males
lines(density(summary$median[which(summary$count_1 >= 12 & summary$count_2 >= 12)]), col = 'green')             # distribution of median estimates for frequently sighted males

# try plotting a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
plot_edges$dyad_id <- as.numeric(plot_edges$dyad)
draws <- left_join(plot_edges, summary, by = 'dyad_id')                      # combine to add node and median value data
which(is.na(draws) == TRUE)                                                  # check fully merged

draws$dyad_id <- reorder(draws$dyad_id, draws$count_dyad)                    # order data frame based on total sightings per pair
ggplot(data = draws, mapping = aes(x = edge_draw)) +                         # plot data frame
  geom_density(colour = 'blue') +                                            # plot density plots of probability distributions output per dyad
  facet_wrap(. ~ dyad_id, nrow = 20) +                                       # split by dyad, allow to vary height depending on dyad
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3) + # add line showing where the median estimate is
  geom_vline(mapping = aes(xintercept = sri), colour = 'red') +              # add line showing where the SRI value is
  theme(strip.text.x = element_blank(), strip.background = element_blank())  # remove facet labels so plots can be bigger

draws$dyad_id <- reorder(draws$dyad_id, draws$sri)                           # order dataframe based on SRI
ggplot(data = draws, mapping = aes(x = edge_draw))+                          # plot dataframe
  geom_density(colour = 'blue')+                                             # plot density plots of probability distirbtuions output per dyad
  facet_wrap(. ~ dyad_id, nrow = 20)+                                        # split by dyad, allow to vary height depending on dyad
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+  # add line showing where the median estimate is
  geom_vline(mapping = aes(xintercept = sri), colour = 'red') +              # add line showing where the SRI value is
  theme(strip.text.x = element_blank(), strip.background = element_blank())  # remove facet labels so plots can be bigger
ggplot(data = draws, mapping = aes(x = edge_draw))+                          # repeat but with free y axes
  geom_density(colour = 'blue')+                                             # plot density plots of probability distirbtuions output per dyad
  facet_wrap(. ~ dyad_id, nrow = 20, scales = 'free_y')+                     # split by dyad, allow to vary height depending on dyad
  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+  # add line showing where the median estimate is
  geom_vline(mapping = aes(xintercept = sri), colour = 'red')+               # add line showing where the SRI value is
  theme(strip.text.x = element_blank(), strip.background = element_blank())  # remove facet labels so plots can be bigger

write_csv(draws, '../data_processed/motnp_sampledyads_conditionalprior.csv') # save output for future reference

### clean environment
rm(draws, plot_edges, n_test, plot_dyads) ; gc()
save.image('motnp_edgeweights_conditionalprior.RData')

### end pdf
dev.off()

### clear workspace
rm(counts_df, n_dyads, counts_ls, n_eles, fit_edges_motnp, edges, plot_dyads, plot_edges, nodes, edge_weights_matrix, edge_samples)

