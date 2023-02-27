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
counts_df_non0 <- read_csv('../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_shortwindows_impossiblepairsremoved.csv')

ate <- readxl::read_excel(path = '../data_raw/Raw_ATE_MaleSightingsCleaned_Fishlock220808.xlsx', sheet = 1) %>% janitor::clean_names()
ate$id <- paste('M',ate$casename, sep = '')
ate$node_id <- as.integer(as.factor(ate$casename))
ate$obs_date <- lubridate::as_date(ate$obs_date)
ate <- separate(ate, obs_time, into = c('wrong_date','correct_time'), remove = F, sep = ' ')
ate$correct_time_hms <- hms::as_hms(ate$correct_time)
ate$corrected_time <- lubridate::hour(ate$correct_time_hms)*60*60 + lubridate::minute(ate$correct_time_hms) + lubridate::second(ate$correct_time_hms)
lu <- function(x) { length(unique(x)) }
ate_nums <- tapply(X = ate$obs_num, INDEX = ate$obs_date, FUN = lu )
ate$obs_num <- ifelse(ate$obs_num == '0','00', ate$obs_num)
ate$obs_num <- ifelse(ate$obs_num == '0a','0A', ate$obs_num)
ate$obs_num <- ifelse(ate$obs_num == '0b','0B', ate$obs_num)
ate$obs_num <- ifelse(ate$obs_num == '1','01', ate$obs_num)
ate$obs_num_std <- NA
for(i in 1:length(ate)){
  date_row <- ate[ate$obs_date == ate$obs_date[i],]
  date_row$obs_num_std <- as.integer(as.factor(sort(date_row$obs_num)))
  ate$obs_num_std[i] <- date_row$obs_num_std[which(date_row$obs_id == ate$obs_id[i])[1]]
}
ate <- ate[,c(1:3,25,26,4,27,28,8,29,9:24)]
head(ate)
table(ate$num_bulls)
ate$num_bulls_recount <- NA
for(i in 1:nrow(ate)){
  x <- ate[ate$obs_id == ate$obs_id[i],]
  ate$num_bulls_recount[i] <- nrow(x)
}
mean(ate$num_bulls_recount) ; sd(ate$num_bulls_recount)
mean(ate$num_bulls) ; sd(ate$num_bulls)
rm(date_row, ate_nums, lu, i)

counts <- as.data.frame(table(ate$id)) ; colnames(counts) <- c('id','count')
median(counts$count)

sightings <- ate[,c('obs_id','obs_date','correct_time_hms','grid_code','obs_type_old')] %>% distinct()
length(unique(sightings$obs_id))         # 24386
sightings$distinct <- NA
for(i in 1:nrow(sightings)) {
  x <- sightings[sightings$obs_id == sightings$obs_id[i], ]
  sightings$distinct[i] <- nrow(x)
  if(is.na(sightings$obs_type_old[i]) == TRUE ) {
    obs_type <- unique(x$obs_type_old[!is.na(x$obs_type_old)])
    sightings$obs_type_old[i] <- ifelse(length(obs_type) == 1, obs_type, 'U')
  }
}
sightings <- distinct(sightings)
table(sightings$obs_type_old)
#sightings <- sightings[c(1:12,15:24176),]

males <- readxl::read_excel('data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()
males$id <- paste0('M',males$casename)
ids <- sort(unique(ate$id))
males <- males %>% dplyr::filter(id %in% ids)

# counts per male per time period
males$count_all <- NA
males$p1 <- NA ; males$p2 <- NA ; males$p3 <- NA ; males$p4 <- NA ; males$p5 <- NA ; males$p6 <- NA ; males$p7 <- NA
males$p8 <- NA ; males$p9 <- NA ; males$p10 <- NA ; males$p11 <- NA ; males$p12 <- NA ; males$p13 <- NA
males$p14 <- NA ; males$p15 <- NA ; males$p16 <- NA ; males$p17 <- NA ; males$p18 <- NA ; males$p19 <- NA
males$p20 <- NA ; males$p21 <- NA ; males$p22 <- NA ; males$p23 <- NA ; males$p24 <- NA ; males$p25 <- NA
males$p26 <- NA ; males$p27 <- NA ; males$p28 <- NA ; males$p29 <- NA ; males$p30 <- NA ; males$p31 <- NA
periods <- seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 32)
for(i in 1:nrow(males)){
  counts <- ate[ate$id == males$id[i],]
  males$count_all[i] <- nrow(counts)
  males$p1[i] <- length(which(counts$obs_date >= periods[1] & counts$obs_date < periods[2]))
  males$p2[i] <- length(which(counts$obs_date >= periods[2] & counts$obs_date < periods[3]))
  males$p3[i] <- length(which(counts$obs_date >= periods[3] & counts$obs_date < periods[4]))
  males$p4[i] <- length(which(counts$obs_date >= periods[4] & counts$obs_date < periods[5]))
  males$p5[i] <- length(which(counts$obs_date >= periods[5] & counts$obs_date < periods[6]))
  males$p6[i] <- length(which(counts$obs_date >= periods[6] & counts$obs_date < periods[7]))
  males$p7[i] <- length(which(counts$obs_date >= periods[7] & counts$obs_date < periods[8]))
  males$p8[i] <- length(which(counts$obs_date >= periods[8] & counts$obs_date < periods[9]))
  males$p9[i] <- length(which(counts$obs_date >= periods[9] & counts$obs_date < periods[10]))
  males$p10[i] <- length(which(counts$obs_date >= periods[10] & counts$obs_date < periods[11]))
  males$p11[i] <- length(which(counts$obs_date >= periods[11] & counts$obs_date < periods[12]))
  males$p12[i] <- length(which(counts$obs_date >= periods[12] & counts$obs_date < periods[13]))
  males$p13[i] <- length(which(counts$obs_date >= periods[13] & counts$obs_date < periods[14]))
  males$p14[i] <- length(which(counts$obs_date >= periods[14] & counts$obs_date < periods[15]))
  males$p15[i] <- length(which(counts$obs_date >= periods[15] & counts$obs_date < periods[16]))
  males$p16[i] <- length(which(counts$obs_date >= periods[16] & counts$obs_date < periods[17]))
  males$p17[i] <- length(which(counts$obs_date >= periods[17] & counts$obs_date < periods[18]))
  males$p18[i] <- length(which(counts$obs_date >= periods[18] & counts$obs_date < periods[19]))
  males$p19[i] <- length(which(counts$obs_date >= periods[19] & counts$obs_date < periods[20]))
  males$p20[i] <- length(which(counts$obs_date >= periods[20] & counts$obs_date < periods[21]))
  males$p21[i] <- length(which(counts$obs_date >= periods[21] & counts$obs_date < periods[22]))
  males$p22[i] <- length(which(counts$obs_date >= periods[22] & counts$obs_date < periods[23]))
  males$p23[i] <- length(which(counts$obs_date >= periods[23] & counts$obs_date < periods[24]))
  males$p24[i] <- length(which(counts$obs_date >= periods[24] & counts$obs_date < periods[25]))
  males$p25[i] <- length(which(counts$obs_date >= periods[25] & counts$obs_date < periods[26]))
  males$p26[i] <- length(which(counts$obs_date >= periods[26] & counts$obs_date < periods[27]))
  males$p27[i] <- length(which(counts$obs_date >= periods[27] & counts$obs_date < periods[28]))
  males$p28[i] <- length(which(counts$obs_date >= periods[28] & counts$obs_date < periods[29]))
  males$p29[i] <- length(which(counts$obs_date >= periods[29] & counts$obs_date < periods[30]))
  males$p30[i] <- length(which(counts$obs_date >= periods[30] & counts$obs_date < periods[31]))
  males$p31[i] <- length(which(counts$obs_date >= periods[31] & counts$obs_date < periods[32]))
}

# time windows
periods <- data.frame(period = 1:31,
                      period_start = seq(from = min(sightings$obs_date),
                                         to = max(sightings$obs_date),
                                         length.out = 32)[1:31],
                      period_end = seq(from = min(sightings$obs_date),
                                       to = max(sightings$obs_date),
                                       length.out = 32)[2:32])

# check out total sightings per individual/dyad per period and see if they are reasonable
table(counts_df_non0$period_count_1) ; hist(counts_df_non0$period_count_1)
table(counts_df_non0$period_count_2) ; hist(counts_df_non0$period_count_2)
table(counts_df_non0$period_count_dyad) ; hist(counts_df_non0$period_count_dyad, breaks = 30)
table(counts_df_non0$event_count) ; hist(counts_df_non0$event_count, breaks = 3)
# many are not seen much, but repeat sightings of pairs do seem to be genuinely rare

rm(ate, counts, sightings, i)

### add time marker
print(paste0('data read in at ', Sys.time()))

#### Period 1 ####
# add pdf output file
pdf(file = '../outputs/anp1_bisonr_edgeweight.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[, c('node_1_males','node_2_males','event_count','count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### set priors
priors <- get_default_priors('binary')
priors$edge <- 'normal(-2.5, 1.5)'    # slightly right skewed so that more density is towards 0 but still very weak and all values are possible
priors$fixed <- 'normal(0, 1.5)'  # slightly less wide that before, but zero centred so allows age to have a positive or negative impact on centrality and edge weight
prior_check(priors, 'binary')

### run edge weight model
anp1_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #partial_pooling = TRUE,
  #mc_cores = 4,
  priors = priors
)

### run diagnostic plots
plot_trace(anp1_edge_weights, par_ids = 2)
plot_predictions(anp1_edge_weights, num_draws = 20, type = "density")
plot_predictions(anp1_edge_weights, num_draws = 20, type = "point")

### extract edge weight summaries
edgelist <- get_edgelist(anp1_edge_weights, ci = 0.9, transform = TRUE)
plot(density(edgelist$median))
summary(anp1_edge_weights)

### add time marker
print(paste0('model completed at ', Sys.time()))

## non-random edge weights ####
### run null model
anp1_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

### compare null model with fitted model -- bisonR model stacking
model_comparison(list(non_random_model = anp1_edge_weights, random_model = anp1_edges_null)) # NETWORK IS RANDOM

### compare null model with fitted model -- Jordan pseudo model averaging
model_averaging <- function(models) {
  loos <- lapply(models, function(model) loo::loo(model$log_lik, r_eff=NA))
  
  results_matrix <- loo::loo_model_weights(loos, method="pseudobma")
  if (!is.null(names(models))) {
    names(results_matrix) <- names(models)
  }
  results_matrix
}
model_averaging(models = list(non_random_model = anp1_edge_weights, random_model = anp1_edges_null))

# save workspace image
save.image('anp1_bisonr_edgescalculated.RData')

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
anp1_ages <- readRDS('../data_processed/anp1_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

### create nodes data frame
nodes <- data.frame(bull = sort(unique(anp1_ages$id)),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  x <- anp1_ages[anp1_ages$id == nodes$bull[i],]
  nodes$age[i] <- mean(x$age)
  if(nodes$bull[i] != 'M99'){
    y <- counts_df[counts_df$id_1 == nodes$bull[i], c('id_1','count_1')]
    nodes$sightings[i] <- y[1,2]
  }
}
nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

### plot network
plot_network_threshold(anp1_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
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
plot_network_threshold2(obj = anp1_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)


### add time marker
print(paste0('network plots completed at ', Sys.time()))

## compare edge weight distributions to simple SRI ####
# load in workspace image with edge weights already calculated
#load('anp1_bisonr_edgescalculated.RData')

### plot against SRI
#head(edgelist)
#colnames(edgelist)[1:2] <- colnames(counts_df_model)[1:2]
#edgelist$node_1_id <- as.integer(edgelist$node_1_id) ; edgelist$node_2_id <- as.integer(edgelist$node_2_id)
#summary <- left_join(edgelist, counts_df_model, by = c('node_1_id','node_2_id'))
#counts <- counts_df[,c('node_1_males','node_2_males','count_1','count_2')]
#colnames(counts)[1:2] <- c('node_1_id','node_2_id')
#summary <- left_join(summary, counts, by = c('node_1_id','node_2_id'))
#summary$sri <- summary$event / (summary$duration)

#plot(density(summary$sri), main = 'SRI vs model output: blue=all,\nred=both seen 8 times, green=both 12 times')
#lines(density(summary$median), col = 'blue')
#lines(density(summary$median[which(summary$count_1 >= 8 & summary$count_2 >= 8)]), col = 'red')
#lines(density(summary$median[which(summary$count_1 >= 12 & summary$count_2 >= 12)]), col = 'green')

# try plotting a subset with facets showing the draw distributions and lines indicating the position of standard SRI calculation
#dyads <- counts_df[,c('dyad_id','node_1_males','node_2_males')]
#colnames(dyads) <- c('dyad_id','node_1_id','node_2_id')
#dyads <- left_join(dyads, counts_df_model, by = c('node_1_id','node_2_id'))
#length(which(is.na(dyads$duration) == TRUE))

#draws <- as.data.frame(anp1_edge_weights$chain) %>% pivot_longer(cols = everything(), names_to = 'dyad_model', values_to = 'edge')
#draws$dyad_id <- rep(counts_df$dyad_id, 4000) ## IS THIS RIGHT??
#draws$weight <- gtools::inv.logit(draws$edge)
#draws$draw <- rep(1:4000,  each = nrow(counts_df_model))
#draws <- left_join(draws, dyads, by = 'dyad_id')
#draws$sri <- draws$event / draws$duration

#set.seed(15)
#subset_draws <- draws[draws$dyad_id %in% sample(draws$dyad_id, 150, replace = F),]
#subset_draws$median <- NA
#for(i in 1:nrow(subset_draws)){
#  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]
#  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), subset_draws$median[i])
#}
#head(subset_draws)
#which(is.na(subset_draws$median) == TRUE)[1]

#subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)
#ggplot(data = subset_draws, mapping = aes(x = weight))+
#  geom_density(colour = 'blue')+
#  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+
#  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
#  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

#write_csv(subset_draws, '../data_processed/anp1_sampledyads_random_binary_vs_sri.csv')

#subset_draws <- draws[draws$sri > 0.2,]
#subset_draws$median <- NA
#for(i in 1:nrow(subset_draws)){
#  x <- subset_draws[subset_draws$dyad_id == subset_draws$dyad_id[i],]
#  subset_draws$median[i] <- ifelse(subset_draws$dyad_id[i] == x$dyad_id[1], median(x$weight), subset_draws$median[i])
#}
#head(subset_draws)
#which(is.na(subset_draws$median) == TRUE)[1]

#subset_draws$dyad_id <- reorder(subset_draws$dyad_id, subset_draws$duration)
#ggplot(data = subset_draws, mapping = aes(x = weight))+
#  geom_density(colour = 'blue')+
#  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15)+
#  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
#  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

#ggplot(data = subset_draws, mapping = aes(x = weight))+
#  geom_density(colour = 'blue')+
#  facet_wrap(. ~ dyad_id, nrow = 10, ncol = 15, scales = 'free_y')+
#  geom_vline(mapping = aes(xintercept = median), colour = 'blue', lty = 3)+
#  geom_vline(mapping = aes(xintercept = sri), colour = 'red')

#write_csv(subset_draws, '../data_processed/anp1_sampledyads_sri0.2_binary_vs_sri.csv')

# clean environment
#rm(draws, dyads, priors, subset_draws, x) ; gc()
#dev.off()

## coefficient of variation of edge weights (aka social differentiation) ####
# add pdf output file
#pdf(file = '../outputs/anp1_bisonr_edgeweight_cv.pdf')

# extract cv for model
global_cv <- extract_metric(anp1_edge_weights, "global_cv")
head(global_cv)
hist(global_cv)

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2%
counts_df_model <- counts_df[,c('node_1_males','node_2_males','event_count','count_dyad','id_1','id_2','dyad_males')]
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ms_chain <- as.data.frame(anp1_edge_weights$chain)
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

ms_edgelist <- bisonR::get_edgelist(anp1_edge_weights)
quantile(ms_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

### save pdf
dev.off()

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

## clean up ####
rm(list= ls()[!(ls() %in% c('anp1_edge_weights','anp1_edges_null',
                            'counts_df','counts_df_model',
                            'gbi_males','m_mat',
                            #'random_nodes',
                            'random_networks'))])

### write over saved workspace image -- only keep the things you might need later that take a long time to run
save.image('anp1_bisonr_edgescalculated.RData')

