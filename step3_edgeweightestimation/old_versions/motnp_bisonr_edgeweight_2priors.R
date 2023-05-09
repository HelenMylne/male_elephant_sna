#### set up ####
# load packages
# library(tidyverse) ; library(dplyr) ; library(cmdstanr) ; library(bisonR) ; library(asnipe) ; library(sna) ; library(raster)
library(tidyverse, lib.loc = '../packages/')     # library(tidyverse)
library(dplyr, lib.loc = '../packages/')         # library(dplyr)
#library(rstan, lib.loc = '../packages/')        # library(rstan)
library(cmdstanr, lib.loc = '../packages/')      # library(cmdstanr)
library(extraDistr, lib.loc = '../packages/')    # library(extraDistr)
library(bisonR, lib.loc = '../packages/')        # library(bisonR)
library(asnipe, lib.loc = '../packages/')        # library(asnipe)
library(sna, lib.loc = '../packages/')           # library(sna)
library(raster, lib.loc = '../packages/')        # library(raster)

# information
sessionInfo()
R.Version()
rstan::stan_version()

# set seed
set.seed(12345)

#### create data lists ####
### import data for aggregated model (binomial) -- counts of positive associations and total sightings
counts_df <- read_csv('../data_processed/motnp_binomialpairwiseevents_malesonly.csv')

### load in ages 
motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

### add time marker
print(paste0('data read in at ', Sys.time()))

#### edge weights -- stronger priors, males only ####
# add pdf output file
pdf(file = '../outputs/motnp_bisonr_edgeweight_strongprior.pdf')

### create data frame for edge weight model
counts_df_model <- counts_df[, c('node_1_males','node_2_males','event_count','count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

### separate priors for zero-inflated dyads (never seen together)
plot(density( LaplacesDemon::invlogit(rnorm(1000, -5, 3)) ), col = 'blue')      # new -- dyads never sighted together
lines(density( LaplacesDemon::invlogit(rnorm(1000, -2.5, 1.5)) ), col = 'red')   # original -- dyads ≥ 1 seen together 

### set priors for dyads seen together at any point
priors_seen <- get_default_priors('binary') # obtain structure for bison model priors
priors_seen$edge <- 'normal(-2.5,1.5)' 
prior_check(priors_seen, 'binary')

### set priors for dyads never seen together at any point
priors_unseen <- get_default_priors('binary') # obtain structure for bison model priors
priors_unseen$edge <- 'normal(-5,3)' 
prior_check(priors_unseen, 'binary')

### create dummy variable indicating if ever seen together or not -- use to select which prior to draw edge weight from
counts_df_model$seen_together <- ifelse(counts_df_model$event > 0, 1, 0)

### run edge weight model with prior for pairs seen together at least once
motnp_fit_seen <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id),   # count of sightings together given count of total sightings as a result of the individuals contained within the dyad
  data = counts_df_model,
  model_type = "binary",
  priors = priors_seen
)

### run edge weight model with prior for pairs never seen together
motnp_fit_unseen <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id),   # count of sightings together given count of total sightings as a result of the individuals contained within the dyad
  data = counts_df_model,
  model_type = "binary",
  priors = priors_unseen
)

### save workspace
save.image('motnp_bisonr_edgescalculated_2priors.RData')
# load('motnp_bisonr_edgescalculated_2priors.RData')

### combine models together
motnp_fit_edges <- motnp_fit_unseen           # new model with new name
edges_unseen <- motnp_fit_edges$edge_samples  # extract edge samples for unseen dyads
edges_seen <- motnp_fit_seen$edge_samples     # extract edge samples for observed dyads
dyads <- motnp_fit_edges$dyad_to_idx %>%      # extract dyad id numbers for comparisons
  as.data.frame()
colnames(dyads) <- c('node_1_id','node_2_id') # rename for joining
dyads$dyad_id_model <- row_number(dyads)      # add dyad id column
dyads <- left_join(dyads, counts_df_model[,c('node_1_id','node_2_id','seen_together')],
                   by = c('node_1_id','node_2_id'))  # join on seen_together data
for(i in 1:nrow(dyads)){
  if(dyads$seen_together[i] == 1){
    edges_unseen[,i] <- edges_seen[,i]        # replace edges for seen dyads with those from other model
  }
}

### make space in workspace
save.image('motnp_bisonr_edgescalculated_2priors.RData')
rm(motnp_fit_unseen, motnp_fit_seen) ; gc()

### run diagnostic plots
plot_trace(motnp_fit_edges, par_ids = 2)                             # trace plot
plot_predictions(motnp_fit_edges, num_draws = 20, type = "density")  # compare predictions to raw data -- predictions are more variable than observations, both overestimating the number of pairs only seen together a few times, but also allowing for together scores higher than observed
plot_predictions(motnp_fit_edges, num_draws = 20, type = "point")    # compare predictions to raw data -- anything below the line is predicting below SRI (e.g. upper end SHOULD be massively below line because we don't want scores of 1 from pairs seen once)

### compare edge weights to prior
edges <- as.data.frame(motnp_fit_edges$chain) %>%               # extract chain of values from model
  pivot_longer(cols = everything(), names_to = 'dyad', values_to = 'draw') %>%  # convert to long format
  mutate(dyad_id = rep(counts_df$dyad_id, 4000),                                # add column to identify dyad per draw
         draw = plogis(draw))                                                   # convert draws to proportion
head(edges)
plot(NULL, xlim = c(0,1), ylim = c(0,50), las = 1,                              # prepare plot window
     main = 'edge distributions',ylab = 'density', xlab = 'edge weight')
plot_seen <- sample(counts_df$dyad_id[which(counts_df$event_count > 0)], 5000, replace = F)      # sample 5000 dyads that were seen together to plot
plot_unseen <- sample(counts_df$dyad_id[which(counts_df$event_count == 0)], 5000, replace = F)   # sample 5000 dyads that were NOT seen together to plot
plot_samples <- c(plot_seen,plot_unseen)                                        # sample 10000 dyads to plot
for(dyad in plot_samples) {                                                     # plot probability density for sampled dyads
  x <- edges[edges$dyad_id == dyad,]
  lines(density(x$draw),
        col = ifelse(dyad %in% plot_seen, 'red','blue'))
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

#priors$edge <- 'normal(-2.5, 1.5)'                                              # set prior value (already set -- doesn't actually change anything, just in as a reminder of what it was)
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

