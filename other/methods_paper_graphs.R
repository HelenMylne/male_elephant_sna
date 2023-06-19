#### script for plotting graphs required for methods paper ####
library(cmdstanr)
library(tidyverse)
library(LaplacesDemon)
library(bisonR) # library(bisonR, lib.loc = '../packages/')

set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0/')

n_chains <- 4

# define plot theme
theme_set(theme_bw())

#### Figure 1a: 6 facet plot of 10 sightings vs 2 sightings, never together vs half together vs always together ####
# define sequence over which to plot
x <- seq(0, 1, length = 100)

# identify probability of each value of x depending on shape of distributions
never2  <- dbeta(x = x, shape1 = 1, shape2 = 5)
never10 <- dbeta(x = x, shape1 = 1, shape2 = 15)
sometimes2  <- dbeta(x = x, shape1 = 2, shape2 = 2)
sometimes10 <- dbeta(x = x, shape1 = 10, shape2 = 10)
always2  <- dbeta(x = x, shape1 = 5, shape2 = 1)
always10 <- dbeta(x = x, shape1 = 15, shape2 = 1)

# tidy
data <- data.frame(x = x,
                   never2 = never2, never10 = never10,
                   sometimes2 = sometimes2, sometimes10 = sometimes10,
                   always2 = always2, always10 = always10) %>% 
  pivot_longer(names_to = 'category', values_to = 'dbeta',
               cols = c("never2", "never10",
                        "sometimes2", "sometimes10",
                        "always2", "always10")) %>% 
  mutate(sightings = ifelse(category == 'never2' |
                              category == 'sometimes2' |
                              category == 'always2',
                            2, 10),
         together = ifelse(category == 'never2' | category == 'never10',
                           'never',
                           ifelse(category == 'always2' | category == 'always10',
                                  'always', 'sometimes')),
         labels_together = factor(ifelse(together == 'never',
                                               'never seen together',
                                               ifelse(together == 'sometimes',
                                                      'seen together half the time',
                                                      'always seen together')),
                                        levels = c('never seen together',
                                                   'seen together half the time',
                                                   'always seen together')),
         labels_sightings = ifelse(sightings == 2, 'seen twice','seen 10 times'),
         sri = ifelse(together == 'never', 0,
                      ifelse(together == 'sometimes',0.5, 1))
         )

(g <- ggplot(data)+
  geom_line(aes(x = x, y = dbeta), lwd = 1.2, colour = rgb(0,0.5,0.6,1))+
  geom_vline(aes(xintercept = sri), lwd = 1.2)+
  scale_x_continuous(name = 'edge weight')+
  scale_y_continuous(name = 'probability density', limits = c(0,3.5))+
  facet_grid(labels_sightings ~ labels_together)) # , scales = 'free_y'

#### Figure 2: beta(2,2) prior ####
# calculate probability of all values of x assuming beta(2,2) prior shape
beta2.2 <- dbeta(x = x, shape1 = 2, shape2 = 2)
data <- data.frame(x = x, dbeta = beta2.2)

# plot
(g <- ggplot(data)+
    geom_line(aes(x = x, y = dbeta), lwd = 2, colour = rgb(0,0.5,0.6,1))+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'probability density'))

#### Figure 3: posterior distributions for F52 and M40 using beta(2,2) prior ####
# read in MOTNP data
motnp <- read_csv('../data_processed/motnp_bayesian_binomialpairwiseevents.csv') %>% 
  mutate(count_dyad = (count_1 + count_2) - event_count,
         apart = count_dyad - event_count)
counts_ls <- list(n_dyads = nrow(motnp),
                  apart = motnp$apart,
                  together = motnp$event_count)

# compile edge model
edge_binary <- cmdstan_model("other/methods_paper_graphs/beta_2.2_nopooling.stan")

# fit model
fit_edges_motnp <- edge_binary$sample(
  data = counts_ls, 
  chains = n_chains, 
  parallel_chains = n_chains)

# extract edges
edge_samples <- fit_edges_motnp$draws() %>% as.data.frame()
edge_samples <- edge_samples[,(n_chains+1):ncol(edge_samples)]

motnp$dyad <- paste0(motnp$id_1,'_',motnp$id_2)
colnames(edge_samples) <- rep(motnp$dyad, each = n_chains)

# select well sampled elephants/dyads of interest
hist(motnp$count_dyad[motnp$id_1 == 'F52' | motnp$id_2 == 'F52'])
hist(motnp$count_dyad[motnp$id_1 == 'M40' | motnp$id_2 == 'M40'])
edges_well_sampled <- edge_samples[,c('F52_F60','F52_U17','F52_F98',
                                      'F52_F8','F52_M15','F52_M1',
                                      'F60_M40','M40_U17','F98_M40',
                                      'F8_M40','M15_M40','M1_M40')] %>% 
  pivot_longer(cols = everything(), names_to = 'dyad', values_to = 'edge_sample') %>% 
  mutate(relationship = ifelse(dyad == 'F52_F60' |
                                 dyad == 'F52_U17' |
                                 dyad == 'F52_F98', 'related','unrelated'),
         partner_sex = case_when(dyad == 'F52_F60' ~ 'female',
                                 dyad == 'F52_U17' ~ 'unknown',
                                 dyad == 'F52_F98' ~ 'female',
                                 dyad == 'F52_F8' ~ 'female',
                                 dyad == 'F52_M15' ~ 'male',
                                 dyad == 'F52_M1' ~ 'male',
                                 dyad == 'F60_M40' ~ 'female',
                                 dyad == 'M40_U17' ~ 'unknown',
                                 dyad == 'F98_M40' ~ 'female',
                                 dyad == 'F8_M40' ~ 'female',
                                 dyad == 'M15_M40' ~ 'male',
                                 dyad == 'M1_M40' ~ 'male'),
         focal = case_when(dyad == 'F52_F60' ~ 'female',
                           dyad == 'F52_U17' ~ 'female',
                           dyad == 'F52_F98' ~ 'female',
                           dyad == 'F52_F8' ~ 'female',
                           dyad == 'F52_M15' ~ 'female',
                           dyad == 'F52_M1' ~ 'female',
                           dyad == 'F60_M40' ~ 'male',
                           dyad == 'M40_U17' ~ 'male',
                           dyad == 'F98_M40' ~ 'male',
                           dyad == 'F8_M40' ~ 'male',
                           dyad == 'M15_M40' ~ 'male',
                           dyad == 'M1_M40' ~ 'male'))
edges_well_sampled$relationship <- factor(edges_well_sampled$relationship,
                                          levels = c('unrelated','related'))

ggplot(edges_well_sampled)+
  geom_density(aes(x = edge_sample, group = dyad,
                   lty = relationship, colour = partner_sex),
               lwd = 1, key_glyph = draw_key_path)+
  scale_x_continuous(name = 'edge weight', limits = c(0,1))+
  scale_y_continuous(name = 'probability density')+
  scale_colour_viridis_d()+
  #facet_grid(relationship ~ focal)
  facet_grid(. ~ focal)+
  theme(legend.position = 'bottom')+
  guides(colour = guide_legend(title.position="top", title.hjust = 0.5),
         lty = guide_legend(title.position="top", title.hjust = 0.5))+
  labs(colour = 'sex of dyad partner',
       lty = 'relationship with dyad partner')

# select poorly sampled elephants/dyads of interest
edges_poor_sampled <- edge_samples[,c('F101_F132','F103_F132','F132_U10',
                                      'F132_M105','F132_M110','F132_M136',
                                      'F101_M161','F103_M161','M161_U10',
                                      'M105_M161','M110_M161','M136_M161')] %>% 
  pivot_longer(cols = everything(), names_to = 'dyad', values_to = 'edge_sample') %>% 
  mutate(relationship = 'unrelated',
         partner_sex = case_when(dyad == 'F101_F132' ~ 'female',
                                 dyad == 'F103_F132' ~ 'female',
                                 dyad == 'F132_M105' ~ 'male',
                                 dyad == 'F132_M110' ~ 'male',
                                 dyad == 'F132_M136' ~ 'male',
                                 dyad == 'F132_U10' ~ 'unknown',
                                 dyad == 'F101_M161' ~ 'female',
                                 dyad == 'F103_M161' ~ 'female',
                                 dyad == 'M105_M161' ~ 'male',
                                 dyad == 'M110_M161' ~ 'male',
                                 dyad == 'M136_M161' ~ 'male',
                                 dyad == 'M161_U10' ~ 'unknown'),
         focal = case_when(dyad == 'F101_F132' ~ 'female',
                           dyad == 'F103_F132' ~ 'female',
                           dyad == 'F132_M105' ~ 'female',
                           dyad == 'F132_M110' ~ 'female',
                           dyad == 'F132_M136' ~ 'female',
                           dyad == 'F132_U10' ~ 'female',
                           dyad == 'F101_M161' ~ 'male',
                           dyad == 'F103_M161' ~ 'male',
                           dyad == 'M105_M161' ~ 'male',
                           dyad == 'M110_M161' ~ 'male',
                           dyad == 'M136_M161' ~ 'male',
                           dyad == 'M161_U10' ~ 'male'))
edges_poor_sampled$relationship <- factor(edges_poor_sampled$relationship,
                                          levels = c('related','unrelated'))

ggplot(edges_poor_sampled)+
  geom_density(aes(x = edge_sample, group = dyad,
                   lty = relationship, colour = partner_sex),
               lwd = 1, key_glyph = draw_key_path)+
  scale_x_continuous(name = 'edge weight', limits = c(0,1))+
  scale_y_continuous(name = 'probability density')+
  scale_colour_viridis_d(limits = c('female', 'male','unknown'))+
  facet_grid(. ~ focal)+
  theme(legend.position = 'bottom')+
  guides(colour = guide_legend(title.position = "top", title.hjust = 0.5),
         lty = guide_legend(title.position = "top", title.hjust = 0.5))+
  labs(lty = 'relationship with dyad partner',
       colour = 'sex of dyad partner')

edges_well_sampled$sampling <- 'well sampled'
edges_poor_sampled$sampling <- 'poorly sampled'
edges_sampled <- rbind(edges_well_sampled, edges_poor_sampled)

ggplot(edges_sampled)+
  geom_density(aes(x = edge_sample, group = dyad,
                   lty = relationship, colour = partner_sex),
               lwd = 1, key_glyph = draw_key_path)+
  scale_x_continuous(name = 'edge weight', limits = c(0,1))+
  scale_y_continuous(name = 'probability density')+
  scale_colour_viridis_d()+
  facet_grid(sampling ~ focal)+
  theme(legend.position = 'bottom')+
  guides(colour = guide_legend(title.position="top", title.hjust = 0.5),
         lty = guide_legend(title.position="top", title.hjust = 0.5))+
  labs(colour = 'sex of dyad partner',
       lty = 'relationship with dyad partner')

ggplot(edges_sampled)+
  geom_density(aes(x = edge_sample, group = dyad,
                   lty = relationship, colour = partner_sex),
               lwd = 1, key_glyph = draw_key_path)+
  scale_x_continuous(name = 'edge weight', limits = c(0,1))+
  scale_y_continuous(name = 'probability density')+
  scale_colour_viridis_d()+
  facet_grid(sampling ~ focal, scales = 'free_y')+
  theme(legend.position = 'bottom')+
  guides(colour = guide_legend(title.position="top", title.hjust = 0.5),
         lty = guide_legend(title.position="top", title.hjust = 0.5))+
  labs(colour = 'sex of dyad partner',
       lty = 'relationship with dyad partner')

#### Figure 4: posterior distributions for F52 and M40 using beta(2,2) prior with partial pooling ####
# compile edge model
edge_binary <- cmdstan_model("other/methods_paper_graphs/beta_2.2_withpooling.stan")

# fit model
fit_edges_motnp <- edge_binary$sample(
  data = counts_ls, 
  chains = n_chains, 
  parallel_chains = n_chains)

# extract edges
edge_samples <- fit_edges_motnp$draws() %>% as.data.frame()
edge_samples <- edge_samples[,(n_chains+1):ncol(edge_samples)]

motnp$dyad <- paste0(motnp$id_1,'_',motnp$id_2)
colnames(edge_samples) <- rep(motnp$dyad, each = n_chains)

# select well sampled elephants/dyads of interest
hist(motnp$count_dyad[motnp$id_1 == 'F52' | motnp$id_2 == 'F52'])
hist(motnp$count_dyad[motnp$id_1 == 'M40' | motnp$id_2 == 'M40'])
edges_well_sampled <- edge_samples[,c('F52_F60','F52_U17','F52_F98',
                                      'F52_F8','F52_M15','F52_M1',
                                      'F60_M40','M40_U17','F98_M40',
                                      'F8_M40','M15_M40','M1_M40')] %>% 
  pivot_longer(cols = everything(), names_to = 'dyad', values_to = 'edge_sample') %>% 
  mutate(relationship = ifelse(dyad == 'F52_F60' |
                                 dyad == 'F52_U17' |
                                 dyad == 'F52_F98', 'related','unrelated'),
         partner_sex = case_when(dyad == 'F52_F60' ~ 'female',
                                 dyad == 'F52_U17' ~ 'unknown',
                                 dyad == 'F52_F98' ~ 'female',
                                 dyad == 'F52_F8' ~ 'female',
                                 dyad == 'F52_M15' ~ 'male',
                                 dyad == 'F52_M1' ~ 'male',
                                 dyad == 'F60_M40' ~ 'female',
                                 dyad == 'M40_U17' ~ 'unknown',
                                 dyad == 'F98_M40' ~ 'female',
                                 dyad == 'F8_M40' ~ 'female',
                                 dyad == 'M15_M40' ~ 'male',
                                 dyad == 'M1_M40' ~ 'male'),
         focal = case_when(dyad == 'F52_F60' ~ 'female',
                           dyad == 'F52_U17' ~ 'female',
                           dyad == 'F52_F98' ~ 'female',
                           dyad == 'F52_F8' ~ 'female',
                           dyad == 'F52_M15' ~ 'female',
                           dyad == 'F52_M1' ~ 'female',
                           dyad == 'F60_M40' ~ 'male',
                           dyad == 'M40_U17' ~ 'male',
                           dyad == 'F98_M40' ~ 'male',
                           dyad == 'F8_M40' ~ 'male',
                           dyad == 'M15_M40' ~ 'male',
                           dyad == 'M1_M40' ~ 'male'))
edges_well_sampled$relationship <- factor(edges_well_sampled$relationship,
                                          levels = c('unrelated','related'))

ggplot(edges_well_sampled)+
  geom_density(aes(x = edge_sample, group = dyad,
                   lty = relationship, colour = partner_sex),
               lwd = 1, key_glyph = draw_key_path)+
  scale_x_continuous(name = 'edge weight', limits = c(0,1))+
  scale_y_continuous(name = 'probability density')+
  scale_colour_viridis_d()+
  #facet_grid(relationship ~ focal)
  facet_grid(. ~ focal)+
  theme(legend.position = 'bottom')+
  guides(colour = guide_legend(title.position="top", title.hjust = 0.5),
         lty = guide_legend(title.position="top", title.hjust = 0.5))+
  labs(colour = 'sex of dyad partner',
       lty = 'relationship with dyad partner')

# select poorly sampled elephants/dyads of interest
edges_poor_sampled <- edge_samples[,c('F101_F132','F103_F132','F132_U10',
                                      'F132_M105','F132_M110','F132_M136',
                                      'F101_M161','F103_M161','M161_U10',
                                      'M105_M161','M110_M161','M136_M161')] %>% 
  pivot_longer(cols = everything(), names_to = 'dyad', values_to = 'edge_sample') %>% 
  mutate(relationship = 'unrelated',
         partner_sex = case_when(dyad == 'F101_F132' ~ 'female',
                                 dyad == 'F103_F132' ~ 'female',
                                 dyad == 'F132_M105' ~ 'male',
                                 dyad == 'F132_M110' ~ 'male',
                                 dyad == 'F132_M136' ~ 'male',
                                 dyad == 'F132_U10' ~ 'unknown',
                                 dyad == 'F101_M161' ~ 'female',
                                 dyad == 'F103_M161' ~ 'female',
                                 dyad == 'M105_M161' ~ 'male',
                                 dyad == 'M110_M161' ~ 'male',
                                 dyad == 'M136_M161' ~ 'male',
                                 dyad == 'M161_U10' ~ 'unknown'),
         focal = case_when(dyad == 'F101_F132' ~ 'female',
                           dyad == 'F103_F132' ~ 'female',
                           dyad == 'F132_M105' ~ 'female',
                           dyad == 'F132_M110' ~ 'female',
                           dyad == 'F132_M136' ~ 'female',
                           dyad == 'F132_U10' ~ 'female',
                           dyad == 'F101_M161' ~ 'male',
                           dyad == 'F103_M161' ~ 'male',
                           dyad == 'M105_M161' ~ 'male',
                           dyad == 'M110_M161' ~ 'male',
                           dyad == 'M136_M161' ~ 'male',
                           dyad == 'M161_U10' ~ 'male'))
edges_poor_sampled$relationship <- factor(edges_poor_sampled$relationship,
                                          levels = c('related','unrelated'))

ggplot(edges_poor_sampled)+
  geom_density(aes(x = edge_sample, group = dyad,
                   lty = relationship, colour = partner_sex),
               lwd = 1, key_glyph = draw_key_path)+
  scale_x_continuous(name = 'edge weight', limits = c(0,1))+
  scale_y_continuous(name = 'probability density')+
  scale_colour_viridis_d(limits = c('female', 'male','unknown'))+
  facet_grid(. ~ focal)+
  theme(legend.position = 'bottom')+
  guides(colour = guide_legend(title.position = "top", title.hjust = 0.5),
         lty = guide_legend(title.position = "top", title.hjust = 0.5))+
  labs(lty = 'relationship with dyad partner',
       colour = 'sex of dyad partner')

edges_well_sampled$sampling <- 'well sampled'
edges_poor_sampled$sampling <- 'poorly sampled'
edges_sampled <- rbind(edges_well_sampled, edges_poor_sampled)

ggplot(edges_sampled)+
  geom_density(aes(x = edge_sample, group = dyad,
                   lty = relationship, colour = partner_sex),
               lwd = 1, key_glyph = draw_key_path)+
  scale_x_continuous(name = 'edge weight', limits = c(0,1))+
  scale_y_continuous(name = 'probability density')+
  scale_colour_viridis_d()+
  facet_grid(sampling ~ focal)+
  theme(legend.position = 'bottom')+
  guides(colour = guide_legend(title.position="top", title.hjust = 0.5),
         lty = guide_legend(title.position="top", title.hjust = 0.5))+
  labs(colour = 'sex of dyad partner',
       lty = 'relationship with dyad partner')

ggplot(edges_sampled)+
  geom_density(aes(x = edge_sample, group = dyad,
                   lty = relationship, colour = partner_sex),
               lwd = 1, key_glyph = draw_key_path)+
  scale_x_continuous(name = 'edge weight', limits = c(0,1))+
  scale_y_continuous(name = 'probability density')+
  scale_colour_viridis_d()+
  facet_grid(sampling ~ focal, scales = 'free_y')+
  theme(legend.position = 'bottom')+
  guides(colour = guide_legend(title.position="top", title.hjust = 0.5),
         lty = guide_legend(title.position="top", title.hjust = 0.5))+
  labs(colour = 'sex of dyad partner',
       lty = 'relationship with dyad partner')

save.image('other/methods_paper_graphs/partialpooling.RData')

#### Figure 5: logit_gaussian(-1,2) prior ####
# calculate probability of all values of x assuming logit(gaussian(-1,2)) prior shape
logit_gaussian_weak <- dnorm(x = logit(x), mean = -1, sd = 2)
data <- data.frame(x = x, dnorm = logit_gaussian_weak)

# plot
(g <- ggplot(data)+
    geom_line(aes(x = x, y = dnorm), lwd = 2, colour = rgb(0,0.5,0.6,1))+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'probability density'))

#### Figure 6: median edge weight vs dyad sightings ####
load('other/methods_paper_graphs/partialpooling.RData')

### load new data
counts_df <- read_csv('../data_processed/motnp_binomialpairwiseevents_malesonly.csv')
counts_ls <- counts_df[,c('dyad_males','node_1','node_2','event_count','count_dyad')]
colnames(counts_ls) <- c('dyad_id','node_1_id','node_2_id','event','duration')

### set priors
priors <- get_default_priors('binary') # obtain structure for bison model priors
priors$edge <- 'normal(-1,2)' 
prior_check(priors, 'binary')

### run edge weight model
motnp_fit_edges <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id),   # count of sightings together given count of total sightings as a result of the individuals contained within the dyad
  data = counts_ls,
  model_type = "binary",
  priors = priors
)

### get median edge weights
colnames(counts_ls)[2:3] <- c('node_1','node_2')
edgelist <- get_edgelist(motnp_fit_edges)
edgelist$node_1 <- as.numeric(edgelist$node_1)
edgelist$node_2 <- as.numeric(edgelist$node_2)
edgelist <-  left_join(edgelist, counts_ls, by = c('node_1', 'node_2'))

### plot median edge against count_dyad
ggplot(edgelist)+
  geom_point(aes(x = duration, y = median), shape = 19,
             colour = rgb(68/255, 1/255, 84/255, 0.05))+
  geom_smooth(aes(x = duration, y = median),
              #method = 'glm', method.args = list(family = 'binomial'),
              #formula = cbind(event,duration) ~ duration,
              lwd = 2, colour = rgb(0,0.5,0.6,1))+
  scale_x_continuous(name = 'total sightings per dyad')+
  scale_y_continuous(name = 'median of edge weight distribution estimate')


### get full draw list
edges <- extract_metric(motnp_fit_edges, 'edge_weight') %>% 
  plogis() %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'dyad_v', values_to = 'edge_draw')
edges$dyad_id <- rep(counts_ls$dyad_id, 1000)
edges <- left_join(edges, edgelist, by = 'dyad_id')

### clean up environment before running to minimise size of global environment
rm(bad, edgelist, edges_sampled, fit_edges_motnp, motnp_fit_edges,
   data, edge_binary, edge_samples, edges_poor_sampled, edges_well_sampled, g,
   always10, always2, never10, never2, sometimes10, sometimes2,
   beta2.2, conditional1, conditional2, logit_gaussian_strong, logit_gaussian_weak, x) ; gc()

### plot full draws against count_dyad
ggplot(edges)+
  geom_point(aes(x = duration, y = edge_draw),
             colour = rgb(0,0.5,0.6,0.05))+
  #geom_point(aes(x = duration, y = median),
  #           colour = rgb(68/255, 1/255, 84/255, 0.1))+
  scale_x_continuous(name = 'total sightings per dyad')+
  scale_y_continuous(name = 'dyad edge weight')

### clean up environment
rm(edges)

#### Figure 7: logit_gaussian(-2.5,1.5) prior ####
# calculate probability of all values of x assuming logit(gaussian(-2.5,1.5)) prior shape
logit_gaussian_strong <- dnorm(x = logit(x), mean = -2.5, sd = 1.5)
data <- data.frame(x = x, dnorm = logit_gaussian_strong)

# plot
(g <- ggplot(data)+
    geom_line(aes(x = x, y = dnorm), lwd = 2, colour = rgb(0,0.5,0.6,1))+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'probability density'))

#### Figure 8a: explanation of why effect of sightings is exacerbated in eigenvector measures beyond edge weights ####
# build node list
nodes <- data.frame(node = 1:14,
                    count = seq(1,27, by = 2))

# build edge list
edgelist <- data.frame(node = rep(1:14, each = 14),
                       node_2 = rep(1:14, 14),
                       weight = NA) %>% 
  filter(node < node_2) %>% 
  mutate(dyad = 1:91) %>% 
  left_join(nodes, by = 'node') %>% 
  rename(node_1 = node, count_1 = count, node = node_2) %>% 
  left_join(nodes, by = 'node') %>% 
  rename(node_2 = node, count_2 = count)

# simulate weights
plot(NULL, xlim = c(0,1), xlab = 'edge_weight', ylim = c(0,10), ylab = 'density', las = 1)
for(i in 1:nrow(edgelist)){
  lines(density(rbeta(1000, 1, sum(edgelist$count_1[i], edgelist$count_2[i])/4)),
        col = rgb(0,0,1,0.2))
}
for(i in 1:nrow(edgelist)){
  edgelist$weight[i] <- list(
    rbeta(1000, 1, sum(edgelist$count_1[i], edgelist$count_2[i])/4)
    )
}

edges_unlist <- data.frame(node = rep(rep(1:14, each = 14),1000),
                           node_2 = rep(rep(1:14, 14), 1000),
                           weight = NA) %>% 
  filter(node < node_2) %>% 
  mutate(dyad = rep(1:91, 1000)) %>% 
  left_join(nodes, by = 'node') %>% 
  rename(node_1 = node, count_1 = count, node = node_2) %>% 
  left_join(nodes, by = 'node') %>% 
  rename(node_2 = node, count_2 = count)
for(i in 1:nrow(edgelist)){
  for(j in 1:1000){
    x <- unlist(edgelist$weight[i])
    edges_unlist$weight[edges_unlist$node_1 == edgelist$node_1[i] & 
                          edges_unlist$node_2 == edgelist$node_2[i]] <- x
  }
}

# plot edge weights
ggplot()+
  geom_density(data = edges_unlist[edges_unlist$node_1 != 1 &
                                     edges_unlist$node_2 != 14,],
               aes(x = weight, group = dyad), lwd = 1,
               colour = rgb(33/255,145/255,140/255))+
  geom_density(data = edges_unlist[edges_unlist$node_1 == 1,],
               aes(x = weight, group = dyad), lwd = 1,
               colour = rgb(68/255,1/255,84/255))+
  geom_density(data = edges_unlist[edges_unlist$node_2 == 14,],
               aes(x = weight, group = dyad), lwd = 1,
               colour = rgb(253/255,231/255,37/255))+
  scale_x_continuous(name = 'edge weight')+
  scale_y_continuous(name = 'probability density')

# develop actual eigenvector distributions based on synthesised edge weights
library(igraph)
nodes_eigen <- data.frame(node = rep(1:14, each = 1000),
                          count = rep(seq(1,27, by = 2), each = 1000),
                          eigen = NA,
                          draw = rep(1:1000, 14))
edges_unlist$draw <- rep(1:1000, each = 91)

# build eigenvector centrality list
for(draw in 1:1000){
  # extract edge weights into adjacency matrix
  #edgelist$weight <- edges[draw,]
  edges <- edges_unlist[edges_unlist$draw == draw,]
  adj_mat <- matrix(NA, nrow = 14, ncol = 14,
                    dimnames = list(1:14, 1:14))
  for(i in 1:nrow(adj_mat)){
    for(j in 1:ncol(adj_mat)){
      if(i == j) {
        adj_mat[i,j] <- 0
      } else {
        if(i < j){
          adj_mat[i,j] <- edges$weight[edgelist$node_1 == i & edgelist$node_2 == j]
        } else {
          adj_mat[i,j] <- edges$weight[edgelist$node_1 == j & edgelist$node_2 == i]
        }
        
      }
    }
  }
  # build network
  net <- igraph::graph_from_adjacency_matrix(adjmatrix = adj_mat,
                                             weighted = T, mode = 'undirected')
  # extract eigen centralities
  eigen <- igraph::eigen_centrality(net, directed = F, weights = E(net)$weight)
  eigen <- eigen$vector
  # store eigen centralities
  nodes_eigen$eigen[nodes_eigen$draw == draw] <- eigen
}

# plot
ggplot(nodes_eigen)+
  geom_density(aes(x = eigen, group = as.factor(node),
                   colour = as.factor(node)), lwd = 1)+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()

ggplot()+
  geom_density(data = nodes_eigen[nodes_eigen$node != 1 &
                                    nodes_eigen$node != 14,],
               aes(x = eigen, group = as.factor(node)),
               lwd = 1, colour = rgb(33/255, 145/255, 140/255))+
  geom_density(data = nodes_eigen[nodes_eigen$node == 1,],
               aes(x = eigen),
               lwd = 1.5, colour = rgb(68/255, 1/255, 84/255))+
  geom_density(data = nodes_eigen[nodes_eigen$node == 14,],
               aes(x = eigen),
               lwd = 1.5, colour = rgb(253/255, 231/255, 37/255))+
  theme(legend.position = 'none')+
  scale_x_continuous(name = 'eigenvector centrality')+
  scale_y_continuous(name = 'probability density')

#### Figure 8b: eigenvector vs sighting count for logit_gaussian(-2.5,1.5) ####
load('motnp_bisonr_edgescalculated_strongprior.RData')

# create nodes data frame
nodes <- data.frame(node_1_males = unique(c(counts_df$node_1_males, counts_df$node_2_males))) %>% 
  left_join(distinct(counts_df[,c('node_1_males','age_category_1','count_1')]), by = 'node_1_males') %>% 
  mutate(node_2_males = node_1_males) %>% 
  left_join(distinct(counts_df[,c('node_2_males','age_category_2','count_2')]), by = 'node_2_males') %>% 
  mutate(node = node_1_males,
         age = ifelse(is.na(age_category_1) == FALSE, age_category_1, age_category_2),
         count = ifelse(is.na(count_1) == FALSE, count_1, count_2)) %>% 
  select(node, age, count)

# extract mean of posterior for eigenvector
mean_eigen_values <- extract_metric(motnp_edge_weights_strongpriors, "node_eigen") %>%
  as.data.frame() %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = everything()) %>%
  rename(bison_node_eigen=value) %>%
  mutate(node = nodes$node,
         id = nodes$id)

# extract posterior draws for eigenvector
eigen_values <- extract_metric(motnp_edge_weights_strongpriors, "node_eigen") %>%
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'name', values_to = 'eigen') %>% 
  left_join(mean_eigen_values, by = 'name') %>% 
  rename(mean_eigen = bison_node_eigen) %>% 
  select(-name) %>% 
  mutate(draw = rep(1:1000, each = length(unique(node))))

# combine
nodes <- left_join(nodes, mean_eigen_values, by = 'node')
nodes <- nodes %>% select(-name)

# plot eigenvector against sightings per individual
nodes %>% ggplot()+
  geom_point(aes(x = count, y = bison_node_eigen),
             colour = '#21918c')+
  theme(legend.position = 'none')+
  scale_x_continuous(name = 'number of observations per elephant')+
  scale_y_continuous(name = 'eigenvector centrality from model\n(mean per elephant)'#, limits = c(0,1)
                     )+
  geom_smooth(aes(x = count, y = bison_node_eigen),
              colour = '#440154')

eigen_values <- left_join(eigen_values, nodes, by = 'node')
eigen_values %>% ggplot()+
  geom_point(aes(x = count, y = eigen), shape = 19,
             colour = '#21918c', alpha = 0.05)+
  geom_point(aes(x = count, y = mean_eigen), shape = 19,
             colour = '#fde725', alpha = 0.05)+
  geom_smooth(aes(x = count, y = mean_eigen), colour = '#440154')+
  scale_x_continuous(name = 'number of observations per elephant')+
  scale_y_continuous(name = 'eigenvector centrality from model\n(all draws per elephant)')+
  theme(legend.position = 'none')

#### Figure 9: conditional prior ####
# calculate probability of all values of x assuming conditional prior shape
conditional1 <- dbeta(x = x, shape1 = 0.7, shape2 = 10)
conditional2 <- dbeta(x = x, shape1 = 1, shape2 = 5)
data <- data.frame(x = x,
                   dbeta1 = conditional1,
                   dbeta2 = conditional2)
ggplot(data)+
    geom_line(aes(x = x, y = dbeta1), lwd = 2, colour = rgb(0,0.5,0.6,1))+
    geom_line(aes(x = x, y = dbeta2), lwd = 2, colour = rgb(68/255,1/255,84/255,1))+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'probability density')

#### Figure 10a: posterior distributions using conditional prior ####
load('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')

if(length(which(cdf_1$event_count >= 1)) >= 200){ n_test <- 200 } else { n_test <- length(which(cdf_1$event_count >= 1)) } # set number of dyads to sample: will be equal numbers of dyads that have and have not ever been seen together, so if there are 200 or more dyads that have been seen together at some point then cap at 200, otherwise use as many as there are
plot_dyads <- c(sample(cdf_1$dyad_id[cdf_1$event_count >= 1], size = n_test, replace = F), # randomly sample dyads that have been seen together
                sample(cdf_1$dyad_id[cdf_1$event_count == 0], size = n_test, replace = F)) # randomly sample dyads never seen together
plot_edges <- edges[edges$dyad %in% plot_dyads,]                    # subset edge weight data using randomly drawn dyads
plot_edges$seen_together <- NA ; for(i in 1:length(plot_dyads)){    # set up for loop
  plot_edges$seen_together[plot_edges$dyad == plot_dyads[i]] <- ifelse(cdf_1$event_count[cdf_1$dyad_id == plot_dyads[i]] > 0, 1, 0) # the value of seen_together is 1 if the dyad has ever been seen in the same group, and 0 if they have not
}

### density plots
plot(NULL, xlim = c(0,1), ylim = c(0,30), las = 1, xlab = 'edge weight', ylab = 'density')         # set up plot window
for(i in 1:length(plot_dyads)){                                                                    # plot randomly sampled dyads
  x <- plot_edges[plot_edges$dyad == plot_dyads[i],]                                               # select data to plot
  lines(density(x$edge_draw), col = ifelse(x$seen_together == 1, rgb(0,0,1,0.1), rgb(1,0,0,0.1)))  # draw edge weight probability plot. blue = seen together at least once, red = never seen together
}

# MOTNP
rm(list = ls())
load('motnp_edgeweights_conditionalpriors.RData')

### density plots
plot(NULL, xlim = c(0,1), ylim = c(0,30), las = 1, xlab = 'edge weight', ylab = 'density')         # set up plot window
for(i in 1:length(plot_dyads)){                                                                    # plot randomly sampled dyads
  x <- plot_edges[plot_edges$dyad == plot_dyads[i],]                                               # select data to plot
  lines(density(x$edge_draw), col = ifelse(x$seen_together == 1, rgb(0,0,1,0.1), rgb(1,0,0,0.1)))  # draw edge weight probability plot. blue = seen together at least once, red = never seen together
}

ggplot(plot_edges)+
  geom_density(aes(x = edge_draw, colour = as.factor(seen_together),
                   group = as.factor(dyad)),
               linewidth = 1)+
  scale_colour_viridis_d(alpha = 0.5)+
  scale_x_continuous(name = 'edge weight')+
  scale_y_continuous(name = 'probability density')+
  theme(legend.position = 'none')

plot_dyads <- c(sample(plot_edges$dyad[plot_edges$seen_together == 1],
                       size = 100, replace = F),
                sample(plot_edges$dyad[plot_edges$seen_together == 0],
                         size = 100, replace = F))
plot_edges <- edges[edges$dyad %in% plot_dyads,]                    # subset edge weight data using randomly drawn dyads
plot_edges$seen_together <- NA ; for(i in 1:length(plot_dyads)){    # set up for loop
  plot_edges$seen_together[plot_edges$dyad == plot_dyads[i]] <- ifelse(counts_df$event_count[counts_df$dyad_id == plot_dyads[i]] > 0, 1, 0) # the value of seen_together is 1 if the dyad has ever been seen in the same group, and 0 if they have not
}

ggplot(plot_edges)+
  geom_density(aes(x = edge_draw, colour = as.factor(seen_together),
                   group = as.factor(dyad)),
               linewidth = 1)+
  scale_colour_viridis_d(alpha = 0.5)+
  scale_x_continuous(name = 'edge weight')+
  scale_y_continuous(name = 'probability density')+
  theme(legend.position = 'none')

#### Figure 10b: eigenvector vs sighting count for conditional prior -- NOT RUN YET (copied directly from anp_nodalregression script) ####
### extract eigenvector centralities
n_eles <- length(unique(c(counts_df$id_1, counts_df$id_2)))
eigen <- array(data = NA, dim = c(n_eles, 4, n_samples*n_chains),
               dimnames = list(nodes$node, c('node','age','sightings','eigenvector'),
                               1:(n_samples*n_chains)))
edges$chain <- ifelse(edges$chain == 'chain1', 1, ifelse(edges$chain == 'chain2', 2, ifelse(edges$chain == 'chain3', 3, 4)))
edges$draw_id <- edges$position + (edges$chain-1) * 1000
colnames(edges)[3] <- 'dyad_id'
edges <- left_join(edges, counts_df[,c('dyad_id','node_1','node_2','id_1','id_2')], by = 'dyad_id')
for(draw in 1:(n_samples*n_chains)){
  eigen[,1,draw] <- nodes$node
  eigen[,2,draw] <- nodes$age
  eigen[,3,draw] <- nodes$sightings
  edges_drawn <- edges[edges$draw_id == draw,]
  adj_mat <- matrix(data = NA, nrow = n_eles, ncol = n_eles, dimnames = list(nodes$id, nodes$id))
  for(i in 1:n_eles){
    for(j in 1:n_eles){
      if(i == j) {
        adj_mat[i,j] <- 0
      } else {
        adj_mat[i,j] <- edges_drawn$edge_draw[which(edges_drawn$id_1 == rownames(adj_mat)[i] &
                                                      edges_drawn$id_2 == colnames(adj_mat)[j] |
                                                      edges_drawn$id_1 == rownames(adj_mat)[j] &
                                                      edges_drawn$id_2 == colnames(adj_mat)[i])]
      }
    }
  }
  network <- igraph::graph_from_adjacency_matrix(adjmatrix = adj_mat, diag = FALSE, mode = 'undirected', weighted = TRUE)
  eigen_values <- as.matrix(igraph::eigen_centrality(network, directed = FALSE)$vector)
  eigen[,4,draw] <- eigen_values[,1]
}

## check eigenvector against sightings
plot(NULL, xlim = c(0,max(nodes$sightings)), ylim = c(0,1),
     las = 1, xlab = 'sightings', ylab = 'eigenvector')
for(i in 1:n_samples){
  points(eigen[,4,i] ~ eigen[,3,i], pch = 19, col = rgb(0.5,0,1,0.1))
}
for(i in 1:n_eles){
  x <- eigen[i,,]
  points(mean(x[4,]) ~ x[3,1], pch = 19, col = 'yellow')
}

# ggplot
nodes$eigen_means <- NA
for(i in 1:n_eles){
  nodes$eigen_means[i] <- mean(eigen[i,4,1:400]) # gave up waiting for it to complete running all 4000 so only have the first 400 currently -- can run it again with all values included if and when it comes to final version
}
eigen_df <- as.data.frame(eigen[,4,]) %>% 
  pivot_longer(everything(), names_to = 'draw', values_to = 'eigen') %>% 
  mutate(node = rep(nodes$node, each = 4000)) %>% 
  left_join(nodes, by = 'node')

eigen_df$draw <- as.numeric(eigen_df$draw)
eigen_df <- eigen_df[eigen_df$draw < 401,] # remove once run all

ggplot(eigen_df)+
  geom_point(aes(y = eigen, x = sightings), colour = '#21918c', alpha = 0.05)+
  geom_point(aes(y = eigen_means, x = sightings), colour = '#fde725')+
  geom_smooth(aes(y = eigen_means, x = sightings), colour = '#440154')+
  scale_x_continuous(name = 'sightings')+
  scale_y_continuous(name = 'eigenvector centrality')

