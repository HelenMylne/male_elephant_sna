#### script for plotting graphs required for methods paper ####
library(cmdstanr)
library(tidyverse)
library(LaplacesDemon)
library(bisonR) # library(bisonR, lib.loc = '../packages/')
library(patchwork)
library(igraph)

#set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0/')

n_chains <- 4

# define plot theme
theme_set(theme_bw(base_size = 12))

#### Figure 1: 6 facet plot of 10 sightings vs 2 sightings, never together vs half together vs always together ####
# define sequence over which to plot
x <- seq(0, 1, length = 100)

# identify probability of each value of x depending on shape of distributions -- I'VE CURRENTLY GONE WITH VALUES THAT LOOK GOOD IN THESE DISTIRBUTIONS, BUT SHOULD THESE BE PARTICULAR VALUES? I DON'T THINK THEY SHOULD BE BUT WANT TO CHECK
rare2  <- dbeta(x = x, shape1 = 1.2, shape2 = 3)
rare10 <- dbeta(x = x, shape1 = 1.9, shape2 = 10)
sometimes2  <- dbeta(x = x, shape1 = 2, shape2 = 2)
sometimes10 <- dbeta(x = x, shape1 = 10, shape2 = 10)
usual2  <- dbeta(x = x, shape1 = 3, shape2 = 1.2)
usual10 <- dbeta(x = x, shape1 = 10, shape2 = 1.9)

# tidy
distributions <- data.frame(x = x,
                            rare2 = rare2, rare10 = rare10,
                            sometimes2 = sometimes2, sometimes10 = sometimes10,
                            usual2 = usual2, usual10 = usual10) %>% 
  pivot_longer(names_to = 'category', values_to = 'dbeta',
               cols = c("rare2", "rare10",
                        "sometimes2", "sometimes10",
                        "usual2", "usual10")) %>% 
  mutate(sightings = ifelse(category == 'rare2' |
                              category == 'sometimes2' |
                              category == 'usual2',
                            2, 10),
         together = ifelse(category == 'rare2' | category == 'rare10',
                           'rare',
                           ifelse(category == 'usual2' | category == 'usual10',
                                  'usual', 'sometimes')),
         labels_together = factor(ifelse(together == 'rare',
                                         'together 10% of the time',
                                         ifelse(together == 'sometimes',
                                                'together 50% of the time',
                                                'together 90% of the time')),
                                  levels = c('together 10% of the time',
                                             'together 50% of the time',
                                             'together 90% of the time')),
         labels_sightings = ifelse(sightings == 2, 'seen twice','seen 10 times')
  )

sri <- data.frame(x = rep(seq(0, 1, length.out = 11), 6),
                  category = rep(unique(distributions$category), each = 11)) %>% 
  mutate(sightings = ifelse(category == 'rare2' |
                              category == 'sometimes2' |
                              category == 'usual2',
                            2, 10),
         observed = x*sightings,
         together = ifelse(category == 'rare2' | category == 'rare10',
                           'rare',
                           ifelse(category == 'usual2' | category == 'usual10',
                                  'usual', 'sometimes')),
         true_value = ifelse(together == 'rare', 0.1,
                             ifelse(together == 'sometimes', 0.5, 0.9)),
         labels_together = factor(ifelse(together == 'rare',
                                         'together 10% of the time',
                                         ifelse(together == 'sometimes',
                                                'together 50% of the time',
                                                'together 90% of the time')),
                                  levels = c('together 10% of the time',
                                             'together 50% of the time',
                                             'together 90% of the time')),
         labels_sightings = ifelse(sightings == 2, 'seen twice','seen 10 times')) %>% 
  mutate(possible = ifelse(sightings == 10, 'yes',
                           ifelse(x == 0 | x == 0.5 | x == 1, 'yes', 'no'))) %>% 
  filter(possible == 'yes') %>% 
  select(-possible) %>% 
  mutate(sri = dbinom(x = observed, size = sightings, prob = true_value),
         plot_y = 0) %>%
  pivot_longer(cols = c('sri', 'plot_y'), names_to = 'pairing', values_to = 'y_value') %>% 
  mutate(pairing = as.integer(as.factor(paste0(x, category))))

true_values <- data.frame(x = rep(c(0.1,0.5,0.9), each = 2),
                          together = rep(unique(distributions$together), each = 2),
                          sightings = rep(c(2,10), 3)) %>% 
  mutate(labels_together = factor(ifelse(together == 'rare',
                                         'together 10% of the time',
                                         ifelse(together == 'sometimes',
                                                'together 50% of the time',
                                                'together 90% of the time')),
                                  levels = c('together 10% of the time',
                                             'together 50% of the time',
                                             'together 90% of the time')),
         labels_sightings = ifelse(sightings == 2, 'seen twice','seen 10 times'))

(figure1 <- ggplot()+
    geom_vline(data = true_values, aes(xintercept = x), linewidth = 1.2, colour = '#fde725')+
    geom_line(data = distributions, aes(x = x, y = dbeta/max(dbeta)), linewidth = 1.2, colour = '#21918c')+
    geom_line(data = sri, aes(x = x, y = y_value, group = pairing), linewidth = 1.2, colour = '#440154')+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density', limits = c(0,1))+
    facet_grid(labels_sightings ~ labels_together) )
ggsave(filename = 'figure1_example_sri_vs_bison.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure1, device = 'png', width = 4200, height = 2400, units = 'px')

# clean up
rm(distributions, sri, true_values, rare10, rare2, sometimes10, sometimes2, usual10, usual2) ; gc()

#### Figure 2: 12 panel plots of SRI/default BISoN/right-skewed BISoN for prior/posterior/mean_edge_vs_sightings/eigenvector_vs_sightingss ####
## panel a: SRI, Prior distribution ####
# calculate probability of all values of x assuming completely flat prior
flat_prior <- dunif(x = x, min = 0, max = 1)
data <- data.frame(x = x, density = flat_prior)

# plot
(figure2a <- ggplot(data)+
    geom_line(aes(x = x, y = density), linewidth = 1.2, colour = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density', limits = c(0,1.05), expand = c(0,0))
)
ggsave(filename = 'figure2a_prior_sri.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2a, device = 'png', width = 1400, height = 800, units = 'px')

## panel b: default, Prior distribution ####
# calculate probability of all values of x assuming default bisonR prior shape
default_prior <- get_default_priors('binary') # default edge prior = normal(0, 2.5)
default_prior <- c(0, 2.5)
default_prior <- dnorm(x = logit(x), mean = default_prior[1], sd = default_prior[2])
data <- data.frame(x = x, density = default_prior)

# plot
(figure2b <- ggplot(data)+
    geom_line(aes(x = x, y = density), linewidth = 1.2, colour = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density',
                       limits = c(0,0.2),
                       #limits = c(0,1),
                       expand = c(0,0))
  )
ggsave(filename = 'figure2b_prior_default.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2b, device = 'png', width = 1400, height = 800, units = 'px')

## panel c: right-skewed, Prior distribution ####
# calculate probability of all values of x assuming right skewed prior shape
skewed_prior <- c(-2.5, 1.5)
skewed_prior <- dnorm(x = logit(x), mean = skewed_prior[1], sd = skewed_prior[2])
data <- data.frame(x = x, density = skewed_prior)

# plot
(figure2c <- ggplot(data)+
    geom_line(aes(x = x, y = density), linewidth = 1.2, colour = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density',
                       limits = c(0,0.3),
                       #limits = c(0,1),
                       expand = c(0,0))
  )
ggsave(filename = 'figure2c_prior_rightskewed.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2c, device = 'png', width = 1400, height = 800, units = 'px')

# clean up
rm(list = ls()[!ls() %in% c('figure2a','figure2b','figure2c')]) ; gc()

## panel d: SRI, Posterior distribution ####
# read in MOTNP data
load('motnp_edgeweights_conditionalprior.RData')
rm(edge_binary, edge_samples, edgelist, edges, fit_edges_motnp)

# calculate sri
counts_df$sri <- counts_df$event_count / counts_df$count_dyad

# calculate percentages of 0s and 1s using SRI
( length(which(counts_df$sri==0))/length(counts_df$sri) ) * 100
( length(which(counts_df$sri==1))/length(counts_df$sri) ) * 100
length(unique(c(counts_df$id_1, counts_df$id_2)))

# calculate percentages of 0s and 1s using SRI, assuming a 5 sighting threshold per elephant
subset5 <- counts_df %>% filter(count_1 >= 5 & count_2 >= 5)
( length(which(subset5$sri==0))/length(subset5$sri) ) * 100
( length(which(subset5$sri==1))/length(subset5$sri) ) * 100
length(unique(c(subset5$id_1, subset5$id_2)))

# calculate percentages of 0s and 1s using SRI, assuming a 10 sighting threshold per elephant
subset10 <- counts_df %>% filter(count_1 >= 10 & count_2 >= 10)
( length(which(subset10$sri==0))/length(subset10$sri) ) * 100
length(unique(c(subset10$id_1, subset10$id_2)))
100-(( length(unique(c(subset10$id_1, subset10$id_2))) / length(unique(c(counts_df$id_1, counts_df$id_2))) )*100)

# histogram posterior distributions
(figure2d <- ggplot(counts_df)+
    geom_histogram(aes(x = sri), bins = 100, fill = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'SRI value')+
    scale_y_continuous(name = 'number of dyads',
                       limits = c(0,20000),
                       expand = c(0,0))
  )
(figure2d <- ggplot(counts_df)+
    geom_histogram(aes(x = sri), bins = 100, fill = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'SRI value')+
    scale_y_continuous(name = 'number of dyads',
                       #limits = c(0,1000),
                       expand = c(0,0))+
    coord_cartesian(ylim = c(0,1000))
    )

# bar plot
(figure2d <- ggplot()+
    geom_bar(data = counts_df, aes(x = round(sri, 3)), fill = rgb(33/255, 145/255, 140/255),
             colour = rgb(33/255, 145/255, 140/255), linewidth = 0.2)+
    #geom_bar(data = counts_df[counts_df$count_dyad > 10,], aes(x = round(sri, 3)), fill = '#fde725', colour = 'transparent')+
    scale_x_continuous(name = 'SRI value')+
    scale_y_continuous(name = 'number of dyads',
                       #limits = c(0,1000),
                       expand = c(0,0))+
    annotate('text', x = 0.2, y = 175,
             label = paste0(length(which(counts_df$sri == 0)), ' dyads with \nSRI = 0'),
             colour = rgb(33/255, 145/255, 140/255))+
    coord_cartesian(ylim = c(0,200))
)
ggsave(filename = 'figure2d_posterior_sri.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2d, device = 'png', width = 1400, height = 800, units = 'px')

# redo annotation on figure 2d so that size isn't too large for multi plot
(figure2d <- ggplot()+
    geom_bar(data = counts_df, aes(x = round(sri, 3)), fill = rgb(33/255, 145/255, 140/255),
             colour = rgb(33/255, 145/255, 140/255), linewidth = 0.2)+
    #geom_bar(data = counts_df[counts_df$count_dyad > 10,], aes(x = round(sri, 3)), fill = '#fde725', colour = 'transparent')+
    scale_x_continuous(name = 'SRI value')+
    scale_y_continuous(name = 'number of dyads',
                       #limits = c(0,1000),
                       expand = c(0,0))+
    annotate('text', x = 0.23, y = 175,
             label = paste0(length(which(counts_df$sri == 0)), ' dyads with \nSRI = 0'),
             colour = rgb(33/255, 145/255, 140/255),
             size = unit(2, 'pt'))+
    coord_cartesian(ylim = c(0,200))
)

## panel e: default normal(0,2.5), Posterior distribution ####
# compile edge model
edge_binary_generic <- cmdstan_model('other/methods_paper/edge_binary_gaussian.stan')

# load data
counts_ls_default <- list(n_dyads = nrow(counts_df), 
                          dyad_ids = counts_df$dyad_id,
                          together = counts_df$event_count,
                          count_dyad = counts_df$count_dyad,
                          prior_mean = 0,
                          prior_stdev = 2.5)

# fit model
fit_edges_default <- edge_binary_generic$sample(
  data = counts_ls_default, 
  chains = n_chains, 
  parallel_chains = n_chains)

# extract edges
edge_samples <- fit_edges_default$draws() %>% as.data.frame()
edge_samples <- edge_samples[,(n_chains+1):ncol(edge_samples)]
edge_samples <- edge_samples[,seq(1,ncol(edge_samples)-3, by = 4)]
edge_samples <- edge_samples %>% 
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'weight')

# plot
(figure2e <- ggplot(data = edge_samples)+
    geom_density(aes(group = parameter, x = weight), colour = rgb(33/255, 145/255, 140/255, 0.1))+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density')
  )

#set.seed(15) ; plot_dyads <- sample(1:nrow(counts_df), size = 1000, replace = F)
#plot_dyad_ids <- unique(edge_samples$parameter)[plot_dyads]
#(figure2e <- ggplot(data = edge_samples[edge_samples$parameter %in% plot_dyad_ids,])+
#    geom_density(aes(group = parameter, x = weight), colour = rgb(33/255, 145/255, 140/255, 0.1))+
#    scale_x_continuous(name = 'edge weight')+
#    scale_y_continuous(name = 'density')
#)
ggsave(filename = 'figure2e_posterior_default.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2e, device = 'png', width = 1400, height = 800, units = 'px')

# clean workspace
save.image('../outputs/sparse_network_methods_figures/model_run_default.RData')

# clean up
rm(list = ls()[!ls() %in% c('figure2a','figure2b','figure2c','figure2d','figure2e',
                            'edge_binary_generic','counts_df','n_chains')]) ; gc()

## panel f: right-skewed normal(-2.5, 1.5), Posterior distribution ####
# compile edge model
# respecify priors for skewed model
counts_ls_skewed <- list(n_dyads = nrow(counts_df), 
                         dyad_ids = counts_df$dyad_id,
                         together = counts_df$event_count,
                         count_dyad = counts_df$count_dyad,
                         prior_mean = -2.5,
                         prior_stdev = 1.5)

# fit model
fit_edges_skewed <- edge_binary_generic$sample(
  data = counts_ls_skewed, 
  chains = n_chains, 
  parallel_chains = n_chains)

# extract edges
edge_samples <- fit_edges_skewed$draws() %>% as.data.frame()
edge_samples <- edge_samples[,(n_chains+1):ncol(edge_samples)]
edge_samples <- edge_samples[,seq(1,ncol(edge_samples)-3, by = 4)]
edge_samples <- edge_samples %>% 
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'weight')

# plot
(figure2f <- ggplot(data = edge_samples)+
    geom_density(aes(group = parameter, x = weight), colour = rgb(33/255, 145/255, 140/255, 0.1))+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density')
)
#(figure2f <- ggplot(data = edge_samples[edge_samples$parameter %in% plot_dyad_ids,])+
#    geom_density(aes(group = parameter, x = weight), colour = rgb(33/255, 145/255, 140/255, 0.1))+
#    scale_x_continuous(name = 'edge weight')+
#    scale_y_continuous(name = 'density')
#)
ggsave(filename = 'figure2f_posterior_rightskewed.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2f, device = 'png', width = 1400, height = 800, units = 'px')

# clean up and save
rm(data, default_prior, flat_prior, elephants) ; gc()
save.image('../outputs/sparse_network_methods_figures/model_run_skewed.RData')

# clean up
rm(list = ls()[!ls() %in% c('figure2a','figure2b','figure2c','figure2d','figure2e','figure2f')]) ; gc()

## merge first six panels ####
# run patchwork and save
(figure2a + figure2b + figure2c) / (figure2d + figure2e + figure2f) +
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'figure2_tophalf.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2100, height = 1200, units = 'px')

# clean up
#rm(figure2a, figure2b, figure2c, figure2d, figure2e, figure2f) ; gc()

## panel g: SRI, mean edge weight vs dyad sighting count ####
load('motnp_edgeweights_conditionalprior.RData')
rm(edge_binary, edge_samples, edgelist, edges, fit_edges_motnp)

# calculate sri
counts_df$sri <- counts_df$event_count / counts_df$count_dyad

# plot
(figure2g.1 <- ggplot(counts_df, aes(x = count_dyad, y = sri))+
   geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
   scale_x_continuous(name = 'total dyad sightings')+
   scale_y_continuous(name = 'SRI weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2g_edgesightings_sri_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2g.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure2g.2 <- ggplot(counts_df, aes(x = count_dyad, y = sri))+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'SRI weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2g_edgesightings_sri_withline_changealpha.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2g.2, device = 'png', width = 1400, height = 800, units = 'px')

counts_df$sri_rank <- as.integer(as.factor(counts_df$sri))
p7 <- as.data.frame(table(counts_df$count_dyad, counts_df$sri)) %>% 
  mutate(count_dyad = as.numeric(Var1),
         sri_rank = as.numeric(Var2)) %>% 
  filter(Freq > 0) %>% 
  left_join(distinct(counts_df[,c('sri','sri_rank')]), by = 'sri_rank')
(figure2g.3 <- ggplot()+
    geom_point(data = p7, aes(x = count_dyad, y = sri, size = Freq),
               colour = rgb(33/255, 145/255, 140/255))+
    geom_smooth(data = counts_df, aes(x = count_dyad, y = sri),
                colour = rgb(68/255, 1/255, 84/255))+
    theme(legend.position = c(0.8,0.7),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8))+
    scale_size_continuous(range = c(0.2, 2))+
    labs(size = 'number of dyads')+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'SRI weight')
)
ggsave(filename = 'figure2g_edgesightings_sri_withline_changesize.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2g.3, device = 'png', width = 1400, height = 800, units = 'px')

(figure2g.3 <- ggplot()+
    geom_point(data = p7, aes(x = count_dyad, y = sri, size = Freq),
               colour = rgb(33/255, 145/255, 140/255, 0.5))+
    geom_smooth(data = counts_df, aes(x = count_dyad, y = sri),
                colour = rgb(68/255, 1/255, 84/255))+
    theme(legend.position = c(0.7,0.7),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          legend.key.height = unit(2.5, 'mm'),
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 4))+
    scale_size_continuous(range = c(0.2, 2))+
    labs(size = 'number of dyads')+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'SRI weight')
)

# clean up 
rm(list = ls()[!ls() %in% c('figure2a', 'figure2b', 'figure2c',
                            'figure2d', 'figure2e', 'figure2f',
                            'figure2g.1','figure2g.2', 'figure2g.3')]) ; gc()

## panel h: default, mean edge weight vs dyad sighting count ####
load('../outputs/sparse_network_methods_figures/model_run_default.RData')
means <- edge_samples %>% 
  group_by(parameter) %>% 
  mutate(mean = mean(weight)) %>% 
  select(parameter, mean) %>% 
  distinct()

means$dyad_males <- 1:nrow(means)

means <- means %>% 
  left_join(counts_df[,c('dyad_id','node_1_males','node_2_males','event_count','count_dyad','count_1','count_2','dyad_males')],
            by = 'dyad_males')
rm(edge_samples) ; gc()

(figure2h.1 <- ggplot(means, aes(x = count_dyad, y = mean))+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mean weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2h_edgesightings_default_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2h.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure2h.2 <- ggplot(means, aes(x = count_dyad, y = mean))+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mean weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2h_edgesightings_default_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2h.2, device = 'png', width = 1400, height = 800, units = 'px')

# clean up
rm(list = ls()[!ls() %in% c('figure2a', 'figure2b', 'figure2c',
                            'figure2d', 'figure2e', 'figure2f',
                            'figure2g.1','figure2g.2','figure2g.3',
                            'figure2h.1','figure2h.2')]) ; gc()

## panel i: right-skewed, mean edge weight vs dyad sighting count ####
load('../outputs/sparse_network_methods_figures/model_run_skewed.RData')
means <- edge_samples %>% 
  group_by(parameter) %>% 
  mutate(mean = mean(weight)) %>% 
  select(parameter, mean) %>% 
  distinct()

means$dyad_males <- 1:nrow(means)

means <- means %>% 
  left_join(counts_df[,c('dyad_id','node_1_males','node_2_males','event_count','count_dyad','count_1','count_2','dyad_males')],
            by = 'dyad_males', multiple = 'all')
rm(edge_samples, fit_edges_skewed, edge_binary_generic) ; gc()

(figure2i.1 <- ggplot(means, aes(x = count_dyad, y = mean))+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mean weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2i_edgesightings_rightskewed_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2i.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure2i.2 <- ggplot(means, aes(x = count_dyad, y = mean))+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mean weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2i_edgesightings_rightskewed_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2i.2, device = 'png', width = 1400, height = 800, units = 'px')

# clean up
rm(list = ls()[!ls() %in% c('counts_df',
                            'figure2a', 'figure2b', 'figure2c',
                            'figure2d', 'figure2e', 'figure2f',
                            'figure2g.1','figure2g.2','figure2g.3',
                            'figure2h.1','figure2h.2',
                            'figure2i.1','figure2i.2')]) ; gc()

## panel j: SRI, eigenvector centrality vs total number of dyads where together = 0 ####
counts_df$sri <- counts_df$event_count / counts_df$count_dyad
elephants <- unique(c(counts_df$id_1, counts_df$id_2))
sri_adjmat <- matrix(NA, nrow = length(elephants), ncol = length(elephants), dimnames = list(x = elephants, y = elephants))
for(i in 1:nrow(sri_adjmat)){
  for(j in 1:ncol(sri_adjmat)){
    if(i <= j){
      sri_adjmat[i,j] <- ifelse(i == j, 0,
                                counts_df$sri[which(counts_df$id_1 == rownames(sri_adjmat)[i] &
                                                      counts_df$id_2 == colnames(sri_adjmat)[j])])
    } else {
      sri_adjmat[i,j] <- sri_adjmat[j,i]
    }
  }
}
sri_net <- igraph::graph_from_adjacency_matrix(sri_adjmat, weighted = T, mode = 'undirected', diag = FALSE)
eigen_sri <- igraph::eigen_centrality(sri_net)$vector %>% 
  as.data.frame()
colnames(eigen_sri) <- 'eigen'
eigen_sri$id <- rownames(eigen_sri)

together0 <- counts_df %>% filter(event_count == 0)
eigen_sri$together0 <- NA ; for(i in 1:nrow(eigen_sri)) {
  x <- together0 %>%
    filter(id_1 == eigen_sri$id[i] | id_2 == eigen_sri$id[i])
  eigen_sri$together0[i] <- nrow(x)
}

(figure2j.1 <- ggplot(eigen_sri, aes(x = together0, y = eigen))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2j_eigentogether0_sri_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2j.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure2j.2 <- ggplot(eigen_sri, aes(x = together0, y = eigen))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2j_eigentogether0_sri_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2j.2, device = 'png', width = 1400, height = 800, units = 'px')

## panel m: SRI, eigenvector centrality vs individual sighting count ####
counts <- data.frame(id = elephants, count = NA)
for(i in 1:nrow(counts)){
  counts$count[i] <- ifelse(i == nrow(counts),
                            unique(counts_df$count_2[which(counts_df$id_2 == counts$id[i])]),
                            unique(counts_df$count_1[which(counts_df$id_1 == counts$id[i])]))
}

eigen_sri <- eigen_sri %>% 
  left_join(counts, by = 'id')

(figure2m.1 <- ggplot(eigen_sri, aes(x = count, y = eigen))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2m_eigensightings_sri_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2m.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure2m.2 <- ggplot(eigen_sri, aes(x = count, y = eigen))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2m_eigensightings_sri_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2m.2, device = 'png', width = 1400, height = 800, units = 'px')

# clean up
rm(list = ls()[!ls() %in% c('figure2a', 'figure2b', 'figure2c',
                            'figure2d', 'figure2e', 'figure2f',
                            'figure2g.1','figure2g.2','figure2g.3',
                            'figure2h.1','figure2h.2',
                            'figure2i.1','figure2i.2',
                            'figure2j.1','figure2j.2',
                            'figure2m.1','figure2m.2', 'eigen_sri')]) ; gc()

## panel k: default, eigenvector centrality vs total number of dyads where together = 0 ####
load('../outputs/sparse_network_methods_figures/model_run_default.RData')
rm(edge_binary_generic, make_edgelist, plot_network_threshold) ; gc()

# extract edge samples as matrix
edge_samples <- fit_edges_default$draws() %>% as.data.frame()
edge_samples <- edge_samples[,(n_chains+1):ncol(edge_samples)]
edge_samples <- edge_samples[,seq(1,ncol(edge_samples)-3, by = 4)]
colnames(edge_samples) <- counts_df$dyad_id

elephants <- sort(unique(c(counts_df$id_1, counts_df$id_2)))
default_adjarr <- array(NA, dim = c(length(elephants), length(elephants), 1000),
                        dimnames = list(elephants, elephants, NULL))
N <- nrow(counts_df)

for (i in 1:N) {
  dyad_row <- counts_df[i, ]
  default_adjarr[dyad_row$id_1, dyad_row$id_2, ] <- edge_samples[, i]
}
default_adjarr[,,1]

# make matrix to save eigenvector values to
eigen_default <- matrix(NA, nrow = 1000, ncol = length(elephants),
                        dimnames = list(1:1000, elephants))

# calculate eigenvector scores
for(i in 1:nrow(eigen_default)){
  default_net <- igraph::graph_from_adjacency_matrix(default_adjarr[,,i], weighted = T, mode = 'undirected', diag = FALSE)
  eigen <- as.data.frame(igraph::eigen_centrality(default_net)$vector)
  eigen_default[i,] <- eigen$`igraph::eigen_centrality(default_net)$vector`
}
rm(default_net, dyad_row, eigen); gc()
eigen_default[,1]
eigen_default[1,]

eigen_default <- eigen_default %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'eigen') %>% 
  group_by(id) %>% 
  mutate(mean = mean(eigen)) %>% 
  left_join(eigen_sri[,c('id', 'together0', 'count')], by = 'id')

means <- eigen_default %>% 
  select(id, mean, together0, count) %>% 
  distinct()

(figure2k.1 <- ggplot()+
    geom_point(data = eigen_default, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = means, aes(x = together0, y = mean, colour = count),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2k_eigentogether0_default_allpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2k.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure2k.2 <- ggplot()+
    geom_point(data = eigen_default, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = means, aes(x = together0, y = mean),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(data = eigen_default, aes(x = together0, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2k_eigentogether0_default_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2k.2, device = 'png', width = 1400, height = 800, units = 'px')

(figure2k.3 <- ggplot(means)+
    geom_point(aes(x = together0, y = mean, colour = count),
               size = 0.5, shape = 19)+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))+
    labs(colour = 'number of\nsightings')
)
ggsave(filename = 'figure2k_eigentogether0_default_meanpoints_noline_coloursightings.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2k.3, device = 'png', width = 1400, height = 800, units = 'px')

(figure2k.4 <- ggplot(means)+
    geom_point(aes(x = together0, y = mean, colour = count),
               size = 0.5, shape = 19)+
    geom_smooth(aes(x = together0, y = mean),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))+
    labs(colour = 'number of\nsightings')
)
ggsave(filename = 'figure2k_eigentogether0_default_meanpoints_withline_colourssightings.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2k.4, device = 'png', width = 1400, height = 800, units = 'px')

## panel n: default, eigenvector centrality vs individual sighting count ####
(figure2n.1 <- ggplot()+
    geom_point(data = eigen_default, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = means, aes(x = count, y = mean),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
  )
ggsave(filename = 'figure2n_eigensightings_default_allpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2n.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure2n.2 <- ggplot()+
    geom_point(data = eigen_default, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = means, aes(x = count, y = mean),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(data = eigen_default, aes(x = count, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2n_eigensightings_default_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2n.2, device = 'png', width = 1400, height = 800, units = 'px')

(figure2n.3 <- ggplot(data = means, aes(y = mean,
                                        x = count))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean eigenvector centrality',
                       limits = c(0, 1),
                       expand = c(0.02,0.02))
)
ggsave(filename = 'figure2n_eigensightings_default_meanpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2n.3, device = 'png', width = 1400, height = 800, units = 'px')

(figure2n.4 <- ggplot(data = means, aes(x = count, y = mean))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean eigenvector centrality',
                       limits = c(0,1),
                       expand = c(0.02,0.02))
)
ggsave(filename = 'figure2n_eigensightings_default_meanpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2n.4, device = 'png', width = 1400, height = 800, units = 'px')

# clean up
rm(list = ls()[!ls() %in% c('figure2a','figure2b','figure2c',
                            'figure2d','figure2e','figure2f',
                            'figure2g.1','figure2g.2','figure2g.3',
                            'figure2h.1','figure2h.2',
                            'figure2i.1','figure2i.2',
                            'figure2j.1','figure2j.2',
                            'figure2k.1','figure2k.2',
                            'figure2m.1','figure2m.2',
                            'figure2n.1','figure2n.2','figure2n.3','figure2n.4',
                            'eigen_sri')]) ; gc()

## panel l: right-skewed, eigenvector centrality vs total number of dyads where together = 0 ####
load('../outputs/sparse_network_methods_figures/model_run_skewed.RData')

# extract edge samples as matrix
edge_samples <- fit_edges_skewed$draws() %>% as.data.frame()
edge_samples <- edge_samples[,(n_chains+1):ncol(edge_samples)]
edge_samples <- edge_samples[,seq(1,ncol(edge_samples)-3, by = 4)]
colnames(edge_samples) <- counts_ls_skewed$dyad_ids

elephants <- sort(unique(c(counts_df$id_1, counts_df$id_2)))
skewed_adjarr <- array(NA, dim = c(length(elephants), length(elephants), 1000),
                       dimnames = list(elephants, elephants, NULL))
N <- nrow(counts_df)

for (i in 1:N) {
  dyad_row <- counts_df[i, ]
  skewed_adjarr[dyad_row$id_1, dyad_row$id_2, ] <- edge_samples[, i]
}
skewed_adjarr[,,1]

# make matrix to save eigenvector values to
eigen_skewed <- matrix(NA, nrow = 1000, ncol = length(elephants),
                       dimnames = list(1:1000, elephants))

# calculate eigenvector scores
for(i in 1:nrow(eigen_skewed)){
  skewed_net <- igraph::graph_from_adjacency_matrix(skewed_adjarr[,,i], weighted = T, mode = 'undirected', diag = FALSE)
  eigen <- as.data.frame(igraph::eigen_centrality(skewed_net)$vector)
  eigen_skewed[i,] <- eigen$`igraph::eigen_centrality(skewed_net)$vector`
}
eigen_skewed[,1]
eigen_skewed[1,]

eigen_skewed <- eigen_skewed %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'eigen') %>% 
  group_by(id) %>% 
  mutate(mean = mean(eigen)) %>% 
  left_join(eigen_sri[,c('id', 'together0', 'count')], by = 'id')

means <- eigen_skewed %>% 
  select(id, mean, together0, count) %>% 
  distinct()

(figure2l.1 <- ggplot()+
    geom_point(data = eigen_skewed, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = means, aes(x = together0, y = mean),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2l_eigentogether0_skewed_allpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2l.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure2l.2 <- ggplot()+
    geom_point(data = eigen_skewed, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = means, aes(x = together0, y = mean),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(data = eigen_skewed, aes(x = together0, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2l_eigentogether0_skewed_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2l.2, device = 'png', width = 1400, height = 800, units = 'px')

(figure2l.3 <- ggplot(means)+
    geom_point(aes(x = together0, y = mean, colour = count),
               size = 0.5, shape = 19)+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))+
    labs(colour = 'number of\nsightings')
)
ggsave(filename = 'figure2l_eigentogether0_skewed_meanpoints_noline_coloursightings.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2l.3, device = 'png', width = 1400, height = 800, units = 'px')

(figure2l.4 <- ggplot(means)+
    geom_point(aes(x = together0, y = mean, colour = count),
               size = 0.5, shape = 19)+
    geom_smooth(aes(x = together0, y = mean),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))+
    labs(colour = 'number of\nsightings')
)
ggsave(filename = 'figure2l_eigentogether0_skewed_meanpoints_withline_colourssightings.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2l.4, device = 'png', width = 1400, height = 800, units = 'px')

## panel o: right-skewed, eigenvector centrality vs individual sighting count ####
(figure2o.1 <- ggplot()+
    geom_point(data = eigen_skewed, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = means, aes(x = count, y = mean),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2o_eigensightings_skewed_allpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2o.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure2o.2 <- ggplot()+
    geom_point(data = eigen_skewed, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = means, aes(x = count, y = mean),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(data = eigen_skewed, aes(x = count, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2o_eigensightings_skewed_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2o.2, device = 'png', width = 1400, height = 800, units = 'px')

(figure2o.3 <- ggplot(data = means, aes(y = mean,
                                        #y = standardised,
                                        x = count))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean eigenvector centrality',
                       limits = c(0,1),
                       expand = c(0.02,0.02))
)
ggsave(filename = 'figure2o_eigensightings_skewed_meanpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2o.3, device = 'png', width = 1400, height = 800, units = 'px')

(figure2o.4 <- ggplot(data = means, aes(x = count, y = mean))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean eigenvector centrality',
                       limits = c(0,1),
                       expand = c(0,0))
)
ggsave(filename = 'figure2o_eigensightings_skewed_meanpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2o.4, device = 'png', width = 1400, height = 800, units = 'px')

# clean up
rm(list = ls()[!ls() %in% c('figure2a', 'figure2b', 'figure2c',
                            'figure2d', 'figure2e', 'figure2f',
                            'figure2g.1','figure2g.2','figure2g.3',
                            'figure2h.1','figure2h.2',
                            'figure2i.1','figure2i.2',
                            'figure2j.1','figure2j.2',
                            'figure2k.1','figure2k.2',
                            'figure2l.1','figure2l.2',
                            'figure2m.1','figure2m.2',
                            'figure2n.1','figure2n.2','figure2n.3','figure2n.4',
                            'figure2o.1','figure2o.2','figure2o.3','figure2o.4')]) ; gc()

# save for plotting
save.image('../outputs/sparse_network_methods_figures/figure2_plots.RData')

## merge second six panels ####
# run patchwork for plots with no lines
(figure2g.1 + figure2h.1 + figure2i.1) / (figure2j.1 + figure2k.1 + figure2l.1) +
  plot_annotation(tag_levels = list(c('g','h','i','j','k','l')))
ggsave(filename = 'figure2_bottomhalf_nolines.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2100, height = 1200, units = 'px')

# run patchwork for plots with lines
(figure2g.2 + figure2h.2 + figure2i.2) / (figure2j.2 + figure2k.2 + figure2l.2) +
  plot_annotation(tag_levels = list(c('g','h','i','j','k','l')))
ggsave(filename = 'figure2_bottomhalf_withlines.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2100, height = 1200, units = 'px')

## merge final 3 panels ####
# run patchwork for plots with no lines
(figure2m.1 + figure2n.1 + figure2o.1) +
  plot_annotation(tag_levels = list(c('m','n','o')))
ggsave(filename = 'figure2_eigensightings_nolines.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2100, height = 1200, units = 'px')

# run patchwork for plots with lines
(figure2m.2 + figure2n.2 + figure2o.2) +
  plot_annotation(tag_levels = list(c('m','n','o')))
ggsave(filename = 'figure2_eigensightings_withlines.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2100, height = 1200, units = 'px')

## merge all fifteen panels ####
# load('../outputs/sparse_network_methods_figures/figure2_plots.RData')
# run patchwork for plots with no lines
(figure2a + figure2b + figure2c) / (figure2d + figure2e + figure2f) / (figure2g.1 + figure2h.1 + figure2i.1) / (figure2j.1 + figure2k.1 + figure2l.1) / (figure2m.1 + figure2n.1 + figure2o.1) +
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'figure2_nolines.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2100, height = 2400, units = 'px')

# run patchwork for plots with lines
(figure2a + figure2b + figure2c) / (figure2d + figure2e + figure2f) / (figure2g.2 + figure2h.2 + figure2i.2) / (figure2j.2 + figure2k.2 + figure2l.2) / (figure2m.2 + figure2n.2 + figure2o.2)  +
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'figure2_withlines.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2100, height = 2400, units = 'px')

# clean up
rm(list = ls()) ; gc()

#### Figure 3: 4 panel plot of conditional BISoN for prior/posterior/mean_edge_vs_sightings/eigenvector_vs_sightings. if 0: beta(0.7,10), if 1: beta(1,5) ####
## panel a: prior ####
# define sequence over which to plot
x <- seq(0, 1, length = 100)

# calculate probability of all values of x assuming completely flat prior
conditional1 <- dbeta(x = x, shape1 = 0.7, shape2 = 10)
conditional2 <- dbeta(x = x, shape1 = 1, shape2 = 5)
data <- data.frame(x = x,
                   density1 = conditional1,
                   density2 = conditional2)
# plot
(figure3a <- ggplot(data)+
    geom_line(aes(x = x, y = density1), linewidth = 1.2, colour = rgb(33/255, 145/255, 140/255))+
    geom_line(aes(x = x, y = density2), linewidth = 1.2, colour = rgb(68/255,1/255,84/255,1))+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density', expand = c(0,0))+
    coord_cartesian(ylim = c(0,15))
)
ggsave(filename = 'figure3a_prior_conditional.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3a, device = 'png', width = 1400, height = 800, units = 'px')

# clean up
rm(list = ls()[!ls() %in% c('figure3a')]) ; gc()

## panel b: posterior ####
load('motnp_edgeweights_conditionalprior.RData')
counts <- counts_df %>% 
  select(dyad_id, count_dyad, event_count) %>% 
  rename(dyad = dyad_id)
edges <- edges %>% 
  filter(chain == 'chain1') %>% 
  left_join(counts, by = 'dyad')
edges$together <- ifelse(edges$event_count == 0, 'never together', 'together at least once')
rm(edge_binary, edge_samples, edgelist, motnp_ages, nodes, summary, x, n_dyads, n_samples, make_edgelist, plot_network_threshold, i) ; gc()

# plot
(figure3b <- ggplot(data = edges)+
    geom_density(aes(group = dyad, x = edge_draw, colour = together), show.legend = FALSE)+
    stat_density(aes(group = dyad, x = edge_draw, colour = together),
                 geom = "line", position = "identity", linewidth = 0)+
    scale_x_continuous(name = 'edge weight', limits = c(-0.05, 1.05), expand = c(0,0))+
    scale_y_continuous(name = 'density', limits = c(-2, 72), expand = c(0,0))+
    scale_colour_manual(values = c(rgb(33/255, 145/255, 140/255, 0.05), rgb(68/255, 1/255, 84/255, 0.05)),
                        aesthetics = 'colour')+
    #scale_alpha(0.05)+
    theme(legend.position = c(0.8,0.7),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8))+
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))+
    labs(colour = 'sightings together')
)
ggsave(filename = 'figure3b_posterior_default.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3b, device = 'png', width = 1400, height = 800, units = 'px')

# override legend size for big plot
(figure3b <- ggplot(data = edges)+
    geom_density(aes(group = dyad, x = edge_draw, colour = together), show.legend = FALSE)+
    stat_density(aes(group = dyad, x = edge_draw, colour = together),
                 geom = "line", position = "identity", linewidth = 0)+
    scale_x_continuous(name = 'edge weight', limits = c(-0.05, 1.05), expand = c(0,0))+
    scale_y_continuous(name = 'density', limits = c(-2, 72), expand = c(0,0))+
    scale_colour_manual(values = c(rgb(33/255, 145/255, 140/255, 0.05), rgb(68/255, 1/255, 84/255, 0.05)),
                        aesthetics = 'colour')+
    #scale_alpha(0.05)+
    theme(legend.position = c(0.6,0.7),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          legend.key.height = unit(2.5, 'mm'),
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 4))+
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))+
    labs(colour = 'sightings together')
)

#set.seed(15) ; plot_dyads <- sample(1:nrow(counts), size = 1000, replace = F)
#plot_dyad_ids <- unique(edges$dyad)[plot_dyads]

#edges_subset <- edges[edges$dyad %in% plot_dyad_ids,]
#edges_subset <- edges_subset[edges_subset$chain == 'chain1',]
#(figure3b <- ggplot(data = edges_subset)+
#    geom_density(aes(group = dyad, x = edge_draw, colour = together), show.legend = FALSE)+
#    stat_density(aes(group = dyad, x = edge_draw, colour = together),
#                 geom = "line", position = "identity", linewidth = 0)+
#    scale_x_continuous(name = 'edge weight', limits = c(-0.05, 1.05), expand = c(0,0))+
#    scale_y_continuous(name = 'density', limits = c(-2, 72), expand = c(0,0))+
#    scale_colour_manual(values = c(rgb(33/255, 145/255, 140/255, 0.05), rgb(68/255, 1/255, 84/255, 0.05)),
#                        aesthetics = 'colour')+
#    theme(legend.position = c(0.8,0.7),
#          legend.background = element_rect(fill = 'white', colour = 'black'),
#          legend.key.height = unit(4, 'mm'),
#          legend.title = element_text(size = 10),
#          legend.text = element_text(size = 8))+
#    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))+
#    labs(colour = 'sightings together')
#)

## panel c: mean edge weight vs dyad sightings ####
means <- edges %>% 
  group_by(dyad) %>% 
  mutate(mean = mean(edge_draw)) %>% 
  select(dyad, mean) %>% 
  distinct() %>% 
  rename(dyad_id = dyad) %>% 
  left_join(counts_df[,c('dyad_id','node_1_males','node_2_males','event_count','count_dyad','count_1','count_2')],
            by = 'dyad_id')

(figure3c.1 <- ggplot(means, aes(x = count_dyad, y = mean))+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mean weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure3c_edgesightings_default_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3c.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure3c.2 <- ggplot(means, aes(x = count_dyad, y = mean))+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mean weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure3c_edgesightings_default_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3c.2, device = 'png', width = 1400, height = 800, units = 'px')

(figure3c.4 <- ggplot(means, aes(x = count_dyad, y = mean))+
    geom_point(colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mean weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure3c_edgesightings_default_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3c.4, device = 'png', width = 1400, height = 800, units = 'px')

## panel d: eigenvector vs total number of dyads where together = 0 ####
edge_samples <- fit_edges_motnp$draws() %>% as.data.frame()
edge_samples <- edge_samples[,(n_chains+1):ncol(edge_samples)]
edge_samples <- edge_samples[,seq(1,ncol(edge_samples)-3, by = 4)]
colnames(edge_samples) <- counts_df$dyad_id

elephants <- sort(unique(c(counts_df$id_1, counts_df$id_2)))
n_eles <- length(unique(c(counts_df$id_1, counts_df$id_2)))
motnp_adjarr <- array(NA, dim = c(length(elephants), length(elephants), 1000),
                      dimnames = list(elephants, elephants, NULL))
N <- nrow(counts_df)

for (i in 1:N) {
  dyad_row <- counts_df[i, ]
  motnp_adjarr[dyad_row$id_1, dyad_row$id_2, ] <- edge_samples[, i]
  motnp_adjarr[dyad_row$id_2, dyad_row$id_1, ] <- edge_samples[, i]
}
motnp_adjarr[,,1]

# make matrix to save eigenvector values to
eigen_motnp <- matrix(NA, nrow = 1000, ncol = length(elephants),
                      dimnames = list(1:1000, elephants))

# calculate eigenvector scores
for(i in 1:nrow(eigen_motnp)){
  motnp_net <- igraph::graph_from_adjacency_matrix(adjmatrix = motnp_adjarr[,,i], diag = FALSE, mode = 'undirected', weighted = T)
  eigen <- as.data.frame(igraph::eigen_centrality(motnp_net)$vector)
  eigen_motnp[i,] <- eigen$`igraph::eigen_centrality(motnp_net)$vector`
}
eigen_motnp[,1]
eigen_motnp[1,]

eigen_motnp <- eigen_motnp %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'eigen') %>% 
  group_by(id) %>% 
  mutate(mean = mean(eigen))

counts <- data.frame(id = elephants, count = NA, together0 = NA)
together0 <- counts_df %>% filter(event_count == 0)
for(i in 1:nrow(counts)){
  counts$count[i] <- ifelse(i == nrow(counts),
                            unique(counts_df$count_2[which(counts_df$id_2 == counts$id[i])]),
                            unique(counts_df$count_1[which(counts_df$id_1 == counts$id[i])]))
  x <- together0 %>%
    filter(id_1 == counts$id[i] | id_2 == counts$id[i])
  counts$together0[i] <- nrow(x)
}

eigen_motnp <- eigen_motnp %>% 
  left_join(counts, by = 'id')

means <- eigen_motnp %>% 
  select(id, mean, count, together0) %>% 
  distinct()

(figure3d.1 <- ggplot()+
    geom_point(data = eigen_motnp, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01))+
    geom_point(data = means, aes(x = together0, y = mean),
               colour = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure3d_eigentogether_motnp_allpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3d.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure3d.2 <- ggplot()+
    geom_point(data = eigen_motnp, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01))+
    geom_point(data = means, aes(x = together0, y = mean),
               colour = rgb(33/255, 145/255, 140/255))+
    geom_smooth(data = eigen_motnp, aes(x = together0, y = eigen),
                colour = rgb(68/255,1/255,84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure3d_eigentogether_motnp_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3d.2, device = 'png', width = 1400, height = 800, units = 'px')

(figure3d.3 <- ggplot(means)+
    geom_point(aes(x = together0, y = mean, colour = count),
               size = 0.5, shape = 19)+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))+
    labs(colour = 'number of\nsightings')+
    theme(legend.position = c(0.5,0.18),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          legend.direction = 'horizontal',
          legend.key.height = unit(2.5, 'mm'),
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 4))
)
ggsave(filename = 'figure3d_eigentogether_motnp_meanpoints_noline_coloursightings.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3d.3, device = 'png', width = 1400, height = 800, units = 'px')

(figure3d.4 <- ggplot(means)+
    geom_point(aes(x = together0, y = mean, colour = count),
               size = 0.5, shape = 19)+
    geom_smooth(aes(x = together0, y = mean),
                colour = rgb(68/255,1/255,84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))+
    labs(colour = 'number of\nsightings')+
    theme(legend.position = c(0.5,0.2),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          legend.direction = 'horizontal',
          legend.key.height = unit(2.5, 'mm'),
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 4))
)
ggsave(filename = 'figure3d_eigentogether_motnp_meanpoints_withline_coloursightings.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3d.4, device = 'png', width = 1400, height = 800, units = 'px')

## panel e: eigenvector vs individual sightings ####
(figure3e.1 <- ggplot()+
    geom_point(data = eigen_motnp, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01))+
    geom_point(data = means, aes(x = count, y = mean),
               colour = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure3e_eigensightings_motnp_allpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3e.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure3e.2 <- ggplot()+
    geom_point(data = eigen_motnp, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01))+
    geom_point(data = means, aes(x = count, y = mean),
               colour = rgb(33/255, 145/255, 140/255))+
    geom_smooth(data = eigen_motnp, aes(x = count, y = eigen),
                colour = rgb(68/255,1/255,84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure3e_eigensightings_motnp_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3e.2, device = 'png', width = 1400, height = 800, units = 'px')

(figure3e.3 <- ggplot(data = means, aes(y = mean,
                                        x = count))+
    geom_point(colour = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean centrality',
                       limits = c(0, 1),
                       expand = c(0.02,0.02))
)
ggsave(filename = 'figure3e_eigensightings_motnp_meanpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3e.3, device = 'png', width = 1400, height = 800, units = 'px')

(figure3e.4 <- ggplot(data = means, aes(x = count, y = mean))+
    geom_point(colour = rgb(33/255, 145/255, 140/255))+
    geom_smooth(colour = rgb(68/255,1/255,84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean centrality',
                       limits = c(0,1),
                       expand = c(0.02,0.02))
)
ggsave(filename = 'figure3e_eigensightings_motnp_meanpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3e.4, device = 'png', width = 1400, height = 800, units = 'px')

# save image for future
save.image('../outputs/sparse_network_methods_figures/figure3.RData')

## merge panels ####
#load('../outputs/sparse_network_methods_figures/figure3.RData')

# run patchwork and save
(figure3a + figure3b) / (figure3c.2 + figure3d.2) +
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'figure3_withlines_4panels.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2100, height = 1200, units = 'px')

(figure3a + figure3b) / (figure3c.1 + figure3d.1) +
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'figure3_nolines_4panels.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2100, height = 1200, units = 'px')

(figure3a + figure3b + figure3c.2) / (figure3e.2 + figure3d.4) +
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'figure3_withlines_5panels.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2100, height = 1200, units = 'px')

(figure3a + figure3b + figure3c.1) / (figure3e.1 + figure3d.3) +
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'figure3_nolines_5panels.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2100, height = 1200, units = 'px')

# clean up
rm(list = ls()) ; gc()

