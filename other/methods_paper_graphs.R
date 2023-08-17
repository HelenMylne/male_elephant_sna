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

################## NEW GRAPHS #################
# define sequence over which to plot
x <- seq(0, 1, length = 100)

#### Figure 1: 6 facet plot of 10 sightings vs 2 sightings, never together vs half together vs always together ####
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
    scale_y_continuous(name = 'probability density', limits = c(0,1))+
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
    scale_y_continuous(name = 'probability density', limits = c(0,1.05), expand = c(0,0))
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
    scale_y_continuous(name = 'probability density',
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
    scale_y_continuous(name = 'probability density',
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
    scale_y_continuous(name = 'probability density')
  )

#set.seed(15) ; plot_dyads <- sample(1:nrow(counts_df), size = 1000, replace = F)
#plot_dyad_ids <- unique(edge_samples$parameter)[plot_dyads]
#(figure2e <- ggplot(data = edge_samples[edge_samples$parameter %in% plot_dyad_ids,])+
#    geom_density(aes(group = parameter, x = weight), colour = rgb(33/255, 145/255, 140/255, 0.1))+
#    scale_x_continuous(name = 'edge weight')+
#    scale_y_continuous(name = 'probability density')
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
    scale_y_continuous(name = 'probability density')
)
#(figure2f <- ggplot(data = edge_samples[edge_samples$parameter %in% plot_dyad_ids,])+
#    geom_density(aes(group = parameter, x = weight), colour = rgb(33/255, 145/255, 140/255, 0.1))+
#    scale_x_continuous(name = 'edge weight')+
#    scale_y_continuous(name = 'probability density')
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
   #geom_smooth(colour = rgb(68/255, 1/255, 84/255, 0.1))+
   scale_x_continuous(name = 'total dyad sightings')+
   scale_y_continuous(name = 'SRI edge weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2g_edgesightings_sri_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2g.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure2g.1 <- ggplot(counts_df, aes(x = count_dyad, y = sri))+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'SRI edge weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2g_edgesightings_sri_withline_changealpha.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2g.1, device = 'png', width = 1400, height = 800, units = 'px')

counts_df$sri_rank <- as.integer(as.factor(counts_df$sri))
p7 <- as.data.frame(table(counts_df$count_dyad, counts_df$sri)) %>% 
  mutate(count_dyad = as.numeric(Var1),
         sri_rank = as.numeric(Var2)) %>% 
  filter(Freq > 0) %>% 
  left_join(distinct(counts_df[,c('sri','sri_rank')]), by = 'sri_rank')
(figure2g.2 <- ggplot()+
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
    scale_y_continuous(name = 'SRI edge weight')
)
ggsave(filename = 'figure2g_edgesightings_sri_withline_changesize.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2g.2, device = 'png', width = 1400, height = 800, units = 'px')

(figure2g.2 <- ggplot()+
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
    scale_y_continuous(name = 'SRI edge weight')
)

# clean up 
rm(list = ls()[!ls() %in% c('figure2a', 'figure2b', 'figure2c',
                            'figure2d', 'figure2e', 'figure2f',
                            'figure2g.1','figure2g.2')]) ; gc()

## panel h: default, mean edge weight vs dyad sighting count ####
load('../outputs/sparse_network_methods_figures/model_run_default.RData')
means <- edge_samples %>% 
  group_by(parameter) %>% 
  mutate(mean = mean(weight)) %>% 
  select(parameter, mean) %>% 
  distinct()

means$dyad_males <- as.integer(as.factor(means$parameter))

means <- means %>% 
  left_join(counts_df[,c('dyad_id','node_1_males','node_2_males','event_count','count_dyad','count_1','count_2','dyad_males')],
            by = 'dyad_males')
rm(edge_samples, fit_edges_motnp) ; gc()

(figure2h.1 <- ggplot(means, aes(x = count_dyad, y = mean))+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mean edge weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2h_edgesightings_default_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2h.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure2h.2 <- ggplot(means, aes(x = count_dyad, y = mean))+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mean edge weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2h_edgesightings_default_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2h.2, device = 'png', width = 1400, height = 800, units = 'px')

# clean up
rm(list = ls()[!ls() %in% c('figure2a', 'figure2b', 'figure2c',
                            'figure2d', 'figure2e', 'figure2f',
                            'figure2g.1','figure2g.2',
                            'figure2h.1','figure2h.2')]) ; gc()

## panel i: right-skewed, mean edge weight vs dyad sighting count ####
load('../outputs/sparse_network_methods_figures/model_run_skewed.RData')
means <- edge_samples %>% 
  group_by(parameter) %>% 
  mutate(mean = mean(weight)) %>% 
  select(parameter, mean) %>% 
  distinct()

means$dyad_males <- as.integer(as.factor(means$parameter))

means <- means %>% 
  left_join(counts_df[,c('dyad_id','node_1_males','node_2_males','event_count','count_dyad','count_1','count_2','dyad_males')],
            by = 'dyad_males', multiple = 'all')
rm(edge_samples, fit_edges_skewed, edge_binary_skewed, edge_binary_generic, figure2d, figure2e, figure2f) ; gc()

(figure2i.1 <- ggplot(means, aes(x = count_dyad, y = mean))+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mean edge weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2i_edgesightings_rightskewed_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2i.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure2i.2 <- ggplot(means, aes(x = count_dyad, y = mean))+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mean edge weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2i_edgesightings_rightskewed_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2i.2, device = 'png', width = 1400, height = 800, units = 'px')

# clean up
rm(list = ls()[!ls() %in% c('counts_df',
                            'figure2a', 'figure2b', 'figure2c',
                            'figure2d', 'figure2e', 'figure2f',
                            'figure2g.1','figure2g.2',
                            'figure2h.1','figure2h.2',
                            'figure2i.1','figure2i.2')]) ; gc()

## panel j0: SRI, eigenvector centrality vs individual sighting count ####
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

counts <- data.frame(id = elephants, count = NA)
for(i in 1:nrow(counts)){
  counts$count[i] <- ifelse(i == nrow(counts),
                            unique(counts_df$count_2[which(counts_df$id_2 == counts$id[i])]),
                            unique(counts_df$count_1[which(counts_df$id_1 == counts$id[i])]))
}

eigen_sri <- eigen_sri %>% 
  left_join(counts, by = 'id')

(figure2j.1 <- ggplot(eigen_sri, aes(x = count, y = eigen))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2j_eigensightings_sri_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2j.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure2j.2 <- ggplot(eigen_sri, aes(x = count, y = eigen))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2j_eigensightings_sri_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2j.2, device = 'png', width = 1400, height = 800, units = 'px')

# clean up
rm(list = ls()[!ls() %in% c('figure2a', 'figure2b', 'figure2c',
                            'figure2d', 'figure2e', 'figure2f',
                            'figure2g.1','figure2g.2',
                            'figure2h.1','figure2h.2',
                            'figure2i.1','figure2i.2',
                            'figure2j.1','figure2j.2')]) ; gc()

## panel k: default, eigenvector centrality vs individual sighting count ####
load('../outputs/sparse_network_methods_figures/model_run_default.RData')
rm(counts_ls, edge_binary_generic, edge_binary_skewed, figure2d, figure2e, make_edgelist, plot_network_threshold) ; gc()

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
eigen_default[,1]
eigen_default[1,]

eigen_default <- eigen_default %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'eigen') %>% 
  group_by(id) %>% 
  mutate(mean = mean(eigen))

counts <- data.frame(id = elephants, count = NA)
for(i in 1:nrow(counts)){
  counts$count[i] <- ifelse(i == nrow(counts),
                            unique(counts_df$count_2[which(counts_df$id_2 == counts$id[i])]),
                            unique(counts_df$count_1[which(counts_df$id_1 == counts$id[i])]))
}

eigen_default <- eigen_default %>% 
  left_join(counts, by = 'id')

means <- eigen_default %>% 
  select(id, mean, count) %>% 
  distinct()

(figure2k.1 <- ggplot()+
    geom_point(data = eigen_default, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01))+
    geom_point(data = means, aes(x = count, y = mean),
               colour = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
  )
ggsave(filename = 'figure2k_eigensightings_default_allpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2k.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure2k.2 <- ggplot()+
    geom_point(data = eigen_default, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01))+
    geom_point(data = means, aes(x = count, y = mean),
               colour = rgb(33/255, 145/255, 140/255))+
    geom_smooth(data = eigen_default, aes(x = count, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2k_eigensightings_default_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2k.2, device = 'png', width = 1400, height = 800, units = 'px')

(figure2k.3 <- ggplot(data = means, aes(y = mean,
                                        x = count))+
    geom_point(colour = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean eigenvector centrality',
                       limits = c(0, 1),
                       expand = c(0.02,0.02))
)
ggsave(filename = 'figure2k_eigensightings_default_meanpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2k.3, device = 'png', width = 1400, height = 800, units = 'px')

(figure2k.4 <- ggplot(data = means, aes(x = count, y = mean))+
    geom_point(colour = rgb(33/255, 145/255, 140/255))+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean eigenvector centrality',
                       limits = c(0,1),
                       expand = c(0.02,0.02))
)
ggsave(filename = 'figure2k_eigensightings_default_meanpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2k.4, device = 'png', width = 1400, height = 800, units = 'px')

# clean up
rm(list = ls()[!ls() %in% c('figure2a', 'figure2b', 'figure2c',
                            'figure2d', 'figure2e', 'figure2f',
                            'figure2g.1','figure2g.2',
                            'figure2h.1','figure2h.2',
                            'figure2i.1','figure2i.2',
                            'figure2j.1','figure2j.2',
                            'figure2k.1', 'figure2k.2', 'figure2k.3', 'figure2k.4')]) ; gc()

## panel l: right-skewed, eigenvector centrality vs individual sighting count ####
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
  mutate(mean = mean(eigen))

counts <- data.frame(id = elephants, count = NA)
for(i in 1:nrow(counts)){
  counts$count[i] <- ifelse(i == nrow(counts),
                            unique(counts_df$count_2[which(counts_df$id_2 == counts$id[i])]),
                            unique(counts_df$count_1[which(counts_df$id_1 == counts$id[i])]))
}

eigen_skewed <- eigen_skewed %>% 
  left_join(counts, by = 'id')

means <- eigen_skewed %>% 
  select(id, mean, count) %>% 
  distinct()# %>% 
  #mutate(standardised = mean / max(mean))

(figure2l.1 <- ggplot()+
    geom_point(data = eigen_skewed, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01))+
    geom_point(data = means, aes(x = count, y = mean),
               colour = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2l_eigensightings_skewed_allpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2l.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure2l.2 <- ggplot()+
    geom_point(data = eigen_skewed, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01))+
    geom_point(data = means, aes(x = count, y = mean),
               colour = rgb(33/255, 145/255, 140/255))+
    geom_smooth(data = eigen_skewed, aes(x = count, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure2l_eigensightings_skewed_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2l.2, device = 'png', width = 1400, height = 800, units = 'px')

(figure2l.3 <- ggplot(data = means, aes(y = mean,
                                        #y = standardised,
                                        x = count))+
    geom_point(colour = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean eigenvector centrality',
                       limits = c(0,1),
                       expand = c(0.02,0.02))
)
ggsave(filename = 'figure2l_eigensightings_skewed_meanpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2l.3, device = 'png', width = 1400, height = 800, units = 'px')

(figure2l.4 <- ggplot(data = means, aes(x = count, y = mean))+
    geom_point(colour = rgb(33/255, 145/255, 140/255))+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean eigenvector centrality',
                       limits = c(0,1),
                       expand = c(0,0))
)
ggsave(filename = 'figure2l_eigensightings_skewed_meanpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure2l.4, device = 'png', width = 1400, height = 800, units = 'px')

# clean up
rm(list = ls()[!ls() %in% c('figure2a', 'figure2b', 'figure2c',
                            'figure2d', 'figure2e', 'figure2f',
                            'figure2g.1','figure2g.2',
                            'figure2h.1','figure2h.2',
                            'figure2i.1','figure2i.2',
                            'figure2j.1','figure2j.2',
                            'figure2k.1', 'figure2k.2', 'figure2k.3', 'figure2k.4',
                            'figure2l.1', 'figure2l.2', 'figure2l.3', 'figure2l.4')]) ; gc()

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

## merge all twelve panels ####
# load('../outputs/sparse_network_methods_figures/figure2_plots.RData')
# run patchwork for plots with no lines
(figure2a + figure2b + figure2c) / (figure2d + figure2e + figure2f) / (figure2g.1 + figure2h.1 + figure2i.1) / (figure2j.1 + figure2k.1 + figure2l.1) +
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'figure2_nolines.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2100, height = 2400, units = 'px')

# run patchwork for plots with lines
(figure2a + figure2b + figure2c) / (figure2d + figure2e + figure2f) / (figure2g.1 + figure2h.2 + figure2i.2) / (figure2j.2 + figure2k.2 + figure2l.2) +
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'figure2_withlines.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2100, height = 2400, units = 'px')

# clean up
rm(list = ls()) ; gc()

#### Figure 3: 4 panel plot of conditional BISoN for prior/posterior/mean_edge_vs_sightings/eigenvector_vs_sightings. if 0: beta(0.7,10), if 1: beta(1,5) ####
## panel 1: prior ####
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
    scale_y_continuous(name = 'probability density', expand = c(0,0))+
    coord_cartesian(ylim = c(0,15))
)
ggsave(filename = 'figure3a_prior_conditional.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3a, device = 'png', width = 1400, height = 800, units = 'px')

# clean up
rm(list = ls()[!ls() %in% c('figure3a')]) ; gc()

## panel 2: posterior ####
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
    scale_y_continuous(name = 'probability density', limits = c(-2, 72), expand = c(0,0))+
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
    scale_y_continuous(name = 'probability density', limits = c(-2, 72), expand = c(0,0))+
    scale_colour_manual(values = c(rgb(33/255, 145/255, 140/255, 0.05), rgb(68/255, 1/255, 84/255, 0.05)),
                        aesthetics = 'colour')+
    #scale_alpha(0.05)+
    theme(legend.position = c(0.7,0.7),
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
#    scale_y_continuous(name = 'probability density', limits = c(-2, 72), expand = c(0,0))+
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

## panel 3: mean edge weight vs dyad sightings ####
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
    scale_y_continuous(name = 'mean edge weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure3c_edgesightings_default_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3c.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure3c.2 <- ggplot(means, aes(x = count_dyad, y = mean))+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mean edge weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure3c_edgesightings_default_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3c.2, device = 'png', width = 1400, height = 800, units = 'px')

(figure3c.4 <- ggplot(means, aes(x = count_dyad, y = mean))+
    geom_point(colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mean edge weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure3c_edgesightings_default_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3c.4, device = 'png', width = 1400, height = 800, units = 'px')

## panel 4: eigenvector vs individual sightings ####
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

counts <- data.frame(id = elephants, count = NA)
for(i in 1:nrow(counts)){
  counts$count[i] <- ifelse(i == nrow(counts),
                            unique(counts_df$count_2[which(counts_df$id_2 == counts$id[i])]),
                            unique(counts_df$count_1[which(counts_df$id_1 == counts$id[i])]))
}

eigen_motnp <- eigen_motnp %>% 
  left_join(counts, by = 'id')

means <- eigen_motnp %>% 
  select(id, mean, count) %>% 
  distinct()

(figure3d.1 <- ggplot()+
    geom_point(data = eigen_motnp, aes(x = count, y = eigen),
               colour = rgb(33/255, 145/255, 140/255, 0.01))+
    geom_point(data = means, aes(x = count, y = mean),
               colour = rgb(253/255, 231/255, 37/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure3d_eigensightings_motnp_allpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3d.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure3d.2 <- ggplot()+
    geom_point(data = eigen_motnp, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01))+
    geom_point(data = means, aes(x = count, y = mean),
               colour = rgb(33/255, 145/255, 140/255))+
    geom_smooth(data = eigen_motnp, aes(x = count, y = eigen),
                colour = rgb(68/255,1/255,84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'figure3d_eigensightings_motnp_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3d.2, device = 'png', width = 1400, height = 800, units = 'px')

(figure3d.3 <- ggplot(data = means, aes(y = mean,
                                        x = count))+
    geom_point(colour = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean centrality',
                       limits = c(0, 1),
                       expand = c(0.02,0.02))
)
ggsave(filename = 'figure3d_eigensightings_motnp_meanpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3d.3, device = 'png', width = 1400, height = 800, units = 'px')

(figure3d.4 <- ggplot(data = means, aes(x = count, y = mean))+
    geom_point(colour = rgb(33/255, 145/255, 140/255))+
    geom_smooth(colour = rgb(68/255,1/255,84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean centrality',
                       limits = c(0,1),
                       expand = c(0.02,0.02))
)
ggsave(filename = 'figure3d_eigensightings_motnp_meanpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure3d.4, device = 'png', width = 1400, height = 800, units = 'px')

# save image for future
save.image('../outputs/sparse_network_methods_figures/figure3.RData')

## merge panels ####
# run patchwork and save
(figure3a + figure3b) / (figure3c.2 + figure3d.2) +
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'figure3_withlines.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2100, height = 1200, units = 'px')

(figure3a + figure3b) / (figure3c.1 + figure3d.1) +
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'figure3_nolines.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2100, height = 1200, units = 'px')

# clean up
rm(list = ls()) ; gc()

############### ORIGINAL GRAPHS ###############
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
  points(mean(x[4,], na.rm = T) ~ x[3,1], pch = 19, col = 'yellow')
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

