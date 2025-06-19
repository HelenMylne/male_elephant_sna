######## information ########
# run BISoN models using simulated data

######## set up ########
# library(cmdstanr) ; library(tidyverse) ; library(LaplacesDemon) ; library(patchwork) ; library(igraph)
library(cmdstanr, lib.loc = '../packages/')
library(tidyverse, lib.loc = '../packages/')
library(LaplacesDemon, lib.loc = '../packages/')
library(patchwork, lib.loc = '../packages/')
library(igraph, lib.loc = '../packages/')

## set stan parameters
set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0/')
n_chains <- 4

## define plot theme
theme_set(theme_bw(base_size = 12))
colours <- c('#fde725','#5ec962','#21918c','#3b528b','#440154')

######## load data ########
load('methods_paper/simulation_outputs/create_simulated_data.RData')

######## BISoN unconditional: uniform ########
#### run model ####
## compile edge model
edge_binary_uniform <- cmdstan_model('methods_paper/models/edge_binary_uniform.stan')

## create data list for uniform model
counts_ls_uniform <- list(n_dyads = nrow(sim_dyads),
                          dyad_ids = sim_dyads$dyad_id,
                          together = sim_dyads$together,
                          count_dyad = sim_dyads$total_sightings)

## fit model
fit_edges_uniform <- edge_binary_uniform$sample(
  data = counts_ls_uniform,
  chains = n_chains,
  parallel_chains = n_chains)
save.image('methods_paper/simulation_outputs/sim_uniform_run.RData')

## extract edges
edges <- fit_edges_uniform$draws() %>% as.data.frame()
edges <- edges[,(n_chains+1):ncol(edges)]        # remove lp__ columns
edges1 <- edges[,seq(1,ncol(edges)-3, by = 4)]   # select only chain 1
edge_names <- data.frame(name = colnames(edges1)) %>%
  separate(name, into = c('chain','weight'), sep = '.edge_')
colnames(edges1) <- edge_names$weight
edges2 <- edges[,seq(2,ncol(edges)-2, by = 4)] ; colnames(edges2) <- edge_names$weight  # select only chain 2
edges3 <- edges[,seq(3,ncol(edges)-1, by = 4)] ; colnames(edges3) <- edge_names$weight   # select only chain 3
edges4 <- edges[,seq(4,ncol(edges)-0, by = 4)] ; colnames(edges4) <- edge_names$weight   # select only chain 4
edges <- rbind(edges1, edges2, edges3, edges4)
n_samples <- nrow(edges)
edge_samples1 <- edges1 %>%
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'weight')
rm(edges1, edges2, edges3, edges4) ; gc()

#### plot edges ####
## plot all edges
(figure_uniform_posterior <- ggplot(data = edge_samples1)+
   geom_density(aes(group = parameter, x = weight),
                colour = colours[3], alpha = 0.1)+
   scale_x_continuous(name = 'edge weight')+
   scale_y_continuous(name = 'density')
)
ggsave(filename = 'sim_posterior_uniform_alldyads.png',
       path = 'methods_paper/simulation_outputs/',
       plot = figure_uniform_posterior, device = 'png', width = 1600, height = 700, units = 'px')

## sample edges
edge_samples1$dyad_id <- rep(sim_dyads$dyad_id, n_samples/4)
edge_samples1 <- edge_samples1 %>%
  left_join(sim_dyads[,c('dyad_id','node1','node2','sightings1','sightings2','total_sightings','together','sri')],
            by = 'dyad_id') %>%
  mutate(never = ifelse(together == 0, 'never together', 'together at least once'))

set.seed(15) ; plot_dyads <- sample(1:nrow(sim_dyads), size = 100, replace = F)
plot_dyad_ids <- unique(edge_samples1$parameter)[plot_dyads]
edges_subset <- edge_samples1[edge_samples1$parameter %in% plot_dyad_ids,]

## plot sample
(figure_uniform_posterior <- ggplot(data = edges_subset)+
    geom_density(aes(group = parameter, x = weight, colour = never),
                 show.legend = F, linewidth = 0.2)+
    stat_density(aes(group = parameter, x = weight, colour = never),
                 geom = "line", position = "identity", linewidth = 0)+
    scale_x_continuous(name = 'edge weight', limits = c(0,1))+
    scale_y_continuous(name = 'density', limits = c(0,40))+
    scale_colour_manual(values = colours[c(3,5)],
                        aesthetics = 'colour')+
    theme(legend.position = 'bottom',
          legend.background = element_blank(),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 8))+
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))
)
ggsave(filename = 'sim_posterior_uniform_sampledyads.png',
       path = 'methods_paper/simulation_outputs/',
       plot = figure_uniform_posterior, device = 'png',
       width = 1000, height = 1000, units = 'px')

## clean up and save
save.image('methods_paper/simulation_outputs/sim_uniform_run.RData')

#### mean edge weight vs dyad sighting count ####
#load('methods_paper/simulation_outputs/sim_uniform_run.RData')
averages <- data.frame(dyad = colnames(edges),
                       mean = apply(edges, 2, mean),
                       median = apply(edges, 2, median),
                       uniform_lwr = apply(edges, 2, quantile, prob = 0.025),
                       uniform_upr = apply(edges, 2, quantile, prob = 0.975))
averages$dyad_id <- sim_dyads$dyad_id
averages <- averages %>%
  left_join(sim_dyads[,c('dyad_id','node1','node2','together','never','total_sightings','sightings1','sightings2')],
            by = 'dyad_id')

(figure_uniform_edgesightings.2 <- ggplot()+
    geom_point(data = averages,
               aes(x = total_sightings, y = median),
               colour = rgb(33/255, 145/255, 140/255, 0.1),
               size = 0.5,
               shape = 19)+
    geom_smooth(data = averages,
                aes(x = total_sightings, y = median, linetype = as.factor(never)),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'median weight', limits = c(-0.02,1.02), expand = c(0,0))+
    scale_linetype_manual(name = NULL,
                          values = c(1,6))+
    theme(legend.position = 'bottom',
          legend.key.height = unit(4, 'mm'),
          legend.text = element_text(size = 8)))
# ggsave(filename = 'edgesightingsuniform_twolines.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = last_plot(), device = 'png', width = 700, height = 700, units = 'px')

(figure_uniform_edgesightings.3 <- ggplot()+
    geom_point(data = edge_samples1,
               aes(x = total_sightings, y = weight),
               colour = colours[1],
               size = 0.5, shape = 19)+
    geom_point(data = averages,
               aes(x = total_sightings, y = median),
               colour = colours[2],
               size = 1, shape = 19)+
    geom_smooth(data = averages,
                aes(x = total_sightings, y = median,
                    colour = as.factor(never))
    )+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'edge weight',
                       limits = c(-0.02,1.02),
                       expand = c(0,0))+
    scale_colour_manual(values = colours[c(3,5)])+
    theme(legend.position = 'bottom',
          legend.background = element_blank(),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 8))
)
ggsave(filename = 'sim_edgesightingsuniform_allpoints.png',
       path = 'methods_paper/simulation_outputs/',
       plot = last_plot(), device = 'png', width = 700, height = 700, units = 'px')

sim_dyads <- sim_dyads %>% 
  left_join(averages[,c('dyad_id','mean','median','uniform_lwr','uniform_upr')], by = 'dyad_id') %>% 
  rename(uniform_mu = mean,
         uniform_md = median)

save.image('methods_paper/simulation_outputs/plots_bisonuniform.RData')
rm(list = ls()[!ls() %in% c('sim_dyads','n_chains','plot_dyads','plot_dyad_ids','colours')]) ; gc()

######## BISoN unconditional: default normal(0,2.5) ########
#### run model ####
## compile edge model
edge_binary_default <- cmdstan_model('methods_paper/models/edge_binary_gaussian.stan')

## create data list for default model
counts_ls_default <- list(n_dyads = nrow(sim_dyads),
                          dyad_ids = sim_dyads$dyad_id,
                          together = sim_dyads$together,
                          count_dyad = sim_dyads$total_sightings)

## fit model
fit_edges_default <- edge_binary_default$sample(
  data = counts_ls_default,
  chains = n_chains,
  parallel_chains = n_chains)
save.image('methods_paper/simulation_outputs/sim_default_run.RData')

## extract edges
edges <- fit_edges_default$draws() %>% as.data.frame()
edges <- edges[,(n_chains+1):ncol(edges)]        # remove lp__ columns
edges1 <- edges[,seq(1,ncol(edges)-3, by = 4)]   # select only chain 1
edge_names <- data.frame(name = colnames(edges1)) %>%
  separate(name, into = c('chain','weight'), sep = '.edge_')
colnames(edges1) <- edge_names$weight
edges2 <- edges[,seq(2,ncol(edges)-2, by = 4)] ; colnames(edges2) <- edge_names$weight  # select only chain 2
edges3 <- edges[,seq(3,ncol(edges)-1, by = 4)] ; colnames(edges3) <- edge_names$weight  # select only chain 3
edges4 <- edges[,seq(4,ncol(edges)-0, by = 4)] ; colnames(edges4) <- edge_names$weight  # select only chain 4
edges <- rbind(edges1, edges2, edges3, edges4)
n_samples <- nrow(edges)
edge_samples1 <- edges1 %>%
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'weight')
rm(edges1, edges2, edges3, edges4) ; gc()

#### plot edges ####
## plot all edges
(figure_default_posterior <- ggplot(data = edge_samples1)+
   geom_density(aes(group = parameter, x = weight),
                colour = colours[3], alpha = 0.1)+
   scale_x_continuous(name = 'edge weight')+
   scale_y_continuous(name = 'density')
)
ggsave(filename = 'sim_posterior_default_alldyads.png',
       path = 'methods_paper/simulation_outputs/',
       plot = figure_default_posterior, device = 'png', width = 1600, height = 700, units = 'px')

## sample edges
edge_samples1$dyad_id <- rep(sim_dyads$dyad_id, n_samples/4)
edge_samples1 <- edge_samples1 %>%
  left_join(sim_dyads[,c('dyad_id','node1','node2','sightings1','sightings2','total_sightings','together','sri')],
            by = 'dyad_id') %>%
  mutate(never = ifelse(together == 0, 'never together', 'together at least once'))

set.seed(15) ; plot_dyads <- sample(1:nrow(sim_dyads), size = 100, replace = F)
plot_dyad_ids <- unique(edge_samples1$parameter)[plot_dyads]
edges_subset <- edge_samples1[edge_samples1$parameter %in% plot_dyad_ids,]

## plot sample
(figure_default_posterior <- ggplot(data = edges_subset)+
    geom_density(aes(group = parameter, x = weight, colour = never),
                 show.legend = F, linewidth = 0.2)+
    stat_density(aes(group = parameter, x = weight, colour = never),
                 geom = "line", position = "identity", linewidth = 0)+
    scale_x_continuous(name = 'edge weight', limits = c(0,1))+
    scale_y_continuous(name = 'density', limits = c(0,40))+
    scale_colour_manual(values = colours[c(3,5)],
                        aesthetics = 'colour')+
    theme(legend.position = 'bottom',
          legend.background = element_blank(),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 8))+
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))
)
ggsave(filename = 'sim_posterior_default_sampledyads.png',
       path = 'methods_paper/simulation_outputs/',
       plot = figure_default_posterior, device = 'png',
       width = 1000, height = 1000, units = 'px')

## clean up and save
save.image('methods_paper/simulation_outputs/sim_default_run.RData')

#### mean edge weight vs dyad sighting count ####
#load('methods_paper/simulation_outputs/sim_default_run.RData')
averages <- data.frame(dyad = colnames(edges),
                       mean = apply(edges, 2, mean),
                       median = apply(edges, 2, median),
                       default_lwr = apply(edges, 2, quantile, prob = 0.025),
                       default_upr = apply(edges, 2, quantile, prob = 0.975))
averages$dyad_id <- sim_dyads$dyad_id
averages <- averages %>%
  left_join(sim_dyads[,c('dyad_id','node1','node2','together','never','total_sightings','sightings1','sightings2')],
            by = 'dyad_id')

(figure_default_edgesightings.2 <- ggplot()+
    geom_point(data = averages,
               aes(x = total_sightings, y = median),
               colour = rgb(33/255, 145/255, 140/255, 0.1),
               size = 0.5,
               shape = 19)+
    geom_smooth(data = averages,
                aes(x = total_sightings, y = median, linetype = as.factor(never)),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'median weight', limits = c(-0.02,1.02), expand = c(0,0))+
    scale_linetype_manual(name = NULL,
                          values = c(1,6))+
    theme(legend.position = 'bottom',
          legend.key.height = unit(4, 'mm'),
          legend.text = element_text(size = 8)))
# ggsave(filename = 'edgesightingsdefault_twolines.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = last_plot(), device = 'png', width = 700, height = 700, units = 'px')

(figure_default_edgesightings.3 <- ggplot()+
    geom_point(data = edge_samples1,
               aes(x = total_sightings, y = weight),
               colour = colours[1],
               size = 0.5, shape = 19)+
    geom_point(data = averages,
               aes(x = total_sightings, y = median),
               colour = colours[2],
               size = 1, shape = 19)+
    geom_smooth(data = averages,
                aes(x = total_sightings, y = median,
                    colour = as.factor(never))
    )+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'edge weight',
                       limits = c(-0.02,1.02),
                       expand = c(0,0))+
    scale_colour_manual(values = colours[c(3,5)])+
    theme(legend.position = 'bottom',
          legend.background = element_blank(),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 8))
)
ggsave(filename = 'sim_edgesightingsdefault_allpoints.png',
       path = 'methods_paper/simulation_outputs/',
       plot = last_plot(), device = 'png', width = 700, height = 700, units = 'px')

sim_dyads <- sim_dyads %>% 
  left_join(averages[,c('dyad_id','mean','median','default_lwr','default_upr')], by = 'dyad_id') %>% 
  rename(default_mu = mean,
         default_md = median)

save.image('methods_paper/simulation_outputs/plots_bisondefault.RData')
rm(list = ls()[!ls() %in% c('sim_dyads','n_chains','plot_dyads','plot_dyad_ids','colours')]) ; gc()

######## BISoN unconditional: skewed beta(0.7,5) ########
#### run model ####
## compile edge model
edge_binary_skewed <- cmdstan_model('methods_paper/models/edge_binary_skewed.stan')

## create data list for skewed model
counts_ls_skewed <- list(n_dyads = nrow(sim_dyads),
                         dyad_ids = sim_dyads$dyad_id,
                         together = sim_dyads$together,
                         count_dyad = sim_dyads$total_sightings)

## fit model
fit_edges_skewed <- edge_binary_skewed$sample(
  data = counts_ls_skewed,
  chains = n_chains,
  parallel_chains = n_chains)
save.image('methods_paper/simulation_outputs/sim_skewed_run.RData')

## extract edges
edges <- fit_edges_skewed$draws() %>% as.data.frame()
edges <- edges[,(n_chains+1):ncol(edges)]        # remove lp__ columns
edges1 <- edges[,seq(1,ncol(edges)-3, by = 4)]   # select only chain 1
edge_names <- data.frame(name = colnames(edges1)) %>%
  separate(name, into = c('chain','weight'), sep = '.edge_')
colnames(edges1) <- edge_names$weight
edges2 <- edges[,seq(2,ncol(edges)-2, by = 4)] ; colnames(edges2) <- edge_names$weight  # select only chain 2
edges3 <- edges[,seq(3,ncol(edges)-1, by = 4)] ; colnames(edges3) <- edge_names$weight  # select only chain 3
edges4 <- edges[,seq(4,ncol(edges)-0, by = 4)] ; colnames(edges4) <- edge_names$weight  # select only chain 4
edges <- rbind(edges1, edges2, edges3, edges4)
n_samples <- nrow(edges)
edge_samples1 <- edges1 %>%
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'weight')
rm(edges1, edges2, edges3, edges4) ; gc()

#### plot edges ####
## plot all edges
(figure_skewed_posterior <- ggplot(data = edge_samples1)+
   geom_density(aes(group = parameter, x = weight),
                colour = colours[3], alpha = 0.1)+
   scale_x_continuous(name = 'edge weight')+
   scale_y_continuous(name = 'density')
)
ggsave(filename = 'sim_posterior_skewed_alldyads.png',
       path = 'methods_paper/simulation_outputs/',
       plot = figure_skewed_posterior, device = 'png', width = 1600, height = 700, units = 'px')

## sample edges
edge_samples1$dyad_id <- rep(sim_dyads$dyad_id, n_samples/4)
edge_samples1 <- edge_samples1 %>%
  left_join(sim_dyads[,c('dyad_id','node1','node2','sightings1','sightings2','total_sightings','together','sri')],
            by = 'dyad_id') %>%
  mutate(never = ifelse(together == 0, 'never together', 'together at least once'))

set.seed(15) ; plot_dyads <- sample(1:nrow(sim_dyads), size = 100, replace = F)
plot_dyad_ids <- unique(edge_samples1$parameter)[plot_dyads]
edges_subset <- edge_samples1[edge_samples1$parameter %in% plot_dyad_ids,]

## plot sample
(figure_skewed_posterior <- ggplot(data = edges_subset)+
    geom_density(aes(group = parameter, x = weight, colour = never),
                 show.legend = F, linewidth = 0.2)+
    stat_density(aes(group = parameter, x = weight, colour = never),
                 geom = "line", position = "identity", linewidth = 0)+
    scale_x_continuous(name = 'edge weight', limits = c(0,1))+
    scale_y_continuous(name = 'density', limits = c(0,40))+
    scale_colour_manual(values = colours[c(3,5)],
                        aesthetics = 'colour')+
    theme(legend.position = 'bottom',
          legend.background = element_blank(),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 8))+
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))
)
ggsave(filename = 'sim_posterior_skewed_sampledyads.png',
       path = 'methods_paper/simulation_outputs/',
       plot = figure_skewed_posterior, device = 'png',
       width = 1000, height = 1000, units = 'px')

## clean up and save
save.image('methods_paper/simulation_outputs/sim_skewed_run.RData')

#### mean edge weight vs dyad sighting count ####
#load('methods_paper/simulation_outputs/sim_skewed_run.RData')
averages <- data.frame(dyad = colnames(edges),
                       mean = apply(edges, 2, mean),
                       median = apply(edges, 2, median),
                       skewed_lwr = apply(edges, 2, quantile, prob = 0.025),
                       skewed_upr = apply(edges, 2, quantile, prob = 0.975))
averages$dyad_id <- sim_dyads$dyad_id
averages <- averages %>%
  left_join(sim_dyads[,c('dyad_id','node1','node2','together','never','total_sightings','sightings1','sightings2')],
            by = 'dyad_id')

(figure_skewed_edgesightings.2 <- ggplot()+
    geom_point(data = averages,
               aes(x = total_sightings, y = median),
               colour = rgb(33/255, 145/255, 140/255, 0.1),
               size = 0.5,
               shape = 19)+
    geom_smooth(data = averages,
                aes(x = total_sightings, y = median, linetype = as.factor(never)),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'median weight', limits = c(-0.02,1.02), expand = c(0,0))+
    scale_linetype_manual(name = NULL,
                          values = c(1,6))+
    theme(legend.position = 'bottom',
          legend.key.height = unit(4, 'mm'),
          legend.text = element_text(size = 8)))
# ggsave(filename = 'edgesightingsskewed_twolines.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = last_plot(), device = 'png', width = 700, height = 700, units = 'px')

figure_skewed_edgesightings.3 <- ggplot()+
  geom_point(data = edge_samples1,
             aes(x = total_sightings, y = weight),
             colour = colours[1],
             size = 0.5, shape = 19)+
  geom_point(data = averages,
             aes(x = total_sightings, y = median),
             colour = colours[2],
             size = 1, shape = 19)+
  geom_smooth(data = averages,
              aes(x = total_sightings, y = median,
                  colour = as.factor(never))
  )+
  scale_x_continuous(name = 'total dyad sightings')+
  scale_y_continuous(name = 'edge weight',
                     limits = c(-0.02,1.02),
                     expand = c(0,0))+
  scale_colour_manual(values = colours[c(3,5)])+
  theme(legend.position = 'bottom',
        legend.background = element_blank(),
        legend.key.height = unit(4, 'mm'),
        legend.title = element_blank(),
        legend.text = element_text(size = 8))
figure_skewed_edgesightings.3
ggsave(filename = 'sim_edgesightingsskewed_allpoints.png',
       path = 'methods_paper/simulation_outputs/',
       plot = figure_skewed_edgesightings.3, device = 'png', width = 700, height = 700, units = 'px')

sim_dyads <- sim_dyads %>% 
  left_join(averages[,c('dyad_id','mean','median','skewed_lwr','skewed_upr')], by = 'dyad_id') %>% 
  rename(skewed_mu = mean,
         skewed_md = median)

save.image('methods_paper/simulation_outputs/plots_bisonskewed.RData')
# rm(list = ls()[!ls() %in% c('sim_dyads','n_chains','plot_dyads','plot_dyad_ids','colours')]) ; gc()

######## BISoN unconditional: merge plots ########
## load all data -- change name of sim_dyads or it will remove the information for default and skewed by overwriting with uniform only
#load('methods_paper/simulation_outputs/plots_bisonskewed.RData')
all_data <- sim_dyads
rm(list = ls()[!ls() %in% c('all_data','n_chains','plot_dyads','plot_dyad_ids','colours',
                            'figure_skewed_posterior','figure_skewed_edgesightings.3')]) ; gc()
load('methods_paper/simulation_outputs/plots_bisondefault.RData')
rm(list = ls()[!ls() %in% c('all_data','n_chains','plot_dyads','plot_dyad_ids','colours',
                            'figure_default_posterior','figure_default_edgesightings.3',
                            'figure_skewed_posterior','figure_skewed_edgesightings.3')]) ; gc()
load('methods_paper/simulation_outputs/plots_bisonuniform.RData')
sim_dyads <- all_data
rm(list = ls()[!ls() %in% c('sim_dyads','n_chains','plot_dyads','plot_dyad_ids','colours',
                            'figure_uniform_posterior','figure_uniform_edgesightings.3',
                            'figure_default_posterior','figure_default_edgesightings.3',
                            'figure_skewed_posterior','figure_skewed_edgesightings.3')]) ; gc()
save.image('methods_paper/simulation_outputs/plots_unconditional.RData') # load('methods_paper/simulation_outputs/plots_unconditional.RData')

(figure_uniform_posterior <- figure_uniform_posterior+
    scale_y_continuous(limits = c(0, 15)))
(figure_default_posterior <- figure_default_posterior+
    scale_y_continuous(limits = c(0, 15)))
(figure_skewed_posterior <- figure_skewed_posterior+
    scale_y_continuous(limits = c(0, 15)))

## combine posterior distributions
(combined_top <- (figure_uniform_posterior + 
                    figure_default_posterior + 
                    figure_skewed_posterior) +
    plot_annotation(tag_levels = 'a') +
    plot_layout(guides = 'collect') &
    theme(legend.position = 'bottom'))
ggsave(filename = 'sim_posterior_unconditional.png',
       path = 'methods_paper/simulation_outputs/',
       plot = last_plot(), device = 'png', width = 2700, height = 700, units = 'px')

## combine edges vs sightings
combined_bottom <- (figure_uniform_edgesightings.3 + figure_default_edgesightings.3 + figure_skewed_edgesightings.3)+
  plot_annotation(tag_levels = list(c('d','e','f')))
(combined_bottom <- combined_bottom + 
    plot_layout(guides = 'collect') & theme(legend.position = 'bottom'))
ggsave(filename = 'sim_edgesightingsunconditional.png',
       path = 'methods_paper/simulation_outputs/',
       plot = last_plot(), device = 'png', width = 2700, height = 700, units = 'px')

## combine everything into one single plot
# (figure_uniform_posterior + figure_default_posterior + figure_skewed_posterior) /
#   (figure_uniform_edgesightings.3 + figure_default_edgesightings.3 + figure_skewed_edgesightings.3)+
#   plot_annotation(tag_levels = 'a')
(combined_top / combined_bottom)
ggsave(filename = 'sim_outputs_unconditional.png',
       path = 'methods_paper/simulation_outputs/',
       plot = last_plot(), device = 'png', width = 2700, height = 1600, units = 'px')

save.image('methods_paper/simulation_outputs/plots_unconditional_combined.RData')
rm(list = ls()[!ls() %in% c('sim_dyads','n_chains','plot_dyads','plot_dyad_ids','colours')]) ; gc()

######## BISoN conditional ########
#### run model ####
## compile edge model
edge_binary_conditional <- cmdstan_model('methods_paper/models/edge_binary_conditional.stan')

## create data list for conditional model
counts_ls_conditional <- list(n_dyads = nrow(sim_dyads),
                              dyad_ids = sim_dyads$dyad_id,
                              together = sim_dyads$together,
                              count_dyad = sim_dyads$total_sightings)

## fit model
fit_edges_conditional <- edge_binary_conditional$sample(
  data = counts_ls_conditional,
  chains = n_chains,
  parallel_chains = n_chains)
save.image('methods_paper/simulation_outputs/sim_conditional_run.RData')

## extract edges
edges <- fit_edges_conditional$draws() %>% as.data.frame()
edges <- edges[,(n_chains+1):ncol(edges)]        # remove lp__ columns
edges1 <- edges[,seq(1,ncol(edges)-3, by = 4)]   # select only chain 1
edge_names <- data.frame(name = colnames(edges1)) %>%
  separate(name, into = c('chain','weight'), sep = '.edge_')
colnames(edges1) <- edge_names$weight
edges2 <- edges[,seq(2,ncol(edges)-2, by = 4)] ; colnames(edges2) <- edge_names$weight  # select only chain 2
edges3 <- edges[,seq(3,ncol(edges)-1, by = 4)] ; colnames(edges3) <- edge_names$weight  # select only chain 3
edges4 <- edges[,seq(4,ncol(edges)-0, by = 4)] ; colnames(edges4) <- edge_names$weight  # select only chain 4
edges <- rbind(edges1, edges2, edges3, edges4)
n_samples <- nrow(edges)
edge_samples1 <- edges1 %>%
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'weight')
rm(edges1, edges2, edges3, edges4) ; gc()

#### plot edges ####
## plot all edges
(posterior_conditional.1 <- ggplot(data = edge_samples1)+
   geom_density(aes(group = parameter, x = weight),
                colour = colours[3], alpha = 0.1)+
   scale_x_continuous(name = 'edge weight')+
   scale_y_continuous(name = 'density')
)
ggsave(filename = 'sim_posterior_conditional_alldyads.png',
       path = 'methods_paper/simulation_outputs/',
       plot = posterior_conditional.1, device = 'png', width = 1600, height = 700, units = 'px')

## sample edges
edge_samples1$dyad_id <- rep(sim_dyads$dyad_id, n_samples/4)
edge_samples1 <- edge_samples1 %>%
  left_join(sim_dyads[,c('dyad_id','node1','node2','sightings1','sightings2','total_sightings','together','sri')],
            by = 'dyad_id') %>%
  mutate(never = ifelse(together == 0, 'never together', 'together at least once'))

set.seed(15) ; plot_dyads <- sample(1:nrow(sim_dyads), size = 100, replace = F)
plot_dyad_ids <- unique(edge_samples1$parameter)[plot_dyads]
edges_subset <- edge_samples1[edge_samples1$parameter %in% plot_dyad_ids,]

## plot sample
(posterior_conditional.2 <- ggplot(data = edges_subset)+
    geom_density(aes(group = parameter, x = weight, colour = never),
                 show.legend = F, linewidth = 0.2)+
    stat_density(aes(group = parameter, x = weight, colour = never),
                 geom = "line", position = "identity", linewidth = 0)+
    scale_x_continuous(name = 'edge weight', limits = c(0,1))+
    scale_y_continuous(name = 'density', limits = c(0,40))+
    scale_colour_manual(values = colours[c(3,5)],
                        aesthetics = 'colour')+
    theme(legend.position = 'bottom',
          legend.background = element_blank(),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 8))+
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))
)
ggsave(filename = 'sim_posterior_conditional_sampledyads.png',
       path = 'methods_paper/simulation_outputs/',
       plot = posterior_conditional.2, device = 'png',
       width = 1000, height = 1000, units = 'px')

## clean up and save
save.image('methods_paper/simulation_outputs/sim_conditional_run.RData')

#### mean edge weight vs dyad sighting count ####
#load('methods_paper/simulation_outputs/sim_conditional_run.RData')
averages <- data.frame(dyad = colnames(edges),
                       mean = apply(edges, 2, mean),
                       median = apply(edges, 2, median),
                       conditional_lwr = apply(edges, 2, quantile, prob = 0.025),
                       conditional_upr = apply(edges, 2, quantile, prob = 0.975))
averages$dyad_id <- sim_dyads$dyad_id
averages <- averages %>%
  left_join(sim_dyads[,c('dyad_id','node1','node2','together','never','total_sightings','sightings1','sightings2')],
            by = 'dyad_id')

(edgesightingsconditional.2 <- ggplot()+
    geom_point(data = averages,
               aes(x = total_sightings, y = median),
               colour = rgb(33/255, 145/255, 140/255, 0.1),
               size = 0.5,
               shape = 19)+
    geom_smooth(data = averages,
                aes(x = total_sightings, y = median, linetype = as.factor(never)),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'median weight', limits = c(-0.02,1.02), expand = c(0,0))+
    scale_linetype_manual(name = NULL,
                          values = c(1,6))+
    theme(legend.position = 'bottom',
          legend.key.height = unit(4, 'mm'),
          legend.text = element_text(size = 8)))
# ggsave(filename = 'edgesightingsconditional_twolines.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = last_plot(), device = 'png', width = 700, height = 700, units = 'px')

(edgesightingsconditional.3 <- ggplot()+
    geom_point(data = edge_samples1,
               aes(x = total_sightings, y = weight),
               colour = colours[1],
               size = 0.5, shape = 19)+
    geom_point(data = averages,
               aes(x = total_sightings, y = median),
               colour = colours[2],
               size = 1, shape = 19)+
    geom_smooth(data = averages,
                aes(x = total_sightings, y = median,
                    colour = as.factor(never))
    )+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'edge weight',
                       limits = c(-0.02,1.02),
                       expand = c(0,0))+
    scale_colour_manual(values = colours[c(3,5)])+
    theme(legend.position = 'bottom',
          legend.background = element_blank(),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 8))
)
ggsave(filename = 'sim_edgesightingsconditional_allpoints.png',
       path = 'methods_paper/simulation_outputs/',
       plot = last_plot(), device = 'png', width = 700, height = 700, units = 'px')

sim_dyads <- sim_dyads %>% 
  left_join(averages[,c('dyad_id','mean','median','conditional_lwr','conditional_upr')], by = 'dyad_id') %>% 
  rename(conditional_mu = mean,
         conditional_md = median)

save.image('methods_paper/simulation_outputs/plots_bisonconditional.RData')

#### merge ########
# load('methods_paper/simulation_outputs/plots_bisonconditional.RData')
rm(list = ls()[ !ls() %in% c('sim_dyads','colours',
                             'posterior_conditional.1','posterior_conditional.2',
                             'edgesightingsconditional.2','edgesightingsconditional.3',
                             'n_chains','n_samples','plot_dyad_ids','plot_dyads')]) ; gc()

# run patchwork and save
post_nolegend <- posterior_conditional.2 + theme(legend.position = 'none')
post_nolegend

combined <- (post_nolegend + edgesightingsconditional.3) +
  plot_annotation(tag_levels = 'a')
(combined + plot_layout(guides = 'collect') & theme(legend.position = 'bottom'))
ggsave(filename = 'outputs_conditional.png',
       path = 'methods_paper/simulation_outputs/',
       plot = last_plot(), device = 'png', width = 1700, height = 1000, units = 'px')

save.image('methods_paper/simulation_outputs/plots_bisonconditional.RData')
rm(list = ls()[ !ls() %in% c('sim_dyads','colours', 'n_chains','n_samples','plot_dyad_ids','plot_dyads')]) ; gc()

######## plot true vs modelled ########
## create new data frame which includes average values for all methods
plot_results <- sim_dyads %>% 
  pivot_longer(cols = c('sri',
                        'uniform_mu','uniform_md',#'uniform_lwr','uniform_upr',
                        'default_mu','default_md',#'default_lwr','default_upr',
                        'skewed_mu','skewed_md',#'skewed_lwr','skewed_upr',
                        'conditional_mu','conditional_md'#,'conditional_lwr','conditional_upr'
  ),
  names_to = 'measure', values_to = 'edge_estimate') %>% 
  separate(measure, into = c('method','average'), remove = T) %>% 
  mutate(average = ifelse(method == 'sri', 'sri', 
                          ifelse(average == 'mu', 'mean', 'median')))

## plot using mean as point estimate
(mu1 <- plot_results %>% 
    filter(average == 'mean' | average == 'sri') %>% 
    ggplot()+
    geom_point(aes(x = true_edge, y = edge_estimate, colour = method),
               alpha = 0.1)+
    geom_smooth(aes(x = true_edge, y = edge_estimate, colour = method), method = 'lm')+
    geom_abline(slope = 1, intercept = 0)+
    labs(x = 'true edge weight',
         y = 'estimated edge weight')+
    scale_colour_viridis_d()+
    facet_wrap(. ~ method)
)

## add error bars
(mu2 <- plot_results %>% 
    filter(average != 'median') %>% 
    ggplot()+
    geom_errorbar(data = plot_results[plot_results$average == 'mean' & plot_results$method == 'uniform',],
                  aes(x = true_edge, ymin = uniform_lwr, ymax = uniform_upr),
                  alpha = 0.1, colour = colours[1])+
    geom_errorbar(data = plot_results[plot_results$average == 'mean' & plot_results$method == 'default',],
                  aes(x = true_edge, ymin = default_lwr, ymax = default_upr),
                  alpha = 0.1, colour = colours[4])+
    geom_errorbar(data = plot_results[plot_results$average == 'mean' & plot_results$method == 'skewed',],
                  aes(x = true_edge, ymin = skewed_lwr, ymax = skewed_upr),
                  alpha = 0.1, colour = colours[3])+
    geom_errorbar(data = plot_results[plot_results$average == 'mean' & plot_results$method == 'conditional',],
                  aes(x = true_edge, ymin = conditional_lwr, ymax = conditional_upr),
                  alpha = 0.1, colour = colours[5])+
    geom_point(aes(x = true_edge, y = edge_estimate),
               alpha = 0.1, shape = 1)+
    geom_smooth(aes(x = true_edge, y = edge_estimate, group = method), method = 'lm', colour = 'black')+
    geom_abline(slope = 1, intercept = 0, linetype = 2)+
    labs(x = 'true edge weight',
         y = 'estimated edge weight')+
    scale_colour_viridis_d()+
    facet_wrap(. ~ method)
)
ggsave(plot = mu2, device = 'png', width = 3500, height = 2200, unit = 'px',
       filename = 'sim_compare_outputs_to_true.png',
       path = 'methods_paper/simulation_outputs/')

