#### information #####
# script to graph random intercepts of node and window from ANP nodal regression

#### set up #####
## load other packages
# library(StanHeaders) ; library(rstan) ; library(LaplacesDemon) ; library(MASS) ; library(tidyverse)
library(StanHeaders, lib.loc = '../packages/')
library(rstan, lib.loc = '../packages/')
library(tidyverse, lib.loc = '../packages/')
library(LaplacesDemon, lib.loc = '../packages/')
library(MASS, lib.loc = '../packages/')

## set theme for plots
theme_set(theme_bw())

## set seed
set.seed(12345)

#### LONG:  calculate proportion of variation explained by window random effect ####
load('anp_nodalregression/anp_long_nodal.RData')

## calculate values
mean_window_logit_long <- apply(rand_window, 2, mean)
stdv_window_logit_long <- apply(rand_window, 2, sd)

## convert to data frame
mean_window_logit_long <- mean_window_logit_long %>%
  as.data.frame() %>%
  mutate(eles_per_window = as.vector(unlist(model_data_list[grep(pattern = 'num_nodes_window',
                                                                 x = names(model_data_list))]))) %>%
  mutate(window = 1:7) %>%
  relocate(window) %>%
  mutate(stdv_effect = as.vector(stdv_window_logit_long))
colnames(mean_window_logit_long)[2] <- 'mean_effect'

## plot
mean_window_logit_long %>%
  ggplot()+
  geom_errorbar(aes(x = eles_per_window,
                    ymax = mean_effect + stdv_effect,
                    ymin = mean_effect - stdv_effect,
                    colour = as.factor(window)))+
  scale_colour_viridis_d()+
  geom_point(aes(x = eles_per_window,
                 #colour = as.factor(window),
                 y = mean_effect),
             pch = 21,
             colour = 'black',
             fill = 'white')+
  geom_smooth(aes(x = eles_per_window,
                  y = mean_effect),
              colour = 'black')+
  labs(x = 'number of elephants in time window',
       y = 'mean effect of window ID (logit scale)')+
  guides(colour = guide_legend(ncol = 2, nrow = 18,
                               title = 'window ID',
                               keyheight = 1))+
  theme(legend.position = 'right',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 14))
ggsave(plot = last_plot(), device = 'png',
       filename = 'anplong_window_random_effect.png',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')
ggsave(plot = last_plot(), device = 'svg',
       filename = 'anplong_window_random_effect.svg',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')

#### LONG:  calculate proportion of variation explained by node random effect ####
## calculate values
mean_node_logit_long <- apply(rand_node, 2, mean)
stdv_node_logit_long <- apply(rand_node, 2, sd)

## convert to data frame
counts_long <- nodes_all %>%
  mutate(observed = ifelse(sightings == 0, 0, 1)) %>%
  group_by(node) %>%
  summarise(total_sightings = sum(sightings),
            total_windows = sum(observed))
mean_node_logit_long <- mean_node_logit_long %>%
  as.data.frame() %>%
  mutate(node = sort(unique(nodes_all$node))) %>%
  relocate(node) %>%
  left_join(counts_long, by = 'node') %>%
  mutate(stdv_effect = as.vector(stdv_node_logit_long))
colnames(mean_node_logit_long)[2] <- 'mean_effect'

## plot
mean_node_logit_long %>%
  ggplot()+
  geom_errorbar(aes(x = total_sightings,#total_windows,
                    colour = node,
                    ymax = mean_effect + stdv_effect,
                    ymin = mean_effect - stdv_effect))+
  scale_colour_viridis_c(breaks = c(1,seq(25,700, by = 25)),
                         name = 'node ID')+
  geom_point(aes(x = total_sightings,#total_windows,
                 #colour = node,
                 y = mean_effect),
             pch = 21,
             colour = 'black',
             fill = 'white')+
  geom_smooth(aes(x = total_sightings,#total_windows,
                  y = mean_effect),
              colour = 'black')+
  labs(x = 'total number observations per individual',#'number of windows in which individual was observed',
       y = 'mean effect of node ID (logit scale)')+
  guides(colour = guide_legend(ncol = 2,
                               keyheight = 1))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))
ggsave(plot = last_plot(), device = 'png',
       filename = 'anplong_node_random_effect.png',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')
ggsave(plot = last_plot(), device = 'svg',
       filename = 'anplong_node_random_effect.svg',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')

## plot both
random_effects_long <- mean_node_logit_long %>%
  dplyr::select(-node, -total_sightings, -total_windows) %>%
  mutate(type = 'node')
window_effect_long <- mean_window_logit_long %>%
  dplyr::select(-window, -eles_per_window) %>%
  mutate(type = 'window')
random_effects_long <- rbind(random_effects_long,window_effect_long)
rm(window_effect_long)
ggplot()+
  geom_violin(data = random_effects_long,
              aes(x = type,
                  y = mean_effect))+
  geom_jitter(data = random_effects_long,
              aes(x = type,
                  y = mean_effect))

#### SHORT: calculate proportion of variation explained by window random effect ####
rm(list = ls()[!ls() %in% c('mean_node_logit_long','mean_window_logit_long')])
mean_node_logit_long <- mean_node_logit_long %>%
  mutate(window_length = 'long')
mean_window_logit_long <- mean_window_logit_long %>%
  mutate(window_length = 'long')
load('anp_nodalregression/anp_short_nodal.RData')

## calculate values
mean_window_logit_short <- apply(rand_window, 2, mean)
stdv_window_logit_short <- apply(rand_window, 2, sd)

## convert to data frame
mean_window_logit_short <- mean_window_logit_short %>%
  as.data.frame() %>%
  mutate(eles_per_window = as.vector(unlist(model_data_list[grep(pattern = 'num_nodes_window',
                                                                 x = names(model_data_list))]))) %>%
  mutate(window = 1:36) %>%
  relocate(window) %>%
  mutate(stdv_effect = as.vector(stdv_window_logit_short))
colnames(mean_window_logit_short)[2] <- 'mean_effect'

## plot
mean_window_logit_short %>%
  ggplot()+
  geom_errorbar(aes(x = eles_per_window,
                    ymax = mean_effect + stdv_effect,
                    ymin = mean_effect - stdv_effect,
                    colour = as.factor(window)))+
  scale_colour_viridis_d()+
  geom_point(aes(x = eles_per_window,
                 #colour = as.factor(window),
                 y = mean_effect),
             pch = 21,
             colour = 'black',
             fill = 'white')+
  geom_smooth(aes(x = eles_per_window,
                  y = mean_effect),
              colour = 'black')+
  labs(x = 'number of elephants in time window',
       y = 'mean effect of window ID (logit scale)')+
  guides(colour = guide_legend(ncol = 2, nrow = 18,
                               title = 'window ID',
                               keyheight = 1))+
  theme(legend.position = 'right',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 14))
ggsave(plot = last_plot(), device = 'png',
       filename = 'anpshort_window_random_effect.png',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')
ggsave(plot = last_plot(), device = 'svg',
       filename = 'anpshort_window_random_effect.svg',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')

#### SHORT: calculate proportion of variation explained by node random effect ####
## calculate values
mean_node_logit_short <- apply(rand_node, 2, mean)
stdv_node_logit_short <- apply(rand_node, 2, sd)

## convert to data frame
counts_short <- nodes_all %>%
  mutate(observed = ifelse(sightings == 0, 0, 1)) %>%
  group_by(node) %>%
  summarise(total_sightings = sum(sightings),
            total_windows = sum(observed))
mean_node_logit_short <- mean_node_logit_short %>%
  as.data.frame() %>%
  mutate(node = sort(unique(nodes_all$node))) %>%
  relocate(node) %>%
  left_join(counts_short, by = 'node') %>%
  mutate(stdv_effect = as.vector(stdv_node_logit_short))
colnames(mean_node_logit_short)[2] <- 'mean_effect'

## plot
mean_node_logit_short %>%
  ggplot()+
  geom_errorbar(aes(x = total_sightings,#total_windows,
                    colour = node,
                    ymax = mean_effect + stdv_effect,
                    ymin = mean_effect - stdv_effect))+
  scale_colour_viridis_c(breaks = c(1,seq(25,700, by = 25)),
                         name = 'node ID')+
  geom_point(aes(x = total_sightings,#total_windows,
                 #colour = node,
                 y = mean_effect),
             pch = 21,
             colour = 'black',
             fill = 'white')+
  geom_smooth(aes(x = total_sightings,#total_windows,
                  y = mean_effect),
              colour = 'black')+
  labs(x = 'total number observations per individual',#'number of windows in which individual was observed',
       y = 'mean effect of node ID (logit scale)')+
  guides(colour = guide_legend(ncol = 2,
                               keyheight = 1))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))
ggsave(plot = last_plot(), device = 'png',
       filename = 'anpshort_node_random_effect.png',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')
ggsave(plot = last_plot(), device = 'svg',
       filename = 'anpshort_node_random_effect.svg',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')

## plot both
random_effects_short <- mean_node_logit_short %>%
  dplyr::select(-node, -total_sightings, -total_windows) %>%
  mutate(type = 'node')
window_effect_short <- mean_window_logit_short %>%
  dplyr::select(-window, -eles_per_window) %>%
  mutate(type = 'window')
random_effects_short <- rbind(random_effects_short,window_effect_short)
rm(window_effect_short)
ggplot()+
  geom_violin(data = random_effects_short,
              aes(x = type,
                  y = mean_effect))+
  geom_jitter(data = random_effects_short,
              aes(x = type,
                  y = mean_effect))

#### plot together ####
rm(list = ls()[!ls() %in% c('mean_node_logit_long','mean_window_logit_long',
                            'mean_node_logit_short','mean_window_logit_short')])
mean_node_logit_short <- mean_node_logit_short %>%
  mutate(window_length = 'short')
mean_window_logit_short <- mean_window_logit_short %>%
  mutate(window_length = 'short')
save.image('anp_nodalregression/plot_random_effects.RData')

## combine data frames
node <- rbind(mean_node_logit_long, mean_node_logit_short)
window <- rbind(mean_window_logit_long, mean_window_logit_short)

## plot
node %>%
  ggplot()+
  geom_errorbar(aes(x = total_sightings,#total_windows,
                    colour = node,
                    ymax = mean_effect + stdv_effect,
                    ymin = mean_effect - stdv_effect))+
  scale_colour_viridis_c(breaks = c(1,seq(25,700, by = 25)),
                         name = 'node ID')+
  geom_point(aes(x = total_sightings,#total_windows,
                 #colour = node,
                 y = mean_effect),
             pch = 21,
             colour = 'black',
             fill = 'white')+
  geom_smooth(aes(x = total_sightings,#total_windows,
                  y = mean_effect),
              colour = 'black')+
  labs(x = 'total number observations per individual',#'number of windows in which individual was observed',
       y = 'mean effect of node ID (logit scale)')+
  guides(colour = guide_legend(ncol = 2,
                               keyheight = 1))+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  facet_wrap(. ~ factor(window_length, levels = c('short','long')))
ggsave(plot = last_plot(), device = 'png',
       filename = 'anp_node_random_effect_longshort.png',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')
ggsave(plot = last_plot(), device = 'svg',
       filename = 'anp_node_random_effect_longshort.svg',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')

## plot
window %>%
  ggplot()+
  geom_errorbar(aes(x = eles_per_window,
                    ymax = mean_effect + stdv_effect,
                    ymin = mean_effect - stdv_effect,
                    colour = as.factor(window)))+
  scale_colour_viridis_d()+
  geom_point(aes(x = eles_per_window,
                 #colour = as.factor(window),
                 y = mean_effect),
             pch = 21,
             colour = 'black',
             fill = 'white')+
  geom_smooth(aes(x = eles_per_window,
                  y = mean_effect),
              colour = 'black')+
  labs(x = 'number of elephants in time window',
       y = 'mean effect of window ID (logit scale)')+
  guides(colour = guide_legend(ncol = 2, nrow = 18,
                               title = 'window ID',
                               keyheight = 1))+
  theme(legend.position = 'right',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 14))+
  facet_wrap(. ~ factor(window_length,
                        levels = c('short','long')),
             scales = 'free_x')
ggsave(plot = last_plot(), device = 'png',
       filename = 'anp_window_random_effect_longshort.png',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')
ggsave(plot = last_plot(), device = 'svg',
       filename = 'anp_window_random_effect_longshort.svg',
       path = '../outputs/step4_nodalregression/',
       height = 1800, width = 2400, units = 'px')
