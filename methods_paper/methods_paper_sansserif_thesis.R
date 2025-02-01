## take graphs from methods paper and plot each using png instead of png/svg
#### Set up                                             ####
## load packages needed -- library(cmdstanr) ; library(tidyverse) ; library(lubridate) ; library(janitor) ; library(hms) ; library(readxl) ; library(data.table) ; library(spatsoc) ; library(igraph) ; library(patchwork)
library(cmdstanr, lib.loc = '../packages/')      # run models
library(tidyverse, lib.loc = '../packages')      # data manipulation
library(lubridate, lib.loc = '../packages')      # sort date columns out
library(janitor, lib.loc = '../packages')        # clean up data frame names
library(hms, lib.loc = '../packages')            # sort time columns out
library(readxl, lib.loc = '../packages')         # read in .xlxs files
library(data.table, lib.loc = '../packages')     # convert data frames to data tables
library(spatsoc, lib.loc = '../packages')        # create social network data
library(igraph, lib.loc = '../packages/')        # extract eigenvector centralities
library(patchwork, lib.loc = '../packages/')     # plot multiple graphs together
library(LaplacesDemon, lib.loc = '../packages/') # logit and invlogit functions

### set stan path
set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0')

### set seed
set.seed(12345)

### define plot theme
theme_set(theme_bw(base_size = 12))

### set height and width
width <- 750      # min = 789, max = 2250
height <- 700     # max = 2625
resolution <- 300 # min = 300, max = 600, but changes depending on dimensions

#-------- 8) Create SRI plots for paper                 ##########
load('methods_paper/r_workspaces/sri.RData')

#### SRI Edge weights: edges and edge_vs_sightings
edges_sri <- edges_sri+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
edgesightings_sri.2 <- edgesightings_sri.2 +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

(combined <- (edges_sri + edgesightings_sri.2)+
    plot_annotation(tag_levels = 'a'))
ggsave(filename = 'outputs_sri_sansserif.png',
       path = 'methods_paper/outputs/',
       plot = last_plot(), device = 'png',
       width = width*3,
       height = height,
       dpi = resolution,
       units = 'px')

## print progress marker
print('SRI edges plotted')

#### SRI Eigenvector Centrality: eigenvector vs sighting count and non-associations
eigen0_sri.2 <- eigen0_sri.2 +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
eigensightings_sri <- eigensightings_sri +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

(eigen0_sri.2 + eigensightings_sri)+
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'eigen_outputs_sri_sansserif.png',
       path = 'methods_paper/outputs/',
       plot = last_plot(), device = 'png',
       width = (width+50)*2,
       height = height,
       dpi = resolution,
       units = 'px')

## print progress marker
print('SRI eigenvector plotted')

#-------- 9) Create Unconditional BISoN plots for paper ####
#### Unconditional Priors                               ####
x <- seq(0, 1, length = 100)

#### Uniform, uniform(0,1)
# calculate probability of all values of x assuming completely flat prior
flat_prior <- dunif(x = x, min = 0, max = 1)
data <- data.frame(x = x, density = flat_prior)

# plot
(priors_sri <- ggplot(data)+
    geom_line(aes(x = x, y = density),
              linewidth = 1.2,
              colour = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density', limits = c(0,1.05), expand = c(0,0)) +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))
  
)

#### Default, normal(0,2.5)
# calculate probability of all values of x assuming default bisonR prior shape
default_prior <- c(0, 2.5) # obtained using: default_prior <- bisonR::get_default_priors('binary')
default_prior <- dnorm(x = logit(x), mean = default_prior[1], sd = default_prior[2])
data <- data.frame(x = x, density = default_prior)

# plot
(priors_default <- ggplot(data)+
    geom_line(aes(x = x, y = density),
              linewidth = 1.2,
              colour = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density',
                       limits = c(0,0.2),
                       #limits = c(0,1),
                       expand = c(0,0)) +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))
)

#### Skewed, beta(0.7,5)
# calculate probability of all values of x assuming right skewed prior shape
skewed_prior <- c(1, 5)
skewed_prior <- dbeta(x = x, shape1 = skewed_prior[1], shape2 = skewed_prior[2])
data <- data.frame(x = x, density = skewed_prior)

# plot
(priors_skewed <- ggplot(data)+
    geom_line(aes(x = x, y = density),
              linewidth = 1.2,
              colour = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density',
                       #limits = c(0,0.3),
                       #limits = c(0,1),
                       expand = c(0,0)) +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))
)

#### Merge
(priors_sri + priors_default + priors_skewed) +
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'priors_unconditional_sansserif.png',
       path = 'methods_paper/outputs/',
       plot = last_plot(), device = 'png',
       width = width*3,
       height = height,
       dpi = resolution,
       units = 'px')

## print progress marker
print('priors plotted')

rm(list = ls()[! ls() %in% c('width','height','resolution')]) ; gc()

#### Unconditional Outputs                              ####
load('methods_paper/r_workspaces/plots_unconditional.RData')
figure_uniform_posterior <- figure_uniform_posterior +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(colour = '')
figure_default_posterior <- figure_default_posterior +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(colour = '')
figure_skewed_posterior <- figure_skewed_posterior +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(colour = '')

(combined_top <- (figure_uniform_posterior +
                    figure_default_posterior +
                    figure_skewed_posterior) +
    labs(colour = '') +
    plot_annotation(tag_levels = 'a') +
    plot_layout(guides = 'collect') &
    theme(legend.position = 'bottom'))
# ggsave(filename = 'posterior_unconditional_noinset_sansserif.png',
#        path = 'methods_paper/outputs/',
#        plot = combined_top, device = 'png',
#        width = width*3,
#        height = height+100,
#        dpi = resolution,
#        units = 'px')

figure_uniform_edgesightings.3 <- figure_uniform_edgesightings.3 +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(colour = '')
figure_default_edgesightings.3 <- figure_default_edgesightings.3 +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(colour = '')
figure_skewed_edgesightings.3 <- figure_skewed_edgesightings.3 +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  labs(colour = '')

combined_bottom <- (figure_uniform_edgesightings.3 + figure_default_edgesightings.3 + figure_skewed_edgesightings.3)+
  plot_annotation(tag_levels = list(c('d','e','f')))
(combined_bottom <- combined_bottom +
    plot_layout(guides = 'collect') & theme(legend.position = 'bottom'))
# ggsave(filename = 'edgesightings_unconditional_noinset_sansserif.png',
#        path = 'methods_paper/outputs/',
#        plot = combined_bottom, device = 'png',
#        width = width*3,
#        height = height+100,
#        dpi = resolution,
#        units = 'px')

(combined_top / combined_bottom) +
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'outputs_unconditional_noinset_sansserif.png',
       path = 'methods_paper/outputs/',
       plot = last_plot(), device = 'png',
       width = width*3,
       height = (height+100)*2,
       dpi = resolution,
       units = 'px')

# uniform_posterior_inset <- wrap_elements(figure_uniform_posterior +
#                                            inset_element(inset_sri,
#                                                          left = 0.5, right = 0.99,
#                                                          bottom = 0.5, top = 0.99) +
#                                            theme(legend.position = 'none'))
# default_posterior_inset <- wrap_elements(figure_default_posterior +
#                                            inset_element(inset_default,
#                                                          left = 0.5, right = 0.99,
#                                                          bottom = 0.5, top = 0.99) +
#                                            theme(legend.position = 'bottom',
#                                                  legend.title = element_blank()))
# skewed_posterior_inset <- wrap_elements(figure_skewed_posterior +
#                                           inset_element(inset_skewed,
#                                                         left = 0.5, right = 0.99,
#                                                         bottom = 0.5, top = 0.99) +
#                                           theme(legend.position = 'none'))
#
# combined_top <- (uniform_posterior_inset + default_posterior_inset + skewed_posterior_inset)+
#   labs(colour = '') +
#   plot_annotation(tag_levels = 'a') +
#   plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
# (combined_top <- combined_top +
#     theme(legend.position = 'bottom',
#           legend.title = element_blank()) +
#     plot_layout(guides = 'collect')
# )
# ggsave(filename = 'posterior_unconditional_inset_sansserif.png',
#        path = 'methods_paper/outputs/',
#        plot = combined_top, device = 'png',
#        width = width*3,
#        height = height+100,
#        dpi = resolution,
#        units = 'px')
#
# (combined_top / combined_bottom)
# ggsave(filename = 'outputs_unconditional_inset_sansserif.png',
#        path = 'methods_paper/outputs/',
#        plot = last_plot(), device = 'png',
#        width = width*3,
#        height = (height+100)*2,
#        dpi = resolution,
#        units = 'px')

rm(list = ls()[! ls() %in% c('width','height','resolution')]) ; gc()

#### Unconditional Eigenvector Centrality               ####
load('methods_paper/r_workspaces/eigen_checks_unconditional.RData')
rm(list = ls()[!ls() %in% c('eigen0_uniform.2','eigensightings_uniform.2',
                            'eigen0_default.2','eigensightings_default.2',
                            'eigen0_skewed.2', 'eigensightings_skewed.2',
                            'inset_sri','inset_default','inset_skewed',
                            'eigen_sri','elephants','colours','N',
                            'width','height','resolution')]) ; gc()

## plot without insets included
eigen0_uniform.2 <- eigen0_uniform.2 +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
eigen0_default.2 <- eigen0_default.2 +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
eigen0_skewed.2 <- eigen0_skewed.2 +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

eigensightings_uniform.2 <- eigensightings_uniform.2 +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
eigensightings_default.2 <- eigensightings_default.2 +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
eigensightings_skewed.2 <- eigensightings_skewed.2 +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

(eigen0_uniform.2 + eigen0_default.2 + eigen0_skewed.2) /
  (eigensightings_uniform.2 + eigensightings_default.2 + eigensightings_skewed.2)+
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'eigen_outputs_unconditional_noinset_sansserif.png',
       path = 'methods_paper/outputs/',
       plot = last_plot(), device = 'png',
       width = width*3,
       height = (height+100)*2,
       dpi = resolution,
       units = 'px')

# ## set up insets
# uniform_eigen0_inset <- wrap_elements(eigen0_uniform.2 +
#                                         inset_element(inset_sri,
#                                                       left = 0.01, right = 0.25,
#                                                       bottom = 0.75, top = 0.99))
# default_eigen0_inset <- wrap_elements(eigen0_default.2 +
#                                         inset_element(inset_default,
#                                                       left = 0.01, right = 0.25,
#                                                       bottom = 0.75, top = 0.99))
# skewed_eigen0_inset <- wrap_elements(eigen0_skewed.2 +
#                                        inset_element(inset_skewed,
#                                                      left = 0.01, right = 0.25,
#                                                      bottom = 0.75, top = 0.99))
#
# ## plot with priors inset into top row
# (uniform_eigen0_inset + default_eigen0_inset + skewed_eigen0_inset) /
#   (eigensightings_uniform.2 + eigensightings_default.2 + eigensightings_skewed.2)+
#   plot_annotation(tag_levels = 'a')
# ggsave(filename = 'eigen_outputs_unconditional_inset_sansserif.png',
#        path = 'methods_paper/outputs/',
#        plot = last_plot(), device = 'png',
#        width = width*3,
#        height = (height+100)*2,
#        dpi = resolution,
#        units = 'px')

## print progress marker
print('unconditional eigenvector plotted')

rm(list = ls()[! ls() %in% c('width','height','resolution')]) ; gc()

#-------- 10) Create Conditional BISoN plots for paper  ####
#### Conditional Prior                                  ####
# define sequence over which to plot
x <- seq(0, 1, length = 100)

# calculate probability of all values of x assuming completely flat prior
conditional1 <- dbeta(x = x, shape1 = 0.7, shape2 = 10)
conditional2 <- dbeta(x = x, shape1 = 1, shape2 = 5)
data <- data.frame(x = x,
                   density1 = conditional1,
                   density2 = conditional2) %>%
  pivot_longer(cols = c(density1, density2),
               names_to = 'density',
               values_to = 'y') %>%
  mutate(together = ifelse(density == 'density1', 'never together',
                           'together at least once'))

# plot
(prior_conditional <- ggplot(data)+
    geom_line(aes(x = x, y = y, colour = together),
              linetype = 1,
              linewidth = 1.2)+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density',
                       expand = c(0,0),
                       limits = c(0,15))+
    scale_colour_manual(name = NULL,
                        values = c(rgb(33/255, 145/255, 140/255),
                                   rgb(68/255,1/255,84/255,1)))+
    theme(legend.position = 'right',
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))
)
ggsave(filename = 'prior_conditional_sansserif.png',
       path = 'methods_paper/outputs/',
       plot = prior_conditional, device = 'png',
       width = (width)*2,
       height = height,
       dpi = resolution,
       units = 'px')

## print progress marker
print('conditional prior plotted')

#### Conditional Outputs                                ####
load('methods_paper/r_workspaces/plots_conditional.RData')
rm(list = ls()[ !ls() %in% c('counts_df','n_chains',
                             'eigen_sri','elephants','colours','N',
                             'figure_conditional_posterior',
                             'figure_conditional_edgesightings.1',
                             'figure_conditional_edgesightings.2',
                             'figure_conditional_edgesightings.3',
                             'inset_conditional','fit_edges_conditional',
                             'width','height','resolution')]) ; gc()

figure_conditional_posterior <- figure_conditional_posterior +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        # legend.title = element_blank(),
        # legend.position = 'none',
        legend.text = element_text(size = 10))+
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))
figure_conditional_edgesightings.3 <- figure_conditional_edgesightings.3 +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))+
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

## combine plots without prior inset
(combined <- (figure_conditional_posterior + figure_conditional_edgesightings.3) +
    plot_annotation(tag_levels = 'a') &
    theme(legend.position = 'bottom',
          legend.text = element_text(size = 10),
          legend.title = element_blank()) &
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)))
ggsave(filename = 'outputs_conditional_noinset_sansserif.png',
       path = 'methods_paper/outputs/',
       plot = last_plot(), device = 'png',
       width = (width+50)*2,
       height = height+100,
       dpi = resolution,
       units = 'px')

# ## add prior as inset
# conditional_posterior_inset <- wrap_elements(figure_conditional_posterior +
#                                                inset_element(inset_conditional,
#                                                              left = 0.5, right = 0.99,
#                                                              bottom = 0.5, top = 0.99)) +
#   theme(legend.position = 'none')
# combined <- (conditional_posterior_inset + figure_conditional_edgesightings.3) +
#   plot_annotation(tag_levels = 'a')
# (combined + plot_layout(guides = 'collect') &
#     theme(legend.position = 'bottom'))
# ggsave(filename = 'outputs_conditional_inset_sansserif.png',
#        path = 'methods_paper/outputs/',
#        plot = last_plot(), device = 'png',
#        width = (width+50)*2,
#        height = height,
#        dpi = resolution,
#        units = 'px')

# clean up
rm(list = ls()[! ls() %in% c('width','height','resolution')]) ; gc()

## print progress marker
print('conditional posterior plotted')

#### Conditional Eigenvector Centrality                 ####
load('methods_paper/r_workspaces/eigen_checks_conditional.RData')
rm(list = ls()[!ls() %in% c('eigen0_sri.2','eigensightings_sri.2',
                            'eigen0_uniform.2','eigen0_default.2','eigen0_skewed.2',
                            'eigensightings_uniform.2','eigensightings_default.2','eigensightings_skewed.2',
                            'eigen0_conditional.2','eigensightings_conditional.2',
                            'inset_sri','inset_default','inset_skewed','inset_conditional',
                            'eigen_sri','elephants','colours','N',
                            'width','height','resolution')]) ; gc()

# conditional
eigen0_conditional.2 <- eigen0_conditional.2 +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
eigensightings_conditional.2 <- eigensightings_conditional.2 +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

(eigen0_conditional.2 + eigensightings_conditional.2)+
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'eigen_outputs_conditional_noinset_sansserif.png',
       path = 'methods_paper/outputs/',
       plot = last_plot(), device = 'png',
       width = (width+50)*2,
       height = height,
       dpi = resolution,
       units = 'px')

# conditional_eigensightings_inset <- wrap_elements(eigensightings_conditional.2 +
#                                                     inset_element(inset_conditional,
#                                                                   left = 0.75, right = 0.99,
#                                                                   bottom = 0.75, top = 0.99))
#
# (conditional_eigensightings_inset + eigen0_conditional.2)+
#   plot_annotation(tag_levels = 'a')
# ggsave(filename = 'eigen_outputs_conditional_inset_sansserif.png',
#        path = 'methods_paper/outputs/',
#        plot = last_plot(), device = 'png',
#        width = (width+50)*2,
#        height = height,
#        dpi = resolution,
#        units = 'px')

## print progress marker
print('conditional eigenvector plotted')

#-------- 11) Supplementary material: simulation        ####
rm(list = ls()[! ls() %in% c('width','height','resolution')]) ; gc()

# define sequence over which to plot
x <- seq(0, 1, length = 100)

# identify probability of each value of x depending on shape of distributions
distributions <- data.frame(x = x,
                            rare2 = dbeta(x = x, shape1 = 0.1, shape2 = 1.9),
                            rare10 = dbeta(x = x, shape1 = 1, shape2 = 9),
                            sometimes2 = dbeta(x = x, shape1 = 1, shape2 = 1),
                            sometimes10 = dbeta(x = x, shape1 = 10, shape2 = 10),
                            usual2 = dbeta(x = x, shape1 = 1.9, shape2 = 0.1),
                            usual10 = dbeta(x = x, shape1 = 9, shape2 = 1)) %>%
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

(simulation <- ggplot()+
    geom_vline(data = true_values, aes(xintercept = x,
                                       colour = 'True edge weight'),
               linewidth = 1.2)+
    geom_line(data = distributions, aes(x = x,
                                        y = dbeta/max(dbeta[dbeta != 'Inf']),
                                        colour = 'Bayesian posterior distribution'),
              linewidth = 1.2,
              key_glyph = 'rect')+
    geom_line(data = sri, aes(x = x, y = y_value, group = pairing,
                              colour = 'SRI probability'),
              linewidth = 1.2)+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density', limits = c(0,1))+
    facet_grid(labels_sightings ~ labels_together)+
    scale_color_manual(name='Edge weight',
                       breaks=c('True edge weight',
                                'SRI probability',
                                'Bayesian posterior distribution'),
                       values=c('True edge weight'='#fde725',
                                'SRI probability'='#440154',
                                'Bayesian posterior distribution'='#21918c'))+
    theme(legend.position = 'bottom',
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))
)
ggsave(filename = 'simulation_sansserif.png',
       path = 'methods_paper/outputs/',
       plot = simulation,
       device = 'png',
       width = width*3,
       height = (height+100)*2,
       dpi = resolution,
       units = 'px')

#-------- 12) Supplementary material: mixture model     ####
rm(list = ls()[ !ls() %in% c('counts_df','n_chains',
                             'eigen_sri','elephants','colours','N',
                             'figure_conditional_posterior',
                             'figure_conditional_edgesightings.1',
                             'figure_conditional_edgesightings.2',
                             'figure_conditional_edgesightings.3',
                             'inset_conditional','fit_edges_conditional',
                             'width','height','resolution')]) ; gc()
load('methods_paper/r_workspaces/mixture_model.RData')

#### Mixture Outputs                                    ####
edges$ever_together <- ifelse(edges$event_count == 0, 'never together', 'together at least once')
(edges_mixture <- ggplot(edges)+
    geom_bar(aes(x = round(prdctn_invlogit, 3),
                 fill = ever_together,
                 colour = ever_together),
             # fill = rgb(33/255, 145/255, 140/255),
             # colour = rgb(33/255, 145/255, 140/255),
             linewidth = 0.2)+
    scale_x_continuous(name = 'mixture model weight',
                       breaks = c(0.00,0.10,0.20))+
    scale_y_continuous(name = 'number of dyads',
                       expand = c(0,0))+
    scale_colour_viridis_d(end = 0.5, direction = -1)+
    scale_fill_viridis_d(end = 0.5, direction = -1)+
    coord_cartesian(ylim = c(0,5300))+
    theme(legend.position = 'bottom')+
    labs(fill = NULL, colour = NULL)+
    guides(colour = guide_legend(nrow = 2),
           fill = guide_legend(nrow = 2))
)

(edgesightings_mixture <- ggplot()+
    geom_point(data = edges,
               aes(x = count_dyad, y = prdctn_invlogit),
               colour = rgb(94/255, 201/255, 98/255, 0.2),#rgb(33/255, 145/255, 140/255, 0.1),
               size = 1, shape = 19)+
    geom_smooth(data = edges[edges$ever_together == 'together at least once',],
                aes(x = count_dyad, y = prdctn_invlogit,
                    colour = ever_together)
    )+
    geom_smooth(data = edges[edges$ever_together == 'never together',],
                aes(x = count_dyad, y = prdctn_invlogit,
                    colour = ever_together)
    )+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mixture weight', limits = c(-0.02,1.02), expand = c(0,0))+
    scale_colour_manual(values = colours,
                        breaks = c('never together','together at least once'),
                        name = NULL)+
    theme(legend.position = 'right',#c(0.56,0.76),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8))
)

(edges_mixture + edgesightings_mixture)+
  plot_annotation(tag_levels = 'a') &
  theme(legend.position = 'bottom',
        legend.text = element_text(size = 10),
        legend.title = element_blank()) &
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))
ggsave(filename = 'outputs_mixturemodel.tiff',
       path = 'methods_paper/outputs/sans_serif/',
       plot = last_plot(), device = 'tiff',
       width = (width+50)*2,
       height = height+100,
       dpi = resolution,
       units = 'px')

#### Mixture Eigenvector Centrality                     ####
eigen0_mix.2 <- eigen0_mix.2 +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
eigensightings_mix.2 <- eigensightings_mix.2 +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

(eigen0_mix.2 + eigensightings_mix.2) +
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'eigen_outputs_mixture.png',
       path = 'methods_paper/outputs/',
       plot = last_plot(), device = 'png',
       width = (width+50)*2,
       height = height,
       dpi = resolution,
       units = 'px')
