#### script for plotting graphs required for methods paper ####
# library(cmdstanr) ; library(tidyverse) ; library(LaplacesDemon) ; library(patchwork) ; library(igraph) ; library(bisonR)
library(cmdstanr, lib.loc = '../packages/')
library(tidyverse, lib.loc = '../packages/')
library(bisonR, lib.loc = '../packages/')
library(LaplacesDemon, lib.loc = '../packages/')
library(patchwork, lib.loc = '../packages/')
library(igraph, lib.loc = '../packages/')

set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0/')

n_chains <- 4

# define plot theme
theme_set(theme_bw(base_size = 12))

# #### unconditional priors: 3 panel plot of uniform/default BISoN/right-skewed BISoN for priors ####
# x <- seq(0, 1, length = 100)
# 
# ## SRI, prior distribution ####
# # calculate probability of all values of x assuming completely flat prior
# flat_prior <- dunif(x = x, min = 0, max = 1)
# data <- data.frame(x = x, density = flat_prior)
# 
# # plot
# (priors_sri <- ggplot(data)+
#     geom_line(aes(x = x, y = density), linewidth = 1.2, colour = rgb(33/255, 145/255, 140/255))+
#     scale_x_continuous(name = 'edge weight')+
#     scale_y_continuous(name = 'density', limits = c(0,1.05), expand = c(0,0))
# )
# ggsave(filename = 'priors_sri.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = priors_sri, device = 'png', width = 700, height = 800, units = 'px')
# 
# ## default, prior distribution ####
# # calculate probability of all values of x assuming default bisonR prior shape
# default_prior <- get_default_priors('binary') # default edge prior = normal(0, 2.5)
# default_prior <- c(0, 2.5)
# default_prior <- dnorm(x = logit(x), mean = default_prior[1], sd = default_prior[2])
# data <- data.frame(x = x, density = default_prior)
# 
# # plot
# (priors_default <- ggplot(data)+
#     geom_line(aes(x = x, y = density), linewidth = 1.2, colour = rgb(33/255, 145/255, 140/255))+
#     scale_x_continuous(name = 'edge weight')+
#     scale_y_continuous(name = 'density',
#                        limits = c(0,0.2),
#                        #limits = c(0,1),
#                        expand = c(0,0))
# )
# ggsave(filename = 'priors_default.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = priors_default, device = 'png', width = 700, height = 800, units = 'px')
# 
# ## skewed, prior distribution ####
# # calculate probability of all values of x assuming right skewed prior shape
# skewed_prior <- c(1, 5)
# skewed_prior <- dbeta(x = x, shape1 = skewed_prior[1], shape2 = skewed_prior[2])
# data <- data.frame(x = x, density = skewed_prior)
# 
# # plot
# (priors_skewed <- ggplot(data)+
#     geom_line(aes(x = x, y = density), linewidth = 1.2, colour = rgb(33/255, 145/255, 140/255))+
#     scale_x_continuous(name = 'edge weight')+
#     scale_y_continuous(name = 'density',
#                        #limits = c(0,0.3),
#                        #limits = c(0,1),
#                        expand = c(0,0))
# )
# ggsave(filename = 'priors_skewed.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = priors_skewed, device = 'png', width = 700, height = 800, units = 'px')
# 
# # clean up
# rm(list = ls()[!ls() %in% c('priors_sri','priors_default','priors_skewed')]) ; gc()
# 
# ## merge ####
# (priors_sri + priors_default + priors_skewed) +
#   plot_annotation(tag_levels = 'a')
# ggsave(filename = 'priors_unconditional.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = last_plot(), device = 'png', width = 2100, height = 800, units = 'px')
# 
# #### conditional prior: 1 panel plot of conditional prior ####
# # define sequence over which to plot
# x <- seq(0, 1, length = 100)
# 
# # calculate probability of all values of x assuming completely flat prior
# conditional1 <- dbeta(x = x, shape1 = 0.7, shape2 = 10)
# conditional2 <- dbeta(x = x, shape1 = 1, shape2 = 5)
# data <- data.frame(x = x,
#                    density1 = conditional1,
#                    density2 = conditional2)
# # plot
# (prior_conditional <- ggplot(data)+
#     geom_line(aes(x = x, y = density1), linewidth = 1.2, linetype = 1, colour = rgb(33/255, 145/255, 140/255))+
#     geom_line(aes(x = x, y = density2), linewidth = 1.2, linetype = 2, colour = rgb(68/255,1/255,84/255,1))+
#     scale_x_continuous(name = 'edge weight')+
#     scale_y_continuous(name = 'density', expand = c(0,0))+
#     coord_cartesian(ylim = c(0,15))
# )
# ggsave(filename = 'prior_conditional.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = prior_conditional, device = 'png', width = 700, height = 800, units = 'px')
# 
# # clean up
# rm(list = ls()) ; gc()
# 
# #### SRI outputs: 2 panel plot of SRI edges and edge_vs_sightings ####
# ## SRI edge bar graph ####
# # read in MOTNP data
# load('motnp_edgeweights_conditionalprior.RData')
# rm(edge_binary, edge_samples, edgelist, edges, fit_edges_motnp)
# 
# # calculate sri
# counts_df$sri <- counts_df$event_count / counts_df$count_dyad
# 
# # calculate percentages of 0s and 1s using SRI
# ( length(which(counts_df$sri==0))/length(counts_df$sri) ) * 100
# ( length(which(counts_df$sri==1))/length(counts_df$sri) ) * 100
# length(unique(c(counts_df$id_1, counts_df$id_2)))
# 
# # calculate percentages of 0s and 1s using SRI, assuming a 5 sighting threshold per elephant
# subset5 <- counts_df %>% filter(count_1 >= 5 & count_2 >= 5)
# ( length(which(subset5$sri==0))/length(subset5$sri) ) * 100
# ( length(which(subset5$sri==1))/length(subset5$sri) ) * 100
# length(unique(c(subset5$id_1, subset5$id_2)))
# 
# # calculate percentages of 0s and 1s using SRI, assuming a 10 sighting threshold per elephant
# subset10 <- counts_df %>% filter(count_1 >= 10 & count_2 >= 10)
# ( length(which(subset10$sri==0))/length(subset10$sri) ) * 100
# length(unique(c(subset10$id_1, subset10$id_2)))
# 100-(( length(unique(c(subset10$id_1, subset10$id_2))) / length(unique(c(counts_df$id_1, counts_df$id_2))) )*100)
# 
# # bar plot
# (edges_sri <- ggplot()+
#     geom_bar(data = counts_df, aes(x = round(sri, 3)), fill = rgb(33/255, 145/255, 140/255),
#              colour = rgb(33/255, 145/255, 140/255), linewidth = 0.2)+
#     #geom_bar(data = counts_df[counts_df$count_dyad > 10,], aes(x = round(sri, 3)), fill = '#fde725', colour = 'transparent')+
#     scale_x_continuous(name = 'SRI value')+
#     scale_y_continuous(name = 'number of dyads',
#                        #limits = c(0,1000),
#                        expand = c(0,0))+
#     annotate('text', x = 0.42, y = 165,
#              label = paste0(length(which(counts_df$sri == 0)), ' dyads with \nSRI = 0'),
#              size = unit(4, 'pt'),
#              colour = rgb(33/255, 145/255, 140/255))+
#     coord_cartesian(ylim = c(0,200))
# )
# ggsave(filename = 'edges_sri.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = edges_sri, device = 'png', width = 700, height = 800, units = 'px')
# 
# # redo annotation on figure 2d so that size isn't too large for multi plot
# (edges_sri <- ggplot()+
#     geom_bar(data = counts_df, aes(x = round(sri, 3)), fill = rgb(33/255, 145/255, 140/255),
#              colour = rgb(33/255, 145/255, 140/255), linewidth = 0.2)+
#     scale_x_continuous(name = 'SRI value')+
#     scale_y_continuous(name = 'number of dyads',
#                        expand = c(0,0))+
#     annotate('text', x = 0.5, y = 165,
#              label = paste0(length(which(counts_df$sri == 0)), ' dyads with \nSRI = 0'),
#              size = unit(4, 'pt'),
#              colour = rgb(33/255, 145/255, 140/255))+
#     coord_cartesian(ylim = c(0,200))
# )
# 
# ## SRI edges vs sighting count ####
# (edgesightings_sri.1 <- ggplot(counts_df, aes(x = count_dyad, y = sri))+
#    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
#    scale_x_continuous(name = 'total dyad sightings')+
#    scale_y_continuous(name = 'SRI weight', limits = c(-0.02,1.02), expand = c(0,0))
# )
# ggsave(filename = 'edgesightings_sri_noline.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = edgesightings_sri.1, device = 'png', width = 700, height = 800, units = 'px')
# 
# (edgesightings_sri.2 <- ggplot(counts_df, aes(x = count_dyad, y = sri))+
#     geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
#     geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
#     scale_x_continuous(name = 'total dyad sightings')+
#     scale_y_continuous(name = 'SRI weight', limits = c(-0.02,1.02), expand = c(0,0))
# )
# ggsave(filename = 'edgesightings_sri_withline_changealpha.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = edgesightings_sri.2, device = 'png', width = 700, height = 800, units = 'px')
# 
# counts_df$sri_rank <- as.integer(as.factor(counts_df$sri))
# p7 <- as.data.frame(table(counts_df$count_dyad, counts_df$sri)) %>%
#   mutate(count_dyad = as.numeric(Var1),
#          sri_rank = as.numeric(Var2)) %>%
#   filter(Freq > 0) %>%
#   left_join(distinct(counts_df[,c('sri','sri_rank')]), by = 'sri_rank')
# (edgesightings_sri.3 <- ggplot()+
#     geom_point(data = p7, aes(x = count_dyad, y = sri, size = Freq),
#                colour = rgb(33/255, 145/255, 140/255))+
#     geom_smooth(data = counts_df, aes(x = count_dyad, y = sri),
#                 colour = rgb(68/255, 1/255, 84/255))+
#     theme(legend.position = c(0.6,0.7),
#           legend.background = element_rect(fill = 'white', colour = 'black'),
#           legend.key.height = unit(4, 'mm'),
#           legend.title = element_text(size = 10),
#           legend.text = element_text(size = 8))+
#     scale_size_continuous(range = c(0.2, 2))+
#     labs(size = 'number of dyads')+
#     scale_x_continuous(name = 'total dyad sightings')+
#     scale_y_continuous(name = 'SRI weight')
# )
# ggsave(filename = 'edgesightings_sri_withline_changesize.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = edgesightings_sri.3, device = 'png', width = 700, height = 800, units = 'px')
# 
# # # clean up
# rm(list = ls()[!ls() %in% c('edges_sri','edgesightings_sri.2', 'counts_df')]) ; gc()
# 
# ## merge ####
# (edges_sri + edgesightings_sri.2)+
#   plot_annotation(tag_levels = 'a')
# ggsave(filename = 'outputs_sri.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = last_plot(), device = 'png', width = 1400, height = 800, units = 'px')
# 
# #### unconditional outputs: 6 panel plot of uniform/default/skewed for posterior and edge_vs_sightings ####
# ## uniform, posterior distribution ####
# rm(list = ls()) ; gc()
# load('motnp_edgeweights_conditionalprior.RData')
# rm(list = ls()[! ls() %in% c('counts_df','n_chains','eigen_sri')])
# 
# # calculate sri
# counts_df$sri <- counts_df$event_count / counts_df$count_dyad
# print('sri calculated in counts_df')
# 
# # compile edge model
# edge_binary_uniform <- cmdstan_model('other/methods_paper/edge_binary_uniform.stan')
# 
# # create data list for uniform model
# counts_ls_uniform <- list(n_dyads = nrow(counts_df),
#                           dyad_ids = counts_df$dyad_id,
#                           together = counts_df$event_count,
#                           count_dyad = counts_df$count_dyad)
# 
# # fit model
# fit_edges_uniform <- edge_binary_uniform$sample(
#   data = counts_ls_uniform,
#   chains = n_chains,
#   parallel_chains = n_chains)
# 
# # extract edges
# edges <- fit_edges_uniform$draws() %>% as.data.frame()
# edges <- edges[,(n_chains+1):ncol(edges)]        # remove lp__ columns
# edges1 <- edges[,seq(1,ncol(edges)-3, by = 4)]   # select only chain 1
# edge_names <- data.frame(name = colnames(edges1)) %>%
#   separate(name, into = c('chain','weight'), sep = '.edge_')
# colnames(edges1) <- edge_names$weight
# edges2 <- edges[,seq(2,ncol(edges)-2, by = 4)] ; colnames(edges2) <- edge_names$weight  # select only chain 2
# edges3 <- edges[,seq(3,ncol(edges)-1, by = 4)] ; colnames(edges3) <- edge_names$weight   # select only chain 3
# edges4 <- edges[,seq(4,ncol(edges)-0, by = 4)] ; colnames(edges4) <- edge_names$weight   # select only chain 4
# edges <- rbind(edges1, edges2, edges3, edges4)
# n_samples <- nrow(edges)
# edge_samples1 <- edges1 %>%
#   pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'weight')
# rm(edges1, edges2, edges3, edges4) ; gc()
# 
# # plot
# (figure_uniform_posterior <- ggplot(data = edge_samples1)+
#     geom_density(aes(group = parameter, x = weight),
#                  colour = rgb(33/255, 145/255, 140/255, 0.1))+
#     scale_x_continuous(name = 'edge weight')+
#     scale_y_continuous(name = 'density')
# )
# ggsave(filename = 'posterior_uniform_alldyads.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = figure_uniform_posterior, device = 'png', width = 1400, height = 800, units = 'px')
# 
# edge_samples1$dyad_males <- rep(1:nrow(counts_df), n_samples/4)
# edge_samples1 <- edge_samples1 %>%
#   left_join(counts_df[,c('dyad_males','id_1','id_2','count_1','count_2','count_dyad','event_count','sri')],
#             by = 'dyad_males') %>%
#   mutate(together0 = ifelse(event_count == 0, 'never together', 'together at least once'))
# 
# set.seed(15) ; plot_dyads <- sample(1:nrow(counts_df), size = 100, replace = F)
# plot_dyad_ids <- unique(edge_samples1$parameter)[plot_dyads]
# edges_subset <- edge_samples1[edge_samples1$parameter %in% plot_dyad_ids,]
# (figure_uniform_posterior <- ggplot(data = edges_subset)+
#     geom_density(aes(group = parameter, x = weight, colour = together0), show.legend = F)+
#     stat_density(aes(group = parameter, x = weight, colour = together0),
#                  geom = "line", position = "identity", linewidth = 0)+
#     scale_x_continuous(name = 'edge weight', limits = c(0,1))+
#     scale_y_continuous(name = 'density', limits = c(0,40))+
#     scale_colour_manual(values = c(rgb(33/255, 145/255, 140/255, 0.5), rgb(68/255, 1/255, 84/255, 0.5)),
#                         aesthetics = 'colour')+
#     theme(legend.position = c(0.5,0.75),
#           legend.background = element_rect(fill = 'white', colour = 'black'),
#           legend.key.height = unit(4, 'mm'),
#           legend.title = element_text(size = 10),
#           legend.text = element_text(size = 8))+
#     guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))+
#     labs(colour = 'sightings together')
# )
# ggsave(filename = 'posterior_uniform_sampledyads.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = figure_uniform_posterior, device = 'png', width = 1400, height = 800, units = 'px')
# 
# # clean up and save
# save.image('../outputs/sparse_network_methods_figures/model_run_uniform.RData')
# 
# ## uniform, mean edge weight vs dyad sighting count ####
# load('../outputs/sparse_network_methods_figures/model_run_uniform.RData')
# averages <- data.frame(dyad = colnames(edges),
#                        mean = apply(edges, 2, mean),
#                        median = apply(edges, 2, median))
# averages$dyad_males <- 1:nrow(averages)
# 
# averages <- averages %>%
#   left_join(counts_df[,c('dyad_id','node_1_males','node_2_males','event_count','count_dyad','count_1','count_2','dyad_males')],
#             by = 'dyad_males')
# 
# (figure_uniform_edgesightings.1 <- ggplot(averages, aes(x = count_dyad, y = median))+
#     geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1),
#                size = 0.5,
#                shape = 19)+
#     scale_x_continuous(name = 'total dyad sightings')+
#     scale_y_continuous(name = 'median weight', limits = c(-0.02,1.02), expand = c(0,0))
# )
# ggsave(filename = 'edgesightings_uniform_noline.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = figure_uniform_edgesightings.1, device = 'png', width = 1400, height = 800, units = 'px')
# 
# (figure_uniform_edgesightings.2 <- ggplot(averages, aes(x = count_dyad, y = median))+
#     geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1),
#                size = 0.5,
#                shape = 19)+
#     geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
#     scale_x_continuous(name = 'total dyad sightings')+
#     scale_y_continuous(name = 'median weight', limits = c(-0.02,1.02), expand = c(0,0))
# )
# ggsave(filename = 'edgesightings_uniform_withline.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = figure_uniform_edgesightings.2, device = 'png', width = 1400, height = 800, units = 'px')
# 
# (figure_uniform_edgesightings.3 <- ggplot()+
#     geom_point(data = edge_samples1, aes(x = count_dyad, y = weight),
#                colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
#     geom_point(data = averages, aes(x = count_dyad, y = median),
#                colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
#     geom_smooth(data = averages, aes(x = count_dyad, y = median), colour = rgb(68/255, 1/255, 84/255))+
#     scale_x_continuous(name = 'total dyad sightings')+
#     scale_y_continuous(name = 'edge weight', limits = c(-0.02,1.02), expand = c(0,0))
# )
# ggsave(filename = 'edgesightings_uniform_allpoints_withline.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = figure_uniform_edgesightings.3, device = 'png', width = 1400, height = 800, units = 'px')
# 
# save.image('../outputs/sparse_network_methods_figures/plots_bisonuniform.RData')
# rm(list = ls()) ; gc()
# 
# ## default normal(0,2.5), posterior distribution ####
# # compile edge model
# edge_binary_default <- cmdstan_model('other/methods_paper/edge_binary_gaussian.stan')
# 
# # load data
# counts_ls_default <- list(n_dyads = nrow(counts_df),
#                           dyad_ids = counts_df$dyad_id,
#                           together = counts_df$event_count,
#                           count_dyad = counts_df$count_dyad)
# 
# # fit model
# fit_edges_default <- edge_binary_default$sample(
#   data = counts_ls_default,
#   chains = n_chains,
#   parallel_chains = n_chains)
# 
# # extract edges
# edges <- fit_edges_default$draws() %>% as.data.frame()
# edges <- edges[,(n_chains+1):ncol(edges)]        # remove lp__ columns
# edges1 <- edges[,seq(1,ncol(edges)-3, by = 4)]   # select only chain 1
# edge_names <- data.frame(name = colnames(edges1)) %>%
#   separate(name, into = c('chain','weight'), sep = '.edge_')
# colnames(edges1) <- edge_names$weight
# edges2 <- edges[,seq(2,ncol(edges)-2, by = 4)] ; colnames(edges2) <- edge_names$weight  # select only chain 2
# edges3 <- edges[,seq(3,ncol(edges)-1, by = 4)] ; colnames(edges3) <- edge_names$weight   # select only chain 3
# edges4 <- edges[,seq(4,ncol(edges)-0, by = 4)] ; colnames(edges4) <- edge_names$weight   # select only chain 4
# edges <- rbind(edges1, edges2, edges3, edges4)
# n_samples <- nrow(edges)
# edge_samples1 <- edges1 %>%
#   pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'weight')
# rm(edges1, edges2, edges3, edges4) ; gc()
# 
# # plot
# (figure_default_posterior <- ggplot(data = edge_samples1)+
#     geom_density(aes(group = parameter, x = weight),
#                  colour = rgb(33/255, 145/255, 140/255, 0.1))+
#     scale_x_continuous(name = 'edge weight')+
#     scale_y_continuous(name = 'density')
# )
# ggsave(filename = 'posterior_default_alldyads.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = figure_default_posterior, device = 'png', width = 1400, height = 800, units = 'px')
# 
# edge_samples1$dyad_males <- rep(1:nrow(counts_df), n_samples/4)
# edge_samples1 <- edge_samples1 %>%
#   left_join(counts_df[,c('dyad_males','id_1','id_2','count_1','count_2','count_dyad','event_count','sri')],
#             by = 'dyad_males') %>%
#   mutate(together0 = ifelse(event_count == 0, 'never together', 'together at least once'))
# 
# set.seed(15) ; plot_dyads <- sample(1:nrow(counts_df), size = 100, replace = F)
# plot_dyad_ids <- unique(edge_samples1$parameter)[plot_dyads]
# edges_subset <- edge_samples1[edge_samples1$parameter %in% plot_dyad_ids,]
# (figure_default_posterior <- ggplot(data = edges_subset)+
#     geom_density(aes(group = parameter, x = weight, colour = together0), show.legend = F)+#, colour = rgb(33/255, 145/255, 140/255, 0.1))+
#     stat_density(aes(group = parameter, x = weight, colour = together0),
#                  geom = "line", position = "identity", linewidth = 0)+
#     scale_x_continuous(name = 'edge weight', limits = c(0,1))+
#     scale_y_continuous(name = 'density', limits = c(0,40))+
#     scale_colour_manual(values = c(rgb(33/255, 145/255, 140/255, 0.5),
#                                    rgb(68/255, 1/255, 84/255, 0.5)),
#                         aesthetics = 'colour')+
#     theme(legend.position = 'none')
# )
# 
# figure_default_posterior +
#   theme(legend.position = c(0.5,0.75),
#         legend.background = element_rect(fill = 'white', colour = 'black'),
#         legend.key.height = unit(4, 'mm'),
#         legend.title = element_text(size = 10),
#         legend.text = element_text(size = 8))+
#   guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))+
#   labs(colour = 'sightings together')
# ggsave(filename = 'posterior_default_sampledyads.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = last_plot(), device = 'png', width = 1400, height = 800, units = 'px')
# 
# # clean workspace
# save.image('../outputs/sparse_network_methods_figures/model_run_default.RData')
# 
## default, mean edge weight vs dyad sighting count ####
load('../outputs/sparse_network_methods_figures/model_run_default.RData')
averages <- data.frame(dyad = colnames(edges),
                       mean = apply(edges, 2, mean),
                       median = apply(edges, 2, median))
averages$dyad_males <- 1:nrow(averages)

averages <- averages %>%
  left_join(counts_df[,c('dyad_id','node_1_males','node_2_males','event_count','count_dyad','count_1','count_2','dyad_males')],
            by = 'dyad_males')

(figure_default_edgesightings.1 <- ggplot(averages, aes(x = count_dyad, y = median))+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mean weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'edgesightings_default_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure_default_edgesightings.1, device = 'png', width = 1400, height = 800, units = 'px')

(figure_default_edgesightings.2 <- ggplot(averages, aes(x = count_dyad, y = median))+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mean weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'edgesightings_default_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure_default_edgesightings.2, device = 'png', width = 1400, height = 800, units = 'px')

(figure_default_edgesightings.3 <- ggplot()+
    geom_point(data = edge_samples1, aes(x = count_dyad, y = weight),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = averages, aes(x = count_dyad, y = median),
               colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    geom_smooth(data = averages, aes(x = count_dyad, y = median), colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'edge weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'edgesightings_default_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure_default_edgesightings.3, device = 'png', width = 1400, height = 800, units = 'px')

save.image('../outputs/sparse_network_methods_figures/plots_bisondefault.RData')
rm(list = ls()) ; gc()

# ## skewed beta(0.7, 5), posterior distribution ####
# # compile edge model
# edge_binary_skewed <- cmdstan_model('other/methods_paper/edge_binary_skewed.stan')
# 
# # respecify priors for skewed model
# counts_ls_skewed <- list(n_dyads = nrow(counts_df),
#                          dyad_ids = counts_df$dyad_id,
#                          together = counts_df$event_count,
#                          count_dyad = counts_df$count_dyad)
# 
# # fit model
# fit_edges_skewed <- edge_binary_skewed$sample(
#   data = counts_ls_skewed,
#   chains = n_chains,
#   parallel_chains = n_chains)
# 
# # extract edges
# edges <- fit_edges_skewed$draws() %>% as.data.frame()
# edges <- edges[,(n_chains+1):ncol(edges)]        # remove lp__ columns
# edges1 <- edges[,seq(1,ncol(edges)-3, by = 4)]   # select only chain 1
# edge_names <- data.frame(name = colnames(edges1)) %>%
#   separate(name, into = c('chain','weight'), sep = '.edge_')
# colnames(edges1) <- edge_names$weight
# edges2 <- edges[,seq(2,ncol(edges)-2, by = 4)] ; colnames(edges2) <- edge_names$weight  # select only chain 2
# edges3 <- edges[,seq(3,ncol(edges)-1, by = 4)] ; colnames(edges3) <- edge_names$weight   # select only chain 3
# edges4 <- edges[,seq(4,ncol(edges)-0, by = 4)] ; colnames(edges4) <- edge_names$weight   # select only chain 4
# edges <- rbind(edges1, edges2, edges3, edges4)
# n_samples <- nrow(edges)
# edge_samples1 <- edges1 %>%
#   pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'weight')
# 
# # plot
# (figure_skewed_posterior <- ggplot(data = edge_samples1)+
#     geom_density(aes(group = parameter, x = weight), colour = rgb(33/255, 145/255, 140/255, 0.1))+
#     scale_x_continuous(name = 'edge weight')+
#     scale_y_continuous(name = 'density')
# )
# ggsave(filename = 'posterior_skewed_alldyads.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = figure_skewed_posterior, device = 'png', width = 1400, height = 800, units = 'px')
# 
# edge_samples1$dyad_males <- rep(1:nrow(counts_df), n_samples/4)
# edge_samples1 <- edge_samples1 %>%
#   left_join(counts_df[,c('dyad_males','id_1','id_2','count_1','count_2','count_dyad','event_count','sri')],
#             by = 'dyad_males') %>%
#   mutate(together0 = ifelse(event_count == 0, 'never together', 'together at least once'))
# 
# set.seed(15) ; plot_dyads <- sample(1:nrow(counts_df), size = 100, replace = F)
# plot_dyad_ids <- unique(edge_samples1$parameter)[plot_dyads]
# edges_subset <- edge_samples1[edge_samples1$parameter %in% plot_dyad_ids,]
# (figure_skewed_posterior <- ggplot(data = edges_subset)+
#     geom_density(aes(group = parameter, x = weight, colour = together0), show.legend = F)+#, colour = rgb(33/255, 145/255, 140/255, 0.1))+
#     stat_density(aes(group = parameter, x = weight, colour = together0),
#                  geom = "line", position = "identity", linewidth = 0)+
#     scale_x_continuous(name = 'edge weight', limits = c(0,1))+
#     scale_y_continuous(name = 'density', limits = c(0,40))+
#     scale_colour_manual(values = c(rgb(33/255, 145/255, 140/255, 0.5),
#                                    rgb(68/255, 1/255, 84/255, 0.5)),
#                         aesthetics = 'colour')+
#     theme(legend.position = 'none')
# )
# 
# figure_skewed_posterior + 
#   theme(legend.position = c(0.5,0.75),
#         legend.background = element_rect(fill = 'white', colour = 'black'),
#         legend.key.height = unit(4, 'mm'),
#         legend.title = element_text(size = 10),
#         legend.text = element_text(size = 8))+
#   guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))+
#   labs(colour = 'sightings together')
# ggsave(filename = 'posterior_skewed_sampledyads.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = last_plot(), device = 'png', width = 1400, height = 800, units = 'px')
# 
# # clean up and save
# save.image('../outputs/sparse_network_methods_figures/model_run_skewed.RData')
# 
# ## skewed, mean edge weight vs dyad sighting count ####
# load('../outputs/sparse_network_methods_figures/model_run_skewed.RData')
# averages <- data.frame(dyad = colnames(edges),
#                        mean = apply(edges, 2, mean),
#                        median = apply(edges, 2, median))
# averages$dyad_males <- 1:nrow(averages)
# 
# averages <- averages %>% 
#   left_join(counts_df[,c('dyad_id','node_1_males','node_2_males','event_count','count_dyad','count_1','count_2','dyad_males')],
#             by = 'dyad_males')
# 
# (figure_skewed_edgesightings.1 <- ggplot(averages, aes(x = count_dyad, y = median))+
#     geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1),
#                size = 0.5,
#                shape = 19)+
#     scale_x_continuous(name = 'total dyad sightings')+
#     scale_y_continuous(name = 'median weight',
#                        limits = c(-0.02,1.02), expand = c(0,0))
# )
# ggsave(filename = 'edgesightings_skewed_noline.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = figure_skewed_edgesightings.1, device = 'png',
#        width = 1400, height = 800, units = 'px')
# 
# (figure_skewed_edgesightings.2 <- ggplot(averages,
#                                          aes(x = count_dyad, y = median))+
#     geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1),
#                size = 0.5,
#                shape = 19)+
#     geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
#     scale_x_continuous(name = 'total dyad sightings')+
#     scale_y_continuous(name = 'median weight',
#                        limits = c(-0.02,1.02), expand = c(0,0))
# )
# ggsave(filename = 'edgesightings_skewed_withline.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = figure_skewed_edgesightings.2, device = 'png',
#        width = 1400, height = 800, units = 'px')
# 
# (figure_skewed_edgesightings.3 <- ggplot()+
#     geom_point(data = edge_samples1,
#                aes(x = count_dyad, y = weight),
#                colour = rgb(253/255, 231/255, 37/255, 0.01),
#                size = 0.5, shape = 19)+
#     geom_point(data = averages, aes(x = count_dyad, y = median),
#                colour = rgb(33/255, 145/255, 140/255, 0.1),
#                size = 0.5, shape = 19)+
#     geom_smooth(data = averages, aes(x = count_dyad, y = median),
#                 colour = rgb(68/255, 1/255, 84/255))+
#     scale_x_continuous(name = 'total dyad sightings')+
#     scale_y_continuous(name = 'edge weight',
#                        limits = c(-0.02,1.02), expand = c(0,0))
# )
# ggsave(filename = 'edgesightings_skewed_allpoints_withline.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = figure_skewed_edgesightings.3, device = 'png',
#        width = 1400, height = 800, units = 'px')
# 
# ## save outputs
# save.image('../outputs/sparse_network_methods_figures/plots_bisonskewed.RData')
# 
# ## merge ####
# #rm(list = ls()) ; gc()
# #load('../outputs/sparse_network_methods_figures/plots_bisonskewed.RData')
# rm(list = ls()[!ls() %in% c('figure_skewed_posterior','figure_skewed_edgesightings.3')]) ; gc()
# load('../outputs/sparse_network_methods_figures/plots_bisondefault.RData')
# rm(list = ls()[!ls() %in% c('figure_default_posterior','figure_default_edgesightings.3',
#                             'figure_skewed_posterior','figure_skewed_edgesightings.3')]) ; gc()
# load('../outputs/sparse_network_methods_figures/plots_bisonuniform.RData')
# rm(list = ls()[!ls() %in% c('figure_uniform_posterior','figure_uniform_edgesightings.3',
#                             'figure_default_posterior','figure_default_edgesightings.3',
#                             'figure_skewed_posterior','figure_skewed_edgesightings.3')]) ; gc()
# save.image('../outputs/sparse_network_methods_figures/plots_unconditional.RData')
# 
# # load('../outputs/sparse_network_methods_figures/plots_unconditional.RData')
# (figure_uniform_posterior + figure_default_posterior + figure_skewed_posterior) /
#   (figure_uniform_edgesightings.3 + figure_default_edgesightings.3 + figure_skewed_edgesightings.3)+ 
#   plot_annotation(tag_levels = 'a')
# ggsave(filename = 'outputs_unconditional.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = last_plot(), device = 'png', width = 2700, height = 1600, units = 'px')
# 
# (figure_uniform_posterior + figure_default_posterior + figure_skewed_posterior)+
#   plot_annotation(tag_levels = 'a')
# ggsave(filename = 'posterior_unconditional.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = last_plot(), device = 'png', width = 2100, height = 800, units = 'px')
# 
# (figure_uniform_edgesightings.3 + figure_default_edgesightings.3 + figure_skewed_edgesightings.3)+ 
#   plot_annotation(tag_levels = list(c('d','e','f')))
# ggsave(filename = 'edgesightings_unconditional.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = last_plot(), device = 'png', width = 2100, height = 800, units = 'px')
# 
# rm(list = ls()) ; gc()
# 
#### conditional outputs: 2 panel plot of conditional for posterior and edge_vs_sightings ####
# ## conditional, posterior distribution ####
# load('motnp_edgeweights_conditionalprior.RData')
# counts <- counts_df %>% 
#   select(dyad_id, count_dyad, event_count) %>% 
#   rename(dyad = dyad_id)
# edges <- edges %>% 
#   filter(chain == 'chain1') %>% 
#   left_join(counts, by = 'dyad')
# edges$together <- ifelse(edges$event_count == 0, 'never together', 'together at least once')
# rm(edge_binary, edge_samples, edgelist, motnp_ages, nodes, summary, x, n_dyads, n_samples, make_edgelist, plot_network_threshold, i) ; gc()
# 
# # plot
# (posterior_conditional.1 <- ggplot(data = edges)+
#     geom_density(aes(group = dyad, x = edge_draw, colour = together), show.legend = FALSE)+
#     stat_density(aes(group = dyad, x = edge_draw, colour = together),
#                  geom = "line", position = "identity", linewidth = 0)+
#     scale_x_continuous(name = 'edge weight', limits = c(-0.05, 1.05), expand = c(0,0))+
#     scale_y_continuous(name = 'density', limits = c(-2, 72), expand = c(0,0))+
#     scale_colour_manual(values = c(rgb(33/255, 145/255, 140/255, 0.05),
#                                    rgb(68/255, 1/255, 84/255, 0.05)),
#                         aesthetics = 'colour')+
#     theme(legend.position = c(0.6,0.7),
#           legend.background = element_rect(fill = 'white', colour = 'black'),
#           legend.key.height = unit(4, 'mm'),
#           legend.title = element_text(size = 10),
#           legend.text = element_text(size = 8))+
#     guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))+
#     labs(colour = 'sightings together')
# )
# ggsave(filename = 'posterior_conditional_alldyads.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = posterior_conditional.1, device = 'png', width = 1400, height = 800, units = 'px')
# 
# set.seed(15) ; plot_dyads <- sample(1:nrow(counts), size = 100, replace = F)
# plot_dyad_ids <- unique(edges$dyad)[plot_dyads]
# edges_subset <- edges[edges$dyad %in% plot_dyad_ids,]
# edges_subset <- edges_subset[edges_subset$chain == 'chain1',]
# (posterior_conditional.2 <- ggplot(data = edges_subset)+
#     geom_density(aes(group = dyad, x = edge_draw, colour = together), show.legend = FALSE)+
#     stat_density(aes(group = dyad, x = edge_draw, colour = together),
#                  geom = "line", position = "identity", linewidth = 0)+
#     scale_x_continuous(name = 'edge weight', limits = c(-0.05, 1.05), expand = c(0,0))+
#     scale_y_continuous(name = 'density', limits = c(-2, 72), expand = c(0,0))+
#     scale_colour_manual(values = c(rgb(33/255, 145/255, 140/255, 0.5), rgb(68/255, 1/255, 84/255, 0.5)),
#                         aesthetics = 'colour')+
#     theme(legend.position = c(0.5,0.75),
#           legend.background = element_rect(fill = 'white', colour = 'black'),
#           legend.key.height = unit(4, 'mm'),
#           legend.title = element_text(size = 10),
#           legend.text = element_text(size = 8))+
#     guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))+
#     labs(colour = 'sightings together')
# )
# ggsave(filename = 'posterior_conditional_sampledyads.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = posterior_conditional.2, device = 'png', width = 1400, height = 800, units = 'px')
# 
# ## conditional, mean edge weight vs dyad sightings ####
# averages_conditional <- edges %>% 
#   group_by(dyad) %>% 
#   mutate(mean = mean(edge_draw),
#          median = median(edge_draw)) %>% 
#   select(dyad, mean, median) %>% 
#   distinct() %>% 
#   rename(dyad_id = dyad) %>% 
#   left_join(counts_df[,c('dyad_id','node_1_males','node_2_males','event_count','count_dyad','count_1','count_2')],
#             by = 'dyad_id')
# 
# (edgesightings_conditional.1 <- ggplot(averages_conditional,
#                                        aes(x = count_dyad, y = median))+
#     geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
#     scale_x_continuous(name = 'total dyad sightings')+
#     scale_y_continuous(name = 'mean weight', limits = c(-0.02,1.02), expand = c(0,0))
# )
# ggsave(filename = 'edgesightings_conditional_noline.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = edgesightings_conditional.1, device = 'png', width = 1400, height = 800, units = 'px')
# 
# (edgesightings_conditional.2 <- ggplot(averages_conditional,
#                                        aes(x = count_dyad, y = median))+
#     geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
#     geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
#     scale_x_continuous(name = 'total dyad sightings')+
#     scale_y_continuous(name = 'mean weight', limits = c(-0.02,1.02), expand = c(0,0))
# )
# ggsave(filename = 'edgesightings_conditional_withline.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = edgesightings_conditional.2, device = 'png', width = 1400, height = 800, units = 'px')
# 
# (edgesightings_conditional.3 <- ggplot()+
#     geom_point(data = edges, aes(x = count_dyad, y = edge_draw),
#                colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
#     geom_point(data = averages_conditional,
#                aes(x = count_dyad, y = median),
#                colour = rgb(33/255, 145/255,140/255, 0.1), size = 0.5, shape = 19)+
#     geom_smooth(data = averages_conditional,
#                 aes(x = count_dyad, y = median),
#                 colour = rgb(68/255, 1/255, 84/255))+
#     scale_x_continuous(name = 'total dyad sightings')+
#     scale_y_continuous(name = 'edge weight', limits = c(-0.02,1.02), expand = c(0,0))
# )
# ggsave(filename = 'edgesightings_conditional_allpoints_withline.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = edgesightings_conditional.3, device = 'png', width = 1400, height = 800, units = 'px')
# 
# ## save outputs
# save.image('../outputs/sparse_network_methods_figures/plots_conditional.RData')
# 
# ## merge ####
# #load('../outputs/sparse_network_methods_figures/plots_conditional.RData')
# rm(list = ls()[ !ls() %in% c('posterior_conditional.2','edgesightings_conditional.3')]) ; gc()
# 
# # run patchwork and save
# (posterior_conditional.2 + edgesightings_conditional.3) +
#   plot_annotation(tag_levels = 'a')
# ggsave(filename = 'outputs_conditional_narrow.png',
#        path = '../outputs/sparse_network_methods_figures/',
#        plot = last_plot(), device = 'png', width = 1600, height = 800, units = 'px')
# 
# # clean up
# rm(list = ls()) ; gc()
# 
# #### Supplementary material: simulation -- 6 facet plot of 10 sightings vs 2 sightings, never together vs half together vs always together ####
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

(figure1 <- ggplot()+
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
                       breaks=c('True edge weight', 'SRI probability', 'Bayesian posterior distribution'),
                       values=c('True edge weight'='#fde725',
                                'SRI probability'='#440154',
                                'Bayesian posterior distribution'='#21918c'))+
    theme(legend.position = 'bottom')
)
ggsave(filename = 'figure1_example_sri_vs_bison.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure1, device = 'png', width = 2100, height = 1600, units = 'px')

# clean up
rm(list = ls()) ; gc()

#### Supplementary material: eigenvector vs individual sighting count and total dyads together = 0 ####
## SRI, eigenvector centrality vs individual sighting count ####
load('../outputs/sparse_network_methods_figures/model_run_uniform.RData')

# calculate SRI and identify individuals
counts_df$sri <- counts_df$event_count / counts_df$count_dyad
elephants <- unique(c(counts_df$id_1, counts_df$id_2))

# calculate adjacency matrix
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

# create network
sri_net <- igraph::graph_from_adjacency_matrix(sri_adjmat, weighted = T, mode = 'undirected', diag = FALSE)

# calculate eigenvector centrality
eigen_sri <- igraph::eigen_centrality(sri_net)$vector %>%
  as.data.frame()
colnames(eigen_sri) <- 'eigen'
eigen_sri$id <- rownames(eigen_sri)

# create data frame showing number of sightings per individual with other elephants
together0 <- counts_df %>% filter(event_count == 0)
eigen_sri$together0 <- NA ; for(i in 1:nrow(eigen_sri)) {
  x <- together0 %>%
    filter(id_1 == eigen_sri$id[i] | id_2 == eigen_sri$id[i])
  eigen_sri$together0[i] <- nrow(x)
}

# calculate number of sightings per elephant
counts <- data.frame(id = elephants, count = NA)
for(i in 1:nrow(counts)){
  counts$count[i] <- ifelse(i == nrow(counts),
                            unique(counts_df$count_2[which(counts_df$id_2 == counts$id[i])]),
                            unique(counts_df$count_1[which(counts_df$id_1 == counts$id[i])]))
}
eigen_sri <- eigen_sri %>%
  left_join(counts, by = 'id')

# plot
(eigensightings_sri.1 <- ggplot(eigen_sri, aes(x = count, y = eigen))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigensightings_sri_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_sri.1, device = 'png', width = 1400, height = 800, units = 'px')

(eigensightings_sri.2 <- ggplot(eigen_sri, aes(x = count, y = eigen))+
    geom_point(colour = rgb(33/255, 145/255, 140/255),
               size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigensightings_sri_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_sri.2, device = 'png', width = 1400, height = 800, units = 'px')

## SRI, eigenvector centrality vs total number of dyads where together = 0 ####
(eigen0_sri.1 <- ggplot(eigen_sri, aes(x = together0, y = eigen))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigentogether0_sri_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_sri.1, device = 'png', width = 1400, height = 800, units = 'px')

(eigen0_sri.2 <- ggplot(eigen_sri, aes(x = together0, y = eigen))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigentogether0_sri_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_sri.2, device = 'png', width = 1400, height = 800, units = 'px')

## change line breaks
(eigen0_sri.2 <- ggplot(eigen_sri, aes(x = together0, y = eigen))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0',
                       breaks = c(100,150,200))+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

# save workspace
save.image('../outputs/sparse_network_methods_figures/eigen_checks.RData')

## uniform, eigenvector centrality vs individual sighting count ####
uniform_adjarr <- array(NA, dim = c(length(elephants), length(elephants), 1000),
                       dimnames = list(elephants, elephants, NULL))
N <- nrow(counts_df)

for (i in 1:N) {
  dyad_row <- counts_df[i, ]
  uniform_adjarr[dyad_row$id_1, dyad_row$id_2, ] <- edges[1:1000, i]
}
#uniform_adjarr[,,1]

# make matrix to save eigenvector values to
eigen_uniform <- matrix(NA, nrow = 1000, ncol = length(elephants),
                       dimnames = list(1:1000, elephants))

# calculate eigenvector scores
for(i in 1:nrow(eigen_uniform)){
  uniform_net <- igraph::graph_from_adjacency_matrix(uniform_adjarr[,,i], weighted = T, mode = 'undirected', diag = FALSE)
  eigen <- as.data.frame(igraph::eigen_centrality(uniform_net)$vector)
  eigen_uniform[i,] <- eigen$`igraph::eigen_centrality(uniform_net)$vector`
}
#eigen_uniform[,1]
#eigen_uniform[1,]

# make into long format data frame
eigen_uniform <- eigen_uniform %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'eigen') %>%
  group_by(id) %>%
  mutate(mean = mean(eigen),
         median = median(eigen)) %>%
  left_join(eigen_sri[,c('id', 'together0', 'count')], by = 'id')

# create data frame of average eigenvector centrality scores
averages_uniform <- eigen_uniform %>%
  select(id, mean, median, together0, count) %>%
  distinct()

# plot
(eigensightings_uniform.1 <- ggplot()+
    geom_point(data = eigen_uniform, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = averages_uniform, aes(x = count, y = median),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigensightings_uniform_allpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_uniform.1, device = 'png', width = 1400, height = 800, units = 'px')

(eigensightings_uniform.2 <- ggplot()+
    geom_point(data = eigen_uniform, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = averages_uniform, aes(x = count, y = median),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(data = eigen_uniform, aes(x = count, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigensightings_uniform_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_uniform.2, device = 'png', width = 1400, height = 800, units = 'px')

(eigensightings_uniform.3 <- ggplot(data = averages_uniform, aes(y = median,
                                        #y = standardised,
                                        x = count))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean eigenvector centrality',
                       limits = c(0,1),
                       expand = c(0.02,0.02))
)
ggsave(filename = 'eigensightings_uniform_medianpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_uniform.3, device = 'png', width = 1400, height = 800, units = 'px')

(eigensightings_uniform.4 <- ggplot(data = averages_uniform, aes(x = count, y = median))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean eigenvector centrality',
                       limits = c(0,1),
                       expand = c(0,0))
)
ggsave(filename = 'eigensightings_uniform_medianpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_uniform.4, device = 'png', width = 1400, height = 800, units = 'px')

## uniform, eigenvector centrality vs total number of dyads where together = 0 ####
(eigen0_uniform.1 <- ggplot()+
    geom_point(data = eigen_uniform, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = averages_uniform, aes(x = together0, y = median),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigentogether0_uniform_allpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_uniform.1, device = 'png', width = 1400, height = 800, units = 'px')

(eigen0_uniform.2 <- ggplot()+
    geom_point(data = eigen_uniform, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = averages_uniform, aes(x = together0, y = median),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(data = eigen_uniform, aes(x = together0, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigentogether0_uniform_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_uniform.2, device = 'png', width = 1400, height = 800, units = 'px')

## alter axis labels
(eigen0_uniform.2 <- ggplot()+
    geom_point(data = eigen_uniform, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = averages_uniform, aes(x = together0, y = median),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(data = eigen_uniform, aes(x = together0, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0',
                       breaks = c(100,150,200))+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

(eigen0_uniform.3 <- ggplot(averages_uniform)+
    geom_point(aes(x = together0, y = median, colour = count),
               size = 0.5, shape = 19)+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))+
    labs(colour = 'number of\nsightings')
)
ggsave(filename = 'eigentogether0_uniform_medianpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_uniform.3, device = 'png', width = 1400, height = 800, units = 'px')

(eigen0_uniform.4 <- ggplot(averages_uniform)+
    geom_point(aes(x = together0, y = median, colour = count),
               size = 0.5, shape = 19)+
    geom_smooth(aes(x = together0, y = median),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))+
    labs(colour = 'number of\nsightings')
)
ggsave(filename = 'eigentogether0_uniform_medianpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_uniform.4, device = 'png', width = 1400, height = 800, units = 'px')

# save workspace
save.image('../outputs/sparse_network_methods_figures/eigen_checks.RData')

## default, eigenvector centrality vs individual sighting count ####
load('../outputs/sparse_network_methods_figures/model_run_default.RData')

# get adjacency array for default prior
default_adjarr <- array(NA, dim = c(length(elephants), length(elephants), 1000),
                        dimnames = list(elephants, elephants, NULL))
for (i in 1:N) {
  dyad_row <- counts_df[i, ]
  default_adjarr[dyad_row$id_1, dyad_row$id_2, ] <- edges[1:1000, i]
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
#eigen_default[,1]
#eigen_default[1,]

# make long format data frame
eigen_default <- eigen_default %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'eigen') %>%
  group_by(id) %>%
  mutate(mean = mean(eigen),
         median = median(eigen)) %>%
  left_join(eigen_sri[,c('id', 'together0', 'count')], by = 'id')

# calculate averages
averages_default <- eigen_default %>%
  select(id, mean, median, together0, count) %>%
  distinct()

# plot
(eigensightings_default.1 <- ggplot()+
    geom_point(data = eigen_default, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5,
               shape = 19)+
    geom_point(data = averages_default, aes(x = count, y = median),
               colour = rgb(33/255, 145/255, 140/255),
               size = 0.5,
               shape = 19)+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
  )
ggsave(filename = 'eigensightings_default_allpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_default.1, device = 'png', width = 1400, height = 800, units = 'px')

(eigensightings_default.2 <- ggplot()+
    geom_point(data = eigen_default, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5,
               shape = 19)+
    geom_point(data = averages_default, aes(x = count, y = median),
               colour = rgb(33/255, 145/255, 140/255),
               size = 0.5,
               shape = 19)+
    geom_smooth(data = eigen_default, aes(x = count, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigensightings_default_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_default.2, device = 'png', width = 1400, height = 800, units = 'px')

(eigensightings_default.3 <- ggplot(data = averages_default, aes(y = median,
                                        x = count))+
    geom_point(colour = rgb(33/255, 145/255, 140/255),
               size = 0.5,
               shape = 19)+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'median eigenvector centrality',
                       limits = c(0, 1),
                       expand = c(0.02,0.02))
)
ggsave(filename = 'eigensightings_default_medianpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_default.3, device = 'png', width = 1400, height = 800, units = 'px')

(eigensightings_default.4 <- ggplot(data = averages_default, aes(x = count, y = median))+
    geom_point(colour = rgb(33/255, 145/255, 140/255),
               size = 0.5,
               shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'median eigenvector centrality',
                       limits = c(0,1),
                       expand = c(0.02,0.02))
)
ggsave(filename = 'eigensightings_default_medianpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_default.4, device = 'png', width = 1400, height = 800, units = 'px')

## default, eigenvector centrality vs total number of dyads where together = 0 ####
(eigen0_default.1 <- ggplot()+
    geom_point(data = eigen_default, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5,
               shape = 19)+
    geom_point(data = averages_default, aes(x = together0, y = median, colour = count),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigentogether0_default_allpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_default.1, device = 'png', width = 1400, height = 800, units = 'px')

(eigen0_default.2 <- ggplot()+
    geom_point(data = eigen_default, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5,
               shape = 19)+
    geom_point(data = averages_default, aes(x = together0, y = median),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(data = eigen_default, aes(x = together0, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigentogether0_default_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_default.2, device = 'png', width = 1400, height = 800, units = 'px')

## alter axis labels
(eigen0_default.2 <- ggplot()+
    geom_point(data = eigen_default, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5,
               shape = 19)+
    geom_point(data = averages_default, aes(x = together0, y = median),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(data = eigen_default, aes(x = together0, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0',
                       breaks = c(100,150,200))+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

(eigen0_default.3 <- ggplot(averages_default)+
    geom_point(aes(x = together0, y = median, colour = count),
               size = 0.5,
               shape = 19)+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))+
    labs(colour = 'number of\nsightings')
)
ggsave(filename = 'eigentogether0_default_medianpoints_noline_coloursightings.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_default.3, device = 'png', width = 1400, height = 800, units = 'px')

(eigen0_default.4 <- ggplot(averages_default)+
    geom_point(aes(x = together0, y = median, colour = count),
               size = 0.5,
               shape = 19)+
    geom_smooth(aes(x = together0, y = median),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))+
    labs(colour = 'number of\nsightings')
)
ggsave(filename = 'eigentogether0_default_medianpoints_withline_colourssightings.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_default.4, device = 'png', width = 1400, height = 800, units = 'px')

# save workspace
save.image('../outputs/sparse_network_methods_figures/eigen_checks.RData')

## skewed, eigenvector centrality vs individual sighting count ####
load('../outputs/sparse_network_methods_figures/model_run_skewed.RData')

# get adjacency array for skewed prior
skewed_adjarr <- array(NA, dim = c(length(elephants), length(elephants), 1000),
                       dimnames = list(elephants, elephants, NULL))
for (i in 1:N) {
  dyad_row <- counts_df[i, ]
  skewed_adjarr[dyad_row$id_1, dyad_row$id_2, ] <- edges[1:1000, i]
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
#eigen_skewed[,1]
#eigen_skewed[1,]

# make long format data frame
eigen_skewed <- eigen_skewed %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'eigen') %>%
  group_by(id) %>%
  mutate(mean = mean(eigen),
         median = median(eigen)) %>%
  left_join(eigen_sri[,c('id', 'together0', 'count')], by = 'id')

# calculate averages
averages_skewed <- eigen_skewed %>%
  select(id, mean, median, together0, count) %>%
  distinct()

# plot
(eigensightings_skewed.1 <- ggplot()+
    geom_point(data = eigen_skewed, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = averages_skewed, aes(x = count, y = median),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigensightings_skewed_allpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_skewed.1, device = 'png', width = 1400, height = 800, units = 'px')

(eigensightings_skewed.2 <- ggplot()+
    geom_point(data = eigen_skewed, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = averages_skewed, aes(x = count, y = median),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(data = eigen_skewed, aes(x = count, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigensightings_skewed_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_skewed.2, device = 'png', width = 1400, height = 800, units = 'px')

(eigensightings_skewed.3 <- ggplot(data = averages_skewed,
                      aes(y = median,
                          #y = standardised,
                          x = count))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean eigenvector centrality',
                       limits = c(0,1),
                       expand = c(0.02,0.02))
)
ggsave(filename = 'eigensightings_skewed_medianpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_skewed.3, device = 'png', width = 1400, height = 800, units = 'px')

(eigensightings_skewed.4 <- ggplot(data = averages_skewed, aes(x = count, y = median))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean eigenvector centrality',
                       limits = c(0,1),
                       expand = c(0,0))
)
ggsave(filename = 'eigensightings_skewed_medianpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_skewed.4, device = 'png', width = 1400, height = 800, units = 'px')

## skewed, eigenvector centrality vs total number of dyads where together = 0 ####
(eigen0_skewed.1 <- ggplot()+
    geom_point(data = eigen_skewed, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = averages_skewed, aes(x = together0, y = median),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigentogether0_skewed_allpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_skewed.1, device = 'png', width = 1400, height = 800, units = 'px')

(eigen0_skewed.2 <- ggplot()+
    geom_point(data = eigen_skewed, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = averages_skewed, aes(x = together0, y = median),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(data = eigen_skewed, aes(x = together0, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigentogether0_skewed_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_skewed.2, device = 'png', width = 1400, height = 800, units = 'px')

## alter axis breaks
(eigen0_skewed.2 <- ggplot()+
    geom_point(data = eigen_skewed, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = averages_skewed, aes(x = together0, y = median),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(data = eigen_skewed, aes(x = together0, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0',
                       breaks = c(100,150,200))+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

(eigen0_skewed.3 <- ggplot(averages_skewed)+
    geom_point(aes(x = together0, y = median, colour = count),
               size = 0.5, shape = 19)+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))+
    labs(colour = 'number of\nsightings')
)
ggsave(filename = 'eigentogether0_skewed_medianpoints_noline_coloursightings.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_skewed.3, device = 'png', width = 1400, height = 800, units = 'px')

(eigen0_skewed.4 <- ggplot(averages_skewed)+
    geom_point(aes(x = together0, y = median, colour = count),
               size = 0.5, shape = 19)+
    geom_smooth(aes(x = together0, y = median),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))+
    labs(colour = 'number of\nsightings')
)
ggsave(filename = 'eigentogether0_skewed_medianpoints_withline_colourssightings.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_skewed.4, device = 'png', width = 1400, height = 800, units = 'px')

# save workspace
save.image('../outputs/sparse_network_methods_figures/eigen_checks.RData')

## conditional, eigenvector vs total number of dyads where together = 0 ####
load('motnp_edgeweights_conditionalprior.RData')

# get edge samples
edge_samples <- fit_edges_motnp$draws() %>% as.data.frame()
edge_samples <- edge_samples[,(n_chains+1):ncol(edge_samples)]
edge_samples <- edge_samples[,seq(1,ncol(edge_samples)-3, by = 4)]
colnames(edge_samples) <- counts_df$dyad_id

# get adjacency array for default prior
conditional_adjarr <- array(NA, dim = c(length(elephants), length(elephants), 1000),
                      dimnames = list(elephants, elephants, NULL))
for (i in 1:N) {
  dyad_row <- counts_df[i, ]
  conditional_adjarr[dyad_row$id_1, dyad_row$id_2, ] <- edge_samples[, i]
  conditional_adjarr[dyad_row$id_2, dyad_row$id_1, ] <- edge_samples[, i]
}
conditional_adjarr[,,1]

# make matrix to save eigenvector values to
eigen_conditional <- matrix(NA, nrow = 1000, ncol = length(elephants),
                      dimnames = list(1:1000, elephants))

# calculate eigenvector scores
for(i in 1:nrow(eigen_conditional)){
  conditional_net <- igraph::graph_from_adjacency_matrix(adjmatrix = conditional_adjarr[,,i], diag = FALSE, mode = 'undirected', weighted = T)
  eigen <- as.data.frame(igraph::eigen_centrality(conditional_net)$vector)
  eigen_conditional[i,] <- eigen$`igraph::eigen_centrality(conditional_net)$vector`
}
#eigen_conditional[,1]
#eigen_conditional[1,]

# make long format data frame
eigen_conditional <- eigen_conditional %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'eigen') %>%
  group_by(id) %>%
  mutate(mean = mean(eigen),
         median = median(eigen)) %>%
  left_join(eigen_sri[,c('id', 'together0', 'count')], by = 'id')

# calculate averages
averages_conditional <- eigen_conditional %>%
  select(id, mean, median, count, together0) %>%
  distinct()

# plot
(eigen0_conditional.1 <- ggplot()+
    geom_point(data = eigen_conditional, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5)+
    geom_point(data = averages_conditional, aes(x = together0, y = median),
               colour = rgb(33/255, 145/255, 140/255),
               size = 0.5)+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigentogether_conditional_allpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_conditional.1, device = 'png', width = 1400, height = 800, units = 'px')

(eigen0_conditional.2 <- ggplot()+
    geom_point(data = eigen_conditional, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5)+
    geom_point(data = averages_conditional, aes(x = together0, y = median),
               colour = rgb(33/255, 145/255, 140/255),
               size = 0.5)+
    geom_smooth(data = eigen_conditional, aes(x = together0, y = eigen),
                colour = rgb(68/255,1/255,84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigentogether_conditional_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_conditional.2, device = 'png', width = 1400, height = 800, units = 'px')

## alter axis labels
(eigen0_conditional.2 <- ggplot()+
    geom_point(data = eigen_conditional, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5)+
    geom_point(data = averages_conditional, aes(x = together0, y = median),
               colour = rgb(33/255, 145/255, 140/255),
               size = 0.5)+
    geom_smooth(data = eigen_conditional, aes(x = together0, y = eigen),
                colour = rgb(68/255,1/255,84/255))+
    scale_x_continuous(name = 'dyads together = 0',
                       breaks = c(100,150,200))+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

(eigen0_conditional.3 <- ggplot(averages_conditional)+
    geom_point(aes(x = together0, y = median, colour = count),
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
ggsave(filename = 'eigentogether_conditional_medianpoints_noline_coloursightings.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_conditional.3, device = 'png', width = 1400, height = 800, units = 'px')

(eigen0_conditional.4 <- ggplot(averages_conditional)+
    geom_point(aes(x = together0, y = median, colour = count),
               size = 0.5, shape = 19)+
    geom_smooth(aes(x = together0, y = median),
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
ggsave(filename = 'eigentogether_conditional_medianpoints_withline_coloursightings.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_conditional.4, device = 'png', width = 1400, height = 800, units = 'px')

## conditional, eigenvector vs individual sightings ####
(eigensightings_conditional.1 <- ggplot()+
   geom_point(data = eigen_conditional, aes(x = count, y = eigen),
              colour = rgb(253/255, 231/255, 37/255, 0.01))+
   geom_point(data = averages_conditional, aes(x = count, y = median),
              colour = rgb(33/255, 145/255, 140/255))+
   scale_x_continuous(name = 'total node sightings')+
   scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigensightings_conditional_allpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_conditional.1, device = 'png', width = 1400, height = 800, units = 'px')

(eigensightings_conditional.2 <- ggplot()+
    geom_point(data = eigen_conditional, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01))+
    geom_point(data = averages_conditional, aes(x = count, y = median),
               colour = rgb(33/255, 145/255, 140/255))+
    geom_smooth(data = eigen_conditional, aes(x = count, y = eigen),
                colour = rgb(68/255,1/255,84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigensightings_conditional_allpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_conditional.2, device = 'png', width = 1400, height = 800, units = 'px')
(eigensightings_conditional.2 <- ggplot()+
    geom_point(data = eigen_conditional, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5)+
    geom_point(data = averages_conditional, aes(x = count, y = median),
               colour = rgb(33/255, 145/255, 140/255),
               size = 0.5)+
    geom_smooth(data = eigen_conditional, aes(x = count, y = eigen),
                colour = rgb(68/255,1/255,84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

(eigensightings_conditional.3 <- ggplot(data = averages_conditional,
                                        aes(y = median, x = count))+
    geom_point(colour = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean centrality',
                       limits = c(0, 1),
                       expand = c(0.02,0.02))
)
ggsave(filename = 'eigensightings_conditional_medianpoints_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_conditional.3, device = 'png', width = 1400, height = 800, units = 'px')

(eigensightings_conditional.4 <- ggplot(data = averages_conditional,
                                        aes(x = count, y = median))+
    geom_point(colour = rgb(33/255, 145/255, 140/255))+
    geom_smooth(colour = rgb(68/255,1/255,84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean centrality',
                       limits = c(0,1),
                       expand = c(0.02,0.02))
)
ggsave(filename = 'eigensightings_conditional_medianpoints_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_conditional.4, device = 'png', width = 1400, height = 800, units = 'px')

# save workspace
save.image('../outputs/sparse_network_methods_figures/eigen_checks.RData')

## merge ####
#load('../outputs/sparse_network_methods_figures/eigen_checks.RData')
rm(list = ls()[!ls() %in% c('eigen0_sri.2','eigensightings_sri.2',
                            'eigen0_uniform.2','eigen0_default.2','eigen0_skewed.2',
                            'eigensightings_uniform.2','eigensightings_default.2','eigensightings_skewed.2',
                            'eigen0_conditional.2','eigensightings_conditional.2')]) ; gc()

# SRI
(eigen0_sri.2 + eigensightings_sri.2)+
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'eigen_outputs_sri.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 1600, height = 800, units = 'px')
rm(eigen0_sri.2,eigensightings_sri.2) ; gc()

# unconditional
(eigen0_uniform.2 + eigen0_default.2 + eigen0_skewed.2) / 
  (eigensightings_uniform.2 + eigensightings_default.2 + eigensightings_skewed.2)+
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'eigen_outputs_unconditional.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2400, height = 1600, units = 'px')
rm(eigen0_uniform.2,eigen0_default.2,eigen0_skewed.2,
   eigensightings_uniform.2,eigensightings_default.2,eigensightings_skewed.2) ; gc()

# conditional
(eigen0_conditional.2 + eigensightings_conditional.2)+
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'eigen_outputs_conditional.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 1600, height = 800, units = 'px')
