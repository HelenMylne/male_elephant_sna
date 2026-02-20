######## information ########
# run Bayesian mixture model using first simulated then empirical data

######## set up ########
# library(rstan) ; library(StanHeaders)
# library(tidyverse) ; library(LaplacesDemon) ; library(patchwork) ; library(igraph)
# library(cmdstanr)
# library(cmdstanr, lib.loc = '../packages/')
library(StanHeaders, lib.loc = '../packages/')
library(rstan, lib.loc = '../packages/')
library(tidyverse, lib.loc = '../packages/')
library(LaplacesDemon, lib.loc = '../packages/')
library(patchwork, lib.loc = '../packages/')
library(igraph, lib.loc = '../packages/')

## set stan parameters
set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0/')
n_chains <- 4
n_draws <- 1000

## define plot theme
theme_set(theme_bw(base_size = 12))

#### load model ####
# mod <- cmdstan_model('methods_paper/models/edge_mixture_bayesian3.stan')
mod <- stan_model('methods_paper/models/edge_mixture_bayesian3.stan')

# ######## simulation ########
# #### load data ####
# load('methods_paper/simulation_outputs/create_simulated_data.RData')
# sim_dyads <- sim_dyads %>% 
#   mutate(ever = ifelse(together == 0, 0, 1))
# 
# #### create model data list ####
# sim_data <- list(n_dyads = n_dyads,
#                  ever_together = sim_dyads$ever,
#                  count_together = sim_dyads$together,
#                  total_sightings = sim_dyads$total_sightings,
#                  dyad_id = sim_dyads$dyad_id)
# 
# #### run model ####
# fit <- mod$sample(data = sim_data, seed = 12345,
#                   chains = n_chains, parallel_chains = n_chains,
#                   iter_warmup = n_draws, iter_sampling = n_draws)
# 
# save.image('methods_paper/simulation_outputs/fit_bayesian_mixture3.RData')
# 
# ## summary
# fit$summary()
# 
# #### extract edges ####
# ## extract
# draws <- as.data.frame(fit$draws())
# 
# ## put all chains into single column: separate out chains
# draws1 <- draws[,seq(1, ncol(draws), by = 4)]
# draws2 <- draws[,seq(2, ncol(draws), by = 4)]
# draws3 <- draws[,seq(3, ncol(draws), by = 4)]
# draws4 <- draws[,seq(4, ncol(draws), by = 4)]
# 
# ## bind 4 chains into single data frame
# colnames(draws2) <- colnames(draws3) <- colnames(draws4) <- colnames(draws1)
# draws <- rbind(draws1, draws2, draws3, draws4)
# rm(draws1, draws2, draws3, draws4) ; gc()
# 
# ## rename without "1." on start
# params <- data.frame(param = colnames(draws)) %>% 
#   mutate(param = str_remove(pattern = '1.', string = param))
# colnames(draws) <- params$param
# 
# ## add columns for draw/chain
# draws <- draws %>% 
#   mutate(draw = rep(1:n_draws, n_chains),
#          chain = rep(1:n_chains, each = n_draws))
# 
# ## check know where all parameters are coming from
# global <- params[params$param %in% c('lp__','logit_together','average_edge','prob_together'),]
# dyad_effects <- params[grep(x = params$param, pattern = 'dyad_effect'),]
# hypothetical <- params[grep(x = params$param, pattern = 'hypothetical_together'),]
# logit_edges <- params[grep(x = params$param, pattern = 'logit_edge'),]
# edge_weights <- params[grep(x = params$param, pattern = 'edge_weight'),]
# nrow(params) == length(global) + length(dyad_effects) + length(hypothetical) + length(logit_edges) + length(edge_weights)
# global <- global[2:4]
# 
# #### check fit ####
# ## traceplot global parameters
# draws %>% 
#   select(all_of(global), draw, chain) %>% 
#   pivot_longer(cols = all_of(global),
#                names_to = 'param', values_to = 'estimate') %>% 
#   ggplot()+
#   geom_line(aes(x = draw, y = estimate, colour = as.factor(chain)))+
#   facet_wrap(param ~ ., scales = 'free')
# 
# ## sample dyads for plotting dyad-level parameters
# plot_dyads <- sample(1:n_dyads, 12, F)
# 
# ## dyad_effects
# d <- draws %>% 
#   select(all_of(dyad_effects), draw, chain) %>% 
#   select(all_of(plot_dyads), draw, chain)
# d %>% 
#   pivot_longer(cols = 1:(ncol(d)-2),
#                names_to = 'param', values_to = 'estimate') %>% 
#   ggplot()+
#   geom_line(aes(x = draw, y = estimate, colour = as.factor(chain)))+
#   facet_wrap(param ~ ., scales = 'free')
# 
# ## hypothetical
# d <- draws %>% 
#   select(all_of(hypothetical), draw, chain) %>% 
#   select(all_of(plot_dyads), draw, chain)
# d %>% 
#   pivot_longer(cols = 1:(ncol(d)-2),
#                names_to = 'param', values_to = 'estimate') %>% 
#   ggplot()+
#   geom_line(aes(x = draw, y = estimate, colour = as.factor(chain)))+
#   facet_wrap(param ~ ., scales = 'free')
# 
# ## logit_edges
# d <- draws %>% 
#   select(all_of(logit_edges), draw, chain) %>% 
#   select(all_of(plot_dyads), draw, chain)
# d %>% 
#   pivot_longer(cols = 1:(ncol(d)-2),
#                names_to = 'param', values_to = 'estimate') %>% 
#   ggplot()+
#   geom_line(aes(x = draw, y = estimate, colour = as.factor(chain)))+
#   facet_wrap(param ~ ., scales = 'free')
# 
# ## edge_weights
# d <- draws %>% 
#   select(all_of(edge_weights), draw, chain) %>% 
#   select(all_of(plot_dyads), draw, chain)
# d %>% 
#   pivot_longer(cols = 1:(ncol(d)-2),
#                names_to = 'param', values_to = 'estimate') %>% 
#   ggplot()+
#   geom_line(aes(x = draw, y = estimate, colour = as.factor(chain)))+
#   facet_wrap(param ~ ., scales = 'free')
# rm(global,dyad_effects,hypothetical,logit_edges,d) ; gc()
# 
# #### plot edges -- NOT ACTUALLY SURE WHICH ONE EVEN IS THE EDGE WEIGHT IN THE END!! ####
# ## extract estimates
# ew <- draws %>% 
#   select(all_of(edge_weights))
# 
# # ## add probability together??
# # for(i in 1:n_dyads){
# #   ew[,i] <- ew[,i] + draws$prob_together
# # }
# 
# ## calculate average values
# sim_dyads <- sim_dyads %>% 
#   mutate(mean_edge = apply(ew, 2, mean),
#          mid_edge = apply(ew, 2, median),
#          stdv_edge = apply(ew, 2, sd))
# 
# ## plot average estimate against true value
# ggplot(sim_dyads)+
#   geom_point(aes(x = true_edge, y = mean_edge, colour = as.factor(ever)))+
#   geom_abline(slope = 1, intercept = 0)+
#   labs(x = 'true (simulated) dyadic edge weight',
#        y = 'mean edge estimate',
#        colour = 'ever seen\ntogether')
# ggplot(sim_dyads)+
#   geom_point(aes(x = true_edge, y = mid_edge, colour = as.factor(ever)))+
#   geom_abline(slope = 1, intercept = 0)+
#   labs(x = 'true (simulated) dyadic edge weight',
#        y = 'mean edge estimate',
#        colour = 'ever seen\ntogether')
# 
# ## plot uncertainty (sd) against sightings
# ggplot(sim_dyads)+
#   geom_point(aes(x = total_sightings, y = stdv_edge, colour = as.factor(ever)))+
#   labs(x = 'total sightings of dyad',
#        y = 'uncertainty (SD) in edge estimate',
#        colour = 'ever seen\ntogether')
# 
# 
######## empirical data ########
#### create data list -- ONCE IT'S ALL CHECKED, COME BACK AND REMOVE THE FILTERING DOWN TO RANDOM SET OF ELEPHANTS ####
load('motnp_edgeweights_conditionalprior.RData')
rm(list = ls()[! ls() %in% c('counts_df','nodes','n_chains','n_draws','n_dyads','mod')]) ; gc()

## create data list
counts_dl <- list(n_dyads = n_dyads,
                  count_together = counts_df$event_count,
                  total_sightings = counts_df$count_dyad)

#### fit model ####
# fit <- mod$sample(data = counts_dl, seed = 12345,
#                   chains = n_chains, parallel_chains = n_chains,
#                   iter_warmup = n_draws, iter_sampling = n_draws)
fit <- rstan::sampling(mod, data = counts_dl, seed = 12345,
                  chains = n_chains, cores = n_chains,
                  iter = n_draws, warmup = n_draws/2)

## save output
save.image('methods_paper/motnp_fit_bayesian_mixture.RData') # load('methods_paper/motnp_fit_bayesian_mixture.RData')

## summary
# fit$summary()
fit

#### extract draws ####
## extract
draws <- as.data.frame(rstan::extract(fit))
# draws <- as.data.frame(fit$draws())
# 
# ## put all chains into single column: separate out chains
# draws1 <- draws[,seq(1, ncol(draws), by = 4)]
# draws2 <- draws[,seq(2, ncol(draws), by = 4)]
# draws3 <- draws[,seq(3, ncol(draws), by = 4)]
# draws4 <- draws[,seq(4, ncol(draws), by = 4)]
# 
# ## bind 4 chains into single data frame
# colnames(draws2) <- colnames(draws3) <- colnames(draws4) <- colnames(draws1)
# draws <- rbind(draws1, draws2, draws3, draws4)
# rm(draws1, draws2, draws3, draws4) ; gc()

## rename without "1." on start
params <- data.frame(param = colnames(draws)) #%>%
  # mutate(param = str_remove(pattern = '1.', string = param))
colnames(draws) <- params$param

## add columns for draw/chain
draws <- draws %>%
  # mutate(draw = rep(1:(n_draws), n_chains),
  #        chain = rep(1:n_chains, each = n_draws))
  mutate(draw = rep(1:(n_draws/2), n_chains),
         chain = rep(1:n_chains, each = n_draws/2))

#### check fit ####
## function traceplot
trace <- function(draws_dataframe, parameter_type, dyads = NULL){
  d <- draws_dataframe %>% 
    select(all_of(parameter_type), draw, chain)
  if(is.null(dyads) == FALSE){
    d <- d %>% 
      select(all_of(dyads), draw, chain)
  }
  d <- d %>% 
    pivot_longer(cols = 1:(ncol(d)-2),
                 names_to = 'param', values_to = 'estimate')
  g <- ggplot(d)+
    geom_line(aes(x = draw, y = estimate, colour = as.factor(chain)))+
    facet_wrap(param ~ ., scales = 'free')
  return(g)
}

## split into parameter types for checking
global <- params[params$param %in% c('lp__','average_edge','prob_0'),]
dyad_effects <- params[grep(x = params$param, pattern = 'dyad_effect'),]
edge_weights <- params[grep(x = params$param, pattern = 'edge_weight'),]
nrow(params) == length(global) + length(dyad_effects) + length(edge_weights)
global <- global[global != 'lp__']

## traceplot global parameters
trace(draws_dataframe = draws, parameter_type = global)

## sample dyads for plotting dyad-level parameters
plot_dyads <- sample(1:n_dyads, 12, F)

## dyad-level parameters
trace(draws, dyad_effects, plot_dyads)
trace(draws, edge_weights, plot_dyads)
rm(global,dyad_effects,logit_edges) ; gc()

#### plot edges ####
## extract estimates
ew <- draws %>% 
  select(all_of(edge_weights))
ew_mat <- as.matrix(ew)

## calculate average values
counts_df <- counts_df %>% 
  mutate(ever_together = ifelse(event_count == 0, 0, 1),
         together_chr = ifelse(event_count == 0, 'never together', 'together at least once')) %>% 
  mutate(sri = event_count / count_dyad,
         mean_edge = apply(ew_mat, 2, mean),
         stdv_edge = apply(ew_mat, 2, sd),
         mid_edge = apply(ew_mat, 2, median),
         lwr_edge = apply(ew_mat, 2, quantile, prob = 0.025),
         upr_edge = apply(ew_mat, 2, quantile, prob = 0.975))

## plot average estimate against SRI
bkgrnd <- data.frame(x = c(0, 0.25, 1),
                     y = c(0, 0.25, 0.25))
ggplot()+
  geom_area(data = bkgrnd,
            aes(x = x, y = y),
            fill = rgb(0,0,1,0.25))+
  geom_point(data = counts_df,
             aes(x = sri, y = mean_edge, colour = together_chr))+
  geom_abline(slope = 1, intercept = 0)+
  labs(x = 'SRI edge weight',
       y = 'mean edge estimate',
       colour = 'ever seen\ntogether')+
  annotate('text', x = 0.1, y = 0.23, label = 'zero-inflated higher')+
  annotate('text', x = 0.7, y = 0.23, label = 'SRI higher')+
  # annotate('text', x = 0.15, y = 0.21, label = 'zero-inflated higher', angle = 82.5)+
  # annotate('text', x = 0.26, y = 0.21, label = 'SRI higher', angle = 82.5)+
  scale_x_continuous(expand = c(0,0), limits = c(-0.01,1))+
  scale_y_continuous(expand = c(0,0), limits = c(0,0.25))+
  scale_colour_viridis_d()+
  theme(panel.background = element_rect(fill = rgb(0,0,1,0.1)))
ggplot()+
  geom_area(data = bkgrnd,
            aes(x = x, y = y),
            fill = rgb(0,0,1,0.25))+
  geom_point(data = counts_df,
             aes(x = sri, y = mid_edge, colour = together_chr))+
  geom_abline(slope = 1, intercept = 0)+
  labs(x = 'SRI edge weight',
       y = 'median edge estimate',
       colour = 'ever seen\ntogether')+
  annotate('text', x = 0.1, y = 0.23, label = 'zero-inflated higher')+
  annotate('text', x = 0.7, y = 0.23, label = 'SRI higher')+
  # annotate('text', x = 0.15, y = 0.21, label = 'zero-inflated higher', angle = 82.5)+
  # annotate('text', x = 0.26, y = 0.21, label = 'SRI higher', angle = 82.5)+
  scale_x_continuous(expand = c(0,0), limits = c(-0.01, 1))+
  scale_y_continuous(expand = c(0,0), limits = c(0,0.25))+
  scale_colour_viridis_d()+
  theme(panel.background = element_rect(fill = rgb(0,0,1,0.1)))

## plot uncertainty (sd) against sightings
ggplot(counts_df)+
  geom_point(aes(x = count_dyad, y = stdv_edge, colour = together_chr))+
  labs(x = 'total sightings of dyad',
       y = 'uncertainty (SD) in edge estimate',
       colour = 'ever seen\ntogether')

## plot edge weight against sightings
ggplot(counts_df)+
  geom_point(aes(x = count_dyad, y = mean_edge, colour = together_chr))+
  labs(x = 'total sightings of dyad',
       y = 'uncertainty (SD) in edge estimate',
       colour = 'ever seen\ntogether')

#### full model plots ####
ew_long <- ew 
colnames(ew_long) <- counts_df$dyad_id
plot_draws <- sample(1:1000, 50, FALSE)
plot_draws <- c(plot_draws, plot_draws+1000, plot_draws+2000, plot_draws+3000)
ew_long <- ew_long[plot_draws,] %>% 
  pivot_longer(cols = everything(), names_to = 'dyad_id', values_to = 'draw') %>% 
  mutate(dyad_id = as.numeric(dyad_id)) %>% 
  left_join(counts_df[,c('dyad_id','node_1','node_2',
                         'event_count','count_dyad','ever_together',
                         'age_category_1','age_category_2',
                         'sri','mean_edge','mid_edge','stdv_edge','lwr_edge','upr_edge')],
            by = 'dyad_id')

(edge_obs <- ggplot()+
    geom_point(data = ew_long,
               aes(x = count_dyad, y = draw),
               colour = '#fde725', alpha = 0.1)+
    geom_point(data = counts_df,
               aes(x = count_dyad, y = mean_edge),
               colour = '#5ec962')+
    geom_smooth(data = counts_df,
                aes(x = count_dyad, y = mean_edge, colour = together_chr))+
    scale_colour_viridis_d(begin = 0, end = 0.5)+
    labs(x = 'total dyad sightings',
         y = 'edge weight',
         colour = 'ever seen together')+
    scale_y_continuous(name = 'edge weight', limits = c(0,1), expand = c(0,0))+
    theme(legend.position = 'bottom'))
ggsave(filename = 'edgesightings_bayesmix.png',
       # path = '../outputs/sparse_network_methods_figures/',
       path = 'methods_paper/outputs/new_version/',
       plot = edge_obs, device = 'png',
       width = 1600, height = 700, units = 'px')

ew_plot <- ew 
colnames(ew_plot) <- counts_df$dyad_id
plot_dyads <- as.character(sample(counts_df$dyad_id, 100))
(edge_dist <- ew_plot %>% 
    select(all_of(plot_dyads)) %>% 
    pivot_longer(cols = everything(), names_to = 'dyad_id', values_to = 'draw') %>% 
    mutate(dyad_id = as.numeric(dyad_id)) %>% 
    left_join(counts_df[,c('dyad_id','node_1','node_2',
                           'event_count','count_dyad','ever_together','together_chr',
                           'age_category_1','age_category_2')],
              by = 'dyad_id') %>% 
    ggplot()+
    geom_density(aes(x = draw, group = dyad_id, colour = together_chr))+
    scale_colour_viridis_d(end = 0.5, direction = -1)+
    labs(x = 'edge weight',
         y = 'density',
         colour = 'ever seen together')+
    theme(legend.position = 'bottom'))
ggsave(filename = 'edgedistributions_bayesmix.png',
       # path = '../outputs/sparse_network_methods_figures/',
       path = 'methods_paper/outputs/new_version/',
       plot = edge_dist, device = 'png',
       width = 1600, height = 700, units = 'px')

## save output
save.image('methods_paper/motnp_fit_bayesian_mixture.RData') # load('methods_paper/motnp_fit_bayesian_mixture.RData')

#### eigenvector output checks ####
## get adjacency array
adjarr <- array(NA, dim = c(nrow(nodes), nrow(nodes), nrow(ew)),
                dimnames = list(nodes$id, nodes$id, NULL))
for (i in 1:n_dyads) {
  dyad_row <- counts_df[i, ]
  adjarr[dyad_row$id_1, dyad_row$id_2, ] <- ew[, i]
  adjarr[dyad_row$id_2, dyad_row$id_1, ] <- ew[, i]
}
adjarr[,,1]

## make matrix to save eigenvector values to
eigen <- matrix(NA, nrow = 1000, ncol = nrow(nodes),
                dimnames = list(1:1000, nodes$id))

## calculate eigenvector scores
for(i in 1:nrow(eigen)){
  net <- igraph::graph_from_adjacency_matrix(adjmatrix = adjarr[,,i], diag = FALSE, mode = 'undirected', weighted = T)
  e <- as.data.frame(igraph::eigen_centrality(net)$vector)
  eigen[i,] <- e$`igraph::eigen_centrality(net)$vector`
}
#eigen[,1]
#eigen[1,]

## calculate total dyads associated
nodes$nonzero <- NA
for(i in 1:nrow(nodes)){
  cdf <- counts_df %>%
    filter(node_1 == nodes$node[i] | node_2 == nodes$node[i])
  nodes$nonzero[i] <- length(which(cdf$event_count > 0))
}

## make long format data frame
eigen <- eigen %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'eigen') %>%
  group_by(id) %>%
  mutate(mean = mean(eigen),
         median = median(eigen)) %>% 
  left_join(nodes[,c('id','nonzero','sightings')], by = 'id')

## calculate averages
averages <- eigen %>%
  select(id, mean, median, nonzero,sightings) %>%
  distinct()

## plot eigenvector vs total individuals with which elephant was observed associating
(eigen_vs_nonzero <- ggplot()+
    geom_point(data = eigen, aes(x = nonzero, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5)+
    geom_point(data = averages, aes(x = nonzero, y = median),
               colour = rgb(33/255, 145/255, 140/255),
               size = 0.5)+
    geom_smooth(data = eigen, aes(x = nonzero, y = eigen),
                colour = rgb(68/255,1/255,84/255))+
    scale_x_continuous(name = 'total observed association partners',
                       breaks = seq(0, 150, 25))+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigen1_bayesmix.png',
       # path = '../outputs/sparse_network_methods_figures/',
       path = 'methods_paper/outputs/new_version/',
       plot = last_plot(), device = 'png',
       width = 1000, height = 700, units = 'px')

## plot eigenvector vs total sightings per elephant
(eigen_vs_sightings <- ggplot()+
    geom_point(data = eigen, aes(x = sightings, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5)+
    geom_point(data = averages, aes(x = sightings, y = median),
               colour = rgb(33/255, 145/255, 140/255),
               size = 0.5)+
    geom_smooth(data = eigen, aes(x = sightings, y = eigen),
                colour = rgb(68/255,1/255,84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigensightings_bayesmix.png',
       # path = '../outputs/sparse_network_methods_figures/',
       path = 'methods_paper/outputs/new_version/',
       plot = last_plot(), device = 'png',
       width = 1000, height = 700, units = 'px')

#### combine plots ####
library(patchwork)

## edge weights
(edge_dist + edge_obs)+
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'edges_bayesmix.png',
       # path = '../outputs/sparse_network_methods_figures/',
       path = 'methods_paper/outputs/new_version/',
       plot = last_plot(), device = 'png', width = 2000, height = 1000, units = 'px')

## eigenvector
(eigen_vs_nonzero + eigen_vs_sightings)+
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'eigen_bayesmix.png',
       # path = '../outputs/sparse_network_methods_figures/',
       path = 'methods_paper/outputs/new_version/',
       plot = last_plot(), device = 'png', width = 2000, height = 700, units = 'px')

# save workspace
rm(cdf, dyad_row, e, net) ; gc()
save.image('methods_paper/bayesian_zeroinf_motnpoutputs.RData')
