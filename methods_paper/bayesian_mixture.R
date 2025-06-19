######## information ########
# run Bayesian mixture model using first simulated then empirical data

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
n_draws <- 1000

## define plot theme
theme_set(theme_bw(base_size = 12))

######## simulation ########
#### load data ####
load('methods_paper/simulation_outputs/create_simulated_data.RData')
sim_dyads <- sim_dyads %>% 
  mutate(ever = ifelse(together == 0, 0, 1))

#### load model ####
mod <- cmdstan_model('methods_paper/models/edge_mixture_bayesian3.stan')

#### create model data list ####
sim_data <- list(n_dyads = n_dyads,
                 ever_together = sim_dyads$ever,
                 count_together = sim_dyads$together,
                 total_sightings = sim_dyads$total_sightings,
                 dyad_id = sim_dyads$dyad_id)

#### run model ####
fit <- mod$sample(data = sim_data, seed = 12345,
                  chains = n_chains, parallel_chains = n_chains,
                  iter_warmup = n_draws, iter_sampling = n_draws)

save.image('methods_paper/simulation_outputs/fit_bayesian_mixture3.RData')

## summary
fit$summary()

#### extract edges ####
## extract
draws <- as.data.frame(fit$draws())

## put all chains into single column: separate out chains
draws1 <- draws[,seq(1, ncol(draws), by = 4)]
draws2 <- draws[,seq(2, ncol(draws), by = 4)]
draws3 <- draws[,seq(3, ncol(draws), by = 4)]
draws4 <- draws[,seq(4, ncol(draws), by = 4)]

## bind 4 chains into single data frame
colnames(draws2) <- colnames(draws3) <- colnames(draws4) <- colnames(draws1)
draws <- rbind(draws1, draws2, draws3, draws4)
rm(draws1, draws2, draws3, draws4) ; gc()

## rename without "1." on start
params <- data.frame(param = colnames(draws)) %>% 
  mutate(param = str_remove(pattern = '1.', string = param))
colnames(draws) <- params$param

## add columns for draw/chain
draws <- draws %>% 
  mutate(draw = rep(1:n_draws, n_chains),
         chain = rep(1:n_chains, each = n_draws))

## check know where all parameters are coming from
global <- params[params$param %in% c('lp__','logit_together','average_edge','prob_together'),]
dyad_effects <- params[grep(x = params$param, pattern = 'dyad_effect'),]
hypothetical <- params[grep(x = params$param, pattern = 'hypothetical_together'),]
logit_edges <- params[grep(x = params$param, pattern = 'logit_edge'),]
edge_weights <- params[grep(x = params$param, pattern = 'edge_weight'),]
nrow(params) == length(global) + length(dyad_effects) + length(hypothetical) + length(logit_edges) + length(edge_weights)
global <- global[2:4]

#### check fit ####
## traceplot global parameters
draws %>% 
  select(all_of(global), draw, chain) %>% 
  pivot_longer(cols = all_of(global),
               names_to = 'param', values_to = 'estimate') %>% 
  ggplot()+
  geom_line(aes(x = draw, y = estimate, colour = as.factor(chain)))+
  facet_wrap(param ~ ., scales = 'free')

## sample dyads for plotting dyad-level parameters
plot_dyads <- sample(1:n_dyads, 12, F)

## dyad_effects
d <- draws %>% 
  select(all_of(dyad_effects), draw, chain) %>% 
  select(all_of(plot_dyads), draw, chain)
d %>% 
  pivot_longer(cols = 1:(ncol(d)-2),
               names_to = 'param', values_to = 'estimate') %>% 
  ggplot()+
  geom_line(aes(x = draw, y = estimate, colour = as.factor(chain)))+
  facet_wrap(param ~ ., scales = 'free')

## hypothetical
d <- draws %>% 
  select(all_of(hypothetical), draw, chain) %>% 
  select(all_of(plot_dyads), draw, chain)
d %>% 
  pivot_longer(cols = 1:(ncol(d)-2),
               names_to = 'param', values_to = 'estimate') %>% 
  ggplot()+
  geom_line(aes(x = draw, y = estimate, colour = as.factor(chain)))+
  facet_wrap(param ~ ., scales = 'free')

## logit_edges
d <- draws %>% 
  select(all_of(logit_edges), draw, chain) %>% 
  select(all_of(plot_dyads), draw, chain)
d %>% 
  pivot_longer(cols = 1:(ncol(d)-2),
               names_to = 'param', values_to = 'estimate') %>% 
  ggplot()+
  geom_line(aes(x = draw, y = estimate, colour = as.factor(chain)))+
  facet_wrap(param ~ ., scales = 'free')

## edge_weights
d <- draws %>% 
  select(all_of(edge_weights), draw, chain) %>% 
  select(all_of(plot_dyads), draw, chain)
d %>% 
  pivot_longer(cols = 1:(ncol(d)-2),
               names_to = 'param', values_to = 'estimate') %>% 
  ggplot()+
  geom_line(aes(x = draw, y = estimate, colour = as.factor(chain)))+
  facet_wrap(param ~ ., scales = 'free')
rm(global,dyad_effects,hypothetical,logit_edges,d) ; gc()

#### plot edges -- NOT ACTUALLY SURE WHICH ONE EVEN IS THE EDGE WEIGHT IN THE END!! ####
## extract estimates
ew <- draws %>% 
  select(all_of(edge_weights))

# ## add probability together??
# for(i in 1:n_dyads){
#   ew[,i] <- ew[,i] + draws$prob_together
# }

## calculate average values
sim_dyads <- sim_dyads %>% 
  mutate(mean_edge = apply(ew, 2, mean),
         mid_edge = apply(ew, 2, median),
         stdv_edge = apply(ew, 2, sd))

## plot average estimate against true value
ggplot(sim_dyads)+
  geom_point(aes(x = true_edge, y = mean_edge, colour = as.factor(ever)))+
  geom_abline(slope = 1, intercept = 0)+
  labs(x = 'true (simulated) dyadic edge weight',
       y = 'mean edge estimate',
       colour = 'ever seen\ntogether')
ggplot(sim_dyads)+
  geom_point(aes(x = true_edge, y = mid_edge, colour = as.factor(ever)))+
  geom_abline(slope = 1, intercept = 0)+
  labs(x = 'true (simulated) dyadic edge weight',
       y = 'mean edge estimate',
       colour = 'ever seen\ntogether')

## plot uncertainty (sd) against sightings
ggplot(sim_dyads)+
  geom_point(aes(x = total_sightings, y = stdv_edge, colour = as.factor(ever)))+
  labs(x = 'total sightings of dyad',
       y = 'uncertainty (SD) in edge estimate',
       colour = 'ever seen\ntogether')

