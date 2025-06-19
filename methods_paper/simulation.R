######## information ########
# script to simulate dyadic network data for estimating social relationships using two-step BISoN

######## set up ########
# library(tidyverse) ; library(LaplacesDemon) ; library(patchwork) ; library(igraph)
library(tidyverse, lib.loc = '../packages/')
library(LaplacesDemon, lib.loc = '../packages/')
library(patchwork, lib.loc = '../packages/')
library(igraph, lib.loc = '../packages/')

## define plot theme
theme_set(theme_bw(base_size = 12))
colours <- c('#fde725','#5ec962','#21918c','#3b528b','#440154')

######## set parameters for simulation ########
## set seed
set.seed(5)

## set population parameters
n_nodes <- 100
n_dyads <- cumsum(1:n_nodes)[n_nodes-1]

## set probability parameters
p_assoc <- 0.1  # individual edge weight (make some more popular than others)
sd_assoc<- 0.5  # variation around p_assoc
p_apart <- 0.4  # proportion of dyads that truly never associate
p_sight <- 4    # average rate of observation per elephant

######## create nodes data ########
## create data frame
sim_nodes <- data.frame(node = 1:n_nodes,
                        sightings = rnbinom(n = n_nodes, prob = 0.3, size = p_sight),
                        social_score = invlogit(rnorm(n_nodes, p_assoc, 0.5))) %>% 
  mutate(age = rweibull(n = n_nodes, shape = 1.2, scale = 9) + 10) %>% # ages starting at 10 years (approximate age of dispersal)
  mutate(node1 = node, node2 = node)

## plot ages and sightings
ggplot(sim_nodes)+
  geom_density(aes(x = age))
ggplot(sim_nodes)+
  geom_bar(aes(x = as.factor(sightings)))+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 25), minor_breaks = 0:25)+
  labs(x = 'sightings', y = 'count of nodes')

######## create dyadic data ########
## create data frame
sim_dyads <- expand.grid(node1 = sim_nodes$node,
                         node2 = sim_nodes$node) %>% 
  filter(node1 < node2) %>% 
  left_join(sim_nodes[,c('node1', 'sightings','social_score')], by = 'node1') %>% 
  rename(sightings1 = sightings,
         soc1 = social_score) %>% 
  left_join(sim_nodes[,c('node2', 'sightings','social_score')], by = 'node2') %>% 
  rename(sightings2 = sightings,
         soc2 = social_score) %>% 
  mutate(dyad_id = as.integer(as.factor(paste0(node1,'_',node2)))) %>% 
  relocate(dyad_id, .before = 'node1') %>% 
  relocate(sightings2, .after = 'sightings1')

######## simulate true edge weights ########
## set true dyadic edge weights
sim_dyads <- sim_dyads %>% 
  mutate(true_edge = rnorm(n = n_dyads, 
                           # mean = logit(soc1 * soc2), sd = 0.1)) %>% 
                           mean = logit(p_assoc), sd = 1.2)) %>% 
  mutate(true_edge = invlogit(true_edge))
ggplot(sim_dyads)+
  geom_density(aes(x = true_edge))

## set some which truly do not associate
dyads_apart <- sample(x = 1:n_dyads, size = n_dyads*p_apart, replace = F)
sim_dyads$true_edge[dyads_apart] <- 0
ggplot(sim_dyads)+
  geom_density(aes(x = true_edge))

######## simulate raw observation data ########
sim_dyads$together <- NA
for(i in 1:nrow(sim_dyads)){
  sim_dyads$together[i] <- rbinom(size = min(sim_dyads$sightings1[i],sim_dyads$sightings2[i]),
                                  n = 1,
                                  prob = sim_dyads$true_edge[i])
}
sim_dyads <- sim_dyads %>% 
  mutate(total_sightings = (sightings1 + sightings2) - together,
         apart = (sightings1 + sightings2) - together*2)

## plot edges vs observations (should be a slight negative correlation: high edges are unlikely anyway due to simulation, then some may be removed and made to be zero, strong association = lower total number of sightings together because more often observed at the same time)
ggplot(sim_dyads)+
  geom_point(aes(x = true_edge, y = total_sightings),
             colour = colours[5], alpha = 0.05)+
  geom_smooth(aes(x = true_edge, y = total_sightings),
              method = 'lm',
              colour = colours[3])+
  scale_y_continuous(breaks = seq(0,max(sim_dyads$total_sightings), by = 5))+
  labs(y = 'total observations',
       x = 'true edge weight')

## plot edges vs observations together (should be a positive correlation)
ggplot(sim_dyads)+
  geom_point(aes(x = true_edge, y = together, size = total_sightings),
             colour = colours[5], alpha = 0.05)+
  geom_smooth(aes(x = true_edge, y = together),
              method = 'lm',
              colour = colours[3])+
  scale_y_continuous(breaks = 0:max(sim_dyads$total_sightings))+
  labs(size = 'total observations',
       x = 'true edge weight',
       y = 'sightings together')

## save simulation
save.image('methods_paper/simulation_outputs/create_simulated_data.RData')
write_csv(sim_dyads, 'methods_paper/simulation_outputs/sim_data.csv')

######## SRI ########################
#### SRI edge bar graph ####
# calculate sri
sim_dyads$sri <- sim_dyads$together / sim_dyads$total_sightings

# calculate percentages of 0s and 1s using SRI
( length(which(sim_dyads$sri==0))/length(sim_dyads$sri) ) * 100
( length(which(sim_dyads$sri==1))/length(sim_dyads$sri) ) * 100

# calculate percentages of 0s and 1s using SRI, assuming a 3 sighting threshold per node
subset3 <- sim_dyads %>% filter(sightings1 >= 3 & sightings2 >= 3)
( length(which(subset3$sri==0))/length(subset3$sri) ) * 100
( length(which(subset3$sri==1))/length(subset3$sri) ) * 100

# calculate percentages of 0s and 1s using SRI, assuming a 5 sighting threshold per node
subset5 <- sim_dyads %>% filter(sightings1 >= 5 & sightings2 >= 5)
( length(which(subset5$sri==0))/length(subset5$sri) ) * 100

# redo annotation on so that size isn't too large for multi plot
(edges_sri <- ggplot()+
    geom_bar(data = sim_dyads, aes(x = round(sri, 3)), fill = rgb(33/255, 145/255, 140/255),
             colour = rgb(33/255, 145/255, 140/255), linewidth = 0.2)+
    scale_x_continuous(name = 'SRI value')+
    scale_y_continuous(name = 'number of dyads',
                       expand = c(0,0))+
    annotate('text', x = 0.5, y = 165,
             label = paste0(length(which(sim_dyads$sri == 0)), ' dyads with \nSRI = 0'),
             size = unit(4, 'pt'),
             colour = rgb(33/255, 145/255, 140/255))+
    coord_cartesian(ylim = c(0,200))
)

#### SRI edges vs sighting count ####
sim_dyads$never <- ifelse(sim_dyads$together == 0,
                             'never together', 'together at least once')
(edgesightingssri <- ggplot()+
    geom_point(data = sim_dyads,
               aes(x = total_sightings, y = sri),
               colour = colours[2],
               size = 1, shape = 19)+
    geom_smooth(data = sim_dyads[sim_dyads$never == 'together at least once',],
                aes(x = total_sightings, y = sri, colour = never),
                # colour = colours[3]
    )+
    geom_smooth(data = data.frame(x = 2:max(sim_dyads$total_sightings[sim_dyads$together == 0]),
                                  y = 0, never = 'never together'),
                aes(x = x, y = y, colour = never)
                # colour = colours[5]
    )+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'SRI weight', limits = c(-0.02,1.02), expand = c(0,0))+
    scale_colour_manual(values = colours[c(3,5)],
                        breaks = c('never together','together at least once'),
                        name = NULL
    )+
    theme(legend.position = 'right',
          legend.key.height = unit(4, 'mm'),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8)))

## clean up
rm(subset3, subset5, dyads_apart, i) ; gc()

#### merge ####
(combined <- (edges_sri + edgesightingssri)+
   plot_annotation(tag_levels = 'a'))
ggsave(filename = 'sim_sri.png',
       path = 'methods_paper/simulation_outputs/',
       plot = last_plot(), device = 'png', width = 2100, height = 700, units = 'px')

