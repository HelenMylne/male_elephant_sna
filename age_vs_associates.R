######## Information ####
# Analysis script for paper on statistical methods for social network analysis using gambit-of-the-group data on male elephant population of Mosi-Oa-Tunya National Park, Zambia.
# Data collected: 19th May 2016-21st December 2017.
# Collected by: Mr David Youldon, Mr Dabwiso Sakala, Miss Helen Mylne and African Lion and Environmental Research Trust during this time
# Data supplied by: African Lion and Environmental Research Trust (ALERT)

######## Set up ####
library(tidyverse, lib.loc = '../packages')      # library(tidyverse)

# define plot theme
theme_set(theme_bw(base_size = 12))

######## plot age versus total SRI = 0 #######
#### motnp ####
load('methods_paper/r_workspaces/model_run_conditional.RData')

nodes$together0 <- NA
for(i in 1:nrow(nodes)){
  n1 <- length(which(counts_df$event_count == 0 & counts_df$id_1 == nodes$id[i]))
  n2 <- length(which(counts_df$event_count == 0 & counts_df$id_2 == nodes$id[i]))
  nodes$together0[i] <- n1 + n2
}

nodes$together1 <- (nrow(motnp) - 1) - nodes$together0

motnp <- nodes
ggplot(data = motnp)+
  geom_boxplot(aes(x = age, y = together1),
               notch = T)+
  geom_point(aes(x = age, y = together1, size = sightings),
             colour = rgb(0,0,1,0.2))+
  labs(y = 'number of partners')

rm(list = ls()[! ls() %in% c('motnp')]) ; gc()

#### short ####
load('anp_edgecalculations/excluding_under10s/anpshort1_edgeweights_conditionalprior.RData')

short <- nodes %>% 
  mutate(period = 1)
rm(list = ls()[! ls() %in% c('motnp','counts_df','short')]) ; gc()

for(window in 2:max(counts_df$period)){
  cdf <- counts_df %>% 
    filter(period == window)
  nodes <- cdf %>% 
    dplyr::select(id_1, node_1, age_start_1, period_count_1, period) %>% 
    rename(id = id_1,
           node = node_1,
           age = age_start_1,
           sightings = period_count_1) %>% 
    distinct()
  short <- rbind(short, nodes)
}
rm(list = ls()[! ls() %in% c('motnp','counts_df','short')]) ; gc()

short$together0 <- NA
for(i in 1:nrow(short)){
  n1 <- length(which(counts_df$event_count == 0 & counts_df$id_1 == short$id[i] & counts_df$period == short$period[i]))
  n2 <- length(which(counts_df$event_count == 0 & counts_df$id_2 == short$id[i] & counts_df$period == short$period[i]))
  short$together0[i] <- n1 + n2
}

short$together1 <- NA
for(i in 1:nrow(short)){
  x <- length(which(short$period == short$period[i]))
  short$together1[i] <- (x-1) - short$together0[i]
}

(g <- ggplot(data = short, mapping = aes(x = age, y = together1))+
  geom_point(colour = '#21918c', alpha = 0.2)+
  geom_smooth(method = 'loess', colour = '#440154')+
  labs(y = 'number of partners') )
g + facet_wrap(period ~ .)

#### long ####
load('anp_edgecalculations/excluding_under10s/anplong1_edgeweights_conditionalprior.RData')

long <- nodes %>% 
  mutate(period = 1)
rm(list = ls()[! ls() %in% c('motnp','counts_df','short','long')]) ; gc()

for(window in 2:max(counts_df$period)){
  cdf <- counts_df %>% 
    filter(period == window)
  nodes <- cdf %>% 
    dplyr::select(id_1, node_1, age_start_1, count_period_1, period) %>% 
    rename(id = id_1,
           node = node_1,
           age = age_start_1,
           sightings = count_period_1) %>% 
    distinct()
  long <- rbind(long, nodes)
}

long$together0 <- NA
for(i in 1:nrow(long)){
  n1 <- length(which(counts_df$event_count == 0 & counts_df$id_1 == long$id[i] & counts_df$period == long$period[i]))
  n2 <- length(which(counts_df$event_count == 0 & counts_df$id_2 == long$id[i] & counts_df$period == long$period[i]))
  long$together0[i] <- n1 + n2
}

long$together1 <- NA
for(i in 1:nrow(long)){
  x <- length(which(long$period == long$period[i]))
  long$together1[i] <- (x-1) - long$together0[i]
}

rm(list = ls()[! ls() %in% c('motnp','short','long')]) ; gc()

(g <- ggplot(data = long, mapping = aes(x = age, y = together1))+
    geom_point(colour = '#21918c', alpha = 0.2)+
    geom_smooth(method = 'loess', colour = '#440154')+
    labs(y = 'number of partners') )
g + facet_wrap(period ~ .)

short$length <- 'short'
long$length <- 'long'
anp <- rbind(short, long)

ggplot(data = anp,
       mapping = aes(x = age,
                     y = together0,
                     colour = length))+
  geom_point(alpha = 0.2)+
  geom_smooth(method = 'loess')+
  scale_colour_viridis_d()

ggplot(data = anp,
       mapping = aes(x = age, y = together1, colour = length))+
  geom_point(alpha = 0.2)+
  geom_smooth(method = 'loess')+
  scale_colour_viridis_d()+
  labs(y = 'number of partners')
