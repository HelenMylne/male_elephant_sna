#### information ####
# script to just see if there is evidence for sociality affecting lifespan in ANP -- would be handy to join up 2 analyses into a single paper

#### set up ####
library(tidyverse)
theme_set(theme_bw())

#### create data ####
load('anp_nodalregression/anp_short_nodal_modelprep.RData')
nodes <- nodes_all %>% 
  group_by(id, node, node_random) %>% 
  summarise(eigen = mean(mean_eigen),
            sightings = mean(sightings))

lifespan <- readxl::read_xlsx('../data_raw/Raw_ATE_LifeHistoryData_Fishlock220808.xlsx') %>% 
  janitor::clean_names() %>% 
  mutate(id = paste0('M',casename),
         age_death = dyr - byr) %>% 
  dplyr::select(id, age_death, byr, dyr) %>% 
  distinct()

#### plot lifespan versus mortality ####
lifespan_dead <- lifespan %>% 
  filter(age_death > 0)

nodes_dead <- nodes %>% 
  filter(id %in% lifespan_dead$id)

lifespan_dead <- lifespan_dead %>% 
  filter(id %in% nodes_dead$id)

nodes_dead <- nodes_dead %>% 
  left_join(lifespan_dead, by = 'id')

nodes_dead %>% 
  mutate(byr_estimated = ifelse(byr >= 1972, 'measured', 'estimated')) %>% 
  ggplot()+
  geom_point(aes(x = eigen, y = age_death, colour = byr_estimated))

nodes_alive <- nodes %>% 
  anti_join(nodes_dead)

lifespan_alive <- lifespan %>% 
  anti_join(lifespan_dead) 

nodes_alive <- nodes_alive %>% 
  left_join(lifespan_alive, by = 'id')

nodes_alive %>% 
  mutate(current_age = 2019 - byr) %>% 
  ggplot()+
  geom_point(aes(x = eigen, y = current_age))

nodes_alive$alive <- 'yes'
nodes_dead$alive <- 'no'

nodes <- rbind(nodes_alive, nodes_dead) %>% 
  mutate(final_age = ifelse(age_death < 0, 2019-byr, age_death))

nodes %>% 
  mutate(byr_estimated = ifelse(byr >= 1972, 'measured', 'estimated')) %>% 
  filter( ! (alive == 'yes' & byr_estimated == 'estimated')) %>% 
  ggplot()+
  geom_point(aes(x = eigen, y = final_age,
                 colour = alive, shape = byr_estimated, size = sightings),
             alpha = 0.5)+
  geom_smooth(aes(x = eigen, y = final_age,
                  colour = alive, linetype = byr_estimated),
              method = 'lm')+
  scale_colour_viridis_d(end = 0.5)

nodes %>% 
  mutate(byr_estimated = ifelse(byr >= 1972, 'measured', 'estimated')) %>% 
  filter( ! (alive == 'yes' & byr_estimated == 'estimated')) %>% 
  mutate(alive_byr = paste0(alive,'_',byr_estimated)) %>% 
  ggplot()+
  geom_point(aes(x = eigen, y = final_age,
                 colour = sightings,
                 shape = alive_byr, fill = alive_byr))+
  geom_smooth(aes(x = eigen, y = final_age,
                  linetype = alive_byr),
              method = 'loess')+
  scale_colour_viridis_c()+
  scale_shape_manual(values = c(21, 21, 24),
                     breaks = c('yes_measured', 'no_measured', 'no_estimated'))+
  scale_fill_manual(values = c('black', 'transparent', 'transparent'),
                    breaks = c('yes_measured', 'no_measured', 'no_estimated'))

#### plot annual mortality risk vs centrality ####
rm(list = ls()) ; gc()
n <- 1000
a_0 <- rnorm(n, -5.13, 0.720)
a_1 <- rnorm(n, 3.00, 0.105)
c   <- rnorm(n, 0.0258, 0.00669)
b_0 <- rnorm(n, -5.07, 0.567)
b_1 <- rnorm(n, 0.0948, 0.0178)

go_bt <- function(age){
  mortality <- exp(b_0 + b_1 * age) + exp(a_0 - a_1 * age) + c
  return(mortality)
}

age <- 1:40
mortality <- matrix(NA, nrow = n, ncol = length(age))
for( j in 1:ncol(mortality) ){
  mortality[,j] <- go_bt(age[j])
}

mort <- as.data.frame(mortality)
colnames(mort) <- age
mort <- mort %>% 
  pivot_longer(cols = everything(), names_to = 'age', values_to = 'mortality') %>% 
  mutate(age = as.numeric(age))

ggplot(mort)+
  geom_point(aes(x = age, y = mortality),
             colour = 'blue', alpha = 0.1)

mort2 <- data.frame(age = age,
                    mean = apply(mortality, 2, mean),
                    median = apply(mortality, 2, median),
                    upr = apply(mortality, 2, quantile, prob = 0.975),
                    lwr = apply(mortality, 2, quantile, prob = 0.025))

ggplot(mort2)+
  geom_ribbon(aes(x = age, ymin = lwr, ymax = upr),
              alpha = 0.2)+
  geom_line(aes(x = age, y = mean),
             colour = 'blue', lwd = 1.2)

## read in elephant data
load('anp_nodalregression/anp_short_nodal_modelprep.RData')
rm(cents_all, covs_all, nodes) ; gc()

mortality <- matrix(NA, nrow = n, ncol = nrow(nodes_all))
for( j in 1:ncol(mortality) ){
  mortality[,j] <- go_bt(nodes_all$age[j])
}

nodes_all$mort_mu  <- apply(mortality, 2, mean)
nodes_all$mort_mid <- apply(mortality, 2, median)
nodes_all$mort_upr <- apply(mortality, 2, quantile, prob = 0.975)
nodes_all$mort_lwr <- apply(mortality, 2, quantile, prob = 0.025)

ggplot(nodes_all)+
  geom_point(aes(x = age, y = mort_mu),
             colour = 'blue', alpha = 0.1)


ggplot(nodes_all)+
  geom_point(aes(x = mean_eigen,
                 y = mort_mu),
             colour = '#1F968BFF', 
             alpha = 0.2)+
  labs(x = 'risk of mortality (mean of posterior)',
       y = 'logit eigenvector centrality (mean of posterior)')

ggplot(nodes_all)+
  geom_point(aes(x = mean_eigen,
                 y = mort_mu),
             colour = '#1F968BFF', 
             alpha = 0.2)+
  labs(x = 'risk of mortality (mean of posterior)',
       y = 'logit eigenvector centrality (mean of posterior)')
ggsave(plot = last_plot(),
       filename = 'mortality_vs_eigen.png',
       path = '../outputs/',
       device = 'png', height = 1200, width = 1500, unit = 'px')























