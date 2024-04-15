#### Information ####
# script to run basic frequentist analysis of edge weights using mixture model
# to respond to Colin's comments on methods paper

#### setup ####
## load packages
library(tidyverse)
library(glmmTMB)
library(LaplacesDemon)
library(patchwork)

## set seed
set.seed(1)

## set plotting aesthetics
theme_set(theme_bw(base_size = 12))

#### load data ####
load('motnp_edgeweights_conditionalprior.RData')
rm(list = ls()[! ls() %in% c('counts_df')])

#### run model ####
## fit a zero-inflated binomial model using glmmTMB
mixture_fit_joint_random <- glmmTMB(cbind(event_count, apart) ~ as.factor(dyad_males),
                                    ziformula = ~1,
                                    data = counts_df,
                                    family = binomial(link = 'logit'))

## summary of the model
summary(mixture_fit_joint_random)

## view data used for model
head(mixture_fit_joint_random$frame)

## extract fit
mixture_fit_joint_random$fit

## save edges
edges <- mixture_fit_joint_random$fit$parfull %>% 
  as.data.frame()
colnames(edges) <- 'edge_weight'
rownames(edges) <- c('beta','betazi',counts_df$dyad_males,'theta')
edges <- edges %>% 
  mutate(parameter = c('beta','betazi',counts_df$dyad_males,'theta'))
params <- edges %>% 
  filter(parameter %in% c('beta','betazi','theta'))
edges <- edges %>% 
  anti_join(params) %>% 
  rename(dyad_males = parameter) %>% 
  mutate(dyad_males = as.integer(dyad_males)) %>% 
  left_join(counts_df, by = 'dyad_males')

## extract standard error
mixture_fit_joint_random$sdr

## extract predictions
predictions <- predict(object = mixture_fit_joint_random, newdata = counts_df)
hist(invlogit(predictions), breaks = 50) # all low

## save predictions
predictions <- predictions %>% 
  as.data.frame() %>% 
  mutate(dyad_males = counts_df$dyad_males)
colnames(predictions)[1] <- 'prediction'
edges <- edges %>% 
  left_join(predictions, by = 'dyad_males') %>% 
  relocate(prediction, .after = 'edge_weight')

#### plot outputs ####
edges <- edges %>% 
  mutate(weight_invlogit = invlogit(edge_weight),
         prdctn_invlogit = invlogit(prediction))

## basic plots
plot(edges$prediction ~ edges$edge_weight)      # exactly match, just at different scale values
plot(edges$weight_invlogit ~ edges$count_dyad)  # looks decent, just too high
plot(edges$prdctn_invlogit ~ edges$count_dyad)  # looks similar to other plots

## nice plot edges
(edges_mixture <- ggplot(edges)+
    geom_bar(aes(x = round(prdctn_invlogit, 3)),
             fill = rgb(33/255, 145/255, 140/255),
             colour = rgb(33/255, 145/255, 140/255), linewidth = 0.2)+
    scale_x_continuous(name = 'mixture model weight')+
    scale_y_continuous(name = 'number of dyads',
                       expand = c(0,0))+
    # annotate('text', x = 0.5, y = 165,
    #          label = paste0(length(which(edges$prdctn_invlogit == 0)),
    #                         ' dyads with \nSRI = 0'),
    #          size = unit(4, 'pt'),
    #          colour = rgb(33/255, 145/255, 140/255))+
    coord_cartesian(ylim = c(0,5300))
)

## nice plot edgesightings
(edgesightings_mixture1 <- ggplot(edges,
       aes(x = count_dyad, y = prdctn_invlogit))+
  geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1),
             size = 0.5, shape = 19)+
  geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
  scale_x_continuous(name = 'total dyad sightings')+
  scale_y_continuous(name = 'mixture model weight', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'edgesightings_mixturemodel_withline_changealpha.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = edgesightings_mixture1, device = 'png', width = 700, height = 700, units = 'px')

edges$hack_linetype_all <- 'all values'
edges$hack_linetype_together <- 'together at least once'
(edgesightings_mixture2 <- ggplot()+
    geom_point(data = edges,
               aes(x = count_dyad, y = prdctn_invlogit),
               colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    geom_smooth(data = edges[edges$event_count > 0,],
                aes(x = count_dyad, y = prdctn_invlogit,
                    linetype = hack_linetype_together),
                linewidth = 0.8,
                colour = rgb(68/255, 1/255, 84/255))+
    geom_smooth(data = edges,
                aes(x = count_dyad, y = prdctn_invlogit,
                    linetype = hack_linetype_all),
                linewidth = 0.8,
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mixture weight', limits = c(-0.02,1.02), expand = c(0,0))+
    scale_linetype_manual(values = c(1,6),
                          #breaks = c('never together','together at least once'),
                          breaks = c('all values','together at least once'),
                          name = 'sightings together')+
    theme(legend.position = c(0.56,0.76),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8))
)
ggsave(filename = 'edgesightings_mixture_twolines_changealpha.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = edgesightings_mixture2, device = 'png', width = 700, height = 700, units = 'px')

## merge
(edges_mixture + edgesightings_mixture2)+
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'outputs_mixturemodel.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 1600, height = 700, units = 'px')

#### I think these are wrong ####
## fit a zero-inflated Poisson model using glmmTMB
mixture_fit_joint <- glmmTMB(cbind(event_count, apart) ~ 1,
                             ziformula = ~1, data = counts_df, family = binomial(link = 'logit'))
mixture_fit_together <- glmmTMB(event_count ~ count_dyad,
                                ziformula = ~1, data = counts_df, family = poisson)
mixture_fit_together_random <- glmmTMB(event_count ~ count_dyad,
                                       ziformula = ~1, data = counts_df, family = poisson)

## summary of the model
summary(mixture_fit_joint)
summary(mixture_fit_together)

## view data used for model
head(mixture_fit_joint$frame)
View(mixture_fit_together$frame) # takes data from both variables

## extract fit
mixture_fit_joint$fit
mixture_fit_together$fit

## extract standard error
mixture_fit_joint$sdr
mixture_fit_together$sdr

## extract edges
predictions <- predict(object = mixture_fit_joint, newdata = counts_df)
hist(invlogit(predictions), breaks = 50) # all the same
predictions <- predict(object = mixture_fit_together, newdata = counts_df)
hist(predictions)
