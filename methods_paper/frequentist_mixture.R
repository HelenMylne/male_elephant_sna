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
mixture_fit_joint_random <- glmmTMB(cbind(event_count, apart) ~ (1|dyad_males),
                                    ziformula = ~1,
                                    data = counts_df,
                                    family = binomial(link = 'logit'))

## summary of the model
(summary <- summary(mixture_fit_joint_random))
invlogit(summary$coefficients$cond[1]) # estimate
invlogit(summary$coefficients$cond[1] + summary$coefficients$cond[2]) # upper
invlogit(summary$coefficients$cond[1] - summary$coefficients$cond[2]) # lower

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
width <- 750      # min = 789, max = 2250
height <- 700     # max = 2625
resolution <- 300 # min = 300, max = 600, but changes depending on dimensions

colours <- c('#21918c','#440154')
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
    scale_x_continuous(name = 'mixture model weight',
                       breaks = c(0.00,0.10,0.20))+
    scale_y_continuous(name = 'number of dyads',
                       expand = c(0,0))+
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

edges$together <- ifelse(counts_df$event_count == 0,
                             'never together', 'together at least once')

(edgesightings_mixture2 <- ggplot()+
    geom_point(data = edges,
               aes(x = count_dyad, y = prdctn_invlogit),
               colour = rgb(94/255, 201/255, 98/255, 0.2),#rgb(33/255, 145/255, 140/255, 0.1),
               size = 1, shape = 19)+
    geom_smooth(data = edges[edges$together == 'together at least once',],
                aes(x = count_dyad, y = prdctn_invlogit,
                    colour = together)#linetype = hack_linetype_together),
                #colour = rgb(68/255, 1/255, 84/255)
                )+
    # geom_smooth(data = edges,
    #             aes(x = count_dyad, y = prdctn_invlogit,
    #                 linetype = hack_linetype_all),
    #             linewidth = 0.8,
    #             colour = rgb(68/255, 1/255, 84/255))+
    geom_smooth(data = edges[edges$together == 'never together',],
                aes(x = count_dyad, y = prdctn_invlogit,
                    colour = together)
                # data = data.frame(x = 2:max(counts_df$count_dyad[counts_df$event_count == 0]),
                #                   y = 0, together = 'never together'),
                # aes(x = x, y = y, colour = together)#,
    )+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mixture weight', limits = c(-0.02,1.02), expand = c(0,0))+
    scale_colour_manual(values = colours,
                        breaks = c('never together','together at least once'),
                        name = NULL)+
    # scale_linetype_manual(values = c(1,6),
    #                       #breaks = c('never together','together at least once'),
    #                       breaks = c('all values','together at least once'),
    #                       name = 'sightings together')+
    theme(legend.position = 'right',#c(0.56,0.76),
          #legend.background = element_rect(fill = 'white', colour = 'black'),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8))
)
ggsave(filename = 'edgesightings_mixture_twolines_changealpha.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = edgesightings_mixture2, device = 'png', width = 700, height = 700, units = 'px')

## merge
edges_mixture <- edges_mixture +
  theme(text = element_text(family = 'serif'),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        # legend.title = element_blank(),
        # legend.position = 'none',
        legend.text = element_text(size = 10))+
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))
edgesightings_mixture2 <- edgesightings_mixture2 +
  theme(text = element_text(family = 'serif'),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))+
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

(edges_mixture + edgesightings_mixture2)+
  plot_annotation(tag_levels = 'a') &
    theme(legend.position = 'bottom',
          text = element_text(family = 'serif'),
          legend.text = element_text(size = 10),
          legend.title = element_blank()) &
    guides(colour = guide_legend(nrow = 2, byrow = TRUE))
ggsave(filename = 'outputs_mixturemodel.tiff',
       path = 'methods_paper/outputs/',
       plot = last_plot(), device = 'tiff',
       width = (width+50)*2,
       height = height+100,
       dpi = resolution,
       units = 'px')

#### calculate eigenvector ####
# identify elephants
elephants <- unique(c(counts_df$id_1, counts_df$id_2))

# calculate adjacency matrix
mix_adjmat <- matrix(NA, nrow = length(elephants), ncol = length(elephants),
                     dimnames = list(x = elephants, y = elephants))
for(i in 1:nrow(mix_adjmat)){
  for(j in 1:ncol(mix_adjmat)){
    if(i <= j){
      mix_adjmat[i,j] <- ifelse(i == j, 0,
                                edges$prediction[which(counts_df$id_1 == rownames(mix_adjmat)[i] &
                                                      counts_df$id_2 == colnames(mix_adjmat)[j])])
    } else {
      mix_adjmat[i,j] <- mix_adjmat[j,i]
    }
  }
}
mix_adjmat <- invlogit(mix_adjmat)

# create network
mix_net <- igraph::graph_from_adjacency_matrix(mix_adjmat, weighted = T, mode = 'undirected', diag = FALSE)

# calculate eigenvector centrality
eigen_mix <- igraph::eigen_centrality(mix_net)$vector %>%
  as.data.frame()
colnames(eigen_mix) <- 'eigen'
eigen_mix$id <- rownames(eigen_mix)

# create data frame showing number of sightings per individual with other elephants
together0 <- counts_df %>% filter(event_count == 0)
eigen_mix$together0 <- NA ; for(i in 1:nrow(eigen_mix)) {
  x <- together0 %>%
    filter(id_1 == eigen_mix$id[i] | id_2 == eigen_mix$id[i])
  eigen_mix$together0[i] <- nrow(x)
}

# calculate number of sightings per elephant
counts <- data.frame(id = elephants, count = NA)
for(i in 1:nrow(counts)){
  counts$count[i] <- ifelse(i == nrow(counts),
                            unique(counts_df$count_2[which(counts_df$id_2 == counts$id[i])]),
                            unique(counts_df$count_1[which(counts_df$id_1 == counts$id[i])]))
}
eigen_mix <- eigen_mix %>%
  left_join(counts, by = 'id')

# plot
(eigensightings_mix.1 <- ggplot(eigen_mix, aes(x = count, y = eigen))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigensightings_mix_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_mix.1, device = 'png', width = 1600, height = 700, units = 'px')

(eigensightings_mix.2 <- ggplot(eigen_mix, aes(x = count, y = eigen))+
    geom_point(colour = rgb(33/255, 145/255, 140/255),
               size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigensightings_mix_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigensightings_mix.2, device = 'png', width = 1600, height = 700, units = 'px')

(eigen0_mix.1 <- ggplot(eigen_mix, aes(x = together0, y = eigen))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigentogether0_mix_noline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_mix.1, device = 'png', width = 1600, height = 700, units = 'px')

(eigen0_mix.2 <- ggplot(eigen_mix, aes(x = together0, y = eigen))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)
ggsave(filename = 'eigentogether0_mix_withline.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = eigen0_mix.2, device = 'png', width = 1600, height = 700, units = 'px')

## change line breaks
(eigen0_mix.2 <- ggplot(eigen_mix, aes(x = together0, y = eigen))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0',
                       breaks = c(100,150,200))+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

# save workspace
save.image('../outputs/sparse_network_methods_figures/mixture_model.RData')

# merge
eigen0_mix.2 <- eigen0_mix.2 +
  theme(text = element_text(family = 'serif'),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
eigensightings_mix.2 <- eigensightings_mix.2 +
  theme(text = element_text(family = 'serif'),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

(eigen0_mix.2 + eigensightings_mix.2) +
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'eigen_outputs_mixture.tiff',
       path = 'methods_paper/outputs/',
       plot = last_plot(), device = 'tiff',
       width = (width+50)*2,
       height = height,
       dpi = resolution,
       units = 'px')
# ggsave(filename = 'eigen_outputs_mixture.png',
#        path = 'methods_paper/outputs/',
#        plot = last_plot(), device = 'png', width = 1600, height = 700, units = 'px')

# #### I think these are wrong ####
# ## fit a zero-inflated Poisson model using glmmTMB
# mixture_fit_joint <- glmmTMB(cbind(event_count, apart) ~ 1,
#                              ziformula = ~1, data = counts_df, family = binomial(link = 'logit'))
# mixture_fit_together <- glmmTMB(event_count ~ count_dyad,
#                                 ziformula = ~1, data = counts_df, family = poisson)
# mixture_fit_together_random <- glmmTMB(event_count ~ count_dyad,
#                                        ziformula = ~1, data = counts_df, family = poisson)
# 
# ## summary of the model
# summary(mixture_fit_joint)
# summary(mixture_fit_together)
# 
# ## view data used for model
# head(mixture_fit_joint$frame)
# View(mixture_fit_together$frame) # takes data from both variables
# 
# ## extract fit
# mixture_fit_joint$fit
# mixture_fit_together$fit
# 
# ## extract standard error
# mixture_fit_joint$sdr
# mixture_fit_together$sdr
# 
# ## extract edges
# predictions <- predict(object = mixture_fit_joint, newdata = counts_df)
# hist(invlogit(predictions), breaks = 50) # all the same
# predictions <- predict(object = mixture_fit_together, newdata = counts_df)
# hist(predictions)
