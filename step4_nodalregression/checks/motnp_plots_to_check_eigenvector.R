#### information and set up ####
# script to check the values being produced by bisonR for eigenvector centrality because they don't match those produced by igraph and sna for the SRI data, and the centrality appears extremely closely linked to the number of sightings per individual

library(tidyverse) ; library(bisonR) ; library(igraph)
theme_set(theme_classic())

#### load data ####
#load('motnp_bisonr_edgescalculated_widebinomconj_v2.RData')
load('motnp_nodalregression_meanage.RData')

### read in data (produced in posterior check v2 of motnp_nodalregression_meanage.R)
nodes <- read_csv('../data_processed/motnp_eigenvector_estimates.csv') %>% 
  select(-model) # column was a single draw from the posterior, not means

# extract mean of posterior for eigenvector
mean_eigen_values <- extract_metric(motnp_fit_edges, "node_eigen") %>%
  as.data.frame() %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = everything()) %>%
  rename(bison_node_eigen=value) %>%
  mutate(node = nodes$node,
         id = nodes$id)

# extract posterior draws for eigenvector
eigen_values <- extract_metric(motnp_fit_edges, "node_eigen") %>%
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'name', values_to = 'eigen') %>% 
  left_join(mean_eigen_values, by = 'name') %>% 
  rename(mean_eigen = bison_node_eigen) %>% 
  select(-name) %>% 
  mutate(draw = rep(1:1000, each = length(unique(id))))

# add mean model values to nodes data frame
nodes <- left_join(nodes, mean_eigen_values[,2:3], by = 'node')
colnames(nodes)[8] <- 'model'

# clean up environment
rm(counts_df_model, df_nodal, mean_eigen_summary, motnp_ages, motnp_edges_null_strongpriors, priors, age, beta, beta_mu, beta_sigma, i, intercept, mean_age) ; gc()
#rm(counts_df_model, gbi_males, m_mat, random_networks, node_ages) ; gc()

#### plot 1) eigenvector from igraph vs eigenvector from bisonR ####
## v1 -- SRI for igraph = basic SRI
ggplot()+
  geom_point(data = nodes, aes(x = igraph, y = model,  colour = age_cat, size = count))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d(direction = -1)+
  scale_x_continuous('eigenvector from putting SRI through igraph')+
  scale_y_continuous('eigenvector from extract_metric(bison_model)')+
  geom_line(aes(x = seq(0.1,1,0.1), y = seq(0.1,1,0.1)))

## v2 -- SRI for igraph = mean of bisonR distribution
edges <- extract_metric(motnp_fit_edges, 'edge_weight') %>% 
  plogis()
n_eles <- length(unique(c(counts_df$id_1, counts_df$id_2)))
edge_means <- data.frame(dyad_id = counts_df$dyad_id,
                         id_1 = counts_df$id_1,
                         id_2 = counts_df$id_2,
                         mean = apply(edges, 2, mean))
sri_means <- matrix(data = NA, nrow = n_eles, ncol = n_eles,
                    dimnames = list(unique(c(counts_df$id_1, counts_df$id_2)),
                                    unique(c(counts_df$id_1, counts_df$id_2))))
for(i in rownames(sri_means)){
  for(j in colnames(sri_means)){
    if(i == j){
      sri_means[i,j] <- 0
    }
    else {
    sri_means[i,j] <- edge_means$mean[i == edge_means$id_1 & j == edge_means$id_2 |
                                        i == edge_means$id_2 & j == edge_means$id_1]
    }
  }
}
g <- graph_from_adjacency_matrix(sri_means, 'undirected', weighted = TRUE)

nodes$igraph_meanedge <- eigen_centrality(g)$vector
ggplot()+
  geom_point(data = nodes, aes(x = igraph_meanedge, y = model,
                               colour = age_cat, size = count))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d(direction = -1)+
  scale_x_continuous('mean eigen centrality using igraph per bisonR draw')+
  scale_y_continuous('mean eigen centrality using bisonR extract_metric()')+
  geom_line(aes(x = seq(0.3,1,0.1), y = seq(0.3,1,0.1)))

## v3 -- SRI for igraph = loop through posterior draws of edge weight
n_draws <- 20 # nrow(edges)
sri_draws <- array(data = NA, dim = c(n_eles, n_eles, n_draws),
                   dimnames = list(unique(c(counts_df$id_1, counts_df$id_2)),
                                   unique(c(counts_df$id_1, counts_df$id_2)),
                                   NULL))
nodes_draws <- matrix(data = NA, ncol = n_eles, nrow = n_draws,
                      dimnames = list(NULL,unique(c(counts_df$id_1, counts_df$id_2))))
for(k in 1:8){
  for(i in rownames(sri_draws)){
    for(j in colnames(sri_draws)){
      if(i == j){
        sri_draws[i,j,k] <- 0
      }
      else {
        sri_draws[i,j,k] <- edges[k, which(counts_df$id_1 == i & counts_df$id_2 == j |
                                              counts_df$id_2 == i & counts_df$id_1 == j)]
      }
    }
  }
  g <- graph_from_adjacency_matrix(sri_draws[,,k], 'undirected', weighted = TRUE)
  #nodes_draws$eigen[nodes_draws$eigen == k] <- eigen_centrality(g)$vector
  nodes_draws[k,] <- eigen_centrality(g)$vector
}
# x <- which(is.na(nodes_draws[,ncol(nodes_draws)]) == TRUE)[1] ; nodes_draws <- nodes_draws[1:x-1,]    # if you have to stop loop above early, go from wherever it reached

# write_csv(as.data.frame(nodes_draws), 'igraph_model_compare.csv') # keeps crashing so save it just in case and then you can just carry on
# nodes_draws <- read_csv('igraph_model_compare.csv')
# rm(mean_eigen_values, mean_motnp_ages, g, motnp_edge_weights_strongpriors,post_eigen) ; gc()

nodes_draws_long <- nodes_draws %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(),
               names_to = 'id', values_to = 'eigenvector') %>% 
  mutate(draw = rep(1:nrow(nodes_draws), each = n_eles))

nodes_draws_long <- left_join(nodes_draws_long,
                              eigen_values[1:(n_draws*n_eles),c('eigen','id','draw')],
                              by = c('id','draw'))
colnames(nodes_draws_long) <- c('id','igraph','draw','model')
nodes_draws_long <- left_join(nodes_draws_long,
                              nodes[,c('id','age_cat','count')],
                              by = 'id')

ggplot()+
  geom_point(data = nodes_draws_long, aes(x = igraph, y = model,
                                          colour = count))+
  #theme(legend.position = 'none')+
  geom_line(aes(x = seq(0.3,1,0.1), y = seq(0.3,1,0.1)))

#### plot 2) observations per individual vs mean group size observed ####
nodes <- read_csv('../data_processed/motnp_eigenvector_estimates.csv') %>% 
  select(-model) # column was a single draw from the posterior, not means

eles_long <- read_csv('../data_processed/motnp_eles_long.csv')
nodes$mean_group_size <- NA
for(i in 1:nrow(nodes)){
  x <- eles_long[eles_long$elephant == nodes$id[i],]
  x <- x[!is.na(x$encounter),]
  nodes$mean_group_size[i] <- mean(x$total_elephants_numeric)
}

ggplot(nodes, aes(x = count, y = mean_group_size))+
  geom_point(aes(#size = igraph,
                 colour = age_cat),
             size = 2)+
  theme(legend.position = 'none')+
  scale_colour_viridis_d(direction = -1)+
  geom_smooth()+
  scale_x_continuous(name = 'total observations per elephant')+
  scale_y_continuous(name = 'mean group size')

# males only
males <- eles_long[eles_long$elephant %in% nodes$id,]
females <- anti_join(eles_long, males)

males$males_idd <- NA
for(i in 1:nrow(males)){
  x <- males[males$encounter == males$encounter[i],]
  x <- x[!is.na(x$encounter),]
  males$males_idd[i] <- nrow(x)
}

nodes$mean_males_per_group <- NA
for(i in 1:nrow(nodes)){
  x <- males[males$elephant == nodes$id[i],]
  x <- x[!is.na(x$encounter),]
  nodes$mean_males_per_group[i] <- mean(x$males_idd)
}

ggplot(nodes, aes(x = count, y = mean_males_per_group))+
  geom_point(aes(#size = igraph,
    colour = age_cat),
    size = 2)+
  theme(legend.position = 'none')+
  scale_colour_viridis_d(direction = -1)+
  geom_smooth()+
  scale_x_continuous(name = 'total observations per elephant')+
  scale_y_continuous(name = 'mean group size')

#### plot 3a) observations per individual vs eigenvector ####
## v1 -- igraph (SRI version)
nodes %>% ggplot()+
  geom_point(aes(x = count, y = igraph, colour = age_cat))+
  scale_colour_viridis_d(direction = -1)+
  theme(legend.position = 'none')+
  scale_x_continuous(name = 'number of observations')+
  scale_y_continuous(name = 'eigenvector centrality from SRI')+
  geom_smooth(aes(x = count, y = igraph))

## v2 -- bisonr
nodes %>% ggplot()+
  geom_point(aes(x = count, y = model, colour = as.factor(mean_age)))+
  scale_colour_viridis_d(direction = -1)+
  theme(legend.position = 'none')+
  scale_x_continuous(name = 'number of observations per elephant')+
  scale_y_continuous(name = 'eigenvector centrality from model\n(mean per elephant)')+
  geom_smooth(aes(x = count, y = model))
eigen_values <- left_join(eigen_values, nodes[,c(1:3,5)], by = c('node','id'))
eigen_values %>% ggplot()+
  geom_point(aes(x = count, y = eigen, colour = age_cat), shape = 19)+
  scale_colour_viridis_d(direction = -1, alpha = 0.2)+
  theme(legend.position = 'none')+
  scale_x_continuous(name = 'number of observations per elephant')+
  scale_y_continuous(name = 'eigenvector centrality from model\n(all draws per elephant)')+
  geom_smooth(aes(x = count, y = mean_eigen))

#### plot 3b) observations per individual vs edge weight ####
## bisonr
edge_means <- left_join(edge_means, counts_df[,c('dyad_id','event_count','count_dyad')],
                        by = 'dyad_id')
edge_means %>% ggplot()+
  geom_point(aes(x = count_dyad, y = mean), colour = alpha('red',0.2))+
  scale_x_continuous(name = 'number of observations of dyad')+
  scale_y_continuous(name = 'mean bisonR-estimated edge weight')+
  geom_smooth(aes(x = count_dyad, y = mean))

edges_long <- edges %>% as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = 'dyad_id', values_to = 'edge_weight')
edges_long$dyad_id <- as.numeric(edges_long$dyad_id)
edges_long <- edges_long %>% 
  left_join(counts_df[,c('dyad_id','event_count','count_dyad')], by = 'dyad_id')
edges_long %>% ggplot()+
  geom_point(aes(x = count_dyad, y = edge_weight), colour = alpha('red',0.02))+
  scale_x_continuous(name = 'number of observations')+
  scale_y_continuous(name = 'bisonR-estimated edge weight')

#### plot 4) low sampling + high group size, bisonR vs SRI ####
## v1 -- plot all and colour by sampling and group size: edge weight
edges <- extract_metric(motnp_fit_edges, 'edge_weight') %>% plogis()
colnames(edges) <- counts_df$dyad_id
counts_df$sri <- counts_df$event_count / counts_df$count_dyad
hist(counts_df$sri)

nodes$mean_edge_weight <- NA
nodes$mean_sri <- NA
for(i in 1:nrow(nodes)){
  x <- counts_df[counts_df$id_1 == nodes$id[i] | 
                   counts_df$id_2 == nodes$id[i],]
  nodes$mean_sri[i] <- mean(x$sri)
  x <- edges[,colnames(edges) %in% x$dyad_id]
  nodes$mean_edge_weight[i] <- mean(x)
}

nodes %>% ggplot()+
  geom_point(aes(x = mean_sri, y = mean_edge_weight,
                 colour = mean_group_size, size = count))+
  scale_colour_viridis_c(option = 'B')+
  theme(legend.position = 'none')+
  geom_line(data = data.frame(x = seq(0.01,0.08,0.01), y = seq(0.01,0.08,0.01)),
            aes(x = x, y = y))+
  scale_x_continuous(name = 'mean edge weight from SRI')+
  scale_y_continuous(name = 'mean egde weight from bisonR')

## v2 -- plot all and colour by bisonR vs size: edge weight
nodes %>% ggplot()+
  geom_point(aes(y = mean_males_per_group, x = count,
                 colour = mean_edge_weight, size = mean_sri))+
  scale_colour_viridis_c(option = 'B')+
  #theme(legend.position = 'none')+
  scale_y_continuous(name = 'average males per group sighting')+
  scale_x_continuous(name = 'number of observations')

## v3 -- subset data for low sampling, high group size: edge weight
nodes %>% filter(count < 10 & mean_males_per_group > 10) %>% 
  ggplot()+
  geom_point(aes(x = mean_sri, y = mean_edge_weight,
                 colour = mean_males_per_group, size = count))+
  scale_colour_viridis_c(option = 'B')+
  #theme(legend.position = 'none')+
  #geom_line(data = data.frame(x = seq(0,0.07,0.01), y = seq(0,0.07,0.01)),
  #aes(x = x, y = y))+
  scale_x_continuous(name = 'eigenvector centrality from SRI')+
  scale_y_continuous(name = 'eigenvector centrality from bisonR')

## v4 -- plot all and colour by sampling and group size: eigenvector
ggplot()+
  geom_point(data = nodes, aes(x = igraph, y = model,
                               colour = mean_males_per_group, size = count))+
  scale_colour_viridis_c(option = 'B')+
  theme(legend.position = 'none')+
  scale_x_continuous(name = 'eigenvector centrality from SRI')+
  scale_y_continuous(name = 'eigenvector centrality from bisonR')+
  geom_line(aes(x = seq(0.25,1,0.05), y = seq(0.25,1,0.05)))

## v5 -- plot all and colour by bisonR vs size: eigenvector
ggplot()+
  geom_point(data = nodes, aes(y = mean_males_per_group, x = count,
                               colour = model, size = igraph))+
  scale_colour_viridis_c(option = 'B')+
  theme(legend.position = 'none')+
  scale_y_continuous(name = 'average males per group sighting')+
  scale_x_continuous(name = 'number of observations')

## v6 -- subset data for low sampling, high group size: eigenvector
nodes %>% filter(count < 10 & mean_males_per_group > 10) %>% 
  ggplot()+
  geom_point(aes(x = igraph, y = model,
                 colour = mean_males_per_group, size = count))+
  scale_colour_viridis_c(option = 'B')+
  #theme(legend.position = 'none')+
  geom_line(data = data.frame(x = seq(0,1,0.1), y = seq(0,1,0.1)), aes(x = x, y = y))+
  scale_x_continuous(name = 'eigenvector centrality from SRI')+
  scale_y_continuous(name = 'eigenvector centrality from bisonR')

#### plot 5) egocentric networks of individuals with extreme eigenvector centrality scores ####
extreme_model <- nodes[nodes$model > quantile(nodes$model, 0.98),]
extreme_igrph <- nodes[nodes$igraph > quantile(nodes$igraph, 0.98),]
xtm <- rbind(extreme_model, extreme_igrph)
rm(extreme_model, extreme_igrph) ; gc()

### generate 2 full networks, one for SRI and one for bisonR
# SRI version
counts_df$sri <- counts_df$event_count / counts_df$count_dyad
n_eles <- length(unique(c(counts_df$id_1, counts_df$id_2)))
sri_mat <- matrix(NA, nrow = n_eles, ncol = n_eles,
                  dimnames = list(unique(c(counts_df$id_1, counts_df$id_2)),
                                  unique(c(counts_df$id_1, counts_df$id_2))))
for(i in rownames(sri_mat)){
  for(j in colnames(sri_mat)){
    if(i == j) { sri_mat[i,j] <- 0 }
    else {
    sri_mat[i,j] <- counts_df$sri[counts_df$id_1 == i & counts_df$id_2 == j |
                                    counts_df$id_1 == j & counts_df$id_2 == i]
    }
  }
}
g_sri <- graph_from_adjacency_matrix(sri_mat, mode = 'undirected', weighted = TRUE)
coords <- igraph::layout_in_circle(g_sri)
id_xtm <- which(nodes$id == xtm$id[1])
coords[id_xtm, ] <- 0
igraph::plot.igraph(g_sri, layout = coords,
                    vertex.label = nodes$id,
                    vertex.label.color = ifelse(nodes$id == xtm$id[1],'black','transparent'),#ifelse(nodes$age_years < 20, 'black', 'white'),
                    vertex.color = ifelse(nodes$age_years < 15, '#FDE725FF',
                                          ifelse(nodes$age_years < 20, '#55C667FF',
                                                 ifelse(nodes$age_years < 30, '#1F968BFF', 
                                                        ifelse(nodes$age_years < 40, '#39568CFF', '#440154FF')))), 
                    vertex.size = 1,#nodes$count,
                    frame.color = NA, frame.width = 0,
                    edge.color = NA, edge.arrow.size = 0, edge.width = 0)
igraph::plot.igraph(g_sri, layout = coords, add = TRUE,
                    vertex.label = NA, vertex.color = 'transparent', vertex.size = 0, 
                    frame.color = NA, frame.width = 0,
                    edge.color = 'black', edge.arrow.size = 0,
                    edge.width = E(g_sri)$weight)

### create egocentric network for each individual in xtm
sri_mat_egocentric <- sri_mat
for(i in rownames(sri_mat_egocentric)){
  for(j in colnames(sri_mat_egocentric)){
    if(i == xtm$id[1] | j == xtm$id[1]){
      sri_mat_egocentric[i,j] <- sri_mat_egocentric[i,j]
    } else { sri_mat_egocentric[i,j] <- 0 }
  }
}
g_ego <- graph_from_adjacency_matrix(sri_mat_egocentric, mode = 'undirected', weighted = T)
coords <- igraph::layout_in_circle(g_ego)
id_xtm <- which(colnames(sri_mat_egocentric) == xtm$id[1])
coords[id_xtm, ] <- 0
igraph::plot.igraph(g_ego, layout = coords,
                    vertex.label = nodes$id,
                    vertex.label.color = ifelse(nodes$id == xtm$id[1],'black','transparent'),#ifelse(nodes$age_years < 20, 'black', 'white'),
                    vertex.color = ifelse(nodes$age_years < 15, '#FDE725FF',
                                          ifelse(nodes$age_years < 20, '#55C667FF',
                                                 ifelse(nodes$age_years < 30, '#1F968BFF', 
                                                        ifelse(nodes$age_years < 40, '#39568CFF', '#440154FF')))), 
                    vertex.size = 1,#nodes$count,
                    frame.color = NA, frame.width = 0,
                    edge.color = NA, edge.arrow.size = 0, edge.width = 0)
igraph::plot.igraph(g_ego, layout = coords, add = TRUE,
                    vertex.label = NA, vertex.color = 'transparent', vertex.size = 0, 
                    frame.color = NA, frame.width = 0,
                    edge.color = 'black', edge.arrow.size = 0,
                    edge.width = E(g_ego)$weight*5)

## repeat using same individual but now with bisonR values --> compare valeus and connections between different methods
counts_df$mod_mean <- NA
weights <- extract_metric(motnp_edge_weights_strongpriors, 'edge_weight') %>% plogis()
for(i in 1:nrow(counts_df)){
  counts_df$mod_mean[i] <- mean(weights[,i])
}
mod_mat <- matrix(NA, nrow = n_eles, ncol = n_eles,
                  dimnames = list(unique(c(counts_df$id_1, counts_df$id_2)),
                                  unique(c(counts_df$id_1, counts_df$id_2))))
for(i in rownames(mod_mat)){
  for(j in colnames(mod_mat)){
    if(i == j) { mod_mat[i,j] <- 0 }
    else {
      mod_mat[i,j] <- counts_df$mod_mean[counts_df$id_1 == i & counts_df$id_2 == j |
                                           counts_df$id_1 == j & counts_df$id_2 == i]
    }
  }
}
mod_mat_egocentric <- mod_mat
for(i in rownames(mod_mat_egocentric)){
  for(j in colnames(mod_mat_egocentric)){
    if(i == xtm$id[1] | j == xtm$id[1]){
      mod_mat_egocentric[i,j] <- mod_mat_egocentric[i,j]
    } else { mod_mat_egocentric[i,j] <- 0 }
  }
}
g_ego_mod <- graph_from_adjacency_matrix(mod_mat_egocentric, mode = 'undirected', weighted = T)
igraph::plot.igraph(g_ego_mod, layout = coords,
                    vertex.label = nodes$id,
                    vertex.label.color = ifelse(nodes$id == xtm$id[1],'black','transparent'),#ifelse(nodes$age_years < 20, 'black', 'white'),
                    vertex.color = ifelse(nodes$age_years < 15, '#FDE725FF',
                                          ifelse(nodes$age_years < 20, '#55C667FF',
                                                 ifelse(nodes$age_years < 30, '#1F968BFF', 
                                                        ifelse(nodes$age_years < 40, '#39568CFF', '#440154FF')))), 
                    vertex.size = 1,#nodes$count,
                    frame.color = NA, frame.width = 0,
                    edge.color = NA, edge.arrow.size = 0, edge.width = 0)
igraph::plot.igraph(g_ego_mod, layout = coords, add = TRUE,
                    vertex.label = NA, vertex.color = 'transparent', vertex.size = 0, 
                    frame.color = NA, frame.width = 0,
                    edge.color = 'black', edge.arrow.size = 0,
                    edge.width = E(g_ego_mod)$weight*5)

### create egocentric network for each individual in xtm
pdf('../outputs/plots_to_check_eigenvector_measures/egocentric_extreme_sri.pdf')
for(elephant in xtm$id){
  sri_mat_egocentric <- sri_mat
  for(i in rownames(sri_mat_egocentric)){
    for(j in colnames(sri_mat_egocentric)){
      if(i == elephant | j == elephant){
        sri_mat_egocentric[i,j] <- sri_mat_egocentric[i,j]
      } else { sri_mat_egocentric[i,j] <- 0 }
    }
  }
  g_ego <- graph_from_adjacency_matrix(sri_mat_egocentric, mode = 'undirected', weighted = T)
  coords <- igraph::layout_in_circle(g_ego)
  id_xtm <- which(colnames(sri_mat_egocentric) == elephant)
  coords[id_xtm, ] <- 0
  igraph::plot.igraph(g_ego, layout = coords,
                      vertex.label = nodes$id,
                      vertex.label.color = ifelse(nodes$id == elephant,'black','transparent'),#ifelse(nodes$age_years < 20, 'black', 'white'),
                      vertex.color = ifelse(nodes$age_years < 15, '#FDE725FF',
                                            ifelse(nodes$age_years < 20, '#55C667FF',
                                                   ifelse(nodes$age_years < 30, '#1F968BFF', 
                                                          ifelse(nodes$age_years < 40, '#39568CFF', '#440154FF')))), 
                      vertex.size = 1,#nodes$count,
                      frame.color = NA, frame.width = 0,
                      edge.color = NA, edge.arrow.size = 0, edge.width = 0)
  igraph::plot.igraph(g_ego, layout = coords, add = TRUE,
                      vertex.label = NA, vertex.color = 'transparent', vertex.size = 0, 
                      frame.color = NA, frame.width = 0,
                      edge.color = 'black', edge.arrow.size = 0,
                      edge.width = E(g_ego)$weight*5)
  
  ## repeat using same individual but now with bisonR values --> compare valeus and connections between different methods
  counts_df$mod_mean <- NA
  weights <- extract_metric(motnp_edge_weights_strongpriors, 'edge_weight') %>% plogis()
  for(i in 1:nrow(counts_df)){
    counts_df$mod_mean[i] <- mean(weights[,i])
  }
  mod_mat <- matrix(NA, nrow = n_eles, ncol = n_eles,
                    dimnames = list(unique(c(counts_df$id_1, counts_df$id_2)),
                                    unique(c(counts_df$id_1, counts_df$id_2))))
  for(i in rownames(mod_mat)){
    for(j in colnames(mod_mat)){
      if(i == j) { mod_mat[i,j] <- 0 }
      else {
        mod_mat[i,j] <- counts_df$mod_mean[counts_df$id_1 == i & counts_df$id_2 == j |
                                             counts_df$id_1 == j & counts_df$id_2 == i]
      }
    }
  }
  mod_mat_egocentric <- mod_mat
  for(i in rownames(mod_mat_egocentric)){
    for(j in colnames(mod_mat_egocentric)){
      if(i == elephant | j == elephant){
        mod_mat_egocentric[i,j] <- mod_mat_egocentric[i,j]
      } else { mod_mat_egocentric[i,j] <- 0 }
    }
  }
  g_ego_mod <- graph_from_adjacency_matrix(mod_mat_egocentric, mode = 'undirected', weighted = T)
  igraph::plot.igraph(g_ego_mod, layout = coords,
                      vertex.label = nodes$id,
                      vertex.label.color = ifelse(nodes$id == elephant,'black','transparent'),#ifelse(nodes$age_years < 20, 'black', 'white'),
                      vertex.color = ifelse(nodes$age_years < 15, '#FDE725FF',
                                            ifelse(nodes$age_years < 20, '#55C667FF',
                                                   ifelse(nodes$age_years < 30, '#1F968BFF', 
                                                          ifelse(nodes$age_years < 40, '#39568CFF', '#440154FF')))), 
                      vertex.size = 1,#nodes$count,
                      frame.color = NA, frame.width = 0,
                      edge.color = NA, edge.arrow.size = 0, edge.width = 0)
  igraph::plot.igraph(g_ego_mod, layout = coords, add = TRUE,
                      vertex.label = NA, vertex.color = 'transparent', vertex.size = 0, 
                      frame.color = NA, frame.width = 0,
                      edge.color = 'black', edge.arrow.size = 0,
                      edge.width = E(g_ego_mod)$weight*5)
  
}
dev.off()











sri_mat_egocentric <- sri_mat_egocentric[which(rowSums(sri_mat_egocentric) > 0),
                                         which(colSums(sri_mat_egocentric) > 0)]
g_ego <- graph_from_adjacency_matrix(sri_mat_egocentric, mode = 'undirected', weighted = T)
coords <- igraph::layout_in_circle(g_ego)
id_xtm <- which(colnames(sri_mat_egocentric) == xtm$id[1])
coords[id_xtm, ] <- 0
igraph::plot.igraph(g_ego, layout = coords,
                    vertex.label = colnames(sri_mat_egocentric),
                    vertex.label.color = ifelse(nodes$age_years < 20, 'black', 'white'), #ifelse(nodes$id == xtm$id[1],'black','transparent'),
                    vertex.color = ifelse(nodes$age_years < 15, '#FDE725FF',
                                          ifelse(nodes$age_years < 20, '#55C667FF',
                                                 ifelse(nodes$age_years < 30, '#1F968BFF', 
                                                        ifelse(nodes$age_years < 40, '#39568CFF', '#440154FF')))), 
                    vertex.size = 20,#nodes$count,
                    frame.color = NA, frame.width = 0,
                    edge.color = NA, edge.arrow.size = 0, edge.width = 0)
igraph::plot.igraph(g_ego, layout = coords, add = TRUE,
                    vertex.label = NA, vertex.color = 'transparent', vertex.size = 0, 
                    frame.color = NA, frame.width = 0,
                    edge.color = 'black', edge.arrow.size = 0,
                    edge.width = E(g_ego)$weight)










net <- make_ego_graph(g_sri, order = 1, nodes = xtm$id[1], mindist = 0)
is_weighted(net[[1]])
plot(net[[1]])

pdf('../outputs/plots_to_check_eigenvector_measures/egocentric_extreme_sri.pdf')
for(i in 1:nrow(xtm)){
  plot_network_threshold(ego_graphs_sri[[i]])
}
dev.off()




# bison version
edges <- extract_metric(motnp_edge_weights_strongpriors, 'edge_weight') %>% 
  plogis()
model_means <- data.frame(dyad_id = counts_df$dyad_id,
                          id_1 = counts_df$id_1,
                          id_2 = counts_df$id_2,
                          mean = apply(edges, 2, mean))
edge_means <- matrix(data = NA, nrow = n_eles, ncol = n_eles,
                     dimnames = list(unique(c(counts_df$id_1, counts_df$id_2)),
                                     unique(c(counts_df$id_1, counts_df$id_2))))
for(i in rownames(edge_means)){
  for(j in colnames(edge_means)){
    if(i == j){
      edge_means[i,j] <- 0
    }
    else {
      edge_means[i,j] <- model_means$mean[i == model_means$id_1 & j == model_means$id_2 |
                                            i == model_means$id_2 & j == model_means$id_1]
    }
  }
}
g_mod <- graph_from_adjacency_matrix(edge_means, 'undirected', weighted = TRUE)



