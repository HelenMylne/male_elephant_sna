#### set up ####
# library(tidyverse) ; library(dplyr) ; library(rstan) ; library(rethinking) ; library(igraph) ; library(cmdstanr)
library(cmdstanr)
library(Rcpp)
library(tidyverse, lib.loc = '../packages/')
library(dplyr, lib.loc = '../packages/')
library(igraph, lib.loc = '../packages/')
library(janitor, lib.loc = '../packages/')
library(lubridate, lib.loc = '../packages/')
library(readxl, lib.loc = '../packages/')

# set stan path
#set_cmdstan_path('/Users/helen/.cmdstanr/cmdstan-8.2')

# set seed
set.seed(12345)

# make pdf
pdf('experimenting_with_stan_motnp.pdf')

#### import data frame ####
counts_df <- read_csv('../data_processed/motnp_binomialpairwiseevents_malesonly.csv')

#### run using more complex model where I actually specify the prior properly ####
### compile Stan model
edge_binary <- cmdstan_model("models/edge_binary_basic.stan")
edge_binary

### create data list -- can contain no NA values in any column, even if column is not specified in model
num_dyads <- nrow(counts_df)
num_eles <- length(unique(c(counts_df$id_1, counts_df$id_2)))
counts_ls <- list(
  n_dyads = num_dyads,                 # total number of times one or other of the dyad was observed
  dyad_ids = counts_df$dyad_males,
  together = counts_df$event_count,    # count number of sightings seen together
  count_dyad = counts_df$count_dyad
  )

### Fit model
num_chains <- 4
num_samples <- 1000
fit_weight_motnp <- edge_binary$sample(
  data = counts_ls, 
  chains = num_chains, 
  parallel_chains = num_chains)

# Extract event predictions from the fitted model
event_pred <- rstan::extract(fit_weight_motnp)$event_pred
num_iterations <- dim(event_pred)[1]

# Plot the density of the observed event counts
plot(density(counts_df$event_count), main="", xlab="Dyadic event counts", las=1)

save.image('motnp_edges_fit_experimenting_with_stan.RData')
#load('motnp_edges_fit_experimenting_with_stan.RData')

# Plot the densities of the predicted event counts, repeat for 10 samples
df_copy <- counts_df
for (i in 1:20) {
  df_copy$event <- event_pred[sample(1:num_iterations, size=1), ]
  df_agg_copy <- df_copy %>% 
    group_by(node_1, node_2) %>%
    summarise(event_count=sum(event))
  lines(density(df_agg_copy$event_count), col=rgb(0, 0, 1, 0.5))
}

# Extract edge weights
logit_edge_samples <- rstan::extract(fit_weight_motnp)$logit_edge # Logit scale edge weights -- SHOULD THIS BE LOGIT_EDGE_WEIGHT ??
edge_samples <- plogis(logit_edge_samples) # (0, 1) scale edge weights


# Create edge list
dyad_name <- do.call(paste, c(counts_df[c("node_1", "node_2")], sep=" <-> "))
edge_lower <- apply(edge_samples, 2, function(x) quantile(x, probs=0.025))
edge_upper <- apply(edge_samples, 2, function(x) quantile(x, probs=0.975))
edge_median <- apply(edge_samples, 2, function(x) quantile(x, probs=0.5))
edge_list <- cbind(
  "median"=round(edge_median, 3), 
  "2.5%"=round(edge_lower, 3), 
  "97.5%"=round(edge_upper, 3)
)
rownames(edge_list) <- dyad_name
edge_list

### check model
fit_edges_motnp$summary()
fit_edges_motnp$draws()

# Extract posterior samples
posterior_samples <- fit_edges_motnp$draws()

# Convert the array to a matrix -- save for eigenvector centralities
edge_weights_matrix <- posterior_samples[,,2:(nrow(counts_df)+1)]

# convert matrix to data frame -- save samples and plot outputs
edges <- as.data.frame(edge_weights_matrix[,,1])
colnames(edges) <- c('chain1','chain2','chain3','chain4')
edges <- pivot_longer(edges, everything(), values_to = 'edge_draw', names_to = 'chain')
edges$dyad <- counts_ls$dyad_ids[1]
edges$position <- rep(1:num_samples, each = num_chains)
for(i in 2:num_dyads){
  x <- as.data.frame(edge_weights_matrix[,,i])
  colnames(x) <- c('chain1','chain2','chain3','chain4')
  x <- pivot_longer(x, everything(), values_to = 'edge_draw', names_to = 'chain')
  x$dyad <- counts_ls$dyad_ids[i]
  x$position <- rep(1:num_samples, each = num_chains)
  edges <- rbind(edges, x)
}

### save data 
saveRDS(edges, '../data_processed/motnp_edgedistributions_beta1.5.RDS')
#edges <- readRDS('../data_processed/motnp_edgedistributions_beta1.5.RDS')


#### check outputs -- NETWORK PLOTS USING MY ADAPTED BISONR CODE NOT CURRENTLY WORKING, BUT COME BACK TO THAT ####
# Assign random set of columns to check
if(length(which(counts_df$event_count >= 1)) >= 200){ n_test <- 200 } else { n_test <- length(which(counts_df$event_count >= 1)) }
plot_dyads <- c(sample(counts_df$dyad_males[counts_df$event_count >= 1], size = n_test, replace = F),
                sample(counts_df$dyad_males[counts_df$event_count == 0], size = n_test, replace = F))
plot_edges <- edges[edges$dyad %in% plot_dyads,]
plot_edges$seen_together <- NA
for(i in 1:length(plot_dyads)){
  plot_edges$seen_together[plot_edges$dyad == plot_dyads[i]] <- ifelse(counts_df$event_count[counts_df$dyad_males == plot_dyads[i]] > 0, 1, 0)
}

### build traceplots
ggplot(data = plot_edges[plot_edges$seen_together == 1,], aes(y = edge_draw, x = position, colour = chain))+
  geom_line()+
  facet_wrap(. ~ dyad)+
  theme_classic()+
  theme(legend.position = 'none', strip.background = element_blank(), strip.text = element_blank())
ggplot(data = plot_edges[plot_edges$seen_together == 0,], aes(y = edge_draw, x = position, colour = chain))+
  geom_line()+
  facet_wrap(. ~ dyad)+
  theme_classic()+
  theme(legend.position = 'none', strip.background = element_blank(), strip.text = element_blank())

### density plots
plot(NULL, xlim = c(0,1), ylim = c(0,50), las = 1, xlab = 'edge weight', ylab = 'density')
for(i in 1:length(plot_dyads)){
  x <- plot_edges[plot_edges$dyad == plot_dyads[i],]
  lines(density(x$edge_draw), col = ifelse(x$seen_together == 1, rgb(0,0,1,0.1), rgb(1,0,0,0.1)))
}

### end pdf
dev.off()

### save image
save.image('motnp_edgeweights_stanmodel_beta1.5.RData')

### plot network ####
# Create adjacency array
eles <- unique(c(counts_df$id_1, counts_df$id_2))
num_eles <- length(eles)
adj_tensor <- array(0, c(num_eles, num_eles, num_samples*num_chains),
                    dimnames = list(eles, eles, NULL))
logit_adj_tensor <- array(0, c(num_eles, num_eles, num_samples*num_chains),
                          dimnames = list(eles, eles, NULL))
counts_df$dyad_id_model <- 1:counts_ls$n_dyads
for (dyad_id in 1:counts_ls$n_dyads) {
  dyad_row <- counts_df[counts_df$dyad_id_model == dyad_id, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- edge_samples[, dyad_id]
  adj_tensor[dyad_row$id_2, dyad_row$id_1, ] <- edge_samples[, dyad_id]
  logit_adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- logit_edge_samples[, dyad_id]
  logit_adj_tensor[dyad_row$id_2, dyad_row$id_1, ] <- logit_edge_samples[, dyad_id]
}
adj_tensor[,,1] # Print the first sample of the posterior distribution over adjacency matrices

# plot network
# Calculate lower, median, and upper quantiles of edge weights. Lower and upper give credible intervals.
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]

# Calculate width of credible intervals.
adj_range <- adj_upper - adj_lower
adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one form the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid, mode="undirected", weighted=TRUE)
g_range <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Plot the median graph first and then the standardised width graph to show uncertainty over edges.
coords <- igraph::layout_nicely(g_mid)
plot(g_mid, edge.width=3 * E(g_mid)$weight, edge.color="black",  layout=coords)
plot(g_mid, edge.width=20 * E(g_range)$weight, edge.color=rgb(0, 0, 0, 0.25), 
     vertex.label=eles,vertex.label.color="black", layout=coords, add=TRUE)

# save output
dev.off()
save.image('experimenting_with_stan_motnp.RData')
