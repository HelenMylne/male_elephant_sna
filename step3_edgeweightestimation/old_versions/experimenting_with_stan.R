#### set up ####
# library(tidyverse) ; library(dplyr) ; library(rstan) ; library(rethinking) ; library(igraph) ; library(cmdstanr)
library(tidyverse, lib.loc = 'packages/')
library(dplyr, lib.loc = 'packages/')
#library(rstan, lib.loc = 'packages/')
library(rethinking, lib.loc = 'packages/')
library(igraph, lib.loc = 'packages/')
library(cmdstanr, lib.loc = 'packages/')

# set stan path
#set_cmdstan_path('/Users/helen/.cmdstanr/cmdstan-8.2')

# set seed
set.seed(12345)

#### import data frame ####
counts_df <- read_csv('../data_processed/motnp_binomialpairwiseevents_malesonly.csv')

#### run using more complex model where I actually specify the prior properly ####
### compile Stan model
edge_binary <- cmdstan_model("models/edge_binary_basic.stan")
edge_binary

### create data list -- can contain no NA values in any column, even if column is not specified in model
counts_ls <- list(
  n_dyads = nrow(counts_df),          # total number of times one or other of the dyad was observed
  dyad_ids = counts_df$dyad_males,
  together = counts_df$event_count,    # count number of sightings seen together
  #apart = counts_df$apart,          # count number of sightings seen apart
  #seen_together = ifelse(counts_df$event_count > 0, 1, 0),
  #mu = ifelse(counts_df$event_count > 0, -1.5, -3),
  #sigma = ifelse(counts_df$event_count > 0, 0.8, 1),
  count_dyad = counts_df$count_dyad
  )

### Fit model
num_samples <- 1000
num_chains <- 4
num_dyads <- nrow(counts_df)
fit_weight_motnp <- edge_binary$sample(
  data = counts_ls, 
  chains = num_chains, 
  parallel_chains = num_chains)

# Extract event predictions from the fitted model
event_pred <- rstan::extract(fit_weight_motnp)$event_pred
num_iterations <- dim(event_pred)[1]

# Plot the density of the observed event counts
plot(density(counts_df$event_count), main="", xlab="Dyadic event counts", las=1)

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

# Extract posterior samples
posterior_samples <- fit_weight_motnp$draws()

# Convert the array to a matrix -- save for eigenvector centralities
total_samples <- n_chains*n_samples
edge_weights_matrix <- posterior_samples[,,2:(length(posterior_samples)/(total_samples))]

# convert matrix to data frame -- save samples and plot outputs
edges <- as.data.frame(edge_weights_matrix[,,1])
colnames(edges) <- c('chain1','chain2','chain3','chain4')
edges <- pivot_longer(edges, everything(), values_to = 'edge_draw', names_to = 'chain')
edges$dyad <- counts_ls$dyad_ids[1]
edges$position <- rep(1:n_samples, each = n_chains)
for(i in 2:num_dyads){
  x <- as.data.frame(edge_weights_matrix[,,i])
  colnames(x) <- c('chain1','chain2','chain3','chain4')
  x <- pivot_longer(x, everything(), values_to = 'edge_draw', names_to = 'chain')
  x$dyad <- counts_ls$dyad_ids[i]
  x$position <- rep(1:n_samples, each = n_chains)
  edges <- rbind(edges, x)
}

### save data 
saveRDS(edges, '../data_processed/motnp_edgedistributions.RDS')

# Plot the densities of the association strengths
sample_seen <- sample(counts_df$dyad_id[counts_df$event_count > 0], 500, replace = F)
sample_unsn <- sample(counts_df$dyad_id[counts_df$event_count ==0], 500, replace = F)
counts_df$sri <- counts_df$event_count / counts_df$count_dyad








colnames(edge_samples) <- counts_df$dyad_id









plot(NULL, main="", xlab="edge weight", ylab="density", las=1, xlim=c(0,1), ylim=c(0,50))
for (i in 1:length(sample_seen)) {
  lines(density(edge_samples[,colnames(edge_samples) == sample_seen[i]]), col=rgb(0, 0, 1, 0.05))
  lines(density(edge_samples[,colnames(edge_samples) == sample_unsn[i]]), col=rgb(1, 0, 0, 0.05))
}
lines(density(counts_df$sri[counts_df$event_count > 0]), lwd = 2, col = 'green')
lines(density(counts_df$sri), lwd = 2)

### check traces
par(mfrow = c(8,8))
for (i in 1:32) {
  plot(edge_samples[,colnames(edge_samples) == sample_seen[i]], type = 'l', col=rgb(0, 0, 1, 0.05))
  plot(edge_samples[,colnames(edge_samples) == sample_unsn[i]], type = 'l', col=rgb(1, 0, 0, 0.05))
}
par(mfrow = c(1,1))

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
