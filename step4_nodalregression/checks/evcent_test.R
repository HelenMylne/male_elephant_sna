# check sna and igraph produce the same things just without standardising in sna
library(tidyverse)
library(igraph)
library(sna)

## load data and remove additional data
load('anp_nodalregression/anpshort1_nodalregression_conditionaledge.RData')
rm(counts_df, adj_mat, edges, ele_obs, obs, summary, make_edgelist, plot_network_threshold_anp, dyad_row, eigen_values, network, eigen, i, n_windows, periods) ; gc()

## build adjacency matrix
ncol(edge_samples) # 1326
cdf_1$node_1_id <- as.integer(as.factor(cdf_1$node_1))
cdf_1$node_2_id <- as.integer(as.factor(cdf_1$node_2))+1
cdf_1$dyad_rank <- as.integer(as.factor(cdf_1$dyad_id))

adj_mat <- matrix(NA, nrow = n_eles, ncol = n_eles,
                  dimnames = list(ele_ids, ele_ids))
colnames(edge_samples) <- cdf_1$dyad_rank
for(i in 1:n_eles){
  for(j in 1:n_eles){
    if(i < j) {       # add eigen scores for top triangle as normal
      dyad_rank <- cdf_1$dyad_rank[which(cdf_1$id_1 == colnames(adj_mat)[i] &
                                           cdf_1$id_2 == rownames(adj_mat)[j] )]
      adj_mat[i,j] <- edge_samples[1,dyad_rank]
    } else {
      if(i > j) {     # also fill in bottom half -- sna::evcent() doesn't work for undirected graph if matrix is not symmetrical
        dyad_rank <- cdf_1$dyad_rank[which(cdf_1$id_1 == colnames(adj_mat)[j] &
                                             cdf_1$id_2 == rownames(adj_mat)[i] )]
        adj_mat[i,j] <- edge_samples[1,dyad_rank]
      } else {
        adj_mat[i,j] <- 0
      }
    }
  }
}

## extract eigenvector using igraph
g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE)
eigen <- igraph::eigen_centrality(g)$vector %>% 
  as.data.frame()
colnames(eigen) <- 'igraph'
eigen$id <- rownames(eigen)

## extract eigenvector using SNA
eigen$sna <- sna::evcent(adj_mat, gmode = 'graph', diag = F)
hist(eigen$sna)

## compare SNA and igraph outputs
eigen <- eigen %>% 
  relocate(id) %>% 
  mutate(sna_std = sna / max(sna)) %>% 
  mutate(check = ifelse(round(igraph,2) == round(sna_std,2),
                        'match','different'))
table(eigen$check) # all are the same, therefore I can use sna::evcent instead of igraph to obtain eigenvector centralities, which avoids the issues of logit(1) = Inf caused by igraph always giving one elephant a centrality of 1
