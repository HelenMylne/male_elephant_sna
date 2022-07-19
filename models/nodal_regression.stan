data {
  int num_nodes;                               // Number of nodes
  int num_node_types;                          // length(unique(node_age))
  vector[num_nodes] centrality_mu;             // Means of Gaussian approximation of logit edge weights
  matrix[num_nodes, num_nodes] centrality_cov; // Covariance matrix of Gaussian approximation of logit edge weights
  int node_age[num_nodes];                     // Node ages (point estimate)
}

parameters {
  vector[num_node_types] beta_nodetype;
  real<lower=0> sigma;
}

transformed parameters {
  vector[num_nodes] predictor;
  predictor = beta_nodetype[node_age];         // should this be [num_node_types]?
}

model {
  centrality_mu ~ multi_normal(predictor, centrality_cov + diag_matrix(rep_vector(sigma, num_nodes)));
  beta_nodetype ~ normal(0, 1);
  sigma ~ normal(0, 1);
}
