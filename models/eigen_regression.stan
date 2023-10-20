data {
  int num_nodes;                    // Number of nodes
  vector[num_nodes] nodes;          // Node IDs
  vector[num_nodes] centrality_mu;  // Means of centrality estimates
  matrix[num_nodes, num_nodes] centrality_cov;  // standard deviations of centrality estimates
  vector[num_nodes] node_age;       // Node ages (point estimate)
}

parameters {
  real beta_age;
  real<lower=0> sigma;
}

transformed parameters {
  // linear model
  vector[num_nodes] predictor;
  predictor = beta_age*node_age;
}

model {
  // priors
  beta_age ~ normal(0, 0.1);
  sigma ~ normal(0, 1); // exponential(1)??
  
  // likelihood
  centrality_mu ~ multi_normal(predictor, centrality_cov + diag_matrix(rep_vector(sigma, num_nodes)));
  
}
