data {
  int num_nodes;                    // Number of nodes
  vector[num_nodes] nodes;          // Node IDs
  vector[num_nodes] centrality_mu;  // Means of centrality estimates (0-1 bounded then standardised)
  matrix[num_nodes, num_nodes] centrality_cov;  // standard deviations of centrality estimates
  vector[num_nodes] node_age;       // Node ages (point estimate)
}

parameters {
  real beta_age;
  real<lower=0> sigma;
  real intercept;
}

transformed parameters {
  // linear model
  vector[num_nodes] predictor;
  predictor = intercept + beta_age*node_age;
}

model {
  // priors
  beta_age ~ normal(0, 0.5);
  sigma ~ exponential(2);
  intercept ~ normal(0,5);
  
  // likelihood
  centrality_mu ~ multi_normal(predictor, centrality_cov + diag_matrix(rep_vector(sigma, num_nodes)));
}
