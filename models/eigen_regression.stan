data {
  int num_nodes;                    // Number of nodes
  vector[num_nodes] nodes;          // Node IDs
  vector[num_nodes] centrality_mu;  // Means of centrality estimates
  vector[num_nodes] centrality_sd;  // standard deviations of centrality estimates
  vector[num_nodes] node_age;       // Node ages (point estimate)
}

parameters {
  real beta_age;
  real<lower=0> sigma;
}

transformed parameters {
  vector[num_nodes] predictor;
  predictor = beta_age*node_age;
}

model {
  // priors
  beta_age ~ normal(0, 0.1);
  
  // likelihood
  logit(centrality_mu) ~ normal(predictor, centrality_sd);

  
}
