data {
  int num_nodes;                               // Number of nodes
  //int num_node_types;                        // length(unique(node_age))
  vector[num_nodes] centrality_mu;             // Means of Gaussian approximation of logit edge weights
  matrix[num_nodes, num_nodes] centrality_cov; // Covariance matrix of Gaussian approximation of logit edge weights
  vector[num_nodes] node_age;                  // Node ages (point estimate)
  vector[num_nodes] node_age2;                  // Node ages (point estimate)
}

parameters {
  real beta_age;
  real beta_age2;
  real intercept;
  real<lower=0> sigma;
}

transformed parameters {
  vector[num_nodes] predictor;
  predictor = intercept + beta_age*node_age + beta_age2*node_age2;
}

model {
  centrality_mu ~ multi_normal(predictor, centrality_cov + diag_matrix(rep_vector(square(sigma), num_nodes)));
  intercept ~ normal(0, 1);
  beta_age ~ normal(0, 0.02);
  beta_age2 ~ normal(0, 0.003);
  sigma ~ normal(0, 1);
}
