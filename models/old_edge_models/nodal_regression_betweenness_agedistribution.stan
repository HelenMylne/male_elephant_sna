data {
  int num_nodes;                               // Number of nodes
  //int num_node_types;                        // length(unique(node_age))
  vector[num_nodes] centrality_mu;             // Means of Gaussian approximation of logit edge weights
  matrix[num_nodes, num_nodes] centrality_cov; // Covariance matrix of Gaussian approximation of logit edge weights
  matrix[num_nodes, 8000] node_age;            // Node ages (draws from posterior distribution)
  //matrix[num_nodes, 8000] node_age2;           // Node ages squared (same as above, all estimates squared)
}

parameters {
  real beta_age;
  //real beta_age2;
  real intercept;
  real<lower=0> sigma;
}

transformed parameters {
  vector[num_nodes] predictor;
  for(i in 1:num_nodes){
    predictor = intercept + beta_age*node_age[,i];// + beta_age2*node_age2[,i];
  }
}

model {
  centrality_mu ~ multi_normal(predictor, centrality_cov + diag_matrix(rep_vector(square(sigma), num_nodes)));
  intercept ~ normal(2, 1.4);
  beta_age ~ normal(0, 0.02);
  //beta_age2 ~ normal(0, 0.001);
  sigma ~ normal(0, 1);
}
