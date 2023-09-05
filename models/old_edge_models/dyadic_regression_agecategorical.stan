data {
  int num_dyads;                               // Number of edges
  int num_nodes;                               // Number of nodes
  vector[num_dyads] logit_edge_mu;             // Means of Gaussian approximation of logit edge weights
  matrix[num_dyads, num_dyads] logit_edge_cov; // Covariance matrix of Gaussian approximation of logit edge weights
  int dyad_type[num_dyads];                    // Dyad types of edges
  //int node_1_id[num_dyads];                    // Node 1 IDs for multimembership terms
  //int node_2_id[num_dyads];                    // Node 2 IDs for multimembership terms
}

parameters {
  vector[3] beta_dyadtype;
  //vector[num_nodes] mm_nodes;
  real<lower=0> sigma;
  //real<lower=0> sigma_mm;
}

transformed parameters {
  vector[num_dyads] predictor;
  predictor = beta_dyadtype[dyad_type];// + mm_nodes[node_1_id] + mm_nodes[node_2_id];
}

model {
  logit_edge_mu ~ multi_normal(predictor, logit_edge_cov + diag_matrix(rep_vector(sigma, num_dyads)));
  beta_dyadtype ~ normal(0, 1);
  //mm_nodes ~ normal(0, sigma_mm);
  sigma ~ normal(0, 1);
  //sigma_mm ~ normal(0, 1);
}
