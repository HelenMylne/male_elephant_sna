data {
  int num_dyads; // Number of edges
  int num_nodes; // Number of nodes
  vector[num_dyads] logit_edge_mu;             // Means of Gaussian approximation of logit edge weights
  matrix[num_dyads, num_dyads] logit_edge_cov; // Covariance matrix of Gaussian approximation of logit edge weights
  vector[num_dyads] age_diff; // Dyad difference between mean ages (not mean difference between ages)
  vector[num_dyads] age_1;    // Node 1 mean age
  vector[num_dyads] age_2;    // Node 2 mean age
  //int node_1_id[num_dyads];   // Node 1 IDs for multimembership terms
  //int node_2_id[num_dyads];   // Node 2 IDs for multimembership terms
}

parameters {
  real beta_agediff;
  real beta_age;
  //vector[num_nodes] mm_nodes;
  real<lower=0> sigma;
  //real<lower=0> sigma_mm;
}

transformed parameters {
  vector[num_dyads] age_effect;
  age_effect = beta_agediff*age_diff + beta_age*age_1 + beta_age*age_2;// + mm_nodes[node_1_id] + mm_nodes[node_2_id];
}

model {
  logit_edge_mu ~ multi_normal(age_effect, logit_edge_cov + diag_matrix(rep_vector(sigma, num_dyads)));
  beta_agediff ~ normal(0, 1);
  beta_age ~ normal(0, 1);
  //mm_nodes ~ normal(0, sigma_mm);
  sigma ~ normal(0, 1);
  //sigma_mm ~ normal(0, 1);
}
