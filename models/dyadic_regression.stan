data {
  int num_dyads;                               // Number of edges
  int num_nodes;                               // Number of nodes
  vector[num_dyads] logit_edge_mu;             // Means of Gaussian approximation of logit edge weights
  matrix[num_dyads, num_dyads] logit_edge_cov; // Covariance matrix of Gaussian approximation of logit edge weights
  array[num_dyads] real age_min;               // age of younger dyad member
  array[num_dyads] real<lower=0> age_diff;     // difference in age between dyad members
  array[num_dyads] int node_1;                 // Node 1 IDs for multimembership terms
  array[num_dyads] int node_2;                 // Node 2 IDs for multimembership terms
}

parameters {
  real beta_age_diff;
  real beta_age_min;
  vector[num_nodes] mm_nodes;
  real<lower=0> sigma;
  real<lower=0> sigma_mm;
}

transformed parameters {
  //  Cholesky factor of the covariance matrix
  matrix[num_dyads, num_dyads] L_cov;
  L_cov = cholesky_decompose(logit_edge_cov + diag_matrix(rep_vector(sigma, num_dyads)));
  
  // regression equation
  vector[num_dyads] predictor;
  for (i in 1:num_dyads) {
    predictor[i] =  beta_age_min * age_min[i] + beta_age_diff * age_diff[i] + mm_nodes[node_1[i]] + mm_nodes[node_2[i]];
  }
}

model {
  // priors
  beta_age_diff ~ normal(0,0.1);
  beta_age_min ~ normal(0,0.1);
  mm_nodes ~ normal(0, sigma_mm);
  sigma ~ exponential(1);
  sigma_mm ~ exponential(1);

  // likelihood using Cholesky decomposition
  logit_edge_mu ~ multi_normal_cholesky(predictor, L_cov);
}
