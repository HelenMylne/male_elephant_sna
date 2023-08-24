data {
  int num_dyads;                               // Number of edges
  int num_nodes;                               // Number of nodes
  vector[num_dyads] logit_edge_mu;             // Means of Gaussian approximation of logit edge weights
  matrix[num_dyads, num_dyads] logit_edge_cov; // Covariance matrix of Gaussian approximation of logit edge weights
  real age_min[num_dyads];                     // age of younger dyad member
  real<lower=0> age_diff[num_dyads];           // difference in age between dyad members
  int node_1[num_dyads];                       // Node 1 IDs for multimembership terms
  int node_2[num_dyads];                       // Node 2 IDs for multimembership terms
  //real jitter;                                 // jitter to add to the diag of the cov matrix for numerical stability -- TRY TO REMOVE WHEN YOU PUT SIGMA BACK IN PLACE
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
