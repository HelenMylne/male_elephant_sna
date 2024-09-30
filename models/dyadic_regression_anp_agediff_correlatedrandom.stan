data {
  // global data size
  int num_data;                            // Number of edges in total data
  int num_dyads;                           // Number of unique dyads across all windows
  int num_nodes;                           // Number of unique nodes across all windows
  int num_windows;                         // Number of time windows
  
  // Gaussian approximation of edge weights
  vector[num_data] logit_edge_mu;          // Means of Gaussian approximation of logit edge weights
  vector[num_data] logit_edge_sd;          // SD of Gaussian approximation of logit edge weights
  
  // explanatory variable
  array[num_data] real age_diff;           // age difference between dyad members
  
  // random effects
  array[num_data] int window;              // window ID for random effect
  
  // multimembership terms
  array[num_data] int node_1;              // Node 1 IDs for multimembership terms
  array[num_data] int node_2;              // Node 2 IDs for multimembership terms
}
parameters {
  // exposure slope
  real beta_age_diff;
  
  // global variance
  real tau_sigma_raw;                      // Unconstrained real value for global scale
  
  // multimembership effects
  vector[num_nodes] mm_nodes;
  real mu_mm;
  vector[num_nodes] rand_mm;
  real<lower=0> tau_mm;
  vector[num_nodes] sigma_mm;              // Unconstrained real values for node-specific scales
  
  // random effect of time window
  real mu_window;
  vector[num_windows] rand_window;
  real<lower=0> sigma_window;
  
  // random effect correlations
  vector<lower=0>[2] sigma_rand;           // random effects standard deviations
  cholesky_factor_corr[2] chol_corr_rand;  // Cholesky factor of correlation matrix
  matrix[num_nodes,num_windows] mtx_rand;  // random effect matrix
}
transformed parameters {
  // random effects
  vector[num_windows] window_random_effect;
  window_random_effect = mu_window + rand_window * sigma_window;
  
  // multimembership effects
  vector[num_nodes] node_mean;
  node_mean = mu_mm + rand_mm * tau_mm;
  
  // global sigma
  real tau_sigma;
  tau_sigma = exp(tau_sigma_raw);          // Global scale is positive
  
  // node-level sigma
  vector[num_nodes] node_sigma;
  node_sigma = exp(sigma_mm);              // Node-specific deviations are positive

  // correlation matrix of random effects
  matrix[num_nodes,num_windows] mtx_corr_rand;
  mtx_corr_rand = diag_pre_multiply(sigma_rand, chol_corr_rand) * mtx_rand;
  
  // linear model
  vector[num_data] predictor;
  for(i in 1:num_data) {
    predictor[i] = beta_age_diff * age_diff[i] + mm_nodes[node_1[i]] + mm_nodes[node_2[i]] + window_random_effect[window[i]] + mtx_corr_rand[node_1[i],window[i]] + mtx_corr_rand[node_2[i],window[i]];
  }
}
model {
  // slope priors
  beta_age_diff ~ normal(0,1);

  // variance
  tau_sigma_raw ~ normal(0,2);

  // multimembership priors
  mm_nodes ~ normal(node_mean, node_sigma);
  mu_mm ~ normal(0,0.5);
  rand_mm ~ normal(0,1);
  tau_mm ~ normal(0,0.5);//tau_mm ~ cauchy(0,0.5);
  sigma_mm ~ normal(0,1);

  // random effect of window
  mu_window ~ normal(0,0.4);//1);
  rand_window ~ normal(0,1);
  sigma_window ~ exponential(1);
  
  // correlation matrix of random effects
  chol_corr_rand ~ lkj_corr_cholesky(1.5); // LKJ prior for the correlation matrix
  to_vector(mtx_rand) ~ normal(0,2);
  sigma_rand ~ exponential(1);

  // likelihood
  logit_edge_mu ~ normal(predictor, logit_edge_sd + rep_vector(tau_sigma, num_data));
  
}

// generated quantities {
//   matrix[num_nodes, num_nodes] Omega;
//   Omega = chol_corr_rand * chol_corr_rand';           // reverse Cholesky to return correlation matrix between window and ID
//   
//   // // list of resources for building model with correlated random effects
//   // https://mlisi.xyz/post/bayesian-multilevel-models-r-stan/
// }


