data {
  // global data size
  int num_data;                            // Number of edges in total data
  int num_dyads;                           // Number of unique dyads across all windows
  int num_nodes;                           // Number of unique nodes across all windows
  int num_windows;                         // Number of time windows
  int max_age_diff;  // maximum possible age difference between dyad
  
  // Gaussian approximation of edge weights
  vector[num_data] logit_edge_mu;          // Means of Gaussian approximation of logit edge weights
  vector[num_data] logit_edge_sd;          // SD of Gaussian approximation of logit edge weights
  
  // explanatory variable
  array[num_data] int age_diff_plus1;           // age difference between dyad members
  
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
  
  // random effect correlations
  vector[max_age_diff] diff;
  matrix[max_age_diff,num_windows] z_window;
  vector[max_age_diff] tau_window;
  cholesky_factor_corr[max_age_diff] rho_window;
  
  matrix[max_age_diff,num_nodes] z_nodes;
  vector[max_age_diff] tau_nodes;
  cholesky_factor_corr[max_age_diff] rho_nodes;
  
}
transformed parameters {
  // global sigma
  real tau_sigma;
  tau_sigma = exp(tau_sigma_raw);          // Global scale is positive
  
  // correlation matrix of random effects
  vector<lower=0>[max_age_diff] sigma_window;
  sigma_window = exp(tau_window);
  vector<lower=0>[max_age_diff] sigma_nodes;
  sigma_nodes = exp(tau_nodes);
  
  matrix[num_windows,max_age_diff] window_effect;
  window_effect = diag_pre_multiply(sigma_window, rho_window) * z_window;
  matrix[num_nodes,max_age_diff] mm_nodes;
  mm_nodes = diag_pre_multiply(sigma_nodes, rho_nodes) * z_nodes;
  
  // linear model
  vector[num_data] predictor;
  for(i in 1:num_data) {
    predictor[i] = beta_age_diff * age_diff_plus1[i] + mm_nodes[node_1[i],age_diff_plus1[i]] + mm_nodes[node_2[i],age_diff_plus1[i]] + window_effect[window[i],age_diff_plus1[i]];
  }
}
model {
  // slope priors
  beta_age_diff ~ normal(0,1);

  // variance
  tau_sigma_raw ~ normal(0,2);

  // correlation matrix of random effects
  rho_window ~ lkj_corr_cholesky( 2 );
  tau_window ~ normal(0,1);
  rho_nodes ~ lkj_corr_cholesky( 2 );
  tau_nodes ~ normal(0,1);
  
  to_vector(z_window) ~ normal(0,1);
  to_vector(z_nodes) ~ normal(0,1);

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




  // // multimembership effects
  // vector[num_nodes] mm_nodes;
  // real mu_mm;
  // vector[num_nodes] rand_mm;
  // real<lower=0> tau_mm;
  // vector[num_nodes] sigma_mm;              // Unconstrained real values for node-specific scales
  
  // // random effect of time window
  // real mu_window;
  // vector[num_windows] rand_window;
  // real<lower=0> sigma_window;
  



  // // random effects
  // vector[num_windows] window_random_effect;
  // window_random_effect = mu_window + rand_window * sigma_window;
  
  // // multimembership effects
  // vector[num_nodes] node_mean;
  // node_mean = mu_mm + rand_mm * tau_mm;
  
  // // node-level sigma
  // vector[num_nodes] node_sigma;
  // node_sigma = exp(sigma_mm);              // Node-specific deviations are positive




  // // multimembership priors
  // mm_nodes ~ normal(node_mean, node_sigma);
  // mu_mm ~ normal(0,0.5);
  // rand_mm ~ normal(0,1);
  // tau_mm ~ normal(0,0.5);//tau_mm ~ cauchy(0,0.5);
  // sigma_mm ~ normal(0,1);

  // // random effect of window
  // mu_window ~ normal(0,0.4);//1);
  // rand_window ~ normal(0,1);
  // sigma_window ~ exponential(1);
  


