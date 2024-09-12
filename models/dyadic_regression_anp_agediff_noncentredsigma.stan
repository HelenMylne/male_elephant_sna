data {
  // global data size
  int num_data;                       // Number of edges in total data
  int num_dyads;                      // Number of unique dyads across all windows
  int num_nodes;                      // Number of unique nodes across all windows
  int num_windows;                    // Number of time windows
  
  // Gaussian approximation of edge weights
  vector[num_data] logit_edge_mu;     // Means of Gaussian approximation of logit edge weights
  vector[num_data] logit_edge_sd;     // SD of Gaussian approximation of logit edge weights
  
  // explanatory variable
  array[num_data] real age_diff;      // age difference between dyad members
  
  // random effects
  array[num_data] int window;         // window ID for random effect
  
  // multimembership terms
  array[num_data] int node_1;         // Node 1 IDs for multimembership terms
  array[num_data] int node_2;         // Node 2 IDs for multimembership terms
  
}

parameters {
  // intercept
  real intercept;
  
  // exposure slope
  real beta_age_diff;
  
  // global variance
  real tau_sigma_raw;            // Unconstrained real value for global scale
  // real<lower=0> sigma;
  
  // // multimembership effects
  // real mu_mm;
  // vector[num_nodes] rand_mm;
  // real<lower=0> sigma_mm;
  
  // multimembership effects
  vector[num_nodes] mm_nodes;
  real mu_mm;
  vector[num_nodes] rand_mm;
  real<lower=0> tau_mm;
  // vector<lower=0>[num_nodes] node_sigma;
  // real theta_node;
  vector[num_nodes] raw_sigma;   // Unconstrained real values for node-specific scales
  
  // random effect of time window
  real mu_window;
  vector[num_windows] rand_window;
  real<lower=0> sigma_window;
  
}

transformed parameters {
  // random effects
  vector[num_windows] window_random_effect;
  window_random_effect = mu_window + rand_window * sigma_window;
  
  // // multimembership effects
  // vector[num_nodes] mm_nodes;
  // mm_nodes = mu_mm + rand_mm * sigma_mm;
  
  // multimembership effects
  vector[num_nodes] node_mean;
  node_mean = mu_mm + rand_mm * tau_mm;
  
  // global sigma
  real tau_sigma;
  tau_sigma = exp(tau_sigma_raw);      // Global scale is positive
  
  // node-level sigma
  vector[num_nodes] node_sigma;
  node_sigma = exp(raw_sigma);         // Node-specific deviations are positive

  
  // linear model
  vector[num_data] predictor;
  for(i in 1:num_data) {
    predictor[i] = intercept + beta_age_diff * age_diff[i] + mm_nodes[node_1[i]] + mm_nodes[node_2[i]] + window_random_effect[window[i]];
  }
}

model {
  // intercept prior
  intercept ~ normal(0,0.5);//2);

  // slope priors
  beta_age_diff ~ normal(0,1);

  // variance
  tau_sigma_raw ~ normal(0, 2);        //sigma ~ exponential(2); --> Prior on global scale before transformation

  // // multimembership priors
  // mu_mm ~ normal(0,0.5);
  // rand_mm ~ normal(0,1);
  // sigma_mm ~ exponential(2);//normal(0,1)T[0,];
  
  // multimembership priors
  mm_nodes ~ normal(node_mean, node_sigma);
  mu_mm ~ normal(0,0.5);
  rand_mm ~ normal(0,1);
  tau_mm ~ cauchy(0,2.5);
  // node_sigma ~ exponential(theta_node);
  // theta_node ~ normal(0,2);
  raw_sigma ~ normal(0, 1);            // Non-centered parameterization

  // random effect of window
  mu_window ~ normal(0,0.4);//1);
  rand_window ~ normal(0,1);
  sigma_window ~ exponential(2);//1);

  // likelihood
  logit_edge_mu ~ normal(predictor, logit_edge_sd + rep_vector(tau_sigma, num_data));//sigma, num_data));
  
}
