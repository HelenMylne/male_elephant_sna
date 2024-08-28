data {
  // global data size
  int num_data;                                   // Number of edges in total data
  int num_nodes;                                  // Number of unique nodes across all windows
  int num_windows;                                // Number of time windows
  
  // Gaussian approximation of edge weights
  vector[num_data] logit_edge_mu;                 // Means of Gaussian approximation of logit edge weights
  vector[num_data] logit_edge_sd;                 // SD of Gaussian approximation of logit edge weights
  
  // explanatory variable
  array[num_data] real age_mean;                  // mean age of dyad members
  
  // time window
  array[num_data] int window;                     // window ID
  
  // multimembership terms
  array[num_data] int node_1;                     // Node 1 IDs for multimembership terms
  array[num_data] int node_2;                     // Node 2 IDs for multimembership terms
  
}

parameters {
  // intercept
  real intercept;
  
  // exposure slope
  real beta_age_mean;
  
  // variance
  real<lower=0> sigma;
  
  // multimembership effects
  real mu_mm;
  vector[num_nodes] rand_mm;
  real<lower=0> sigma_mm;
  
  // random effect of time window
  vector[num_windows] rand_window;
  real<lower=0> sigma_window;
}

transformed parameters {
  // random effects
  vector[num_windows] window_random_effect;
  window_random_effect = rand_window * sigma_window;

  // multimembership effects
  vector[num_nodes] mm_nodes;
  mm_nodes = mu_mm + rand_mm * sigma_mm;
  
  // linear model
  vector[num_data] predictor;
  for(i in 1:num_data) {
    predictor[i] = intercept + beta_age_mean * age_mean[i] + mm_nodes[node_1[i]] + mm_nodes[node_2[i]] + window_random_effect[window[i]];
  }
}

model {
  // intercept prior
  intercept ~ normal(0,2);
  
  // slope priors
  beta_age_mean ~ normal(0,1);
  
  // variance
  sigma ~ exponential(2);
  
  // multimembership priors
  mu_mm ~ normal(0,1);
  rand_mm ~ normal(0,1);
  sigma_mm ~ normal(0,1)T[0,];
  
  // random effect of window
  rand_window ~ normal(0,1);
  sigma_window ~ exponential(1);
  
  // likelihood
  logit_edge_mu ~ normal(predictor, logit_edge_sd + rep_vector(sigma, num_data));
}
