data {
  // global data size
  int num_dyads;                       // Number of edges
  int num_nodes;                       // Number of nodes

  // Gaussian approximation of edge weights
  vector[num_dyads] logit_edge_mu;     // Means of logit edge weights
  vector[num_dyads] logit_edge_sd;     // Standard deviation of logit edge weights
  
  // multimembership terms
  array[num_dyads] int node_1;         // Node 1 IDs for multimembership terms
  array[num_dyads] int node_2;         // Node 2 IDs for multimembership terms
}

parameters {
  // intercept
  real intercept;
  
  // global variance
  real tau_sigma_raw;                  // Unconstrained real value for global scale
  // real<lower=0> sigma;
  
  // multimembership effects
  vector[num_nodes] mm_nodes;
  real mu_mm;
  vector[num_nodes] rand_mm;
  real<lower=0> tau_mm;
  vector[num_nodes] raw_sigma;         // Unconstrained real values for node-specific scales
  
}

transformed parameters {
  // multimembership effects
  vector[num_nodes] node_mean;
  node_mean = mu_mm + rand_mm * tau_mm;
  
  // global sigma
  real tau_sigma;
  tau_sigma = exp(tau_sigma_raw);      // Global scale is positive
  
  // node-level sigma
  vector[num_nodes] node_sigma;
  node_sigma = exp(raw_sigma);         // Node-specific deviations are positive

  // regression equation
  vector[num_dyads] predictor;
  for (i in 1:num_dyads) {
    predictor[i] = intercept + mm_nodes[node_1[i]] + mm_nodes[node_2[i]];
  }
}

model {
  // intercept prior
  intercept ~ normal(0,2);
  
  // variance
  tau_sigma_raw ~ normal(0,2);

  // multimembership priors
  mm_nodes ~ normal(node_mean, node_sigma);
  mu_mm ~ normal(0,0.5);
  rand_mm ~ normal(0,1);
  tau_mm ~ normal(0,0.5);//tau_mm ~ cauchy(0,0.5);
  raw_sigma ~ normal(0,0.5);            // Non-centered parameterization

  // likelihood
  logit_edge_mu ~ normal(predictor, logit_edge_sd + rep_vector(tau_sigma, num_dyads));
}
