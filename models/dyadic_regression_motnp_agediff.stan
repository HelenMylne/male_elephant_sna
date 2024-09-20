data {
  // global data size
  int num_dyads;                       // Number of edges
  int num_nodes;                       // Number of nodes
  int num_age_diffs;                   // Number of unique age difference categories
  int length_dirichlet;                // Number of unique age difference categories + 1
  
  // Gaussian approximation of edge weights
  vector[num_dyads] logit_edge_mu;     // Means of logit edge weights
  vector[num_dyads] logit_edge_sd;     // Standard deviation of logit edge weights
  
  // explanatory variable
  array[num_dyads] int age_diff;       // ordinal categorical where 1 = same age, 2 = 1 category different, 3 = 2 categories different... etc.
  
  // multimembership terms
  array[num_dyads] int node_1;         // Node 1 IDs for multimembership terms
  array[num_dyads] int node_2;         // Node 2 IDs for multimembership terms
  
  // prior values for Dirichlet
  vector[num_age_diffs] prior_diff;      // Dirichlet prior values
}

parameters {
  // // intercept
  // real intercept;
  
  // exposure slope
  real beta_age_diff;
  
  // global variance
  real tau_sigma_raw;                  // Unconstrained real value for global scale
  // real<lower=0> sigma;
  
  // multimembership effects
  vector[num_nodes] mm_nodes;
  real mu_mm;
  vector[num_nodes] rand_mm;
  real<lower=0> tau_mm;
  vector[num_nodes] raw_sigma;         // Unconstrained real values for node-specific scales
  
  // difference between age categories
  simplex[num_age_diffs] delta_diff;
}

transformed parameters {
  // create prior for cumulative probability of each age category
  vector[length_dirichlet] delta_j_diff;
  delta_j_diff = append_row(0, delta_diff);
  
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
    predictor[i] = beta_age_diff * sum(delta_j_diff[1:age_diff[i]]) + mm_nodes[node_1[i]] + mm_nodes[node_2[i]]; // + intercept;
  }
}

model {
  // // intercept prior
  // intercept ~ normal(0,2);
  
  // age priors
  beta_age_diff ~ normal(0,3);
  delta_diff ~ dirichlet(prior_diff);
  
  // variance
  tau_sigma_raw ~ normal(0,1);

  // multimembership priors
  mm_nodes ~ normal(node_mean, node_sigma);
  mu_mm ~ normal(0,1);//normal(logit(0.1),1);
  rand_mm ~ normal(0,1);
  tau_mm ~ normal(0,1);
  raw_sigma ~ normal(0,0.5);            // Non-centered parameterization

  // likelihood
  logit_edge_mu ~ normal(predictor, logit_edge_sd + rep_vector(tau_sigma, num_dyads));
}
