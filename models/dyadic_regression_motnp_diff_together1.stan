data {
  // global data size
  int num_dyads;                       // Number of edges
  int num_nodes;                       // Number of nodes
  int num_age_diffs;                   // Number of unique age difference categories
  int length_dirichlet;                // Number of unique age difference categories + 1
  
  // explanatory variables
  array[num_dyads] int age_diff;       // ordinal categorical where 1 = same age, 2 = 1 category different, 3 = 2 categories different... etc.
  
  // observed variable
  array[num_dyads] int together;       // observed values of were they ever seen together
  
  // multimembership terms
  array[num_dyads] int node_1;         // Node 1 IDs for multimembership terms
  array[num_dyads] int node_2;         // Node 2 IDs for multimembership terms
  
  // prior values for Dirichlet
  vector[num_age_diffs] prior_diff;    // Dirichlet prior values
}

parameters {
  // exposure slopes
  real beta_diff_together;             // effect of age on Pr(ever together)
  
  // global variance
  real tau_sigma_raw;                  // Unconstrained real value for global scale
  
  // multimembership effects
  vector[num_nodes] mm_nodes_togt;
  real mu_mm_togt;
  vector[num_nodes] rand_mm_togt;
  real<lower=0> tau_mm_togt;
  vector[num_nodes] raw_sigma_togt;         // Unconstrained real values for node-specific scales
  
  // difference between age categories
  simplex[num_age_diffs] delta_diff_together;
}

transformed parameters {
  // prior cumulative probability of each age category: effect on Pr(ever together)
  vector[length_dirichlet] dj_diff_together;
  dj_diff_together = append_row(0, delta_diff_together);
  
  // multimembership effects
  vector[num_nodes] node_mean_togt;
  node_mean_togt = mu_mm_togt + rand_mm_togt * tau_mm_togt;
  
  // global sigma
  real tau_sigma;
  tau_sigma = exp(tau_sigma_raw);      // Global scale is positive
  
  // node-level sigma
  vector[num_nodes] node_sigma_togt;
  node_sigma_togt = exp(raw_sigma_togt);         // Node-specific deviations are positive

  // regression equation for Pr(ever together)
  vector[num_dyads] predictor_together;
  for (i in 1:num_dyads) {
    predictor_together[i] = beta_diff_together * sum(dj_diff_together[1:age_diff[i]]) + mm_nodes_togt[node_1[i]] + mm_nodes_togt[node_2[i]];
  }
}

model {
  // age priors
  beta_diff_together ~ normal(0,2);
  
  // difference between age categories
  delta_diff_together ~ dirichlet(prior_diff);
  
  // variance
  tau_sigma_raw ~ normal(0,0.25);

  // multimembership priors
  mm_nodes_togt ~ normal(node_mean_togt, node_sigma_togt);
  mu_mm_togt ~ normal(logit(0.1),1);
  rand_mm_togt ~ normal(0,1);
  tau_mm_togt ~ normal(0,0.5);
  raw_sigma_togt ~ normal(0,0.5);            // Non-centered parameterization

  // likelihoods
  together ~ bernoulli_logit(predictor_together);
}
