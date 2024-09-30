data {
  // global data size
  int num_dyads;                       // Number of edges
  int num_nodes;                       // Number of nodes
  int num_age_diffs;                   // Number of unique age difference categories
  int length_dirichlet;                // Number of unique age difference categories + 1
  // int num_samples;
  
  // Gaussian approximation of edge weights
  vector[num_dyads] logit_edge_mu;     // Means of logit edge weights
  vector[num_dyads] logit_edge_sd;     // Standard deviation of logit edge weights
  
  // explanatory variables
  array[num_dyads] int age_diff;       // ordinal categorical where 1 = same age, 2 = 1 category different, 3 = 2 categories different... etc.
  // array[num_samples,num_dyads] int together_imputed; // predicted together vs not together from previous model
  vector[num_dyads] pr_together;       // mean predicted probability together vs not together from previous model
  
  // multimembership terms
  array[num_dyads] int node_1;         // Node 1 IDs for multimembership terms
  array[num_dyads] int node_2;         // Node 2 IDs for multimembership terms
  
  // prior values for Dirichlet
  vector[num_age_diffs] prior_diff;    // Dirichlet prior values
}

parameters {
  // exposure slopes
  // vector[2] beta_together_weight;      // effect of together 1/0 on overall edge weight
  real beta_together_weight;           // effect of Pr(ever together) on overall edge weight
  real beta_diff_weight;               // effect of age on overall edge weight
  
  // global variance
  real tau_sigma_raw_weight;           // Unconstrained real value for global scale
  
  // multimembership effects
  vector[num_nodes] mm_nodes_weight;
  real mu_mm_weight;
  vector[num_nodes] rand_mm_weight;
  real<lower=0> tau_mm_weight;
  vector[num_nodes] raw_sigma_weight;
  
  // difference between age categories
  simplex[num_age_diffs] delta_diff_weight;
}

transformed parameters {
  // prior cumulative probability of each age category: effect on edge weight
  vector[length_dirichlet] dj_diff_weight;
  dj_diff_weight = append_row(0, delta_diff_weight);
  
  // multimembership effects
  vector[num_nodes] node_mean_weight;
  node_mean_weight = mu_mm_weight + rand_mm_weight * tau_mm_weight;
  
  // global sigma
  real tau_sigma_weight;
  tau_sigma_weight = exp(tau_sigma_raw_weight);      // Global scale is positive
  
  // node-level sigma
  vector[num_nodes] node_sigma_weight;
  node_sigma_weight = exp(raw_sigma_weight);         // Node-specific deviations are positive

  // regression equation for Pr(ever together)
  vector[num_dyads] predictor_edge;
  for (i in 1:num_dyads) {
    predictor_edge[i] = beta_diff_weight * sum(dj_diff_weight[1:age_diff[i]]) + mm_nodes_weight[node_1[i]] + mm_nodes_weight[node_2[i]] + beta_together_weight * pr_together[i]; // + beta_together_weight[together_imputed[i]]
  
  }
}

model {
  // age priors
  beta_together_weight ~ normal(1,3);
  beta_diff_weight ~ normal(0,1);
  
  // difference between age categories
  delta_diff_weight ~ dirichlet(prior_diff);
  
  // variance
  tau_sigma_raw_weight ~ normal(0,1);

  // multimembership priors
  mm_nodes_weight ~ normal(node_mean_weight, node_sigma_weight);
  mu_mm_weight ~ normal(logit(0.1),1);
  rand_mm_weight ~ normal(0,1);
  tau_mm_weight ~ normal(0,1);
  raw_sigma_weight ~ normal(0,0.5);            // Non-centered parameterization

  // likelihoods
  logit_edge_mu ~ normal(predictor_edge, logit_edge_sd + rep_vector(tau_sigma_weight, num_dyads));
}

