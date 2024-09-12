data {
  // global data size
  int num_dyads;                      // Number of edges
  int num_nodes;                      // Number of nodes
  int num_age_cat;                    // Number of unique age categories
  int length_dirichlet;               // Number of unique age categories + 1
  
  // Gaussian approximation of edge weights
  vector[num_dyads] logit_edge_mu;    // Means of logit edge weights
  vector[num_dyads] logit_edge_sd;    // Standard deviation of logit edge weights
  
  // explanatory variable
  array[num_dyads] int age_min_cat;   // age of younger dyad member
  
  // multimembership terms
  array[num_dyads] int node_1;        // Node 1 IDs for multimembership terms
  array[num_dyads] int node_2;        // Node 2 IDs for multimembership terms
  
  // prior values for Dirichlet
  vector[num_age_cat] prior_min;      // Dirichlet prior values
}

parameters {
  // intercept
  real intercept;
  
  // exposure slope
  real beta_age_min;
  
  // global variance
  real tau_sigma_raw;            // Unconstrained real value for global scale
  // real<lower=0> sigma;
  
  // multimembership effects
  vector[num_nodes] mm_nodes;
  real mu_mm;
  vector[num_nodes] rand_mm;
  real<lower=0> tau_mm;
  // vector<lower=0>[num_nodes] node_sigma;
  // real theta_node;
  vector[num_nodes] raw_sigma;   // Unconstrained real values for node-specific scales
  
  // difference between age categories
  simplex[num_age_cat] delta_min;
}

transformed parameters {
  // create prior for cumulative probability of each age category
  vector[length_dirichlet] delta_j_min;
  delta_j_min = append_row(0, delta_min);
  
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

  // regression equation
  vector[num_dyads] predictor;
  for (i in 1:num_dyads) {
    predictor[i] = intercept + beta_age_min * sum(delta_j_min[1:age_min_cat[i]]) + mm_nodes[node_1[i]] + mm_nodes[node_2[i]];
  }
}

model {
  // intercept prior
  intercept ~ normal(0,2);
  
  // age priors
  beta_age_min ~ normal(0,5);
  delta_min ~ dirichlet(prior_min);
  
  // variance
  tau_sigma_raw ~ normal(0,1);        //sigma ~ exponential(2); --> Prior on global scale before transformation

  // multimembership priors
  mm_nodes ~ normal(node_mean, node_sigma);
  mu_mm ~ normal(0,0.5);
  rand_mm ~ normal(0,1);
  tau_mm ~ cauchy(0,0.5);
  // node_sigma ~ exponential(theta_node);
  // theta_node ~ normal(0,2);
  raw_sigma ~ normal(0,0.5);            // Non-centered parameterization

  // likelihood
  logit_edge_mu ~ normal(predictor, logit_edge_sd + rep_vector(tau_sigma, num_dyads));//sigma, num_dyads));
}
