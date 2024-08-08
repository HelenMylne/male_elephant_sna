data {
  int num_dyads;                      // Number of edges
  int num_nodes;                      // Number of nodes
  int num_age_cat;                    // Number of unique age categories
  int length_dirichlet;               // Number of unique age categories + 1
  vector[num_dyads] logit_edge_mu;    // Means of logit edge weights
  vector[num_dyads] logit_edge_sd;    // Standard deviation of logit edge weights
  array[num_dyads] int age_min_cat;   // age of younger dyad member
  array[num_dyads] int node_1;        // Node 1 IDs for multimembership terms
  array[num_dyads] int node_2;        // Node 2 IDs for multimembership terms
  vector[num_age_cat] prior_min;      // Dirichlet prior values
}

parameters {
  // intercept
  real intercept;
  // exposure slopes
  real beta_age_min;
  // variance
  real<lower=0> sigma;
  // multimembership effects
  vector[num_nodes] mm_nodes;
  real<lower=0> sigma_mm;
  // difference between age categories
  simplex[num_age_cat] delta_min;
}

transformed parameters {
  // create prior for cumulative probability of each age category
  vector[length_dirichlet] delta_j_min;
  delta_j_min = append_row(0, delta_min);
  
  // regression equation
  vector[num_dyads] predictor;
  for (i in 1:num_dyads) {
    predictor[i] = intercept + beta_age_min * sum(delta_j_min[1:age_min_cat[i]]) + mm_nodes[node_1[i]] + mm_nodes[node_2[i]];
  }
}

model {
  // priors
  intercept ~ normal(0,2);
  beta_age_min ~ normal(0,2);
  sigma ~ exponential(1);
  mm_nodes ~ normal(0, sigma_mm);
  sigma_mm ~ exponential(2);
  delta_min ~ dirichlet(prior_min);

  // likelihood
  logit_edge_mu ~ normal(predictor, logit_edge_sd + rep_vector(sigma, num_dyads));
  
}
