data {
  int num_data;                     // Number of rows in total data
  int num_nodes;                    // Number of unique nodes across all windows
  int num_windows;                  // Number of time windows
  vector[num_data] centrality_mu;  // Means of centrality estimates (0-1 bounded then standardised)
  matrix[num_data, num_data] centrality_cov;  // standard deviations of centrality estimates
  array[num_data] real node_age;       // Node ages (point estimate)
  array[num_data] int window;                  // ID of time window (random effect)
  array[num_data] int nodes;          // Node IDs
}

parameters {
  // intercept
  real intercept;
  // exposure slope
  real beta_age;
  //variance
  real<lower=0> sigma;
  // random effects
  vector[num_windows] rand_window;
  vector[num_dyads] rand_dyad;
}

transformed parameters {
  // linear model
  vector[num_data] predictor;
  for(i in 1:num_data) {
    predictor[i] = intercept + beta_age * node_age[i] + rand_window[window[i]] + rand_dyad[nodes[i]];
  }
}

 model {
  // priors
  beta_age ~ normal(0, 0.8);
  sigma ~ exponential(2);
  intercept ~ normal(0,0.8);
  rand_window ~ normal(0,1);
  rand_dyad ~ normal(0,1);

  // likelihood
  centrality_mu ~ multi_normal(predictor, centrality_cov + diag_matrix(rep_vector(sigma, num_nodes)));
}
