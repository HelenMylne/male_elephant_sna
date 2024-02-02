data {
  int num_data;                          // Number of rows in total data
  int num_nodes;                         // Number of unique nodes across all windows
  //int num_windows;                       // Number of time windows
  int num_nodes_window1; // Number of rows in data in time window 1
  vector[num_nodes_window1] centrality_mu_1;      // Means of centrality estimates (0-1 bounded then standardised), time window 1
  matrix[num_nodes_window1, num_nodes_window1] centrality_cov_1;  // standard deviations of centrality estimates, time window 1
  array[num_data] real node_age;         // Node ages (point estimate)
  //array[num_data] int window;            // ID of time window (random effect)
  //array[num_data] int nodes;             // Node IDs
  array[num_nodes_window1] int nodes_window1;             // Node IDs
}

parameters {
  // intercept
  real intercept;
  // exposure slope
  real beta_age;
  //variance
  real<lower=0> sigma;
  // random effects
  //vector[num_windows] rand_window;
  //vector[num_nodes] rand_node;
}

transformed parameters {
  // linear model
  vector[num_nodes_window1] predictor_window1;
  for(i in 1:num_nodes_window1) {
    predictor_window1[i] = intercept + beta_age * node_age[i];// + rand_window[1] + rand_node[nodes_window1[i]];
  }
}

 model {
  // priors
  beta_age ~ normal(-1, 0.1);    // beta_age ~ normal(0, 0.8);
  sigma ~ exponential(2);
  intercept ~ normal(2,0.1);     // intercept ~ normal(0,0.8);
  rand_window ~ normal(0,1);     // rand_window ~ normal(0,1);
  rand_node ~ normal(0,0.2);     // rand_node ~ normal(0,1);

  // likelihood
  centrality_mu_1 ~ multi_normal(predictor_window1, centrality_cov_1 + diag_matrix(rep_vector(sigma, num_nodes_window1)));
}

