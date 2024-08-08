data {
  int num_data;                          // Number of rows in total data
  int num_nodes;                         // Number of unique nodes across all windows
  int num_windows;                       // Number of time windows
  int num_nodes_window1; // Number of rows in data in time window 1
  int num_nodes_window2; // Number of rows in data in time window 2
  int num_nodes_window3; // Number of rows in data in time window 2
  array[num_windows] int num_nodes_prev_windows; // Number of rows in data in time window 2
  vector[num_nodes_window1] centrality_mu_1;      // Means of centrality estimates (0-1 bounded then standardised), time window 1
  vector[num_nodes_window2] centrality_mu_2;      // Means of centrality estimates (0-1 bounded then standardised), time window 2
  vector[num_nodes_window3] centrality_mu_3;      // Means of centrality estimates (0-1 bounded then standardised), time window 3
  matrix[num_nodes_window1, num_nodes_window1] centrality_cov_1;  // standard deviations of centrality estimates, time window 1
  matrix[num_nodes_window2, num_nodes_window2] centrality_cov_2;  // standard deviations of centrality estimates, time window 2
  matrix[num_nodes_window3, num_nodes_window3] centrality_cov_3;  // standard deviations of centrality estimates, time window 3
  array[num_data] real node_age;         // Node ages (point estimate)
  array[num_nodes_window1] int nodes_window1;             // Node IDs
  array[num_nodes_window2] int nodes_window2;             // Node IDs
  array[num_nodes_window3] int nodes_window3;             // Node IDs
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
  vector[num_nodes] rand_node;
}

transformed parameters {
  // linear model
  vector[num_nodes_window1] predictor_window1;
  for(i in 1:num_nodes_window1) {
    predictor_window1[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[1])] + rand_window[1] + rand_node[nodes_window1[i]];
  }
  vector[num_nodes_window2] predictor_window2;
  for(i in 1:num_nodes_window2) {
    predictor_window2[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[2])] + rand_window[2] + rand_node[nodes_window2[i]];
  }
  vector[num_nodes_window3] predictor_window3;
  for(i in 1:num_nodes_window3) {
    predictor_window3[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[3])] + rand_window[3] + rand_node[nodes_window3[i]];
  }
}

 model {
  // priors
  beta_age ~ normal(0, 0.8);
  sigma ~ exponential(2);
  intercept ~ normal(-2,1); //intercept ~ normal(0,0.8);
  rand_window ~ normal(0,sigma_window);
  sigma_window ~ exponential(1);
  rand_node ~ normal(0,sigma_node);
  sigma_node ~ exponential(1);
  
  // likelihood
  centrality_mu_1 ~ multi_normal(predictor_window1, centrality_cov_1 + diag_matrix(rep_vector(sigma, num_nodes_window1)));
  centrality_mu_2 ~ multi_normal(predictor_window2, centrality_cov_2 + diag_matrix(rep_vector(sigma, num_nodes_window2)));
  centrality_mu_3 ~ multi_normal(predictor_window3, centrality_cov_3 + diag_matrix(rep_vector(sigma, num_nodes_window3)));
}

