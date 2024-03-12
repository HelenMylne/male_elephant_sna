data {
  // global data size
  int num_data;                                   // Number of rows in total data
  int num_nodes;                                  // Number of unique nodes across all windows
  int num_windows;                                // Number of time windows
  
  // age categories
  int num_age_cat;                                // Number of unique age categories
  int length_dirichlet;                           // Number of unique age categories + 1
  vector[num_age_cat] prior_age;                  // Dirichlet prior values

  // per time window data size
  int num_nodes_window1;                          // Number of rows in data in time window 1
  int num_nodes_window2;                          // Number of rows in data in time window 2
  int num_nodes_window3;                          // Number of rows in data in time window 3
  int num_nodes_window4;                          // Number of rows in data in time window 4
  int num_nodes_window5;                          // Number of rows in data in time window 5
  
  // number of nodes in all preceding time windows for node_age indexing
  array[num_windows] int num_nodes_prev_windows;  // Number of rows in data in time window 2
  
  // centrality means per time window
  vector[num_nodes_window1] centrality_mu_1;      // Means of centrality estimates (0-1 bounded then standardised), time window 1
  vector[num_nodes_window2] centrality_mu_2;      // Means of centrality estimates (0-1 bounded then standardised), time window 2
  vector[num_nodes_window3] centrality_mu_3;      // Means of centrality estimates (0-1 bounded then standardised), time window 3
  vector[num_nodes_window4] centrality_mu_4;      // Means of centrality estimates (0-1 bounded then standardised), time window 4
  vector[num_nodes_window5] centrality_mu_5;      // Means of centrality estimates (0-1 bounded then standardised), time window 5
  
  // covariance matrices per time window
  matrix[num_nodes_window1, num_nodes_window1] centrality_cov_1;  // standard deviations of centrality estimates, time window 1
  matrix[num_nodes_window2, num_nodes_window2] centrality_cov_2;  // standard deviations of centrality estimates, time window 2
  matrix[num_nodes_window3, num_nodes_window3] centrality_cov_3;  // standard deviations of centrality estimates, time window 3
  matrix[num_nodes_window4, num_nodes_window4] centrality_cov_4;  // standard deviations of centrality estimates, time window 4
  matrix[num_nodes_window5, num_nodes_window5] centrality_cov_5;  // standard deviations of centrality estimates, time window 5

  // node IDs for all time windows
  array[num_nodes_window1] int nodes_window1;     // Node IDs, time window 1
  array[num_nodes_window2] int nodes_window2;     // Node IDs, time window 2
  array[num_nodes_window3] int nodes_window3;     // Node IDs, time window 3
  array[num_nodes_window4] int nodes_window4;     // Node IDs, time window 4
  array[num_nodes_window5] int nodes_window5;     // Node IDs, time window 5

  // node ages for all individuals
  array[num_data] int node_age;         // Node ages (categorical)
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
  // difference between age categories
  simplex[num_age_cat] delta;
}

transformed parameters {
  // create prior for cumulative probability of each age category
  vector[length_dirichlet] delta_j;
  delta_j = append_row(0, delta);
  
  // random effects
  vector[num_windows] window_random_effect;
  window_random_effect = rand_window * sigma_window;
  vector[num_nodes] node_random_effect;
  node_random_effect = rand_node * sigma_node;
  
  // linear model
  vector[num_nodes_window1] predictor_window1;
  for(i in 1:num_nodes_window1) {
    predictor_window1[i] = intercept + beta_age * sum(delta_j[1:node_age[(i + num_nodes_prev_windows[1])]]) + rand_window[1] + rand_node[nodes_window1[i]];
  }
  vector[num_nodes_window2] predictor_window2;
  for(i in 1:num_nodes_window2) {
    predictor_window2[i] = intercept + beta_age * sum(delta_j[1:node_age[(i + num_nodes_prev_windows[2])]]) + rand_window[2] + rand_node[nodes_window2[i]];
  }
  vector[num_nodes_window3] predictor_window3;
  for(i in 1:num_nodes_window3) {
    predictor_window3[i] = intercept + beta_age * sum(delta_j[1:node_age[(i + num_nodes_prev_windows[3])]]) + rand_window[3] + rand_node[nodes_window3[i]];
  }
  vector[num_nodes_window4] predictor_window4;
  for(i in 1:num_nodes_window4) {
    predictor_window4[i] = intercept + beta_age * sum(delta_j[1:node_age[(i + num_nodes_prev_windows[4])]]) + rand_window[4] + rand_node[nodes_window4[i]];
  }
  vector[num_nodes_window5] predictor_window5;
  for(i in 1:num_nodes_window5) {
    predictor_window5[i] = intercept + beta_age * sum(delta_j[1:node_age[(i + num_nodes_prev_windows[5])]]) + rand_window[5] + rand_node[nodes_window5[i]];
  }
}

 model {
  // priors
  intercept ~ normal(logit(0.05),2); //0.05 best estimate of intercept for Amboseli population eigen ~ age (edge weights from half weight index) from Chiyo 2011 -- estimated intercept from figure. Logit because on 0-1 scale 
  delta ~ dirichlet(prior_age);
  beta_age  ~ normal(0,0.8);
  sigma ~ exponential(2);
  // random effects
  rand_window ~ normal(0,1); // rand_window ~ normal(0,sigma_window);
  sigma_window ~ exponential(2);
  rand_node ~ normal(0,1);   // rand_node ~ normal(0,sigma_node);
  sigma_node ~ exponential(2);
  
  // likelihood
  centrality_mu_1 ~ multi_normal(predictor_window1, centrality_cov_1 + diag_matrix(rep_vector(sigma, num_nodes_window1)));
  centrality_mu_2 ~ multi_normal(predictor_window2, centrality_cov_2 + diag_matrix(rep_vector(sigma, num_nodes_window2)));
  centrality_mu_3 ~ multi_normal(predictor_window3, centrality_cov_3 + diag_matrix(rep_vector(sigma, num_nodes_window3)));
  centrality_mu_4 ~ multi_normal(predictor_window4, centrality_cov_4 + diag_matrix(rep_vector(sigma, num_nodes_window4)));
  centrality_mu_5 ~ multi_normal(predictor_window5, centrality_cov_5 + diag_matrix(rep_vector(sigma, num_nodes_window5)));
}

