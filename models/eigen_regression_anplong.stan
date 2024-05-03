data {
  // global data size
  int num_data;                          // Number of rows in total data
  int num_nodes;                         // Number of unique nodes across all windows
  int num_windows;                       // Number of time windows
  
  // per time window data size
  int num_nodes_window1;  // Number of rows in data in time window 1
  int num_nodes_window2;  // Number of rows in data in time window 2
  int num_nodes_window3;  // Number of rows in data in time window 3
  int num_nodes_window4;  // Number of rows in data in time window 4
  int num_nodes_window5;  // Number of rows in data in time window 5
  int num_nodes_window6;  // Number of rows in data in time window 6
  int num_nodes_window7;  // Number of rows in data in time window 7

  // number of nodes in all preceding time windows for node_age indexing
  array[num_windows] int num_nodes_prev_windows; // Number of rows in data in time window 2
  
  // centrality means per time window
  vector[num_nodes_window1] centrality_mu_1;      // Means of centrality estimates (0-1 bounded then standardised), time window 1
  vector[num_nodes_window2] centrality_mu_2;      // Means of centrality estimates (0-1 bounded then standardised), time window 2
  vector[num_nodes_window3] centrality_mu_3;      // Means of centrality estimates (0-1 bounded then standardised), time window 3
  vector[num_nodes_window4] centrality_mu_4;      // Means of centrality estimates (0-1 bounded then standardised), time window 4
  vector[num_nodes_window5] centrality_mu_5;      // Means of centrality estimates (0-1 bounded then standardised), time window 5
  vector[num_nodes_window6] centrality_mu_6;      // Means of centrality estimates (0-1 bounded then standardised), time window 6
  vector[num_nodes_window7] centrality_mu_7;      // Means of centrality estimates (0-1 bounded then standardised), time window 7
  
  // covariance matrices per time window
  matrix[num_nodes_window1, num_nodes_window1] centrality_cov_1;  // standard deviations of centrality estimates, time window 1
  matrix[num_nodes_window2, num_nodes_window2] centrality_cov_2;  // standard deviations of centrality estimates, time window 2
  matrix[num_nodes_window3, num_nodes_window3] centrality_cov_3;  // standard deviations of centrality estimates, time window 3
  matrix[num_nodes_window4, num_nodes_window4] centrality_cov_4;  // standard deviations of centrality estimates, time window 4
  matrix[num_nodes_window5, num_nodes_window5] centrality_cov_5;  // standard deviations of centrality estimates, time window 5
  matrix[num_nodes_window6, num_nodes_window6] centrality_cov_6;  // standard deviations of centrality estimates, time window 6
  matrix[num_nodes_window7, num_nodes_window7] centrality_cov_7;  // standard deviations of centrality estimates, time window 7
  
  // node IDs for all time windows
  array[num_nodes_window1] int nodes_window1;             // Node IDs, time window 1
  array[num_nodes_window2] int nodes_window2;             // Node IDs, time window 2
  array[num_nodes_window3] int nodes_window3;             // Node IDs, time window 3
  array[num_nodes_window4] int nodes_window4;             // Node IDs, time window 4
  array[num_nodes_window5] int nodes_window5;             // Node IDs, time window 5
  array[num_nodes_window6] int nodes_window6;             // Node IDs, time window 6
  array[num_nodes_window7] int nodes_window7;             // Node IDs, time window 7

  // node ages for all individuals
  array[num_data] real node_age;         // Node ages (years)
  
}

parameters {
  // intercept
  real intercept;
  // exposure slope
  real beta_age;
  // variance and degrees of freedom
  real<lower=0> sigma;
  // real<lower=1,upper=5> nu;
  // random effects
  vector[num_windows] rand_window;
  real<lower=0> sigma_window;
  vector[num_nodes] rand_node;
  real<lower=0> sigma_node;
}

transformed parameters {
  // random effects
  vector[num_windows] window_random_effect;
  window_random_effect = rand_window * sigma_window;
  vector[num_nodes] node_random_effect;
  node_random_effect = rand_node * sigma_node;
  
  // linear model
  vector[num_nodes_window1] predictor_window1;
  for(i in 1:num_nodes_window1) {
    predictor_window1[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[1])] + window_random_effect[1] + node_random_effect[nodes_window1[i]];
  }
  vector[num_nodes_window2] predictor_window2;
  for(i in 1:num_nodes_window2) {
    predictor_window2[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[2])] + window_random_effect[2] + node_random_effect[nodes_window2[i]];
  }
  vector[num_nodes_window3] predictor_window3;
  for(i in 1:num_nodes_window3) {
    predictor_window3[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[3])] + window_random_effect[3] + node_random_effect[nodes_window3[i]];
  }
  vector[num_nodes_window4] predictor_window4;
  for(i in 1:num_nodes_window4) {
    predictor_window4[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[4])] + window_random_effect[4] + node_random_effect[nodes_window4[i]];
  }
  vector[num_nodes_window5] predictor_window5;
  for(i in 1:num_nodes_window5) {
    predictor_window5[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[5])] + window_random_effect[5] + node_random_effect[nodes_window5[i]];
  }
  vector[num_nodes_window6] predictor_window6;
  for(i in 1:num_nodes_window6) {
    predictor_window6[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[6])] + window_random_effect[6] + node_random_effect[nodes_window6[i]];
  }
  vector[num_nodes_window7] predictor_window7;
  for(i in 1:num_nodes_window7) {
    predictor_window7[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[7])] + window_random_effect[7] + node_random_effect[nodes_window7[i]];
  }
}

model {
  // priors
  intercept ~ normal(logit(0.05),2); // 0.05 best estimate of intercept for Amboseli population eigen ~ age (edge weights from half weight index) from Chiyo 2011 -- estimated intercept from figure. Logit because on 0-1 scale 
  beta_age ~ normal(0,0.8);
  sigma ~ exponential(2);
  // nu ~ normal(3,1);
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
  centrality_mu_6 ~ multi_normal(predictor_window6, centrality_cov_6 + diag_matrix(rep_vector(sigma, num_nodes_window6)));
  centrality_mu_7 ~ multi_normal(predictor_window7, centrality_cov_7 + diag_matrix(rep_vector(sigma, num_nodes_window7)));
  // centrality_mu_1 ~ multi_student_t(nu, predictor_window1, centrality_cov_1 + diag_matrix(rep_vector(sigma, num_nodes_window1)));
  // centrality_mu_2 ~ multi_student_t(nu, predictor_window2, centrality_cov_2 + diag_matrix(rep_vector(sigma, num_nodes_window2)));
  // centrality_mu_3 ~ multi_student_t(nu, predictor_window3, centrality_cov_3 + diag_matrix(rep_vector(sigma, num_nodes_window3)));
  // centrality_mu_4 ~ multi_student_t(nu, predictor_window4, centrality_cov_4 + diag_matrix(rep_vector(sigma, num_nodes_window4)));
  // centrality_mu_5 ~ multi_student_t(nu, predictor_window5, centrality_cov_5 + diag_matrix(rep_vector(sigma, num_nodes_window5)));
  // centrality_mu_6 ~ multi_student_t(nu, predictor_window6, centrality_cov_6 + diag_matrix(rep_vector(sigma, num_nodes_window6)));
  // centrality_mu_7 ~ multi_student_t(nu, predictor_window7, centrality_cov_7 + diag_matrix(rep_vector(sigma, num_nodes_window7)));

}

