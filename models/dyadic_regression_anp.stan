data {
  // global data size
  int num_data;                                   // Number of edges in total data
  int num_dyads;                                  // Number of unique dyads across all windows
  int num_nodes;                                  // Number of unique nodes across all windows
  int num_windows;                                // Number of time windows
  
  // per time window data size
  int num_dyads_window1;                          // Number of rows in data in time window 1
  int num_dyads_window2;                          // Number of rows in data in time window 2
  int num_dyads_window3;                          // Number of rows in data in time window 3
  // int num_dyads_window4;                          // Number of rows in data in time window 4
  // int num_dyads_window5;                          // Number of rows in data in time window 5
  // int num_dyads_window6;                          // Number of rows in data in time window 6
  // int num_dyads_window7;                          // Number of rows in data in time window 7
  
  // number of dyads in all preceding time windows for age indexing
  array[num_windows] int num_dyads_prev_windows;  // Number of rows in data in time window -1
  
  // edge means per time window
  //vector[num_data] logit_edge_mu;                 // Means of Gaussian approximation of logit edge weights
  vector[num_dyads_window1] logit_edge_mu_1;      // Mean normal approximation logit edge weights, window 1
  vector[num_dyads_window2] logit_edge_mu_2;      // Mean normal approximation logit edge weights, window 2
  vector[num_dyads_window3] logit_edge_mu_3;      // Mean normal approximation logit edge weights, window 3
  // vector[num_dyads_window4] logit_edge_mu_4;      // Mean normal approximation logit edge weights, window 4
  // vector[num_dyads_window5] logit_edge_mu_5;      // Mean normal approximation logit edge weights, window 5
  // vector[num_dyads_window6] logit_edge_mu_6;      // Mean normal approximation logit edge weights, window 6
  // vector[num_dyads_window7] logit_edge_mu_7;      // Mean normal approximation logit edge weights, window 7
  
  // covariance matrices per time window
  //matrix[num_data, num_data] logit_edge_cov;     // Covariance matrix of Gaussian approximation of logit edge weights
  matrix[num_dyads_window1, num_dyads_window1] logit_edge_cov_1;  // SD logit edge weights, window 1
  matrix[num_dyads_window2, num_dyads_window2] logit_edge_cov_2;  // SD logit edge weights, window 2
  matrix[num_dyads_window3, num_dyads_window3] logit_edge_cov_3;  // SD logit edge weights, window 3
  // matrix[num_dyads_window4, num_dyads_window4] logit_edge_cov_4;  // SD logit edge weights, window 4
  // matrix[num_dyads_window5, num_dyads_window5] logit_edge_cov_5;  // SD logit edge weights, window 5
  // matrix[num_dyads_window6, num_dyads_window6] logit_edge_cov_6;  // SD logit edge weights, window 6
  // matrix[num_dyads_window7, num_dyads_window7] logit_edge_cov_7;  // SD logit edge weights, window 7
  
  // explanatory variables
  array[num_data] real age_min;                // age of younger dyad member
  array[num_data] real age_max;                // age of older dyad member
  
  // multimembership terms
  array[num_data] int node_1;                  // Node 1 IDs for multimembership terms
  array[num_data] int node_2;                  // Node 2 IDs for multimembership terms
  
  // dyad IDs for all time windows
  array[num_dyads_window1] int dyads_window1;             // Dyad IDs, window 1
  array[num_dyads_window2] int dyads_window2;             // Dyad IDs, window 2
  array[num_dyads_window3] int dyads_window3;             // Dyad IDs, window 3
  // array[num_dyads_window4] int dyads_window4;             // Dyad IDs, window 4
  // array[num_dyads_window5] int dyads_window5;             // Dyad IDs, window 5
  // array[num_dyads_window6] int dyads_window6;             // Dyad IDs, window 6
  // array[num_dyads_window7] int dyads_window7;             // Dyad IDs, window 7

}

parameters {
  // intercept
  real intercept;
  // exposure slopes
  real beta_age_max;
  real beta_age_min;
  // Cholesky variance
  real<lower=0> sigma;
  // multimembership effects
  vector[num_nodes] mm_nodes;
  real<lower=0> sigma_mm;
  // random effect of time window
  vector[num_windows] rand_window;
  real sigma_window;
  // random effect of dyad ID
  vector[num_dyads] rand_dyad;
  real sigma_dyad;
}

transformed parameters {
  // //  Cholesky factor of the covariance matrix
  // matrix[num_data, num_data] L_cov;
  // L_cov = cholesky_decompose(logit_edge_cov + diag_matrix(rep_vector(sigma, num_data)));
  
  // Cholesky factor of each covariance matrix
  matrix[num_dyads_window1, num_dyads_window1] L_cov_1;
  L_cov_1 = cholesky_decompose(logit_edge_cov_1 + diag_matrix(rep_vector(sigma, num_dyads_window1)));
  matrix[num_dyads_window2, num_dyads_window2] L_cov_2;
  L_cov_2 = cholesky_decompose(logit_edge_cov_2 + diag_matrix(rep_vector(sigma, num_dyads_window2)));
  matrix[num_dyads_window3, num_dyads_window3] L_cov_3;
  L_cov_3 = cholesky_decompose(logit_edge_cov_3 + diag_matrix(rep_vector(sigma, num_dyads_window3)));
  // matrix[num_dyads_window4, num_dyads_window4] L_cov_4;
  // L_cov_4 = cholesky_decompose(logit_edge_cov_4 + diag_matrix(rep_vector(sigma, num_dyads_window4)));
  // matrix[num_dyads_window5, num_dyads_window5] L_cov_5;
  // L_cov_5 = cholesky_decompose(logit_edge_cov_5 + diag_matrix(rep_vector(sigma, num_dyads_window5)));
  // matrix[num_dyads_window6, num_dyads_window6] L_cov_6;
  // L_cov_6 = cholesky_decompose(logit_edge_cov_6 + diag_matrix(rep_vector(sigma, num_dyads_window6)));
  // matrix[num_dyads_window7, num_dyads_window7] L_cov_7;
  // L_cov_7 = cholesky_decompose(logit_edge_cov_7 + diag_matrix(rep_vector(sigma, num_dyads_window7)));
  
  // random effects
  vector[num_windows] window_random_effect;
  window_random_effect = rand_window * sigma_window;
  vector[num_dyads] dyad_random_effect;
  dyad_random_effect = rand_dyad * sigma_dyad;
  vector[num_nodes] mm_nodes;
  mm_nodes = rand_mm * sigma_mm;
  
  // // regression equation
  // vector[num_data] predictor;
  // for (i in 1:num_data) {
  //   predictor[i] = intercept + beta_age_min * age_min[i] + beta_age_max * age_max[i] + mm_nodes[node_1[i]] + mm_nodes[node_2[i]] + rand_window[window[i]] + rand_dyad[dyad_id[i]];
  // }

  // linear model
  vector[num_dyads_window1] predictor_window1;
  for(i in 1:num_dyads_window1) {
    predictor_window1[i] = intercept + beta_age_min * age_min[(i + num_dyads_prev_windows[1])] + beta_age_max * age_max[(i + num_dyads_prev_windows[1])] + mm_nodes[node_1[(i + num_dyads_prev_windows[1])]] + mm_nodes[node_2[(i + num_dyads_prev_windows[1])]] + window_random_effect[1] + dyad_random_effect[dyads_window1[i]];
  }
  vector[num_dyads_window2] predictor_window2;
  for(i in 1:num_dyads_window2) {
    predictor_window2[i] = intercept + beta_age_min * age_min[(i + num_dyads_prev_windows[2])] + beta_age_max * age_max[(i + num_dyads_prev_windows[2])] + mm_nodes[node_1[(i + num_dyads_prev_windows[2])]] + mm_nodes[node_2[(i + num_dyads_prev_windows[2])]] + window_random_effect[2] + dyad_random_effect[dyads_window2[i]];
  }
  vector[num_dyads_window3] predictor_window3;
  for(i in 1:num_dyads_window3) {
    predictor_window3[i] = intercept + beta_age_min * age_min[(i + num_dyads_prev_windows[3])] + beta_age_max * age_max[(i + num_dyads_prev_windows[3])] + mm_nodes[node_1[(i + num_dyads_prev_windows[3])]] + mm_nodes[node_2[(i + num_dyads_prev_windows[3])]] + window_random_effect[3] + dyad_random_effect[dyads_window3[i]];
  }
  // vector[num_dyads_window4] predictor_window4;
  // for(i in 1:num_dyads_window4) {
  //   predictor_window4[i] = intercept + beta_age_min * age_min[(i + num_dyads_prev_windows[4])] + beta_age_max * age_max[(i + num_dyads_prev_windows[4])] + mm_nodes[node_1[(i + num_dyads_prev_windows[4])]] + mm_nodes[node_2[(i + num_dyads_prev_windows[4])]] + window_random_effect[4] + dyad_random_effect[dyads_window4[i]];
  // }
  // vector[num_dyads_window5] predictor_window5;
  // for(i in 1:num_dyads_window5) {
  //   predictor_window5[i] = intercept + beta_age_min * age_min[(i + num_dyads_prev_windows[5])] + beta_age_max * age_max[(i + num_dyads_prev_windows[5])] + mm_nodes[node_1[(i + num_dyads_prev_windows[5])]] + mm_nodes[node_2[(i + num_dyads_prev_windows[5])]] + window_random_effect[5] + dyad_random_effect[dyads_window5[i]];
  // }
  // vector[num_dyads_window6] predictor_window6;
  // for(i in 1:num_dyads_window6) {
  //   predictor_window6[i] = intercept + beta_age_min * age_min[(i + num_dyads_prev_windows[6])] + beta_age_max * age_max[(i + num_dyads_prev_windows[6])] + mm_nodes[node_1[(i + num_dyads_prev_windows[6])]] + mm_nodes[node_2[(i + num_dyads_prev_windows[6])]] + window_random_effect[6] + dyad_random_effect[dyads_window6[i]];
  // }
  // vector[num_dyads_window7] predictor_window7;
  // for(i in 1:num_dyads_window7) {
  //   predictor_window7[i] = intercept + beta_age_min * age_min[(i + num_dyads_prev_windows[7])] + beta_age_max * age_max[(i + num_dyads_prev_windows[7])] + mm_nodes[node_1[(i + num_dyads_prev_windows[7])]] + mm_nodes[node_2[(i + num_dyads_prev_windows[7])]] + window_random_effect[7] + dyad_random_effect[dyads_window7[i]];
  // }
  
}

model {
  // intercept prior
  intercept ~ normal(0,1);
  // slope priors
  beta_age_max ~ normal(0,1);
  beta_age_min ~ normal(0,1);
  // Cholesky variance
  sigma ~ exponential(2);
  // multimembership priors
  rand_mm ~ normal(0,1);
  sigma_mm ~ exponential(2);
  // random effect of window 
  rand_window ~ normal(0,1);
  sigma_window ~ exponential(2);
  // random effect of dyad ID
  rand_window ~ normal(0,1);
  sigma_dyad ~ exponential(2);
  
  // likelihood
  logit_edge_mu_1 ~ multi_normal(predictor_window1, L_cov_1);
  logit_edge_mu_2 ~ multi_normal(predictor_window2, L_cov_2);
  logit_edge_mu_3 ~ multi_normal(predictor_window3, L_cov_3);
  // logit_edge_mu_4 ~ multi_normal(predictor_window4, L_cov_4);
  // logit_edge_mu_5 ~ multi_normal(predictor_window5, L_cov_5);
  // logit_edge_mu_6 ~ multi_normal(predictor_window6, L_cov_6);
  // logit_edge_mu_7 ~ multi_normal(predictor_window7, L_cov_7);
}
