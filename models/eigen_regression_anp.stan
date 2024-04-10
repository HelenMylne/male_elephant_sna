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
  int num_nodes_window8;  // Number of rows in data in time window 8
  int num_nodes_window9;  // Number of rows in data in time window 9
  int num_nodes_window10; // Number of rows in data in time window 10
  int num_nodes_window11;  // Number of rows in data in time window 11
  int num_nodes_window12;  // Number of rows in data in time window 12
  int num_nodes_window13;  // Number of rows in data in time window 13
  int num_nodes_window14;  // Number of rows in data in time window 14
  int num_nodes_window15;  // Number of rows in data in time window 15
  int num_nodes_window16;  // Number of rows in data in time window 16
  int num_nodes_window17;  // Number of rows in data in time window 17
  int num_nodes_window18;  // Number of rows in data in time window 18
  int num_nodes_window19;  // Number of rows in data in time window 19
  int num_nodes_window20; // Number of rows in data in time window 10
  int num_nodes_window21;  // Number of rows in data in time window 21
  int num_nodes_window22;  // Number of rows in data in time window 22
  int num_nodes_window23;  // Number of rows in data in time window 23
  int num_nodes_window24;  // Number of rows in data in time window 24
  int num_nodes_window25;  // Number of rows in data in time window 25
  int num_nodes_window26;  // Number of rows in data in time window 26
  int num_nodes_window27;  // Number of rows in data in time window 27
  int num_nodes_window28;  // Number of rows in data in time window 28
  int num_nodes_window29;  // Number of rows in data in time window 29
  int num_nodes_window30; // Number of rows in data in time window 30
  int num_nodes_window31;  // Number of rows in data in time window 31
  int num_nodes_window32;  // Number of rows in data in time window 32
  int num_nodes_window33;  // Number of rows in data in time window 33
  int num_nodes_window34;  // Number of rows in data in time window 34
  int num_nodes_window35;  // Number of rows in data in time window 35
  int num_nodes_window36;  // Number of rows in data in time window 36
  
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
  vector[num_nodes_window8] centrality_mu_8;      // Means of centrality estimates (0-1 bounded then standardised), time window 8
  vector[num_nodes_window9] centrality_mu_9;      // Means of centrality estimates (0-1 bounded then standardised), time window 9
  vector[num_nodes_window10] centrality_mu_10;      // Means of centrality estimates (0-1 bounded then standardised), time window 10
  vector[num_nodes_window11] centrality_mu_11;      // Means of centrality estimates (0-1 bounded then standardised), time window 11
  vector[num_nodes_window12] centrality_mu_12;      // Means of centrality estimates (0-1 bounded then standardised), time window 12
  vector[num_nodes_window13] centrality_mu_13;      // Means of centrality estimates (0-1 bounded then standardised), time window 13
  vector[num_nodes_window14] centrality_mu_14;      // Means of centrality estimates (0-1 bounded then standardised), time window 14
  vector[num_nodes_window15] centrality_mu_15;      // Means of centrality estimates (0-1 bounded then standardised), time window 15
  vector[num_nodes_window16] centrality_mu_16;      // Means of centrality estimates (0-1 bounded then standardised), time window 16
  vector[num_nodes_window17] centrality_mu_17;      // Means of centrality estimates (0-1 bounded then standardised), time window 17
  vector[num_nodes_window18] centrality_mu_18;      // Means of centrality estimates (0-1 bounded then standardised), time window 18
  vector[num_nodes_window19] centrality_mu_19;      // Means of centrality estimates (0-1 bounded then standardised), time window 19
  vector[num_nodes_window20] centrality_mu_20;      // Means of centrality estimates (0-1 bounded then standardised), time window 20
  vector[num_nodes_window21] centrality_mu_21;      // Means of centrality estimates (0-1 bounded then standardised), time window 21
  vector[num_nodes_window22] centrality_mu_22;      // Means of centrality estimates (0-1 bounded then standardised), time window 22
  vector[num_nodes_window23] centrality_mu_23;      // Means of centrality estimates (0-1 bounded then standardised), time window 23
  vector[num_nodes_window24] centrality_mu_24;      // Means of centrality estimates (0-1 bounded then standardised), time window 24
  vector[num_nodes_window25] centrality_mu_25;      // Means of centrality estimates (0-1 bounded then standardised), time window 25
  vector[num_nodes_window26] centrality_mu_26;      // Means of centrality estimates (0-1 bounded then standardised), time window 26
  vector[num_nodes_window27] centrality_mu_27;      // Means of centrality estimates (0-1 bounded then standardised), time window 27
  vector[num_nodes_window28] centrality_mu_28;      // Means of centrality estimates (0-1 bounded then standardised), time window 28
  vector[num_nodes_window29] centrality_mu_29;      // Means of centrality estimates (0-1 bounded then standardised), time window 29
  vector[num_nodes_window30] centrality_mu_30;      // Means of centrality estimates (0-1 bounded then standardised), time window 30
  vector[num_nodes_window31] centrality_mu_31;      // Means of centrality estimates (0-1 bounded then standardised), time window 31
  vector[num_nodes_window32] centrality_mu_32;      // Means of centrality estimates (0-1 bounded then standardised), time window 32
  vector[num_nodes_window33] centrality_mu_33;      // Means of centrality estimates (0-1 bounded then standardised), time window 33
  vector[num_nodes_window34] centrality_mu_34;      // Means of centrality estimates (0-1 bounded then standardised), time window 34
  vector[num_nodes_window35] centrality_mu_35;      // Means of centrality estimates (0-1 bounded then standardised), time window 35
  vector[num_nodes_window36] centrality_mu_36;      // Means of centrality estimates (0-1 bounded then standardised), time window 36
  
  // covariance matrices per time window
  matrix[num_nodes_window1, num_nodes_window1] centrality_cov_1;  // standard deviations of centrality estimates, time window 1
  matrix[num_nodes_window2, num_nodes_window2] centrality_cov_2;  // standard deviations of centrality estimates, time window 2
  matrix[num_nodes_window3, num_nodes_window3] centrality_cov_3;  // standard deviations of centrality estimates, time window 3
  matrix[num_nodes_window4, num_nodes_window4] centrality_cov_4;  // standard deviations of centrality estimates, time window 4
  matrix[num_nodes_window5, num_nodes_window5] centrality_cov_5;  // standard deviations of centrality estimates, time window 5
  matrix[num_nodes_window6, num_nodes_window6] centrality_cov_6;  // standard deviations of centrality estimates, time window 6
  matrix[num_nodes_window7, num_nodes_window7] centrality_cov_7;  // standard deviations of centrality estimates, time window 7
  matrix[num_nodes_window8, num_nodes_window8] centrality_cov_8;  // standard deviations of centrality estimates, time window 8
  matrix[num_nodes_window9, num_nodes_window9] centrality_cov_9;  // standard deviations of centrality estimates, time window 9
  matrix[num_nodes_window10, num_nodes_window10] centrality_cov_10;  // standard deviations of centrality estimates, time window 10
  matrix[num_nodes_window11, num_nodes_window11] centrality_cov_11;  // standard deviations of centrality estimates, time window 11
  matrix[num_nodes_window12, num_nodes_window12] centrality_cov_12;  // standard deviations of centrality estimates, time window 12
  matrix[num_nodes_window13, num_nodes_window13] centrality_cov_13;  // standard deviations of centrality estimates, time window 13
  matrix[num_nodes_window14, num_nodes_window14] centrality_cov_14;  // standard deviations of centrality estimates, time window 14
  matrix[num_nodes_window15, num_nodes_window15] centrality_cov_15;  // standard deviations of centrality estimates, time window 15
  matrix[num_nodes_window16, num_nodes_window16] centrality_cov_16;  // standard deviations of centrality estimates, time window 16
  matrix[num_nodes_window17, num_nodes_window17] centrality_cov_17;  // standard deviations of centrality estimates, time window 17
  matrix[num_nodes_window18, num_nodes_window18] centrality_cov_18;  // standard deviations of centrality estimates, time window 18
  matrix[num_nodes_window19, num_nodes_window19] centrality_cov_19;  // standard deviations of centrality estimates, time window 19
  matrix[num_nodes_window20, num_nodes_window20] centrality_cov_20;  // standard deviations of centrality estimates, time window 20
  matrix[num_nodes_window21, num_nodes_window21] centrality_cov_21;  // standard deviations of centrality estimates, time window 21
  matrix[num_nodes_window22, num_nodes_window22] centrality_cov_22;  // standard deviations of centrality estimates, time window 22
  matrix[num_nodes_window23, num_nodes_window23] centrality_cov_23;  // standard deviations of centrality estimates, time window 23
  matrix[num_nodes_window24, num_nodes_window24] centrality_cov_24;  // standard deviations of centrality estimates, time window 24
  matrix[num_nodes_window25, num_nodes_window25] centrality_cov_25;  // standard deviations of centrality estimates, time window 25
  matrix[num_nodes_window26, num_nodes_window26] centrality_cov_26;  // standard deviations of centrality estimates, time window 26
  matrix[num_nodes_window27, num_nodes_window27] centrality_cov_27;  // standard deviations of centrality estimates, time window 27
  matrix[num_nodes_window28, num_nodes_window28] centrality_cov_28;  // standard deviations of centrality estimates, time window 28
  matrix[num_nodes_window29, num_nodes_window29] centrality_cov_29;  // standard deviations of centrality estimates, time window 29
  matrix[num_nodes_window30, num_nodes_window30] centrality_cov_30;  // standard deviations of centrality estimates, time window 30
  matrix[num_nodes_window31, num_nodes_window31] centrality_cov_31;  // standard deviations of centrality estimates, time window 31
  matrix[num_nodes_window32, num_nodes_window32] centrality_cov_32;  // standard deviations of centrality estimates, time window 32
  matrix[num_nodes_window33, num_nodes_window33] centrality_cov_33;  // standard deviations of centrality estimates, time window 33
  matrix[num_nodes_window34, num_nodes_window34] centrality_cov_34;  // standard deviations of centrality estimates, time window 34
  matrix[num_nodes_window35, num_nodes_window35] centrality_cov_35;  // standard deviations of centrality estimates, time window 35
  matrix[num_nodes_window36, num_nodes_window36] centrality_cov_36;  // standard deviations of centrality estimates, time window 36
  
  // node IDs for all time windows
  array[num_nodes_window1] int nodes_window1;             // Node IDs, time window 1
  array[num_nodes_window2] int nodes_window2;             // Node IDs, time window 2
  array[num_nodes_window3] int nodes_window3;             // Node IDs, time window 3
  array[num_nodes_window4] int nodes_window4;             // Node IDs, time window 4
  array[num_nodes_window5] int nodes_window5;             // Node IDs, time window 5
  array[num_nodes_window6] int nodes_window6;             // Node IDs, time window 6
  array[num_nodes_window7] int nodes_window7;             // Node IDs, time window 7
  array[num_nodes_window8] int nodes_window8;             // Node IDs, time window 8
  array[num_nodes_window9] int nodes_window9;             // Node IDs, time window 9
  array[num_nodes_window10] int nodes_window10;             // Node IDs, time window 10
  array[num_nodes_window11] int nodes_window11;             // Node IDs, time window 11
  array[num_nodes_window12] int nodes_window12;             // Node IDs, time window 12
  array[num_nodes_window13] int nodes_window13;             // Node IDs, time window 13
  array[num_nodes_window14] int nodes_window14;             // Node IDs, time window 14
  array[num_nodes_window15] int nodes_window15;             // Node IDs, time window 15
  array[num_nodes_window16] int nodes_window16;             // Node IDs, time window 16
  array[num_nodes_window17] int nodes_window17;             // Node IDs, time window 17
  array[num_nodes_window18] int nodes_window18;             // Node IDs, time window 18
  array[num_nodes_window19] int nodes_window19;             // Node IDs, time window 19
  array[num_nodes_window20] int nodes_window20;             // Node IDs, time window 20
  array[num_nodes_window21] int nodes_window21;             // Node IDs, time window 21
  array[num_nodes_window22] int nodes_window22;             // Node IDs, time window 22
  array[num_nodes_window23] int nodes_window23;             // Node IDs, time window 23
  array[num_nodes_window24] int nodes_window24;             // Node IDs, time window 24
  array[num_nodes_window25] int nodes_window25;             // Node IDs, time window 25
  array[num_nodes_window26] int nodes_window26;             // Node IDs, time window 26
  array[num_nodes_window27] int nodes_window27;             // Node IDs, time window 27
  array[num_nodes_window28] int nodes_window28;             // Node IDs, time window 28
  array[num_nodes_window29] int nodes_window29;             // Node IDs, time window 29
  array[num_nodes_window30] int nodes_window30;             // Node IDs, time window 30
  array[num_nodes_window31] int nodes_window31;             // Node IDs, time window 31
  array[num_nodes_window32] int nodes_window32;             // Node IDs, time window 32
  array[num_nodes_window33] int nodes_window33;             // Node IDs, time window 33
  array[num_nodes_window34] int nodes_window34;             // Node IDs, time window 34
  array[num_nodes_window35] int nodes_window35;             // Node IDs, time window 35
  array[num_nodes_window36] int nodes_window36;             // Node IDs, time window 36

  // node ages for all individuals
  array[num_data] real node_age;         // Node ages (years)
  
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
  vector[num_nodes_window8] predictor_window8;
  for(i in 1:num_nodes_window8) {
    predictor_window8[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[8])] + window_random_effect[8] + node_random_effect[nodes_window8[i]];
  }
  vector[num_nodes_window9] predictor_window9;
  for(i in 1:num_nodes_window9) {
    predictor_window9[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[9])] + window_random_effect[9] + node_random_effect[nodes_window9[i]];
  }
  vector[num_nodes_window10] predictor_window10;
  for(i in 1:num_nodes_window10) {
    predictor_window10[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[10])] + window_random_effect[10] + node_random_effect[nodes_window10[i]];
  }
  vector[num_nodes_window11] predictor_window11;
  for(i in 1:num_nodes_window11) {
    predictor_window11[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[11])] + window_random_effect[11] + node_random_effect[nodes_window11[i]];
  }
  vector[num_nodes_window12] predictor_window12;
  for(i in 1:num_nodes_window12) {
    predictor_window12[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[12])] + window_random_effect[12] + node_random_effect[nodes_window12[i]];
  }
  vector[num_nodes_window13] predictor_window13;
  for(i in 1:num_nodes_window13) {
    predictor_window13[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[13])] + window_random_effect[13] + node_random_effect[nodes_window13[i]];
  }
  vector[num_nodes_window14] predictor_window14;
  for(i in 1:num_nodes_window14) {
    predictor_window14[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[14])] + window_random_effect[14] + node_random_effect[nodes_window14[i]];
  }
  vector[num_nodes_window15] predictor_window15;
  for(i in 1:num_nodes_window15) {
    predictor_window15[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[15])] + window_random_effect[15] + node_random_effect[nodes_window15[i]];
  }
  vector[num_nodes_window16] predictor_window16;
  for(i in 1:num_nodes_window16) {
    predictor_window16[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[16])] + window_random_effect[16] + node_random_effect[nodes_window16[i]];
  }
  vector[num_nodes_window17] predictor_window17;
  for(i in 1:num_nodes_window17) {
    predictor_window17[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[17])] + window_random_effect[17] + node_random_effect[nodes_window17[i]];
  }
  vector[num_nodes_window18] predictor_window18;
  for(i in 1:num_nodes_window18) {
    predictor_window18[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[18])] + window_random_effect[18] + node_random_effect[nodes_window18[i]];
  }
  vector[num_nodes_window19] predictor_window19;
  for(i in 1:num_nodes_window19) {
    predictor_window19[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[19])] + window_random_effect[19] + node_random_effect[nodes_window19[i]];
  }
  vector[num_nodes_window20] predictor_window20;
  for(i in 1:num_nodes_window20) {
    predictor_window20[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[20])] + window_random_effect[20] + node_random_effect[nodes_window20[i]];
  }
  vector[num_nodes_window21] predictor_window21;
  for(i in 1:num_nodes_window21) {
    predictor_window21[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[21])] + window_random_effect[21] + node_random_effect[nodes_window21[i]];
  }
  vector[num_nodes_window22] predictor_window22;
  for(i in 1:num_nodes_window22) {
    predictor_window22[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[22])] + window_random_effect[22] + node_random_effect[nodes_window22[i]];
  }
  vector[num_nodes_window23] predictor_window23;
  for(i in 1:num_nodes_window23) {
    predictor_window23[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[23])] + window_random_effect[23] + node_random_effect[nodes_window23[i]];
  }
  vector[num_nodes_window24] predictor_window24;
  for(i in 1:num_nodes_window24) {
    predictor_window24[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[24])] + window_random_effect[24] + node_random_effect[nodes_window24[i]];
  }
  vector[num_nodes_window25] predictor_window25;
  for(i in 1:num_nodes_window25) {
    predictor_window25[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[25])] + window_random_effect[25] + node_random_effect[nodes_window25[i]];
  }
  vector[num_nodes_window26] predictor_window26;
  for(i in 1:num_nodes_window26) {
    predictor_window26[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[26])] + window_random_effect[26] + node_random_effect[nodes_window26[i]];
  }
  vector[num_nodes_window27] predictor_window27;
  for(i in 1:num_nodes_window27) {
    predictor_window27[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[27])] + window_random_effect[27] + node_random_effect[nodes_window27[i]];
  }
  vector[num_nodes_window28] predictor_window28;
  for(i in 1:num_nodes_window28) {
    predictor_window28[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[28])] + window_random_effect[28] + node_random_effect[nodes_window28[i]];
  }
  vector[num_nodes_window29] predictor_window29;
  for(i in 1:num_nodes_window29) {
    predictor_window29[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[29])] + window_random_effect[29] + node_random_effect[nodes_window29[i]];
  }
  vector[num_nodes_window30] predictor_window30;
  for(i in 1:num_nodes_window30) {
    predictor_window30[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[30])] + window_random_effect[30] + node_random_effect[nodes_window30[i]];
  }
  vector[num_nodes_window31] predictor_window31;
  for(i in 1:num_nodes_window31) {
    predictor_window31[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[31])] + window_random_effect[31] + node_random_effect[nodes_window31[i]];
  }
  vector[num_nodes_window32] predictor_window32;
  for(i in 1:num_nodes_window32) {
    predictor_window32[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[32])] + window_random_effect[32] + node_random_effect[nodes_window32[i]];
  }
  vector[num_nodes_window33] predictor_window33;
  for(i in 1:num_nodes_window33) {
    predictor_window33[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[33])] + window_random_effect[33] + node_random_effect[nodes_window33[i]];
  }
  vector[num_nodes_window34] predictor_window34;
  for(i in 1:num_nodes_window34) {
    predictor_window34[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[34])] + window_random_effect[34] + node_random_effect[nodes_window34[i]];
  }
  vector[num_nodes_window35] predictor_window35;
  for(i in 1:num_nodes_window35) {
    predictor_window35[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[35])] + window_random_effect[35] + node_random_effect[nodes_window35[i]];
  }
  vector[num_nodes_window36] predictor_window36;
  for(i in 1:num_nodes_window36) {
    predictor_window36[i] = intercept + beta_age * node_age[(i + num_nodes_prev_windows[36])] + window_random_effect[36] + node_random_effect[nodes_window36[i]];
  }
  
}

model {
  // priors
  beta_age ~ normal(0,0.8);
  sigma ~ exponential(2);
  intercept ~ normal(logit(0.05),2); // 0.05 best estimate of intercept for Amboseli population eigen ~ age (edge weights from half weight index) from Chiyo 2011 -- estimated intercept from figure. Logit because on 0-1 scale 
  rand_window ~ normal(0,1); // rand_window ~ normal(0,sigma_window);
  sigma_window ~ exponential(2);
  rand_node ~ normal(0,1);   // rand_node ~ normal(0,sigma_node);
  sigma_node ~ exponential(2);
  
  // likelihood
  // centrality_mu_1 ~ multi_normal(predictor_window1, centrality_cov_1 + diag_matrix(rep_vector(sigma, num_nodes_window1)));
  // centrality_mu_2 ~ multi_normal(predictor_window2, centrality_cov_2 + diag_matrix(rep_vector(sigma, num_nodes_window2)));
  // centrality_mu_3 ~ multi_normal(predictor_window3, centrality_cov_3 + diag_matrix(rep_vector(sigma, num_nodes_window3)));
  // centrality_mu_4 ~ multi_normal(predictor_window4, centrality_cov_4 + diag_matrix(rep_vector(sigma, num_nodes_window4)));
  // centrality_mu_5 ~ multi_normal(predictor_window5, centrality_cov_5 + diag_matrix(rep_vector(sigma, num_nodes_window5)));
  // centrality_mu_6 ~ multi_normal(predictor_window6, centrality_cov_6 + diag_matrix(rep_vector(sigma, num_nodes_window6)));
  // centrality_mu_7 ~ multi_normal(predictor_window7, centrality_cov_7 + diag_matrix(rep_vector(sigma, num_nodes_window7)));
  // centrality_mu_8 ~ multi_normal(predictor_window8, centrality_cov_8 + diag_matrix(rep_vector(sigma, num_nodes_window8)));
  // centrality_mu_9 ~ multi_normal(predictor_window9, centrality_cov_9 + diag_matrix(rep_vector(sigma, num_nodes_window9)));
  // centrality_mu_10 ~ multi_normal(predictor_window10, centrality_cov_10 + diag_matrix(rep_vector(sigma, num_nodes_window10)));
  // centrality_mu_11 ~ multi_normal(predictor_window11, centrality_cov_11 + diag_matrix(rep_vector(sigma, num_nodes_window11)));
  // centrality_mu_12 ~ multi_normal(predictor_window12, centrality_cov_12 + diag_matrix(rep_vector(sigma, num_nodes_window12)));
  // centrality_mu_13 ~ multi_normal(predictor_window13, centrality_cov_13 + diag_matrix(rep_vector(sigma, num_nodes_window13)));
  // centrality_mu_14 ~ multi_normal(predictor_window14, centrality_cov_14 + diag_matrix(rep_vector(sigma, num_nodes_window14)));
  // centrality_mu_15 ~ multi_normal(predictor_window15, centrality_cov_15 + diag_matrix(rep_vector(sigma, num_nodes_window15)));
  // centrality_mu_16 ~ multi_normal(predictor_window16, centrality_cov_16 + diag_matrix(rep_vector(sigma, num_nodes_window16)));
  // centrality_mu_17 ~ multi_normal(predictor_window17, centrality_cov_17 + diag_matrix(rep_vector(sigma, num_nodes_window17)));
  // centrality_mu_18 ~ multi_normal(predictor_window18, centrality_cov_18 + diag_matrix(rep_vector(sigma, num_nodes_window18)));
  // centrality_mu_19 ~ multi_normal(predictor_window19, centrality_cov_19 + diag_matrix(rep_vector(sigma, num_nodes_window19)));
  // centrality_mu_20 ~ multi_normal(predictor_window20, centrality_cov_20 + diag_matrix(rep_vector(sigma, num_nodes_window20)));
  // centrality_mu_21 ~ multi_normal(predictor_window21, centrality_cov_21 + diag_matrix(rep_vector(sigma, num_nodes_window21)));
  // centrality_mu_22 ~ multi_normal(predictor_window22, centrality_cov_22 + diag_matrix(rep_vector(sigma, num_nodes_window22)));
  // centrality_mu_23 ~ multi_normal(predictor_window23, centrality_cov_23 + diag_matrix(rep_vector(sigma, num_nodes_window23)));
  // centrality_mu_24 ~ multi_normal(predictor_window24, centrality_cov_24 + diag_matrix(rep_vector(sigma, num_nodes_window24)));
  // centrality_mu_25 ~ multi_normal(predictor_window25, centrality_cov_25 + diag_matrix(rep_vector(sigma, num_nodes_window25)));
  // centrality_mu_26 ~ multi_normal(predictor_window26, centrality_cov_26 + diag_matrix(rep_vector(sigma, num_nodes_window26)));
  // centrality_mu_27 ~ multi_normal(predictor_window27, centrality_cov_27 + diag_matrix(rep_vector(sigma, num_nodes_window27)));
  // centrality_mu_28 ~ multi_normal(predictor_window28, centrality_cov_28 + diag_matrix(rep_vector(sigma, num_nodes_window28)));
  // centrality_mu_29 ~ multi_normal(predictor_window29, centrality_cov_29 + diag_matrix(rep_vector(sigma, num_nodes_window29)));
  // centrality_mu_30 ~ multi_normal(predictor_window30, centrality_cov_30 + diag_matrix(rep_vector(sigma, num_nodes_window30)));
  // centrality_mu_31 ~ multi_normal(predictor_window31, centrality_cov_31 + diag_matrix(rep_vector(sigma, num_nodes_window31)));
  // centrality_mu_32 ~ multi_normal(predictor_window32, centrality_cov_32 + diag_matrix(rep_vector(sigma, num_nodes_window32)));
  // centrality_mu_33 ~ multi_normal(predictor_window33, centrality_cov_33 + diag_matrix(rep_vector(sigma, num_nodes_window33)));
  // centrality_mu_34 ~ multi_normal(predictor_window34, centrality_cov_34 + diag_matrix(rep_vector(sigma, num_nodes_window34)));
  // centrality_mu_35 ~ multi_normal(predictor_window35, centrality_cov_35 + diag_matrix(rep_vector(sigma, num_nodes_window35)));
  // centrality_mu_36 ~ multi_normal(predictor_window36, centrality_cov_36 + diag_matrix(rep_vector(sigma, num_nodes_window36)));
  centrality_mu_1 ~ multi_student_t(num_nodes_window1, predictor_window1, centrality_cov_1 + diag_matrix(rep_vector(sigma, num_nodes_window1)));
  centrality_mu_2 ~ multi_student_t(num_nodes_window2, predictor_window2, centrality_cov_2 + diag_matrix(rep_vector(sigma, num_nodes_window2)));
  centrality_mu_3 ~ multi_student_t(num_nodes_window3, predictor_window3, centrality_cov_3 + diag_matrix(rep_vector(sigma, num_nodes_window3)));
  centrality_mu_4 ~ multi_student_t(num_nodes_window4, predictor_window4, centrality_cov_4 + diag_matrix(rep_vector(sigma, num_nodes_window4)));
  centrality_mu_5 ~ multi_student_t(num_nodes_window5, predictor_window5, centrality_cov_5 + diag_matrix(rep_vector(sigma, num_nodes_window5)));
  centrality_mu_6 ~ multi_student_t(num_nodes_window6, predictor_window6, centrality_cov_6 + diag_matrix(rep_vector(sigma, num_nodes_window6)));
  centrality_mu_7 ~ multi_student_t(num_nodes_window7, predictor_window7, centrality_cov_7 + diag_matrix(rep_vector(sigma, num_nodes_window7)));
  centrality_mu_8 ~ multi_student_t(num_nodes_window8, predictor_window8, centrality_cov_8 + diag_matrix(rep_vector(sigma, num_nodes_window8)));
  centrality_mu_9 ~ multi_student_t(num_nodes_window9, predictor_window9, centrality_cov_9 + diag_matrix(rep_vector(sigma, num_nodes_window9)));
  centrality_mu_10 ~ multi_student_t(num_nodes_window10, predictor_window10, centrality_cov_10 + diag_matrix(rep_vector(sigma, num_nodes_window10)));
  centrality_mu_11 ~ multi_student_t(num_nodes_window11, predictor_window11, centrality_cov_11 + diag_matrix(rep_vector(sigma, num_nodes_window11)));
  centrality_mu_12 ~ multi_student_t(num_nodes_window12, predictor_window12, centrality_cov_12 + diag_matrix(rep_vector(sigma, num_nodes_window12)));
  centrality_mu_13 ~ multi_student_t(num_nodes_window13, predictor_window13, centrality_cov_13 + diag_matrix(rep_vector(sigma, num_nodes_window13)));
  centrality_mu_14 ~ multi_student_t(num_nodes_window14, predictor_window14, centrality_cov_14 + diag_matrix(rep_vector(sigma, num_nodes_window14)));
  centrality_mu_15 ~ multi_student_t(num_nodes_window15, predictor_window15, centrality_cov_15 + diag_matrix(rep_vector(sigma, num_nodes_window15)));
  centrality_mu_16 ~ multi_student_t(num_nodes_window16, predictor_window16, centrality_cov_16 + diag_matrix(rep_vector(sigma, num_nodes_window16)));
  centrality_mu_17 ~ multi_student_t(num_nodes_window17, predictor_window17, centrality_cov_17 + diag_matrix(rep_vector(sigma, num_nodes_window17)));
  centrality_mu_18 ~ multi_student_t(num_nodes_window18, predictor_window18, centrality_cov_18 + diag_matrix(rep_vector(sigma, num_nodes_window18)));
  centrality_mu_19 ~ multi_student_t(num_nodes_window19, predictor_window19, centrality_cov_19 + diag_matrix(rep_vector(sigma, num_nodes_window19)));
  centrality_mu_20 ~ multi_student_t(num_nodes_window20, predictor_window20, centrality_cov_20 + diag_matrix(rep_vector(sigma, num_nodes_window20)));
  centrality_mu_21 ~ multi_student_t(num_nodes_window21, predictor_window21, centrality_cov_21 + diag_matrix(rep_vector(sigma, num_nodes_window21)));
  centrality_mu_22 ~ multi_student_t(num_nodes_window22, predictor_window22, centrality_cov_22 + diag_matrix(rep_vector(sigma, num_nodes_window22)));
  centrality_mu_23 ~ multi_student_t(num_nodes_window23, predictor_window23, centrality_cov_23 + diag_matrix(rep_vector(sigma, num_nodes_window23)));
  centrality_mu_24 ~ multi_student_t(num_nodes_window24, predictor_window24, centrality_cov_24 + diag_matrix(rep_vector(sigma, num_nodes_window24)));
  centrality_mu_25 ~ multi_student_t(num_nodes_window25, predictor_window25, centrality_cov_25 + diag_matrix(rep_vector(sigma, num_nodes_window25)));
  centrality_mu_26 ~ multi_student_t(num_nodes_window26, predictor_window26, centrality_cov_26 + diag_matrix(rep_vector(sigma, num_nodes_window26)));
  centrality_mu_27 ~ multi_student_t(num_nodes_window27, predictor_window27, centrality_cov_27 + diag_matrix(rep_vector(sigma, num_nodes_window27)));
  centrality_mu_28 ~ multi_student_t(num_nodes_window28, predictor_window28, centrality_cov_28 + diag_matrix(rep_vector(sigma, num_nodes_window28)));
  centrality_mu_29 ~ multi_student_t(num_nodes_window29, predictor_window29, centrality_cov_29 + diag_matrix(rep_vector(sigma, num_nodes_window29)));
  centrality_mu_30 ~ multi_student_t(num_nodes_window30, predictor_window30, centrality_cov_30 + diag_matrix(rep_vector(sigma, num_nodes_window30)));
  centrality_mu_31 ~ multi_student_t(num_nodes_window31, predictor_window31, centrality_cov_31 + diag_matrix(rep_vector(sigma, num_nodes_window31)));
  centrality_mu_32 ~ multi_student_t(num_nodes_window32, predictor_window32, centrality_cov_32 + diag_matrix(rep_vector(sigma, num_nodes_window32)));
  centrality_mu_33 ~ multi_student_t(num_nodes_window33, predictor_window33, centrality_cov_33 + diag_matrix(rep_vector(sigma, num_nodes_window33)));
  centrality_mu_34 ~ multi_student_t(num_nodes_window34, predictor_window34, centrality_cov_34 + diag_matrix(rep_vector(sigma, num_nodes_window34)));
  centrality_mu_35 ~ multi_student_t(num_nodes_window35, predictor_window35, centrality_cov_35 + diag_matrix(rep_vector(sigma, num_nodes_window35)));
  centrality_mu_36 ~ multi_student_t(num_nodes_window36, predictor_window36, centrality_cov_36 + diag_matrix(rep_vector(sigma, num_nodes_window36)));

}

