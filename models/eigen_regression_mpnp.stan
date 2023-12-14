data {
  //int num_data;                               // Numer of total data points (same node in two windows = 2)
  int num_nodes;                                // Number of unique nodes across all windows (same node in two windows = 1)
  int num_age_cat;                              // Number of unique age categories
  int length_dirichlet;                         // Number of unique age categories + 1
  //int num_windows;                              // Number of time windows
  vector[num_nodes] centrality_mu;              // Means of centrality estimates (0-1 bounded then standardised)
  matrix[num_nodes, num_nodes] centrality_cov;   // standard deviations of centrality estimates
  array[num_nodes] int node_age;                // Node ages (factor)
  //array[num_data] int window;                   // Window ID per data point
  //array[num_data] int node;                     // Node ID per data point (same individual in 2 windwos gets same value)
  vector[num_age_cat] prior_age;               // Dirichlet prior values
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
  //real mu_window;
  //real sigma_window;
  //vector[num_nodes] rand_node;
  //real mu_node;
  //real sigma_node;
  // difference between age categories
  simplex[num_age_cat] delta;
}

transformed parameters {
  // create prior for cumulative probability of each age category
  vector[length_dirichlet] delta_j;
  delta_j = append_row(0, delta);
  // linear model
  vector[num_nodes] predictor;
  for(i in 1:num_nodes) {
    predictor[i] = intercept + beta_age * sum(delta_j[1:node_age[i]]);// + rand_window[window[i]] + rand_node[node[i]];
  }
}

 model {
  // priors
  intercept ~ normal(0,0.8);
  delta ~ dirichlet(prior_age);
  beta_age  ~ normal(0,1);
  sigma ~ exponential(2);
  //rand_window ~ normal(mu_window,sigma_window);
  //mu_window ~ normal(0,1);
  //sigma_window ~ exponential(1);
  //rand_node ~ normal(mu_node,sigma_node);
  //mu_node ~ normal(0,1);
  //sigma_node ~ exponential(1);

  // likelihood
  centrality_mu ~ multi_normal(predictor, centrality_cov + diag_matrix(rep_vector(sigma, num_nodes)));
}
