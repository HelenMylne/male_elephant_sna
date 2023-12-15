data {
  int num_nodes;                                // Number of unique nodes across all windows (same node in two windows = 1)
  int num_age_cat;                              // Number of unique age categories
  int length_dirichlet;                         // Number of unique age categories + 1
  vector[num_nodes] centrality_mu;              // Means of centrality estimates (0-1 bounded then standardised)
  matrix[num_nodes, num_nodes] centrality_cov;  // standard deviations of centrality estimates
  array[num_nodes] int node_age;                // Node ages (factor)
  vector[num_age_cat] prior_age;                // Dirichlet prior values
}

parameters {
  // intercept
  real intercept;
  // exposure slope
  real beta_age;
  //variance
  real<lower=0> sigma;
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
    predictor[i] = intercept + beta_age * sum(delta_j[1:node_age[i]]);
  }
}

 model {
  // priors
  intercept ~ normal(0,1);
  delta ~ dirichlet(prior_age);
  beta_age  ~ normal(0,1);
  sigma ~ exponential(2);

  // likelihood
  centrality_mu ~ multi_normal(predictor, centrality_cov + diag_matrix(rep_vector(sigma, num_nodes)));
}
