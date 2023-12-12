data {
  int num_nodes;                    // Number of nodes
  int num_age_cat;                  // Number of unique age categories
  int length_dirichlet;             // number of unique age categories + 1
  vector[num_nodes] nodes;          // Node IDs
  vector[num_nodes] centrality_mu;  // Means of centrality estimates (0-1 bounded then standardised)
  matrix[num_nodes, num_nodes] centrality_cov;  // standard deviations of centrality estimates
  int node_age[num_nodes];       // Node ages (point estimate)
  vector[num_age_cat] prior_age;    // dirichlet prior values (0.5,0.5,0.5,0.5,0.5)
  //vector[num_age_cat + 1] delta_j;  // append zero to dirichlet prior (0, 0.5,0.5,0.5,0.5,0.5)
}

parameters {
  real beta_age;
  real<lower=0> sigma;
  real intercept;
  simplex[num_age_cat] delta;
}

transformed parameters {
  // create prior for cumulative probability of each age category
  vector[length_dirichlet] delta_j;
  delta_j = append_row(0, delta);
  // linear model
  vector[num_nodes] predictor;
  for(i in 1:num_nodes){
    predictor[i] = intercept + beta_age*sum(delta_j[1:node_age[i]]);
  }
}

model {
  // priors
  delta ~ dirichlet(prior_age);
  beta_age ~ normal(0,1);
  sigma ~ exponential(1);
  intercept ~ normal(0,0.8);
  
  // likelihood
  centrality_mu ~ multi_normal(predictor, centrality_cov + diag_matrix(rep_vector(sigma, num_nodes)));
}
