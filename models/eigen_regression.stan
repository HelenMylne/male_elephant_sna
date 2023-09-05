data {
  int num_nodes;                    // Number of nodes
  vector[num_nodes] nodes;          // Node IDs
  vector[num_nodes] centrality_mu;  // Means of centrality estimates
  vector[num_nodes] centrality_sd;  // standard deviations of centrality estimates
  vector[num_nodes] node_age;       // Node ages (point estimate)
  //vector[num_nodes] node_age2;      // Node ages squared (point estimate)
}

parameters {
  real beta_age;
  //real beta_age2;
  real<lower=0> sigma;
  //real intercept;
}

transformed parameters {
  vector[num_nodes] predictor;
  predictor = beta_age*node_age;// + beta_age2*node_age2; //intercept + 
}

model {
  // priors
  beta_age ~ normal(0, 0.02);
  //beta_age2 ~ normal(0, 0.003);
  sigma ~ exponential(1);
  //intercept ~ normal(0, 0.05);
  
  // likelihood
  centrality_mu ~ normal(predictor, centrality_sd);
}
