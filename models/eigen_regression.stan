data {
  int num_nodes;                    // Number of nodes
  vector[num_nodes] nodes;          // Node IDs
  vector[num_nodes] centrality_mu;  // Means of centrality estimates
  vector[num_nodes] centrality_sd;  // standard deviations of centrality estimates
  vector[num_nodes] node_age;       // Node ages (point estimate)
}

parameters {
  real beta_age;
  real<lower=0> sigma;
  real centrality[num_nodes];
}

transformed parameters {
  vector[num_nodes] predictor;
  predictor = beta_age*node_age;
}

model {
  // priors
  beta_age ~ normal(0, 0.1);
  
  // prior for sigma -- is sigma actually necessary here?? for dyadic it was part of the covariance matrix, but here there is no covariance matrix involved so I think it only needs the standard deviation of the centrality?
  //sigma ~ exponential(1);
  
  // likelihood
  logit(centrality_mu) ~ normal(predictor, centrality_sd);//sigma);
  // for( i in 1:num_nodes ){
  //  inv_logit(centrality[i]) ~ normal(centrality_mu[i], centrality_sd[i]);
  //  centrality_mu[i] ~ normal(predictor[i], sigma);
  // }
  
  
}
