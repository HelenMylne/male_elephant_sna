data {
  int num_dyads; // number of edges
  int num_nodes; // number of nodes
  vector[num_dyads] edge_mu;             // means of Gaussian approximation of logit edge weights
  matrix[num_dyads, num_dyads] edge_cov; // covariance matrix of Gaussian approximation of logit edge weights
  vector[num_dyads] age_min;    // age of younger member of dyad
  vector[num_dyads] age_max;    // age of  older  member of dyad
  //int node_1[num_dyads];   // node 1 for multimembership
  //int node_2[num_dyads];   // node 2 for multimembership
}

parameters {
  // age
  real b_min;
  real b_max;
  real b_int;

  // multivariate normal
  real<lower=0> sigma;
  
  // multimembership
  //vector[num_nodes] mm_nodes;
  //real<lower=0> sigma_mm;
}

transformed parameters {
  vector[num_dyads] age_effect;
  for(i in 1:num_dyads){
    age_effect[i] = b_min*age_min[i] + b_max*age_max[i] + b_int*(age_min[i]*age_max[i]);// + mm_nodes[node_1] + mm_nodes[node_2]
  }
}

model {
  // model
  edge_mu ~ multi_normal(age_effect, edge_cov + diag_matrix(rep_vector(sigma, num_dyads)));
  
  // age priors
  b_min ~ normal(0,0.1);
  b_max ~ normal(0,0.1);
  b_int ~ normal(0,0.05);
  
  // multivariate normal prior
  sigma ~ exponential(1); // PREVIOUSLY NORMAL(0,1)
  
  // multimembership priors
  //mm_nodes ~ normal(0, sigma_mm);
  //sigma_mm ~ exponential(1);
}

