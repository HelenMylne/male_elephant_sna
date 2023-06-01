data {
  int<lower=0> n_dyads;                     // Number of dyads
  array[n_dyads] int<lower=0> dyad_ids;   // Dyad ID corresponding to each data point
  array[n_dyads] int together;            // Total sightings of dyad in which they were together
  array[n_dyads] int count_dyad;          // Total sightings of dyad
  array[n_dyads] int node_1;              // node ids for multimembership random effects
  array[n_dyads] int node_2;              // node ids for multimembership random effects
  array[n_dyads] int age_min;             // lower node age of dyad
  array[n_dyads] int age_max;             // upper node age of dyad
}

parameters {
  // edge weight
  vector<lower=0, upper=1>[n_dyads] edge_weight;      // edge weights for each dyad
  real<lower=0, upper=1> theta;                       // mixture prior multiplication factor
  
  // likelihood
  real<lower=0> sigma_weight;             // Error standard deviation
  
  // age effects
  real intercept;                         // intercept for model
  real b_min;                             // slope of effect of minimum age
  real b_max;                             // slope of effect of maximum age
  //real b_int;                             // slope of effect of interaction between minimum and maximum age
  
  // multimembership effects
  vector[n_dyads] mm;                     // multimembership effects
  real<lower=0> sigma_mm;                 // standard deviation of multimembership effects
}

model {
    // estimate probability of association from times together out of times observed
    together ~ binomial(count_dyad, edge_weight);
    // mixture prior for edge weight
    for(i in 1:n_dyads){
      target += log_mix(theta,
                        beta_lpdf(edge_weight[i] | 0.7, 10),
                        beta_lpdf(edge_weight[i] | 1, 5));
    }
    // prior for theta = flat prior 0-1
    theta ~ beta(1,1);
  
  // model
  for(i in 1:n_dyads) {
    logit(edge_weight[i]) ~ normal(intercept + b_min*age_min[i] + b_max*age_max[i] + mm[node_1[i]] + mm[node_2[i]], sigma_weight);  // b_int*age_min*age_max + 
  }
  
  // priors for likelihood
  sigma_weight ~ exponential(1);     // Cauchy prior for error standard deviation sigma
  
  // priors for age effects
  intercept ~ normal(-2,1);
  b_min ~ normal(0,0.1);
  b_max ~ normal(0,0.1);
  //b_int ~ normal(0,0.05);
  
  // priors for multimembership effects
  mm ~ normal(0, sigma_mm);
  sigma_mm ~ exponential(1);

}

