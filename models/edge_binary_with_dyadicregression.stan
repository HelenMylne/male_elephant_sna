data {
  int<lower=0> n_dyads;                     // Number of dyads
  array[n_dyads] int<lower=0> dyad_ids;   // Dyad ID corresponding to each data point
  array[n_dyads] int together;            // Total sightings of dyad in which they were together
  array[n_dyads] int count_dyad;          // Total sightings of dyad
  array[n_dyads] int node_1;              // node ids for multimembership random effects
  array[n_dyads] int node_2;              // node ids for multimembership random effects
  array[n_dyads] int age_min;             // lower node age of dyad
  array[n_dyads] int age_max;             // upper node age of dyad
  //int<lower=0> dyad_ids[n_dyads];           // Dyad ID corresponding to each data point
  //int together[n_dyads];                    // Total sightings of dyad in which they were together
  //int count_dyad[n_dyads];                  // Total sightings of dyad
  //int node_1[n_dyads];                      // node ids for multimembership random effects
  //int node_2[n_dyads];                      // node ids for multimembership random effects
  //int age_min[n_dyads];                     // lower node age of dyad
  //int age_max[n_dyads];                     // upper node age of dyad
  //int<lower=0> n_samples[n_dyads];          // number of samples to run in model
}

parameters {
  // edge weight
  vector<lower=0, upper=1>[n_dyads] edge_weight;      // edge weights for each dyad
  
  // likelihood
  real<lower=0> sigma_weight;             // Error standard deviation
  
  // age effects
  real intercept;                         // intercept for model
  real b_min;                             // slope of effect of minimum age
  real b_max;                             // slope of effect of maximum age
  real b_int;                             // slope of effect of interaction between minimum and maximum age
  
  // multimembership effects
  vector[n_dyads] mm;                     // multimembership effects
  real<lower=0> sigma_mm;                 // standard deviation of multimembership effects
}

model {
    // priors for edge weight (conditional -- based on if ever seen together)
    for (i in 1:n_dyads) {
        if (together[i] == 0)
            edge_weight[i] ~ beta(0.7, 10);
        else
            edge_weight[i] ~ beta(1, 5);
        }
    // produce edge model
    together ~ binomial(count_dyad, edge_weight);   // estimate edge weight from sightings together out of total sightings
  
  logit(edge_weight) ~ normal(intercept + b_min*age_min + b_max*age_max + mm[node_1] + mm[node_2], sigma_weight);  // b_int*age_min*age_max + 
  
  // priors for likelihood
  sigma_weight ~ exponential(1);     // Cauchy prior for error standard deviation sigma
  
  // priors for age effects
  intercept ~ normal(0,1);
  b_min ~ normal(0,0.05);
  b_max ~ normal(0,0.05);
  //b_int ~ normal(0,0.05);
  
  // priors for multimembership effects
  mm ~ normal(0, sigma_mm);
  sigma_mm ~ exponential(1);

}
