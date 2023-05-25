data {
  int<lower=0> n_dyads;                   // Number of dyads
  array[n_dyads] int<lower=0> dyad_ids;   // Dyad ID corresponding to each data point
  array[n_dyads] int together;            // Total sightings of dyad in which they were together
  array[n_dyads] int count_dyad;          // Total sightings of dyad
  array[n_dyads] int node_1;              // node ids for multimembership random effects
  array[n_dyads] int node_2;              // node ids for multimembership random effects
  array[n_dyads] int age_min;             // lower node age of dyad
  array[n_dyads] int age_max;             // upper node age of dyad
}

parameters {
  vector<lower=0, upper=1>[n_dyads] edge_weight;      // edge weights for each dyad
  vector[n_dyads] mm;                                 // multimembership effects
  real<lower=0> sigma_mm;                             // standard deviation of multimembership effects
}

model {
    for (i in 1:n_dyads) {
        // Conditional priors
        if (together[i] == 0)
            edge_weight[i] ~ beta(0.7, 10);
        else
            edge_weight[i] ~ beta(1, 5);
        }
    // produce edge model
    together ~ binomial(count_dyad, edge_weight);
  // take values from edge model and use in dyadic regression
  for (d in 1:n_dyads) {
    i = node_1[d];
    j = node_2[d];
    edge_weight ~ age_min + age_max + age_min*age_max + mm[i] + mm[j];
    //edge_weight ~ b_min*age_min + b_max*age_max + b_int*age_min*age_max + mm[node_1[d]] + mm[node_2[d]];
    //b_min ~ normal(0,1);
    //b_max ~ normal(0,1);
    //b_int ~ normal(0,1);
  }
  mm ~ gaussian(0, sigma_mm);
  sigma_mm ~ exponential(1);
}

