data {
  int<lower=0> n_dyads;                   // Number of dyads
  array[n_dyads] int<lower=0> dyad_ids;   // Dyad ID corresponding to each data point
  array[n_dyads] int together;            // Total sightings of dyad in which they were together
  array[n_dyads] int count_dyad;          // Total sightings of dyad
  
}

parameters {
  vector<lower=0, upper=1>[n_dyads] edge_weight;      // edge weights for each dyad
  real<lower=0, upper=1> theta;
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
}
