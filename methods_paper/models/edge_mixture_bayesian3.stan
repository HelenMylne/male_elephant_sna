data {
  // number of data points
  int<lower=0> n_dyads;
  
  // observation data
  array[n_dyads] int<lower=0> count_together;
  array[n_dyads] int<lower=0> total_sightings;
  
}

parameters {
  // probability of structural zero
  real<lower=0, upper=1> prob_0;
  
  // probability of associating (population-level average and dyad-level deviation)
  array[n_dyads] real dyad_effect;
  real average_edge;
  
}

transformed parameters {
 
  // edge weight
  array[n_dyads] real edge_weight;
  edge_weight = inv_logit(average_edge + dyad_effect);
  
}

model {
  
  for(dyad in 1:n_dyads){
    
    if (count_together[dyad] == 0) {
    // step 1: probability of a structural zero
    target += log_mix(prob_0,
                      0,
                      binomial_lpmf(0 | total_sightings[dyad], edge_weight[dyad]));
    } else {
    // step 2: time spent associating given that they associate
      target += log1m(prob_0) + binomial_lpmf(count_together[dyad] | total_sightings[dyad], edge_weight[dyad] );
    }

  }
  
  // priors on probabilities
  prob_0 ~ beta(2,2);
  dyad_effect ~ normal(0,1);
  average_edge ~ normal(0,1);
}
