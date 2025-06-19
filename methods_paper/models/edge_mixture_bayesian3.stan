data {
  // number of data points
  int<lower=0> n_dyads;
  
  // observation data
  array[n_dyads] int ever_together;
  array[n_dyads] int count_together;
  array[n_dyads] int total_sightings;
  
  // dyad ID
  array[n_dyads] int dyad_id;
  
}

parameters {
  // probability of associating at all
  real logit_together;
  
  // time spent associating given that they do associate
  array[n_dyads] real dyad_effect;
  real average_edge;
  
  // hypothetical number of interactions dyad WOULD have IF socially connected -- this feels like it should be an integer value if I've understood Dan correctly about what it was representing, but parameters can't be integers. Also threw errors if I allowed it to be outside range of [0,1]. Not entirely sure what this is actually estimating!
  array[n_dyads] real<lower=0,upper=1> hypothetical_together;
  
}

transformed parameters {
  // probability of associating at all
  real prob_together;
  prob_together = inv_logit(logit_together);
  
  // time spent associating given that they do associate
  array[n_dyads] real logit_edge;
  array[n_dyads] real edge_weight;
  for(i in 1:n_dyads){
    logit_edge[i] = average_edge + dyad_effect[i];
    edge_weight[i] = inv_logit(logit_edge[i]);
  }
  
}

model {
  
  for(dyad in 1:n_dyads){
    
    // step 1: probability of ever associating
    ever_together[dyad] ~ bernoulli(prob_together);
  
    // step 2: time spent associating given that they associate
    if (count_together[dyad] == 0) {
      target += log_sum_exp(log(prob_together),
                            log1m(prob_together)
                              + binomial_lpmf(count_together[dyad] | total_sightings[dyad], edge_weight[dyad] ));
    } else {
      target += log1m(prob_together)
                  + binomial_lpmf(count_together[dyad] | total_sightings[dyad], edge_weight[dyad] );
  }

  }
  
  // priors on probabilities
  logit_together ~ normal(0,5);
  dyad_effect ~ normal(0,1);
  average_edge ~ normal(0,1);
}
