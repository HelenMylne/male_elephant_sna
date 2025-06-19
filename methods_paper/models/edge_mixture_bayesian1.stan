functions {
  // Zero-inflated Binomial log-PMF
  real zero_inflated_binomial_lpmf(int y,       // observed count   (0 … n)
                                   int n,       // number of trials (≥ 0)
                                   real psi,    // structural–zero probability   (0–1)
                                   real p) {    // success probability in Binomial (0–1)
    // same as log( psi * I[y==0] + (1-psi) * Binom(y | n, p) )
    if (y == 0) {
      return log_mix(psi,
                     0,                         // log(1) for point-mass at 0
                     binomial_lpmf(0 | n, p));
    } else {
      return log1m(psi) +
             binomial_lpmf(y | n, p);
    }
  }
}

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
    count_together[dyad] ~ zero_inflated_binomial(total_sightings[dyad], hypothetical_together[dyad], edge_weight[dyad]); /*hypothetical_together ~ binomial(total_sightings, edge_weight);*/ // not convinced I've understood this function correctly??

  }
  
  // priors on probabilities
  logit_together ~ normal(0,5);
  dyad_effect ~ normal(0,1);
  average_edge ~ normal(0,1);
}
