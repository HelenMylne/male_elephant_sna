data{
  int n_dyads;                          // number of dyads in model
  int<lower=0> together[n_dyads];                // number of interactions
  //vector[n_dyads] together;           // number of sightings of dyad members where they were together
  int<lower=0> count_dyad[n_dyads];              // number of opportunities to interact
  //vector[n_dyads] count_dyad;         // total number of sightings
  int<lower=0,upper=1> seen_together[n_dyads];   // condition for choosing prior
  //vector[n_dyads] seen_together;      // binary vairable: 1 indicates together > 0, 0 indicates together = 0
}

parameters {
  vector[n_dyads] weight;  // probability of association
}

model {
  for(i in 1:n_dyads){
  // prior
    if (seen_together[i] == 0) {
      inv_logit(weight) ~ beta(0.7, 10);
    } else {
      inv_logit(weight) ~ beta(2, 8);
    }
  
  // likelihood
  together[i] ~ binomial(count_dyad[i], inv_logit(weight[i])); // binomial likelihood with beta prior
  }
}








data{
  int n_dyads;                          // number of dyads in model
  int<lower=0> together[n_dyads];                // number of interactions
  //vector[n_dyads] together;           // number of sightings of dyad members where they were together
  int<lower=0> count_dyad[n_dyads];              // number of opportunities to interact
  //vector[n_dyads] count_dyad;         // total number of sightings
  int<lower=0,upper=1> seen_together[n_dyads];   // condition for choosing prior
  //vector[n_dyads] seen_together;      // binary vairable: 1 indicates together > 0, 0 indicates together = 0
}

parameters {
  real<lower=0,upper=1> weight;  // probability of association
}

model {
  //for(i in 1:n_dyads){
  // prior
  //  if (seen_together[i] == 0) {
  //    weight ~ beta(0.7, 10);
  //  } else {
  //    weight ~ beta(2, 8);
  //  }
  
  // likelihood
  //together[i] ~ binomial(count_dyad[i], weight); // binomial likelihood with beta prior
  //}

  // prior
  if (seen_together == 0) {
      weight ~ beta(0.7, 10);
    } else {
      weight ~ beta(2, 8);
    }
  
  // likelihood
  together ~ binomial(count_dyad, weight); // binomial likelihood with beta prior
}


