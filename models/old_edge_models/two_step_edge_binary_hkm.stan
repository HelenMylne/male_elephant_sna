functions {
  real conditional_prior_lpdf(real mu, int seen_together){
      if (seen_together == 1) {
      return beta_lpdf(mu | 2, 8);
    } else {
      return beta_lpdf(mu | 0.7, 10);
    }
  }
}

data {
  // for whole model, refer to https://github.com/JHart96/bison_examples/blob/main/models/binary_model.stan
  int<lower=0> n_dyads;             // Number of dyads
  int<lower=0> dyad_ids[n_dyads];   // Dyad ID corresponding to each data point
  int together[n_dyads];            // Total sightings of dyad in which they were together
  int count_dyad[n_dyads];          // Total sightings of dyad
  int<lower=0, upper=1> seen_together[n_dyads];       // Was dyad ever seen together? 1 = Yes, 0 = No.
  //int unsn_together[n_dyads];       // Was dyad ever seen together? 0 = Yes, 1 = No.
}

parameters {
  vector[n_dyads] logit_edge_weight; // Logit edge weights for each dyad that have been seen together at least once
}

transformed parameters {
  vector[n_dyads] logit_pn = logit_edge_weight[dyad_ids]; // Logit probability of a social event for each observation.
  //vector[n_dyads] mu = seen_together[dyad_ids];
}

model {
    // Main model
    together ~ binomial_logit(count_dyad, logit_pn);

    // Priors
    //logit_edge_weight ~ normal(mu, sigma);
    logit_edge_weight ~ conditional_prior_lpdf( mu, seen_together ) ; 
    
    //if(seen_together == 1){
    //  logit_edge_weight ~ normal(-1.5, 0.8);     // or beta(0.7, 10) for never together
    //} else { logit_edge_weight ~ normal(-3, 1) } // or beta(2,8) for sometimes together
}

generated quantities {
  int event_pred[n_dyads] = binomial_rng(count_dyad, inv_logit(logit_pn));
}
