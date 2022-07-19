data{
    int<lower=0> n_dyads;       // Number of dyads
    vector[n_dyads] together;   // number of times dyad observed together
    vector[n_dyads] apart;      // number of times dyad observed in total
    vector[n_dyads] age_diff;   // age difference between dyad pair
}
parameters {
    vector<lower=0,upper=1>[n_dyads] events; // edge weight to estimate from together vs apart  -- measurement of association
    vector<lower=0,upper=1>[n_dyads] weight; // edge weight to estimate regression coefficients -- probability of association
    //real beta_0;                           // intercept of age variable -- not variable by individual as want to see overall population trend
    real beta_1;                             // slope of age -- no boundaries as could be positive or negative
    real<lower=0> sigma;                     // variance in weights
    //vector[n_dyads] mu;                    // weight mean
    //vector[n_eles] u_i;                    // elephant-specific probability of associating with others of different age -- see comments below
    //vector[n_eles] u_j;                    // elephant-specific probability of associating with others of different age -- see comments below
}

model {
  vector[n_dyads] mu;
  
  // priors
  beta_1 ~ normal(-0.5,2);               // weak prior focussed slightly below zero -- expect a wider age gap to induce lower rates of interaction but allow for young males to prefer interacting with older ones.
  sigma ~ exponential(1);
  //u_i[n_eles] ~ normal(0,1); -- see comment below
  //u_j[n_eles] ~ normal(0,1); -- see comment below

  // dyadic regression
  for( i in 1:n_dyads) {
    mu[i] = events[i] + beta_1*age_diff[i]; // + u_i[n_eles[i]] + u_j[n_eles[i]];  -- should be something included here about individual gregariousness, but I can't work out how to do this so that the order of which is id_1 and which is id_2 is not important.
    weight ~ normal(mu, sigma);          // previously was complaining that "normal" was not an option and needed to specify normal_lpdf(), now telling me I shouldn't include the _lpdf and should just put normal()
  }
  
  // edge weight
  events ~ beta( 1 + together, 1 + apart );  // weight values
}

