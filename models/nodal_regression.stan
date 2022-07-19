data{
    int<lower=0> n_dyads;       // Number of dyads
    int<lower=0> n_eles;        // this isn't in Jordan's model, but the metric must surely need it
    vector[n_dyads] together;   // number of times dyad observed together
    vector[n_dyads] apart;      // number of times dyad observed in total
    vector[n_eles] age_cat_id;  // ages of individuals
}
parameters {
    vector<lower=0,upper=1>[n_dyads] weight; // edge weight to estimate from together vs apart
    vector[n_eles] betweenness;              // betweenness estimate for each individual
    real beta_0;                             // intercept of age variable -- should this be variable by individual or not?
    real beta_1;                             // slope of age -- no boundaries as could be positive or negative
    real<lower=0> sigma;                     // variance in betweenness
}
transformed parameters {
    real mu;                                 // betweenness mean
}
model {
  // edge weight
  weight ~ beta( 2 + together, 2 + apart );                    // weight values
  
  // nodal regression
  betweenness ~ normal_lpdf(mu, sigma); // is normal right here? Taken from Hart et al. 2021 PAPER (not github), which doesn't mention a specific centrality metric. Also check the _lpdf thing -- no idea if this is right or not
  mu = beta_0[n_eles] + beta_1[n_eles]*age_cat_id;
  
  // priors
  beta_0 ~ weight;                        // 
  beta_1 ~ normal_rng(0,1)                                     // weak prior focussed around 0: Chiyo >0, Goldenberg ≤0
  sigma ~ exponential(1);
}

// betweenness centrality = sum( total shortest paths from node A to node B that do NOT pass through C / total shortest paths from node A to node B )




// Jordan Hart count_model.stan -- this is NOT a nodal regression, just a poisson model of interaction strengths
data {
  int<lower=0> N; // Number of data points
  int<lower=0> M; // Number of dyads
  int<lower=0> L; // Number of locations
  int<lower=0> dyad_ids[N]; // Dyad ID corresponding to each data point
  int<lower=0> event_count[N]; // Outcome corresponding to each data point (presence/absence)
  int<lower=0> location_ids[N]; // Location ID corresponding to each data point
}

parameters {
  vector[M] log_p; // Logit edge weights for each dyad.
  vector[L] beta_loc; // Parameters for location effects.
  real<lower=0> loc_sigma; // Hyperparameter for location effect adaptive prior standard deviation.
}

transformed parameters {
  vector[N] log_pn = log_p[dyad_ids] + beta_loc[location_ids]; // Logit probability of a social event for each observation.
}

model {
  // Main model
  event_count ~ poisson(exp(log_pn));

  // Adaptive prior over location effects
  beta_loc ~ normal(0, loc_sigma);

  // Priors
  log_p ~ normal(0, 2.5);
  loc_sigma ~ normal(0, 1);
}

generated quantities {
  int event_pred[N] = poisson_rng(exp(log_pn));
}
