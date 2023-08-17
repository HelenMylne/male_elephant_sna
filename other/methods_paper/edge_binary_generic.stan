data {
  int<lower=0> n_dyads;                   // Number of dyads
  array[n_dyads] int<lower=0> dyad_ids;   // Dyad ID corresponding to each data point
  array[n_dyads] int together;            // Total sightings of dyad in which they were together
  array[n_dyads] int count_dyad;          // Total sightings of dyad
  int prior_param1_0;
  int prior_param2_0;
  int prior_param1_1;
  int prior_param2_1;
  int model_type; // 1 = logit gaussian, 2 = beta
}

parameters {
  vector<lower=0, upper=1>[n_dyads] edge_weight;      // edge weights for each dyad.
}

model {
    for (i in 1:n_dyads) {
        // Conditional priors
        if (together[i] == 0)
            edge_weight[i] ~ beta(prior_param1_0, prior_param2_0);
        else
            edge_weight[i] ~ beta(prior_param1_1, prior_param2_1);
        }

    together ~ binomial(count_dyad, edge_weight);

}
