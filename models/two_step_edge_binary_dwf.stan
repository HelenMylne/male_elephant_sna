data {
  // for whole model, refer to https://github.com/JHart96/bison_examples/blob/main/models/binary_model.stan
  int<lower=0> n_dyads;                   // Number of dyads
  array[n_dyads] int<lower=0> dyad_ids;   // Dyad ID corresponding to each data point
  array[n_dyads] int together;            // Total sightings of dyad in which they were together
  array[n_dyads] int count_dyad;          // Total sightings of dyad
}

parameters {
  vector<lower=0, upper=1>[n_dyads] edge_weight;      // edge weights for each dyad.
}

transformed parameters {
  vector[n_dyads] pn = edge_weight[dyad_ids];         // Probability of a social event to the correct dyad
}

model {
    for (i in 1:n_dyads) {
        // Conditional priors
        if (together[i] == 0)
            edge_weight[i] ~ beta(0.7, 10);
        else
            edge_weight[i] ~ beta(2, 8);
    }

    together ~ binomial(count_dyad, pn);
}

generated quantities {
  array[n_dyads] int event_pred = binomial_rng(count_dyad, pn);
}

