
data {
  int<lower=0> N;                   // Number of dyads
  array[N] int together;            // Outcome: total sightings of dyad in which they were together
  array[N] int count_dyad;          // Total sightings of dyad
}

parameters {
  vector<lower=0, upper=1>[N] edge_weight;      // edge weights for each dyad.
}

model {
  for (i in 1:N) {
    // Conditional prior on edge weight
    if (together[i] == 0)
        edge_weight[i] ~ beta(0.7, 10);
    else
        edge_weight[i] ~ beta(2, 8);
  }

  together ~ binomial(count_dyad, edge_weight);
}
