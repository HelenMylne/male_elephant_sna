data{
    int n_dyads;
    vector[n_dyads] apart;
    vector[n_dyads] together;
}
parameters {
    vector<lower=0,upper=1>[n_dyads] weight;
    real<lower=0> a;
    real<lower=0> b;
}
model {
  weight ~ beta( a + together, b + apart );
  a ~ normal(2,0.5);
  b ~ normal(2,0.5);
}
