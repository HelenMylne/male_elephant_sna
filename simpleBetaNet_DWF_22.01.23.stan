data{
    int n_dyads;
    vector[n_dyads] apart;
    vector[n_dyads] together;
}
parameters {
    vector<lower=0,upper=1>[n_dyads] weight;            // Outcome variable
}
model {
  for (n in 1:n_dyads){
    weight[n] ~ beta( 1 + together[n], 1 +  apart[n] );
  }
}
