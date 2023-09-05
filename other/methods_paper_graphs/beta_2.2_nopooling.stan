data{
    int n_dyads;
    vector[n_dyads] apart;
    vector[n_dyads] together;
}
parameters {
    vector<lower=0,upper=1>[n_dyads] weight; 
}
model {
  weight ~ beta( 2 + together, 2 + apart );
}
