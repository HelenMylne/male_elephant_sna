data{
    int n_dyads;                    // number of dyads in model
    vector[n_dyads] apart;          // number of sightings of dyad members where they weren't together
    vector[n_dyads] together;       // number of sightings of dyad members where they were together
    //vector[n_dyads] count_dyad;     // total number of sightings
    //vector[n_dyads] seen_together;  // binary variable: 1 indicates together > 0, 0 indicates together = 0
    vector[n_dyads] alpha;          // start off by just coding which alpha parameter to use in data
    vector[n_dyads] beta;           // start off by just coding which beta  parameter to use in data
}
parameters {
    vector<lower=0,upper=1>[n_dyads] weight;
    // real sigma; // partial pooling parameter
}
model {
  weight ~ beta( together + alpha , apart + beta );
}



