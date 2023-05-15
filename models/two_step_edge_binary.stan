data{
    int n_dyads;                    // number of dyads in model
    vector[n_dyads] apart;          // number of sightings of dyad members where they weren't together
    vector[n_dyads] together;       // number of sightings of dyad members where they were together
    vector[n_dyads] count_dyad;     // total number of sightings
    vector[n_dyads] seen_together;  // binary vairable: 1 indicates together > 0, 0 indicates together = 0
    vector[n_dyads] alpha;          // start off by just coding which alpha parameter to use in data
    vector[n_dyads] beta;           // start off by just coding which beta  parameter to use in data
}
parameters {
    vector<lower=0,upper=1>[n_dyads] weight;
    //vector[n_dyads] p; // probability of being observed together
    // real<lower=0,upper=1> sigma; // partial pooling parameter -- add this in once you've got the rest of it working
}
transformed parameters {
    vector<lower=0,upper=1>[n_dyads] p;
    //vector<lower=0,upper=1>[n_dyads] alpha;
    //vector<lower=0,upper=1>[n_dyads] beta;
   // alpha = if(seen_together == 1) { 2 } else { 0.7 } ;
   // beta  = if(seen_together == 1) { 8 } else { 10 };
}
model {
  //weight ~ binomial(count_dyad, p); // binomial won't accept having weight or count_dyad as vectors, must be integer
  
  for( i in 1:n_dyads) {
    weight[i] ~ binomial(count_dyad[i], p);   // "real return type required for probability function"
    p ~ beta(alpha[i], beta[i]);
  }
   // p ~ beta( 2 + together, 8 + apart )*seen_together + beta( 0.7 + together, 10 + apart )*(1-seen_together); // can't recognise beta() as probability function when not the first thing after the ~
  
//ME MESSING AROUND 5:
  //if (seen_together == 0) {
  //    a = 0.7;
  //    b = 10;
  //  } else {
  //    a = 2;
  //    b = 8;
  //  }
  //weight ~ beta( together + a , apart + b ) ;

//ME MESSING AROUND 4:
  //weight ~ beta_lpdf( together | (seen_together * count_dyad) ) ; // can't cope with *, removing "* count_dyad" gives "no matches for: available argument signature for beta_lpdf:"
  //  if (seen_together == 0) {
  //   return beta_lpdf(mu | 0.7, 10);
  //  } else {
  //    return beta_lpdf(mu | 2, 8);
  //  }
  //}

//ME MESSING AROUND 3 -- adapting Dan's code from a couple of weeks ago:
  //real conditional_prior_lpdf(real mu, int x) { ;  // parser expected ";"??
  //  if (seen_together == 0) {
  //   return beta_lpdf(mu | 0.7, 10);
  //  } else {
  //    return beta_lpdf(mu | 2, 8);
  //  }
  //}

//ME MESSING AROUND 2:
  //for(i in 1:n_dyads){
  //   if(seen_together[i] == 1){
  //     p ~ beta( 2 + together, 8 + apart );
  //   } else {
  //     p ~ beta( 0.7 + together, 10 + apart );
  //   }
  //   weight ~ binomial(count_dyad, p);  // "real return type required for probability function" 
  // }

//ME MESSING AROUND 1:
  //weight ~ binomial(count_dyad, p);   // "real return type required for probability function"
  //p ~ beta( 2 + together, 8 + apart )*seen_together + beta( 0.7 + together, 10 + apart )*(1-seen_together); // can't recognise beta() as probability function when not the first thing after the ~

//ORIGINAL AND WORKING:
  // weight ~ beta( 2 + together, 2 + apart );
  

}



