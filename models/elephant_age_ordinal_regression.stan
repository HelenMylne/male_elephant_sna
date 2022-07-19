functions{
//gompertz_bathtub <- function(a0, a1, c, b0, b1, age){
//  out <- exp(a0 - a1*age) + c + exp(b0 + b1*age)

//  vector gompertz_bathtub(real a0, real a1, real c, real b0, real b1, vector age){
//    vector[num_elements(age)] prob;
//    for(i in 1:num_elements(age)){
//      prob[i] = ( exp(b0 + b1*age[i]) + c + exp(a0 - a1*age[i]) );
//    }
//    return(prob);
//  }

 real gompertz_bathtub_lpdf(real a0, real a1, real c, real b0, real b1, real age){
    real prob;
    prob = ( exp(b0 + b1*age) + c + exp(a0 - a1*age) );
    return(prob);
  }
}

data {
    int N; // number of individuals
    int K; // number of age categories
    int<lower=1, upper=K> age_category_index[N]; // age categories
}

parameters {
  vector<lower=0>[N] observed_age_std;
  real<lower=0> true_age;
  real<lower=0> sigma_age;
  real a0;
  real a1;
  real c;
  real b0;
  real b1;
  //real independence_age;
}

transformed parameters {
  ordered[K-1] thresholds;
  vector<lower=0>[N] observed_age;
  // Thresholds for age classes
  thresholds[1] = 5;
  thresholds[2] = 10;
  thresholds[3] = 15;
  thresholds[4] = 20;
  thresholds[5] = 25;
  thresholds[6] = 40;
  // Non-centred age
  observed_age = true_age + sigma_age*observed_age_std;// + independence_age; // if one sigma for all
}

model {
  for(i in 1:N) {
    age_category_index[i] ~ ordered_logistic(observed_age[i], thresholds);
  }
  sigma_age ~ exponential(0.5);
  true_age ~ gompertz_bathtub(a0, a1, c, b0, b1);
   // estimates and standard errors from basta model of Amboseli population
  a0 ~ normal(-5.13, 0.72);
  a1 ~ normal( 3.00, 0.10);
  c  ~ normal( 0.026, 0.006);
  b0 ~ normal(-5.08, 0.56);
  b1 ~ normal( 0.09, 0.018);
  // literature value for age of becoming independent = 10-18, peak at 14
  //independence_age ~ normal(14, 2);
}
