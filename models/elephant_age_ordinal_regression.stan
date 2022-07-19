functions {
  real gompertz_bathtub_lpdf(vector age, real a0, real a1, real c, real b0, real b1) {
    // 1 - gompertz to give survival to an age
    return( sum(log(1.0 - (exp(b0 + b1*age) + c + exp(a0 - a1*age)))) );
  }
}

data {
    int N; // number of individuals
    int K; // number of age categories
    int<lower=1, upper=K> age_category_index[N]; // age categories
}

parameters {
  vector<lower=0>[N] observed_age_std;
  vector<lower=0>[N] true_age;
  // real<lower=0> sigma_age[N];
  real<lower=0> sigma_age;
  real a0_std;
  real a1_std;
  real b0_std;
  real b1_std;
  real c_std;
}

transformed parameters {
  ordered[K-1] thresholds;
  vector<lower=0>[N] observed_age;
  real a0;
  real a1;
  real b0;
  real b1;
  real c;

  // Thresholds for age classes
  thresholds[1] = 15;
  thresholds[2] = 30;
  thresholds[3] = 45;
  thresholds[4] = 60;
  thresholds[5] = 75;
  // Non-centred age. The same as observed_age ~ normal(true_age,sigma_age)
  observed_age = true_age + sigma_age*observed_age_std; // if one sigma for all
  
  // In a loop to have a different sigma for each individual
  //for(i in 1:N) {
  //  observed_age[i] = true_age[i] + sigma_age[i]*observed_age_std[i];
  //}
  // Non-centred shape. 

  // non-centred estimates (i.e. a0 = -5.13 + 0.72*a0_std is the same as a0 = normal(-5.13, 0.72))
  // estimates and standard errors from basta model of Amboseli population
  a0 = -5.13 + 0.72*a0_std;
  a1 = 3.0 + 0.1*a1_std;
  c  = 0.026 + 0.006*c_std;
  b0 = -5.08 + 0.56*b0_std;
  b1 = 0.09 + 0.018*b1_std;
}

model {
  for(i in 1:N) {
    age_category_index[i] ~ ordered_logistic(observed_age[i], thresholds);
  }
  observed_age_std ~ std_normal();
  sigma_age ~ exponential(0.5);
  true_age ~ gompertz_bathtub(a0, a1, c, b0, b1);
  a0_std ~ std_normal();
  a1_std ~ std_normal();
  b0_std ~ std_normal();
  b1_std ~ std_normal();
  c_std ~ std_normal();
}




//functions {
  // 1 - gompertz to give survival to an age
  // a is the asymtote
  // b shifts x axis
  // c + number, 'growth' rate
  // t number of years (first parameter is response)
//  real gompertz_lpdf(vector t, real a, real b, real c) {
//    return sum(log(1-(a * exp(-b*exp(-c*t)))));
//  }
//}
