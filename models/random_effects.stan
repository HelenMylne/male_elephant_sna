// data {
//   int<lower=0> num_data;                            // num individuals
//   int<lower=1> num_exposures;                       // num ind predictors -- 1
//   int<lower=1> num_windows;                         // num groups
//   int<lower=1> L;                                   // num group predictors
//   int<lower=1,upper=num_windows> window[num_data];  // group for individual
//   matrix[num_data, num_exposures] age;              // individual predictors
//   row_vector[L] b_window[num_windows];              // group predictors
//   vector[num_data] y;                               // outcomes
// }
// parameters {
//   corr_matrix[num_exposures] omega_window;          // prior correlation
//   vector<lower=0>[num_exposures] sigma_window;      // prior scale
//   matrix[L, num_exposures] mu_window;               // group coeffs
//   vector[num_exposures] beta[num_windows];          // indiv coeffs by group
//   real<lower=0> sigma;                              // prediction error scale
// }
// model {
//   sigma_window ~ cauchy(0, 2.5);
//   omega_window ~ lkj_corr(2);
//   to_vector(mu_window) ~ normal(0, 5);
//   {
//     row_vector[num_exposures] b_mu_window[J];
//     for (j in 1:num_windows)
//       b_mu_window[j] = b_window[j] * mu_window;
//     beta ~ multi_normal(b_mu_window, quad_form_diag(omega_window, sigma_window));
//   }
//   for (n in 1:num_data)
//     y[n] ~ normal(age[n] * beta[window[n]], sigma);
// }


data {
  int<lower=0> num_data;                                 // num individuals
  int<lower=1> num_exposures;                            // num ind predictors -- 1
  int<lower=1> num_windows;                              // num groups -- random 1
  int<lower=1> num_dyads;                                // num groups -- random 2
  int<lower=1> num_random_effects;                       // num group predictors -- 2
  int<lower=1,upper=num_windows> window[num_data];       // group for individual
  int<lower=1,upper=num_dyads> dyad_id[num_data];        // group for individual
  matrix[num_data, num_exposures] age;                   // individual predictors
  row_vector[num_random_effects] b_random[num_windows*num_dyads];  // group predictors
  vector[num_data] y;                                    // outcomes
}
parameters {
  corr_matrix[num_exposures] omega_window;               // prior correlation
  vector<lower=0>[num_exposures] sigma_window;           // prior scale
  matrix[num_random_effects, num_exposures] mu_window;   // group coeffs
  vector[num_exposures] beta[num_windows];               // indiv coeffs by group
  real<lower=0> sigma;                                   // prediction error scale
}
model {
  sigma_random ~ cauchy(0, 2.5);
  omega_random ~ lkj_corr(2);
  to_vector(mu_window) ~ normal(0, 5);
  {
    row_vector[num_exposures] b_mu_window[num_windows];
    row_vector[num_exposures] b_mu_dyad[num_dyads];
    for (window in 1:num_windows) {
      for(dyad in 1:num_dyads){
        b_mu_window[window] = b_window[window] * mu_window;
        b_mu_dyad = b_dyad[dyad] * mu_dyad;
      }
    }
    beta ~ multi_normal(b_mu_random, quad_form_diag(omega_random, sigma_random));
  }
  for (n in 1:num_data)
    y[n] ~ normal(age[n] * beta[window[n]], sigma);
}
