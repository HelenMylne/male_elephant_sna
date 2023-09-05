data {
  int<lower=1> n_nodes;   // total number of observations
  vector[n_nodes] eigen;  // response variable
  int<lower=1> K;         // number of population-level effects
  matrix[n_nodes, K] age; // population-level design matrix
}
transformed data {
  int K_centre = K - 1;
  matrix[n_nodes, K_centre] age_centre; // centered version of age without an intercept
  vector[K_centre] means_age;           // column means of age before centering
  for (i in 2 : K) {                    // centre age values
    means_age[i - 1] = mean(age[ : , i]);
    age_centre[ : , i - 1] = age[ : , i] - means_age[i - 1];
  }
}
parameters {
  vector[K_centre] b;  // population-level effects
  real Intercept;      // temporary intercept for centered predictors
  real<lower=0> sigma; // dispersion parameter
}
transformed parameters {
  real lprior = 0;     // prior contributions to the log posterior
  lprior += student_t_lpdf(Intercept | 3, 0.7, 2.5);
  lprior += student_t_lpdf(sigma | 3, 0, 2.5)
            - 1 * student_t_lccdf(0 | 3, 0, 2.5);
}
model {
  // likelihood including constants
  target += normal_id_glm_lpdf(eigen | age_centre, Intercept, b, sigma);
  // priors including constants
  target += lprior;
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_age, b);
}
