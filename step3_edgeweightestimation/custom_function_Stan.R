library(brms)

# Dummy data
data <- data.frame(y = rbeta(100, 2, 2),
                   x = factor(sample(0:1, 100, replace = TRUE)),
                   z = rnorm(100))

# dan code ####
# Define the prior function
myprior <- custom_prior("
  real conditional_prior_lpdf(real mu, int x) {
    if (x == 1) {
      return beta_lpdf(mu | 2, 2);
    } else {
      return beta_lpdf(mu | 1, 1);
    }
  }
")

# Fit the model
fit <- brm(
  y ~ z * x,
  data = data,
  family = beta(),
  prior = prior(conditional_prior(b | x), class = "b", coef = "zx1"),
  custom_prior = myprior
)

# try again ####
myprior_hkm <- prior_string("
real conditional_prior(real mu, int x) {
    if (x == 1) {
      return normal(mu | -2.5, 1.5);
    } else {
      return normal(mu | -5, 3);
    }
  }
")

# Fit the model
fit <- brm(
  y ~ z*x,
  data = data,
  family = gaussian(),
  #prior = prior(conditional_prior(b | x), class = "b", coef = "x1"),
  prior = prior(normal(ifelse(data$x == 1, -2.5, -5), ifelse(data$x == 1, 1.5, 3)), class = "b", coef = "x1"),
  custom_prior = myprior_hkm
)

make_stancode(y ~ z * x,
              data = data,
              family = gaussian(),
              prior = prior(conditional_prior(b | x), class = "b", coef = "x1"),
              custom_prior = myprior_hkm)
'// generated with brms 2.18.0
functions {
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> sigma;  // dispersion parameter
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += student_t_lpdf(Intercept | 3, 0.5, 2.5);
  lprior += student_t_lpdf(sigma | 3, 0, 2.5)
  - 1 * student_t_lccdf(0 | 3, 0, 2.5);
}
model {
  // likelihood including constants
  if (!prior_only) {
    target += normal_id_glm_lpdf(Y | Xc, Intercept, b, sigma);
  }
  // priors including constants
  target += lprior;
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}'





"data {
int<lower=0> num_rows; // Number of data points
int<lower=0> num_edges; // Number of edge weights to estimate
int<lower=0> num_fixed; // Number of fixed effect parameters
int<lower=0> num_random; // Number of random effect parameters
int<lower=0> num_random_groups; // Number of random effect groups

array[num_rows] int event; // Outcome for each data point (presence/absence)
array[num_rows] int divisor; // Duration of each observation
array[num_rows] int dyad_ids; // Dyad IDs of each observation for indexing edge weights
matrix[num_rows, num_fixed] design_fixed; // Design matrix for fixed effects
matrix[num_rows, num_random] design_random; // Design matrix for random effects.
array[num_random] int<lower=0, upper=num_random_groups> random_group_index; // Index for groupings for random effects

real prior_edge_mu_0; // Prior mean for edge weight when never seen together
real<lower=0> prior_edge_sigma_0; // Prior standard deviation for edge weight when never seen together
real prior_edge_mu_1; // Prior mean for edge weight when seen together at least once
real<lower=0> prior_edge_sigma_1; // Prior standard deviation for edge weight when seen together at least once

real prior_fixed_mu; // Prior mean for fixed effects
real<lower=0> prior_fixed_sigma; // Prior standard deviation for fixed effects
real prior_random_mean_mu; // Prior mean on centralisation of random effects
real<lower=0> prior_random_mean_sigma; // Prior standard deviation on centralisation of random effects
real<lower=0> prior_random_std_sigma; // Prior standard deviation on dispersion of random effects
real<lower=0> prior_zero_prob_alpha; // Prior alpha on zero inflation
real<lower=0> prior_zero_prob_beta; // Prior beta on zero inflation

int<lower=0, upper=1> priors_only; // Whether to sample from only the priors
int<lower=0, upper=1> partial_pooling; // Whether to pool edge weight estimates
int<lower=0, upper=1> zero_inflated; // Whether to use zero-inflated edge model
}

parameters {
vector[num_edges] edge_weight; // Parameters for edge weights.
vector[num_fixed] beta_fixed; // Parameters for fixed effects.
vector[num_random] beta_random; // Parameters for random effects.
vector[num_random_groups] random_group_mu; // Hyperpriors for random effects (mean).
vector<lower=0>[num_random_groups] random_group_sigma; // Hyperpriors for random effects (std. dev.).
vector<lower=0>[partial_pooling] edge_sigma; // Random effect for edge weight pooling.                             
vector<lower=0>[zero_inflated] zero_prob; // Zero inflated parameter for probability of zeroes.                    
}                                                                                                                    

transformed parameters {                                                                                             
vector[num_rows] predictor;                                                                                        
predictor = rep_vector(0, num_rows);                                                                               
if (num_edges > 0) {                                                                                               
  predictor += edge_weight[dyad_ids];                                                                              
}                                                                                                                  
if (num_fixed > 0) {                                                                                               
  predictor += design_fixed * beta_fixed;                                                                          
}                                                                                                                  
if (num_random > 0) {                                                                                              
  predictor += design_random * beta_random;                                                                        
}                                                                                                                  
}                                                                                                                    

model {                                                                                                              
if (!priors_only) {                                                                                                
  // Main model                                                                                                    
  if (zero_inflated == 0) {                                                                                        
    event ~ binomial(divisor, inv_logit(predictor));                                                               
  } else {                                                                                                         
    for (i in 1:num_rows) {                                                                                        
      if (event[i] == 0) {                                                                                         
        target += log_sum_exp(                                                                                     
          bernoulli_lpmf(1 | zero_prob[1]),                                                                        
          bernoulli_lpmf(0 | zero_prob[1]) + binomial_lpmf(event[i] | divisor[i], inv_logit(predictor[i]))         
        );                                                                                                         
      } else {                                                                                                     
        target += bernoulli_lpmf(0 | zero_prob[1]) + binomial_lpmf(event[i] | divisor[i], inv_logit(predictor[i]));
      }                                                                                                            
    }                                                                                                              
    zero_prob ~ beta(prior_zero_prob_alpha, prior_zero_prob_beta);                                                 
  }                                                                                                                
}                                                                                                                  
// Priors
  if (num_edges > 0) {
    for (i in 1:num_edges) {
      log_lik[i] = binomial_lpmf(event[i] | divisor[i], inv_logit(predictor[i]));
      # add for loop to tell it which prior to use -- for (i in 1:num_edges){ if x == 1 use prior1, otherwise use prior2 }
      if (partial_pooling == 0) {
        if (seen_together == 1){
          edge_weight ~ normal(prior_edge_mu_1, prior_edge_sigma_1);
        }
        else { edge_weight ~ normal(prior_edge_mu_0, prior_edge_sigma_0) };
      } else {
        edge_weight ~ normal(prior_edge_mu, edge_sigma[1]);
        edge_sigma ~ normal(0, prior_edge_sigma);
      }
    }
  }
  
  if (num_fixed > 0) {
  beta_fixed ~ normal(prior_fixed_mu, prior_fixed_sigma);
}

  if (num_random > 0) {
  beta_random ~ normal(random_group_mu[random_group_index], random_group_sigma[random_group_index]);
  // Hyperpriors
  random_group_mu ~ normal(prior_random_mean_mu, prior_random_mean_sigma);
  random_group_sigma ~ normal(0, prior_random_std_sigma);
}
}

generated quantities {
array[num_rows] int event_pred;
event_pred = binomial_rng(divisor, inv_logit(predictor));
vector[num_rows] log_lik;
for (i in 1:num_rows) {
   log_lik[i] = binomial_lpmf(event[i] | divisor[i], inv_logit(predictor[i]));
}
}"

