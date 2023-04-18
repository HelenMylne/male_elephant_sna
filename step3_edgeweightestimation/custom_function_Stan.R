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
myprior_hkm <- "
real conditional_prior_lpdf(real mu, int x) {
    if (x == 1) {
      return normal_lpdf(mu | -2.5, 1.5);
    } else {
      return normal_lpdf(mu | -5, 3);
    }
  }
"

# Fit the model
fit <- brm(
  y ~ z * x,
  data = data,
  family = gaussian(),
  prior = prior(conditional_prior(b | x), class = "b", coef = "x1"),
  custom_prior = myprior_hkm
)

