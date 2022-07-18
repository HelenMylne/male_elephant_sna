library(tidyverse)
# The final common type of network analysis we will cover here is nodal regression, where a regression is performed to analyse the relationship between a nodal network metric(such as centrality) and nodal traits (such as age and sex).  These analyses are usually used to assess how network position depends on various biological factors, but can also be used where network position is a predictor.  Since node metrics are derivative measures of the network, uncertainty in edge weights should ideally propagate through to node metrics, and on to coefficients in regression analyses, giving an accurate estimate of the total uncertainty in inferred parameters.  The core of the nodal regression is similar to dyadic regression, taking the form:
#X(n)i j∼Poisson(λ(n)i jD(n)i j)
#logÄλ(n)i jä=ωi j+L(n)i j
#Mi(ω)∼Normal(β0+β1xi,σ2)
# where β0 is the intercept parameter, β1 is the slope parameter, and σ is the standard deviation of the residuals. Mi(ω) denotes the node metric estimates when applied to the edge weights estimated in the top part of the model. Calculating node metrics within the model may present a practical challenge when using standard model fitting software, as many node metrics can not be described purely in terms of closed-form equations (Freeman, 1978).  In this case, splitting up the model along the dashed line becomes an important step, as it allows the edge weights to be estimated using a standard piece of software, and leaves the regression part of the model to be fit using a sampler that supports numeric estimates of network centralities. As part of the code we have made available, we have written a custom Metropolis-Hastings sampler that allows numeric estimates of node metrics to be used when fitting the model. Importantly, our sampler samples from the joint distribution of edge weights, rather than sampling edge weights independently. This ensures that the effect of any structure in the observation data (such as location effects) is maintained and propagated through to the node metric estimates, and subsequently to regression coefficient estimates.

##### import data #####
cent_mat <- read_csv('data_processed/motnp_bayesian_nodalregression_centralitymatrix_22.03.04.csv')
nodes <- read_csv('data_processed/motnp_elenodes_centrality_22.03.03.csv')
centrality_matrix <- data.matrix(cent_mat)
head(centrality_matrix)

##### including just sex ####
# loglik function from BISoN github, adapted to have 3 categorical levels for sex. Retained 0±2.5 priors to allow variation positive and negative.
loglik <- function(params, Y, X, index) {
  # Define parameters
  intercept <- params[1]
  beta_females <- params[2]
  beta_males <- params[3]
  beta_unks <- params[4]
  sigma <- exp(params[5]) # Exponential keeps underlying value unconstrained, which is much easier for the sampler.
  
  # Sample data according to index
  y <- Y[index %% dim(Y)[1] + 1, ]
  
  # Define model
  target <- 0
  target <- target + sum(dnorm(y, mean = intercept + beta_females*X[,1] + beta_males*X[,2] + beta_unks*X[,3], sd = sigma, log = TRUE))                                 # Main model
  target <- target + dnorm(intercept,    mean = 0, sd = 2.5, log = TRUE)  # Prior on intercept
  target <- target + dnorm(beta_females, mean = 0, sd = 2.5, log = TRUE)  # Prior on female coefficient
  target <- target + dnorm(beta_males,   mean = 0, sd = 2.5, log = TRUE)  # Prior on male coefficient
  target <- target + dnorm(beta_unks,    mean = 0, sd = 2.5, log = TRUE)  # Prior on unknown coefficient
  target <- target + dexp(sigma, 1, log = TRUE)                           # Prior on sigma
  
  return(target)
}

# The predictor matrix is simply a matrix with 3 columns and 463 rows, corresponding to whether each of the nodes is a female (column 1), a male (column 2), or an unknown (column 3). -- Adapted Jordan's code to add an additional column for unknowns.
predictor_matrix <- matrix(0, nrow = length(nodes$id), ncol = 3)
colnames(predictor_matrix) <- c('female', 'male','unknown')
predictor_matrix[which(nodes$sex == 'F'), 1] <- 1
predictor_matrix[which(nodes$sex == 'M'), 2] <- 1
predictor_matrix[which(nodes$sex == 'U'), 3] <- 1
predictor_matrix

# Since network strength is strictly positive, a Gaussian error is not a reasonable model for the data. The Gaussian family model is much easier to implement as well as interpret than many other models, so we will standardise the centralities by taking z-scores. -- UNCHANGED FROM JORDAN'S CODE
centrality_matrix_std <- (centrality_matrix - apply(centrality_matrix, 1, mean))/apply(centrality_matrix, 1, sd)

plot(centrality_matrix_std[,1])
plot(centrality_matrix_std[,10])
plot(centrality_matrix_std[,100])
plot(centrality_matrix_std[,150])
plot(centrality_matrix_std[,200])
plot(centrality_matrix_std[,250])
plot(centrality_matrix_std[,300])
plot(centrality_matrix_std[,350])
plot(centrality_matrix_std[,400])
plot(centrality_matrix_std[,463])

par(mfrow = c(1,1))

plot(density(centrality_matrix_std[,1]), ylim=c(0, 10), xlim = c(-2.5,2.5), main="", xlab="Standardised network strength")
for (i in 1:463) {
  lines(density(centrality_matrix_std[,i]), col = ifelse(nodes$age_class[i] == 'Adult', rgb(0, 0, 1, 0.25),
                                                         ifelse(nodes$age_class[i] == 'Pubescent', rgb(1, 0, 0, 0.25),
                                                                rgb(0,1,0,0.25))))
}

plot(density(centrality_matrix_std[1, ]), ylim=c(0, 0.7), main="", xlab="Standardised network strength")
sample_ids <- sample(1:1000, size=200)
for (i in sample_ids) {
  lines(density(centrality_matrix_std[i, ]), col=rgb(0, 0, 0, 0.25))
}

# Now we’re in a position to fit the model. To do this, we define the target function, which is simply a function that maps candidate parameters and a network centrality index to the log-likelihood of that function for the given sample of the centrality posterior. This means the target function can be written as a function of the data centrality_matrix_std and predictor_matrix. -- UNCHANGED FROM JORDAN'S CODE
target <- function(params, index) loglik(params, centrality_matrix_std, predictor_matrix, index)

# define metropolis sampler function -- UNCHANGED FROM JORDAN'S CODE
metropolis <- function(target, initial, iterations=10000, warmup=2000, thin=4, chain_id=1, refresh=2000) {
  k <- length(initial)
  chain <- matrix(0, iterations + warmup, k)
  ll_obs <- NULL
  pred <- NULL
  chain[1, ] <- initial
  acceptances <- c(1)
  prop_dispersion <- 2.38
  for (i in 2:(iterations + warmup)) {
    current <- chain[i - 1, ]
    # Single Component Adaptive Metropolis (SCAM)
    if (i <= 100) {
      candidate <- rnorm(k, mean=current, sd=5)
    } else {
      prop_var <- prop_dispersion * (current_var + 0.05)
      candidate <- rnorm(k, mean=current, sd=sqrt(prop_var))
    }
    current_lk <- target(current, i)
    candidate_lk <- target(candidate, i)
    A <- exp(candidate_lk - current_lk)
    if (runif(1) < A) { # Accept
      chain[i, ] <- candidate
      acceptances[length(acceptances) + 1] <- 1
    } else {
      acceptances[length(acceptances) + 1] <- 0
      chain[i, ] <- current
    }
    # Update current mean and current var for next iteration.
    if (i == 100) {
      current_mean <- colMeans(chain[1:i, ])
      current_var <- apply(chain[1:i, ], 2, var)
    } else if (i >= 100) {
      prev_mean <- current_mean
      prev_var <- current_var
      current_mean <- (i - 1)/i * prev_mean + (1/i) * chain[i, ]
      current_var <- (i - 2)/(i - 1) * prev_var + prev_mean^2 + 1/(i - 1) * (chain[i, ])^2 - i/(i - 1) * current_mean^2
    }
    # Adjust acceptance rate for batch.
    if (i %% 100 == 0) {
      acc_batch <- mean(acceptances[(i - 99):i])
      if (acc_batch < 0.23) {
        prop_dispersion <- prop_dispersion/exp(sqrt(1/(i/100)))
      } else {
        prop_dispersion <- prop_dispersion * exp(sqrt(1/(i/100)))
      }
    }
    # Print out progress
    if (refresh != 0 && i %% refresh == 0) {
      if (i < warmup) {
        cat(paste0("Chain: ", chain_id, " | Iteration: ", i, "/", warmup + iterations, " (Warmup)\n"))
      } else {
        cat(paste0("Chain: ", chain_id, " | Iteration: ", i, "/", warmup + iterations, " (Sampling)\n"))
      }
    }
  }
  # close(pb)
  cat(paste0("Acceptance Rate: ", mean(acceptances), "\n"))
  return(chain[seq(warmup, iterations + warmup - 1, thin), ])
}

# The function metropolis can now be used to fit the model using the provided target function, an initial set of parameters, and some additional MCMC options.
chain <- metropolis(target, c(0, 0, 0, 0, 0), iterations=100000, thin=100, refresh=10000)

colnames(chain) <- c("intercept", "beta_female", "beta_male", "beta_unknown", "sigma")
head(chain)

par(mfrow=c(3,2))
for (i in 1:5) {
  plot(chain[, i], type="l") # these don't look great apart from the sigma one, and 3 sexes are almost identical
}
par(mfrow=c(1,1))

plot(density(centrality_matrix_std[1, ]), ylim=c(0, 0.7), main="", xlab="Standardised network strength")
sample_ids <- sample(1:1000, size=200)
for (i in sample_ids) {
  pred <- rnorm(463, mean = chain[i,"intercept"] + chain[i,"beta_female"]*predictor_matrix[,1] + chain[i,"beta_male"]*predictor_matrix[,2] + chain[i,"beta_unknown"]*predictor_matrix[,3], sd = exp(chain[i,"sigma"]))
  lines(density(centrality_matrix_std[i, ]), col=rgb(0, 0, 0, 0.25)) # has a double peak but I don't really know why
  lines(density(pred), col=rgb(0, 0, 1, 0.25))                       # single peak between those of centrality_matrix_std
}
