# prior predictive check
priors$fixed
prior_check(priors, 'binary')

age <- 1:60
plot(NULL, xlim = c(0,60), ylim = c(0,1), las = 1,
     xlab = 'age', ylab = 'eigenvector centrality')
for(i in 1:100){
  intercept <- rbeta(1,1,1)
  beta <- rnorm(1, 0, 0.005)
  lines(x = age, y = intercept + age*beta, col = 'blue')
  lines(x = age, y = intercept + age*plogis(beta), col = 'red') # can't ever be a negative effect so that can't be right
}

mean_age <- mean(motnp_ages$age)
plot(NULL, xlim = c(10,60), ylim = c(0,1), las = 1,
     xlab = 'age', ylab = 'eigenvector centrality')
for(i in 1:100){
  intercept <- rbeta(1,2,2)   # this isn't right but I don't think there is a prior for the intercept?? I've gone for a symmetrical one here that in itself explores most of the parameter space and allows it to see whether some of the lines are steep enough to go from top to bottom, but on the assumption that when combined, they will explore only the space relevant to their starting position
  beta <- rnorm(1, 0, 0.005)
  lines(x = age, y = intercept + (age - mean_age)*beta, col = rgb(0,0,1,0.5)) # vast majority come out somewhere sensible, and those that don't would if they started at a different value for age 10 so that comes down to my ability to work out what the intercept prior should actually be -- think this is a good prior for the slope
}

