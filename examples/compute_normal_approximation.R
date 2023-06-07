# import data
data <- readRDS("examples/binary.RData")
df <- data$df
df_agg <- data$df_agg
logit_edge_samples <- data$logit_edge_samples

# calculate edge means and standard deviations
logit_edge_mu <- apply(logit_edge_samples, 2, mean)
logit_edge_cov <- cov(logit_edge_samples)

# compute multivariate Gaussian approximation
edge_samples <- MASS::mvrnorm(1e5, logit_edge_mu, logit_edge_cov)

# plot edge weights against normal approximation for edge 1
par(mfrow=c(1, 2))
plot(density(logit_edge_samples[, 1]), lwd=2, main="Estimated logit edge weight vs normal approximation", xlab="Logit edge weight")
lines(density(edge_samples[, 1]), col=rgb(0, 0, 1, 0.5), lwd=2)

# plot edge 1 against edge 2 to check for covariance
plot(edge_samples[, 1], edge_samples[, 2], col=rgb(0, 0, 1, 0.05), main="Covariance between edges 1 & 2", xlab="Edge 1 samples", ylab="Edge 2 samples")

