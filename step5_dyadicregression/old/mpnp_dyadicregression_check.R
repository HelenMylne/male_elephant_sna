logit_edge_draws_cov <- readRDS('../data_processed/step5_dyadicregression/mpnp_logit_edgedraws_cov.RDS')
print('matrix imported')

(total_size <- nrow(logit_edge_draws_cov) * ncol(logit_edge_draws_cov))
(na_size <- length(which(is.na(logit_edge_draws_cov))))
na_size == total_size
print('check NAs done -- if TRUE, Viking didn't manage any covariance calculations')

edge_samples <- readRDS('../data_processed/step3_edgeweightestimation/mpnplong_edgedistributions_conditionalprior_wideformat.RDS')
print('edges imported')

edge_samples <- as.matrix(edge_samples)
edge_samples <- edge_samples[,1:2]
print('edges cut down to only first two dyads')

cov(edge_samples)
print('covariance calculated for first two dyads only')




