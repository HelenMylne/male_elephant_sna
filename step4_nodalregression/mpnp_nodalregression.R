#### information ####
# Makgadikgadi Pans National Park -- regression to test effect of individual age on eigenvector centrality

#### set up ####
#options(future.globals.maxSize = 10000*(1024^2))   # running model with full age distributions = 5.01GB of globals, which exceeds default maximum allowed size. Set to allow up to 10 GB to allow model to run

# library(rstan) ; library(igraph) ; library(tidyverse) ; library(LaplacesDemon) ; library(MASS)
library(StanHeaders, lib.loc = '../packages/')
library(rstan, lib.loc = '../packages/')
library(sna, lib.loc = '../packages/')             # library(sna)
#library(igraph, lib.loc = '../packages/')          # library(igraph)
library(tidyverse, lib.loc = '../packages/')       # library(tidyverse)
library(LaplacesDemon, lib.loc = '../packages/')   # library(LaplacesDemon)
library(MASS, lib.loc = '../packages/')            # library(MASS)
#library(cmdstanr, lib.loc = '../packages/')        # library(cmdstanr)
#library(brms, lib.loc = '../packages/')            # library(brms)
#library(bisonR, lib.loc = '../packages/')          # library(bisonR)
#library(janitor, lib.loc = '../packages/')         # library(janitor)

## set cmdstan path
#set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0/')

## set seed for reproducibility
set.seed(12345)

## set up pdf
#pdf('../outputs/step4_nodalregression/mpnplong_nodalregression_modelprep.pdf')

## set up theme for plots
theme_set(theme_classic())

# ######## long window ########
# ## load data and remove additional data
# load('mpnp_edgecalculations/mpnplong_edgeweights_conditionalprior.RData')
# rm(eles, edgelist, summary, edge_binary, fit_edges_mpnp1, mean_ages, missing_age, make_edgelist, plot_network_threshold_mpnp) ; gc()
# 
# #### filter down to only elephants with known age categories ####
# ## ages
# hist(nodes$age, breaks = 50)
# nodes$age_cat <- ifelse(nodes$age < 10, 1,
#                         ifelse(nodes$age < 16, 2,
#                                ifelse(nodes$age < 21, 3,
#                                       ifelse(nodes$age < 26, 4,
#                                              ifelse(nodes$age < 36, 5, 6)))))
# ## nodes data frame
# nodes <- nodes %>%
#   filter(! is.na(age_cat))
# n_age_cat <- length(unique(nodes$age_cat))
# ele_ids <- nodes$id
# n_eles <- length(ele_ids)
# 
# ## dyads data frame
# counts_df <- counts_df %>%
#   filter(id_1 %in% ele_ids) %>%
#   filter(id_2 %in% ele_ids)
# n_dyads <- nrow(counts_df)
# 
# #### extract centralities ####
# ## clean up nodes
# ncol(edge_samples) ; nrow(counts_df)
# nodes <- nodes %>%
#   mutate(node_rank = as.integer(as.factor(node)))
# node_join <- nodes %>%
#   dplyr::select(node_rank, node, id) %>%
#   rename(node_1 = node, id_1 = id, node_rank_1 = node_rank) %>%
#   mutate(node_2 = node_1, id_2 = id_1, node_rank_2 = node_rank_1)
# 
# ## add rank IDs to counts_df
# counts_df <- counts_df %>%
#   left_join(node_join[,c('node_rank_1','node_1','id_1')],
#             by = c('node_1', 'id_1')) %>%
#   left_join(node_join[,c('node_rank_2','node_2','id_2')],
#             by = c('node_2', 'id_2'))
# 
# ## build adjacency tensor
# adj_tensor <- array(NA, c(n_eles, n_eles, n_samples*n_chains),
#                     dimnames = list(ele_ids, ele_ids, NULL))
# for(i in 1:n_dyads) {
#   dyad_row <- counts_df[i,]
#   adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- edge_samples[,i]
#   adj_tensor[dyad_row$id_2, dyad_row$id_1, ] <- edge_samples[,i]
# }
# adj_tensor[,,1]
# 
# ## calculate centrality and store posterior samples in a matrix
# centrality_samples_invlogit <- matrix(0, n_chains*n_samples, n_eles,
#                                       dimnames = list(NULL,ele_ids))
# for (draw in 1:(n_chains*n_samples)) {
#   centrality_samples_invlogit[draw, ] <- sna::evcent(dat = adj_tensor[,,draw],
#                                                      gmode = 'graph', diag = F)
# }
# head(centrality_samples_invlogit)      # unstandardised eigenvector centrality
# 
# ## convert to logit scale
# centrality_samples <- logit(centrality_samples_invlogit)
# head(centrality_samples)
# 
# ## save progress
# save.image('mpnp_nodalregression/mpnplong_nodalregression.RData')
# 
# ## visualise centralities
# df_wide <- data.frame(centrality_samples)
# colnames(df_wide) <- 1:n_eles
# df_long <- pivot_longer(df_wide, cols = everything(),
#                         names_to = "node_rank", values_to = "centrality") %>%
#   mutate(node_rank = as.integer(node_rank)) %>%
#   left_join(nodes[,c('node_rank','age')], by = 'node_rank')
# df_long %>%
#   filter(node_rank <= 50) %>%
#   mutate(nodes_reordered = fct_reorder(.f = as.factor(node_rank), .x = age, .desc = T)) %>%
#   ggplot(aes(x = centrality, fill = age)) +
#   geom_density(linewidth = 0.4) +
#   facet_grid(rows = vars(as.factor(nodes_reordered)), scales = "free") +
#   labs(x = "Eigenvector centrality (standardised)") +
#   theme_void() +
#   theme(strip.text.y = element_text(size = 12),
#         axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
#         axis.title.x = element_text(size = 12),
#         plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
# 
# #### compute normal approximation ####
# ## check covariance
# plot(centrality_samples[, 1], centrality_samples[, 2],
#      xlab = 'standardised eigenvectors ID1',
#      ylab = 'standardised eigenvectors ID2',
#      las = 1, pch = 19, col = rgb(0,0,1,0.2))
# plot(centrality_samples[, which(nodes$age_cat == 1)[1]],
#      centrality_samples[, which(nodes$age_cat == 6)[1]],
#      xlab = 'standardised eigenvectors ID1',
#      ylab = 'standardised eigenvectors ID2',
#      las = 1, pch = 19, col = rgb(0,0,1,0.2))
# 
# ## check variance
# par(mfrow = c(5,5), mai = c(0.2,0.2,0.2,0.2))
# to_plot <- sample(1:ncol(centrality_samples), 25, replace = F)
# for(i in 1:length(to_plot)){
#   hist(centrality_samples[,to_plot[i]], main = '')
# }
# 
# ## compute normal approximation
# centrality_mu <- apply(centrality_samples, 2, mean)
# centrality_cov <- cov(centrality_samples)
# centrality_samples_sim <- MASS::mvrnorm(1e5, centrality_mu, centrality_cov)
# for(i in 1:length(to_plot)){
#   plot(density(centrality_samples[, i]), lwd = 2, main = "", xlab = "")
#   lines(density(centrality_samples_sim[, i]), col = rgb(0,0,1,0.5), lwd = 2)
# }
# par(mfrow = c(1,1), mai = c(1,1,1,1))
# rm(g, edge_samples, adj_tensor, i) ; gc()
# 
# save.image('mpnp_nodalregression/mpnplong_nodalregression.RData')
# 
# #### prior predictive check ####
# ## set values
# n <- 100
# n_age_cat <- length(unique(nodes$age_cat))
# beta_age <- rnorm(n, 0, 0.8)
# intercept  <- rnorm(n, logit(0.05), 2)
# age_dirichlet <- rdirichlet(n, c(1,1,1,1,1))
# 
# ## plot
# plot(NULL, las = 1, xlab = 'age category', ylab = 'logit eigenvector',
#      ylim = c(min(centrality_mu)-5, max(centrality_mu)+5), xlim = c(1,n_age_cat))
# abline(h = min(centrality_mu), lty = 2) ; abline(h = max(centrality_mu), lty = 2)
# x <- 1:n_age_cat
# for(i in 1:n){
#   y <- rep(NA, length(x))
#   for(j in 1:length(x)){
#     y[j] <- intercept[i] + beta_age[i]*sum(age_dirichlet[i,][1:x[j]])
#   }
#   lines(x = x, y = y, col = rgb(0,0,1,0.4))
# }
# rm(n, beta_age, intercept, age_dirichlet, sigma, x, y, df_plot, df_wide) ; gc()
# 
# ## save and reset plotting
# save.image('mpnp_nodalregression/mpnplong_nodalregression.RData')
# dev.off()
# pdf('../outputs/step4_nodalregression/mpnplong_nodalregression_modelchecks.pdf')
# 
# #### run model ####
# ## create data list
# eigen_list <- list(num_nodes = n_eles,
#                    num_age_cat = n_age_cat,
#                    length_dirichlet = n_age_cat+1,
#                    centrality_mu = centrality_mu,
#                    centrality_cov = centrality_cov,
#                    node_age = as.integer(nodes$age_cat),
#                    prior_age = rep(1, n_age_cat))
# 
# ## check inputs
# boxplot(centrality_mu ~ nodes$age_cat, notch = T)
# 
# ## load model
# nodal_regression <- stan_model('models/eigen_regression_motnp.stan') # same structure as MOTNP, but now using a student-t distribution instead of a normal because the data are of a slightly different structure
# 
# ## run model
# fit_mpnp_eigen <- sampling(nodal_regression,
#                            data = eigen_list,
#                            cores = n_chains,
#                            chains = n_chains)
# 
# ## save output
# rm(edge_samples, adj_tensor, i, to_plot, draw,summary,centrality_samples_sim, dyad_row) ; gc()
# save.image('mpnp_nodalregression/mpnplong_nodalregression.RData')
# 
# #### check outputs ####
# # load('mpnp_nodalregression/mpnplong_nodalregression.RData')
# ## traceplot linear effect size
# traceplot(fit_mpnp_eigen, pars = c('intercept','beta_age','sigma','nu','predictor[1]','predictor[2]','predictor[3]','predictor[4]','predictor[5]','predictor[6]','predictor[7]','predictor[8]'))
# traceplot(fit_mpnp_eigen, pars = c('delta_j[1]','delta_j[2]','delta_j[3]',
#                                    'delta_j[4]','delta_j[5]','delta_j[6]'))
# 
# ## summarise
# (summary <- as.data.frame(round(summary(fit_mpnp_eigen)$summary, 3)))
# summary$parameter <- rownames(summary)
# par(mfrow = c(2,1))
# hist(summary$Rhat, breaks = 50)
# hist(summary$n_eff, breaks = 50)
# par(mfrow = c(1,1))
# 
# ## extract posterior
# params <- rstan::extract(fit_mpnp_eigen)
# 
# ## plot posterior age effect
# plot(density(params$beta_age),
#      main = "Posterior age effect estimate")
# abline(v = 0, lty = 2)
# 
# ## compare parameters
# # par(mfrow = c(2,2))
# beta_diff_12 <- params$delta[, 1] - params$delta[, 2]
# beta_diff_23 <- params$delta[, 2] - params$delta[, 3]
# beta_diff_34 <- params$delta[, 3] - params$delta[, 4]
# beta_diff_45 <- params$delta[, 4] - params$delta[, 5]
# beta_diff_56 <- params$delta[, 5] - params$delta[, 6]
# # plot(density(beta_diff_12), main="Posterior difference:\n<10 and 10-15") ; abline(v=0, lty=2)
# # plot(density(beta_diff_12), main="Posterior difference:\n10-15 and 15-20") ; abline(v=0, lty=2)
# # plot(density(beta_diff_23), main="Posterior difference:\n15-20 and 21-25") ; abline(v=0, lty=2)
# # plot(density(beta_diff_34), main="Posterior difference:\n21-25 and 26-35") ; abline(v=0, lty=2)
# # plot(density(beta_diff_45), main="Posterior difference:\n26-35 and >35") ; abline(v=0, lty=2)
# # par(mfrow = c(1,1))
# plot(density(beta_diff_12), col = 'green',
#      main="Posterior differences: green=1-2, turquoise=2-3,\nblue=3-4, purple=4-5, red=5-6",
#      ylim = c(0, 5), las = 1, xlim = c(-1,1))
# lines(density(beta_diff_23), col = 'turquoise')
# lines(density(beta_diff_34), col = 'blue')
# lines(density(beta_diff_45), col = 'purple')
# lines(density(beta_diff_56), col = 'red')
# abline(v=0, lty=2)
# 
# #### posterior predictive check ####
# plot(density(centrality_samples[1, ]), main="Posterior predictive density of responses:\nblack = data, blue = predicted", col=rgb(0, 0, 0, 0.25), ylim=c(0, 10))
# for (i in 1:100) {
#   j <- sample(1:(n_chains*n_samples), 1)
#   lines(density(centrality_samples[j, ]), col=rgb(0, 0, 0, 0.25))
#   mu <- rep(NA, n_eles)
#   for(k in 1:n_eles){
#     mu[k] <- params$beta_age[j]*sum(params$delta_j[j, 1:nodes$age_cat[k]]) + params$intercept[j]
#   }
#   sigma <- centrality_cov + diag(rep(params$sigma[j], n_eles))
#   nu <- params$nu[j]
#   lines(density(LaplacesDemon::rmvt(n = 1, mu = mu, S = sigma, df = nu)), col=rgb(0, 0, 1, 0.25))
#   #lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
# }
# 
# ## save and reset plotting
# save.image('mpnp_nodalregression/mpnplong_nodalregression.RData')
# dev.off()
# print('model checks complete')
# pdf('../outputs/step4_nodalregression/mpnplong_nodalregression_modelpredictions.pdf')
# 
# #### predict from model ####
# ## create functions for predictions
# predict_mean_centrality <- function(params, pred_data){
#   mu_matrix <- matrix(data = NA,
#                       nrow = length(params$beta_age),
#                       ncol = nrow(pred_data),
#                       dimnames = list(1:length(params$beta_age),
#                                       pred_data$node))
#   for(i in 1:nrow(mu_matrix)){
#     for(j in 1:ncol(mu_matrix)){
#       mu_matrix[i,j] <- params$intercept[i] + params$beta_age[i] * sum(params$delta_j[i,(1:pred_data$age_cat[j])])
#     }
#   }
#   return(mu_matrix)
# }
# predict_full_centrality <- function(params, mu_matrix, cov_matrix){
#   full_matrix <- mu_matrix
#   for(i in 1:nrow(full_matrix)){
#     full_matrix[i,] <- MASS::mvrnorm(n = 1, mu = mu_matrix[i,],
#                                      Sigma = cov_matrix + diag(rep(params$sigma[i], ncol(full_matrix))))
#   }
#   return(full_matrix)
# }
# 
# ## predicted means
# sim_mean <- predict_mean_centrality(params = params, pred_data = nodes)
# nodes$mu_avg_predict <- apply(sim_mean, 2, mean)
# nodes$mu_lwr_predict <- apply(sim_mean, 2, quantile, probs = 0.025)
# nodes$mu_upr_predict <- apply(sim_mean, 2, quantile, probs = 0.975)
# 
# ## full predictions
# sim_full <- predict_full_centrality(params = params, mu_matrix = sim_mean, cov_matrix = centrality_cov)
# nodes$full_lwr_predict <- apply(sim_full, 2, quantile, probs = 0.025)
# nodes$full_upr_predict <- apply(sim_full, 2, quantile, probs = 0.975)
# 
# ## save output
# save.image('mpnp_nodalregression/mpnplong_nodalregression.RData')
# print('predictions made')
# 
# ## compare to raw data
# nodes <- nodes %>%
#   mutate(mu_raw = centrality_mu,
#          mu_raw_invlogit = LaplacesDemon::invlogit(centrality_mu),
#          age_cat_chr = ifelse(age_cat == 1, '<10',
#                               ifelse(age_cat == 2, '10-15',
#                                      ifelse(age_cat == 3, '16-20',
#                                             ifelse(age_cat == 4, '21-25',
#                                                    ifelse(age_cat == 5, '26-35','36+')))))) %>%
#   relocate(age_cat_chr, .after = age_cat)
# ggplot(nodes)+
#   geom_ribbon(aes(x = age_cat, ymin = full_lwr_predict, ymax = full_upr_predict),
#               alpha = 0.3, fill = 'purple')+
#   geom_ribbon(aes(x = age_cat, ymin = mu_lwr_predict, ymax = mu_upr_predict),
#               alpha = 0.3, fill = 'blue')+
#   geom_jitter(aes(x = age_cat, y = mu_raw), width = 0.2)+
#   geom_line(aes(x = age_cat, y = mu_avg_predict))+
#   theme_classic()
# 
# ## summarise predictions -- all the same for all nodes in the same category because age is the only predictor
# nodes_full <- sim_full %>%
#   as.data.frame() %>%
#   pivot_longer(cols = everything(), names_to = 'node', values_to = 'node_pred') %>%
#   mutate(draw = rep(1:n_eles, each = nrow(sim_full)),
#          node = as.numeric(node),
#          invlogit_pred = LaplacesDemon::invlogit(node_pred)) %>%
#   left_join(nodes[,c('id','node','age_cat_chr','age_cat')],
#             by = 'node') %>%
#   group_by(id) %>%
#   mutate(pred_full_lwr = quantile(node_pred, probs = 0.025),
#          pred_full_mean = mean(node_pred),
#          pred_full_mid = quantile(node_pred, probs = 0.5),
#          pred_full_upr = quantile(node_pred, probs = 0.975),
#          pred_full_lwr_invlogit = quantile(invlogit_pred, probs = 0.025),
#          pred_full_mean_invlogit = mean(invlogit_pred),
#          pred_full_mid_invlogit = quantile(invlogit_pred, probs = 0.5),
#          pred_full_upr_invlogit = quantile(invlogit_pred, probs = 0.975))
# 
# ## plot simulations
# ggplot()+
#   geom_violin(data = nodes_full,
#               aes(x = age_cat_chr, y = node_pred,
#                   fill = factor(age_cat, levels = 1:6)),
#               alpha = 0.5)+
#   geom_boxplot(data = nodes_full,
#                aes(x = age_cat_chr, y = pred_full_mid,
#                    fill = factor(age_cat, levels = 1:6)))+
#   geom_jitter(data = nodes,
#               aes(x = age_cat, y = mu_raw,
#                   fill = factor(age_cat, levels = 1:6),
#                   size = sightings),
#               width = 0.2, pch = 21)+
#   scale_fill_viridis_d()+
#   #scale_colour_viridis_d()+
#   labs(fill = 'age category',
#        #colour = 'age category',
#        x = 'age category',
#        y = 'predicted eigenvector centrality')
# ggsave(file = '../outputs/step4_nodalregression/mpnplong_nodal_violin_logit.png', device = 'png',
#        plot = last_plot(), width = 2100, height = 1600, units = 'px')
# 
# ggplot()+
#   geom_violin(data = nodes_full,
#               aes(x = age_cat_chr, y = invlogit_pred,
#                   fill = factor(age_cat, levels = 1:6)),
#               alpha = 0.5)+
#   geom_boxplot(data = nodes_full,
#                aes(x = age_cat_chr, y = pred_full_mid_invlogit,
#                    fill = factor(age_cat, levels = 1:6)))+
#   geom_jitter(data = nodes,
#               aes(x = age_cat_chr, y = mu_raw_invlogit,
#                   fill = factor(age_cat, levels = 1:6),
#                   size = sightings),
#               width = 0.2, pch = 21)+
#   scale_fill_viridis_d()+
#   #scale_colour_viridis_d()+
#   labs(fill = 'age category',
#        #colour = 'age category',
#        x = 'age category',
#        y = 'predicted eigenvector centrality')
# ggsave(file = '../outputs/step4_nodalregression/mpnplong_nodal_violin.png', device = 'png',
#        plot = last_plot(), width = 2100, height = 1600, units = 'px')
# 
# ## save output
# save.image('mpnp_nodalregression/mpnplong_nodalregression.RData')
# 
# #### extract contrasts from predictions ####
# ## predict for age categroy + 1
# pred_data_new <- nodes[order(nodes$age_cat),] %>%
#   mutate(age_cat_original = age_cat,
#          age_cat = ifelse(age_cat == 6, 1, age_cat + 1))
# pred_mu_new <- predict_mean_centrality(params = params,
#                                        pred_data = pred_data_new)
# pred_full_new <- predict_full_centrality(params = params,
#                                          mu_matrix = pred_mu_new,
#                                          cov_matrix = centrality_cov)
# 
# ## calculate contrast -- all changes
# contrast <- pred_full_new - sim_full
# mean(contrast)
# quantile(contrast, prob = c(0.025, 0.975))
# 
# ## calculate contrast -- 1 vs 2
# contrast12 <- contrast[,which(nodes$age_cat == 1)]
# mean(contrast12)
# quantile(contrast12, prob = c(0.025, 0.975))
# 
# ## calculate contrast -- 2 vs 3
# contrast23 <- contrast[,which(nodes$age_cat == 2)]
# mean(contrast23)
# quantile(contrast23, prob = c(0.025, 0.975))
# 
# ## calculate contrast -- 3 vs 4
# contrast34 <- contrast[,which(nodes$age_cat == 3)]
# mean(contrast34)
# quantile(contrast34, prob = c(0.025, 0.975))
# 
# ## calculate contrast -- 4 vs 5
# contrast45 <- contrast[,which(nodes$age_cat == 4)]
# mean(contrast45)
# quantile(contrast45, prob = c(0.025, 0.975))
# 
# ## calculate contrast -- 5 vs 6
# contrast56 <- contrast[,which(nodes$age_cat == 5)]
# mean(contrast56)
# quantile(contrast56, prob = c(0.025, 0.975))
# 
# ## calculate contrast -- 6 vs 1
# contrast61 <- contrast[,which(nodes$age_cat == 6)]
# mean(contrast61)
# quantile(contrast61, prob = c(0.026, 0.976))
# 
# ## plot contrasts
# contrasts <- data.frame(contrast = c('1 vs 2', '2 vs 3', '3 vs 4', '4 vs 5', '5 vs 6'),
#                         mean = c(mean(contrast12),mean(contrast23),
#                                  mean(contrast34),mean(contrast45),
#                                  mean(contrast56)),
#                         median = c(median(contrast12),median(contrast23),
#                                    median(contrast34),median(contrast45),
#                                    median(contrast56)),
#                         upr = c(quantile(contrast12, prob = c(0.975)),
#                                 quantile(contrast23, prob = c(0.975)),
#                                 quantile(contrast34, prob = c(0.975)),
#                                 quantile(contrast45, prob = c(0.975)),
#                                 quantile(contrast56, prob = c(0.975))),
#                         lwr = c(quantile(contrast12, prob = c(0.025)),
#                                 quantile(contrast23, prob = c(0.025)),
#                                 quantile(contrast34, prob = c(0.025)),
#                                 quantile(contrast45, prob = c(0.025)),
#                                 quantile(contrast56, prob = c(0.025)))) %>%
#   pivot_longer(cols = c('mean', 'upr', 'lwr'),
#                names_to = 'stat', values_to = 'value')
# ggplot()+
#   geom_line(data = contrasts[contrasts$stat != 'mean',],
#             aes(y = value, x = contrast))+
#   geom_hline(aes(yintercept = 0), linetype = 2)+
#   geom_point(data = contrasts[contrasts$stat == 'mean',],
#              aes(y = value, x = contrast))+
#   scale_x_discrete(#expand = c(0.05,0.05),
#     name = 'age category comparison')+
#   scale_y_continuous(name = 'difference in prediction',
#                      limits = c(-1,1))+
#   theme_bw()+
#   coord_flip()
# 
# ## save
# save.image('mpnp_nodalregression/mpnplong_nodalregression.RData')
# 
# #### plot nicely ####
# load('mpnp_nodalregression/mpnplong_nodalregression.RData')
# pdf('../outputs/step4_nodalregression/mpnplong_nodalregression_niceplots.pdf')
# 
# ## clean data frame
# clean <- data.frame(age_cat = 1:6,
#                     mean = NA,
#                     predict_mu_lwr = NA, predict_mu_upr = NA,
#                     predict_full_lwr = NA, predict_full_upr = NA,
#                     mean_invlogit = NA,
#                     predict_mu_lwr_invlogit = NA, predict_mu_upr_invlogit = NA,
#                     predict_full_lwr_invlogit = NA, predict_full_upr_invlogit = NA)
# 
# ## predict means only
# predict_clean_mean <- predict_mean_centrality(params, clean)
# clean$mean <- mean(predict_clean_mean)
# clean$predict_mu_lwr <- apply(predict_clean_mean, 2, quantile, prob = 0.025)
# clean$predict_mu_upr <- apply(predict_clean_mean, 2, quantile, prob = 0.975)
# 
# ## convert to invlogit scale
# predict_clean_mean_invlogit <- LaplacesDemon::invlogit(predict_clean_mean)
# clean$mean_invlogit <- mean(predict_clean_mean_invlogit)
# clean$predict_mu_lwr_invlogit <- apply(predict_clean_mean_invlogit, 2, quantile, prob = 0.025)
# clean$predict_mu_upr_invlogit <- apply(predict_clean_mean_invlogit, 2, quantile, prob = 0.975)
# 
# ## average fulls for whole data frame (can't just predict from clean data frame because dimensions don't fit with covariance matrix)
# for(i in 1:n_age_cat){
#   sim_age <- sim_full[,which(nodes$age_cat == i)]
#   sim_age_invlogit <- LaplacesDemon::invlogit(sim_age)
#   clean$predict_full_lwr[clean$age_cat == i] <- mean(apply(sim_age, 2,
#                                                                quantile, prob = 0.025))
#   clean$predict_full_upr[clean$age_cat == i] <- mean(apply(sim_age, 2,
#                                                                quantile, prob = 0.975))
#   clean$predict_full_lwr_invlogit[clean$age_cat == i] <- mean(apply(sim_age_invlogit, 2,
#                                                                         quantile, prob = 0.025))
#   clean$predict_full_upr_invlogit[clean$age_cat == i] <- mean(apply(sim_age_invlogit, 2,
#                                                                         quantile, prob = 0.975))
# }
# 
# ## add node data
# df_long <- df_long %>%
#   mutate(centrality_invlogit = LaplacesDemon::invlogit(centrality)) %>% 
#   left_join(nodes, by = 'node_rank')
# 
# ## plot on logit scale
# (clean_plot <- ggplot()+
#     geom_ribbon(data = clean,
#                 aes(x = as.numeric(age_cat),
#                     ymin = predict_full_lwr, ymax = predict_full_upr),    # shade simulations
#                 colour = 'transparent', fill = rgb(0,0,0,0.1))+
#     geom_ribbon(data = clean,
#                 aes(x = as.numeric(age_cat),
#                     ymin = predict_mu_lwr, ymax = predict_mu_upr),  # shade mean distribution
#                 colour = 'transparent',
#                 fill = rgb(33/255, 145/255, 140/255, 0.5))+
#     geom_line(data = clean,
#               aes(x = as.numeric(age_cat),
#                   y = mean),                            # mean line
#               colour = rgb(33/255, 145/255, 140/255),
#               linewidth = 1)+
#     scale_x_continuous('age category')+
#     scale_y_continuous('eigenvector centrality')+
#     theme(axis.text = element_text(size = 18),
#           axis.title = element_text(size = 22),
#           legend.text = element_text(size = 18),
#           legend.title = element_text(size = 22)))
# 
# clean_plot +
#   geom_point(data = nodes,
#              aes(x = as.numeric(age_cat),
#                  y = mu_raw,
#                  size = sightings),
#              colour = rgb(68/255, 1/255, 84/255))
# ggsave(filename = '../outputs/step4_nodalregression/mpnplong_nodalregression_line_meanpoints_logit.png', device = 'png',
#        plot = last_plot(), width = 2800, height = 1600, units = 'px')
# 
# clean_plot +
#   geom_point(data = df_long,
#              aes(x = as.numeric(age_cat),
#                  y = centrality),
#              colour = rgb(253/255, 231/255, 37/255, 0.01)) +
#   geom_point(data = nodes,
#              aes(x = as.numeric(age_cat),
#                  y = mu_raw,
#                  size = sightings),
#              colour = rgb(68/255, 1/255, 84/255))
# ggsave(filename = '../outputs/step4_nodalregression/mpnplong_nodalregression_line_allpoints_logit.png', device = 'png',
#        plot = last_plot(), width = 2800, height = 1600, units = 'px')
# 
# ## plot on invlogit scale
# (clean_plot <- ggplot()+
#     geom_ribbon(data = clean,
#                 aes(x = as.numeric(age_cat),
#                     ymin = predict_full_lwr_invlogit,
#                     ymax = predict_full_upr_invlogit),           # shade simulations
#                 colour = 'transparent', fill = rgb(0,0,0,0.1))+
#     geom_ribbon(data = clean,
#                 aes(x = as.numeric(age_cat),
#                     ymin = predict_mu_lwr_invlogit,
#                     ymax = predict_mu_upr_invlogit),             # shade mean distribution
#                 colour = 'transparent',
#                 fill = rgb(33/255, 145/255, 140/255, 0.5))+
#     geom_line(data = clean,
#               aes(x = as.numeric(age_cat),
#                   y = mean_invlogit),                            # mean line
#               colour = rgb(33/255, 145/255, 140/255),
#               linewidth = 1)+
#     scale_x_continuous('age category')+
#     scale_y_continuous('eigenvector centrality')+
#     theme(axis.text = element_text(size = 18),
#           axis.title = element_text(size = 22),
#           legend.text = element_text(size = 18),
#           legend.title = element_text(size = 22)))
# 
# clean_plot +
#   geom_point(data = nodes,
#              aes(x = as.numeric(age_cat),
#                  y = mu_raw_invlogit,
#                  size = sightings),
#              colour = rgb(68/255, 1/255, 84/255))
# ggsave(filename = '../outputs/step4_nodalregression/mpnplong_nodalregression_line_meanpoints_invlogit.png',
#        device = 'png', plot = last_plot(), width = 2800, height = 1600, units = 'px')
# 
# clean_plot +
#   geom_point(data = df_long,
#              aes(x = as.numeric(age_cat),
#                  y = centrality_invlogit),
#              colour = rgb(253/255, 231/255, 37/255, 0.01)) +
#   geom_point(data = nodes,
#              aes(x = as.numeric(age_cat),
#                  y = mu_raw_invlogit,
#                  size = sightings),
#              colour = rgb(68/255, 1/255, 84/255))
# ggsave(filename = '../outputs/step4_nodalregression/mpnplong_nodalregression_line_allpoints_invlogit.png', device = 'png',
#        plot = last_plot(), width = 2800, height = 1600, units = 'px')
# 
# ## save
# save.image('mpnp_nodalregression/mpnplong_nodalregression.RData')
# dev.off()
# 
######## short window ########
rm(list = ls()) ; gc()
pdf('../outputs/step4_nodalregression/mpnpshort_nodalregression_modelprep.pdf')

## set seed for reproducibility
set.seed(12345)

## load data and remove additional data
load('mpnp_edgecalculations/mpnpshort5_edgeweights_conditionalprior.RData')
rm(eles, edges, edgelist, summary, edge_binary, fit_edges_mpnp, mean_ages, make_edgelist, plot_network_threshold_mpnp) ; gc()

nodes_all <- nodes %>% 
  mutate(window = 5)
counts_df_all <- counts_df
edge_samples_all <- list() ; edge_samples_all[[5]] <- edge_samples
num_nodes_all <- list() ; num_nodes_all[[5]] <- n_dyads

for(time_window in 4:1){
  # load next window
  load(paste0('mpnp_edgecalculations/mpnpshort',time_window,'_edgeweights_conditionalprior.RData'))
  rm(eles, edges, edgelist, summary, edge_binary, fit_edges_mpnp, make_edgelist, plot_network_threshold_mpnp, counts_ls) ; gc()
  
  # find age data for ones that don't have it already
  if(length(colnames(nodes_all)) != length(colnames(nodes))+1){
    ages <- readRDS(paste0('../data_processed/step2_ageestimation/mpnp',time_window,'_ageestimates_mcmcoutput.rds'))
    
    nodes <- nodes %>% 
      filter(id %in% colnames(ages)) %>% 
      mutate(age = NA)
    for(i in 1:nrow(nodes)){
      nodes$age[i] <- mean(ages[,which(colnames(ages) == nodes$id[i])])
    }
  }
  
  # join data
  nodes <- nodes %>% 
    mutate(window = time_window)
  nodes_all <- rbind(nodes, nodes_all)
  counts_df_all <- rbind(counts_df, counts_df_all)
  edge_samples_all[[time_window]] <- edge_samples
  num_nodes_all[[time_window]] <- n_dyads
}

cdf_all <- counts_df_all

rm(counts_df_all, edge_samples, nodes, n_dyads, time_window) ; gc()

#### filter down to only elephants with known age categories ####
## ages
hist(nodes_all$age, breaks = 50)
nodes_all$age_cat <- ifelse(nodes_all$age < 16, 1,
                        ifelse(nodes_all$age < 21, 2,
                               ifelse(nodes_all$age < 26, 3,
                                      ifelse(nodes_all$age < 36, 4, 5))))
## nodes data frame
nodes_all <- nodes_all %>%
  filter(! is.na(age_cat))
n_age_cat <- length(unique(nodes_all$age_cat))

## randomise node IDs
node_random <- nodes_all %>% 
  dplyr::select(id, node) %>% 
  distinct()
node_random$node_rand <- sample(1:nrow(node_random), replace = F)
nodes_all <- nodes_all %>% 
  left_join(node_random, by = c('id','node')) %>% 
  mutate(node_window = paste0(node_rand,'_',window))

## dyads data frame
ele_ids <- node_random$id
counts_df <- counts_df %>%
  filter(id_1 %in% ele_ids) %>%
  filter(id_2 %in% ele_ids)

## define population parameters
n_eles <- length(ele_ids)
n_dyads <- nrow(counts_df)

## see if ages have similar average across time windows -- not especially, but not horrendous
ggplot(nodes_all)+
  geom_bar(aes(x = as.factor(age_cat)))+
  facet_wrap(. ~ window)

print('population information collated')

#### extract centralities ####
## create function to extract centrality
extract_eigen_centrality <- function(nodes_df, dyads_df, edgeweight_matrix, logit = TRUE, window){
  ## calculate data size parameters
  num_nodes <- nrow(nodes_df)
  num_dyads <- nrow(dyads_df)
  num_samples <- nrow(edgeweight_matrix)
  
  ## build adjacency tensor
  dyads_df$node_1_id <- as.integer(as.factor(dyads_df$node_1_randomised))
  dyads_df$node_2_id <- as.integer(as.factor(dyads_df$node_2_randomised))+1
  adj_tensor <- array(0, c(num_samples, num_nodes, num_nodes),
                      dimnames = list(NULL, nodes_df$node_window, nodes_df$node_window))
  
  ## fill adjacency tensor
  for(dyad_id in 1:num_dyads) {
    dyad_row <- dyads_df[dyad_id, ]
    adj_tensor[, dyad_row$node_1_id, dyad_row$node_2_id] <- edgeweight_matrix[, dyad_id]
    adj_tensor[, dyad_row$node_2_id, dyad_row$node_1_id] <- edgeweight_matrix[, dyad_id]
  }
  
  ## calculate centrality and store posterior samples in a matrix
  centrality_samples_invlogit <- matrix(0, num_samples, num_nodes,
                                        dimnames = list(NULL, nodes_df$node_window))
  for(i in 1:(num_samples)){
    centrality_samples_invlogit[i, ] <- sna::evcent(adj_tensor[i,,], gmode = 'graph')
  }
  
  ## convert to logit scale
  if(logit == TRUE) {
    centrality_samples <- logit(centrality_samples_invlogit)
    return(centrality_samples)
  } else {
    return(centrality_samples_invlogit)
  }
}

## select nodes data frame for window 1 -- do this from nodes_all, not the nodes data frame already in the saved workspace, so now working from randomised node ID, not original ID
nodes_join <- node_random %>%
  #filter(window == 5) %>%
  mutate(node_1 = node,
         node_2 = node,
         node_1_randomised = node_rand,
         node_2_randomised = node_rand) %>% 
  distinct()

## randomise node IDs in dyads data frame
cdf_all <- cdf_all %>%
  left_join(nodes_join[,c('node_1','node_1_randomised')], by = 'node_1') %>%
  left_join(nodes_join[,c('node_2','node_2_randomised')], by = 'node_2')

## extract centrality for window 5
cents_all <- extract_eigen_centrality(nodes_df = nodes_all[nodes_all$window == 5,],
                                      dyads_df = cdf_all[cdf_all$period == 5,],
                                      edgeweight_matrix = edge_samples_all[[5]],
                                      window = 5)

## extract mean and covariance
cent_mu5 <- apply(cents_all, 2, mean)
cent_cov5 <- cov(cents_all)

## add mean estimate per ID to nodes data frame
mean_df <- as.data.frame(cent_mu5) %>%
  rename(mean_eigen = cent_mu5)
mean_df$node_window <- rownames(mean_df)
nodes_all <- nodes_all %>%
  left_join(mean_df, by = 'node_window')

## prep for combining everything to a single data frame
covs_all <- list()
covs_all[[5]] <- cent_cov5

## clean up
rm(cent_mu5, cent_cov5, mean_df) ; gc()

## for loop
for(time_window in (n_windows-1):1){
  ## import workspace image
  load(paste0('mpnp_edgecalculations/mpnpshort',time_window,'_edgeweights_conditionalprior.RData'))
  
  if('edge_weights_matrix' %in% ls()) {
    edge_samples <- edge_weights_matrix
  }
  
  rm(list = ls()[! ls() %in% c('counts_df','cdf_all','cents_all','covs_all','edge_samples','edge_samples_all','ele_ids','extract_eigen_centrality','n_age_cat','n_chains','n_dyads','n_eles','n_samples','n_windows','node_random','nodes_all','nodes_join','num_nodes_all','time_window')]) ; gc()
  
  ## select randomised nodes data frame for window
  nodes <- nodes_all %>%
    filter(window == time_window) %>%
    mutate(node_1 = node, node_2 = node,
           node_1_randomised = node_rand, node_2_randomised = node_rand) %>%
    dplyr::select(-mean_eigen)
  nodes_all <- nodes_all %>%
    filter(window != time_window)
  
  ## randomise node IDs in dyads data frame
  cdf <- counts_df %>%
    left_join(nodes[,c('node_1','node_1_randomised')], by = 'node_1') %>%
    left_join(nodes[,c('node_2','node_2_randomised')], by = 'node_2')
  cdf <- cdf %>% 
    filter(!is.na(node_1_randomised)) %>% 
    filter(!is.na(node_2_randomised))
  
  ## extract centrality
  cent_new <- extract_eigen_centrality(nodes_df = nodes,
                                       dyads_df = cdf,
                                       edgeweight_matrix = edge_samples,
                                       window = time_window)
  
  ## extract mean and covariance
  cent_mu <- apply(cent_new, 2, mean)
  cent_cov <- cov(cent_new)
  
  ## add mean estimate per ID to nodes data frame
  mean_df <- as.data.frame(cent_mu) %>%
    rename(mean_eigen = cent_mu)
  mean_df$node_window <- rownames(mean_df)
  nodes <- nodes %>%
    dplyr::select(-node_1, -node_2, -node_1_randomised, -node_2_randomised) %>%
    left_join(mean_df, by = 'node_window')
  
  ## combine everything to a single data frame
  nodes_all <- rbind(nodes, nodes_all)
  cents_all <- cbind(cent_new, cents_all)
  covs_all[[time_window]] <- cent_cov
  
  ## add progress marker
  print(time_window)
}

## save workspace
save.image('mpnp_nodalregression/mpnpshort_nodalregression.RData')

## visualise centralities
df_wide <- data.frame(cents_all)
#colnames(df_wide) <- 1:n_eles # THIS IS WRONG
df_long <- pivot_longer(df_wide, cols = everything(),
                        names_to = "x_node_window", values_to = "centrality") %>%
  separate(x_node_window, into = c('x','node_window'), sep = 1) %>% 
  dplyr::select(-x) %>% 
  #mutate(node_rank = as.integer(node_rank)) %>%
  left_join(nodes_all[,c('node_window','node_rand','age_cat')], by = 'node_window')
df_long %>%
  filter(node_rand <= 50) %>%
  #mutate(nodes_reordered = fct_reorder(.f = as.factor(node_rand), .x = age, .desc = T)) %>%
  ggplot(aes(x = centrality, fill = as.factor(age_cat), group = node_window)) +
  geom_density(linewidth = 0.4) +
  scale_fill_viridis_d()+
  facet_grid(rows = vars(as.factor(node_rand)), scales = "free") +
  labs(x = "Eigenvector centrality (standardised)") +
  theme_void() +
  theme(strip.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

## plot full unstandardised distribution
plot(cents_all[1,which(nodes_all$window == 1)] ~ nodes_all$age[which(nodes_all$window == 1)],
     col = 'red', pch = 19, ylim = c(-8,0))      # plot simulated values against age for window 1 (first row of simulated centralities only)
points(cents_all[1,which(nodes_all$window == 2)] ~ nodes_all$age[which(nodes_all$window == 2)],
       col = 'blue', pch = 19)                     # plot simulated values against age for window 2 (first row of simulated centralities only)
points(cents_all[1,which(nodes_all$window == 3)] ~ nodes_all$age[which(nodes_all$window == 3)],
       col = 'green', pch = 19)                    # plot simulated values against age for window 3 (first row of simulated centralities only)
points(cents_all[1,which(nodes_all$window == 4)] ~ nodes_all$age[which(nodes_all$window == 4)],
       col = 'purple', pch = 19)                   # plot simulated values against age for window 4 (first row of simulated centralities only)
points(cents_all[1,which(nodes_all$window == 5)] ~ nodes_all$age[which(nodes_all$window == 5)],
       col = 'yellow', pch = 19)                   # plot simulated values against age for window 5 (first row of simulated centralities only)
print('centralities calculated and plotted')

rm(cdf, cent_cov, cent_new, counts_df, df_long, df_wide, edge_samples, mean_df, node_random, nodes, nodes_join) ; gc()

#### compute normal approximation ####
# compute normal approximation by window -- calculate means per node
sim_cent_mu <- list()
for(time_window in 1:n_windows){
  sim_cent_mu[[time_window]] <- apply(cents_all[,which(nodes_all$window == time_window)], 2, mean)
}

# compute normal approximation by window -- calculate covariance matrix
sim_cent_cov <- list()
for(time_window in 1:n_windows){
  sim_cent_cov[[time_window]] <- cov(cents_all[,which(nodes_all$window == time_window)])
}

## check normal approximation -- simulate from combined mean and covariance, plot curve against relative node ID
par(mfrow = c(round(n_windows/2, 0),2))
for(time_window in 1:n_windows){
  sim_cent_samples <- MASS::mvrnorm(1e5, sim_cent_mu[[time_window]], sim_cent_cov[[time_window]])     # simulate from multivariate normal
  node_id_sample <- nodes_all$node_window[sample(which(nodes_all$window == time_window),1)]
  plot(density(cents_all[,node_id_sample]), lwd = 2, las = 1,                         # plot true density curve
       main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
  lines(density(sim_cent_samples[, node_id_sample]), col = rgb(0,0,1,0.5), lwd = 2)      # overlay normal approximation
}
par(mfrow = c(1,1))
rm(sim_cent_samples, node_id_sample) ; gc()

colours <- RColorBrewer::brewer.pal(n = 11, name = "Spectral")

for(time_window in 1:n_windows){
  node_id_sample <- nodes_all$node_window[sample(which(nodes_all$window == time_window),1)]
  plot(density(cents_all[,node_id_sample]), lwd = 2, las = 1, ylim = c(0,2), xlim = c(-5,0),
       main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
  for(nu in 1:11){
    sim_cent_samples <- LaplacesDemon::rmvt(1e5, mu = sim_cent_mu[[time_window]], 
                                            S = sim_cent_cov[[time_window]], df = seq(2,5,length.out=11)[nu])
    lines(density(sim_cent_samples[, node_id_sample]), 
          col = colours[nu], lwd = 2)      # overlay normal approximation
  }
}
par(mfrow = c(1,1))
rm(sim_cent_samples, node_id_sample) ; gc()


print('normal approximation completed')

#### prior predictive check ####
n <- 100
beta_age <- rnorm(n, 0, 0.8)
intercept  <- rnorm(n, logit(0.05), 2)
age_dirichlet <- rdirichlet(n, rep(1, n_age_cat))
cent_min <- min(c(sim_cent_mu[[1]], sim_cent_mu[[2]], sim_cent_mu[[3]],
                  sim_cent_mu[[4]], sim_cent_mu[[5]]))
cent_max <- max(c(sim_cent_mu[[1]], sim_cent_mu[[2]], sim_cent_mu[[3]],
                  sim_cent_mu[[4]], sim_cent_mu[[5]]))
plot(NULL, las = 1, xlab = 'age category', ylab = 'logit(eigenvector)',
     ylim = c(cent_min-2, cent_max+2),
     xlim = c(1, n_age_cat))
abline(h = cent_min, lty = 2)
abline(h = cent_max, lty = 2)
x <- min(nodes_all$age_cat):max(nodes_all$age_cat)
for(i in 1:n){
  y <- rep(NA, length(x))
  for(j in 1:length(x)){
    y[j] <- intercept[i] + beta_age[i]*sum(age_dirichlet[i,][1:x[j]])
  }
  lines(x = x, y = y, col = rgb(0,0,1,0.4))
}
rm(n, beta_age, intercept, age_dirichlet, x, y) ; gc()

print('prior predictive check completed')

## save output
dev.off()
pdf('../outputs/step4_nodalregression/mpnpshort_nodalregression_modelchecks.pdf')

#### run model ####
## create data
eigen_list <- list(
  # global data size
  num_data = nrow(nodes_all),
  num_nodes = length(unique(nodes_all$id)),
  num_windows = n_windows,
  num_age_cat = n_age_cat,
  length_dirichlet = n_age_cat + 1,
  # prior for age effect
  prior_age = rep(1,n_age_cat),
  # per time window data size
  num_nodes_window1 = num_nodes_all[[1]],
  num_nodes_window2 = num_nodes_all[[2]],
  num_nodes_window3 = num_nodes_all[[3]],
  num_nodes_window4 = num_nodes_all[[4]],
  num_nodes_window5 = num_nodes_all[[5]],
  # number of nodes in all preceding time windows for node age indexing
  num_nodes_prev_windows = c(0, num_nodes_all[[1]], num_nodes_all[[2]],
                             num_nodes_all[[3]], num_nodes_all[[4]]),
  # centrality means per time window
  centrality_mu_1 = sim_cent_mu[[1]],
  centrality_mu_2 = sim_cent_mu[[2]],
  centrality_mu_3 = sim_cent_mu[[3]],
  centrality_mu_4 = sim_cent_mu[[4]],
  centrality_mu_5 = sim_cent_mu[[5]],
  # covariance matrix per time window
  centrality_cov_1 = sim_cent_cov[[1]],
  centrality_cov_2 = sim_cent_cov[[2]],
  centrality_cov_3 = sim_cent_cov[[3]],
  centrality_cov_4 = sim_cent_cov[[4]],
  centrality_cov_5 = sim_cent_cov[[5]],
  # node IDs for all time windows
  nodes_window1 = nodes_all$node_rand[nodes_all$window == 1],
  nodes_window2 = nodes_all$node_rand[nodes_all$window == 2],
  nodes_window3 = nodes_all$node_rand[nodes_all$window == 3],
  nodes_window4 = nodes_all$node_rand[nodes_all$window == 4],
  nodes_window5 = nodes_all$node_rand[nodes_all$window == 5],
  # exposure variable
  node_age = nodes_all$age_cat)

## check inputs
plot(eigen_list$centrality_mu_1 ~ jitter(nodes_all$age_cat[nodes_all$window == 1]),
     pch = 19, col = 'red', las = 1, ylim = -8, 0)
points(eigen_list$centrality_mu_2 ~ jitter(nodes_all$age_cat[nodes_all$window == 2]),
       pch = 19, col = 'blue')
points(eigen_list$centrality_mu_3 ~ jitter(nodes_all$age_cat[nodes_all$window == 3]),
       pch = 19, col = 'green')
points(eigen_list$centrality_mu_4 ~ jitter(nodes_all$age_cat[nodes_all$window == 4]),
       pch = 19, col = 'purple')
points(eigen_list$centrality_mu_5 ~ jitter(nodes_all$age_cat[nodes_all$window == 5]),
       pch = 19, col = 'yellow')

print('model data list created')

## load model
nodal_regression <- cmdstan_model('models/eigen_regression_mpnp.stan')

## run model
n_chains <- 4
n_samples <- 1000
mpnp_short_fit <- sampling(nodal_regression, data = eigen_list, chains = n_chains, cores = n_chains)
# mpnp_short_fit <- nodal_regression$sample(data = eigen_list,
#                                           chains = n_chains, parallel_chains = n_chains,
#                                           iter_warmup = n_samples, iter_sampling = n_samples)
save.image('mpnp_nodalregression/mpnpshort_nodalregression.RData')

#### check outputs ####
## extract model fit
mpnp_short_fit$summary()
summary <- mpnp_short_fit$summary()
par(mfrow = c(3,1))
hist(summary$rhat, breaks = 50)
hist(summary$ess_bulk, breaks = 50)
hist(summary$ess_tail, breaks = 50)
par(mfrow = c(1,1))

## summarise
(summary <- as.data.frame(round(summary(mpnp_short_fit)$summary, 3)))
# summary <- mpnp_short_fit$summary()
summary$parameter <- rownames(summary)
par(mfrow = c(2,1))
hist(summary$Rhat, breaks = 50)
hist(summary$n_eff, breaks = 50)
# par(mfrow = c(3,1))
# hist(summary$rhat, breaks = 50)
# hist(summary$ess_bulk, breaks = 50)
# hist(summary$ess_tail, breaks = 50)
par(mfrow = c(1,1))

## extract posterior
params <- rstan::extract(mpnp_short_fit)
#params <- mpnp_short_fit$draws(format = 'draws_df')

## extract delta and delta_j parameters
delta <- params %>%
  select(`delta[1]`,`delta[2]`,`delta[3]`,`delta[4]`,`delta[5]`,
         `.chain`,`.iteration`,`.draw`)
delta_j <- params %>%
  select(`delta_j[1]`,`delta_j[2]`,`delta_j[3]`,`delta_j[4]`,`delta_j[5]`,`delta_j[6]`,
         `.chain`,`.iteration`,`.draw`)

## separate random effects from global parameters
rand_window <- params %>%
  dplyr::select(grep('rand_window', colnames(params), value=TRUE))

## traceplot all parameters
#traceplot(mpnp_short_fit, pars = c('intercept','beta_age','sigma','predictor[1]','predictor[50]','predictor[100]'))
plot_params <- c('intercept','beta_age','sigma','nu',
                 'rand_node[1]','rand_node[50]','rand_node[100]',
                 'rand_window[1]','rand_window[2]','rand_window[3]',
                 'predictor_window1[1]','predictor_window1[25]','predictor_window1[50]',
                 'predictor_window2[1]','predictor_window2[25]','predictor_window2[50]')
params %>%
  select(all_of(plot_params),`.draw`,`.chain`,`.iteration`) %>%
  pivot_longer(cols = all_of(plot_params), names_to = 'parameter', values_to = 'draw') %>%
  rename(chain = .chain,
         chain_position = .iteration,
         draw_id = .draw) %>%
  #filter(chain == 1) %>% # inspect individual chains -- check for wandery sections that might be hidden by other chains when all plotted together
  ggplot(aes(x = chain_position, y = draw, colour = as.factor(chain)))+
  geom_line()+
  facet_wrap(. ~ parameter, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'none')
rm(plot_params) ; gc()
delta %>%
  pivot_longer(cols = c(`delta[1]`,`delta[2]`,`delta[3]`,`delta[4]`,`delta[5]`),
               names_to = 'parameter', values_to = 'value') %>%
  rename(chain_position = .iteration,
         chain = .chain,
         draw = .draw) %>%
  #filter(chain == 4) %>%
  ggplot(aes(x = chain_position, y = value, colour = as.factor(chain)))+
  geom_line()+
  facet_wrap(. ~ parameter, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'none')
delta_j %>%
  pivot_longer(cols = c(`delta_j[1]`,`delta_j[2]`,`delta_j[3]`,`delta_j[4]`,`delta_j[5]`,`delta_j[6]`),
               names_to = 'parameter', values_to = 'value') %>%
  rename(chain_position = .iteration,
         chain = .chain,
         draw = .draw) %>%
  #filter(chain == 4) %>%
  ggplot(aes(x = chain_position, y = value, colour = as.factor(chain)))+
  geom_line()+
  facet_wrap(. ~ parameter, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'none')

#### posterior predictive check ####
par(mfrow = c(n_windows,1))

## write function
ppcheck <- function(eigen_mat, eigen_df, cent_cov, window, params, rand_node, rand_window){
  plot(density(eigen_mat[1, which(eigen_df$window == window)]), las = 1, ylim = c(0,0.3),
       main = "Posterior predictive check:\nblack = data, blue = predicted",
       col = rgb(0, 0, 0, 0.25))
  n_nodes <- length(which(eigen_df$window == window))
  eigen_data <- data.frame(node_age = eigen_df$age_cat[which(eigen_df$window == window)],
                           window = rep(window, n_nodes),
                           nodes = eigen_df$node_random[eigen_df$window == window],
                           nodes_window = eigen_df$node_random[eigen_df$window == window])
  for (i in 1:100) {
    j <- sample(1:length(params$beta_age), 1)
    lines(density(eigen_mat[j, which(eigen_df$window == window)]), col=rgb(0, 0, 0, 0.25))
    mu <- rep(NA, length(eigen_data$node_age))
    for(k in 1:length(mu)) {
      mu[k] <- params$intercept[j] + params$beta_age[j]*sum(delta_j[j,(1:eigen_data$node_age[k])]) + as.numeric(rand_window[j,eigen_data$window[k]]) + as.numeric(rand_node[j,eigen_data$nodes[k]])
    }
    sigma <- cent_cov + diag(rep(params$sigma[j], n_nodes))
    #lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
    nu <- params$nu[j]
    lines(density(LaplacesDemon::rmvt(n = 1, mu = mu, S = sigma, df = nu)),
          col=rgb(0, 0, 1, 0.25))
  }
}

## check per time window
for(time_window in 1:n_windows){
  ppcheck(eigen_mat = sim_dat, eigen_df = sim,
          cent_cov = sim_cent_cov[[time_window]], window = time_window,
          params = params, rand_node = rand_node, rand_window = rand_window)
}
par(mfrow = c(1,1))

## save output
dev.off()
#pdf('../outputs/step4_nodalregression/mpnpshort_nodalregression_modelpredictions.pdf')

# #### predict from model -- NOT YET CHECKED AND CONVERTED FROM SIMULATION ####
# ## create mean prediction function
# get_mean_predictions <- function(predict_df, delta_j, parameters, include_node = TRUE, include_window = TRUE){
#   ## create empty matrix to fill with predictions
#   mean_matrix <- matrix(NA, nrow = nrow(parameters), ncol = nrow(predict_df),
#                         dimnames = list(NULL, predict_df$node_window))
# 
#   ## populate matrix = mean centrality values per node, predicting for real data
#   for(i in 1:nrow(mean_matrix)){
#     for(j in 1:ncol(mean_matrix)){
#       mean_matrix[i,j] <- parameters$intercept[i] + parameters$beta_age[i] * sum(delta_j[i,(1:predict_df$age_cat[j])])
#     }
#   }
# 
#   if(include_window == TRUE){
#     window_effect <- parameters %>% dplyr::select(grep('rand_window', colnames(parameters), value=TRUE))
#     for(i in 1:nrow(mean_matrix)){
#       mean_matrix[i,] <- mean_matrix[i,] + as.numeric(window_effect[i,predict_df$window])
#     }
#   }
# 
#   if(include_node == TRUE){
#     node_effect <- parameters %>% dplyr::select(grep('rand_node', colnames(parameters), value=TRUE))
#     for(i in 1:nrow(mean_matrix)){
#       mean_matrix[i,] <- mean_matrix[i,] + as.numeric(node_effect[i, predict_df$node_random])
#     }
#   }
# 
#   return(mean_matrix)
# 
# }
# 
# ## get mean predictions
# mu_std <- get_mean_predictions(predict_df = sim, parameters = params, delta_j = delta_j,
#                                include_window = TRUE, include_node = TRUE)
# 
# ## add mean and CI of predicted means to input data frame for comparison
# sim$mu_mean_std <- apply(mu_std, 2, mean)
# sim$mu_lwr_std <- apply(mu_std, 2, quantile, prob = 0.975)
# sim$mu_upr_std <- apply(mu_std, 2, quantile, prob = 0.025)
# 
# ## plot mean of model vs mean of raw data
# ggplot()+
#   geom_point(data = sim, aes(x = mu, y = mu_mean_std, colour = as.factor(window)))+
#   scale_colour_viridis_d()+
#   labs(colour = 'window', x = 'simulated mean (standardised)', y = 'predicted mean (standardised)')+
#   geom_abline(slope = 1, intercept = 0) # add line showing where points would lie if model fit was perfect
# 
# ## put together sigma arrays, separated by time window
# sigma_all <- list()
# for(time_window in 1:n_windows){
#   cent_cov <- eigen_list[[grep('centrality_cov', names(eigen_list))[time_window] ]]
#   n_nodes_window <- eigen_list[[grep('num_nodes_window', names(eigen_list))[time_window] ]]
#   nodes_window <- eigen_list[[grep('nodes_window', names(eigen_list))[time_window + n_windows] ]]
#   sigma_array <- array(NA, dim = c(n_nodes_window,
#                                    n_nodes_window,
#                                    nrow(params)),
#                        dimnames = list(nodes_window,
#                                        nodes_window,
#                                        NULL))
# 
#   for(i in 1:nrow(params)){
#     sigma_array[,,i] <- cent_cov + diag(rep(params$sigma[i], n_nodes_window))
#   }
#   sigma_all[[time_window]] <- sigma_array
# }
# 
# ## create empty matrix to take full set of predicted values per elephant
# predictions_std <- matrix(NA, nrow = nrow(params), ncol = nrow(sim), dimnames = list(NULL, sim$node_window))
# 
# ## populate matrix using mean values in matrix mu_std, and sigma values based on time window
# for(time_window in 1:n_windows){
#   sigma_array <- sigma_all[[time_window]]
#   for(i in 1:nrow(predictions_std)){
#     predictions_std[i,sim$window == time_window] <- MASS::mvrnorm(1, mu_std[i,sim$window == time_window], sigma_array[,,i])
#   }
# }
# 
# ## add CI of predicted data points to input data frame for comparison
# sim$predict_lwr_std <- apply(predictions_std, 2, quantile, prob = 0.025)
# sim$predict_upr_std <- apply(predictions_std, 2, quantile, prob = 0.975)
# 
# ## plot predictions
# ggplot(sim)+
#   geom_ribbon(aes(x = age, ymin = predict_lwr_std, ymax = predict_upr_std, fill = as.factor(window)),
#               alpha = 0.2)+                      # background layer showing the 95% CI of all predictions
#   geom_ribbon(aes(x = age, ymin = mu_lwr_std, ymax = mu_upr_std, fill = as.factor(window)),
#               alpha = 0.4)+                      # mid layer showing the 95% CI of predicted means
#   geom_line(aes(x = age, y = mu_mean_std, colour = as.factor(window)))+  # line showing mean of predicted means
#   geom_point(aes(x = age, y = mu))+              # original data points (standardised centrality, actual age)
#   scale_colour_viridis_d(begin = 0, end = 0.7)+
#   scale_fill_viridis_d(begin = 0, end = 0.7)+
#   facet_wrap(. ~ as.factor(window))+             # separate plots per window
#   theme(legend.position = 'bottom')+
#   labs(colour = 'time window', fill = 'time window',
#        y = 'eigenvector centrality', x = 'age (years)')
# 
# ## save output
# dev.off()
# pdf('step4_nodalregression/checks/simulation_mpnp_extractoriginal.pdf')
# 
# #### extract contrasts from predictions -- TAKE FROM LONG VERSION ONCE YOU'RE SURE IT WORKS ####
# #### plot nicely -- TAKE FROM LONG VERSION ONCE YOU'RE SURE IT WORKS ####
