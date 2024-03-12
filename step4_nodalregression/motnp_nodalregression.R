#### information #####
# 1) to estimate the ages of individuals in other populations, first came up with a survival curve for the Amboseli elephants which have good birth and death data.
# 2) selected the Gompertz bathtub distribution
# 3) use distribution to convert MOTNP ages in categories to probability distributions of true age (first part of this script)
# 4) extract node centrality distribution from draws of edge weight
# 5) use MOTNP age distributions to predict network centrality (third part of this script)

#### set up ####
#options(future.globals.maxSize = 10000*(1024^2))   # running model with full age distributions = 5.01GB of globals, which exceeds default maximum allowed size. Set to allow up to 10 GB to allow model to run

# library(rstan) ; library(sna) ; library(tidyverse) ; library(LaplacesDemon) ; library(MASS)
library(StanHeaders, lib.loc = '../packages/')
library(rstan, lib.loc = '../packages/')
library(sna, lib.loc = '../packages/')             # library(sna)
library(tidyverse, lib.loc = '../packages/')       # library(tidyverse)
library(LaplacesDemon, lib.loc = '../packages/')   # library(LaplacesDemon)
library(MASS, lib.loc = '../packages/')            # library(MASS)
#library(cmdstanr, lib.loc = '../packages/')        # library(cmdstanr)
#library(igraph, lib.loc = '../packages/')          # library(igraph)
#library(brms, lib.loc = '../packages/')            # library(brms)
#library(Rcpp, lib.loc = '../packages/')            # library(Rcpp)
#library(ggdist, lib.loc = '../packages/')          # library(ggdist)
#library(posterior, lib.loc = '../packages/')       # library(posterior)
#library(bayesplot, lib.loc = '../packages/')       # library(bayesplot)
#library(bisonR, lib.loc = '../packages/')          # library(bisonR)
#library(janitor, lib.loc = '../packages/')         # library(janitor)

## set cmdstan path
#set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0/')

## set seed for reproducibility
set.seed(12345)

## set up pdf
pdf('../outputs/motnp_nodalregression.pdf')

## set up theme for plots
theme_set(theme_classic())

#### read in data ####
# load edge weight model and data frames
load('motnp_edgeweights_conditionalprior.RData')
rm(edgelist, x, make_edgelist, plot_network_threshold, i, motnp_ages) ; gc()

# df_nodal <- distinct(counts_df[,c('node_1_males','id_1')])
# colnames(df_nodal) <- c('node_2_males','id_2')
# df_nodal <- rbind(df_nodal, counts_df[nrow(counts_df),c('node_2_males','id_2')])
# colnames(df_nodal) <- c('node','id')
# 
# motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>% 
#   dplyr::select(sort(unique(c(counts_df$id_1, counts_df$id_2)))) %>% 
#   pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')
# motnp_ages <- left_join(motnp_ages, df_nodal, by = 'id')
# motnp_ages$draw <- rep(1:8000, length(unique(motnp_ages$id)))
# 
# mean_motnp_ages <- df_nodal
# mean_motnp_ages$age <- NA
# for(i in 1:nrow(mean_motnp_ages)){
#   x <- motnp_ages[motnp_ages$id == mean_motnp_ages$id[i],]
#   mean_motnp_ages$age[i] <- mean(x$age)
#   rm(x)
# }

## convert age data back to categories
nodes <- nodes %>% 
  rename(id_1 = id) %>% 
  left_join(distinct(counts_df[,c('id_1','age_category_1','age_cat_id_1')]),
            by = 'id_1') %>% 
  rename(id = id_1, age_cat_chr = age_category_1, age_cat_num = age_cat_id_1)
nodes$age_cat_chr[nrow(nodes)] <- counts_df$age_category_2[nrow(counts_df)]
nodes$age_cat_num[nrow(nodes)] <- counts_df$age_cat_id_2[nrow(counts_df)]
nodes$age_cat_num <- as.numeric(nodes$age_cat_num)
nodes$age_cat_fct <- as.factor(nodes$age_cat_num - 2)

## clean up data frame
ele_ids <- nodes$id
n_eles <- length(ele_ids)
edges$chain <- ifelse(edges$chain == 'chain1', 1,
                      ifelse(edges$chain == 'chain2', 2,
                             ifelse(edges$chain == 'chain3', 3, 4)))
edges$draw_id <- edges$position + (edges$chain-1) * 1000
edges <- edges %>% 
  rename(dyad_id = dyad) %>% 
  left_join(counts_df[,c('dyad_id','node_1','node_2','id_1','id_2',
                         'count_1','count_2','event_count','count_dyad')],
            by = 'dyad_id')

## build adjacency tensor
adj_tensor <- array(NA, c(n_eles, n_eles, n_samples*n_chains),
                    dimnames = list(ele_ids, ele_ids, NULL))
for(i in 1:n_dyads) {
  dyad_row <- counts_df[i,]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- edge_samples[,i]
  adj_tensor[dyad_row$id_2, dyad_row$id_1, ] <- edge_samples[,i]
}
adj_tensor[,,1]

#### extract centralities ####
# ## calculate centrality and store posterior samples in a matrix
# centrality_samples <- matrix(0, n_chains*n_samples, n_eles, dimnames = list(NULL,nodes$id))
# #centrality_samples_std <- matrix(0, n_chains*n_samples, n_eles, dimnames = list(NULL,nodes$id))
# for (draw in 1:(n_chains*n_samples)) {
#   g <- graph_from_adjacency_matrix(adj_tensor[,,draw], mode = "undirected",
#                                    diag = FALSE, # TRUE in Jordan's example and in ANP -- shouldn't make a difference as NA anyway
#                                    weighted = TRUE)
#   centrality_samples[draw, ] <- eigen_centrality(g)$vector
#   #centrality_samples_std[draw, ] <- (centrality_samples[draw, ] - mean(centrality_samples[draw, ]))/sd(centrality_samples[draw, ])
# }
# head(centrality_samples)      # unstandardised eigenvector centrality
# #head(centrality_samples_std)  #  standardised  eigenvector centrality

## calculate centrality and store posterior samples in a matrix
centrality_samples_invlogit <- matrix(0, n_chains*n_samples, n_eles,
                                      dimnames = list(NULL,nodes$id))
for (draw in 1:(n_chains*n_samples)) {
  centrality_samples_invlogit[draw, ] <- sna::evcent(dat = adj_tensor[,,draw],
                                                     gmode = 'graph', diag = F)
}
head(centrality_samples_invlogit)      # unstandardised eigenvector centrality

## convert to logit scale
centrality_samples <- logit(centrality_samples_invlogit)

## visualise centralities
nodes$node_rank <- as.integer(as.factor(nodes$node))
df_wide <- data.frame(centrality_samples)
colnames(df_wide) <- 1:n_eles
df_long <- pivot_longer(df_wide, cols = everything(),
                        names_to = "node_rank", values_to = "centrality") %>% 
  mutate(node_rank = as.integer(node_rank)) %>% 
  left_join(nodes[,c('node_rank','age_cat_num')], by = 'node_rank')
df_long %>% 
  filter(node_rank <= 30) %>% 
  mutate(nodes_reordered = fct_reorder(.f = as.factor(node_rank), .x = age_cat_num, .desc = T)) %>% 
  ggplot(aes(x = centrality, fill = as.factor(age_cat_num))) +
  geom_density(linewidth = 0.4) +
  scale_fill_viridis_d() +
  facet_grid(rows = vars(as.factor(nodes_reordered)), scales = "free") +
  labs(x = "Eigenvector centrality",
       fill = 'age cat') + 
  theme_void() + 
  theme(strip.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

#### compute normal approximation ####
## check covariance
plot(centrality_samples[, 1], centrality_samples[, 2],
     xlab = 'standardised eigenvectors ID1',
     ylab = 'standardised eigenvectors ID2',
     las = 1, pch = 19, col = rgb(0,0,1,0.2))
plot(centrality_samples[, which(nodes$age_cat_num == 3)[1]],
     centrality_samples[, which(nodes$age_cat_num == 7)[1]],
     xlab = 'standardised eigenvectors ID1',
     ylab = 'standardised eigenvectors ID2',
     las = 1, pch = 19, col = rgb(0,0,1,0.2))

## check variance
par(mfrow = c(5,5), mai = c(0.2,0.2,0.2,0.2))
to_plot <- sample(1:ncol(centrality_samples), 25, replace = F)
for(i in 1:length(to_plot)){
  hist(centrality_samples[,to_plot[i]], main = '')
}

## compute normal approximation
centrality_mu <- apply(centrality_samples, 2, mean)
centrality_cov <- cov(centrality_samples)
centrality_samples_sim <- MASS::mvrnorm(1e5, centrality_mu, centrality_cov)
for(i in 1:length(to_plot)){
  plot(density(centrality_samples[, i]), lwd = 2, main = "", xlab = "")
  lines(density(centrality_samples_sim[, i]), col = rgb(0,0,1,0.5), lwd = 2)
  
}
par(mfrow = c(1,1), mai = c(1,1,1,1))

save.image('motnp_nodalregression.RData')

#### prior predictive check ####
## set values
n <- 100
n_age_cat <- length(unique(nodes$age_cat_fct))
beta_age <- rnorm(n, 0, 0.8)
intercept  <- rnorm(n, logit(0.05), 2)
age_dirichlet <- rdirichlet(n, c(1,1,1,1,1))

## plot
plot(NULL, las = 1, xlab = 'age category', ylab = 'logit eigenvector',
     ylim = c(min(centrality_mu)-5, max(centrality_mu)+5), xlim = c(1,n_age_cat))
abline(h = min(centrality_mu), lty = 2) ; abline(h = max(centrality_mu), lty = 2)
x <- 1:n_age_cat
for(i in 1:n){
  y <- rep(NA, length(x))
  for(j in 1:length(x)){
    y[j] <- intercept[i] + beta_age[i]*sum(age_dirichlet[i,][1:x[j]])
  }
  lines(x = x, y = y, col = rgb(0,0,1,0.4))
}
rm(n, beta_age, intercept, age_dirichlet, sigma, x, y, df_plot, df_wide) ; gc()

#### run model -- using ordered categorical exposure ####
#load('motnp_nodalregression.RData')

## create data list
eigen_list <- list(num_nodes = n_eles,
                   num_age_cat = n_age_cat,
                   length_dirichlet = n_age_cat+1,
                   centrality_mu = centrality_mu,
                   centrality_cov = centrality_cov,
                   node_age = as.integer(nodes$age_cat_fct),
                   prior_age = rep(1, n_age_cat))

## check inputs
boxplot(centrality_mu ~ nodes$age_cat_fct, notch = T)

## load model
nodal_regression <- stan_model('models/eigen_regression_motnp.stan')

## run model
fit_motnp_eigen <- sampling(nodal_regression,
                            data = eigen_list,
                            cores = n_chains,
                            chains = n_chains)

## save output (get it saved, then clean it up once it hasn't crashed, then save the cleaner version!)
save.image('motnp_nodalregression.RData')
rm(edge_samples, adj_tensor, i, to_plot, draw,summary,centrality_samples_sim, dyad_row) ; gc()
save.image('motnp_nodalregression.RData')

#### check outputs ####
# load('motnp_nodalregression.RData')
## traceplot linear effect size
traceplot(fit_motnp_eigen, pars = c('intercept','beta_age','sigma','predictor[1]','predictor[2]','predictor[3]','predictor[4]','predictor[5]','predictor[6]','predictor[7]','predictor[8]','predictor[9]'))
traceplot(fit_motnp_eigen, pars = c('delta_j[1]','delta_j[2]','delta_j[3]',
                                    'delta_j[4]','delta_j[5]','delta_j[6]'))

## summarise
(summary <- as.data.frame(round(summary(fit_motnp_eigen)$summary, 3)))
summary$parameter <- rownames(summary)
par(mfrow = c(2,1))
hist(summary$Rhat, breaks = 50)
hist(summary$n_eff, breaks = 50)
par(mfrow = c(1,1))

## extract posterior
params <- rstan::extract(fit_motnp_eigen)

## plot posterior age effect
plot(density(params$beta_age),
     main = "Posterior age effect estimate")
abline(v = 0, lty = 2)

#### posterior predictive check ####
plot(density(centrality_samples[1, ]), main="Posterior predictive density of responses:\nblack = data, blue = predicted", col=rgb(0, 0, 0, 0.25), ylim=c(0, 3))
for (i in 1:100) {
  j <- sample(1:(n_chains*n_samples), 1)
  lines(density(centrality_samples[j, ]), col=rgb(0, 0, 0, 0.25))
  mu <- rep(NA, n_eles)
  for(k in 1:n_eles){
    mu[k] <- params$beta_age[j]*sum(params$delta_j[j,1:nodes$age_cat_fct[k]]) + params$intercept[j]
  }
  sigma <- centrality_cov + diag(rep(params$sigma[j], n_eles))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}

#### predict from model -- raw data ####
# load('motnp_nodalregression.RData')
## create functions for predictions
predict_mean_centrality <- function(params, pred_data){
  mu_matrix <- matrix(data = NA,
                      nrow = length(params$beta_age),
                      ncol = nrow(pred_data),
                      dimnames = list(1:length(params$beta_age),
                                      pred_data$node))
  for(i in 1:nrow(mu_matrix)){
    for(j in 1:ncol(mu_matrix)){
      mu_matrix[i,j] <- params$intercept[i] + params$beta_age[i] * sum(params$delta_j[i,(1:pred_data$age_cat_fct[j])])
    }
  }
  return(mu_matrix)
}
predict_full_centrality <- function(params, mu_matrix, cov_matrix){
  full_matrix <- mu_matrix
  for(i in 1:nrow(full_matrix)){
    full_matrix[i,] <- MASS::mvrnorm(n = 1, mu = mu_matrix[i,],
                                     Sigma = cov_matrix + diag(rep(params$sigma[i], ncol(full_matrix))))
  }
  return(full_matrix)
}

## predicted means
sim_mean <- predict_mean_centrality(params = params, pred_data = nodes)
nodes$mu_avg_predict <- apply(sim_mean, 2, mean)
nodes$mu_lwr_predict <- apply(sim_mean, 2, quantile, probs = 0.025)
nodes$mu_upr_predict <- apply(sim_mean, 2, quantile, probs = 0.975)

## full predictions
sim_full <- predict_full_centrality(params = params, mu_matrix = sim_mean, cov_matrix = centrality_cov)
nodes$full_lwr_predict <- apply(sim_full, 2, quantile, probs = 0.025)
nodes$full_upr_predict <- apply(sim_full, 2, quantile, probs = 0.975)

## compare to raw data
nodes <- nodes %>% 
  mutate(mu_raw = centrality_mu,
         mu_raw_invlogit = LaplacesDemon::invlogit(centrality_mu))
ggplot(nodes)+
  geom_ribbon(aes(x = age_cat_num, ymin = full_lwr_predict, ymax = full_upr_predict),
              alpha = 0.3, fill = 'purple')+
  geom_ribbon(aes(x = age_cat_num, ymin = mu_lwr_predict, ymax = mu_upr_predict),
              alpha = 0.3, fill = 'blue')+
  geom_jitter(aes(x = age_cat_num, y = mu_raw), width = 0.2)+
  geom_line(aes(x = age_cat_num, y = mu_avg_predict))+
  theme_classic()

## save image
save.image('motnp_nodalregression.RData')

## summarise predictions -- all the same for all nodes in the same category because age is the only predictor
nodes_full <- sim_full %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'node', values_to = 'node_pred') %>% 
  mutate(draw = rep(1:n_eles, each = nrow(sim_full)),
         node = as.numeric(node),
         invlogit_pred = LaplacesDemon::invlogit(node_pred)) %>% 
  left_join(nodes[,c('id','node','age_cat_chr','age_cat_fct')],
            by = 'node') %>% 
  group_by(id) %>% 
  mutate(pred_full_lwr = quantile(node_pred, probs = 0.025),
         pred_full_mean = mean(node_pred),
         pred_full_mid = quantile(node_pred, probs = 0.5),
         pred_full_upr = quantile(node_pred, probs = 0.975),
         pred_full_lwr_invlogit = quantile(invlogit_pred, probs = 0.025),
         pred_full_mean_invlogit = mean(invlogit_pred),
         pred_full_mid_invlogit = quantile(v, probs = 0.5),
         pred_full_upr_invlogit = quantile(invlogit_pred, probs = 0.975))

ggplot()+
  geom_violin(data = nodes_full,
              aes(x = age_cat_fct, y = node_pred,
                  fill = factor(age_cat_chr, 
                                levels = c('10-15', '15-19', '20-25',
                                           '25-40','40+'))),
              alpha = 0.5)+
  geom_boxplot(data = nodes_full,
               aes(x = age_cat_fct, y = pred_full_mid, 
                   fill = age_cat_chr))+
  geom_jitter(data = nodes,
              aes(x = age_cat_fct, y = mu_raw, 
                  fill = age_cat_chr, size = sightings),
              width = 0.2, pch = 21)+
  scale_fill_viridis_d()+
  #scale_colour_viridis_d()+
  labs(fill = 'age category',
       #colour = 'age category',
       x = 'age category',
       y = 'predicted eigenvector centrality')
ggsave(file = '../outputs/step4_nodalregression/motnp_nodal_violin_logit.png', device = 'png',
       plot = last_plot(), width = 2100, height = 1600, units = 'px')

ggplot()+
  geom_violin(data = nodes_full,
              aes(x = age_cat_fct, y = invlogit_pred,
                  fill = factor(age_cat_chr,
                                levels = c('10-15', '15-19', '20-25',
                                           '25-40','40+'))),
              alpha = 0.5)+
  geom_boxplot(data = nodes_full,
               aes(x = age_cat_fct, y = pred_full_mid_invlogit,
                   fill = age_cat_chr))+
  geom_jitter(data = nodes,
              aes(x = age_cat_fct, y = mu_raw_invlogit,
                  fill = age_cat_chr, size = sightings),
              width = 0.2, pch = 21)+
  scale_fill_viridis_d()+
  #scale_colour_viridis_d()+
  labs(fill = 'age category',
       #colour = 'age category',
       x = 'age category',
       y = 'predicted eigenvector centrality')
ggsave(file = '../outputs/step4_nodalregression/motnp_nodal_violin.png', device = 'png',
       plot = last_plot(), width = 2100, height = 1600, units = 'px')

## save output
save.image('motnp_nodalregression.RData')

#### extract contrasts from predictions ####
## extract raw values
par(mfrow = c(2,2))
beta_diff_12 <- params$delta[, 1] - params$delta[, 2]
plot(density(beta_diff_12), main="Posterior difference:\n10-15 and 15-20") ; abline(v=0, lty=2)
beta_diff_23 <- params$delta[, 2] - params$delta[, 3]
plot(density(beta_diff_23), main="Posterior difference:\n15-20 and 21-25") ; abline(v=0, lty=2)
beta_diff_34 <- params$delta[, 3] - params$delta[, 4]
plot(density(beta_diff_34), main="Posterior difference:\n21-25 and 26-40") ; abline(v=0, lty=2)
beta_diff_45 <- params$delta[, 4] - params$delta[, 5]
plot(density(beta_diff_45), main="Posterior difference:\n26-40 and 40+") ; abline(v=0, lty=2)

## plot
par(mfrow = c(1,1))
plot(density(beta_diff_12), main="Posterior differences")
lines(density(beta_diff_23), col = 'red')
lines(density(beta_diff_34), col = 'blue')
lines(density(beta_diff_45), col = 'purple')
abline(v=0, lty=2)

## predict for age categroy + 1
pred_data_new <- nodes[order(nodes$age_cat_num),] %>% 
    mutate(age_cat_num_original = age_cat_num,
           age_cat_fct_original = age_cat_fct,
           age_cat_num = ifelse(age_cat_num == 7, 3, age_cat_num + 1),
           age_cat_fct = as.integer(as.factor(age_cat_num)))
pred_mu_new <- predict_mean_centrality(params = params,
                                       pred_data = pred_data_new)
pred_full_new <- predict_full_centrality(params = params,
                                         mu_matrix = pred_mu_new,
                                         cov_matrix = centrality_cov)

## calculate contrast -- all changes
contrast <- pred_full_new - sim_full
mean(contrast)
quantile(contrast, prob = c(0.025, 0.975))

## calculate contrast -- 1 vs 2
contrast12 <- contrast[,which(nodes$age_cat_fct == 1)]
mean(contrast12)
quantile(contrast12, prob = c(0.025, 0.975))

## calculate contrast -- 2 vs 3
contrast23 <- contrast[,which(nodes$age_cat_fct == 2)]
mean(contrast23)
quantile(contrast23, prob = c(0.025, 0.975))

## calculate contrast -- 3 vs 4
contrast34 <- contrast[,which(nodes$age_cat_fct == 3)]
mean(contrast34)
quantile(contrast34, prob = c(0.025, 0.975))

## calculate contrast -- 4 vs 5
contrast45 <- contrast[,which(nodes$age_cat_fct == 4)]
mean(contrast45)
quantile(contrast45, prob = c(0.025, 0.975))

## calculate contrast -- 5 vs 1
contrast51 <- contrast[,which(nodes$age_cat_fct == 5)]
mean(contrast51)
quantile(contrast51, prob = c(0.025, 0.975))

## plot contrasts
contrasts <- data.frame(contrast = c('1 vs 2', '2 vs 3', '3 vs 4', '4 vs 5'),
                        mean = c(mean(contrast12),mean(contrast23),
                                 mean(contrast34),mean(contrast45)),
                        upr = c(quantile(contrast12, prob = c(0.975)),
                                quantile(contrast23, prob = c(0.975)),
                                quantile(contrast34, prob = c(0.975)),
                                quantile(contrast45, prob = c(0.975))),
                        lwr = c(quantile(contrast12, prob = c(0.025)),
                                quantile(contrast23, prob = c(0.025)),
                                quantile(contrast34, prob = c(0.025)),
                                quantile(contrast45, prob = c(0.025)))) %>% 
  pivot_longer(cols = c('mean', 'upr', 'lwr'),
               names_to = 'stat', values_to = 'value')
ggplot()+
  geom_line(data = contrasts[contrasts$stat != 'mean',],
            aes(y = value, x = contrast))+
  geom_hline(aes(yintercept = 0), linetype = 2)+
  geom_point(data = contrasts[contrasts$stat == 'mean',],
             aes(y = value, x = contrast))+
  scale_x_discrete(#expand = c(0.05,0.05),
                   name = 'age category comparison')+
  scale_y_continuous(name = 'difference in prediction',
                     limits = c(-1,1))+
  theme_bw()+
  coord_flip()

## save
save.image('motnp_nodalregression.RData')

#### plot nicely ####
## clean data frame
clean <- data.frame(age_cat_fct = 1:5,
                    mean = NA,
                    predict_mu_lwr = NA, predict_mu_upr = NA,
                    predict_full_lwr = NA, predict_full_upr = NA,
                    mean_invlogit = NA,
                    predict_mu_lwr_invlogit = NA, predict_mu_upr_invlogit = NA,
                    predict_full_lwr_invlogit = NA, predict_full_upr_invlogit = NA)

## predict means only
predict_clean_mean <- predict_mean_centrality(params, clean)
clean$mean <- mean(predict_clean_mean)
clean$predict_mu_lwr <- apply(predict_clean_mean, 2, quantile, prob = 0.025)
clean$predict_mu_upr <- apply(predict_clean_mean, 2, quantile, prob = 0.975)

## convert to invlogit scale
predict_clean_mean_invlogit <- LaplacesDemon::invlogit(predict_clean_mean)
clean$mean_invlogit <- mean(predict_clean_mean_invlogit)
clean$predict_mu_lwr_invlogit <- apply(predict_clean_mean_invlogit, 2, quantile, prob = 0.025)
clean$predict_mu_upr_invlogit <- apply(predict_clean_mean_invlogit, 2, quantile, prob = 0.975)

## average fulls for whole data frame (can't just predict from clean data frame because dimensions don't fit with covariance matrix)
for(i in 1:n_age_cat){
  sim_age <- sim_full[,which(nodes$age_cat_fct == i)]
  sim_age_invlogit <- LaplacesDemon::invlogit(sim_age)
  clean$predict_full_lwr[clean$age_cat_fct == i] <- mean(apply(sim_age, 2,
                                                               quantile, prob = 0.025))
  clean$predict_full_upr[clean$age_cat_fct == i] <- mean(apply(sim_age, 2,
                                                               quantile, prob = 0.975))
  clean$predict_full_lwr_invlogit[clean$age_cat_fct == i] <- mean(apply(sim_age_invlogit, 2,
                                                                        quantile, prob = 0.025))
  clean$predict_full_upr_invlogit[clean$age_cat_fct == i] <- mean(apply(sim_age_invlogit, 2,
                                                                        quantile, prob = 0.975))
}

## add node data
df_long <- df_long %>%
  mutate(centrality_invlogit = LaplacesDemon::invlogit(centrality)) %>% 
  left_join(nodes, by = 'node_rank') %>% 
  dplyr::select(-age_cat_num.x) %>% 
  rename(age_cat_num = age_cat_num.y)

## plot on logit scale
(clean_plot <- ggplot()+
  geom_ribbon(data = clean,
              aes(x = as.numeric(age_cat_fct),
                  ymin = predict_full_lwr, ymax = predict_full_upr),    # shade simulations
              colour = 'transparent', fill = rgb(0,0,0,0.1))+
  geom_ribbon(data = clean,
              aes(x = as.numeric(age_cat_fct),
                  ymin = predict_mu_lwr, ymax = predict_mu_upr),  # shade mean distribution
              colour = 'transparent',
              fill = rgb(33/255, 145/255, 140/255, 0.5))+
  geom_line(data = clean,
            aes(x = as.numeric(age_cat_fct),
                y = mean),                            # mean line
            colour = rgb(33/255, 145/255, 140/255),
            linewidth = 1)+
  scale_x_continuous('age category')+
  scale_y_continuous('eigenvector centrality')+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22)))

clean_plot +
  geom_point(data = nodes,
           aes(x = as.numeric(age_cat_fct),
               y = mu_raw,
               size = sightings),
           colour = rgb(68/255, 1/255, 84/255))
ggsave(filename = '../outputs/step4_nodalregression/motnp_nodalregression_line_meanpoints_logit.png', device = 'png',
       plot = last_plot(), width = 2800, height = 1600, units = 'px')

clean_plot +
  geom_point(data = df_long,
             aes(x = as.numeric(age_cat_fct),
                 y = centrality),
             colour = rgb(253/255, 231/255, 37/255, 0.01)) +
  geom_point(data = nodes,
             aes(x = as.numeric(age_cat_fct),
                 y = mu_raw,
                 size = sightings),
             colour = rgb(68/255, 1/255, 84/255))
ggsave(filename = '../outputs/step4_nodalregression/motnp_nodalregression_line_allpoints_logit.png', device = 'png',
       plot = last_plot(), width = 2800, height = 1600, units = 'px')

## plot on invlogit scale
(clean_plot <- ggplot()+
    geom_ribbon(data = clean,
                aes(x = as.numeric(age_cat_fct),
                    ymin = predict_full_lwr_invlogit,
                    ymax = predict_full_upr_invlogit),           # shade simulations
                colour = 'transparent', fill = rgb(0,0,0,0.1))+
    geom_ribbon(data = clean,
                aes(x = as.numeric(age_cat_fct),
                    ymin = predict_mu_lwr_invlogit,
                    ymax = predict_mu_upr_invlogit),             # shade mean distribution
                colour = 'transparent',
                fill = rgb(33/255, 145/255, 140/255, 0.5))+
    geom_line(data = clean,
              aes(x = as.numeric(age_cat_fct),
                  y = mean_invlogit),                            # mean line
              colour = rgb(33/255, 145/255, 140/255),
              linewidth = 1)+
    scale_x_continuous('age category')+
    scale_y_continuous('eigenvector centrality')+
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 22),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 22)))

clean_plot +
  geom_point(data = nodes,
             aes(x = as.numeric(age_cat_fct),
                 y = mu_raw_invlogit,
                 size = sightings),
             colour = rgb(68/255, 1/255, 84/255))
ggsave(filename = '../outputs/step4_nodalregression/motnp_nodalregression_line_meanpoints_invlogit.png',
       device = 'png', plot = last_plot(), width = 2800, height = 1600, units = 'px')

clean_plot +
  geom_point(data = df_long,
             aes(x = as.numeric(age_cat_fct),
                 y = centrality_invlogit),
             colour = rgb(253/255, 231/255, 37/255, 0.01)) +
  geom_point(data = nodes,
             aes(x = as.numeric(age_cat_fct),
                 y = mu_raw_invlogit,
                 size = sightings),
             colour = rgb(68/255, 1/255, 84/255))
ggsave(filename = '../outputs/step4_nodalregression/motnp_nodalregression_line_allpoints_invlogit.png', device = 'png',
       plot = last_plot(), width = 2800, height = 1600, units = 'px')

## save
save.image('motnp_nodalregression.RData')
dev.off()
