##### information #####
# 1) to estimate the ages of individuals in other populations, first came up with a survival curve for the Amboseli elephants which have good birth and death data.
# 2) selected the Gompertz bathtub distribution
# 3) use distribution to convert MOTNP ages in categories to probability distributions of true age (first part of this script)
# 4) extract node centrality distribution from draws of edge weight
# 5) use MOTNP age distributions to predict network centrality (third part of this script)

##### set up ####
#options(future.globals.maxSize = 10000*(1024^2))   # running model with full age distributions = 5.01GB of globals, which exceeds default maximum allowed size. Set to allow up to 10 GB to allow model to run

# library(rstan) ; library(igraph) ; library(tidyverse) ; library(LaplacesDemon) ; library(MASS)
library(rstan, lib.loc = '../packages/')
library(igraph, lib.loc = '../packages/')          # library(igraph)
library(tidyverse, lib.loc = '../packages/')       # library(tidyverse)
library(LaplacesDemon, lib.loc = '../packages/')   # library(LaplacesDemon)
library(MASS, lib.loc = '../packages/')            # library(MASS)
#library(cmdstanr, lib.loc = '../packages/')        # library(cmdstanr)
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

##### prior predictive check ####
## simulate
age <- 10:60
mu_mean <- 0
mu_stdv <- 0.1
mu <- rnorm(100, mu_mean, mu_stdv)
mu <- sort(mu)
mu_mtrx <- matrix(NA, nrow = length(mu), ncol = length(age))
for(i in 1:nrow(mu_mtrx)){
  mu_mtrx[i,] <- mu[i]*age
}
sigma_range <- rexp(25, 2)
sigma_range <- sort(sigma_range)
par(mfrow = c(5,5), mai = c(0.2,0.2,0.2,0.2))
for(j in 1:25){
  plot(NULL, xlim = c(10,60), ylim = c(0,1), las = 1,
       xlab = '', ylab = '')
  sigma <- diag(rep(sigma_range[j], length(age)))
  for(i in 1:length(mu)){
    predictor <- mu_mtrx[i,]
    y <- MASS::mvrnorm(1, predictor, sigma)
    lines(x = age, y = invlogit(y), col = rgb(0,0,1,1))
  }
}
par(mfrow = c(1,1))
rm(list = ls()) ; gc()

##### read in data ####
# load edge weight model and data frames
load('motnp_edgeweights_conditionalprior.RData')
rm(edgelist, x, make_edgelist, plot_network_threshold, i) ; gc()

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

## clean up data frame
ele_ids <- unique(c(counts_df$id_1, counts_df$id_2))
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
}
adj_tensor[,,1]

### extract centralities ####
## calculate centrality and store posterior samples in a matrix
centrality_samples <- matrix(0, n_chains*n_samples, n_eles, dimnames = list(NULL,nodes$id))
centrality_samples_std <- matrix(0, n_chains*n_samples, n_eles, dimnames = list(NULL,nodes$id))
for (draw in 1:(n_chains*n_samples)) {
  g <- graph_from_adjacency_matrix(adj_tensor[,,draw], mode = "undirected",
                                   diag = FALSE, # TRUE in Jordan's example and in ANP -- shouldn't make a difference as NA anyway
                                   weighted = TRUE)
  centrality_samples[draw, ] <- eigen_centrality(g)$vector
  centrality_samples_std[draw, ] <- (centrality_samples[draw, ] - mean(centrality_samples[draw, ]))/sd(centrality_samples[draw, ])
}
head(centrality_samples)      # unstandardised eigenvector centrality
head(centrality_samples_std)  #  standardised  eigenvector centrality

## visualise centralities
nodes$node_rank <- as.integer(as.factor(nodes$node))
df_wide <- data.frame(centrality_samples_std)
colnames(df_wide) <- 1:n_eles
df_long <- pivot_longer(df_wide, cols = everything(),
                        names_to = "node_rank", values_to = "centrality") %>% 
  mutate(node_rank = as.integer(node_rank)) %>% 
  left_join(nodes[,c('node_rank','age')], by = 'node_rank')
df_long %>% 
  filter(node_rank <= 30) %>% 
  mutate(nodes_reordered = fct_reorder(.f = as.factor(node_rank), .x = age, .desc = T)) %>% 
  ggplot(aes(x = centrality, fill = age)) +
  geom_density(linewidth = 0.4) +
  facet_grid(rows = vars(as.factor(nodes_reordered)), scales = "free") +
  labs(x = "Eigenvector centrality (standardised)") + 
  theme_void() + 
  theme(strip.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

### compute normal approximation ####
## check covariance
plot(centrality_samples_std[, 1], centrality_samples_std[, 2],
     xlab = 'standardised eigenvectors ID1',
     ylab = 'standardised eigenvectors ID2',
     las = 1, pch = 19, col = rgb(0,0,1,0.2))

## check variance
par(mfrow = c(5,5), mai = c(0.1,0.1,0.1,0.1))
to_plot <- sample(1:ncol(centrality_samples), 25, replace = F)
for(i in 1:length(to_plot)){
  hist(centrality_samples[,to_plot[i]], main = '')
}
par(mfrow = c(1,1), mai = c(1,1,1,1))

## compute normal approximation
centrality_mu <- apply(centrality_samples_std, 2, mean)
centrality_cov <- cov(centrality_samples_std)

centrality_samples_sim <- MASS::mvrnorm(1e5, centrality_mu, centrality_cov)

plot(density(centrality_samples_std[, 1]), lwd = 2, main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
lines(density(centrality_samples_sim[, 1]), col = rgb(0,0,1,0.5), lwd = 2)

save.image('motnp_nodalregression_conditionaledge.RData')

### run model -- CURRENTLY RUNNING USING NORMAL APPROXIMATION OF AGE, WILL NEED TO COME BACK TO THIS AND WORK OUT HOW TO DO IT PROPERLY ONCE FIXED AGE MODEL ####
load('motnp_nodalregression_conditionaledge.RData')

## extract age distributions (using a normal for now)
motnp_ages <- readRDS('../data_processed/step2_ageestimation/motnp_ageestimates_mcmcoutput.rds')
motnp_ages <- motnp_ages[,which(colnames(motnp_ages) %in% ele_ids)]
age_mu <- apply(motnp_ages, 2, mean)
age_sd <- apply(motnp_ages, 2, sd)

## create data list
eigen_list <- list(num_nodes = n_eles,
                   nodes = nodes$node_rank,
                   centrality_mu = centrality_mu,
                   centrality_cov = centrality_cov,
                   age_mu = age_mu,
                   age_sd = age_sd)

# load model
nodal_regression <- stan_model('models/eigen_regression_agedist.stan')

## run model
fit_motnp_eigen <- sampling(nodal_regression,
                            data = eigen_list,
                            cores = n_chains,
                            chains = n_chains)

## save output (get it saved, then clean it up once it hasn't crashed, then save the cleaner version!)
save.image('motnp_nodalregression_conditionaledge_rstan.RData')
rm(g, edge_samples, adj_tensor, i, to_plot, draw,summary,centrality_samples_sim, dyad_row,motnp_ages) ; gc()
save.image('motnp_nodalregression_conditionaledge_rstan.RData')

### posterior check ####
# load('motnp_nodalregression_conditionaledge_rstan.RData')
## traceplot linear effect size
traceplot(fit_motnp_eigen, pars = c('beta_age','sigma','predictor[1]','predictor[2]','predictor[3]','predictor[4]','predictor[5]','predictor[6]','predictor[7]','predictor[8]','predictor[9]','predictor[10]'))

## posterior predictive check
params <- rstan::extract(fit_motnp_eigen)
plot(density(centrality_samples_std[1, ]), main="Posterior predictive density of responses (standardised centrality)", col=rgb(0, 0, 0, 0.25), ylim=c(0, 1))
for (i in 1:100) {
  j <- sample(1:(n_chains*n_samples), 1)
  lines(density(centrality_samples_std[j, ]), col=rgb(0, 0, 0, 0.25))
  mu <- params$beta_age[j]*params$node_age[j,]
  sigma <- centrality_cov + diag(rep(params$sigma[j], n_eles))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}

## interpret model
plot(density(params$beta_age), # don't need to do the contrast because continuous -- is there an equivalent using the do operator? by standardising have we already dealt with the logit transformation?
     main = "Posterior difference between node types")
abline(v = 0, lty = 2)

### predict from model ####
(summary <- as.data.frame(round(summary(fit_motnp_eigen)$summary[1:2, c(1, 4, 8)], 3)))
summary$parameter <- rownames(summary)

## calculate mean predictions for model
mod_mu <- data.frame(age = floor(min(nodes$age)):ceiling(max(nodes$age))) %>% 
  mutate(lwr = age*summary$`2.5%`[1],
         mid = age*summary$mean[1],
         upr = age*summary$`97.5%`[1])

## plot mean predictions
plot(mid ~ age, data = mod_mu, type = 'l', las = 1, col = 'blue', lwd = 2,
     ylim = c(min(df_long$centrality), max(df_long$centrality)),
     xlab = 'age (years)', ylab = 'eigenvector centrality (standardised)')
lines(lwr ~ age, data = mod_mu, lty = 2)
lines(upr ~ age, data = mod_mu, lty = 2)

## simulate full predictions for model
sim <- matrix(NA, nrow = n_chains*n_samples, ncol = nrow(mod_mu),
              dimnames = list(1:(n_chains*n_samples), mod_mu$age))
for(i in 1:nrow(sim)){
  for(j in 1:ncol(sim)){
    #sim[i,j] <- params$beta_age[i]*x$age[j]
    sim[i,j] <- MASS::mvrnorm(n = 1, mu = params$beta_age[i]*mod_mu$age[j],
                              Sigma = params$sigma[i])
  }
}

## plot simulations
sim_df <- as.data.frame(sim) %>% 
  pivot_longer(cols = everything(), names_to = 'age', values_to = 'eigen_sim')
points(sim_df$eigen_sim ~ sim_df$age, col = rgb(0,0,0,0.01), pch = 19, cex = 0.5)

## summarise simulations
sim_summary <- data.frame(age = as.numeric(unique(sim_df$age)),
                          lwr = NA, mid = NA, upr = NA)
for(i in 1:nrow(sim_summary)){
  x <- sim_df %>% filter(age == sim_summary$age[i])
  sim_summary$lwr[i] <- quantile(x$eigen_sim, 0.025)
  sim_summary$mid[i] <- quantile(x$eigen_sim, 0.5)
  sim_summary$upr[i] <- quantile(x$eigen_sim, 0.975)
}

## plot raw with model output
ages <- data.frame(age_mean = age_mu,
                   age_stdv = age_sd)
ages$id <- rownames(ages)
nodes <- nodes %>% 
  left_join(ages, by = 'id') %>% 
  mutate(age_lwr = age_mean - age_stdv,
         age_upr = age_mean + age_stdv) %>% 
  left_join(distinct(df_long[,c('node_rank','mean_eigen')]))
length(which((nodes$age == nodes$age_mean) == FALSE))

df_long <- df_long %>%
  group_by(node_rank) %>% 
  mutate(mean_eigen = mean(centrality)) %>% 
  ungroup() %>% 
  left_join(nodes[,c('node_rank','age_stdv','age_lwr','age_upr','sightings')], by = 'node_rank') %>% 
  rename(age_mean = age)

ggplot()+
  geom_ribbon(data = sim_summary, aes(x = age, ymin = lwr, ymax = upr),          # shade simulations
              colour = 'transparent', fill = rgb(0,0,0,0.1))+
  geom_ribbon(data = mod_mu, aes(x = age, ymin = lwr, ymax = upr),               # shade mean distribution
              colour = 'transparent', fill = rgb(33/255, 145/255, 140/255, 0.5))+
  geom_point(data = df_long, aes(x = age_mean, y = centrality),                  # all eigenvector draws
             colour = rgb(253/255, 231/255, 37/255, 0.01))+
  geom_point(data = nodes, aes(x = age_mean, y = mean_eigen, size = sightings),  # mean eigenvector
             colour = rgb(68/255, 1/255, 84/255))+
  geom_line(data = mod_mu, aes(x = age, y = mid),                                # mean line
            colour = rgb(33/255, 145/255, 140/255), linewidth = 1)+
  geom_errorbar(data = nodes, aes(xmin = age_lwr, xmax = age_upr,                # age distribution
                                    y = mean_eigen, group = node_rank),
            colour = rgb(68/255, 1/255, 84/255), linewidth = 0.5, width = 0)+
  scale_x_continuous('age (years)')+
  scale_y_continuous('eigenvector centrality (standardised)')+
  theme_classic()+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22))
ggsave(filename = '../outputs/motnp_nodalregression.png', device = 'png',
       plot = last_plot(), width = 16.6, height = 11.6)

## save output
save.image('motnp_nodalregression_conditionaledge_rstan.RData')
dev.off()
rm(list = ls()[!ls() %in% 'nodal_regression'])

# below here everything is old ----------------
### set priors ####
# define PDF output
pdf('../outputs/motnp_nodalregression_plots_meanage.pdf')

# prior predictive check
priors <- get_default_priors('binary')
priors$fixed
prior_check(priors, 'binary')

age <- 1:60
beta_mu <- 0
beta_sigma <- 0.005

#plot(NULL, xlim = c(0,60), ylim = c(0,1), las = 1,
#     xlab = 'age', ylab = 'eigenvector centrality')
#for(i in 1:100){
#  intercept <- rbeta(1,1,1)
#  beta <- rnorm(1, beta_mu, beta_sigma)
#  lines(x = age, y = intercept + age*beta, col = 'blue')
#  lines(x = age, y = intercept + age*plogis(beta), col = 'red') # can't ever be a negative effect so that can't be right
#}

mean_age <- mean(motnp_ages$age)
plot(NULL, xlim = c(10,60), ylim = c(0,1), las = 1,
     xlab = 'age', ylab = 'eigenvector centrality')
for(i in 1:100){
  intercept <- rbeta(1,2,2)   # this isn't right but I don't think there is a prior for the intercept?? I've gone for a symmetrical one here that in itself explores most of the parameter space and allows it to see whether some of the lines are steep enough to go from top to bottom, but on the assumption that when combined, they will explore only the space relevant to their starting position
  beta <- rnorm(1, beta_mu, beta_sigma)
  lines(x = age, y = intercept + (age - mean_age)*beta, col = rgb(0,0,1,0.5)) # vast majority come out somewhere sensible, and those that don't would if they started at a different value for age 10 so that comes down to my ability to work out what the intercept prior should actually be -- think this is a good prior for the slope
}

motnp_edge_weights_strongpriors$model_data$prior_fixed_mu    <- beta_mu
motnp_edge_weights_strongpriors$model_data$prior_fixed_sigma <- beta_sigma

### run model -- mean age value only ####
mean_motnp_eigen <- bison_brm(
  bison(node_eigen(node)) ~ age,
  motnp_edge_weights_strongpriors,
  mean_motnp_ages,
  chains = 4,
  iter = 10000,
  thin = 2
)
summary(mean_motnp_eigen) # FIT 100 IMPUTED MODELS, 4 CHAINS PER MODEL, EACH 1000 DRAWS LONG (+1000 DRAWS WARM UP). WARNING AT END OF MODEL RUN THAT CHAINS <3 DRAWS LONG AS ACTUAL CHAIN IS ONLY 1 WARMUP AND 1 SAMPLE. ONLY IMPUTED CHAINS FOLLOW THE SPECIFIED ITERATIONS AND THINNING

## save output
save.image('motnp_nodalregression_meanage.RData')

# posterior check v1 -- some is good but forgot to draw from posterior so most is rubbish ####
#load('motnp_nodalregression_meanage.RData')
#rm(df_nodal, motnp_ages) ; gc()

summary(mean_motnp_eigen)

# plot
mean_eigen_values <- mean_motnp_eigen$data
plot(mean_eigen_values$bison_node_eigen ~ mean_eigen_values$age, 
     las = 1, pch = 19, col = rgb(0,0,1,0.2),
     xlab = 'mean age estimate', ylab = 'eigenvector centrality',
     main = 'effect of age on eigenvector centrality')

### plot rhat values to check convergence
hist(mean_motnp_eigen$rhats[,2], las = 1, main = 'Rhat values for 100 imputed model runs', xlab = 'Rhat')

### plot traces to check mixing
post_eigen <- as.data.frame(as_draws_df(mean_motnp_eigen)) %>% janitor::clean_names()
plot(data = post_eigen[post_eigen$chain == 1,], b_age ~ draw, type = 'l', xlim = c(0,10000))
lines(data = post_eigen[post_eigen$chain == 2,], b_age ~ draw, col = 'red')
lines(data = post_eigen[post_eigen$chain == 3,], b_age ~ draw, col = 'blue')
lines(data = post_eigen[post_eigen$chain == 4,], b_age ~ draw, col = 'green')
lines(data = post_eigen[post_eigen$chain == 5,], b_age ~ draw, col = 'purple')
unique(post_eigen$chain)

post_eigen$imputation <- rep(1:100, each = 10000)
post_eigen$chain_imp  <- as.factor(rep(1:4, each = length(post_eigen$chain[post_eigen$chain == 1])))
ggplot(post_eigen, aes(x = iteration, y = b_age, colour = chain_imp))+
  geom_line()+
  facet_wrap(imputation ~ .)   # looks good, well mixed

### plot slopes
hist(post_eigen$b_age)
plot(conditional_effects(mean_motnp_eigen), points = TRUE,
     #ylab = 'mean eigenvector centrality',
     theme = theme_classic(base_size = 14))  # older elephants have lower network centrality than younger

### boxplot original categories against eigenvector
# recreate categories
mean_eigen_values$age_cat <- ifelse(mean_eigen_values$age < 15, '10-15',
                                    ifelse(mean_eigen_values$age < 20, '16-20',
                                           ifelse(mean_eigen_values$age < 25, '21-25',
                                                  ifelse(mean_eigen_values$age < 40, '25-40', '40+'))))
mean_eigen_values$age_cat <- factor(mean_eigen_values$age_cat,
                                    levels = c('10-15','16-20','21-25','25-40','40+'))

# combine age and count data
data <- left_join(mean_eigen_values, mean_motnp_ages, by = 'age')
which(is.na(data$id))
length(unique(data$id))
counts1 <- counts_df[,c('id_1','count_1')] %>% distinct 
colnames(counts1) <- c('id','count')
counts2 <- counts_df[,c('id_2','count_2')] %>% distinct
colnames(counts2) <- c('id','count')
counts <- rbind(counts1, counts2) %>% distinct()
rm(counts1, counts2) ; gc()
data <- left_join(data, counts, by = 'id')

# plot boxplot with points sized by count number
ggplot(data = data, aes(x = age_cat, y = bison_node_eigen,
                        fill = factor(mean_eigen_values$age_cat,
                                      levels = c('40+','25-40','21-25','16-20','10-15'))))+
  geom_boxplot(notch = T)+
  geom_jitter(width = 0.2, #shape = 1,
              mapping = aes(size = count))+
  scale_x_discrete(name = 'age category')+
  scale_y_continuous(name = 'mean eigenvector centrality')+
  scale_fill_viridis_d()+
  theme_classic()+
  theme(legend.position = 'none', axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))

# compare empirical distribution to posterior predictive distribution
y <- mean_eigen_values$bison_node_eigen                    # extract eigenvector centralities
yrep <- posterior_predict(mean_motnp_eigen, draws = 500)   # make predictions of eigenvector centrality
dim(yrep)
ppc_dens_overlay(y, yrep[1:1000, ])                        # plot 1000 predictions over empirical distribution (is it ok that this slightly extends beyond x = 1??)
ppc_hist(y, yrep[1:55, ])                                  # compare 55 predictions to empirical distribution
ppc_pit_ecdf(y, yrep[1:50,])                               # no idea....
ppc_ecdf_overlay(y, yrep[1:50,]) + xaxis_text()            # also no idea...

## so... 
# conditional effects plot shows a slight negative effect of age on centrality
# boxplot indicates that older pubescents and young adults are the most central
# effects are weak: all categories show nearly full spread -- BIGGEST EFFECT SEEM TO BE NUMBER OF INDIVIDUAL SIGHTINGS: MORE SIGHTINGS = LOWER CENTRALITY
# model very good at predicting

# posterior check v2 -- run using SRI and see how it compares ####
#library(igraph) ; library(tidyverse) ; library(sna) ; library(brms)
#load('motnp_bisonr_edgescalculated_strongprior.RData')
#load('motnp_nodalregression_meanage.RData')
#rm(df_nodal, mean_eigen_summary, priors, age, beta, beta_mu, beta_sigma, i, intercept, mean_age) ; gc()

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
edges <- counts_df_model[,c('node_1_id','node_2_id','sri')]
colnames(edges)[1:2] <- c('from','to')
nodes <- data.frame(node_1_males = unique(c(counts_df$node_1_males, counts_df$node_2_males))) %>% 
  left_join(distinct(counts_df[,c('node_1_males','age_category_1','count_1')]), by = 'node_1_males') %>% 
  mutate(node_2_males = node_1_males) %>% 
  left_join(distinct(counts_df[,c('node_2_males','age_category_2','count_2')]), by = 'node_2_males') %>% 
  mutate(node = node_1_males,
         age = ifelse(is.na(age_category_1) == FALSE, age_category_1, age_category_2),
         count = ifelse(is.na(count_1) == FALSE, count_1, count_2)) %>% 
  select(node, age, count)
nodes$model <- mean_eigen_values$bison_node_eigen

mean_eigen_values <- left_join(mean_eigen_values, mean_motnp_ages, by = 'age')
length(which(is.na(mean_eigen_values$bison_node_eigen) == TRUE))
nodes <- left_join(nodes, mean_eigen_values, by = 'node')
nodes$model == nodes$bison_node_eigen
colnames(nodes)[c(2,6)] <- c('age_cat','age_years')
nodes <- nodes %>% select(-bison_node_eigen)
nodes <- nodes[,c(1,6,2,5,3,4)]

## check same number of individuals in every category -- yes there are
nrow(nodes) == motnp_edge_weights_strongpriors$num_nodes
nodes_age <- table(nodes$age)
barplot(nodes_age)
model_age <- table(round(mean_motnp_ages$age,0))
barplot(model_age)
model_age <- as.data.frame(model_age)
nodes_age <- as.data.frame(nodes_age)
model_age$Freq == nodes_age$Freq
rm(motnp_edge_weights_strongpriors, motnp_edges_null_strongpriors, nodes_age, mean_eigen_values, mean_motnp_ages) ; gc()

## create adjacency matrix
adj_mat <- matrix(0, nrow = nrow(nodes), ncol = nrow(nodes))
for (i in 1:nrow(edges)) {
  dyad_row <- edges[i, ]
  adj_mat[dyad_row$from, dyad_row$to] <- dyad_row$sri
  adj_mat[dyad_row$to, dyad_row$from] <- dyad_row$sri
}

# Generate igraph object
g <- graph_from_adjacency_matrix(adj_mat, 
                                 mode="undirected", weighted=TRUE)

# check igraph object
#coords <- layout_nicely(g)
#plot(g,
#     edge.width = E(g)$weight,
#     vertex.size = nodes$count,
#     vertex.color = as.factor(nodes$age),
#     edge.color = rgb(0, 0, 0, 0.25),
#     layout = coords)

# calculate centrality
igraph_cent <- eigen_centrality(g, directed = F) # how do I know what order these are calculating?? is it node order, or order nodes appear in E(g)??

# check order using randomised vertex list and see if they come out the same with igraph
new_order <- sample(nodes$node, size = length(nodes$node), replace = F)
adj_mat_randomised <- matrix(0, nrow = nrow(nodes), ncol = nrow(nodes),
                             dimnames = list(x = new_order, y = new_order))
for (i in 1:nrow(edges)) {
  dyad_row <- edges[i, ]
  adj_mat_randomised[which(rownames(adj_mat_randomised) == dyad_row$from[1]), 
                     which(colnames(adj_mat_randomised) == dyad_row$to[1])] <- dyad_row$sri
  adj_mat_randomised[which(rownames(adj_mat_randomised) == dyad_row$to[1]), 
                     which(colnames(adj_mat_randomised) == dyad_row$from[1])] <- dyad_row$sri
}
g_randomised <- graph_from_adjacency_matrix(adj_mat_randomised, 
                                            mode="undirected", weighted=TRUE)
c_randomised <- eigen_centrality(g_randomised, directed = F)
c_randomised$vector == igraph_cent$vector # no -- they are not automatically reverting to node name order
V(g_randomised)
nodes$igraph <- NA # ; nodes$random <- NA
for(i in 1:nrow(nodes)){
  #random_position <- which(new_order == nodes$node[i])
  #nodes$random[i] <- c_randomised$vector[random_position]
  nodes$igraph[i] <- igraph_cent$vector[i]
}
round(nodes$igraph,5) == round(nodes$random,5) # yes -- calculates in order that igraph object was produced
nodes <- nodes %>% select(-random)

# check if individual with eigenvalue = 1 is the same in both data sets (spoiler -- no.)
nodes$node[which(nodes$model == max(nodes$model))]   # 135
nodes$node[which(nodes$igraph == max(nodes$igraph))] # 187

# plot boxplot
ggplot(data = nodes, aes(x = age_cat, y = igraph,
                         fill = factor(age_cat, levels = c('40+','25-40','20-25','15-19','10-15'))))+
  geom_boxplot(notch = T)+
  geom_jitter(width = 0.2, shape = 1,
              aes(size = count))+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_fill_viridis_d()+
  scale_x_discrete(name = 'age category')+
  scale_y_continuous(name = 'eigenvector centrality')

ggplot(data = nodes)+
  geom_point(aes(x = count, y = model,
                 colour = factor(age_cat, levels = c('40+','25-40','20-25','15-19','10-15'))),
             shape = 19)+
  geom_smooth(aes(x = count, y = model), colour = 'red')+
  geom_point(aes(x = count, y = igraph,
                 colour = factor(age_cat, levels = c('40+','25-40','20-25','15-19','10-15'))),
             shape = 4)+
  geom_smooth(aes(x = count, y = igraph))+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_viridis_d()+
  scale_x_continuous(name = 'count')+
  scale_y_continuous(name = 'eigenvector centrality')

## repeat using only individuals sighted 10 or more times
nodes10 <- nodes[nodes$count >= 10,]
edges10 <- edges[edges$from %in% nodes10$node & edges$to %in% nodes10$node,]
adj_mat <- matrix(0, nrow = nrow(nodes10), ncol = nrow(nodes10))
colnames(adj_mat) <- nodes10$node
rownames(adj_mat) <- nodes10$node

for (i in 1:nrow(edges10)) {
  dyad_row <- edges10[i, ]
  adj_mat[which(rownames(adj_mat) == dyad_row$from),
          which(colnames(adj_mat) == dyad_row$to)] <- dyad_row$sri
}

# Generate igraph object
g <- graph_from_adjacency_matrix(adj_mat, 
                                 mode="undirected", weighted=TRUE)

c <- eigen_centrality(g, directed = F)
nodes10$eigen <- c$vector
nodes10$age2 <- factor(nodes10$age,
                       levels = c('40+','25-40',
                                  '20-25','15-19',
                                  '10-15'))

ggplot(data = nodes10, aes(x = age, y = eigen,
                           fill = age2))+
  geom_boxplot()+
  geom_jitter(width = 0.2, shape = 1,
              aes(size = count))+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_fill_viridis_d()+
  scale_x_discrete(name = 'age category')+
  scale_y_continuous(name = 'eigenvector centrality')

ggplot(data = nodes10, aes(x = count, y = eigen))+
  geom_point(aes(colour = age2))+
  geom_smooth()+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18))+
  scale_color_viridis_d()+
  scale_x_continuous(name = 'count')+
  scale_y_continuous(name = 'eigenvector centrality')

# try using sna function to measure and check for any similarity
sna_cent <- evcent(dat = adj_mat, gmode = 'graph', ignore.eval = F)
nodes$sna_cent <- sna_cent

# plot 3 different centrality measures to see if there is any agreement at all -- no there's not
plot(nodes$model, col = 'red', pch = 19, ylim = c(0,1))
points(nodes$igraph, col = 'blue', pch = 4)
points(nodes$sna_cent, col = 'purple', pch = 2)

# test is sna output values match to either once you convert to relative values not absolute
max_sna <- max(sna_cent)
sna_relative <- sna_cent/max_sna
points(sna_relative, col = 'black', pch = 2) # YES -- SNA PACKAGE AND IGRAPH AGREE ON VALUES, JUST THAT ONE IS ABSOLUTE AND ONE IS RELATIVE. ALSO IN SAME ORDER SO AS LONG AS IGRAPH ONES ARE IN RIGHT ORDER IN NODES, SO ARE SNA_CENT

ggplot(nodes, aes(x = age_cat, y = sna_cent,
                  fill = factor(age_cat, levels = c('40+','25-40','20-25','15-19','10-15'))))+
  geom_boxplot(notch = T)+
  geom_jitter(width = 0.2, aes(size = count), colour = rgb(0,0,0,0.5))+
  theme_classic()+
  theme(legend.position = 'none', axis.title = element_text(size = 18), axis.text = element_text(size = 14))+
  scale_fill_viridis_d()+
  scale_x_discrete(name = 'age category')+
  scale_y_continuous(name = 'eigenvector centrality')

write_csv(nodes, '../data_processed/motnp_eigenvector_estimates.csv')

# posterior check v3 -- take draws from posterior and calculate values from that ####
#library(bisonR) ; library(brms) ; library(tidyverse) ; library(rethinking)
#load('motnp_nodalregression_meanage.RData')
##rm(biologylibs, environlibs, homedrive, homelibs, homelibsprofile, mathlibs, psychlibs, rlibs, Rversion)

# extract mean of posterior for eigenvector
mean_eigen_values <- extract_metric(motnp_edge_weights_strongpriors, "node_eigen") %>%
  as.data.frame() %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = everything()) %>%
  rename(bison_node_eigen=value) %>%
  mutate(node = df_nodal$node,
         id = df_nodal$id)

# extract posterior draws for eigenvector
eigen_values <- extract_metric(motnp_edge_weights_strongpriors, "node_eigen") %>%
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'name', values_to = 'eigen') %>% 
  left_join(mean_eigen_values, by = 'name') %>% 
  rename(mean_eigen = bison_node_eigen) %>% 
  select(-name)

# combine with mean age values to create data frame for plotting scatter plots
eigen_values <- left_join(eigen_values, mean_motnp_ages, by = c('node','id'))

# create categorical age variable for plotting boxplots
eigen_values$age_cat <- ifelse(eigen_values$age < 15, '10-15',
                               ifelse(eigen_values$age < 20, '15-20',
                                      ifelse(eigen_values$age < 25, '20-25',
                                             ifelse(eigen_values$age < 40, '25-40', '40+'))))
rm(mean_eigen_values, mean_motnp_ages, df_nodal, counts_df, counts_df_model, motnp_ages, motnp_edges_null_strongpriors, priors, age, beta, beta_mu, beta_sigma, i, intercept, mean_age) ; gc()

# plot eigenvector distributions
ggplot(eigen_values)+
  geom_density(aes(x = eigen, group = node), colour = rgb(0,0,1,0.4))+
  theme_classic()

# plot eigenvector vs age
#ggplot(data = eigen_values#, aes(x = jitter(age))
#       )+
#  geom_point(aes(x = age, y = mean_eigen))+
#  geom_violin(aes(x = age, y = eigen, group = node), fill = rgb(0,0,1,0.2))+
#  theme_classic()+
#  theme(legend.position = 'none')+
#  scale_x_continuous(name = 'mean age')+
#  scale_y_continuous(name = 'eigenvector centrality')

# plot eigenvector vs age category
#library(ggridges)
#ggplot(eigen_values, aes(x = eigen, group = node,
#                         y = age_cat, fill = age_cat)) +
#  geom_density_ridges() +
#  theme_ridges() + 
#  scale_fill_viridis_d(alpha = 0.4)+
#  theme(legend.position = "none") +
#  xlab("eigenvector centrality") +
#  ylab("age category")
#ggplot(eigen_values, aes(x = eigen, y = factor(age),
#                         group = node,
#                         fill = age_cat)) +
#  geom_density_ridges(scale = 2) +
#  theme_ridges() + 
#  scale_fill_viridis_d()+
#  theme(legend.position = "none",
#        axis.text.x = element_blank()) +
#  xlab("eigenvector centrality") +
#  ylab("age category")+
#  coord_flip()

# extract posterior intercept and mean slope
post_eigen <- as.data.frame(as_draws_df(mean_motnp_eigen)) %>% janitor::clean_names()

# compute posterior means
age <- seq(10,55,1)
mu <- matrix(nrow = nrow(post_eigen), ncol = length(age))
for(i in 1:nrow(mu)){
  for(j in 1:ncol(mu)){
    mu[i,j] <- post_eigen$b_intercept[i] + post_eigen$b_age[i]*age[j]
  }
}
mu_mean <- apply(mu, 2, mean)
mu_hpdi <- apply(mu, 2, HPDI, prob = 0.95)

# compute posterior predictions -- I think this does the same as my loop above but not totally sure. The outputs are similar, but not the same.
predict <- as.data.frame(posterior_predict(mean_motnp_eigen,
                                           newdata = as.data.frame(age)))
predict_mean <- apply(predict, 2, mean)
predict_hpdi <- apply(predict, 2, HPDI, prob = 0.95)

# simulate from posterior
sim_eigen <- matrix(ncol = length(age), nrow = 1000)
for(i in 1:nrow(sim_eigen)){
  for(j in 1:ncol(sim_eigen)){
    a <- post_eigen$b_intercept[i] + post_eigen$b_age[i]*age[j]
    b <- post_eigen$sigma[i]
    #sim_eigen[i,j] <- rbeta(1, shape1 = a/(a+b), shape2 = 1 - (a/(a+b))) # I THINK this one is wrong as this would be shape 1 = mean of beta distribution and shape 2 = spread, where instead I want the actual shape values, but I've left it here just in case I'm wrong and need to change it back! produce fairly similar values
    sim_eigen[i,j] <- rbeta(1, shape1 = a, shape2 = b)
  }
}
hist(sim_eigen) # far too many getting a centrality score of 1
sim_pi <- apply(sim_eigen, 2, HPDI, prob = 0.95)

# plot
plot(eigen ~ age, data = eigen_values, col = rgb(0,0,1,0.01), # raw sightings
     pch = 19, las = 1, xlim = c(10,55), ylim = c(0,1),
     xlab = 'mean age', ylab = 'eigenvector centrality')
lines(age, mu_mean, lwd = 2, col = 'red')                     # add mean line
shade(mu_hpdi, age)                                           # add mean uncertainty
shade(sim_pi, age)                                            # add predictions

mean_line <- data.frame(age = age, mean = mu_mean)
mean_shade <- data.frame(age = age, lb = mu_hpdi[1,], ub = mu_hpdi[2,])
pred_shade <- data.frame(age = age, lb = sim_pi[1,], ub = sim_pi[2,])

ggplot()+
  geom_ribbon(data = mean_shade, aes(x = age, ymin = lb, ymax = ub),
              fill = rgb(0,0,0,0.15))+
  geom_ribbon(data = pred_shade, aes(x = age, ymin = lb, ymax = ub),
              fill = rgb(0,0,0,0.15))+
  geom_point(data = eigen_values, aes(x = age, y = eigen,
                                      colour = as.factor(age)),
             pch = 19, alpha = 0.01)+
  geom_line(data = mean_line, aes(x = age, y = mean),
            linewidth = 1.5, colour = 'black')+
  theme_classic()+
  theme(legend.position = 'none',
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))+
  scale_colour_viridis_d(direction = -1)+
  scale_x_continuous(name = 'mean of age distribution')+
  scale_y_continuous(name = 'eigenvector centrality\n(1000 draws/elephant)')

ggplot()+
  geom_ribbon(data = mean_shade, aes(x = age, ymin = lb, ymax = ub),
              fill = rgb(0,0,0,0.15))+
  #geom_ribbon(data = pred_shade, aes(x = age, ymin = lb, ymax = ub),
  #            fill = rgb(0,0,0,0.15))+
  geom_point(data = eigen_values, aes(x = age, y = eigen,
                                      colour = as.factor(age)),
             pch = 19, alpha = 0.01, size = 2)+
  geom_line(data = mean_line, aes(x = age, y = mean),
            linewidth = 2.5, colour = 'black')+
  theme_classic()+
  theme(legend.position = 'none',
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 24))+
  scale_colour_viridis_d(direction = -1)+
  scale_x_continuous(name = 'mean of age distribution')+
  scale_y_continuous(name = 'eigenvector centrality\n(1000 draws/elephant)')

## save workspace image
save.image('motnp_nodalregression_posteriorplots.RData')

### run model -- full age distribution ####
#pdf('../outputs/motnp_nodalregression_eigen_agedistribution.pdf')

## define priors
motnp_edge_weights_strongpriors$model_data$prior_fixed_sigma <- 0.005

## reduce age data to manageable number of draws
motnp_ages$draw <- rep(1:8000, each = length(unique(motnp_ages$id)))
test_ages <- motnp_ages[motnp_ages$draw %in% sample(motnp_ages$draw, 5, replace = F),]

## run model
x <- bison_mice(edgemodel_list = list(motnp_edge_weights_strongpriors),
                data_list = test_ages,
                param_names = unique(c(counts_df_model$node_1_id, counts_df_model$node_2_id)),
                target_name = 'node',
                metric_name = 'eigen',
                num_draws = 10,
                z_score = F)

edgemodel_list <- list(motnp_edge_weights_strongpriors)
data_list_long <- test_ages %>% 
  rename(id_1 = id) %>% 
  left_join(distinct(counts_df[,c('node_1_males','id_1')]), by = 'id_1') %>% 
  rename(id_2 = id_1) %>% 
  left_join(distinct(counts_df[nrow(counts_df),c('node_2_males','id_2')]), by = 'id_2') %>% 
  mutate(id = ifelse(is.na(id_2), id_1, id_2),
         node = ifelse(is.na(node_2_males), node_1_males, node_2_males)) %>% 
  filter(is.na(node) == F) %>% 
  select(node, age)
#data_list <- pivot_wider(data = data_list_long,
#                         names_from = 'node', 
#                         values_from = 'age') %>% 
#  unnest(cols = everything())

data_list <- data.frame()

param_names <- node  #unique(c(counts_df$node_1_males, counts_df$node_2_males))
target_name <- 'node'
metric_name <- 'eigen'
num_draws <- 10
z_score <- F

motnp_ages_list <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>% 
  dplyr::select(sort(unique(c(counts_df$id_1, counts_df$id_2)))) %>% 
  as.list()
for(i in 1:length(motnp_ages_list)){
  id <- names(motnp_ages_list[i])
  motnp_ages_list[[i]] <- as.data.frame(motnp_ages_list[[i]])
  colnames(id)
}

#mean_motnp_eigen <- bison_mice(
  edgemodel_list = list(motnp_edge_weights_strongpriors)#,
  #data_list = as.data.frame(motnp_ages_list)#,
  #data_list = list(df_nodal$node, motnp_ages_list)
  data_list = array(data = NA, dim = c(213, 8000))#,
                    #dimnames = list(df_nodal$id, 'age_value',1:8000))
for(i in 1:nrow(motnp_ages)){
  for(j in 1:ncol(motnp_ages)){
  data_list[i,,j] <- motnp_ages[,i]
  }
}
  param_names = df_nodal$node#, c(counts_df_model$id_1,counts_df_model$id_2), 'node_id' ,df_nodal$id , 
  target_name = 'node'#,
  metric_name = 'eigen'#,
  num_draws = 10#,
  z_score = T
#)

#### bison_mice code ####
if (class(edgemodel_list)[1] != "list") edgemodel_list <- list(edgemodel_list)
#if (class(data_list)[1] != "list") data_list <- list(data_list)

new_bison_term <- paste0("bison_", target_name, "_", metric_name)

posterior_samples_list <- list()

# actual bison code
if (target_name == "node") {
  node_ids <- sapply(
    dplyr::pull(data_list[[i]], param_names[1]),                   # doesn't work because param_names = nodes, data_list[[1]] = ages
    function(x) which(names(edgemodel_list[[i]]$node_to_idx) == x) # produces a list per elephant of 5 zeroes because age doesn't match to node ids
  )
  posterior_samples <- extract_metric(edgemodel_list[[i]], paste0(target_name, "_", metric_name), num_draws)#[, node_ids]
}

# me messing around until something does something without an error
if (target_name == "node") {
  # node_ids <- 1:213
  node_ids <- sapply(
    dplyr::pull(as.data.frame(data_list[[1]]), param_names[1]),
    function(x) which(names(edgemodel_list[[1]]$node_to_idx) == x)
  )
  posterior_samples <- extract_metric(edgemodel_list[[1]], paste0(target_name, "_", metric_name), num_draws)#[, node_ids]
}

# Make sure posterior samples are in matrix format (for num_draws=1).
posterior_samples_list[[i]] <- matrix(posterior_samples, nrow=num_draws)

# Generate list of dataframes.
imputed_dataframes <- lapply(1:num_draws, function(i) {
  # For each posterior draw, combine dataframes from different edge models and data
  new_data_list <- lapply(2,#1:length(edgemodel_list),
                          function(j) {
    new_data <- data_list[[j]]
    if (z_score) {
      new_data[new_bison_term] <- (posterior_samples_list[[j]][i, ] - mean(posterior_samples_list[[j]][i, ]))/sd(posterior_samples_list[[j]][i, ])
    } else {
      new_data[new_bison_term] <- posterior_samples_list[[j]][i, ]
    }
    new_data["bison_network"] <- j
    new_data
  })
  new_data <- dplyr::bind_rows(new_data_list)
  new_data$bison_network <- factor(new_data$bison_network)
  new_data
})

imputed_dataframes <- dplyr::bind_rows(imputed_dataframes, .id=".imp")
imputed_dataframes <- dplyr::mutate(imputed_dataframes, .imp=as.integer(.imp))

original_dataframes <- lapply(1:length(edgemodel_list), function(j) {
  new_data <- data_list[[j]]
  new_data["bison_network"] <- j
  new_data$bison_network <- factor(new_data$bison_network)
  new_data
})
original_dataframes <- dplyr::bind_rows(original_dataframes)
original_dataframes[".imp"] <- 0

combined_dataframes <- dplyr::bind_rows(list(original_dataframes, imputed_dataframes))

mice::as.mids(as.data.frame(combined_dataframes))



### extract eigenvector centralities -- old method ####

## create array for eigen values to be saved into
eigen <- array(data = NA, dim = c(n_eles, 4, n_samples*n_chains),
               dimnames = list(nodes$node,
                               c('node','age','sightings','eigenvector'),
                               1:(n_samples*n_chains)))
eigen[,1,] <- nodes$node
eigen[,2,] <- nodes$age
eigen[,3,] <- nodes$sightings

## fill array
for(draw in 1:(n_samples*n_chains)){
  network <- graph_from_adjacency_matrix(adjmatrix = adj_tensor[,,draw],
                                         diag = FALSE, mode = 'undirected', weighted = TRUE)
  eigen_values <- as.matrix(igraph::eigen_centrality(network, directed = FALSE)$vector)
  eigen[,4,draw] <- eigen_values[,1]
}

## save workspace for future
rm(dyad_row, edge_binary, edgelist, eigen_values, network, adj_tensor, draw, i) ; gc()
save.image('motnp_nodalregression_conditionaledge.RData')

## check eigenvector against sightings
plot(NULL, xlim = c(0,max(nodes$sightings)), ylim = c(0,1),
     las = 1, xlab = 'sightings', ylab = 'eigenvector',
     main = 'time window = 1')
for(i in 1:n_samples){
  points(eigen[,4,i] ~ eigen[,3,i], pch = 19, col = rgb(0.5,0,1,0.1))
}
for(i in 1:n_eles){
  x <- eigen[i,,]
  points(mean(x[4,]) ~ x[3,1], pch = 19, col = 'yellow')
}

## check covariance of eigenvector
plot(eigen[1,4,] ~ eigen[2,4,], col = rgb(0,0,1,0.2), las = 1, pch = 19)
plot(eigen[3,4,] ~ eigen[4,4,], col = rgb(0,0,1,0.2), las = 1, pch = 19)
plot(eigen[5,4,] ~ eigen[6,4,], col = rgb(0,0,1,0.2), las = 1, pch = 19)















