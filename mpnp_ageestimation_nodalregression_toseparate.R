#### information #####
# To estimate the ages of individuals in other populations, first come up with a survival curve for the Amboseli elephants which have good birth and death data. -- this identified the best curve as being the Gompertz-Bathtub.
# Next fit this to the MPNP dataset

############ First half of script: Obtain age estimates ############
#### load packages ####
library(tidyverse)
library(cmdstanr)
library(ggdist)
library(posterior)
library(bayesplot)

#### load model -- adapted from MOTNP to use MPNP thresholds ####
# read in Stan model to estimate ages based on Gompertz bathtub distribution from ANP
latent_age_ordinal_model <- cmdstan_model("models/age_estimation/mpnp_elephant_latent_age_ordinal_regression_hkm_22.07.07.stan")
# Age_Range_ID
# 1	>1
# 2	1-4
# 3	5-9
# 4	10-15
# 5	16-20
# 6	21-25
# 7	26-35
# 8	36+
#10	UK

#### load MPNP data ####
#mpnp1_males <- read_delim('data_processed/mpnp_elenodes_22.03.08.csv', delim = ' ') %>% 
#  filter(sex == 'M')
#unique(mpnp1_males$age_class) # AGE CLASS RECORDED HERE IS MODAL CATEGORY ACROSS ALL SIGHTINGS PER INDIVIDUAL
#rm(mpnp1_males)
#mpnp_groups <- readxl::read_excel('data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211214.xlsx')
#colnames(mpnp_groups) <- mpnp_groups[2,]
#mpnp_groups <- mpnp_groups[3:nrow(mpnp_groups),c(1:23,57)] %>% janitor::clean_names()
mpnp_long <- readxl::read_excel('data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211019.xlsx') %>%
  janitor::clean_names()

# identify period 1 only
periods <- seq(from = min(mpnp_long$date, na.rm = T),
               to = max(mpnp_long$date, na.rm = T),
               length.out = 7)
periods[7] <- periods[7]+1 # set to be one higher than the final higher otherwise it takes the last date and creates a whole new period
#as.numeric(max(mpnp_groups$date, na.rm = T)) - as.numeric(min(mpnp_groups$date, na.rm = T))
#mpnp1_groups <- mpnp_groups[mpnp_groups$date < periods[2],]

# remove non-B-numbered elephants
mpnp_long <- mpnp_long %>% 
  filter(elephant_id != '-') %>% 
  separate(elephant_id, into = c('BTF','num'), sep = 1, remove = F) %>% 
  filter(BTF == 'B') %>% 
  select(-BTF, -num)

# select only period 1
mpnp1_long <- mpnp_long[mpnp_long$date < periods[2],]
max(mpnp_long$date, na.rm = T) - min(mpnp_long$date, na.rm = T)

# get range of ages
mpnp1_long$age_range_id <- as.numeric(mpnp1_long$age_range_id)
mpnp1_long$age_min <- NA
mpnp1_long$age_max <- NA
mpnp1_long$age_mid <- NA
mpnp1_long$age_mean <- NA
mpnp1_long$age_range <- NA
for(i in 1:nrow(mpnp1_long)){
  ele <- mpnp1_long[mpnp1_long$elephant_id == mpnp1_long$elephant_id[i],]
  ele <- ele[!is.na(ele$age_range_id),]
  ele <- ele[ele$age_range_id != 10,]
  mpnp1_long$age_min[i] <- ifelse(nrow(ele) == 0, 10, min(ele$age_range_id))
  mpnp1_long$age_max[i] <- ifelse(nrow(ele) == 0, 10, max(ele$age_range_id, na.rm = T))
  mpnp1_long$age_mid[i] <- ifelse(nrow(ele) == 0, 10, median(ele$age_range_id, na.rm = T))
  mpnp1_long$age_mean[i] <- ifelse(nrow(ele) == 0, 10, round(mean(ele$age_range_id, na.rm = T),0))
  mpnp1_long$age_range[i] <- ifelse(nrow(ele) == 0, 0, max(ele$age_range_id, na.rm = T) - min(ele$age_range_id, na.rm = T))
  rm(ele)
}

summary(mpnp1_long$age_range)
to.check <- mpnp1_long[mpnp1_long$age_range > 0,]

mpnp1_long$age_average <- ifelse(mpnp1_long$age_mid - mpnp1_long$age_mean == 0, mpnp1_long$age_mid, 999)
to.check <- mpnp1_long[mpnp1_long$age_average > 10,c(3,7,32:37)]
length(unique(to.check$elephant_id)) # 32

# DISCUSS WITH DAN AND COLIN HOW BEST TO USE THESE












#### create data list ####
N_mpnp1 <- nrow(mpnp1_males)
K <- 9
mpnp1_ls <- list(
  N = N_mpnp1,
  K = K,
  age_category_index = mpnp1_males$age_cat_id)
hist(mpnp1_ls$age_category_index)
#hist(elephants_ls$age_category_index)

#### fit model to MPNP data ####
# Fit model with cmdstanr
age_mpnp1_fit <- latent_age_ordinal_model$sample(
  data = mpnp1_ls, 
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 2000
)

# Examine the estimates. We can plot the estimated ages against the biologist assigned ages
age_est_mat <- age_mpnp1_fit$summary()[(N_mpnp1+2):(N_mpnp1*2+1), ]
summary(age_est_mat)
hist(age_est_mat$mean)
hist(age_est_mat$rhat, breaks = 20)

plot_data <- data.frame(age = ifelse(mpnp1_ls$age == 1, 1, 
                                     ifelse(mpnp1_ls$age == 2, 3,
                                            ifelse(mpnp1_ls$age == 3, 7,
                                                   ifelse(mpnp1_ls$age == 4, 12,
                                                          ifelse(mpnp1_ls$age == 5, 18,
                                                                 ifelse(mpnp1_ls$age == 6, 22, 
                                                                        ifelse(mpnp1_ls$age == 7, 30, 45))))))),
                        model_age = age_est_mat$mean) # Mean modelled age

plot_data %>%
  ggplot(aes(x=factor(age), y=model_age)) +
  geom_point(size=4,col = 'blue', alpha=0.6) +
  geom_vline(xintercept = c(5, 10, 15, 20, 25, 40, 60), linetype = "dashed", alpha = 0.6) +
  geom_hline(yintercept = c(5, 10, 15, 20, 25, 40, 60), linetype = "dashed", alpha = 0.6) +
  #geom_abline(slope = 1, intercept = 0)+
  scale_y_continuous(limits = c(0,60))+
  theme_minimal() + 
  xlab("Assigned age") + ylab("Modelled age")

# posterior predictive plot using draws from distribution to show uncertainty around mean age
true_ages <- age_mpnp1_fit$draws("true_age", format="df")
#mcmc_dens(true_ages)
true_ages <- true_ages[,1:N_mpnp1]

df <- as.data.frame(do.call(rbind, true_ages)) %>%
  mutate(age_cat = mpnp1_ls$age) %>% relocate(age_cat) %>%
  mutate(ID = mpnp1_males$id) %>% relocate(ID)

df <- df %>% pivot_longer(cols = 3:102) %>% select(-name)

df$true_age <- ifelse(df$age_cat == 1, 3, 
                      ifelse(df$age_cat == 2, 8,
                             ifelse(df$age_cat == 3, 12,
                                    ifelse(df$age_cat == 4, 18,
                                           ifelse(df$age_cat == 5, 22, 
                                                  ifelse(df$age_cat == 6, 32, 45))))))

df %>% ggplot(aes(x=true_age, y=value, group=factor(ID))) +
  geom_point(size=2,col = 'blue', alpha=0.1) +
  #stat_halfeye() +
  geom_vline(xintercept=c(5,10,15,20,25,40,60), linetype="dashed", alpha=0.6) +
  geom_hline(yintercept=c(5,10,15,20,25,40,60), linetype="dashed", alpha=0.6) +
  theme_bw() + 
  xlab("Assigned age") + ylab("Modelled age")

############ Second half of script: Use in nodal regression ############
rm(elephants_ls, age_estimation_fit, latent_age_ordinal_model, plot_data)
#### load packages ####
library(rstan)
library(igraph)
library(tidyverse)
library(LaplacesDemon)

#### read in MPNP data ########
### clean environment ###
rm(df_long_mpnp1, df_wide_mpnp1, )

### import data for aggregated model (binomial)
df_agg_mpnp1 <- read_delim('data_processed/mpnp1_bayesian_allpairwiseevents_splitbygrouptype_22.01.13.csv', delim = ',') %>% 
  filter(dem_class_1 == 'AM' | dem_class_1 == 'PM') %>% 
  filter(dem_class_2 == 'AM' | dem_class_2 == 'PM')
df_agg_mpnp1$sex_1 <- 'M'
df_agg_mpnp1$age_cat_id_1 <- ifelse(df_agg_mpnp1$age_category_1 == '9-10', 2,
                                    ifelse(df_agg_mpnp1$age_category_1 == '10-15', 3,
                                           ifelse(df_agg_mpnp1$age_category_1 == '15-19', 4,
                                                  ifelse(df_agg_mpnp1$age_category_1 == '20-25', 5,
                                                         ifelse(df_agg_mpnp1$age_category_1 == '25-40', 6,
                                                                ifelse(df_agg_mpnp1$age_category_1 == '40+', 7,
                                                                       df_agg_mpnp1$age_category_1))))))
df_agg_mpnp1$age_class_1 <- ifelse(df_agg_mpnp1$age_cat_id_1 == 2, 'Juvenile',
                                   ifelse(df_agg_mpnp1$age_cat_id_1 > 4, 'Adult','Pubescent'))
df_agg_mpnp1$age_cat_id_2 <- ifelse(df_agg_mpnp1$age_category_2 == '9-10', 2,
                                    ifelse(df_agg_mpnp1$age_category_2 == '10-15', 3,
                                           ifelse(df_agg_mpnp1$age_category_2 == '15-19', 4,
                                                  ifelse(df_agg_mpnp1$age_category_2 == '20-25', 5,
                                                         ifelse(df_agg_mpnp1$age_category_2 == '25-40', 6,
                                                                ifelse(df_agg_mpnp1$age_category_2 == '40+', 7,
                                                                       df_agg_mpnp1$age_category_2))))))
df_agg_mpnp1$age_class_2 <- ifelse(df_agg_mpnp1$age_cat_id_2 == 2, 'Juvenile',
                                   ifelse(df_agg_mpnp1$age_cat_id_2 > 4, 'Adult','Pubescent'))
df_agg_mpnp1$dem_class_1 <- ifelse(df_agg_mpnp1$age_class_1 == 'Adult', 'AM',
                                   ifelse(df_agg_mpnp1$age_class_1 == 'Pubescent', 'PM', 'JM'))
df_agg_mpnp1$dem_class_2 <- ifelse(df_agg_mpnp1$age_class_2 == 'Adult', 'AM',
                                   ifelse(df_agg_mpnp1$age_class_2 == 'Pubescent', 'PM', 'JM'))
df_agg_mpnp1$dem_type <- ifelse(df_agg_mpnp1$age_cat_id_1 >= df_agg_mpnp1$age_cat_id_2,
                                paste(df_agg_mpnp1$dem_class_1, df_agg_mpnp1$dem_class_2, sep = '_'),
                                paste(df_agg_mpnp1$dem_class_2, df_agg_mpnp1$dem_class_1, sep = '_'))
df_agg_mpnp1$count_dyad <- (df_agg_mpnp1$count_1 + df_agg_mpnp1$count_2) - df_agg_mpnp1$all_events  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.
df_agg_mpnp1$node_1_nogaps <- as.integer(as.factor(df_agg_mpnp1$node_1))
df_agg_mpnp1$node_2_nogaps <- as.integer(as.factor(df_agg_mpnp1$node_2))+1
df_agg_mpnp1$dyad_id_nogaps <- as.integer(as.factor(df_agg_mpnp1$dyad))

### load the edge weights
mpnp1 <- readRDS('data_processed/mpnp1_bayesian_edgedistributions_a2.b2_22.02.07.rds') %>% 
  select(-`1.lp__`)
mpnp1 <- mpnp1[, which(colnames(mpnp1) %in% df_agg_mpnp1$dyad)]
logit_edge_samples_mpnp1 <- logit(mpnp1)

#### calculate posterior centralities ####
mpnp1 <- as.matrix(logit_edge_samples_mpnp1)
edge_samples_mpnp1 <- plogis(mpnp1)

# Build adjacency tensor
num_nodes <- length(unique(c(df_agg_mpnp1$id_1, df_agg_mpnp1$id_2)))
num_samples <- 4000
adj_tensor <- array(0, c(num_samples, num_nodes, num_nodes),
                    dimnames = list(NULL,
                                    unique(c(df_agg_mpnp1$id_1, df_agg_mpnp1$id_2)),
                                    unique(c(df_agg_mpnp1$id_1, df_agg_mpnp1$id_2))))
for (dyad_id in 1:nrow(df_agg_mpnp1)) {
  dyad_row <- df_agg_mpnp1[df_agg_mpnp1$dyad_id_nogaps == dyad_id, ]
  adj_tensor[, dyad_row$id_1, dyad_row$id_2] <- edge_samples_mpnp1[, dyad_id]
}
rm(dyad_row)

# Calculate centrality and store posterior samples in a matrix
eigen_samples_mpnp1 <- matrix(0, num_samples, num_nodes)
eigen_samples_mpnp1_std <- matrix(0, num_samples, num_nodes)
for (i in 1:num_samples) {
  g <- graph_from_adjacency_matrix(adj_tensor[i, , ], mode="undirected", weighted=TRUE)
  eigen_samples_mpnp1[i, ] <- eigen_centrality(g)$vector
  eigen_samples_mpnp1_std[i, ] <- (eigen_samples_mpnp1[i, ] - mean(eigen_samples_mpnp1[i, ]))/sd(eigen_samples_mpnp1[i, ])
}
head(eigen_samples_mpnp1)     # Unstandardised eigenvector centrality
head(eigen_samples_mpnp1_std) # Standardised eigenvector centrality
rm(g)

# Computing these matrices is essential for model fitting, but it’s hard to understand differences in centralities intuitively from these raw values. Instead we can visualise the centralities of each node using the following code, using tidyverse as it simplifies the process:

df_wide_mpnp1 <- data.frame(eigen_samples_mpnp1_std)
colnames(df_wide_mpnp1) <- 1:num_nodes
df_long_mpnp1 <- pivot_longer(df_wide_mpnp1, cols = 1:num_nodes, names_to="node_id", values_to="eigenvector_centrality")
ggplot(df_long_mpnp1, aes(x = eigenvector_centrality)) +
  geom_density(fill="#387780", alpha=0.7, size=0.4) +
  #facet_grid(rows=vars(as.factor(node_id)), scales="free") +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        #strip.text.y = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

nodes1 <- df_agg_mpnp1 %>% select(id_1, age_class_1, age_cat_id_1) %>% distinct()
nodes2 <- df_agg_mpnp1 %>% select(id_2, age_class_2, age_cat_id_2) %>% distinct()
colnames(nodes1) <- c('id', 'age_class', 'age_cat_id')
colnames(nodes2) <- c('id', 'age_class', 'age_cat_id')
node_ages_mpnp1 <- rbind(nodes1, nodes2) %>% distinct()
rm(nodes1, nodes2)

node_ages_mpnp1$node_id <- as.integer(as.factor(node_ages_mpnp1$id))
df_long_mpnp1$node_id <- as.integer(df_long_mpnp1$node_id)
df_long_mpnp1 <- left_join(df_long_mpnp1, node_ages_mpnp1, by = 'node_id')

ggplot(plot_eles, aes(x = eigenvector_centrality)) +
  geom_density(aes(fill = AGE), alpha=0.7, size=0.4) +
  facet_grid(rows=vars(as.factor(node_id)), scales="free") +
  labs(x="Eigenvector centrality (standardised)") + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 0, size=12, debug = FALSE),
        axis.title.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

#### compute normal approximation ####
# The posterior centralities can now be characterised by the multivariate normal distribution as an approximation to their true posterior distributions. To highlight the importance of sampling from a multivariate normal rather than a univariate normal, we can plot centrality samples from two nodes against each other to see how they co-vary.
plot(eigen_samples_mpnp1_std[, 1], eigen_samples_mpnp1_std[, 2])

# This joint distribution contains important information that we would ideally like to include in the model. The multivariate normal will allow us to do this. Parameterising the multivariate normal approximation is a relatively simple process and can be done by calculating the sample mean and sample covariance matrix of the posterior centrality samples:
eigen_mu_mpnp1 <- apply(eigen_samples_mpnp1_std, 2, mean) 
eigen_cov_mpnp1 <- cov(eigen_samples_mpnp1_std)

# These quantities will be given to the Stan model as data to model joint posteriors of centrality in the regression. We can run a few quick plots to see how well the approximation is working and the covariance for one of the nodes:
eigen_samples_mpnp1_sim <- MASS::mvrnorm(1e5, eigen_mu_mpnp1, eigen_cov_mpnp1)
plot(density(eigen_samples_mpnp1_std[, 1]), lwd=2, main="Estimated standardised centrality vs normal approximation", xlab="Logit edge weight")
lines(density(eigen_samples_mpnp1_sim[, 1]), col=rgb(0, 0, 1, 0.5), lwd=2)

# If we’re happy with the approximation, we can now get ready to fit the model.

#### define the model ####
model_nodal_cont <- stan_model("models/nodal_regression_hkm_continuousage_22.07.06.stan")

#### prior predictive check -- WORK OUT HOW TO DO THIS ####













#### fit the model ####
colnames(true_ages) <- mpnp1_males$id # check this shouldn't be as.integer(as.factor(mpnp1_males$id)), but 90% sure this is correct -- goes into the model in order so should come out in same order?
mpnp1_ap <- mpnp1_males[mpnp1_males$id %in% node_ages_mpnp1$id, ]
true_ages_ap <- true_ages[, colnames(true_ages) %in% node_ages_mpnp1$id]

node_ages_mpnp1$age_mean <- NA
for(i in 1:nrow(node_ages_mpnp1)){
  age_ele <- unlist(true_ages_ap[,i])
  node_ages_mpnp1$age_mean[i] <- mean(age_ele, na.rm = T)
}
hist(node_ages_mpnp1$age_mean)

model_data_mpnp1 <- list(
  num_nodes = num_nodes,                          # Number of dyads
  centrality_mu  = eigen_mu_mpnp1,                # Sample means of logit edge weights
  centrality_cov = eigen_cov_mpnp1,               # Sample covariance of logit edge weights
  node_age = node_ages_mpnp1$age_mean,            # Age of individual -- these qre in the wrong order? should be M1, M10, M100, M101...?
  node_age2 = (node_ages_mpnp1$age_mean)^2        # Age squared for quadratic term
)
str(model_data_mpnp1)

plot(model_data_mpnp1$centrality_mu ~ model_data_mpnp1$node_age, las = 1, pch = 19, col = rgb(0,0,1,0.3))

fit_nodal_cont <- sampling(model_nodal_cont, data = model_data_mpnp1, cores = 4, chains = 4)

#### diagnostics ####
traceplot(fit_nodal_cont)

### posterior predictive checks -- I think this is probably a sign that I've written the model wrong because I can't generate an mvnorm from this output -> params$beta_age is 1-dimensional so can't select j'th network and also ages
params <- rstan::extract(fit_nodal_cont)
plot(density(eigen_samples_mpnp1_std[1, ]), main = "Posterior predictive density (standardised centrality)",
     col = rgb(0, 0, 0, 0.25), ylim = c(0, 1), las = 1)
for (i in 1:100) {
  j <- sample(1:num_samples, 1)                                        # select a network to plot
  lines(density(eigen_samples_mpnp1_std[j, ]), col=rgb(0, 0, 0, 0.25)) # plot centrality density for network j (black)
  mu <- params$beta_age[j]#, model_data_mpnp1$node_age]                  # extract age slope parameter for network j
  sigma <- eigen_cov_mpnp1 + diag(rep(params$sigma[j], num_nodes))     # extract variance parameter for network j
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))  # plot predicted centralities for network j (blue)
}

#### interpreting the model ####
# Calculate the 95% credible intervals of the model parameters
round(summary(fit_nodal_cont)$summary[1:3, c(1:4, 8)], 2) # basically no effect whatsoever
