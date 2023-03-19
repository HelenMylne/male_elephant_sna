#### Information ####
# 1) to estimate the ages of individuals in other populations, first came up with a survival curve for the Amboseli elephants which have good birth and death data.
# 2) selected the Gompertz bathtub distribution
# 3) use distribution to convert MOTNP ages in categories to probability distributions of true age (first part of this script)
# 4) extract node centrality distribution from draws of edge weight
# 5) use MOTNP age distributions to predict network centrality (third part of this script)

#### set up ####
options(future.globals.maxSize = 10000*(1024^2))   # running model with full age distributions = 5.01GB of globals, which exceeds default maximum allowed size. Set to allow up to 10 GB to allow model to run

library(tidyverse, lib.loc = '../packages/')       # library(tidyverse)
library(cmdstanr, lib.loc = '../packages/')        # library(cmdstanr)
library(brms, lib.loc = '../packages/')            # library(brms)
library(Rcpp, lib.loc = '../packages/')            # library(Rcpp)
library(ggdist, lib.loc = '../packages/')          # library(ggdist)
library(posterior, lib.loc = '../packages/')       # library(posterior)
library(bayesplot, lib.loc = '../packages/')       # library(bayesplot)
#library(rstan, lib.loc = '../packages/')           # library(rstan)
library(igraph, lib.loc = '../packages/')          # library(igraph)
library(LaplacesDemon, lib.loc = '../packages/')   # library(LaplacesDemon)
library(bisonR, lib.loc = '../packages/')          # library(bisonR)
library(janitor, lib.loc = '../packages/')         # library(janitor)

# load edge weight model and data frames
load('motnp_bisonr_edgescalculated_strongprior.RData')
#rm(counts_df_model, edgelist, females_df, motnp_edges_null_strongpriors) ; gc()

#### read in data ####
df_nodal <- distinct(counts_df[,c('node_1_males','id_1')])
colnames(df_nodal) <- c('node_2_males','id_2')
df_nodal <- rbind(df_nodal, counts_df[nrow(counts_df),c('node_2_males','id_2')])
colnames(df_nodal) <- c('node','id')

motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>% 
  select(sort(unique(c(counts_df$id_1, counts_df$id_2)))) %>% 
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')
motnp_ages <- left_join(motnp_ages, df_nodal, by = 'id')
motnp_ages$draw <- rep(1:8000, length(unique(motnp_ages$id)))

mean_motnp_ages <- df_nodal
mean_motnp_ages$age <- NA
for(i in 1:nrow(mean_motnp_ages)){
  x <- motnp_ages[motnp_ages$id == mean_motnp_ages$id[i],]
  mean_motnp_ages$age[i] <- mean(x$age)
  rm(x)
}

# define PDF output
pdf('../outputs/motnp_nodalregression_plots_meanage.pdf')

#### set priors ####
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

#### run model ####
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

#### posterior check v1 -- some is good but forgot to draw from posterior so most is rubbish ####
#load('motnp_nodalregression_meanage.RData')
#rm(df_nodal, motnp_ages) ; gc()

summary(mean_motnp_eigen)

# plot
mean_eigen_values <- mean_motnp_eigen$data
plot(mean_eigen_values$bison_node_eigen ~ mean_eigen_values$age, 
     las = 1, pch = 19, col = rgb(0,0,1,0.2),
     xlab = 'mean age estimate', ylab = 'eigenvector centrality',
     main = 'effect of age on eigenvector centrality')

mean_eigen_summary <- mean_motnp_eigen$fit

hist(mean_motnp_eigen$rhats[,2], las = 1, main = 'Rhat values for 100 imputed model runs', xlab = 'Rhat')

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

hist(post_eigen$b_age)         # natural scale?
hist(plogis(post_eigen$b_age)) # logit scale?

#plot(mean_motnp_eigen)

plot(conditional_effects(mean_motnp_eigen), points = TRUE,
     #ylab = 'mean eigenvector centrality',
     theme = theme_classic(base_size = 14),
     ) # older elephants have lower network centrality than younger

mean_eigen_values$age_cat <- ifelse(mean_eigen_values$age < 15, '10-15',
                                    ifelse(mean_eigen_values$age < 20, '16-20',
                                           ifelse(mean_eigen_values$age < 25, '21-25',
                                                  ifelse(mean_eigen_values$age < 40, '25-40', '40+'))))
mean_eigen_values$age_cat <- factor(mean_eigen_values$age_cat,
                                    levels = c('10-15','16-20','21-25','25-40','40+'))
mean_eigen_values$age_cat2 <- factor(mean_eigen_values$age_cat,
                                    levels = c('40+','25-40','21-25','16-20','10-15'))

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

ggplot(data = data, aes(x = age_cat, y = bison_node_eigen,
                        fill = age_cat2))+
  geom_boxplot(notch = T)+
  geom_jitter(width = 0.2, shape = 1,
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

#### posterior check v2 -- run using SRI and see how it compares ####
library(igraph) ; library(tidyverse) ; library(sna) ; library(brms)
load('motnp_bisonr_edgescalculated_strongprior.RData')
load('motnp_nodalregression_meanage.RData')
rm(df_nodal, mean_eigen_summary, priors, age, beta, beta_mu, beta_sigma, i, intercept, mean_age) ; gc()

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

#### posterior check v3 -- take draws from posterior and calculate values from that ####
library(bisonR) ; library(brms) ; library(tidyverse) ; library(rethinking)
load('motnp_nodalregression_meanage.RData')
#rm(biologylibs, environlibs, homedrive, homelibs, homelibsprofile, mathlibs, psychlibs, rlibs, Rversion)

# extract posterior draws for eigenvector
mean_eigen_values <- extract_metric(motnp_edge_weights_strongpriors, "node_eigen") %>%
  as.data.frame() %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = everything()) %>%
  rename(bison_node_eigen=value) %>%
  mutate(node = df_nodal$node,
         id = df_nodal$id)
eigen_values <- extract_metric(motnp_edge_weights_strongpriors, "node_eigen") %>%
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'name', values_to = 'eigen') %>% 
  left_join(mean_eigen_values, by = 'name') %>% 
  rename(mean_eigen = bison_node_eigen) %>% 
  select(-name)
eigen_values <- left_join(eigen_values, mean_motnp_ages, by = c('node','id'))
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
ggplot(data = eigen_values#, aes(x = jitter(age))
       )+
  geom_point(aes(x = age, y = mean_eigen))+
  geom_violin(aes(x = age, y = eigen, group = node), fill = rgb(0,0,1,0.2))+
  theme_classic()+
  theme(legend.position = 'none')+
  scale_x_continuous(name = 'mean age')+
  scale_y_continuous(name = 'eigenvector centrality')

# plot eigenvector vs age category
library(ggridges)
ggplot(eigen_values, aes(x = eigen, group = node,
                         y = age_cat, fill = age_cat)) +
  geom_density_ridges() +
  theme_ridges() + 
  scale_fill_viridis_d(alpha = 0.4)+
  theme(legend.position = "none") +
  xlab("eigenvector centrality") +
  ylab("age category")
ggplot(eigen_values, aes(x = eigen, y = factor(age),
                         group = node,
                         fill = age_cat)) +
  geom_density_ridges(scale = 2) +
  theme_ridges() + 
  scale_fill_viridis_d()+
  theme(legend.position = "none",
        axis.text.x = element_blank()) +
  xlab("eigenvector centrality") +
  ylab("age category")+
  coord_flip()

# extract posterior intercept and mean slope
post_eigen <- as.data.frame(as_draws_df(mean_motnp_eigen)) %>% janitor::clean_names()
#age_mean <- mean(post_eigen$b_age)
#age_pi <- rethinking::PI(post_eigen$b_age)
#intercept_mean <- mean(post_eigen$b_intercept)
#rethinking::PI(post_eigen$b_intercept)

# compute posterior predictions
#pp <- as.data.frame(posterior_predict(mean_motnp_eigen))
#colnames(pp) <- eigen_values$id[1:213]
#pp <- pivot_longer(pp, everything(), names_to = 'id', values_to = 'eigen_predict')
#rethinking::HPDI(pp$eigen_predict, prob = 0.95)

# plot
#raw_eigen <- extract_metric(motnp_edge_weights_strongpriors, "node_eigen") %>%
#  as.data.frame()
age <- seq(10,55,1)
mu <- matrix(nrow = nrow(post_eigen)/100, ncol = length(age))
for(i in 1:nrow(mu)){
  for(j in 1:ncol(mu)){
    mu[i,j] <- post_eigen$b_intercept[i] + post_eigen$b_age[i]*age[j]
  }
}
mu_mean <- apply(mu, 2, mean)
mu_hpdi <- apply(mu, 2, HPDI, prob = 0.95)

#predict <- as.data.frame(posterior_predict(mean_motnp_eigen))
#colnames(predict) <- eigen_values$id[1:213]
#pp_hpdi <- apply(predict, 2, HPDI, prob = 0.95)

sim_eigen <- matrix(ncol = length(age), nrow = 1000)
for(i in 1:nrow(sim_eigen)){
  for(j in 1:ncol(sim_eigen)){
    a <- post_eigen$b_intercept[i] + post_eigen$b_age[i]*age[j]
    b <- post_eigen$sigma[i]
    #shapes <- simstudy::betaGetShapes(a, b)
    #sim_eigen[i,j] <- rbeta(1, shape1 = shapes$shape1, shape2 = shapes$shape2)
    sim_eigen[i,j] <- rbeta(1, shape1 = a, shape2 = b)
  }
}
hist(sim_eigen)
sim_pi <- apply(sim_eigen, 2, HPDI, prob = 0.95)

plot(eigen ~ age, data = eigen_values, col = rgb(0,0,1,0.01),
     pch = 19, las = 1, xlim = c(10,55), ylim = c(0,1),
     xlab = 'mean age', ylab = 'eigenvector centrality')
lines(age, mu_mean, lwd = 2, col = 'red')
shade(mu_hpdi, age)
shade(sim_pi, age)

