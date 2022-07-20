#### Bayesian analysis of EfA data ####
# Script to process association data from Makgadikgadi Pans National Park, Botswana
# Data collected: 2012-2021 by Elephants for Africa
# Data supplied by: Kate Evans

#### Set up ####
# load packages
library(tidyverse)
library(dplyr)
#library(rstan)
#library(rethinking)
library(igraph)
library(dagitty)
library(cmdstanr)

# information
sessionInfo()
R.Version()
#rstan::stan_version()
#packageVersion('')
#citation('')

# set stan path
set_cmdstan_path('/Users/helen/.cmdstanr/cmdstan-2.28.2')

# load model
mod_2.2 <- cmdstan_model("models/edge_weight_estimation/simpleBetaNet_HKM_2.2_22.02.03.stan")
mod_2.2

# set seed
set.seed(12345)

################ 1) Draw DAGS ################
# plot with full names
binom <- dagitty::dagitty("dag{
                         age_dyad [exposure];
                         sex_dyad [exposure];
                         weight [outcome];
                         relatedness_dyad [unobserved];
                         age_dyad -> weight <- sex_dyad;
                         weight <- relatedness_dyad;
                         }")
dagitty::coordinates(binom) <- list(x = c(age_dyad = 0, sex_dyad = 1, weight = 2, relatedness = 3),
                                    y = c(age_dyad = 0, sex_dyad = 2, weight = 1, relatedness = 0))
drawdag(binom)

# plot with letters
binom <- dagitty::dagitty("dag{
                         A [exposure];
                         S [exposure];
                         W [outcome];
                         R [unobserved];
                         A -> W <- S;
                         W <- R;
                         }")
dagitty::coordinates(binom) <- list(x = c(A = 0.5, S = 1, W = 1.0, R = 1.5),
                                    y = c(A = 0.5, S = 2, W = 1.2, R = 0.5))
drawdag(binom, radius = 6, cex = 1.6)

# clear environment and reset plot window
rm(binom)
dev.off()
################ 2) Create data lists ################
counts_df1 <- read_csv('data_processed/mpnp_period1_pairwiseevents_22.04.23.csv')
str(counts_df1)

################ Run model on real standardised data -- period 1 ################
### create data list
counts_df_period1 <- counts_df1[counts_df1$period == 1,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period1),          # total number of times one or other of the dyad was observed
  together = counts_df_period1$event_count,    # count number of sightings seen together
  apart    = counts_df_period1$apart,          # count number of sightings seen apart
  period   = counts_df_period1$period)         # which period it's within

### Fit model
weight_mpnp1_2.2_period1 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_mpnp1_2.2_period1
#variable        mean     median    sd   mad         q5        q95 rhat ess_bulk ess_tail
#lp__      -132933.05 -132930.00 97.62 97.85 -133098.00 -132776.00 1.00     1394     1933
#weight[1]       0.05       0.05  0.04  0.03       0.01       0.12 1.00     3484     1794
#weight[2]       0.05       0.04  0.04  0.03       0.01       0.12 1.00     3514     1903
#weight[3]       0.08       0.07  0.05  0.05       0.01       0.17 1.00     4739     2340
#weight[4]       0.08       0.07  0.05  0.05       0.01       0.18 1.00     4905     2312
#weight[5]       0.18       0.17  0.09  0.09       0.05       0.35 1.00     4790     2435
#weight[6]       0.07       0.06  0.05  0.04       0.01       0.16 1.00     4968     2548
#weight[7]       0.20       0.19  0.09  0.09       0.08       0.36 1.00     5330     2555
#weight[8]       0.09       0.08  0.06  0.05       0.02       0.19 1.00     4911     2139
#weight[9]       0.06       0.05  0.04  0.04       0.01       0.14 1.00     4463     2009
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_mpnp1_2.2_period1$output_files()[1])
draws1_mpnp1_1 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_mpnp1_2.2_period1$output_files()[2])
draws2_mpnp1_1 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_mpnp1_2.2_period1$output_files()[3])
draws3_mpnp1_1 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_mpnp1_2.2_period1$output_files()[4])
draws4_mpnp1_1 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_mpnp1_1 <- rbind(draws1_mpnp1_1, draws2_mpnp1_1, draws3_mpnp1_1, draws4_mpnp1_1)

colnames(draws_mpnp1_1)[2:ncol(draws_mpnp1_1)] <- counts_df_period1$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_mpnp1_1), size = 30, replace = F)

### save data 
write_csv(draws_mpnp1_1, 'data_processed/mpnp1_bayesian_edgedistributions_a2.b2_period1_22.04.15.csv')
rm(draws1_mpnp1_1, draws2_mpnp1_1, draws3_mpnp1_1, draws4_mpnp1_1)

### build traceplots -- period 1 -- very wide ####
plot(draws_mpnp1_1[,plot_cols[1]], type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_mpnp1_1[,plot_cols[2]], col = 'tan')
lines(draws_mpnp1_1[,plot_cols[3]], col = 'orange')
lines(draws_mpnp1_1[,plot_cols[4]], col = 'green')
lines(draws_mpnp1_1[,plot_cols[5]], col = 'chocolate')
lines(draws_mpnp1_1[,plot_cols[6]], col = 'blue')
lines(draws_mpnp1_1[,plot_cols[7]], col = 'red')
lines(draws_mpnp1_1[,plot_cols[8]], col = 'seagreen')
lines(draws_mpnp1_1[,plot_cols[9]], col = 'purple')
lines(draws_mpnp1_1[,plot_cols[10]],col = 'magenta')
lines(draws_mpnp1_1[,plot_cols[11]],col = 'black')
lines(draws_mpnp1_1[,plot_cols[12]], col = 'tan')
lines(draws_mpnp1_1[,plot_cols[13]], col = 'orange')
lines(draws_mpnp1_1[,plot_cols[14]], col = 'green')
lines(draws_mpnp1_1[,plot_cols[15]], col = 'chocolate')
lines(draws_mpnp1_1[,plot_cols[16]], col = 'blue')
lines(draws_mpnp1_1[,plot_cols[17]], col = 'red')
lines(draws_mpnp1_1[,plot_cols[18]], col = 'seagreen')
lines(draws_mpnp1_1[,plot_cols[19]], col = 'purple')
lines(draws_mpnp1_1[,plot_cols[20]],col = 'magenta')
lines(draws_mpnp1_1[,plot_cols[21]],col = 'black')
lines(draws_mpnp1_1[,plot_cols[22]], col = 'tan')
lines(draws_mpnp1_1[,plot_cols[23]], col = 'orange')
lines(draws_mpnp1_1[,plot_cols[24]], col = 'green')
lines(draws_mpnp1_1[,plot_cols[25]], col = 'chocolate')
lines(draws_mpnp1_1[,plot_cols[26]], col = 'blue')
lines(draws_mpnp1_1[,plot_cols[27]], col = 'red')
lines(draws_mpnp1_1[,plot_cols[28]], col = 'seagreen')
lines(draws_mpnp1_1[,plot_cols[29]], col = 'purple')
lines(draws_mpnp1_1[,plot_cols[30]],col = 'magenta')

### density plots -- period 1 ####
dens(draws_mpnp1_1[,2], ylim = c(0,10),xlim = c(0,1), las = 1)
for(i in 1:30){  # looks like a poisson distribution - peaks just above 0 and then exponentiol decline. almost nothing above 0.5
  dens(add = T, draws_mpnp1_1[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 1 ####
# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_mpnp1_1[,2:ncol(draws_mpnp1_1)]),
                        min = rep(NA, ncol(draws_mpnp1_1)-1),
                        max = rep(NA, ncol(draws_mpnp1_1)-1),
                        mean = rep(NA, ncol(draws_mpnp1_1)-1),
                        median = rep(NA, ncol(draws_mpnp1_1)-1),
                        sd = rep(NA, ncol(draws_mpnp1_1)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_mpnp1_1[,i+1])
  summaries$max[i]    <- max(draws_mpnp1_1[,i+1])
  summaries$mean[i]   <- mean(draws_mpnp1_1[,i+1])
  summaries$median[i] <- median(draws_mpnp1_1[,i+1])
  summaries$sd[i]     <- sd(draws_mpnp1_1[,i+1])
}

# organise dem_class -- report age category based on age at start of period
plot_data_mpnp1_1 <- left_join(x = summaries, y = counts_df_period1, by = 'dyad')
head(plot_data_mpnp1_1)
plot_data_mpnp1_1$age_cat_1 <- ifelse(plot_data_mpnp1_1$age_start_1 < 6, 'C',
                                     ifelse(plot_data_mpnp1_1$age_start_1 < 11, 'J',
                                            ifelse(plot_data_mpnp1_1$age_start_1 > 19, 'A','P')))
plot_data_mpnp1_1$age_cat_2 <- ifelse(plot_data_mpnp1_1$age_start_2 < 6, 'C',
                                     ifelse(plot_data_mpnp1_1$age_start_2 < 11, 'J',
                                            ifelse(plot_data_mpnp1_1$age_start_2 > 19, 'A','P')))
plot_data_mpnp1_1$age_catnum_1 <- ifelse(plot_data_mpnp1_1$age_start_1 < 6, 1,
                                        ifelse(plot_data_mpnp1_1$age_start_1 < 11, 2,
                                               ifelse(plot_data_mpnp1_1$age_start_1 < 16, 3,
                                                      ifelse(plot_data_mpnp1_1$age_start_1 < 20, 4,
                                                             ifelse(plot_data_mpnp1_1$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_mpnp1_1$age_start_1 < 40, 6, 7))))))
plot_data_mpnp1_1$age_catnum_2 <- ifelse(plot_data_mpnp1_1$age_start_2 < 6, 1,
                                        ifelse(plot_data_mpnp1_1$age_start_2 < 11, 2,
                                               ifelse(plot_data_mpnp1_1$age_start_2 < 16, 3,
                                                      ifelse(plot_data_mpnp1_1$age_start_2 < 20, 4,
                                                             ifelse(plot_data_mpnp1_1$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_mpnp1_1$age_start_2 < 40, 6, 7))))))

plot_data_mpnp1_1$age_dyad <- ifelse(plot_data_mpnp1_1$age_catnum_1 >= plot_data_mpnp1_1$age_catnum_2,
                                    paste(plot_data_mpnp1_1$age_cat_1, plot_data_mpnp1_1$age_cat_2, sep = ''),
                                    paste(plot_data_mpnp1_1$age_cat_2, plot_data_mpnp1_1$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_mpnp1_1, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4, fill = c('blue','purple','blue','grey','purple','blue'))+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_mpnp1_1, aes(x = age_diff, y = mean))+
    geom_jitter(alpha = 0.2)+
    theme_classic()+
    theme(legend.position = 'none')+
    scale_x_continuous('age difference between dyad members',
                       expand = c(0.02,0))+
    scale_y_continuous('mean edge weight',
                       expand = c(0.02,0),
                       limits = c(0,1)))

# values for reporting
summary(summaries) # generally about the same as MOTNP, but higher SD
quantile(summaries$median, seq(0,1,length.out = 101))
quantile(summaries$mean,   seq(0,1,length.out = 101))
hist(summaries$median, breaks = 100, xlim = c(0,1), las = 1,
     xlab = 'median estimated dyad strength', main = '', col = 'skyblue')

m_sum <- plot_data_mpnp1_1[plot_data_mpnp1_1$age_dyad == 'AA' | plot_data_mpnp1_1$age_dyad == 'AP' | plot_data_mpnp1_1$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 1 ################
head(summaries)
length(unique(plot_data_mpnp1_1$id_1))+1 # number of individuals = 55

par(mai = c(0.1,0.1,0.1,0.1))

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period1$id_1))+1,
                         NROW(unique(counts_df_period1$id_2))+1,
                         NROW(draws_mpnp1_1)),
                    dimnames = list(c(unique(counts_df_period1$id_1),
                                      unique(counts_df_period1$id_2)[length(unique(counts_df_period1$id_2))]),
                                    c(unique(counts_df_period1$id_1),
                                      unique(counts_df_period1$id_2)[length(unique(counts_df_period1$id_2))]),
                                    NULL))
N <- nrow(counts_df_period1)

for (i in 1:N) {
  dyad_row <- counts_df_period1[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_mpnp1_1[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids1 <- c(unique(counts_df_period1$id_1), unique(counts_df_period1$id_2)[length(unique(counts_df_period1$id_2))])
males1 <- males[,c(21,6,9,22:53)]
males1 <- males1 %>% dplyr::filter(id %in% ids1)
males1

# create variables for different degrees of node connectedness
males1$degree_0.1 <- NA
males1$degree_0.2 <- NA
males1$degree_0.3 <- NA
males1$degree_0.4 <- NA
males1$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males1)){
  rows <- summaries[summaries$id_1 == males1$id[i] | summaries$id_2 == males1$id[i],]
  males1$degree_0.1[i] <- length(which(rows$median > 0.1))
  males1$degree_0.2[i] <- length(which(rows$median > 0.2))
  males1$degree_0.3[i] <- length(which(rows$median > 0.3))
  males1$degree_0.4[i] <- length(which(rows$median > 0.4))
  males1$degree_0.5[i] <- length(which(rows$median > 0.5))
}

which(males1$degree_0.1 < males1$degree_0.2)
which(males1$degree_0.2 < males1$degree_0.3)
which(males1$degree_0.3 < males1$degree_0.4)
which(males1$degree_0.4 < males1$degree_0.5)

# age variable
males1$age <- lubridate::year(periods$period_start[periods$period == 1]) - males1$byr
summary(males1$age)
males1$age_class <- ifelse(males1$age < 10, 2,
                            ifelse(males1$age < 15, 3,
                                   ifelse(males1$age < 20, 4,
                                          ifelse(males1$age < 25, 5,
                                                 ifelse(males1$age < 40, 6, 7)))))

# Plot whole network
coords <- igraph::layout_nicely(g_mid)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = rgb(0, 0, 0, 0.25), 
     vertex.label = NA,
     vertex.size = 5,
     layout = coords)
plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     edge.color = 'black',
     vertex.size = 8,
     vertex.label = males1$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males1$age_class == 7,'seagreen4',
                          ifelse(males1$age_class == 6,'seagreen3',
                                 ifelse(males1$age_class == 5,'seagreen2',
                                        ifelse(males1$age_class == 4,'steelblue3',
                                               ifelse(males1$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(adj_mid < 0.2,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(adj_mid < 0.2,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males1$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males1$age_class == 7,'seagreen4',
                          ifelse(males1$age_class == 6,'seagreen3',
                                 ifelse(males1$age_class == 5,'seagreen2',
                                        ifelse(males1$age_class == 4,'steelblue3',
                                               ifelse(males1$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.2 -- period 1 ####
g_mid_0.2 <- delete.vertices(graph = g_mid, v = males1$id[which(males1$degree_0.2 == 0)])
g_rng_0.2 <- delete.vertices(graph = g_rng, v = males1$id[which(males1$degree_0.2 == 0)])

set.seed(3)
coords_0.2 <- layout_nicely(g_mid_0.2)
plot(g_mid_0.2,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.2)$weight > 0.2, E(g_rng_0.2)$weight*3, 0),
     vertex.size = 6,
     vertex.label = NA,
     layout = coords_0.2)
plot(g_mid_0.2,
     edge.width = ifelse(E(g_mid_0.2)$weight > 0.2, E(g_mid_0.2)$weight*3, 0),
     edge.color = 'black',
     vertex.size = 6,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males1[which(males1$degree_0.2 > 0),]$age_class == 7,'seagreen4',
                          ifelse(males1[which(males1$degree_0.2 > 0),]$age_class == 6,'seagreen3',
                                 ifelse(males1[which(males1$degree_0.2 > 0),]$age_class == 5,'seagreen2',
                                        ifelse(males1[which(males1$degree_0.2 > 0),]$age_class == 4,'steelblue3',
                                               ifelse(males1[which(males1$degree_0.2 > 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.2, add = TRUE)

par(mai = c(0.1,0.4,0.1,0.2))
plot(g_mid_0.2,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.2)$weight > 0.2, E(g_rng_0.2)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.2)
plot(g_mid_0.2,
     edge.width = ifelse(E(g_mid_0.2)$weight > 0.2, E(g_mid_0.2)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males1$p1[which(males1$degree_0.2 != 0)],
     vertex.label = males1$id[which(males1$degree_0.2 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males1[which(males1$degree_0.2 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males1[which(males1$degree_0.2 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males1[which(males1$degree_0.2 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males1[which(males1$degree_0.2 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males1[which(males1$degree_0.2 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.2, add = TRUE)

### save summary data -- period 1 ####
dyad_period_weights <- counts_df1
summaries$period <- 1
summaries <- summaries[,c(1,4:9)]
dyad_period_weights <- left_join(x = dyad_period_weights, y = summaries,
                                 by = c('dyad','period'))
rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.2,
   counts_df_period1, draws_mpnp1_1, dyad_row, g_mid,g_mid_0.2, g_rng, g_rng_0.2, males1, plot_data_mpnp1_1, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids1, N, plot_cols)


################ Graphs of MPNP1 -- THESE TWO SCRIPTS DO NOT MATCH UP: DIFFERENT INPUT DATA ################
### elephant data
counts_df1 <- read_csv('data_processed/mpnp_period1_pairwiseevents_22.05.30.csv')

### model data
draws_mpnp1 <- readRDS('data_processed/mpnp1_bayesian_edgedistributions_a2.b2_period1_22.05.30.rds')

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_mpnp1), size = 30, replace = F)

### build traceplots -- period 1 -- very wide ####
plot(draws_mpnp1[,plot_cols[1]], type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_mpnp1[,plot_cols[2]], col = 'tan')
lines(draws_mpnp1[,plot_cols[3]], col = 'orange')
lines(draws_mpnp1[,plot_cols[4]], col = 'green')
lines(draws_mpnp1[,plot_cols[5]], col = 'chocolate')
lines(draws_mpnp1[,plot_cols[6]], col = 'blue')
lines(draws_mpnp1[,plot_cols[7]], col = 'red')
lines(draws_mpnp1[,plot_cols[8]], col = 'seagreen')
lines(draws_mpnp1[,plot_cols[9]], col = 'purple')
lines(draws_mpnp1[,plot_cols[10]],col = 'magenta')
lines(draws_mpnp1[,plot_cols[11]],col = 'black')
lines(draws_mpnp1[,plot_cols[12]], col = 'tan')
lines(draws_mpnp1[,plot_cols[13]], col = 'orange')
lines(draws_mpnp1[,plot_cols[14]], col = 'green')
lines(draws_mpnp1[,plot_cols[15]], col = 'chocolate')
lines(draws_mpnp1[,plot_cols[16]], col = 'blue')
lines(draws_mpnp1[,plot_cols[17]], col = 'red')
lines(draws_mpnp1[,plot_cols[18]], col = 'seagreen')
lines(draws_mpnp1[,plot_cols[19]], col = 'purple')
lines(draws_mpnp1[,plot_cols[20]],col = 'magenta')
lines(draws_mpnp1[,plot_cols[21]],col = 'black')
lines(draws_mpnp1[,plot_cols[22]], col = 'tan')
lines(draws_mpnp1[,plot_cols[23]], col = 'orange')
lines(draws_mpnp1[,plot_cols[24]], col = 'green')
lines(draws_mpnp1[,plot_cols[25]], col = 'chocolate')
lines(draws_mpnp1[,plot_cols[26]], col = 'blue')
lines(draws_mpnp1[,plot_cols[27]], col = 'red')
lines(draws_mpnp1[,plot_cols[28]], col = 'seagreen')
lines(draws_mpnp1[,plot_cols[29]], col = 'purple')
lines(draws_mpnp1[,plot_cols[30]],col = 'magenta')

### density plots -- period 1 ####
dens(draws_mpnp1[,2], ylim = c(0,10),xlim = c(0,1), las = 1)
for(i in 1:30){  # looks like a poisson distribution - peaks just above 0 and then exponential decline. almost nothing above 0.5
  dens(add = T, draws_mpnp1[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

# print progress stamp
print(paste0('Chain checking for period 1 completed at ', Sys.time()))

### summarise ####
summaries <- data.frame(dyad = colnames(draws_mpnp1[,2:ncol(draws_mpnp1)]),
                        min = rep(NA, ncol(draws_mpnp1)-1),
                        max = rep(NA, ncol(draws_mpnp1)-1),
                        mean = rep(NA, ncol(draws_mpnp1)-1),
                        median = rep(NA, ncol(draws_mpnp1)-1),
                        sd = rep(NA, ncol(draws_mpnp1)-1),
                        q2.5 = rep(NA, ncol(draws_mpnp1)-1),
                        q97.5 = rep(NA, ncol(draws_mpnp1)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_mpnp1[,i+1])
  summaries$max[i]    <- max(draws_mpnp1[,i+1])
  summaries$mean[i]   <- mean(draws_mpnp1[,i+1])
  summaries$median[i] <- median(draws_mpnp1[,i+1])
  summaries$sd[i]     <- sd(draws_mpnp1[,i+1])
  summaries$q2.5[i]   <- quantile(draws_mpnp1[,i+1], 0.025)
  summaries$q97.5[i]  <- quantile(draws_mpnp1[,i+1], 0.975)
  if(i %% 1000 == 0) { print(paste0('Summary row ', i, ' completed at: ', Sys.time() )) }
}

# print progress stamp
print(paste0('Summaries for period 1 completed at ', Sys.time()))

### create network plots -- period 1 ################
head(summaries)
length(unique(summaries$id_1))+1 # number of individuals = 55

par(mai = c(0.1,0.1,0.1,0.1))

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df1$id_1))+1,
                         NROW(unique(counts_df1$id_2))+1,
                         NROW(draws_mpnp1)),
                    dimnames = list(c(unique(counts_df1$id_1),
                                      unique(counts_df1$id_2)[length(unique(counts_df1$id_2))]),
                                    c(unique(counts_df1$id_1),
                                      unique(counts_df1$id_2)[length(unique(counts_df1$id_2))]),
                                    NULL))
N <- nrow(counts_df1)

for (i in 1:N) {
  dyad_row <- counts_df1[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_mpnp1[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
males1 <- data.frame(id = sort(unique(c(counts_df1$id_1, counts_df1$id_2))),
                     age_median = NA,
                     count_period = NA)
for(i in 1:nrow(males1)){
  male_id1 <- counts_df1[counts_df1$id_1 == males1$id[i],c('id_1','age_median_1','count_period_1')]
  male_id2 <- counts_df1[counts_df1$id_2 == males1$id[i],c('id_2','age_median_2','count_period_2')]
  colnames(male_id1) <- c('id','age_median','count_period')
  colnames(male_id2) <- c('id','age_median','count_period')
  male <- rbind(male_id1, male_id2) %>% distinct()
  males1$count_period[i] <- male$count_period[1]
  if(nrow(male) == 1){
    males1$age_median[i] <- ifelse(is.na(male$age_median[1]) == TRUE, 10, male$age_median[1])
  }
  else {
    males1$age_median[i] <- median(male$age_median)
  }
  rm(male_id1, male_id2, male)
}

# create variables for different degrees of node connectedness
males1$degree_0.1 <- NA
males1$degree_0.2 <- NA
males1$degree_0.3 <- NA
males1$degree_0.4 <- NA
males1$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males1)){
  rows <- summaries[summaries$id_1 == males1$id[i] | summaries$id_2 == males1$id[i],]
  males1$degree_0.1[i] <- length(which(rows$median > 0.1))
  males1$degree_0.2[i] <- length(which(rows$median > 0.2))
  males1$degree_0.3[i] <- length(which(rows$median > 0.3))
  males1$degree_0.4[i] <- length(which(rows$median > 0.4))
  males1$degree_0.5[i] <- length(which(rows$median > 0.5))
}

which(males1$degree_0.1 < males1$degree_0.2)
which(males1$degree_0.2 < males1$degree_0.3)
which(males1$degree_0.3 < males1$degree_0.4)
which(males1$degree_0.4 < males1$degree_0.5)

# age variable
males1$age_class <- ifelse(males1$age_median == 10, 10,
                           ifelse(males1$age_median < 4, 2,
                                  ifelse(males1$age_median < 5, 3,
                                         ifelse(males1$age_median < 6, 4, 
                                                ifelse(males1$age_median < 7, 5, 
                                                       ifelse(males1$age_median < 8, 6, 7))))))

# print progress stamp
print(paste0('Networks ready to plot for period 1 completed at ', Sys.time()))

# Plot whole network
coords <- igraph::layout_nicely(g_mid)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = rgb(0, 0, 0, 0.25), 
     vertex.label = NA,
     vertex.size = 5,
     layout = coords)
plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     edge.color = 'black',
     vertex.size = 8,
     vertex.label = males1$elephant,
     vertex.label.dist = 0,
     vertex.label.color = ifelse(males1$age_class == 10,'white','black'),
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males1$age_class == 7,'seagreen4',
                          ifelse(males1$age_class == 6,'seagreen3',
                                 ifelse(males1$age_class == 5,'seagreen2',
                                        ifelse(males1$age_class == 4,'steelblue3',
                                               ifelse(males1$age_class == 3,'steelblue1',
                                                      ifelse(males1$age_class == 2,'yellow',
                                                             'black')))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(adj_mid < 0.2,'transparent',rgb(0,0,0,0.25)),
     vertex.label = NA,
     vertex.size = 5,
     layout = coords)
plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     edge.color = ifelse(adj_mid < 0.2,'transparent','black'),
     vertex.size = 8,
     vertex.label = males1$elephant,
     vertex.label.dist = 0,
     vertex.label.color = ifelse(males1$age_class == 10,'white','black'),
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males1$age_class == 7,'seagreen4',
                          ifelse(males1$age_class == 6,'seagreen3',
                                 ifelse(males1$age_class == 5,'seagreen2',
                                        ifelse(males1$age_class == 4,'steelblue3',
                                               ifelse(males1$age_class == 3,'steelblue1',
                                                      ifelse(males1$age_class == 2,'yellow',
                                                             'black')))))),
     layout = coords, add = TRUE)

# print progress stamp
print(paste0('Full population plots for period 1 completed at ', Sys.time()))

### only elephants degree > 0.3 -- period 1 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males1$id[which(males1$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males1$id[which(males1$degree_0.3 == 0)])

males1_0.3 <- males1[males1$degree_0.3 > 0,]

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*10, 0),
     edge.color = rgb(0,0,0,0.25),
     vertex.size = 6,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*10, 0),
     edge.color = 'black',
     vertex.size = 6,
     vertex.label = males1_0.3$elephant,
     vertex.label.color = ifelse(males1_0.3$age_class == 10,'white','black'),
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(males1_0.3$age_class == 7,'seagreen4',
                           ifelse(males1_0.3$age_class == 6,'seagreen3',
                                  ifelse(males1_0.3$age_class == 5,'seagreen2',
                                         ifelse(males1_0.3$age_class == 4,'steelblue3',
                                                ifelse(males1_0.3$age_class == 3,'steelblue1',
                                                       ifelse(males1_0.3$age_class == 2,'yellow',
                                                              'black')))))),
     layout = coords_0.3, add = TRUE)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*10, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*10, 0),
     edge.color = 'black',
     vertex.size = males1_0.3$count_period*8,
     vertex.label = males1_0.3$id,
     vertex.label.family = 'Helvetica',
     vertex.label.color = ifelse(males1_0.3$age_class == 10,'white','black'),
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(males1_0.3$age_class == 7,'seagreen4',
                           ifelse(males1_0.3$age_class == 6,'seagreen3',
                                  ifelse(males1_0.3$age_class == 5,'seagreen2',
                                         ifelse(males1_0.3$age_class == 4,'steelblue3',
                                                ifelse(males1_0.3$age_class == 3,'steelblue1',
                                                       ifelse(males1_0.3$age_class == 2,'yellow',
                                                              'black')))))),
     layout = coords_0.3, add = TRUE)

# print progress stamp
print(paste0('All network plots for period 1 completed at ', Sys.time()))

### save summary data -- period 1 ####
dyad_period_weights <- counts_df1
summaries$period <- 1
summaries <- summaries[,c(1,4:11)]
dyad_period_weights <- left_join(x = dyad_period_weights, y = summaries,
                                 by = c('dyad','period'))
rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df1, draws_mpnp1, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males1, rows, summaries, adj_quantiles, adj_tensor, i, N, plot_cols)

# print progress stamp
print(paste0('Time period 1 completed at ', Sys.time()))

# save summary data
write_csv(dyad_period_weights, 'data_processed/mpnp_dyad_weightdistributions_2.2_period1_22.05.31.csv')

# save graphs
dev.off()
