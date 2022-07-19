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

################ 6) Run model on real standardised data -- all ################
### create data list
counts_ls1 <- list(
  n_dyads  = nrow(counts_df1),          # total number of times one or other of the dyad was observed
  together = counts_df1$together   ,    # count number of sightings seen together
  apart    = counts_df1$apart,          # count number of sightings seen apart
  period   = counts_df1$period)         # which period it's within

### Fit model
weight_mpnp1_2.2 <- mod_2.2$sample(
  data = counts_ls1, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_mpnp1_2.2
# variable        mean      median     sd    mad          q5         q95 rhat ess_bulk ess_tail
#lp__      -3068878.46 -3068870.00 554.23 548.56 -3069780.00 -3067960.00 1.00     1334     2343
#weight[1]        0.06        0.06   0.04   0.04        0.01        0.15 1.00     4577     1902
#weight[2]        0.06        0.05   0.04   0.04        0.01        0.15 1.00     4262     2326
#weight[3]        0.06        0.05   0.04   0.04        0.01        0.15 1.00     4455     2159
#weight[4]        0.06        0.05   0.04   0.04        0.01        0.15 1.00     4756     2696
#weight[5]        0.05        0.04   0.03   0.03        0.01        0.11 1.00     3622     1682
#weight[6]        0.05        0.04   0.03   0.03        0.01        0.11 1.00     5156     2216
#weight[7]        0.04        0.04   0.03   0.03        0.01        0.10 1.00     4752     2229
#weight[8]        0.05        0.04   0.03   0.03        0.01        0.11 1.00     3381     1475
#weight[9]        0.05        0.04   0.03   0.03        0.01        0.11 1.00     4109     1927
rm(counts_ls1)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_mpnp1_2.2$output_files()[1])
draws1_mpnp1_2.2 <- as.data.frame(output1$post_warmup_draws)
rm(output1)
write_csv('data_processed/mpnp1_bayesian_edgedistributions_a2.b2_chain1_22.04.06.csv', x = draws1_mpnp1_2.2)

output2 <- read_cmdstan_csv(weight_mpnp1_2.2$output_files()[2])
draws2_mpnp1_2.2 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_mpnp1_2.2$output_files()[3])
draws3_mpnp1_2.2 <- as.data.frame(output3$post_warmup_draws)
r(output3)

output4 <- read_cmdstan_csv(weight_mpnp1_2.2$output_files()[4])
draws4_mpnp1_2.2 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_mpnp1_2.2 <- rbind(draws1_mpnp1_2.2, draws2_mpnp1_2.2, draws3_mpnp1_2.2, draws4_mpnp1_2.2)

colnames(draws_mpnp1_2.2)[2:531532] <- counts_df1$dyad

### save data 
write_csv(draws_mpnp1_2.2, 'data_processed/mpnp1_bayesian_edgedistributions_a2.b2_22.04.05.csv')

################ check single chain ################
draws1_mpnp1_2.2 <- read_csv('data_processed/mpnp1_bayesian_edgedistributions_a2.b2_chain1_22.04.06.csv')
plot_cols <- sample(2:531532, size = 30, replace = F)

### build traceplots -- some quite wide, generally not bad ####
plot(draws1_mpnp1_2.2[,plot_cols[1]], type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws1_mpnp1_2.2[,plot_cols[2]], col = 'tan')
lines(draws1_mpnp1_2.2[,plot_cols[3]], col = 'orange')
lines(draws1_mpnp1_2.2[,plot_cols[4]], col = 'green')
lines(draws1_mpnp1_2.2[,plot_cols[5]], col = 'chocolate')
lines(draws1_mpnp1_2.2[,plot_cols[6]], col = 'blue')
lines(draws1_mpnp1_2.2[,plot_cols[7]], col = 'red')
lines(draws1_mpnp1_2.2[,plot_cols[8]], col = 'seagreen')
lines(draws1_mpnp1_2.2[,plot_cols[9]], col = 'purple')
lines(draws1_mpnp1_2.2[,plot_cols[10]],col = 'magenta')
lines(draws1_mpnp1_2.2[,plot_cols[11]],col = 'black')
lines(draws1_mpnp1_2.2[,plot_cols[12]], col = 'tan')
lines(draws1_mpnp1_2.2[,plot_cols[13]], col = 'orange')
lines(draws1_mpnp1_2.2[,plot_cols[14]], col = 'green')
lines(draws1_mpnp1_2.2[,plot_cols[15]], col = 'chocolate')
lines(draws1_mpnp1_2.2[,plot_cols[16]], col = 'blue')
lines(draws1_mpnp1_2.2[,plot_cols[17]], col = 'red')
lines(draws1_mpnp1_2.2[,plot_cols[18]], col = 'seagreen')
lines(draws1_mpnp1_2.2[,plot_cols[19]], col = 'purple')
lines(draws1_mpnp1_2.2[,plot_cols[20]],col = 'magenta')
lines(draws1_mpnp1_2.2[,plot_cols[21]],col = 'black')
lines(draws1_mpnp1_2.2[,plot_cols[22]], col = 'tan')
lines(draws1_mpnp1_2.2[,plot_cols[23]], col = 'orange')
lines(draws1_mpnp1_2.2[,plot_cols[24]], col = 'green')
lines(draws1_mpnp1_2.2[,plot_cols[25]], col = 'chocolate')
lines(draws1_mpnp1_2.2[,plot_cols[26]], col = 'blue')
lines(draws1_mpnp1_2.2[,plot_cols[27]], col = 'red')
lines(draws1_mpnp1_2.2[,plot_cols[28]], col = 'seagreen')
lines(draws1_mpnp1_2.2[,plot_cols[29]], col = 'purple')
lines(draws1_mpnp1_2.2[,plot_cols[30]],col = 'magenta')

### density plots -- most are very wide, high uncertainty####
dens(draws1_mpnp1_2.2[,2], ylim = c(0,100),xlim = c(0,1), las = 1)
for(i in 3:8000){  # looks like a poisson distribution - peaks just above 0 and then exponentiol decline. almost nothing above 0.5
  dens(add = T, draws1_mpnp1_2.2[,i], col = col.alpha('blue', alpha = 0.01))
}



### summarise data ####
# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws1_mpnp1_2.2[,2:531532]),
                        min = rep(NA, ncol(draws1_mpnp1_2.2)-1),
                        max = rep(NA, ncol(draws1_mpnp1_2.2)-1),
                        mean = rep(NA, ncol(draws1_mpnp1_2.2)-1),
                        median = rep(NA, ncol(draws1_mpnp1_2.2)-1),
                        sd = rep(NA, ncol(draws1_mpnp1_2.2)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws1_mpnp1_2.2[,i+1])
  summaries$max[i]    <- max(draws1_mpnp1_2.2[,i+1])
  summaries$mean[i]   <- mean(draws1_mpnp1_2.2[,i+1])
  summaries$median[i] <- median(draws1_mpnp1_2.2[,i+1])
  summaries$sd[i]     <- sd(draws1_mpnp1_2.2[,i+1])
}

summary(summaries) # generally higher than MOTNP but SD much higher

# organise dem_class -- report age category based on age at start of period
plot_data_mpnp1_2.2 <- left_join(x = summaries, y = counts_df1, by = 'dyad')
head(plot_data_mpnp1_2.2)
plot_data_mpnp1_2.2$age_cat_1 <- ifelse(plot_data_mpnp1_2.2$age_start_1 < 6, 'C',
                                     ifelse(plot_data_mpnp1_2.2$age_start_1 < 11, 'J',
                                            ifelse(plot_data_mpnp1_2.2$age_start_1 > 19, 'A','P')))
plot_data_mpnp1_2.2$age_cat_2 <- ifelse(plot_data_mpnp1_2.2$age_start_2 < 6, 'C',
                                     ifelse(plot_data_mpnp1_2.2$age_start_2 < 11, 'J',
                                            ifelse(plot_data_mpnp1_2.2$age_start_2 > 19, 'A','P')))
plot_data_mpnp1_2.2$age_catnum_1 <- ifelse(plot_data_mpnp1_2.2$age_start_1 < 6, 1,
                                        ifelse(plot_data_mpnp1_2.2$age_start_1 < 11, 2,
                                               ifelse(plot_data_mpnp1_2.2$age_start_1 < 16, 3,
                                                      ifelse(plot_data_mpnp1_2.2$age_start_1 < 20, 4,
                                                             ifelse(plot_data_mpnp1_2.2$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_mpnp1_2.2$age_start_1 < 40, 6, 7))))))
plot_data_mpnp1_2.2$age_catnum_2 <- ifelse(plot_data_mpnp1_2.2$age_start_2 < 6, 1,
                                        ifelse(plot_data_mpnp1_2.2$age_start_2 < 11, 2,
                                               ifelse(plot_data_mpnp1_2.2$age_start_2 < 16, 3,
                                                      ifelse(plot_data_mpnp1_2.2$age_start_2 < 20, 4,
                                                             ifelse(plot_data_mpnp1_2.2$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_mpnp1_2.2$age_start_2 < 40, 6, 7))))))

plot_data_mpnp1_2.2$age_dyad <- ifelse(plot_data_mpnp1_2.2$age_catnum_1 >= plot_data_mpnp1_2.2$age_catnum_2,
                                    paste(plot_data_mpnp1_2.2$age_cat_1, plot_data_mpnp1_2.2$age_cat_2, sep = ''),
                                    paste(plot_data_mpnp1_2.2$age_cat_2, plot_data_mpnp1_2.2$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_mpnp1_2.2, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4, fill = c('blue','purple','blue','grey','purple','blue'))+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_mpnp1_2.2, aes(x = age_diff, y = mean))+
    geom_jitter(alpha = 0.2)+
    theme_classic()+
    theme(legend.position = 'none')+
    scale_x_continuous('age difference between dyad members',
                       expand = c(0.02,0))+
    scale_y_continuous('mean edge weight',
                       expand = c(0.02,0),
                       limits = c(0,1)))

# values for reporting
summary(summaries)
quantile(summaries$median, seq(0,1,length.out = 101))
quantile(summaries$mean,   seq(0,1,length.out = 101))
hist(summaries$median, breaks = 100, xlim = c(0,1), las = 1,
     xlab = 'median estimated dyad strength', main = '', col = 'skyblue')

m_sum <- plot_data_mpnp1_2.2[plot_data_mpnp1_2.2$age_dyad == 'AA' | plot_data_mpnp1_2.2$age_dyad == 'AP' | plot_data_mpnp1_2.2$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(counts_ls)
rm(draws1_mpnp1_2.2, draws2_mpnp1_2.2, draws3_mpnp1_2.2, draws4_mpnp1_2.2, tidy_draws1_2.2)
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- single chain ################
head(summaries)
length(unique(plot_data_mpnp1_2.2$id_1))+1 # number of individuals = 55

par(mai = c(0.1,0.1,0.1,0.1))

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df1$id_1))+1,
                         NROW(unique(counts_df1$id_2))+1,
                         NROW(draws1_mpnp1_2.2)),
                    dimnames = list(c(unique(counts_df1$id_1),'M119'),
                                    c('M2',unique(counts_df1$id_2)),
                                    NULL))
N <- nrow(counts_df1)

for (i in 1:N) {
  dyad_row <- counts_df1[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws1_mpnp1_2.2[, i+1]
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
males <- males[,c(21,6,9,22:53)]
View(males) # all ones born before start of study

# create variables for different degrees of node connectedness
males$degree_0.1 <- NA
males$degree_0.2 <- NA
males$degree_0.3 <- NA
males$degree_0.4 <- NA
males$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males)){
  rows <- summaries[summaries$id_1 == males$id[i] | summaries$id_2 == males$id[i],]
  males$degree_0.1[i] <- length(which(rows$median > 0.1))
  males$degree_0.2[i] <- length(which(rows$median > 0.2))
  males$degree_0.3[i] <- length(which(rows$median > 0.3))
  males$degree_0.4[i] <- length(which(rows$median > 0.4))
  males$degree_0.5[i] <- length(which(rows$median > 0.5))
}

which(males$degree_0.1 < males$degree_0.2)
which(males$degree_0.2 < males$degree_0.3)
which(males$degree_0.3 < males$degree_0.4)
which(males$degree_0.4 < males$degree_0.5)

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
     vertex.label = males$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     #vertex.color= ifelse(males$age_class == 'Adult','seagreen1',
     #                    ifelse(males$age_class == 'Pubescent','skyblue',
     #                           ifelse(males$age_class == 'Juvenile', 'yellow','magenta'))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(adj_mid < 0.5,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(adj_mid < 0.5,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     #vertex.color= ifelse(males$age_class == 'Adult','seagreen1',
     #                     ifelse(males$age_class == 'Pubescent','skyblue',
     #                            ifelse(males$age_class == 'Juvenile', 'yellow','magenta'))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males$id[which(males$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males$id[which(males$degree_0.3 == 0)])

set.seed(2)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 6,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = 6,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     #vertex.color = ifelse(nodes[which(males$degree_0.3 != 0),]$age_class == 'Adult',
     #                      'seagreen1',
     #                      ifelse(nodes[which(males$degree_0.3 != 0),]$age_class == 'Pubescent','skyblue',
     #                             ifelse(nodes[which(males$degree_0.3 != 0),]$age_class == 'Juvenile','yellow','magenta'))),
     layout = coords_0.3, add = TRUE)

plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = log(males$count_all[which(males$degree_0.3 != 0)]),
     vertex.color = 'green',
     vertex.label = males$id[which(males$degree_0.3 != 0)],
     vertex.label.color = 'magenta',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     #vertex.color = ifelse(males[which(males$degree_0.3 != 0),]$age_class == 'Adult',
     #                      'seagreen1',
     #                      ifelse(males[which(males$degree_0.3 != 0),]$age_class == 'Pubescent','skyblue',
     #                             ifelse(males[which(males$degree_0.3 != 0),]$age_class == 'Juvenile','yellow','magenta'))),
     layout = coords_0.3, add = TRUE)

################ repeat for all chains ################
# draws_mpnp1_2.2 <- read_csv('data_processed/mpnp1_bayesian_edgedistributions_a2.b2_22.04.05.csv') %>% 
#   data.matrix()
plot_cols <- sample(2:531532, size = 20, replace = F)
### build traceplots -- most are very wide, high uncertainty ####
plot(draws_mpnp1_2.2[,plot_cols[1]], type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_mpnp1_2.2[,plot_cols[2]], col = 'tan')
lines(draws_mpnp1_2.2[,plot_cols[3]], col = 'orange')
lines(draws_mpnp1_2.2[,plot_cols[4]], col = 'green')
lines(draws_mpnp1_2.2[,plot_cols[5]], col = 'chocolate')
lines(draws_mpnp1_2.2[,plot_cols[6]], col = 'blue')
lines(draws_mpnp1_2.2[,plot_cols[7]], col = 'red')
lines(draws_mpnp1_2.2[,plot_cols[8]], col = 'seagreen')
lines(draws_mpnp1_2.2[,plot_cols[9]], col = 'purple')
lines(draws_mpnp1_2.2[,plot_cols[10]],col = 'magenta')
lines(draws_mpnp1_2.2[,plot_cols[11]],col = 'black')
lines(draws_mpnp1_2.2[,plot_cols[12]], col = 'tan')
lines(draws_mpnp1_2.2[,plot_cols[13]], col = 'orange')
lines(draws_mpnp1_2.2[,plot_cols[14]], col = 'green')
lines(draws_mpnp1_2.2[,plot_cols[15]], col = 'chocolate')
lines(draws_mpnp1_2.2[,plot_cols[16]], col = 'blue')
lines(draws_mpnp1_2.2[,plot_cols[17]], col = 'red')
lines(draws_mpnp1_2.2[,plot_cols[18]], col = 'seagreen')
lines(draws_mpnp1_2.2[,plot_cols[19]], col = 'purple')
lines(draws_mpnp1_2.2[,plot_cols[20]],col = 'magenta')

### density plots -- most are very wide, high uncertainty ####
dens(draws_mpnp1_2.2[,2], ylim = c(0,10),xlim = c(0,1), las = 1)
for(i in 3:531532){
  dens(add = T, draws_mpnp1_2.2[,i], col = col.alpha('blue', alpha = 0.1))
}

### summarise data ####
# summarise -- look for any anomalies in draw values or chain variation
summaries <- data.frame(dyad = colnames(draws_mpnp1_2.2[,2:531532]),
                        min = rep(NA, ncol(draws_mpnp1_2.2)-1),
                        max = rep(NA, ncol(draws_mpnp1_2.2)-1),
                        mean = rep(NA, ncol(draws_mpnp1_2.2)-1),
                        median = rep(NA, ncol(draws_mpnp1_2.2)-1),
                        sd = rep(NA, ncol(draws_mpnp1_2.2)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_mpnp1_2.2[,i+1])
  summaries$max[i]    <- max(draws_mpnp1_2.2[,i+1])
  summaries$mean[i]   <- mean(draws_mpnp1_2.2[,i+1])
  summaries$median[i] <- median(draws_mpnp1_2.2[,i+1])
  summaries$sd[i]     <- sd(draws_mpnp1_2.2[,i+1])
}

summary(summaries) # generally higher than MOTNP but SD much higher -- higher mean/median doesn't really make sense given that no pair was ever seen more than twice together in any 2-year period

# organise dem_class -- report age category based on age at start of period
plot_data_mpnp1_2.2 <- left_join(x = summaries, y = counts_df1, by = 'dyad')
head(plot_data_mpnp1_2.2)
plot_data_mpnp1_2.2$age_cat_1 <- ifelse(plot_data_mpnp1_2.2$age_start_1 < 6, 'C',
                                     ifelse(plot_data_mpnp1_2.2$age_start_1 < 11, 'J',
                                            ifelse(plot_data_mpnp1_2.2$age_start_1 > 19, 'A','P')))
plot_data_mpnp1_2.2$age_cat_2 <- ifelse(plot_data_mpnp1_2.2$age_start_2 < 6, 'C',
                                     ifelse(plot_data_mpnp1_2.2$age_start_2 < 11, 'J',
                                            ifelse(plot_data_mpnp1_2.2$age_start_2 > 19, 'A','P')))
plot_data_mpnp1_2.2$age_catnum_1 <- ifelse(plot_data_mpnp1_2.2$age_start_1 < 6, 1,
                                        ifelse(plot_data_mpnp1_2.2$age_start_1 < 11, 2,
                                               ifelse(plot_data_mpnp1_2.2$age_start_1 < 16, 3,
                                                      ifelse(plot_data_mpnp1_2.2$age_start_1 < 20, 4,
                                                             ifelse(plot_data_mpnp1_2.2$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_mpnp1_2.2$age_start_1 < 40, 6, 7))))))
plot_data_mpnp1_2.2$age_catnum_2 <- ifelse(plot_data_mpnp1_2.2$age_start_2 < 6, 1,
                                        ifelse(plot_data_mpnp1_2.2$age_start_2 < 11, 2,
                                               ifelse(plot_data_mpnp1_2.2$age_start_2 < 16, 3,
                                                      ifelse(plot_data_mpnp1_2.2$age_start_2 < 20, 4,
                                                             ifelse(plot_data_mpnp1_2.2$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_mpnp1_2.2$age_start_2 < 40, 6, 7))))))

plot_data_mpnp1_2.2$age_dyad <- ifelse(plot_data_mpnp1_2.2$age_catnum_1 >= plot_data_mpnp1_2.2$age_catnum_2,
                                    paste(plot_data_mpnp1_2.2$age_cat_1, plot_data_mpnp1_2.2$age_cat_2, sep = ''),
                                    paste(plot_data_mpnp1_2.2$age_cat_2, plot_data_mpnp1_2.2$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_mpnp1_2.2, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4, fill = c('blue','purple','blue','grey','purple','blue'))+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_mpnp1_2.2, aes(x = age_diff, y = mean))+
    geom_jitter(alpha = 0.2)+
    theme_classic()+
    theme(legend.position = 'none')+
    scale_x_continuous('age difference between dyad members',
                       expand = c(0.02,0))+
    scale_y_continuous('mean edge weight',
                       expand = c(0.02,0),
                       limits = c(0,1)))

# values for reporting
summary(summaries)
quantile(summaries$median, seq(0,1,length.out = 101))
quantile(summaries$mean,   seq(0,1,length.out = 101))
hist(summaries$median, breaks = 100, xlim = c(0,1), las = 1,
     xlab = 'median estimated dyad strength', main = '', col = 'skyblue')

m_sum <- plot_data_mpnp1_2.2[plot_data_mpnp1_2.2$age_dyad == 'AA' | plot_data_mpnp1_2.2$age_dyad == 'AP' | plot_data_mpnp1_2.2$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(counts_ls)
rm(draws1_mpnp1_2.2, draws2_mpnp1_2.2, draws3_mpnp1_2.2, draws4_mpnp1_2.2, tidy_draws_2.2)
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots ################
head(summaries)
length(unique(plot_data_mpnp1_2.2$id_1))+1 # number of individuals = 55

par(mai = c(0.1,0.1,0.1,0.1))

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df1$id_1))+1,
                         NROW(unique(counts_df1$id_2))+1,
                         NROW(draws_mpnp1_2.2)),
                    dimnames = list(c(unique(counts_df1$id_1),'M119'),
                                    c('M2',unique(counts_df1$id_2)),
                                    NULL))
N <- nrow(counts_df1)

for (i in 1:N) {
  dyad_row <- counts_df1[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_mpnp1_2.2[, i+1]
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
males <- males[,c(21,6,9,22:53)]
View(males) # all ones born before start of study

# create variables for different degrees of node connectedness
males$degree_0.1 <- NA
males$degree_0.2 <- NA
males$degree_0.3 <- NA
males$degree_0.4 <- NA
males$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males)){
  rows <- summaries[summaries$id_1 == males$id[i] | summaries$id_2 == males$id[i],]
  males$degree_0.1[i] <- length(which(rows$median > 0.1))
  males$degree_0.2[i] <- length(which(rows$median > 0.2))
  males$degree_0.3[i] <- length(which(rows$median > 0.3))
  males$degree_0.4[i] <- length(which(rows$median > 0.4))
  males$degree_0.5[i] <- length(which(rows$median > 0.5))
}

which(males$degree_0.1 < males$degree_0.2)
which(males$degree_0.2 < males$degree_0.3)
which(males$degree_0.3 < males$degree_0.4)
which(males$degree_0.4 < males$degree_0.5)

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
     vertex.label = males$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     #vertex.color= ifelse(males$age_class == 'Adult','seagreen1',
     #                    ifelse(males$age_class == 'Pubescent','skyblue',
     #                           ifelse(males$age_class == 'Juvenile', 'yellow','magenta'))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(adj_mid < 0.5,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(adj_mid < 0.5,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     #vertex.color= ifelse(males$age_class == 'Adult','seagreen1',
     #                     ifelse(males$age_class == 'Pubescent','skyblue',
     #                            ifelse(males$age_class == 'Juvenile', 'yellow','magenta'))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males$id[which(males$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males$id[which(males$degree_0.3 == 0)])

set.seed(2)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 6,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = 6,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     #vertex.color = ifelse(nodes[which(males$degree_0.3 != 0),]$age_class == 'Adult',
     #                      'seagreen1',
     #                      ifelse(nodes[which(males$degree_0.3 != 0),]$age_class == 'Pubescent','skyblue',
     #                             ifelse(nodes[which(males$degree_0.3 != 0),]$age_class == 'Juvenile','yellow','magenta'))),
     layout = coords_0.3, add = TRUE)

plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = log(males$count_all[which(males$degree_0.3 != 0)]),
     vertex.color = 'green',
     vertex.label = males$id[which(males$degree_0.3 != 0)],
     vertex.label.color = 'magenta',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     #vertex.color = ifelse(males[which(males$degree_0.3 != 0),]$age_class == 'Adult',
     #                      'seagreen1',
     #                      ifelse(males[which(males$degree_0.3 != 0),]$age_class == 'Pubescent','skyblue',
     #                             ifelse(males[which(males$degree_0.3 != 0),]$age_class == 'Juvenile','yellow','magenta'))),
     layout = coords_0.3, add = TRUE)



################ 7.1) Run model on real standardised data -- period 1 ################
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

################ 7.15) Run model on real standardised data -- period 15 ################
### create data list
counts_df_period15 <- counts_df1[counts_df1$period == 15,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period15),          # total number of times one or other of the dyad was observed
  together = counts_df_period15$event_count,    # count number of sightings seen together
  apart    = counts_df_period15$apart,          # count number of sightings seen apart
  period   = counts_df_period15$period)         # which period it's within

### Fit model
weight_mpnp1_2.2_period15 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_mpnp1_2.2_period15
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
output1 <- read_cmdstan_csv(weight_mpnp1_2.2_period15$output_files()[1])
draws1_mpnp1_15 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_mpnp1_2.2_period15$output_files()[2])
draws2_mpnp1_15 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_mpnp1_2.2_period15$output_files()[3])
draws3_mpnp1_15 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_mpnp1_2.2_period15$output_files()[4])
draws4_mpnp1_15 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_mpnp1_15 <- rbind(draws1_mpnp1_15, draws2_mpnp1_15, draws3_mpnp1_15, draws4_mpnp1_15)

colnames(draws_mpnp1_15)[2:ncol(draws_mpnp1_15)] <- counts_df_period15$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_mpnp1_15), size = 30, replace = F)

### save data 
write_csv(draws_mpnp1_15, 'data_processed/mpnp1_bayesian_edgedistributions_a2.b2_period15_22.04.15.csv')
rm(draws1_mpnp1_15, draws2_mpnp1_15, draws3_mpnp1_15, draws4_mpnp1_15)

### build traceplots -- period 15 ####
# reset plot window
dev.off()

# traceplot
plot(draws_mpnp1_15[,plot_cols[1]], type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_mpnp1_15[,plot_cols[2]], col = 'tan')
lines(draws_mpnp1_15[,plot_cols[3]], col = 'orange')
lines(draws_mpnp1_15[,plot_cols[4]], col = 'green')
lines(draws_mpnp1_15[,plot_cols[5]], col = 'chocolate')
lines(draws_mpnp1_15[,plot_cols[6]], col = 'blue')
lines(draws_mpnp1_15[,plot_cols[7]], col = 'red')
lines(draws_mpnp1_15[,plot_cols[8]], col = 'seagreen')
lines(draws_mpnp1_15[,plot_cols[9]], col = 'purple')
lines(draws_mpnp1_15[,plot_cols[10]],col = 'magenta')
lines(draws_mpnp1_15[,plot_cols[11]],col = 'black')
lines(draws_mpnp1_15[,plot_cols[12]], col = 'tan')
lines(draws_mpnp1_15[,plot_cols[13]], col = 'orange')
lines(draws_mpnp1_15[,plot_cols[14]], col = 'green')
lines(draws_mpnp1_15[,plot_cols[15]], col = 'chocolate')
lines(draws_mpnp1_15[,plot_cols[16]], col = 'blue')
lines(draws_mpnp1_15[,plot_cols[17]], col = 'red')
lines(draws_mpnp1_15[,plot_cols[18]], col = 'seagreen')
lines(draws_mpnp1_15[,plot_cols[19]], col = 'purple')
lines(draws_mpnp1_15[,plot_cols[20]],col = 'magenta')
lines(draws_mpnp1_15[,plot_cols[21]],col = 'black')
lines(draws_mpnp1_15[,plot_cols[22]], col = 'tan')
lines(draws_mpnp1_15[,plot_cols[23]], col = 'orange')
lines(draws_mpnp1_15[,plot_cols[24]], col = 'green')
lines(draws_mpnp1_15[,plot_cols[25]], col = 'chocolate')
lines(draws_mpnp1_15[,plot_cols[26]], col = 'blue')
lines(draws_mpnp1_15[,plot_cols[27]], col = 'red')
lines(draws_mpnp1_15[,plot_cols[28]], col = 'seagreen')
lines(draws_mpnp1_15[,plot_cols[29]], col = 'purple')
lines(draws_mpnp1_15[,plot_cols[30]],col = 'magenta')

### density plots -- period 15 ####
dens(draws_mpnp1_15[,2], ylim = c(0,30),xlim = c(0,1), las = 1)
for(i in 1:30){  # looks like a poisson distribution - peaks just above 0 and then exponentiol decline. almost nothing above 0.5
  dens(add = T, draws_mpnp1_15[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 15 ####
# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_mpnp1_15[,2:ncol(draws_mpnp1_15)]),
                        min = rep(NA, ncol(draws_mpnp1_15)-1),
                        max = rep(NA, ncol(draws_mpnp1_15)-1),
                        mean = rep(NA, ncol(draws_mpnp1_15)-1),
                        median = rep(NA, ncol(draws_mpnp1_15)-1),
                        sd = rep(NA, ncol(draws_mpnp1_15)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_mpnp1_15[,i+1])
  summaries$max[i]    <- max(draws_mpnp1_15[,i+1])
  summaries$mean[i]   <- mean(draws_mpnp1_15[,i+1])
  summaries$median[i] <- median(draws_mpnp1_15[,i+1])
  summaries$sd[i]     <- sd(draws_mpnp1_15[,i+1])
}

# organise dem_class -- report age category based on age at start of period
plot_data_mpnp1_15 <- left_join(x = summaries, y = counts_df_period15, by = 'dyad')
head(plot_data_mpnp1_15)
plot_data_mpnp1_15$age_cat_1 <- ifelse(plot_data_mpnp1_15$age_start_1 < 6, 'C',
                                     ifelse(plot_data_mpnp1_15$age_start_1 < 11, 'J',
                                            ifelse(plot_data_mpnp1_15$age_start_1 > 19, 'A','P')))
plot_data_mpnp1_15$age_cat_2 <- ifelse(plot_data_mpnp1_15$age_start_2 < 6, 'C',
                                     ifelse(plot_data_mpnp1_15$age_start_2 < 11, 'J',
                                            ifelse(plot_data_mpnp1_15$age_start_2 > 19, 'A','P')))
plot_data_mpnp1_15$age_catnum_1 <- ifelse(plot_data_mpnp1_15$age_start_1 < 6, 1,
                                        ifelse(plot_data_mpnp1_15$age_start_1 < 11, 2,
                                               ifelse(plot_data_mpnp1_15$age_start_1 < 16, 3,
                                                      ifelse(plot_data_mpnp1_15$age_start_1 < 20, 4,
                                                             ifelse(plot_data_mpnp1_15$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_mpnp1_15$age_start_1 < 40, 6, 7))))))
plot_data_mpnp1_15$age_catnum_2 <- ifelse(plot_data_mpnp1_15$age_start_2 < 6, 1,
                                        ifelse(plot_data_mpnp1_15$age_start_2 < 11, 2,
                                               ifelse(plot_data_mpnp1_15$age_start_2 < 16, 3,
                                                      ifelse(plot_data_mpnp1_15$age_start_2 < 20, 4,
                                                             ifelse(plot_data_mpnp1_15$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_mpnp1_15$age_start_2 < 40, 6, 7))))))

plot_data_mpnp1_15$age_dyad <- ifelse(plot_data_mpnp1_15$age_catnum_1 >= plot_data_mpnp1_15$age_catnum_2,
                                    paste(plot_data_mpnp1_15$age_cat_1, plot_data_mpnp1_15$age_cat_2, sep = ''),
                                    paste(plot_data_mpnp1_15$age_cat_2, plot_data_mpnp1_15$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_mpnp1_15, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4, fill = c('blue','purple','blue','grey','purple','blue'))+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_mpnp1_15, aes(x = age_diff, y = mean))+
    geom_jitter(alpha = 0.2)+
    theme_classic()+
    theme(legend.position = 'none')+
    scale_x_continuous('age difference between dyad members',
                       expand = c(0.02,0))+
    scale_y_continuous('mean edge weight',
                       expand = c(0.02,0),
                       limits = c(0,1)))

# values for reporting
summary(summaries) # generally about the same as MOTNP
quantile(summaries$median, seq(0,1,length.out = 101))
quantile(summaries$mean,   seq(0,1,length.out = 101))
hist(summaries$median, breaks = 100, xlim = c(0,1), las = 1,
     xlab = 'median estimated dyad strength', main = '', col = 'skyblue')

m_sum <- plot_data_mpnp1_15[plot_data_mpnp1_15$age_dyad == 'AA' | plot_data_mpnp1_15$age_dyad == 'AP' | plot_data_mpnp1_15$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))
hist(m_sum$median, breaks = 100, xlim = c(0,1), las = 1,
     xlab = 'median estimated dyad strength', main = '', col = 'skyblue')

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 15 ################
head(summaries)
length(unique(plot_data_mpnp1_15$id_1))+1 # number of individuals = 181

par(mai = c(0.1,0.1,0.1,0.1))

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period15$id_1))+1,
                         NROW(unique(counts_df_period15$id_2))+1,
                         NROW(draws_mpnp1_15)),
                    dimnames = list(c(unique(counts_df_period15$id_1),
                                      unique(counts_df_period15$id_2)[length(unique(counts_df_period15$id_2))]),
                                    c(unique(counts_df_period15$id_1),
                                      unique(counts_df_period15$id_2)[length(unique(counts_df_period15$id_2))]),
                                    NULL))
N <- nrow(counts_df_period15)

for (i in 1:N) {
  dyad_row <- counts_df_period15[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_mpnp1_15[, i+1]
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
ids15 <- c(unique(counts_df_period15$id_1), unique(counts_df_period15$id_2)[length(unique(counts_df_period15$id_2))])
males15 <- males[,c(21,6,9,22:53)]
males15 <- males15 %>% dplyr::filter(id %in% ids15)
males15

# create variables for different degrees of node connectedness
males15$degree_0.1 <- NA
males15$degree_0.2 <- NA
males15$degree_0.3 <- NA
males15$degree_0.4 <- NA
males15$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males15)){
  rows <- summaries[summaries$id_1 == males15$id[i] | summaries$id_2 == males15$id[i],]
  males15$degree_0.1[i] <- length(which(rows$median > 0.1))
  males15$degree_0.2[i] <- length(which(rows$median > 0.2))
  males15$degree_0.3[i] <- length(which(rows$median > 0.3))
  males15$degree_0.4[i] <- length(which(rows$median > 0.4))
  males15$degree_0.5[i] <- length(which(rows$median > 0.5))
}

which(males15$degree_0.1 < males15$degree_0.2)
which(males15$degree_0.2 < males15$degree_0.3)
which(males15$degree_0.3 < males15$degree_0.4)
which(males15$degree_0.4 < males15$degree_0.5)

# age variable
males15$age <- lubridate::year(periods$period_start[periods$period == 15]) - males15$byr
summary(males15$age)
males15$age_class <- ifelse(males15$age < 10, 2,
                            ifelse(males15$age < 15, 3,
                                   ifelse(males15$age < 20, 4,
                                          ifelse(males15$age < 25, 5,
                                                 ifelse(males15$age < 40, 6, 7)))))

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
     vertex.label = males15$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males15$age_class == 7,'seagreen4',
                          ifelse(males15$age_class == 6,'seagreen3',
                                 ifelse(males15$age_class == 5,'seagreen2',
                                        ifelse(males15$age_class == 4,'steelblue3',
                                               ifelse(males15$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(adj_mid < 0.15,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(adj_mid < 0.15,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males15$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males15$age_class == 7,'seagreen4',
                          ifelse(males15$age_class == 6,'seagreen3',
                                 ifelse(males15$age_class == 5,'seagreen2',
                                        ifelse(males15$age_class == 4,'steelblue3',
                                               ifelse(males15$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.2 -- period 15 ####
g_mid_0.2 <- delete.vertices(graph = g_mid, v = males15$id[which(males15$degree_0.2 == 0)])
g_rng_0.2 <- delete.vertices(graph = g_rng, v = males15$id[which(males15$degree_0.2 == 0)])

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
     vertex.color= ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 3,'steelblue1',
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
     vertex.size = males15$p15[which(males15$degree_0.2 != 0)],
     vertex.label = males15$id[which(males15$degree_0.2 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.2, add = TRUE)

### save summary data -- period 15 ####
summaries$period <- 15
summaries <- summaries[,c(1,4:9)]
dyad_period_weights <- left_join(x = dyad_period_weights, y = summaries,
                                 by = c('dyad','period'))
colnames(dyad_period_weights)[36:40] <- c('min','max','mean','median','sd')
for(i in 1:nrow(dyad_period_weights)){
  dyad_period_weights$min[i] <- ifelse(is.na(dyad_period_weights$min[i]) == TRUE,
                                       ifelse(is.na(dyad_period_weights$min.x[i]) == FALSE,
                                              dyad_period_weights$min.x[i], NA),
                                       dyad_period_weights$min[i])
  dyad_period_weights$max[i] <- ifelse(is.na(dyad_period_weights$max[i]) == TRUE,
                                       ifelse(is.na(dyad_period_weights$max.x[i]) == FALSE,
                                              dyad_period_weights$max.x[i], NA),
                                       dyad_period_weights$max[i])
  dyad_period_weights$mean[i] <- ifelse(is.na(dyad_period_weights$mean[i]) == TRUE,
                                        ifelse(is.na(dyad_period_weights$mean.x[i]) == FALSE,
                                               dyad_period_weights$mean.x[i], NA),
                                        dyad_period_weights$mean[i])
  dyad_period_weights$median[i] <- ifelse(is.na(dyad_period_weights$median[i]) == TRUE,
                                          ifelse(is.na(dyad_period_weights$median.x[i]) == FALSE,
                                                 dyad_period_weights$median.x[i], NA),
                                          dyad_period_weights$median[i])
  dyad_period_weights$sd[i] <- ifelse(is.na(dyad_period_weights$sd[i]) == TRUE,
                                      ifelse(is.na(dyad_period_weights$sd.x[i]) == FALSE,
                                             dyad_period_weights$sd.x[i], NA),
                                      dyad_period_weights$sd[i])
}
dyad_period_weights <- dyad_period_weights[,c(1:30,36:40)]

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.2,
   counts_df_period15, draws_mpnp1_15, dyad_row, g_mid,g_mid_0.2, g_rng, g_rng_0.2, males15, plot_data_mpnp1_15, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids15, N, plot_cols)

################ 7.28) Run model on real standardised data -- period 28 ################
### create data list
counts_df_period28 <- counts_df1[counts_df1$period == 28,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period28),          # total number of times one or other of the dyad was observed
  together = counts_df_period28$event_count,    # count number of sightings seen together
  apart    = counts_df_period28$apart,          # count number of sightings seen apart
  period   = counts_df_period28$period)         # which period it's within

### Fit model
weight_mpnp1_2.2_period28 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_mpnp1_2.2_period28
# variable       mean     median     sd    mad         q5        q95 rhat ess_bulk ess_tail
#lp__      -124881.32 -124878.00 109.57 111.19 -125063.00 -124701.00 1.00     1361     2334
#weight[1]       0.14       0.13   0.09   0.09       0.03       0.33 1.00     3988     2130
#weight[2]       0.12       0.11   0.07   0.07       0.02       0.26 1.00     3463     1692
#weight[3]       0.11       0.10   0.06   0.06       0.03       0.22 1.00     4186     2399
#weight[4]       0.10       0.09   0.07   0.06       0.02       0.23 1.00     3647     2072
#weight[5]       0.15       0.13   0.09   0.09       0.03       0.32 1.00     3699     2219
#weight[6]       0.18       0.16   0.09   0.09       0.05       0.34 1.00     4891     2856
#weight[7]       0.22       0.20   0.11   0.11       0.06       0.42 1.00     3636     2532
#weight[8]       0.20       0.19   0.10   0.10       0.06       0.39 1.00     4324     2476
#weight[9]       0.09       0.08   0.06   0.06       0.02       0.21 1.00     3740     1912
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_mpnp1_2.2_period28$output_files()[1])
draws1_mpnp1_28 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_mpnp1_2.2_period28$output_files()[2])
draws2_mpnp1_28 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_mpnp1_2.2_period28$output_files()[3])
draws3_mpnp1_28 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_mpnp1_2.2_period28$output_files()[4])
draws4_mpnp1_28 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_mpnp1_28 <- rbind(draws1_mpnp1_28, draws2_mpnp1_28, draws3_mpnp1_28, draws4_mpnp1_28)

colnames(draws_mpnp1_28)[2:ncol(draws_mpnp1_28)] <- counts_df_period28$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_mpnp1_28), size = 30, replace = F)

### save data 
write_csv(draws_mpnp1_28, 'data_processed/mpnp1_bayesian_edgedistributions_a2.b2_period28_22.04.15.csv')
rm(draws1_mpnp1_28, draws2_mpnp1_28, draws3_mpnp1_28, draws4_mpnp1_28)

### build traceplots  -- period 28 ####
# reset plot window
dev.off()

# traceplot
plot(draws_mpnp1_28[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_mpnp1_28[,plot_cols[2]], col = 'tan')
lines(draws_mpnp1_28[,plot_cols[3]], col = 'orange')
lines(draws_mpnp1_28[,plot_cols[4]], col = 'green')
lines(draws_mpnp1_28[,plot_cols[5]], col = 'chocolate')
lines(draws_mpnp1_28[,plot_cols[6]], col = 'blue')
lines(draws_mpnp1_28[,plot_cols[7]], col = 'red')
lines(draws_mpnp1_28[,plot_cols[8]], col = 'seagreen')
lines(draws_mpnp1_28[,plot_cols[9]], col = 'purple')
lines(draws_mpnp1_28[,plot_cols[10]],col = 'magenta')
lines(draws_mpnp1_28[,plot_cols[11]],col = 'black')
lines(draws_mpnp1_28[,plot_cols[12]], col = 'tan')
lines(draws_mpnp1_28[,plot_cols[13]], col = 'orange')
lines(draws_mpnp1_28[,plot_cols[14]], col = 'green')
lines(draws_mpnp1_28[,plot_cols[15]], col = 'chocolate')
lines(draws_mpnp1_28[,plot_cols[16]], col = 'blue')
lines(draws_mpnp1_28[,plot_cols[17]], col = 'red')
lines(draws_mpnp1_28[,plot_cols[18]], col = 'seagreen')
lines(draws_mpnp1_28[,plot_cols[19]], col = 'purple')
lines(draws_mpnp1_28[,plot_cols[20]],col = 'magenta')
lines(draws_mpnp1_28[,plot_cols[21]],col = 'black')
lines(draws_mpnp1_28[,plot_cols[22]], col = 'tan')
lines(draws_mpnp1_28[,plot_cols[23]], col = 'orange')
lines(draws_mpnp1_28[,plot_cols[24]], col = 'green')
lines(draws_mpnp1_28[,plot_cols[25]], col = 'chocolate')
lines(draws_mpnp1_28[,plot_cols[26]], col = 'blue')
lines(draws_mpnp1_28[,plot_cols[27]], col = 'red')
lines(draws_mpnp1_28[,plot_cols[28]], col = 'seagreen')
lines(draws_mpnp1_28[,plot_cols[29]], col = 'purple')
lines(draws_mpnp1_28[,plot_cols[30]],col = 'magenta')

### density plots  -- period 28 ####
dens(draws_mpnp1_28[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:30){  # looks like a poisson distribution - peaks just above 0 and then exponentiol decline. almost nothing above 0.5
  dens(add = T, draws_mpnp1_28[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 28 ####
# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_mpnp1_28[,2:ncol(draws_mpnp1_28)]),
                        min = rep(NA, ncol(draws_mpnp1_28)-1),
                        max = rep(NA, ncol(draws_mpnp1_28)-1),
                        mean = rep(NA, ncol(draws_mpnp1_28)-1),
                        median = rep(NA, ncol(draws_mpnp1_28)-1),
                        sd = rep(NA, ncol(draws_mpnp1_28)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_mpnp1_28[,i+1])
  summaries$max[i]    <- max(draws_mpnp1_28[,i+1])
  summaries$mean[i]   <- mean(draws_mpnp1_28[,i+1])
  summaries$median[i] <- median(draws_mpnp1_28[,i+1])
  summaries$sd[i]     <- sd(draws_mpnp1_28[,i+1])
}

# organise dem_class -- report age category based on age at start of period
plot_data_mpnp1_28 <- left_join(x = summaries, y = counts_df_period28, by = 'dyad')
head(plot_data_mpnp1_28)
plot_data_mpnp1_28$age_cat_1 <- ifelse(plot_data_mpnp1_28$age_start_1 < 6, 'C',
                                     ifelse(plot_data_mpnp1_28$age_start_1 < 11, 'J',
                                            ifelse(plot_data_mpnp1_28$age_start_1 > 19, 'A','P')))
plot_data_mpnp1_28$age_cat_2 <- ifelse(plot_data_mpnp1_28$age_start_2 < 6, 'C',
                                     ifelse(plot_data_mpnp1_28$age_start_2 < 11, 'J',
                                            ifelse(plot_data_mpnp1_28$age_start_2 > 19, 'A','P')))
plot_data_mpnp1_28$age_catnum_1 <- ifelse(plot_data_mpnp1_28$age_start_1 < 6, 1,
                                        ifelse(plot_data_mpnp1_28$age_start_1 < 11, 2,
                                               ifelse(plot_data_mpnp1_28$age_start_1 < 16, 3,
                                                      ifelse(plot_data_mpnp1_28$age_start_1 < 20, 4,
                                                             ifelse(plot_data_mpnp1_28$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_mpnp1_28$age_start_1 < 40, 6, 7))))))
plot_data_mpnp1_28$age_catnum_2 <- ifelse(plot_data_mpnp1_28$age_start_2 < 6, 1,
                                        ifelse(plot_data_mpnp1_28$age_start_2 < 11, 2,
                                               ifelse(plot_data_mpnp1_28$age_start_2 < 16, 3,
                                                      ifelse(plot_data_mpnp1_28$age_start_2 < 20, 4,
                                                             ifelse(plot_data_mpnp1_28$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_mpnp1_28$age_start_2 < 40, 6, 7))))))

plot_data_mpnp1_28$age_dyad <- ifelse(plot_data_mpnp1_28$age_catnum_1 >= plot_data_mpnp1_28$age_catnum_2,
                                    paste(plot_data_mpnp1_28$age_cat_1, plot_data_mpnp1_28$age_cat_2, sep = ''),
                                    paste(plot_data_mpnp1_28$age_cat_2, plot_data_mpnp1_28$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_mpnp1_28, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4, fill = c('blue','purple','blue','grey','purple','blue'))+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_mpnp1_28, aes(x = age_diff, y = mean))+
    geom_jitter(alpha = 0.2)+
    theme_classic()+
    theme(legend.position = 'none')+
    scale_x_continuous('age difference between dyad members',
                       expand = c(0.02,0))+
    scale_y_continuous('mean edge weight',
                       expand = c(0.02,0),
                       limits = c(0,1)))

# values for reporting
summary(summaries) # generally about the same as MOTNP
quantile(summaries$median, seq(0,1,length.out = 101))
quantile(summaries$mean,   seq(0,1,length.out = 101))
hist(summaries$median, breaks = 100, xlim = c(0,1), las = 1,
     xlab = 'median estimated dyad strength', main = '', col = 'skyblue')

m_sum <- plot_data_mpnp1_28[plot_data_mpnp1_28$age_dyad == 'AA' | plot_data_mpnp1_28$age_dyad == 'AP' | plot_data_mpnp1_28$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 28 ################
head(summaries)
length(unique(plot_data_mpnp1_28$id_1))+1 # number of individuals = 198

par(mai = c(0.1,0.1,0.1,0.1))

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period28$id_1))+1,
                         NROW(unique(counts_df_period28$id_2))+1,
                         NROW(draws_mpnp1_28)),
                    dimnames = list(c(unique(counts_df_period28$id_1),
                                      unique(counts_df_period28$id_2)[length(unique(counts_df_period28$id_2))]),
                                    c(unique(counts_df_period28$id_1),
                                      unique(counts_df_period28$id_2)[length(unique(counts_df_period28$id_2))]),
                                    NULL))
N <- nrow(counts_df_period28)

for (i in 1:N) {
  dyad_row <- counts_df_period28[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_mpnp1_28[, i+1]
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
ids28 <- c(unique(counts_df_period28$id_1), unique(counts_df_period28$id_2)[length(unique(counts_df_period28$id_2))])
males28 <- males[,c(21,6,9,22:53)]
males28 <- males28 %>% dplyr::filter(id %in% ids28)
males28

# create variables for different degrees of node connectedness
males28$degree_0.1 <- NA
males28$degree_0.2 <- NA
males28$degree_0.3 <- NA
males28$degree_0.4 <- NA
males28$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males28)){
  rows <- summaries[summaries$id_1 == males28$id[i] | summaries$id_2 == males28$id[i],]
  males28$degree_0.1[i] <- length(which(rows$median > 0.1))
  males28$degree_0.2[i] <- length(which(rows$median > 0.2))
  males28$degree_0.3[i] <- length(which(rows$median > 0.3))
  males28$degree_0.4[i] <- length(which(rows$median > 0.4))
  males28$degree_0.5[i] <- length(which(rows$median > 0.5))
}

which(males28$degree_0.1 < males28$degree_0.2)
which(males28$degree_0.2 < males28$degree_0.3)
which(males28$degree_0.3 < males28$degree_0.4)
which(males28$degree_0.4 < males28$degree_0.5)

# age variable
males28$age <- lubridate::year(periods$period_start[periods$period == 28]) - males28$byr
summary(males28$age)
males28$age_class <- ifelse(males28$age < 10, 2,
                            ifelse(males28$age < 15, 3,
                                   ifelse(males28$age < 20, 4,
                                          ifelse(males28$age < 25, 5,
                                                 ifelse(males28$age < 40, 6, 7)))))

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
     vertex.label = males28$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males28$age_class == 7,'seagreen4',
                          ifelse(males28$age_class == 6,'seagreen3',
                                 ifelse(males28$age_class == 5,'seagreen2',
                                        ifelse(males28$age_class == 4,'steelblue3',
                                               ifelse(males28$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(adj_mid < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(adj_mid < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males28$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males28$age_class == 7,'seagreen4',
                          ifelse(males28$age_class == 6,'seagreen3',
                                 ifelse(males28$age_class == 5,'seagreen2',
                                        ifelse(males28$age_class == 4,'steelblue3',
                                               ifelse(males28$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.2 -- period 28 ####
g_mid_0.2 <- delete.vertices(graph = g_mid, v = males28$id[which(males28$degree_0.2 == 0)])
g_rng_0.2 <- delete.vertices(graph = g_rng, v = males28$id[which(males28$degree_0.2 == 0)])

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
     vertex.color= ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 3,'steelblue1',
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
     vertex.size = males28$p28[which(males28$degree_0.2 != 0)],
     vertex.label = males28$id[which(males28$degree_0.2 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.2, add = TRUE)

### save summary data -- period 28 ####
summaries$period <- 28
summaries <- summaries[,c(1,4:9)]
dyad_period_weights <- left_join(x = dyad_period_weights, y = summaries,
                                 by = c('dyad','period'))
colnames(dyad_period_weights)[36:40] <- c('min','max','mean','median','sd')
for(i in 1:nrow(dyad_period_weights)){
  dyad_period_weights$min[i] <- ifelse(is.na(dyad_period_weights$min[i]) == TRUE,
                                       ifelse(is.na(dyad_period_weights$min.x[i]) == FALSE,
                                              dyad_period_weights$min.x[i], NA),
                                       dyad_period_weights$min[i])
  dyad_period_weights$max[i] <- ifelse(is.na(dyad_period_weights$max[i]) == TRUE,
                                       ifelse(is.na(dyad_period_weights$max.x[i]) == FALSE,
                                              dyad_period_weights$max.x[i], NA),
                                       dyad_period_weights$max[i])
  dyad_period_weights$mean[i] <- ifelse(is.na(dyad_period_weights$mean[i]) == TRUE,
                                        ifelse(is.na(dyad_period_weights$mean.x[i]) == FALSE,
                                               dyad_period_weights$mean.x[i], NA),
                                        dyad_period_weights$mean[i])
  dyad_period_weights$median[i] <- ifelse(is.na(dyad_period_weights$median[i]) == TRUE,
                                          ifelse(is.na(dyad_period_weights$median.x[i]) == FALSE,
                                                 dyad_period_weights$median.x[i], NA),
                                          dyad_period_weights$median[i])
  dyad_period_weights$sd[i] <- ifelse(is.na(dyad_period_weights$sd[i]) == TRUE,
                                      ifelse(is.na(dyad_period_weights$sd.x[i]) == FALSE,
                                             dyad_period_weights$sd.x[i], NA),
                                      dyad_period_weights$sd[i])
}
dyad_period_weights <- dyad_period_weights[,c(1:30,36:40)]

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.2,
   counts_df_period28, draws_mpnp1_28, dyad_row, g_mid,g_mid_0.2, g_rng, g_rng_0.2, males28, plot_data_mpnp1_28, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids28, N, plot_cols)
