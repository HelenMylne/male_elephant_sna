#### Bayesian analysis of ATE data ####
# Script to process association data from Amboseli National Park, Kenya
# Data collected: 1972-2021 by Amboseli Trust for Elephants
# Data supplied by: Vicki Fishlock (March 2022) and Phyllis Lee (December 2021)
# Data input: raw data provided by ATE processed in scripts 22.03.20_anp_dataprocessing1.R and 22.03.22_anp_dataprocessing2.R

#### Set up ####
# set working directory for working on research cluster
#setwd(EleSNA_researchcluster/)

# install packages
# install.packages('tidyverse')
# install.packages('dplyr')
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages(c("StanHeaders","rstan"),type="source")
# install.packages("remotes") ; remotes::install_github("stan-dev/cmdstanr") ## OR USE # install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages('rethinking')
# install.packages('igraph')
# install.packages('dagitty')
# install.packages('janitor')
# install.packages('lubridate')
# install.packages('hms')
# install.packages('readxl')

# load packages
library(tidyverse, lib.loc = 'packages/')   # library(tidyverse)
library(dplyr, lib.loc = 'packages/')       # library(dplyr)
#library(rstan, lib.loc = 'packages/')      # library(rstan)
#library(cmdstanr, lib.loc = 'packages/')    # library(cmdstanr)
# check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
# install_cmdstan(cores = 2)
#library(rethinking, lib.loc = 'packages/')  # library(rethinking)
library(igraph, lib.loc = 'packages/')      # library(igraph)
#library(dagitty, lib.loc = 'packages/')     # library(dagitty)
library(janitor, lib.loc = 'packages/')     # library(janitor)
library(lubridate, lib.loc = 'packages/')   # library(lubridate)
library(hms, lib.loc = 'packages/')         # library(hms)
#library(readxl, lib.loc = 'packages/')      # library(readxl)

# set seed
set.seed(12345)

# create file of output graphs
pdf('data_processed/anp_edgeweights_2.2_allperiods_networkplots_22.07.04.pdf', width = 10, height = 10)

#### read in data ################
counts_df_non0 <- read_csv('data_processed/anp_bayesian_pairwiseevents_22.04.06.csv')

ate <- readxl::read_excel(path = 'data_raw/Raw_ATE_Sightings_Fishlock220224.xlsx', sheet = 1) %>% janitor::clean_names()
ate$id <- paste('M',ate$casename, sep = '')
ate$node_id <- as.integer(as.factor(ate$casename))
ate$obs_date <- lubridate::as_date(ate$obs_date)
ate <- separate(ate, obs_time, into = c('wrong_date','correct_time'), remove = F, sep = ' ')
ate$correct_time_hms <- hms::as_hms(ate$correct_time)
ate$corrected_time <- lubridate::hour(ate$correct_time_hms)*60*60 + lubridate::minute(ate$correct_time_hms) + lubridate::second(ate$correct_time_hms)
lu <- function(x) { length(unique(x)) }
ate_nums <- tapply(X = ate$obs_num, INDEX = ate$obs_date, FUN = lu )
ate$obs_num <- ifelse(ate$obs_num == '0','00', ate$obs_num)
ate$obs_num <- ifelse(ate$obs_num == '0a','0A', ate$obs_num)
ate$obs_num <- ifelse(ate$obs_num == '0b','0B', ate$obs_num)
ate$obs_num <- ifelse(ate$obs_num == '1','01', ate$obs_num)
ate$obs_num_std <- NA
for(i in 1:length(ate)){
  date_row <- ate[ate$obs_date == ate$obs_date[i],]
  date_row$obs_num_std <- as.integer(as.factor(sort(date_row$obs_num)))
  ate$obs_num_std[i] <- date_row$obs_num_std[which(date_row$obs_id == ate$obs_id[i])[1]]
}
ate <- ate[,c(1:3,25,26,4,27,28,8,29,9:24)]
head(ate)
rm(date_row, ate_nums, lu, i)

sightings <- ate[,c('obs_id','obs_date','correct_time_hms','obs_num_std','grid_code')] %>% distinct()
sightings <- sightings[c(1:12,15:24176),]

periods <- seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 32)

males <- readxl::read_excel('data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()
males$id <- paste0('M',males$casename)
ids <- sort(unique(ate$id))
males <- males %>% dplyr::filter(id %in% ids)

# counts per male per time period
males$count_all <- NA
males$p1 <- NA ; males$p2 <- NA ; males$p3 <- NA ; males$p4 <- NA ; males$p5 <- NA ; males$p6 <- NA ; males$p7 <- NA
males$p8 <- NA ; males$p9 <- NA ; males$p10 <- NA ; males$p11 <- NA ; males$p12 <- NA ; males$p13 <- NA
males$p14 <- NA ; males$p15 <- NA ; males$p16 <- NA ; males$p17 <- NA ; males$p18 <- NA ; males$p19 <- NA
males$p20 <- NA ; males$p21 <- NA ; males$p22 <- NA ; males$p23 <- NA ; males$p24 <- NA ; males$p25 <- NA
males$p26 <- NA ; males$p27 <- NA ; males$p28 <- NA ; males$p29 <- NA ; males$p30 <- NA ; males$p31 <- NA
for(i in 1:nrow(males)){
  counts <- ate[ate$id == males$id[i],]
  males$count_all[i] <- nrow(counts)
  males$p1[i] <- length(which(counts$obs_date >= periods[1] & counts$obs_date < periods[2]))
  males$p2[i] <- length(which(counts$obs_date >= periods[2] & counts$obs_date < periods[3]))
  males$p3[i] <- length(which(counts$obs_date >= periods[3] & counts$obs_date < periods[4]))
  males$p4[i] <- length(which(counts$obs_date >= periods[4] & counts$obs_date < periods[5]))
  males$p5[i] <- length(which(counts$obs_date >= periods[5] & counts$obs_date < periods[6]))
  males$p6[i] <- length(which(counts$obs_date >= periods[6] & counts$obs_date < periods[7]))
  males$p7[i] <- length(which(counts$obs_date >= periods[7] & counts$obs_date < periods[8]))
  males$p8[i] <- length(which(counts$obs_date >= periods[8] & counts$obs_date < periods[9]))
  males$p9[i] <- length(which(counts$obs_date >= periods[9] & counts$obs_date < periods[10]))
  males$p10[i] <- length(which(counts$obs_date >= periods[10] & counts$obs_date < periods[11]))
  males$p11[i] <- length(which(counts$obs_date >= periods[11] & counts$obs_date < periods[12]))
  males$p12[i] <- length(which(counts$obs_date >= periods[12] & counts$obs_date < periods[13]))
  males$p13[i] <- length(which(counts$obs_date >= periods[13] & counts$obs_date < periods[14]))
  males$p14[i] <- length(which(counts$obs_date >= periods[14] & counts$obs_date < periods[15]))
  males$p15[i] <- length(which(counts$obs_date >= periods[15] & counts$obs_date < periods[16]))
  males$p16[i] <- length(which(counts$obs_date >= periods[16] & counts$obs_date < periods[17]))
  males$p17[i] <- length(which(counts$obs_date >= periods[17] & counts$obs_date < periods[18]))
  males$p18[i] <- length(which(counts$obs_date >= periods[18] & counts$obs_date < periods[19]))
  males$p19[i] <- length(which(counts$obs_date >= periods[19] & counts$obs_date < periods[20]))
  males$p20[i] <- length(which(counts$obs_date >= periods[20] & counts$obs_date < periods[21]))
  males$p21[i] <- length(which(counts$obs_date >= periods[21] & counts$obs_date < periods[22]))
  males$p22[i] <- length(which(counts$obs_date >= periods[22] & counts$obs_date < periods[23]))
  males$p23[i] <- length(which(counts$obs_date >= periods[23] & counts$obs_date < periods[24]))
  males$p24[i] <- length(which(counts$obs_date >= periods[24] & counts$obs_date < periods[25]))
  males$p25[i] <- length(which(counts$obs_date >= periods[25] & counts$obs_date < periods[26]))
  males$p26[i] <- length(which(counts$obs_date >= periods[26] & counts$obs_date < periods[27]))
  males$p27[i] <- length(which(counts$obs_date >= periods[27] & counts$obs_date < periods[28]))
  males$p28[i] <- length(which(counts$obs_date >= periods[28] & counts$obs_date < periods[29]))
  males$p29[i] <- length(which(counts$obs_date >= periods[29] & counts$obs_date < periods[30]))
  males$p30[i] <- length(which(counts$obs_date >= periods[30] & counts$obs_date < periods[31]))
  males$p31[i] <- length(which(counts$obs_date >= periods[31] & counts$obs_date < periods[32]))
}
males <- males[,c(21,6,9,22:53)]

periods <- data.frame(period = 1:31,
                      period_start = periods[1:31],
                      period_end = periods[2:32])

################ Period 1 ################
print(paste0('Period 1 started at ', Sys.time()))

draws_anp_1 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period1_22.04.15.csv')
saveRDS(draws_anp_1, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period1_22.04.15.rds')
draws_anp_1 <- data.matrix(draws_anp_1)

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_1[,2:ncol(draws_anp_1)]),
                        min = rep(NA, ncol(draws_anp_1)-1),
                        max = rep(NA, ncol(draws_anp_1)-1),
                        mean = rep(NA, ncol(draws_anp_1)-1),
                        median = rep(NA, ncol(draws_anp_1)-1),
                        sd = rep(NA, ncol(draws_anp_1)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_1[,i+1])
  summaries$max[i]    <- max(draws_anp_1[,i+1])
  summaries$mean[i]   <- mean(draws_anp_1[,i+1])
  summaries$median[i] <- median(draws_anp_1[,i+1])
  summaries$sd[i]     <- sd(draws_anp_1[,i+1])
}

### create network plots
counts_df_period1 <- counts_df_non0[counts_df_non0$period == 1,]

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period1$id_1))+1,
                         NROW(unique(counts_df_period1$id_2))+1,
                         NROW(draws_anp_1)),
                    dimnames = list(c(unique(counts_df_period1$id_1),
                                      unique(counts_df_period1$id_2)[length(unique(counts_df_period1$id_2))]),
                                    c(unique(counts_df_period1$id_1),
                                      unique(counts_df_period1$id_2)[length(unique(counts_df_period1$id_2))]),
                                    NULL))
N <- nrow(counts_df_period1)

for (i in 1:N) {
  dyad_row <- counts_df_period1[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_1[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids1 <- c(unique(counts_df_period1$id_1), unique(counts_df_period1$id_2)[length(unique(counts_df_period1$id_2))])
males1 <- males %>% dplyr::filter(id %in% ids1)

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

# age variable
males1$age <- lubridate::year(periods$period_start[periods$period == 1]) - males1$byr
summary(males1$age)
males1$age_class <- ifelse(males1$age < 10, 2,
                           ifelse(males1$age < 15, 3,
                                  ifelse(males1$age < 20, 4,
                                         ifelse(males1$age < 25, 5,
                                                ifelse(males1$age < 40, 6, 7)))))

g_mid_0.2 <- delete.vertices(graph = g_mid, v = males1$id[which(males1$degree_0.2 == 0)])
g_rng_0.2 <- delete.vertices(graph = g_rng, v = males1$id[which(males1$degree_0.2 == 0)])

set.seed(3)
coords_0.2 <- layout_nicely(g_mid_0.2)

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

# clear environment
rm(summaries, draws_anp_1, counts_df_period1, males1, g_mid, g_rng, g_mid_0.2, g_rng_0.2, ids1)

################ Period 2 ################
print(paste0('Period 2 started at ', Sys.time()))

draws_anp_2 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period2_22.04.27.csv')
saveRDS(draws_anp_2, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period2_22.04.27.rds')
draws_anp_2 <- data.matrix(draws_anp_2)

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_2[,2:ncol(draws_anp_2)]),
                        min = rep(NA, ncol(draws_anp_2)-1),
                        max = rep(NA, ncol(draws_anp_2)-1),
                        mean = rep(NA, ncol(draws_anp_2)-1),
                        median = rep(NA, ncol(draws_anp_2)-1),
                        sd = rep(NA, ncol(draws_anp_2)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_2[,i+1])
  summaries$max[i]    <- max(draws_anp_2[,i+1])
  summaries$mean[i]   <- mean(draws_anp_2[,i+1])
  summaries$median[i] <- median(draws_anp_2[,i+1])
  summaries$sd[i]     <- sd(draws_anp_2[,i+1])
}

counts_df_period2 <- counts_df_non0[counts_df_non0$period == 2,]

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period2$id_1))+1,
                         NROW(unique(counts_df_period2$id_2))+1,
                         NROW(draws_anp_2)),
                    dimnames = list(c(unique(counts_df_period2$id_1),
                                      unique(counts_df_period2$id_2)[length(unique(counts_df_period2$id_2))]),
                                    c(unique(counts_df_period2$id_1),
                                      unique(counts_df_period2$id_2)[length(unique(counts_df_period2$id_2))]),
                                    NULL))
N <- nrow(counts_df_period2)

for (i in 1:N) {
  dyad_row <- counts_df_period2[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_2[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids2 <- c(unique(counts_df_period2$id_1), unique(counts_df_period2$id_2)[length(unique(counts_df_period2$id_2))])
males2 <- males %>% dplyr::filter(id %in% ids2)

# create variables for different degrees of node connectedness
males2$degree_0.1 <- NA
males2$degree_0.2 <- NA
males2$degree_0.3 <- NA
males2$degree_0.4 <- NA
males2$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males2)){
  rows <- summaries[summaries$id_1 == males2$id[i] | summaries$id_2 == males2$id[i],]
  males2$degree_0.1[i] <- length(which(rows$median > 0.1))
  males2$degree_0.2[i] <- length(which(rows$median > 0.2))
  males2$degree_0.3[i] <- length(which(rows$median > 0.3))
  males2$degree_0.4[i] <- length(which(rows$median > 0.4))
  males2$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males2$age <- lubridate::year(periods$period_start[periods$period == 2]) - males2$byr
summary(males2$age)
males2$age_class <- ifelse(males2$age < 5, 1,
                           ifelse(males2$age < 10, 2,
                                  ifelse(males2$age < 15, 3,
                                         ifelse(males2$age < 20, 4,
                                                ifelse(males2$age < 25, 5,
                                                       ifelse(males2$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males2$id[which(males2$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males2$id[which(males2$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males2$p2[which(males2$degree_0.3 != 0)],
     vertex.label = males2$id[which(males2$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males2[which(males2$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males2[which(males2$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males2[which(males2$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males2[which(males2$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males2[which(males2$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 3 ################
print(paste0('Period 3 started at ', Sys.time()))

draws_anp_3 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period3_22.04.27.csv')
saveRDS(draws_anp_3, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period3_22.04.27.rds')
draws_anp_3 <- data.matrix(draws_anp_3)

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_3[,2:ncol(draws_anp_3)]),
                        min = rep(NA, ncol(draws_anp_3)-1),
                        max = rep(NA, ncol(draws_anp_3)-1),
                        mean = rep(NA, ncol(draws_anp_3)-1),
                        median = rep(NA, ncol(draws_anp_3)-1),
                        sd = rep(NA, ncol(draws_anp_3)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_3[,i+1])
  summaries$max[i]    <- max(draws_anp_3[,i+1])
  summaries$mean[i]   <- mean(draws_anp_3[,i+1])
  summaries$median[i] <- median(draws_anp_3[,i+1])
  summaries$sd[i]     <- sd(draws_anp_3[,i+1])
}

counts_df_period3 <- counts_df_non0[counts_df_non0$period == 3,]

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period3$id_1))+1,
                         NROW(unique(counts_df_period3$id_2))+1,
                         NROW(draws_anp_3)),
                    dimnames = list(c(unique(counts_df_period3$id_1),
                                      unique(counts_df_period3$id_2)[length(unique(counts_df_period3$id_2))]),
                                    c(unique(counts_df_period3$id_1),
                                      unique(counts_df_period3$id_2)[length(unique(counts_df_period3$id_2))]),
                                    NULL))
N <- nrow(counts_df_period3)

for (i in 1:N) {
  dyad_row <- counts_df_period3[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_3[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids3 <- c(unique(counts_df_period3$id_1), unique(counts_df_period3$id_2)[length(unique(counts_df_period3$id_2))])
males3 <- males %>% dplyr::filter(id %in% ids3)

# create variables for different degrees of node connectedness
males3$degree_0.1 <- NA
males3$degree_0.2 <- NA
males3$degree_0.3 <- NA
males3$degree_0.4 <- NA
males3$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males3)){
  rows <- summaries[summaries$id_1 == males3$id[i] | summaries$id_2 == males3$id[i],]
  males3$degree_0.1[i] <- length(which(rows$median > 0.1))
  males3$degree_0.2[i] <- length(which(rows$median > 0.2))
  males3$degree_0.3[i] <- length(which(rows$median > 0.3))
  males3$degree_0.4[i] <- length(which(rows$median > 0.4))
  males3$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males3$age <- lubridate::year(periods$period_start[periods$period == 3]) - males3$byr
summary(males3$age)
males3$age_class <- ifelse(males3$age < 5, 1,
                           ifelse(males3$age < 10, 3,
                                  ifelse(males3$age < 15, 3,
                                         ifelse(males3$age < 20, 4,
                                                ifelse(males3$age < 25, 5,
                                                       ifelse(males3$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males3$id[which(males3$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males3$id[which(males3$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males3$p3[which(males3$degree_0.3 != 0)],
     vertex.label = males3$id[which(males3$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males3[which(males3$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males3[which(males3$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males3[which(males3$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males3[which(males3$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males3[which(males3$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 4 ################
print(paste0('Period 4 started at ', Sys.time()))

draws_anp_4 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period4_22.04.27.csv')
saveRDS(draws_anp_4, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period4_22.04.27.rds')
draws_anp_4 <- data.matrix(draws_anp_4)

### create data list
counts_df_period4 <- counts_df_non0[counts_df_non0$period == 4,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_4[,2:ncol(draws_anp_4)]),
                        min = rep(NA, ncol(draws_anp_4)-1),
                        max = rep(NA, ncol(draws_anp_4)-1),
                        mean = rep(NA, ncol(draws_anp_4)-1),
                        median = rep(NA, ncol(draws_anp_4)-1),
                        sd = rep(NA, ncol(draws_anp_4)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_4[,i+1])
  summaries$max[i]    <- max(draws_anp_4[,i+1])
  summaries$mean[i]   <- mean(draws_anp_4[,i+1])
  summaries$median[i] <- median(draws_anp_4[,i+1])
  summaries$sd[i]     <- sd(draws_anp_4[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period4$id_1))+1,
                         NROW(unique(counts_df_period4$id_2))+1,
                         NROW(draws_anp_4)),
                    dimnames = list(c(unique(counts_df_period4$id_1),
                                      unique(counts_df_period4$id_2)[length(unique(counts_df_period4$id_2))]),
                                    c(unique(counts_df_period4$id_1),
                                      unique(counts_df_period4$id_2)[length(unique(counts_df_period4$id_2))]),
                                    NULL))
N <- nrow(counts_df_period4)

for (i in 1:N) {
  dyad_row <- counts_df_period4[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_4[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids4 <- c(unique(counts_df_period4$id_1), unique(counts_df_period4$id_2)[length(unique(counts_df_period4$id_2))])
males4 <- males %>% dplyr::filter(id %in% ids4)

# create variables for different degrees of node connectedness
males4$degree_0.1 <- NA
males4$degree_0.2 <- NA
males4$degree_0.3 <- NA
males4$degree_0.4 <- NA
males4$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males4)){
  rows <- summaries[summaries$id_1 == males4$id[i] | summaries$id_2 == males4$id[i],]
  males4$degree_0.1[i] <- length(which(rows$median > 0.1))
  males4$degree_0.2[i] <- length(which(rows$median > 0.2))
  males4$degree_0.3[i] <- length(which(rows$median > 0.3))
  males4$degree_0.4[i] <- length(which(rows$median > 0.4))
  males4$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males4$age <- lubridate::year(periods$period_start[periods$period == 4]) - males4$byr
summary(males4$age)
males4$age_class <- ifelse(males4$age < 5, 1,
                           ifelse(males4$age < 10, 2,
                                  ifelse(males4$age < 15, 3,
                                         ifelse(males4$age < 20, 4,
                                                ifelse(males4$age < 25, 5,
                                                       ifelse(males4$age < 40, 6, 7))))))

g_mid_0.2 <- delete.vertices(graph = g_mid, v = males4$id[which(males4$degree_0.2 == 0)])
g_rng_0.2 <- delete.vertices(graph = g_rng, v = males4$id[which(males4$degree_0.2 == 0)])

set.seed(3)
coords_0.2 <- layout_nicely(g_mid_0.2)

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
     vertex.size = males4$p4[which(males4$degree_0.2 != 0)],
     vertex.label = males4$id[which(males4$degree_0.2 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males4[which(males4$degree_0.2 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males4[which(males4$degree_0.2 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males4[which(males4$degree_0.2 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males4[which(males4$degree_0.2 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males4[which(males4$degree_0.2 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.2, add = TRUE)

################ Period 5 ################
print(paste0('Period 5 started at ', Sys.time()))

draws_anp_5 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period5_22.04.27.csv')
saveRDS(draws_anp_5, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period5_22.04.27.rds')
draws_anp_5 <- data.matrix(draws_anp_5)

### create data list
counts_df_period5 <- counts_df_non0[counts_df_non0$period == 5,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_5[,2:ncol(draws_anp_5)]),
                        min = rep(NA, ncol(draws_anp_5)-1),
                        max = rep(NA, ncol(draws_anp_5)-1),
                        mean = rep(NA, ncol(draws_anp_5)-1),
                        median = rep(NA, ncol(draws_anp_5)-1),
                        sd = rep(NA, ncol(draws_anp_5)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_5[,i+1])
  summaries$max[i]    <- max(draws_anp_5[,i+1])
  summaries$mean[i]   <- mean(draws_anp_5[,i+1])
  summaries$median[i] <- median(draws_anp_5[,i+1])
  summaries$sd[i]     <- sd(draws_anp_5[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period5$id_1))+1,
                         NROW(unique(counts_df_period5$id_2))+1,
                         NROW(draws_anp_5)),
                    dimnames = list(c(unique(counts_df_period5$id_1),
                                      unique(counts_df_period5$id_2)[length(unique(counts_df_period5$id_2))]),
                                    c(unique(counts_df_period5$id_1),
                                      unique(counts_df_period5$id_2)[length(unique(counts_df_period5$id_2))]),
                                    NULL))
N <- nrow(counts_df_period5)

for (i in 1:N) {
  dyad_row <- counts_df_period5[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_5[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids5 <- c(unique(counts_df_period5$id_1), unique(counts_df_period5$id_2)[length(unique(counts_df_period5$id_2))])
males5 <- males %>% dplyr::filter(id %in% ids5)

# create variables for different degrees of node connectedness
males5$degree_0.1 <- NA
males5$degree_0.2 <- NA
males5$degree_0.3 <- NA
males5$degree_0.4 <- NA
males5$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males5)){
  rows <- summaries[summaries$id_1 == males5$id[i] | summaries$id_2 == males5$id[i],]
  males5$degree_0.1[i] <- length(which(rows$median > 0.1))
  males5$degree_0.2[i] <- length(which(rows$median > 0.2))
  males5$degree_0.3[i] <- length(which(rows$median > 0.3))
  males5$degree_0.4[i] <- length(which(rows$median > 0.4))
  males5$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males5$age <- lubridate::year(periods$period_start[periods$period == 5]) - males5$byr
summary(males5$age)
males5$age_class <- ifelse(males5$age < 5, 1,
                           ifelse(males5$age < 10, 2,
                                  ifelse(males5$age < 15, 3,
                                         ifelse(males5$age < 20, 4,
                                                ifelse(males5$age < 25, 5,
                                                       ifelse(males5$age < 40, 6, 7))))))

g_mid_0.2 <- delete.vertices(graph = g_mid, v = males5$id[which(males5$degree_0.2 == 0)])
g_rng_0.2 <- delete.vertices(graph = g_rng, v = males5$id[which(males5$degree_0.2 == 0)])

set.seed(3)
coords_0.2 <- layout_nicely(g_mid_0.2)

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
     vertex.size = males5$p5[which(males5$degree_0.2 != 0)],
     vertex.label = males5$id[which(males5$degree_0.2 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males5[which(males5$degree_0.2 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males5[which(males5$degree_0.2 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males5[which(males5$degree_0.2 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males5[which(males5$degree_0.2 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males5[which(males5$degree_0.2 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.2, add = TRUE)

################ Period 6 ################
print(paste0('Period 6 started at ', Sys.time()))

draws_anp_6 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period6_22.04.27.csv')
saveRDS(draws_anp_6, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period6_22.04.27.rds')
draws_anp_6 <- data.matrix(draws_anp_6)

### create data list
counts_df_period6 <- counts_df_non0[counts_df_non0$period == 6,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_6[,2:ncol(draws_anp_6)]),
                        min = rep(NA, ncol(draws_anp_6)-1),
                        max = rep(NA, ncol(draws_anp_6)-1),
                        mean = rep(NA, ncol(draws_anp_6)-1),
                        median = rep(NA, ncol(draws_anp_6)-1),
                        sd = rep(NA, ncol(draws_anp_6)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_6[,i+1])
  summaries$max[i]    <- max(draws_anp_6[,i+1])
  summaries$mean[i]   <- mean(draws_anp_6[,i+1])
  summaries$median[i] <- median(draws_anp_6[,i+1])
  summaries$sd[i]     <- sd(draws_anp_6[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period6$id_1))+1,
                         NROW(unique(counts_df_period6$id_2))+1,
                         NROW(draws_anp_6)),
                    dimnames = list(c(unique(counts_df_period6$id_1),
                                      unique(counts_df_period6$id_2)[length(unique(counts_df_period6$id_2))]),
                                    c(unique(counts_df_period6$id_1),
                                      unique(counts_df_period6$id_2)[length(unique(counts_df_period6$id_2))]),
                                    NULL))
N <- nrow(counts_df_period6)

for (i in 1:N) {
  dyad_row <- counts_df_period6[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_6[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids6 <- c(unique(counts_df_period6$id_1), unique(counts_df_period6$id_2)[length(unique(counts_df_period6$id_2))])
males6 <- males %>% dplyr::filter(id %in% ids6)

# create variables for different degrees of node connectedness
males6$degree_0.1 <- NA
males6$degree_0.2 <- NA
males6$degree_0.3 <- NA
males6$degree_0.4 <- NA
males6$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males6)){
  rows <- summaries[summaries$id_1 == males6$id[i] | summaries$id_2 == males6$id[i],]
  males6$degree_0.1[i] <- length(which(rows$median > 0.1))
  males6$degree_0.2[i] <- length(which(rows$median > 0.2))
  males6$degree_0.3[i] <- length(which(rows$median > 0.3))
  males6$degree_0.4[i] <- length(which(rows$median > 0.4))
  males6$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males6$age <- lubridate::year(periods$period_start[periods$period == 6]) - males6$byr
summary(males6$age)
males6$age_class <- ifelse(males6$age < 5, 1,
                           ifelse(males6$age < 10, 2,
                                  ifelse(males6$age < 15, 3,
                                         ifelse(males6$age < 20, 4,
                                                ifelse(males6$age < 25, 5,
                                                       ifelse(males6$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males6$id[which(males6$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males6$id[which(males6$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males6$p6[which(males6$degree_0.3 != 0)],
     vertex.label = males6$id[which(males6$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males6[which(males6$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males6[which(males6$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males6[which(males6$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males6[which(males6$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males6[which(males6$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 7 ################
print(paste0('Period 7 started at ', Sys.time()))

draws_anp_7 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period7_22.04.27.csv')
saveRDS(draws_anp_7, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period7_22.04.27.rds')
draws_anp_7 <- data.matrix(draws_anp_7)

### create data list
counts_df_period7 <- counts_df_non0[counts_df_non0$period == 7,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_7[,2:ncol(draws_anp_7)]),
                        min = rep(NA, ncol(draws_anp_7)-1),
                        max = rep(NA, ncol(draws_anp_7)-1),
                        mean = rep(NA, ncol(draws_anp_7)-1),
                        median = rep(NA, ncol(draws_anp_7)-1),
                        sd = rep(NA, ncol(draws_anp_7)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_7[,i+1])
  summaries$max[i]    <- max(draws_anp_7[,i+1])
  summaries$mean[i]   <- mean(draws_anp_7[,i+1])
  summaries$median[i] <- median(draws_anp_7[,i+1])
  summaries$sd[i]     <- sd(draws_anp_7[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period7$id_1))+1,
                         NROW(unique(counts_df_period7$id_2))+1,
                         NROW(draws_anp_7)),
                    dimnames = list(c(unique(counts_df_period7$id_1),
                                      unique(counts_df_period7$id_2)[length(unique(counts_df_period7$id_2))]),
                                    c(unique(counts_df_period7$id_1),
                                      unique(counts_df_period7$id_2)[length(unique(counts_df_period7$id_2))]),
                                    NULL))
N <- nrow(counts_df_period7)

for (i in 1:N) {
  dyad_row <- counts_df_period7[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_7[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids7 <- c(unique(counts_df_period7$id_1), unique(counts_df_period7$id_2)[length(unique(counts_df_period7$id_2))])
males7 <- males %>% dplyr::filter(id %in% ids7)

# create variables for different degrees of node connectedness
males7$degree_0.1 <- NA
males7$degree_0.2 <- NA
males7$degree_0.3 <- NA
males7$degree_0.4 <- NA
males7$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males7)){
  rows <- summaries[summaries$id_1 == males7$id[i] | summaries$id_2 == males7$id[i],]
  males7$degree_0.1[i] <- length(which(rows$median > 0.1))
  males7$degree_0.2[i] <- length(which(rows$median > 0.2))
  males7$degree_0.3[i] <- length(which(rows$median > 0.3))
  males7$degree_0.4[i] <- length(which(rows$median > 0.4))
  males7$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males7$age <- lubridate::year(periods$period_start[periods$period == 7]) - males7$byr
summary(males7$age)
males7$age_class <- ifelse(males7$age < 5, 1,
                           ifelse(males7$age < 10, 2,
                                  ifelse(males7$age < 15, 3,
                                         ifelse(males7$age < 20, 4,
                                                ifelse(males7$age < 25, 5,
                                                       ifelse(males7$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males7$id[which(males7$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males7$id[which(males7$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males7$p7[which(males7$degree_0.3 != 0)],
     vertex.label = males7$id[which(males7$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males7[which(males7$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males7[which(males7$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males7[which(males7$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males7[which(males7$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males7[which(males7$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 8 ################
print(paste0('Period 8 started at ', Sys.time()))

draws_anp_8 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period8_22.04.27.csv')
saveRDS(draws_anp_8, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period8_22.04.27.rds')
draws_anp_8 <- data.matrix(draws_anp_8)

### create data list
counts_df_period8 <- counts_df_non0[counts_df_non0$period == 8,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_8[,2:ncol(draws_anp_8)]),
                        min = rep(NA, ncol(draws_anp_8)-1),
                        max = rep(NA, ncol(draws_anp_8)-1),
                        mean = rep(NA, ncol(draws_anp_8)-1),
                        median = rep(NA, ncol(draws_anp_8)-1),
                        sd = rep(NA, ncol(draws_anp_8)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_8[,i+1])
  summaries$max[i]    <- max(draws_anp_8[,i+1])
  summaries$mean[i]   <- mean(draws_anp_8[,i+1])
  summaries$median[i] <- median(draws_anp_8[,i+1])
  summaries$sd[i]     <- sd(draws_anp_8[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period8$id_1))+1,
                         NROW(unique(counts_df_period8$id_2))+1,
                         NROW(draws_anp_8)),
                    dimnames = list(c(unique(counts_df_period8$id_1),
                                      unique(counts_df_period8$id_2)[length(unique(counts_df_period8$id_2))]),
                                    c(unique(counts_df_period8$id_1),
                                      unique(counts_df_period8$id_2)[length(unique(counts_df_period8$id_2))]),
                                    NULL))
N <- nrow(counts_df_period8)

for (i in 1:N) {
  dyad_row <- counts_df_period8[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_8[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids8 <- c(unique(counts_df_period8$id_1), unique(counts_df_period8$id_2)[length(unique(counts_df_period8$id_2))])
males8 <- males %>% dplyr::filter(id %in% ids8)

# create variables for different degrees of node connectedness
males8$degree_0.1 <- NA
males8$degree_0.2 <- NA
males8$degree_0.3 <- NA
males8$degree_0.4 <- NA
males8$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males8)){
  rows <- summaries[summaries$id_1 == males8$id[i] | summaries$id_2 == males8$id[i],]
  males8$degree_0.1[i] <- length(which(rows$median > 0.1))
  males8$degree_0.2[i] <- length(which(rows$median > 0.2))
  males8$degree_0.3[i] <- length(which(rows$median > 0.3))
  males8$degree_0.4[i] <- length(which(rows$median > 0.4))
  males8$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males8$age <- lubridate::year(periods$period_start[periods$period == 8]) - males8$byr
summary(males8$age)
males8$age_class <- ifelse(males8$age < 5, 1,
                           ifelse(males8$age < 10, 2,
                                  ifelse(males8$age < 15, 3,
                                         ifelse(males8$age < 20, 4,
                                                ifelse(males8$age < 25, 5,
                                                       ifelse(males8$age < 40, 6, 7))))))

g_mid_0.2 <- delete.vertices(graph = g_mid, v = males8$id[which(males8$degree_0.2 == 0)])
g_rng_0.2 <- delete.vertices(graph = g_rng, v = males8$id[which(males8$degree_0.2 == 0)])

set.seed(3)
coords_0.2 <- layout_nicely(g_mid_0.2)

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
     vertex.size = males8$p8[which(males8$degree_0.2 != 0)],
     vertex.label = males8$id[which(males8$degree_0.2 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males8[which(males8$degree_0.2 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males8[which(males8$degree_0.2 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males8[which(males8$degree_0.2 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males8[which(males8$degree_0.2 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males8[which(males8$degree_0.2 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.2, add = TRUE)

################ Period 9 ################
print(paste0('Period 9 started at ', Sys.time()))

draws_anp_9 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period9_22.04.27.csv')
saveRDS(draws_anp_9, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period9_22.04.27.rds')
draws_anp_9 <- data.matrix(draws_anp_9)

### create data list
counts_df_period9 <- counts_df_non0[counts_df_non0$period == 9,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_9[,2:ncol(draws_anp_9)]),
                        min = rep(NA, ncol(draws_anp_9)-1),
                        max = rep(NA, ncol(draws_anp_9)-1),
                        mean = rep(NA, ncol(draws_anp_9)-1),
                        median = rep(NA, ncol(draws_anp_9)-1),
                        sd = rep(NA, ncol(draws_anp_9)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_9[,i+1])
  summaries$max[i]    <- max(draws_anp_9[,i+1])
  summaries$mean[i]   <- mean(draws_anp_9[,i+1])
  summaries$median[i] <- median(draws_anp_9[,i+1])
  summaries$sd[i]     <- sd(draws_anp_9[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period9$id_1))+1,
                         NROW(unique(counts_df_period9$id_2))+1,
                         NROW(draws_anp_9)),
                    dimnames = list(c(unique(counts_df_period9$id_1),
                                      unique(counts_df_period9$id_2)[length(unique(counts_df_period9$id_2))]),
                                    c(unique(counts_df_period9$id_1),
                                      unique(counts_df_period9$id_2)[length(unique(counts_df_period9$id_2))]),
                                    NULL))
N <- nrow(counts_df_period9)

for (i in 1:N) {
  dyad_row <- counts_df_period9[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_9[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids9 <- c(unique(counts_df_period9$id_1), unique(counts_df_period9$id_2)[length(unique(counts_df_period9$id_2))])
males9 <- males %>% dplyr::filter(id %in% ids9)

# create variables for different degrees of node connectedness
males9$degree_0.1 <- NA
males9$degree_0.2 <- NA
males9$degree_0.3 <- NA
males9$degree_0.4 <- NA
males9$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males9)){
  rows <- summaries[summaries$id_1 == males9$id[i] | summaries$id_2 == males9$id[i],]
  males9$degree_0.1[i] <- length(which(rows$median > 0.1))
  males9$degree_0.2[i] <- length(which(rows$median > 0.2))
  males9$degree_0.3[i] <- length(which(rows$median > 0.3))
  males9$degree_0.4[i] <- length(which(rows$median > 0.4))
  males9$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males9$age <- lubridate::year(periods$period_start[periods$period == 9]) - males9$byr
summary(males9$age)
males9$age_class <- ifelse(males9$age < 5, 1,
                           ifelse(males9$age < 10, 2,
                                  ifelse(males9$age < 15, 3,
                                         ifelse(males9$age < 20, 4,
                                                ifelse(males9$age < 25, 5,
                                                       ifelse(males9$age < 40, 6, 7))))))

g_mid_0.2 <- delete.vertices(graph = g_mid, v = males9$id[which(males9$degree_0.2 == 0)])
g_rng_0.2 <- delete.vertices(graph = g_rng, v = males9$id[which(males9$degree_0.2 == 0)])

set.seed(3)
coords_0.2 <- layout_nicely(g_mid_0.2)

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
     vertex.size = males9$p9[which(males9$degree_0.2 != 0)],
     vertex.label = males9$id[which(males9$degree_0.2 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males9[which(males9$degree_0.2 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males9[which(males9$degree_0.2 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males9[which(males9$degree_0.2 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males9[which(males9$degree_0.2 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males9[which(males9$degree_0.2 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.2, add = TRUE)

################ Period 10 ################
print(paste0('Period 10 started at ', Sys.time()))

draws_anp_10 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period10_22.04.27.csv')
saveRDS(draws_anp_10, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period10_22.04.27.rds')
draws_anp_10 <- data.matrix(draws_anp_10)

### create data list
counts_df_period10 <- counts_df_non0[counts_df_non0$period == 10,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_10[,2:ncol(draws_anp_10)]),
                        min = rep(NA, ncol(draws_anp_10)-1),
                        max = rep(NA, ncol(draws_anp_10)-1),
                        mean = rep(NA, ncol(draws_anp_10)-1),
                        median = rep(NA, ncol(draws_anp_10)-1),
                        sd = rep(NA, ncol(draws_anp_10)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_10[,i+1])
  summaries$max[i]    <- max(draws_anp_10[,i+1])
  summaries$mean[i]   <- mean(draws_anp_10[,i+1])
  summaries$median[i] <- median(draws_anp_10[,i+1])
  summaries$sd[i]     <- sd(draws_anp_10[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period10$id_1))+1,
                         NROW(unique(counts_df_period10$id_2))+1,
                         NROW(draws_anp_10)),
                    dimnames = list(c(unique(counts_df_period10$id_1),
                                      unique(counts_df_period10$id_2)[length(unique(counts_df_period10$id_2))]),
                                    c(unique(counts_df_period10$id_1),
                                      unique(counts_df_period10$id_2)[length(unique(counts_df_period10$id_2))]),
                                    NULL))
N <- nrow(counts_df_period10)

for (i in 1:N) {
  dyad_row <- counts_df_period10[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_10[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids10 <- c(unique(counts_df_period10$id_1), unique(counts_df_period10$id_2)[length(unique(counts_df_period10$id_2))])
males10 <- males %>% dplyr::filter(id %in% ids10)

# create variables for different degrees of node connectedness
males10$degree_0.1 <- NA
males10$degree_0.2 <- NA
males10$degree_0.3 <- NA
males10$degree_0.4 <- NA
males10$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males10)){
  rows <- summaries[summaries$id_1 == males10$id[i] | summaries$id_2 == males10$id[i],]
  males10$degree_0.1[i] <- length(which(rows$median > 0.1))
  males10$degree_0.2[i] <- length(which(rows$median > 0.2))
  males10$degree_0.3[i] <- length(which(rows$median > 0.3))
  males10$degree_0.4[i] <- length(which(rows$median > 0.4))
  males10$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males10$age <- lubridate::year(periods$period_start[periods$period == 10]) - males10$byr
summary(males10$age)
males10$age_class <- ifelse(males10$age < 5, 1,
                            ifelse(males10$age < 10, 2,
                                   ifelse(males10$age < 15, 3,
                                          ifelse(males10$age < 20, 4,
                                                 ifelse(males10$age < 25, 5,
                                                        ifelse(males10$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males10$id[which(males10$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males10$id[which(males10$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males10$p10[which(males10$degree_0.3 != 0)],
     vertex.label = males10$id[which(males10$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males10[which(males10$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males10[which(males10$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males10[which(males10$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males10[which(males10$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males10[which(males10$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 11 ################
print(paste0('Period 11 started at ', Sys.time()))

draws_anp_11 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period11_22.04.27.csv')
saveRDS(draws_anp_11, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period11_22.04.27.rds')
draws_anp_11 <- data.matrix(draws_anp_11)

### create data list
counts_df_period11 <- counts_df_non0[counts_df_non0$period == 11,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_11[,2:ncol(draws_anp_11)]),
                        min = rep(NA, ncol(draws_anp_11)-1),
                        max = rep(NA, ncol(draws_anp_11)-1),
                        mean = rep(NA, ncol(draws_anp_11)-1),
                        median = rep(NA, ncol(draws_anp_11)-1),
                        sd = rep(NA, ncol(draws_anp_11)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_11[,i+1])
  summaries$max[i]    <- max(draws_anp_11[,i+1])
  summaries$mean[i]   <- mean(draws_anp_11[,i+1])
  summaries$median[i] <- median(draws_anp_11[,i+1])
  summaries$sd[i]     <- sd(draws_anp_11[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period11$id_1))+1,
                         NROW(unique(counts_df_period11$id_2))+1,
                         NROW(draws_anp_11)),
                    dimnames = list(c(unique(counts_df_period11$id_1),
                                      unique(counts_df_period11$id_2)[length(unique(counts_df_period11$id_2))]),
                                    c(unique(counts_df_period11$id_1),
                                      unique(counts_df_period11$id_2)[length(unique(counts_df_period11$id_2))]),
                                    NULL))
N <- nrow(counts_df_period11)

for (i in 1:N) {
  dyad_row <- counts_df_period11[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_11[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids11 <- c(unique(counts_df_period11$id_1), unique(counts_df_period11$id_2)[length(unique(counts_df_period11$id_2))])
males11 <- males %>% dplyr::filter(id %in% ids11)

# create variables for different degrees of node connectedness
males11$degree_0.1 <- NA
males11$degree_0.2 <- NA
males11$degree_0.3 <- NA
males11$degree_0.4 <- NA
males11$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males11)){
  rows <- summaries[summaries$id_1 == males11$id[i] | summaries$id_2 == males11$id[i],]
  males11$degree_0.1[i] <- length(which(rows$median > 0.1))
  males11$degree_0.2[i] <- length(which(rows$median > 0.2))
  males11$degree_0.3[i] <- length(which(rows$median > 0.3))
  males11$degree_0.4[i] <- length(which(rows$median > 0.4))
  males11$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males11$age <- lubridate::year(periods$period_start[periods$period == 11]) - males11$byr
summary(males11$age)
males11$age_class <- ifelse(males11$age < 5, 1,
                            ifelse(males11$age < 10, 2,
                                   ifelse(males11$age < 15, 3,
                                          ifelse(males11$age < 20, 4,
                                                 ifelse(males11$age < 25, 5,
                                                        ifelse(males11$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males11$id[which(males11$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males11$id[which(males11$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males11$p11[which(males11$degree_0.3 != 0)],
     vertex.label = males11$id[which(males11$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males11[which(males11$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males11[which(males11$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males11[which(males11$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males11[which(males11$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males11[which(males11$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 12 ################
print(paste0('Period 12 started at ', Sys.time()))

draws_anp_12 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period12_22.04.27.csv')
saveRDS(draws_anp_12, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period12_22.04.27.rds')
draws_anp_12 <- data.matrix(draws_anp_12)

### create data list
counts_df_period12 <- counts_df_non0[counts_df_non0$period == 12,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_12[,2:ncol(draws_anp_12)]),
                        min = rep(NA, ncol(draws_anp_12)-1),
                        max = rep(NA, ncol(draws_anp_12)-1),
                        mean = rep(NA, ncol(draws_anp_12)-1),
                        median = rep(NA, ncol(draws_anp_12)-1),
                        sd = rep(NA, ncol(draws_anp_12)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_12[,i+1])
  summaries$max[i]    <- max(draws_anp_12[,i+1])
  summaries$mean[i]   <- mean(draws_anp_12[,i+1])
  summaries$median[i] <- median(draws_anp_12[,i+1])
  summaries$sd[i]     <- sd(draws_anp_12[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period12$id_1))+1,
                         NROW(unique(counts_df_period12$id_2))+1,
                         NROW(draws_anp_12)),
                    dimnames = list(c(unique(counts_df_period12$id_1),
                                      unique(counts_df_period12$id_2)[length(unique(counts_df_period12$id_2))]),
                                    c(unique(counts_df_period12$id_1),
                                      unique(counts_df_period12$id_2)[length(unique(counts_df_period12$id_2))]),
                                    NULL))
N <- nrow(counts_df_period12)

for (i in 1:N) {
  dyad_row <- counts_df_period12[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_12[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids12 <- c(unique(counts_df_period12$id_1), unique(counts_df_period12$id_2)[length(unique(counts_df_period12$id_2))])
males12 <- males %>% dplyr::filter(id %in% ids12)

# create variables for different degrees of node connectedness
males12$degree_0.1 <- NA
males12$degree_0.2 <- NA
males12$degree_0.3 <- NA
males12$degree_0.4 <- NA
males12$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males12)){
  rows <- summaries[summaries$id_1 == males12$id[i] | summaries$id_2 == males12$id[i],]
  males12$degree_0.1[i] <- length(which(rows$median > 0.1))
  males12$degree_0.2[i] <- length(which(rows$median > 0.2))
  males12$degree_0.3[i] <- length(which(rows$median > 0.3))
  males12$degree_0.4[i] <- length(which(rows$median > 0.4))
  males12$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males12$age <- lubridate::year(periods$period_start[periods$period == 12]) - males12$byr
summary(males12$age)
males12$age_class <- ifelse(males12$age < 5, 1,
                            ifelse(males12$age < 10, 2,
                                   ifelse(males12$age < 15, 3,
                                          ifelse(males12$age < 20, 4,
                                                 ifelse(males12$age < 25, 5,
                                                        ifelse(males12$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males12$id[which(males12$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males12$id[which(males12$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males12$p12[which(males12$degree_0.3 != 0)],
     vertex.label = males12$id[which(males12$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males12[which(males12$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males12[which(males12$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males12[which(males12$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males12[which(males12$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males12[which(males12$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 13 ################
print(paste0('Period 13 started at ', Sys.time()))

draws_anp_13 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period13_22.04.27.csv')
saveRDS(draws_anp_13, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period13_22.04.27.rds')
draws_anp_13 <- data.matrix(draws_anp_13)

### create data list
counts_df_period13 <- counts_df_non0[counts_df_non0$period == 13,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_13[,2:ncol(draws_anp_13)]),
                        min = rep(NA, ncol(draws_anp_13)-1),
                        max = rep(NA, ncol(draws_anp_13)-1),
                        mean = rep(NA, ncol(draws_anp_13)-1),
                        median = rep(NA, ncol(draws_anp_13)-1),
                        sd = rep(NA, ncol(draws_anp_13)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_13[,i+1])
  summaries$max[i]    <- max(draws_anp_13[,i+1])
  summaries$mean[i]   <- mean(draws_anp_13[,i+1])
  summaries$median[i] <- median(draws_anp_13[,i+1])
  summaries$sd[i]     <- sd(draws_anp_13[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period13$id_1))+1,
                         NROW(unique(counts_df_period13$id_2))+1,
                         NROW(draws_anp_13)),
                    dimnames = list(c(unique(counts_df_period13$id_1),
                                      unique(counts_df_period13$id_2)[length(unique(counts_df_period13$id_2))]),
                                    c(unique(counts_df_period13$id_1),
                                      unique(counts_df_period13$id_2)[length(unique(counts_df_period13$id_2))]),
                                    NULL))
N <- nrow(counts_df_period13)

for (i in 1:N) {
  dyad_row <- counts_df_period13[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_13[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids13 <- c(unique(counts_df_period13$id_1), unique(counts_df_period13$id_2)[length(unique(counts_df_period13$id_2))])
males13 <- males %>% dplyr::filter(id %in% ids13)

# create variables for different degrees of node connectedness
males13$degree_0.1 <- NA
males13$degree_0.2 <- NA
males13$degree_0.3 <- NA
males13$degree_0.4 <- NA
males13$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males13)){
  rows <- summaries[summaries$id_1 == males13$id[i] | summaries$id_2 == males13$id[i],]
  males13$degree_0.1[i] <- length(which(rows$median > 0.1))
  males13$degree_0.2[i] <- length(which(rows$median > 0.2))
  males13$degree_0.3[i] <- length(which(rows$median > 0.3))
  males13$degree_0.4[i] <- length(which(rows$median > 0.4))
  males13$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males13$age <- lubridate::year(periods$period_start[periods$period == 13]) - males13$byr
summary(males13$age)
males13$age_class <- ifelse(males13$age < 5, 1,
                            ifelse(males13$age < 10, 2,
                                   ifelse(males13$age < 15, 3,
                                          ifelse(males13$age < 20, 4,
                                                 ifelse(males13$age < 25, 5,
                                                        ifelse(males13$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males13$id[which(males13$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males13$id[which(males13$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males13$p13[which(males13$degree_0.3 != 0)],
     vertex.label = males13$id[which(males13$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males13[which(males13$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males13[which(males13$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males13[which(males13$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males13[which(males13$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males13[which(males13$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 14 ################
print(paste0('Period 14 started at ', Sys.time()))

draws_anp_14 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period14_22.04.27.csv')
saveRDS(draws_anp_14, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period14_22.04.27.rds')
draws_anp_14 <- data.matrix(draws_anp_14)

### create data list
counts_df_period14 <- counts_df_non0[counts_df_non0$period == 14,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_14[,2:ncol(draws_anp_14)]),
                        min = rep(NA, ncol(draws_anp_14)-1),
                        max = rep(NA, ncol(draws_anp_14)-1),
                        mean = rep(NA, ncol(draws_anp_14)-1),
                        median = rep(NA, ncol(draws_anp_14)-1),
                        sd = rep(NA, ncol(draws_anp_14)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_14[,i+1])
  summaries$max[i]    <- max(draws_anp_14[,i+1])
  summaries$mean[i]   <- mean(draws_anp_14[,i+1])
  summaries$median[i] <- median(draws_anp_14[,i+1])
  summaries$sd[i]     <- sd(draws_anp_14[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period14$id_1))+1,
                         NROW(unique(counts_df_period14$id_2))+1,
                         NROW(draws_anp_14)),
                    dimnames = list(c(unique(counts_df_period14$id_1),
                                      unique(counts_df_period14$id_2)[length(unique(counts_df_period14$id_2))]),
                                    c(unique(counts_df_period14$id_1),
                                      unique(counts_df_period14$id_2)[length(unique(counts_df_period14$id_2))]),
                                    NULL))
N <- nrow(counts_df_period14)

for (i in 1:N) {
  dyad_row <- counts_df_period14[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_14[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids14 <- c(unique(counts_df_period14$id_1), unique(counts_df_period14$id_2)[length(unique(counts_df_period14$id_2))])
males14 <- males %>% dplyr::filter(id %in% ids14)

# create variables for different degrees of node connectedness
males14$degree_0.1 <- NA
males14$degree_0.2 <- NA
males14$degree_0.3 <- NA
males14$degree_0.4 <- NA
males14$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males14)){
  rows <- summaries[summaries$id_1 == males14$id[i] | summaries$id_2 == males14$id[i],]
  males14$degree_0.1[i] <- length(which(rows$median > 0.1))
  males14$degree_0.2[i] <- length(which(rows$median > 0.2))
  males14$degree_0.3[i] <- length(which(rows$median > 0.3))
  males14$degree_0.4[i] <- length(which(rows$median > 0.4))
  males14$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males14$age <- lubridate::year(periods$period_start[periods$period == 14]) - males14$byr
summary(males14$age)
males14$age_class <- ifelse(males14$age < 5, 1,
                            ifelse(males14$age < 10, 2,
                                   ifelse(males14$age < 15, 3,
                                          ifelse(males14$age < 20, 4,
                                                 ifelse(males14$age < 25, 5,
                                                        ifelse(males14$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males14$id[which(males14$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males14$id[which(males14$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males14$p14[which(males14$degree_0.3 != 0)],
     vertex.label = males14$id[which(males14$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males14[which(males14$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males14[which(males14$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males14[which(males14$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males14[which(males14$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males14[which(males14$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 15 ################
print(paste0('Period 15 started at ', Sys.time()))

draws_anp_15 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period15_22.04.27.csv')
saveRDS(draws_anp_15, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period15_22.04.27.rds')
draws_anp_15 <- data.matrix(draws_anp_15)

### create data list
counts_df_period15 <- counts_df_non0[counts_df_non0$period == 15,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_15[,2:ncol(draws_anp_15)]),
                        min = rep(NA, ncol(draws_anp_15)-1),
                        max = rep(NA, ncol(draws_anp_15)-1),
                        mean = rep(NA, ncol(draws_anp_15)-1),
                        median = rep(NA, ncol(draws_anp_15)-1),
                        sd = rep(NA, ncol(draws_anp_15)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_15[,i+1])
  summaries$max[i]    <- max(draws_anp_15[,i+1])
  summaries$mean[i]   <- mean(draws_anp_15[,i+1])
  summaries$median[i] <- median(draws_anp_15[,i+1])
  summaries$sd[i]     <- sd(draws_anp_15[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period15$id_1))+1,
                         NROW(unique(counts_df_period15$id_2))+1,
                         NROW(draws_anp_15)),
                    dimnames = list(c(unique(counts_df_period15$id_1),
                                      unique(counts_df_period15$id_2)[length(unique(counts_df_period15$id_2))]),
                                    c(unique(counts_df_period15$id_1),
                                      unique(counts_df_period15$id_2)[length(unique(counts_df_period15$id_2))]),
                                    NULL))
N <- nrow(counts_df_period15)

for (i in 1:N) {
  dyad_row <- counts_df_period15[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_15[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids15 <- c(unique(counts_df_period15$id_1), unique(counts_df_period15$id_2)[length(unique(counts_df_period15$id_2))])
males15 <- males %>% dplyr::filter(id %in% ids15)

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

# age variable
males15$age <- lubridate::year(periods$period_start[periods$period == 15]) - males15$byr
summary(males15$age)
males15$age_class <- ifelse(males15$age < 5, 1,
                            ifelse(males15$age < 10, 2,
                                   ifelse(males15$age < 15, 3,
                                          ifelse(males15$age < 20, 4,
                                                 ifelse(males15$age < 25, 5,
                                                        ifelse(males15$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males15$id[which(males15$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males15$id[which(males15$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males15$p15[which(males15$degree_0.3 != 0)],
     vertex.label = males15$id[which(males15$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males15[which(males15$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males15[which(males15$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males15[which(males15$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males15[which(males15$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males15[which(males15$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 16 ################
print(paste0('Period 16 started at ', Sys.time()))

draws_anp_16 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period16_22.04.27.csv')
saveRDS(draws_anp_16, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period16_22.04.27.rds')
draws_anp_16 <- data.matrix(draws_anp_16)

### create data list
counts_df_period16 <- counts_df_non0[counts_df_non0$period == 16,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_16[,2:ncol(draws_anp_16)]),
                        min = rep(NA, ncol(draws_anp_16)-1),
                        max = rep(NA, ncol(draws_anp_16)-1),
                        mean = rep(NA, ncol(draws_anp_16)-1),
                        median = rep(NA, ncol(draws_anp_16)-1),
                        sd = rep(NA, ncol(draws_anp_16)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_16[,i+1])
  summaries$max[i]    <- max(draws_anp_16[,i+1])
  summaries$mean[i]   <- mean(draws_anp_16[,i+1])
  summaries$median[i] <- median(draws_anp_16[,i+1])
  summaries$sd[i]     <- sd(draws_anp_16[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period16$id_1))+1,
                         NROW(unique(counts_df_period16$id_2))+1,
                         NROW(draws_anp_16)),
                    dimnames = list(c(unique(counts_df_period16$id_1),
                                      unique(counts_df_period16$id_2)[length(unique(counts_df_period16$id_2))]),
                                    c(unique(counts_df_period16$id_1),
                                      unique(counts_df_period16$id_2)[length(unique(counts_df_period16$id_2))]),
                                    NULL))
N <- nrow(counts_df_period16)

for (i in 1:N) {
  dyad_row <- counts_df_period16[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_16[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids16 <- c(unique(counts_df_period16$id_1), unique(counts_df_period16$id_2)[length(unique(counts_df_period16$id_2))])
males16 <- males %>% dplyr::filter(id %in% ids16)

# create variables for different degrees of node connectedness
males16$degree_0.1 <- NA
males16$degree_0.2 <- NA
males16$degree_0.3 <- NA
males16$degree_0.4 <- NA
males16$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males16)){
  rows <- summaries[summaries$id_1 == males16$id[i] | summaries$id_2 == males16$id[i],]
  males16$degree_0.1[i] <- length(which(rows$median > 0.1))
  males16$degree_0.2[i] <- length(which(rows$median > 0.2))
  males16$degree_0.3[i] <- length(which(rows$median > 0.3))
  males16$degree_0.4[i] <- length(which(rows$median > 0.4))
  males16$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males16$age <- lubridate::year(periods$period_start[periods$period == 16]) - males16$byr
summary(males16$age)
males16$age_class <- ifelse(males16$age < 5, 1,
                            ifelse(males16$age < 10, 2,
                                   ifelse(males16$age < 15, 3,
                                          ifelse(males16$age < 20, 4,
                                                 ifelse(males16$age < 25, 5,
                                                        ifelse(males16$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males16$id[which(males16$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males16$id[which(males16$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males16$p16[which(males16$degree_0.3 != 0)],
     vertex.label = males16$id[which(males16$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males16[which(males16$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males16[which(males16$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males16[which(males16$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males16[which(males16$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males16[which(males16$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 17 ################
print(paste0('Period 17 started at ', Sys.time()))

draws_anp_17 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period17_22.04.27.csv')
saveRDS(draws_anp_17, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period17_22.04.27.rds')
draws_anp_17 <- data.matrix(draws_anp_17)

### create data list
counts_df_period17 <- counts_df_non0[counts_df_non0$period == 17,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_17[,2:ncol(draws_anp_17)]),
                        min = rep(NA, ncol(draws_anp_17)-1),
                        max = rep(NA, ncol(draws_anp_17)-1),
                        mean = rep(NA, ncol(draws_anp_17)-1),
                        median = rep(NA, ncol(draws_anp_17)-1),
                        sd = rep(NA, ncol(draws_anp_17)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_17[,i+1])
  summaries$max[i]    <- max(draws_anp_17[,i+1])
  summaries$mean[i]   <- mean(draws_anp_17[,i+1])
  summaries$median[i] <- median(draws_anp_17[,i+1])
  summaries$sd[i]     <- sd(draws_anp_17[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period17$id_1))+1,
                         NROW(unique(counts_df_period17$id_2))+1,
                         NROW(draws_anp_17)),
                    dimnames = list(c(unique(counts_df_period17$id_1),
                                      unique(counts_df_period17$id_2)[length(unique(counts_df_period17$id_2))]),
                                    c(unique(counts_df_period17$id_1),
                                      unique(counts_df_period17$id_2)[length(unique(counts_df_period17$id_2))]),
                                    NULL))
N <- nrow(counts_df_period17)

for (i in 1:N) {
  dyad_row <- counts_df_period17[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_17[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids17 <- c(unique(counts_df_period17$id_1), unique(counts_df_period17$id_2)[length(unique(counts_df_period17$id_2))])
males17 <- males %>% dplyr::filter(id %in% ids17)

# create variables for different degrees of node connectedness
males17$degree_0.1 <- NA
males17$degree_0.2 <- NA
males17$degree_0.3 <- NA
males17$degree_0.4 <- NA
males17$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males17)){
  rows <- summaries[summaries$id_1 == males17$id[i] | summaries$id_2 == males17$id[i],]
  males17$degree_0.1[i] <- length(which(rows$median > 0.1))
  males17$degree_0.2[i] <- length(which(rows$median > 0.2))
  males17$degree_0.3[i] <- length(which(rows$median > 0.3))
  males17$degree_0.4[i] <- length(which(rows$median > 0.4))
  males17$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males17$age <- lubridate::year(periods$period_start[periods$period == 17]) - males17$byr
summary(males17$age)
males17$age_class <- ifelse(males17$age < 5, 1,
                            ifelse(males17$age < 10, 2,
                                   ifelse(males17$age < 15, 3,
                                          ifelse(males17$age < 20, 4,
                                                 ifelse(males17$age < 25, 5,
                                                        ifelse(males17$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males17$id[which(males17$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males17$id[which(males17$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males17$p17[which(males17$degree_0.3 != 0)],
     vertex.label = males17$id[which(males17$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males17[which(males17$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males17[which(males17$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males17[which(males17$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males17[which(males17$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males17[which(males17$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 18 ################
print(paste0('Period 18 started at ', Sys.time()))

draws_anp_18 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period18_22.04.27.csv')
saveRDS(draws_anp_18, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period18_22.04.27.rds')
draws_anp_18 <- data.matrix(draws_anp_18)

### create data list
counts_df_period18 <- counts_df_non0[counts_df_non0$period == 18,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_18[,2:ncol(draws_anp_18)]),
                        min = rep(NA, ncol(draws_anp_18)-1),
                        max = rep(NA, ncol(draws_anp_18)-1),
                        mean = rep(NA, ncol(draws_anp_18)-1),
                        median = rep(NA, ncol(draws_anp_18)-1),
                        sd = rep(NA, ncol(draws_anp_18)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_18[,i+1])
  summaries$max[i]    <- max(draws_anp_18[,i+1])
  summaries$mean[i]   <- mean(draws_anp_18[,i+1])
  summaries$median[i] <- median(draws_anp_18[,i+1])
  summaries$sd[i]     <- sd(draws_anp_18[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period18$id_1))+1,
                         NROW(unique(counts_df_period18$id_2))+1,
                         NROW(draws_anp_18)),
                    dimnames = list(c(unique(counts_df_period18$id_1),
                                      unique(counts_df_period18$id_2)[length(unique(counts_df_period18$id_2))]),
                                    c(unique(counts_df_period18$id_1),
                                      unique(counts_df_period18$id_2)[length(unique(counts_df_period18$id_2))]),
                                    NULL))
N <- nrow(counts_df_period18)

for (i in 1:N) {
  dyad_row <- counts_df_period18[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_18[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids18 <- c(unique(counts_df_period18$id_1), unique(counts_df_period18$id_2)[length(unique(counts_df_period18$id_2))])
males18 <- males %>% dplyr::filter(id %in% ids18)

# create variables for different degrees of node connectedness
males18$degree_0.1 <- NA
males18$degree_0.2 <- NA
males18$degree_0.3 <- NA
males18$degree_0.4 <- NA
males18$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males18)){
  rows <- summaries[summaries$id_1 == males18$id[i] | summaries$id_2 == males18$id[i],]
  males18$degree_0.1[i] <- length(which(rows$median > 0.1))
  males18$degree_0.2[i] <- length(which(rows$median > 0.2))
  males18$degree_0.3[i] <- length(which(rows$median > 0.3))
  males18$degree_0.4[i] <- length(which(rows$median > 0.4))
  males18$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males18$age <- lubridate::year(periods$period_start[periods$period == 18]) - males18$byr
summary(males18$age)
males18$age_class <- ifelse(males18$age < 5, 1,
                            ifelse(males18$age < 10, 2,
                                   ifelse(males18$age < 15, 3,
                                          ifelse(males18$age < 20, 4,
                                                 ifelse(males18$age < 25, 5,
                                                        ifelse(males18$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males18$id[which(males18$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males18$id[which(males18$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males18$p18[which(males18$degree_0.3 != 0)],
     vertex.label = males18$id[which(males18$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males18[which(males18$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males18[which(males18$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males18[which(males18$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males18[which(males18$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males18[which(males18$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 19 ################
print(paste0('Period 19 started at ', Sys.time()))

draws_anp_19 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period19_22.04.27.csv')
saveRDS(draws_anp_19, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period19_22.04.27.rds')
draws_anp_19 <- data.matrix(draws_anp_19)

### create data list
counts_df_period19 <- counts_df_non0[counts_df_non0$period == 19,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_19[,2:ncol(draws_anp_19)]),
                        min = rep(NA, ncol(draws_anp_19)-1),
                        max = rep(NA, ncol(draws_anp_19)-1),
                        mean = rep(NA, ncol(draws_anp_19)-1),
                        median = rep(NA, ncol(draws_anp_19)-1),
                        sd = rep(NA, ncol(draws_anp_19)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_19[,i+1])
  summaries$max[i]    <- max(draws_anp_19[,i+1])
  summaries$mean[i]   <- mean(draws_anp_19[,i+1])
  summaries$median[i] <- median(draws_anp_19[,i+1])
  summaries$sd[i]     <- sd(draws_anp_19[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period19$id_1))+1,
                         NROW(unique(counts_df_period19$id_2))+1,
                         NROW(draws_anp_19)),
                    dimnames = list(c(unique(counts_df_period19$id_1),
                                      unique(counts_df_period19$id_2)[length(unique(counts_df_period19$id_2))]),
                                    c(unique(counts_df_period19$id_1),
                                      unique(counts_df_period19$id_2)[length(unique(counts_df_period19$id_2))]),
                                    NULL))
N <- nrow(counts_df_period19)

for (i in 1:N) {
  dyad_row <- counts_df_period19[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_19[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids19 <- c(unique(counts_df_period19$id_1), unique(counts_df_period19$id_2)[length(unique(counts_df_period19$id_2))])
males19 <- males %>% dplyr::filter(id %in% ids19)

# create variables for different degrees of node connectedness
males19$degree_0.1 <- NA
males19$degree_0.2 <- NA
males19$degree_0.3 <- NA
males19$degree_0.4 <- NA
males19$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males19)){
  rows <- summaries[summaries$id_1 == males19$id[i] | summaries$id_2 == males19$id[i],]
  males19$degree_0.1[i] <- length(which(rows$median > 0.1))
  males19$degree_0.2[i] <- length(which(rows$median > 0.2))
  males19$degree_0.3[i] <- length(which(rows$median > 0.3))
  males19$degree_0.4[i] <- length(which(rows$median > 0.4))
  males19$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males19$age <- lubridate::year(periods$period_start[periods$period == 19]) - males19$byr
summary(males19$age)
males19$age_class <- ifelse(males19$age < 5, 1,
                            ifelse(males19$age < 10, 2,
                                   ifelse(males19$age < 15, 3,
                                          ifelse(males19$age < 20, 4,
                                                 ifelse(males19$age < 25, 5,
                                                        ifelse(males19$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males19$id[which(males19$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males19$id[which(males19$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males19$p19[which(males19$degree_0.3 != 0)],
     vertex.label = males19$id[which(males19$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males19[which(males19$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males19[which(males19$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males19[which(males19$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males19[which(males19$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males19[which(males19$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 20 ################
print(paste0('Period 20 started at ', Sys.time()))

draws_anp_20 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period20_22.04.27.csv')
saveRDS(draws_anp_20, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period20_22.04.27.rds')
draws_anp_20 <- data.matrix(draws_anp_20)

### create data list
counts_df_period20 <- counts_df_non0[counts_df_non0$period == 20,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_20[,2:ncol(draws_anp_20)]),
                        min = rep(NA, ncol(draws_anp_20)-1),
                        max = rep(NA, ncol(draws_anp_20)-1),
                        mean = rep(NA, ncol(draws_anp_20)-1),
                        median = rep(NA, ncol(draws_anp_20)-1),
                        sd = rep(NA, ncol(draws_anp_20)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_20[,i+1])
  summaries$max[i]    <- max(draws_anp_20[,i+1])
  summaries$mean[i]   <- mean(draws_anp_20[,i+1])
  summaries$median[i] <- median(draws_anp_20[,i+1])
  summaries$sd[i]     <- sd(draws_anp_20[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period20$id_1))+1,
                         NROW(unique(counts_df_period20$id_2))+1,
                         NROW(draws_anp_20)),
                    dimnames = list(c(unique(counts_df_period20$id_1),
                                      unique(counts_df_period20$id_2)[length(unique(counts_df_period20$id_2))]),
                                    c(unique(counts_df_period20$id_1),
                                      unique(counts_df_period20$id_2)[length(unique(counts_df_period20$id_2))]),
                                    NULL))
N <- nrow(counts_df_period20)

for (i in 1:N) {
  dyad_row <- counts_df_period20[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_20[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids20 <- c(unique(counts_df_period20$id_1), unique(counts_df_period20$id_2)[length(unique(counts_df_period20$id_2))])
males20 <- males %>% dplyr::filter(id %in% ids20)

# create variables for different degrees of node connectedness
males20$degree_0.1 <- NA
males20$degree_0.2 <- NA
males20$degree_0.3 <- NA
males20$degree_0.4 <- NA
males20$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males20)){
  rows <- summaries[summaries$id_1 == males20$id[i] | summaries$id_2 == males20$id[i],]
  males20$degree_0.1[i] <- length(which(rows$median > 0.1))
  males20$degree_0.2[i] <- length(which(rows$median > 0.2))
  males20$degree_0.3[i] <- length(which(rows$median > 0.3))
  males20$degree_0.4[i] <- length(which(rows$median > 0.4))
  males20$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males20$age <- lubridate::year(periods$period_start[periods$period == 20]) - males20$byr
summary(males20$age)
males20$age_class <- ifelse(males20$age < 5, 1,
                            ifelse(males20$age < 10, 2,
                                   ifelse(males20$age < 15, 3,
                                          ifelse(males20$age < 20, 4,
                                                 ifelse(males20$age < 25, 5,
                                                        ifelse(males20$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males20$id[which(males20$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males20$id[which(males20$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males20$p20[which(males20$degree_0.3 != 0)],
     vertex.label = males20$id[which(males20$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males20[which(males20$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males20[which(males20$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males20[which(males20$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males20[which(males20$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males20[which(males20$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 21 ################
print(paste0('Period 21 started at ', Sys.time()))

draws_anp_21 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period21_22.04.27.csv')
saveRDS(draws_anp_21, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period21_22.04.27.rds')
draws_anp_21 <- data.matrix(draws_anp_21)

### create data list
counts_df_period21 <- counts_df_non0[counts_df_non0$period == 21,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_21[,2:ncol(draws_anp_21)]),
                        min = rep(NA, ncol(draws_anp_21)-1),
                        max = rep(NA, ncol(draws_anp_21)-1),
                        mean = rep(NA, ncol(draws_anp_21)-1),
                        median = rep(NA, ncol(draws_anp_21)-1),
                        sd = rep(NA, ncol(draws_anp_21)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_21[,i+1])
  summaries$max[i]    <- max(draws_anp_21[,i+1])
  summaries$mean[i]   <- mean(draws_anp_21[,i+1])
  summaries$median[i] <- median(draws_anp_21[,i+1])
  summaries$sd[i]     <- sd(draws_anp_21[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period21$id_1))+1,
                         NROW(unique(counts_df_period21$id_2))+1,
                         NROW(draws_anp_21)),
                    dimnames = list(c(unique(counts_df_period21$id_1),
                                      unique(counts_df_period21$id_2)[length(unique(counts_df_period21$id_2))]),
                                    c(unique(counts_df_period21$id_1),
                                      unique(counts_df_period21$id_2)[length(unique(counts_df_period21$id_2))]),
                                    NULL))
N <- nrow(counts_df_period21)

for (i in 1:N) {
  dyad_row <- counts_df_period21[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_21[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids21 <- c(unique(counts_df_period21$id_1), unique(counts_df_period21$id_2)[length(unique(counts_df_period21$id_2))])
males21 <- males %>% dplyr::filter(id %in% ids21)

# create variables for different degrees of node connectedness
males21$degree_0.1 <- NA
males21$degree_0.2 <- NA
males21$degree_0.3 <- NA
males21$degree_0.4 <- NA
males21$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males21)){
  rows <- summaries[summaries$id_1 == males21$id[i] | summaries$id_2 == males21$id[i],]
  males21$degree_0.1[i] <- length(which(rows$median > 0.1))
  males21$degree_0.2[i] <- length(which(rows$median > 0.2))
  males21$degree_0.3[i] <- length(which(rows$median > 0.3))
  males21$degree_0.4[i] <- length(which(rows$median > 0.4))
  males21$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males21$age <- lubridate::year(periods$period_start[periods$period == 21]) - males21$byr
summary(males21$age)
males21$age_class <- ifelse(males21$age < 5, 1,
                            ifelse(males21$age < 10, 2,
                                   ifelse(males21$age < 15, 3,
                                          ifelse(males21$age < 20, 4,
                                                 ifelse(males21$age < 25, 5,
                                                        ifelse(males21$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males21$id[which(males21$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males21$id[which(males21$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males21$p21[which(males21$degree_0.3 != 0)],
     vertex.label = males21$id[which(males21$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males21[which(males21$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males21[which(males21$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males21[which(males21$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males21[which(males21$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males21[which(males21$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 22 ################
print(paste0('Period 22 started at ', Sys.time()))

draws_anp_22 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period22_22.04.27.csv')
saveRDS(draws_anp_22, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period22_22.04.27.rds')
draws_anp_22 <- data.matrix(draws_anp_22)

### create data list
counts_df_period22 <- counts_df_non0[counts_df_non0$period == 22,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_22[,2:ncol(draws_anp_22)]),
                        min = rep(NA, ncol(draws_anp_22)-1),
                        max = rep(NA, ncol(draws_anp_22)-1),
                        mean = rep(NA, ncol(draws_anp_22)-1),
                        median = rep(NA, ncol(draws_anp_22)-1),
                        sd = rep(NA, ncol(draws_anp_22)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_22[,i+1])
  summaries$max[i]    <- max(draws_anp_22[,i+1])
  summaries$mean[i]   <- mean(draws_anp_22[,i+1])
  summaries$median[i] <- median(draws_anp_22[,i+1])
  summaries$sd[i]     <- sd(draws_anp_22[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period22$id_1))+1,
                         NROW(unique(counts_df_period22$id_2))+1,
                         NROW(draws_anp_22)),
                    dimnames = list(c(unique(counts_df_period22$id_1),
                                      unique(counts_df_period22$id_2)[length(unique(counts_df_period22$id_2))]),
                                    c(unique(counts_df_period22$id_1),
                                      unique(counts_df_period22$id_2)[length(unique(counts_df_period22$id_2))]),
                                    NULL))
N <- nrow(counts_df_period22)

for (i in 1:N) {
  dyad_row <- counts_df_period22[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_22[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids22 <- c(unique(counts_df_period22$id_1), unique(counts_df_period22$id_2)[length(unique(counts_df_period22$id_2))])
males22 <- males %>% dplyr::filter(id %in% ids22)

# create variables for different degrees of node connectedness
males22$degree_0.1 <- NA
males22$degree_0.2 <- NA
males22$degree_0.3 <- NA
males22$degree_0.4 <- NA
males22$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males22)){
  rows <- summaries[summaries$id_1 == males22$id[i] | summaries$id_2 == males22$id[i],]
  males22$degree_0.1[i] <- length(which(rows$median > 0.1))
  males22$degree_0.2[i] <- length(which(rows$median > 0.2))
  males22$degree_0.3[i] <- length(which(rows$median > 0.3))
  males22$degree_0.4[i] <- length(which(rows$median > 0.4))
  males22$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males22$age <- lubridate::year(periods$period_start[periods$period == 22]) - males22$byr
summary(males22$age)
males22$age_class <- ifelse(males22$age < 5, 1,
                            ifelse(males22$age < 10, 2,
                                   ifelse(males22$age < 15, 3,
                                          ifelse(males22$age < 20, 4,
                                                 ifelse(males22$age < 25, 5,
                                                        ifelse(males22$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males22$id[which(males22$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males22$id[which(males22$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males22$p22[which(males22$degree_0.3 != 0)],
     vertex.label = males22$id[which(males22$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males22[which(males22$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males22[which(males22$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males22[which(males22$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males22[which(males22$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males22[which(males22$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 23 ################
print(paste0('Period 23 started at ', Sys.time()))

draws_anp_23 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period23_22.04.27.csv')
saveRDS(draws_anp_23, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period23_22.04.27.rds')
draws_anp_23 <- data.matrix(draws_anp_23)

### create data list
counts_df_period23 <- counts_df_non0[counts_df_non0$period == 23,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_23[,2:ncol(draws_anp_23)]),
                        min = rep(NA, ncol(draws_anp_23)-1),
                        max = rep(NA, ncol(draws_anp_23)-1),
                        mean = rep(NA, ncol(draws_anp_23)-1),
                        median = rep(NA, ncol(draws_anp_23)-1),
                        sd = rep(NA, ncol(draws_anp_23)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_23[,i+1])
  summaries$max[i]    <- max(draws_anp_23[,i+1])
  summaries$mean[i]   <- mean(draws_anp_23[,i+1])
  summaries$median[i] <- median(draws_anp_23[,i+1])
  summaries$sd[i]     <- sd(draws_anp_23[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period23$id_1))+1,
                         NROW(unique(counts_df_period23$id_2))+1,
                         NROW(draws_anp_23)),
                    dimnames = list(c(unique(counts_df_period23$id_1),
                                      unique(counts_df_period23$id_2)[length(unique(counts_df_period23$id_2))]),
                                    c(unique(counts_df_period23$id_1),
                                      unique(counts_df_period23$id_2)[length(unique(counts_df_period23$id_2))]),
                                    NULL))
N <- nrow(counts_df_period23)

for (i in 1:N) {
  dyad_row <- counts_df_period23[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_23[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids23 <- c(unique(counts_df_period23$id_1), unique(counts_df_period23$id_2)[length(unique(counts_df_period23$id_2))])
males23 <- males %>% dplyr::filter(id %in% ids23)

# create variables for different degrees of node connectedness
males23$degree_0.1 <- NA
males23$degree_0.2 <- NA
males23$degree_0.3 <- NA
males23$degree_0.4 <- NA
males23$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males23)){
  rows <- summaries[summaries$id_1 == males23$id[i] | summaries$id_2 == males23$id[i],]
  males23$degree_0.1[i] <- length(which(rows$median > 0.1))
  males23$degree_0.2[i] <- length(which(rows$median > 0.2))
  males23$degree_0.3[i] <- length(which(rows$median > 0.3))
  males23$degree_0.4[i] <- length(which(rows$median > 0.4))
  males23$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males23$age <- lubridate::year(periods$period_start[periods$period == 23]) - males23$byr
summary(males23$age)
males23$age_class <- ifelse(males23$age < 5, 1,
                            ifelse(males23$age < 10, 2,
                                   ifelse(males23$age < 15, 3,
                                          ifelse(males23$age < 20, 4,
                                                 ifelse(males23$age < 25, 5,
                                                        ifelse(males23$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males23$id[which(males23$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males23$id[which(males23$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males23$p23[which(males23$degree_0.3 != 0)],
     vertex.label = males23$id[which(males23$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males23[which(males23$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males23[which(males23$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males23[which(males23$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males23[which(males23$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males23[which(males23$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 24 ################
print(paste0('Period 24 started at ', Sys.time()))

draws_anp_24 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period24_22.04.27.csv')
saveRDS(draws_anp_24, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period24_22.04.27.rds')
draws_anp_24 <- data.matrix(draws_anp_24)

### create data list
counts_df_period24 <- counts_df_non0[counts_df_non0$period == 24,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_24[,2:ncol(draws_anp_24)]),
                        min = rep(NA, ncol(draws_anp_24)-1),
                        max = rep(NA, ncol(draws_anp_24)-1),
                        mean = rep(NA, ncol(draws_anp_24)-1),
                        median = rep(NA, ncol(draws_anp_24)-1),
                        sd = rep(NA, ncol(draws_anp_24)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_24[,i+1])
  summaries$max[i]    <- max(draws_anp_24[,i+1])
  summaries$mean[i]   <- mean(draws_anp_24[,i+1])
  summaries$median[i] <- median(draws_anp_24[,i+1])
  summaries$sd[i]     <- sd(draws_anp_24[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period24$id_1))+1,
                         NROW(unique(counts_df_period24$id_2))+1,
                         NROW(draws_anp_24)),
                    dimnames = list(c(unique(counts_df_period24$id_1),
                                      unique(counts_df_period24$id_2)[length(unique(counts_df_period24$id_2))]),
                                    c(unique(counts_df_period24$id_1),
                                      unique(counts_df_period24$id_2)[length(unique(counts_df_period24$id_2))]),
                                    NULL))
N <- nrow(counts_df_period24)

for (i in 1:N) {
  dyad_row <- counts_df_period24[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_24[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids24 <- c(unique(counts_df_period24$id_1), unique(counts_df_period24$id_2)[length(unique(counts_df_period24$id_2))])
males24 <- males %>% dplyr::filter(id %in% ids24)

# create variables for different degrees of node connectedness
males24$degree_0.1 <- NA
males24$degree_0.2 <- NA
males24$degree_0.3 <- NA
males24$degree_0.4 <- NA
males24$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males24)){
  rows <- summaries[summaries$id_1 == males24$id[i] | summaries$id_2 == males24$id[i],]
  males24$degree_0.1[i] <- length(which(rows$median > 0.1))
  males24$degree_0.2[i] <- length(which(rows$median > 0.2))
  males24$degree_0.3[i] <- length(which(rows$median > 0.3))
  males24$degree_0.4[i] <- length(which(rows$median > 0.4))
  males24$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males24$age <- lubridate::year(periods$period_start[periods$period == 24]) - males24$byr
summary(males24$age)
males24$age_class <- ifelse(males24$age < 5, 1,
                            ifelse(males24$age < 10, 2,
                                   ifelse(males24$age < 15, 3,
                                          ifelse(males24$age < 20, 4,
                                                 ifelse(males24$age < 25, 5,
                                                        ifelse(males24$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males24$id[which(males24$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males24$id[which(males24$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males24$p24[which(males24$degree_0.3 != 0)],
     vertex.label = males24$id[which(males24$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males24[which(males24$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males24[which(males24$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males24[which(males24$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males24[which(males24$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males24[which(males24$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 25 ################
print(paste0('Period 25 started at ', Sys.time()))

draws_anp_25 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period25_22.04.27.csv')
saveRDS(draws_anp_25, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period25_22.04.27.rds')
draws_anp_25 <- data.matrix(draws_anp_25)

### create data list
counts_df_period25 <- counts_df_non0[counts_df_non0$period == 25,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_25[,2:ncol(draws_anp_25)]),
                        min = rep(NA, ncol(draws_anp_25)-1),
                        max = rep(NA, ncol(draws_anp_25)-1),
                        mean = rep(NA, ncol(draws_anp_25)-1),
                        median = rep(NA, ncol(draws_anp_25)-1),
                        sd = rep(NA, ncol(draws_anp_25)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_25[,i+1])
  summaries$max[i]    <- max(draws_anp_25[,i+1])
  summaries$mean[i]   <- mean(draws_anp_25[,i+1])
  summaries$median[i] <- median(draws_anp_25[,i+1])
  summaries$sd[i]     <- sd(draws_anp_25[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period25$id_1))+1,
                         NROW(unique(counts_df_period25$id_2))+1,
                         NROW(draws_anp_25)),
                    dimnames = list(c(unique(counts_df_period25$id_1),
                                      unique(counts_df_period25$id_2)[length(unique(counts_df_period25$id_2))]),
                                    c(unique(counts_df_period25$id_1),
                                      unique(counts_df_period25$id_2)[length(unique(counts_df_period25$id_2))]),
                                    NULL))
N <- nrow(counts_df_period25)

for (i in 1:N) {
  dyad_row <- counts_df_period25[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_25[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids25 <- c(unique(counts_df_period25$id_1), unique(counts_df_period25$id_2)[length(unique(counts_df_period25$id_2))])
males25 <- males %>% dplyr::filter(id %in% ids25)

# create variables for different degrees of node connectedness
males25$degree_0.1 <- NA
males25$degree_0.2 <- NA
males25$degree_0.3 <- NA
males25$degree_0.4 <- NA
males25$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males25)){
  rows <- summaries[summaries$id_1 == males25$id[i] | summaries$id_2 == males25$id[i],]
  males25$degree_0.1[i] <- length(which(rows$median > 0.1))
  males25$degree_0.2[i] <- length(which(rows$median > 0.2))
  males25$degree_0.3[i] <- length(which(rows$median > 0.3))
  males25$degree_0.4[i] <- length(which(rows$median > 0.4))
  males25$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males25$age <- lubridate::year(periods$period_start[periods$period == 25]) - males25$byr
summary(males25$age)
males25$age_class <- ifelse(males25$age < 5, 1,
                            ifelse(males25$age < 10, 2,
                                   ifelse(males25$age < 15, 3,
                                          ifelse(males25$age < 20, 4,
                                                 ifelse(males25$age < 25, 5,
                                                        ifelse(males25$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males25$id[which(males25$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males25$id[which(males25$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males25$p25[which(males25$degree_0.3 != 0)],
     vertex.label = males25$id[which(males25$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males25[which(males25$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males25[which(males25$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males25[which(males25$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males25[which(males25$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males25[which(males25$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 26 ################
print(paste0('Period 26 started at ', Sys.time()))

draws_anp_26 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period26_22.04.27.csv')
saveRDS(draws_anp_26, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period26_22.04.27.rds')
draws_anp_26 <- data.matrix(draws_anp_26)

### create data list
counts_df_period26 <- counts_df_non0[counts_df_non0$period == 26,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_26[,2:ncol(draws_anp_26)]),
                        min = rep(NA, ncol(draws_anp_26)-1),
                        max = rep(NA, ncol(draws_anp_26)-1),
                        mean = rep(NA, ncol(draws_anp_26)-1),
                        median = rep(NA, ncol(draws_anp_26)-1),
                        sd = rep(NA, ncol(draws_anp_26)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_26[,i+1])
  summaries$max[i]    <- max(draws_anp_26[,i+1])
  summaries$mean[i]   <- mean(draws_anp_26[,i+1])
  summaries$median[i] <- median(draws_anp_26[,i+1])
  summaries$sd[i]     <- sd(draws_anp_26[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period26$id_1))+1,
                         NROW(unique(counts_df_period26$id_2))+1,
                         NROW(draws_anp_26)),
                    dimnames = list(c(unique(counts_df_period26$id_1),
                                      unique(counts_df_period26$id_2)[length(unique(counts_df_period26$id_2))]),
                                    c(unique(counts_df_period26$id_1),
                                      unique(counts_df_period26$id_2)[length(unique(counts_df_period26$id_2))]),
                                    NULL))
N <- nrow(counts_df_period26)

for (i in 1:N) {
  dyad_row <- counts_df_period26[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_26[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids26 <- c(unique(counts_df_period26$id_1), unique(counts_df_period26$id_2)[length(unique(counts_df_period26$id_2))])
males26 <- males %>% dplyr::filter(id %in% ids26)

# create variables for different degrees of node connectedness
males26$degree_0.1 <- NA
males26$degree_0.2 <- NA
males26$degree_0.3 <- NA
males26$degree_0.4 <- NA
males26$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males26)){
  rows <- summaries[summaries$id_1 == males26$id[i] | summaries$id_2 == males26$id[i],]
  males26$degree_0.1[i] <- length(which(rows$median > 0.1))
  males26$degree_0.2[i] <- length(which(rows$median > 0.2))
  males26$degree_0.3[i] <- length(which(rows$median > 0.3))
  males26$degree_0.4[i] <- length(which(rows$median > 0.4))
  males26$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males26$age <- lubridate::year(periods$period_start[periods$period == 26]) - males26$byr
summary(males26$age)
males26$age_class <- ifelse(males26$age < 5, 1,
                            ifelse(males26$age < 10, 2,
                                   ifelse(males26$age < 15, 3,
                                          ifelse(males26$age < 20, 4,
                                                 ifelse(males26$age < 25, 5,
                                                        ifelse(males26$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males26$id[which(males26$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males26$id[which(males26$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males26$p26[which(males26$degree_0.3 != 0)],
     vertex.label = males26$id[which(males26$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males26[which(males26$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males26[which(males26$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males26[which(males26$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males26[which(males26$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males26[which(males26$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 27 ################
print(paste0('Period 27 started at ', Sys.time()))

draws_anp_27 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period27_22.04.27.csv')
saveRDS(draws_anp_27, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period27_22.04.27.rds')
draws_anp_27 <- data.matrix(draws_anp_27)

### create data list
counts_df_period27 <- counts_df_non0[counts_df_non0$period == 27,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_27[,2:ncol(draws_anp_27)]),
                        min = rep(NA, ncol(draws_anp_27)-1),
                        max = rep(NA, ncol(draws_anp_27)-1),
                        mean = rep(NA, ncol(draws_anp_27)-1),
                        median = rep(NA, ncol(draws_anp_27)-1),
                        sd = rep(NA, ncol(draws_anp_27)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_27[,i+1])
  summaries$max[i]    <- max(draws_anp_27[,i+1])
  summaries$mean[i]   <- mean(draws_anp_27[,i+1])
  summaries$median[i] <- median(draws_anp_27[,i+1])
  summaries$sd[i]     <- sd(draws_anp_27[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period27$id_1))+1,
                         NROW(unique(counts_df_period27$id_2))+1,
                         NROW(draws_anp_27)),
                    dimnames = list(c(unique(counts_df_period27$id_1),
                                      unique(counts_df_period27$id_2)[length(unique(counts_df_period27$id_2))]),
                                    c(unique(counts_df_period27$id_1),
                                      unique(counts_df_period27$id_2)[length(unique(counts_df_period27$id_2))]),
                                    NULL))
N <- nrow(counts_df_period27)

for (i in 1:N) {
  dyad_row <- counts_df_period27[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_27[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids27 <- c(unique(counts_df_period27$id_1), unique(counts_df_period27$id_2)[length(unique(counts_df_period27$id_2))])
males27 <- males %>% dplyr::filter(id %in% ids27)

# create variables for different degrees of node connectedness
males27$degree_0.1 <- NA
males27$degree_0.2 <- NA
males27$degree_0.3 <- NA
males27$degree_0.4 <- NA
males27$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males27)){
  rows <- summaries[summaries$id_1 == males27$id[i] | summaries$id_2 == males27$id[i],]
  males27$degree_0.1[i] <- length(which(rows$median > 0.1))
  males27$degree_0.2[i] <- length(which(rows$median > 0.2))
  males27$degree_0.3[i] <- length(which(rows$median > 0.3))
  males27$degree_0.4[i] <- length(which(rows$median > 0.4))
  males27$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males27$age <- lubridate::year(periods$period_start[periods$period == 27]) - males27$byr
summary(males27$age)
males27$age_class <- ifelse(males27$age < 5, 1,
                            ifelse(males27$age < 10, 2,
                                   ifelse(males27$age < 15, 3,
                                          ifelse(males27$age < 20, 4,
                                                 ifelse(males27$age < 25, 5,
                                                        ifelse(males27$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males27$id[which(males27$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males27$id[which(males27$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males27$p27[which(males27$degree_0.3 != 0)],
     vertex.label = males27$id[which(males27$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males27[which(males27$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males27[which(males27$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males27[which(males27$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males27[which(males27$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males27[which(males27$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 28 ################
print(paste0('Period 28 started at ', Sys.time()))

draws_anp_28 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period28_22.04.27.csv')
saveRDS(draws_anp_28, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period28_22.04.27.rds')
draws_anp_28 <- data.matrix(draws_anp_28)

### create data list
counts_df_period28 <- counts_df_non0[counts_df_non0$period == 28,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_28[,2:ncol(draws_anp_28)]),
                        min = rep(NA, ncol(draws_anp_28)-1),
                        max = rep(NA, ncol(draws_anp_28)-1),
                        mean = rep(NA, ncol(draws_anp_28)-1),
                        median = rep(NA, ncol(draws_anp_28)-1),
                        sd = rep(NA, ncol(draws_anp_28)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_28[,i+1])
  summaries$max[i]    <- max(draws_anp_28[,i+1])
  summaries$mean[i]   <- mean(draws_anp_28[,i+1])
  summaries$median[i] <- median(draws_anp_28[,i+1])
  summaries$sd[i]     <- sd(draws_anp_28[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period28$id_1))+1,
                         NROW(unique(counts_df_period28$id_2))+1,
                         NROW(draws_anp_28)),
                    dimnames = list(c(unique(counts_df_period28$id_1),
                                      unique(counts_df_period28$id_2)[length(unique(counts_df_period28$id_2))]),
                                    c(unique(counts_df_period28$id_1),
                                      unique(counts_df_period28$id_2)[length(unique(counts_df_period28$id_2))]),
                                    NULL))
N <- nrow(counts_df_period28)

for (i in 1:N) {
  dyad_row <- counts_df_period28[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_28[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids28 <- c(unique(counts_df_period28$id_1), unique(counts_df_period28$id_2)[length(unique(counts_df_period28$id_2))])
males28 <- males %>% dplyr::filter(id %in% ids28)

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

# age variable
males28$age <- lubridate::year(periods$period_start[periods$period == 28]) - males28$byr
summary(males28$age)
males28$age_class <- ifelse(males28$age < 5, 1,
                            ifelse(males28$age < 10, 2,
                                   ifelse(males28$age < 15, 3,
                                          ifelse(males28$age < 20, 4,
                                                 ifelse(males28$age < 25, 5,
                                                        ifelse(males28$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males28$id[which(males28$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males28$id[which(males28$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males28$p28[which(males28$degree_0.3 != 0)],
     vertex.label = males28$id[which(males28$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males28[which(males28$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males28[which(males28$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males28[which(males28$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males28[which(males28$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males28[which(males28$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)


rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.2,
   counts_df_period28, draws_anp_28, dyad_row, g_mid,g_mid_0.2, g_rng, g_rng_0.2, males28, plot_data_anp_28, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids28, N, plot_cols)

################ Period 29 ################
print(paste0('Period 29 started at ', Sys.time()))

draws_anp_29 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period29_22.07.04.csv')
saveRDS(draws_anp_29, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period29_22.07.04.rds')
draws_anp_29 <- data.matrix(draws_anp_29)

### create data list
counts_df_period29 <- counts_df_non0[counts_df_non0$period == 29,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_29[,2:ncol(draws_anp_29)]),
                        min = rep(NA, ncol(draws_anp_29)-1),
                        max = rep(NA, ncol(draws_anp_29)-1),
                        mean = rep(NA, ncol(draws_anp_29)-1),
                        median = rep(NA, ncol(draws_anp_29)-1),
                        sd = rep(NA, ncol(draws_anp_29)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_29[,i+1])
  summaries$max[i]    <- max(draws_anp_29[,i+1])
  summaries$mean[i]   <- mean(draws_anp_29[,i+1])
  summaries$median[i] <- median(draws_anp_29[,i+1])
  summaries$sd[i]     <- sd(draws_anp_29[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period29$id_1))+1,
                         NROW(unique(counts_df_period29$id_2))+1,
                         NROW(draws_anp_29)),
                    dimnames = list(c(unique(counts_df_period29$id_1),
                                      unique(counts_df_period29$id_2)[length(unique(counts_df_period29$id_2))]),
                                    c(unique(counts_df_period29$id_1),
                                      unique(counts_df_period29$id_2)[length(unique(counts_df_period29$id_2))]),
                                    NULL))
N <- nrow(counts_df_period29)

for (i in 1:N) {
  dyad_row <- counts_df_period29[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_29[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids29 <- c(unique(counts_df_period29$id_1), unique(counts_df_period29$id_2)[length(unique(counts_df_period29$id_2))])
males29 <- males %>% dplyr::filter(id %in% ids29)

# create variables for different degrees of node connectedness
males29$degree_0.1 <- NA
males29$degree_0.2 <- NA
males29$degree_0.3 <- NA
males29$degree_0.4 <- NA
males29$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males29)){
  rows <- summaries[summaries$id_1 == males29$id[i] | summaries$id_2 == males29$id[i],]
  males29$degree_0.1[i] <- length(which(rows$median > 0.1))
  males29$degree_0.2[i] <- length(which(rows$median > 0.2))
  males29$degree_0.3[i] <- length(which(rows$median > 0.3))
  males29$degree_0.4[i] <- length(which(rows$median > 0.4))
  males29$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males29$age <- lubridate::year(periods$period_start[periods$period == 29]) - males29$byr
summary(males29$age)
males29$age_class <- ifelse(males29$age < 5, 1,
                            ifelse(males29$age < 10, 2,
                                   ifelse(males29$age < 15, 3,
                                          ifelse(males29$age < 20, 4,
                                                 ifelse(males29$age < 25, 5,
                                                        ifelse(males29$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males29$id[which(males29$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males29$id[which(males29$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males29$p28[which(males29$degree_0.3 != 0)],
     vertex.label = males29$id[which(males29$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males29[which(males29$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males29[which(males29$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males29[which(males29$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males29[which(males29$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males29[which(males29$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 30 ################
print(paste0('Period 30 started at ', Sys.time()))

draws_anp_30 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period30_22.07.04.csv')
saveRDS(draws_anp_30, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period30_22.07.04.rds')
draws_anp_30 <- data.matrix(draws_anp_30)

### create data list
counts_df_period30 <- counts_df_non0[counts_df_non0$period == 30,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_30[,2:ncol(draws_anp_30)]),
                        min = rep(NA, ncol(draws_anp_30)-1),
                        max = rep(NA, ncol(draws_anp_30)-1),
                        mean = rep(NA, ncol(draws_anp_30)-1),
                        median = rep(NA, ncol(draws_anp_30)-1),
                        sd = rep(NA, ncol(draws_anp_30)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_30[,i+1])
  summaries$max[i]    <- max(draws_anp_30[,i+1])
  summaries$mean[i]   <- mean(draws_anp_30[,i+1])
  summaries$median[i] <- median(draws_anp_30[,i+1])
  summaries$sd[i]     <- sd(draws_anp_30[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period30$id_1))+1,
                         NROW(unique(counts_df_period30$id_2))+1,
                         NROW(draws_anp_30)),
                    dimnames = list(c(unique(counts_df_period30$id_1),
                                      unique(counts_df_period30$id_2)[length(unique(counts_df_period30$id_2))]),
                                    c(unique(counts_df_period30$id_1),
                                      unique(counts_df_period30$id_2)[length(unique(counts_df_period30$id_2))]),
                                    NULL))
N <- nrow(counts_df_period30)

for (i in 1:N) {
  dyad_row <- counts_df_period30[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_30[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids30 <- c(unique(counts_df_period30$id_1), unique(counts_df_period30$id_2)[length(unique(counts_df_period30$id_2))])
males30 <- males %>% dplyr::filter(id %in% ids30)

# create variables for different degrees of node connectedness
males30$degree_0.1 <- NA
males30$degree_0.2 <- NA
males30$degree_0.3 <- NA
males30$degree_0.4 <- NA
males30$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males30)){
  rows <- summaries[summaries$id_1 == males30$id[i] | summaries$id_2 == males30$id[i],]
  males30$degree_0.1[i] <- length(which(rows$median > 0.1))
  males30$degree_0.2[i] <- length(which(rows$median > 0.2))
  males30$degree_0.3[i] <- length(which(rows$median > 0.3))
  males30$degree_0.4[i] <- length(which(rows$median > 0.4))
  males30$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males30$age <- lubridate::year(periods$period_start[periods$period == 30]) - males30$byr
summary(males30$age)
males30$age_class <- ifelse(males30$age < 5, 1,
                            ifelse(males30$age < 10, 2,
                                   ifelse(males30$age < 15, 3,
                                          ifelse(males30$age < 20, 4,
                                                 ifelse(males30$age < 25, 5,
                                                        ifelse(males30$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males30$id[which(males30$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males30$id[which(males30$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males30$p28[which(males30$degree_0.3 != 0)],
     vertex.label = males30$id[which(males30$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males30[which(males30$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males30[which(males30$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males30[which(males30$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males30[which(males30$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males30[which(males30$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

################ Period 31 ################
print(paste0('Period 31 started at ', Sys.time()))

draws_anp_31 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_period31_22.07.04.csv')
saveRDS(draws_anp_31, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period31_22.07.04.rds')
draws_anp_31 <- data.matrix(draws_anp_31)

### create data list
counts_df_period31 <- counts_df_non0[counts_df_non0$period == 31,]

# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_31[,2:ncol(draws_anp_31)]),
                        min = rep(NA, ncol(draws_anp_31)-1),
                        max = rep(NA, ncol(draws_anp_31)-1),
                        mean = rep(NA, ncol(draws_anp_31)-1),
                        median = rep(NA, ncol(draws_anp_31)-1),
                        sd = rep(NA, ncol(draws_anp_31)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_31[,i+1])
  summaries$max[i]    <- max(draws_anp_31[,i+1])
  summaries$mean[i]   <- mean(draws_anp_31[,i+1])
  summaries$median[i] <- median(draws_anp_31[,i+1])
  summaries$sd[i]     <- sd(draws_anp_31[,i+1])
}

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period31$id_1))+1,
                         NROW(unique(counts_df_period31$id_2))+1,
                         NROW(draws_anp_31)),
                    dimnames = list(c(unique(counts_df_period31$id_1),
                                      unique(counts_df_period31$id_2)[length(unique(counts_df_period31$id_2))]),
                                    c(unique(counts_df_period31$id_1),
                                      unique(counts_df_period31$id_2)[length(unique(counts_df_period31$id_2))]),
                                    NULL))
N <- nrow(counts_df_period31)

for (i in 1:N) {
  dyad_row <- counts_df_period31[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_31[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid   <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids31 <- c(unique(counts_df_period31$id_1), unique(counts_df_period31$id_2)[length(unique(counts_df_period31$id_2))])
males31 <- males %>% dplyr::filter(id %in% ids31)

# create variables for different degrees of node connectedness
males31$degree_0.1 <- NA
males31$degree_0.2 <- NA
males31$degree_0.3 <- NA
males31$degree_0.4 <- NA
males31$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males31)){
  rows <- summaries[summaries$id_1 == males31$id[i] | summaries$id_2 == males31$id[i],]
  males31$degree_0.1[i] <- length(which(rows$median > 0.1))
  males31$degree_0.2[i] <- length(which(rows$median > 0.2))
  males31$degree_0.3[i] <- length(which(rows$median > 0.3))
  males31$degree_0.4[i] <- length(which(rows$median > 0.4))
  males31$degree_0.5[i] <- length(which(rows$median > 0.5))
}

# age variable
males31$age <- lubridate::year(periods$period_start[periods$period == 31]) - males31$byr
summary(males31$age)
males31$age_class <- ifelse(males31$age < 5, 1,
                            ifelse(males31$age < 10, 2,
                                   ifelse(males31$age < 15, 3,
                                          ifelse(males31$age < 20, 4,
                                                 ifelse(males31$age < 25, 5,
                                                        ifelse(males31$age < 40, 6, 7))))))

g_mid_0.3 <- delete.vertices(graph = g_mid, v = males31$id[which(males31$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males31$id[which(males31$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)

par(mai = c(0.1,0.4,0.1,0.3))
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males31$p28[which(males31$degree_0.3 != 0)],
     vertex.label = males31$id[which(males31$degree_0.3 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males31[which(males31$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males31[which(males31$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males31[which(males31$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males31[which(males31$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males31[which(males31$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

#### Save outputs to file #####
dev.off()
