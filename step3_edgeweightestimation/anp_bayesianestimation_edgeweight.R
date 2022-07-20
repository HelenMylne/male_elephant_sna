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
library(cmdstanr, lib.loc = 'packages/')    # library(cmdstanr)
# check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
# install_cmdstan(cores = 2)
library(rethinking, lib.loc = 'packages/')  # library(rethinking)
library(igraph, lib.loc = 'packages/')      # library(igraph)
library(dagitty, lib.loc = 'packages/')     # library(dagitty)
library(janitor, lib.loc = 'packages/')     # library(janitor)
library(lubridate, lib.loc = 'packages/')   # library(lubridate)
library(hms, lib.loc = 'packages/')         # library(hms)
library(readxl, lib.loc = 'packages/')      # library(readxl)

# information
sessionInfo()
R.Version()
rstan::stan_version()
#packageVersion('')
#citation('')

# set stan path
set_cmdstan_path('/home/userfs/h/hkm513/.cmdstan/cmdstan-2.29.2')

# load model
mod_2.2 <- cmdstan_model("models/simpleBetaNet_HKM_2.2_22.02.03.stan")
mod_2.2

# set seed
set.seed(12345)

# create file of output graphs
pdf('data_processed/anp_edgeweights_2.2_period1to5_22.05.06.pdf', width = 10, height = 10)

################ 2) Create data lists ################
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
periods <- seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 32)
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

# time windows
periods <- data.frame(period = 1:31,
                      period_start = seq(from = min(sightings$obs_date),
                                         to = max(sightings$obs_date),
                                         length.out = 32)[1:31],
                      period_end = seq(from = min(sightings$obs_date),
                                       to = max(sightings$obs_date),
                                       length.out = 32)[2:32])

# check out total sightings per individual/dyad per period and see if they are reasonable
table(counts_df_non0$period_count_1) ; hist(counts_df_non0$period_count_1)
table(counts_df_non0$period_count_2) ; hist(counts_df_non0$period_count_2)
table(counts_df_non0$period_count_dyad) ; hist(counts_df_non0$period_count_dyad, breaks = 30)
table(counts_df_non0$event_count) ; hist(counts_df_non0$event_count, breaks = 3)
# many are not seen much, but repeat sightings of pairs do seem to be genuinely rare

rm(ate, counts, sightings, i)

################ Time windows 1:5 ################
################ 7.1) Run model on real standardised data -- period 1 ################
### create data list
counts_df_period1 <- counts_df_non0[counts_df_non0$period == 1,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period1),          # total number of times one or other of the dyad was observed
  together = counts_df_period1$event_count,    # count number of sightings seen together
  apart    = counts_df_period1$apart,          # count number of sightings seen apart
  period   = counts_df_period1$period)         # which period it's within

### Fit model
weight_anp_2.2_period1 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period1
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
output1 <- read_cmdstan_csv(weight_anp_2.2_period1$output_files()[1])
draws1_anp_1 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period1$output_files()[2])
draws2_anp_1 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period1$output_files()[3])
draws3_anp_1 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period1$output_files()[4])
draws4_anp_1 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_1 <- rbind(draws1_anp_1, draws2_anp_1, draws3_anp_1, draws4_anp_1)

colnames(draws_anp_1)[2:ncol(draws_anp_1)] <- counts_df_period1$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_1), size = 30, replace = F)

### save data 
write_csv(draws_anp_1, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period1_22.04.15.csv')
rm(draws1_anp_1, draws2_anp_1, draws3_anp_1, draws4_anp_1)

### build traceplots -- period 1 -- very wide ####
plot(draws_anp_1[,plot_cols[1]], type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_1[,plot_cols[2]], col = 'tan')
lines(draws_anp_1[,plot_cols[3]], col = 'orange')
lines(draws_anp_1[,plot_cols[4]], col = 'green')
lines(draws_anp_1[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_1[,plot_cols[6]], col = 'blue')
lines(draws_anp_1[,plot_cols[7]], col = 'red')
lines(draws_anp_1[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_1[,plot_cols[9]], col = 'purple')
lines(draws_anp_1[,plot_cols[10]],col = 'magenta')
lines(draws_anp_1[,plot_cols[11]],col = 'black')
lines(draws_anp_1[,plot_cols[12]], col = 'tan')
lines(draws_anp_1[,plot_cols[13]], col = 'orange')
lines(draws_anp_1[,plot_cols[14]], col = 'green')
lines(draws_anp_1[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_1[,plot_cols[16]], col = 'blue')
lines(draws_anp_1[,plot_cols[17]], col = 'red')
lines(draws_anp_1[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_1[,plot_cols[19]], col = 'purple')
lines(draws_anp_1[,plot_cols[20]],col = 'magenta')
lines(draws_anp_1[,plot_cols[21]],col = 'black')
lines(draws_anp_1[,plot_cols[22]], col = 'tan')
lines(draws_anp_1[,plot_cols[23]], col = 'orange')
lines(draws_anp_1[,plot_cols[24]], col = 'green')
lines(draws_anp_1[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_1[,plot_cols[26]], col = 'blue')
lines(draws_anp_1[,plot_cols[27]], col = 'red')
lines(draws_anp_1[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_1[,plot_cols[29]], col = 'purple')
lines(draws_anp_1[,plot_cols[30]],col = 'magenta')

### density plots -- period 1 ####
dens(draws_anp_1[,2], ylim = c(0,10),xlim = c(0,1), las = 1)
for(i in 1:30){  # looks like a poisson distribution - peaks just above 0 and then exponentiol decline. almost nothing above 0.5
  dens(add = T, draws_anp_1[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 1 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_1 <- left_join(x = summaries, y = counts_df_period1, by = 'dyad')
head(plot_data_anp_1)
plot_data_anp_1$age_cat_1 <- ifelse(plot_data_anp_1$age_start_1 < 6, 'C',
                                    ifelse(plot_data_anp_1$age_start_1 < 11, 'J',
                                           ifelse(plot_data_anp_1$age_start_1 > 19, 'A','P')))
plot_data_anp_1$age_cat_2 <- ifelse(plot_data_anp_1$age_start_2 < 6, 'C',
                                    ifelse(plot_data_anp_1$age_start_2 < 11, 'J',
                                           ifelse(plot_data_anp_1$age_start_2 > 19, 'A','P')))
plot_data_anp_1$age_catnum_1 <- ifelse(plot_data_anp_1$age_start_1 < 6, 1,
                                       ifelse(plot_data_anp_1$age_start_1 < 11, 2,
                                              ifelse(plot_data_anp_1$age_start_1 < 16, 3,
                                                     ifelse(plot_data_anp_1$age_start_1 < 20, 4,
                                                            ifelse(plot_data_anp_1$age_start_1 < 25, 5,
                                                                   ifelse(plot_data_anp_1$age_start_1 < 40, 6, 7))))))
plot_data_anp_1$age_catnum_2 <- ifelse(plot_data_anp_1$age_start_2 < 6, 1,
                                       ifelse(plot_data_anp_1$age_start_2 < 11, 2,
                                              ifelse(plot_data_anp_1$age_start_2 < 16, 3,
                                                     ifelse(plot_data_anp_1$age_start_2 < 20, 4,
                                                            ifelse(plot_data_anp_1$age_start_2 < 25, 5,
                                                                   ifelse(plot_data_anp_1$age_start_2 < 40, 6, 7))))))

plot_data_anp_1$age_dyad <- ifelse(plot_data_anp_1$age_catnum_1 >= plot_data_anp_1$age_catnum_2,
                                   paste(plot_data_anp_1$age_cat_1, plot_data_anp_1$age_cat_2, sep = ''),
                                   paste(plot_data_anp_1$age_cat_2, plot_data_anp_1$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_1, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4, fill = c('blue','purple','blue','grey','purple','blue'))+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_1, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_1[plot_data_anp_1$age_dyad == 'AA' | plot_data_anp_1$age_dyad == 'AP' | plot_data_anp_1$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 1 ################
head(summaries)
length(unique(plot_data_anp_1$id_1))+1 # number of individuals = 55

par(mai = c(0.1,0.1,0.1,0.1))

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
dyad_period_weights <- counts_df_non0
summaries$period <- 1
summaries <- summaries[,c(1,4:9)]
dyad_period_weights <- left_join(x = dyad_period_weights, y = summaries,
                                 by = c('dyad','period'))
rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.2,
   counts_df_period1, draws_anp_1, dyad_row, g_mid,g_mid_0.2, g_rng, g_rng_0.2, males1, plot_data_anp_1, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids1, N, plot_cols)

################ 7.2) Run model on real standardised data -- period 2 ################
### create data list
counts_df_period2 <- counts_df_non0[counts_df_non0$period == 2,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period2),          # total number of times one or other of the dyad was observed
  together = counts_df_period2$event_count,    # count number of sightings seen together
  apart    = counts_df_period2$apart,          # count number of sightings seen apart
  period   = counts_df_period2$period)         # which period it's within

### Fit model
weight_anp_2.2_period2 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period2
# variable       mean     median     sd    mad         q5        q95 rhat ess_bulk ess_tail
#lp__      -5618.76 -5618.45 25.76 25.55 -5661.68 -5576.80 1.00     1505     1915
#weight[1]     0.50     0.50  0.19  0.21     0.20     0.81 1.00     4889     2734
#weight[2]     0.20     0.18  0.12  0.12     0.04     0.43 1.00     3569     2095
#weight[3]     0.28     0.27  0.15  0.17     0.07     0.57 1.00     4288     1961
#weight[4]     0.25     0.23  0.14  0.15     0.05     0.52 1.00     4435     1842
#weight[5]     0.25     0.23  0.15  0.15     0.05     0.52 1.00     4453     2656
#weight[6]     0.20     0.18  0.12  0.12     0.04     0.43 1.00     4252     2100
#weight[7]     0.29     0.26  0.16  0.16     0.07     0.58 1.00     4738     2564
#weight[8]     0.25     0.23  0.14  0.15     0.06     0.52 1.00     5230     2484
#weight[9]     0.22     0.20  0.13  0.14     0.04     0.48 1.00     4733     2047
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period2$output_files()[1])
draws1_anp_2 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period2$output_files()[2])
draws2_anp_2 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period2$output_files()[3])
draws3_anp_2 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period2$output_files()[4])
draws4_anp_2 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_2 <- rbind(draws1_anp_2, draws2_anp_2, draws3_anp_2, draws4_anp_2)

colnames(draws_anp_2)[2:ncol(draws_anp_2)] <- counts_df_period2$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_2), size = 30, replace = F)

### save data 
write_csv(draws_anp_2, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period2_22.04.27.csv')
rm(draws1_anp_2, draws2_anp_2, draws3_anp_2, draws4_anp_2)

### build traceplots  -- period 2 ####
plot(draws_anp_2[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_2[,plot_cols[2]], col = 'tan')
lines(draws_anp_2[,plot_cols[3]], col = 'orange')
lines(draws_anp_2[,plot_cols[4]], col = 'green')
lines(draws_anp_2[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_2[,plot_cols[6]], col = 'blue')
lines(draws_anp_2[,plot_cols[7]], col = 'red')
lines(draws_anp_2[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_2[,plot_cols[9]], col = 'purple')
lines(draws_anp_2[,plot_cols[10]],col = 'magenta')
lines(draws_anp_2[,plot_cols[11]],col = 'black')
lines(draws_anp_2[,plot_cols[12]], col = 'tan')
lines(draws_anp_2[,plot_cols[13]], col = 'orange')
lines(draws_anp_2[,plot_cols[14]], col = 'green')
lines(draws_anp_2[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_2[,plot_cols[16]], col = 'blue')
lines(draws_anp_2[,plot_cols[17]], col = 'red')
lines(draws_anp_2[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_2[,plot_cols[19]], col = 'purple')
lines(draws_anp_2[,plot_cols[20]],col = 'magenta')
lines(draws_anp_2[,plot_cols[21]],col = 'black')
lines(draws_anp_2[,plot_cols[22]], col = 'tan')
lines(draws_anp_2[,plot_cols[23]], col = 'orange')
lines(draws_anp_2[,plot_cols[24]], col = 'green')
lines(draws_anp_2[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_2[,plot_cols[26]], col = 'blue')
lines(draws_anp_2[,plot_cols[27]], col = 'red')
lines(draws_anp_2[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_2[,plot_cols[29]], col = 'purple')
lines(draws_anp_2[,plot_cols[30]],col = 'magenta')

### density plots  -- period 2 ####
dens(draws_anp_2[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:nrow(counts_df_period2)){  # VERY wide, pretty limited
  dens(add = T, draws_anp_2[,i], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 2 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_2 <- left_join(x = summaries, y = counts_df_period2, by = 'dyad')
head(plot_data_anp_2)
plot_data_anp_2$age_cat_1 <- ifelse(plot_data_anp_2$age_start_1 < 6, 'C',
                                    ifelse(plot_data_anp_2$age_start_1 < 11, 'J',
                                           ifelse(plot_data_anp_2$age_start_1 > 19, 'A','P')))
plot_data_anp_2$age_cat_2 <- ifelse(plot_data_anp_2$age_start_2 < 6, 'C',
                                    ifelse(plot_data_anp_2$age_start_2 < 11, 'J',
                                           ifelse(plot_data_anp_2$age_start_2 > 19, 'A','P')))
plot_data_anp_2$age_catnum_1 <- ifelse(plot_data_anp_2$age_start_1 < 6, 1,
                                       ifelse(plot_data_anp_2$age_start_1 < 11, 2,
                                              ifelse(plot_data_anp_2$age_start_1 < 16, 3,
                                                     ifelse(plot_data_anp_2$age_start_1 < 20, 4,
                                                            ifelse(plot_data_anp_2$age_start_1 < 25, 5,
                                                                   ifelse(plot_data_anp_2$age_start_1 < 40, 6, 7))))))
plot_data_anp_2$age_catnum_2 <- ifelse(plot_data_anp_2$age_start_2 < 6, 1,
                                       ifelse(plot_data_anp_2$age_start_2 < 11, 2,
                                              ifelse(plot_data_anp_2$age_start_2 < 16, 3,
                                                     ifelse(plot_data_anp_2$age_start_2 < 20, 4,
                                                            ifelse(plot_data_anp_2$age_start_2 < 25, 5,
                                                                   ifelse(plot_data_anp_2$age_start_2 < 40, 6, 7))))))

plot_data_anp_2$age_dyad <- ifelse(plot_data_anp_2$age_catnum_1 >= plot_data_anp_2$age_catnum_2,
                                   paste(plot_data_anp_2$age_cat_1, plot_data_anp_2$age_cat_2, sep = ''),
                                   paste(plot_data_anp_2$age_cat_2, plot_data_anp_2$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_2, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4, fill = c('blue','purple','purple','blue','grey','purple','purple','blue'))+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_2, aes(x = age_diff, y = mean))+
    geom_jitter(alpha = 0.2)+
    theme_classic()+
    theme(legend.position = 'none')+
    scale_x_continuous('age difference between dyad members',
                       expand = c(0.02,0))+
    scale_y_continuous('mean edge weight',
                       expand = c(0.02,0),
                       limits = c(0,1)))

# values for reporting
summary(summaries) # generally about the same as MOTNP, higher SD
quantile(summaries$median, seq(0,1,length.out = 101))
quantile(summaries$mean,   seq(0,1,length.out = 101))
hist(summaries$median, breaks = 100, xlim = c(0,1), las = 1,
     xlab = 'median estimated dyad strength', main = '', col = 'skyblue')

m_sum <- plot_data_anp_2[plot_data_anp_2$age_dyad == 'AA' | plot_data_anp_2$age_dyad == 'AP' | plot_data_anp_2$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 2 ################
head(summaries)
length(unique(plot_data_anp_2$id_1))+1 # number of individuals = 198

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids2 <- c(unique(counts_df_period2$id_1), unique(counts_df_period2$id_2)[length(unique(counts_df_period2$id_2))])
males2 <- males[,c(21,6,9,22:53)]
males2 <- males2 %>% dplyr::filter(id %in% ids2)
males2

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

which(males2$degree_0.1 < males2$degree_0.2)
which(males2$degree_0.2 < males2$degree_0.3)
which(males2$degree_0.3 < males2$degree_0.4)
which(males2$degree_0.4 < males2$degree_0.5)

# age variable
males2$age <- lubridate::year(periods$period_start[periods$period == 2]) - males2$byr
summary(males2$age)
males2$age_class <- ifelse(males2$age < 5, 1,
                           ifelse(males2$age < 10, 2,
                                  ifelse(males2$age < 15, 3,
                                         ifelse(males2$age < 20, 4,
                                                ifelse(males2$age < 25, 5,
                                                       ifelse(males2$age < 40, 6, 7))))))

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
     vertex.label = males2$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males2$age_class == 7,'seagreen4',
                          ifelse(males2$age_class == 6,'seagreen3',
                                 ifelse(males2$age_class == 5,'seagreen2',
                                        ifelse(males2$age_class == 4,'steelblue3',
                                               ifelse(males2$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males2$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males2$age_class == 7,'seagreen4',
                          ifelse(males2$age_class == 6,'seagreen3',
                                 ifelse(males2$age_class == 5,'seagreen2',
                                        ifelse(males2$age_class == 4,'steelblue3',
                                               ifelse(males2$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 2 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males2$id[which(males2$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males2$id[which(males2$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males2[which(males2$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males2[which(males2$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males2[which(males2$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males2[which(males2$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males2[which(males2$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 2 ####
summaries$period <- 2
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period2, draws_anp_2, dyad_row, g_mid, g_mid_0.3, g_rng, g_rng_0.3, males2, plot_data_anp_2, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids2, N, plot_cols)

################ 7.3) Run model on real standardised data -- period 3 ################
### create data list
counts_df_period3 <- counts_df_non0[counts_df_non0$period == 3,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period3),          # total number of times one or other of the dyad was observed
  together = counts_df_period3$event_count,    # count number of sightings seen together
  apart    = counts_df_period3$apart,          # count number of sightings seen apart
  period   = counts_df_period3$period)         # which period it's within

### Fit model
weight_anp_2.2_period3 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period3
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#lp__      -31403.22 -31403.50 47.19 48.04 -31480.41 -31325.80 1.00     1159     2298
#weight[1]      0.14      0.13  0.07  0.07      0.04      0.27 1.00     8979     2701
#weight[2]      0.13      0.12  0.06  0.06      0.05      0.24 1.00     9758     2510
#weight[3]      0.10      0.08  0.06  0.06      0.02      0.22 1.00     6814     2318
#weight[4]      0.07      0.06  0.04  0.04      0.02      0.15 1.00     7786     2463
#weight[5]      0.09      0.08  0.06  0.06      0.02      0.21 1.00     8305     2564
#weight[6]      0.11      0.09  0.07  0.06      0.02      0.24 1.00     8379     2454
#weight[7]      0.15      0.14  0.07  0.07      0.05      0.29 1.00     8549     2982
#weight[8]      0.16      0.16  0.06  0.06      0.08      0.27 1.00     8564     2599
#weight[9]      0.14      0.13  0.07  0.07      0.05      0.26 1.00     7938     2374
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period3$output_files()[1])
draws1_anp_3 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period3$output_files()[2])
draws2_anp_3 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period3$output_files()[3])
draws3_anp_3 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period3$output_files()[4])
draws4_anp_3 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_3 <- rbind(draws1_anp_3, draws2_anp_3, draws3_anp_3, draws4_anp_3)

colnames(draws_anp_3)[2:ncol(draws_anp_3)] <- counts_df_period3$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 3:ncol(draws_anp_3), size = 30, replace = F)

### save data 
write_csv(draws_anp_3, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period3_22.04.27.csv')
rm(draws1_anp_3, draws2_anp_3, draws3_anp_3, draws4_anp_3)

### build traceplots  -- period 3 -- much narrower than 2 ####
plot(draws_anp_3[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_3[,plot_cols[2]], col = 'tan')
lines(draws_anp_3[,plot_cols[3]], col = 'orange')
lines(draws_anp_3[,plot_cols[4]], col = 'green')
lines(draws_anp_3[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_3[,plot_cols[6]], col = 'blue')
lines(draws_anp_3[,plot_cols[7]], col = 'red')
lines(draws_anp_3[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_3[,plot_cols[9]], col = 'purple')
lines(draws_anp_3[,plot_cols[10]],col = 'magenta')
lines(draws_anp_3[,plot_cols[11]],col = 'black')
lines(draws_anp_3[,plot_cols[12]], col = 'tan')
lines(draws_anp_3[,plot_cols[13]], col = 'orange')
lines(draws_anp_3[,plot_cols[14]], col = 'green')
lines(draws_anp_3[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_3[,plot_cols[16]], col = 'blue')
lines(draws_anp_3[,plot_cols[17]], col = 'red')
lines(draws_anp_3[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_3[,plot_cols[19]], col = 'purple')
lines(draws_anp_3[,plot_cols[20]],col = 'magenta')
lines(draws_anp_3[,plot_cols[21]],col = 'black')
lines(draws_anp_3[,plot_cols[22]], col = 'tan')
lines(draws_anp_3[,plot_cols[23]], col = 'orange')
lines(draws_anp_3[,plot_cols[24]], col = 'green')
lines(draws_anp_3[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_3[,plot_cols[26]], col = 'blue')
lines(draws_anp_3[,plot_cols[27]], col = 'red')
lines(draws_anp_3[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_3[,plot_cols[29]], col = 'purple')
lines(draws_anp_3[,plot_cols[30]],col = 'magenta')

### density plots  -- period 3 ####
dens(draws_anp_3[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:nrow(counts_df_period3)){
  dens(add = T, draws_anp_3[,i], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 3 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_3 <- left_join(x = summaries, y = counts_df_period3, by = 'dyad')
head(plot_data_anp_3)
plot_data_anp_3$age_cat_1 <- ifelse(plot_data_anp_3$age_start_1 < 6, 'C',
                                    ifelse(plot_data_anp_3$age_start_1 < 11, 'J',
                                           ifelse(plot_data_anp_3$age_start_1 > 19, 'A','P')))
plot_data_anp_3$age_cat_2 <- ifelse(plot_data_anp_3$age_start_2 < 6, 'C',
                                    ifelse(plot_data_anp_3$age_start_2 < 11, 'J',
                                           ifelse(plot_data_anp_3$age_start_2 > 19, 'A','P')))
plot_data_anp_3$age_catnum_1 <- ifelse(plot_data_anp_3$age_start_1 < 6, 1,
                                       ifelse(plot_data_anp_3$age_start_1 < 11, 3,
                                              ifelse(plot_data_anp_3$age_start_1 < 16, 3,
                                                     ifelse(plot_data_anp_3$age_start_1 < 20, 4,
                                                            ifelse(plot_data_anp_3$age_start_1 < 25, 5,
                                                                   ifelse(plot_data_anp_3$age_start_1 < 40, 6, 7))))))
plot_data_anp_3$age_catnum_2 <- ifelse(plot_data_anp_3$age_start_2 < 6, 1,
                                       ifelse(plot_data_anp_3$age_start_2 < 11, 3,
                                              ifelse(plot_data_anp_3$age_start_2 < 16, 3,
                                                     ifelse(plot_data_anp_3$age_start_2 < 20, 4,
                                                            ifelse(plot_data_anp_3$age_start_2 < 25, 5,
                                                                   ifelse(plot_data_anp_3$age_start_2 < 40, 6, 7))))))

plot_data_anp_3$age_dyad <- ifelse(plot_data_anp_3$age_catnum_1 >= plot_data_anp_3$age_catnum_2,
                                   paste(plot_data_anp_3$age_cat_1, plot_data_anp_3$age_cat_2, sep = ''),
                                   paste(plot_data_anp_3$age_cat_2, plot_data_anp_3$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_3, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4, fill = c('blue','purple','purple','blue','grey','purple','purple','blue'))+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_3, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_3[plot_data_anp_3$age_dyad == 'AA' | plot_data_anp_3$age_dyad == 'AP' | plot_data_anp_3$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 3 ################
head(summaries)
length(unique(plot_data_anp_3$id_1))+1 # number of individuals = 90

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids3 <- c(unique(counts_df_period3$id_1), unique(counts_df_period3$id_2)[length(unique(counts_df_period3$id_2))])
males3 <- males[,c(21,6,9,22:53)]
males3 <- males3 %>% dplyr::filter(id %in% ids3)
males3

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

which(males3$degree_0.1 < males3$degree_0.2)
which(males3$degree_0.2 < males3$degree_0.3)
which(males3$degree_0.3 < males3$degree_0.4)
which(males3$degree_0.4 < males3$degree_0.5)

# age variable
males3$age <- lubridate::year(periods$period_start[periods$period == 3]) - males3$byr
summary(males3$age)
males3$age_class <- ifelse(males3$age < 5, 1,
                           ifelse(males3$age < 10, 3,
                                  ifelse(males3$age < 15, 3,
                                         ifelse(males3$age < 20, 4,
                                                ifelse(males3$age < 25, 5,
                                                       ifelse(males3$age < 40, 6, 7))))))

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
     vertex.label = males3$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males3$age_class == 7,'seagreen4',
                          ifelse(males3$age_class == 6,'seagreen3',
                                 ifelse(males3$age_class == 5,'seagreen2',
                                        ifelse(males3$age_class == 4,'steelblue3',
                                               ifelse(males3$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males3$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males3$age_class == 7,'seagreen4',
                          ifelse(males3$age_class == 6,'seagreen3',
                                 ifelse(males3$age_class == 5,'seagreen2',
                                        ifelse(males3$age_class == 4,'steelblue3',
                                               ifelse(males3$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 3 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males3$id[which(males3$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males3$id[which(males3$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males3[which(males3$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males3[which(males3$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males3[which(males3$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males3[which(males3$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males3[which(males3$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 3 ####
summaries$period <- 3
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
   counts_df_period3, draws_anp_3, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males3, plot_data_anp_3, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids3, N, plot_cols)


################ 7.4) Run model on real standardised data -- period 4 ################
### create data list
counts_df_period4 <- counts_df_non0[counts_df_non0$period == 4,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period4),          # total number of times one or other of the dyad was observed
  together = counts_df_period4$event_count,    # count number of sightings seen together
  apart    = counts_df_period4$apart,          # count number of sightings seen apart
  period   = counts_df_period4$period)         # which period it's within

### Fit model
weight_anp_2.2_period4 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period4
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#lp__      -82579.29 -82577.70 56.61 56.78 -82674.02 -82488.10 1.00     1286     2421
#weight[1]      0.06      0.06  0.03  0.03      0.02      0.12 1.00     9213     2523
#weight[2]      0.03      0.03  0.02  0.02      0.01      0.07 1.00     6821     2281
#weight[3]      0.08      0.08  0.03  0.03      0.04      0.13 1.00     7485     2724
#weight[4]      0.06      0.06  0.04  0.03      0.02      0.13 1.00     7905     1884
#weight[5]      0.07      0.07  0.03  0.03      0.02      0.14 1.00     8006     2700
#weight[6]      0.06      0.06  0.03  0.03      0.02      0.12 1.00     7862     2694
#weight[7]      0.06      0.06  0.03  0.03      0.02      0.12 1.00     7811     2493
#weight[8]      0.06      0.06  0.03  0.03      0.02      0.12 1.00     7067     2590
#weight[9]      0.09      0.09  0.04  0.04      0.04      0.15 1.00     7780     2979
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period4$output_files()[1])
draws1_anp_4 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period4$output_files()[2])
draws2_anp_4 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period4$output_files()[3])
draws3_anp_4 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period4$output_files()[4])
draws4_anp_4 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_4 <- rbind(draws1_anp_4, draws2_anp_4, draws3_anp_4, draws4_anp_4)

colnames(draws_anp_4)[2:ncol(draws_anp_4)] <- counts_df_period4$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_4), size = 30, replace = F)

### save data 
write_csv(draws_anp_4, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period4_22.04.27.csv')
rm(draws1_anp_4, draws2_anp_4, draws3_anp_4, draws4_anp_4)

### build traceplots  -- period 4 -- very tight and low ####
plot(draws_anp_4[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_4[,plot_cols[2]], col = 'tan')
lines(draws_anp_4[,plot_cols[3]], col = 'orange')
lines(draws_anp_4[,plot_cols[4]], col = 'green')
lines(draws_anp_4[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_4[,plot_cols[6]], col = 'blue')
lines(draws_anp_4[,plot_cols[7]], col = 'red')
lines(draws_anp_4[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_4[,plot_cols[9]], col = 'purple')
lines(draws_anp_4[,plot_cols[10]],col = 'magenta')
lines(draws_anp_4[,plot_cols[11]],col = 'black')
lines(draws_anp_4[,plot_cols[12]], col = 'tan')
lines(draws_anp_4[,plot_cols[13]], col = 'orange')
lines(draws_anp_4[,plot_cols[14]], col = 'green')
lines(draws_anp_4[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_4[,plot_cols[16]], col = 'blue')
lines(draws_anp_4[,plot_cols[17]], col = 'red')
lines(draws_anp_4[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_4[,plot_cols[19]], col = 'purple')
lines(draws_anp_4[,plot_cols[20]],col = 'magenta')
lines(draws_anp_4[,plot_cols[21]],col = 'black')
lines(draws_anp_4[,plot_cols[22]], col = 'tan')
lines(draws_anp_4[,plot_cols[23]], col = 'orange')
lines(draws_anp_4[,plot_cols[24]], col = 'green')
lines(draws_anp_4[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_4[,plot_cols[26]], col = 'blue')
lines(draws_anp_4[,plot_cols[27]], col = 'red')
lines(draws_anp_4[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_4[,plot_cols[29]], col = 'purple')
lines(draws_anp_4[,plot_cols[30]],col = 'magenta')

### density plots  -- period 4 ####
dens(draws_anp_4[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_4[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 4 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_4 <- left_join(x = summaries, y = counts_df_period4, by = 'dyad')
head(plot_data_anp_4)
plot_data_anp_4$age_cat_1 <- ifelse(plot_data_anp_4$age_start_1 < 6, 'C',
                                    ifelse(plot_data_anp_4$age_start_1 < 11, 'J',
                                           ifelse(plot_data_anp_4$age_start_1 > 19, 'A','P')))
plot_data_anp_4$age_cat_2 <- ifelse(plot_data_anp_4$age_start_2 < 6, 'C',
                                    ifelse(plot_data_anp_4$age_start_2 < 11, 'J',
                                           ifelse(plot_data_anp_4$age_start_2 > 19, 'A','P')))
plot_data_anp_4$age_catnum_1 <- ifelse(plot_data_anp_4$age_start_1 < 6, 1,
                                       ifelse(plot_data_anp_4$age_start_1 < 11, 2,
                                              ifelse(plot_data_anp_4$age_start_1 < 16, 3,
                                                     ifelse(plot_data_anp_4$age_start_1 < 20, 4,
                                                            ifelse(plot_data_anp_4$age_start_1 < 25, 5,
                                                                   ifelse(plot_data_anp_4$age_start_1 < 40, 6, 7))))))
plot_data_anp_4$age_catnum_2 <- ifelse(plot_data_anp_4$age_start_2 < 6, 1,
                                       ifelse(plot_data_anp_4$age_start_2 < 11, 2,
                                              ifelse(plot_data_anp_4$age_start_2 < 16, 3,
                                                     ifelse(plot_data_anp_4$age_start_2 < 20, 4,
                                                            ifelse(plot_data_anp_4$age_start_2 < 25, 5,
                                                                   ifelse(plot_data_anp_4$age_start_2 < 40, 6, 7))))))

plot_data_anp_4$age_dyad <- ifelse(plot_data_anp_4$age_catnum_1 >= plot_data_anp_4$age_catnum_2,
                                   paste(plot_data_anp_4$age_cat_1, plot_data_anp_4$age_cat_2, sep = ''),
                                   paste(plot_data_anp_4$age_cat_2, plot_data_anp_4$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_4, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_4, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_4[plot_data_anp_4$age_dyad == 'AA' | plot_data_anp_4$age_dyad == 'AP' | plot_data_anp_4$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 4 ################
head(summaries)
length(unique(plot_data_anp_4$id_1))+1 # number of individuals = 110

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids4 <- c(unique(counts_df_period4$id_1), unique(counts_df_period4$id_2)[length(unique(counts_df_period4$id_2))])
males4 <- males[,c(21,6,9,22:53)]
males4 <- males4 %>% dplyr::filter(id %in% ids4)
males4

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

which(males4$degree_0.1 < males4$degree_0.2)
which(males4$degree_0.2 < males4$degree_0.3)
which(males4$degree_0.3 < males4$degree_0.4)
which(males4$degree_0.4 < males4$degree_0.5)

# age variable
males4$age <- lubridate::year(periods$period_start[periods$period == 4]) - males4$byr
summary(males4$age)
males4$age_class <- ifelse(males4$age < 5, 1,
                           ifelse(males4$age < 10, 2,
                                  ifelse(males4$age < 15, 3,
                                         ifelse(males4$age < 20, 4,
                                                ifelse(males4$age < 25, 5,
                                                       ifelse(males4$age < 40, 6, 7))))))

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
     vertex.label = males4$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males4$age_class == 7,'seagreen4',
                          ifelse(males4$age_class == 6,'seagreen3',
                                 ifelse(males4$age_class == 5,'seagreen2',
                                        ifelse(males4$age_class == 4,'steelblue3',
                                               ifelse(males4$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males4$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males4$age_class == 7,'seagreen4',
                          ifelse(males4$age_class == 6,'seagreen3',
                                 ifelse(males4$age_class == 5,'seagreen2',
                                        ifelse(males4$age_class == 4,'steelblue3',
                                               ifelse(males4$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.2 -- period 4 ####
g_mid_0.2 <- delete.vertices(graph = g_mid, v = males4$id[which(males4$degree_0.2 == 0)])
g_rng_0.2 <- delete.vertices(graph = g_rng, v = males4$id[which(males4$degree_0.2 == 0)])

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
     vertex.color= ifelse(males4[which(males4$degree_0.2 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males4[which(males4$degree_0.2 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males4[which(males4$degree_0.2 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males4[which(males4$degree_0.2 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males4[which(males4$degree_0.2 != 0),]$age_class == 3,'steelblue1',
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

### save summary data -- period 4 ####
summaries$period <- 4
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
   counts_df_period4, draws_anp_4, dyad_row, g_mid,g_mid_0.2, g_rng, g_rng_0.2, males4, plot_data_anp_4, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids4, N, plot_cols)

################ 7.5) Run model on real standardised data -- period 5 ################
### create data list
counts_df_period5 <- counts_df_non0[counts_df_non0$period == 5,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period5),          # total number of times one or other of the dyad was observed
  together = counts_df_period5$event_count,    # count number of sightings seen together
  apart    = counts_df_period5$apart,          # count number of sightings seen apart
  period   = counts_df_period5$period)         # which period it's within

### Fit model
weight_anp_2.2_period5 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period5
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#lp__      -77204.58 -77204.95 66.51 66.12 -77313.62 -77091.90 1.00     1269     2401
#weight[1]      0.05      0.04  0.03  0.03      0.01      0.11 1.00     6298     2354
#weight[2]      0.04      0.03  0.03  0.02      0.01      0.09 1.00     5475     1898
#weight[3]      0.07      0.06  0.04  0.03      0.02      0.14 1.00     5436     2391
#weight[4]      0.05      0.04  0.03  0.03      0.01      0.11 1.00     5694     2612
#weight[5]      0.05      0.04  0.02  0.02      0.01      0.09 1.00     5260     2244
#weight[6]      0.06      0.06  0.04  0.03      0.02      0.13 1.00     6084     2049
#weight[7]      0.03      0.03  0.02  0.02      0.01      0.07 1.00     5936     2185
#weight[8]      0.05      0.05  0.03  0.03      0.02      0.11 1.00     6178     2009
#weight[9]      0.04      0.03  0.02  0.02      0.01      0.08 1.00     5983     2242

rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period5$output_files()[1])
draws1_anp_5 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period5$output_files()[2])
draws2_anp_5 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period5$output_files()[3])
draws3_anp_5 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period5$output_files()[4])
draws4_anp_5 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_5 <- rbind(draws1_anp_5, draws2_anp_5, draws3_anp_5, draws4_anp_5)

colnames(draws_anp_5)[2:ncol(draws_anp_5)] <- counts_df_period5$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_5), size = 30, replace = F)

### save data 
write_csv(draws_anp_5, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period5_22.04.27.csv')
rm(draws1_anp_5, draws2_anp_5, draws3_anp_5, draws4_anp_5)

### build traceplots  -- period 5 ####
plot(draws_anp_5[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_5[,plot_cols[2]], col = 'tan')
lines(draws_anp_5[,plot_cols[3]], col = 'orange')
lines(draws_anp_5[,plot_cols[4]], col = 'green')
lines(draws_anp_5[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_5[,plot_cols[6]], col = 'blue')
lines(draws_anp_5[,plot_cols[7]], col = 'red')
lines(draws_anp_5[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_5[,plot_cols[9]], col = 'purple')
lines(draws_anp_5[,plot_cols[10]],col = 'magenta')
lines(draws_anp_5[,plot_cols[11]],col = 'black')
lines(draws_anp_5[,plot_cols[12]], col = 'tan')
lines(draws_anp_5[,plot_cols[13]], col = 'orange')
lines(draws_anp_5[,plot_cols[14]], col = 'green')
lines(draws_anp_5[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_5[,plot_cols[16]], col = 'blue')
lines(draws_anp_5[,plot_cols[17]], col = 'red')
lines(draws_anp_5[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_5[,plot_cols[19]], col = 'purple')
lines(draws_anp_5[,plot_cols[20]],col = 'magenta')
lines(draws_anp_5[,plot_cols[21]],col = 'black')
lines(draws_anp_5[,plot_cols[22]], col = 'tan')
lines(draws_anp_5[,plot_cols[23]], col = 'orange')
lines(draws_anp_5[,plot_cols[24]], col = 'green')
lines(draws_anp_5[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_5[,plot_cols[26]], col = 'blue')
lines(draws_anp_5[,plot_cols[27]], col = 'red')
lines(draws_anp_5[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_5[,plot_cols[29]], col = 'purple')
lines(draws_anp_5[,plot_cols[30]],col = 'magenta')

### density plots  -- period 5 ####
dens(draws_anp_5[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_5[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 5 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_5 <- left_join(x = summaries, y = counts_df_period5, by = 'dyad')
head(plot_data_anp_5)
plot_data_anp_5$age_cat_1 <- ifelse(plot_data_anp_5$age_start_1 < 6, 'C',
                                    ifelse(plot_data_anp_5$age_start_1 < 11, 'J',
                                           ifelse(plot_data_anp_5$age_start_1 > 19, 'A','P')))
plot_data_anp_5$age_cat_2 <- ifelse(plot_data_anp_5$age_start_2 < 6, 'C',
                                    ifelse(plot_data_anp_5$age_start_2 < 11, 'J',
                                           ifelse(plot_data_anp_5$age_start_2 > 19, 'A','P')))
plot_data_anp_5$age_catnum_1 <- ifelse(plot_data_anp_5$age_start_1 < 6, 1,
                                       ifelse(plot_data_anp_5$age_start_1 < 11, 2,
                                              ifelse(plot_data_anp_5$age_start_1 < 16, 3,
                                                     ifelse(plot_data_anp_5$age_start_1 < 20, 4,
                                                            ifelse(plot_data_anp_5$age_start_1 < 25, 5,
                                                                   ifelse(plot_data_anp_5$age_start_1 < 40, 6, 7))))))
plot_data_anp_5$age_catnum_2 <- ifelse(plot_data_anp_5$age_start_2 < 6, 1,
                                       ifelse(plot_data_anp_5$age_start_2 < 11, 2,
                                              ifelse(plot_data_anp_5$age_start_2 < 16, 3,
                                                     ifelse(plot_data_anp_5$age_start_2 < 20, 4,
                                                            ifelse(plot_data_anp_5$age_start_2 < 25, 5,
                                                                   ifelse(plot_data_anp_5$age_start_2 < 40, 6, 7))))))

plot_data_anp_5$age_dyad <- ifelse(plot_data_anp_5$age_catnum_1 >= plot_data_anp_5$age_catnum_2,
                                   paste(plot_data_anp_5$age_cat_1, plot_data_anp_5$age_cat_2, sep = ''),
                                   paste(plot_data_anp_5$age_cat_2, plot_data_anp_5$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_5, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_5, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_5[plot_data_anp_5$age_dyad == 'AA' | plot_data_anp_5$age_dyad == 'AP' | plot_data_anp_5$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 5 ################
head(summaries)
length(unique(plot_data_anp_5$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids5 <- c(unique(counts_df_period5$id_1), unique(counts_df_period5$id_2)[length(unique(counts_df_period5$id_2))])
males5 <- males[,c(21,6,9,22:53)]
males5 <- males5 %>% dplyr::filter(id %in% ids5)
males5

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

which(males5$degree_0.1 < males5$degree_0.2)
which(males5$degree_0.2 < males5$degree_0.3)
which(males5$degree_0.3 < males5$degree_0.4)
which(males5$degree_0.4 < males5$degree_0.5)

# age variable
males5$age <- lubridate::year(periods$period_start[periods$period == 5]) - males5$byr
summary(males5$age)
males5$age_class <- ifelse(males5$age < 5, 1,
                           ifelse(males5$age < 10, 2,
                                  ifelse(males5$age < 15, 3,
                                         ifelse(males5$age < 20, 4,
                                                ifelse(males5$age < 25, 5,
                                                       ifelse(males5$age < 40, 6, 7))))))

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
     vertex.label = males5$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males5$age_class == 7,'seagreen4',
                          ifelse(males5$age_class == 6,'seagreen3',
                                 ifelse(males5$age_class == 5,'seagreen2',
                                        ifelse(males5$age_class == 4,'steelblue3',
                                               ifelse(males5$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males5$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males5$age_class == 7,'seagreen4',
                          ifelse(males5$age_class == 6,'seagreen3',
                                 ifelse(males5$age_class == 5,'seagreen2',
                                        ifelse(males5$age_class == 4,'steelblue3',
                                               ifelse(males5$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.2 -- period 5 ####
g_mid_0.2 <- delete.vertices(graph = g_mid, v = males5$id[which(males5$degree_0.2 == 0)])
g_rng_0.2 <- delete.vertices(graph = g_rng, v = males5$id[which(males5$degree_0.2 == 0)])

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
     vertex.color= ifelse(males5[which(males5$degree_0.2 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males5[which(males5$degree_0.2 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males5[which(males5$degree_0.2 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males5[which(males5$degree_0.2 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males5[which(males5$degree_0.2 != 0),]$age_class == 3,'steelblue1',
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

### save summary data -- period 5 ####
summaries$period <- 5
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period5, draws_anp_5, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males5, plot_data_anp_5, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids5, N, plot_cols)

#### Save outputs to file #####
# write out csv file
write_csv(dyad_period_weights, 'data_processed/anp_dyad_weightdistributions_2.2_periods1to5_22.05.06.csv')

# produce PDF of graphs
dev.off()

# add progress stamp
print(paste0('Time windows 1:5 completed at ', Sys.time()))

################ Time windows 6:10 ################
# set seed
set.seed(12345)

# create file of output graphs
pdf('data_processed/anp_edgeweights_2.2_period6to10_22.05.06.pdf', width = 10, height = 10)

# read in previous summary file
#dyad_period_weights <- read_csv('data_processed/anp_dyad_weightdistributions_2.2_periods1to5_22.05.06.csv')

################ 7.6) Run model on real standardised data -- period 6 ################
### create data list
counts_df_period6 <- counts_df_non0[counts_df_non0$period == 6,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period6),          # total number of times one or other of the dyad was observed
  together = counts_df_period6$event_count,    # count number of sightings seen together
  apart    = counts_df_period6$apart,          # count number of sightings seen apart
  period   = counts_df_period6$period)         # which period it's within

### Fit model
weight_anp_2.2_period6 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period6
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#lp__      -91497.34 -91497.10 66.47 67.53 -91607.90 -91390.40 1.00     1445     2249
#weight[1]      0.04      0.04  0.03  0.02      0.01      0.09 1.00     6590     2588
#weight[2]      0.03      0.03  0.02  0.02      0.01      0.07 1.00     5221     2732
#weight[3]      0.13      0.13  0.04  0.03      0.08      0.20 1.00     7163     2879
#weight[4]      0.10      0.10  0.03  0.03      0.05      0.17 1.00     6325     2452
#weight[5]      0.16      0.16  0.04  0.04      0.11      0.23 1.00     7535     2858
#weight[6]      0.12      0.11  0.04  0.04      0.06      0.18 1.00     5362     2364
#weight[7]      0.04      0.04  0.02  0.02      0.01      0.08 1.00     6300     2560
#weight[8]      0.06      0.06  0.03  0.02      0.03      0.11 1.00     6350     2447
#weight[9]      0.06      0.05  0.03  0.03      0.02      0.11 1.00     7642     2693
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period6$output_files()[1])
draws1_anp_6 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period6$output_files()[2])
draws2_anp_6 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period6$output_files()[3])
draws3_anp_6 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period6$output_files()[4])
draws4_anp_6 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_6 <- rbind(draws1_anp_6, draws2_anp_6, draws3_anp_6, draws4_anp_6)

colnames(draws_anp_6)[2:ncol(draws_anp_6)] <- counts_df_period6$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_6), size = 30, replace = F)

### save data 
write_csv(draws_anp_6, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period6_22.04.27.csv')
rm(draws1_anp_6, draws2_anp_6, draws3_anp_6, draws4_anp_6)

### build traceplots  -- period 6 ####
plot(draws_anp_6[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_6[,plot_cols[2]], col = 'tan')
lines(draws_anp_6[,plot_cols[3]], col = 'orange')
lines(draws_anp_6[,plot_cols[4]], col = 'green')
lines(draws_anp_6[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_6[,plot_cols[6]], col = 'blue')
lines(draws_anp_6[,plot_cols[7]], col = 'red')
lines(draws_anp_6[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_6[,plot_cols[9]], col = 'purple')
lines(draws_anp_6[,plot_cols[10]],col = 'magenta')
lines(draws_anp_6[,plot_cols[11]],col = 'black')
lines(draws_anp_6[,plot_cols[12]], col = 'tan')
lines(draws_anp_6[,plot_cols[13]], col = 'orange')
lines(draws_anp_6[,plot_cols[14]], col = 'green')
lines(draws_anp_6[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_6[,plot_cols[16]], col = 'blue')
lines(draws_anp_6[,plot_cols[17]], col = 'red')
lines(draws_anp_6[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_6[,plot_cols[19]], col = 'purple')
lines(draws_anp_6[,plot_cols[20]],col = 'magenta')
lines(draws_anp_6[,plot_cols[21]],col = 'black')
lines(draws_anp_6[,plot_cols[22]], col = 'tan')
lines(draws_anp_6[,plot_cols[23]], col = 'orange')
lines(draws_anp_6[,plot_cols[24]], col = 'green')
lines(draws_anp_6[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_6[,plot_cols[26]], col = 'blue')
lines(draws_anp_6[,plot_cols[27]], col = 'red')
lines(draws_anp_6[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_6[,plot_cols[29]], col = 'purple')
lines(draws_anp_6[,plot_cols[30]],col = 'magenta')

### density plots  -- period 6 ####
dens(draws_anp_6[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_6[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 6 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_6 <- left_join(x = summaries, y = counts_df_period6, by = 'dyad')
head(plot_data_anp_6)
plot_data_anp_6$age_cat_1 <- ifelse(plot_data_anp_6$age_start_1 < 6, 'C',
                                    ifelse(plot_data_anp_6$age_start_1 < 11, 'J',
                                           ifelse(plot_data_anp_6$age_start_1 > 19, 'A','P')))
plot_data_anp_6$age_cat_2 <- ifelse(plot_data_anp_6$age_start_2 < 6, 'C',
                                    ifelse(plot_data_anp_6$age_start_2 < 11, 'J',
                                           ifelse(plot_data_anp_6$age_start_2 > 19, 'A','P')))
plot_data_anp_6$age_catnum_1 <- ifelse(plot_data_anp_6$age_start_1 < 6, 1,
                                       ifelse(plot_data_anp_6$age_start_1 < 11, 2,
                                              ifelse(plot_data_anp_6$age_start_1 < 16, 3,
                                                     ifelse(plot_data_anp_6$age_start_1 < 20, 4,
                                                            ifelse(plot_data_anp_6$age_start_1 < 25, 5,
                                                                   ifelse(plot_data_anp_6$age_start_1 < 40, 6, 7))))))
plot_data_anp_6$age_catnum_2 <- ifelse(plot_data_anp_6$age_start_2 < 6, 1,
                                       ifelse(plot_data_anp_6$age_start_2 < 11, 2,
                                              ifelse(plot_data_anp_6$age_start_2 < 16, 3,
                                                     ifelse(plot_data_anp_6$age_start_2 < 20, 4,
                                                            ifelse(plot_data_anp_6$age_start_2 < 25, 5,
                                                                   ifelse(plot_data_anp_6$age_start_2 < 40, 6, 7))))))

plot_data_anp_6$age_dyad <- ifelse(plot_data_anp_6$age_catnum_1 >= plot_data_anp_6$age_catnum_2,
                                   paste(plot_data_anp_6$age_cat_1, plot_data_anp_6$age_cat_2, sep = ''),
                                   paste(plot_data_anp_6$age_cat_2, plot_data_anp_6$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_6, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_6, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_6[plot_data_anp_6$age_dyad == 'AA' | plot_data_anp_6$age_dyad == 'AP' | plot_data_anp_6$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 6 ################
head(summaries)
length(unique(plot_data_anp_6$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids6 <- c(unique(counts_df_period6$id_1), unique(counts_df_period6$id_2)[length(unique(counts_df_period6$id_2))])
males6 <- males[,c(21,6,9,22:53)]
males6 <- males6 %>% dplyr::filter(id %in% ids6)
males6

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

which(males6$degree_0.1 < males6$degree_0.2)
which(males6$degree_0.2 < males6$degree_0.3)
which(males6$degree_0.3 < males6$degree_0.4)
which(males6$degree_0.4 < males6$degree_0.5)

# age variable
males6$age <- lubridate::year(periods$period_start[periods$period == 6]) - males6$byr
summary(males6$age)
males6$age_class <- ifelse(males6$age < 5, 1,
                           ifelse(males6$age < 10, 2,
                                  ifelse(males6$age < 15, 3,
                                         ifelse(males6$age < 20, 4,
                                                ifelse(males6$age < 25, 5,
                                                       ifelse(males6$age < 40, 6, 7))))))

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
     vertex.label = males6$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males6$age_class == 7,'seagreen4',
                          ifelse(males6$age_class == 6,'seagreen3',
                                 ifelse(males6$age_class == 5,'seagreen2',
                                        ifelse(males6$age_class == 4,'steelblue3',
                                               ifelse(males6$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males6$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males6$age_class == 7,'seagreen4',
                          ifelse(males6$age_class == 6,'seagreen3',
                                 ifelse(males6$age_class == 5,'seagreen2',
                                        ifelse(males6$age_class == 4,'steelblue3',
                                               ifelse(males6$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 6 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males6$id[which(males6$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males6$id[which(males6$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males6[which(males6$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males6[which(males6$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males6[which(males6$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males6[which(males6$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males6[which(males6$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 6 ####
summaries$period <- 6
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period6, draws_anp_6, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males6, plot_data_anp_6, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids6, N, plot_cols)

################ 7.7) Run model on real standardised data -- period 7 ################
### create data list
counts_df_period7 <- counts_df_non0[counts_df_non0$period == 7,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period7),          # total number of times one or other of the dyad was observed
  together = counts_df_period7$event_count,    # count number of sightings seen together
  apart    = counts_df_period7$apart,          # count number of sightings seen apart
  period   = counts_df_period7$period)         # which period it's within

### Fit model
weight_anp_2.2_period7 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period7
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#lp__      -83596.45 -83596.30 83.62 83.62 -83735.70 -83462.10 1.00     1468     2450
#weight[1]      0.22      0.20  0.13  0.13      0.05      0.46 1.00     5920     2729
#weight[2]      0.14      0.12  0.09  0.09      0.03      0.32 1.00     5620     2440
#weight[3]      0.20      0.18  0.12  0.12      0.04      0.42 1.00     4932     2129
#weight[4]      0.16      0.14  0.10  0.09      0.03      0.35 1.00     5649     2415
#weight[5]      0.09      0.07  0.06  0.05      0.02      0.19 1.00     3632     1344
#weight[6]      0.10      0.09  0.07  0.06      0.02      0.24 1.00     4763     2172
#weight[7]      0.06      0.05  0.04  0.04      0.01      0.13 1.00     4769     2192
#weight[8]      0.17      0.15  0.10  0.10      0.03      0.36 1.00     4676     2310
#weight[9]      0.17      0.15  0.10  0.09      0.04      0.35 1.00     5076     2559
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period7$output_files()[1])
draws1_anp_7 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period7$output_files()[2])
draws2_anp_7 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period7$output_files()[3])
draws3_anp_7 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period7$output_files()[4])
draws4_anp_7 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_7 <- rbind(draws1_anp_7, draws2_anp_7, draws3_anp_7, draws4_anp_7)

colnames(draws_anp_7)[2:ncol(draws_anp_7)] <- counts_df_period7$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_7), size = 30, replace = F)

### save data 
write_csv(draws_anp_7, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period7_22.04.27.csv')
rm(draws1_anp_7, draws2_anp_7, draws3_anp_7, draws4_anp_7)

### build traceplots  -- period 7 ####
plot(draws_anp_7[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_7[,plot_cols[2]], col = 'tan')
lines(draws_anp_7[,plot_cols[3]], col = 'orange')
lines(draws_anp_7[,plot_cols[4]], col = 'green')
lines(draws_anp_7[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_7[,plot_cols[6]], col = 'blue')
lines(draws_anp_7[,plot_cols[7]], col = 'red')
lines(draws_anp_7[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_7[,plot_cols[9]], col = 'purple')
lines(draws_anp_7[,plot_cols[10]],col = 'magenta')
lines(draws_anp_7[,plot_cols[11]],col = 'black')
lines(draws_anp_7[,plot_cols[12]], col = 'tan')
lines(draws_anp_7[,plot_cols[13]], col = 'orange')
lines(draws_anp_7[,plot_cols[14]], col = 'green')
lines(draws_anp_7[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_7[,plot_cols[16]], col = 'blue')
lines(draws_anp_7[,plot_cols[17]], col = 'red')
lines(draws_anp_7[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_7[,plot_cols[19]], col = 'purple')
lines(draws_anp_7[,plot_cols[20]],col = 'magenta')
lines(draws_anp_7[,plot_cols[21]],col = 'black')
lines(draws_anp_7[,plot_cols[22]], col = 'tan')
lines(draws_anp_7[,plot_cols[23]], col = 'orange')
lines(draws_anp_7[,plot_cols[24]], col = 'green')
lines(draws_anp_7[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_7[,plot_cols[26]], col = 'blue')
lines(draws_anp_7[,plot_cols[27]], col = 'red')
lines(draws_anp_7[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_7[,plot_cols[29]], col = 'purple')
lines(draws_anp_7[,plot_cols[30]],col = 'magenta')

### density plots  -- period 7 ####
dens(draws_anp_7[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_7[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 7 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_7 <- left_join(x = summaries, y = counts_df_period7, by = 'dyad')
head(plot_data_anp_7)
plot_data_anp_7$age_cat_1 <- ifelse(plot_data_anp_7$age_start_1 < 6, 'C',
                                    ifelse(plot_data_anp_7$age_start_1 < 11, 'J',
                                           ifelse(plot_data_anp_7$age_start_1 > 19, 'A','P')))
plot_data_anp_7$age_cat_2 <- ifelse(plot_data_anp_7$age_start_2 < 6, 'C',
                                    ifelse(plot_data_anp_7$age_start_2 < 11, 'J',
                                           ifelse(plot_data_anp_7$age_start_2 > 19, 'A','P')))
plot_data_anp_7$age_catnum_1 <- ifelse(plot_data_anp_7$age_start_1 < 6, 1,
                                       ifelse(plot_data_anp_7$age_start_1 < 11, 2,
                                              ifelse(plot_data_anp_7$age_start_1 < 16, 3,
                                                     ifelse(plot_data_anp_7$age_start_1 < 20, 4,
                                                            ifelse(plot_data_anp_7$age_start_1 < 25, 5,
                                                                   ifelse(plot_data_anp_7$age_start_1 < 40, 6, 7))))))
plot_data_anp_7$age_catnum_2 <- ifelse(plot_data_anp_7$age_start_2 < 6, 1,
                                       ifelse(plot_data_anp_7$age_start_2 < 11, 2,
                                              ifelse(plot_data_anp_7$age_start_2 < 16, 3,
                                                     ifelse(plot_data_anp_7$age_start_2 < 20, 4,
                                                            ifelse(plot_data_anp_7$age_start_2 < 25, 5,
                                                                   ifelse(plot_data_anp_7$age_start_2 < 40, 6, 7))))))

plot_data_anp_7$age_dyad <- ifelse(plot_data_anp_7$age_catnum_1 >= plot_data_anp_7$age_catnum_2,
                                   paste(plot_data_anp_7$age_cat_1, plot_data_anp_7$age_cat_2, sep = ''),
                                   paste(plot_data_anp_7$age_cat_2, plot_data_anp_7$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_7, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_7, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_7[plot_data_anp_7$age_dyad == 'AA' | plot_data_anp_7$age_dyad == 'AP' | plot_data_anp_7$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 7 ################
head(summaries)
length(unique(plot_data_anp_7$id_1))+1 # number of individuals = 156

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids7 <- c(unique(counts_df_period7$id_1), unique(counts_df_period7$id_2)[length(unique(counts_df_period7$id_2))])
males7 <- males[,c(21,6,9,22:53)]
males7 <- males7 %>% dplyr::filter(id %in% ids7)
males7

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

which(males7$degree_0.1 < males7$degree_0.2)
which(males7$degree_0.2 < males7$degree_0.3)
which(males7$degree_0.3 < males7$degree_0.4)
which(males7$degree_0.4 < males7$degree_0.5)

# age variable
males7$age <- lubridate::year(periods$period_start[periods$period == 7]) - males7$byr
summary(males7$age)
males7$age_class <- ifelse(males7$age < 5, 1,
                           ifelse(males7$age < 10, 2,
                                  ifelse(males7$age < 15, 3,
                                         ifelse(males7$age < 20, 4,
                                                ifelse(males7$age < 25, 5,
                                                       ifelse(males7$age < 40, 6, 7))))))

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
     vertex.label = males7$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males7$age_class == 7,'seagreen4',
                          ifelse(males7$age_class == 6,'seagreen3',
                                 ifelse(males7$age_class == 5,'seagreen2',
                                        ifelse(males7$age_class == 4,'steelblue3',
                                               ifelse(males7$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males7$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males7$age_class == 7,'seagreen4',
                          ifelse(males7$age_class == 6,'seagreen3',
                                 ifelse(males7$age_class == 5,'seagreen2',
                                        ifelse(males7$age_class == 4,'steelblue3',
                                               ifelse(males7$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 7 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males7$id[which(males7$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males7$id[which(males7$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males7[which(males7$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males7[which(males7$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males7[which(males7$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males7[which(males7$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males7[which(males7$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 7 ####
summaries$period <- 7
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period7, draws_anp_7, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males7, plot_data_anp_7, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids7, N, plot_cols)

################ 7.8) Run model on real standardised data -- period 8 ################
### create data list
counts_df_period8 <- counts_df_non0[counts_df_non0$period == 8,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period8),          # total number of times one or other of the dyad was observed
  together = counts_df_period8$event_count,    # count number of sightings seen together
  apart    = counts_df_period8$apart,          # count number of sightings seen apart
  period   = counts_df_period8$period)         # which period it's within

### Fit model
weight_anp_2.2_period8 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period8
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period8$output_files()[1])
draws1_anp_8 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period8$output_files()[2])
draws2_anp_8 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period8$output_files()[3])
draws3_anp_8 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period8$output_files()[4])
draws4_anp_8 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_8 <- rbind(draws1_anp_8, draws2_anp_8, draws3_anp_8, draws4_anp_8)

colnames(draws_anp_8)[2:ncol(draws_anp_8)] <- counts_df_period8$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_8), size = 30, replace = F)

### save data 
write_csv(draws_anp_8, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period8_22.04.27.csv')
rm(draws1_anp_8, draws2_anp_8, draws3_anp_8, draws4_anp_8)

### build traceplots  -- period 8 ####
plot(draws_anp_8[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_8[,plot_cols[2]], col = 'tan')
lines(draws_anp_8[,plot_cols[3]], col = 'orange')
lines(draws_anp_8[,plot_cols[4]], col = 'green')
lines(draws_anp_8[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_8[,plot_cols[6]], col = 'blue')
lines(draws_anp_8[,plot_cols[7]], col = 'red')
lines(draws_anp_8[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_8[,plot_cols[9]], col = 'purple')
lines(draws_anp_8[,plot_cols[10]],col = 'magenta')
lines(draws_anp_8[,plot_cols[11]],col = 'black')
lines(draws_anp_8[,plot_cols[12]], col = 'tan')
lines(draws_anp_8[,plot_cols[13]], col = 'orange')
lines(draws_anp_8[,plot_cols[14]], col = 'green')
lines(draws_anp_8[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_8[,plot_cols[16]], col = 'blue')
lines(draws_anp_8[,plot_cols[17]], col = 'red')
lines(draws_anp_8[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_8[,plot_cols[19]], col = 'purple')
lines(draws_anp_8[,plot_cols[20]],col = 'magenta')
lines(draws_anp_8[,plot_cols[21]],col = 'black')
lines(draws_anp_8[,plot_cols[22]], col = 'tan')
lines(draws_anp_8[,plot_cols[23]], col = 'orange')
lines(draws_anp_8[,plot_cols[24]], col = 'green')
lines(draws_anp_8[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_8[,plot_cols[26]], col = 'blue')
lines(draws_anp_8[,plot_cols[27]], col = 'red')
lines(draws_anp_8[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_8[,plot_cols[29]], col = 'purple')
lines(draws_anp_8[,plot_cols[30]],col = 'magenta')

### density plots  -- period 8 ####
dens(draws_anp_8[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_8[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 8 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_8 <- left_join(x = summaries, y = counts_df_period8, by = 'dyad')
head(plot_data_anp_8)
plot_data_anp_8$age_cat_1 <- ifelse(plot_data_anp_8$age_start_1 < 6, 'C',
                                    ifelse(plot_data_anp_8$age_start_1 < 11, 'J',
                                           ifelse(plot_data_anp_8$age_start_1 > 19, 'A','P')))
plot_data_anp_8$age_cat_2 <- ifelse(plot_data_anp_8$age_start_2 < 6, 'C',
                                    ifelse(plot_data_anp_8$age_start_2 < 11, 'J',
                                           ifelse(plot_data_anp_8$age_start_2 > 19, 'A','P')))
plot_data_anp_8$age_catnum_1 <- ifelse(plot_data_anp_8$age_start_1 < 6, 1,
                                       ifelse(plot_data_anp_8$age_start_1 < 11, 2,
                                              ifelse(plot_data_anp_8$age_start_1 < 16, 3,
                                                     ifelse(plot_data_anp_8$age_start_1 < 20, 4,
                                                            ifelse(plot_data_anp_8$age_start_1 < 25, 5,
                                                                   ifelse(plot_data_anp_8$age_start_1 < 40, 6, 7))))))
plot_data_anp_8$age_catnum_2 <- ifelse(plot_data_anp_8$age_start_2 < 6, 1,
                                       ifelse(plot_data_anp_8$age_start_2 < 11, 2,
                                              ifelse(plot_data_anp_8$age_start_2 < 16, 3,
                                                     ifelse(plot_data_anp_8$age_start_2 < 20, 4,
                                                            ifelse(plot_data_anp_8$age_start_2 < 25, 5,
                                                                   ifelse(plot_data_anp_8$age_start_2 < 40, 6, 7))))))

plot_data_anp_8$age_dyad <- ifelse(plot_data_anp_8$age_catnum_1 >= plot_data_anp_8$age_catnum_2,
                                   paste(plot_data_anp_8$age_cat_1, plot_data_anp_8$age_cat_2, sep = ''),
                                   paste(plot_data_anp_8$age_cat_2, plot_data_anp_8$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_8, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_8, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_8[plot_data_anp_8$age_dyad == 'AA' | plot_data_anp_8$age_dyad == 'AP' | plot_data_anp_8$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 8 ################
head(summaries)
length(unique(plot_data_anp_8$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids8 <- c(unique(counts_df_period8$id_1), unique(counts_df_period8$id_2)[length(unique(counts_df_period8$id_2))])
males8 <- males[,c(21,6,9,22:53)]
males8 <- males8 %>% dplyr::filter(id %in% ids8)
males8

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

which(males8$degree_0.1 < males8$degree_0.2)
which(males8$degree_0.2 < males8$degree_0.3)
which(males8$degree_0.3 < males8$degree_0.4)
which(males8$degree_0.4 < males8$degree_0.5)

# age variable
males8$age <- lubridate::year(periods$period_start[periods$period == 8]) - males8$byr
summary(males8$age)
males8$age_class <- ifelse(males8$age < 5, 1,
                           ifelse(males8$age < 10, 2,
                                  ifelse(males8$age < 15, 3,
                                         ifelse(males8$age < 20, 4,
                                                ifelse(males8$age < 25, 5,
                                                       ifelse(males8$age < 40, 6, 7))))))

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
     vertex.label = males8$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males8$age_class == 7,'seagreen4',
                          ifelse(males8$age_class == 6,'seagreen3',
                                 ifelse(males8$age_class == 5,'seagreen2',
                                        ifelse(males8$age_class == 4,'steelblue3',
                                               ifelse(males8$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males8$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males8$age_class == 7,'seagreen4',
                          ifelse(males8$age_class == 6,'seagreen3',
                                 ifelse(males8$age_class == 5,'seagreen2',
                                        ifelse(males8$age_class == 4,'steelblue3',
                                               ifelse(males8$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.2 -- period 8 ####
g_mid_0.2 <- delete.vertices(graph = g_mid, v = males8$id[which(males8$degree_0.2 == 0)])
g_rng_0.2 <- delete.vertices(graph = g_rng, v = males8$id[which(males8$degree_0.2 == 0)])

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
     vertex.color= ifelse(males8[which(males8$degree_0.2 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males8[which(males8$degree_0.2 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males8[which(males8$degree_0.2 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males8[which(males8$degree_0.2 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males8[which(males8$degree_0.2 != 0),]$age_class == 3,'steelblue1',
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

### save summary data -- period 8 ####
summaries$period <- 8
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period8, draws_anp_8, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males8, plot_data_anp_8, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids8, N, plot_cols)

################ 7.9) Run model on real standardised data -- period 9 ################
### create data list
counts_df_period9 <- counts_df_non0[counts_df_non0$period == 9,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period9),          # total number of times one or other of the dyad was observed
  together = counts_df_period9$event_count,    # count number of sightings seen together
  apart    = counts_df_period9$apart,          # count number of sightings seen apart
  period   = counts_df_period9$period)         # which period it's within

### Fit model
weight_anp_2.2_period9 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period9
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period9$output_files()[1])
draws1_anp_9 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period9$output_files()[2])
draws2_anp_9 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period9$output_files()[3])
draws3_anp_9 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period9$output_files()[4])
draws4_anp_9 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_9 <- rbind(draws1_anp_9, draws2_anp_9, draws3_anp_9, draws4_anp_9)

colnames(draws_anp_9)[2:ncol(draws_anp_9)] <- counts_df_period9$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_9), size = 30, replace = F)

### save data 
write_csv(draws_anp_9, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period9_22.04.27.csv')
rm(draws1_anp_9, draws2_anp_9, draws3_anp_9, draws4_anp_9)

### build traceplots  -- period 9 ####
plot(draws_anp_9[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_9[,plot_cols[2]], col = 'tan')
lines(draws_anp_9[,plot_cols[3]], col = 'orange')
lines(draws_anp_9[,plot_cols[4]], col = 'green')
lines(draws_anp_9[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_9[,plot_cols[6]], col = 'blue')
lines(draws_anp_9[,plot_cols[7]], col = 'red')
lines(draws_anp_9[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_9[,plot_cols[9]], col = 'purple')
lines(draws_anp_9[,plot_cols[10]],col = 'magenta')
lines(draws_anp_9[,plot_cols[11]],col = 'black')
lines(draws_anp_9[,plot_cols[12]], col = 'tan')
lines(draws_anp_9[,plot_cols[13]], col = 'orange')
lines(draws_anp_9[,plot_cols[14]], col = 'green')
lines(draws_anp_9[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_9[,plot_cols[16]], col = 'blue')
lines(draws_anp_9[,plot_cols[17]], col = 'red')
lines(draws_anp_9[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_9[,plot_cols[19]], col = 'purple')
lines(draws_anp_9[,plot_cols[20]],col = 'magenta')
lines(draws_anp_9[,plot_cols[21]],col = 'black')
lines(draws_anp_9[,plot_cols[22]], col = 'tan')
lines(draws_anp_9[,plot_cols[23]], col = 'orange')
lines(draws_anp_9[,plot_cols[24]], col = 'green')
lines(draws_anp_9[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_9[,plot_cols[26]], col = 'blue')
lines(draws_anp_9[,plot_cols[27]], col = 'red')
lines(draws_anp_9[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_9[,plot_cols[29]], col = 'purple')
lines(draws_anp_9[,plot_cols[30]],col = 'magenta')

### density plots  -- period 9 ####
dens(draws_anp_9[,2], ylim = c(0,40),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_9[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 9 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_9 <- left_join(x = summaries, y = counts_df_period9, by = 'dyad')
head(plot_data_anp_9)
plot_data_anp_9$age_cat_1 <- ifelse(plot_data_anp_9$age_start_1 < 6, 'C',
                                    ifelse(plot_data_anp_9$age_start_1 < 11, 'J',
                                           ifelse(plot_data_anp_9$age_start_1 > 19, 'A','P')))
plot_data_anp_9$age_cat_2 <- ifelse(plot_data_anp_9$age_start_2 < 6, 'C',
                                    ifelse(plot_data_anp_9$age_start_2 < 11, 'J',
                                           ifelse(plot_data_anp_9$age_start_2 > 19, 'A','P')))
plot_data_anp_9$age_catnum_1 <- ifelse(plot_data_anp_9$age_start_1 < 6, 1,
                                       ifelse(plot_data_anp_9$age_start_1 < 11, 2,
                                              ifelse(plot_data_anp_9$age_start_1 < 16, 3,
                                                     ifelse(plot_data_anp_9$age_start_1 < 20, 4,
                                                            ifelse(plot_data_anp_9$age_start_1 < 25, 5,
                                                                   ifelse(plot_data_anp_9$age_start_1 < 40, 6, 7))))))
plot_data_anp_9$age_catnum_2 <- ifelse(plot_data_anp_9$age_start_2 < 6, 1,
                                       ifelse(plot_data_anp_9$age_start_2 < 11, 2,
                                              ifelse(plot_data_anp_9$age_start_2 < 16, 3,
                                                     ifelse(plot_data_anp_9$age_start_2 < 20, 4,
                                                            ifelse(plot_data_anp_9$age_start_2 < 25, 5,
                                                                   ifelse(plot_data_anp_9$age_start_2 < 40, 6, 7))))))

plot_data_anp_9$age_dyad <- ifelse(plot_data_anp_9$age_catnum_1 >= plot_data_anp_9$age_catnum_2,
                                   paste(plot_data_anp_9$age_cat_1, plot_data_anp_9$age_cat_2, sep = ''),
                                   paste(plot_data_anp_9$age_cat_2, plot_data_anp_9$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_9, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_9, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_9[plot_data_anp_9$age_dyad == 'AA' | plot_data_anp_9$age_dyad == 'AP' | plot_data_anp_9$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 9 ################
head(summaries)
length(unique(plot_data_anp_9$id_1))+1 # number of individuals = 166

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids9 <- c(unique(counts_df_period9$id_1), unique(counts_df_period9$id_2)[length(unique(counts_df_period9$id_2))])
males9 <- males[,c(21,6,9,22:53)]
males9 <- males9 %>% dplyr::filter(id %in% ids9)
males9

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

which(males9$degree_0.1 < males9$degree_0.2)
which(males9$degree_0.2 < males9$degree_0.3)
which(males9$degree_0.3 < males9$degree_0.4)
which(males9$degree_0.4 < males9$degree_0.5)

# age variable
males9$age <- lubridate::year(periods$period_start[periods$period == 9]) - males9$byr
summary(males9$age)
males9$age_class <- ifelse(males9$age < 5, 1,
                           ifelse(males9$age < 10, 2,
                                  ifelse(males9$age < 15, 3,
                                         ifelse(males9$age < 20, 4,
                                                ifelse(males9$age < 25, 5,
                                                       ifelse(males9$age < 40, 6, 7))))))

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
     vertex.label = males9$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males9$age_class == 7,'seagreen4',
                          ifelse(males9$age_class == 6,'seagreen3',
                                 ifelse(males9$age_class == 5,'seagreen2',
                                        ifelse(males9$age_class == 4,'steelblue3',
                                               ifelse(males9$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males9$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males9$age_class == 7,'seagreen4',
                          ifelse(males9$age_class == 6,'seagreen3',
                                 ifelse(males9$age_class == 5,'seagreen2',
                                        ifelse(males9$age_class == 4,'steelblue3',
                                               ifelse(males9$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.2 -- period 9 ####
g_mid_0.2 <- delete.vertices(graph = g_mid, v = males9$id[which(males9$degree_0.2 == 0)])
g_rng_0.2 <- delete.vertices(graph = g_rng, v = males9$id[which(males9$degree_0.2 == 0)])

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
     vertex.color= ifelse(males9[which(males9$degree_0.2 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males9[which(males9$degree_0.2 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males9[which(males9$degree_0.2 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males9[which(males9$degree_0.2 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males9[which(males9$degree_0.2 != 0),]$age_class == 3,'steelblue1',
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

### save summary data -- period 9 ####
summaries$period <- 9
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period9, draws_anp_9, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males9, plot_data_anp_9, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids9, N, plot_cols)

################ 7.10) Run model on real standardised data -- period 10 ################
### create data list
counts_df_period10 <- counts_df_non0[counts_df_non0$period == 10,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period10),          # total number of times one or other of the dyad was observed
  together = counts_df_period10$event_count,    # count number of sightings seen together
  apart    = counts_df_period10$apart,          # count number of sightings seen apart
  period   = counts_df_period10$period)         # which period it's within

### Fit model
weight_anp_2.2_period10 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period10
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period10$output_files()[1])
draws1_anp_10 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period10$output_files()[2])
draws2_anp_10 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period10$output_files()[3])
draws3_anp_10 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period10$output_files()[4])
draws4_anp_10 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_10 <- rbind(draws1_anp_10, draws2_anp_10, draws3_anp_10, draws4_anp_10)

colnames(draws_anp_10)[2:ncol(draws_anp_10)] <- counts_df_period10$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_10), size = 30, replace = F)

### save data 
write_csv(draws_anp_10, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period10_22.04.27.csv')
rm(draws1_anp_10, draws2_anp_10, draws3_anp_10, draws4_anp_10)

### build traceplots  -- period 10 ####
plot(draws_anp_10[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_10[,plot_cols[2]], col = 'tan')
lines(draws_anp_10[,plot_cols[3]], col = 'orange')
lines(draws_anp_10[,plot_cols[4]], col = 'green')
lines(draws_anp_10[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_10[,plot_cols[6]], col = 'blue')
lines(draws_anp_10[,plot_cols[7]], col = 'red')
lines(draws_anp_10[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_10[,plot_cols[9]], col = 'purple')
lines(draws_anp_10[,plot_cols[10]],col = 'magenta')
lines(draws_anp_10[,plot_cols[11]],col = 'black')
lines(draws_anp_10[,plot_cols[12]], col = 'tan')
lines(draws_anp_10[,plot_cols[13]], col = 'orange')
lines(draws_anp_10[,plot_cols[14]], col = 'green')
lines(draws_anp_10[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_10[,plot_cols[16]], col = 'blue')
lines(draws_anp_10[,plot_cols[17]], col = 'red')
lines(draws_anp_10[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_10[,plot_cols[19]], col = 'purple')
lines(draws_anp_10[,plot_cols[20]],col = 'magenta')
lines(draws_anp_10[,plot_cols[21]],col = 'black')
lines(draws_anp_10[,plot_cols[22]], col = 'tan')
lines(draws_anp_10[,plot_cols[23]], col = 'orange')
lines(draws_anp_10[,plot_cols[24]], col = 'green')
lines(draws_anp_10[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_10[,plot_cols[26]], col = 'blue')
lines(draws_anp_10[,plot_cols[27]], col = 'red')
lines(draws_anp_10[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_10[,plot_cols[29]], col = 'purple')
lines(draws_anp_10[,plot_cols[30]],col = 'magenta')

### density plots  -- period 10 ####
dens(draws_anp_10[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_10[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 10 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_10 <- left_join(x = summaries, y = counts_df_period10, by = 'dyad')
head(plot_data_anp_10)
plot_data_anp_10$age_cat_1 <- ifelse(plot_data_anp_10$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_10$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_10$age_start_1 > 19, 'A','P')))
plot_data_anp_10$age_cat_2 <- ifelse(plot_data_anp_10$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_10$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_10$age_start_2 > 19, 'A','P')))
plot_data_anp_10$age_catnum_1 <- ifelse(plot_data_anp_10$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_10$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_10$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_10$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_10$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_10$age_start_1 < 40, 6, 7))))))
plot_data_anp_10$age_catnum_2 <- ifelse(plot_data_anp_10$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_10$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_10$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_10$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_10$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_10$age_start_2 < 40, 6, 7))))))

plot_data_anp_10$age_dyad <- ifelse(plot_data_anp_10$age_catnum_1 >= plot_data_anp_10$age_catnum_2,
                                    paste(plot_data_anp_10$age_cat_1, plot_data_anp_10$age_cat_2, sep = ''),
                                    paste(plot_data_anp_10$age_cat_2, plot_data_anp_10$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_10, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_10, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_10[plot_data_anp_10$age_dyad == 'AA' | plot_data_anp_10$age_dyad == 'AP' | plot_data_anp_10$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 10 ################
head(summaries)
length(unique(plot_data_anp_10$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids10 <- c(unique(counts_df_period10$id_1), unique(counts_df_period10$id_2)[length(unique(counts_df_period10$id_2))])
males10 <- males[,c(21,6,9,22:53)]
males10 <- males10 %>% dplyr::filter(id %in% ids10)
males10

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

which(males10$degree_0.1 < males10$degree_0.2)
which(males10$degree_0.2 < males10$degree_0.3)
which(males10$degree_0.3 < males10$degree_0.4)
which(males10$degree_0.4 < males10$degree_0.5)

# age variable
males10$age <- lubridate::year(periods$period_start[periods$period == 10]) - males10$byr
summary(males10$age)
males10$age_class <- ifelse(males10$age < 5, 1,
                            ifelse(males10$age < 10, 2,
                                   ifelse(males10$age < 15, 3,
                                          ifelse(males10$age < 20, 4,
                                                 ifelse(males10$age < 25, 5,
                                                        ifelse(males10$age < 40, 6, 7))))))

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
     vertex.label = males10$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males10$age_class == 7,'seagreen4',
                          ifelse(males10$age_class == 6,'seagreen3',
                                 ifelse(males10$age_class == 5,'seagreen2',
                                        ifelse(males10$age_class == 4,'steelblue3',
                                               ifelse(males10$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males10$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males10$age_class == 7,'seagreen4',
                          ifelse(males10$age_class == 6,'seagreen3',
                                 ifelse(males10$age_class == 5,'seagreen2',
                                        ifelse(males10$age_class == 4,'steelblue3',
                                               ifelse(males10$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 10 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males10$id[which(males10$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males10$id[which(males10$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males10[which(males10$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males10[which(males10$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males10[which(males10$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males10[which(males10$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males10[which(males10$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 10 ####
summaries$period <- 10
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period10, draws_anp_10, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males10, plot_data_anp_10, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids10, N, plot_cols)

#### Save outputs to file #####
# write out csv file
write_csv(dyad_period_weights, 'data_processed/anp_dyad_weightdistributions_2.2_periods1to10_22.05.06.csv')

# produce PDF of graphs
dev.off()

# add progress stamp
print(paste0('Time windows 6:10 completed at ', Sys.time()))

################ Time windows 11:15 -- not yet tested individually ################
# set seed
set.seed(12345)

# create file of output graphs
pdf('data_processed/anp_edgeweights_2.2_period11to15_22.05.06.pdf', width = 10, height = 10)

# read in previous summary file
#dyad_period_weights <- read_csv('data_processed/anp_dyad_weightdistributions_2.2_periods6to10_22.05.06.csv')

################ 7.11) Run model on real standardised data -- period 11 ################
### create data list
counts_df_period11 <- counts_df_non0[counts_df_non0$period == 11,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period11),          # total number of times one or other of the dyad was observed
  together = counts_df_period11$event_count,    # count number of sightings seen together
  apart    = counts_df_period11$apart,          # count number of sightings seen apart
  period   = counts_df_period11$period)         # which period it's within

### Fit model
weight_anp_2.2_period11 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period11
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period11$output_files()[1])
draws1_anp_11 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period11$output_files()[2])
draws2_anp_11 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period11$output_files()[3])
draws3_anp_11 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period11$output_files()[4])
draws4_anp_11 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_11 <- rbind(draws1_anp_11, draws2_anp_11, draws3_anp_11, draws4_anp_11)

colnames(draws_anp_11)[2:ncol(draws_anp_11)] <- counts_df_period11$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_11), size = 30, replace = F)

### save data 
write_csv(draws_anp_11, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period11_22.04.27.csv')
rm(draws1_anp_11, draws2_anp_11, draws3_anp_11, draws4_anp_11)

### build traceplots  -- period 11 ####
plot(draws_anp_11[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_11[,plot_cols[2]], col = 'tan')
lines(draws_anp_11[,plot_cols[3]], col = 'orange')
lines(draws_anp_11[,plot_cols[4]], col = 'green')
lines(draws_anp_11[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_11[,plot_cols[6]], col = 'blue')
lines(draws_anp_11[,plot_cols[7]], col = 'red')
lines(draws_anp_11[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_11[,plot_cols[9]], col = 'purple')
lines(draws_anp_11[,plot_cols[10]],col = 'magenta')
lines(draws_anp_11[,plot_cols[11]],col = 'black')
lines(draws_anp_11[,plot_cols[12]], col = 'tan')
lines(draws_anp_11[,plot_cols[13]], col = 'orange')
lines(draws_anp_11[,plot_cols[14]], col = 'green')
lines(draws_anp_11[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_11[,plot_cols[16]], col = 'blue')
lines(draws_anp_11[,plot_cols[17]], col = 'red')
lines(draws_anp_11[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_11[,plot_cols[19]], col = 'purple')
lines(draws_anp_11[,plot_cols[20]],col = 'magenta')
lines(draws_anp_11[,plot_cols[21]],col = 'black')
lines(draws_anp_11[,plot_cols[22]], col = 'tan')
lines(draws_anp_11[,plot_cols[23]], col = 'orange')
lines(draws_anp_11[,plot_cols[24]], col = 'green')
lines(draws_anp_11[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_11[,plot_cols[26]], col = 'blue')
lines(draws_anp_11[,plot_cols[27]], col = 'red')
lines(draws_anp_11[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_11[,plot_cols[29]], col = 'purple')
lines(draws_anp_11[,plot_cols[30]],col = 'magenta')

### density plots  -- period 11 ####
dens(draws_anp_11[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_11[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 11 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_11 <- left_join(x = summaries, y = counts_df_period11, by = 'dyad')
head(plot_data_anp_11)
plot_data_anp_11$age_cat_1 <- ifelse(plot_data_anp_11$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_11$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_11$age_start_1 > 19, 'A','P')))
plot_data_anp_11$age_cat_2 <- ifelse(plot_data_anp_11$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_11$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_11$age_start_2 > 19, 'A','P')))
plot_data_anp_11$age_catnum_1 <- ifelse(plot_data_anp_11$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_11$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_11$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_11$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_11$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_11$age_start_1 < 40, 6, 7))))))
plot_data_anp_11$age_catnum_2 <- ifelse(plot_data_anp_11$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_11$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_11$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_11$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_11$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_11$age_start_2 < 40, 6, 7))))))

plot_data_anp_11$age_dyad <- ifelse(plot_data_anp_11$age_catnum_1 >= plot_data_anp_11$age_catnum_2,
                                    paste(plot_data_anp_11$age_cat_1, plot_data_anp_11$age_cat_2, sep = ''),
                                    paste(plot_data_anp_11$age_cat_2, plot_data_anp_11$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_11, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_11, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_11[plot_data_anp_11$age_dyad == 'AA' | plot_data_anp_11$age_dyad == 'AP' | plot_data_anp_11$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 11 ################
head(summaries)
length(unique(plot_data_anp_11$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids11 <- c(unique(counts_df_period11$id_1), unique(counts_df_period11$id_2)[length(unique(counts_df_period11$id_2))])
males11 <- males[,c(21,6,9,22:53)]
males11 <- males11 %>% dplyr::filter(id %in% ids11)
males11

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

which(males11$degree_0.1 < males11$degree_0.2)
which(males11$degree_0.2 < males11$degree_0.3)
which(males11$degree_0.3 < males11$degree_0.4)
which(males11$degree_0.4 < males11$degree_0.5)

# age variable
males11$age <- lubridate::year(periods$period_start[periods$period == 11]) - males11$byr
summary(males11$age)
males11$age_class <- ifelse(males11$age < 5, 1,
                            ifelse(males11$age < 10, 2,
                                   ifelse(males11$age < 15, 3,
                                          ifelse(males11$age < 20, 4,
                                                 ifelse(males11$age < 25, 5,
                                                        ifelse(males11$age < 40, 6, 7))))))

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
     vertex.label = males11$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males11$age_class == 7,'seagreen4',
                          ifelse(males11$age_class == 6,'seagreen3',
                                 ifelse(males11$age_class == 5,'seagreen2',
                                        ifelse(males11$age_class == 4,'steelblue3',
                                               ifelse(males11$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males11$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males11$age_class == 7,'seagreen4',
                          ifelse(males11$age_class == 6,'seagreen3',
                                 ifelse(males11$age_class == 5,'seagreen2',
                                        ifelse(males11$age_class == 4,'steelblue3',
                                               ifelse(males11$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 11 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males11$id[which(males11$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males11$id[which(males11$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males11[which(males11$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males11[which(males11$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males11[which(males11$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males11[which(males11$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males11[which(males11$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 11 ####
summaries$period <- 11
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period11, draws_anp_11, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males11, plot_data_anp_11, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids11, N, plot_cols)

################ 7.12) Run model on real standardised data -- period 12 ################
### create data list
counts_df_period12 <- counts_df_non0[counts_df_non0$period == 12,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period12),          # total number of times one or other of the dyad was observed
  together = counts_df_period12$event_count,    # count number of sightings seen together
  apart    = counts_df_period12$apart,          # count number of sightings seen apart
  period   = counts_df_period12$period)         # which period it's within

### Fit model
weight_anp_2.2_period12 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period12
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period12$output_files()[1])
draws1_anp_12 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period12$output_files()[2])
draws2_anp_12 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period12$output_files()[3])
draws3_anp_12 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period12$output_files()[4])
draws4_anp_12 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_12 <- rbind(draws1_anp_12, draws2_anp_12, draws3_anp_12, draws4_anp_12)

colnames(draws_anp_12)[2:ncol(draws_anp_12)] <- counts_df_period12$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_12), size = 30, replace = F)

### save data 
write_csv(draws_anp_12, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period12_22.04.27.csv')
rm(draws1_anp_12, draws2_anp_12, draws3_anp_12, draws4_anp_12)

### build traceplots  -- period 12 ####
plot(draws_anp_12[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_12[,plot_cols[2]], col = 'tan')
lines(draws_anp_12[,plot_cols[3]], col = 'orange')
lines(draws_anp_12[,plot_cols[4]], col = 'green')
lines(draws_anp_12[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_12[,plot_cols[6]], col = 'blue')
lines(draws_anp_12[,plot_cols[7]], col = 'red')
lines(draws_anp_12[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_12[,plot_cols[9]], col = 'purple')
lines(draws_anp_12[,plot_cols[10]],col = 'magenta')
lines(draws_anp_12[,plot_cols[11]],col = 'black')
lines(draws_anp_12[,plot_cols[12]], col = 'tan')
lines(draws_anp_12[,plot_cols[13]], col = 'orange')
lines(draws_anp_12[,plot_cols[14]], col = 'green')
lines(draws_anp_12[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_12[,plot_cols[16]], col = 'blue')
lines(draws_anp_12[,plot_cols[17]], col = 'red')
lines(draws_anp_12[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_12[,plot_cols[19]], col = 'purple')
lines(draws_anp_12[,plot_cols[20]],col = 'magenta')
lines(draws_anp_12[,plot_cols[21]],col = 'black')
lines(draws_anp_12[,plot_cols[22]], col = 'tan')
lines(draws_anp_12[,plot_cols[23]], col = 'orange')
lines(draws_anp_12[,plot_cols[24]], col = 'green')
lines(draws_anp_12[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_12[,plot_cols[26]], col = 'blue')
lines(draws_anp_12[,plot_cols[27]], col = 'red')
lines(draws_anp_12[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_12[,plot_cols[29]], col = 'purple')
lines(draws_anp_12[,plot_cols[30]],col = 'magenta')

### density plots  -- period 12 ####
dens(draws_anp_12[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_12[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 12 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_12 <- left_join(x = summaries, y = counts_df_period12, by = 'dyad')
head(plot_data_anp_12)
plot_data_anp_12$age_cat_1 <- ifelse(plot_data_anp_12$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_12$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_12$age_start_1 > 19, 'A','P')))
plot_data_anp_12$age_cat_2 <- ifelse(plot_data_anp_12$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_12$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_12$age_start_2 > 19, 'A','P')))
plot_data_anp_12$age_catnum_1 <- ifelse(plot_data_anp_12$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_12$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_12$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_12$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_12$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_12$age_start_1 < 40, 6, 7))))))
plot_data_anp_12$age_catnum_2 <- ifelse(plot_data_anp_12$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_12$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_12$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_12$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_12$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_12$age_start_2 < 40, 6, 7))))))

plot_data_anp_12$age_dyad <- ifelse(plot_data_anp_12$age_catnum_1 >= plot_data_anp_12$age_catnum_2,
                                    paste(plot_data_anp_12$age_cat_1, plot_data_anp_12$age_cat_2, sep = ''),
                                    paste(plot_data_anp_12$age_cat_2, plot_data_anp_12$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_12, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_12, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_12[plot_data_anp_12$age_dyad == 'AA' | plot_data_anp_12$age_dyad == 'AP' | plot_data_anp_12$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 12 ################
head(summaries)
length(unique(plot_data_anp_12$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids12 <- c(unique(counts_df_period12$id_1), unique(counts_df_period12$id_2)[length(unique(counts_df_period12$id_2))])
males12 <- males[,c(21,6,9,22:53)]
males12 <- males12 %>% dplyr::filter(id %in% ids12)
males12

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

which(males12$degree_0.1 < males12$degree_0.2)
which(males12$degree_0.2 < males12$degree_0.3)
which(males12$degree_0.3 < males12$degree_0.4)
which(males12$degree_0.4 < males12$degree_0.5)

# age variable
males12$age <- lubridate::year(periods$period_start[periods$period == 12]) - males12$byr
summary(males12$age)
males12$age_class <- ifelse(males12$age < 5, 1,
                            ifelse(males12$age < 10, 2,
                                   ifelse(males12$age < 15, 3,
                                          ifelse(males12$age < 20, 4,
                                                 ifelse(males12$age < 25, 5,
                                                        ifelse(males12$age < 40, 6, 7))))))

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
     vertex.label = males12$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males12$age_class == 7,'seagreen4',
                          ifelse(males12$age_class == 6,'seagreen3',
                                 ifelse(males12$age_class == 5,'seagreen2',
                                        ifelse(males12$age_class == 4,'steelblue3',
                                               ifelse(males12$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males12$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males12$age_class == 7,'seagreen4',
                          ifelse(males12$age_class == 6,'seagreen3',
                                 ifelse(males12$age_class == 5,'seagreen2',
                                        ifelse(males12$age_class == 4,'steelblue3',
                                               ifelse(males12$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 12 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males12$id[which(males12$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males12$id[which(males12$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males12[which(males12$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males12[which(males12$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males12[which(males12$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males12[which(males12$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males12[which(males12$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 12 ####
summaries$period <- 12
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period12, draws_anp_12, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males12, plot_data_anp_12, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids12, N, plot_cols)

################ 7.13) Run model on real standardised data -- period 13 ################
### create data list
counts_df_period13 <- counts_df_non0[counts_df_non0$period == 13,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period13),          # total number of times one or other of the dyad was observed
  together = counts_df_period13$event_count,    # count number of sightings seen together
  apart    = counts_df_period13$apart,          # count number of sightings seen apart
  period   = counts_df_period13$period)         # which period it's within

### Fit model
weight_anp_2.2_period13 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period13
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period13$output_files()[1])
draws1_anp_13 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period13$output_files()[2])
draws2_anp_13 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period13$output_files()[3])
draws3_anp_13 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period13$output_files()[4])
draws4_anp_13 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_13 <- rbind(draws1_anp_13, draws2_anp_13, draws3_anp_13, draws4_anp_13)

colnames(draws_anp_13)[2:ncol(draws_anp_13)] <- counts_df_period13$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_13), size = 30, replace = F)

### save data 
write_csv(draws_anp_13, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period13_22.04.27.csv')
rm(draws1_anp_13, draws2_anp_13, draws3_anp_13, draws4_anp_13)

### build traceplots  -- period 13 ####
plot(draws_anp_13[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_13[,plot_cols[2]], col = 'tan')
lines(draws_anp_13[,plot_cols[3]], col = 'orange')
lines(draws_anp_13[,plot_cols[4]], col = 'green')
lines(draws_anp_13[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_13[,plot_cols[6]], col = 'blue')
lines(draws_anp_13[,plot_cols[7]], col = 'red')
lines(draws_anp_13[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_13[,plot_cols[9]], col = 'purple')
lines(draws_anp_13[,plot_cols[10]],col = 'magenta')
lines(draws_anp_13[,plot_cols[11]],col = 'black')
lines(draws_anp_13[,plot_cols[12]], col = 'tan')
lines(draws_anp_13[,plot_cols[13]], col = 'orange')
lines(draws_anp_13[,plot_cols[14]], col = 'green')
lines(draws_anp_13[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_13[,plot_cols[16]], col = 'blue')
lines(draws_anp_13[,plot_cols[17]], col = 'red')
lines(draws_anp_13[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_13[,plot_cols[19]], col = 'purple')
lines(draws_anp_13[,plot_cols[20]],col = 'magenta')
lines(draws_anp_13[,plot_cols[21]],col = 'black')
lines(draws_anp_13[,plot_cols[22]], col = 'tan')
lines(draws_anp_13[,plot_cols[23]], col = 'orange')
lines(draws_anp_13[,plot_cols[24]], col = 'green')
lines(draws_anp_13[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_13[,plot_cols[26]], col = 'blue')
lines(draws_anp_13[,plot_cols[27]], col = 'red')
lines(draws_anp_13[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_13[,plot_cols[29]], col = 'purple')
lines(draws_anp_13[,plot_cols[30]],col = 'magenta')

### density plots  -- period 13 ####
dens(draws_anp_13[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_13[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 13 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_13 <- left_join(x = summaries, y = counts_df_period13, by = 'dyad')
head(plot_data_anp_13)
plot_data_anp_13$age_cat_1 <- ifelse(plot_data_anp_13$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_13$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_13$age_start_1 > 19, 'A','P')))
plot_data_anp_13$age_cat_2 <- ifelse(plot_data_anp_13$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_13$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_13$age_start_2 > 19, 'A','P')))
plot_data_anp_13$age_catnum_1 <- ifelse(plot_data_anp_13$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_13$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_13$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_13$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_13$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_13$age_start_1 < 40, 6, 7))))))
plot_data_anp_13$age_catnum_2 <- ifelse(plot_data_anp_13$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_13$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_13$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_13$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_13$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_13$age_start_2 < 40, 6, 7))))))

plot_data_anp_13$age_dyad <- ifelse(plot_data_anp_13$age_catnum_1 >= plot_data_anp_13$age_catnum_2,
                                    paste(plot_data_anp_13$age_cat_1, plot_data_anp_13$age_cat_2, sep = ''),
                                    paste(plot_data_anp_13$age_cat_2, plot_data_anp_13$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_13, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_13, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_13[plot_data_anp_13$age_dyad == 'AA' | plot_data_anp_13$age_dyad == 'AP' | plot_data_anp_13$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 13 ################
head(summaries)
length(unique(plot_data_anp_13$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids13 <- c(unique(counts_df_period13$id_1), unique(counts_df_period13$id_2)[length(unique(counts_df_period13$id_2))])
males13 <- males[,c(21,6,9,22:53)]
males13 <- males13 %>% dplyr::filter(id %in% ids13)
males13

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

which(males13$degree_0.1 < males13$degree_0.2)
which(males13$degree_0.2 < males13$degree_0.3)
which(males13$degree_0.3 < males13$degree_0.4)
which(males13$degree_0.4 < males13$degree_0.5)

# age variable
males13$age <- lubridate::year(periods$period_start[periods$period == 13]) - males13$byr
summary(males13$age)
males13$age_class <- ifelse(males13$age < 5, 1,
                            ifelse(males13$age < 10, 2,
                                   ifelse(males13$age < 15, 3,
                                          ifelse(males13$age < 20, 4,
                                                 ifelse(males13$age < 25, 5,
                                                        ifelse(males13$age < 40, 6, 7))))))

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
     vertex.label = males13$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males13$age_class == 7,'seagreen4',
                          ifelse(males13$age_class == 6,'seagreen3',
                                 ifelse(males13$age_class == 5,'seagreen2',
                                        ifelse(males13$age_class == 4,'steelblue3',
                                               ifelse(males13$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males13$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males13$age_class == 7,'seagreen4',
                          ifelse(males13$age_class == 6,'seagreen3',
                                 ifelse(males13$age_class == 5,'seagreen2',
                                        ifelse(males13$age_class == 4,'steelblue3',
                                               ifelse(males13$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 13 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males13$id[which(males13$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males13$id[which(males13$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males13[which(males13$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males13[which(males13$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males13[which(males13$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males13[which(males13$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males13[which(males13$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 13 ####
summaries$period <- 13
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period13, draws_anp_13, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males13, plot_data_anp_13, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids13, N, plot_cols)

################ 7.14) Run model on real standardised data -- period 14 ################
### create data list
counts_df_period14 <- counts_df_non0[counts_df_non0$period == 14,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period14),          # total number of times one or other of the dyad was observed
  together = counts_df_period14$event_count,    # count number of sightings seen together
  apart    = counts_df_period14$apart,          # count number of sightings seen apart
  period   = counts_df_period14$period)         # which period it's within

### Fit model
weight_anp_2.2_period14 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period14
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period14$output_files()[1])
draws1_anp_14 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period14$output_files()[2])
draws2_anp_14 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period14$output_files()[3])
draws3_anp_14 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period14$output_files()[4])
draws4_anp_14 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_14 <- rbind(draws1_anp_14, draws2_anp_14, draws3_anp_14, draws4_anp_14)

colnames(draws_anp_14)[2:ncol(draws_anp_14)] <- counts_df_period14$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_14), size = 30, replace = F)

### save data 
write_csv(draws_anp_14, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period14_22.04.27.csv')
rm(draws1_anp_14, draws2_anp_14, draws3_anp_14, draws4_anp_14)

### build traceplots  -- period 14 ####
plot(draws_anp_14[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_14[,plot_cols[2]], col = 'tan')
lines(draws_anp_14[,plot_cols[3]], col = 'orange')
lines(draws_anp_14[,plot_cols[4]], col = 'green')
lines(draws_anp_14[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_14[,plot_cols[6]], col = 'blue')
lines(draws_anp_14[,plot_cols[7]], col = 'red')
lines(draws_anp_14[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_14[,plot_cols[9]], col = 'purple')
lines(draws_anp_14[,plot_cols[10]],col = 'magenta')
lines(draws_anp_14[,plot_cols[11]],col = 'black')
lines(draws_anp_14[,plot_cols[12]], col = 'tan')
lines(draws_anp_14[,plot_cols[13]], col = 'orange')
lines(draws_anp_14[,plot_cols[14]], col = 'green')
lines(draws_anp_14[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_14[,plot_cols[16]], col = 'blue')
lines(draws_anp_14[,plot_cols[17]], col = 'red')
lines(draws_anp_14[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_14[,plot_cols[19]], col = 'purple')
lines(draws_anp_14[,plot_cols[20]],col = 'magenta')
lines(draws_anp_14[,plot_cols[21]],col = 'black')
lines(draws_anp_14[,plot_cols[22]], col = 'tan')
lines(draws_anp_14[,plot_cols[23]], col = 'orange')
lines(draws_anp_14[,plot_cols[24]], col = 'green')
lines(draws_anp_14[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_14[,plot_cols[26]], col = 'blue')
lines(draws_anp_14[,plot_cols[27]], col = 'red')
lines(draws_anp_14[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_14[,plot_cols[29]], col = 'purple')
lines(draws_anp_14[,plot_cols[30]],col = 'magenta')

### density plots  -- period 14 ####
dens(draws_anp_14[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_14[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 14 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_14 <- left_join(x = summaries, y = counts_df_period14, by = 'dyad')
head(plot_data_anp_14)
plot_data_anp_14$age_cat_1 <- ifelse(plot_data_anp_14$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_14$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_14$age_start_1 > 19, 'A','P')))
plot_data_anp_14$age_cat_2 <- ifelse(plot_data_anp_14$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_14$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_14$age_start_2 > 19, 'A','P')))
plot_data_anp_14$age_catnum_1 <- ifelse(plot_data_anp_14$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_14$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_14$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_14$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_14$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_14$age_start_1 < 40, 6, 7))))))
plot_data_anp_14$age_catnum_2 <- ifelse(plot_data_anp_14$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_14$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_14$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_14$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_14$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_14$age_start_2 < 40, 6, 7))))))

plot_data_anp_14$age_dyad <- ifelse(plot_data_anp_14$age_catnum_1 >= plot_data_anp_14$age_catnum_2,
                                    paste(plot_data_anp_14$age_cat_1, plot_data_anp_14$age_cat_2, sep = ''),
                                    paste(plot_data_anp_14$age_cat_2, plot_data_anp_14$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_14, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_14, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_14[plot_data_anp_14$age_dyad == 'AA' | plot_data_anp_14$age_dyad == 'AP' | plot_data_anp_14$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 14 ################
head(summaries)
length(unique(plot_data_anp_14$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids14 <- c(unique(counts_df_period14$id_1), unique(counts_df_period14$id_2)[length(unique(counts_df_period14$id_2))])
males14 <- males[,c(21,6,9,22:53)]
males14 <- males14 %>% dplyr::filter(id %in% ids14)
males14

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

which(males14$degree_0.1 < males14$degree_0.2)
which(males14$degree_0.2 < males14$degree_0.3)
which(males14$degree_0.3 < males14$degree_0.4)
which(males14$degree_0.4 < males14$degree_0.5)

# age variable
males14$age <- lubridate::year(periods$period_start[periods$period == 14]) - males14$byr
summary(males14$age)
males14$age_class <- ifelse(males14$age < 5, 1,
                            ifelse(males14$age < 10, 2,
                                   ifelse(males14$age < 15, 3,
                                          ifelse(males14$age < 20, 4,
                                                 ifelse(males14$age < 25, 5,
                                                        ifelse(males14$age < 40, 6, 7))))))

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
     vertex.label = males14$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males14$age_class == 7,'seagreen4',
                          ifelse(males14$age_class == 6,'seagreen3',
                                 ifelse(males14$age_class == 5,'seagreen2',
                                        ifelse(males14$age_class == 4,'steelblue3',
                                               ifelse(males14$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males14$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males14$age_class == 7,'seagreen4',
                          ifelse(males14$age_class == 6,'seagreen3',
                                 ifelse(males14$age_class == 5,'seagreen2',
                                        ifelse(males14$age_class == 4,'steelblue3',
                                               ifelse(males14$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 14 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males14$id[which(males14$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males14$id[which(males14$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males14[which(males14$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males14[which(males14$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males14[which(males14$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males14[which(males14$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males14[which(males14$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 14 ####
summaries$period <- 14
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period14, draws_anp_14, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males14, plot_data_anp_14, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids14, N, plot_cols)

################ 7.15) Run model on real standardised data -- period 15 ################
### create data list
counts_df_period15 <- counts_df_non0[counts_df_non0$period == 15,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period15),          # total number of times one or other of the dyad was observed
  together = counts_df_period15$event_count,    # count number of sightings seen together
  apart    = counts_df_period15$apart,          # count number of sightings seen apart
  period   = counts_df_period15$period)         # which period it's within

### Fit model
weight_anp_2.2_period15 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period15
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period15$output_files()[1])
draws1_anp_15 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period15$output_files()[2])
draws2_anp_15 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period15$output_files()[3])
draws3_anp_15 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period15$output_files()[4])
draws4_anp_15 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_15 <- rbind(draws1_anp_15, draws2_anp_15, draws3_anp_15, draws4_anp_15)

colnames(draws_anp_15)[2:ncol(draws_anp_15)] <- counts_df_period15$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_15), size = 30, replace = F)

### save data 
write_csv(draws_anp_15, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period15_22.04.27.csv')
rm(draws1_anp_15, draws2_anp_15, draws3_anp_15, draws4_anp_15)

### build traceplots  -- period 15 ####
plot(draws_anp_15[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_15[,plot_cols[2]], col = 'tan')
lines(draws_anp_15[,plot_cols[3]], col = 'orange')
lines(draws_anp_15[,plot_cols[4]], col = 'green')
lines(draws_anp_15[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_15[,plot_cols[6]], col = 'blue')
lines(draws_anp_15[,plot_cols[7]], col = 'red')
lines(draws_anp_15[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_15[,plot_cols[9]], col = 'purple')
lines(draws_anp_15[,plot_cols[10]],col = 'magenta')
lines(draws_anp_15[,plot_cols[11]],col = 'black')
lines(draws_anp_15[,plot_cols[12]], col = 'tan')
lines(draws_anp_15[,plot_cols[13]], col = 'orange')
lines(draws_anp_15[,plot_cols[14]], col = 'green')
lines(draws_anp_15[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_15[,plot_cols[16]], col = 'blue')
lines(draws_anp_15[,plot_cols[17]], col = 'red')
lines(draws_anp_15[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_15[,plot_cols[19]], col = 'purple')
lines(draws_anp_15[,plot_cols[20]],col = 'magenta')
lines(draws_anp_15[,plot_cols[21]],col = 'black')
lines(draws_anp_15[,plot_cols[22]], col = 'tan')
lines(draws_anp_15[,plot_cols[23]], col = 'orange')
lines(draws_anp_15[,plot_cols[24]], col = 'green')
lines(draws_anp_15[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_15[,plot_cols[26]], col = 'blue')
lines(draws_anp_15[,plot_cols[27]], col = 'red')
lines(draws_anp_15[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_15[,plot_cols[29]], col = 'purple')
lines(draws_anp_15[,plot_cols[30]],col = 'magenta')

### density plots  -- period 15 ####
dens(draws_anp_15[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_15[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 15 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_15 <- left_join(x = summaries, y = counts_df_period15, by = 'dyad')
head(plot_data_anp_15)
plot_data_anp_15$age_cat_1 <- ifelse(plot_data_anp_15$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_15$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_15$age_start_1 > 19, 'A','P')))
plot_data_anp_15$age_cat_2 <- ifelse(plot_data_anp_15$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_15$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_15$age_start_2 > 19, 'A','P')))
plot_data_anp_15$age_catnum_1 <- ifelse(plot_data_anp_15$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_15$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_15$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_15$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_15$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_15$age_start_1 < 40, 6, 7))))))
plot_data_anp_15$age_catnum_2 <- ifelse(plot_data_anp_15$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_15$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_15$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_15$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_15$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_15$age_start_2 < 40, 6, 7))))))

plot_data_anp_15$age_dyad <- ifelse(plot_data_anp_15$age_catnum_1 >= plot_data_anp_15$age_catnum_2,
                                    paste(plot_data_anp_15$age_cat_1, plot_data_anp_15$age_cat_2, sep = ''),
                                    paste(plot_data_anp_15$age_cat_2, plot_data_anp_15$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_15, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_15, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_15[plot_data_anp_15$age_dyad == 'AA' | plot_data_anp_15$age_dyad == 'AP' | plot_data_anp_15$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 15 ################
head(summaries)
length(unique(plot_data_anp_15$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
males15$age_class <- ifelse(males15$age < 5, 1,
                            ifelse(males15$age < 10, 2,
                                   ifelse(males15$age < 15, 3,
                                          ifelse(males15$age < 20, 4,
                                                 ifelse(males15$age < 25, 5,
                                                        ifelse(males15$age < 40, 6, 7))))))

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
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
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

### only elephants degree > 0.3 -- period 15 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males15$id[which(males15$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males15$id[which(males15$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males15[which(males15$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males15[which(males15$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males15[which(males15$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males15[which(males15$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males15[which(males15$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period15, draws_anp_15, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males15, plot_data_anp_15, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids15, N, plot_cols)

#### Save outputs to file #####
# write out csv file
write_csv(dyad_period_weights, 'data_processed/anp_dyad_weightdistributions_2.2_periods1to15_22.05.06.csv')

# produce PDF of graphs
dev.off()

################ Time windows 16:20 ################
# set seed
set.seed(12345)

# create file of output graphs
pdf('anp_edgeweights_2.2_period16to20_22.05.06.pdf', width = 10, height = 10)

# read in previous summary file
#dyad_period_weights <- read_csv('data_processed/anp_dyad_weightdistributions_2.2_periods1to15_22.05.06.csv')

################ 7.16) Run model on real standardised data -- period 16 ################
### create data list
counts_df_period16 <- counts_df_non0[counts_df_non0$period == 16,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period16),          # total number of times one or other of the dyad was observed
  together = counts_df_period16$event_count,    # count number of sightings seen together
  apart    = counts_df_period16$apart,          # count number of sightings seen apart
  period   = counts_df_period16$period)         # which period it's within

### Fit model
weight_anp_2.2_period16 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period16
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period16$output_files()[1])
draws1_anp_16 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period16$output_files()[2])
draws2_anp_16 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period16$output_files()[3])
draws3_anp_16 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period16$output_files()[4])
draws4_anp_16 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_16 <- rbind(draws1_anp_16, draws2_anp_16, draws3_anp_16, draws4_anp_16)

colnames(draws_anp_16)[2:ncol(draws_anp_16)] <- counts_df_period16$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_16), size = 30, replace = F)

### save data 
write_csv(draws_anp_16, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period16_22.04.27.csv')
rm(draws1_anp_16, draws2_anp_16, draws3_anp_16, draws4_anp_16)

### build traceplots  -- period 16 ####
plot(draws_anp_16[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_16[,plot_cols[2]], col = 'tan')
lines(draws_anp_16[,plot_cols[3]], col = 'orange')
lines(draws_anp_16[,plot_cols[4]], col = 'green')
lines(draws_anp_16[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_16[,plot_cols[6]], col = 'blue')
lines(draws_anp_16[,plot_cols[7]], col = 'red')
lines(draws_anp_16[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_16[,plot_cols[9]], col = 'purple')
lines(draws_anp_16[,plot_cols[10]],col = 'magenta')
lines(draws_anp_16[,plot_cols[11]],col = 'black')
lines(draws_anp_16[,plot_cols[12]], col = 'tan')
lines(draws_anp_16[,plot_cols[13]], col = 'orange')
lines(draws_anp_16[,plot_cols[14]], col = 'green')
lines(draws_anp_16[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_16[,plot_cols[16]], col = 'blue')
lines(draws_anp_16[,plot_cols[17]], col = 'red')
lines(draws_anp_16[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_16[,plot_cols[19]], col = 'purple')
lines(draws_anp_16[,plot_cols[20]],col = 'magenta')
lines(draws_anp_16[,plot_cols[21]],col = 'black')
lines(draws_anp_16[,plot_cols[22]], col = 'tan')
lines(draws_anp_16[,plot_cols[23]], col = 'orange')
lines(draws_anp_16[,plot_cols[24]], col = 'green')
lines(draws_anp_16[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_16[,plot_cols[26]], col = 'blue')
lines(draws_anp_16[,plot_cols[27]], col = 'red')
lines(draws_anp_16[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_16[,plot_cols[29]], col = 'purple')
lines(draws_anp_16[,plot_cols[30]],col = 'magenta')

### density plots  -- period 16 ####
dens(draws_anp_16[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_16[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 16 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_16 <- left_join(x = summaries, y = counts_df_period16, by = 'dyad')
head(plot_data_anp_16)
plot_data_anp_16$age_cat_1 <- ifelse(plot_data_anp_16$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_16$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_16$age_start_1 > 19, 'A','P')))
plot_data_anp_16$age_cat_2 <- ifelse(plot_data_anp_16$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_16$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_16$age_start_2 > 19, 'A','P')))
plot_data_anp_16$age_catnum_1 <- ifelse(plot_data_anp_16$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_16$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_16$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_16$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_16$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_16$age_start_1 < 40, 6, 7))))))
plot_data_anp_16$age_catnum_2 <- ifelse(plot_data_anp_16$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_16$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_16$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_16$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_16$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_16$age_start_2 < 40, 6, 7))))))

plot_data_anp_16$age_dyad <- ifelse(plot_data_anp_16$age_catnum_1 >= plot_data_anp_16$age_catnum_2,
                                    paste(plot_data_anp_16$age_cat_1, plot_data_anp_16$age_cat_2, sep = ''),
                                    paste(plot_data_anp_16$age_cat_2, plot_data_anp_16$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_16, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_16, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_16[plot_data_anp_16$age_dyad == 'AA' | plot_data_anp_16$age_dyad == 'AP' | plot_data_anp_16$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 16 ################
head(summaries)
length(unique(plot_data_anp_16$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids16 <- c(unique(counts_df_period16$id_1), unique(counts_df_period16$id_2)[length(unique(counts_df_period16$id_2))])
males16 <- males[,c(21,6,9,22:53)]
males16 <- males16 %>% dplyr::filter(id %in% ids16)
males16

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

which(males16$degree_0.1 < males16$degree_0.2)
which(males16$degree_0.2 < males16$degree_0.3)
which(males16$degree_0.3 < males16$degree_0.4)
which(males16$degree_0.4 < males16$degree_0.5)

# age variable
males16$age <- lubridate::year(periods$period_start[periods$period == 16]) - males16$byr
summary(males16$age)
males16$age_class <- ifelse(males16$age < 5, 1,
                            ifelse(males16$age < 10, 2,
                                   ifelse(males16$age < 15, 3,
                                          ifelse(males16$age < 20, 4,
                                                 ifelse(males16$age < 25, 5,
                                                        ifelse(males16$age < 40, 6, 7))))))

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
     vertex.label = males16$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males16$age_class == 7,'seagreen4',
                          ifelse(males16$age_class == 6,'seagreen3',
                                 ifelse(males16$age_class == 5,'seagreen2',
                                        ifelse(males16$age_class == 4,'steelblue3',
                                               ifelse(males16$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males16$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males16$age_class == 7,'seagreen4',
                          ifelse(males16$age_class == 6,'seagreen3',
                                 ifelse(males16$age_class == 5,'seagreen2',
                                        ifelse(males16$age_class == 4,'steelblue3',
                                               ifelse(males16$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 16 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males16$id[which(males16$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males16$id[which(males16$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males16[which(males16$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males16[which(males16$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males16[which(males16$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males16[which(males16$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males16[which(males16$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 16 ####
summaries$period <- 16
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period16, draws_anp_16, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males16, plot_data_anp_16, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids16, N, plot_cols)

################ 7.17) Run model on real standardised data -- period 17 ################
### create data list
counts_df_period17 <- counts_df_non0[counts_df_non0$period == 17,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period17),          # total number of times one or other of the dyad was observed
  together = counts_df_period17$event_count,    # count number of sightings seen together
  apart    = counts_df_period17$apart,          # count number of sightings seen apart
  period   = counts_df_period17$period)         # which period it's within

### Fit model
weight_anp_2.2_period17 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period17
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period17$output_files()[1])
draws1_anp_17 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period17$output_files()[2])
draws2_anp_17 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period17$output_files()[3])
draws3_anp_17 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period17$output_files()[4])
draws4_anp_17 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_17 <- rbind(draws1_anp_17, draws2_anp_17, draws3_anp_17, draws4_anp_17)

colnames(draws_anp_17)[2:ncol(draws_anp_17)] <- counts_df_period17$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_17), size = 30, replace = F)

### save data 
write_csv(draws_anp_17, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period17_22.04.27.csv')
rm(draws1_anp_17, draws2_anp_17, draws3_anp_17, draws4_anp_17)

### build traceplots  -- period 17 ####
plot(draws_anp_17[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_17[,plot_cols[2]], col = 'tan')
lines(draws_anp_17[,plot_cols[3]], col = 'orange')
lines(draws_anp_17[,plot_cols[4]], col = 'green')
lines(draws_anp_17[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_17[,plot_cols[6]], col = 'blue')
lines(draws_anp_17[,plot_cols[7]], col = 'red')
lines(draws_anp_17[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_17[,plot_cols[9]], col = 'purple')
lines(draws_anp_17[,plot_cols[10]],col = 'magenta')
lines(draws_anp_17[,plot_cols[11]],col = 'black')
lines(draws_anp_17[,plot_cols[12]], col = 'tan')
lines(draws_anp_17[,plot_cols[13]], col = 'orange')
lines(draws_anp_17[,plot_cols[14]], col = 'green')
lines(draws_anp_17[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_17[,plot_cols[16]], col = 'blue')
lines(draws_anp_17[,plot_cols[17]], col = 'red')
lines(draws_anp_17[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_17[,plot_cols[19]], col = 'purple')
lines(draws_anp_17[,plot_cols[20]],col = 'magenta')
lines(draws_anp_17[,plot_cols[21]],col = 'black')
lines(draws_anp_17[,plot_cols[22]], col = 'tan')
lines(draws_anp_17[,plot_cols[23]], col = 'orange')
lines(draws_anp_17[,plot_cols[24]], col = 'green')
lines(draws_anp_17[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_17[,plot_cols[26]], col = 'blue')
lines(draws_anp_17[,plot_cols[27]], col = 'red')
lines(draws_anp_17[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_17[,plot_cols[29]], col = 'purple')
lines(draws_anp_17[,plot_cols[30]],col = 'magenta')

### density plots  -- period 17 ####
dens(draws_anp_17[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_17[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 17 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_17 <- left_join(x = summaries, y = counts_df_period17, by = 'dyad')
head(plot_data_anp_17)
plot_data_anp_17$age_cat_1 <- ifelse(plot_data_anp_17$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_17$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_17$age_start_1 > 19, 'A','P')))
plot_data_anp_17$age_cat_2 <- ifelse(plot_data_anp_17$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_17$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_17$age_start_2 > 19, 'A','P')))
plot_data_anp_17$age_catnum_1 <- ifelse(plot_data_anp_17$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_17$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_17$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_17$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_17$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_17$age_start_1 < 40, 6, 7))))))
plot_data_anp_17$age_catnum_2 <- ifelse(plot_data_anp_17$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_17$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_17$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_17$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_17$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_17$age_start_2 < 40, 6, 7))))))

plot_data_anp_17$age_dyad <- ifelse(plot_data_anp_17$age_catnum_1 >= plot_data_anp_17$age_catnum_2,
                                    paste(plot_data_anp_17$age_cat_1, plot_data_anp_17$age_cat_2, sep = ''),
                                    paste(plot_data_anp_17$age_cat_2, plot_data_anp_17$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_17, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_17, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_17[plot_data_anp_17$age_dyad == 'AA' | plot_data_anp_17$age_dyad == 'AP' | plot_data_anp_17$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 17 ################
head(summaries)
length(unique(plot_data_anp_17$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids17 <- c(unique(counts_df_period17$id_1), unique(counts_df_period17$id_2)[length(unique(counts_df_period17$id_2))])
males17 <- males[,c(21,6,9,22:53)]
males17 <- males17 %>% dplyr::filter(id %in% ids17)
males17

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

which(males17$degree_0.1 < males17$degree_0.2)
which(males17$degree_0.2 < males17$degree_0.3)
which(males17$degree_0.3 < males17$degree_0.4)
which(males17$degree_0.4 < males17$degree_0.5)

# age variable
males17$age <- lubridate::year(periods$period_start[periods$period == 17]) - males17$byr
summary(males17$age)
males17$age_class <- ifelse(males17$age < 5, 1,
                            ifelse(males17$age < 10, 2,
                                   ifelse(males17$age < 15, 3,
                                          ifelse(males17$age < 20, 4,
                                                 ifelse(males17$age < 25, 5,
                                                        ifelse(males17$age < 40, 6, 7))))))

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
     vertex.label = males17$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males17$age_class == 7,'seagreen4',
                          ifelse(males17$age_class == 6,'seagreen3',
                                 ifelse(males17$age_class == 5,'seagreen2',
                                        ifelse(males17$age_class == 4,'steelblue3',
                                               ifelse(males17$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males17$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males17$age_class == 7,'seagreen4',
                          ifelse(males17$age_class == 6,'seagreen3',
                                 ifelse(males17$age_class == 5,'seagreen2',
                                        ifelse(males17$age_class == 4,'steelblue3',
                                               ifelse(males17$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 17 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males17$id[which(males17$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males17$id[which(males17$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males17[which(males17$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males17[which(males17$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males17[which(males17$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males17[which(males17$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males17[which(males17$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 17 ####
summaries$period <- 17
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period17, draws_anp_17, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males17, plot_data_anp_17, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids17, N, plot_cols)

################ 7.18) Run model on real standardised data -- period 18 ################
### create data list
counts_df_period18 <- counts_df_non0[counts_df_non0$period == 18,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period18),          # total number of times one or other of the dyad was observed
  together = counts_df_period18$event_count,    # count number of sightings seen together
  apart    = counts_df_period18$apart,          # count number of sightings seen apart
  period   = counts_df_period18$period)         # which period it's within

### Fit model
weight_anp_2.2_period18 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period18
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period18$output_files()[1])
draws1_anp_18 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period18$output_files()[2])
draws2_anp_18 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period18$output_files()[3])
draws3_anp_18 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period18$output_files()[4])
draws4_anp_18 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_18 <- rbind(draws1_anp_18, draws2_anp_18, draws3_anp_18, draws4_anp_18)

colnames(draws_anp_18)[2:ncol(draws_anp_18)] <- counts_df_period18$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_18), size = 30, replace = F)

### save data 
write_csv(draws_anp_18, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period18_22.04.27.csv')
rm(draws1_anp_18, draws2_anp_18, draws3_anp_18, draws4_anp_18)

### build traceplots  -- period 18 ####
plot(draws_anp_18[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_18[,plot_cols[2]], col = 'tan')
lines(draws_anp_18[,plot_cols[3]], col = 'orange')
lines(draws_anp_18[,plot_cols[4]], col = 'green')
lines(draws_anp_18[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_18[,plot_cols[6]], col = 'blue')
lines(draws_anp_18[,plot_cols[7]], col = 'red')
lines(draws_anp_18[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_18[,plot_cols[9]], col = 'purple')
lines(draws_anp_18[,plot_cols[10]],col = 'magenta')
lines(draws_anp_18[,plot_cols[11]],col = 'black')
lines(draws_anp_18[,plot_cols[12]], col = 'tan')
lines(draws_anp_18[,plot_cols[13]], col = 'orange')
lines(draws_anp_18[,plot_cols[14]], col = 'green')
lines(draws_anp_18[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_18[,plot_cols[16]], col = 'blue')
lines(draws_anp_18[,plot_cols[17]], col = 'red')
lines(draws_anp_18[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_18[,plot_cols[19]], col = 'purple')
lines(draws_anp_18[,plot_cols[20]],col = 'magenta')
lines(draws_anp_18[,plot_cols[21]],col = 'black')
lines(draws_anp_18[,plot_cols[22]], col = 'tan')
lines(draws_anp_18[,plot_cols[23]], col = 'orange')
lines(draws_anp_18[,plot_cols[24]], col = 'green')
lines(draws_anp_18[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_18[,plot_cols[26]], col = 'blue')
lines(draws_anp_18[,plot_cols[27]], col = 'red')
lines(draws_anp_18[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_18[,plot_cols[29]], col = 'purple')
lines(draws_anp_18[,plot_cols[30]],col = 'magenta')

### density plots  -- period 18 ####
dens(draws_anp_18[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_18[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 18 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_18 <- left_join(x = summaries, y = counts_df_period18, by = 'dyad')
head(plot_data_anp_18)
plot_data_anp_18$age_cat_1 <- ifelse(plot_data_anp_18$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_18$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_18$age_start_1 > 19, 'A','P')))
plot_data_anp_18$age_cat_2 <- ifelse(plot_data_anp_18$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_18$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_18$age_start_2 > 19, 'A','P')))
plot_data_anp_18$age_catnum_1 <- ifelse(plot_data_anp_18$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_18$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_18$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_18$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_18$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_18$age_start_1 < 40, 6, 7))))))
plot_data_anp_18$age_catnum_2 <- ifelse(plot_data_anp_18$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_18$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_18$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_18$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_18$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_18$age_start_2 < 40, 6, 7))))))

plot_data_anp_18$age_dyad <- ifelse(plot_data_anp_18$age_catnum_1 >= plot_data_anp_18$age_catnum_2,
                                    paste(plot_data_anp_18$age_cat_1, plot_data_anp_18$age_cat_2, sep = ''),
                                    paste(plot_data_anp_18$age_cat_2, plot_data_anp_18$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_18, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_18, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_18[plot_data_anp_18$age_dyad == 'AA' | plot_data_anp_18$age_dyad == 'AP' | plot_data_anp_18$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 18 ################
head(summaries)
length(unique(plot_data_anp_18$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids18 <- c(unique(counts_df_period18$id_1), unique(counts_df_period18$id_2)[length(unique(counts_df_period18$id_2))])
males18 <- males[,c(21,6,9,22:53)]
males18 <- males18 %>% dplyr::filter(id %in% ids18)
males18

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

which(males18$degree_0.1 < males18$degree_0.2)
which(males18$degree_0.2 < males18$degree_0.3)
which(males18$degree_0.3 < males18$degree_0.4)
which(males18$degree_0.4 < males18$degree_0.5)

# age variable
males18$age <- lubridate::year(periods$period_start[periods$period == 18]) - males18$byr
summary(males18$age)
males18$age_class <- ifelse(males18$age < 5, 1,
                            ifelse(males18$age < 10, 2,
                                   ifelse(males18$age < 15, 3,
                                          ifelse(males18$age < 20, 4,
                                                 ifelse(males18$age < 25, 5,
                                                        ifelse(males18$age < 40, 6, 7))))))

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
     vertex.label = males18$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males18$age_class == 7,'seagreen4',
                          ifelse(males18$age_class == 6,'seagreen3',
                                 ifelse(males18$age_class == 5,'seagreen2',
                                        ifelse(males18$age_class == 4,'steelblue3',
                                               ifelse(males18$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males18$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males18$age_class == 7,'seagreen4',
                          ifelse(males18$age_class == 6,'seagreen3',
                                 ifelse(males18$age_class == 5,'seagreen2',
                                        ifelse(males18$age_class == 4,'steelblue3',
                                               ifelse(males18$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 18 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males18$id[which(males18$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males18$id[which(males18$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males18[which(males18$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males18[which(males18$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males18[which(males18$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males18[which(males18$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males18[which(males18$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 18 ####
summaries$period <- 18
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period18, draws_anp_18, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males18, plot_data_anp_18, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids18, N, plot_cols)

################ 7.19) Run model on real standardised data -- period 19 ################
### create data list
counts_df_period19 <- counts_df_non0[counts_df_non0$period == 19,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period19),          # total number of times one or other of the dyad was observed
  together = counts_df_period19$event_count,    # count number of sightings seen together
  apart    = counts_df_period19$apart,          # count number of sightings seen apart
  period   = counts_df_period19$period)         # which period it's within

### Fit model
weight_anp_2.2_period19 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period19
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period19$output_files()[1])
draws1_anp_19 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period19$output_files()[2])
draws2_anp_19 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period19$output_files()[3])
draws3_anp_19 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period19$output_files()[4])
draws4_anp_19 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_19 <- rbind(draws1_anp_19, draws2_anp_19, draws3_anp_19, draws4_anp_19)

colnames(draws_anp_19)[2:ncol(draws_anp_19)] <- counts_df_period19$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_19), size = 30, replace = F)

### save data 
write_csv(draws_anp_19, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period19_22.04.27.csv')
rm(draws1_anp_19, draws2_anp_19, draws3_anp_19, draws4_anp_19)

### build traceplots  -- period 19 ####
plot(draws_anp_19[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_19[,plot_cols[2]], col = 'tan')
lines(draws_anp_19[,plot_cols[3]], col = 'orange')
lines(draws_anp_19[,plot_cols[4]], col = 'green')
lines(draws_anp_19[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_19[,plot_cols[6]], col = 'blue')
lines(draws_anp_19[,plot_cols[7]], col = 'red')
lines(draws_anp_19[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_19[,plot_cols[9]], col = 'purple')
lines(draws_anp_19[,plot_cols[10]],col = 'magenta')
lines(draws_anp_19[,plot_cols[11]],col = 'black')
lines(draws_anp_19[,plot_cols[12]], col = 'tan')
lines(draws_anp_19[,plot_cols[13]], col = 'orange')
lines(draws_anp_19[,plot_cols[14]], col = 'green')
lines(draws_anp_19[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_19[,plot_cols[16]], col = 'blue')
lines(draws_anp_19[,plot_cols[17]], col = 'red')
lines(draws_anp_19[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_19[,plot_cols[19]], col = 'purple')
lines(draws_anp_19[,plot_cols[20]],col = 'magenta')
lines(draws_anp_19[,plot_cols[21]],col = 'black')
lines(draws_anp_19[,plot_cols[22]], col = 'tan')
lines(draws_anp_19[,plot_cols[23]], col = 'orange')
lines(draws_anp_19[,plot_cols[24]], col = 'green')
lines(draws_anp_19[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_19[,plot_cols[26]], col = 'blue')
lines(draws_anp_19[,plot_cols[27]], col = 'red')
lines(draws_anp_19[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_19[,plot_cols[29]], col = 'purple')
lines(draws_anp_19[,plot_cols[30]],col = 'magenta')

### density plots  -- period 19 ####
dens(draws_anp_19[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_19[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 19 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_19 <- left_join(x = summaries, y = counts_df_period19, by = 'dyad')
head(plot_data_anp_19)
plot_data_anp_19$age_cat_1 <- ifelse(plot_data_anp_19$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_19$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_19$age_start_1 > 19, 'A','P')))
plot_data_anp_19$age_cat_2 <- ifelse(plot_data_anp_19$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_19$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_19$age_start_2 > 19, 'A','P')))
plot_data_anp_19$age_catnum_1 <- ifelse(plot_data_anp_19$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_19$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_19$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_19$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_19$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_19$age_start_1 < 40, 6, 7))))))
plot_data_anp_19$age_catnum_2 <- ifelse(plot_data_anp_19$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_19$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_19$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_19$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_19$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_19$age_start_2 < 40, 6, 7))))))

plot_data_anp_19$age_dyad <- ifelse(plot_data_anp_19$age_catnum_1 >= plot_data_anp_19$age_catnum_2,
                                    paste(plot_data_anp_19$age_cat_1, plot_data_anp_19$age_cat_2, sep = ''),
                                    paste(plot_data_anp_19$age_cat_2, plot_data_anp_19$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_19, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_19, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_19[plot_data_anp_19$age_dyad == 'AA' | plot_data_anp_19$age_dyad == 'AP' | plot_data_anp_19$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 19 ################
head(summaries)
length(unique(plot_data_anp_19$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids19 <- c(unique(counts_df_period19$id_1), unique(counts_df_period19$id_2)[length(unique(counts_df_period19$id_2))])
males19 <- males[,c(21,6,9,22:53)]
males19 <- males19 %>% dplyr::filter(id %in% ids19)
males19

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

which(males19$degree_0.1 < males19$degree_0.2)
which(males19$degree_0.2 < males19$degree_0.3)
which(males19$degree_0.3 < males19$degree_0.4)
which(males19$degree_0.4 < males19$degree_0.5)

# age variable
males19$age <- lubridate::year(periods$period_start[periods$period == 19]) - males19$byr
summary(males19$age)
males19$age_class <- ifelse(males19$age < 5, 1,
                            ifelse(males19$age < 10, 2,
                                   ifelse(males19$age < 15, 3,
                                          ifelse(males19$age < 20, 4,
                                                 ifelse(males19$age < 25, 5,
                                                        ifelse(males19$age < 40, 6, 7))))))

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
     vertex.label = males19$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males19$age_class == 7,'seagreen4',
                          ifelse(males19$age_class == 6,'seagreen3',
                                 ifelse(males19$age_class == 5,'seagreen2',
                                        ifelse(males19$age_class == 4,'steelblue3',
                                               ifelse(males19$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males19$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males19$age_class == 7,'seagreen4',
                          ifelse(males19$age_class == 6,'seagreen3',
                                 ifelse(males19$age_class == 5,'seagreen2',
                                        ifelse(males19$age_class == 4,'steelblue3',
                                               ifelse(males19$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 19 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males19$id[which(males19$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males19$id[which(males19$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males19[which(males19$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males19[which(males19$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males19[which(males19$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males19[which(males19$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males19[which(males19$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 19 ####
summaries$period <- 19
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period19, draws_anp_19, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males19, plot_data_anp_19, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids19, N, plot_cols)

################ 7.20) Run model on real standardised data -- period 20 ################
### create data list
counts_df_period20 <- counts_df_non0[counts_df_non0$period == 20,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period20),          # total number of times one or other of the dyad was observed
  together = counts_df_period20$event_count,    # count number of sightings seen together
  apart    = counts_df_period20$apart,          # count number of sightings seen apart
  period   = counts_df_period20$period)         # which period it's within

### Fit model
weight_anp_2.2_period20 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period20
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period20$output_files()[1])
draws1_anp_20 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period20$output_files()[2])
draws2_anp_20 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period20$output_files()[3])
draws3_anp_20 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period20$output_files()[4])
draws4_anp_20 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_20 <- rbind(draws1_anp_20, draws2_anp_20, draws3_anp_20, draws4_anp_20)

colnames(draws_anp_20)[2:ncol(draws_anp_20)] <- counts_df_period20$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_20), size = 30, replace = F)

### save data 
write_csv(draws_anp_20, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period20_22.04.27.csv')
rm(draws1_anp_20, draws2_anp_20, draws3_anp_20, draws4_anp_20)

### build traceplots  -- period 20 ####
plot(draws_anp_20[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_20[,plot_cols[2]], col = 'tan')
lines(draws_anp_20[,plot_cols[3]], col = 'orange')
lines(draws_anp_20[,plot_cols[4]], col = 'green')
lines(draws_anp_20[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_20[,plot_cols[6]], col = 'blue')
lines(draws_anp_20[,plot_cols[7]], col = 'red')
lines(draws_anp_20[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_20[,plot_cols[9]], col = 'purple')
lines(draws_anp_20[,plot_cols[10]],col = 'magenta')
lines(draws_anp_20[,plot_cols[11]],col = 'black')
lines(draws_anp_20[,plot_cols[12]], col = 'tan')
lines(draws_anp_20[,plot_cols[13]], col = 'orange')
lines(draws_anp_20[,plot_cols[14]], col = 'green')
lines(draws_anp_20[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_20[,plot_cols[16]], col = 'blue')
lines(draws_anp_20[,plot_cols[17]], col = 'red')
lines(draws_anp_20[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_20[,plot_cols[19]], col = 'purple')
lines(draws_anp_20[,plot_cols[20]],col = 'magenta')
lines(draws_anp_20[,plot_cols[21]],col = 'black')
lines(draws_anp_20[,plot_cols[22]], col = 'tan')
lines(draws_anp_20[,plot_cols[23]], col = 'orange')
lines(draws_anp_20[,plot_cols[24]], col = 'green')
lines(draws_anp_20[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_20[,plot_cols[26]], col = 'blue')
lines(draws_anp_20[,plot_cols[27]], col = 'red')
lines(draws_anp_20[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_20[,plot_cols[29]], col = 'purple')
lines(draws_anp_20[,plot_cols[30]],col = 'magenta')

### density plots  -- period 20 ####
dens(draws_anp_20[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_20[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 20 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_20 <- left_join(x = summaries, y = counts_df_period20, by = 'dyad')
head(plot_data_anp_20)
plot_data_anp_20$age_cat_1 <- ifelse(plot_data_anp_20$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_20$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_20$age_start_1 > 19, 'A','P')))
plot_data_anp_20$age_cat_2 <- ifelse(plot_data_anp_20$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_20$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_20$age_start_2 > 19, 'A','P')))
plot_data_anp_20$age_catnum_1 <- ifelse(plot_data_anp_20$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_20$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_20$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_20$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_20$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_20$age_start_1 < 40, 6, 7))))))
plot_data_anp_20$age_catnum_2 <- ifelse(plot_data_anp_20$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_20$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_20$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_20$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_20$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_20$age_start_2 < 40, 6, 7))))))

plot_data_anp_20$age_dyad <- ifelse(plot_data_anp_20$age_catnum_1 >= plot_data_anp_20$age_catnum_2,
                                    paste(plot_data_anp_20$age_cat_1, plot_data_anp_20$age_cat_2, sep = ''),
                                    paste(plot_data_anp_20$age_cat_2, plot_data_anp_20$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_20, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_20, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_20[plot_data_anp_20$age_dyad == 'AA' | plot_data_anp_20$age_dyad == 'AP' | plot_data_anp_20$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 20 ################
head(summaries)
length(unique(plot_data_anp_20$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids20 <- c(unique(counts_df_period20$id_1), unique(counts_df_period20$id_2)[length(unique(counts_df_period20$id_2))])
males20 <- males[,c(21,6,9,22:53)]
males20 <- males20 %>% dplyr::filter(id %in% ids20)
males20

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

which(males20$degree_0.1 < males20$degree_0.2)
which(males20$degree_0.2 < males20$degree_0.3)
which(males20$degree_0.3 < males20$degree_0.4)
which(males20$degree_0.4 < males20$degree_0.5)

# age variable
males20$age <- lubridate::year(periods$period_start[periods$period == 20]) - males20$byr
summary(males20$age)
males20$age_class <- ifelse(males20$age < 5, 1,
                            ifelse(males20$age < 10, 2,
                                   ifelse(males20$age < 15, 3,
                                          ifelse(males20$age < 20, 4,
                                                 ifelse(males20$age < 25, 5,
                                                        ifelse(males20$age < 40, 6, 7))))))

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
     vertex.label = males20$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males20$age_class == 7,'seagreen4',
                          ifelse(males20$age_class == 6,'seagreen3',
                                 ifelse(males20$age_class == 5,'seagreen2',
                                        ifelse(males20$age_class == 4,'steelblue3',
                                               ifelse(males20$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males20$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males20$age_class == 7,'seagreen4',
                          ifelse(males20$age_class == 6,'seagreen3',
                                 ifelse(males20$age_class == 5,'seagreen2',
                                        ifelse(males20$age_class == 4,'steelblue3',
                                               ifelse(males20$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 20 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males20$id[which(males20$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males20$id[which(males20$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males20[which(males20$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males20[which(males20$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males20[which(males20$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males20[which(males20$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males20[which(males20$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 20 ####
summaries$period <- 20
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period20, draws_anp_20, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males20, plot_data_anp_20, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids20, N, plot_cols)

#### Save outputs to file #####
# write out csv file
write_csv(dyad_period_weights, 'data_processed/anp_dyad_weightdistributions_2.2_periods1to20_22.05.06.csv')

# produce PDF of graphs
dev.off()

# add progress stamp
print(paste0('Time windows 16:20 completed at ', Sys.time()))

######## periods 21:24 #########
# set seed
set.seed(12345)

# create file of output graphs
pdf('anp_edgeweights_2.2_period21to24_22.05.06.pdf', width = 10, height = 10)

# read in previous summary file
#dyad_period_weights <- read_csv('data_processed/anp_dyad_weightdistributions_2.2_periods1to20_22.05.06.csv')

################ 7.21) Run model on real standardised data -- period 21 ################
### create data list
counts_df_period21 <- counts_df_non0[counts_df_non0$period == 21,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period21),          # total number of times one or other of the dyad was observed
  together = counts_df_period21$event_count,    # count number of sightings seen together
  apart    = counts_df_period21$apart,          # count number of sightings seen apart
  period   = counts_df_period21$period)         # which period it's within

### Fit model
weight_anp_2.2_period21 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period21
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period21$output_files()[1])
draws1_anp_21 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period21$output_files()[2])
draws2_anp_21 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period21$output_files()[3])
draws3_anp_21 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period21$output_files()[4])
draws4_anp_21 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_21 <- rbind(draws1_anp_21, draws2_anp_21, draws3_anp_21, draws4_anp_21)

colnames(draws_anp_21)[2:ncol(draws_anp_21)] <- counts_df_period21$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_21), size = 30, replace = F)

### save data 
write_csv(draws_anp_21, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period21_22.04.27.csv')
rm(draws1_anp_21, draws2_anp_21, draws3_anp_21, draws4_anp_21)

### build traceplots  -- period 21 ####
plot(draws_anp_21[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_21[,plot_cols[2]], col = 'tan')
lines(draws_anp_21[,plot_cols[3]], col = 'orange')
lines(draws_anp_21[,plot_cols[4]], col = 'green')
lines(draws_anp_21[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_21[,plot_cols[6]], col = 'blue')
lines(draws_anp_21[,plot_cols[7]], col = 'red')
lines(draws_anp_21[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_21[,plot_cols[9]], col = 'purple')
lines(draws_anp_21[,plot_cols[10]],col = 'magenta')
lines(draws_anp_21[,plot_cols[11]],col = 'black')
lines(draws_anp_21[,plot_cols[12]], col = 'tan')
lines(draws_anp_21[,plot_cols[13]], col = 'orange')
lines(draws_anp_21[,plot_cols[14]], col = 'green')
lines(draws_anp_21[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_21[,plot_cols[16]], col = 'blue')
lines(draws_anp_21[,plot_cols[17]], col = 'red')
lines(draws_anp_21[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_21[,plot_cols[19]], col = 'purple')
lines(draws_anp_21[,plot_cols[20]],col = 'magenta')
lines(draws_anp_21[,plot_cols[21]],col = 'black')
lines(draws_anp_21[,plot_cols[22]], col = 'tan')
lines(draws_anp_21[,plot_cols[23]], col = 'orange')
lines(draws_anp_21[,plot_cols[24]], col = 'green')
lines(draws_anp_21[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_21[,plot_cols[26]], col = 'blue')
lines(draws_anp_21[,plot_cols[27]], col = 'red')
lines(draws_anp_21[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_21[,plot_cols[29]], col = 'purple')
lines(draws_anp_21[,plot_cols[30]],col = 'magenta')

### density plots  -- period 21 ####
dens(draws_anp_21[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_21[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 21 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_21 <- left_join(x = summaries, y = counts_df_period21, by = 'dyad')
head(plot_data_anp_21)
plot_data_anp_21$age_cat_1 <- ifelse(plot_data_anp_21$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_21$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_21$age_start_1 > 19, 'A','P')))
plot_data_anp_21$age_cat_2 <- ifelse(plot_data_anp_21$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_21$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_21$age_start_2 > 19, 'A','P')))
plot_data_anp_21$age_catnum_1 <- ifelse(plot_data_anp_21$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_21$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_21$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_21$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_21$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_21$age_start_1 < 40, 6, 7))))))
plot_data_anp_21$age_catnum_2 <- ifelse(plot_data_anp_21$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_21$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_21$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_21$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_21$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_21$age_start_2 < 40, 6, 7))))))

plot_data_anp_21$age_dyad <- ifelse(plot_data_anp_21$age_catnum_1 >= plot_data_anp_21$age_catnum_2,
                                    paste(plot_data_anp_21$age_cat_1, plot_data_anp_21$age_cat_2, sep = ''),
                                    paste(plot_data_anp_21$age_cat_2, plot_data_anp_21$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_21, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_21, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_21[plot_data_anp_21$age_dyad == 'AA' | plot_data_anp_21$age_dyad == 'AP' | plot_data_anp_21$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 21 ################
head(summaries)
length(unique(plot_data_anp_21$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids21 <- c(unique(counts_df_period21$id_1), unique(counts_df_period21$id_2)[length(unique(counts_df_period21$id_2))])
males21 <- males[,c(21,6,9,22:53)]
males21 <- males21 %>% dplyr::filter(id %in% ids21)
males21

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

which(males21$degree_0.1 < males21$degree_0.2)
which(males21$degree_0.2 < males21$degree_0.3)
which(males21$degree_0.3 < males21$degree_0.4)
which(males21$degree_0.4 < males21$degree_0.5)

# age variable
males21$age <- lubridate::year(periods$period_start[periods$period == 21]) - males21$byr
summary(males21$age)
males21$age_class <- ifelse(males21$age < 5, 1,
                            ifelse(males21$age < 10, 2,
                                   ifelse(males21$age < 15, 3,
                                          ifelse(males21$age < 20, 4,
                                                 ifelse(males21$age < 25, 5,
                                                        ifelse(males21$age < 40, 6, 7))))))

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
     vertex.label = males21$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males21$age_class == 7,'seagreen4',
                          ifelse(males21$age_class == 6,'seagreen3',
                                 ifelse(males21$age_class == 5,'seagreen2',
                                        ifelse(males21$age_class == 4,'steelblue3',
                                               ifelse(males21$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males21$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males21$age_class == 7,'seagreen4',
                          ifelse(males21$age_class == 6,'seagreen3',
                                 ifelse(males21$age_class == 5,'seagreen2',
                                        ifelse(males21$age_class == 4,'steelblue3',
                                               ifelse(males21$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 21 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males21$id[which(males21$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males21$id[which(males21$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males21[which(males21$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males21[which(males21$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males21[which(males21$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males21[which(males21$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males21[which(males21$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 21 ####
summaries$period <- 21
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period21, draws_anp_21, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males21, plot_data_anp_21, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids21, N, plot_cols)

################ 7.22) Run model on real standardised data -- period 22 ################
### create data list
counts_df_period22 <- counts_df_non0[counts_df_non0$period == 22,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period22),          # total number of times one or other of the dyad was observed
  together = counts_df_period22$event_count,    # count number of sightings seen together
  apart    = counts_df_period22$apart,          # count number of sightings seen apart
  period   = counts_df_period22$period)         # which period it's within

### Fit model
weight_anp_2.2_period22 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period22
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period22$output_files()[1])
draws1_anp_22 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period22$output_files()[2])
draws2_anp_22 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period22$output_files()[3])
draws3_anp_22 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period22$output_files()[4])
draws4_anp_22 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_22 <- rbind(draws1_anp_22, draws2_anp_22, draws3_anp_22, draws4_anp_22)

colnames(draws_anp_22)[2:ncol(draws_anp_22)] <- counts_df_period22$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_22), size = 30, replace = F)

### save data 
write_csv(draws_anp_22, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period22_22.04.27.csv')
rm(draws1_anp_22, draws2_anp_22, draws3_anp_22, draws4_anp_22)

### build traceplots  -- period 22 ####
plot(draws_anp_22[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_22[,plot_cols[2]], col = 'tan')
lines(draws_anp_22[,plot_cols[3]], col = 'orange')
lines(draws_anp_22[,plot_cols[4]], col = 'green')
lines(draws_anp_22[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_22[,plot_cols[6]], col = 'blue')
lines(draws_anp_22[,plot_cols[7]], col = 'red')
lines(draws_anp_22[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_22[,plot_cols[9]], col = 'purple')
lines(draws_anp_22[,plot_cols[10]],col = 'magenta')
lines(draws_anp_22[,plot_cols[11]],col = 'black')
lines(draws_anp_22[,plot_cols[12]], col = 'tan')
lines(draws_anp_22[,plot_cols[13]], col = 'orange')
lines(draws_anp_22[,plot_cols[14]], col = 'green')
lines(draws_anp_22[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_22[,plot_cols[16]], col = 'blue')
lines(draws_anp_22[,plot_cols[17]], col = 'red')
lines(draws_anp_22[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_22[,plot_cols[19]], col = 'purple')
lines(draws_anp_22[,plot_cols[20]],col = 'magenta')
lines(draws_anp_22[,plot_cols[21]],col = 'black')
lines(draws_anp_22[,plot_cols[22]], col = 'tan')
lines(draws_anp_22[,plot_cols[23]], col = 'orange')
lines(draws_anp_22[,plot_cols[24]], col = 'green')
lines(draws_anp_22[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_22[,plot_cols[26]], col = 'blue')
lines(draws_anp_22[,plot_cols[27]], col = 'red')
lines(draws_anp_22[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_22[,plot_cols[29]], col = 'purple')
lines(draws_anp_22[,plot_cols[30]],col = 'magenta')

### density plots  -- period 22 ####
dens(draws_anp_22[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_22[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 22 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_22 <- left_join(x = summaries, y = counts_df_period22, by = 'dyad')
head(plot_data_anp_22)
plot_data_anp_22$age_cat_1 <- ifelse(plot_data_anp_22$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_22$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_22$age_start_1 > 19, 'A','P')))
plot_data_anp_22$age_cat_2 <- ifelse(plot_data_anp_22$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_22$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_22$age_start_2 > 19, 'A','P')))
plot_data_anp_22$age_catnum_1 <- ifelse(plot_data_anp_22$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_22$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_22$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_22$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_22$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_22$age_start_1 < 40, 6, 7))))))
plot_data_anp_22$age_catnum_2 <- ifelse(plot_data_anp_22$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_22$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_22$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_22$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_22$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_22$age_start_2 < 40, 6, 7))))))

plot_data_anp_22$age_dyad <- ifelse(plot_data_anp_22$age_catnum_1 >= plot_data_anp_22$age_catnum_2,
                                    paste(plot_data_anp_22$age_cat_1, plot_data_anp_22$age_cat_2, sep = ''),
                                    paste(plot_data_anp_22$age_cat_2, plot_data_anp_22$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_22, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_22, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_22[plot_data_anp_22$age_dyad == 'AA' | plot_data_anp_22$age_dyad == 'AP' | plot_data_anp_22$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 22 ################
head(summaries)
length(unique(plot_data_anp_22$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids22 <- c(unique(counts_df_period22$id_1), unique(counts_df_period22$id_2)[length(unique(counts_df_period22$id_2))])
males22 <- males[,c(21,6,9,22:53)]
males22 <- males22 %>% dplyr::filter(id %in% ids22)
males22

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

which(males22$degree_0.1 < males22$degree_0.2)
which(males22$degree_0.2 < males22$degree_0.3)
which(males22$degree_0.3 < males22$degree_0.4)
which(males22$degree_0.4 < males22$degree_0.5)

# age variable
males22$age <- lubridate::year(periods$period_start[periods$period == 22]) - males22$byr
summary(males22$age)
males22$age_class <- ifelse(males22$age < 5, 1,
                            ifelse(males22$age < 10, 2,
                                   ifelse(males22$age < 15, 3,
                                          ifelse(males22$age < 20, 4,
                                                 ifelse(males22$age < 25, 5,
                                                        ifelse(males22$age < 40, 6, 7))))))

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
     vertex.label = males22$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males22$age_class == 7,'seagreen4',
                          ifelse(males22$age_class == 6,'seagreen3',
                                 ifelse(males22$age_class == 5,'seagreen2',
                                        ifelse(males22$age_class == 4,'steelblue3',
                                               ifelse(males22$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males22$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males22$age_class == 7,'seagreen4',
                          ifelse(males22$age_class == 6,'seagreen3',
                                 ifelse(males22$age_class == 5,'seagreen2',
                                        ifelse(males22$age_class == 4,'steelblue3',
                                               ifelse(males22$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 22 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males22$id[which(males22$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males22$id[which(males22$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males22[which(males22$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males22[which(males22$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males22[which(males22$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males22[which(males22$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males22[which(males22$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 22 ####
summaries$period <- 22
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period22, draws_anp_22, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males22, plot_data_anp_22, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids22, N, plot_cols)

################ 7.23) Run model on real standardised data -- period 23 ################
### create data list
counts_df_period23 <- counts_df_non0[counts_df_non0$period == 23,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period23),          # total number of times one or other of the dyad was observed
  together = counts_df_period23$event_count,    # count number of sightings seen together
  apart    = counts_df_period23$apart,          # count number of sightings seen apart
  period   = counts_df_period23$period)         # which period it's within

### Fit model
weight_anp_2.2_period23 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period23
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period23$output_files()[1])
draws1_anp_23 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period23$output_files()[2])
draws2_anp_23 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period23$output_files()[3])
draws3_anp_23 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period23$output_files()[4])
draws4_anp_23 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_23 <- rbind(draws1_anp_23, draws2_anp_23, draws3_anp_23, draws4_anp_23)

colnames(draws_anp_23)[2:ncol(draws_anp_23)] <- counts_df_period23$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_23), size = 30, replace = F)

### save data 
write_csv(draws_anp_23, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period23_22.04.27.csv')
rm(draws1_anp_23, draws2_anp_23, draws3_anp_23, draws4_anp_23)

### build traceplots  -- period 23 ####
plot(draws_anp_23[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_23[,plot_cols[2]], col = 'tan')
lines(draws_anp_23[,plot_cols[3]], col = 'orange')
lines(draws_anp_23[,plot_cols[4]], col = 'green')
lines(draws_anp_23[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_23[,plot_cols[6]], col = 'blue')
lines(draws_anp_23[,plot_cols[7]], col = 'red')
lines(draws_anp_23[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_23[,plot_cols[9]], col = 'purple')
lines(draws_anp_23[,plot_cols[10]],col = 'magenta')
lines(draws_anp_23[,plot_cols[11]],col = 'black')
lines(draws_anp_23[,plot_cols[12]], col = 'tan')
lines(draws_anp_23[,plot_cols[13]], col = 'orange')
lines(draws_anp_23[,plot_cols[14]], col = 'green')
lines(draws_anp_23[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_23[,plot_cols[16]], col = 'blue')
lines(draws_anp_23[,plot_cols[17]], col = 'red')
lines(draws_anp_23[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_23[,plot_cols[19]], col = 'purple')
lines(draws_anp_23[,plot_cols[20]],col = 'magenta')
lines(draws_anp_23[,plot_cols[21]],col = 'black')
lines(draws_anp_23[,plot_cols[22]], col = 'tan')
lines(draws_anp_23[,plot_cols[23]], col = 'orange')
lines(draws_anp_23[,plot_cols[24]], col = 'green')
lines(draws_anp_23[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_23[,plot_cols[26]], col = 'blue')
lines(draws_anp_23[,plot_cols[27]], col = 'red')
lines(draws_anp_23[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_23[,plot_cols[29]], col = 'purple')
lines(draws_anp_23[,plot_cols[30]],col = 'magenta')

### density plots  -- period 23 ####
dens(draws_anp_23[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_23[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 23 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_23 <- left_join(x = summaries, y = counts_df_period23, by = 'dyad')
head(plot_data_anp_23)
plot_data_anp_23$age_cat_1 <- ifelse(plot_data_anp_23$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_23$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_23$age_start_1 > 19, 'A','P')))
plot_data_anp_23$age_cat_2 <- ifelse(plot_data_anp_23$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_23$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_23$age_start_2 > 19, 'A','P')))
plot_data_anp_23$age_catnum_1 <- ifelse(plot_data_anp_23$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_23$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_23$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_23$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_23$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_23$age_start_1 < 40, 6, 7))))))
plot_data_anp_23$age_catnum_2 <- ifelse(plot_data_anp_23$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_23$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_23$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_23$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_23$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_23$age_start_2 < 40, 6, 7))))))

plot_data_anp_23$age_dyad <- ifelse(plot_data_anp_23$age_catnum_1 >= plot_data_anp_23$age_catnum_2,
                                    paste(plot_data_anp_23$age_cat_1, plot_data_anp_23$age_cat_2, sep = ''),
                                    paste(plot_data_anp_23$age_cat_2, plot_data_anp_23$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_23, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_23, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_23[plot_data_anp_23$age_dyad == 'AA' | plot_data_anp_23$age_dyad == 'AP' | plot_data_anp_23$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 23 ################
head(summaries)
length(unique(plot_data_anp_23$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids23 <- c(unique(counts_df_period23$id_1), unique(counts_df_period23$id_2)[length(unique(counts_df_period23$id_2))])
males23 <- males[,c(21,6,9,22:53)]
males23 <- males23 %>% dplyr::filter(id %in% ids23)
males23

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

which(males23$degree_0.1 < males23$degree_0.2)
which(males23$degree_0.2 < males23$degree_0.3)
which(males23$degree_0.3 < males23$degree_0.4)
which(males23$degree_0.4 < males23$degree_0.5)

# age variable
males23$age <- lubridate::year(periods$period_start[periods$period == 23]) - males23$byr
summary(males23$age)
males23$age_class <- ifelse(males23$age < 5, 1,
                            ifelse(males23$age < 10, 2,
                                   ifelse(males23$age < 15, 3,
                                          ifelse(males23$age < 20, 4,
                                                 ifelse(males23$age < 25, 5,
                                                        ifelse(males23$age < 40, 6, 7))))))

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
     vertex.label = males23$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males23$age_class == 7,'seagreen4',
                          ifelse(males23$age_class == 6,'seagreen3',
                                 ifelse(males23$age_class == 5,'seagreen2',
                                        ifelse(males23$age_class == 4,'steelblue3',
                                               ifelse(males23$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males23$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males23$age_class == 7,'seagreen4',
                          ifelse(males23$age_class == 6,'seagreen3',
                                 ifelse(males23$age_class == 5,'seagreen2',
                                        ifelse(males23$age_class == 4,'steelblue3',
                                               ifelse(males23$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 23 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males23$id[which(males23$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males23$id[which(males23$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males23[which(males23$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males23[which(males23$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males23[which(males23$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males23[which(males23$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males23[which(males23$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 23 ####
summaries$period <- 23
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period23, draws_anp_23, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males23, plot_data_anp_23, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids23, N, plot_cols)

################ 7.24) Run model on real standardised data -- period 24 ################
### create data list
counts_df_period24 <- counts_df_non0[counts_df_non0$period == 24,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period24),          # total number of times one or other of the dyad was observed
  together = counts_df_period24$event_count,    # count number of sightings seen together
  apart    = counts_df_period24$apart,          # count number of sightings seen apart
  period   = counts_df_period24$period)         # which period it's within

### Fit model
weight_anp_2.2_period24 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period24
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period24$output_files()[1])
draws1_anp_24 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period24$output_files()[2])
draws2_anp_24 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period24$output_files()[3])
draws3_anp_24 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period24$output_files()[4])
draws4_anp_24 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_24 <- rbind(draws1_anp_24, draws2_anp_24, draws3_anp_24, draws4_anp_24)

colnames(draws_anp_24)[2:ncol(draws_anp_24)] <- counts_df_period24$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_24), size = 30, replace = F)

### save data 
write_csv(draws_anp_24, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period24_22.04.27.csv')
rm(draws1_anp_24, draws2_anp_24, draws3_anp_24, draws4_anp_24)

### build traceplots  -- period 24 ####
plot(draws_anp_24[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_24[,plot_cols[2]], col = 'tan')
lines(draws_anp_24[,plot_cols[3]], col = 'orange')
lines(draws_anp_24[,plot_cols[4]], col = 'green')
lines(draws_anp_24[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_24[,plot_cols[6]], col = 'blue')
lines(draws_anp_24[,plot_cols[7]], col = 'red')
lines(draws_anp_24[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_24[,plot_cols[9]], col = 'purple')
lines(draws_anp_24[,plot_cols[10]],col = 'magenta')
lines(draws_anp_24[,plot_cols[11]],col = 'black')
lines(draws_anp_24[,plot_cols[12]], col = 'tan')
lines(draws_anp_24[,plot_cols[13]], col = 'orange')
lines(draws_anp_24[,plot_cols[14]], col = 'green')
lines(draws_anp_24[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_24[,plot_cols[16]], col = 'blue')
lines(draws_anp_24[,plot_cols[17]], col = 'red')
lines(draws_anp_24[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_24[,plot_cols[19]], col = 'purple')
lines(draws_anp_24[,plot_cols[20]],col = 'magenta')
lines(draws_anp_24[,plot_cols[21]],col = 'black')
lines(draws_anp_24[,plot_cols[22]], col = 'tan')
lines(draws_anp_24[,plot_cols[23]], col = 'orange')
lines(draws_anp_24[,plot_cols[24]], col = 'green')
lines(draws_anp_24[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_24[,plot_cols[26]], col = 'blue')
lines(draws_anp_24[,plot_cols[27]], col = 'red')
lines(draws_anp_24[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_24[,plot_cols[29]], col = 'purple')
lines(draws_anp_24[,plot_cols[30]],col = 'magenta')

### density plots  -- period 24 ####
dens(draws_anp_24[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_24[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 24 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_24 <- left_join(x = summaries, y = counts_df_period24, by = 'dyad')
head(plot_data_anp_24)
plot_data_anp_24$age_cat_1 <- ifelse(plot_data_anp_24$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_24$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_24$age_start_1 > 19, 'A','P')))
plot_data_anp_24$age_cat_2 <- ifelse(plot_data_anp_24$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_24$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_24$age_start_2 > 19, 'A','P')))
plot_data_anp_24$age_catnum_1 <- ifelse(plot_data_anp_24$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_24$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_24$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_24$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_24$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_24$age_start_1 < 40, 6, 7))))))
plot_data_anp_24$age_catnum_2 <- ifelse(plot_data_anp_24$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_24$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_24$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_24$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_24$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_24$age_start_2 < 40, 6, 7))))))

plot_data_anp_24$age_dyad <- ifelse(plot_data_anp_24$age_catnum_1 >= plot_data_anp_24$age_catnum_2,
                                    paste(plot_data_anp_24$age_cat_1, plot_data_anp_24$age_cat_2, sep = ''),
                                    paste(plot_data_anp_24$age_cat_2, plot_data_anp_24$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_24, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_24, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_24[plot_data_anp_24$age_dyad == 'AA' | plot_data_anp_24$age_dyad == 'AP' | plot_data_anp_24$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 24 ################
head(summaries)
length(unique(plot_data_anp_24$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids24 <- c(unique(counts_df_period24$id_1), unique(counts_df_period24$id_2)[length(unique(counts_df_period24$id_2))])
males24 <- males[,c(21,6,9,22:53)]
males24 <- males24 %>% dplyr::filter(id %in% ids24)
males24

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

which(males24$degree_0.1 < males24$degree_0.2)
which(males24$degree_0.2 < males24$degree_0.3)
which(males24$degree_0.3 < males24$degree_0.4)
which(males24$degree_0.4 < males24$degree_0.5)

# age variable
males24$age <- lubridate::year(periods$period_start[periods$period == 24]) - males24$byr
summary(males24$age)
males24$age_class <- ifelse(males24$age < 5, 1,
                            ifelse(males24$age < 10, 2,
                                   ifelse(males24$age < 15, 3,
                                          ifelse(males24$age < 20, 4,
                                                 ifelse(males24$age < 25, 5,
                                                        ifelse(males24$age < 40, 6, 7))))))

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
     vertex.label = males24$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males24$age_class == 7,'seagreen4',
                          ifelse(males24$age_class == 6,'seagreen3',
                                 ifelse(males24$age_class == 5,'seagreen2',
                                        ifelse(males24$age_class == 4,'steelblue3',
                                               ifelse(males24$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males24$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males24$age_class == 7,'seagreen4',
                          ifelse(males24$age_class == 6,'seagreen3',
                                 ifelse(males24$age_class == 5,'seagreen2',
                                        ifelse(males24$age_class == 4,'steelblue3',
                                               ifelse(males24$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 24 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males24$id[which(males24$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males24$id[which(males24$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males24[which(males24$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males24[which(males24$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males24[which(males24$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males24[which(males24$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males24[which(males24$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 24 ####
summaries$period <- 24
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period24, draws_anp_24, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males24, plot_data_anp_24, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids24, N, plot_cols)

#### Save outputs to file #####
# write out csv file
write_csv(dyad_period_weights, 'data_processed/anp_dyad_weightdistributions_2.2_periods1to24_22.05.06.csv')

# produce PDF of graphs
dev.off()

# add progress stamp
print(paste0('Time windows 21:24 completed at ', Sys.time()))

##################### time windows 25:28 ########################
# set seed
set.seed(12345)

# create file of output graphs
pdf('anp_edgeweights_2.2_period25to28_22.05.06.pdf', width = 10, height = 10)

# read in previous summary file
#dyad_period_weights <- read_csv('data_processed/anp_dyad_weightdistributions_2.2_periods1to24_22.05.06.csv')

################ 7.25) Run model on real standardised data -- period 25 ################
### create data list
counts_df_period25 <- counts_df_non0[counts_df_non0$period == 25,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period25),          # total number of times one or other of the dyad was observed
  together = counts_df_period25$event_count,    # count number of sightings seen together
  apart    = counts_df_period25$apart,          # count number of sightings seen apart
  period   = counts_df_period25$period)         # which period it's within

### Fit model
weight_anp_2.2_period25 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period25
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period25$output_files()[1])
draws1_anp_25 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period25$output_files()[2])
draws2_anp_25 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period25$output_files()[3])
draws3_anp_25 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period25$output_files()[4])
draws4_anp_25 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_25 <- rbind(draws1_anp_25, draws2_anp_25, draws3_anp_25, draws4_anp_25)

colnames(draws_anp_25)[2:ncol(draws_anp_25)] <- counts_df_period25$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_25), size = 30, replace = F)

### save data 
write_csv(draws_anp_25, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period25_22.04.27.csv')
rm(draws1_anp_25, draws2_anp_25, draws3_anp_25, draws4_anp_25)

### build traceplots  -- period 25 ####
plot(draws_anp_25[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_25[,plot_cols[2]], col = 'tan')
lines(draws_anp_25[,plot_cols[3]], col = 'orange')
lines(draws_anp_25[,plot_cols[4]], col = 'green')
lines(draws_anp_25[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_25[,plot_cols[6]], col = 'blue')
lines(draws_anp_25[,plot_cols[7]], col = 'red')
lines(draws_anp_25[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_25[,plot_cols[9]], col = 'purple')
lines(draws_anp_25[,plot_cols[10]],col = 'magenta')
lines(draws_anp_25[,plot_cols[11]],col = 'black')
lines(draws_anp_25[,plot_cols[12]], col = 'tan')
lines(draws_anp_25[,plot_cols[13]], col = 'orange')
lines(draws_anp_25[,plot_cols[14]], col = 'green')
lines(draws_anp_25[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_25[,plot_cols[16]], col = 'blue')
lines(draws_anp_25[,plot_cols[17]], col = 'red')
lines(draws_anp_25[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_25[,plot_cols[19]], col = 'purple')
lines(draws_anp_25[,plot_cols[20]],col = 'magenta')
lines(draws_anp_25[,plot_cols[21]],col = 'black')
lines(draws_anp_25[,plot_cols[22]], col = 'tan')
lines(draws_anp_25[,plot_cols[23]], col = 'orange')
lines(draws_anp_25[,plot_cols[24]], col = 'green')
lines(draws_anp_25[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_25[,plot_cols[26]], col = 'blue')
lines(draws_anp_25[,plot_cols[27]], col = 'red')
lines(draws_anp_25[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_25[,plot_cols[29]], col = 'purple')
lines(draws_anp_25[,plot_cols[30]],col = 'magenta')

### density plots  -- period 25 ####
dens(draws_anp_25[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_25[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 25 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_25 <- left_join(x = summaries, y = counts_df_period25, by = 'dyad')
head(plot_data_anp_25)
plot_data_anp_25$age_cat_1 <- ifelse(plot_data_anp_25$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_25$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_25$age_start_1 > 19, 'A','P')))
plot_data_anp_25$age_cat_2 <- ifelse(plot_data_anp_25$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_25$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_25$age_start_2 > 19, 'A','P')))
plot_data_anp_25$age_catnum_1 <- ifelse(plot_data_anp_25$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_25$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_25$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_25$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_25$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_25$age_start_1 < 40, 6, 7))))))
plot_data_anp_25$age_catnum_2 <- ifelse(plot_data_anp_25$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_25$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_25$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_25$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_25$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_25$age_start_2 < 40, 6, 7))))))

plot_data_anp_25$age_dyad <- ifelse(plot_data_anp_25$age_catnum_1 >= plot_data_anp_25$age_catnum_2,
                                    paste(plot_data_anp_25$age_cat_1, plot_data_anp_25$age_cat_2, sep = ''),
                                    paste(plot_data_anp_25$age_cat_2, plot_data_anp_25$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_25, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_25, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_25[plot_data_anp_25$age_dyad == 'AA' | plot_data_anp_25$age_dyad == 'AP' | plot_data_anp_25$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 25 ################
head(summaries)
length(unique(plot_data_anp_25$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids25 <- c(unique(counts_df_period25$id_1), unique(counts_df_period25$id_2)[length(unique(counts_df_period25$id_2))])
males25 <- males[,c(21,6,9,22:53)]
males25 <- males25 %>% dplyr::filter(id %in% ids25)
males25

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

which(males25$degree_0.1 < males25$degree_0.2)
which(males25$degree_0.2 < males25$degree_0.3)
which(males25$degree_0.3 < males25$degree_0.4)
which(males25$degree_0.4 < males25$degree_0.5)

# age variable
males25$age <- lubridate::year(periods$period_start[periods$period == 25]) - males25$byr
summary(males25$age)
males25$age_class <- ifelse(males25$age < 5, 1,
                            ifelse(males25$age < 10, 2,
                                   ifelse(males25$age < 15, 3,
                                          ifelse(males25$age < 20, 4,
                                                 ifelse(males25$age < 25, 5,
                                                        ifelse(males25$age < 40, 6, 7))))))

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
     vertex.label = males25$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males25$age_class == 7,'seagreen4',
                          ifelse(males25$age_class == 6,'seagreen3',
                                 ifelse(males25$age_class == 5,'seagreen2',
                                        ifelse(males25$age_class == 4,'steelblue3',
                                               ifelse(males25$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males25$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males25$age_class == 7,'seagreen4',
                          ifelse(males25$age_class == 6,'seagreen3',
                                 ifelse(males25$age_class == 5,'seagreen2',
                                        ifelse(males25$age_class == 4,'steelblue3',
                                               ifelse(males25$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 25 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males25$id[which(males25$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males25$id[which(males25$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males25[which(males25$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males25[which(males25$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males25[which(males25$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males25[which(males25$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males25[which(males25$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 25 ####
summaries$period <- 25
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period25, draws_anp_25, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males25, plot_data_anp_25, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids25, N, plot_cols)

################ 7.26) Run model on real standardised data -- period 26 ################
### create data list
counts_df_period26 <- counts_df_non0[counts_df_non0$period == 26,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period26),          # total number of times one or other of the dyad was observed
  together = counts_df_period26$event_count,    # count number of sightings seen together
  apart    = counts_df_period26$apart,          # count number of sightings seen apart
  period   = counts_df_period26$period)         # which period it's within

### Fit model
weight_anp_2.2_period26 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period26
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period26$output_files()[1])
draws1_anp_26 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period26$output_files()[2])
draws2_anp_26 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period26$output_files()[3])
draws3_anp_26 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period26$output_files()[4])
draws4_anp_26 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_26 <- rbind(draws1_anp_26, draws2_anp_26, draws3_anp_26, draws4_anp_26)

colnames(draws_anp_26)[2:ncol(draws_anp_26)] <- counts_df_period26$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_26), size = 30, replace = F)

### save data 
write_csv(draws_anp_26, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period26_22.04.27.csv')
rm(draws1_anp_26, draws2_anp_26, draws3_anp_26, draws4_anp_26)

### build traceplots  -- period 26 ####
plot(draws_anp_26[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_26[,plot_cols[2]], col = 'tan')
lines(draws_anp_26[,plot_cols[3]], col = 'orange')
lines(draws_anp_26[,plot_cols[4]], col = 'green')
lines(draws_anp_26[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_26[,plot_cols[6]], col = 'blue')
lines(draws_anp_26[,plot_cols[7]], col = 'red')
lines(draws_anp_26[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_26[,plot_cols[9]], col = 'purple')
lines(draws_anp_26[,plot_cols[10]],col = 'magenta')
lines(draws_anp_26[,plot_cols[11]],col = 'black')
lines(draws_anp_26[,plot_cols[12]], col = 'tan')
lines(draws_anp_26[,plot_cols[13]], col = 'orange')
lines(draws_anp_26[,plot_cols[14]], col = 'green')
lines(draws_anp_26[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_26[,plot_cols[16]], col = 'blue')
lines(draws_anp_26[,plot_cols[17]], col = 'red')
lines(draws_anp_26[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_26[,plot_cols[19]], col = 'purple')
lines(draws_anp_26[,plot_cols[20]],col = 'magenta')
lines(draws_anp_26[,plot_cols[21]],col = 'black')
lines(draws_anp_26[,plot_cols[22]], col = 'tan')
lines(draws_anp_26[,plot_cols[23]], col = 'orange')
lines(draws_anp_26[,plot_cols[24]], col = 'green')
lines(draws_anp_26[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_26[,plot_cols[26]], col = 'blue')
lines(draws_anp_26[,plot_cols[27]], col = 'red')
lines(draws_anp_26[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_26[,plot_cols[29]], col = 'purple')
lines(draws_anp_26[,plot_cols[30]],col = 'magenta')

### density plots  -- period 26 ####
dens(draws_anp_26[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_26[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 26 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_26 <- left_join(x = summaries, y = counts_df_period26, by = 'dyad')
head(plot_data_anp_26)
plot_data_anp_26$age_cat_1 <- ifelse(plot_data_anp_26$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_26$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_26$age_start_1 > 19, 'A','P')))
plot_data_anp_26$age_cat_2 <- ifelse(plot_data_anp_26$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_26$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_26$age_start_2 > 19, 'A','P')))
plot_data_anp_26$age_catnum_1 <- ifelse(plot_data_anp_26$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_26$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_26$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_26$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_26$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_26$age_start_1 < 40, 6, 7))))))
plot_data_anp_26$age_catnum_2 <- ifelse(plot_data_anp_26$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_26$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_26$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_26$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_26$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_26$age_start_2 < 40, 6, 7))))))

plot_data_anp_26$age_dyad <- ifelse(plot_data_anp_26$age_catnum_1 >= plot_data_anp_26$age_catnum_2,
                                    paste(plot_data_anp_26$age_cat_1, plot_data_anp_26$age_cat_2, sep = ''),
                                    paste(plot_data_anp_26$age_cat_2, plot_data_anp_26$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_26, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_26, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_26[plot_data_anp_26$age_dyad == 'AA' | plot_data_anp_26$age_dyad == 'AP' | plot_data_anp_26$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 26 ################
head(summaries)
length(unique(plot_data_anp_26$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids26 <- c(unique(counts_df_period26$id_1), unique(counts_df_period26$id_2)[length(unique(counts_df_period26$id_2))])
males26 <- males[,c(21,6,9,22:53)]
males26 <- males26 %>% dplyr::filter(id %in% ids26)
males26

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

which(males26$degree_0.1 < males26$degree_0.2)
which(males26$degree_0.2 < males26$degree_0.3)
which(males26$degree_0.3 < males26$degree_0.4)
which(males26$degree_0.4 < males26$degree_0.5)

# age variable
males26$age <- lubridate::year(periods$period_start[periods$period == 26]) - males26$byr
summary(males26$age)
males26$age_class <- ifelse(males26$age < 5, 1,
                            ifelse(males26$age < 10, 2,
                                   ifelse(males26$age < 15, 3,
                                          ifelse(males26$age < 20, 4,
                                                 ifelse(males26$age < 25, 5,
                                                        ifelse(males26$age < 40, 6, 7))))))

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
     vertex.label = males26$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males26$age_class == 7,'seagreen4',
                          ifelse(males26$age_class == 6,'seagreen3',
                                 ifelse(males26$age_class == 5,'seagreen2',
                                        ifelse(males26$age_class == 4,'steelblue3',
                                               ifelse(males26$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males26$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males26$age_class == 7,'seagreen4',
                          ifelse(males26$age_class == 6,'seagreen3',
                                 ifelse(males26$age_class == 5,'seagreen2',
                                        ifelse(males26$age_class == 4,'steelblue3',
                                               ifelse(males26$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 26 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males26$id[which(males26$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males26$id[which(males26$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males26[which(males26$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males26[which(males26$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males26[which(males26$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males26[which(males26$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males26[which(males26$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 26 ####
summaries$period <- 26
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period26, draws_anp_26, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males26, plot_data_anp_26, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids26, N, plot_cols)

################ 7.27) Run model on real standardised data -- period 27 ################
### create data list
counts_df_period27 <- counts_df_non0[counts_df_non0$period == 27,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period27),          # total number of times one or other of the dyad was observed
  together = counts_df_period27$event_count,    # count number of sightings seen together
  apart    = counts_df_period27$apart,          # count number of sightings seen apart
  period   = counts_df_period27$period)         # which period it's within

### Fit model
weight_anp_2.2_period27 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period27
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period27$output_files()[1])
draws1_anp_27 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period27$output_files()[2])
draws2_anp_27 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period27$output_files()[3])
draws3_anp_27 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period27$output_files()[4])
draws4_anp_27 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_27 <- rbind(draws1_anp_27, draws2_anp_27, draws3_anp_27, draws4_anp_27)

colnames(draws_anp_27)[2:ncol(draws_anp_27)] <- counts_df_period27$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_27), size = 30, replace = F)

### save data 
write_csv(draws_anp_27, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period27_22.04.27.csv')
rm(draws1_anp_27, draws2_anp_27, draws3_anp_27, draws4_anp_27)

### build traceplots  -- period 27 ####
plot(draws_anp_27[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_27[,plot_cols[2]], col = 'tan')
lines(draws_anp_27[,plot_cols[3]], col = 'orange')
lines(draws_anp_27[,plot_cols[4]], col = 'green')
lines(draws_anp_27[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_27[,plot_cols[6]], col = 'blue')
lines(draws_anp_27[,plot_cols[7]], col = 'red')
lines(draws_anp_27[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_27[,plot_cols[9]], col = 'purple')
lines(draws_anp_27[,plot_cols[10]],col = 'magenta')
lines(draws_anp_27[,plot_cols[11]],col = 'black')
lines(draws_anp_27[,plot_cols[12]], col = 'tan')
lines(draws_anp_27[,plot_cols[13]], col = 'orange')
lines(draws_anp_27[,plot_cols[14]], col = 'green')
lines(draws_anp_27[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_27[,plot_cols[16]], col = 'blue')
lines(draws_anp_27[,plot_cols[17]], col = 'red')
lines(draws_anp_27[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_27[,plot_cols[19]], col = 'purple')
lines(draws_anp_27[,plot_cols[20]],col = 'magenta')
lines(draws_anp_27[,plot_cols[21]],col = 'black')
lines(draws_anp_27[,plot_cols[22]], col = 'tan')
lines(draws_anp_27[,plot_cols[23]], col = 'orange')
lines(draws_anp_27[,plot_cols[24]], col = 'green')
lines(draws_anp_27[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_27[,plot_cols[26]], col = 'blue')
lines(draws_anp_27[,plot_cols[27]], col = 'red')
lines(draws_anp_27[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_27[,plot_cols[29]], col = 'purple')
lines(draws_anp_27[,plot_cols[30]],col = 'magenta')

### density plots  -- period 27 ####
dens(draws_anp_27[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_27[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 27 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_27 <- left_join(x = summaries, y = counts_df_period27, by = 'dyad')
head(plot_data_anp_27)
plot_data_anp_27$age_cat_1 <- ifelse(plot_data_anp_27$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_27$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_27$age_start_1 > 19, 'A','P')))
plot_data_anp_27$age_cat_2 <- ifelse(plot_data_anp_27$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_27$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_27$age_start_2 > 19, 'A','P')))
plot_data_anp_27$age_catnum_1 <- ifelse(plot_data_anp_27$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_27$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_27$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_27$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_27$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_27$age_start_1 < 40, 6, 7))))))
plot_data_anp_27$age_catnum_2 <- ifelse(plot_data_anp_27$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_27$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_27$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_27$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_27$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_27$age_start_2 < 40, 6, 7))))))

plot_data_anp_27$age_dyad <- ifelse(plot_data_anp_27$age_catnum_1 >= plot_data_anp_27$age_catnum_2,
                                    paste(plot_data_anp_27$age_cat_1, plot_data_anp_27$age_cat_2, sep = ''),
                                    paste(plot_data_anp_27$age_cat_2, plot_data_anp_27$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_27, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_27, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_27[plot_data_anp_27$age_dyad == 'AA' | plot_data_anp_27$age_dyad == 'AP' | plot_data_anp_27$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 27 ################
head(summaries)
length(unique(plot_data_anp_27$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids27 <- c(unique(counts_df_period27$id_1), unique(counts_df_period27$id_2)[length(unique(counts_df_period27$id_2))])
males27 <- males[,c(21,6,9,22:53)]
males27 <- males27 %>% dplyr::filter(id %in% ids27)
males27

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

which(males27$degree_0.1 < males27$degree_0.2)
which(males27$degree_0.2 < males27$degree_0.3)
which(males27$degree_0.3 < males27$degree_0.4)
which(males27$degree_0.4 < males27$degree_0.5)

# age variable
males27$age <- lubridate::year(periods$period_start[periods$period == 27]) - males27$byr
summary(males27$age)
males27$age_class <- ifelse(males27$age < 5, 1,
                            ifelse(males27$age < 10, 2,
                                   ifelse(males27$age < 15, 3,
                                          ifelse(males27$age < 20, 4,
                                                 ifelse(males27$age < 25, 5,
                                                        ifelse(males27$age < 40, 6, 7))))))

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
     vertex.label = males27$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males27$age_class == 7,'seagreen4',
                          ifelse(males27$age_class == 6,'seagreen3',
                                 ifelse(males27$age_class == 5,'seagreen2',
                                        ifelse(males27$age_class == 4,'steelblue3',
                                               ifelse(males27$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males27$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males27$age_class == 7,'seagreen4',
                          ifelse(males27$age_class == 6,'seagreen3',
                                 ifelse(males27$age_class == 5,'seagreen2',
                                        ifelse(males27$age_class == 4,'steelblue3',
                                               ifelse(males27$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 -- period 27 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males27$id[which(males27$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males27$id[which(males27$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males27[which(males27$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males27[which(males27$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males27[which(males27$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males27[which(males27$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males27[which(males27$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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

### save summary data -- period 27 ####
summaries$period <- 27
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

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.3,
   counts_df_period27, draws_anp_27, dyad_row, g_mid,g_mid_0.3, g_rng, g_rng_0.3, males27, plot_data_anp_27, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids27, N, plot_cols)

################ 7.28) Run model on real standardised data -- period 28 ################
### create data list
counts_df_period28 <- counts_df_non0[counts_df_non0$period == 28,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period28),          # total number of times one or other of the dyad was observed
  together = counts_df_period28$event_count,    # count number of sightings seen together
  apart    = counts_df_period28$apart,          # count number of sightings seen apart
  period   = counts_df_period28$period)         # which period it's within

### Fit model
weight_anp_2.2_period28 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period28
#variable       mean    median    sd   mad        q5       q95 rhat ess_bulk ess_tail
#
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period28$output_files()[1])
draws1_anp_28 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period28$output_files()[2])
draws2_anp_28 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period28$output_files()[3])
draws3_anp_28 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period28$output_files()[4])
draws4_anp_28 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_28 <- rbind(draws1_anp_28, draws2_anp_28, draws3_anp_28, draws4_anp_28)

colnames(draws_anp_28)[2:ncol(draws_anp_28)] <- counts_df_period28$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_28), size = 30, replace = F)

### save data 
write_csv(draws_anp_28, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period28_22.04.27.csv')
rm(draws1_anp_28, draws2_anp_28, draws3_anp_28, draws4_anp_28)

### build traceplots  -- period 28 ####
plot(draws_anp_28[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_28[,plot_cols[2]], col = 'tan')
lines(draws_anp_28[,plot_cols[3]], col = 'orange')
lines(draws_anp_28[,plot_cols[4]], col = 'green')
lines(draws_anp_28[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_28[,plot_cols[6]], col = 'blue')
lines(draws_anp_28[,plot_cols[7]], col = 'red')
lines(draws_anp_28[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_28[,plot_cols[9]], col = 'purple')
lines(draws_anp_28[,plot_cols[10]],col = 'magenta')
lines(draws_anp_28[,plot_cols[11]],col = 'black')
lines(draws_anp_28[,plot_cols[12]], col = 'tan')
lines(draws_anp_28[,plot_cols[13]], col = 'orange')
lines(draws_anp_28[,plot_cols[14]], col = 'green')
lines(draws_anp_28[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_28[,plot_cols[16]], col = 'blue')
lines(draws_anp_28[,plot_cols[17]], col = 'red')
lines(draws_anp_28[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_28[,plot_cols[19]], col = 'purple')
lines(draws_anp_28[,plot_cols[20]],col = 'magenta')
lines(draws_anp_28[,plot_cols[21]],col = 'black')
lines(draws_anp_28[,plot_cols[22]], col = 'tan')
lines(draws_anp_28[,plot_cols[23]], col = 'orange')
lines(draws_anp_28[,plot_cols[24]], col = 'green')
lines(draws_anp_28[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_28[,plot_cols[26]], col = 'blue')
lines(draws_anp_28[,plot_cols[27]], col = 'red')
lines(draws_anp_28[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_28[,plot_cols[29]], col = 'purple')
lines(draws_anp_28[,plot_cols[30]],col = 'magenta')

### density plots  -- period 28 ####
dens(draws_anp_28[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:length(plot_cols)){
  dens(add = T, draws_anp_28[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 28 ####
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

# organise dem_class -- report age category based on age at start of period
plot_data_anp_28 <- left_join(x = summaries, y = counts_df_period28, by = 'dyad')
head(plot_data_anp_28)
plot_data_anp_28$age_cat_1 <- ifelse(plot_data_anp_28$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_28$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_28$age_start_1 > 19, 'A','P')))
plot_data_anp_28$age_cat_2 <- ifelse(plot_data_anp_28$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_28$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_28$age_start_2 > 19, 'A','P')))
plot_data_anp_28$age_catnum_1 <- ifelse(plot_data_anp_28$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_28$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_28$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_28$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_28$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_28$age_start_1 < 40, 6, 7))))))
plot_data_anp_28$age_catnum_2 <- ifelse(plot_data_anp_28$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_28$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_28$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_28$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_28$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_28$age_start_2 < 40, 6, 7))))))

plot_data_anp_28$age_dyad <- ifelse(plot_data_anp_28$age_catnum_1 >= plot_data_anp_28$age_catnum_2,
                                    paste(plot_data_anp_28$age_cat_1, plot_data_anp_28$age_cat_2, sep = ''),
                                    paste(plot_data_anp_28$age_cat_2, plot_data_anp_28$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_28, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 #fill = c('blue','purple','purple','blue','grey','purple','purple','blue')
    )+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_28, aes(x = age_diff, y = mean))+
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

m_sum <- plot_data_anp_28[plot_data_anp_28$age_dyad == 'AA' | plot_data_anp_28$age_dyad == 'AP' | plot_data_anp_28$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 28 ################
head(summaries)
length(unique(plot_data_anp_28$id_1))+1 # number of individuals = 

par(mai = c(0.1,0.1,0.1,0.1))

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
males28$age_class <- ifelse(males28$age < 5, 1,
                            ifelse(males28$age < 10, 2,
                                   ifelse(males28$age < 15, 3,
                                          ifelse(males28$age < 20, 4,
                                                 ifelse(males28$age < 25, 5,
                                                        ifelse(males28$age < 40, 6, 7))))))

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
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(E(g_mid)$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
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

### only elephants degree > 0.3 -- period 28 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males28$id[which(males28$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males28$id[which(males28$degree_0.3 == 0)])

set.seed(3)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.35),
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
     vertex.color= ifelse(males28[which(males28$degree_0.3 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males28[which(males28$degree_0.3 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males28[which(males28$degree_0.3 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males28[which(males28$degree_0.3 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males28[which(males28$degree_0.3 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.3, add = TRUE)

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
   counts_df_period28, draws_anp_28, dyad_row, g_mid,g_mid_0.2, g_rng, g_rng_0.2, males28, plot_data_anp_28, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids28, N, plot_cols)

#### Save outputs to file #####
# write out csv file
write_csv(dyad_period_weights, 'data_processed/anp_dyad_weightdistributions_2.2_period1to28_22.05.06.csv')

# produce PDF of graphs
dev.off()

# add progress stamp
print(paste0('Time windows 25:28 completed at ', Sys.time()))
