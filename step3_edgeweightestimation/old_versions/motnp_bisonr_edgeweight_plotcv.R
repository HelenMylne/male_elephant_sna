#### set up ####
# load packages
library(tidyverse) #, lib.loc = '../packages/') # library(tidyverse)
library(bisonR, lib.loc = '../packages/')    # library(bisonR)
library(asnipe, lib.loc = '../packages/')    # library(asnipe)
library(sna, lib.loc = '../packages/')       # library(sna)
library(raster, lib.loc = '../packages/')    # library(raster)
library(spatsoc, lib.loc = '../packages/')   # library(spatsoc)

# information
sessionInfo()
R.Version()
rstan::stan_version()

# set seed
set.seed(12345)

# set cmdstanr path
#set_cmdstan_path('H:/rlibs/4.2.1/')

## coefficient of variation of edge weights (aka social differentiation) ####
load('motnp_bisonr_edgescalculated_strongprior.RData')

# add pdf output file
pdf(file = '../outputs/motnp_bisonr_edgeweight_strongprior_cv2.pdf')

# extract cv for model
global_cv_strongprior <- extract_metric(motnp_edge_weights_strongpriors, "global_cv")
head(global_cv_strongprior)
hist(global_cv_strongprior)

# calculate SRI for all dyads
counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
summary(counts_df_model$sri)

# calculate CV of SRI for all dyads
raster::cv(counts_df_model$sri)  # very high, but lower than gbi_matrix
#raster::cv(m$sri[m$sri > 0])     # still massive even when I remove the 0s -- zero inflation is real

### create SRI matrix
# generate matrix
N <- length(unique(c(counts_df_model$id_1, counts_df_model$id_2)))
ids <- unique(c(counts_df_model$id_1, counts_df_model$id_2))
m_mat <- diag(nrow = N)
colnames(m_mat) <- ids
rownames(m_mat) <- ids

# populate matrix with SRI values
for( i in 1:N ) {
  for( j in 1:N ) {
    if( i >= j ) {
      m_mat[i,j] <- m_mat[j,i]
    }
    else {
      id1 <- colnames(m_mat)[i]
      id2 <- rownames(m_mat)[j]
      m_mat[i,j] <- counts_df_model$sri[which(counts_df_model$id_1 == id1 & counts_df_model$id_2 == id2)]
    }
  }
}

# calculate CV for matrix
raster::cv(m_mat) # still way too high
which(is.nan(m_mat) == TRUE)
sd(m_mat)/mean(m_mat) # matches raster version -- check it is doing what I think it is doing

# create gbi_matrix
eles <- read_csv(file = '../data_processed/motnp_eles_long.csv')
eles$location <- paste(eles$gps_s, eles$gps_e, sep = '_')              # make single variable for unique locations
eles <- eles[,c(1,16,2,3,17,4,5,14,7,8,10,13)]                         # rearrange variables
nodes_data <- read_csv(file = '../data_processed/motnp_elenodes.csv' ) # read in node data
colnames(nodes_data)
str(nodes_data)
eles_asnipe <- eles[,c(3,4,2,5)]                                       # date, time, elephant, location
eles_asnipe$Date <- as.integer(eles_asnipe$date)                       # make date numeric
eles_asnipe$Date <- 1+eles_asnipe$Date - min(eles_asnipe$Date)         # start from 1, not 1st January 1970
eles_asnipe$Time <- ifelse(eles_asnipe$time > 1, NA, eles_asnipe$time) # time = proportion of day so anything >1 has to be wrong
eles_asnipe$Time <- eles_asnipe$Time*(24*60*60)                        # convert time values to seconds through day
which(is.na(eles_asnipe$Time))                                         # 161 698 1122 1469 1770
eles_asnipe[c(161,698,1122,1469,1770),]                                # all 1 sighting of B7+M44
eles_asnipe <- eles_asnipe[,c(5,6,3,4)]                                # create data frame to produce gbi matrix from
colnames(eles_asnipe) <- c('Date','Time','ID','Location')              # rename variables for get_gbi
eles_asnipe$ID <- as.character(eles_asnipe$ID)                         # correct data type
str(eles_asnipe)
eles_asnipe$d_pad <- str_pad(eles_asnipe$Date, 3, pad = '0')  # 0-pad dates
eles_asnipe$encounter <- paste(eles_asnipe$d_pad, eles_asnipe$Time, eles_asnipe$Location, sep = '_') # unique value for each sighting
eles_asnipe$group <- as.integer(as.factor(eles_asnipe$encounter)) # unique factor for every sighting
max(eles_asnipe$group)                         # 574 -- agrees with number of different sightings for which elephants were identified
eles_dt <- eles_asnipe[,c(3,7)]                # create data table for gbi matrix
eles_dt <- data.table::setDT(eles_dt)          # create data table for gbi matrix
gbi_matrix <- spatsoc::get_gbi(DT = eles_dt, group = 'group', id = 'ID')  # create gbi matrix
gbi_males <- gbi_matrix[,colnames(gbi_matrix) %in% ids]
gbi_males <- gbi_males[rowSums(gbi_males) > 0,] # remove male-only sightings

# set up permutations
N_networks <- 1000000

# run network permutations
random_networks <- asnipe::network_permutation(association_data = gbi_males,
                                               association_matrix = m_mat,
                                               permutations = N_networks)
print('random networks generated')

cv_random_networks <- rep(0,N_networks)
for (i in c(1:N_networks)) {
  net_rand <- random_networks[i,,]
  cv_random_networks[i] <- raster::cv(net_rand)
  if(i %% 1000 == 0) {print(i)}
}
print(paste0('network permutations for entire network completed at ', Sys.time()))

# compare permuted networks to actual network
hist(cv_random_networks, xlim = c(min(cv_random_networks)-50,max(cv_random_networks)+50),
     main = 'permutations for network of all male elephants')
abline(v = raster::cv(m_mat), col = 'red')
text(round(raster::cv(m_mat),3), col = 'red', x = median(cv_random_networks)+50, y = 2000)

plot_cv <- as.data.frame(cv_random_networks)
cv_network <- raster::cv(m_mat)
cv_random_networks <- plot_cv$cv_random_networks
ggplot(data = plot_cv)+
  geom_histogram(aes(cv_random_networks),
                 fill = rgb(253/255, 231/255, 37/255),
                 colour = 'black',
                 bins = 100)+
  theme_classic()+
  scale_x_continuous(name = 'coefficient of variation',
                     expand = c(0,0),
                     limits = c(min(cv_random_networks)-10,
                                cv_network+10))+
  scale_y_continuous(name = 'frequency',
                     expand = c(0,0))+
  geom_vline(xintercept = cv_network, linewidth = 1.5,
             colour = rgb(68/255, 1/255, 84/255))+
  #geom_text(x = cv_network - 50, y = 700,
  #          colour = rgb(68/255, 1/255, 84/255),
  #          label = 'coefficient of/nvariation for/nmeasured network')+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))

### write out outputs for future reference
write_csv(plot_cv, '../data_processed/motnp_networkpermutations_cv_strongprior.csv')

### plot chain output and look for highest values -- random networks indicating only 2% likelihood of non-random model being best = look to see what value of edges are top 2% ####
counts_df_model <- counts_df[,c('node_1_males','node_2_males','event_count','count_dyad','id_1','id_2','dyad_males')]
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','id_1','id_2','dyad_id')
ms_chain <- as.data.frame(motnp_edge_weights_strongpriors$chain)
colnames(ms_chain) <- counts_df_model$dyad_id
ms_chain <- pivot_longer(ms_chain, cols = everything(), names_to = 'dyad_id', values_to = 'draw')
ms_chain$chain_position <- rep(1:length(unique(ms_chain$dyad_id)), each = 4000)
ms_chain$draw <- LaplacesDemon::invlogit(ms_chain$draw)
ms_chain$mean <- NA
hist(ms_chain$draw)
plot(NULL, xlim = c(0,1), ylim = c(0,30), main = 'distribution of edges for male network using strong prior', xlab = 'edge weight', ylab = 'density', las = 1)
for(i in sample(unique(ms_chain$dyad_id), size = 1000, replace = F)){
  x <- ms_chain[ms_chain$dyad_id == i,]
  ms_chain$mean <- ifelse(ms_chain$dyad_id[i] == x$dyad_id[1], mean(x$draw), ms_chain$mean)
  lines(density(x$draw), col = rgb(0,0,1,0.1))
}
lines(density(ms_chain$draw), lwd = 2)

quantile(ms_chain$draw, 0.98)

ms_edgelist <- bisonR::get_edgelist(motnp_edge_weights_strongpriors)
quantile(ms_edgelist$median, 0.98)

counts_df_model$sri <- counts_df_model$event / counts_df_model$duration
quantile(counts_df_model$sri, 0.98)

## clean up ####
### save pdf
dev.off()

### add time marker
print(paste0('strong-priored model completed at ', Sys.time()))

rm(list= ls()[!(ls() %in% c('motnp_edge_weights_strongpriors','motnp_edges_null_strongpriors',
                            'counts_df','counts_df_model',
                            'gbi_males','m_mat',
                            #'random_nodes',
                            'random_networks'))])

### write over saved workspace image -- only keep the things you might need later that take a long time to run
save.image('motnp_bisonr_edgescalculated_strongprior.RData')

