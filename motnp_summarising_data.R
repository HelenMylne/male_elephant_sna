### Basic values for report
library(tidyverse)
library(igraph)
### Number of elephants ####
# read in data ####
counts_df <- read_delim('data_processed/motnp_bayesian_allpairwiseevents_splitbygrouptype_22.01.13.csv', delim = ',')

# correct sex_1, which has loaded in as a logical vector not a character/factor
unique(counts_df$sex_1) # FALSE or NA
sex_1 <- data.frame(sex_1 = counts_df$id_1)
sex_1 <- sex_1 %>% separate(sex_1, c("sex", "number"), sep = 1, remove = FALSE) ; unique(sex_1$sex) # F, M, U
counts_df$sex_1 <- as.character(sex_1$sex) ; rm(sex_1)
counts_df$sex_1 <- as.character(counts_df$sex_1)
str(counts_df)  # sex_1 still comes up as logical?

# create variable for age difference
unique(counts_df$age_category_1) # "0-3","1-2","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-15","15-19","20-25","20-35","25-40","35-50","40+","50+",NA 
counts_df$age_category_1 <- ifelse(counts_df$age_category_1 == '1-2','0-3',counts_df$age_category_1)
counts_df$age_cat_id_1 <- ifelse(counts_df$age_category_1 == '0-3', 1,
                                 ifelse(counts_df$age_category_1 == '3-4', 1,
                                        ifelse(counts_df$age_category_1 == '4-5', 1,
                                               ifelse(counts_df$age_category_1 == '5-6', 2,
                                                      ifelse(counts_df$age_category_1 == '6-7', 2,
                                                             ifelse(counts_df$age_category_1 == '7-8', 2,
                                                                    ifelse(counts_df$age_category_1 == '8-9', 2,
                                                                           ifelse(counts_df$age_category_1 == '9-10', 2,
                                                                                  ifelse(counts_df$age_category_1 == '10-15', 3,
                                                                                         ifelse(counts_df$age_category_1 == '15-19', 4,
                                                                                                ifelse(counts_df$age_category_1 == '20-25', 5,
                                                                                                       ifelse(counts_df$age_category_1 == '20-35', 5,
                                                                                                              ifelse(counts_df$age_category_1 == '25-40', 6,
                                                                                                                     ifelse(counts_df$age_category_1 == '35-50', 6,
                                                                                                                            ifelse(counts_df$age_category_1 == '40+', 7,
                                                                                                                                   ifelse(counts_df$age_category_1 == '50+', 7, counts_df$age_category_1))))))))))))))))
counts_df[is.na(counts_df$age_cat_id_1),]   # U8 doesn't have an age but not a problem -- won't be part of the main analysis. Is a calf so age_cat_id = 1
counts_df$age_cat_id_1[which(is.na(counts_df$age_cat_id_1))] <- 1

unique(counts_df$age_category_1[counts_df$age_class_1 == 'Calf'])      # shouldn't include any ages over 4-5
unique(counts_df$age_category_1[counts_df$age_class_1 == 'Juvenile'])  # shouldn't include any ages under 5-6
unique(counts_df$age_category_1[counts_df$age_class_1 == 'Pubescent']) # shouldn't include any ages under 9-10 or over 15-19
unique(counts_df$age_category_1[counts_df$age_class_1 == 'Adult'])     # shouldn't include any ages under 20-25

counts_df$age_class_1 <- ifelse(counts_df$age_cat_id_1 == 1, 'Calf',
                                ifelse(counts_df$age_cat_id_1 == 2, 'Juvenile',
                                       ifelse(counts_df$age_cat_id_1 > 4, 'Adult','Pubescent')))

unique(counts_df$age_category_2) # "0-3","1-2","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-15","15-19","20-25","20-35","25-40","35-50","40+","50+",NA 
counts_df$age_category_2 <- ifelse(counts_df$age_category_2 == '1-2','0-3',counts_df$age_category_2)
counts_df$age_cat_id_2 <- ifelse(counts_df$age_category_2 == '0-3', 1,
                                 ifelse(counts_df$age_category_2 == '3-4', 1,
                                        ifelse(counts_df$age_category_2 == '4-5', 1,
                                               ifelse(counts_df$age_category_2 == '5-6', 2,
                                                      ifelse(counts_df$age_category_2 == '6-7', 2,
                                                             ifelse(counts_df$age_category_2 == '7-8', 2,
                                                                    ifelse(counts_df$age_category_2 == '8-9', 2,
                                                                           ifelse(counts_df$age_category_2 == '9-10', 2,
                                                                                  ifelse(counts_df$age_category_2 == '10-15', 3,
                                                                                         ifelse(counts_df$age_category_2 == '15-19', 4,
                                                                                                ifelse(counts_df$age_category_2 == '20-25', 5,
                                                                                                       ifelse(counts_df$age_category_2 == '20-35', 5,
                                                                                                              ifelse(counts_df$age_category_2 == '25-40', 6,
                                                                                                                     ifelse(counts_df$age_category_2 == '35-50', 6,
                                                                                                                            ifelse(counts_df$age_category_2 == '40+', 7,
                                                                                                                                   ifelse(counts_df$age_category_2 == '50+', 7, counts_df$age_category_2))))))))))))))))
counts_df[is.na(counts_df$age_cat_id_2),]   # U8 doesn't have an age but not a problem -- won't be part of the main analysis. Is a calf so age_cat_id = 1
counts_df$age_cat_id_2[which(is.na(counts_df$age_cat_id_2))] <- 1

counts_df$age_class_2 <- ifelse(counts_df$age_cat_id_2 == 1, 'Calf',
                                ifelse(counts_df$age_cat_id_2 == 2, 'Juvenile',
                                       ifelse(counts_df$age_cat_id_2 > 4, 'Adult','Pubescent')))

unique(counts_df$age_category_2[counts_df$age_class_2 == 'Calf'])      # shouldn't include any ages over 4-5
unique(counts_df$age_category_2[counts_df$age_class_2 == 'Juvenile'])  # shouldn't include any ages under 5-6
unique(counts_df$age_category_2[counts_df$age_class_2 == 'Pubescent']) # shouldn't include any ages under 9-10 or over 15-19
unique(counts_df$age_category_2[counts_df$age_class_2 == 'Adult'])     # shouldn't include any ages under 20-25

### add column for age difference between dyad
counts_df$age_diff <- abs(as.numeric(counts_df$age_cat_id_1) - as.numeric(counts_df$age_cat_id_2))

### correct dem_class with corrected age classes
counts_df$dem_class_1 <- ifelse(counts_df$age_class_1 == 'Adult', paste('A',counts_df$sex_1, sep = ''),
                                ifelse(counts_df$age_class_1 == 'Pubescent', paste('P',counts_df$sex_1, sep = ''),
                                       ifelse(counts_df$age_class_1 == 'Juvenile', paste('J',counts_df$sex_1, sep = ''),
                                              paste('C',counts_df$sex_1, sep = ''))))
counts_df$dem_class_2 <- ifelse(counts_df$age_class_2 == 'Adult', paste('A',counts_df$sex_2, sep = ''),
                                ifelse(counts_df$age_class_2 == 'Pubescent', paste('P',counts_df$sex_2, sep = ''),
                                       ifelse(counts_df$age_class_2 == 'Juvenile', paste('J',counts_df$sex_2, sep = ''),
                                              paste('C',counts_df$sex_2, sep = ''))))

### add column for total number of sightings per pair
counts_df$count_dyad <- (counts_df$count_1 + counts_df$count_2) - counts_df$all_events  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.

### add column for total number of sightings per pair where they were NOT together
counts_df$apart <- counts_df$count_dyad - counts_df$all_events

### ele_nodes -- ages corrected in counts_df so use this for most things
eles <- read_csv('data_processed/motnp_elenodes_22.01.13.csv')
str(eles)
eles[463,]
U9 <- c('U9',463,NA,'Calf','0-3','U',2,'CU',1)

### observations
sightings <- read_csv('data_processed/motnp_eles_long_22.01.13.csv')

# actually produce the values I need ####
nodes <- counts_df[c(2,4,12,14,16,18,20,22,28)]
colnames(nodes) <- c("id", "node", "name", "age_class", "age_category", "sex", "count", "dem_class", "age_cat_id")
nodes <- distinct(nodes)
nodes <- rbind(nodes,U9)

# number of elephants in each demographic group
table(nodes$age_class, nodes$sex)

# sightings per individual
m <- nodes[nodes$dem_class == 'AM' | nodes$dem_class == 'PM',]
str(m)
m$count <- as.numeric(m$count)
summary(m$count)
sd(m$count)
rethinking::dens(m$count)

f <- nodes[nodes$dem_class == 'AF' | nodes$dem_class == 'PF',]
f$count <- as.numeric(f$count)
summary(f$count)
sd(f$count)
rethinking::dens(f$count)

# observations
s <- sightings[,c(1:7,10,13:14)] %>% 
  distinct()
table(s$total_elephants_numeric)
nrow(s)-table(s$total_elephants_numeric)[1] # 485
rethinking::dens(s$total_elephants_numeric)
quantile(s$total_elephants_numeric, 0.7)

loners <- s[s$total_elephants_numeric == 1,]
table(loners$herd_type)
summary(s$total_elephants_numeric)

groups <- s[s$total_elephants_numeric > 1,]
table(groups$herd_type)

u <- sightings[sightings$herd_type == 'UK',]
table(u$elephant, u$encounter)

max(sightings$total_elephants_numeric)

### simulate graphs for planning ####
sim_eles <- data.frame(dyad   = 1:1200,
                       pop    = rep(c('motnp','mpnp','anp'), each = 400),
                       k_pop  = rep(c(3,2,4), each = 400),
                       k_dyad = rep(NA, 1200))
sim_eles$k_dyad[1:400]    <- sample(1:3, 400, prob = c(0.6,0.3,0.1), replace = T)
sim_eles$k_dyad[401:800]  <- sample(1:2, 400, prob = c(0.9,0.1), replace = T)
sim_eles$k_dyad[801:1200] <- sample(1:4, 400, prob = c(0.6,0.2,0.15,0.05), replace = T)

table(sim_eles$k_dyad, sim_eles$pop)

sim_eles$edge <- NA
for(i in 1:nrow(sim_eles)){
  sim_eles$edge[i] <- rbeta(n = 1,
                            shape1 = ifelse(sim_eles$k_dyad[i] == 1, 2,
                                            ifelse(sim_eles$k_dyad[i] == 2, 4,
                                                   ifelse(sim_eles$k_dyad[i] == 3, 8, 16))),
                            shape2 = ifelse(sim_eles$k_dyad[i] == 4, 2,
                                            ifelse(sim_eles$k_dyad[i] == 3, 4,
                                                   ifelse(sim_eles$k_dyad[i] == 2, 8, 16))))
}
sim_eles$k_dyad <- as.character(sim_eles$k_dyad)
ggplot(sim_eles, aes(x = pop, y = edge, fill = k_dyad))+
  geom_boxplot(notch = T)+
  scale_fill_viridis_d()+
  labs(fill = 'k', x = 'population', y = 'edge weight')+
  theme_light()

### 
N <- 1200
sim_eles <- data.frame(id  = 1:N,
                       pop = rep(c('motnp','mpnp','anp'), each = N/3),
                       age = sample(3:7, replace = T, size = N),
                       bet = rep(NA, N),
                       deg = rep(NA, N),
                       eig = rep(NA, N))
for(i in 1:N){
  sim_eles$bet[i] <- ifelse(sim_eles$pop[i] == 'anp',
                            rbeta(1, shape1 = (sim_eles$age[i]-2)*4, shape2 = sim_eles$age[i]),
                            ifelse(sim_eles$pop[i] == 'motnp',
                                   rbeta(1, shape1 = sim_eles$age[i], shape2 = sim_eles$age[i]*2),
                                   rbeta(1, shape1 = (sim_eles$age[i]+6)*2, shape2 = sim_eles$age[i]*2)))
  sim_eles$deg[i] <- ifelse(sim_eles$pop[i] == 'anp',
                            rpois(n = 1, lambda = sim_eles$age[i]*2),
                            ifelse(sim_eles$pop[i] == 'motnp',
                                   rpois(n = 1, lambda = sim_eles$age[i]),
                                   rpois(n = 1, lambda = 3)))
  sim_eles$eig[i] <- ifelse(sim_eles$pop[i] == 'anp',
                            rbeta(1, shape1 = sim_eles$age[i], shape2 = sim_eles$age[i]),
                            ifelse(sim_eles$pop[i] == 'motnp',
                                   rbeta(1, shape1 = sim_eles$age[i]*3, shape2 = sim_eles$age[i]*2),
                                   rbeta(1, shape1 = sim_eles$age[i], shape2 = sim_eles$age[i]*2)))
}

ggplot(sim_eles, aes(x = age, col = pop, y = bet))+
  geom_jitter()+
  geom_smooth()+
  scale_y_continuous(limits = c(0,1))

ggplot(sim_eles, aes(x = age, col = pop, y = deg))+
  geom_jitter()+
  geom_smooth()

ggplot(sim_eles, aes(x = age, col = pop, y = eig))+
  geom_jitter()+
  geom_smooth()+
  scale_y_continuous(limits = c(0,1))

sim_eles2 <- pivot_longer(sim_eles, cols = c(bet, deg, eig),
                          values_to = 'centrality', names_to = 'measure')
sim_eles2$cent <- ifelse(sim_eles2$measure == 'bet','betweenness',
                         ifelse(sim_eles2$measure == 'deg','degree','eigenvector'))
ggplot(sim_eles2, aes(x = age, col = pop, y = centrality))+
  geom_jitter()+
  geom_smooth()+
  facet_wrap(.~cent, scale = 'free')+
  scale_color_viridis_d(name="population",
                     labels=c("ANP","MOTNP","MPNP"))+
  theme_light()
  
### edge weight ####
draws <- read_csv('data_processed/motnp_bayesian_edgedistributions_a2.b2_22.02.07.csv')
draws <- draws[,2:106954]
weight_summary <- data.frame(dyad = colnames(draws),
                             mean = apply(draws, 2, mean),
                             median = apply(draws, 2, median),
                             min = apply(draws, 2, min),
                             max = apply(draws, 2, max),
                             stdev = apply(draws, 2, sd))
median(weight_summary$median)

100*length(which(weight_summary$mean < 0.1)) / length(weight_summary$mean)

#Only XX dyads comprising males of independent age had mean values over 0.3, indicating a potential for spending over 30% of their time together (Fig 1c), but these individuals all had low sighting frequencies (range) so high uncertainty estimations around their mean values (standard deviations).
males <- draws[,which(counts_df$dem_type == 'AM_AM' | counts_df$dem_type == 'AM_PM' | counts_df$dem_type == 'PM_PM')]
male_summary <- data.frame(dyad = colnames(males),
                           mean = apply(males, 2, mean),
                           median = apply(males, 2, median),
                           min = apply(males, 2, min),
                           max = apply(males, 2, max),
                           stdev = apply(males, 2, sd),
                           upr = apply(males, 2, quantile, 0.975),
                           lwr = apply(males, 2, quantile, 0.025))
male_summary$range <- male_summary$upr - male_summary$lwr
100*length(which(male_summary$mean > 0.3)) / length(male_summary$mean)

male_counts <- counts_df[counts_df$dem_type == 'AM_AM' | counts_df$dem_type == 'AM_PM' | counts_df$dem_type == 'PM_PM',]
m0.3 <- male_counts[which(male_summary$mean > 0.3),]
males0.3 <- data.frame(id = sort(unique(c(unique(m0.3$id_1),unique(m0.3$id_2)))))
eles <- read_csv('data_processed/motnp_elenodes_22.01.13.csv')
males0.3 <- left_join(x = males0.3, y = eles, by = 'id')
table(males0.3$count)

summary(male_summary$stdev[which(male_summary$mean > 0.3)])
hist(male_summary$stdev[which(male_summary$mean > 0.3)], breaks = 100)
quantile(male_summary$stdev[which(male_summary$mean > 0.3)], 0.001)

summary(male_summary$stdev[which(male_summary$mean < 0.3)])
quantile(male_summary$stdev[which(male_summary$mean < 0.3)], 0.73)

male_summary <- left_join(male_summary, counts_df, by = 'dyad')
head(male_summary)

unique(male_summary$age_category_1)
unique(male_summary$age_category_2)
male_summary$age_type2 <- ifelse(male_summary$age_cat_id_1 > male_summary$age_cat_id_2,
                                 paste(male_summary$age_cat_id_2, male_summary$age_cat_id_1, sep = '_'),
                                 paste(male_summary$age_cat_id_1, male_summary$age_cat_id_2, sep = '_'))
male_plot_median <- male_summary[,c('id_1','id_2','median','age_type2')]
male_plot_range  <- male_summary[,c('id_1','id_2','range','age_type2')]
colnames(male_plot_median)[3:4] <- c('weight','type')
colnames(male_plot_range)[3:4]  <- c('weight','type')

unique(male_plot_median$type)
male_plot_median <- male_plot_median[male_plot_median$type != '2_2' &
                                       male_plot_median$type != '2_3' &
                                       male_plot_median$type != '2_4' &
                                       male_plot_median$type != '2_5' &
                                       male_plot_median$type != '2_6' &
                                       male_plot_median$type != '2_7',]
male_plot_range <- male_plot_range[male_plot_range$type != '2_2' &
                                       male_plot_range$type != '2_3' &
                                       male_plot_range$type != '2_4' &
                                       male_plot_range$type != '2_5' &
                                       male_plot_range$type != '2_6' &
                                       male_plot_range$type != '2_7',]

### plots for TAP report ####
# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- igraph::graph_from_data_frame(male_plot_median, directed = F)
g_rng <- igraph::graph_from_data_frame(male_plot_range,  directed = F)

male_nodes1 <- male_counts[,c(2,12,14,16,18,20,22)]
male_nodes2 <- male_counts[,c(3,13,15,17,19,21,23)]
colnames(male_nodes1) <- c('id','name','age_class','age_category','sex','count','dem_class')
colnames(male_nodes2) <- c('id','name','age_class','age_category','sex','count','dem_class')
male_nodes <- rbind(male_nodes1, male_nodes2) %>% 
  distinct()
male_nodes <- male_nodes[male_nodes$age_class != 'Juvenile',]

# Plot all
coords <- igraph::layout_nicely(g_mid)
plot(g_mid,
     edge.width = male_plot_median$weight*5,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = 'black',
     layout = coords)
plot(g_mid,
     edge.width = male_plot_range$weight*5,
     edge.color = rgb(0, 0, 0, 0.25), 
     vertex.size = 8,
     vertex.label = male_nodes$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(male_nodes$age_class == 'Adult','seagreen1',
                          ifelse(male_nodes$age_class == 'Pubescent','skyblue','yellow')),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = male_plot_median$weight,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(male_plot_median$weight < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = male_plot_range$weight,
     edge.color = ifelse(male_plot_median$weight < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = male_nodes$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(male_nodes$age_class == 'Adult','seagreen1',
                          ifelse(male_nodes$age_class == 'Pubescent','skyblue','yellow')),
     layout = coords, add = TRUE)

# create variables for different degrees of node connectedness
male_nodes$degree_0.3 <- NA
for(i in 1:NROW(male_nodes)){
  male_rows <- male_plot_median[male_plot_median$id_1 == male_nodes$id[i] |
                                  male_plot_median$id_2 == male_nodes$id[i],]
  male_nodes$degree_0.3[i] <- length(which(male_rows$weight > 0.3))
}

#  plot network with reduced nodes -- only those with degree values > 0.3
g_mid_0.3 <- igraph::delete.vertices(graph = g_mid,
                             v = male_nodes$id[which(male_nodes$degree_0.3 == 0)])
g_rng_0.3 <- igraph::delete.vertices(graph = g_rng,
                             v = male_nodes$id[which(male_nodes$degree_0.3 == 0)])

coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight, 0),
     edge.color = 'black',
     vertex.size = 7,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight, 0),
     vertex.size = 7,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(male_nodes[which(male_nodes$degree_0.3 != 0),]$age_class == 'Adult',
                           'seagreen1',
                           ifelse(male_nodes[which(male_nodes$degree_0.3 != 0),]$age_class == 'Pubescent', 'skyblue','yellow')),
     layout = coords_0.3, add = TRUE)

write_csv(male_nodes, 'Random/gephi_males_22.2.22.csv')
write_csv(male_plot_median, 'Random/gephi_maleweights_22.2.22.csv')
write_csv(male_plot_range, 'Random/gephi_maleuncertainties_22.2.22.csv')

