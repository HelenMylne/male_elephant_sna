# script to look at number of sightings alone per individual vs age
#### load packages ####
library(tidyverse)
library(patchwork)

#### import data ####
## import counts_df from edge weights
load('motnp_edgeweights_conditionalprior.RData')
rm(list = ls()[! ls() %in% c('counts_df')])

## load original nodes data
nodes <- read_csv('../data_processed/step1_dataprocessing/motnp_elenodes.csv') %>% 
  filter(dem_class == 'AM' | dem_class == 'PM') %>% 
  filter(age_category != '9-10') %>% 
  filter(name != 'Richard') %>% 
  filter(name != 'Gabriel')

## load original sightings data
sightings <- read_csv('../data_processed/step1_dataprocessing/motnp_eles_long.csv') %>% 
  select (-number) %>% 
  filter(elephant %in% nodes$id)

#### calculate basic values to report ####
## count number of sightings of different types
length(unique(sightings$encounter))                                  # 481
length(which(sightings$total_id_hkm == 1))                           # 120
length(which(sightings$total_elephants_numeric == 1))                # 89
length(unique(sightings$encounter[which(sightings$total_elephants_numeric > 1 &
                                          sightings$type == 'MO')])) # 164
length(unique(sightings$encounter[which(sightings$total_elephants_numeric > 1 &
                                          sightings$type == 'MX')])) # 198
length(unique(sightings$encounter[which(sightings$total_elephants_numeric > 1 &
                                          sightings$type == 'BH')])) # 29
length(unique(sightings$encounter[which(sightings$total_elephants_numeric > 1 &
                                          sightings$type == 'UK')])) # 5

## correct encounter 50 which is listed as 1 elephant but 5 are identified (and 8 were originally but 3 deleted as <10 yrs)
solo <- sightings %>% 
  filter(total_elephants_numeric == 1)
table(solo$encounter) # 50 has 5 sightings
sightings <- sightings %>% 
  mutate(total_elephants_numeric = ifelse(encounter == 50, 5, total_elephants_numeric)) %>% 
  mutate(perc_id_hkm = ifelse(encounter == 50, NA, perc_id_hkm))

#### look at which elephants were seen alone ####
solo <- sightings %>% 
  filter(total_elephants_numeric == 1)
length(unique(solo$elephant)) # 41 / 84 are unique
table(solo$elephant)

## look at ages of males seen alone
solo %>% 
  rename(id = elephant) %>% 
  left_join(nodes, by = 'id') %>% 
  mutate(age_category = factor(age_category,
                               levels = c('10-15','15-19','20-25',
                                          '25-40','40+'))) %>% 
  ggplot(aes(x = reorder(id, -as.numeric(age_category)),
             fill = age_category))+
  geom_bar()+
  scale_fill_viridis_d()+
  coord_flip()+
  labs(colour = 'age category',
       x = 'individual ID',
       y = 'sightings alone')+
  theme_bw()+
  theme(legend.position = c(0.8,0.84))

## compare to ages to total population
alone <- solo %>% 
  rename(id = elephant) %>% 
  left_join(nodes, by = 'id') %>% 
  mutate(age_category = factor(age_category,
                               levels = c('10-15','15-19','20-25',
                                          '25-40','40+'))) %>% 
  ggplot(aes(x = age_category,
             fill = age_category))+
  geom_bar()+
  scale_fill_viridis_d(begin = 0.2)+
  scale_x_discrete(drop = F)+
  labs(fill = 'age category',
       x = 'age category',
       y = 'sightings of lone individuals')+
  theme_bw()+
  theme(legend.position = 'none')
total <- sightings %>% 
  rename(id = elephant) %>% 
  left_join(nodes, by = 'id') %>% 
  mutate(age_category = factor(age_category,
                               levels = c('10-15','15-19','20-25',
                                          '25-40','40+'))) %>% 
  ggplot(aes(x = age_category,
             fill = age_category))+
  geom_bar()+
  scale_fill_viridis_d()+
  labs(fill = 'age category',
       x = 'age category',
       y = 'total sightings')+
  theme_bw()
(alone + total)+
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides = 'collect')

#### count sightings alone ####
nodes$count_alone <- NA
for(i in 1:nrow(nodes)){
  nodes$count_alone[i] <- length(which(solo$elephant == nodes$id[i]))
}

### plot individual sightings alone vs age
ggplot(nodes, aes(x = count, y = count_alone,
                  colour = age_category))+
  geom_point()+
  scale_colour_viridis_d()+
  labs(colour = 'age category',
       x = 'total sightings per elephant',
       y = 'sightings alone per elephant')+
  theme_bw()+
  geom_abline(linetype = 2)
nodes$prop_alone <- nodes$count_alone / nodes$count

ggplot(nodes, aes(x = age_category, y = prop_alone,
                  fill = age_category))+
  geom_violin(alpha = 0.6)+
  geom_jitter(pch = 21)+
  scale_fill_viridis_d()+
  labs(fill = 'age category',
       x = 'age category',
       y = 'proportion of sightings alone per elephant')+
  theme_bw()

#### repeat from counts_df instead of raw ####
length(unique(c(counts_df$id_1, counts_df$id_2)))

if(counts_df$count_2[nrow(counts_df)] != 1){
  eles <- counts_df %>% 
    dplyr::select(id_1, node_1, node_1_males, age_category_1, count_1, age_cat_id_1) %>% 
    distinct() %>% 
    filter(count_1 == 1)
} else {
  print('id2 only seen once -- COME BACK TO THIS BIT AND REDO IT')
}

colnames(eles) <- c('id','node','node_males','age_category','count','age_cat_id')

eles$seen_alone <- NA
for(i in 1:nrow(eles)){
  x <- counts_df %>% 
    filter(id_1 == eles$id[i] | id_2 == eles$id[i])
  eles$seen_alone[i] <- ifelse(sum(x$event_count) == 0, 'yes','no')
}
table(eles$seen_alone)

#### compare outputs -- different values for different methods because sightings data also includes female counts (total_elephants_numeric != total males) ####
nodes$id[nodes$count_alone == 1 & nodes$count == 1] # "M195" "M201" "M214"
eles$id[eles$seen_alone == 'yes']     # "M128" "M169" "M195" "M201" "M214"

nodes$count[nodes$id == 'M128']       # 1
nodes$count[nodes$id == 'M169']       # 1
nodes$count_alone[nodes$id == 'M128'] # 0
nodes$count_alone[nodes$id == 'M169'] # 0

encounter128 <- sightings$encounter[sightings$elephant == 'M128']
encounter169 <- sightings$encounter[sightings$elephant == 'M169']

sightings$total_elephants_numeric[sightings$encounter == encounter128] # 30
sightings$total_elephants_numeric[sightings$encounter == encounter169] # 5
sightings$total_id_hkm[sightings$encounter == encounter128]            # 7
sightings$total_id_hkm[sightings$encounter == encounter169]            # 2

all_sightings <- read_csv('../data_processed/step1_dataprocessing/motnp_eles_long.csv')
all_sightings$elephant[all_sightings$encounter == encounter128]
all_sightings$elephant[all_sightings$encounter == encounter169]










#
###############
true_ages <- readRDS('data_processed/motnp_ageestimates_mcmcoutput_22.07.13.rds')
true_ages <- true_ages[, which(colnames(true_ages) %in% nodes$id)]

age_long <- pivot_longer(data = true_ages, cols = everything(),
                         names_to = 'elephant', values_to = 'age_draw')

true_ages <- as.matrix(true_ages)
age_mean <- data.frame(elephant = colnames(true_ages),
                       mean_age = NA,
                       stdv_age = NA)
for(i in 1:nrow(age_mean)){
  age_mean$mean_age[i] <- mean(true_ages[,i])
  age_mean$stdv_age[i] <- sd(true_ages[,i])
}

#ggplot(age_long, aes(x = elephant, y = age_draw))+
#  geom_point(shape = 1, size = 2, colour = rgb(0,0,1,0.2))

age_long_solo <- age_long[age_long$elephant %in% solo$elephant,]
ggplot(age_long_solo, aes(x = elephant, y = age_draw))+
  geom_point()

# ages vs group size
s <- sightings[,c('encounter','elephant','total_elephants_numeric')]
age_gs_full <- full_join(age_long, s, by = 'elephant')
age_gs_mean <- full_join(age_mean, s, by = 'elephant')
age_gs_mean$lwr_age <- age_gs_mean$mean_age - age_gs_mean$stdv_age
age_gs_mean$upr_age <- age_gs_mean$mean_age + age_gs_mean$stdv_age
age_gs_mean$max_gs <- NA
age_gs_mean$mean_gs <- NA
age_gs_mean$stdv_gs <- NA
for(i in 1:nrow(age_gs_mean)){
  ele <- age_gs_mean[age_gs_mean$elephant == age_gs_mean$elephant[i],]
  age_gs_mean$max_gs[i] <- ifelse(age_gs_mean$total_elephants_numeric[i] == max(ele$total_elephants_numeric), 
                                  'max','not_max')
  age_gs_mean$mean_gs[i] <- mean(ele$total_elephants_numeric, na.rm = T)
  age_gs_mean$stdv_gs[i] <- sd(ele$total_elephants_numeric, na.rm = T)
  rm(ele)
}

plot(age_gs_mean$total_elephants_numeric ~ age_gs_mean$mean_age,
     pch = 19, col = rgb(0,0,1,0.3), las = 1,
     xlab = 'mean age estimate',
     ylab = 'total elephants in group')
for(i in 1:nrow(age_gs_mean)){
  lines(x = c(age_gs_mean$lwr_age[i], age_gs_mean$upr_age[i]),
        y = c(age_gs_mean$total_elephants_numeric[i], age_gs_mean$total_elephants_numeric[i]),
        col = rgb(0,0,1,0.3))
}

ggplot(age_gs_mean[age_gs_mean$max_gs == 'max',], aes(x = mean_age, y = total_elephants_numeric))+
  geom_point(colour = rgb(0,0,1,0.3))

plot(age_gs_mean$mean_gs ~ jitter(x = age_gs_mean$mean_age, amount = 2),
     pch = 19, col = rgb(0,0,1,0.3), las = 1,
     xlab = 'mean age estimate',
     ylab = 'total elephants in group')
for(i in 1:nrow(age_gs_mean)){
  lines(x = c(age_gs_mean$lwr_age[i], age_gs_mean$upr_age[i]),
        y = c(age_gs_mean$total_elephants_numeric[i], age_gs_mean$total_elephants_numeric[i]),
        col = rgb(0,0,1,0.3))
}

plot_means <- age_gs_mean[,c('elephant','mean_age','stdv_age','mean_gs','stdv_gs')] %>% distinct()
ggplot(plot_means)+
  geom_errorbar(aes(x = jitter(mean_age, amount = 2), ymin = mean_gs - stdv_gs, ymax = mean_gs + stdv_gs),
                colour = rgb(0.5,0,1,0.5))+
  labs(x = 'mean age estimate (jittered by 2)', y = 'mean group size')+
  theme_classic()

plot_means$age_cat <- ifelse(plot_means$mean_age < 10, 1,
                             ifelse(plot_means$mean_age < 15, 2,
                                    ifelse(plot_means$mean_age < 20, 3,
                                           ifelse(plot_means$mean_age < 25, 4,
                                                  ifelse(plot_means$mean_age < 40, 5, 6)))))
plot_means$index <- NA
for(i in 1:nrow(plot_means)){
  ele <- plot_means[plot_means$age_cat == plot_means$age_cat[i],]
  plot_means$index[i] <- ifelse(is.na(ele$index[1]) == TRUE, 1, max(ele$index, na.rm = T) + 1)
  rm(ele)
}
ggplot(plot_means)+
  geom_point(aes(x = index, y = mean_gs), colour = rgb(0.5,0,1,0.5))+
  geom_errorbar(aes(x = index, ymin = mean_gs - stdv_gs, ymax = mean_gs + stdv_gs),
                colour = rgb(0.5,0,1,0.5))+
  labs(x = 'individual', y = 'mean group size')+
  theme_classic()+
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+
  facet_wrap(. ~ age_cat, scales = 'free_x')











