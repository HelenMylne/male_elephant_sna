#### information ####
# script to show how many elephants on average within an equivalent time window (to compare population sizes), and produce descriptive plots for supplementary materials of number of individuals for each window, age distributions and anything that is important to the study.

#### set up ####
library(tidyverse)
theme_set(theme_bw())

#### eles per age category ####
## motnp ####
mot <- read_csv('../data_processed/step1_dataprocessing/motnp_elenodes.csv') %>% 
  filter(sex == 'M') %>% 
  filter(age_class %in% c('Adult','Pubescent')) %>% 
  filter(age_category != '9-10') %>% 
  filter(name != 'Gabriel') %>% 
  filter(name != 'Richard') %>% 
  mutate(age_category = ifelse(age_category == '10-15', '10-14',
                               ifelse(age_category == '25-40', '26-39',
                                      ifelse(age_category == '40+', '≥40', age_category)))) %>% 
  mutate(age_category = factor(age_category, levels = c('10-14','15-19','20-25','26-39','≥40')))

counts <- as.data.frame(table(mot$age_category)) %>% 
  rename(age_category = Var1,
         count = Freq)
percentages <- as.data.frame(prop.table(table(mot$age_category))*100) %>% 
  rename(age_category = Var1,
         percent = Freq) %>% 
  mutate(percent = round(percent, 1)) %>% 
  mutate(percent = paste0(percent, '%')) %>% 
  left_join(counts, by = 'age_category')
ggplot()+
  geom_bar(data = mot,
           aes(x = age_category,
               fill = age_category))+
  scale_y_continuous(name = 'number of elephants',
                     expand = c(0,0),
                     limits = c(0,72),
                     breaks = c(0,10,20,30,40,50,60,70))+
  scale_x_discrete(name = 'age category (years)')+
  scale_fill_viridis_d(direction = -1)+
  labs(fill = 'age category (years)')+
  geom_text(aes(x = age_category, y = 18, label = percent), data = percentages)+
  geom_text(aes(x = age_category, y = 23, label = count), data = percentages)+
  theme(legend.position = 'bottom')
ggsave('motnp_descriptive.png',
       plot = last_plot(),
       device = 'png',
       path = '../outputs/descriptive_stats/',
       width = 1600, height = 1200, units = 'px')

## anp short ####
rm(list = ls()) ; gc()
load('anp_edgecalculations/excluding_under10s/anpshort1_edgeweights_conditionalprior.RData')

anp_dyads <- counts_df
anp1 <- anp_dyads %>%
  dplyr::select(id_1, node_1, period, age_start_1) %>%
  distinct() %>%
  rename(id = id_1, node = node_1, age_start = age_start_1)
anp2 <- anp_dyads %>%
  dplyr::select(id_2, node_2, period, age_start_2) %>%
  distinct() %>%
  rename(id = id_2, node = node_2, age_start = age_start_2)
anp <- rbind(anp1, anp2) %>%
  distinct() %>%
  mutate(age_category = ifelse(age_start < 15, '10-14',
                               ifelse(age_start < 21, '15-19',
                                      ifelse(age_start < 26, '20-25',
                                             ifelse(age_start < 41, '26-39', '≥40'))))) %>%
  mutate(age_category = factor(age_category, levels = c('10-14','15-19','20-25','26-39','≥40')))
rm(list = ls()[! ls() %in% c('anp','anp_dyads')])

counts <- as.data.frame(table(anp$age_category,anp$period)) %>% 
  rename(age_category = Var1,
         period = Var2,
         count = Freq) %>% 
  mutate(period = as.numeric(period))

percentages <- as.data.frame(prop.table(table(anp$age_category,anp$period))*100) %>% 
  rename(age_category = Var1,
         period = Var2,
         percent = Freq) %>% 
  mutate(period = as.numeric(period)) %>% 
  mutate(percent = round(percent, 2)) %>% 
  mutate(percent_long = paste0(percent, '%')) %>% 
  left_join(counts, by = c('age_category','period'))

totals <- counts %>% 
  group_by(period) %>% 
  summarise(n = sum(count)) %>% 
  ungroup()
percentages <- counts %>% 
  left_join(totals, by = 'period') %>% 
  mutate(percent = round((count/n)*100, 1)) %>% 
  mutate(percent_long = paste0(percent,'%'))

percentages2 <- percentages %>% 
  mutate(age_plot = ifelse(age_category == '10-14', 12,
                           ifelse(age_category == '15-19', 18,
                                  ifelse(age_category == '20-25', 25,
                                         ifelse(age_category == '26-39', 33,
                                                50))))) %>%
  mutate(percent_short = paste0(round(percent,0),'%'))

periods <- anp_dyads %>% 
  dplyr::select(period, period_start, period_end) %>% 
  distinct() %>% 
  mutate(start_my = format(period_start, format = "%b %Y"),
         end_my = format(period_end, format = "%b %Y")) %>% 
  mutate(start_end = paste0(start_my, ' - ', end_my))

anp <- anp %>% 
  left_join(periods, by = 'period') %>% 
  mutate(age_cat_mid = ifelse(age_category == '10-14',12,
                              ifelse(age_category == '15-19',17,
                                     ifelse(age_category == '20-25',22.5,
                                            ifelse(age_category == '26-39',32.5,50)))),
         width = ifelse(age_category == '10-14', 5,
                        ifelse(age_category == '15-19', 5,
                               ifelse(age_category == '20-25', 6,
                                      ifelse(age_category == '26-39', 14, 20)))))
percentages2 <- percentages2 %>% 
  left_join(periods[,c('period','end_my','start_end')], by = 'period') %>% 
  mutate(age_cat_mid = ifelse(age_category == '10-14',12,
                              ifelse(age_category == '15-19',17,
                                     ifelse(age_category == '20-25',22.5,
                                            ifelse(age_category == '26-39',32.5,50)))),
         width = ifelse(age_category == '10-14', 5,
                        ifelse(age_category == '15-19', 5,
                               ifelse(age_category == '20-25', 6,
                                      ifelse(age_category == '26-39', 14, 20)))))

(ggplot()+
  geom_bar(data = anp,
           aes(x = age_category,
               # y = count,
               fill = age_category))+
  facet_wrap(. ~ reorder(start_end, period), ncol = 4)+
  scale_y_continuous(name = 'number of elephants',
                     expand = c(0,0),
                     # limits = c(0,72),
                     # breaks = c(0,10,20,30,40,50,60,70)
                     )+
  scale_x_discrete(name = 'age category (years)')+
  scale_fill_viridis_d(direction = -1)+
  labs(fill = 'age category (years)')+
  geom_text(aes(x = age_category, y = 18, label = percent_short),
            data = percentages2, size = 2)+
  geom_text(aes(x = age_category, y = 30, label = count),
            data = percentages2, size = 3)+
  theme(legend.position = 'bottom',
        strip.text = element_text(size = 8)))
ggsave('anp_short_descriptive_bar.png',
       plot = last_plot(),
       device = 'png',
       path = '../outputs/descriptive_stats/',
       width = 4800, height = 2400, units = 'px')

ggplot()+
  geom_histogram(data = anp,
                 aes(x = age_start,
                     fill = age_category),
                 binwidth = 1)+
  facet_wrap(. ~ reorder(start_end, period), ncol = 4)+
  scale_y_continuous(name = 'number of elephants',
                     expand = c(0,0),
                     limits = c(0,25))+
  scale_x_continuous(name = 'age (years)')+
  scale_fill_viridis_d(direction = -1)+
  labs(fill = 'age category (years)')+
  geom_text(aes(x = age_plot, y = 16, label = percent_short),
            data = percentages2, size = 2)+
  geom_text(aes(x = age_plot, y = 20, label = count),
            data = percentages2, size = 3)+
  theme(legend.position = 'bottom',
        strip.text = element_text(size = 8))
ggsave('anp_short_descriptive_hist.png',
       plot = last_plot(),
       device = 'png',
       path = '../outputs/descriptive_stats/',
       width = 2400, height = 3200, units = 'px')

ggplot()+
    geom_histogram(data = anp[anp$period < 19,],
                   aes(x = age_start,
                       fill = age_category),
                   binwidth = 1)+
    facet_wrap(. ~ reorder(start_end, period), ncol = 3)+
    scale_y_continuous(name = 'number of elephants',
                       expand = c(0,0),
                       limits = c(0,25))+
    scale_x_continuous(name = 'age (years)')+
    scale_fill_viridis_d(direction = -1)+
    labs(fill = 'age category (years)')+
    geom_text(data = percentages2[percentages2$period < 19,],
              aes(x = age_plot, y = 16, label = percent_short),
              size = 3)+
    geom_text(data = percentages2[percentages2$period < 19,],
              aes(x = age_plot, y = 20, label = count),
              size = 4)+
    theme(legend.position = 'bottom',
          strip.text = element_text(size = 10))
ggsave('anp_short_descriptive_hist_firsthalf.png',
       plot = last_plot(),
       device = 'png',
       path = '../outputs/descriptive_stats/',
       width = 2400, height = 3200, units = 'px')

ggplot()+
  geom_histogram(data = anp[anp$period > 18,],
                 aes(x = age_start,
                     fill = age_category),
                 binwidth = 1)+
  facet_wrap(. ~ reorder(start_end, period), ncol = 3)+
  scale_y_continuous(name = 'number of elephants',
                     expand = c(0,0),
                     limits = c(0,25))+
  scale_x_continuous(name = 'age (years)')+
  scale_fill_viridis_d(direction = -1)+
  labs(fill = 'age category (years)')+
  geom_text(data = percentages2[percentages2$period > 18,],
            aes(x = age_plot, y = 16, label = percent_short),
            size = 3)+
  geom_text(data = percentages2[percentages2$period > 18,],
            aes(x = age_plot, y = 20, label = count),
            size = 4)+
  theme(legend.position = 'bottom',
        strip.text = element_text(size = 10))
ggsave('anp_short_descriptive_hist_secondhalf.png',
       plot = last_plot(),
       device = 'png',
       path = '../outputs/descriptive_stats/',
       width = 2400, height = 3200, units = 'px')

## anp long ####
rm(list = ls()) ; gc()
load('anp_edgecalculations/excluding_under10s/anplong1_edgeweights_conditionalprior.RData')

anp_dyads <- counts_df
anp1 <- anp_dyads %>% 
  dplyr::select(id_1, node_1, period, age_start_1) %>% 
  distinct() %>% 
  rename(id = id_1, node = node_1, age_start = age_start_1)
anp2 <- anp_dyads %>% 
  dplyr::select(id_2, node_2, period, age_start_2) %>% 
  distinct() %>% 
  rename(id = id_2, node = node_2, age_start = age_start_2)
anp <- rbind(anp1, anp2) %>% 
  distinct() %>% 
  mutate(age_category = ifelse(age_start < 15, '10-14',
                               ifelse(age_start < 21, '15-19',
                                      ifelse(age_start < 26, '20-25',
                                             ifelse(age_start < 41, '26-39', '≥40'))))) %>% 
  mutate(age_category = factor(age_category, levels = c('10-14','15-19','20-25','26-39','≥40')))
rm(list = ls()[! ls() %in% c('anp','anp_dyads')])

counts <- as.data.frame(table(anp$age_category,anp$period)) %>% 
  rename(age_category = Var1,
         period = Var2,
         count = Freq) %>% 
  mutate(period = as.numeric(period))

percentages <- as.data.frame(prop.table(table(anp$age_category,anp$period))*100) %>% 
  rename(age_category = Var1,
         period = Var2,
         percent = Freq) %>% 
  mutate(period = as.numeric(period)) %>% 
  mutate(percent = round(percent, 2)) %>% 
  mutate(percent_long = paste0(percent, '%')) %>% 
  left_join(counts, by = c('age_category','period'))

totals <- counts %>% 
  group_by(period) %>% 
  summarise(n = sum(count)) %>% 
  ungroup()
percentages <- counts %>% 
  left_join(totals, by = 'period') %>% 
  mutate(percent = round((count/n)*100, 1)) %>% 
  mutate(percent_long = paste0(percent,'%'),
         percent_short = paste0(round(percent,0),'%'),
         age_plot = ifelse(age_category == '10-14', 12,
                           ifelse(age_category == '15-19', 18,
                                  ifelse(age_category == '20-25', 25,
                                         ifelse(age_category == '26-39', 33,
                                                50)))))
periods <- anp_dyads %>% 
  dplyr::select(period, period_start, period_end) %>% 
  distinct() %>% 
  mutate(start_my = format(period_start, format = "%b %Y"),
         end_my = format(period_end, format = "%b %Y")) %>% 
  mutate(start_end = paste0(start_my, ' - ', end_my))
anp <- anp %>% 
  left_join(periods, by = 'period') %>% 
  mutate(age_cat_mid = ifelse(age_category == '10-14',12,
                              ifelse(age_category == '15-19',17,
                                     ifelse(age_category == '20-25',22.5,
                                            ifelse(age_category == '26-39',32.5,50)))),
         width = ifelse(age_category == '10-14', 5,
                        ifelse(age_category == '15-19', 5,
                               ifelse(age_category == '20-25', 6,
                                      ifelse(age_category == '26-39', 14, 20)))))
percentages <- percentages %>% 
  left_join(periods[,c('period','end_my','start_end')], by = 'period') %>% 
  mutate(age_cat_mid = ifelse(age_category == '10-14',12,
                              ifelse(age_category == '15-19',17,
                                     ifelse(age_category == '20-25',22.5,
                                            ifelse(age_category == '26-39',32.5,50)))),
         width = ifelse(age_category == '10-14', 5,
                        ifelse(age_category == '15-19', 5,
                               ifelse(age_category == '20-25', 6,
                                      ifelse(age_category == '26-39', 14, 20))))) %>% 
  mutate(age_plot = ifelse(age_category == '10-14', 12,
                           ifelse(age_category == '15-19', 17,
                                  ifelse(age_category == '20-25', 22.5,
                                         ifelse(age_category == '26-39', 32.5,
                                                50))))) %>%
  mutate(percent_short = paste0(round(percent,0),'%'))

ggplot()+
  geom_bar(data = anp,
           aes(x = age_category,
               # y = count,
               fill = age_category))+
  facet_wrap(. ~ period)+
  scale_y_continuous(name = 'number of elephants',
                     expand = c(0,0),
                     # limits = c(0,72),
                     # breaks = c(0,10,20,30,40,50,60,70)
  )+
  scale_x_discrete(name = 'age category (years)')+
  scale_fill_viridis_d(direction = -1)+
  labs(fill = 'age category (years)')+
  geom_text(aes(x = age_category, y = 18, label = percent_long),
            data = percentages, size = 3)+
  geom_text(aes(x = age_category, y = 30, label = count),
            data = percentages, size = 4)+
  theme(legend.position = 'bottom')
ggsave('anp_long_descriptive_bar.png',
       plot = last_plot(),
       device = 'png',
       path = '../outputs/descriptive_stats/',
       width = 4800, height = 2400, units = 'px')

ggplot()+
  geom_histogram(data = anp,
                 aes(x = age_start,
                     fill = age_category),
                 binwidth = 1)+
  geom_text(aes(x = age_plot, y = 25, label = percent_short),
            data = percentages, size = 3)+
  geom_text(aes(x = age_plot, y = 30, label = count),
            data = percentages, size = 4)+
  facet_wrap(. ~ reorder(start_end, period), ncol = 3)+
  scale_y_continuous(name = 'number of elephants',
                     expand = c(0,0),
                     limits = c(0,35))+
  scale_x_continuous(name = 'age (years)')+
  scale_fill_viridis_d(direction = -1)+
  labs(fill = 'age category (years)')+
  theme(legend.position = 'bottom',
        strip.text = element_text(size = 10))
ggsave('anp_long_descriptive_hist.png',
       plot = last_plot(),
       device = 'png',
       path = '../outputs/descriptive_stats/',
       width = 2400, height = 1750, units = 'px')

#### eles total ####
rm(list = ls()) ; gc()

## motnp ####
mot <- read_csv('../data_processed/step1_dataprocessing/motnp_elenodes.csv') %>% 
  filter(sex == 'M') %>% 
  filter(age_class %in% c('Adult','Pubescent')) %>% 
  filter(age_category != '9-10') %>% 
  filter(name != 'Gabriel') %>% 
  filter(name != 'Richard') %>% 
  mutate(age_category = ifelse(age_category == '10-15', '10-14',
                               ifelse(age_category == '25-40', '26-39',
                                      ifelse(age_category == '40+', '≥40', age_category)))) %>% 
  mutate(age_category = factor(age_category, levels = c('10-14','15-19','20-25','26-39','≥40')))
length(unique(mot$id_no))

## anp short ####
load('anp_edgecalculations/excluding_under10s/anpshort1_edgeweights_conditionalprior.RData')
anp1 <- counts_df %>%
  dplyr::select(id_1, node_1, period, age_start_1) %>%
  distinct() %>%
  rename(id = id_1, node = node_1, age_start = age_start_1)
anp2 <- counts_df %>%
  dplyr::select(id_2, node_2, period, age_start_2) %>%
  distinct() %>%
  rename(id = id_2, node = node_2, age_start = age_start_2)
anp_nodes_short <- rbind(anp1, anp2) %>%
  distinct() %>%
  mutate(age_category = ifelse(age_start < 15, '10-14',
                               ifelse(age_start < 21, '15-19',
                                      ifelse(age_start < 26, '20-25',
                                             ifelse(age_start < 41, '26-39', '≥40'))))) %>%
  mutate(age_category = factor(age_category, levels = c('10-14','15-19','20-25','26-39','≥40')))
rm(list = ls()[! ls() %in% c('anp_nodes_short','mot','counts_df')])

short_nodes <- anp_nodes_short %>% 
  group_by(period) %>% 
  summarise(n = length(unique(id)))

short_dyads <- counts_df %>% 
  group_by(period) %>% 
  summarise(n = length(unique(dyad_id)))

## anp long ####
load('anp_edgecalculations/excluding_under10s/anplong1_edgeweights_conditionalprior.RData')
anp1 <- counts_df %>% 
  dplyr::select(id_1, node_1, period, age_start_1) %>% 
  distinct() %>% 
  rename(id = id_1, node = node_1, age_start = age_start_1)
anp2 <- counts_df %>% 
  dplyr::select(id_2, node_2, period, age_start_2) %>% 
  distinct() %>% 
  rename(id = id_2, node = node_2, age_start = age_start_2)
anp_nodes_long <- rbind(anp1, anp2) %>% 
  distinct() %>% 
  mutate(age_category = ifelse(age_start < 15, '10-14',
                               ifelse(age_start < 21, '15-19',
                                      ifelse(age_start < 26, '20-25',
                                             ifelse(age_start < 41, '26-39', '≥40'))))) %>% 
  mutate(age_category = factor(age_category, levels = c('10-14','15-19','20-25','26-39','≥40')))
rm(list = ls()[! ls() %in% c('short','mot','short_nodes','short_dyads','anp_nodes_long','counts_df')])

long_nodes <- anp_nodes_long %>% 
  group_by(period) %>% 
  summarise(nodes = length(unique(id)))

long_dyads <- counts_df %>% 
  group_by(period) %>% 
  summarise(nodes = length(unique(dyad_id)))

#### sightings per elephant ####
## anp short ####
load('anp_edgecalculations/excluding_under10s/anpshort1_edgeweights_conditionalprior.RData')

anp1 <- counts_df %>% 
  dplyr::select(id_1, node_1, period, period_count_1, age_start_1) %>% 
  rename(id = id_1, node = node_1, count = period_count_1, age = age_start_1) %>% 
  distinct()
anp2 <- counts_df %>% 
  dplyr::select(id_2, node_2, period, period_count_2, age_start_2) %>% 
  rename(id = id_2, node = node_2, count = period_count_2, age = age_start_2) %>% 
  distinct()
anp <- rbind(anp1, anp2) %>% 
  distinct()
rm(list = ls()[! ls() %in% c('anp','counts_df')]) ; gc()

counts <- anp %>% 
  group_by(period) %>% 
  summarise(mid = median(count),
            min = min(count),
            max = max(count))

## anp long ####
rm(list = ls()) ; gc()
load('anp_edgecalculations/excluding_under10s/anplong1_edgeweights_conditionalprior.RData')

anp1 <- counts_df %>% 
  dplyr::select(id_1, node_1, period, count_period_1, age_start_1) %>% 
  rename(id = id_1, node = node_1, count = count_period_1, age = age_start_1) %>% 
  distinct()
anp2 <- counts_df %>% 
  dplyr::select(id_2, node_2, period, count_period_2, age_start_2) %>% 
  rename(id = id_2, node = node_2, count = count_period_2, age = age_start_2) %>% 
  distinct()
anp <- rbind(anp1, anp2) %>% 
  distinct()
rm(list = ls()[! ls() %in% c('anp','counts_df')]) ; gc()

counts <- anp %>% 
  group_by(period) %>% 
  summarise(mid = median(count),
            min = min(count),
            max = max(count))
