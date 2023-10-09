### information ####
# script to recreate model results produced by Chiyo 2011 and see how they relate to our results
# Chiyo data: all-male sightings from June-December 2005, 2006 and 2007. Only chose elephants with a minimum of 15 sightings in an all-male group and whose frequent associates were sighted at least 15 times. Total sample = 47 individuals (13% of the male population 10 years and older): sightings mean = 45, median = 39, mode = 46, max = 107.

### set up ####
library(tidyverse)

### import data and filter ####
anp <- read_csv('../data_processed/step1_dataprocessing/anp_sightings_rawcombined.csv') %>%
  filter(obs_date > '2005-05-30' & obs_date < '2008-01-01') %>%
  separate(obs_date, sep = '-', into = c('year','month','day'), remove = F) %>%
  mutate(month = as.numeric(month)) %>%
  filter(month > 5)
# anp <- readxl::read_xlsx('../data_raw/Raw_ATE_Sightings_Fishlock220224.xlsx') %>% 
#   janitor::clean_names() %>%
#   filter(obs_date > '2005-06-01' & obs_date < '2008-01-01') %>% 
#   separate(obs_date, sep = '-', into = c('year','month','day'), remove = F) %>% 
#   mutate(month = as.numeric(month)) %>% 
#   filter(month > 5) %>% 
#   mutate(id = paste0('M',casename))

length(unique(anp$casename))

counts <- data.frame(id = unique(anp$id), count = NA)
for( i in 1:nrow(counts)){
  x <- anp %>% filter(id == counts$id[i])
  counts$count[i] <- nrow(x)
}
hist(counts$count)
length(which(counts$count >= 15))












