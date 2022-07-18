## Bayesian analysis of ATE data
#### Information ####
# Script to process association data from Amboseli National Park, Kenya.
# Based on script 22.01.13_ALERT_bayesian.R

# Data collected: 
# Collected by: Dr Phyllis Lee, Dr Joyce Poole, Dr Vicki Fishlock and other staff at ATE during the collection period
# Data supplied by: Dr Phyllis Lee and Dr Vicki Fishlock, 21st January 2022.

#### Set up ####
library(tidyverse)  # data manipulation
library(lubridate)  # sort date columns out
library(zoo)        # sort date columns out
library(asnipe)     # generating networks

################ 1) Create basic data frames from raw data ################
#### Elephants ####
males <- readxl::read_excel("data_raw/Raw_ATE_Males_Lee220121.xlsx", sheet = 1)
colnames(males) <- c('natal_name','case','name','family',
                     'birth_month','birth_year','birth_accuracy',
                     'death_month','death_year','death_accuracy',
                     'death_cause','death_cause_accuracy',
                     'indep_month','indep_year','indep_accuracy_value','indep_accuracy_word',
                     'musth_month','musth_year','musth_known','age_first_musth')
head(males)
