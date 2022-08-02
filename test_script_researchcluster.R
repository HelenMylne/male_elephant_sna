# file to test whether sourcing a script directly from GitHub to research cluster will work

library(tidyverse, lib.loc = 'packages/')
print('tidyverse loaded')

read_csv('data_processed/motnp_id.csv')
