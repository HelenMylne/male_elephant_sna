# information: script to incorporate additional changes made by Vicki Fishlock to ANP data into original sightings

library(tidyverse)

# script to incorporate and correct data from Amboseli
old <- readxl::read_xlsx('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/Raw_ATE_Sightings_Fishlock220224.xlsx') %>% 
  janitor::clean_names() # read in original data
new <- readxl::read_xlsx('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/Raw_ATE_Corrections_Fishlock220622.xlsx', sheet = 2) %>% 
  janitor::clean_names() # read in data with notes for what needs updating

# create variable with unique value for every individual observation - will use to join changes to original
old$obs_casename <- paste0(old$obs_id, old$casename)
new$obs_casename <- paste0(new$obs_id, new$casename)

# check which ones to delete immediately
table(new$change_to)
delete <- new[new$change_to == 'delete' | new$change_to == 'delete - is UNKID female group',] # rows to delete
delete$obs_id              # check obs_id of rows to delete
old$row_num <- 1:nrow(old) # give all original data a unique row number

# remove data in deleted rows
updated <- anti_join(x = old, y = delete, by = "obs_casename")
rm(delete)
new <- new[new$change_to != 'delete' & new$change_to != 'delete - is UNKID female group',] # remove rows already corrected

# check which ones need alteration
corrections <- new[new$change_to != '261 - as other pop',]
to_correct <- updated[updated$obs_casename %in% corrections$obs_casename,] # identify rows in original data require correction
which(sort(corrections$obs_casename) != sort(to_correct$obs_casename))     # identify where corrections needed
corrections <- corrections[corrections$obs_casename != "48678217",]

# correct "to_correct" data with "corrections"
barplot(table(to_correct$obs_id)) # all are 1 -- no overlap, can just merge by obs_id
corrections <- corrections[,c('obs_id','change_to')]
to_correct <- left_join(to_correct, corrections, by = 'obs_id')
to_correct <- to_correct[,c(1,25,3:24)]
colnames(to_correct)[2] <- 'casename'

# check remaining corrections
which(sort(to_correct$casename) != sort(corrections$change_to))

# join corrected data into original
updated2 <- anti_join(x = updated, y = to_correct, by = 'row_num')
rm(updated)
updated <- rbind(updated2, to_correct)
rm(updated2)
to_correct

# remove corrected data from list of data to be removed
new <- anti_join(new, corrections, by = 'obs_id')
rm(corrections, to_correct)

# additional corrections
updated[updated$obs_id == new$obs_id[1],]
updated$casename[which(updated$obs_id == new$obs_id[1] & updated$casename == new$casename[1])] <- '261'

updated[updated$obs_id == new$obs_id[2],]
updated[updated$obs_date == new$obs_date[2],]
# this sighting does not exist in original data. As it is a single bull, ignore it.

# In addition, the following males from obsID 36407 should have their group details changed to obsID 18084: 229, 230, 231, 251, 259, 263, 335
sort(unique(updated$casename[updated$obs_id == '36407'])) # all others are below 229
updated$obs_id <- ifelse(updated$obs_id == '36407' & as.numeric(updated$casename) > 228,
                         '18084', updated$obs_id)

write_csv(updated, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_sightings_updated.csv')

