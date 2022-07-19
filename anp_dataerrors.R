sightings <- readxl::read_excel('data_raw/Raw_ATE_Sightings_Fishlock220224.xlsx') %>% janitor::clean_names()
males <- readxl::read_excel('data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()
males <- males[males$casename %in% sightings$casename,]

males$first <- NA ; males$last <- NA
for(i in 1:nrow(males)){
  casename <- sightings[sightings$casename == males$casename[i],]
  males$first[i] <- lubridate::year(casename$obs_date[1])
  males$last[i] <- lubridate::year(casename$obs_date[nrow(casename)])
}

m <- males[,c(2,6,9,21,22)]
m$dyr_living <- ifelse(m$dyr < 100, 3000, m$dyr) # give living elephants impossibly high death year so don't get counted as dying before their last sighting

m$error <- ifelse(m$byr > m$first, 'before_birth',
                  ifelse(m$dyr_living < m$last, 'after_death','good'))
table(m$error)

before <- m[m$error == 'before_birth',]
after <- m[m$error == 'after_death',]

m$all_errors <- ifelse(m$dyr_living < m$first, 'all_after_death','ignore')
table(m$all_errors)
