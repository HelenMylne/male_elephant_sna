### summarise data from different sources

#### MOTNP ####
summary <- read_csv('data_processed/motnp_6thJan2022/motnp_elenodes_22.01.06.csv')
summary <- summary[summary$id != 'F157' & summary$id != 'F158' & summary$id != 'F176' &
                     summary$id != 'M125' & summary$id != 'M13' & summary$id != 'M138' &
                     summary$id != 'M21' & summary$id != 'M223' & summary$id != 'M227',] # remove all individuals that can't be in analysis due to identity queries
summary$age_cat_id <- ifelse(summary$age_category == '0-3', 1, 
                             ifelse(summary$age_category == '1-2', 1,
                               ifelse(summary$age_category == '3-4', 1,
                                 ifelse(summary$age_category == '4-5', 1,
                                   ifelse(summary$age_category == '5-6', 2,
                                     ifelse(summary$age_category == '6-7', 2,
                                       ifelse(summary$age_category == '7-8', 2,
                                         ifelse(summary$age_category == '8-9', 2,
                                           ifelse(summary$age_category == '9-10', 2,
                                             ifelse(summary$age_category == '10-15', 3,
                                               ifelse(summary$age_category == '15-19', 4,
                                                 ifelse(summary$age_category == '20-25', 5,
                                                   ifelse(summary$age_category == '20-35', 5,
                                                     ifelse(summary$age_category == '25-40', 6,
                                                       ifelse(summary$age_category == '35-50', 6,
                                                         ifelse(summary$age_category == '40+', 7,
                                                           ifelse(summary$age_category == '50+', 7,
                                                                  summary$age_category)))))))))))))))))
table(summary$sex, summary$age_category)
table(summary$sex, summary$age_cat_id)
independents <- summary[summary$age_cat_id > 2,]
prop.table(table(independents$sex))*100

m <- independents[independents$sex == 'M',] ; m <- m[!is.na(m$id_no),]
f <- independents[independents$sex == 'F',] ; f <- f[!is.na(f$id_no),]

mean(m$count)
mean(f$count)
sd(m$count)
sd(f$count)
min(m$count)
min(f$count)
max(m$count)
max(f$count)

eles_long <- read_csv('data_processed/motnp_6thJan2022/motnp_eles_long_22.01.06.csv')
gs <- eles_long$total_elephants_numeric
which(is.na(gs))
mean(gs)
median(gs)
sd(gs)
min(gs)
max(gs)

obs <- read_delim('data_processed/motnp_6thJan2022/motnp_recording_sessions.csv', delim = ';')
sum(obs$mins)/60
min(obs$date)
max(obs$date)
max(obs$date) - min(obs$date) # 581 days

sightings <- read_csv('data_processed/motnp_encounters_22.01.13.csv')
head(sightings)
summary(sightings$total_elephants_numeric, na.rm = T)
quantile(sightings$total_elephants_numeric, seq(0,1,length.out = 101), na.rm = T)

# edge weight
draws_motnp1.1 <- read_csv('data_processed/motnp_bayesian_edgedistributions_a1.b1_22.03.03.csv') %>%
  data.matrix()
summaries <- data.frame(dyad = colnames(draws_motnp1.1[,2:106954]),
                        min = rep(NA, ncol(draws_motnp1.1)-1),
                        max = rep(NA, ncol(draws_motnp1.1)-1),
                        mean = rep(NA, ncol(draws_motnp1.1)-1),
                        median = rep(NA, ncol(draws_motnp1.1)-1),
                        sd = rep(NA, ncol(draws_motnp1.1)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_motnp1.1[,i+1])
  summaries$max[i]    <- max(draws_motnp1.1[,i+1])
  summaries$mean[i]   <- mean(draws_motnp1.1[,i+1])
  summaries$median[i] <- median(draws_motnp1.1[,i+1])
  summaries$sd[i]     <- sd(draws_motnp1.1[,i+1])
}

# nodes data
ele_nodes <- read_csv('data_processed/motnp_elenodes_22.01.13.csv')
nodes <- data.frame(id = sort(unique(ele_nodes$id)))
nodes <- left_join(nodes, ele_nodes, by = 'id')
nodes$sex       <- as.factor(nodes$sex)
nodes$age_class <- as.factor(nodes$age_class)
nodes$dem_class <- as.factor(nodes$dem_class)
str(nodes)

# correct age and demographic classes for nodes
unique(nodes$age_category) # "50+"   "10-15" "35-50" "20-35" "15-19" "8-9"   "9-10"  "4-5"   "5-6"   "6-7"   "0-3"   "7-8"   "20-25" "25-40" "40+"   "3-4" "1-2"   NA
nodes$age_category <- ifelse(nodes$age_category == '1-2','0-3',nodes$age_category)
nodes$age_cat_id <- ifelse(nodes$age_category == '0-3', 1,
                           ifelse(nodes$age_category == '3-4', 1,
                                  ifelse(nodes$age_category == '4-5', 1,
                                         ifelse(nodes$age_category == '5-6', 2,
                                                ifelse(nodes$age_category == '6-7', 2,
                                                       ifelse(nodes$age_category == '7-8', 2,
                                                              ifelse(nodes$age_category == '8-9', 2,
                                                                     ifelse(nodes$age_category == '9-10', 2,
                                                                            ifelse(nodes$age_category == '10-15', 3,
                                                                                   ifelse(nodes$age_category == '15-19', 4,
                                                                                          ifelse(nodes$age_category == '20-25', 5,
                                                                                                 ifelse(nodes$age_category == '20-35', 5,
                                                                                                        ifelse(nodes$age_category == '25-40', 6,
                                                                                                               ifelse(nodes$age_category == '35-50', 6,
                                                                                                                      ifelse(nodes$age_category == '40+', 7,
                                                                                                                             ifelse(nodes$age_category == '50+', 7, nodes$age_category))))))))))))))))
nodes[is.na(nodes$age_category),]   # U8 doesn't have an age but not a problem -- won't be part of the main analysis. Is a calf so age_cat_id = 1
nodes$age_cat_id[which(is.na(nodes$age_cat_id))] <- 1

unique(nodes$age_category[nodes$age_class == 'Calf'])      # shouldn't include any ages over 4-5
unique(nodes$age_category[nodes$age_class == 'Juvenile'])  # shouldn't include any ages under 5-6
unique(nodes$age_category[nodes$age_class == 'Pubescent']) # shouldn't include any ages under 9-10 or over 15-19
unique(nodes$age_category[nodes$age_class == 'Adult'])     # shouldn't include any ages under 20-25

nodes$age_class <- ifelse(nodes$age_cat_id == 1, 'Calf',
                          ifelse(nodes$age_cat_id == 2, 'Juvenile',
                                 ifelse(nodes$age_cat_id > 4, 'Adult','Pubescent')))

nodes$dem_class <- ifelse(nodes$age_class == 'Adult', paste('A',nodes$sex, sep = ''),
                          ifelse(nodes$age_class == 'Pubescent', paste('P',nodes$sex, sep = ''),
                                 ifelse(nodes$age_class == 'Juvenile', paste('J',nodes$sex, sep = ''),
                                        paste('C',nodes$sex, sep = ''))))

rm(ele_nodes)

# create variables for different degrees of node connectedness
nodes$degree_0.1 <- NA
nodes$degree_0.2 <- NA
nodes$degree_0.3 <- NA
nodes$degree_0.4 <- NA
nodes$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(nodes)){
  rows <- summaries[summaries$id_1 == nodes$id[i] | summaries$id_2 == nodes$id[i],]
  nodes$degree_0.1[i] <- length(which(rows$median > 0.1))
  nodes$degree_0.2[i] <- length(which(rows$median > 0.2))
  nodes$degree_0.3[i] <- length(which(rows$median > 0.3))
  nodes$degree_0.4[i] <- length(which(rows$median > 0.4))
  nodes$degree_0.5[i] <- length(which(rows$median > 0.5))
}

which(nodes$degree_0.1 < nodes$degree_0.2)
which(nodes$degree_0.2 < nodes$degree_0.3)
which(nodes$degree_0.3 < nodes$degree_0.4)
which(nodes$degree_0.4 < nodes$degree_0.5)

# create age_sex variable
nodes$age_sex <- paste(nodes$age_cat_id ,nodes$sex, sep = '')
nodes$family[which(nodes$family == "M26 M87 Dead 1 Jul 17 Poached")] <- 'M26 M87'
nodes$family[which(nodes$family == 'n/a')] <- NA
nodes <- separate(nodes, family, into = c('family.1','family.2','family.3'), sep = ' ')

### dyads ####
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

### correct dem_class with corrected age classes
counts_df$dem_class_1 <- ifelse(counts_df$age_class_1 == 'Adult', paste('A',counts_df$sex_1, sep = ''),
                                ifelse(counts_df$age_class_1 == 'Pubescent', paste('P',counts_df$sex_1, sep = ''),
                                       ifelse(counts_df$age_class_1 == 'Juvenile', paste('J',counts_df$sex_1, sep = ''),
                                              paste('C',counts_df$sex_1, sep = ''))))
counts_df$dem_class_2 <- ifelse(counts_df$age_class_2 == 'Adult', paste('A',counts_df$sex_2, sep = ''),
                                ifelse(counts_df$age_class_2 == 'Pubescent', paste('P',counts_df$sex_2, sep = ''),
                                       ifelse(counts_df$age_class_2 == 'Juvenile', paste('J',counts_df$sex_2, sep = ''),
                                              paste('C',counts_df$sex_2, sep = ''))))

### correct dem_type with new dem_class
counts_df$age_class_id_1 <- ifelse(counts_df$age_class_1 == 'Adult',4,
                                   ifelse(counts_df$age_class_1 == 'Pubescent',3,
                                          ifelse(counts_df$age_class_1 == 'Juvenile',2,1)))
counts_df$age_class_id_2 <- ifelse(counts_df$age_class_2 == 'Adult',4,
                                   ifelse(counts_df$age_class_2 == 'Pubescent',3,
                                          ifelse(counts_df$age_class_2 == 'Juvenile',2,1)))
counts_df$dem_type <- ifelse(counts_df$age_class_id_1 > counts_df$age_class_id_2,
                             paste(counts_df$dem_class_1, counts_df$dem_class_2, sep = '_'),
                             ifelse(counts_df$age_class_id_1 < counts_df$age_class_id_2,
                                    paste(counts_df$dem_class_2, counts_df$dem_class_1, sep = '_'),
                                    paste(counts_df$dem_class_1, counts_df$dem_class_2, sep = '_')))
sort(unique(counts_df$dem_type))

### add column for age difference between dyad
counts_df$age_diff <- abs(as.numeric(counts_df$age_cat_id_1) - as.numeric(counts_df$age_cat_id_2))

### add column for total number of sightings per pair
counts_df$count_dyad <- (counts_df$count_1 + counts_df$count_2) - counts_df$all_events  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.

### add column for total number of sightings per pair where they were NOT together
counts_df$apart <- counts_df$count_dyad - counts_df$all_events

# reassign dyad numbers to remove gaps
counts_df$node_1_nogaps <- as.integer(as.factor(counts_df$node_1))
counts_df$node_2_nogaps <- as.integer(as.factor(counts_df$node_2))+1

# add in age_sex variable
age_sex_family <- nodes[,c('id','age_sex','family.1','family.2','family.3')]
colnames(age_sex_family) <- c('id_1','age_sex_1','family_1.1','family_1.2','family_1.3')
counts_df <- left_join(x = counts_df, y = age_sex_family, by = 'id_1')
colnames(age_sex_family) <- c('id_2','age_sex_2','family_2.1','family_2.2','family_2.3')
counts_df <- left_join(x = counts_df, y = age_sex_family, by = 'id_2')
rm(age_sex_family)
counts_df$dem_dyad <- paste(counts_df$age_sex_1, counts_df$age_sex_2, sep = '_')
barplot(table(counts_df$dem_dyad), las = 1, horiz = T) # some have only a very few pairs

# add in variable for family relationships
for(i in 1:nrow(counts_df)){
  counts_df$family_1.1[i] <- ifelse(is.na(counts_df$family_1.1[i]) == TRUE, 'NA', counts_df$family_1.1[i])
  counts_df$family_1.2[i] <- ifelse(is.na(counts_df$family_1.2[i]) == TRUE, 'NA', counts_df$family_1.2[i])
  counts_df$family_1.3[i] <- ifelse(is.na(counts_df$family_1.3[i]) == TRUE, 'NA', counts_df$family_1.3[i])
  counts_df$family_2.1[i] <- ifelse(is.na(counts_df$family_2.1[i]) == TRUE, 'NA', counts_df$family_2.1[i])
  counts_df$family_2.2[i] <- ifelse(is.na(counts_df$family_2.2[i]) == TRUE, 'NA', counts_df$family_2.2[i])
  counts_df$family_2.3[i] <- ifelse(is.na(counts_df$family_2.3[i]) == TRUE, 'NA', counts_df$family_2.2[i])
}

counts_df$family1 <- ifelse(counts_df$family_1.1 == counts_df$id_2, 'family', 
                            ifelse(counts_df$family_2.1 == counts_df$id_1, 'family', 'not_family'))
counts_df$family2 <- ifelse(counts_df$family_1.2 == counts_df$id_2, 'family',  
                            ifelse(counts_df$family_2.2 == counts_df$id_1, 'family', 'not_family'))
counts_df$family3 <- ifelse(counts_df$family_1.3 == counts_df$id_2, 'family', 
                            ifelse(counts_df$family_2.3 == counts_df$id_1, 'family', 'not_family'))

counts_df$family <- ifelse(counts_df$family1 == 'family' | counts_df$family2 == 'family' | counts_df$family3 == 'family',
                           'mother_calf','not_family')
table(counts_df$family)

family <- counts_df[counts_df$family == 'mother_calf',]
unique(family$dem_type)

counts_df$family <- ifelse(counts_df$family == 'mother_calf' & counts_df$dem_type == 'AF_AF','family',
                           ifelse(counts_df$family == 'mother_calf' & counts_df$dem_type == 'AF_PF','family',
                                  ifelse(counts_df$family == 'mother_calf' & counts_df$dem_type == 'AF_PM','family',
                                         ifelse(counts_df$family == 'mother_calf' & counts_df$dem_type == 'JM_CM','family',
                                                ifelse(counts_df$family == 'mother_calf' & counts_df$dem_type == 'CM_CU','family',
                                                       ifelse(counts_df$family == 'mother_calf' & counts_df$dem_type == 'JM_CU','family',
                                                              ifelse(counts_df$family == 'mother_calf' & counts_df$dem_type == 'PM_CM','family',
                                                                     ifelse(counts_df$family == 'mother_calf' & counts_df$dem_type == 'JU_CU','family',counts_df$family))))))))

table(counts_df$family)

# add some summary data
counts_df <- left_join(counts_df, summaries, by = 'dyad')
counts_df <- counts_df[,c(1:6,12:36,39:43)]
colnames(counts_df)[2:3] <- c('id_1','id_2')

head(counts_df)
hist(counts_df$median)
hist(counts_df$mean)
hist(counts_df$sd)

med0.1 <- length(counts_df$median[which(counts_df$median > 0.1)])
100*(med0.1/length(counts_df$median)) # 26.4%

m_ind <- counts_df[counts_df$dem_type == 'AM_AM' | counts_df$dem_type == 'AM_PM' | counts_df$dem_type == 'PM_PM',]
med_m0.3 <- length(m_ind$median[which(m_ind$median > 0.3)])
100*(med_m0.3/length(m_ind$median))
unique(m_ind$sex_1) ; unique(m_ind$sex_2)
unique(m_ind$age_category_1) ; unique(m_ind$age_category_2)

length(which(counts_df$sd[which(counts_df$median > 0.3)] > 0.1)) / length(counts_df$sd[which(counts_df$median > 0.3)])
length(which(counts_df$sd[which(counts_df$median < 0.3)] > 0.1)) / length(counts_df$sd[which(counts_df$median < 0.3)])

length(which(m_ind$sd[which(m_ind$median > 0.3)] > 0.1)) / length(m_ind$sd[which(m_ind$median > 0.3)])
length(which(m_ind$sd[which(m_ind$median < 0.3)] > 0.1)) / length(m_ind$sd[which(m_ind$median < 0.3)])






