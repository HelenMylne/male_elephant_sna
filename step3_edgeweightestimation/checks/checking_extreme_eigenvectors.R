library(tidyverse)
eles <- read_csv('../data_processed/motnp_eles_long.csv')

## extreme eigenvector values from SRI
m49 <- eles[eles$elephant == 'M49',]
m108 <- eles[eles$elephant == 'M108',]
m153 <- eles[eles$elephant == 'M153',]
m242 <- eles[eles$elephant == 'M242',]
m244 <- eles[eles$elephant == 'M244',]

## extreme eigenvector values in bisonR
m23 <- eles[eles$elephant == 'M23',]
m47 <- eles[eles$elephant == 'M47',]
m145 <- eles[eles$elephant == 'M145',]
m188 <- eles[eles$elephant == 'M188',]
m205 <- eles[eles$elephant == 'M205',]

## combine
xtm_sri <- rbind(m49, m108, m153, m242, m244)
xtm_mod <- rbind(m23, m47,  m145, m188, m205)
xtm_sri$eigen_type <- 'sri'
xtm_mod$eigen_type <- 'bison_model'
xtm <- rbind(xtm_sri, xtm_mod)

## group size
summary(xtm_sri$total_elephants_numeric)
summary(xtm_mod$total_elephants_numeric)
boxplot(xtm$total_elephants_numeric ~ xtm$eigen_type) # none of the SRI extremes seen in the mega group, but overall much larger groups
ggplot(xtm, aes(x = eigen_type, y = total_elephants_numeric, fill = eigen_type))+
  geom_boxplot()+
  geom_jitter(width = 0.3)+
  theme_classic()+
  scale_x_discrete(name = 'model producing extreme eigenvector score')+
  scale_y_continuous(name = 'group size')

summary(m23$total_elephants_numeric)
summary(m47$total_elephants_numeric)
summary(m145$total_elephants_numeric)
summary(m188$total_elephants_numeric)
summary(m205$total_elephants_numeric)

## sighting count
table(xtm_sri$elephant)
table(xtm_mod$elephant)
barplot(table(xtm$elephant, xtm$eigen_type), beside = T) # 
ggplot(xtm, aes(x = elephant, fill = eigen_type))+
  geom_bar()

## check encounter 203
load('motnp_bisonr_edgescalculated_strongprior.RData')
motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds')
s203 <- eles[eles$encounter == 203,] %>% 
  filter(elephant %in% counts_df$id_1 | elephant %in% counts_df$id_2)

m108_cdf <- counts_df[counts_df$id_1 == 'M108' | counts_df$id_2 == 'M108',] %>% 
  filter(event_count > 0) %>% 
  mutate(sri = event_count / count_dyad)
summary(m108_cdf$sri)










