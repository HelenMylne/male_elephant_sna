# EfA data model prep
### Information ####
# Data collected by Elephants for Africa (2012-2021)
# Data supplied by Dr Kate Evans
### Set up ####
library(tidyverse)
library(rstan)
library(rethinking)
library(cmdstanr)
library(lubridate)

### read in processed data files ####
# check all fine ####
aa <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1.250_22.03.10.csv')
ab <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings251.387_22.03.08.csv')
ac <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings388.400_22.03.08.csv')
ad <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings401.450_22.03.08.csv')
ae <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings451.500_22.03.08.csv')
af <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings501.550_22.03.08.csv')
ag <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings551.600_22.03.08.csv')
ah <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings601.650_22.03.08.csv')
ai <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings651.700_22.03.08.csv')
aj <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings701.750_22.03.08.csv')
ak <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings751.800_22.03.08.csv')
al <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings801.850_22.03.08.csv')
am <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings851.900_22.03.08.csv')
an <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings901.950_22.03.08.csv')
ao <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings951.1000_22.03.08.csv')
ap <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1001.1050_22.03.08.csv')
aq <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1051.1100_22.03.08.csv')
ar <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1101.1150_22.03.08.csv')
as <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1151.1200_22.03.08.csv')
at <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1201.1250_22.03.08.csv')
au <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1251.1300_22.03.08.csv')
av <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1301.1350_22.03.08.csv')
aw <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1351.1400_22.03.08.csv')
ax <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1401.1450_22.03.08.csv')
ay <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1451.1500_22.03.08.csv')
az <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1501.1550_22.03.08.csv')

ba <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1551.1600_22.03.08.csv')
bb <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1601.1650_22.03.08.csv')
bc <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1651.1700_22.03.08.csv')
bd <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1701.1750_22.03.08.csv')
be <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1751.1800_22.03.08.csv')
bf <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1801.1850_22.03.08.csv')
bg <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1851.1900_22.03.08.csv')
bh <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1901.1950_22.03.08.csv')
bi <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1951.2000_22.03.08.csv')
bj <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2001.2050_22.03.08.csv')
bk <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2051.2100_22.03.08.csv')
bl <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2101.2150_22.03.08.csv')
bm <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2151.2200_22.03.08.csv')
bn <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2201.2250_22.03.08.csv')
bo <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2251.2300_22.03.08.csv')
bp <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2301.2350_22.03.08.csv')
bq <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2351.2400_22.03.08.csv')
br <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2401.2450_22.03.08.csv')
bs <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2451.2500_22.03.08.csv')
bt <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2501.2550_22.03.08.csv')
bu <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2551.2600_22.03.08.csv')
bv <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2601.2650_22.03.08.csv')
bw <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2651.2700_22.03.08.csv')
bx <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2701.2750_22.03.08.csv')
by <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2751.2800_22.03.08.csv')
bz <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2801.2850_22.03.08.csv')

ca <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2851.2900_22.03.08.csv')
cb <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2901.2950_22.03.08.csv')
cc <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2951.3000_22.03.08.csv')
cd <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3001.3050_22.03.08.csv')
ce <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3051.3100_22.03.08.csv')
cf <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3101.3150_22.03.08.csv')
cg <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3151.3200_22.03.08.csv')
ch <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3201.3250_22.03.08.csv')
ci <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3251.3300_22.03.08.csv')
cj <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3301.3350_22.03.08.csv')
ck <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3351.3400_22.03.08.csv')
cl <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3401.3450_22.03.08.csv')
cm <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3451.3500_22.03.08.csv')
cn <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3501.3550_22.03.08.csv')
co <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3551.3600_22.03.08.csv')
cp <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3601.3650_22.03.08.csv')
cq <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3651.3700_22.03.08.csv')
cr <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3701.3750_22.03.08.csv')
cs <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3751.3765_22.03.08.csv')

# merge 
all <- rbind(aa,ab,ac,ad,ae,af,ag,ah,ai,aj,ak,al,am,an,ao,ap,aq,ar,as,at,au,av,aw,ax,ay,az,
             ba,bb,bc,bd,be,bf,bg,bh,bi,bj,bk,bl,bm,bn,bo,bp,bq,br,bs,bt,bu,bv,bw,bx,by,bz,
             ca,cb,cc,cd,ce,cf,cg,ch,ci,cj,ck,cl,cm,cn,co,cp,cq,cr,cs)
length(unique(all$obs_id)) # 3765 -- correct number!
distinct <- distinct(all)
nrow(all) - nrow(distinct) # 13043 duplicates -- this is all of the times that a pair were seen together?? all others must be 0 or single sightings
rm(list = ls())

# short version ####
aa <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1.250_22.03.10.csv') %>% distinct()
ab <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings251.387_22.03.08.csv') %>% distinct()
ac <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings388.400_22.03.08.csv') %>% distinct()
ad <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings401.450_22.03.08.csv') %>% distinct()
ae <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings451.500_22.03.08.csv') %>% distinct()
af <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings501.550_22.03.08.csv') %>% distinct()
ag <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings551.600_22.03.08.csv') %>% distinct()
ah <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings601.650_22.03.08.csv') %>% distinct()
ai <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings651.700_22.03.08.csv') %>% distinct()
aj <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings701.750_22.03.08.csv') %>% distinct()
ak <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings751.800_22.03.08.csv') %>% distinct()
al <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings801.850_22.03.08.csv') %>% distinct()
am <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings851.900_22.03.08.csv') %>% distinct()
an <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings901.950_22.03.08.csv') %>% distinct()
ao <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings951.1000_22.03.08.csv') %>% distinct()
ap <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1001.1050_22.03.08.csv') %>% distinct()
aq <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1051.1100_22.03.08.csv') %>% distinct()
ar <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1101.1150_22.03.08.csv') %>% distinct()
as <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1151.1200_22.03.08.csv') %>% distinct()
at <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1201.1250_22.03.08.csv') %>% distinct()
au <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1251.1300_22.03.08.csv') %>% distinct()
av <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1301.1350_22.03.08.csv') %>% distinct()
aw <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1351.1400_22.03.08.csv') %>% distinct()
ax <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1401.1450_22.03.08.csv') %>% distinct()
ay <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1451.1500_22.03.08.csv') %>% distinct()
az <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1501.1550_22.03.08.csv') %>% distinct()

ba <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1551.1600_22.03.08.csv') %>% distinct()
bb <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1601.1650_22.03.08.csv') %>% distinct()
bc <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1651.1700_22.03.08.csv') %>% distinct()
bd <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1701.1750_22.03.08.csv') %>% distinct()
be <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1751.1800_22.03.08.csv') %>% distinct()
bf <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1801.1850_22.03.08.csv') %>% distinct()
bg <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1851.1900_22.03.08.csv') %>% distinct()
bh <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1901.1950_22.03.08.csv') %>% distinct()
bi <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1951.2000_22.03.08.csv') %>% distinct()
bj <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2001.2050_22.03.08.csv') %>% distinct()
bk <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2051.2100_22.03.08.csv') %>% distinct()
bl <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2101.2150_22.03.08.csv') %>% distinct()
bm <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2151.2200_22.03.08.csv') %>% distinct()
bn <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2201.2250_22.03.08.csv') %>% distinct()
bo <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2251.2300_22.03.08.csv') %>% distinct()
bp <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2301.2350_22.03.08.csv') %>% distinct()
bq <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2351.2400_22.03.08.csv') %>% distinct()
br <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2401.2450_22.03.08.csv') %>% distinct()
bs <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2451.2500_22.03.08.csv') %>% distinct()
bt <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2501.2550_22.03.08.csv') %>% distinct()
bu <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2551.2600_22.03.08.csv') %>% distinct()
bv <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2601.2650_22.03.08.csv') %>% distinct()
bw <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2651.2700_22.03.08.csv') %>% distinct()
bx <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2701.2750_22.03.08.csv') %>% distinct()
by <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2751.2800_22.03.08.csv') %>% distinct()
bz <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2801.2850_22.03.08.csv') %>% distinct()

ca <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2851.2900_22.03.08.csv') %>% distinct()
cb <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2901.2950_22.03.08.csv') %>% distinct()
cc <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2951.3000_22.03.08.csv') %>% distinct()
cd <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3001.3050_22.03.08.csv') %>% distinct()
ce <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3051.3100_22.03.08.csv') %>% distinct()
cf <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3101.3150_22.03.08.csv') %>% distinct()
cg <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3151.3200_22.03.08.csv') %>% distinct()
ch <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3201.3250_22.03.08.csv') %>% distinct()
ci <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3251.3300_22.03.08.csv') %>% distinct()
cj <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3301.3350_22.03.08.csv') %>% distinct()
ck <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3351.3400_22.03.08.csv') %>% distinct()
cl <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3401.3450_22.03.08.csv') %>% distinct()
cm <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3451.3500_22.03.08.csv') %>% distinct()
cn <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3501.3550_22.03.08.csv') %>% distinct()
co <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3551.3600_22.03.08.csv') %>% distinct()
cp <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3601.3650_22.03.08.csv') %>% distinct()
cq <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3651.3700_22.03.08.csv') %>% distinct()
cr <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3701.3750_22.03.08.csv') %>% distinct()
cs <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3751.3765_22.03.08.csv') %>% distinct()

# merge
all <- rbind(aa,ab,ac,ad,ae,af,ag,ah,ai,aj,ak,al,am,an,ao,ap,aq,ar,as,at,au,av,aw,ax,ay,az,
             ba,bb,bc,bd,be,bf,bg,bh,bi,bj,bk,bl,bm,bn,bo,bp,bq,br,bs,bt,bu,bv,bw,bx,by,bz,
             ca,cb,cc,cd,ce,cf,cg,ch,ci,cj,ck,cl,cm,cn,co,cp,cq,cr,cs)
length(unique(all$obs_id))          # 3765

# clean environment
rm(aa,ab,ac,ad,ae,af,ag,ah,ai,aj,ak,al,am,an,ao,ap,aq,ar,as,at,au,av,aw,ax,ay,az,ba,bb,bc,bd,be,bf,bg,bh,bi,bj,bk,bl,bm,bn,bo,bp,bq,br,bs,bt,bu,bv,bw,bx,by,bz,ca,cb,cc,cd,ce,cf,cg,ch,ci,cj,ck,cl,cm,cn,co,cp,cq,cr,cs)

### read in sightings data ####
#mpnp <- read_csv('data_processed/mpnp_eles_long_22.03.08.csv') %>% select(-total_id, -perc_id)
#mpnp$obs_id_intfact <- as.integer(as.factor(mpnp$encounter))
#max(mpnp$obs_id_intfact, na.rm = T) # 3735

# shrink dataset to only unique sightings
#observations <- mpnp[,c('encounter','obs_id_intfact','date','location')] %>% distinct()
#nodes <- mpnp[,c('encounter','obs_id_intfact','elephant','sex','age_range')] %>% distinct()
#rm(mpnp)

#### Look at how data were produced to match up obs_id to actual encounters ####
# sightings data
s <- readxl::read_excel('data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211214.xlsx')
str(s)
colnames(s)[c(1:23,57)] <- s[2,c(1:23,57)]
colnames(s)[24:56] <- c('CM','CF','CU','CM','CF','CU','JM','JF','JU','YPM','YPF','YPU','OPM','OPF','OPU',
                        'YAM','YAF','YAU','MAM','MAF','MAU','OAM','OAF','OAU','UM','UF','UU','SM','SF','SU',
                        'AM','AF','AU')
s <- s[3:nrow(s),]
s <- janitor::clean_names(s)

# individual data
efa <- readxl::read_excel('data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211019.xlsx')
str(efa)
efa$time_cat <- lubridate::hour(efa$Time)
efa <- separate(efa, Time, into = c('wrong_date','time'), sep = ' ')
efa$time <- hms::as_hms(efa$time)
efa <- efa[,c(1:5,7,33,8:32)]
efa$Age_Range_ID <- as.factor(efa$Age_Range_ID)
efa$Activity_ID <- as.factor(efa$Activity_ID)
efa$Distance_To_Observer <- as.numeric(efa$Distance_To_Observer)
efa$Physical_Condition_ID <- as.factor(efa$Physical_Condition_ID)
efa <- efa[,c(1:15,18:27,31:32)]
efa$Date <- lubridate::as_date(efa$Date)
efa <- janitor::clean_names(efa)

str(efa)
efa_long <- data.frame(encounter = efa$elephant_sighting_id,
                       date = efa$date,
                       time = efa$time,
                       gps_s = NA,
                       gps_e = NA,
                       total_elephants = NA,
                       total_id = NA,
                       perc_id = NA,
                       type = ifelse(efa$sex_id == 1, 'MO', 'BH/MX'),
                       elephant = ifelse(efa$elephant_id == '-', NA, efa$elephant_id),
                       sex = ifelse(efa$sex_id == 1, 'M','F/U'),
                       age_range = efa$age_range_id)

for(i in 1:length(efa_long$total_elephants)) {
  efa_long$total_elephants[i] <- length(which(efa$elephant_sighting_id == efa_long$encounter[i]))
}

efa_id <- efa_long[!is.na(efa_long$elephant),]
for(i in 1:length(efa_id$total_id)) {
  efa_id$total_id[i] <- length(which(efa_id$encounter == efa_long$encounter[i]))
}

efa_id$perc_id <- 100*(efa_id$total_id/efa_id$total_elephants)

colnames(s)[1] <- 'encounter'
s$encounter <- as.numeric(s$encounter)
efa_gps <- left_join(x = efa_id, y = s, by = 'encounter') # first few rows are blank because 2 datasets formed at slightly different times -- group sightings only goes as far as 29th January 2021, individual sightings go until 24th September 2021
efa_gps <- efa_gps[,c(1,13,10,2,3,18,19,6:9,11,12)]
colnames(efa_gps)[4] <- c('date')

efa_gps$location <- paste(efa_gps$latitude, efa_gps$longitude, sep = '_')

length(unique(efa_long$encounter)) # 5162
length(unique(efa_id$encounter))   # 3736
length(unique(efa_gps$encounter))  # 3736

efa_gps$longitude[which(efa_gps$longitude < 20)] <- NA
efa_gps$longitude[which(efa_gps$longitude < 21)] # 20.512689999999
efa_gps$longitude[which(efa_gps$latitude < -20)] # 24.750029999999999 24.707709999999999 24.707709999999999 24.768260000000001 24.764379999999999 24.765647999999999 24.528880000000001 24.528880000000001

### create group-by-individual matrix
eles_asnipe <- efa_gps[,c(1,4,5,3,14)]  # encounter, date, time, elephant, location
eles_asnipe$Date <- as.integer(eles_asnipe$date)
eles_asnipe$Date <- 1+eles_asnipe$Date - min(eles_asnipe$Date)         # start from 1, not 1st January 1970
eles_asnipe$Time <- hour(eles_asnipe$time)*60*60 + minute(eles_asnipe$time)*60 + second(eles_asnipe$time) # convert time values to seconds through day

eles_asnipe <- eles_asnipe[,c(1,6,7,4,5)]
colnames(eles_asnipe)[2:5] <- c('Date','Time','ID','Location')
eles_asnipe$ID <- as.character(eles_asnipe$ID)
eles_asnipe$group <- with(eles_asnipe, paste(Date, Time, Location, sep = '_'))

length(unique(eles_asnipe$encounter)) # 3736
length(unique(eles_asnipe$group))     # 3765

#eles_asnipe$encounter_id <- as.integer(as.factor(eles_asnipe$encounter))
#eles_asnipe$group_id <- as.integer(as.factor(eles_asnipe$group))
#eles_asnipe$newencounter <- c('new',rep(NA,nrow(eles_asnipe)-1))
#eles_asnipe$newgroup <- c('new',rep(NA,nrow(eles_asnipe)-1))
#for(i in 2:nrow(eles_asnipe)){
#  eles_asnipe$newencounter[i] <- ifelse(eles_asnipe$encounter_id[i] == eles_asnipe$encounter_id[i-1],'old','new')
#  eles_asnipe$newgroup[i]     <- ifelse(eles_asnipe$group_id[i] == eles_asnipe$group_id[i-1],'old','new')
#}
#which(eles_asnipe$newgroup != eles_asnipe$newencounter) # 203  235  263  264  544  545  699  960 1021 1112 1475 1521 2021 2142 2330 2460 2465 2901 2977 3097 3274 3372 3374 3661 3935 4120 4322 4580 4637 4758 4862 5073 5240 5301 5546 5661 5707 6702 6703 6704 6955 6956 6960 7039 7040 7047 7233 7238 7472 7492 7498 7499 7500 7566 7567 7568 7569 7570 7571 7575 7576 7578 7593 7594 7616 7617 7633 7741 7743 7774 7775 7782 7783 7792 7812 7985 8001 8002 8176 -- 79 places where group doesn't equal encounter, only 29 additional obs_id over sightings
# 203 = 2 encounters same group.  235  263  264  544  545  699  960 1021 1112 1475 1521 2021 2142 2330 2460 2465 2901 2977 3097 3274 3372 3374 3661 3935 4120 4322 4580 4637 4758 4862 5073 5240 5301 5546 5661 5707 6702 6703 6704 6955 6956 6960 7039 7040 7047 7233 7238 7472 7492 7498 7499 7500 7566 7567 7568 7569 7570 7571 7575 7576 7578 7593 7594 7616 7617 7633 7741 7743 7774 7775 7782 7783 7792 7812 7985 8001 8002 8176
#which(eles_asnipe$newgroup == 'new' & eles_asnipe$newencounter == 'old')
#which(eles_asnipe$newgroup == 'old' & eles_asnipe$newencounter == 'new')

# get_gbi generates a group by individual matrix. The function accepts a data.table with individual identifiers and a group column. The group by individual matrix can then be used to build a network using asnipe::get_network.
eles_asnipe$d_pad <- str_pad(eles_asnipe$Date, 3, pad = '0')
eles_asnipe$encounter <- paste(eles_asnipe$d_pad, eles_asnipe$Time, eles_asnipe$Location, sep = '_')
eles_asnipe$group <- as.integer(as.factor(eles_asnipe$encounter))
length(unique(eles_asnipe$encounter)) # 3765
max(eles_asnipe$group)                # 3765
eles_asnipe2 <- eles_asnipe[,c(4,6)]
eles_asnipe2 <- data.table::setDT(eles_asnipe2)
gbi_matrix <- spatsoc::get_gbi(DT = eles_asnipe2, group = 'group', id = 'ID') # 3765 rows
gbi_df <- as.data.frame(gbi_matrix)
gbi_df$identified <- rowSums(gbi_df)
summary(gbi_df$identified)

# ROW 1 IN gbi_matrix = obs_id 1 IN FINAL OUTPUT
unique(all$node_1[all$obs_id_intfact == 1])

obs_id <- rownames(gbi_matrix)
asnipe <- unique(sort(as.integer(as.factor(eles_asnipe$group))))
which(obs_id != asnipe) # none -- both 3765 long
asnipe[3765]


### SO: Create a new variable using eles_asnipe group -- do a for loop that adds 1 to the value every time a group changes and retains it whenever it is the same -- this will create a variable which is unique for 3765 encounters but which goes in the same order as the obs_id. 
### use observations dataframe which is unique by GROUP and then add number in order of group
eles_asnipe$in_order_date <- c('yes',rep(NA,nrow(eles_asnipe)-1))
eles_asnipe$in_order_time <- c('yes',rep(NA,nrow(eles_asnipe)-1))
for(i in 2:nrow(eles_asnipe)){
  eles_asnipe$in_order_date[i] <- ifelse(eles_asnipe$Date[i] > eles_asnipe$Date[i-1], 'no', 'yes')
  eles_asnipe$in_order_time[i] <- ifelse(eles_asnipe$Date[i] == eles_asnipe$Date[i-1], 
                                         ifelse(eles_asnipe$Time[i] >= eles_asnipe$Time[i-1], 'yes', 'no'),
                                         'yes') # THIS DOESN'T WORK BECAUSE THE SIGHTINGS GO IN REVERSE ORDER BY DATE, BUT WITHIN DATE GO IN CHRONOLOGICAL ORDER
}
table(eles_asnipe$in_order_date) # all dates are in order
table(eles_asnipe$in_order_time) # 122 times that are on the same day are not in order
which(eles_asnipe$in_order_time == 'no') # again, doesn't appear to be any logical order to when they do follow on and when they don't!

eles_asnipe3 <- eles_asnipe %>% arrange(desc(encounter)) %>% arrange(desc(Time)) %>% arrange(desc(Date)) # convert to an order where group is genuinely going in order
which(eles_asnipe$encounter != eles_asnipe3$encounter)                      # check to see which ones have changed -- row 9 = day 3421 = now in descending order as with date
eles_asnipe3$in_order_date <- c('yes',rep(NA,nrow(eles_asnipe3)-1))
eles_asnipe3$in_order_time <- c('yes',rep(NA,nrow(eles_asnipe3)-1))
for(i in 2:nrow(eles_asnipe3)){
  eles_asnipe3$in_order_date[i] <- ifelse(eles_asnipe3$Date[i] > eles_asnipe3$Date[i-1], 'no', 'yes')
  eles_asnipe3$in_order_time[i] <- ifelse(eles_asnipe3$Date[i] == eles_asnipe3$Date[i-1], 
                                         ifelse(eles_asnipe3$Time[i] <= eles_asnipe3$Time[i-1], 'yes', 'no'),
                                         'yes')
}
table(eles_asnipe3$in_order_date) ; table(eles_asnipe3$in_order_time)

N <- length(unique(eles_asnipe$group)) # 3765
length(unique(eles_asnipe$encounter)) - N
eles_asnipe3$group_id <- c(N, rep(NA, nrow(eles_asnipe3)-1))
for(i in 2:nrow(eles_asnipe3)){
  eles_asnipe3$group_id[i] <- ifelse(eles_asnipe3$group[i] == eles_asnipe3$group[i-1],
                                     eles_asnipe3$group_id[i-1],
                                     eles_asnipe3$group_id[i-1]-1)
}
summary(eles_asnipe3$group_id)

observations <- eles_asnipe3[,c('Date','Time','Location','group')] %>% distinct() # 3765 observations
observations$group_id <- nrow(observations):1

eles_asnipe4 <- left_join(x = eles_asnipe3, y = observations, by = 'group_id')
which(eles_asnipe4$group.x != eles_asnipe4$group.y)  # this successfully allocates the group to the id number. THIS WILL BE SUFFICIENTLY ACCURATE TO ASSIGN SIGHTINGS TO 2-YEAR TIME WINDOWS. IT WILL LIKELY NOT BE 100% ACCURATE IF THERE ARE PLACES THAT INDIVIDUAL SIGHTINGS HAVE BEEN SWITCHED AROUND (e.g. same date and time but different location). THEREFORE IT IS NOT TO BE USED FOR ANY KIND OF BERNOULLI ANALYSIS WITHOUT FURTHER CHECKING, BUT CAN BE USED FOR AGGREGATING THE MAIN DATAFRAME.

## SO... to identify which sightings in main data frame are actually represented in binomial data (all): merge 'observations' data with Binomial social interactions
rm(efa, efa_gps, efa_id, efa_long, eles_asnipe, eles_asnipe3, eles_asnipe4, gbi_matrix, s, i)
### merge sightings information into dyad data #####
colnames(observations)[5] <- 'obs_id'
all <- left_join(x = all, y = observations, by = 'obs_id')
head(all,10)
tail(all,10)

### convert to Bernoulli model data format -- can't actually use Bernoulli as would require too much computing power ####
all$dyad <- paste(all$node_1, all$node_2, sep = '_')
all$dyad_id <- as.integer(as.factor(all$dyad))
all$location_id <- as.integer(as.factor(all$location))
head(all)

colnames(all)[5:7] <- c('date','time','location')

max(all$date) - min(all$date) # 3430 days. Time period for ALERT data = 581 days --> split into 6 time windows (7 boundaries)
periods <- seq(from = min(all$date), to = max(all$date), length.out = 7)
events <- all[all$social_event == 1,] ; events <- events[!is.na(events$node_1),] # cuts down to 55490933 dyad pairs from 13043 -- assigning to time window doesn't take forever
events$period <- NA
for(i in 1:nrow(events)){
  events$period[i] <- which(periods <= events$date[i])[length(which(periods <= events$date[i]))] # take last value in vector
  if(i %% 1000 == 0) {print(i)}
}
range(events$obs_id)

periods ; View(events[c(sample(x = 1:nrow(events), size = 20, replace = F)),]) # visual check that periods have come out right

# check elephants all match up
length(unique(all$node_1))
length(unique(all$node_2))
length(unique(events$node_1))
length(unique(events$node_2))
#rm(all)

### convert to Binomial model data format -- aggregate all sightings of each dyad together into a count ####
## all time:
df_agg <- events %>%
  group_by(node_1, node_2) %>%
  summarise(event_count=sum(social_event),
            dyad_id=cur_group_id())
length(unique(df_agg$node_1)) ; length(unique(df_agg$node_2))
head(df_agg)

## by time window:
df_split <- events %>%
  group_by(node_1, node_2, period) %>%
  summarise(event_count=sum(social_event),
            dyad_id=cur_group_id())
head(df_split)

# move this bit somewhere else when cleaning!#####
eles <- read_csv('data_processed/mpnp_eles_long_22.03.08.csv') %>% 
  select(elephant, sex, age_range) %>% 
  distinct()
nrow(eles) - length(unique(eles$elephant)) # 472 elephants age/sex reclassified
unique(eles$sex)
eles <- eles[,c(1,3)] %>% distinct()
nrow(eles) - length(unique(eles$elephant)) # 450 elephants age reclassified

table(eles$age_range) # no 1, 9 is meant to be UK and 10 is not classified, but no 9s and 148 10s I'm thinking 10 is unknown age
eles$age_range_NA <- ifelse(eles$age_range == 10, NA, eles$age_range)

eles$age_unsure <- NA ; eles$age_min <- NA ; eles$age_max <- NA ; eles$age_maxmin <- NA ; eles$age_median <- NA
for (i in 1:nrow(eles)) {
  individual <- eles[eles$elephant == eles$elephant[i],]
  summary <- summary(individual$age_range_NA, na.rm = T)
  eles$age_unsure[i] <- nrow(individual)
  eles$age_min[i] <- summary[1]
  eles$age_max[i] <- summary[6]
  eles$age_maxmin[i] <- eles$age_max[i] - eles$age_min[i]
  eles$age_median[i] <- summary[3]
}
summary(eles$age_median)
rm(individual, summary)
#####
eles <- read_csv('data_processed/mpnp_eles_long_22.03.08.csv') %>% 
  select(elephant) %>% 
  distinct()
eles$node_1 <- as.integer(as.factor(eles$elephant))
colnames(eles)[1] <- c('id_1')
df <- left_join(df_split, eles, by = 'node_1') %>% distinct()
colnames(eles) <- c('id_2','node_2')
df <- left_join(df, eles, by = 'node_2') %>% distinct()
head(df) ; tail(df)

df$dyad <- paste(df$id_1, df$id_2, sep = '_')
df$dyad_id_period <- df$dyad_id              # every dyad has it's own ID number, including if same dyad in a different time window
df$dyad_id <- as.integer(as.factor(df$dyad)) # every dyad has it's own ID number, but same dyad in different windows share ID number

### create dyad row for all pairs per period ####
dyads <- data.frame(id_1 = rep(sort(eles$id_2), each = nrow(eles)),
                    id_2 = rep(sort(eles$id_2), nrow(eles)))

colnames(eles) <- c('id_1','node_1')
dyads <- left_join(dyads, eles, by = 'id_1')
colnames(eles) <- c('id_2','node_2')
dyads <- left_join(dyads, eles, by = 'id_2')
dyads <- dyads[dyads$node_1 < dyads$node_2,]

dyads <- data.frame(id_1 = rep(dyads$id_1, length(unique(df$period))),
                    id_2 = rep(dyads$id_2, length(unique(df$period))),
                    node_1 = rep(dyads$node_1, length(unique(df$period))),
                    node_2 = rep(dyads$node_2, length(unique(df$period))),
                    period = rep(sort(unique(df$period)), each = nrow(dyads)))
head(df) ; head(dyads)
















rm(all, df_agg, df_split, eles, events, observations, i, N, periods) # may need to clear environment a little before running next step





data <- left_join(x = dyads, y = df, by = c('id_1','id_2','period'))
data <- data[,c(1:5,8:10)]
colnames(data)[3:4] <- c('node_1','node_2')
data$event_count <- ifelse(is.na(data$event_count) == TRUE, 0, data$event_count)
table(data$event_count)
data$dyad <- paste(data$id_1, data$id_2, sep = '_')
data$dyad_id <- as.integer(as.factor(data$dyad))
head(data, 20)

periods <- data.frame(period_start = seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 32)[1:31],
                      period = 1:31)
data <- left_join(x = data, y = periods, by = 'period')
head(data)

### add data about nodes ####
males <- readxl::read_excel('data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()
males$id <- paste('M',males$casename, sep = '')
males <- males[,c(21,2,5,6,8,9,14,18)]
colnames(males) ; colnames(data)

colnames(males)[1] <- 'id_1'
data <- left_join(x = data, y = males, by = 'id_1')
colnames(data)[c(10:16)] <- c('id_no_1','bmo_1','byr_1','dmo_1','dyr_1','indyr_1','musthyr_1')
colnames(males)[1] <- 'id_2'
data <- left_join(x = data, y = males, by = 'id_2')
colnames(data)[c(17:23)] <- c('id_no_2','bmo_2','byr_2','dmo_2','dyr_2','indyr_2','musthyr_2')
data <- data[,c(8,7,1:5,9,6,10,17,11,18,12,19,13,20,14,21,15,22,16,23)]
head(data)

### write to file ####
readr::write_delim(data, 'data_processed/anp_bayesian_pairwiseevents_22.03.28.csv', delim = ',')

## clean environment
rm(list = ls())
