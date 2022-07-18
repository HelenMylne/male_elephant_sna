# ATE data model prep
### Information ####
# Data collected by Amboseli Trust for Elephants (ATE) 1972-2021
# Data supplied by Vicki Fishlock, 24th February 2022
### Set up ####
library(tidyverse)
library(rstan)
library(rethinking)
library(cmdstanr)
library(lubridate)

### read in processed data files ####
aa <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings1.250.csv')
ab <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings251.500.csv')
ac <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings501.750.csv')
ad <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings751.999.csv')
ae <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings1000.2000.csv')
af <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings2001.2250.csv')
ag <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings2251.2500.csv')
ah <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings2501.2750.csv')
ai <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings2751.2999.csv')
aj <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings3000.csv')
ak <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings3001.3250.csv')
al <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings3251.3500.csv')
am <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings3501.3750.csv')
an <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings3751.4000.csv')
ao <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings4001.4250.csv')
ap <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings4251.4500.csv')
aq <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings4501.4750.csv')
ar <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings4751.5000.csv')
as <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings5001.5250.csv')
at <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings5251.5500.csv')
au <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings5501.5750.csv')
av <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings5751.6000.csv')
aw <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings6001.6250.csv')
ax <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings6251.6500.csv')
ay <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings6501.6750.csv')
az <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings6751.7000.csv')

ba <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings7001.7250.csv')
bb <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings7251.7500.csv')
bc <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings7501.7750.csv')
bd <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings7751.8000.csv')
be <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings8001.8250.csv')
bf <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings8251.8500.csv')
bg <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings8501.8750.csv')
bh <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings8751.9000.csv')
bi <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings9001.9250.csv')
bj <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings9251.9500.csv')
bk <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings9501.9750.csv')
bl <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings9751.10000.csv')
bm <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings10001.10250.csv')
bn <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings10251.10500.csv')
bo <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings10501.10750.csv')
bp <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings10751.11000.csv')
bq <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings11001.11250.csv')
br <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings11251.11500.csv')
bs <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings11501.11750.csv')
bt <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings11751.12000.csv')
bu <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings12001.12250.csv')
bv <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings12251.12500.csv')
bw <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings12501.12750.csv')
bx <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings12751.13000.csv')
by <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings13001.13250.csv')
bz <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings13251.13500.csv')

ca <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings13501.13750.csv')
cb <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings13751.14000.csv')
cc <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings14001.14250.csv')
cd <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings14251.14500.csv')
ce <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings14501.14750.csv')
cf <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings14751.15000.csv')
cg <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings15001.15250.csv')
ch <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings15251.15500.csv')
ci <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings15501.15750.csv')
cj <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings15751.16000.csv')
ck <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings16001.16250.csv')
cl <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings16251.16500.csv')
cm <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings16501.16750.csv')
cn <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings16751.17000.csv')
co <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings17001.17250.csv')
cp <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings17251.17500.csv')
cq <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings17501.17750.csv')
cr <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings17751.18000.csv')
cs <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings18001.18250.csv')
ct <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings18251.18500.csv')
cu <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings18501.18750.csv')
cv <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings18751.19000.csv')
cw <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings19001.19250.csv')
cx <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings19251.19500.csv')
cy <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings19501.19750.csv')
cz <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings19751.20000.csv')

da <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings20001.20250.csv')
db <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings20251.20500.csv')
dc <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings20501.20750.csv')
dd <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings20751.21000.csv')
de <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings21001.21250.csv')
df <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings21251.21500.csv')
dg <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings21501.21750.csv')
dh <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings21751.22000.csv')
di <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings22001.22250.csv')
dj <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings22251.22500.csv')
dk <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings22501.22750.csv')
dl <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings22751.23000.csv')
dm <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings23001.23250.csv')
dn <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings23251.23500.csv')
do <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings23501.23750.csv')
dp <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings23751.24000.csv')
dq <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings24001.24174.csv')

all <- rbind(aa,ab,ac,ad,ae,af,ag,ah,ai,aj,ak,al,am,an,ao,ap,aq,ar,as,at,au,av,aw,ax,ay,az,
             ba,bb,bc,bd,be,bf,bg,bh,bi,bj,bk,bl,bm,bn,bo,bp,bq,br,bs,bt,bu,bv,bw,bx,by,bz,
             ca,cb,cc,cd,ce,cf,cg,ch,ci,cj,ck,cl,cm,cn,co,cp,cq,cr,cs,ct,cu,cv,cw,cx,cy,cz,
             da,db,dc,dd,de,df,dg,dh,di,dj,dk,dl,dm,dn,do,dp,dq)
all_distinct <- distinct(all)
nrow(all) - nrow(all_distinct) # 114209 = number of duplicates caused by pairs of elephants observed together in a sighting

# testing - just ignore this bit ####
aa <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings1.250.csv') %>% distinct()
ab <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings251.500.csv') %>% distinct()
ac <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings501.750.csv') %>% distinct()
ad <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings751.999.csv') %>% distinct()
ae <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings1000.2000.csv') %>% distinct()
af <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings2001.2250.csv') %>% distinct()
ag <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings2251.2500.csv') %>% distinct()
ah <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings2501.2750.csv') %>% distinct()
ai <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings2751.2999.csv') %>% distinct()
aj <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings3000.csv') %>% distinct()
ak <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings3001.3250.csv') %>% distinct()
al <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings3251.3500.csv') %>% distinct()
am <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings3501.3750.csv') %>% distinct()
an <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings3751.4000.csv') %>% distinct()
ao <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings4001.4250.csv') %>% distinct()
ap <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings4251.4500.csv') %>% distinct()
aq <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings4501.4750.csv') %>% distinct()
ar <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings4751.5000.csv') %>% distinct()
as <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings5001.5250.csv') %>% distinct()
at <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings5251.5500.csv') %>% distinct()
au <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings5501.5750.csv') %>% distinct()
av <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings5751.6000.csv') %>% distinct()
aw <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings6001.6250.csv') %>% distinct()
ax <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings6251.6500.csv') %>% distinct()
ay <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings6501.6750.csv') %>% distinct()
az <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings6751.7000.csv') %>% distinct()

ba <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings7001.7250.csv') %>% distinct()
bb <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings7251.7500.csv') %>% distinct()
bc <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings7501.7750.csv') %>% distinct()
bd <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings7751.8000.csv') %>% distinct()
be <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings8001.8250.csv') %>% distinct()
bf <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings8251.8500.csv') %>% distinct()
bg <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings8501.8750.csv') %>% distinct()
bh <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings8751.9000.csv') %>% distinct()
bi <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings9001.9250.csv') %>% distinct()
bj <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings9251.9500.csv') %>% distinct()
bk <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings9501.9750.csv') %>% distinct()
bl <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings9751.10000.csv') %>% distinct()
bm <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings10001.10250.csv') %>% distinct()
bn <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings10251.10500.csv') %>% distinct()
bo <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings10501.10750.csv') %>% distinct()
bp <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings10751.11000.csv') %>% distinct()
bq <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings11001.11250.csv') %>% distinct()
br <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings11251.11500.csv') %>% distinct()
bs <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings11501.11750.csv') %>% distinct()
bt <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings11751.12000.csv') %>% distinct()
bu <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings12001.12250.csv') %>% distinct()
bv <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings12251.12500.csv') %>% distinct()
bw <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings12501.12750.csv') %>% distinct()
bx <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings12751.13000.csv') %>% distinct()
by <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings13001.13250.csv') %>% distinct()
bz <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings13251.13500.csv') %>% distinct()

ca <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings13501.13750.csv') %>% distinct()
cb <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings13751.14000.csv') %>% distinct()
cc <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings14001.14250.csv') %>% distinct()
cd <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings14251.14500.csv') %>% distinct()
ce <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings14501.14750.csv') %>% distinct()
cf <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings14751.15000.csv') %>% distinct()
cg <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings15001.15250.csv') %>% distinct()
ch <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings15251.15500.csv') %>% distinct()
ci <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings15501.15750.csv') %>% distinct()
cj <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings15751.16000.csv') %>% distinct()
ck <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings16001.16250.csv') %>% distinct()
cl <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings16251.16500.csv') %>% distinct()
cm <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings16501.16750.csv') %>% distinct()
cn <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings16751.17000.csv') %>% distinct()
co <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings17001.17250.csv') %>% distinct()
cp <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings17251.17500.csv') %>% distinct()
cq <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings17501.17750.csv') %>% distinct()
cr <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings17751.18000.csv') %>% distinct()
cs <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings18001.18250.csv') %>% distinct()
ct <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings18251.18500.csv') %>% distinct()
cu <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings18501.18750.csv') %>% distinct()
cv <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings18751.19000.csv') %>% distinct()
cw <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings19001.19250.csv') %>% distinct()
cx <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings19251.19500.csv') %>% distinct()
cy <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings19501.19750.csv') %>% distinct()
cz <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings19751.20000.csv') %>% distinct()

da <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings20001.20250.csv') %>% distinct()
db <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings20251.20500.csv') %>% distinct()
dc <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings20501.20750.csv') %>% distinct()
dd <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings20751.21000.csv') %>% distinct()
de <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings21001.21250.csv') %>% distinct()
df <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings21251.21500.csv') %>% distinct()
dg <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings21501.21750.csv') %>% distinct()
dh <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings21751.22000.csv') %>% distinct()
di <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings22001.22250.csv') %>% distinct()
dj <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings22251.22500.csv') %>% distinct()
dk <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings22501.22750.csv') %>% distinct()
dl <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings22751.23000.csv') %>% distinct()
dm <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings23001.23250.csv') %>% distinct()
dn <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings23251.23500.csv') %>% distinct()
do <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings23501.23750.csv') %>% distinct()
dp <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings23751.24000.csv') %>% distinct()
dq <- read_csv('data_processed/anp_pairwiseevents/anp_bayesian_allpairwiseevents_22.03.03_sightings24001.24174.csv') %>% distinct()

test <- rbind(aa,ab,ac,ad,ae,af,ag,ah,ai,aj,ak,al,am,an,ao,ap,aq,ar,as,at,au,av,aw,ax,ay,az,
              ba,bb,bc,bd,be,bf,bg,bh,bi,bj,bk,bl,bm,bn,bo,bp,bq,br,bs,bt,bu,bv,bw,bx,by,bz,
              ca,cb,cc,cd,ce,cf,cg,ch,ci,cj,ck,cl,cm,cn,co,cp,cq,cr,cs,ct,cu,cv,cw,cx,cy,cz,
              da,db,dc,dd,de,df,dg,dh,di,dj,dk,dl,dm,dn,do,dp,dq)

nrow(all_distinct) - nrow(test) # 0 -- no overlap between sections of data prep

# clean environment ####
rm(test,all,aa,ab,ac,ad,ae,af,ag,ah,ai,aj,ak,al,am,an,ao,ap,aq,ar,as,at,au,av,aw,ax,ay,az,ba,bb,bc,bd,be,bf,bg,bh,bi,bj,bk,bl,bm,bn,bo,bp,bq,br,bs,bt,bu,bv,bw,bx,by,bz,ca,cb,cc,cd,ce,cf,cg,ch,ci,cj,ck,cl,cm,cn,co,cp,cq,cr,cs,ct,cu,cv,cw,cx,cy,cz,da,db,dc,dd,de,df,dg,dh,di,dj,dk,dl,dm,dn,do,dp,dq)

### add sightings information to data frame ####
# read in sightings data ####
ate <- readxl::read_excel(path = 'data_raw/Raw_ATE_Sightings_Fishlock220224.xlsx', sheet = 1)
ate <- janitor::clean_names(ate)
ate$id <- paste('M',ate$casename, sep = '')
ate$node_id <- as.integer(as.factor(ate$casename))
ate$obs_date <- lubridate::as_date(ate$obs_date)
ate <- separate(ate, obs_time, into = c('wrong_date','correct_time'), remove = F, sep = ' ')
ate$correct_time_hms <- hms::as_hms(ate$correct_time)
ate$corrected_time <- lubridate::hour(ate$correct_time_hms)*60*60 + lubridate::minute(ate$correct_time_hms) + lubridate::second(ate$correct_time_hms)
lu <- function(x) { length(unique(x)) }
ate_nums <- tapply(X = ate$obs_num, INDEX = ate$obs_date, FUN = lu )
ate$obs_num <- ifelse(ate$obs_num == '0','00', ate$obs_num)
ate$obs_num <- ifelse(ate$obs_num == '0a','0A', ate$obs_num)
ate$obs_num <- ifelse(ate$obs_num == '0b','0B', ate$obs_num)
ate$obs_num <- ifelse(ate$obs_num == '1','01', ate$obs_num)
ate$obs_num_std <- NA
for(i in 1:length(ate)){
  date_row <- ate[ate$obs_date == ate$obs_date[i],]
  date_row$obs_num_std <- as.integer(as.factor(sort(date_row$obs_num)))
  ate$obs_num_std[i] <- date_row$obs_num_std[which(date_row$obs_id == ate$obs_id[i])[1]]
}
ate <- ate[,c(1:3,25,26,4,27,28,8,29,9:24)]
head(ate)

# merge sightings information into data - messy ####
ate_asnipe <- ate[,c(4,1)] ; ate_asnipe <- data.table::setDT(ate_asnipe)
gbi_matrix <- spatsoc::get_gbi(DT = ate_asnipe, group = 'obs_id', id = 'id')

sightings <- ate[,c('obs_id','obs_date','correct_time_hms',"corrected_time","obs_num","obs_num_std","grid_code","hab_code_1","hab_code_2","bull_q_r","obs_type","grp_q_c","grp_size","bull_q_c","num_bulls","bulls_1_2","bulls_3_5","musth_male","oestrus_fem","act_code")] %>% distinct()
which(sightings$obs_id == 12584) ; which(sightings$obs_id == 33413) ; which(sightings$obs_id == 35241) ; which(sightings$obs_id == 40007) # none of these match the row numbers in gbi_matrix (also doesn't work with ate dataframe as that has multiple per encounter)
unique(all_distinct$obs_id[which(all_distinct$node_1 == 1 & all_distinct$social_event == 1)]) # 3584, 7057, 8464 and 12745 -- output matches gbi_matrix, just need to work out how gbi_matrix got it's values for the row numbers

nrow(gbi_matrix)                 # 24174
length(unique(sightings$obs_id)) # 24174
summary(sightings$obs_id)        # not 1:24174
which(gbi_matrix[,1] == 1)       # M1 seen 4 times: obs_id = 12584, 33413, 35241 and 40007, row numbers = 3584, 7057, 8464 and 12745
gbi_matrix[3584,1]


sightings$obs_id_intfact <- as.integer(as.factor(sightings$obs_id))
sightings$obs_id[which(sightings$obs_id_intfact == 3584)]
sightings$obs_id[which(sightings$obs_id_intfact == 7057)]
sightings$obs_id[which(sightings$obs_id_intfact == 8464)]
sightings$obs_id[which(sightings$obs_id_intfact == 12745)]

ate <- left_join(x = ate, y = sightings[,c('obs_id','obs_id_intfact')], by = 'obs_id')

# all_distinct$obs_id == 1 --> elephant 408 only
# all_distinct$obs_id == 2 --> elephants 48 & 39
ate$obs_id[which(ate$node_id == 408)]
ate$node_id[which(ate$obs_id_intfact == 1)] # sighting 5 = elephant 408
ate$obs_id[which(ate$node_id == 39)]
ate$obs_id[which(ate$node_id == 48)]
ate$node_id[which(ate$obs_id_intfact == 2)] # sighting 8 = elephants 39&48

colnames(all_distinct)[4] <- 'obs_id_intfact'
data <- left_join(x = all_distinct, y = sightings, by = 'obs_id_intfact') %>% distinct()
nrow(data) - nrow(all_distinct) # 5528 additional lines but where from??

which(data$node_1 != all_distinct$node_1) # 5584807 onwards
View(data[5584806:5584808,])
View(all_distinct[5584806:5584808,])

d <- sort(unique(data$obs_id_intfact))
a <- sort(unique(all_distinct$obs_id_intfact))
s <- sort(unique(sightings$obs_id_intfact))
s2 <- sort(sightings$obs_id_intfact)
which(d != a)
which(d != s)
which(s != a)
which(s != s2)
sightings[2775:2777,]
View(all_distinct[all_distinct$obs_id_intfact == 2192,])
sightings[sightings$obs_id_intfact == 2192,]

d <- table(data$obs_id_intfact)
a <- table(all_distinct$obs_id_intfact)
which(d != a)     # 2775 2776
d[2775] ; a[2775] # d = 2770 dyads, a = 1385 dyads
d[2776] ; a[2776] # d = 8286 dyads, a = 4143 dyads
(2770-1385) + (8286-4143) # 5528 = number of additional lines that appeared when merging data together
View(data[data$obs_id_intfact == 2775 | data$obs_id_intfact == 2776,])
View(sightings[sightings$obs_id_intfact == 2775 | sightings$obs_id_intfact == 2776,]) # duplicate line not removed by distinct() because of an NA in one column in half the rows

which(is.na(sightings$obs_num_std) == TRUE & sightings$obs_id_intfact == 2775) # 13
which(is.na(sightings$obs_num_std) == TRUE & sightings$obs_id_intfact == 2776) # 14

nrow(sightings) # 24176
sightings <- sightings[c(1:12,15:24176),]
data <- left_join(x = all_distinct, y = sightings, by = 'obs_id_intfact') %>% distinct()

# clean version ####
sightings <- ate[,c('obs_id','obs_date','correct_time_hms','obs_num_std','grid_code')] %>% distinct()

sightings$obs_id_intfact <- as.integer(as.factor(sightings$obs_id))

colnames(all_distinct)[4] <- 'obs_id_intfact'

which(is.na(sightings$obs_num_std) == TRUE & sightings$obs_id_intfact == 2775) # 13
which(is.na(sightings$obs_num_std) == TRUE & sightings$obs_id_intfact == 2776) # 14
nrow(sightings) # 24176 -- remove duplicates in rows 13 and 14
sightings <- sightings[c(1:12,15:24176),]

data <- left_join(x = all_distinct, y = sightings, by = 'obs_id_intfact') %>% distinct()

## clean environment
rm(ate_nums, ate2, date_row, i, ate_asnipe, d, a, s)
#rm(all_distinct, sightings)

### convert to Bernoulli model data format -- can't actually use Bernoulli as would require too much computing power ####
eles <- ate[,c('id','node_id')]
colnames(eles) <- c('id_1','node_1')
data <- left_join(data, eles, by = 'node_1')
colnames(eles) <- c('id_2','node_2')
data <- left_join(data, eles, by = 'node_2')

data$dyad <- paste(data$id_1, data$id_2, sep = '_')
data$dyad <- paste(data$id_1, data$id_2, sep = '_')
data$dyad_id <- as.integer(as.factor(data$dyad))
data$location_id <- as.integer(as.factor(data$location))
data <- dplyr::distinct(data)
head(data)

### convert to Binomial model data format -- aggregate all sightings of each dyad together into a count
## all time:
df_agg <- data %>%
  group_by(id_1, id_2) %>%
  summarise(event_count=sum(social_event), dyad_id=cur_group_id()) %>%
  mutate(node_1_id=as.integer(as.factor(id_1)), node_2_id=as.integer(as.factor(id_2)))
length(df_agg$id_1) == cumsum(1:471)[471] # check have correct number of dyads -- number will be the (n-1)th value of the triangular number sequence in which n = total number of elephants in analysis (472). If TRUE, correct number of pairs.
head(df_agg) ; tail(df_agg)
##   id_1  id_2  event_count dyad_id node_1_id node_2_id
##   <chr> <chr>       <dbl>   <int>     <int>     <int>
## 1 F1    F10             0       1         1         1  -- both F1 and F10 have registered as elephant number 1
## 2 F1    F100            0       2         1         2
## 3 F1    F101            0       3         1         3
## 4 F1    F102            0       4         1         4
## 5 F1    F103            0       5         1         5
## 6 F1    F104            0       6         1         6  -- F10, F100-F104 never interacted with F1
## | |     |               |       |         |         |
## | |     |               |       |         |         |
## | |     |               |       |         |         |
## 1 U67   U7              0  111151         1         1  -- U67 and U7 both registering as elephant number 1 (like F1 and F10 above)
## 2 U67   U8              0  111152         1         2
## 3 U67   U9              0  111153         1         3
## 4 U7    U8              1  111154         1         1
## 5 U7    U9              1  111155         1         2
## 6 U8    U9              2  111156         1         1
# All good except node_1_id and node_2_id reset to 1 for every new value of id_1, so node_1_id contains nothing but "1" in every cell, and node_2_id counts all values when F1 is node_1, all-1 for F10, all-2 for F100..., only U8 and U9 when U7 is node_1, and only U9 when U8 is node_1

### correct values in node_1_id and node_2_id using factor values.
df_agg <- df_agg[,c(1:4)]
df_agg$node_1 <- as.integer(as.factor(df_agg$id_1))
df_agg$node_2 <- as.integer(as.factor(df_agg$id_2))+1 # add 1 so starts at 2 and F1 is "1"
head(df_agg,10) ; tail(df_agg,10)

## per 2 year window, as in MOTNP:






















### add data about nodes
colnames(nodes)
nodes <- nodes[,c(1,3:5,9,11:13)]
nodes$id_1 <- nodes$id ; nodes$id_2 <- nodes$id
colnames(nodes) ; colnames(df_agg)
dyads <- left_join(x = df_agg, y = nodes, by = 'id_1')
colnames(dyads)[c(2,7:15)] <- c('id_2','id_pad_1','name_1','age_class_1','age_category_1','sex_1','id1_deletecolumn','count_1','dem_class_1','deletecolumn1')
dyads <- left_join(x = dyads, y = nodes, by = 'id_2')
colnames(dyads)[c(1,16:24)] <- c('id_1','id_pad_2','name_2','age_class_2','age_category_2','sex_2','id2_deletecolumn','count_2','dem_class_2','deletecolumn2')
dyads <- dyads[,c(4,1,2,5,6,3,7,16,8,17,9,18,10,19,11,20,13,22,14,23)]
head(dyads)

### remove any elephants from whom their is disagreement in the different data frames regarding their names or ID numbers
unique(gbi_distinct$id_1) # 155 females, 250 males, 65 unknowns
unique(nodes$id_1)        # 151 females, 245 males, 66 unknowns
length(which(is.na(dyads$id_pad_1)))                 # 2639 entries where elephants have no information
unique(dyads$id_1[which(is.na(dyads$id_pad_1))])     # "F157" "F158" "F176" "M125" "M13"  "M138" "M21"  "M223" "M227"
length(which(is.na(dyads$id_pad_2)))                 # 1600 entries where elephants have no information
unique(dyads$node_2[which(is.na(dyads$id_pad_1))])   # 412 elephants -- all individuals that come after F157 in the sequence
# all of these should be removed -- these are sightings of elephants which were deleted previously from the ele_nodes and ele_links data frames due to very high uncertainty in their identity, and so their sightings are unreliable.

dyads <- dyads[dyads$id_1 != "F157" & dyads$id_1 != "F158" & dyads$id_1 != "F176" & 
                 dyads$id_1 != "M125" & dyads$id_1 != "M13" & dyads$id_1 !=  "M138" & 
                 dyads$id_1 != "M21" & dyads$id_1 != "M223" & dyads$id_1 != "M227", ]
dyads <- dyads[dyads$id_2 != "F157" & dyads$id_2 != "F158" & dyads$id_2 != "F176" & 
                 dyads$id_2 != "M125" & dyads$id_2 != "M13" & dyads$id_2 !=  "M138" & 
                 dyads$id_2 != "M21" & dyads$id_2 != "M223" & dyads$id_2 != "M227", ]
length(which(is.na(dyads$id_pad_1)))                 # 0 entries where elephants have no information
length(which(is.na(dyads$id_pad_2)))                 # 0 entries where elephants have no information

### write csv
readr::write_delim(dyads, 'data_processed/motnp_bayesian_trimmedpairwiseevents_22.01.10.csv', delim = ',')

