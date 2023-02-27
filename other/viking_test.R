#### information #####
# 1) to estimate the ages of individuals in other populations, first came up with a survival curve for the Amboseli elephants which have good birth and death data.
# 2) selected the Gompertz bathtub distribution
# 3) use distribution to convert MOTNP ages in categories to probability distributions of true age (first part of this script)
# 4) extract node centrality distribution from draws of edge weight
# 5) use MOTNP age distributions to predict network centrality (third part of this script)

#### set up ####
options(future.globals.maxSize = 10000*(1024^2))   # running model with full age distributions = 5.01GB of globals, which exceeds default maximum allowed size. Set to allow up to 10 GB to allow model to run

library(tidyverse, lib.loc = '../packages/')       # library(tidyverse)
print('successfully loaded tidyverse')
library(igraph, lib.loc = '../packages/')          # library(igraph)
print('successfully loaded igraph')
library(ggdist, lib.loc = '../packages/')          # library(ggdist)
print('successfully loaded ggdist')
library(LaplacesDemon, lib.loc = '../packages/')   # library(LaplacesDemon)
print('successfully loaded LaplacesDemon')
library(posterior, lib.loc = '../packages/')       # library(posterior)
print('successfully loaded posterior')
library(Rcpp, lib.loc = '../packages/')            # library(Rcpp)
print('successfully loaded Rcpp')
library(bayesplot, lib.loc = '../packages/')       # library(bayesplot)
print('successfully loaded bayesplot')
library(cmdstanr, lib.loc = '../packages/')        # library(cmdstanr)
library(brms, lib.loc = '../packages/')            # library(brms)
#library(rstan, lib.loc = '../packages/')           # library(rstan)
library(bisonR, lib.loc = '../packages/')          # library(bisonR)

# load edge weight model and data frames
load('motnp_bisonr_edgescalculated_strongprior.RData') # currently all saved versions of this use prior N(0,1.5) for fixed, but running regression at the moment with normal(0,1) instead

print('script run')