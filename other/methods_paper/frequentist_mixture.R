#### Information ####
# script to run basic frequentist analysis of edge weights using mixture model
# to respond to Colin's comments on methods paper

#### setup ####
## load packages
library(tidyverse)
library(glmmTMB)
library(LaplacesDemon)

## set seed
set.seed(1)

#### load data ####
load('motnp_edgeweights_conditionalprior.RData')
rm(list = ls()[! ls() %in% c('counts_df')])

###### run model #####
## fit a zero-inflated Poisson model using glmmTMB
mixture_fit_together <- glmmTMB(event_count ~ count_dyad,
                                ziformula = ~1, data = counts_df, family = poisson)
mixture_fit_together_random <- glmmTMB(event_count ~ count_dyad,
                                       ziformula = ~1, data = counts_df, family = poisson)

## summary of the model
summary(mixture_fit_together)

## view data used for model
View(mixture_fit_together$frame) # takes data from both variables

## extract fit
mixture_fit_together$fit

## extract standard error
mixture_fit_together$sdr

## extract edges
predictions <- predict(object = mixture_fit_together, newdata = counts_df)
hist(predictions)


#### joint, as Dan originally mentioned ####
# Fit a zero-inflated Poisson model using glmmTMB
mixture_fit_joint <- glmmTMB(cbind(event_count, apart) ~ 1,
                             ziformula = ~1, data = counts_df, family = binomial(link = 'logit'))
mixture_fit_joint_random <- glmmTMB(cbind(event_count, apart) ~ 1 + (1|node_1_males) + (1|node_2_males),
                                    ziformula = ~1, data = counts_df, family = binomial(link = 'logit'))

# Summary of the model
summary(mixture_fit_joint)
summary(mixture_fit_joint_random)

# view data used for model
View(mixture_fit_joint$frame)
View(mixture_fit_joint_random$frame)

# extract fit
mixture_fit_joint$fit
mixture_fit_joint_random$fit

# extract standard error
mixture_fit_joint$sdr
mixture_fit_joint_random$sdr

# extract edges
predictions <- predict(object = mixture_fit_joint, newdata = counts_df)
hist(invlogit(predictions), breaks = 50) # all the same
predictions <- predict(object = mixture_fit_joint_random, newdata = counts_df)
hist(invlogit(predictions), breaks = 50) # all low



