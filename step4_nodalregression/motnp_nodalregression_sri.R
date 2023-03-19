### set up ####
library(tidyverse)
library(brms)

### read in data (produced in posterior check of motnp_nodalregression_meanage.R)
nodes <- read_csv('../data_processed/motnp_eigenvector_estimates.csv')
nodes <- nodes %>% select(-model) # this column is wrong

### compare eigenvectors from bisonR model to those from igraph
#plot(nodes$model ~ nodes$igraph, las = 1, pch = 19, xlim = c(0,1), ylim = c(0,1),
#     xlab = 'SRI eigenvector centrality', ylab = 'bisonR eigenvector centrality')
#abline(a = 0, b = 1)

## set priors
priors <- get_prior(bf(sna_cent ~ age_years,
                       phi ~ age_years,
                       zi ~ 1),
                    data = nodes,
                    family = zero_inflated_beta())
priors$prior[which(priors$coef == 'age_years')] <- 'normal(0, 0.01)'
priors$prior[which(priors$class == 'Intercept')] <- 'beta(1,10)'
priors$prior[which(priors$class == 'b' & priors$prior == '(flat)')] <- 'normal(0, 1)'

fit_prior <-  brm(bf(sna_cent ~ age_years,
                     phi ~ age_years,
                     zi ~ 1),
                  data = nodes,
                  prior = priors,
                  family = zero_inflated_beta(),
                  iter = 10000,
                  sample_prior = 'only')
#save.image('motnp_nodalregression_sri_ppcheck.RData')
gc()
pp_check(fit_prior, ndraws = 100)

fit <-  brm(bf(sna_cent ~ age_years,
               phi ~ age_years,
               zi ~ 1),
            data = nodes,
            prior = priors,
            family = zero_inflated_beta(),
            iter = 10000)

## run model
save.image('motnp_nodalregression_sri.RData')

# load('motnp_nodalregression_sri.RData') ; library(tidyverse) ; library(brms)
summary(fit)
plot(fit)
pp_check(fit, ndraws = 100) # not even slightly similar

# plot eigenvector vs age
ggplot(data = nodes)+
  geom_point(aes(x = age_years, y = sna_cent))+
  theme_classic()+
  scale_x_continuous(name = 'mean age')+
  scale_y_continuous(name = 'eigenvector centrality')

# plot eigenvector vs age category
ggplot(nodes, aes(y = sna_cent, x = age_cat,
                  fill = factor(age_cat, levels = c('40+','25-40','20-25','15-19','10-15')))) +
  geom_boxplot() +
  theme_classic() + 
  scale_fill_viridis_d()+
  theme(legend.position = "none") +
  ylab("eigenvector centrality") +
  xlab("age category")

# extract posterior intercept and mean slope
post <- as_draws(fit)
N <- length(unlist(post$`1`[1]))
post <- data.frame(draw = 1:N,
                   intercept = c(unlist(post$`1`[1]),
                                 unlist(post$`2`[1]),
                                 unlist(post$`3`[1]),
                                 unlist(post$`4`[1])),
                   slope = c(unlist(post$`1`[2]),
                             unlist(post$`2`[2]),
                             unlist(post$`3`[2]),
                             unlist(post$`4`[2])))

# plot
age <- seq(10,55,1)
mu <- matrix(nrow = nrow(post), ncol = length(age))
for(i in 1:nrow(mu)){
  for(j in 1:ncol(mu)){
    mu[i,j] <- post$intercept[i] + post$slope[i]*age[j]
  }
}
mu_mean <- apply(mu, 2, mean)
mu_hpdi <- apply(mu, 2, HPDI, prob = 0.95)

# simulate eigenvector values from posterior
sim_data <- nodes[sample(1:nrow(nodes), size = length(age), replace = F),c(1,2,4)]
sim_data$age_years <- age
sim_eigen <- posterior_predict(object = fit, newdata = sim_data)
sim_pi <- apply(sim_eigen, 2, HPDI, prob = 0.5)
sim_pi

# plot
plot(sna_cent ~ age_years, data = nodes, col = rgb(0,0,1,0.1),
     pch = 19, las = 1, xlim = c(10,55), ylim = c(0,1),
     xlab = 'mean age', ylab = 'eigenvector centrality')
lines(age, mu_mean, lwd = 2, col = 'red')
shade(mu_hpdi, age)
shade(sim_pi, age)


summary(fit)
