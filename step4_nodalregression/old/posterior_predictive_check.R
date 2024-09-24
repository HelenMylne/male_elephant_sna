library(bisonR) ; library(brms) ; library(tidyverse) ; library(rethinking)
load('motnp_nodalregression_meanage.RData')

# extract mean of posterior for eigenvector
mean_eigen_values <- extract_metric(motnp_edge_weights_strongpriors, "node_eigen") %>%
  as.data.frame() %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = everything()) %>%
  rename(bison_node_eigen=value) %>%
  mutate(node = df_nodal$node,
         id = df_nodal$id)

# extract posterior draws for eigenvector
eigen_values <- extract_metric(motnp_edge_weights_strongpriors, "node_eigen") %>%
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'name', values_to = 'eigen') %>% 
  left_join(mean_eigen_values, by = 'name') %>% 
  rename(mean_eigen = bison_node_eigen) %>% 
  select(-name)

# combine with mean age values to create data frame for plotting scatter plots
eigen_values <- left_join(eigen_values, mean_motnp_ages, by = c('node','id'))

# extract posterior intercept and mean slope
post_eigen <- as.data.frame(as_draws_df(mean_motnp_eigen)) %>% janitor::clean_names()

# compute posterior means
age <- seq(10,55,1)
mu <- matrix(nrow = nrow(post_eigen), ncol = length(age))
for(i in 1:nrow(mu)){
  for(j in 1:ncol(mu)){
    mu[i,j] <- post_eigen$b_intercept[i] + post_eigen$b_age[i]*age[j]  # mean = intercept + slope*age
  }
}
mu_mean <- apply(mu, 2, mean)               # mean per age
mu_hpdi <- apply(mu, 2, HPDI, prob = 0.95)  # CI per age

# simulate from posterior -- not sure I've got this right, but the output graph looks OK
sim_eigen <- matrix(ncol = length(age), nrow = 1000)
for(i in 1:nrow(sim_eigen)){
  for(j in 1:ncol(sim_eigen)){
    a <- post_eigen$b_intercept[i] + post_eigen$b_age[i]*age[j]
    b <- post_eigen$sigma[i]
    sim_eigen[i,j] <- rbeta(1, shape1 = a, shape2 = b)
  }
}
hist(sim_eigen) # far too many getting a centrality score of 1
sim_pi <- apply(sim_eigen, 2, HPDI, prob = 0.95)

# compute posterior predictions -- this appears to produce something too narrow given the data?
predict <- as.data.frame(posterior_predict(mean_motnp_eigen,
                                           newdata = as.data.frame(age)))
predict_mean <- apply(predict, 2, mean)
predict_hpdi <- apply(predict, 2, HPDI, prob = 0.95)

# plot
mean_line <- data.frame(age = age, mean = mu_mean)
mean_shade <- data.frame(age = age, lb = mu_hpdi[1,], ub = mu_hpdi[2,])
pred_shade <- data.frame(age = age, lb = predict_hpdi[1,], ub = predict_hpdi[2,])
sim_shade <- data.frame(age = age, lb = sim_pi[1,], ub = sim_pi[2,])

ggplot()+
  geom_ribbon(data = mean_shade, aes(x = age, ymin = lb, ymax = ub),   # credible interval of the mean
              fill = rgb(0,0,0,0.15))+
  geom_ribbon(data = sim_shade, aes(x = age, ymin = lb, ymax = ub),    # credible interval rbeta(interval+slope*age, sigma)
              fill = rgb(0,0,0,0.15))+
  geom_ribbon(data = pred_shade, aes(x = age, ymin = lb, ymax = ub),   # credible interval from posterior_predict()
              fill = rgb(0,0,0,0.15))+
  geom_point(data = eigen_values, aes(x = age, y = eigen,              # 1000 draws per individual
                                      colour = as.factor(age2)),
             pch = 19, alpha = 0.01)+
  geom_line(data = mean_line, aes(x = age, y = mean),                  # mean posterior slope
            linewidth = 1.5, colour = 'black')+
  theme_classic()+
  theme(legend.position = 'none',
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))+
  scale_colour_viridis_d(direction = -1)+                              # colour it nicely -- direction = -1 reverses so yellow is youngest
  scale_x_continuous(name = 'mean of age distribution')+
  scale_y_continuous(name = 'eigenvector centrality\n(1000 draws/elephant)')

eigen_values$age2 <- ifelse(eigen_values$age < 15, '10-15', 
                            ifelse(eigen_values$age < 20, '15-20', 
                                   ifelse(eigen_values$age < 25, '20-25', 
                                          ifelse(eigen_values$age < 40, '25-40', '40+'))))
eigen_values$age_colours <- ifelse(eigen_values$age < 15, '#f0e328', 
                                   ifelse(eigen_values$age < 20, '#5ec962', 
                                          ifelse(eigen_values$age < 25, '#31688e', 
                                                 ifelse(eigen_values$age < 40, '#0d0887', '#440154'))))
colours <- eigen_values$age_colours
ggplot()+ # poster version
  geom_ribbon(data = mean_shade, aes(x = age, ymin = lb, ymax = ub),   # credible interval of the mean
              fill = rgb(0,0,0,0.3))+
  geom_point(data = eigen_values, aes(x = age, y = eigen),              # 1000 draws per individual
                                      #colour = age2),
             pch = 19, alpha = 0.01,
             colour = colours)+
  geom_line(data = mean_line, aes(x = age, y = mean),                  # mean posterior slope
            linewidth = 1.5, colour = 'black')+
  theme_classic()+
  theme(legend.position = 'none',
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 20))+
  #scale_colour_viridis_d(direction = -1)+                              # colour it nicely -- direction = -1 reverses so yellow is youngest
  scale_x_continuous(name = 'mean of age distribution')+
  scale_y_continuous(name = 'eigenvector centrality\n(1000 draws/elephant)')
