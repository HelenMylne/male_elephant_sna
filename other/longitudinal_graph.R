age <- 1:50
beta_age <- rnorm(5, 0, 0.005)
start_age <- runif(5,0.2,0.8)

centrality <- data.frame(ele = as.factor(rep(1:5, each = 50)),
                         age = rep(age,5),
                         b_a = rep(beta_age, each = 50),
                         s_a = rep(start_age, each = 50),
                         eig = NA)
centrality$eig <- centrality$s_a + centrality$age*centrality$b_a

ggplot(centrality, aes(x = age, y = eig, colour = ele))+
  geom_smooth()+
  scale_y_continuous(limits = c(0,1))+
  theme_classic()

ba1 <- rnorm(50, beta_age[1], 0.001)
ba2 <- rnorm(50, beta_age[2], 0.001)
ba3 <- rnorm(50, beta_age[3], 0.001)
ba4 <- rnorm(50, beta_age[4], 0.001)
ba5 <- rnorm(50, beta_age[5], 0.001)

start_age <- runif(5,0.2,0.8)
centrality <- data.frame(ele = as.factor(rep(1:5, each = 50)),
                         age = rep(age,5),
                         b_a = c(ba1, ba2, ba3, ba4, ba5),
                         s_a = rep(start_age, each = 50),
                         eig = NA)
centrality$eig <- centrality$s_a + centrality$age*centrality$b_a

ggplot(centrality, aes(x = age, y = eig, colour = ele))+
  geom_smooth()+
  scale_y_continuous(limits = c(0,1), name = 'eigenvector centrality')+
  theme_classic()+
  theme(legend.position = 'none', axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))+
  scale_color_viridis_d()


