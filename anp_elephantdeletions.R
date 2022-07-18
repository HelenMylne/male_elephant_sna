### Recreate model from this paper: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009792
library(rethinking)
library(tidyverse)
library(igraph)

### Simulate population ####
# threeclangroups[i.e.,K1(bondgroupsB1,B2,B3andB4);K2(B5,B6andB7);andK3(B8)]
# Doesn't specify actual number of elephants
# 10 core groups, 8 bond groups, 3 clan groups: K1 = [B1 = 1 group of 10, B2 = 1 group of 6, B3 = 2 groups of 9 and 10, B4 = 1 group of 4], K2 = [B5 = 2 groups of 6 and 10, B6 = 1 group of 11, B7 = 1 group of 9], K3 = [B8 = 1 group of 8].
# 4 age classes (oldest is matriarch)
# AI WCG = 0.9, BCG = 0.7, BBG = 0.4, BKG = 0.1 -- NOT SPECIFIED IN PAPER, I JUST MADE THESE UP
N <- sum(10, 6, 9, 10, 4, 6, 10, 11, 9, 8)
sim_eles <- data.frame(id_1 = 1:N,
                       id_2 = 1:N,
                       core = c(rep('C1',10),rep('C2',6),rep('C3',9),rep('C4',10),rep('C5',4),
                                rep('C6',6),rep('C7',10),rep('C8',11),rep('C9',9),rep('C10',8)),
                       bond = c(rep('B1',10),rep('B2',6),rep('B3',19),rep('B4',4),rep('B5',16),
                                rep('B6',11),rep('B7',9),rep('B8',8)),
                       clan = c(rep('K1',39),rep('K2',36),rep('K3',8)),
                       age = c(4,rep(NA,N-1)))
for(i in 2:N){
  sim_eles$age[i] <- ifelse(sim_eles$core[i] == sim_eles$core[i-1], sample(x = 1:3, size = 1), 4)
}

sim_dyads <- data.frame(dyad = rep(NA,N^2),
                        dyad_id = rep(NA,N^2),
                        id_1 = rep(sim_eles$id_1, each = N),
                        id_2 = rep(sim_eles$id_1, N),
                        ai = rep(NA,N^2))
sim_dyads$dyad <- ifelse(sim_dyads$id_1 > sim_dyads$id_2,
                         paste(sim_dyads$id_2, sim_dyads$id_1, sep = '_'),
                         paste(sim_dyads$id_1, sim_dyads$id_2, sep = '_'))
sim_dyads <- sim_dyads[which(sim_dyads$id_1 != sim_dyads$id_2),]
sim_dyads <- sim_dyads[which(sim_dyads$id_1 < sim_dyads$id_2),]
sim_dyads$dyad_id <- as.integer(as.factor(sim_dyads$dyad))
sim_dyads <- left_join(x = sim_dyads, y = sim_eles, by = 'id_1')
colnames(sim_dyads)[c(4,7:10)] <- c('id_2','core_1','bond_1','clan_1','age_1')
sim_dyads <- left_join(x = sim_dyads, y = sim_eles, by = 'id_2')
colnames(sim_dyads)[c(3,12:15)] <- c('id_1','core_2','bond_2','clan_2','age_2')
sim_dyads <- sim_dyads[,c(1:5,7,12,8,13,9,14,10,15)]

sim_dyads$together_100 <- NA ; sim_dyads$apart_100 <- NA
sim_dyads$together_200 <- NA ; sim_dyads$apart_200 <- NA
sim_dyads$together_300 <- NA ; sim_dyads$apart_300 <- NA
sim_dyads$together_400 <- NA ; sim_dyads$apart_400 <- NA
sim_dyads$together_500 <- NA ; sim_dyads$apart_500 <- NA

set.seed(1)
for(i in 1:nrow(sim_dyads)){
  sim_dyads$ai[i] <- ifelse(sim_dyads$core_1[i] == sim_dyads$core_2[i], rbeta2(1,0.9,50),
                            ifelse(sim_dyads$bond_1[i] == sim_dyads$bond_2[i], rbeta2(1,0.7,50),
                                   ifelse(sim_dyads$clan_1[i] == sim_dyads$clan_2[i], rbeta2(1,0.4,50),
                                          rbeta2(1,0.1,50))))
  sim_dyads$together_100[i] <- rbinom(n = 1, size = 100, prob = sim_dyads$ai[i])
  sim_dyads$apart_100[i] <- 100-sim_dyads$together_100[i]
  sim_dyads$together_200[i] <- sim_dyads$together_100[i] + rbinom(1, 100, prob = sim_dyads$ai[i])
  sim_dyads$apart_200[i] <- 200-sim_dyads$together_200[i]
  sim_dyads$together_300[i] <- sim_dyads$together_200[i] + rbinom(1, 100, prob = sim_dyads$ai[i])
  sim_dyads$apart_300[i] <- 300-sim_dyads$together_300[i]
  sim_dyads$together_400[i] <- sim_dyads$together_300[i] + rbinom(1, 100, prob = sim_dyads$ai[i])
  sim_dyads$apart_400[i] <- 400-sim_dyads$together_400[i]
  sim_dyads$together_500[i] <- sim_dyads$together_400[i] + rbinom(1, 100, prob = sim_dyads$ai[i])
  sim_dyads$apart_500[i] <- 500-sim_dyads$together_500[i]
}

### CALCULATE BETWEENNESS ####
# frequentist = make igraph object and use betweenness(graph, v = V(graph))
sim_dyads$type <- ifelse(sim_dyads$core_1 == sim_dyads$core_2, 'wcg',
                         ifelse(sim_dyads$bond_1 == sim_dyads$bond_2, 'bcg',
                                ifelse(sim_dyads$clan_1 == sim_dyads$clan_2, 'bbg','bclan')))
links_100 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = sim_dyads$together_100,
                        type = sim_dyads$type)
sim_graph_100 <- igraph::graph_from_data_frame(d = links_100, vertices = sim_eles, directed = F)
coords <- layout_with_fr(sim_graph_100)
plot(sim_graph_100, edge.width = E(sim_graph_100)$weight/100, edge.color = 'black',
     layout = coords)

links_200 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = sim_dyads$together_200,
                        type = sim_dyads$type)
sim_graph_200 <- igraph::graph_from_data_frame(d = links_200, vertices = sim_eles, directed = F)
coords <- layout_with_fr(sim_graph_200)
plot(sim_graph_200, edge.width = E(sim_graph_200)$weight/200, edge.color = 'black',
     layout = coords)

links_300 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = sim_dyads$together_300,
                        type = sim_dyads$type)
sim_graph_300 <- igraph::graph_from_data_frame(d = links_300, vertices = sim_eles, directed = F)
coords <- layout_with_fr(sim_graph_300)
plot(sim_graph_300, edge.width = E(sim_graph_300)$weight/300, edge.color = 'black',
     layout = coords)

links_400 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = sim_dyads$together_400,
                        type = sim_dyads$type)
sim_graph_400 <- igraph::graph_from_data_frame(d = links_400, vertices = sim_eles, directed = F)
coords <- layout_with_fr(sim_graph_400)
plot(sim_graph_400, edge.width = E(sim_graph_400)$weight/400, edge.color = 'black',
     layout = coords)

links_500 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = sim_dyads$together_500,
                        type = sim_dyads$type)
sim_graph_500 <- igraph::graph_from_data_frame(d = links_500, vertices = sim_eles, directed = F)
coords <- layout_with_fr(sim_graph_500)
plot(sim_graph_500, edge.width = E(sim_graph_500)$weight/500, edge.color = 'black',
     layout = coords)

#betweenness(sim_graph_500, v = V(sim_graph_500), directed = F)

btwn <- data.frame(id_1 = sim_eles$id_1,
                   #bet_100 = betweenness(sim_graph_100, v = V(sim_graph_100), directed = F),
                   bet_200 = betweenness(sim_graph_200, v = V(sim_graph_200), directed = F),
                   bet_300 = betweenness(sim_graph_300, v = V(sim_graph_300), directed = F),
                   bet_400 = betweenness(sim_graph_400, v = V(sim_graph_400), directed = F),
                   bet_500 = betweenness(sim_graph_500, v = V(sim_graph_500), directed = F))

btwn$order <- as.integer(as.factor(btwn$bet_500))
table(btwn$order)

### Make deletions: oldest 20% ####
table(sim_eles$age)  # 10 matriarchs
N*0.2                # 16.6 individuals -- delete all matriarchs and a random 6 or 7 others from class 3
N*0.04               # 3.32 -- delete 3 or 4 at a time, random order

# is it delete the elephants and then resimulate the sightings before calculating betweenness, or just delete them and then calculate immediately?
# set up matrix of which elephants to delete in each set
oldest_delete <- matrix(data = NA, nrow = 17, ncol = 100)
oldest_delete[1:10,] <- sim_eles$id_1[which(sim_eles$age == 4)]
for(i in 1:100){
  oldest_delete[11:17,i] <- sample(size = 7, sim_eles$id_1[which(sim_eles$age == 3)])
}
oldest_delete_randomised <- matrix(data = NA, nrow = 17, ncol = 100)
for(i in 1:100){
  oldest_delete_randomised[,i] <- sample(size = nrow(oldest_delete_randomised),
                                     x = oldest_delete[,i],
                                     replace = F)
}

# set up second matrix of sim_dyads values for all dyads present and NA for deleted elephants
oldest_deleted_4 <- matrix(data = NA, nrow = nrow(sim_dyads), ncol = 100)
rownames(oldest_deleted_4) <- sim_dyads$dyad
# for loop debugging ####
oldest_deleted_4[1,1] <- (sim_dyads$together_500[1])
oldest_deleted_4[2,1] <- (sim_dyads$together_500[2])

oldest_deleted_4[3,1] <- ifelse(sim_dyads$id_1[3] == 26 |
                                  sim_dyads$id_2[3] == 26,
                                NA, sim_dyads$together_500[3])
oldest_deleted_4[3,1] <- ifelse(sim_dyads$id_1[3] == 1 |
                                  sim_dyads$id_2[3] == 1,
                                NA, sim_dyads$together_500[3])
oldest_deleted_4[3,1] <- ifelse(sim_dyads$id_1[3] == oldest_delete_randomised[1,1] |
                                  sim_dyads$id_2[3] == oldest_delete_randomised[1,1],
                                NA, sim_dyads$together_500[3])
oldest_deleted_4[3,1] <- ifelse(sim_dyads$id_1[3] == oldest_delete_randomised[6,1] |
                                  sim_dyads$id_2[3] == oldest_delete_randomised[6,1],
                                NA, sim_dyads$together_500[3])

for(i in 1:nrow(sim_dyads)){
  oldest_deleted_4[i,1] <- ifelse(sim_dyads$id_1[i] == 1 |
                                    sim_dyads$id_2[i] == 1,
                                  NA, sim_dyads$together_500[i])
}
length(which(is.na(oldest_deleted_4[,1])))


for(i in 1:nrow(sim_dyads)){
  for(j in 1:100){
    for(k in 1:4){
      oldest_deleted_4[i,j] <- ifelse(sim_dyads$id_1[i] == oldest_delete_randomised[k,j] |
                                        sim_dyads$id_2[i] == oldest_delete_randomised[k,j],
                                      NA, sim_dyads$together_500[i])
    }
  }
}
length(which(is.na(oldest_deleted_4[,1])))

for(i in 1:nrow(sim_dyads)){
  for(j in 1:100){
    oldest_deleted_4[i,j] <- ifelse(sim_dyads$id_1[i] == 5 | sim_dyads$id_1[i] == 56 , NA,
                                    ifelse(sim_dyads$id_2[i] == 5 | sim_dyads$id_1[i] == 56,
                                           NA, sim_dyads$together_500[i]))
  }
}
length(which(is.na(oldest_deleted_4[,1])))

for(i in 1:nrow(sim_dyads)){
  for(j in 1:100){
    oldest_deleted_4[i,j] <- ifelse(sim_dyads$id_1[i] == 5 | 
                                      sim_dyads$id_1[i] == oldest_delete_randomised[2,j] , NA,
                                    ifelse(sim_dyads$id_2[i] == 5 |
                                             sim_dyads$id_1[i] == oldest_delete_randomised[2,j],
                                           NA, sim_dyads$together_500[i]))
  }
}
length(which(is.na(oldest_deleted_4[,1])))










for(i in 1:nrow(sim_dyads)){
  for(j in 1:100){
    for(k in 1:4){
      oldest_deleted_4[i,j] <- ifelse(sim_dyads$id_1[i] == oldest_delete_randomised[k,j], NA,
                                      ifelse(sim_dyads$id_2[i] == oldest_delete_randomised[k,j],
                                             NA, sim_dyads$together_500[i]))
    }
  }
}
length(which(is.na(oldest_deleted_4[,1])))






oldest_deleted_4 <-matrix(data = c(sim_dyads$dyad,
                                   rep(sim_dyads$together_500, 100)),
                          nrow = nrow(sim_dyads), ncol = 101)
for(i in 1:nrow(sim_dyads)){
  for(j in 1:N){
    na_dyads <- sim_dyads[sim_dyads$id_1 == j | sim_dyads$id_2 == j,]
    oldest_deleted_4[i,2] <- na_dyads$together_500[which(na_dyads$dyad == oldest_deleted_4[[1]])]
  }
}
length(which(is.na(oldest_deleted_4[,1])))

for(i in 1:nrow(sim_dyads)){
  for(j in 1:100){
    for(k in 1:4){
      oldest_deleted_4[i,j] <- ifelse(sim_dyads$id_1[i] == oldest_delete_randomised[k,j] |
                                        sim_dyads$id_2[i] == oldest_delete_randomised[k,j],
                                      NA, sim_dyads$together_500[i])
    }
  }
}
length(which(is.na(oldest_deleted_4[,1])))

oldest_deleted_4 <- data.frame(deletion1 = rep(NA, nrow(sim_dyads)))
for( i in 1:nrow(sim_dyads)){
  oldest_deleted_4$deletion1[i] <- ifelse(sim_dyads$id_1[i] == oldest_delete_randomised[1,1] |
                                            sim_dyads$id_2[i] == oldest_delete_randomised[1,1],
                                          NA,sim_dyads$together_500[i])
}
length(which(is.na(oldest_deleted_4$deletion1)))

for( i in 1:nrow(sim_dyads)){
  oldest_deleted_4$deletion1[i] <- ifelse(sim_dyads$id_1[i] == oldest_delete_randomised[1,1] |
                                            sim_dyads$id_2[i] == oldest_delete_randomised[1,1] |
                                            sim_dyads$id_1[i] == oldest_delete_randomised[2,1] |
                                            sim_dyads$id_2[i] == oldest_delete_randomised[2,1] |
                                            sim_dyads$id_1[i] == oldest_delete_randomised[3,1] |
                                            sim_dyads$id_2[i] == oldest_delete_randomised[3,1] |
                                            sim_dyads$id_1[i] == oldest_delete_randomised[4,1] |
                                            sim_dyads$id_2[i] == oldest_delete_randomised[4,1],
                                          NA,sim_dyads$together_500[i])
}
length(which(is.na(oldest_deleted_4$deletion1)))

oldest_deleted_4 <- matrix(data = NA, nrow = nrow(sim_dyads), ncol = 100)
for( i in 1:nrow(sim_dyads)){
  for( j in 1:100){
    oldest_deleted_4[i,j] <- ifelse(sim_dyads$id_1[i] == oldest_delete_randomised[1,j] |
                                      sim_dyads$id_2[i] == oldest_delete_randomised[1,j] |
                                      sim_dyads$id_1[i] == oldest_delete_randomised[2,j] |
                                      sim_dyads$id_2[i] == oldest_delete_randomised[2,j] |
                                      sim_dyads$id_1[i] == oldest_delete_randomised[3,j] |
                                      sim_dyads$id_2[i] == oldest_delete_randomised[3,j] |
                                      sim_dyads$id_1[i] == oldest_delete_randomised[4,j] |
                                      sim_dyads$id_2[i] == oldest_delete_randomised[4,j],
                                    NA,sim_dyads$together_500[i])
  }
}
length(which(is.na(oldest_deleted_4[,1])))
length(which(is.na(oldest_deleted_4[,2])))
length(which(is.na(oldest_deleted_4[,3])))
length(which(is.na(oldest_deleted_4[,4])))
length(which(is.na(oldest_deleted_4[,5])))



# OK here we go! ####
oldest_deleted_4 <- matrix(data = NA, nrow = nrow(sim_dyads), ncol = 100)
for( i in 1:nrow(sim_dyads)){
  for( j in 1:100){
    oldest_deleted_4[i,j] <- ifelse(sim_dyads$id_1[i] == oldest_delete_randomised[1,j] |
                                      sim_dyads$id_2[i] == oldest_delete_randomised[1,j] |
                                      sim_dyads$id_1[i] == oldest_delete_randomised[2,j] |
                                      sim_dyads$id_2[i] == oldest_delete_randomised[2,j] |
                                      sim_dyads$id_1[i] == oldest_delete_randomised[3,j] |
                                      sim_dyads$id_2[i] == oldest_delete_randomised[3,j] |
                                      sim_dyads$id_1[i] == oldest_delete_randomised[4,j] |
                                      sim_dyads$id_2[i] == oldest_delete_randomised[4,j],
                                    NA,sim_dyads$together_500[i])
  }
}
length(which(is.na(oldest_deleted_4[,1])))

oldest_deleted_8 <- matrix(data = NA, nrow = nrow(sim_dyads), ncol = 100)
for( i in 1:nrow(sim_dyads)){
  for( j in 1:100){
    oldest_deleted_8[i,j] <- ifelse(sim_dyads$id_1[i] == oldest_delete_randomised[1,j] |
                                      sim_dyads$id_2[i] == oldest_delete_randomised[1,j] |
                                      sim_dyads$id_1[i] == oldest_delete_randomised[2,j] |
                                      sim_dyads$id_2[i] == oldest_delete_randomised[2,j] |
                                      sim_dyads$id_1[i] == oldest_delete_randomised[3,j] |
                                      sim_dyads$id_2[i] == oldest_delete_randomised[3,j] |
                                      sim_dyads$id_1[i] == oldest_delete_randomised[4,j] |
                                      sim_dyads$id_2[i] == oldest_delete_randomised[4,j] |
                                      sim_dyads$id_1[i] == oldest_delete_randomised[5,j] |
                                      sim_dyads$id_2[i] == oldest_delete_randomised[5,j] |
                                      sim_dyads$id_1[i] == oldest_delete_randomised[6,j] |
                                      sim_dyads$id_2[i] == oldest_delete_randomised[6,j] |
                                      sim_dyads$id_1[i] == oldest_delete_randomised[7,j] |
                                      sim_dyads$id_2[i] == oldest_delete_randomised[7,j],
                                    NA,sim_dyads$together_500[i])
  }
}
length(which(is.na(oldest_deleted_8[,1])))

oldest_deleted_12 <- matrix(data = NA, nrow = nrow(sim_dyads), ncol = 100)
for( i in 1:nrow(sim_dyads)){
  for( j in 1:100){
    oldest_deleted_12[i,j] <- ifelse(sim_dyads$id_1[i] == oldest_delete_randomised[1,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[1,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[2,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[2,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[3,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[3,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[4,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[4,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[5,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[5,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[6,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[6,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[7,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[7,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[8,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[8,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[9,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[9,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[10,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[10,j],
                                     NA,sim_dyads$together_500[i])
  }
}
length(which(is.na(oldest_deleted_12[,1])))

oldest_deleted_16 <- matrix(data = NA, nrow = nrow(sim_dyads), ncol = 100)
for( i in 1:nrow(sim_dyads)){
  for( j in 1:100){
    oldest_deleted_16[i,j] <- ifelse(sim_dyads$id_1[i] == oldest_delete_randomised[1,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[1,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[2,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[2,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[3,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[3,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[4,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[4,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[5,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[5,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[6,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[6,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[7,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[7,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[8,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[8,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[9,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[9,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[10,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[10,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[11,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[11,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[12,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[12,j] | 
                                       sim_dyads$id_1[i] == oldest_delete_randomised[13,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[13,j],
                                     NA,sim_dyads$together_500[i])
  }
}
length(which(is.na(oldest_deleted_16[,1])))


oldest_deleted_20 <- matrix(data = NA, nrow = nrow(sim_dyads), ncol = 100)
for( i in 1:nrow(sim_dyads)){
  for( j in 1:100){
    oldest_deleted_20[i,j] <- ifelse(sim_dyads$id_1[i] == oldest_delete_randomised[1,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[1,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[2,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[2,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[3,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[3,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[4,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[4,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[5,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[5,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[6,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[6,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[7,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[7,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[8,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[8,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[9,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[9,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[10,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[10,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[11,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[11,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[12,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[12,j] | 
                                       sim_dyads$id_1[i] == oldest_delete_randomised[13,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[13,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[14,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[14,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[15,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[15,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[16,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[16,j] |
                                       sim_dyads$id_1[i] == oldest_delete_randomised[17,j] |
                                       sim_dyads$id_2[i] == oldest_delete_randomised[17,j],
                                     NA,sim_dyads$together_500[i])
  }
}
length(which(is.na(oldest_deleted_20[,1])))

### Make deletions: most central 20% ####
btwn <- btwn[,c('id_1','bet_500','order')]
sim_eles <- left_join(x = sim_eles, y = btwn, by = 'id_1')
table(sim_eles$order)

f <- max(btwn$order)
central_delete <- matrix(data = NA, nrow = 17, ncol = 100)
for(i in 1:100){
  central_delete[1:17,i] <- sample(sim_eles$id_1[which(sim_eles$order > f-17)],
                                   size = 17, replace = F)
}

# set up second matrix of sim_dyads values for all dyads present and NA for deleted elephants
central_deleted_4 <- matrix(data = NA, nrow = nrow(sim_dyads), ncol = 100)
for( i in 1:nrow(sim_dyads)){
  for( j in 1:100){
    central_deleted_4[i,j] <- ifelse(sim_dyads$id_1[i] == central_delete[1,j] |
                                      sim_dyads$id_2[i] == central_delete[1,j] |
                                      sim_dyads$id_1[i] == central_delete[2,j] |
                                      sim_dyads$id_2[i] == central_delete[2,j] |
                                      sim_dyads$id_1[i] == central_delete[3,j] |
                                      sim_dyads$id_2[i] == central_delete[3,j] |
                                      sim_dyads$id_1[i] == central_delete[4,j] |
                                      sim_dyads$id_2[i] == central_delete[4,j],
                                    NA,sim_dyads$together_500[i])
  }
}
length(which(is.na(central_deleted_4[,1])))

central_deleted_8 <- matrix(data = NA, nrow = nrow(sim_dyads), ncol = 100)
for( i in 1:nrow(sim_dyads)){
  for( j in 1:100){
    central_deleted_8[i,j] <- ifelse(sim_dyads$id_1[i] == central_delete[1,j] |
                                      sim_dyads$id_2[i] == central_delete[1,j] |
                                      sim_dyads$id_1[i] == central_delete[2,j] |
                                      sim_dyads$id_2[i] == central_delete[2,j] |
                                      sim_dyads$id_1[i] == central_delete[3,j] |
                                      sim_dyads$id_2[i] == central_delete[3,j] |
                                      sim_dyads$id_1[i] == central_delete[4,j] |
                                      sim_dyads$id_2[i] == central_delete[4,j] |
                                      sim_dyads$id_1[i] == central_delete[5,j] |
                                      sim_dyads$id_2[i] == central_delete[5,j] |
                                      sim_dyads$id_1[i] == central_delete[6,j] |
                                      sim_dyads$id_2[i] == central_delete[6,j] |
                                      sim_dyads$id_1[i] == central_delete[7,j] |
                                      sim_dyads$id_2[i] == central_delete[7,j],
                                    NA,sim_dyads$together_500[i])
  }
}
length(which(is.na(central_deleted_8[,1])))

central_deleted_12 <- matrix(data = NA, nrow = nrow(sim_dyads), ncol = 100)
for( i in 1:nrow(sim_dyads)){
  for( j in 1:100){
    central_deleted_12[i,j] <- ifelse(sim_dyads$id_1[i] == central_delete[1,j] |
                                       sim_dyads$id_2[i] == central_delete[1,j] |
                                       sim_dyads$id_1[i] == central_delete[2,j] |
                                       sim_dyads$id_2[i] == central_delete[2,j] |
                                       sim_dyads$id_1[i] == central_delete[3,j] |
                                       sim_dyads$id_2[i] == central_delete[3,j] |
                                       sim_dyads$id_1[i] == central_delete[4,j] |
                                       sim_dyads$id_2[i] == central_delete[4,j] |
                                       sim_dyads$id_1[i] == central_delete[5,j] |
                                       sim_dyads$id_2[i] == central_delete[5,j] |
                                       sim_dyads$id_1[i] == central_delete[6,j] |
                                       sim_dyads$id_2[i] == central_delete[6,j] |
                                       sim_dyads$id_1[i] == central_delete[7,j] |
                                       sim_dyads$id_2[i] == central_delete[7,j] |
                                       sim_dyads$id_1[i] == central_delete[8,j] |
                                       sim_dyads$id_2[i] == central_delete[8,j] |
                                       sim_dyads$id_1[i] == central_delete[9,j] |
                                       sim_dyads$id_2[i] == central_delete[9,j] |
                                       sim_dyads$id_1[i] == central_delete[10,j] |
                                       sim_dyads$id_2[i] == central_delete[10,j],
                                     NA,sim_dyads$together_500[i])
  }
}
length(which(is.na(central_deleted_12[,1])))

central_deleted_16 <- matrix(data = NA, nrow = nrow(sim_dyads), ncol = 100)
for( i in 1:nrow(sim_dyads)){
  for( j in 1:100){
    central_deleted_16[i,j] <- ifelse(sim_dyads$id_1[i] == central_delete[1,j] |
                                       sim_dyads$id_2[i] == central_delete[1,j] |
                                       sim_dyads$id_1[i] == central_delete[2,j] |
                                       sim_dyads$id_2[i] == central_delete[2,j] |
                                       sim_dyads$id_1[i] == central_delete[3,j] |
                                       sim_dyads$id_2[i] == central_delete[3,j] |
                                       sim_dyads$id_1[i] == central_delete[4,j] |
                                       sim_dyads$id_2[i] == central_delete[4,j] |
                                       sim_dyads$id_1[i] == central_delete[5,j] |
                                       sim_dyads$id_2[i] == central_delete[5,j] |
                                       sim_dyads$id_1[i] == central_delete[6,j] |
                                       sim_dyads$id_2[i] == central_delete[6,j] |
                                       sim_dyads$id_1[i] == central_delete[7,j] |
                                       sim_dyads$id_2[i] == central_delete[7,j] |
                                       sim_dyads$id_1[i] == central_delete[8,j] |
                                       sim_dyads$id_2[i] == central_delete[8,j] |
                                       sim_dyads$id_1[i] == central_delete[9,j] |
                                       sim_dyads$id_2[i] == central_delete[9,j] |
                                       sim_dyads$id_1[i] == central_delete[10,j] |
                                       sim_dyads$id_2[i] == central_delete[10,j] |
                                       sim_dyads$id_1[i] == central_delete[11,j] |
                                       sim_dyads$id_2[i] == central_delete[11,j] |
                                       sim_dyads$id_1[i] == central_delete[12,j] |
                                       sim_dyads$id_2[i] == central_delete[12,j] | 
                                       sim_dyads$id_1[i] == central_delete[13,j] |
                                       sim_dyads$id_2[i] == central_delete[13,j],
                                     NA,sim_dyads$together_500[i])
  }
}
length(which(is.na(central_deleted_16[,1])))


central_deleted_20 <- matrix(data = NA, nrow = nrow(sim_dyads), ncol = 100)
for( i in 1:nrow(sim_dyads)){
  for( j in 1:100){
    central_deleted_20[i,j] <- ifelse(sim_dyads$id_1[i] == central_delete[1,j] |
                                       sim_dyads$id_2[i] == central_delete[1,j] |
                                       sim_dyads$id_1[i] == central_delete[2,j] |
                                       sim_dyads$id_2[i] == central_delete[2,j] |
                                       sim_dyads$id_1[i] == central_delete[3,j] |
                                       sim_dyads$id_2[i] == central_delete[3,j] |
                                       sim_dyads$id_1[i] == central_delete[4,j] |
                                       sim_dyads$id_2[i] == central_delete[4,j] |
                                       sim_dyads$id_1[i] == central_delete[5,j] |
                                       sim_dyads$id_2[i] == central_delete[5,j] |
                                       sim_dyads$id_1[i] == central_delete[6,j] |
                                       sim_dyads$id_2[i] == central_delete[6,j] |
                                       sim_dyads$id_1[i] == central_delete[7,j] |
                                       sim_dyads$id_2[i] == central_delete[7,j] |
                                       sim_dyads$id_1[i] == central_delete[8,j] |
                                       sim_dyads$id_2[i] == central_delete[8,j] |
                                       sim_dyads$id_1[i] == central_delete[9,j] |
                                       sim_dyads$id_2[i] == central_delete[9,j] |
                                       sim_dyads$id_1[i] == central_delete[10,j] |
                                       sim_dyads$id_2[i] == central_delete[10,j] |
                                       sim_dyads$id_1[i] == central_delete[11,j] |
                                       sim_dyads$id_2[i] == central_delete[11,j] |
                                       sim_dyads$id_1[i] == central_delete[12,j] |
                                       sim_dyads$id_2[i] == central_delete[12,j] | 
                                       sim_dyads$id_1[i] == central_delete[13,j] |
                                       sim_dyads$id_2[i] == central_delete[13,j] |
                                       sim_dyads$id_1[i] == central_delete[14,j] |
                                       sim_dyads$id_2[i] == central_delete[14,j] |
                                       sim_dyads$id_1[i] == central_delete[15,j] |
                                       sim_dyads$id_2[i] == central_delete[15,j] |
                                       sim_dyads$id_1[i] == central_delete[16,j] |
                                       sim_dyads$id_2[i] == central_delete[16,j] |
                                       sim_dyads$id_1[i] == central_delete[17,j] |
                                       sim_dyads$id_2[i] == central_delete[17,j],
                                     NA,sim_dyads$together_500[i])
  }
}
length(which(is.na(central_deleted_20[,1])))

### Make deletions: random 20% ####
# set up matrix of which elephants to delete in each set
random_delete <- matrix(data = NA, nrow = 17, ncol = 100)
for(i in 1:100){
  random_delete[1:17,i] <- sample(size = 17, sim_eles$id_1, replace = F)
}

# set up second matrix of sim_dyads values for all dyads present and NA for deleted elephants
random_deleted_4 <- matrix(data = NA, nrow = nrow(sim_dyads), ncol = 100)
for( i in 1:nrow(sim_dyads)){
  for( j in 1:100){
    random_deleted_4[i,j] <- ifelse(sim_dyads$id_1[i] == random_delete[1,j] |
                                      sim_dyads$id_2[i] == random_delete[1,j] |
                                      sim_dyads$id_1[i] == random_delete[2,j] |
                                      sim_dyads$id_2[i] == random_delete[2,j] |
                                      sim_dyads$id_1[i] == random_delete[3,j] |
                                      sim_dyads$id_2[i] == random_delete[3,j] |
                                      sim_dyads$id_1[i] == random_delete[4,j] |
                                      sim_dyads$id_2[i] == random_delete[4,j],
                                    NA,sim_dyads$together_500[i])
  }
}
length(which(is.na(random_deleted_4[,1])))

random_deleted_8 <- matrix(data = NA, nrow = nrow(sim_dyads), ncol = 100)
for( i in 1:nrow(sim_dyads)){
  for( j in 1:100){
    random_deleted_8[i,j] <- ifelse(sim_dyads$id_1[i] == random_delete[1,j] |
                                      sim_dyads$id_2[i] == random_delete[1,j] |
                                      sim_dyads$id_1[i] == random_delete[2,j] |
                                      sim_dyads$id_2[i] == random_delete[2,j] |
                                      sim_dyads$id_1[i] == random_delete[3,j] |
                                      sim_dyads$id_2[i] == random_delete[3,j] |
                                      sim_dyads$id_1[i] == random_delete[4,j] |
                                      sim_dyads$id_2[i] == random_delete[4,j] |
                                      sim_dyads$id_1[i] == random_delete[5,j] |
                                      sim_dyads$id_2[i] == random_delete[5,j] |
                                      sim_dyads$id_1[i] == random_delete[6,j] |
                                      sim_dyads$id_2[i] == random_delete[6,j] |
                                      sim_dyads$id_1[i] == random_delete[7,j] |
                                      sim_dyads$id_2[i] == random_delete[7,j],
                                    NA,sim_dyads$together_500[i])
  }
}
length(which(is.na(random_deleted_8[,1])))

random_deleted_12 <- matrix(data = NA, nrow = nrow(sim_dyads), ncol = 100)
for( i in 1:nrow(sim_dyads)){
  for( j in 1:100){
    random_deleted_12[i,j] <- ifelse(sim_dyads$id_1[i] == random_delete[1,j] |
                                       sim_dyads$id_2[i] == random_delete[1,j] |
                                       sim_dyads$id_1[i] == random_delete[2,j] |
                                       sim_dyads$id_2[i] == random_delete[2,j] |
                                       sim_dyads$id_1[i] == random_delete[3,j] |
                                       sim_dyads$id_2[i] == random_delete[3,j] |
                                       sim_dyads$id_1[i] == random_delete[4,j] |
                                       sim_dyads$id_2[i] == random_delete[4,j] |
                                       sim_dyads$id_1[i] == random_delete[5,j] |
                                       sim_dyads$id_2[i] == random_delete[5,j] |
                                       sim_dyads$id_1[i] == random_delete[6,j] |
                                       sim_dyads$id_2[i] == random_delete[6,j] |
                                       sim_dyads$id_1[i] == random_delete[7,j] |
                                       sim_dyads$id_2[i] == random_delete[7,j] |
                                       sim_dyads$id_1[i] == random_delete[8,j] |
                                       sim_dyads$id_2[i] == random_delete[8,j] |
                                       sim_dyads$id_1[i] == random_delete[9,j] |
                                       sim_dyads$id_2[i] == random_delete[9,j] |
                                       sim_dyads$id_1[i] == random_delete[10,j] |
                                       sim_dyads$id_2[i] == random_delete[10,j],
                                     NA,sim_dyads$together_500[i])
  }
}
length(which(is.na(random_deleted_12[,1])))

random_deleted_16 <- matrix(data = NA, nrow = nrow(sim_dyads), ncol = 100)
for( i in 1:nrow(sim_dyads)){
  for( j in 1:100){
    random_deleted_16[i,j] <- ifelse(sim_dyads$id_1[i] == random_delete[1,j] |
                                       sim_dyads$id_2[i] == random_delete[1,j] |
                                       sim_dyads$id_1[i] == random_delete[2,j] |
                                       sim_dyads$id_2[i] == random_delete[2,j] |
                                       sim_dyads$id_1[i] == random_delete[3,j] |
                                       sim_dyads$id_2[i] == random_delete[3,j] |
                                       sim_dyads$id_1[i] == random_delete[4,j] |
                                       sim_dyads$id_2[i] == random_delete[4,j] |
                                       sim_dyads$id_1[i] == random_delete[5,j] |
                                       sim_dyads$id_2[i] == random_delete[5,j] |
                                       sim_dyads$id_1[i] == random_delete[6,j] |
                                       sim_dyads$id_2[i] == random_delete[6,j] |
                                       sim_dyads$id_1[i] == random_delete[7,j] |
                                       sim_dyads$id_2[i] == random_delete[7,j] |
                                       sim_dyads$id_1[i] == random_delete[8,j] |
                                       sim_dyads$id_2[i] == random_delete[8,j] |
                                       sim_dyads$id_1[i] == random_delete[9,j] |
                                       sim_dyads$id_2[i] == random_delete[9,j] |
                                       sim_dyads$id_1[i] == random_delete[10,j] |
                                       sim_dyads$id_2[i] == random_delete[10,j] |
                                       sim_dyads$id_1[i] == random_delete[11,j] |
                                       sim_dyads$id_2[i] == random_delete[11,j] |
                                       sim_dyads$id_1[i] == random_delete[12,j] |
                                       sim_dyads$id_2[i] == random_delete[12,j] | 
                                       sim_dyads$id_1[i] == random_delete[13,j] |
                                       sim_dyads$id_2[i] == random_delete[13,j],
                                     NA,sim_dyads$together_500[i])
  }
}
length(which(is.na(random_deleted_16[,1])))


random_deleted_20 <- matrix(data = NA, nrow = nrow(sim_dyads), ncol = 100)
for( i in 1:nrow(sim_dyads)){
  for( j in 1:100){
    random_deleted_20[i,j] <- ifelse(sim_dyads$id_1[i] == random_delete[1,j] |
                                       sim_dyads$id_2[i] == random_delete[1,j] |
                                       sim_dyads$id_1[i] == random_delete[2,j] |
                                       sim_dyads$id_2[i] == random_delete[2,j] |
                                       sim_dyads$id_1[i] == random_delete[3,j] |
                                       sim_dyads$id_2[i] == random_delete[3,j] |
                                       sim_dyads$id_1[i] == random_delete[4,j] |
                                       sim_dyads$id_2[i] == random_delete[4,j] |
                                       sim_dyads$id_1[i] == random_delete[5,j] |
                                       sim_dyads$id_2[i] == random_delete[5,j] |
                                       sim_dyads$id_1[i] == random_delete[6,j] |
                                       sim_dyads$id_2[i] == random_delete[6,j] |
                                       sim_dyads$id_1[i] == random_delete[7,j] |
                                       sim_dyads$id_2[i] == random_delete[7,j] |
                                       sim_dyads$id_1[i] == random_delete[8,j] |
                                       sim_dyads$id_2[i] == random_delete[8,j] |
                                       sim_dyads$id_1[i] == random_delete[9,j] |
                                       sim_dyads$id_2[i] == random_delete[9,j] |
                                       sim_dyads$id_1[i] == random_delete[10,j] |
                                       sim_dyads$id_2[i] == random_delete[10,j] |
                                       sim_dyads$id_1[i] == random_delete[11,j] |
                                       sim_dyads$id_2[i] == random_delete[11,j] |
                                       sim_dyads$id_1[i] == random_delete[12,j] |
                                       sim_dyads$id_2[i] == random_delete[12,j] | 
                                       sim_dyads$id_1[i] == random_delete[13,j] |
                                       sim_dyads$id_2[i] == random_delete[13,j] |
                                       sim_dyads$id_1[i] == random_delete[14,j] |
                                       sim_dyads$id_2[i] == random_delete[14,j] |
                                       sim_dyads$id_1[i] == random_delete[15,j] |
                                       sim_dyads$id_2[i] == random_delete[15,j] |
                                       sim_dyads$id_1[i] == random_delete[16,j] |
                                       sim_dyads$id_2[i] == random_delete[16,j] |
                                       sim_dyads$id_1[i] == random_delete[17,j] |
                                       sim_dyads$id_2[i] == random_delete[17,j],
                                     NA,sim_dyads$together_500[i])
  }
}
length(which(is.na(random_deleted_20[,1])))

### calculate betweenness after deletion -- oldest ####
# first ####
old_4.1 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = oldest_deleted_4[,1],
                      type = sim_dyads$type)
old_4.1 <- old_4.1[!is.na(old_4.1$weight),]
length(unique(old_4.1$from))
old_4.1_g <- igraph::graph_from_data_frame(d = old_4.1, vertices = sim_eles, directed = F)

old_8.1 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = oldest_deleted_8[,1],
                      type = sim_dyads$type)
old_8.1 <- old_8.1[!is.na(old_8.1$weight),]
length(unique(old_8.1$from))
old_8.1_g <- igraph::graph_from_data_frame(d = old_8.1, vertices = sim_eles, directed = F)

old_12.1 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = oldest_deleted_12[,1],
                      type = sim_dyads$type)
old_12.1 <- old_12.1[!is.na(old_12.1$weight),]
length(unique(old_12.1$from))
old_12.1_g <- igraph::graph_from_data_frame(d = old_12.1, vertices = sim_eles, directed = F)

old_16.1 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = oldest_deleted_16[,1],
                      type = sim_dyads$type)
old_16.1 <- old_16.1[!is.na(old_16.1$weight),]
length(unique(old_16.1$from))
old_16.1_g <- igraph::graph_from_data_frame(d = old_16.1, vertices = sim_eles, directed = F)

old_20.1 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = oldest_deleted_20[,1],
                      type = sim_dyads$type)
old_20.1 <- old_20.1[!is.na(old_20.1$weight),]
length(unique(old_20.1$from))
old_20.1_g <- igraph::graph_from_data_frame(d = old_20.1, vertices = sim_eles, directed = F)

btwn_old.1 <- data.frame(id_1 = sim_eles$id_1,
                         old_4.1 = betweenness(old_4.1_g, v = V(old_4.1_g), directed = F),
                         old_8.1 = betweenness(old_8.1_g, v = V(old_8.1_g), directed = F),
                         old_12.1 = betweenness(old_12.1_g, v = V(old_12.1_g), directed = F),
                         old_16.1 = betweenness(old_16.1_g, v = V(old_16.1_g), directed = F),
                         old_20.1 = betweenness(old_20.1_g, v = V(old_20.1_g), directed = F))

colnames(btwn_old.1) <- c('id','4','8','12','16','20')
btwn_old.1_long <- pivot_longer(btwn_old.1, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_old.1_long$deletion <- as.numeric(btwn_old.1_long$deletion)
btwn_old.1_long$id <- as.factor(btwn_old.1_long$id)

ggplot(data = btwn_old.1_long, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'percentage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')

btwn_old.1_na <- btwn_old.1_long
btwn_old.1_na$betweenness <- ifelse(btwn_old.1_na$betweenness == 0, NA, btwn_old.1_na$betweenness)
sim_eles$id <- as.factor(sim_eles$id_1)
btwn_old.1_na <- left_join(x = btwn_old.1_na, y = sim_eles, by = 'id')
btwn_old.1_na$age_cat <- factor(ifelse(btwn_old.1_na$age == 1, 'young',
                                ifelse(btwn_old.1_na$age == 2, 'prime',
                                ifelse(btwn_old.1_na$age == 3, 'mature','matriarch'))),
                                levels = c('young','prime','mature','matriarch'))

ggplot(data = btwn_old.1_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  facet_wrap(.~age_cat)+
  scale_x_continuous(expand = c(0,0), 'percentage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, fewer increase than decrease

# second ####
old_4.2 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = oldest_deleted_4[,2],
                      type = sim_dyads$type)
old_4.2 <- old_4.2[!is.na(old_4.2$weight),]
length(unique(old_4.2$from))
old_4.2_g <- igraph::graph_from_data_frame(d = old_4.2, vertices = sim_eles, directed = F)

old_8.2 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = oldest_deleted_8[,2],
                      type = sim_dyads$type)
old_8.2 <- old_8.2[!is.na(old_8.2$weight),]
length(unique(old_8.2$from))
old_8.2_g <- igraph::graph_from_data_frame(d = old_8.2, vertices = sim_eles, directed = F)

old_12.2 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_12[,2],
                       type = sim_dyads$type)
old_12.2 <- old_12.2[!is.na(old_12.2$weight),]
length(unique(old_12.2$from))
old_12.2_g <- igraph::graph_from_data_frame(d = old_12.2, vertices = sim_eles, directed = F)

old_16.2 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_16[,2],
                       type = sim_dyads$type)
old_16.2 <- old_16.2[!is.na(old_16.2$weight),]
length(unique(old_16.2$from))
old_16.2_g <- igraph::graph_from_data_frame(d = old_16.2, vertices = sim_eles, directed = F)

old_20.2 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_20[,2],
                       type = sim_dyads$type)
old_20.2 <- old_20.2[!is.na(old_20.2$weight),]
length(unique(old_20.2$from))
old_20.2_g <- igraph::graph_from_data_frame(d = old_20.2, vertices = sim_eles, directed = F)

btwn_old.2 <- data.frame(id_1 = sim_eles$id_1,
                         old_4.2 = betweenness(old_4.2_g, v = V(old_4.2_g), directed = F),
                         old_8.2 = betweenness(old_8.2_g, v = V(old_8.2_g), directed = F),
                         old_12.2 = betweenness(old_12.2_g, v = V(old_12.2_g), directed = F),
                         old_16.2 = betweenness(old_16.2_g, v = V(old_16.2_g), directed = F),
                         old_20.2 = betweenness(old_20.2_g, v = V(old_20.2_g), directed = F))

colnames(btwn_old.2) <- c('id','4','8','12','16','20')
btwn_old.2_long <- pivot_longer(btwn_old.2, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_old.2_long$deletion <- as.numeric(btwn_old.2_long$deletion)
btwn_old.2_long$id <- as.factor(btwn_old.2_long$id)
btwn_old.2_na <- btwn_old.2_long
btwn_old.2_na$betweenness <- ifelse(btwn_old.2_na$betweenness == 0, NA, btwn_old.2_na$betweenness)

ggplot(data = btwn_old.2_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'percentage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, fewer increase than decrease

# third ####
old_4.3 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = oldest_deleted_4[,3],
                      type = sim_dyads$type)
old_4.3 <- old_4.3[!is.na(old_4.3$weight),]
length(unique(old_4.3$from))
old_4.3_g <- igraph::graph_from_data_frame(d = old_4.3, vertices = sim_eles, directed = F)

old_8.3 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = oldest_deleted_8[,3],
                      type = sim_dyads$type)
old_8.3 <- old_8.3[!is.na(old_8.3$weight),]
length(unique(old_8.3$from))
old_8.3_g <- igraph::graph_from_data_frame(d = old_8.3, vertices = sim_eles, directed = F)

old_12.3 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_12[,3],
                       type = sim_dyads$type)
old_12.3 <- old_12.3[!is.na(old_12.3$weight),]
length(unique(old_12.3$from))
old_12.3_g <- igraph::graph_from_data_frame(d = old_12.3, vertices = sim_eles, directed = F)

old_16.3 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_16[,3],
                       type = sim_dyads$type)
old_16.3 <- old_16.3[!is.na(old_16.3$weight),]
length(unique(old_16.3$from))
old_16.3_g <- igraph::graph_from_data_frame(d = old_16.3, vertices = sim_eles, directed = F)

old_20.3 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_20[,3],
                       type = sim_dyads$type)
old_20.3 <- old_20.3[!is.na(old_20.3$weight),]
length(unique(old_20.3$from))
old_20.3_g <- igraph::graph_from_data_frame(d = old_20.3, vertices = sim_eles, directed = F)

btwn_old.3 <- data.frame(id_1 = sim_eles$id_1,
                         old_4.3 = betweenness(old_4.3_g, v = V(old_4.3_g), directed = F),
                         old_8.3 = betweenness(old_8.3_g, v = V(old_8.3_g), directed = F),
                         old_12.3 = betweenness(old_12.3_g, v = V(old_12.3_g), directed = F),
                         old_16.3 = betweenness(old_16.3_g, v = V(old_16.3_g), directed = F),
                         old_20.3 = betweenness(old_20.3_g, v = V(old_20.3_g), directed = F))

colnames(btwn_old.3) <- c('id','4','8','12','16','20')
btwn_old.3_long <- pivot_longer(btwn_old.3, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_old.3_long$deletion <- as.numeric(btwn_old.3_long$deletion)
btwn_old.3_long$id <- as.factor(btwn_old.3_long$id)
btwn_old.3_na <- btwn_old.3_long
btwn_old.3_na$betweenness <- ifelse(btwn_old.3_na$betweenness == 0, NA, btwn_old.3_na$betweenness)

ggplot(data = btwn_old.3_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'percentage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, fewer increase than decrease

# fourth ####
old_4.4 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = oldest_deleted_4[,4],
                      type = sim_dyads$type)
old_4.4 <- old_4.4[!is.na(old_4.4$weight),]
length(unique(old_4.4$from))
old_4.4_g <- igraph::graph_from_data_frame(d = old_4.4, vertices = sim_eles, directed = F)

old_8.4 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = oldest_deleted_8[,4],
                      type = sim_dyads$type)
old_8.4 <- old_8.4[!is.na(old_8.4$weight),]
length(unique(old_8.4$from))
old_8.4_g <- igraph::graph_from_data_frame(d = old_8.4, vertices = sim_eles, directed = F)

old_12.4 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_12[,4],
                       type = sim_dyads$type)
old_12.4 <- old_12.4[!is.na(old_12.4$weight),]
length(unique(old_12.4$from))
old_12.4_g <- igraph::graph_from_data_frame(d = old_12.4, vertices = sim_eles, directed = F)

old_16.4 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_16[,4],
                       type = sim_dyads$type)
old_16.4 <- old_16.4[!is.na(old_16.4$weight),]
length(unique(old_16.4$from))
old_16.4_g <- igraph::graph_from_data_frame(d = old_16.4, vertices = sim_eles, directed = F)

old_20.4 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_20[,4],
                       type = sim_dyads$type)
old_20.4 <- old_20.4[!is.na(old_20.4$weight),]
length(unique(old_20.4$from))
old_20.4_g <- igraph::graph_from_data_frame(d = old_20.4, vertices = sim_eles, directed = F)

btwn_old.4 <- data.frame(id_1 = sim_eles$id_1,
                         old_4.4 = betweenness(old_4.4_g, v = V(old_4.4_g), directed = F),
                         old_8.4 = betweenness(old_8.4_g, v = V(old_8.4_g), directed = F),
                         old_12.4 = betweenness(old_12.4_g, v = V(old_12.4_g), directed = F),
                         old_16.4 = betweenness(old_16.4_g, v = V(old_16.4_g), directed = F),
                         old_20.4 = betweenness(old_20.4_g, v = V(old_20.4_g), directed = F))

colnames(btwn_old.4) <- c('id','4','8','12','16','20')
btwn_old.4_long <- pivot_longer(btwn_old.4, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_old.4_long$deletion <- as.numeric(btwn_old.4_long$deletion)
btwn_old.4_long$id <- as.factor(btwn_old.4_long$id)
btwn_old.4_na <- btwn_old.4_long
btwn_old.4_na$betweenness <- ifelse(btwn_old.4_na$betweenness == 0, NA, btwn_old.4_na$betweenness)

ggplot(data = btwn_old.4_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'percentage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, fewer increase than decrease

# fifth ####
old_4.5 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = oldest_deleted_4[,5],
                      type = sim_dyads$type)
old_4.5 <- old_4.5[!is.na(old_4.5$weight),]
length(unique(old_4.5$from))
old_4.5_g <- igraph::graph_from_data_frame(d = old_4.5, vertices = sim_eles, directed = F)

old_8.5 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = oldest_deleted_8[,5],
                      type = sim_dyads$type)
old_8.5 <- old_8.5[!is.na(old_8.5$weight),]
length(unique(old_8.5$from))
old_8.5_g <- igraph::graph_from_data_frame(d = old_8.5, vertices = sim_eles, directed = F)

old_12.5 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_12[,5],
                       type = sim_dyads$type)
old_12.5 <- old_12.5[!is.na(old_12.5$weight),]
length(unique(old_12.5$from))
old_12.5_g <- igraph::graph_from_data_frame(d = old_12.5, vertices = sim_eles, directed = F)

old_16.5 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_16[,5],
                       type = sim_dyads$type)
old_16.5 <- old_16.5[!is.na(old_16.5$weight),]
length(unique(old_16.5$from))
old_16.5_g <- igraph::graph_from_data_frame(d = old_16.5, vertices = sim_eles, directed = F)

old_20.5 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_20[,5],
                       type = sim_dyads$type)
old_20.5 <- old_20.5[!is.na(old_20.5$weight),]
length(unique(old_20.5$from))
old_20.5_g <- igraph::graph_from_data_frame(d = old_20.5, vertices = sim_eles, directed = F)

btwn_old.5 <- data.frame(id_1 = sim_eles$id_1,
                         old_4.5 = betweenness(old_4.5_g, v = V(old_4.5_g), directed = F),
                         old_8.5 = betweenness(old_8.5_g, v = V(old_8.5_g), directed = F),
                         old_12.5 = betweenness(old_12.5_g, v = V(old_12.5_g), directed = F),
                         old_16.5 = betweenness(old_16.5_g, v = V(old_16.5_g), directed = F),
                         old_20.5 = betweenness(old_20.5_g, v = V(old_20.5_g), directed = F))

colnames(btwn_old.5) <- c('id','4','8','12','16','20')
btwn_old.5_long <- pivot_longer(btwn_old.5, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_old.5_long$deletion <- as.numeric(btwn_old.5_long$deletion)
btwn_old.5_long$id <- as.factor(btwn_old.5_long$id)
btwn_old.5_na <- btwn_old.5_long
btwn_old.5_na$betweenness <- ifelse(btwn_old.5_na$betweenness == 0, NA, btwn_old.5_na$betweenness)

ggplot(data = btwn_old.5_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'percentage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, fewer increase than decrease

# sixth ####
old_4.6 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = oldest_deleted_4[,6],
                      type = sim_dyads$type)
old_4.6 <- old_4.6[!is.na(old_4.6$weight),]
length(unique(old_4.6$from))
old_4.6_g <- igraph::graph_from_data_frame(d = old_4.6, vertices = sim_eles, directed = F)

old_8.6 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = oldest_deleted_8[,6],
                      type = sim_dyads$type)
old_8.6 <- old_8.6[!is.na(old_8.6$weight),]
length(unique(old_8.6$from))
old_8.6_g <- igraph::graph_from_data_frame(d = old_8.6, vertices = sim_eles, directed = F)

old_12.6 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_12[,6],
                       type = sim_dyads$type)
old_12.6 <- old_12.6[!is.na(old_12.6$weight),]
length(unique(old_12.6$from))
old_12.6_g <- igraph::graph_from_data_frame(d = old_12.6, vertices = sim_eles, directed = F)

old_16.6 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_16[,6],
                       type = sim_dyads$type)
old_16.6 <- old_16.6[!is.na(old_16.6$weight),]
length(unique(old_16.6$from))
old_16.6_g <- igraph::graph_from_data_frame(d = old_16.6, vertices = sim_eles, directed = F)

old_20.6 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_20[,6],
                       type = sim_dyads$type)
old_20.6 <- old_20.6[!is.na(old_20.6$weight),]
length(unique(old_20.6$from))
old_20.6_g <- igraph::graph_from_data_frame(d = old_20.6, vertices = sim_eles, directed = F)

btwn_old.6 <- data.frame(id_1 = sim_eles$id_1,
                         old_4.6 = betweenness(old_4.6_g, v = V(old_4.6_g), directed = F),
                         old_8.6 = betweenness(old_8.6_g, v = V(old_8.6_g), directed = F),
                         old_12.6 = betweenness(old_12.6_g, v = V(old_12.6_g), directed = F),
                         old_16.6 = betweenness(old_16.6_g, v = V(old_16.6_g), directed = F),
                         old_20.6 = betweenness(old_20.6_g, v = V(old_20.6_g), directed = F))

colnames(btwn_old.6) <- c('id','4','8','12','16','20')
btwn_old.6_long <- pivot_longer(btwn_old.6, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_old.6_long$deletion <- as.numeric(btwn_old.6_long$deletion)
btwn_old.6_long$id <- as.factor(btwn_old.6_long$id)
btwn_old.6_na <- btwn_old.6_long
btwn_old.6_na$betweenness <- ifelse(btwn_old.6_na$betweenness == 0, NA, btwn_old.6_na$betweenness)

ggplot(data = btwn_old.6_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'percentage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, fewer increase than decrease

# seventh ####
old_4.7 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = oldest_deleted_4[,7],
                      type = sim_dyads$type)
old_4.7 <- old_4.7[!is.na(old_4.7$weight),]
length(unique(old_4.7$from))
old_4.7_g <- igraph::graph_from_data_frame(d = old_4.7, vertices = sim_eles, directed = F)

old_8.7 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = oldest_deleted_8[,7],
                      type = sim_dyads$type)
old_8.7 <- old_8.7[!is.na(old_8.7$weight),]
length(unique(old_8.7$from))
old_8.7_g <- igraph::graph_from_data_frame(d = old_8.7, vertices = sim_eles, directed = F)

old_12.7 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_12[,7],
                       type = sim_dyads$type)
old_12.7 <- old_12.7[!is.na(old_12.7$weight),]
length(unique(old_12.7$from))
old_12.7_g <- igraph::graph_from_data_frame(d = old_12.7, vertices = sim_eles, directed = F)

old_16.7 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_16[,7],
                       type = sim_dyads$type)
old_16.7 <- old_16.7[!is.na(old_16.7$weight),]
length(unique(old_16.7$from))
old_16.7_g <- igraph::graph_from_data_frame(d = old_16.7, vertices = sim_eles, directed = F)

old_20.7 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_20[,7],
                       type = sim_dyads$type)
old_20.7 <- old_20.7[!is.na(old_20.7$weight),]
length(unique(old_20.7$from))
old_20.7_g <- igraph::graph_from_data_frame(d = old_20.7, vertices = sim_eles, directed = F)

btwn_old.7 <- data.frame(id_1 = sim_eles$id_1,
                         old_4.7 = betweenness(old_4.7_g, v = V(old_4.7_g), directed = F),
                         old_8.7 = betweenness(old_8.7_g, v = V(old_8.7_g), directed = F),
                         old_12.7 = betweenness(old_12.7_g, v = V(old_12.7_g), directed = F),
                         old_16.7 = betweenness(old_16.7_g, v = V(old_16.7_g), directed = F),
                         old_20.7 = betweenness(old_20.7_g, v = V(old_20.7_g), directed = F))

colnames(btwn_old.7) <- c('id','4','8','12','16','20')
btwn_old.7_long <- pivot_longer(btwn_old.7, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_old.7_long$deletion <- as.numeric(btwn_old.7_long$deletion)
btwn_old.7_long$id <- as.factor(btwn_old.7_long$id)
btwn_old.7_na <- btwn_old.7_long
btwn_old.7_na$betweenness <- ifelse(btwn_old.7_na$betweenness == 0, NA, btwn_old.7_na$betweenness)

ggplot(data = btwn_old.7_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'percentage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, fewer increase than decrease

# eighth ####
old_4.8 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = oldest_deleted_4[,8],
                      type = sim_dyads$type)
old_4.8 <- old_4.8[!is.na(old_4.8$weight),]
length(unique(old_4.8$from))
old_4.8_g <- igraph::graph_from_data_frame(d = old_4.8, vertices = sim_eles, directed = F)

old_8.8 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = oldest_deleted_8[,8],
                      type = sim_dyads$type)
old_8.8 <- old_8.8[!is.na(old_8.8$weight),]
length(unique(old_8.8$from))
old_8.8_g <- igraph::graph_from_data_frame(d = old_8.8, vertices = sim_eles, directed = F)

old_12.8 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_12[,8],
                       type = sim_dyads$type)
old_12.8 <- old_12.8[!is.na(old_12.8$weight),]
length(unique(old_12.8$from))
old_12.8_g <- igraph::graph_from_data_frame(d = old_12.8, vertices = sim_eles, directed = F)

old_16.8 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_16[,8],
                       type = sim_dyads$type)
old_16.8 <- old_16.8[!is.na(old_16.8$weight),]
length(unique(old_16.8$from))
old_16.8_g <- igraph::graph_from_data_frame(d = old_16.8, vertices = sim_eles, directed = F)

old_20.8 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = oldest_deleted_20[,8],
                       type = sim_dyads$type)
old_20.8 <- old_20.8[!is.na(old_20.8$weight),]
length(unique(old_20.8$from))
old_20.8_g <- igraph::graph_from_data_frame(d = old_20.8, vertices = sim_eles, directed = F)

btwn_old.8 <- data.frame(id_1 = sim_eles$id_1,
                         old_4.8 = betweenness(old_4.8_g, v = V(old_4.8_g), directed = F),
                         old_8.8 = betweenness(old_8.8_g, v = V(old_8.8_g), directed = F),
                         old_12.8 = betweenness(old_12.8_g, v = V(old_12.8_g), directed = F),
                         old_16.8 = betweenness(old_16.8_g, v = V(old_16.8_g), directed = F),
                         old_20.8 = betweenness(old_20.8_g, v = V(old_20.8_g), directed = F))

colnames(btwn_old.8) <- c('id','4','8','12','16','20')
btwn_old.8_long <- pivot_longer(btwn_old.8, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_old.8_long$deletion <- as.numeric(btwn_old.8_long$deletion)
btwn_old.8_long$id <- as.factor(btwn_old.8_long$id)
btwn_old.8_na <- btwn_old.8_long
btwn_old.8_na$betweenness <- ifelse(btwn_old.8_na$betweenness == 0, NA, btwn_old.8_na$betweenness)

ggplot(data = btwn_old.8_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'percentage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, fewer increase than decrease

# now work out how to do it for all of them at once -- ADD IN A ROW FOR INITIAL BETWEENNESS BEFORE ANY DELETIONS ####
old_4.1 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted_4[,1])],
                      to = sim_dyads$id_2[!is.na(oldest_deleted_4[,1])],
                      weight = oldest_deleted_4[!is.na(oldest_deleted_4[,1]),1],
                      type = sim_dyads$type[!is.na(oldest_deleted_4[,1])])
old_4.1_g <- igraph::graph_from_data_frame(d = old_4.1, vertices = sim_eles, directed = F)

old_8.1 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted_8[,1])],
                      to = sim_dyads$id_2[!is.na(oldest_deleted_8[,1])],
                      weight = oldest_deleted_8[!is.na(oldest_deleted_8[,1]),1],
                      type = sim_dyads$type[!is.na(oldest_deleted_8[,1])])
old_8.1_g <- igraph::graph_from_data_frame(d = old_8.1, vertices = sim_eles, directed = F)

old_12.1 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted_12[,1])],
                       to = sim_dyads$id_2[!is.na(oldest_deleted_12[,1])],
                       weight = oldest_deleted_12[!is.na(oldest_deleted_12[,1]),1],
                       type = sim_dyads$type[!is.na(oldest_deleted_12[,1])])
old_12.1_g <- igraph::graph_from_data_frame(d = old_12.1, vertices = sim_eles, directed = F)

old_16.1 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted_16[,1])],
                       to = sim_dyads$id_2[!is.na(oldest_deleted_16[,1])],
                       weight = oldest_deleted_16[!is.na(oldest_deleted_16[,1]),1],
                       type = sim_dyads$type[!is.na(oldest_deleted_16[,1])])
old_16.1_g <- igraph::graph_from_data_frame(d = old_16.1, vertices = sim_eles, directed = F)

old_20.1 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted_20[,1])],
                       to = sim_dyads$id_2[!is.na(oldest_deleted_20[,1])],
                       weight = oldest_deleted_20[!is.na(oldest_deleted_20[,1]),1],
                       type = sim_dyads$type[!is.na(oldest_deleted_20[,1])])
old_20.1_g <- igraph::graph_from_data_frame(d = old_20.1, vertices = sim_eles, directed = F)

btwn_old.1 <- data.frame(id_1 = rep(sim_eles$id_1,5),
                         deletion = rep(c(4,8,12,16,20), each = length(sim_eles$id_1)),
                         betweenness = c(
                           betweenness(old_4.1_g, v = V(old_4.1_g), directed = F),
                           betweenness(old_8.1_g, v = V(old_8.1_g), directed = F),
                           betweenness(old_12.1_g, v = V(old_12.1_g), directed = F),
                           betweenness(old_16.1_g, v = V(old_16.1_g), directed = F),
                           betweenness(old_20.1_g, v = V(old_20.1_g), directed = F)))
btwn_old.1$id <- as.factor(btwn_old.1$id)
btwn_old.1$betweenness <- ifelse(btwn_old.1$betweenness == 0, NA, btwn_old.1$betweenness)
btwn_old.1 <- left_join(x = btwn_old.1, y = sim_eles, by = 'id_1')
btwn_old.1$age_cat <- factor(ifelse(btwn_old.1$age == 1, 'young',
                                    ifelse(btwn_old.1$age == 2, 'prime',
                                           ifelse(btwn_old.1$age == 3, 'mature','matriarch'))),
                             levels = c('young','prime','mature','matriarch'))

ggplot(data = btwn_old.1, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  facet_wrap(.~age_cat)+
  scale_x_continuous(expand = c(0,0), 'percentage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')


btwn_old <- data.frame(id_1 = rep(rep(sim_eles$id_1,5),2),
                       deletion = rep(rep(c(4,8,12,16,20), each = length(sim_eles$id_1)),2),
                       run = rep(c(1:2),each = length(sim_eles$id_1)*5),
                       betweenness = rep(NA,length(sim_eles$id_1)*5*2))

for(i in 1:2){
  old_4 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted_4[,i])],
                      to = sim_dyads$id_2[!is.na(oldest_deleted_4[,i])],
                      weight = oldest_deleted_4[!is.na(oldest_deleted_4[,i]),i],
                      type = sim_dyads$type[!is.na(oldest_deleted_4[,i])])
  old_4_g <- igraph::graph_from_data_frame(d = old_4, vertices = sim_eles, directed = F)
  
  old_8 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted_8[,i])],
                      to = sim_dyads$id_2[!is.na(oldest_deleted_8[,i])],
                      weight = oldest_deleted_8[!is.na(oldest_deleted_8[,i]),i],
                      type = sim_dyads$type[!is.na(oldest_deleted_8[,i])])
  old_8_g <- igraph::graph_from_data_frame(d = old_8, vertices = sim_eles, directed = F)
  
  old_12 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted_12[,i])],
                       to = sim_dyads$id_2[!is.na(oldest_deleted_12[,i])],
                       weight = oldest_deleted_12[!is.na(oldest_deleted_12[,i]),i],
                       type = sim_dyads$type[!is.na(oldest_deleted_12[,i])])
  old_12_g <- igraph::graph_from_data_frame(d = old_12, vertices = sim_eles, directed = F)
  
  old_16 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted_16[,i])],
                       to = sim_dyads$id_2[!is.na(oldest_deleted_16[,i])],
                       weight = oldest_deleted_16[!is.na(oldest_deleted_16[,i]),i],
                       type = sim_dyads$type[!is.na(oldest_deleted_16[,i])])
  old_16_g <- igraph::graph_from_data_frame(d = old_16, vertices = sim_eles, directed = F)
  
  old_20 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted_20[,i])],
                       to = sim_dyads$id_2[!is.na(oldest_deleted_20[,i])],
                       weight = oldest_deleted_20[!is.na(oldest_deleted_20[,i]),i],
                       type = sim_dyads$type[!is.na(oldest_deleted_20[,i])])
  old_20_g <- igraph::graph_from_data_frame(d = old_20, vertices = sim_eles, directed = F)
  
  for(j in 1:length(btwn_old$betweenness)){
    for(k in 1:N){
      btwn_old$betweenness[j] <- ifelse(btwn_old$deletion[j] == 4 & btwn_old$id_1[j] == k & btwn_old$run[j] == i,
                                        betweenness(old_4_g, v = V(old_4_g), directed = F)[k],
                                  ifelse(btwn_old$deletion[j] == 8 & btwn_old$id_1[j] == k & btwn_old$run[j] == i,
                                         betweenness(old_8_g, v = V(old_8_g), directed = F)[k],
                                   ifelse(btwn_old$deletion[j] == 12 & btwn_old$id_1[j] == k & btwn_old$run[j] == i,
                                          betweenness(old_12_g, v = V(old_12_g), directed = F)[k],
                                    ifelse(btwn_old$deletion[j] == 16 & btwn_old$id_1[j] == k & btwn_old$run[j] == i,
                                           betweenness(old_16_g, v = V(old_16_g), directed = F)[k],
                                     ifelse(btwn_old$deletion[j] == 20 & btwn_old$id_1[j] == k & btwn_old$run[j] == i,
                                            betweenness(old_20_g, v = V(old_20_g), directed = F)[k],
                                            btwn_old$betweenness[j])))))
    }
  }
}

# this loop definitely runs and works for all (though takes a long time)
btwn_old <- data.frame(id_1 = rep(rep(sim_eles$id_1,5),100),
                       deletion = rep(rep(c(4,8,12,16,20), each = length(sim_eles$id_1)),100),
                       run = rep(c(1:100),each = length(sim_eles$id_1)*5),
                       betweenness = rep(NA,length(sim_eles$id_1)*5*100))
for(i in 1:100){
  old_4 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted_4[,i])],
                      to = sim_dyads$id_2[!is.na(oldest_deleted_4[,i])],
                      weight = oldest_deleted_4[!is.na(oldest_deleted_4[,i]),i],
                      type = sim_dyads$type[!is.na(oldest_deleted_4[,i])])
  old_4_g <- igraph::graph_from_data_frame(d = old_4, vertices = sim_eles, directed = F)
  
  old_8 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted_8[,i])],
                      to = sim_dyads$id_2[!is.na(oldest_deleted_8[,i])],
                      weight = oldest_deleted_8[!is.na(oldest_deleted_8[,i]),i],
                      type = sim_dyads$type[!is.na(oldest_deleted_8[,i])])
  old_8_g <- igraph::graph_from_data_frame(d = old_8, vertices = sim_eles, directed = F)
  
  old_12 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted_12[,i])],
                       to = sim_dyads$id_2[!is.na(oldest_deleted_12[,i])],
                       weight = oldest_deleted_12[!is.na(oldest_deleted_12[,i]),i],
                       type = sim_dyads$type[!is.na(oldest_deleted_12[,i])])
  old_12_g <- igraph::graph_from_data_frame(d = old_12, vertices = sim_eles, directed = F)
  
  old_16 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted_16[,i])],
                       to = sim_dyads$id_2[!is.na(oldest_deleted_16[,i])],
                       weight = oldest_deleted_16[!is.na(oldest_deleted_16[,i]),i],
                       type = sim_dyads$type[!is.na(oldest_deleted_16[,i])])
  old_16_g <- igraph::graph_from_data_frame(d = old_16, vertices = sim_eles, directed = F)
  
  old_20 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted_20[,i])],
                       to = sim_dyads$id_2[!is.na(oldest_deleted_20[,i])],
                       weight = oldest_deleted_20[!is.na(oldest_deleted_20[,i]),i],
                       type = sim_dyads$type[!is.na(oldest_deleted_20[,i])])
  old_20_g <- igraph::graph_from_data_frame(d = old_20, vertices = sim_eles, directed = F)
  
  for(j in 1:length(btwn_old$betweenness)){
    for(k in 1:N){
      btwn_old$betweenness[j] <- ifelse(btwn_old$deletion[j] == 4 & btwn_old$id_1[j] == k & btwn_old$run[j] == i,
                                        betweenness(old_4_g, v = V(old_4_g), directed = F)[k],
                                        ifelse(btwn_old$deletion[j] == 8 & btwn_old$id_1[j] == k & btwn_old$run[j] == i,
                                               betweenness(old_8_g, v = V(old_8_g), directed = F)[k],
                                               ifelse(btwn_old$deletion[j] == 12 & btwn_old$id_1[j] == k & btwn_old$run[j] == i,
                                                      betweenness(old_12_g, v = V(old_12_g), directed = F)[k],
                                                      ifelse(btwn_old$deletion[j] == 16 & btwn_old$id_1[j] == k & btwn_old$run[j] == i,
                                                             betweenness(old_16_g, v = V(old_16_g), directed = F)[k],
                                                             ifelse(btwn_old$deletion[j] == 20 & btwn_old$id_1[j] == k & btwn_old$run[j] == i,
                                                                    betweenness(old_20_g, v = V(old_20_g), directed = F)[k],
                                                                    btwn_old$betweenness[j])))))
    }
  }
}

# improve by having only one input (but haven't tested this because it takes too long to run)
btwn_old. <- data.frame(id_1 = rep(rep(sim_eles$id_1,5),100),
                        deletion = rep(rep(c(4,8,12,16,20), each = length(sim_eles$id_1)),100),
                        run = rep(c(1:100),each = length(sim_eles$id_1)*5),
                        betweenness = rep(NA,length(sim_eles$id_1)*5*100))

oldest_deleted <- array(data = 0, dim = c(3403,100,5))
oldest_deleted[,,1] <- oldest_deleted_4
oldest_deleted[,,2] <- oldest_deleted_8
oldest_deleted[,,3] <- oldest_deleted_12
oldest_deleted[,,4] <- oldest_deleted_16
oldest_deleted[,,5] <- oldest_deleted_20

for(i in 1:100){
  old_4 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted[,i,1])],
                      to = sim_dyads$id_2[!is.na(oldest_deleted[,i,1])],
                      weight = oldest_deleted_4[!is.na(oldest_deleted[,i,1]),i],
                      type = sim_dyads$type[!is.na(oldest_deleted[,i,1])])
  old_4_g <- igraph::graph_from_data_frame(d = old_4, vertices = sim_eles, directed = F)
  
  old_8 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted[,i,2])],
                      to = sim_dyads$id_2[!is.na(oldest_deleted[,i,2])],
                      weight = oldest_deleted_8[!is.na(oldest_deleted[,i,2]),i],
                      type = sim_dyads$type[!is.na(oldest_deleted[,i,2])])
  old_8_g <- igraph::graph_from_data_frame(d = old_8, vertices = sim_eles, directed = F)
  
  old_12 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted[,i,3])],
                       to = sim_dyads$id_2[!is.na(oldest_deleted[,i,3])],
                       weight = oldest_deleted_12[!is.na(oldest_deleted[,i,3]),i],
                       type = sim_dyads$type[!is.na(oldest_deleted[,i,3])])
  old_12_g <- igraph::graph_from_data_frame(d = old_12, vertices = sim_eles, directed = F)
  
  old_16 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted[,i,4])],
                       to = sim_dyads$id_2[!is.na(oldest_deleted[,i,4])],
                       weight = oldest_deleted_16[!is.na(oldest_deleted[,i,4]),i],
                       type = sim_dyads$type[!is.na(oldest_deleted[,i,4])])
  old_16_g <- igraph::graph_from_data_frame(d = old_16, vertices = sim_eles, directed = F)
  
  old_20 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted[,i,5])],
                       to = sim_dyads$id_2[!is.na(oldest_deleted[,i,5])],
                       weight = oldest_deleted_20[!is.na(oldest_deleted[,i,5]),i],
                       type = sim_dyads$type[!is.na(oldest_deleted[,i,5])])
  old_20_g <- igraph::graph_from_data_frame(d = old_20, vertices = sim_eles, directed = F)
  
  for(j in 1:length(btwn_old$betweenness)){
    for(k in 1:N){
      btwn_old$betweenness[j] <- ifelse(btwn_old$deletion[j] == 4 & btwn_old$id_1[j] == k & btwn_old$run[j] == i,
                                        betweenness(old_4_g, v = V(old_4_g), directed = F)[k],
                                        ifelse(btwn_old$deletion[j] == 8 & btwn_old$id_1[j] == k & btwn_old$run[j] == i,
                                               betweenness(old_8_g, v = V(old_8_g), directed = F)[k],
                                               ifelse(btwn_old$deletion[j] == 12 & btwn_old$id_1[j] == k & btwn_old$run[j] == i,
                                                      betweenness(old_12_g, v = V(old_12_g), directed = F)[k],
                                                      ifelse(btwn_old$deletion[j] == 16 & btwn_old$id_1[j] == k & btwn_old$run[j] == i,
                                                             betweenness(old_16_g, v = V(old_16_g), directed = F)[k],
                                                             ifelse(btwn_old$deletion[j] == 20 & btwn_old$id_1[j] == k & btwn_old$run[j] == i,
                                                                    betweenness(old_20_g, v = V(old_20_g), directed = F)[k],
                                                                    btwn_old$betweenness[j])))))
    }
  }
}

# plot
btwn_old$id <- factor(btwn_old$id_1)
ggplot(btwn_old, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  facet_wrap(.~run)+
  theme_classic()+
  theme(legend.position = 'none')

btwn_old_na <- btwn_old
btwn_old_na$betweenness <- ifelse(btwn_old_na$betweenness == 0, NA, btwn_old_na$betweenness)
ggplot(btwn_old_na[btwn_old_na$run < 26,],
       aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  facet_wrap(.~run)+
  theme_classic()+
  theme(legend.position = 'none')
ggplot(btwn_old_na[btwn_old_na$run > 25 & btwn_old_na$run < 51,],
       aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  facet_wrap(.~run)+
  theme_classic()+
  theme(legend.position = 'none')
ggplot(btwn_old_na[btwn_old_na$run > 50 & btwn_old_na$run < 76,],
       aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  facet_wrap(.~run)+
  theme_classic()+
  theme(legend.position = 'none')
ggplot(btwn_old_na[btwn_old_na$run > 75,],
       aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  facet_wrap(.~run)+
  theme_classic()+
  theme(legend.position = 'none')

btwn_old_all <- left_join(x = btwn_old_na, y = sim_eles, by = 'id_1')
btwn_old_all$age_cat <- ifelse(btwn_old_all$age == 1, 'young',
                               ifelse(btwn_old_all$age == 2, 'prime',
                                      ifelse(btwn_old_all$age == 3, 'mature', 'matriarch')))
ggplot(btwn_old_all[btwn_old_all$run < 11,],
       aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  facet_grid(age_cat ~ run)+
  theme_light()+
  theme(legend.position = 'none')
ggplot(btwn_old_all[btwn_old_all$run > 10 & btwn_old_all$run < 21,],
       aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  facet_grid(age_cat ~ run)+
  theme_light()+
  theme(legend.position = 'none')
ggplot(btwn_old_all[btwn_old_all$run > 20 & btwn_old_all$run < 31,],
       aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  facet_grid(age_cat ~ run)+
  theme_light()+
  theme(legend.position = 'none')
ggplot(btwn_old_all[btwn_old_all$run > 30 & btwn_old_all$run < 41,],
       aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  facet_grid(age_cat ~ run)+
  theme_light()+
  theme(legend.position = 'none')
ggplot(btwn_old_all[btwn_old_all$run > 40 & btwn_old_all$run < 51,],
       aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  facet_grid(age_cat ~ run)+
  theme_light()+
  theme(legend.position = 'none')
ggplot(btwn_old_all[btwn_old_all$run > 50 & btwn_old_all$run < 61,],
       aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  facet_grid(age_cat ~ run)+
  theme_light()+
  theme(legend.position = 'none')
ggplot(btwn_old_all[btwn_old_all$run > 60 & btwn_old_all$run < 71,],
       aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  facet_grid(age_cat ~ run)+
  theme_light()+
  theme(legend.position = 'none')
ggplot(btwn_old_all[btwn_old_all$run > 70 & btwn_old_all$run < 81,],
       aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  facet_grid(age_cat ~ run)+
  theme_light()+
  theme(legend.position = 'none')
ggplot(btwn_old_all[btwn_old_all$run > 80 & btwn_old_all$run < 91,],
       aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  facet_grid(age_cat ~ run)+
  theme_light()+
  theme(legend.position = 'none')
ggplot(btwn_old_all[btwn_old_all$run > 90,],
       aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  facet_grid(age_cat ~ run)+
  theme_light()+
  theme(legend.position = 'none')

# modularity ####
metrics_old <- data.frame(run = rep(1:100,5),
                          deletion = rep(c(4,8,12,16,20),
                                         each = 100),
                          modularity = rep(NA, 5*100),
                          clustering = rep(NA, 5*100),
                          diameter = rep(NA, 5*100),
                          efficiency = rep(NA, 5*100))

for(i in 1:100){
  old_4 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted[,i,1])],
                      to = sim_dyads$id_2[!is.na(oldest_deleted[,i,1])],
                      weight = oldest_deleted_4[!is.na(oldest_deleted[,i,1]),i],
                      type = sim_dyads$type[!is.na(oldest_deleted[,i,1])])
  old_4_g <- igraph::graph_from_data_frame(d = old_4, vertices = sim_eles, directed = F)
  
  old_8 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted[,i,2])],
                      to = sim_dyads$id_2[!is.na(oldest_deleted[,i,2])],
                      weight = oldest_deleted_8[!is.na(oldest_deleted[,i,2]),i],
                      type = sim_dyads$type[!is.na(oldest_deleted[,i,2])])
  old_8_g <- igraph::graph_from_data_frame(d = old_8, vertices = sim_eles, directed = F)
  
  old_12 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted[,i,3])],
                       to = sim_dyads$id_2[!is.na(oldest_deleted[,i,3])],
                       weight = oldest_deleted_12[!is.na(oldest_deleted[,i,3]),i],
                       type = sim_dyads$type[!is.na(oldest_deleted[,i,3])])
  old_12_g <- igraph::graph_from_data_frame(d = old_12, vertices = sim_eles, directed = F)
  
  old_16 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted[,i,4])],
                       to = sim_dyads$id_2[!is.na(oldest_deleted[,i,4])],
                       weight = oldest_deleted_16[!is.na(oldest_deleted[,i,4]),i],
                       type = sim_dyads$type[!is.na(oldest_deleted[,i,4])])
  old_16_g <- igraph::graph_from_data_frame(d = old_16, vertices = sim_eles, directed = F)
  
  old_20 <- data.frame(from = sim_dyads$id_1[!is.na(oldest_deleted[,i,5])],
                       to = sim_dyads$id_2[!is.na(oldest_deleted[,i,5])],
                       weight = oldest_deleted_20[!is.na(oldest_deleted[,i,5]),i],
                       type = sim_dyads$type[!is.na(oldest_deleted[,i,5])])
  old_20_g <- igraph::graph_from_data_frame(d = old_20, vertices = sim_eles, directed = F)
  
  for(j in 1:length(metrics_old$modularity)){
      metrics_old$modularity[j] <- ifelse(metrics_old$deletion[j] == 4 & metrics_old$run[j] == i,
                                        modularity(old_4_g, membership = V(old_4_g)),
                                        ifelse(metrics_old$deletion[j] == 8 & metrics_old$run[j] == i,
                                               modularity(old_8_g, membership = V(old_8_g)),
                                               ifelse(metrics_old$deletion[j] == 12 & metrics_old$run[j] == i,
                                                      modularity(old_12_g, membership = V(old_12_g)),
                                                      ifelse(metrics_old$deletion[j] == 16 & metrics_old$run[j] == i,
                                                             modularity(old_16_g, membership = V(old_16_g)),
                                                             ifelse(metrics_old$deletion[j] == 20 & metrics_old$run[j] == i,
                                                                    modularity(old_20_g, membership = V(old_20_g)),
                                                                    metrics_old$modularity[j])))))
      metrics_old$clustering[j] <- ifelse(metrics_old$deletion[j] == 4 & metrics_old$run[j] == i,
                                          transitivity(old_4_g, type = 'undirected', vids = V(old_4_g)),
                                          ifelse(metrics_old$deletion[j] == 8 & metrics_old$run[j] == i,
                                                 transitivity(old_8_g, type = 'undirected', vids = V(old_8_g)),
                                                 ifelse(metrics_old$deletion[j] == 12 & metrics_old$run[j] == i,
                                                        transitivity(old_12_g, type = 'undirected', vids = V(old_12_g)),
                                                        ifelse(metrics_old$deletion[j] == 16 & metrics_old$run[j] == i,
                                                               transitivity(old_16_g, type = 'undirected', vids = V(old_16_g)),
                                                               ifelse(metrics_old$deletion[j] == 20 & metrics_old$run[j] == i,
                                                                      transitivity(old_20_g, type = 'undirected', vids = V(old_20_g), isolates = 'zero'),
                                                                      metrics_old$clustering[j])))))
      metrics_old$diameter[j] <- ifelse(metrics_old$deletion[j] == 4 & metrics_old$run[j] == i,
                                          diameter(old_4_g, directed = F, unconnected = T),
                                          ifelse(metrics_old$deletion[j] == 8 & metrics_old$run[j] == i,
                                                 diameter(old_8_g, directed = F, unconnected = T),
                                                 ifelse(metrics_old$deletion[j] == 12 & metrics_old$run[j] == i,
                                                        diameter(old_12_g, directed = F, unconnected = T),
                                                        ifelse(metrics_old$deletion[j] == 16 & metrics_old$run[j] == i,
                                                               diameter(old_16_g, directed = F, unconnected = T),
                                                               ifelse(metrics_old$deletion[j] == 20 & metrics_old$run[j] == i,
                                                                      diameter(old_20_g, directed = F, unconnected = T),
                                                                      metrics_old$diameter[j])))))
      metrics_old$efficiency[j] <- ifelse(metrics_old$deletion[j] == 4 & metrics_old$run[j] == i,
                                          brainGraph::efficiency(old_4_g, type = 'global'),
                                          ifelse(metrics_old$deletion[j] == 8 & metrics_old$run[j] == i,
                                                 brainGraph::efficiency(old_8_g, type = 'global'),
                                                 ifelse(metrics_old$deletion[j] == 12 & metrics_old$run[j] == i,
                                                        brainGraph::efficiency(old_12_g, type = 'global'),
                                                        ifelse(metrics_old$deletion[j] == 16 & metrics_old$run[j] == i,
                                                               brainGraph::efficiency(old_16_g, type = 'global'),
                                                               ifelse(metrics_old$deletion[j] == 20 & metrics_old$run[j] == i,
                                                                      brainGraph::efficiency(old_20_g, type = 'global'),
                                                                      metrics_old$efficiency[j])))))
  }
}

plot(metrics_old$modularity ~ metrics_old$deletion)  # 1 y value for all sharing x
plot(metrics_old$clustering ~ metrics_old$deletion)  # all y = 1
plot(metrics_old$diameter   ~ metrics_old$deletion,  # decent
      las = 1, pch = 19, col = rgb(0,0,1,0.1))
plot(metrics_old$efficiency ~ metrics_old$deletion,  # some variation in y within x, but not much
     las = 1, pch = 19, col = rgb(0,0,1,0.1))










###calculate betweenness after deletion -- most central ####
# first ####
cent_4.1 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = central_deleted_4[,1],
                      type = sim_dyads$type)
cent_4.1 <- cent_4.1[!is.na(cent_4.1$weight),]
length(unique(cent_4.1$from))
cent_4.1_g <- igraph::graph_from_data_frame(d = cent_4.1, vertices = sim_eles, directed = F)

cent_8.1 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = central_deleted_8[,1],
                      type = sim_dyads$type)
cent_8.1 <- cent_8.1[!is.na(cent_8.1$weight),]
length(unique(cent_8.1$from))
cent_8.1_g <- igraph::graph_from_data_frame(d = cent_8.1, vertices = sim_eles, directed = F)

cent_12.1 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_12[,1],
                       type = sim_dyads$type)
cent_12.1 <- cent_12.1[!is.na(cent_12.1$weight),]
length(unique(cent_12.1$from))
cent_12.1_g <- igraph::graph_from_data_frame(d = cent_12.1, vertices = sim_eles, directed = F)

cent_16.1 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_16[,1],
                       type = sim_dyads$type)
cent_16.1 <- cent_16.1[!is.na(cent_16.1$weight),]
length(unique(cent_16.1$from))
cent_16.1_g <- igraph::graph_from_data_frame(d = cent_16.1, vertices = sim_eles, directed = F)

cent_20.1 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_20[,1],
                       type = sim_dyads$type)
cent_20.1 <- cent_20.1[!is.na(cent_20.1$weight),]
length(unique(cent_20.1$from))
cent_20.1_g <- igraph::graph_from_data_frame(d = cent_20.1, vertices = sim_eles, directed = F)

btwn_cent.1 <- data.frame(id_1 = sim_eles$id_1,
                         cent_4.1 = betweenness(cent_4.1_g, v = V(cent_4.1_g), directed = F),
                         cent_8.1 = betweenness(cent_8.1_g, v = V(cent_8.1_g), directed = F),
                         cent_12.1 = betweenness(cent_12.1_g, v = V(cent_12.1_g), directed = F),
                         cent_16.1 = betweenness(cent_16.1_g, v = V(cent_16.1_g), directed = F),
                         cent_20.1 = betweenness(cent_20.1_g, v = V(cent_20.1_g), directed = F))

colnames(btwn_cent.1) <- c('id','4','8','12','16','20')
btwn_cent.1_long <- pivot_longer(btwn_cent.1, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_cent.1_long$deletion <- as.numeric(btwn_cent.1_long$deletion)
btwn_cent.1_long$id <- as.factor(btwn_cent.1_long$id)

ggplot(data = btwn_cent.1_long, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'percentage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')

btwn_cent.1_na <- btwn_cent.1_long
btwn_cent.1_na$betweenness <- ifelse(btwn_cent.1_na$betweenness == 0, NA, btwn_cent.1_na$betweenness)

ggplot(data = btwn_cent.1_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'percentage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# Much more changing, a lot more that increase this time

# second ####
cent_4.2 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = central_deleted_4[,2],
                      type = sim_dyads$type)
cent_4.2 <- cent_4.2[!is.na(cent_4.2$weight),]
length(unique(cent_4.2$from))
cent_4.2_g <- igraph::graph_from_data_frame(d = cent_4.2, vertices = sim_eles, directed = F)

cent_8.2 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = central_deleted_8[,2],
                      type = sim_dyads$type)
cent_8.2 <- cent_8.2[!is.na(cent_8.2$weight),]
length(unique(cent_8.2$from))
cent_8.2_g <- igraph::graph_from_data_frame(d = cent_8.2, vertices = sim_eles, directed = F)

cent_12.2 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_12[,2],
                       type = sim_dyads$type)
cent_12.2 <- cent_12.2[!is.na(cent_12.2$weight),]
length(unique(cent_12.2$from))
cent_12.2_g <- igraph::graph_from_data_frame(d = cent_12.2, vertices = sim_eles, directed = F)

cent_16.2 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_16[,2],
                       type = sim_dyads$type)
cent_16.2 <- cent_16.2[!is.na(cent_16.2$weight),]
length(unique(cent_16.2$from))
cent_16.2_g <- igraph::graph_from_data_frame(d = cent_16.2, vertices = sim_eles, directed = F)

cent_20.2 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_20[,2],
                       type = sim_dyads$type)
cent_20.2 <- cent_20.2[!is.na(cent_20.2$weight),]
length(unique(cent_20.2$from))
cent_20.2_g <- igraph::graph_from_data_frame(d = cent_20.2, vertices = sim_eles, directed = F)

btwn_cent.2 <- data.frame(id_1 = sim_eles$id_1,
                         cent_4.2 = betweenness(cent_4.2_g, v = V(cent_4.2_g), directed = F),
                         cent_8.2 = betweenness(cent_8.2_g, v = V(cent_8.2_g), directed = F),
                         cent_12.2 = betweenness(cent_12.2_g, v = V(cent_12.2_g), directed = F),
                         cent_16.2 = betweenness(cent_16.2_g, v = V(cent_16.2_g), directed = F),
                         cent_20.2 = betweenness(cent_20.2_g, v = V(cent_20.2_g), directed = F))

colnames(btwn_cent.2) <- c('id','4','8','12','16','20')
btwn_cent.2_long <- pivot_longer(btwn_cent.2, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_cent.2_long$deletion <- as.numeric(btwn_cent.2_long$deletion)
btwn_cent.2_long$id <- as.factor(btwn_cent.2_long$id)
btwn_cent.2_na <- btwn_cent.2_long
btwn_cent.2_na$betweenness <- ifelse(btwn_cent.2_na$betweenness == 0, NA, btwn_cent.2_na$betweenness)

ggplot(data = btwn_cent.2_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'percentage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# more like the age-based deletions than the first set of centrality based ones, but still has more changes and more steep increases in centrality

# third ####
cent_4.3 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = central_deleted_4[,3],
                      type = sim_dyads$type)
cent_4.3 <- cent_4.3[!is.na(cent_4.3$weight),]
length(unique(cent_4.3$from))
cent_4.3_g <- igraph::graph_from_data_frame(d = cent_4.3, vertices = sim_eles, directed = F)

cent_8.3 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = central_deleted_8[,3],
                      type = sim_dyads$type)
cent_8.3 <- cent_8.3[!is.na(cent_8.3$weight),]
length(unique(cent_8.3$from))
cent_8.3_g <- igraph::graph_from_data_frame(d = cent_8.3, vertices = sim_eles, directed = F)

cent_12.3 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_12[,3],
                       type = sim_dyads$type)
cent_12.3 <- cent_12.3[!is.na(cent_12.3$weight),]
length(unique(cent_12.3$from))
cent_12.3_g <- igraph::graph_from_data_frame(d = cent_12.3, vertices = sim_eles, directed = F)

cent_16.3 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_16[,3],
                       type = sim_dyads$type)
cent_16.3 <- cent_16.3[!is.na(cent_16.3$weight),]
length(unique(cent_16.3$from))
cent_16.3_g <- igraph::graph_from_data_frame(d = cent_16.3, vertices = sim_eles, directed = F)

cent_20.3 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_20[,3],
                       type = sim_dyads$type)
cent_20.3 <- cent_20.3[!is.na(cent_20.3$weight),]
length(unique(cent_20.3$from))
cent_20.3_g <- igraph::graph_from_data_frame(d = cent_20.3, vertices = sim_eles, directed = F)

btwn_cent.3 <- data.frame(id_1 = sim_eles$id_1,
                         cent_4.3 = betweenness(cent_4.3_g, v = V(cent_4.3_g), directed = F),
                         cent_8.3 = betweenness(cent_8.3_g, v = V(cent_8.3_g), directed = F),
                         cent_12.3 = betweenness(cent_12.3_g, v = V(cent_12.3_g), directed = F),
                         cent_16.3 = betweenness(cent_16.3_g, v = V(cent_16.3_g), directed = F),
                         cent_20.3 = betweenness(cent_20.3_g, v = V(cent_20.3_g), directed = F))

colnames(btwn_cent.3) <- c('id','4','8','12','16','20')
btwn_cent.3_long <- pivot_longer(btwn_cent.3, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_cent.3_long$deletion <- as.numeric(btwn_cent.3_long$deletion)
btwn_cent.3_long$id <- as.factor(btwn_cent.3_long$id)
btwn_cent.3_na <- btwn_cent.3_long
btwn_cent.3_na$betweenness <- ifelse(btwn_cent.3_na$betweenness == 0, NA, btwn_cent.3_na$betweenness)

ggplot(data = btwn_cent.3_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'percentage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, some sharp inclines and declines

# fourth ####
cent_4.4 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = central_deleted_4[,4],
                      type = sim_dyads$type)
cent_4.4 <- cent_4.4[!is.na(cent_4.4$weight),]
length(unique(cent_4.4$from))
cent_4.4_g <- igraph::graph_from_data_frame(d = cent_4.4, vertices = sim_eles, directed = F)

cent_8.4 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = central_deleted_8[,4],
                      type = sim_dyads$type)
cent_8.4 <- cent_8.4[!is.na(cent_8.4$weight),]
length(unique(cent_8.4$from))
cent_8.4_g <- igraph::graph_from_data_frame(d = cent_8.4, vertices = sim_eles, directed = F)

cent_12.4 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_12[,4],
                       type = sim_dyads$type)
cent_12.4 <- cent_12.4[!is.na(cent_12.4$weight),]
length(unique(cent_12.4$from))
cent_12.4_g <- igraph::graph_from_data_frame(d = cent_12.4, vertices = sim_eles, directed = F)

cent_16.4 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_16[,4],
                       type = sim_dyads$type)
cent_16.4 <- cent_16.4[!is.na(cent_16.4$weight),]
length(unique(cent_16.4$from))
cent_16.4_g <- igraph::graph_from_data_frame(d = cent_16.4, vertices = sim_eles, directed = F)

cent_20.4 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_20[,4],
                       type = sim_dyads$type)
cent_20.4 <- cent_20.4[!is.na(cent_20.4$weight),]
length(unique(cent_20.4$from))
cent_20.4_g <- igraph::graph_from_data_frame(d = cent_20.4, vertices = sim_eles, directed = F)

btwn_cent.4 <- data.frame(id_1 = sim_eles$id_1,
                         cent_4.4 = betweenness(cent_4.4_g, v = V(cent_4.4_g), directed = F),
                         cent_8.4 = betweenness(cent_8.4_g, v = V(cent_8.4_g), directed = F),
                         cent_12.4 = betweenness(cent_12.4_g, v = V(cent_12.4_g), directed = F),
                         cent_16.4 = betweenness(cent_16.4_g, v = V(cent_16.4_g), directed = F),
                         cent_20.4 = betweenness(cent_20.4_g, v = V(cent_20.4_g), directed = F))

colnames(btwn_cent.4) <- c('id','4','8','12','16','20')
btwn_cent.4_long <- pivot_longer(btwn_cent.4, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_cent.4_long$deletion <- as.numeric(btwn_cent.4_long$deletion)
btwn_cent.4_long$id <- as.factor(btwn_cent.4_long$id)
btwn_cent.4_na <- btwn_cent.4_long
btwn_cent.4_na$betweenness <- ifelse(btwn_cent.4_na$betweenness == 0, NA, btwn_cent.4_na$betweenness)

ggplot(data = btwn_cent.4_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'percentage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, fewer increase than decrease

# fifth ####
cent_4.5 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = central_deleted_4[,5],
                      type = sim_dyads$type)
cent_4.5 <- cent_4.5[!is.na(cent_4.5$weight),]
length(unique(cent_4.5$from))
cent_4.5_g <- igraph::graph_from_data_frame(d = cent_4.5, vertices = sim_eles, directed = F)

cent_8.5 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = central_deleted_8[,5],
                      type = sim_dyads$type)
cent_8.5 <- cent_8.5[!is.na(cent_8.5$weight),]
length(unique(cent_8.5$from))
cent_8.5_g <- igraph::graph_from_data_frame(d = cent_8.5, vertices = sim_eles, directed = F)

cent_12.5 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_12[,5],
                       type = sim_dyads$type)
cent_12.5 <- cent_12.5[!is.na(cent_12.5$weight),]
length(unique(cent_12.5$from))
cent_12.5_g <- igraph::graph_from_data_frame(d = cent_12.5, vertices = sim_eles, directed = F)

cent_16.5 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_16[,5],
                       type = sim_dyads$type)
cent_16.5 <- cent_16.5[!is.na(cent_16.5$weight),]
length(unique(cent_16.5$from))
cent_16.5_g <- igraph::graph_from_data_frame(d = cent_16.5, vertices = sim_eles, directed = F)

cent_20.5 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_20[,5],
                       type = sim_dyads$type)
cent_20.5 <- cent_20.5[!is.na(cent_20.5$weight),]
length(unique(cent_20.5$from))
cent_20.5_g <- igraph::graph_from_data_frame(d = cent_20.5, vertices = sim_eles, directed = F)

btwn_cent.5 <- data.frame(id_1 = sim_eles$id_1,
                         cent_4.5 = betweenness(cent_4.5_g, v = V(cent_4.5_g), directed = F),
                         cent_8.5 = betweenness(cent_8.5_g, v = V(cent_8.5_g), directed = F),
                         cent_12.5 = betweenness(cent_12.5_g, v = V(cent_12.5_g), directed = F),
                         cent_16.5 = betweenness(cent_16.5_g, v = V(cent_16.5_g), directed = F),
                         cent_20.5 = betweenness(cent_20.5_g, v = V(cent_20.5_g), directed = F))

colnames(btwn_cent.5) <- c('id','4','8','12','16','20')
btwn_cent.5_long <- pivot_longer(btwn_cent.5, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_cent.5_long$deletion <- as.numeric(btwn_cent.5_long$deletion)
btwn_cent.5_long$id <- as.factor(btwn_cent.5_long$id)
btwn_cent.5_na <- btwn_cent.5_long
btwn_cent.5_na$betweenness <- ifelse(btwn_cent.5_na$betweenness == 0, NA, btwn_cent.5_na$betweenness)

ggplot(data = btwn_cent.5_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'percentage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, fewer increase than decrease

# sixth ####
cent_4.6 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = central_deleted_4[,6],
                      type = sim_dyads$type)
cent_4.6 <- cent_4.6[!is.na(cent_4.6$weight),]
length(unique(cent_4.6$from))
cent_4.6_g <- igraph::graph_from_data_frame(d = cent_4.6, vertices = sim_eles, directed = F)

cent_8.6 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = central_deleted_8[,6],
                      type = sim_dyads$type)
cent_8.6 <- cent_8.6[!is.na(cent_8.6$weight),]
length(unique(cent_8.6$from))
cent_8.6_g <- igraph::graph_from_data_frame(d = cent_8.6, vertices = sim_eles, directed = F)

cent_12.6 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_12[,6],
                       type = sim_dyads$type)
cent_12.6 <- cent_12.6[!is.na(cent_12.6$weight),]
length(unique(cent_12.6$from))
cent_12.6_g <- igraph::graph_from_data_frame(d = cent_12.6, vertices = sim_eles, directed = F)

cent_16.6 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_16[,6],
                       type = sim_dyads$type)
cent_16.6 <- cent_16.6[!is.na(cent_16.6$weight),]
length(unique(cent_16.6$from))
cent_16.6_g <- igraph::graph_from_data_frame(d = cent_16.6, vertices = sim_eles, directed = F)

cent_20.6 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_20[,6],
                       type = sim_dyads$type)
cent_20.6 <- cent_20.6[!is.na(cent_20.6$weight),]
length(unique(cent_20.6$from))
cent_20.6_g <- igraph::graph_from_data_frame(d = cent_20.6, vertices = sim_eles, directed = F)

btwn_cent.6 <- data.frame(id_1 = sim_eles$id_1,
                         cent_4.6 = betweenness(cent_4.6_g, v = V(cent_4.6_g), directed = F),
                         cent_8.6 = betweenness(cent_8.6_g, v = V(cent_8.6_g), directed = F),
                         cent_12.6 = betweenness(cent_12.6_g, v = V(cent_12.6_g), directed = F),
                         cent_16.6 = betweenness(cent_16.6_g, v = V(cent_16.6_g), directed = F),
                         cent_20.6 = betweenness(cent_20.6_g, v = V(cent_20.6_g), directed = F))

colnames(btwn_cent.6) <- c('id','4','8','12','16','20')
btwn_cent.6_long <- pivot_longer(btwn_cent.6, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_cent.6_long$deletion <- as.numeric(btwn_cent.6_long$deletion)
btwn_cent.6_long$id <- as.factor(btwn_cent.6_long$id)
btwn_cent.6_na <- btwn_cent.6_long
btwn_cent.6_na$betweenness <- ifelse(btwn_cent.6_na$betweenness == 0, NA, btwn_cent.6_na$betweenness)

ggplot(data = btwn_cent.6_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'percentage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, fewer increase than decrease

# seventh ####
cent_4.7 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = central_deleted_4[,7],
                      type = sim_dyads$type)
cent_4.7 <- cent_4.7[!is.na(cent_4.7$weight),]
length(unique(cent_4.7$from))
cent_4.7_g <- igraph::graph_from_data_frame(d = cent_4.7, vertices = sim_eles, directed = F)

cent_8.7 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = central_deleted_8[,7],
                      type = sim_dyads$type)
cent_8.7 <- cent_8.7[!is.na(cent_8.7$weight),]
length(unique(cent_8.7$from))
cent_8.7_g <- igraph::graph_from_data_frame(d = cent_8.7, vertices = sim_eles, directed = F)

cent_12.7 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_12[,7],
                       type = sim_dyads$type)
cent_12.7 <- cent_12.7[!is.na(cent_12.7$weight),]
length(unique(cent_12.7$from))
cent_12.7_g <- igraph::graph_from_data_frame(d = cent_12.7, vertices = sim_eles, directed = F)

cent_16.7 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_16[,7],
                       type = sim_dyads$type)
cent_16.7 <- cent_16.7[!is.na(cent_16.7$weight),]
length(unique(cent_16.7$from))
cent_16.7_g <- igraph::graph_from_data_frame(d = cent_16.7, vertices = sim_eles, directed = F)

cent_20.7 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_20[,7],
                       type = sim_dyads$type)
cent_20.7 <- cent_20.7[!is.na(cent_20.7$weight),]
length(unique(cent_20.7$from))
cent_20.7_g <- igraph::graph_from_data_frame(d = cent_20.7, vertices = sim_eles, directed = F)

btwn_cent.7 <- data.frame(id_1 = sim_eles$id_1,
                         cent_4.7 = betweenness(cent_4.7_g, v = V(cent_4.7_g), directed = F),
                         cent_8.7 = betweenness(cent_8.7_g, v = V(cent_8.7_g), directed = F),
                         cent_12.7 = betweenness(cent_12.7_g, v = V(cent_12.7_g), directed = F),
                         cent_16.7 = betweenness(cent_16.7_g, v = V(cent_16.7_g), directed = F),
                         cent_20.7 = betweenness(cent_20.7_g, v = V(cent_20.7_g), directed = F))

colnames(btwn_cent.7) <- c('id','4','8','12','16','20')
btwn_cent.7_long <- pivot_longer(btwn_cent.7, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_cent.7_long$deletion <- as.numeric(btwn_cent.7_long$deletion)
btwn_cent.7_long$id <- as.factor(btwn_cent.7_long$id)
btwn_cent.7_na <- btwn_cent.7_long
btwn_cent.7_na$betweenness <- ifelse(btwn_cent.7_na$betweenness == 0, NA, btwn_cent.7_na$betweenness)

ggplot(data = btwn_cent.7_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'percentage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, fewer increase than decrease

# eighth ####
cent_4.8 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = central_deleted_4[,8],
                      type = sim_dyads$type)
cent_4.8 <- cent_4.8[!is.na(cent_4.8$weight),]
length(unique(cent_4.8$from))
cent_4.8_g <- igraph::graph_from_data_frame(d = cent_4.8, vertices = sim_eles, directed = F)

cent_8.8 <- data.frame(from = sim_dyads$id_1,
                      to = sim_dyads$id_2,
                      weight = central_deleted_8[,8],
                      type = sim_dyads$type)
cent_8.8 <- cent_8.8[!is.na(cent_8.8$weight),]
length(unique(cent_8.8$from))
cent_8.8_g <- igraph::graph_from_data_frame(d = cent_8.8, vertices = sim_eles, directed = F)

cent_12.8 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_12[,8],
                       type = sim_dyads$type)
cent_12.8 <- cent_12.8[!is.na(cent_12.8$weight),]
length(unique(cent_12.8$from))
cent_12.8_g <- igraph::graph_from_data_frame(d = cent_12.8, vertices = sim_eles, directed = F)

cent_16.8 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_16[,8],
                       type = sim_dyads$type)
cent_16.8 <- cent_16.8[!is.na(cent_16.8$weight),]
length(unique(cent_16.8$from))
cent_16.8_g <- igraph::graph_from_data_frame(d = cent_16.8, vertices = sim_eles, directed = F)

cent_20.8 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = central_deleted_20[,8],
                       type = sim_dyads$type)
cent_20.8 <- cent_20.8[!is.na(cent_20.8$weight),]
length(unique(cent_20.8$from))
cent_20.8_g <- igraph::graph_from_data_frame(d = cent_20.8, vertices = sim_eles, directed = F)

btwn_cent.8 <- data.frame(id_1 = sim_eles$id_1,
                         cent_4.8 = betweenness(cent_4.8_g, v = V(cent_4.8_g), directed = F),
                         cent_8.8 = betweenness(cent_8.8_g, v = V(cent_8.8_g), directed = F),
                         cent_12.8 = betweenness(cent_12.8_g, v = V(cent_12.8_g), directed = F),
                         cent_16.8 = betweenness(cent_16.8_g, v = V(cent_16.8_g), directed = F),
                         cent_20.8 = betweenness(cent_20.8_g, v = V(cent_20.8_g), directed = F))

colnames(btwn_cent.8) <- c('id','4','8','12','16','20')
btwn_cent.8_long <- pivot_longer(btwn_cent.8, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_cent.8_long$deletion <- as.numeric(btwn_cent.8_long$deletion)
btwn_cent.8_long$id <- as.factor(btwn_cent.8_long$id)
btwn_cent.8_na <- btwn_cent.8_long
btwn_cent.8_na$betweenness <- ifelse(btwn_cent.8_na$betweenness == 0, NA, btwn_cent.8_na$betweenness)

ggplot(data = btwn_cent.8_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'percentage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, fewer increase than decrease



###calculate betweenness after deletion -- random ####
# first ####
rand_4.1 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = random_deleted_4[,1],
                       type = sim_dyads$type)
rand_4.1 <- rand_4.1[!is.na(rand_4.1$weight),]
length(unique(rand_4.1$from))
rand_4.1_g <- igraph::graph_from_data_frame(d = rand_4.1, vertices = sim_eles, directed = F)

rand_8.1 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = random_deleted_8[,1],
                       type = sim_dyads$type)
rand_8.1 <- rand_8.1[!is.na(rand_8.1$weight),]
length(unique(rand_8.1$from))
rand_8.1_g <- igraph::graph_from_data_frame(d = rand_8.1, vertices = sim_eles, directed = F)

rand_12.1 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_12[,1],
                        type = sim_dyads$type)
rand_12.1 <- rand_12.1[!is.na(rand_12.1$weight),]
length(unique(rand_12.1$from))
rand_12.1_g <- igraph::graph_from_data_frame(d = rand_12.1, vertices = sim_eles, directed = F)

rand_16.1 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_16[,1],
                        type = sim_dyads$type)
rand_16.1 <- rand_16.1[!is.na(rand_16.1$weight),]
length(unique(rand_16.1$from))
rand_16.1_g <- igraph::graph_from_data_frame(d = rand_16.1, vertices = sim_eles, directed = F)

rand_20.1 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_20[,1],
                        type = sim_dyads$type)
rand_20.1 <- rand_20.1[!is.na(rand_20.1$weight),]
length(unique(rand_20.1$from))
rand_20.1_g <- igraph::graph_from_data_frame(d = rand_20.1, vertices = sim_eles, directed = F)

btwn_rand.1 <- data.frame(id_1 = sim_eles$id_1,
                          rand_4.1 = betweenness(rand_4.1_g, v = V(rand_4.1_g), directed = F),
                          rand_8.1 = betweenness(rand_8.1_g, v = V(rand_8.1_g), directed = F),
                          rand_12.1 = betweenness(rand_12.1_g, v = V(rand_12.1_g), directed = F),
                          rand_16.1 = betweenness(rand_16.1_g, v = V(rand_16.1_g), directed = F),
                          rand_20.1 = betweenness(rand_20.1_g, v = V(rand_20.1_g), directed = F))

colnames(btwn_rand.1) <- c('id','4','8','12','16','20')
btwn_rand.1_long <- pivot_longer(btwn_rand.1, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_rand.1_long$deletion <- as.numeric(btwn_rand.1_long$deletion)
btwn_rand.1_long$id <- as.factor(btwn_rand.1_long$id)

ggplot(data = btwn_rand.1_long, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'perrandage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')

btwn_rand.1_na <- btwn_rand.1_long
btwn_rand.1_na$betweenness <- ifelse(btwn_rand.1_na$betweenness == 0, NA, btwn_rand.1_na$betweenness)

ggplot(data = btwn_rand.1_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'perrandage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# Much more changing, a lot more that increase this time

# second ####
rand_4.2 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = random_deleted_4[,2],
                       type = sim_dyads$type)
rand_4.2 <- rand_4.2[!is.na(rand_4.2$weight),]
length(unique(rand_4.2$from))
rand_4.2_g <- igraph::graph_from_data_frame(d = rand_4.2, vertices = sim_eles, directed = F)

rand_8.2 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = random_deleted_8[,2],
                       type = sim_dyads$type)
rand_8.2 <- rand_8.2[!is.na(rand_8.2$weight),]
length(unique(rand_8.2$from))
rand_8.2_g <- igraph::graph_from_data_frame(d = rand_8.2, vertices = sim_eles, directed = F)

rand_12.2 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_12[,2],
                        type = sim_dyads$type)
rand_12.2 <- rand_12.2[!is.na(rand_12.2$weight),]
length(unique(rand_12.2$from))
rand_12.2_g <- igraph::graph_from_data_frame(d = rand_12.2, vertices = sim_eles, directed = F)

rand_16.2 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_16[,2],
                        type = sim_dyads$type)
rand_16.2 <- rand_16.2[!is.na(rand_16.2$weight),]
length(unique(rand_16.2$from))
rand_16.2_g <- igraph::graph_from_data_frame(d = rand_16.2, vertices = sim_eles, directed = F)

rand_20.2 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_20[,2],
                        type = sim_dyads$type)
rand_20.2 <- rand_20.2[!is.na(rand_20.2$weight),]
length(unique(rand_20.2$from))
rand_20.2_g <- igraph::graph_from_data_frame(d = rand_20.2, vertices = sim_eles, directed = F)

btwn_rand.2 <- data.frame(id_1 = sim_eles$id_1,
                          rand_4.2 = betweenness(rand_4.2_g, v = V(rand_4.2_g), directed = F),
                          rand_8.2 = betweenness(rand_8.2_g, v = V(rand_8.2_g), directed = F),
                          rand_12.2 = betweenness(rand_12.2_g, v = V(rand_12.2_g), directed = F),
                          rand_16.2 = betweenness(rand_16.2_g, v = V(rand_16.2_g), directed = F),
                          rand_20.2 = betweenness(rand_20.2_g, v = V(rand_20.2_g), directed = F))

colnames(btwn_rand.2) <- c('id','4','8','12','16','20')
btwn_rand.2_long <- pivot_longer(btwn_rand.2, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_rand.2_long$deletion <- as.numeric(btwn_rand.2_long$deletion)
btwn_rand.2_long$id <- as.factor(btwn_rand.2_long$id)
btwn_rand.2_na <- btwn_rand.2_long
btwn_rand.2_na$betweenness <- ifelse(btwn_rand.2_na$betweenness == 0, NA, btwn_rand.2_na$betweenness)

ggplot(data = btwn_rand.2_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'perrandage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# more like the age-based deletions than the first set of randomity based ones, but still has more changes and more steep increases in randomity

# third ####
rand_4.3 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = random_deleted_4[,3],
                       type = sim_dyads$type)
rand_4.3 <- rand_4.3[!is.na(rand_4.3$weight),]
length(unique(rand_4.3$from))
rand_4.3_g <- igraph::graph_from_data_frame(d = rand_4.3, vertices = sim_eles, directed = F)

rand_8.3 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = random_deleted_8[,3],
                       type = sim_dyads$type)
rand_8.3 <- rand_8.3[!is.na(rand_8.3$weight),]
length(unique(rand_8.3$from))
rand_8.3_g <- igraph::graph_from_data_frame(d = rand_8.3, vertices = sim_eles, directed = F)

rand_12.3 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_12[,3],
                        type = sim_dyads$type)
rand_12.3 <- rand_12.3[!is.na(rand_12.3$weight),]
length(unique(rand_12.3$from))
rand_12.3_g <- igraph::graph_from_data_frame(d = rand_12.3, vertices = sim_eles, directed = F)

rand_16.3 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_16[,3],
                        type = sim_dyads$type)
rand_16.3 <- rand_16.3[!is.na(rand_16.3$weight),]
length(unique(rand_16.3$from))
rand_16.3_g <- igraph::graph_from_data_frame(d = rand_16.3, vertices = sim_eles, directed = F)

rand_20.3 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_20[,3],
                        type = sim_dyads$type)
rand_20.3 <- rand_20.3[!is.na(rand_20.3$weight),]
length(unique(rand_20.3$from))
rand_20.3_g <- igraph::graph_from_data_frame(d = rand_20.3, vertices = sim_eles, directed = F)

btwn_rand.3 <- data.frame(id_1 = sim_eles$id_1,
                          rand_4.3 = betweenness(rand_4.3_g, v = V(rand_4.3_g), directed = F),
                          rand_8.3 = betweenness(rand_8.3_g, v = V(rand_8.3_g), directed = F),
                          rand_12.3 = betweenness(rand_12.3_g, v = V(rand_12.3_g), directed = F),
                          rand_16.3 = betweenness(rand_16.3_g, v = V(rand_16.3_g), directed = F),
                          rand_20.3 = betweenness(rand_20.3_g, v = V(rand_20.3_g), directed = F))

colnames(btwn_rand.3) <- c('id','4','8','12','16','20')
btwn_rand.3_long <- pivot_longer(btwn_rand.3, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_rand.3_long$deletion <- as.numeric(btwn_rand.3_long$deletion)
btwn_rand.3_long$id <- as.factor(btwn_rand.3_long$id)
btwn_rand.3_na <- btwn_rand.3_long
btwn_rand.3_na$betweenness <- ifelse(btwn_rand.3_na$betweenness == 0, NA, btwn_rand.3_na$betweenness)

ggplot(data = btwn_rand.3_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'perrandage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, some sharp inclines and declines

# fourth ####
rand_4.4 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = random_deleted_4[,4],
                       type = sim_dyads$type)
rand_4.4 <- rand_4.4[!is.na(rand_4.4$weight),]
length(unique(rand_4.4$from))
rand_4.4_g <- igraph::graph_from_data_frame(d = rand_4.4, vertices = sim_eles, directed = F)

rand_8.4 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = random_deleted_8[,4],
                       type = sim_dyads$type)
rand_8.4 <- rand_8.4[!is.na(rand_8.4$weight),]
length(unique(rand_8.4$from))
rand_8.4_g <- igraph::graph_from_data_frame(d = rand_8.4, vertices = sim_eles, directed = F)

rand_12.4 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_12[,4],
                        type = sim_dyads$type)
rand_12.4 <- rand_12.4[!is.na(rand_12.4$weight),]
length(unique(rand_12.4$from))
rand_12.4_g <- igraph::graph_from_data_frame(d = rand_12.4, vertices = sim_eles, directed = F)

rand_16.4 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_16[,4],
                        type = sim_dyads$type)
rand_16.4 <- rand_16.4[!is.na(rand_16.4$weight),]
length(unique(rand_16.4$from))
rand_16.4_g <- igraph::graph_from_data_frame(d = rand_16.4, vertices = sim_eles, directed = F)

rand_20.4 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_20[,4],
                        type = sim_dyads$type)
rand_20.4 <- rand_20.4[!is.na(rand_20.4$weight),]
length(unique(rand_20.4$from))
rand_20.4_g <- igraph::graph_from_data_frame(d = rand_20.4, vertices = sim_eles, directed = F)

btwn_rand.4 <- data.frame(id_1 = sim_eles$id_1,
                          rand_4.4 = betweenness(rand_4.4_g, v = V(rand_4.4_g), directed = F),
                          rand_8.4 = betweenness(rand_8.4_g, v = V(rand_8.4_g), directed = F),
                          rand_12.4 = betweenness(rand_12.4_g, v = V(rand_12.4_g), directed = F),
                          rand_16.4 = betweenness(rand_16.4_g, v = V(rand_16.4_g), directed = F),
                          rand_20.4 = betweenness(rand_20.4_g, v = V(rand_20.4_g), directed = F))

colnames(btwn_rand.4) <- c('id','4','8','12','16','20')
btwn_rand.4_long <- pivot_longer(btwn_rand.4, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_rand.4_long$deletion <- as.numeric(btwn_rand.4_long$deletion)
btwn_rand.4_long$id <- as.factor(btwn_rand.4_long$id)
btwn_rand.4_na <- btwn_rand.4_long
btwn_rand.4_na$betweenness <- ifelse(btwn_rand.4_na$betweenness == 0, NA, btwn_rand.4_na$betweenness)

ggplot(data = btwn_rand.4_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'perrandage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, fewer increase than decrease

# fifth ####
rand_4.5 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = random_deleted_4[,5],
                       type = sim_dyads$type)
rand_4.5 <- rand_4.5[!is.na(rand_4.5$weight),]
length(unique(rand_4.5$from))
rand_4.5_g <- igraph::graph_from_data_frame(d = rand_4.5, vertices = sim_eles, directed = F)

rand_8.5 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = random_deleted_8[,5],
                       type = sim_dyads$type)
rand_8.5 <- rand_8.5[!is.na(rand_8.5$weight),]
length(unique(rand_8.5$from))
rand_8.5_g <- igraph::graph_from_data_frame(d = rand_8.5, vertices = sim_eles, directed = F)

rand_12.5 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_12[,5],
                        type = sim_dyads$type)
rand_12.5 <- rand_12.5[!is.na(rand_12.5$weight),]
length(unique(rand_12.5$from))
rand_12.5_g <- igraph::graph_from_data_frame(d = rand_12.5, vertices = sim_eles, directed = F)

rand_16.5 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_16[,5],
                        type = sim_dyads$type)
rand_16.5 <- rand_16.5[!is.na(rand_16.5$weight),]
length(unique(rand_16.5$from))
rand_16.5_g <- igraph::graph_from_data_frame(d = rand_16.5, vertices = sim_eles, directed = F)

rand_20.5 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_20[,5],
                        type = sim_dyads$type)
rand_20.5 <- rand_20.5[!is.na(rand_20.5$weight),]
length(unique(rand_20.5$from))
rand_20.5_g <- igraph::graph_from_data_frame(d = rand_20.5, vertices = sim_eles, directed = F)

btwn_rand.5 <- data.frame(id_1 = sim_eles$id_1,
                          rand_4.5 = betweenness(rand_4.5_g, v = V(rand_4.5_g), directed = F),
                          rand_8.5 = betweenness(rand_8.5_g, v = V(rand_8.5_g), directed = F),
                          rand_12.5 = betweenness(rand_12.5_g, v = V(rand_12.5_g), directed = F),
                          rand_16.5 = betweenness(rand_16.5_g, v = V(rand_16.5_g), directed = F),
                          rand_20.5 = betweenness(rand_20.5_g, v = V(rand_20.5_g), directed = F))

colnames(btwn_rand.5) <- c('id','4','8','12','16','20')
btwn_rand.5_long <- pivot_longer(btwn_rand.5, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_rand.5_long$deletion <- as.numeric(btwn_rand.5_long$deletion)
btwn_rand.5_long$id <- as.factor(btwn_rand.5_long$id)
btwn_rand.5_na <- btwn_rand.5_long
btwn_rand.5_na$betweenness <- ifelse(btwn_rand.5_na$betweenness == 0, NA, btwn_rand.5_na$betweenness)

ggplot(data = btwn_rand.5_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'perrandage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, fewer increase than decrease

# sixth ####
rand_4.6 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = random_deleted_4[,6],
                       type = sim_dyads$type)
rand_4.6 <- rand_4.6[!is.na(rand_4.6$weight),]
length(unique(rand_4.6$from))
rand_4.6_g <- igraph::graph_from_data_frame(d = rand_4.6, vertices = sim_eles, directed = F)

rand_8.6 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = random_deleted_8[,6],
                       type = sim_dyads$type)
rand_8.6 <- rand_8.6[!is.na(rand_8.6$weight),]
length(unique(rand_8.6$from))
rand_8.6_g <- igraph::graph_from_data_frame(d = rand_8.6, vertices = sim_eles, directed = F)

rand_12.6 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_12[,6],
                        type = sim_dyads$type)
rand_12.6 <- rand_12.6[!is.na(rand_12.6$weight),]
length(unique(rand_12.6$from))
rand_12.6_g <- igraph::graph_from_data_frame(d = rand_12.6, vertices = sim_eles, directed = F)

rand_16.6 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_16[,6],
                        type = sim_dyads$type)
rand_16.6 <- rand_16.6[!is.na(rand_16.6$weight),]
length(unique(rand_16.6$from))
rand_16.6_g <- igraph::graph_from_data_frame(d = rand_16.6, vertices = sim_eles, directed = F)

rand_20.6 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_20[,6],
                        type = sim_dyads$type)
rand_20.6 <- rand_20.6[!is.na(rand_20.6$weight),]
length(unique(rand_20.6$from))
rand_20.6_g <- igraph::graph_from_data_frame(d = rand_20.6, vertices = sim_eles, directed = F)

btwn_rand.6 <- data.frame(id_1 = sim_eles$id_1,
                          rand_4.6 = betweenness(rand_4.6_g, v = V(rand_4.6_g), directed = F),
                          rand_8.6 = betweenness(rand_8.6_g, v = V(rand_8.6_g), directed = F),
                          rand_12.6 = betweenness(rand_12.6_g, v = V(rand_12.6_g), directed = F),
                          rand_16.6 = betweenness(rand_16.6_g, v = V(rand_16.6_g), directed = F),
                          rand_20.6 = betweenness(rand_20.6_g, v = V(rand_20.6_g), directed = F))

colnames(btwn_rand.6) <- c('id','4','8','12','16','20')
btwn_rand.6_long <- pivot_longer(btwn_rand.6, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_rand.6_long$deletion <- as.numeric(btwn_rand.6_long$deletion)
btwn_rand.6_long$id <- as.factor(btwn_rand.6_long$id)
btwn_rand.6_na <- btwn_rand.6_long
btwn_rand.6_na$betweenness <- ifelse(btwn_rand.6_na$betweenness == 0, NA, btwn_rand.6_na$betweenness)

ggplot(data = btwn_rand.6_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'perrandage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, fewer increase than decrease

# seventh ####
rand_4.7 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = random_deleted_4[,7],
                       type = sim_dyads$type)
rand_4.7 <- rand_4.7[!is.na(rand_4.7$weight),]
length(unique(rand_4.7$from))
rand_4.7_g <- igraph::graph_from_data_frame(d = rand_4.7, vertices = sim_eles, directed = F)

rand_8.7 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = random_deleted_8[,7],
                       type = sim_dyads$type)
rand_8.7 <- rand_8.7[!is.na(rand_8.7$weight),]
length(unique(rand_8.7$from))
rand_8.7_g <- igraph::graph_from_data_frame(d = rand_8.7, vertices = sim_eles, directed = F)

rand_12.7 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_12[,7],
                        type = sim_dyads$type)
rand_12.7 <- rand_12.7[!is.na(rand_12.7$weight),]
length(unique(rand_12.7$from))
rand_12.7_g <- igraph::graph_from_data_frame(d = rand_12.7, vertices = sim_eles, directed = F)

rand_16.7 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_16[,7],
                        type = sim_dyads$type)
rand_16.7 <- rand_16.7[!is.na(rand_16.7$weight),]
length(unique(rand_16.7$from))
rand_16.7_g <- igraph::graph_from_data_frame(d = rand_16.7, vertices = sim_eles, directed = F)

rand_20.7 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_20[,7],
                        type = sim_dyads$type)
rand_20.7 <- rand_20.7[!is.na(rand_20.7$weight),]
length(unique(rand_20.7$from))
rand_20.7_g <- igraph::graph_from_data_frame(d = rand_20.7, vertices = sim_eles, directed = F)

btwn_rand.7 <- data.frame(id_1 = sim_eles$id_1,
                          rand_4.7 = betweenness(rand_4.7_g, v = V(rand_4.7_g), directed = F),
                          rand_8.7 = betweenness(rand_8.7_g, v = V(rand_8.7_g), directed = F),
                          rand_12.7 = betweenness(rand_12.7_g, v = V(rand_12.7_g), directed = F),
                          rand_16.7 = betweenness(rand_16.7_g, v = V(rand_16.7_g), directed = F),
                          rand_20.7 = betweenness(rand_20.7_g, v = V(rand_20.7_g), directed = F))

colnames(btwn_rand.7) <- c('id','4','8','12','16','20')
btwn_rand.7_long <- pivot_longer(btwn_rand.7, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_rand.7_long$deletion <- as.numeric(btwn_rand.7_long$deletion)
btwn_rand.7_long$id <- as.factor(btwn_rand.7_long$id)
btwn_rand.7_na <- btwn_rand.7_long
btwn_rand.7_na$betweenness <- ifelse(btwn_rand.7_na$betweenness == 0, NA, btwn_rand.7_na$betweenness)

ggplot(data = btwn_rand.7_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'perrandage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, fewer increase than decrease

# eighth ####
rand_4.8 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = random_deleted_4[,8],
                       type = sim_dyads$type)
rand_4.8 <- rand_4.8[!is.na(rand_4.8$weight),]
length(unique(rand_4.8$from))
rand_4.8_g <- igraph::graph_from_data_frame(d = rand_4.8, vertices = sim_eles, directed = F)

rand_8.8 <- data.frame(from = sim_dyads$id_1,
                       to = sim_dyads$id_2,
                       weight = random_deleted_8[,8],
                       type = sim_dyads$type)
rand_8.8 <- rand_8.8[!is.na(rand_8.8$weight),]
length(unique(rand_8.8$from))
rand_8.8_g <- igraph::graph_from_data_frame(d = rand_8.8, vertices = sim_eles, directed = F)

rand_12.8 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_12[,8],
                        type = sim_dyads$type)
rand_12.8 <- rand_12.8[!is.na(rand_12.8$weight),]
length(unique(rand_12.8$from))
rand_12.8_g <- igraph::graph_from_data_frame(d = rand_12.8, vertices = sim_eles, directed = F)

rand_16.8 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_16[,8],
                        type = sim_dyads$type)
rand_16.8 <- rand_16.8[!is.na(rand_16.8$weight),]
length(unique(rand_16.8$from))
rand_16.8_g <- igraph::graph_from_data_frame(d = rand_16.8, vertices = sim_eles, directed = F)

rand_20.8 <- data.frame(from = sim_dyads$id_1,
                        to = sim_dyads$id_2,
                        weight = random_deleted_20[,8],
                        type = sim_dyads$type)
rand_20.8 <- rand_20.8[!is.na(rand_20.8$weight),]
length(unique(rand_20.8$from))
rand_20.8_g <- igraph::graph_from_data_frame(d = rand_20.8, vertices = sim_eles, directed = F)

btwn_rand.8 <- data.frame(id_1 = sim_eles$id_1,
                          rand_4.8 = betweenness(rand_4.8_g, v = V(rand_4.8_g), directed = F),
                          rand_8.8 = betweenness(rand_8.8_g, v = V(rand_8.8_g), directed = F),
                          rand_12.8 = betweenness(rand_12.8_g, v = V(rand_12.8_g), directed = F),
                          rand_16.8 = betweenness(rand_16.8_g, v = V(rand_16.8_g), directed = F),
                          rand_20.8 = betweenness(rand_20.8_g, v = V(rand_20.8_g), directed = F))

colnames(btwn_rand.8) <- c('id','4','8','12','16','20')
btwn_rand.8_long <- pivot_longer(btwn_rand.8, cols = 2:6, names_to = 'deletion', values_to = 'betweenness')
btwn_rand.8_long$deletion <- as.numeric(btwn_rand.8_long$deletion)
btwn_rand.8_long$id <- as.factor(btwn_rand.8_long$id)
btwn_rand.8_na <- btwn_rand.8_long
btwn_rand.8_na$betweenness <- ifelse(btwn_rand.8_na$betweenness == 0, NA, btwn_rand.8_na$betweenness)

ggplot(data = btwn_rand.8_na, aes(x = deletion, y = betweenness, col = id))+
  geom_line()+
  geom_point()+
  scale_x_continuous(expand = c(0,0), 'perrandage deleted')+
  scale_y_continuous(expand = c(0,2))+
  theme_light()+
  theme(legend.position = 'none')
# majority remain fairly stable level of betweenness throughout, fewer increase than decrease






