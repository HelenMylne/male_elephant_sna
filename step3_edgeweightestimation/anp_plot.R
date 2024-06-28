#### set up ####
library(igraph)
library(tidyverse)
par(mai = c(0.2,0.2,0.2,0.2))

#### plot ANP short ####
n_windows <- 36
for(time_window in 1:n_windows){
  ## load data
  load(paste0('anp_edgecalculations/anpshort',time_window,'_edgeweights_conditionalprior.RData'))
  if('cdf_1' %in% ls()){
    cdf <- cdf_1
    rm(cdf_1)
  }
  
  ## filter out any elephants <10 years old
  cdf <- cdf %>% 
    filter(age_start_1 > 9) %>% 
    filter(age_start_2 > 9)
  if('edge_weight[1]' %in% colnames(edge_samples)){
    edge_samples <- edge_samples[,cdf$dyad_rank]
  } else {
    edge_samples <- edge_samples[,as.character(cdf$dyad_id)]
  }
  nodes <- nodes %>% 
    filter(age > 9)
  
  ## overwrite plotting function currently saved in workspace
  plot_network_threshold_anp <- function (output, device = c('pdf','png'),
                                          edge_samples, dyad_data, lwd = 2, threshold = 0.3,
                                          label.colour = 'transparent', label.font = 'Helvetica',
                                          node.size = 4, node.colour = 'seagreen1',
                                          link.colour1 = 'black', link.colour2 = rgb(0, 0, 0, 0.3))
  {
    ## create edge list
    dyad_name <- do.call(paste, c(dyad_data[c("node_1", "node_2")], sep=" <-> "))
    edge_lower <- apply(edge_samples, 2, function(x) quantile(x, probs=0.025))
    edge_upper <- apply(edge_samples, 2, function(x) quantile(x, probs=0.975))
    edge_median <- apply(edge_samples, 2, function(x) quantile(x, probs=0.5))
    edge_list <- cbind(
      "median"=round(edge_median, 3),
      "2.5%"=round(edge_lower, 3),
      "97.5%"=round(edge_upper, 3)
    )
    rownames(edge_list) <- dyad_name
    edgelist <- as.data.frame(edge_list)
    edgelist$node_1 <- as.character(dyad_data$node_1)
    edgelist$node_2 <- as.character(dyad_data$node_2)
    edgelist <- edgelist[,c(4:5,1:3)]
    
    ## set threshold
    threshold_edges <- edgelist[edgelist$median >= threshold,]
    if(nrow(threshold_edges) == 0) { stop('No edges above threshold') }
    
    ## create network
    #net_all <- igraph::graph_from_edgelist(as.matrix(edgelist[, 1:2]), directed = F)
    net <- igraph::graph_from_edgelist(as.matrix(threshold_edges[, 1:2]), directed = F)
    
    ## determine node size
    if(is.data.frame(node.size) == TRUE ) {
      nodes_list <- data.frame(node = rep(NA, length(unique(c(threshold_edges$node_1, threshold_edges$node_2)))), #as.numeric(names(net_all[[1]])),
                               sightings = NA)
      for(i in 1:nrow(nodes_list)){
        nodes_all <- rep(NA, 2*nrow(threshold_edges))
        for(a in 1:2){
          for(b in 1:nrow(threshold_edges)){
            nodes_all[a + (b-1)*2] <- threshold_edges[b,a]
          }
        }
        nodes_list$node <- unique(nodes_all)
        nodes_list$sightings[i] <- nodes$sightings[which(nodes$node == nodes_list$node[i])]
      }
      node_sightings <- nodes_list$sightings #log(nodes_list$sightings)*5
    } else { node_sightings <- node.size }
    
    ## determine node colour
    if(is.data.frame(node.colour) == TRUE ) {
      nodes_list <- data.frame(node = rep(NA, length(unique(c(threshold_edges$node_1, threshold_edges$node_2)))), #as.numeric(names(net_all[[1]])),
                               age = NA)
      for(i in 1:nrow(nodes_list)){
        nodes_all <- rep(NA, 2*nrow(threshold_edges))
        for(a in 1:2){
          for(b in 1:nrow(threshold_edges)){
            nodes_all[a + (b-1)*2] <- threshold_edges[b,a]
          }
        }
        nodes_list$node <- unique(nodes_all)
        nodes_list$age[i] <- nodes$age[which(nodes$node == nodes_list$node[i])]
      }
      node_age <- nodes_list$age
    } else { node_age <- node.colour }
    
    ## identify weights for plotting: median and upper bound
    md <- threshold_edges[, 3]
    ub <- threshold_edges[, 5]
    
    ## set up graph layout
    coords <- igraph::layout_nicely(net)
    
    ## create output file
    if(device == 'pdf'){
      pdf(paste0(output,'.pdf'))
    } else {
      png(paste0(output,'.png'))
    }
    
    ## plot
    igraph::plot.igraph(net, layout = coords,
                        vertex.label.color = ifelse(is.null(label.colour) == TRUE,
                                                    ifelse(node_age < 20, 'black', 'white'),
                                                    label.colour),
                        label.family = label.font,
                        vertex.color = ifelse(node_age < 15, '#FDE725FF',
                                              ifelse(node_age < 20, '#55C667FF',
                                                     ifelse(node_age < 25, '#1F968BFF',
                                                            ifelse(node_age < 40, '#39568CFF', '#440154FF')))),
                        vertex.size = node_sightings,
                        frame.color = NA, frame.width = 0,
                        edge.color = NA, edge.arrow.size = 0, edge.width = 0)
    igraph::plot.igraph(net, layout = coords, add = TRUE,
                        vertex.label = NA, vertex.color = 'transparent', vertex.size = 0,
                        frame.color = NA, frame.width = 0,
                        edge.color = link.colour1, edge.arrow.size = 0, edge.width = md * lwd)
    igraph::plot.igraph(net, layout = coords, add = TRUE,
                        vertex.label = NA, vertex.color = 'transparent', vertex.size = 0,
                        frame.color = NA, frame.width = 0,
                        edge.color = link.colour2, edge.arrow.size = 0, edge.width = ub * lwd)
    
    ## save output
    dev.off()
  }
  
  ## plot
  plot_network_threshold_anp(output = paste0('../outputs/step3_edgeweightestimation/anpshort',
                                             time_window,
                                             '_network_0.10'),
                             device = 'pdf',
                             edge_samples = edge_samples, dyad_data = cdf, threshold = 0.10,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  plot_network_threshold_anp(output = paste0('../outputs/step3_edgeweightestimation/anpshort',
                                             time_window,
                                             '_network_0.10'),
                             device = 'png',
                             edge_samples = edge_samples, dyad_data = cdf, threshold = 0.10,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  plot_network_threshold_anp(output = paste0('../outputs/step3_edgeweightestimation/anpshort',
                                             time_window,
                                             '_network_0.15'),
                             device = 'pdf',
                             edge_samples = edge_samples, dyad_data = cdf, threshold = 0.15,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  plot_network_threshold_anp(output = paste0('../outputs/step3_edgeweightestimation/anpshort',
                                             time_window,
                                             '_network_0.15'),
                             device = 'png',
                             edge_samples = edge_samples, dyad_data = cdf, threshold = 0.15,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  
  ## clear workspace
  rm(list = ls()[!ls() %in% c('time_window','n_windows')]) ; gc()
  
  ## print progress marker
  print(paste0('time window ',time_window,' complete'))
}

#### plot ANP long  ####
n_windows <- 7
for(time_window in 1:n_windows){
  ## load data
  load(paste0('anp_edgecalculations/anplong',time_window,'_edgeweights_conditionalprior.RData'))
  if('cdf_1' %in% ls()){
    cdf <- cdf_1
  }
  
  ## filter out any elephants <10 years old
  cdf <- cdf %>% 
    filter(age_start_1 > 9) %>% 
    filter(age_start_2 > 9)
  if('edge_weight[1]' %in% colnames(edge_samples)){
    edge_samples <- edge_samples[,cdf$dyad_rank]
  } else {
    edge_samples <- edge_samples[,as.character(cdf$dyad_id)]
  }
  nodes <- nodes %>% 
    filter(age > 9)
  
  ## overwrite plotting function currently saved in workspace
  plot_network_threshold_anp <- function (output, device = c('pdf','png'),
                                          edge_samples, dyad_data, lwd = 2, threshold = 0.3,
                                          label.colour = 'transparent', label.font = 'Helvetica',
                                          node.size = 4, node.colour = 'seagreen1',
                                          link.colour1 = 'black', link.colour2 = rgb(0, 0, 0, 0.3))
  {
    ## create edge list
    dyad_name <- do.call(paste, c(dyad_data[c("node_1", "node_2")], sep=" <-> "))
    edge_lower <- apply(edge_samples, 2, function(x) quantile(x, probs=0.025))
    edge_upper <- apply(edge_samples, 2, function(x) quantile(x, probs=0.975))
    edge_median <- apply(edge_samples, 2, function(x) quantile(x, probs=0.5))
    edge_list <- cbind(
      "median"=round(edge_median, 3),
      "2.5%"=round(edge_lower, 3),
      "97.5%"=round(edge_upper, 3)
    )
    rownames(edge_list) <- dyad_name
    edgelist <- as.data.frame(edge_list)
    edgelist$node_1 <- as.character(dyad_data$node_1)
    edgelist$node_2 <- as.character(dyad_data$node_2)
    edgelist <- edgelist[,c(4:5,1:3)]
    
    ## set threshold
    threshold_edges <- edgelist[edgelist$median >= threshold,]
    if(nrow(threshold_edges) == 0) { stop('No edges above threshold') }
    
    ## create network
    #net_all <- igraph::graph_from_edgelist(as.matrix(edgelist[, 1:2]), directed = F)
    net <- igraph::graph_from_edgelist(as.matrix(threshold_edges[, 1:2]), directed = F)
    
    ## determine node size
    if(is.data.frame(node.size) == TRUE ) {
      nodes_list <- data.frame(node = rep(NA, length(unique(c(threshold_edges$node_1, threshold_edges$node_2)))), #as.numeric(names(net_all[[1]])),
                               sightings = NA)
      for(i in 1:nrow(nodes_list)){
        nodes_all <- rep(NA, 2*nrow(threshold_edges))
        for(a in 1:2){
          for(b in 1:nrow(threshold_edges)){
            nodes_all[a + (b-1)*2] <- threshold_edges[b,a]
          }
        }
        nodes_list$node <- unique(nodes_all)
        nodes_list$sightings[i] <- nodes$sightings[which(nodes$node == nodes_list$node[i])]
      }
      node_sightings <- nodes_list$sightings #log(nodes_list$sightings)*5
    } else { node_sightings <- node.size }
    
    ## determine node colour
    if(is.data.frame(node.colour) == TRUE ) {
      nodes_list <- data.frame(node = rep(NA, length(unique(c(threshold_edges$node_1, threshold_edges$node_2)))), #as.numeric(names(net_all[[1]])),
                               age = NA)
      for(i in 1:nrow(nodes_list)){
        nodes_all <- rep(NA, 2*nrow(threshold_edges))
        for(a in 1:2){
          for(b in 1:nrow(threshold_edges)){
            nodes_all[a + (b-1)*2] <- threshold_edges[b,a]
          }
        }
        nodes_list$node <- unique(nodes_all)
        nodes_list$age[i] <- nodes$age[which(nodes$node == nodes_list$node[i])]
      }
      node_age <- nodes_list$age
    } else { node_age <- node.colour }
    
    ## identify weights for plotting: median and upper bound
    md <- threshold_edges[, 3]
    ub <- threshold_edges[, 5]
    
    ## set up graph layout
    coords <- igraph::layout_nicely(net)
    
    ## create output file
    if(device == 'pdf'){
      pdf(paste0(output,'.pdf'))
    } else {
      png(paste0(output,'.png'))
    }
    
    ## plot
    igraph::plot.igraph(net, layout = coords,
                        vertex.label.color = ifelse(is.null(label.colour) == TRUE,
                                                    ifelse(node_age < 20, 'black', 'white'),
                                                    label.colour),
                        label.family = label.font,
                        vertex.color = ifelse(node_age < 15, '#FDE725FF',
                                              ifelse(node_age < 20, '#55C667FF',
                                                     ifelse(node_age < 25, '#1F968BFF',
                                                            ifelse(node_age < 40, '#39568CFF', '#440154FF')))),
                        vertex.size = node_sightings,
                        frame.color = NA, frame.width = 0,
                        edge.color = NA, edge.arrow.size = 0, edge.width = 0)
    igraph::plot.igraph(net, layout = coords, add = TRUE,
                        vertex.label = NA, vertex.color = 'transparent', vertex.size = 0,
                        frame.color = NA, frame.width = 0,
                        edge.color = link.colour1, edge.arrow.size = 0, edge.width = md * lwd)
    igraph::plot.igraph(net, layout = coords, add = TRUE,
                        vertex.label = NA, vertex.color = 'transparent', vertex.size = 0,
                        frame.color = NA, frame.width = 0,
                        edge.color = link.colour2, edge.arrow.size = 0, edge.width = ub * lwd)
    
    ## save output
    dev.off()
  }
  
  ## plot
  plot_network_threshold_anp(output = paste0('../outputs/step3_edgeweightestimation/anplong',
                                             time_window,
                                             '_network_0.10'),
                             device = 'pdf',
                             edge_samples = edge_samples, dyad_data = cdf, threshold = 0.10,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  plot_network_threshold_anp(output = paste0('../outputs/step3_edgeweightestimation/anplong',
                                             time_window,
                                             '_network_0.10'),
                             device = 'png',
                             edge_samples = edge_samples, dyad_data = cdf, threshold = 0.10,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  plot_network_threshold_anp(output = paste0('../outputs/step3_edgeweightestimation/anplong',
                                             time_window,
                                             '_network_0.15'),
                             device = 'pdf',
                             edge_samples = edge_samples, dyad_data = cdf, threshold = 0.15,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  plot_network_threshold_anp(output = paste0('../outputs/step3_edgeweightestimation/anplong',
                                             time_window,
                                             '_network_0.15'),
                             device = 'png',
                             edge_samples = edge_samples, dyad_data = cdf, threshold = 0.15,
                             node.size = nodes, node.colour = nodes, lwd = 15)
  
  ## clear workspace
  rm(list = ls()[!ls() %in% c('time_window','n_windows')]) ; gc()
  
  ## print progress marker
  print(time_window)
}

#### plot MOTNP ####
## load data
load('motnp_edgeweights_conditionalprior.RData')

## overwrite plotting function currently saved in workspace
plot_network_threshold_mot <- function (output, device = c('pdf','png'),
                                        edge_samples, dyad_data, lwd = 2, threshold = 0.3,
                                        label.colour = 'transparent', label.font = 'Helvetica',
                                        node.size = 4, node.colour = 'seagreen1',
                                        link.colour1 = 'black', link.colour2 = rgb(0, 0, 0, 0.3))
{
  ## create edge list
  dyad_name <- do.call(paste, c(dyad_data[c("node_1", "node_2")], sep=" <-> "))
  edge_lower <- apply(edge_samples, 2, function(x) quantile(x, probs=0.025))
  edge_upper <- apply(edge_samples, 2, function(x) quantile(x, probs=0.975))
  edge_median <- apply(edge_samples, 2, function(x) quantile(x, probs=0.5))
  edge_list <- cbind(
    "median"=round(edge_median, 3),
    "2.5%"=round(edge_lower, 3),
    "97.5%"=round(edge_upper, 3)
  )
  rownames(edge_list) <- dyad_name
  edgelist <- as.data.frame(edge_list)
  edgelist$node_1 <- as.character(dyad_data$node_1)
  edgelist$node_2 <- as.character(dyad_data$node_2)
  edgelist <- edgelist[,c(4:5,1:3)]

  ## set threshold
  threshold_edges <- edgelist[edgelist$median >= threshold,]
  if(nrow(threshold_edges) == 0) { stop('No edges above threshold') }

  ## create network
  #net_all <- igraph::graph_from_edgelist(as.matrix(edgelist[, 1:2]), directed = F)
  net <- igraph::graph_from_edgelist(as.matrix(threshold_edges[, 1:2]), directed = F)

  ## determine node size
  if(is.data.frame(node.size) == TRUE ) {
    nodes_list <- data.frame(node = rep(NA, length(unique(c(threshold_edges$node_1, threshold_edges$node_2)))), #as.numeric(names(net_all[[1]])),
                             sightings = NA)
    for(i in 1:nrow(nodes_list)){
      nodes_all <- rep(NA, 2*nrow(threshold_edges))
      for(a in 1:2){
        for(b in 1:nrow(threshold_edges)){
          nodes_all[a + (b-1)*2] <- threshold_edges[b,a]
        }
      }
      nodes_list$node <- unique(nodes_all)
      nodes_list$sightings[i] <- nodes$sightings[which(nodes$node == nodes_list$node[i])]
    }
    node_sightings <- nodes_list$sightings #log(nodes_list$sightings)*5
  } else { node_sightings <- node.size }

  ## determine node colour
  if(is.data.frame(node.colour) == TRUE ) {
    nodes_list <- data.frame(node = rep(NA, length(unique(c(threshold_edges$node_1, threshold_edges$node_2)))), #as.numeric(names(net_all[[1]])),
                             age = NA)
    for(i in 1:nrow(nodes_list)){
      nodes_all <- rep(NA, 2*nrow(threshold_edges))
      for(a in 1:2){
        for(b in 1:nrow(threshold_edges)){
          nodes_all[a + (b-1)*2] <- threshold_edges[b,a]
        }
      }
      nodes_list$node <- unique(nodes_all)
      nodes_list$age[i] <- nodes$age[which(nodes$node == nodes_list$node[i])]
    }
    node_age <- nodes_list$age
  } else { node_age <- node.colour }

  ## identify weights for plotting: median and upper bound
  md <- threshold_edges[, 3]
  ub <- threshold_edges[, 5]

  ## set up graph layout
  coords <- igraph::layout_nicely(net)

  ## create output file
  if(device == 'pdf'){
    pdf(paste0(output,'.pdf'))
  } else {
    png(paste0(output,'.png'))
  }

  ## plot
  igraph::plot.igraph(net, layout = coords,
                      vertex.label.color = ifelse(is.null(label.colour) == TRUE,
                                                  ifelse(node_age < 20, 'black', 'white'),
                                                  label.colour),
                      label.family = label.font,
                      vertex.color = ifelse(node_age < 15, '#FDE725FF',
                                            ifelse(node_age < 20, '#55C667FF',
                                                   ifelse(node_age < 25, '#1F968BFF',
                                                          ifelse(node_age < 40, '#39568CFF', '#440154FF')))),
                      vertex.size = node_sightings,
                      frame.color = NA, frame.width = 0,
                      edge.color = NA, edge.arrow.size = 0, edge.width = 0)
  igraph::plot.igraph(net, layout = coords, add = TRUE,
                      vertex.label = NA, vertex.color = 'transparent', vertex.size = 0,
                      frame.color = NA, frame.width = 0,
                      edge.color = link.colour1, edge.arrow.size = 0, edge.width = md * lwd)
  igraph::plot.igraph(net, layout = coords, add = TRUE,
                      vertex.label = NA, vertex.color = 'transparent', vertex.size = 0,
                      frame.color = NA, frame.width = 0,
                      edge.color = link.colour2, edge.arrow.size = 0, edge.width = ub * lwd)

  ## save output
  dev.off()
}

## plot
plot_network_threshold_mot(output = '../outputs/step3_edgeweightestimation/motnp_network_0.10',
                           device = 'pdf',
                           edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.10,
                           node.size = nodes, node.colour = nodes, lwd = 15)
plot_network_threshold_mot(output = '../outputs/step3_edgeweightestimation/motnp_network_0.10',
                           device = 'png',
                           edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.10,
                           node.size = nodes, node.colour = nodes, lwd = 15)
plot_network_threshold_mot(output = '../outputs/step3_edgeweightestimation/motnp_network_0.15',
                           device = 'pdf',
                           edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.15,
                           node.size = nodes, node.colour = nodes, lwd = 15)
plot_network_threshold_mot(output = '../outputs/step3_edgeweightestimation/motnp_network_0.15',
                           device = 'png',
                           edge_samples = edge_samples, dyad_data = counts_df, threshold = 0.15,
                           node.size = nodes, node.colour = nodes, lwd = 15)
