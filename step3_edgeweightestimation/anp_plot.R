library(igraph)

plot_network_threshold_anp <- function (edge_samples, dyad_data, lwd = 2, threshold = 0.3, 
          label.colour = "transparent", label.font = "Helvetica", node.size = 4, 
          node.colour = "seagreen1", link.colour1 = "black", link.colour2 = rgb(0, 0, 0, 0.3)) {
  dyad_name <- do.call(paste, c(dyad_data[c("node_1", "node_2")], 
                                sep = " <-> "))
  edge_lower <- apply(edge_samples, 2, function(x) quantile(x, 
                                                            probs = 0.025))
  edge_upper <- apply(edge_samples, 2, function(x) quantile(x, 
                                                            probs = 0.975))
  edge_median <- apply(edge_samples, 2, function(x) quantile(x, 
                                                             probs = 0.5))
  edge_list <- cbind(median = round(edge_median, 3), `2.5%` = round(edge_lower, 
                                                                    3), `97.5%` = round(edge_upper, 3))
  rownames(edge_list) <- dyad_name
  edgelist <- as.data.frame(edge_list)
  edgelist$node_1 <- as.character(dyad_data$node_1)
  edgelist$node_2 <- as.character(dyad_data$node_2)
  edgelist <- edgelist[, c(4:5, 1:3)]
  threshold_edges <- edgelist[edgelist$median >= threshold, 
  ]
  if (nrow(threshold_edges) == 0) {
    stop("No edges above threshold")
  }
  net <- igraph::graph_from_edgelist(as.matrix(threshold_edges[, 
                                                               1:2]), directed = F)
  if (is.data.frame(node.size) == TRUE) {
    nodes_list <- data.frame(node = rep(NA, length(unique(c(threshold_edges$node_1, 
                                                            threshold_edges$node_2)))), sightings = NA)
    for (i in 1:nrow(nodes_list)) {
      nodes_all <- rep(NA, 2 * nrow(threshold_edges))
      for (a in 1:2) {
        for (b in 1:nrow(threshold_edges)) {
          nodes_all[a + (b - 1) * 2] <- threshold_edges[b, 
                                                        a]
        }
      }
      nodes_list$node <- unique(nodes_all)
      nodes_list$sightings[i] <- nodes$sightings[which(nodes$node == 
                                                         nodes_list$node[i])]
    }
    node_sightings <- nodes_list$sightings
  }
  else {
    node_sightings <- node.size
  }
  if (is.data.frame(node.colour) == TRUE) {
    nodes_list <- data.frame(node = rep(NA, length(unique(c(threshold_edges$node_1, 
                                                            threshold_edges$node_2)))), age = NA)
    for (i in 1:nrow(nodes_list)) {
      nodes_all <- rep(NA, 2 * nrow(threshold_edges))
      for (a in 1:2) {
        for (b in 1:nrow(threshold_edges)) {
          nodes_all[a + (b - 1) * 2] <- threshold_edges[b, 
                                                        a]
        }
      }
      nodes_list$node <- unique(nodes_all)
      nodes_list$age[i] <- nodes$age[which(nodes$node == 
                                             nodes_list$node[i])]
    }
    node_age <- nodes_list$age
  }
  else {
    node_age <- node.colour
  }
  md <- threshold_edges[, 3]
  ub <- threshold_edges[, 5]
  coords <- igraph::layout_nicely(net)
  igraph::plot.igraph(net,
                      layout = coords,
                      vertex.label.color = ifelse(is.null(label.colour) == TRUE,
                                                  ifelse(node_age < 20,
                                                         "black",
                                                         "white"),
                                                  label.colour), 
                      label.family = label.font,
                      vertex.color = ifelse(node_age < 15,
                                            "#FDE725FF",
                                            ifelse(node_age < 20,
                                                   "#55C667FF",
                                                   ifelse(node_age < 30,
                                                          "#1F968BFF",
                                                          ifelse(node_age < 40,
                                                                 "#39568CFF",
                                                                 "#440154FF")))),
                      vertex.size = node_sightings, 
                      frame.color = NA,
                      frame.width = 0,
                      edge.color = NA,
                      edge.arrow.size = 0, 
                      edge.width = 0)
  igraph::plot.igraph(net,
                      layout = coords,
                      add = TRUE,
                      vertex.label = NA, 
                      vertex.color = "transparent",
                      vertex.size = 0,
                      frame.color = NA, 
                      frame.width = 0,
                      edge.color = link.colour1, 
                      edge.arrow.size = 0, 
                      edge.width = md * lwd)
  igraph::plot.igraph(net,
                      layout = coords,
                      add = TRUE,
                      vertex.label = NA, 
                      vertex.color = "transparent",
                      vertex.size = 0,
                      frame.color = NA, 
                      frame.width = 0,
                      edge.color = link.colour2,
                      edge.arrow.size = 0, 
                      edge.width = ub * lwd)
}

par(mai = c(0.2,0.2,0.2,0.2))
pdf('../outputs/anp_edgecalculations_networkplots_conditionalprior_changesize.pdf')

load('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')
plot_network_threshold_anp(edge_samples = edge_samples,
                           dyad_data = cdf_1,
                           threshold = 0.10,
                           node.size = nodes,
                           node.colour = nodes,
                           lwd = 15)
plot_network_threshold_anp(edge_samples = edge_samples,
                           dyad_data = cdf_1,
                           threshold = 0.15,
                           node.size = nodes,
                           node.colour = nodes,
                           lwd = 15)
plot_network_threshold_anp(edge_samples = edge_samples,
                           dyad_data = cdf_1,
                           threshold = 0.20,
                           node.size = nodes,
                           node.colour = nodes,
                           lwd = 15)

for(window in 2:36){
  load(paste0('anp_edgecalculations/anpshort',window,'_edgeweights_conditionalprior.RData'))
  
  for(i in c(0.1, 0.15, 0.2)){
    plot_network_threshold_anp(edge_samples = edge_samples,
                               dyad_data = cdf,
                               threshold = i,
                               node.size = nodes,
                               node.colour = nodes,
                               lwd = 15)
  }
  
  rm(list = ls()[! ls() %in% 'window']) ; gc()
}

dev.off()
