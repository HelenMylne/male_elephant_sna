load('mpnp_edgecalculations/mpnpshort1_bisonr_edgescalculated.RData')

x <- data.frame(object = ls(),
                size = NA)
for (object in ls()[!ls() %in% x]) {
  x$size[object] <- object.size(get(object))
}

print(x)

for (thing in ls()) {
  print(
    object.size(get(thing)),
    units='auto')
}


sizes <- c(5.7,30,11.5,280,52.9,7.1,5.7,16.7,7.2,168,8,280,5.1,74.7,22.4,8.4,17.5,35.2,1.2,280,14.3,56)

units <- c('Mb','Mb','Mb','bytes','Gb','Mb','Mb','Gb','Mb','bytes','Kb','bytes','Kb','Mb','Gb','Gb','Kb','Kb','Kb','bytes','Mb','bytes')

objects <- c('counts','counts_df','counts_df_model','draw98','draws','dyads','edgelist','ew_chain','ew_edgelist','filename','global_cv','median98','model_averaging','mpnp_ages','mpnp_edge_weights','mpnp_edges_null','plot_network_threshold','plot_network_threshold2','priors','sri98','summary','time_window')

x <- cbind(objects, sizes, units)
