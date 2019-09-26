library(tidyverse)
library(sf)
library(igraph)

# Read in network from refine_prms_network.R
reach_net <- readRDS('out/reach_network.rds')

# Use igraph to compute a matrix of distances (in meters)
graph <- reach_net$edges %>%
  st_set_geometry(NULL) %>%
  select(start_pt, end_pt, subseg_id, subseg_length) %>%
  igraph::graph_from_data_frame(directed = TRUE) # df must contain "a symbolic edge list in the first two columns. Additional columns are considered as edge attributes"
dists_complete <- distances(graph, weights = edge.attributes(graph)$subseg_length, mode='all') # symmetric
dists_downstream <- distances(graph, weights = edge.attributes(graph)$subseg_length, mode='out') # shortest paths FROM each vertex.
dists_upstream <- distances(graph, weights = edge.attributes(graph)$subseg_length, mode='in') # shortest paths TO each vertex (flowing upstream only)
dists_updown <- dists_downstream
for(i in 1:nrow(dists_downstream)) {
  for(j in 1:ncol(dists_downstream)) {
    if(is.infinite(dists_downstream[i,j]) & !is.infinite(dists_upstream[i,j])) {
      dists_updown[i,j] <- -dists_upstream[i,j]
    }
  }
}

dist_heatmap <- function(dat, title, save_as) {
  dat_df <- as_tibble(dat) %>%
    mutate(start_pt=rownames(dat)) %>%
    gather('end_pt', 'dist_m', -start_pt) %>%
    mutate(
      start_pt = sapply(strsplit(start_pt, ';'), function(splits) splits[1]),
      end_pt = sapply(strsplit(end_pt, ';'), function(splits) splits[1]))
  point_levels <- unique(c(dat_df$start_pt, dat_df$end_pt))
  point_order <- order(sapply(strsplit(point_levels, '(u|d)(;)*'), function(splits) as.integer(splits[1])))
  point_levels <- point_levels[point_order]
  dat_df <- dat_df %>%
    mutate(
      start_pt = ordered(start_pt, levels=rev(point_levels)),
      end_pt = ordered(end_pt, levels=point_levels))
  g <- ggplot(dat_df, aes(y=start_pt, x=end_pt)) + 
    geom_tile(aes(fill = dist_m), color = NA) +
    scale_fill_gradient('Distance (m)', low = "#ffffff", high = "#396a93", na.value='#192f41') +
    scale_x_discrete(position = 'top') +
    xlab('End Point') + ylab('Start Point') +
    theme(panel.grid = element_blank(),
     axis.ticks = element_blank(),
     axis.text = element_text(size = 9, color='grey50'),
     axis.text.x = element_text(angle = 90, hjust = 0)) +
    coord_equal() +
    ggtitle(title)
  if(!missing(save_as)) {
    ggsave(save_as, plot=g, width=11, height=10)
  }
  return(g)
}
dist_heatmap(dists_downstream[1:100,1:100], 'Downstream', 'out/dists_downstream.png')
dist_heatmap(dists_upstream[1:100,1:100], 'Upstream', 'out/dists_upstream.png')
dist_heatmap(dists_complete[1:100,1:100], 'Complete', 'out/dists_complete.png')
dist_heatmap(dists_updown[1:100,1:100], 'Updown', 'out/dists_updown.png')

round(dists_downstream[1:6,1:6])
round(dists_upstream[1:6,1:6])
round(dists_complete[1:6,1:6])
round(dists_updown[1:6,1:6])
# example: from_pt=2u flows to to_pt=8u with reach length = pull(filter(reach_net$edges, subseg_id=='8_1'), subseg_length) = 17677 m
# dists_downstream['3d','4d'] = 1914
# dists_downstream['4','3'] = Inf
# dists_upstream['3','4'] = Inf
# dists_upstream['4','3'] = 1914
# dists_symmetric['3','4'] = dists_symmetric['4','3'] = 1914


# visualize some more more
plot_dists <- function(start_pt, dist_mat, net, title, save_as) {
  dists_from_start <- dist_mat[start_pt,]/1000
  pt_dists <- mutate(net$vertices, dist_from_start=dists_from_start[point_ids]) %>%
    filter(is.finite(dist_from_start))
  g <- ggplot(net$edges) +
    geom_sf(color='lightgrey') +
    geom_sf(data=pt_dists, aes(color=dist_from_start)) +
    geom_sf(data=filter(net$vertices, point_ids==start_pt), color='red') +
    theme_bw() +
    scale_color_continuous('Distance to\nred point (km)') +
    ggtitle(title)
  if(!missing(save_as)) {
    ggsave(save_as, plot=g, width=5, height=7)
  }
  return(g)
}
pt <- '850d;851d;853u'
plot_dists(pt, dists_downstream, reach_net, 'Downstream', 'out/map_dists_downstream.png')
plot_dists(pt, dists_upstream, reach_net, 'Upstream', 'out/map_dists_upstream.png')
plot_dists(pt, dists_complete, reach_net, 'Complete Network', 'out/map_dists_complete.png')
plot_dists(pt, dists_updown, reach_net, 'Upstream or Downstream', 'out/map_dists_updown.png')

# save distance matrices in numpy format
np <- reticulate::import('numpy')
np$savez_compressed(file='out/dists.npz', downstream=dists_downstream, upstream=dists_upstream, complete=dists_complete, updown=dists_updown)
# R access examples:
# loaded <- np$load('out/dists.npz')
# loaded$files
# loaded$f[['updown']]
