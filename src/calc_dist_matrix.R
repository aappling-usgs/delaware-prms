library(tidyverse)
library(sf)
library(igraph)

# Use igraph to compute a matrix of distances (in meters)
calc_dist_matrices <- function(edges) {
  graph <- edges %>%
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
  return(list(complete=dists_complete, downstream=dists_downstream, upstream=dists_upstream, updown=dists_updown))
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
    if(any(dat < 0)) {
      scale_fill_gradient2('Distance (m)', na.value='#192f41')
    } else {
      scale_fill_gradient('Distance (m)', low = "#ffffff", high = "#396a93", na.value='#192f41')
    } +
    scale_x_discrete(position = 'top') +
    xlab('End Point') + ylab('Start Point') +
    theme(panel.grid = element_blank(),
     axis.ticks = element_blank(),
     axis.text = element_blank()) +
     # axis.text = element_text(size = 9, color='grey50'),
     # axis.text.x = element_text(angle = 90, hjust = 0)) +
    coord_equal() +
    ggtitle(title)
  if(!missing(save_as)) {
    ggsave(save_as, plot=g, width=11, height=10)
  }
  return(g)
}

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

# save distance matrices in numpy format
save_dist_matrices <- function(dist_mat_list, out_file) {
  np <- reticulate::import('numpy')
  np$savez_compressed(
    file=out_file,
    downstream=dist_mat_list$downstream,
    upstream=dist_mat_list$upstream,
    complete=dist_mat_list$complete,
    updown=dist_mat_list$updown)
  # R access examples:
  # loaded <- np$load('out/dists.npz')
  # loaded$files
  # loaded$f[['updown']]
}