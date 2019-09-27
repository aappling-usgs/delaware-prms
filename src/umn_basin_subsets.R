library(tidyverse)
library(sf)
source('src/calc_dist_matrix.R')

# Read in network from refine_prms_network.R
reach_net <- readRDS('out/network_full.rds')
# Read in crosswalk to get site names and observation counts (used for subnetworks)
crosswalk_site_reach <- readRDS('out/crosswalk_site_reach.rds')
filtered_crosswalk <- crosswalk_site_reach %>%
  filter(bird_dist_to_subseg_m < 250, abs(fish_dist_to_outlet_m) < 5000) %>%
  mutate(n_obs_bin = base::cut(n_obs, breaks=c(min(n_obs), 10, 100, 1000, max(n_obs)), right=FALSE, labels=c('1-10', '10-100', '100-1000', '1000+')))

#### Full network ####

# Create, save, and explore full-network distance matrices and maps
dists <- calc_dist_matrices(reach_net$edges)
save_dist_matrices(dists, 'out/dists_full.npz')

dist_heatmap(dists$downstream, 'Downstream', 'out/dists_full_downstream.png')
dist_heatmap(dists$upstream, 'Upstream', 'out/dists_full_upstream.png')
dist_heatmap(dists$complete, 'Complete', 'out/dists_full_complete.png')
dist_heatmap(dists$updown, 'Updown', 'out/dists_full_updown.png')

# matrix subset tables: first 6 points
round(dists$downstream[1:6,1:6])
round(dists$upstream[1:6,1:6])
round(dists$complete[1:6,1:6])
round(dists$updown[1:6,1:6])

# single-reach example: from_pt=3d;4u flows to to_pt=4d;5u;9d along subseg 4_1
filter(reach_net$edges, subseg_id=='4_1') # reach length = 1914 m
dists$downstream['3d;4u','4d;5u;9d'] # 1914
dists$downstream['4d;5u;9d','3d;4u'] # Inf
dists$upstream['3d;4u','4d;5u;9d'] # Inf
dists$upstream['4d;5u;9d','3d;4u'] # 1914
dists$updown['3d;4u','4d;5u;9d'] == -dists$updown['4d;5u;9d','3d;4u']

pt <- '850d;851d;853u'
plot_dists(pt, dists$downstream, reach_net, 'Downstream', 'out/map_dists_full_downstream.png')
plot_dists(pt, dists$upstream, reach_net, 'Upstream', 'out/map_dists_full_upstream.png')
plot_dists(pt, dists$complete, reach_net, 'Complete Network', 'out/map_dists_full_complete.png')
plot_dists(pt, dists$updown, reach_net, 'Upstream or Downstream', 'out/map_dists_full_updown.png')

#### ~100-edge subnetwork ####

make_subnetwork <- function(lower_point, exclude_points, drb_net, dists) {
  
  up_from_lowermost <- names(which(dists$upstream[lower_point,] < Inf))
  if(length(exclude_points) > 0) {
    up_from_uppermost <- unlist(lapply(exclude_points, function(exclude_point) {
      names(which(dists$upstream[exclude_point,] < Inf))
    }))
    subnet_points <- setdiff(up_from_lowermost, up_from_uppermost)
  } else {
    subnet_points <- up_from_lowermost
  }
  subnet_reaches <- drb_net$edges %>% filter(end_pt %in% subnet_points & ifelse(is.na(start_pt), TRUE, start_pt %in% subnet_points))
  subnet_points <- filter(drb_net$vertices, point_ids %in% subnet_points)
  
  return(list(edges=subnet_reaches, vertices=subnet_points, lower_point=lower_point, exclude_points=exclude_points))
}
explore_subnetwork <- function(subnet, crosswalk, drb_net) {
  g <- ggplot(drb_net$edges) + geom_sf(color='lightgray') +
    geom_sf(data=subnet$edges, color='seagreen') +
    geom_sf(data=crosswalk, aes(color=n_obs_bin), size=1) +
    geom_sf(data=filter(drb_net$vertices, point_ids %in% subnet$exclude_points), shape=4, color='red') +
    scale_color_brewer('Number of Observations', palette=3) +
    theme_bw() +
    ggtitle('Filtered by bird and fish distance; showing observation counts')
  print(g)
  
  message(sprintf('%d edges, %d vertices', nrow(subnet$edges), nrow(subnet$vertices)))
  
  obs_count <- subnet$edges %>%
    st_drop_geometry() %>%
    left_join(st_drop_geometry(crosswalk), by='subseg_id') %>%
    group_by(subseg_id) %>%
    summarize(n_sites=length(which(!is.na(site_id))), n_obs=sum(n_obs))
  message(sprintf('%d observed reaches, %d observations total', length(which(obs_count$n_sites > 0)), sum(obs_count$n_obs, na.rm=TRUE)))
}
# subnet_big <- make_subnetwork(lower_seg='2752_1', exclude_segs=c('332_1','341_1'), reach_net) # 165 edges, 164 vertices
subnet1 <- make_subnetwork(exclude_points=c("332u;64d;65d", "330d;341u"), lower_point='2752d;2757u', reach_net, dists) # 165 edges, 164 vertices
explore_subnetwork(subnet1, filtered_crosswalk, reach_net)
# 173 edges, 174 vertices
# 130 observed reaches, 111503 observations total

subnet2 <- make_subnetwork(lower_point='886d;895u', exclude_points=c(), reach_net, dists) # 56 edges, 57 vertices
explore_subnetwork(subnet2, filtered_crosswalk, reach_net)
# 56 edges, 57 vertices
# 40 observed reaches, 14528 observations total

subnet3 <- make_subnetwork(lower_point='285u;288d;289d', exclude_points=c(), reach_net, dists) # 36 edges, 37 vertices
explore_subnetwork(subnet3, filtered_crosswalk, reach_net)
# 36 edges, 37 vertices
# 17 observed reaches, 37999 observations total

subnet4 <- make_subnetwork(lower_point='2748u;606d;613d', exclude_points=c(), reach_net, dists) # 41 edges, 42 vertices
explore_subnetwork(subnet4, filtered_crosswalk, reach_net)
# 41 edges, 42 vertices
# 32 observed reaches, 74733 observations total

#### Save selected subnet ####

selected_subnet <- subnet4
saveRDS(selected_subnet, 'out/network_subset.rds')

dists_subnet <- calc_dist_matrices(selected_subnet$edges)
save_dist_matrices(dists_subnet, 'out/dists_subset.npz')

plot_dists('571d;572d;574u', dists_subnet$downstream, selected_subnet, 'Subnetwork - Downstream', 'out/map_dists_subset_downstream.png')
plot_dists('571d;572d;574u', dists_subnet$upstream, selected_subnet, 'Subnetwork - Upstream', 'out/map_dists_subset_upstream.png')
dist_heatmap(dists_subnet$downstream, 'Downstream', 'out/dists_subset_downstream.png')
dist_heatmap(dists_subnet$upstream, 'Upstream', 'out/dists_subset_upstream.png')
dist_heatmap(dists_subnet$updown, 'Upstream and Downstream', 'out/dists_subset_updown.png')

plot_subnet <- function(subnet, reach_net, crosswalk, out_file) {
  just_beyonds <- reach_net$edges %>% 
    filter(start_pt %in% subnet$vertices$point_ids | end_pt %in% subnet$vertices) %>%
    filter(!subseg_id %in% subnet$edges$subseg_id)
  g <- ggplot(subnet$edges) + geom_sf(color='navy') +
    geom_sf(data=just_beyonds, color='gray') +
    geom_sf(data=subnet$vertices, color='navy', shape=4, size=2) +
    geom_sf(data=filter(crosswalk, subseg_id %in% subnet$edges$subseg_id), aes(color=n_obs_bin), size=1) +
    scale_color_brewer('Number of Observations', palette=3) +
    ggtitle('Filtered by bird and fish distance; showing observation counts')
  ggsave(out_file, g, width=7, height=6)
  return(g)
}
plot_subnet(selected_subnet, reach_net, filtered_crosswalk, 'out/map_sites_subset.png')
