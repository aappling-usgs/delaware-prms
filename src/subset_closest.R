library(raster)
library(dplyr)
library(SearchTrees)
library(sf)

#' @param reaches sf Simple Feature Collection with reach identifiers in
#'   seg_id_nat column
#' @param sites inventory data.frame with MonitoringLocationIdentifier,
#'   latitude, and longitude at a minimum
#' @return data.frame of site info: MonitoringLocationIdentifier, seg_id_nat,
#'   lat, lon, and distance from site lat/lon to centroid of the reach segment
#' @examples
#' reaches <- readRDS('out/reach_network.rds')$edges
#' vertices <- readRDS('out/reach_network.rds')$vertices
#' sites <- drb_sites # build drb_sites in map_sites_to_reaches.R
subset_closest <- function(reaches, sites){
  
  # # narrow the search space by doing a simple "within distance" query for each site
  # system.time({ # 22 secs for all 5024 sites and 459 reaches in DRB
  #   # use the max segment length as an upper bound on acceptable distance between site and reach
  #   neighbors <- sf::st_is_within_distance(sites, reaches, max(reaches$subseg_length))
  # })
  # n_neighbors <- sapply(neighbors, length) # yields between 6 and 114 matches per site
  # summary(n_neighbors)
  
  # system.time({ # 22 secs for all 5024 sites and 459 reaches in DRB
  #   dists <- sf::st_distance(sites, reaches)
  #   # nngeo::st_nn() is also an option - it returns nearest neighbors - but it
  #   # uses sf::st_distance and requires an additional package, so skip it
  # })
  # site_reach_matches <- sites %>%
  #   mutate(
  #     min_dist_id = apply(dists, MARGIN=1, function(dists_to_reach) {
  #       which.min(dists_to_reach)
  #     }),
  #     nearest_subseg = reaches$subseg_id[min_dist_id],
  #     dist_to_subseg = dists[matrix(c(1:n(), min_dist_id), ncol=2, byrow=FALSE)]
  #   ) %>%
  #   select(-min_dist_id) # not useful now that we have nearest_subseg adn dist_to_subseg
  # plot(density(site_reach_matches$dist_to_subseg))
  
  system.time({ # 1.3 seconds. see https://gis.stackexchange.com/questions/288570/find-nearest-point-along-polyline-using-sf-package-in-r
    # Locate nearest feature (by nearest point within that feature), preferring
    # the feature whose downstream point is closest
    nearest_subseg_i <- sf::st_nearest_feature(sites, reaches)
    nearest_subseg <- reaches[nearest_subseg_i, ]
    
    nearest_vertex <- vertices[sf::st_nearest_feature(sites, vertices), ]
    bird_dist_to_vertex <- st_length(st_nearest_points(sites, nearest_vertex, pairwise=TRUE))
    
    # Locate nearest point along nearest reach. disagrees with st_distance
    # method for 2648 sites, but all differences are < 0.000000001 m
    shortest_line <- st_nearest_points(sites, nearest_subseg, pairwise=TRUE)
    sites$bird_dist_to_subseg <- st_length(shortest_line)
    
    sites$subseg_id <- nearest_subseg$subseg_id
    
    # Locate the nearest endpoint and the nearest reach, then:
    #   a) if the nearest reach is closer than the nearest point, use the nearest reach and locate the closest point on that reach
    #   b) if the nearest reach and nearest point are the same distance, use the closest reach upstream of that nearest point
    #   c) if the nearest reach is farther than the nearest point, that's cray-cray so throw an error
    # In situation (a), once you have the nearest reach and the nearest point on that reach:
    #   1) if the nearest point is downstream, then use that point
    #   2) else if the nearest point is upstream, then:
    #     i) if the upstream reach does not fork, match to the upstream point and the reach it drains
    #     ii) else if the upstream reach does fork, then:
    #       A) if the upstream point is >4x closer than the downstream point, match to the upstream point and the closest reach it drains
    #       B) else match to the downstream point and the matched reach
    
  })
  
  system.time({ # 37 seconds
    # Locate the endpoint of the matched reach
    sites$subseg_endpoint <- vertices[match(nearest_subseg$end_pt, vertices$point_ids), ]$point_ids
    # Locate the point snapped to the reach and split the reach at that point to
    # find the distance from the snapped point to the reach's downstream end
    snapped_point <- purrr::map_dfr(seq_len(nrow(sites)), function(i) {
      point_on_line <- st_cast(shortest_line[i], 'POINT')[2] # very close to the reach line
      point_snapped_to_line <- st_snap(point_on_line, nearest_subseg[i,], tol=1e-9) # even closer to the reach line
      
      subseg_splits <- st_split(nearest_subseg[i,], point_snapped_to_line)
      subsubseg_to_outlet <- subseg_splits[st_nearest_feature(subseg_endpoint, subseg_splits)]
      fish_dist_to_outlet <- st_length(subsubseg_to_outlet)
      tibble(sites[i,], )
    })
  })
  
  sl <- st_sf(shortest_line) %>%
    mutate(site_id = sites$site_id)
  sites_2reaches <- filter(srm, nearest_subseg == nearest_subseg2) %>% slice(14)
  dists_1s2r <- filter(sl, site_id %in% sites_2reaches$site_id)
  reaches2_sites <- filter(reaches, subseg_id %in% c(sites_2reaches$nearest_subseg, sites_2reaches$nearest_subseg2))
  ggplot(reaches2_sites) + geom_sf(aes(color=subseg_id)) +
    geom_sf(data=sites_2reaches, aes(color=nearest_subseg)) +
    geom_sf(data=dists_1s2r, aes(color=site_id)) +
    theme_bw()
  
  # # select the median point of each reach (approximately the midpoint).
  # stopifnot(unique(sapply(1:nrow(reaches_sp), function(i) length(reaches_sp[i,]@lines))) == 1)
  # stopifnot(unique(sapply(1:nrow(reaches_sp), function(i) length(reaches_sp[i,]@lines[[1]]@Lines))) == 1)
  # reach_medians <-
  #   sapply(seq_len(nrow(reaches_sp)), function(i) {
  #     apply(reaches_sp[i, ]@lines[[1]]@Lines[[1]]@coords, MARGIN=2, FUN=median)
  #   }) %>%
  #   t() %>%
  #   as_tibble(.name_repair='unique') %>%
  #   rename(lon=...1, lat=...2) %>%
  #   SpatialPoints(proj4string = CRS(proj4string(reaches_sp)))
  # 
  # # use K-Nearest-Neighbors to find the reach whose median is nearest each site
  # tree <- createTree(coordinates(reach_medians))
  # site_reach_matches <- knnLookup(tree, newx = sites_subset$longitude, newy=sites_subset$latitude, k=1)[,1]
  # 
  # # attach the reach information (seg_id_nat and geometry) to the key site
  # # information. site_reach_matches is ordered according to sites_subset and
  # # contains indices into reach_medians (which is ordered just like reaches_sp)
  # crosswalk <- reaches_sp[site_reach_matches,]
  # crosswalk@data <- crosswalk@data %>%
  #   bind_cols(dplyr::select(sites_subset@data, MonitoringLocationIdentifier, site_lat=latitude, site_lon=longitude)) %>%
  #   bind_cols(dplyr::rename(as_tibble(reach_medians[site_reach_matches,]), reach_median_lat=lat, reach_median_lon=lon))
  # 
  # # compute distance from each site to its KNN-nearest reach. projection chosen
  # # using http://www.dmap.co.uk/utmworld.htm and
  # # https://www.spatialreference.org/ref/epsg/wgs-84-utm-zone-18n/
  # reach_medians_proj <- sf::st_as_sf(crosswalk@data, coords=c('reach_median_lon', 'reach_median_lat'), crs=4326) %>%
  #   sf::st_transform(crs="+proj=utm +zone=18 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  # sites_proj <- sf::st_as_sf(crosswalk@data, coords=c('site_lon', 'site_lat'), crs=4326) %>%
  #   sf::st_transform(crs="+proj=utm +zone=18 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  # # dist_to_reach_sf = st_distance(reach_medians_proj, sites_proj,
  # # by_element=TRUE) is also a possibility but takes 21 seconds instead of 0.05
  # # seconds below. The results are the same execept that st_distance confirms
  # # that the units are meters
  # crosswalk@data <- crosswalk@data %>%
  #   mutate(dist_to_reach_km = 0.001 * sqrt( 
  #     (sf::st_coordinates(reach_medians_proj)[,1] - sf::st_coordinates(sites_proj)[,1])^2 +
  #       (sf::st_coordinates(reach_medians_proj)[,2] - sf::st_coordinates(sites_proj)[,2])^2))
  # 
  # # convert to sf because it's 2019 =)
  # crosswalk <- sf::st_as_sf(crosswalk)
  # 
  return(crosswalk)
}
