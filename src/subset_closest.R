library(raster)
library(dplyr)
library(SearchTrees)
library(sf)

#' Match each site point to an edge (reach) and its downstream vertex (outlet)
#' 
#' @details 
#' 
#' Algorithm:
#' 
#' Locate the nearest endpoint and the nearest reach, then:
#'   a) if the nearest reach is closer than the nearest point, use the nearest reach
#'      and locate the closest point on that reach, then see below
#'   b) if the nearest reach and nearest point are the same distance, use the nearest
#'      point and the closest reach upstream of that nearest point
#'   c) if the nearest reach is farther than the nearest point, that's cray-cray so 
#'      throw an error
#' 
#' In situation (a), once you have the nearest reach and the nearest point on that reach:
#'   1) if the nearest point is downstream, then use that point
#'   2) else if the nearest point is upstream, then:
#'     i) if the upstream reach does not fork, match to the upstream point and the reach
#'        it drains
#'     ii) else if the upstream reach does fork, then:
#'       A) if the upstream point is >4x closer than the downstream point, match to the
#'          upstream point and the closest reach it drains
#'       B) else match to the downstream point and the matched reach 
#' 
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
  
  # see https://gis.stackexchange.com/questions/288570/find-nearest-point-along-polyline-using-sf-package-in-r
  # for the guidance I used for implementation
  system.time({ # 1.7 seconds
    # Locate the nearest endpoint=vertex and the nearest reach (by nearest point within that reach)
    nearest_subseg <- reaches[sf::st_nearest_feature(sites, reaches), ]
    bird_dist_to_subseg_m <- st_length(st_nearest_points(sites, nearest_subseg, pairwise=TRUE))
    nearest_vertex <- vertices[sf::st_nearest_feature(sites, vertices), ]
    bird_dist_to_vertex_m <- st_length(st_nearest_points(sites, nearest_vertex, pairwise=TRUE))
    nearest <- tibble(
      site_id = sites$site_id,
      nearest_subseg = nearest_subseg, # creates a cool horizontally nested tibble structure with multiple geometries
      # nearest_subseg_id = nearest_subseg$subseg_id,
      bird_dist_to_subseg_m = bird_dist_to_subseg_m %>% units::drop_units(),
      nearest_vertex = nearest_vertex, # creates a cool horizontally nested tibble structure with multiple geometries
      # nearest_point_ids = nearest_vertex$point_ids,
      bird_dist_to_vertex_m = bird_dist_to_vertex_m %>% units::drop_units()
    ) %>%
      st_set_geometry(st_geometry(sites))
  })
  
  # For each site, use nested conditionals to navigate a set of possibilities
  # for what the best site-reach match is. See notes below
  sites_matched <- nearest %>%
    slice(1:5) %>% # for testing
    split(.$site_id) %>% # rowwise
    purrr:::map_dfr(function(site_sf) {
      message(site_sf$site_id)
      site <- tibble(site_id = site_sf$site_id) # initialize a new non-sf tibble to return
      # equidistance_tolerance_m <- 1
      # if(abs(site_sf$bird_dist_to_subseg_m - site_sf$bird_dist_to_vertex_m) < equidistance_tolerance_m) {
      #   # The site point is as close to a reach vertex as to any other point on
      #   # a reach. This occurs in 292 cases for tol=1m (and 355 cases for
      #   # tol=10m, 281 for 0.1m, 277 for 0.01m, 276 for 0.001 down to 1e-12)
      #   
      #   # use the matched point
      #   # site$point_ids <- site_sf$nearest_point_ids
      #   
      #   # use the nearest reach of those that drain to the matched point
      #   upstream_reaches <- filter(reaches, end_pt == site_sf$nearest_vertex$point_ids)
      #   closest_upstream_subseg <- upstream_reaches[st_nearest_feature(site_sf, upstream_reaches),]
      #   site$subseg_id <- closest_upstream_subseg$subseg_id
      #   
      #   # revise the distances
      #   site$bird_dist_to_subseg_m <- st_distance(site_sf, closest_upstream_subseg)[1,1] %>% units::drop_units()
      #   site$fish_dist_to_outlet_m <- equidistance_tolerance_m # this is close enough, anyway
      #   
      # } else if(site_sf$bird_dist_to_subseg_m < site_sf$bird_dist_to_vertex_m) {
        # The nearest reach is closer than the nearest point. This occurs in
        # 4732 cases for tol=1m
        
        # calculate distances to the downstream and upstream vertices of the
        # matched reach. I tried using st_split but see
        # https://gis.stackexchange.com/questions/288570/find-nearest-point-along-polyline-using-sf-package-in-r
        # -- that often doesn't actually split the line even if you snap the
        # point first
        vertex_upstream <- filter(vertices, point_ids == site_sf$nearest_subseg$start_pt)
        vertex_downstream <- filter(vertices, point_ids == site_sf$nearest_subseg$end_pt)
        subseg_as_points <- st_cast(st_geometry(site_sf$nearest_subseg), 'POINT') # would work poorly if there were a big straight reach with no intermediate points
        point_pos_in_subseg <- st_nearest_feature(site_sf, subseg_as_points)
        stopifnot(st_nearest_feature(vertex_upstream, subseg_as_points) == 1) # confirm that the points are listed upstream to downstream
        fish_dist_upstream_m <- st_combine(subseg_as_points[1:point_pos_in_subseg]) %>%
          st_cast('LINESTRING') %>% st_length() %>% units::drop_units()
        fish_dist_downstream_m <- st_combine(subseg_as_points[point_pos_in_subseg:length(subseg_as_points)]) %>%
          st_cast('LINESTRING') %>% st_length() %>% units::drop_units()
        
        # Decide which reach to use. Because the model predicts values for an
        # the downstream point of each stream reach, we will sometimes want to
        # use the reach upstream of the matched reach (if the site point was
        # very close to the upstream point). So we need some conditionals.
        if(fish_dist_downstream_m < fish_dist_upstream_m) {
          # The nearest point (by fish distance) is downstream, so use the current reach
          # site$point_ids <- vertex_downstream$point_ids
          site$subseg_id <- site_sf$nearest_subseg$subseg_id
          site$bird_dist_to_subseg_m <- site_sf$bird_dist_to_subseg_m
          site$fish_dist_to_outlet_m <- fish_dist_downstream_m
        } else {
          # The nearest point is upstream, so count the reaches immediately
          # upstream to decide what to do
          upstream_subsegs <- filter(reaches, end_pt == vertex_upstream$point_ids)
          if(nrow(upstream_subsegs) == 0) {
            # the current reach is a headwater, so use the downstream point and this reach
            # site$point_ids <- vertex_downstream$point_ids
            site$subseg_id <- site_sf$nearest_subseg$subseg_id
            site$fish_dist_to_outlet_m <- fish_dist_downstream_m
          } else if(nrow(upstream_subsegs) == 1) {
            # The upstream reach exists and does not fork, so match to the
            # upstream point and the reach it drains
            # site$point_ids <- vertex_upstream$point_ids
            site$subseg_id <- upstream_subsegs$subseg_id
            site$fish_dist_to_outlet_m <- -fish_dist_upstream_m # negative distance to indicate upstream
          } else if(nrow(upstream_subsegs) > 1) {
            # The upstream reach does fork, so compare upstream and downstream distances
            if((fish_dist_upstream_m*4) < fish_dist_downstream_m) {
              # The upstream point is >4x closer than the downstream point, so
              # match to the upstream point and the closest reach it drains
              closest_upstream_reach <- upstream_subsegs[st_nearest_feature(site_sf, upstream_subsegs),]
              # site$point_ids <- vertex_upstream$point_ids
              site$subseg_id <- closest_upstream_reach$subseg_id
              site$fish_dist_to_outlet_m <- -fish_dist_upstream_m
            } else {
              # The upstream point isn't close enough to justify using it, so
              # match to the downstream point and the matched reach
              # site$point_ids <- vertex_downstream$point_ids
              site$subseg_id <- site_sf$nearest_subseg$subseg_id
              site$fish_dist_to_outlet_m <- fish_dist_downstream_m
            }
          }
        }
        # regardless of whether we return the current or upstream subseg as the
        # best match, the as-a-bird-flies distance to the river should still be
        # the distance to the initial nearest reach, so that we can use this
        # measure to decide whether we were able to snap the site to the river
        # network successfully
        site$bird_dist_to_subseg_m <- site_sf$bird_dist_to_subseg_m
        
      # } else {
      #   stop('st_nearest_feature should never say a vertex is closer than the nearest reach; I must have misinterpreted the docs')
      # }
      
      return(select(site, everything()))
    })
  
  return(crosswalk)
}
