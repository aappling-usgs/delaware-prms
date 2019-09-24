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
    bird_dist_to_subseg <- st_length(st_nearest_points(sites, nearest_subseg, pairwise=TRUE))
    nearest_vertex <- vertices[sf::st_nearest_feature(sites, vertices), ]
    bird_dist_to_vertex <- st_length(st_nearest_points(sites, nearest_vertex, pairwise=TRUE))
    nearest <- tibble(
      site_id = sites$site_id,
      nearest_subseg_id = nearest_subseg$subseg_id,
      bird_dist_to_subseg = bird_dist_to_subseg,
      nearest_point_ids = nearest_vertex$point_ids,
      bird_dist_to_vertex = bird_dist_to_vertex
    ) %>%
      st_set_geometry(st_geometry(sites))
  })
  
  #' In situation (a), once you have the nearest reach and the nearest point on that reach:
  #'   1) if the nearest point is downstream, then use that point
  #'   2) else if the nearest point is upstream, then:
  #'     i) if the upstream reach does not fork, match to the upstream point and the reach
  #'        it drains
  #'     ii) else if the upstream reach does fork, then:
  #'       A) if the upstream point is >4x closer than the downstream point, match to the
  #'          upstream point and the closest reach it drains
  #'       B) else match to the downstream point and the matched reach 
  sites_matched <- nearest %>%
    mutate(subseg_id = NA, point_ids = NA) %>%
    slice(1:5) %>%
    split(.$site_id) %>%
    map_dfr(function(site) {
      message(site$site_id)
      if(abs(site$bird_dist_to_subseg - site$bird_dist_to_vertex) < units::set_units(1, 'm')) {
        # The site point is as close to a reach vertex as to any other point on
        # a reach. This occurs in 292 cases for tol=1m (and 355 cases for
        # tol=10m, 281 for 0.1m, 277 for 0.01m, 276 for 0.001 down to 1e-12)
        
        # use the matched point
        site$point_ids <- site$nearest_point_ids
        
        # use the nearest reach of those that drain to the matched point
        upstream_reaches <- filter(reaches, end_pt == site$point_ids)
        closest_upstream_subseg <- upstream_reaches[st_nearest_feature(site, upstream_reaches),]
        site$subseg_id <- closest_upstream_subseg$subseg_id
        
        # revise the distances
        site$bird_dist_to_subseg <- st_distance(site, closest_upstream_subseg)
        site$fish_dist_to_outlet <- units::set_units(0, 'm') # or close enough, anyway
        
      } else if(site$bird_dist_to_subseg < site$bird_dist_to_vertex) {
        # The nearest reach is closer than the nearest point. This occurs in
        # 4732 cases for tol=1m
        
        # calculate distances to the downstream and upstream vertices of the
        # matched reach. I tried using st_split but see
        # https://gis.stackexchange.com/questions/288570/find-nearest-point-along-polyline-using-sf-package-in-r
        # -- that often doesn't actually split the line if you have to snap the
        # point first
        nearest_subseg <- filter(reaches, subseg_id == site$nearest_subseg_id)
        vertex_upstream <- filter(vertices, point_ids == nearest_subseg$start_pt)
        vertex_downstream <- filter(vertices, point_ids == nearest_subseg$end_pt)
        subseg_as_points <- st_cast(st_geometry(nearest_subseg), 'POINT') # would work poorly if there were a big straight reach with no intermediate points
        point_pos_in_subseg <- st_nearest_feature(site, subseg_as_points)
        stopifnot(st_nearest_feature(vertex_upstream, subseg_as_points) == 1) # confirm that the points are listed upstream to downstream
        fish_dist_upstream <- st_length(st_cast(st_combine(subseg_as_points[1:point_pos_in_subseg]), 'LINESTRING'))
        fish_dist_downstream <- st_length(st_cast(st_combine(subseg_as_points[point_pos_in_subseg:length(subseg_as_points)]), 'LINESTRING'))
        
        # Decide which reach to use. Because the model predicts values for an
        # the downstream point of each stream reach, we will sometimes want to
        # use the reach upstream of the matched reach (if the site point was
        # very close to the upstream point). So we need some conditionals.
        if(fish_dist_downstream < fish_dist_upstream) {
          # The nearest point (by fish distance) is downstream, so use that
          # point and this reach
          site$point_ids <- vertex_downstream$point_ids
          site$subseg_id <- nearest_subseg$subseg_id
        } else {
          # The nearest point is upstream, so count the reaches immediately
          # upstream to decide what to do
          upstream_subsegs <- filter(reaches, end_pt == vertex_upstream$point_ids)
          if(nrow(upstream_subsegs == 0)) {
            # the current reach is a headwater, so use the downstream point and this reach
            site$point_ids <- vertex_downstream$point_ids
            site$subseg_id <- nearest_subseg$subseg_id
            site$fish_dist_to_outlet <- units::set_units(fish_dist_downstream, m)
          } else if(nrow(upstream_subsegs) == 1) {
            # The upstream reach exists and does not fork, so match to the
            # upstream point and the reach it drains
            site$point_ids <- vertex_upstream$point_ids
            site$subseg_id <- upstream_subsegs$subseg_id
            site$fish_dist_to_outlet <- units::set_units(-fish_dist_upstream, m)
          } else if(nrow(upstream_subsegs) > 1) {
            # The upstream reach does fork, so compare upstream and downstream distances
            if((fish_dist_upstream*4) < fish_dist_downstream) {
              # The upstream point is >4x closer than the downstream point, so
              # match to the upstream point and the closest reach it drains
              closest_upstream_reach <- upstream_subsegs[st_nearest_feature(site, upstream_subsegs),]
              site$point_ids <- vertex_upstream$point_ids
              site$subseg_id <- closest_upstream_reach$subseg_id
              site$fish_dist_to_outlet <- units::set_units(-fish_dist_upstream, m)
            } else {
              # The upstream point isn't close enough to justify using it, so
              # match to the downstream point and the matched reach
              site$point_ids <- vertex_downstream$point_ids
              site$subseg_id <- nearest_subseg$subseg_id
              site$fish_dist_to_outlet <- units::set_units(fish_dist_downstream, m)
            }
          }
        }
        # calculate the as-a-bird-flies distance to the final selected reach
        selected_reach <- filter(reaches, subseg_id == site$subseg_id)
        site$bird_dist_to_subseg <- st_distance(site, selected_reach)
        
      } else {
        browser()
        stop('st_nearest_feature should never say a vertex is closer than the nearest reach; I must have misinterpreted the docs')
      }
    })
  
  return(crosswalk)
}
