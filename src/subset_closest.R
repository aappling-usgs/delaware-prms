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
#' reaches <- sf::read_sf('data/delaware_stream_temp_by_segment/delaware_segments/delaware_segments.shp')
#' sites <- feather::read_feather('../2wp-temp-observations/1_wqp_pull/inout/wqp_inventory.feather')
subset_closest <- function(reaches, sites){
  
  # format site points as spatial and project to wgs84, so we can subset them by
  # the reach layer
  # site_coords <- sites %>%
  #   dplyr::select(longitude, latitude) %>%
  #   as.matrix()
  # rownames(site_coords) <- sites$MonitoringLocationIdentifier
  sites_sp <- SpatialPointsDataFrame(
    coords=dplyr::select(sites, longitude, latitude),
    data=sites,
    proj4string=CRS('+proj=longlat +datum=WGS84'))
  
  # format reaches as spatial, project to wgs84, and set reach ID to be seg_id_nat
  reaches_sp <- sf::as_Spatial(reaches) %>%
    spTransform(CRS("+proj=longlat +datum=WGS84"))
  for (i in 1:nrow(reaches_sp)) {
    reaches_sp@lines[[i]]@ID <- as.character(reaches_sp@data$seg_id_nat[i])
  }
  
  # geospatially subset points by  bbox polygon in wgs84 around the reaches
  reach_bbox <- st_bbox(reaches_sp)
  reach_bbox_poly <- as(raster::extent(reach_bbox[[1]],reach_bbox[[3]],reach_bbox[[2]],reach_bbox[[4]]), "SpatialPolygons")
  crs(reach_bbox_poly) <- "+proj=longlat +datum=WGS84"
  sites_subset <- sites_sp[reach_bbox_poly,]
  
  # select the median point of each reach (approximately the midpoint).
  stopifnot(unique(sapply(1:nrow(reaches_sp), function(i) length(reaches_sp[i,]@lines))) == 1)
  stopifnot(unique(sapply(1:nrow(reaches_sp), function(i) length(reaches_sp[i,]@lines[[1]]@Lines))) == 1)
  reach_medians <-
    sapply(seq_len(nrow(reaches_sp)), function(i) {
      apply(reaches_sp[i, ]@lines[[1]]@Lines[[1]]@coords, MARGIN=2, FUN=median)
    }) %>%
    t() %>%
    as_tibble(.name_repair='unique') %>%
    rename(lon=...1, lat=...2) %>%
    SpatialPoints(proj4string = CRS(proj4string(reaches_sp)))
  
  # use K-Nearest-Neighbors to find the reach whose median is nearest each site
  tree <- createTree(coordinates(reach_medians))
  site_reach_matches <- knnLookup(tree, newx = sites_subset$longitude, newy=sites_subset$latitude, k=1)[,1]

  # attach the reach information (seg_id_nat and geometry) to the key site
  # information. site_reach_matches is ordered according to sites_subset and
  # contains indices into reach_medians (which is ordered just like reaches_sp)
  crosswalk <- reaches_sp[site_reach_matches,]
  crosswalk@data <- crosswalk@data %>%
    bind_cols(dplyr::select(sites_subset@data, MonitoringLocationIdentifier, site_lat=latitude, site_lon=longitude)) %>%
    bind_cols(dplyr::rename(as_tibble(reach_medians[site_reach_matches,]), reach_median_lat=lat, reach_median_lon=lon))
  
  # compute distance from each site to its KNN-nearest reach. projection chosen
  # using http://www.dmap.co.uk/utmworld.htm and
  # https://www.spatialreference.org/ref/epsg/wgs-84-utm-zone-18n/
  reach_medians_proj <- sf::st_as_sf(crosswalk@data, coords=c('reach_median_lon', 'reach_median_lat'), crs=4326) %>%
    sf::st_transform(crs="+proj=utm +zone=18 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  sites_proj <- sf::st_as_sf(crosswalk@data, coords=c('site_lon', 'site_lat'), crs=4326) %>%
    sf::st_transform(crs="+proj=utm +zone=18 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  # dist_to_reach_sf = st_distance(reach_medians_proj, sites_proj,
  # by_element=TRUE) is also a possibility but takes 21 seconds instead of 0.05
  # seconds below. The results are the same execept that st_distance confirms
  # that the units are meters
  crosswalk@data <- crosswalk@data %>%
    mutate(dist_to_reach_km = 0.001 * sqrt( 
      (sf::st_coordinates(reach_medians_proj)[,1] - sf::st_coordinates(sites_proj)[,1])^2 +
        (sf::st_coordinates(reach_medians_proj)[,2] - sf::st_coordinates(sites_proj)[,2])^2))
  
  # convert to sf because it's 2019 =)
  crosswalk <- sf::st_as_sf(crosswalk)
  
  return(crosswalk)
}
