library(geojsonio)
library(geojsonsf)

crosswalk_site_reach <- readRDS('out/crosswalk_site_reach.rds')

delaware_pts <- crosswalk_site_reach %>%
  dplyr::mutate(
    dist_to_reach_km = bird_dist_to_subseg_m/1000,
    nobsBin = base::cut(n_obs, breaks=c(min(n_obs), 10, 100, 1000, max(n_obs)+1), right=FALSE)) %>%
  dplyr::select(
    site_id,
    dist_to_reach_km,
    matched_reach_id = seg_id_nat,
    n_obs,
    nobsBin) %>%
  sf::st_transform(crs = 4326)

delaware_geojson <- geojsonsf::sf_geojson(delaware_pts)
geojsonio::geojson_write(delaware_geojson, file = 'out/delaware_sites_summary.geojson', convert_wgs84 = TRUE)

# From here, email to David and Aaron for now (and get set up to push to S3 eventually)
