library(tidyverse)
library(sf)
library(lwgeom)
library(smoothr)

#### Load site information ####

# Read in the revised network and boundary
drb_net <- readRDS('out/reach_network.rds')
drb_boundary <- sf::read_sf('out/drb_boundary/drb_boundary.shp')

# Read in the inventories of temperature observation sites
wqp_partitions <- feather::read_feather('../2wp-temp-observations/1_wqp_pull/inout/wqp_pull_partitions.feather')
wqp_sites <- feather::read_feather('data/190912-2wp-temp-observations/1_wqp_pull/inout/wqp_inventory.feather') %>%
  left_join(wqp_partitions, by='MonitoringLocationIdentifier') %>%
  filter(ResolvedMonitoringLocationTypeName %in% c('Spring', 'Stream'))
nwis_dv_sites <- feather::read_feather('data/190912-2wp-temp-observations/2_nwis_pull/inout/nwis_dv_inventory.feather')
# nwis_uv_sites <- feather::read_feather('data/190912-2wp-temp-observations/2_nwis_pull/inout/nwis_uv_inventory.feather')
nwis_uv_sites <- feather::read_feather('data/190912-2wp-temp-observations/2_nwis_pull/inout/nwis_uv_inventory_reduced.feather') # exclude those in dv

# Join all inventories into one. there are about 500 sites that are redundant across databases
obs_sites <- bind_rows(
  wqp_sites %>%
    mutate(database = 'wqp', file = file.path('data/190912-2wp-temp-observations/1_wqp_pull/out', sprintf('%s.feather', PullTask))) %>%
    select(database, org = OrganizationIdentifier, site_id = MonitoringLocationIdentifier, lat = latitude, lon = longitude, file),
  nwis_dv_sites %>%
    mutate(database = 'nwis-dv', site_id = paste0('USGS-', site_no), file = 'data/190912-2wp-temp-observations/2_nwis_pull/out/nwis_dv_data.rds') %>%
    select(database, org = agency_cd, site_id, lat = dec_lat_va, lon = dec_long_va, file),
  nwis_uv_sites %>%
    mutate(database = 'nwis-uv', site_id = paste0('USGS-', site_no), file = 'data/190912-2wp-temp-observations/2_nwis_pull/out/nwis_uv_data.rds') %>%
    select(database, org = agency_cd, site_id, lat = dec_lat_va, lon = dec_long_va, file)) %>%
  distinct() %>% # distinct() takes it from 289203 to 288731 rows
  filter(!is.na(lat) & !is.na(lon))  # removes one site, USGS-424317071263001

# Convert to sfc
obs_site_points <- purrr::map2(obs_sites$lon, obs_sites$lat, function(lat, lon) {
  st_point(c(lat, lon), dim='XY')
})
obs_sites <- obs_sites %>%
  st_set_geometry(st_sfc(obs_site_points)) %>%
  st_set_crs(4326) %>%
  st_transform(crs=st_crs(drb_net$vertices))

# Some funkiness with "distinct":
# length(unique(obs_site_points)) # 273226
# nrow(obs_sites) # 288730
# nrow(distinct(obs_sites)) # 273226
# nrow(distinct(st_drop_geometry(obs_sites))) # 288730 -- distinct(sf) must operate on geometry only, not data
# length(unique(obs_sites$site_id)) # 285910 - there are duplicate site_ids even after calling distinct()

# Subset to sites in the Delaware River Basin
drb_sites <- obs_sites[st_intersects(drb_boundary, obs_sites)[[1]], ]

#### Explore the site lists ####

# Explore the national site list
nrow(obs_sites) # 288731
table(obs_sites$database) # 3874=dv, 917=uv, 283940=wqp
multi_db_sites <- obs_sites %>%
  arrange(database) %>%
  group_by(site_id) %>%
  summarize(n_dbs = n(), dbs = paste(database, collapse=','), same_lat = length(unique(lat)) == 1, same_lon = length(unique(lon)) == 1) %>%
  filter(n_dbs > 1)
nrow(multi_db_sites) # 2820 sites appear in multiple databases
nrow(filter(multi_db_sites, !same_lat | !same_lon)) # 1644 sites with mismatches (generally minor) in lat/lon across databases
table(multi_db_sites$dbs) # 2515=nwis-dv,wqp, 305=nwis-uv,wqp: sites are either in dv or uv, not both; any overlap is between nwis and wqp

# Explore the DRB site list
nrow(drb_sites) # 5024
table(drb_sites$database) # 145=dv, 33=uv, 4846=wqp
multi_db_sites <- drb_sites %>%
  arrange(database) %>%
  group_by(site_id) %>%
  summarize(n_dbs = n(), dbs = paste(database, collapse=','), same_lat = length(unique(lat)) == 1, same_lon = length(unique(lon)) == 1) %>%
  filter(n_dbs > 1)
multi_db_sites$lat_lon_dist <- st_cast(multi_db_sites, 'MULTIPOINT') %>% st_cast('LINESTRING') %>% st_length()
nrow(multi_db_sites) # 146 sites appear in multiple databases
nrow(filter(multi_db_sites, !same_lat | !same_lon)) # 84 sites with lat-lon mismatches in lat/lon across databases
max(filter(multi_db_sites, !same_lat | !same_lon) %>% pull(lat_lon_dist)) # lat-lon mismatches are minor
table(multi_db_sites$dbs) # 2515=nwis-dv,wqp, 305=nwis-uv,wqp: sites are either in dv or uv, not both; any overlap is between nwis and wqp

ggplot(drb_net$edges) + geom_sf(color='lightgray') +
  geom_sf(data=filter(drb_sites, database=='wqp'), size=0.5, aes(color=database)) +
  geom_sf(data=filter(drb_sites, database!='wqp'), size=0.5, aes(color=database)) +
  theme_bw()

#### Match reaches to sites ####

source('src/subset_closest.R')
system.time({ # 97 seconds
  crosswalk <- subset_closest(sites=drb_sites, reaches=drb_net$edges, vertices=drb_net$vertices)
  # warns: site USGS-01433005 has diverse coordinates across databases, with bbox diagonal = 6.261 m
  # (the algorithm will have matched to the first version of that site, probably the wqp USGS-NY one)
})

crosswalk_site_reach <- left_join(drb_sites, crosswalk, by='site_id')
saveRDS(crosswalk_site_reach, 'out/crosswalk_site_reach.rds')

