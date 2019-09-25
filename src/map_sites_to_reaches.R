library(tidyverse)
library(sf)
library(lwgeom)
library(smoothr)

#### Load site information ####

# Read in the revised network and boundary
drb_net <- readRDS('out/reach_network.rds')
drb_boundary <- sf::read_sf('out/drb_boundary/drb_boundary.shp')

# Read in the inventories of temperature observation sites
data_snapshot_path <- 'data/190912-2wp-temp-observations'
wqp_partitions <- feather::read_feather(file.path(data_snapshot_path, '1_wqp_pull/inout/wqp_pull_partitions.feather'))
wqp_sites <- feather::read_feather(file.path(data_snapshot_path, '1_wqp_pull/inout/wqp_inventory.feather')) %>%
  dplyr::left_join(wqp_partitions, by='MonitoringLocationIdentifier') %>%
  dplyr::filter(ResolvedMonitoringLocationTypeName %in% c('Spring', 'Stream')) %>%
  dplyr::mutate(source = 'wqp', file = file.path(file.path(data_snapshot_path, '1_wqp_pull/out'), sprintf('%s.feather', PullTask))) 
nwis_dv_sites <- feather::read_feather(file.path(data_snapshot_path, '2_nwis_pull/inout/nwis_dv_inventory.feather')) %>%
  dplyr::filter(site_tp_cd %in% 'ST') %>%
  dplyr::mutate(source = 'nwis-dv', site_id = paste0('USGS-', site_no), file = file.path(data_snapshot_path, '2_nwis_pull/out/nwis_dv_data.rds'))
nwis_uv_sites <- feather::read_feather(file.path(data_snapshot_path, '2_nwis_pull/inout/nwis_uv_inventory_reduced.feather')) %>% # 'reduced' excludes those in dv
  dplyr::filter(site_tp_cd %in% 'ST') %>%
  dplyr::mutate(source = 'nwis-uv', site_id = paste0('USGS-', site_no), file = file.path(data_snapshot_path, '2_nwis_pull/out/nwis_uv_data.rds'))

# Join all inventories into one. there are about 500 sites that are redundant across databases
obs_sites <- dplyr::bind_rows(
  wqp_sites %>%
    dplyr::select(source, org = OrganizationIdentifier, site_id = MonitoringLocationIdentifier, lat = latitude, lon = longitude, n_obs = resultCount, file),
  nwis_dv_sites %>%
    dplyr::select(source, org = agency_cd, site_id, lat = dec_lat_va, lon = dec_long_va, n_obs = count_nu, file),
  nwis_uv_sites %>%
    dplyr::select(source, org = agency_cd, site_id, lat = dec_lat_va, lon = dec_long_va, n_obs = count_nu, file)) %>%
  distinct() # distinct() takes it from 287,847 to 287,807 rows

# Convert to sfc
obs_site_points <- purrr::map2(obs_sites$lon, obs_sites$lat, function(lat, lon) {
  st_point(c(lat, lon), dim='XY')
})
obs_sites <- obs_sites %>%
  st_set_geometry(st_sfc(obs_site_points)) %>%
  st_set_crs(4326) %>%
  st_transform(crs=st_crs(drb_net$vertices))

# Subset to sites in the Delaware River Basin
drb_sites <- obs_sites[st_intersects(drb_boundary, obs_sites)[[1]], ] # 5009 rows

#### Match reaches to sites ####

source('src/subset_closest.R')
system.time({ # 97 seconds
  crosswalk <- subset_closest(sites=drb_sites, reaches=drb_net$edges, vertices=drb_net$vertices)
  # warns: site USGS-01433005 has diverse coordinates across databases, with bbox diagonal = 6.261 m
  # (the algorithm will have matched to the first version of that site, probably the wqp USGS-NY one)
})

# add geospatial info (drb_sites) and seg_id_nat (drb_net$edges)
crosswalk_site_reach <- left_join(drb_sites, crosswalk, by='site_id') %>%
  left_join(st_drop_geometry(select(drb_net$edges, subseg_id, seg_id_nat)), by='subseg_id')
saveRDS(crosswalk_site_reach, 'out/crosswalk_site_reach.rds')

#### Explore the site lists ####

# Some funkiness with "distinct":
length(unique(obs_site_points)) # 272299
nrow(distinct(obs_sites)) # 272299
nrow(obs_sites) # 287807
nrow(distinct(st_drop_geometry(obs_sites))) # 287807 -- distinct(sf) must operate on geometry only, not data
length(unique(obs_sites$site_id)) # 284919 - there are duplicate site_ids even after calling distinct()

# Explore the national site list
nrow(obs_sites) # 287807
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
nrow(drb_sites) # 5009
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

#### Explore the crosswalk ####

nrow(crosswalk_site_reach) # 5009 total sites

proposed_filters <- crosswalk_site_reach %>%
  mutate(
    bird_filtered = bird_dist_to_subseg_m < 250,
    fish_filtered = abs(fish_dist_to_outlet_m) < 5000) %>%
  st_drop_geometry()
table(select(proposed_filters, bird_filtered, fish_filtered))
#              fish_filtered
# bird_filtered FALSE TRUE
#         FALSE  1269 1460
#         TRUE    825 1455
# both filters would drop us to 1455 sites

# Distance from site point to reach

ggplot(crosswalk_site_reach, aes(x=bird_dist_to_subseg_m)) + geom_density(fill='gold') + theme_bw()

ggplot(drb_net$edges) + geom_sf(color='lightgray') +
  geom_sf(data=crosswalk_site_reach, aes(color=bird_dist_to_subseg_m < 250), size=0.5) +
  theme_bw() +
  ggtitle('All sites; showing proposed bird distance filter')

ggplot(drb_net$edges) + geom_sf(color='lightgray') +
  geom_sf(data=filter(crosswalk_site_reach, bird_dist_to_subseg_m < 250) %>%
            mutate(fish_dist_to_outlet_km = fish_dist_to_outlet_m/1000),
          aes(color=bird_dist_to_subseg_m), size=0.7) +
  scale_color_gradient() +
  theme_bw() +
  ggtitle('Filtered by bird distance; showing bird distance')

# Fish distance to outlet

ggplot(mutate(crosswalk_site_reach, fish_dist_to_outlet_km=fish_dist_to_outlet_m/1000), aes(x=fish_dist_to_outlet_km)) + geom_density(fill='slateblue') + theme_bw()

ggplot(drb_net$edges) + geom_sf(color='lightgray') +
  geom_sf(data=filter(crosswalk_site_reach, bird_dist_to_subseg_m < 250) %>%
            mutate(fish_dist_to_outlet_km = fish_dist_to_outlet_m/1000),
          aes(color=fish_dist_to_outlet_km), size=0.7) +
  scale_color_gradient2(low='#ff7f00', mid='#984ea3', high='#4daf4a') +
  theme_bw() +
  ggtitle('Filtered by bird distance; showing fish distance')

ggplot(drb_net$edges) + geom_sf(color='lightgray') +
  geom_sf(data=filter(crosswalk_site_reach, bird_dist_to_subseg_m < 250) %>%
            mutate(fish_dist_to_outlet_km = fish_dist_to_outlet_m/1000),
          aes(color=abs(fish_dist_to_outlet_km) < 5), size=0.7) +
  theme_bw() +
  ggtitle('Filtered by bird distance; showing proposed fish distance filter')

ggplot(drb_net$edges) + geom_sf(color='lightgray') +
  geom_sf(data=filter(crosswalk_site_reach, bird_dist_to_subseg_m < 250, abs(fish_dist_to_outlet_m) < 5000) %>%
            mutate(fish_dist_to_outlet_km = fish_dist_to_outlet_m/1000),
          aes(color=fish_dist_to_outlet_km), size=0.7) +
  scale_color_gradient2(low='#ff7f00', mid='#984ea3', high='#4daf4a') +
  theme_bw() +
  ggtitle('Filtered by bird and fish distance; showing fish distance')

# Number of observations

binned_crosswalk <- filter(crosswalk_site_reach, bird_dist_to_subseg_m < 250, abs(fish_dist_to_outlet_m) < 5000) %>%
  mutate(n_obs_bin = base::cut(n_obs, breaks=c(min(n_obs), 10, 100, 1000, max(n_obs)), right=FALSE, labels=c('1-10', '10-100', '100-1000', '1000+')))
ggplot(drb_net$edges) + geom_sf(color='lightgray') +
  geom_sf(data=binned_crosswalk, aes(color=n_obs_bin), size=1) +
  scale_color_brewer('Number of Observations', palette=3) +
  theme_bw() +
  ggtitle('Filtered by bird and fish distance; showing observation counts')

ggplot(binned_crosswalk, aes(n_obs)) + stat_ecdf() + theme_bw() + scale_x_log10() +
  scale_y_continuous(breaks = seq(0,1,by=0.2), labels=seq(0,1,by=0.2)*nrow(binned_crosswalk)) +
  xlab('Number of observations') + ylab('Number of sites with >=X observations')

# Check on site matching process site by site

site <- crosswalk_site_reach %>% filter(bird_dist_to_subseg_m < 250) %>% slice(300)
nearby_reaches <- st_crop(drb_net$edges, st_bbox(st_buffer(site, dist=units::set_units(30000, m))))
ggplot(st_buffer(site, dist=units::set_units(250, m))) + geom_sf() +
  geom_sf(data=site, aes(color=TRUE)) +
  geom_sf(data=nearby_reaches, aes(color=subseg_id == site$subseg_id)) +
  scale_color_discrete('Matched site and reach') +
  theme_bw()

