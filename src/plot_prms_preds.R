library(readr)
library(dplyr)
library(tidyr)
library(sf)

# read the shape file
del_sf <- sf::read_sf('data/delaware_stream_temp_by_segment/delaware_segments/delaware_segments.shp')
st_bbox(del_sf)

# read the temperature predictions
stream_temps <- read_csv('data/delaware_stream_temp_by_segment/seg_tave_water.csv') %>%
  gather("seg_id_nat", "temperature", -date) %>% 
  mutate(
    seg_id_nat = as.integer(seg_id_nat),
    temperature = ifelse(temperature %in% c(-99.9, -98.9), NA, temperature))

# to get a feel for the predictions, extract and plot the temperature map for one day
example_day <- filter(stream_temps, date == stream_temps$date[200]) %>%
  select(seg_id_nat, temperature)
example_day_sf <- left_join(del_sf, example_day, by='seg_id_nat')
plot(example_day_sf['temperature'])

#### Map seg_id_nat to NHDPlusV2 comid ####

# Translate from seg_id_nat (unique reach IDs specific to the Geospatial Fabric)
# to POI (POI_ID is identical to the NHDPlus COMI_ID attribute).
#   The crosswalk from seg_id_nat to POI is at
# https://www.sciencebase.gov/catalog/item/5362b683e4b0c409c6289bf6 in the 2.3GB
# file named GeospatialFabric_National.gdb.zip.
#   Related metadata (not a perfect match, but a start) are at
# https://www.sciencebase.gov/catalog/item/535eda80e4b08e65d60fc834.
st_layers('data/GeospatialFabric_National.gdb') # get the layer list
gf_points <- sf::read_sf('data/GeospatialFabric_National.gdb', layer='POIs')
gf_points$poi_gage_id %>% unique()
gf_segments <- sf::read_sf('data/GeospatialFabric_National.gdb', layer='nsegmentNationalIdentifier')

names(gf_segments)
del_sf_POI <- del_sf %>%
  left_join(as.data.frame(gf_segments), by='seg_id_nat')
# nrow(filter(del_sf_POI, is.na(POI_ID))) # 0 = everything matched
length(unique(del_sf$seg_id_nat)) # 456
length(unique(del_sf_POI$POI_ID)) # 436 but 20 seg_id_nats share a POI_ID with another seg_id_nat
shared_POIs <- del_sf_POI %>%
  group_by(POI_ID) %>%
  mutate(count = n()) %>%
  filter(count > 1)
one_shared_POI <- sort(unique(shared_POIs$POI_ID))[1]
del_sf_POI %>%
  filter(POI_ID == one_shared_POI) %>%
  select(seg_id_nat) %>%
  plot()
del_sf_POI %>%
  filter(POI_ID == one_shared_POI) %>%
  select(POI_ID) %>%
  plot()
gf_segments %>% 
  filter(POI_ID == one_shared_POI) %>%
  select(seg_id_nat) %>%
  plot()

# Translate from POI to read blodgett's PRMS-NHDV2 crosswalk at
# https://www.sciencebase.gov/catalog/item/5a81e805e4b00f54eb30eb87. POI_ID and
# COMID are both the "from" points in this crosswalk (see the iso_metadata at
# the above link). "COMID is always the NHDPlusV2 COMID associated with the
# HUC_12 from pour point (FPP). POI_ID is the Geofabric (GF) identifier for the
# Point of Interest that can be used as a source of streamflow for the HU12 in
# question."
hu12_network_matches <- readr::read_tsv('data/MappingFrom12Di/hu12_network_matches.tsv')
del_sf_V2 <- left_join(del_sf_POI, hu12_network_matches, by='POI_ID') %>%
  mutate(unmatched = factor(is.na(COMID)))
nrow(filter(del_sf_V2, is.na(POI_ID))) # 0
nrow(filter(del_sf_V2, is.na(HUC_12))) # 220
nrow(filter(del_sf_V2, is.na(COMID))) # 220 non-matches
nrow(filter(del_sf_V2, !is.na(COMID))) # 279 good matches
filter(del_sf_V2, !is.na(COMID)) %>% pull(POI_ID) %>% unique() %>% length() # 226 unique good matches
all.equal(filter(del_sf_V2, is.na(COMID))$POI_ID, filter(del_sf_V2, is.na(HUC_12))$POI_ID) # same ones missing

unmatched_segs <- select(del_sf_V2, seg_id_nat, POI_ID, COMID, HUC_12) %>%
  filter(is.na(COMID))
par(mar=c(1,0,1.2,1))
plot(del_sf_V2['unmatched'])
plot(unmatched_segs['seg_id_nat'])

length(unique(del_sf_POI$POI_ID)) # 436 unique POI_IDs in delaware. every PRMS reach gets a non-NA POI_ID.
length(setdiff(del_sf_POI$POI_ID, hu12_network_matches$POI_ID)) # 210 of the delaware POI_IDs aren't present in hu12_network_matches
nrow(hu12_network_matches) # 83508 rows in the crosswalk
nrow(filter(hu12_network_matches, !is.na(POI_ID))) # 49287 of those rows have a non-NA POI_ID, i.e. are potential amtches with the delaware linework
nrow(filter(hu12_network_matches, !is.na(POI_ID), is.na(COMID))) # 0: every row that has a POI_ID also has a non-NA COMID
length(intersect(del_sf_POI$POI_ID, hu12_network_matches$POI_ID)) # 226 
length(unique(intersect(del_sf_POI$POI_ID, hu12_network_matches$POI_ID))) # 226 unique matches expected (as confirmed above)

# we can't just assume that COMID == POI_ID
hu12_network_matches %>% filter(COMID != POI_ID) # 33580 rows

#### Delaware Basin boundary ####

# Attempt to get a complete basin boundary for the Delaware using NLDI. This may
# also be available in the Geospatial Fabric
del_pour_points <- c(
  '020402070203', '020402050802', '020402050803',
  '020402060105', '020402060602', '020402070102',
  '020402060505', '020402070602', '020402070102') # , 
# stalls or very slow on '020402070405', '020402070505', '020402070304'
# rejects '020402070205', '020402060507' 'Cannot open data source'
# no PRMS matches: '020402070202'
del_subbasins <- lapply(del_pour_points, function(dpp) {
  message(dpp)
  sf::read_sf(sprintf('https://cida.usgs.gov/nldi/huc12pp/%s/basin', dpp)) %>%
    st_transform(crs=st_crs(del_prms))
})
plot(st_geometry(del_prms), col='navyblue', add=F)
del_basin <- del_subbasins[[1]]
plot(st_geometry(del_basin), col=NA, border='red', add=T)
for(dsb in del_subbasins[-1]) {
  del_basin <- st_union(del_basin, dsb)
  plot(st_geometry(dsb), col=NA, border='red', add=T)
}
plot(st_geometry(del_basin), col=NA, border='seagreen1', add=T)


#### Other NHDPlusV2 stuff ####

# Download NHDPlusV2 from
# ftp://www.horizon-systems.com/NHDPlus/NHDPlusV21/Data/NationalData/NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07.7z
#nhdplusv2 <- ...

# here's a V1-V2 crosswalk. do we need this? takes a while to load b/c 2.6 million rows, 10 cols. comes from
# the NHDPlusV2 site (http://www.horizon-systems.com/NHDPlus/V2NationalData.php)
# crosswalk_v1v2 <- foreign::read.dbf('data/NHDPlusV21_NationalData_V1_To_V2_Crosswalk_01/NHDPlusNationalData/NHDPlusV1Network_V2Network_Crosswalk.dbf')

