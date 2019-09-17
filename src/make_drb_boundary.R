library(tidyverse)
library(sf)
library(lwgeom)
library(smoothr)

# Read in the Delaware PRMS-stream_temp stream network and the GF
del_prms <- sf::read_sf('data/delaware_stream_temp_by_segment/delaware_segments/delaware_segments.shp')
sf::st_layers('data/GeospatialFabric_National.gdb') # POIs, one, nhdflowline_en, nhdflowline, regionOutletDA, nhruNationalIdentifier, nsegmentNationalIdentifier
gf_reaches <- sf::read_sf('data/GeospatialFabric_National.gdb', layer='nsegmentNationalIdentifier') %>%
  filter(seg_id_nat %in% del_prms$seg_id_nat) # confirmed: gf_reaches$Shape_Length == round(sf::st_length(gf_reaches))
gf_catchments <- sf::read_sf('data/GeospatialFabric_National.gdb', layer='nhruNationalIdentifier') %>%
  filter(POI_ID %in% gf_reaches$POI_ID) %>%
  lwgeom::st_make_valid()

# Create, plot, and save a boundary of the full Delaware River Basin (DRB)
gf_boundary <- st_union(gf_catchments) %>%
  smoothr::fill_holes(threshold = units::set_units(100, km^2))
ggplot(gf_boundary) + geom_sf() + theme_bw()
dir.create('out/drb_boundary', recursive=FALSE, showWarnings=FALSE)
sf::write_sf(gf_boundary, 'out/drb_boundary/drb_boundary.shp')