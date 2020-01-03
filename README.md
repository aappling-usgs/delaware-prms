## Run order

```
make_drb_boundary -> map_sites_to_reaches ->|
                                            |-> umn_basin_subsets -> pull_all_data -> sent_to_makerspace
                      refine_prms_network ->|
```

## Scripts
```
make_drb_boundary.R
  inputs: 'data/delaware_stream_temp_by_segment/delaware_segments/delaware_segments.shp', 'data/GeospatialFabric_National.gdb'
  output: 'out/drb_boundary/drb_boundary.shp'

map_sites_to_reaches.R
  inputs: 'out/network_full.rds', 'out/drb_boundary/drb_boundary.shp', 'data/190912-2wp-temp-observations'
  sources: 'src/subset_closest.R'
  output: 'out/crosswalk_site_reach.rds'

refine_prms_network.R
  effect: Weed out some funky connections among PRMS reaches
  inputs: 
    - 'data/delaware_stream_temp_by_segment/delaware_segments/delaware_segments.shp'
    - 'data/GeospatialFabric_National.gdb'
  outputs: 'out/network_full.rds'

umn_basin_subsets.R
  effect: Create a network subset and distance matrices for the full and subset network. Also make plots.
  inputs: 'out/network_full.rds', 'out/crosswalk_site_reach.rds'
  sources: 'src/calc_dist_matrix.R'
  outputs:
    - 'out/dists_full.npz'
    - 'out/network_subset.rds'
    - 'out/dists_subset.npz'
    - 'out/map_sites_subset.png'
    - 'dists_*.png'
    - 'map_dists_*.png'
  
pull_all_data.R
  inputs:
    - 'out/crosswalk_site_reach.rds'
    - 'out/network_full.rds'
    - 'out/network_subset.rds'
    - 'data/190926-2wp-temp-observations'
    - 'data/drb_discharge_daily_dv.csv'
  outputs: 'out/obs_flow_full.csv', 'out/obs_flow_subset.csv'

send_to_makerspace.R
  effect: Write the observation site information to geojson
  inputs: 'out/crosswalk_site_reach.rds'
  outputs: 'out/delaware_sites_summary.geojson'
```

## Function definition files
```
calc_dist_matrix.R
  - calc_dist_matrices
  - dist_heatmap
  - plot_dists
  - save_dist_matrices

subset_closest.R
  - subset_closest: Match each site point to a PRMS reach
```
