library(tidyverse)
library(sf)
library(igraph)

#### Calculate distances ####

# Read in the Delaware PRMS-stream_temp stream network and the GF
del_prms <- sf::read_sf('data/delaware_stream_temp_by_segment/delaware_segments/delaware_segments.shp')
#sf::st_layers('data/GeospatialFabric_National.gdb') # POIs, one, nhdflowline_en, nhdflowline, regionOutletDA, nhruNationalIdentifier, nsegmentNationalIdentifier
gf_reaches <- sf::read_sf('data/GeospatialFabric_National.gdb', layer='nsegmentNationalIdentifier') %>%
  filter(seg_id_nat %in% del_prms$seg_id_nat) %>%
  mutate(to_seg = ifelse(tosegment == 0, NA, tosegment)) %>% # tosegment is the downstream segment to which water flows next. replace 0 (none next) with NA
  select(seg_id_nat, seg_id, to_seg, region)
gf_points <- sf::read_sf('data/GeospatialFabric_National.gdb', layer='POIs') %>%
  filter(poi_gage_segment %in% gf_reaches$seg_id & region == '02') %>%
  select(GNIS_ID, GNIS_NAME, REACHCODE, seg_id = poi_gage_segment)

# inspect. lessons:
# (1) everything is in region 02, which we used above in filtering gf_points
gf_reaches %>% pull(region) %>% table()
# (2a) reaches usually have a point with the same id at their downstream end
# (2b) reaches don't necessarily connect at their endpoints
ggplot(filter(gf_reaches, seg_id %in% c(1,2,8,3))) + geom_sf(aes(color=factor(seg_id))) +
  geom_sf(data=filter(gf_points, seg_id %in% c(1,2,8,3)), aes(color=factor(seg_id))) +
  theme_bw()
filter(gf_reaches, seg_id %in% c(1,2,3,8) | to_seg %in% c(1,2,3,8)) %>% select(seg_id, to_seg)
# (3) reaches don't necessarily have gf_points at their endpoints (155,156) and
# can even have mismatches between where the reach ends and where its
# corresponding end point is (157)
ggplot(filter(mutate(gf_reaches, highlight = seg_id %in% c(155,156,157)),highlight)) + geom_sf(aes(color=factor(seg_id))) +
  geom_sf(data=filter(mutate(gf_points, highlight = seg_id %in% c(153,154,155,156,157)), highlight), aes(color=factor(seg_id))) +
  theme_bw()
filter(gf_reaches, seg_id %in% c(155,156,157)) %>% select(seg_id, to_seg)
filter(gf_points, seg_id %in% c(155,156,157))
# (4) reach 357 claims to enter reach 359, but the lines sure don't make that look plausible
ggplot(filter(mutate(gf_reaches, highlight = seg_id %in% c(355,357,359)),seg_id > 300, seg_id < 400)) + geom_sf(aes(color=factor(ifelse(highlight, seg_id, NA)))) + theme_bw()
ggplot(filter(mutate(gf_reaches, highlight = seg_id %in% c(355,357,359)),highlight)) + geom_sf(aes(color=factor(seg_id))) +
  geom_sf(data=filter(mutate(gf_points, highlight = seg_id %in% c(355,357,359)), highlight), aes(color=factor(seg_id))) +
  theme_bw()
filter(gf_reaches, seg_id == 357)
ggplot(filter(gf_reaches, seg_id %in% c(355,357,359)))
ggplot(filter(reaches_bounded, seg_id %in% c(355, 357, 359))) + geom_sf(aes(color=factor(seg_id))) + geom_sf(data=reach_points, aes(color=factor(pt_seg)))

# Modify reaches according to lessons learned
gf_reaches[gf_reaches$seg_id == 357, 'to_seg'] <- NA # there's a huge gap between the two reaches

# Augment gf_reaches with all end points of all reaches. I initially tried to
# find points in gf_points, but at least one doesn't exist (the outlet to 155)
# and I realized what I wanted most was the end points anyway.
reaches_bounded <- gf_reaches %>% mutate(
  end_points = lapply(seg_id, function(segid) {
    reach <- filter(gf_reaches, seg_id==segid)
    end_points <- reach[1,] %>% {st_cast(sf::st_geometry(.), "POINT")} %>% {c(head(.,1),tail(.,1))}
    return(end_points)
  }))
reaches_bounded <- reaches_bounded %>%
  mutate(
    which_end_up = sapply(seg_id, function(segid) {
      reach <- filter(reaches_bounded, seg_id == segid)
      if(!is.na(reach$to_seg)) { # non-river-mouths
        # upstream endpoint is the one furthest from the next reach downstream
        to_reach <- filter(reaches_bounded, seg_id == reach$to_seg)
        which.max(st_distance(reach$end_points[[1]], st_geometry(to_reach)))
      } else { # river mouths
        # upstream endpoint is the one closest to a previous reach upstream
        from_reach <- filter(reaches_bounded, to_seg == reach$seg_id)[1,] # only need one
        which.min(st_distance(reach$end_points[[1]], st_geometry(from_reach)))
      } # turns out that the first point is always the upstream point for the DRB at least. that's nice, but we won't rely on it.
    }),
    up_point = purrr::map2(end_points, which_end_up, function(endpoints, whichendup) { endpoints[[whichendup]] }),
    down_point = purrr::map2(end_points, which_end_up, function(endpoints, whichendup) { endpoints[[c(1,2)[-whichendup]]] })) %>%
  select(-end_points)

# split reaches into segments where necessary so that every pair of points is
# connected by a segment, and no segments pass through and beyond a point
reach_net_info <- lapply(unique(reaches_bounded$seg_id), function(segid) {
    reach <- filter(reaches_bounded, seg_id == segid)
    from_reaches <- filter(reaches_bounded, to_seg == reach$seg_id)
    
    reach_points <- tibble(
      point_raw = c(reach$up_point, from_reaches$down_point, reach$down_point),
      pt_seg = c(segid, from_reaches$seg_id, segid),
      pt_seg_flowend = c('up', rep('down', length(point_raw) - 1)), # 'up' and 'down' are positions relative to flow
      point_snapped = lapply(point_raw, sf::st_snap, reach, tol=1e-9)) %>%
      st_set_geometry('point_snapped') %>%
      st_set_crs(st_crs(reach)) %>%
      mutate(point_id = sprintf('%d%s', pt_seg, substring(pt_seg_flowend, 1, 1)))
    
    # split the reach into segments between neighboring pairs of points
    net_geoms <- lwgeom::st_split(st_geometry(reach), st_geometry(reach_points)) %>%
      st_collection_extract("LINESTRING")
    net_edges <- tibble(
      subseg_seg = rep(segid, length(net_geoms))) %>%
      mutate(tmp_subseg_id = 1:n()) %>% # temporary id until we can rename by segment order
      st_set_geometry(net_geoms) %>%
      mutate(
        left = lapply(geometry, function(geom) {
          st_line_sample(geom, sample=0) %>% st_cast('POINT')
        }),
        right = lapply(geometry, function(geom) {
          st_line_sample(geom, sample=1) %>% st_cast('POINT')
        })) # 'left' and 'right' are positions within the linestring
    net_vertices <- net_edges %>% # one row per edge end (non-unique vertex) but geometry is still the edge
      gather('subseg_lineend', 'point', left, right)
    
    # navigate the set of subsegments by point to figure out their order
    start_point <- filter(reach_points, pt_seg_flowend == 'up')
    # vertex_order will contain indices into net_vertices, ordered from upstream to downstream
    net_vertices_sf <- st_set_crs(do.call(c, net_vertices$point), st_crs(reach)) # format vertex points for computing distances to reach points
    vertex_order <- which.min(st_distance(start_point, net_vertices_sf)) # first element of vertex_order
    while(length(vertex_order) < nrow(net_vertices)) {
      edge_start <- net_vertices[tail(vertex_order,1),]
      edge_end_index <- net_vertices %>%
        mutate(is_edge_end = (tmp_subseg_id == edge_start$tmp_subseg_id & subseg_lineend != edge_start$subseg_lineend)) %>%
        pull(is_edge_end) %>% which()
      vertex_order <- append(vertex_order, edge_end_index)
      if(length(vertex_order) < nrow(net_vertices)) {
        edge_end <- net_vertices[edge_end_index,]
        next_edge_index <- net_vertices %>%
          mutate(is_next_edge = (sapply(point, function(pt) { pt == edge_end$point[[1]] }) & tmp_subseg_id != edge_end$tmp_subseg_id)) %>%
          pull(is_next_edge) %>% which()
        vertex_order <- append(vertex_order, next_edge_index)
      }
    }
    # create a flat table of edges (including the geometry) and their two
    # vertices apiece (one vertex per row), ordered by flow from upstream to
    # downstream
    net_vertices_ordered <- net_vertices[vertex_order, ] %>%
      mutate(
        vertex_order = 1:n(),
        subseg_id = sprintf('%d_%d', subseg_seg, ceiling((1:n())/2))) %>%
      group_by(subseg_id) %>%
      arrange(vertex_order) %>%
      mutate(
        subseg_updown=paste(substring(subseg_lineend,1,1), collapse=''), # lr means left side of linestring is upstream, flows to right side = downstream
        pt_subseg_flowend=ifelse(substring(subseg_updown,1,1) == substring(subseg_lineend,1,1), 'up', 'down')) %>% 
      ungroup() %>%
      select(subseg_id, subseg_seg, subseg_updown, geometry, pt_subseg_flowend, point)
    
    # join with info about how these subsegment endpoints relate to the original network
    net_vertices_ord_sf <- st_set_crs(do.call(c, net_vertices_ordered$point), st_crs(reach)) # format REORDERED vertex points for computing distances to reach points
    net_dists <- st_distance(net_vertices_ord_sf, reach_points) # row = net_vertex_ordered, col=reach_point
    net_info_flat <- net_vertices_ordered %>% 
      mutate(reach_point_matches = lapply(seq_len(nrow(net_dists)), function(i) {
        which(net_dists[i,] < units::set_units(1, 'm')) # find essentially equal points. should be completely equal (dist = 0, but be tolerant)
      })) %>%
      unnest(reach_point_matches) %>%
      mutate(point_id = reach_points$point_id[reach_point_matches]) %>%
      left_join(st_drop_geometry(reach_points), by='point_id') %>%
      select(-reach_point_matches)

    # pull out a table of unique vertices, weeding out duplication in (1)
    # original reach endpoints or contributing reach end points and (2) the end
    # and start of subsequent line subsegments generated in this loop
    net_vertices_final <- net_info_flat %>%
      group_by(subseg_id, pt_subseg_flowend) %>%
      summarize(
        point_geom = unique(point_raw),
        point_ids = paste(sort(point_id), collapse=';')) %>%
      ungroup() %>%
      group_by(point_ids) %>%
      st_drop_geometry() %>%
      summarize(
        point_geom = unique(point_geom),
        starts_subseg = if(any(pt_subseg_flowend=='up')) subseg_id[pt_subseg_flowend=='up'] else as.character(NA),
        ends_subseg = if(any(pt_subseg_flowend=='down')) subseg_id[pt_subseg_flowend=='down'] else as.character(NA)) %>%
      { st_set_geometry(., st_set_crs(do.call(st_sfc, .$point_geom), st_crs(reach))) } %>%
      select(-point_geom)
    
    # pull out a table of unique edges (subsegments)
    net_edges_dups <- net_info_flat %>%
      group_by(subseg_id) %>%
      mutate(
        starts_with_pts = paste(sort(point_id[pt_subseg_flowend=='up']), collapse=';'),
        ends_with_pts = paste(sort(point_id[pt_subseg_flowend=='down']), collapse=';'),
        from_segs = paste(sort(pt_seg[pt_seg != segid & pt_subseg_flowend=='up']), collapse=';')
      ) %>%
      ungroup()
    first_dups <- sapply(unique(net_edges_dups$subseg_id), function(ssid) which(net_edges_dups$subseg_id == ssid)[1])
    net_edges_final <- net_edges_dups %>%
      slice(first_dups) %>%
      mutate(
        subseg_id_num = as.integer(sapply(strsplit(subseg_id, '_'), `[`, 2)),
        to_seg = ifelse(subseg_id_num == max(subseg_id_num), reach$to_seg, NA),
        to_subseg = c(subseg_id[-1], ''),
        subseg_length = st_length(geometry)) %>%
      select(subseg_id, subseg_seg, subseg_updown, subseg_length, starts_with_pts, ends_with_pts, from_segs, to_seg, to_subseg, geometry)
  
    return(list(vertices=net_vertices_final, edges=net_edges_final))  
  })

reach_net_vertices <- do.call(rbind, lapply(reach_net_info, `[[`, 'vertices'))
reach_net_edges <- do.call(rbind, lapply(reach_net_info, `[[`, 'edges'))
for(e in reach_net_edges$subseg_id) {
  reach <- filter(reach_net_edges, subseg_id == e)
  if(reach$to_subseg == '') {
    if(!is.na(reach$to_seg)) {
      to_subseg <- filter(
        reach_net_edges,
        reach$to_seg == subseg_seg,
        grepl(sprintf('(^|;)%d(;|$)',reach$subseg_seg), from_segs))
      if(nrow(to_subseg) != 1) {
        print(to_subseg)
        stop(sprintf("found no to_subseg for %s, to_seg=%d", e, reach$to_seg))
      }
      to_subseg_id <- to_subseg$subseg_id
    } else {
      to_subseg_id <- NA
    }
    reach_net_edges[which(reach_net_edges$subseg_id == e), 'to_subseg'] <- to_subseg_id
  }
}








# Augment reach info with connections to the upstream reaches. In each row of
# reaches_updown (an edge in the graph), the outlet (pour point) of from_seg is the upstream point, seg_id is used as
# the ID of the segment's downstream-most point, and the segment
# identified by seg_id (between from_seg and seg_id) has length seg_length.
reach_key_info <- reaches_bounded %>% select(seg_id, to_seg, which_end_up, up_point, down_point)
reaches_updown <- reach_key_info %>%
  rename(next_seg = to_seg) %>% # to_seg will just be confusing in a minute, but we need it to decide whether something is an outlet, so rename it
  #mutate(to_seg = seg_id) %>% # add a column to clarify that the seg_id and the to_pt share an ID
  left_join(
    reach_key_info %>% st_set_geometry(NULL) %>% select(from_seg = seg_id, seg_id = to_seg),
    by = 'seg_id')

# plot the original reaches by their position
reaches_updown %>%
  mutate(position = ifelse(is.na(next_seg), 'outlet', ifelse(is.na(from_seg), 'headwater', 'middle'))) %>%
  ggplot() + geom_sf(aes(color=position), fill=NA) + theme_bw()


# identify the points at the upstream and downstream ends of each reach
reach_endpoints <- lapply(seq_len(nrow(gf_reaches)-454), function(i) {
  reach <- gf_reaches[i,]
  end_points <- st_cast(sf::st_geometry(reach), "POINT")
  head_point <- head(end_points, 1)
  tail_point <- tail(end_points, 1)

  out_point <- filter(gf_points, poi_gage_segment == reach$seg_id)
  out_snap <- sf::st_snap(out_point, reach, tol=1e-9)
  sf::st_geometry(out_snap) == tail_point
  
  %>% {.[1,length(.)]}
  parts <- st_collection_extract(lwgeom::st_split(sf::st_geometry(reach_3), sf::st_geometry(site1_snap3)),"LINESTRING")
  in_point <- 
  print(reach)
})



# Simplify the vertex information so it makes sense relative to the edge information (reaches_updown)
points_updown <- gf_points %>%
  select(pt_id = poi_gage_segment)
ggplot(reaches_updown) + geom_sf(aes(color=seg_id), size=0.5) + geom_sf(data=points_updown, aes(color=pt_id), size=1) + theme_bw()

# Use igraph to compute a matrix of distances (in meters)
graph <- reaches_updown %>%
  st_set_geometry(NULL) %>%
  select(from_pt, to_pt, seg_id, seg_length) %>%
  igraph::graph_from_data_frame(directed = TRUE) # df must contain "a symbolic edge list in the first two columns. Additional columns are considered as edge attributes"
dists_symmetric <- distances(graph, weights = edge.attributes(graph)$seg_length, mode='all') # symmetric
dists_downstream <- distances(graph, weights = edge.attributes(graph)$seg_length, mode='out') # shortest paths FROM each vertex.
dists_upstream <- distances(graph, weights = edge.attributes(graph)$seg_length, mode='in') # shortest paths TO each vertex (flowing upstream only)

# example: from_pt=2 flows to to_pt=8 with reach length = pull(filter(reaches_updown, seg_id=='8'), seg_length) = 19728 m
# dists_downstream['3','4'] = 1914
# dists_downstream['4','3'] = Inf
# dists_upstream['3','4'] = Inf
# dists_upstream['4','3'] = 1914
# dists_symmetric['3','4'] = dists_symmetric['4','3'] = 1914


ggplot(mutate(reaches_updown, example=seg_id %in% c('8'))) + geom_sf(aes(color=example), fill=NA) + theme_bw()
ggplot(mutate(reaches_updown, example=seg_id %in% c('2'))) + geom_sf(aes(color=example), fill=NA) + theme_bw()
ggplot(mutate(filter(reaches_updown, seg_id %in% c('2','3','8')), example=seg_id %in% c('3'))) + geom_sf(aes(color=example), fill=NA) + theme_bw()

# example: 
reaches_updown %>% filter(position=='headwater')
ggplot(mutate(reaches_updown, example=seg_id %in% c('1','7'))) + geom_sf(aes(color=example), fill=NA) + theme_bw()


site_1 <- filter(points_updown, pt_id == 1)
reach_3 <- filter(gf_reaches, seg_id == 3)
site1_snap3 <- st_snap(site_1, reach_3, tol=1e-9)
sf::st_geometry(reach_3)
parts <- st_collection_extract(lwgeom::st_split(sf::st_geometry(reach_3), sf::st_geometry(site1_snap3)),"LINESTRING")

site_2 <- filter(points_updown, pt_id == 2)
reach_8 <- filter(reaches_updown, seg_id == 8)
site_snap <- st_snap(site_2, reach_8, tol=1e-9)
parts <- st_collection_extract(lwgeom::st_split(sf::st_geometry(reach_8), sf::st_geometry(site_snap)),"LINESTRING")

dists_downstream %>%
  as_tibble() %>%
  rownames_to_column('')
mutate(reaches_updown, )       
